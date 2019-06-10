!===============================================================================
!> **ZEUS-TW MHD method/operators module**
!> The type defined/extended here is intended for instantiation/invokation/use.
!> This particular form uses rho, e, \vec{s} = rho * \vec{v}, \vec{B} as the
!> MHD fluid variables.
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE extras_mod
  USE kinds_mod
  USE variable_mod
  USE mesh_mod
  USE timestep_mod
  USE index_mod
  USE force_mod
  USE zeus_variables, only: dt, ix1x2x3, dtold
  implicit none
  private

  integer, parameter:: KindReal = 8

  ! a derived-type to keep track of the ZEUS boundary types; ZEUS currently
  ! handles setting/resetting of boundaries.
  type:: zeus_boundaries_t
    integer:: iib=6, oib=6, ijb=6, ojb=6, ikb=6, okb=6
  contains
  end type

  type, extends(extras_t), public:: mhd_t
    type(scalar_t):: d, e, gp
    type(vector_t):: s, B
    logical:: mhd = .true.
    type(zeus_boundaries_t):: btypes
    type(force_t):: force
    integer:: ix1x2x3
    logical:: initialized = .false.
  contains
    procedure:: init => init_zeusTW
    procedure:: update => update_zeusTW
    procedure:: gas_pressure
    procedure:: MapVectors
    procedure:: courant_condition => do_nothing ! this is handled internally by the ZEUS solver
  end type

PUBLIC CopyToZeusTW, CopyFromZeusTW
CONTAINS

!===============================================================================
!> Initialise ZEUS-TW module.
!===============================================================================
SUBROUTINE init_zeusTW (self)
  use zeus_variables_shared
  use zeus_initial

  class(mhd_t):: self
  integer:: i, iostat
  logical, save :: first_time=.true.
  real, save:: courant=0.5
  namelist / zeustw_params /  gamma, qcon, qlin, courant, ienergy, dtratini, &
                              iohmic, with_ad_1f, ihall, &
                              ohmic_eta_factor, ohmic_exponent, ohmic_rho_transition , &
                              max_ohmic_steps, courno_ohmic, &
                              with_ad_1f_vion_force, courno_ad_1f, ad_dt_floor, &
                              ad_rhoi0, ad_rhon0, ad_alpha, ad_gamma, &
                              hall_coefficient_0, hall_rhoreference, hall_alpha, &
                              hall_number_of_cycles, hall_max_number_of_cycles, &
                              courno_hall, hall_dt_floor, &
                              withchem, tabulate, gbins, numIons, nbins, tablename
  !-----------------------------------------------------------------------------
  if (first_time) then
    call trace_begin ('mhd_mod%init_zeusTW')

    first_time = .false.

    ! default values
    gamma   = 5.0 / 3.0
    qcon    = 2.0
    qlin    = 0.0
    ienergy = 1
    courant = 0.5

    rewind (io%input); read (io%input, zeustw_params, iostat=iostat)
    if (iostat/=0) call io%namelist_warning ('zeustw_params')
    if (io%master) write (*, zeustw_params)

    ! Initialise ZeusTW
    call initialize_zeustw
    ! Initialise DISPATCH time-stepping.
    call timestep%init
  endif

  self%nv = 12
  self%gamma = gamma
  self%courant = courant
  self%ix1x2x3 = 1

  self%kind = 'zeus_tw_mhd_patch'
  ! don't store the time derivatives, only the MHD quantities;
  ! this also implies that `time_order` in `timestep_params` should be set to 0.
  self%nw = 1
  self%idx%d=1   ; self%idx%e=2
  self%idx%px=3  ; self%idx%py=4; self%idx%pz=5
  self%idx%bx=6  ; self%idx%by=7; self%idx%bz=8
  self%idx%tt=-1 ; self%idx%s=-1; ; self%idx%et=-1; self%idx%dphi=-1
  self%idx%ex=10 ; self%idx%ey=11; self%idx%ez=12
  io%nv = 10 ! do not write EMFs to output
  call self%gpatch_t%init()
  self%idx%phi=9 ! ZEUS always has a slot for the grav. pot. included (even if not used).
  do i=1,3
    self%mesh(i)%h(self%idx%px+i-1) = -0.5 ! momentum/velocity components
    self%mesh(i)%h(self%idx%bx+i-1) = -0.5 ! mag. field components
    self%mesh(i)%h(self%idx%ex:self%idx%ez) = -0.5
    self%mesh(i)%h(self%idx%ex+i-1) = 0.0 ! EMFs; each component is cell-centred in its own direction.
  end do
  idx = self%idx

  if (io%restart < 0) then
    self%t(self%it) = 0.0d0
    self%ix1x2x3 = 1
  end if

  call trace_end
END SUBROUTINE init_zeusTW

!===============================================================================
!===============================================================================
SUBROUTINE update_zeusTW (self)
  use zeus_main

  class(mhd_t) :: self
  integer, save:: itimer=0
  real:: error
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  call CopyToZeusTW (self)
  call trace_begin ('zeusTW (external)', itimer=itimer)
  call ZeusTW (self%istep)  ! run ZEUS
  call trace_end(itimer)
  call CopyFromZeusTW (self)
  !-----------------------------------------------------------------------------
  call self%counter_update
END SUBROUTINE update_zeusTW

!===============================================================================
!===============================================================================
SUBROUTINE CopyToZeusTW (self)
  USE zeus_variables_shared
  USE zeus_grid
  USE zeus_main

  class(mhd_t) :: self
  real(KindReal) :: x1low, x2low, x3low, dx1, dx2, dx3, pvelocity(3)
  logical, save :: first_time=.true.
  !$omp threadprivate (first_time)
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: qz, force, p
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, Ux, Uy, Uz
  integer, save:: itimer1=0
  type(index_t):: idx

  !-----------------------------------------------------------------------------
  ! Remember: ZEUS has one extra zone in each direction to handle staggered variables
  ! (in ZEUS-speak this is "*ep3"); this is one more zone than STAGGER.
  ! Remember: ZEUS keeps track of *both* the internal and total energy (even if it
  ! only uses one at a time).
  qz => self%mem(:,:,:,:,self%it,1)

  ! should be in the native coord. system of the patch
  x1lower = self%llc_nat(1)
  x2lower = self%llc_nat(2)
  x3lower = self%llc_nat(3)
  x1upper = self%llc_nat(1) + self%size(1)
  x2upper = self%llc_nat(2) + self%size(2)
  x3upper = self%llc_nat(3) + self%size(3)

  zniib = self%btypes%iib
  znijb = self%btypes%ijb
  znikb = self%btypes%ikb
  znoib = self%btypes%oib
  znojb = self%btypes%ojb
  znokb = self%btypes%okb

  dx1 = self%mesh(1)%d
  dx2 = self%mesh(2)%d
  dx3 = self%mesh(3)%d
  x1low = self%mesh(1)%r(self%mesh(1)%li) + self%mesh(1)%centre_nat - 0.5d0 * dx1
  x2low = self%mesh(2)%r(self%mesh(2)%li) + self%mesh(2)%centre_nat - 0.5d0 * dx2
  x3low = self%mesh(3)%r(self%mesh(3)%li) + self%mesh(3)%centre_nat - 0.5d0 * dx3

  if (self%mesh(1)%type == 'Cartesian_mesh') then
    geometry = 'xyz'
  else if (self%mesh(1)%type == 'cylindrical_mesh') then
    geometry = 'zrp'
  else if (self%mesh(1)%type == 'spherical_mesh') then
    geometry = 'rtp'
  else
    call mpi%abort("ERROR: Mesh type unknown in CopyToZeus! Abort!")
  end if

  ! assign the patch velocity in its native coord. system; currently
  ! only possible in spherical and cylindrical coords., and only for the
  ! phi-(i.e.3-)direction.
  pvelocity(1) = self%velocity(1)
  pvelocity(2) = self%velocity(2)
  pvelocity(3) = self%velocity(3) + self%grid_vel(3)
  
  if (associated(self%force%selected)) then
    d => qz(:,:,:,self%idx%d)
    p => qz(:,:,:,self%idx%px:self%idx%pz)
    allocate (force(self%gn(1),self%gn(2),self%gn(3),3))
    force = self%force%selected(self%time, d, p, Ux, Uy, Uz, self%mesh)
  else
    force => null()
  end if

  call trace%begin ('CopyToZTW (external)', itimer=itimer1)
  call copytoztw (qz, self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, &
                 self%nv, x1low, x2low, x3low, dx1, dx2, dx3, self%id, &
                 pvelocity, force)
  call trace%end (itimer1)

  if (self%do_pebbles) then
    associate (v1=>self%vgas(:,:,:,1,self%it),v2=>self%vgas(:,:,:,2,self%it),v3=>self%vgas(:,:,:,3,self%it))
    call get_vel (self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, v1, v2, v3)
    end associate
    self%pressure(:,:,:,self%it) = self%mem(:,:,:,self%idx%e,self%it,1) * (self%gamma - 1.0)
  end if

  courno = self%courant
  ! initialise the time step for this patch; do only once for each patch.
  if (.not. self%initialized) then
    call initialize_zeus_dt
    self%dt(self%it) = dt
    self%initialized = .true.
  end if

  dt = self%dt(self%it)
  time = self%time
  ix1x2x3 = self%ix1x2x3
  dtold = dt

  ! copy physical boundary logical flags
  lx1ipb = self%mesh(1)%lower_boundary
  lx1opb = self%mesh(1)%upper_boundary
  lx2ipb = self%mesh(2)%lower_boundary
  lx2opb = self%mesh(2)%upper_boundary
  lx3ipb = self%mesh(3)%lower_boundary
  lx3opb = self%mesh(3)%upper_boundary

  if (associated(self%force%selected)) deallocate (force)
  nullify (qz)

END SUBROUTINE CopyToZeusTW

!===============================================================================
SUBROUTINE CopyFromZeusTW (self)
  use zeus_main
  use zeus_fluxes

  class(mhd_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: qz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  qz => self%mem(:,:,:,:,self%new,1)

  call trace_begin ('CopyFromZTW (external)', itimer=itimer)

  ! Save the EMFs as calculated by the solver.
  ! Note that what is actually stored is `E_i * g_i * dx_i`, where `E_i` are the
  ! actual EMFs (= v x B).
  self%mem(:,:,:,self%idx%ex,self%new,1) = amemf1(:,:,:)
  self%mem(:,:,:,self%idx%ey,self%new,1) = amemf2(:,:,:)
  self%mem(:,:,:,self%idx%ez,self%new,1) = amemf3(:,:,:)

  call copyfromztw (qz, self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, self%nv)
  call trace_end (itimer)

  ! update time step.
  !$omp atomic write
  self%dtime = dt
  self%ix1x2x3 = ix1x2x3

  nullify (qz)

END SUBROUTINE CopyFromZeusTW

!===============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%mem(:,:,:,self%idx%d,self%it,1)
  else
    pg = self%mem(:,:,:,self%idx%e,self%it,1) * (self%gamma - 1.0)
  end if
END FUNCTION gas_pressure

!===============================================================================
!> In `download_mod::interpolate_to` at least, this is invoked by the "source"
!> patch (`self`) and applied to `target`.
!> In the case of momentum components in spherical and cylindrical coords., it
!> is assumed that it is the *angular momentum* that is stored within `in`.
!>
!> target -> use `ii`
!> source/self -> use `jj`
!===============================================================================
SUBROUTINE MapVectors (self, in, target, ii, jj, out)
  USE patch_mod, only: patch_t
  class(mhd_t):: self
  class(patch_t), pointer:: target
  real(kind=KindScalarVar), dimension(:), intent(in):: in
  integer, intent(in):: ii(3), jj(3)
  integer:: sidxs(3), Bidxs(3)
  real(kind=KindScalarVar),intent(out):: out(size(in))
  real(kind=KindReal):: tphi, ttheta, tr, sr, stheta, sphi, sx, sy, sz, dvx, dvy, dvz, Bx, By, Bz
  real(kind=KindScalarVar), dimension(3):: BCart, dvCart
  real(kind=KindReal), dimension(3):: txi
  !.............................................................................
  call trace_begin('mhd_t%MapVectors')

  sidxs = [self%idx%px, self%idx%py, self%idx%pz]
  Bidxs = [self%idx%bx, self%idx%by, self%idx%bz]
  select case (self%mesh_type)
  case (1) ! source: Cartesian
    sx = self%mesh(1)%r(jj(1)) + self%mesh(1)%centre_nat
    sy = self%mesh(2)%r(jj(2)) + self%mesh(2)%centre_nat
    sz = self%mesh(3)%r(jj(3)) + self%mesh(3)%centre_nat
    if (target%mesh_type == mesh_types%Cartesian) then
      ! target: Cartesian
      ! **no conversion necessary**
      out = in
    else if (target%mesh_type == mesh_types%spherical) then
      ! target: spherical
      out(sidxs(1:3)) = MomentaCartToSph(in(sidxs(1:3)), sx, sy, sz)
      out(Bidxs(1:3)) = VectorCartToSph(in(Bidxs(1:3)), sx, sy, sz)
    else if (target%mesh_type == mesh_types%cylindrical) then
      ! target: cylindrical
      out(sidxs(1:3)) = MomentaCartToCyl(in(sidxs(1:3)), sx, sy, sz)
      out(Bidxs(1:3)) = VectorCartToCyl(in(Bidxs(1:3)), sx, sy, sz)
    end if

  case (2) ! spherical
    sr     = self%mesh(1)%r(jj(1)) + self%mesh(1)%centre_nat
    stheta = self%mesh(2)%r(jj(2)) + self%mesh(2)%centre_nat
    sphi   = self%mesh(3)%r(jj(3)) + self%mesh(3)%centre_nat
    if (target%mesh_type == mesh_types%Cartesian) then
      ! target: Cartesian
      out(sidxs(1:3)) = AnyVectorSphToCart(in(sidxs(1:3)), sr, stheta, sphi)
      out(Bidxs(1:3)) = AnyVectorSphToCart(in(Bidxs(1:3)), sr, stheta, sphi)
    else if (target%mesh_type == mesh_types%spherical) then
      ! target: spherical
      ! note that the conversion goes through Cartesian coords. to potentially
      ! account for different origins w.r.t. "target" and "source".
      txi = SphericalToCartesian(sr, stheta, sphi)
      dvCart(1:3) = AnyVectorSphToCart(in(sidxs(1:3)), sr, stheta, sphi)
      BCart(1:3)  = AnyVectorSphToCart(in(Bidxs(1:3)), sr, stheta, sphi)
      out(sidxs(1:3)) = MomentaCartToSph(dvCart, txi(1), txi(2), txi(3))
      out(Bidxs(1:3)) = VectorCartToSph(BCart, txi(1), txi(2), txi(3))
    else if (target%mesh_type == mesh_types%cylindrical) then
      ! target: cylindrical
      ! note that the conversion goes through Cartesian coords. to potentially
      ! account for different origins w.r.t. "target" and "source".
      txi = SphericalToCartesian(sr, stheta, sphi)
      dvCart(1:3) = AnyVectorSphToCart(in(sidxs(1:3)), sr, stheta, sphi)
      BCart(1:3)  = AnyVectorSphToCart(in(Bidxs(1:3)), sr, stheta, sphi)
      out(sidxs(1:3)) = MomentaCartToCyl(dvCart, txi(1), txi(2), txi(3))
      out(Bidxs(1:3)) = VectorCartToCyl(BCart, txi(1), txi(2), txi(3))
    end if

  case (3) ! source: cylindrical
    sz   = self%mesh(1)%r(jj(1)) + self%mesh(1)%centre_nat
    sr   = self%mesh(2)%r(jj(2)) + self%mesh(2)%centre_nat
    sphi = self%mesh(3)%r(jj(3)) + self%mesh(3)%centre_nat
    if (target%mesh_type == mesh_types%Cartesian) then
      ! target: Cartesian
      out(sidxs(1:3)) = AnyVectorCylToCart(in(sidxs(1:3)), sz, sr, sphi)
      out(Bidxs(1:3)) = AnyVectorCylToCart(in(Bidxs(1:3)), sz, sr, sphi)
    else if (target%mesh_type == mesh_types%spherical) then
      ! target: spherical
      ! note that the conversion goes through Cartesian coords. to potentially
      ! account for different origins w.r.t. "target" and "source".
      txi = CylindricalToCartesian(sz, sr, sphi)
      dvCart(1:3) = AnyVectorCylToCart(in(sidxs(1:3)), sz, sr, sphi)
      BCart(1:3)  = AnyVectorCylToCart(in(Bidxs(1:3)), sz, sr, sphi)
      out(sidxs(1:3)) = MomentaCartToSph(dvCart, txi(1), txi(2), txi(3))
      out(Bidxs(1:3)) = VectorCartToSph(BCart, txi(1), txi(2), txi(3))
    else if (target%mesh_type == mesh_types%cylindrical) then
      ! target: cylindrical
      ! note that the conversion goes through Cartesian coords. to potentially
      ! account for different origins w.r.t. "target" and "source".
      txi = CylindricalToCartesian(sz, sr, sphi)
      dvCart(1:3) = AnyVectorCylToCart(in(sidxs(1:3)), sz, sr, sphi)
      BCart(1:3)  = AnyVectorCylToCart(in(Bidxs(1:3)), sz, sr, sphi)
      out(sidxs(1:3)) = MomentaCartToCyl(dvCart, txi(1), txi(2), txi(3))
      out(Bidxs(1:3)) = VectorCartToCyl(BCart, txi(1), txi(2), txi(3))
    end if
  end select
  out(1:sidxs(1)-1) = in(1:sidxs(1)-1)

  call trace_end
CONTAINS

  !=============================================================================
  !> Expected format for `input`: A_x, A_y, A_z
  !=============================================================================
  FUNCTION VectorCartToSph (input, x, y, z) RESULT (vec)
    real(kind=KindReal), intent(in):: x, y, z
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = (x * input(1) + y * input(2) + z * input(3)) &
           / sqrt(x**2 + y**2 + z**2)
    vec(2) = (input(3) * (x**2 + y**2) &
           - z * (x * input(1) + y * input(2))) / sqrt(x**2 + y**2) &
           / sqrt(x**2 + y**2 + z**2)
    vec(3) = -(x * input(2) - y * input(1)) / sqrt(x**2 + y**2)
  END FUNCTION VectorCartToSph

  FUNCTION VectorCartToCyl (input, x, y, z) RESULT(vec)
    real(kind=KindReal), intent(in):: x, y, z
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = input(3)
    vec(2) = (x * input(1) + y * input(2)) / sqrt(x**2 + y**2)
    vec(3) = (x * input(2) - y * input(1)) / sqrt(x**2 + y**2)
  END FUNCTION VectorCartToCyl

  FUNCTION MomentaCartToCyl (input, x, y, z) RESULT (vec)
    real(kind=KindReal), intent(in):: x, y, z
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = input(3)
    vec(2) = target%mesh(1)%h2c(ii(1)) * (x * input(1) + y * input(2)) &
           / sqrt(x**2 + y**2)
    vec(3) = target%mesh(1)%h31c(ii(1)) * target%mesh(2)%h32c(ii(2)) &
           * (x * input(2) - y * input(1)) / (x**2 + y**2)
  END FUNCTION MomentaCartToCyl

  FUNCTION MomentaCartToSph (input, x, y, z) RESULT (vec)
    real(kind=KindReal), intent(in):: x, y, z
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = (x * input(1) + y * input(2) + z * input(3)) / sqrt(x**2 + y**2 + z**2)
    vec(2) = target%mesh(1)%h2c(ii(1)) &
           * (input(3) * (x**2 + y**2) - z * (x * input(1) + y * input(2))) &
           / sqrt(x**2 + y**2) / (x**2 + y**2 + z**2)
    vec(3) = -target%mesh(1)%h31c(ii(1)) * target%mesh(2)%h32c(ii(2)) &
           * (x * input(2) - y * input(1)) / (x**2 + y**2)
  END FUNCTION MomentaCartToSph

  FUNCTION AnyVectorCylToCart (input, Z, R, Phi) RESULT (vec)
    real(kind=KindReal), intent(in):: Z, R, Phi
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = input(2) * cos(Phi) - input(3) * sin(Phi)
    vec(2) = input(2) * sin(Phi) + input(3) * cos(Phi)
    vec(3) = input(1)
  END FUNCTION AnyVectorCylToCart

  FUNCTION AnyVectorSphToCart (input, r, theta, phi) RESULT (vec)
    real(kind=KindReal), intent(in):: r, theta, phi
    real(kind=KindScalarVar), dimension(3):: input, vec
    !...........................................................................
    vec(1) = input(1) * sin(theta) * cos(phi) &
           + input(2) * cos(theta) * cos(phi) &
           - input(3) * sin(phi)
    vec(2) = -input(1) * sin(theta) * sin(phi) &
           - input(2) * cos(theta) * sin(phi) &
           - input(3) * cos(phi)
    vec(3) = -input(1) * cos(theta) &
           + input(2) * sin(theta)
  END FUNCTION AnyVectorSphToCart

END SUBROUTINE MapVectors

!===============================================================================
!> This subroutine intentionally DOES NOTHING.
!===============================================================================
SUBROUTINE do_nothing (self, detailed_timer)
  class(mhd_t):: self
  logical, optional:: detailed_timer
  !.............................................................................
END SUBROUTINE do_nothing

END MODULE
