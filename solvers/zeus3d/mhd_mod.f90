!===============================================================================
!> **ZEUS-3D MHD method/operators module**
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
    type(scalar_t):: d, e, et, gp
    type(vector_t):: s, B
    logical:: mhd = .true.
    type(zeus_boundaries_t):: btypes
    type(force_t):: force
    integer:: ix1x2x3
  contains
    procedure:: init => init_zeus
    procedure:: update => update_zeus
    procedure:: gas_pressure
    procedure:: MapVectors
    procedure:: courant_condition => do_nothing ! this is handled internally by the ZEUS solver
  end type

PUBLIC CopyToZeus, CopyFromZeus
CONTAINS

!===============================================================================
!> Initialise ZEUS module.
!===============================================================================
!---------------------------  NAMELIST PARAMETERS  ---------------------------
!
!  qcon      quadratic artificial viscosity (q) constant
!  qlin      linear    artificial viscosity (q) constant
!  gamma     ratio of specific heats
!  iord      order of interpolation for all variables (default 2)
!            = 1 => first order (donor cell) interpolation
!            =-1 => first order implicit interpolation (compressive s*
!                   only)
!            = 2 => second order (van Leer) interpolation
!            =-2 => second order implicit interpolation (compressive s*
!                   only)
!            = 3 => piecewise parabolic interpolation (PPI)
!            = 4 => second order velocity-corrected interpolation
!            = 5 => second order non-harmonic interpolation
!            = 6 => second order interpolation for non-uniform grid
!  iords     order of interpolation for scalars to override "iord"
!            (default iord)
!  istp      discontinuity steepener (PPI)
!            =0 => do not steepen any variables (default)
!            =1 => steepen all variables
!            =2 => steepen only d and only at contact discontinuities
!  itote     =0 => solve the internal energy equation (positive definite
!                  pressures but energy not conserved).
!            =1 => solve the total energy equation (energy conserved
!                  but pressure not positive definite) (default)
!            See notes below on internal vs. total energy equation.
!  iscydf    =0 => no subcycling on diffusion (default)
!            =1 => subcycle on diffusion
!  iscyqq    =0 => no subcycling on artificial viscosity (default)
!            =1 => subcycle on artificial viscosity
!  ix1x2x3   seed for directional splitting sequence (default=1)
!            no symmetry: starting at ix1x2x3=1 = starting at ix1x2x3=3
!                         (5) with the x1-x2-x3 axes interchanged for
!                         x2-x3-x1 (x3-x1-x2).
!            for ISYM   : starting at ix1x2x3=1 = starting at ix1x2x3=3
!                         with the x1-x2 axes interchanged for x2-x1.
!            for JSYM   : starting at ix1x2x3=1 = starting at ix1x2x3=4
!                         with the x1-x3 axes interchanged for x3-x1.
!            for KSYM   : starting at ix1x2x3=1 = starting at ix1x2x3=4
!                         with the x1-x2 axes interchanged for x2-x1.
!            with any two symmetry flags, it should make no difference
!            what value ix1x2x3 starts at.
!  nuscalar  kinematic viscosity (in units of LV) (default=0.0)
!  diffc1, diffc2  diffusion coefficient is (for the second fluid)
!                  set to diffc1 / B**diffc2
!  isetemf   affects flow-in skin values for emf(perp) and flow-in
!            boundary values for emf(par) for btype=10 only.
!            =0 => SVALEMF* and BVLAEMFS don't overwrite (C)MOC*-
!                  computed flow-in emfs with pre-set flow-in arrays
!                  (default).
!            =1 => SVALEMF* overwrite (C)MOC*-computed flow-in emfs with
!                  preset flow-in arrays (e.g., "emf2iib1", etc.), but
!                  BVALEMFS doesn't.
!            =2 => Both SVALEMF* and BVALEMFS overwrite (C)MOC*-computed
!                  flow-in emfs with preset flow-in arrays.
!  cadvs     =1 => consistent advection (CA) used with scalar transport
!            =0 => not (default)
!  cadvc     =1 => CA used with compressive momentum transport (default)
!            =0 => not
!  cadvt     =1 => CA used with transverse momentum transport (default)
!            =0 => not
!    CA applied to energy transport causes undue ringing in Fig 1b and
!    Fig 2b of Ryu & Jones (1995).  Since 12/07, no deleterious effect
!    of disabling CA from energy transport has been detected.  By
!    contrast, as of 5/13, no test has shown CA should be disabled from
!    momentum transport; indeed, some tests are marginally better with.
!  trnvrsn   = 0 => use legacy CMoC algorithm (MoC, HSMOC, original CT)
!            = 1 => use the FIT algorithm (default)
!            = 2 => same as 1, but with transverse Lorentz terms
!                   rendered in conservative form.
!            =-2 => same as 2 but for legacy transport
!            = 3 => same as 2, but with new source terms in CTRAN*
!                   expressed as a true flux triggering "magnetic
!                   centrifugal" terms in STV*
!            =-3 => same as 3 but for legacy transport
!
!  Toggles for highly experimental "downwinding" algorithm:
!
!  uordq     =1 =>   upwinded interpolations of qq (artificial visc.)
!            =2 => downwinded interpolations of qq
!  avordw    =1 => use 2-point averages for v* for et <-> e
!            =2 => downwind v* for et <-> e (default)
!  voutbc    =1 => outflow BC, v1(is) set to min(0,v1(isp1)), original
!                  Wilson scheme (similar for other boundaries; default)
!            =2 => outflow BC, v1(is) set to min(0,v1(is)) to reduce
!                  negative pressures at boundaries
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!              NOTES ON INTERNAL VS. TOTAL ENERGY EQUATION
!
!  Note that it never makes sense to solve the total energy equation
!  with ISO (isothermal equation of state) turned on, since in this
!  limit, the total energy equation is redundant to the momentum and
!  induction equations.  Thus, if ISO is defined and itote=1, "itote" is
!  reset to 0 and a warning message issued.
!
!  In earlier versions, TWOFLUID also caused "itote" to be reset to 0
!  As of August, 2006, the total energy equation has been recast to
!  include the second energy (and thus the diffusion term acts on "e" as
!  well as "e2" for itote=1), and this restriction has been lifted.
!
!  Because viscous heating is part of the total energy flux density,
!  solving the total energy equation precludes viscous subcycling.
!  If itote=1 and iscyqq=1, iscyqq is automatically reset to 0 and a
!  warning message is issued.
!
!  In keeping with the general rule that the user set the *primitive*
!  variables at initialisation, the user should set the first internal
!  energy density, "e", and any required inflow boundary values (e.g.,
!  "eiib1", etc.) and *not* the total energy density, "et", regardless
!  of whether the internal/total energy equation (itote=0/1) is being
!  solved.  A call to TOTNRG in SETUP after the problem generator has
!  been called sets the total energy array "et" and its inflow boundary
!  arrays to the appropriate values based on the user-settings of the
!  primitive variables.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!===============================================================================
SUBROUTINE init_zeus (self)
  USE zparams
  class(mhd_t):: self
  integer:: i, iostat
  logical, save :: first_time=.true.
  real, save:: courant=0.5
  real(8), save:: gamma
  namelist / zeus_params /  gamma   , qcon    , qlin    , nghst   , courant &
                          , iord    , iords   , istp    , itote   , iscydf  &
                          , iscyqq  , ix1x2x3 , nuscalar, isetemf , trnvrsn &
                          , cadvs   , cadvc   , cadvt   , uordq   , avordw  &
                          , voutbc  , diffc1  , diffc2  , dtinit
  !-----------------------------------------------------------------------------
  if (first_time) then
    call trace_begin ('mhd_mod%init_zeus')

    first_time = .false.
    ! default values for namelist parameters
    call zdefault

    gamma   = 5.0 / 3.0
    qcon    = 2.0
    qlin    = 0.0
    nghst   = 2 ! for the time being, ZEUS only plays nicely when it has 2 ghost cells; FIXME!
    iord    = 2
    iords   = 0
    istp    = 0
    itote   = 1
    iscydf  = 0
    iscyqq  = 0
    isetemf = 0
    ix1x2x3 = 1
    cadvs   = 0
    cadvc   = 1
    cadvt   = 1
    uordq   = 1
    avordw  = 2
    voutbc  = 1
    trnvrsn = 1
    nuscalar = 0.0d0
    diffc1  = 1.0d0
    diffc2  = 0.0d0
    dtinit  = huge(1.0d0)
    rewind (io%input); read (io%input, zeus_params, iostat=iostat)
    if (iostat/=0) call io%namelist_warning ('zeus_params')
    if (io%master) write (*, zeus_params)

    ! ---- Fool-proofing ----
    if (itote .ne. 0) itote = 1

    if (itote .eq. 1 .and. gamma .eq. 1.0) then
      write (*, 2000)
2000  format('ZEUS_MHD_MOD: |---WARNING---> ISO (isothermal EOS) is ' &
            ,'inconsistent with solving the',/ &
            ,'total energy equation; "itote" reset to zero!')
      itote = 0
    endif

    if ( (itote .eq. 1) .and. (iscyqq .gt. 0) ) then
      write (*, 2010)
2010  format('ZEUS_MHD_MOD: |---WARNING---> One cannot subcycle on the ' &
            ,'artificial viscous timestep',/ &
            ,'when solving the total energy equation; "iscyqq" reset to zero.')
      iscyqq = 0
    endif

    ! Order of upwinding stuff
    if ( (iord    .lt.-2) .or. (iord    .gt. 6) ) iord    = 2
    if   (iords   .eq. 0)                         iords   = iord
    if ( (iords   .lt.-2) .or. (iords   .gt. 6) ) iords   = 2
    if ( (istp    .lt. 0) .or. (istp    .gt. 2) ) istp    = 0
    istpcd = istp
    if   (istp    .eq. 2)                         istp    = 0

    ! EMF boundary conditions
    if ( (isetemf .lt. 0) .or. (isetemf .gt. 2) ) isetemf = 0

    ! Consistent advection (CA) and upwinding/downwinding
    if ( (cadvs   .lt. 0) .or. (cadvs   .gt. 1) ) cadvs   = 0
    if ( (cadvc   .lt. 0) .or. (cadvc   .gt. 1) ) cadvc   = 1
    if ( (cadvt   .lt. 0) .or. (cadvt   .gt. 1) ) cadvt   = 1
    if ( (uordq   .lt. 1) .or. (uordq   .gt. 2) ) uordq   = 1
    if ( (avordw  .lt. 1) .or. (avordw  .gt. 2) ) avordw  = 2
    if ( (voutbc  .lt. 1) .or. (voutbc  .gt. 2) ) voutbc  = 1

    ! Transport algorithm
    if ( (trnvrsn .lt.-3) .or. (trnvrsn .gt. 3) .or. (trnvrsn .eq.-1) ) &
      trnvrsn = 1

    ! Initialise DISPATCH time-stepping.
    call timestep%init
  endif

  self%ix1x2x3 = ix1x2x3
  self%nv = 13
  self%gamma = gamma
  gamma_z = gamma
  self%courant = courant

  self%kind = 'zeus_mhd_patch'
  self%etype = 'thermal'
  ! don't store the time derivatives, only the MHD quantities;
  ! this also implies that `time_order` in `timestep_params` should be set to 0.
  self%nw = 1
  self%idx%d=1   ; self%idx%e=2 ; self%idx%et=3
  self%idx%px=4  ; self%idx%py=5; self%idx%pz=6
  self%idx%bx=7  ; self%idx%by=8; self%idx%bz=9
  self%idx%tt=-1 ; self%idx%s=-1; self%idx%dphi=-1
  self%idx%ex=11 ; self%idx%ey=12; self%idx%ez=13
  io%nv = 10 ! do not write EMFs to output
  call self%patch_t%init()
  call self%gpatch_t%init()
  self%idx%phi=10 ! ZEUS always has a slot for the grav. pot. included (even if not used).
  do i=1,3
    self%mesh(i)%h(self%idx%px+i-1) = -0.5 ! momentum/velocity components
    self%mesh(i)%h(self%idx%bx+i-1) = -0.5 ! mag. field components
    self%mesh(i)%h(self%idx%ex:self%idx%ez) = -0.5
    self%mesh(i)%h(self%idx%ex+i-1) = 0.0 ! EMFs; each component is cell-centred in its own direction.
  end do
  idx = self%idx

  if (io%restart < 0) then
    self%dt(self%it) = dtinit
    self%t(self%it) = 0.0d0
  end if

  call trace_end
END SUBROUTINE init_zeus

!===============================================================================
!===============================================================================
SUBROUTINE update_zeus (self)
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
  call CopyToZeus (self)
  call trace_begin ('zeus3d (external)', itimer=itimer)
  call zeus3d (self%istep)  ! run ZEUS
  call trace_end(itimer)
  call CopyFromZeus (self)
  !-----------------------------------------------------------------------------
  call self%counter_update
END SUBROUTINE update_zeus

!===============================================================================
!===============================================================================
SUBROUTINE CopyToZeus (self)
  USE zparams
  USE amnodal
  class(mhd_t) :: self
  real(KindReal) :: x1low, x2low, x3low, dx1, dx2, dx3, pvelocity(3)
  logical, save :: first_time=.true.
  !$omp threadprivate (first_time)
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: qz, force, p
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, Ux, Uy, Uz
  integer, save:: itimer1=0, itimer2=0, itimer3=0, itimer4=0
  external abndyall, bndyflgs, bndyall, get_vel
  ! DISPATCH does not have explicit information about the interface of the ZEUScore
  ! procedure, so let's provide it. This allows for passing pointers to arrays
  ! to ZEUScore.
  interface
    subroutine copyatoz (val, mito, mjto, mkto, nvar, x1low, x2low, x3low, &
                         dx1, dx2, dx3, mptr, pvelocity, force)
       integer:: mito, mjto, mkto, nvar, mptr
       real(kind=8):: x1low, x2low, x3low, dx1, dx2, dx3
       real(kind=8), dimension(3):: pvelocity
       real(kind=8), dimension(:,:,:,:), pointer:: val, force
    end subroutine copyatoz
  end interface
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
  x1low = self%mesh(1)%r(1) + self%mesh(1)%centre_nat - 0.5d0 * dx1
  x2low = self%mesh(2)%r(1) + self%mesh(2)%centre_nat - 0.5d0 * dx2
  x3low = self%mesh(3)%r(1) + self%mesh(3)%centre_nat - 0.5d0 * dx3

  if (self%mesh(1)%type == 'Cartesian_mesh') then
    zcoordsys = 'XYZ'
  else if (self%mesh(1)%type == 'cylindrical_mesh') then
    zcoordsys = 'ZRP'
  else if (self%mesh(1)%type == 'spherical_mesh') then
    zcoordsys = 'RTP'
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

  call trace%begin ('copyatoz (external)', itimer=itimer1)
  call copyatoz (qz, self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, &
                 self%nv, x1low, x2low, x3low, dx1, dx2, dx3, self%id, &
                 pvelocity, force)
  call trace%end (itimer1)
  call trace%begin ('abndyall (external)', itimer=itimer2)
  call abndyall
  call trace%end (itimer2)
  call trace%begin ('bndyflgs (external)', itimer=itimer3)
  call bndyflgs
  call trace%end (itimer3)
  call trace%begin ('bndyall (external)', itimer=itimer4)
  call bndyall
  call trace%end (itimer4)
  if (first_time) first_time = .false.

  if (self%do_pebbles) then
    associate (v1=>self%vgas(:,:,:,1,self%it),v2=>self%vgas(:,:,:,2,self%it),v3=>self%vgas(:,:,:,3,self%it))
    call get_vel (self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, v1, v2, v3)
    end associate
    self%pressure(:,:,:,self%it) = self%mem(:,:,:,self%idx%e,self%it,1) * (self%gamma - 1.0)
  end if

  courno = self%courant
  dt = self%dt(self%it)
  prtime = self%time
  ix1x2x3 = self%ix1x2x3

  ! copy physical boundary logical flags
  lx1ipb = self%mesh(1)%lower_boundary
  lx1opb = self%mesh(1)%upper_boundary
  lx2ipb = self%mesh(2)%lower_boundary
  lx2opb = self%mesh(2)%upper_boundary
  lx3ipb = self%mesh(3)%lower_boundary
  lx3opb = self%mesh(3)%upper_boundary

  if (associated(self%force%selected)) deallocate (force)
  nullify (qz)

END SUBROUTINE CopyToZeus

!===============================================================================
SUBROUTINE CopyFromZeus (self)
  USE zparams
  USE amrfluxes, only: amemf1, amemf2, amemf3

  class(mhd_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: qz
  integer, save:: itimer=0
  interface
    subroutine copyztoa (val, mito, mjto, mkto, nvar)
       integer:: mito, mjto, mkto, nvar
       real(kind=8), dimension(:,:,:,:), pointer:: val
    end subroutine copyztoa
  end interface
  !-----------------------------------------------------------------------------
  qz => self%mem(:,:,:,:,self%new,1)

  call trace_begin ('copyztoa (external)', itimer=itimer)

  ! Save the EMFs as calculated by the solver.
  ! Note that what is actually stored is `E_i * g_i * dx_i`, where `E_i` are the
  ! actual EMFs (= v x B).
  self%mem(:,:,:,self%idx%ex,self%new,1) = amemf1(:,:,:)
  self%mem(:,:,:,self%idx%ey,self%new,1) = amemf2(:,:,:)
  self%mem(:,:,:,self%idx%ez,self%new,1) = amemf3(:,:,:)

  call copyztoa (qz, self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, self%nv)
  call trace_end (itimer)

  ! update time step.
  !$omp atomic write
  self%dtime = dt
  self%dtime = min(self%dtime,io%out_time)
  self%ix1x2x3 = ix1x2x3

  nullify (qz)

END SUBROUTINE CopyFromZeus

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
