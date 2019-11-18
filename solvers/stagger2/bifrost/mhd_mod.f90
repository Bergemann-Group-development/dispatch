!=============================================================================
!> Staggered mesh solver using internal energy per unit volume as the energy variable.
!>
!> Flop counting is enable via annotation at the end of lines. To count flops,
!> use a Linux pipeline such as
!>
!> % grep ops: ../../../solvers/bifrost2/mhd_mod.f90 | sed 's/.*ops://' | sumcol 2
!>
!> where 'sumcol' is any script that sums up values in column 2
!=============================================================================
MODULE mhd_mod
  USE io_mod
  USE aux_mod
  USE bits_mod
  USE omp_mod
  USE math_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mesh_mod
  USE trace_mod
  USE kinds_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE vector_ops
  USE stagger_mod
  USE extras_mod
  USE index_mod
  USE timestep_mod
  USE eos_mod
  USE non_ideal_mod
  USE initial_mod
  !USE hyperdiffusion_mod
  implicit none
  private
  type, extends(extras_t), public:: mhd_t
    real, dimension(:,:,:),   pointer    :: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer    :: p, dpdt, B, dBdt
    real, dimension(:,:,:),   pointer    :: lnd, ee, pg, tt=>null()
    real, dimension(:,:,:),   allocatable:: qrad
    procedure(divB_clean1),   pointer    :: divB_clean => null()
    real:: nu(6)
    real:: cdtd
    real :: ca_max = 3e2
    real :: qmax = 8.0
    logical:: mhd=.true.
    logical:: first_time=.true.
  contains
    procedure:: init
    procedure:: pre_init
    procedure:: post_init
    procedure:: pde
    procedure:: pre_update
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
    procedure:: gas_pressure_sub
    procedure :: smooth_loc
    procedure :: cowboy_hat
    procedure, nopass:: cowboy_hat_1d
    procedure, nopass:: cowboy_hat_2d
    procedure, nopass:: cowboy_hat_3d
    procedure :: symmetric_x
    procedure :: symmetric_y
    procedure :: symmetric_z
    procedure :: symmetric_xyz
    procedure :: periodic_x
    procedure :: periodic_y
    procedure :: periodic_z
    procedure :: max_loc
    procedure :: difxyz_b3d
    procedure :: difxyz_u3d
    procedure :: quenchxyz
    procedure :: quenchx
    procedure :: quenchy
    procedure :: quenchz
  end type
  logical, save:: symmetric_e = .true.
  logical, save:: symmetric_b = .false.
  logical, save:: unsigned=.true.
  logical, save:: loginterpolation = .true.
  logical, save:: do_smooth=.false.
  logical, save:: do_temperature = .false.
  integer, save:: verbose=0
  integer, save:: divb_cleaner=2
  integer, save :: max_r = 3
  integer, save :: smooth_r = 3
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  ! ----------------------------------------------------------------------------
  integer:: i
  real, save:: nu(6), ca_max, csound=1.0, cdtd=0.5, courant=0.2
  real(8), save:: gamma=1.4_8
  logical, save:: first_time=.true.
  logical, save:: mhd=.true.
  character(len=16), save:: eos
  namelist /bifrost_params/ mhd, nu, courant, gamma, csound, cdtd, hardwire, eos, &
    ca_max, unsigned, verbose, divb_cleaner, loginterpolation, symmetric_e, symmetric_b, &
    max_r, smooth_r, do_temperature
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%pre_init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    eos = self%eos
    mhd = self%mhd
    rewind (io%input)
    read (io%input, bifrost_params)
    if (io%master) write (*, bifrost_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  !$omp end critical (input_cr)
  self%eos = eos
  self%mhd = mhd
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%courant = courant
  self%cdtd = cdtd
  !-----------------------------------------------------------------------------
  ! Read IC parameters, which determine if this is HD or MHD.  The IC init is
  ! called again, from extras_t%init, when gamma is known and can be passed on.
  !-----------------------------------------------------------------------------
  self%initial%mhd = self%mhd
  call self%idx%init (self%nv, self%mhd)
  call self%initial%init (self%kind, real(self%gamma))
  self%mhd = self%initial%mhd
  if (self%mhd) then
    self%kind = 'stagger2e_mhd_patch'
  else
    self%kind = 'stagger2e_hd_patch'
    self%idx%bx = -1
    self%idx%by = -1
    self%idx%bz = -1
  end if
  if (self%nv == 0) then
    if (self%mhd) then
      self%nv = 8
    else
      self%nv = 5
    end if
  end if
  ! ----------------------------------------------------------------------------
  ! Call extras pre-init
  ! ----------------------------------------------------------------------------
  call self%extras_t%pre_init
 call trace%end()
END SUBROUTINE pre_init

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  ! ----------------------------------------------------------------------------
  integer:: i
  logical, save:: first_time = .true.
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !-----------------------------------------------------------------------------
  ! Allocate memory and mesh, etc
  !-----------------------------------------------------------------------------
  call self%patch_t%init
  self%type = 'mhd_t'
  do i=1,3
    self%mesh(i)%h(self%idx%px+i-1) = -0.5
    if (self%mhd) self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  if (verbose>1 .and. self%track) then
    do i=1,3
      print '("h:",8f5.2)', self%mesh(i)%h
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Now the setup is complete enough to call gpatch_t%init, for possible ICs
  !-----------------------------------------------------------------------------
  call self%gpatch_t%init
  select case(divb_cleaner)
  case(1)
    self%divb_clean => divb_clean1
  case(2)
    self%divb_clean => divb_clean2
  case default
    call io%abort('mhd_mod:: divb cleaner method not understood')
  end select
  self%unsigned(self%idx%d)  = unsigned
  self%pervolume(self%idx%e) = unsigned
  if (first_time.and.io%master) then
    call stagger_test
    first_time = .false.
  end if
  ! ----------------------------------------------------------------------------
  ! init extras
  ! ----------------------------------------------------------------------------
  call self%extras_t%init
 call trace%end()
END SUBROUTINE init

!> =============================================================================
!> Post init, whatever that means
!> =============================================================================
SUBROUTINE post_init (self)
  class(mhd_t):: self
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! Call extras post-init
  ! ----------------------------------------------------------------------------
  call self%extras_t%post_init
 call trace%end()
END SUBROUTINE post_init

!> =============================================================================
SUBROUTINE pde (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:),   pointer:: d, e, dddt, dedt
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi, qr
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Uxc, Uyc, Uzc, pg
  real, dimension(:,:,:,:), pointer:: p, B, dpdt, dBdt
  real, dimension(:,:,:),   pointer:: Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   allocatable, target:: sij, Ca, fudge, lnd, cs, ee, Q, Tij
  real, dimension(:,:,:),   allocatable:: scr1, scr2, scr3, scr4, scr5, scr6, scr7
  real, dimension(:,:,:,:), allocatable, target:: ld, dd, U, U_center, Uv, Um, fm, J, emf, fd
  logical, save:: first_time=.true.
  !.............................................................................
  real(8):: ds(3), dsmax, dsmin, cu
  real:: u_max, safety, sawtooth_factor
  integer:: itimer=0
  integer:: lb(3), ub(3), k
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  associate (i => self%idx)
  self%d => self%mem(:,:,:,i%d      ,self%it,1); self%dddt => self%mem(:,:,:,i%d      ,self%it,2)
  self%e => self%mem(:,:,:,i%e      ,self%it,1); self%dedt => self%mem(:,:,:,i%e      ,self%it,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%it,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  !> -------------------------------------------------------------------
  d=>self%d; dddt=>self%dddt; e=>self%e; dedt=>self%dedt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  dddt = 0.0
  dedt = 0.0
  dpdt = 0.0
  !> -------------------------------------------------------------------
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%it,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3)
    dbdt = 0.0
  end if
  end associate
  ds   = self%ds
  dsmax = maxval(ds,self%n > 2)
  dsmin = minval(ds,self%n > 2)
  self%u_max = 0.0
  lb = self%mesh%lb
  ub = self%mesh%ub
  call allocate_scalars (self%gn, pg)
  call allocate_vectors_a (self%gn, ld, dd, U, U_center, Uv, Um, fm, fd)
  call allocate_scalars_a (self%gn, cs, ca, lnd, Sij, fudge, ee, Q, Tij)
  call allocate_scalars_a (self%gn, scr1, scr2, scr3, scr4, scr5, scr6, scr7)
  Ux=>U  (:,:,:,1); Uy=>U  (:,:,:,2); Uz=>U  (:,:,:,3)
  Uxc => U_center(:,:,:,1);Uyc => U_center(:,:,:,2);Uzc => U_center(:,:,:,3)
  call allocate_vectors_a (self%gn, J, emf)
  Ex=>emf(:,:,:,1); Ey=>emf(:,:,:,2); Ez=>emf(:,:,:,3)
  Jx=>J(:,:,:,1); Jy=>J(:,:,:,2); Jz=>J(:,:,:,3)
  
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)                                                                  ! ops: 0 0 0 0 3
  ld  = down(lnd)
  dd  = exp(ld)                                                                 ! ops: 0 0 0 0 3
  U   = p/dd                                                                    ! ops: 0 0 3      valid 4:n-2
  U_center = up(U)                                                              ! valid 3:n-3
  !-----------------------------------------------------------------------------
  ! Gas pressure. 
  !-----------------------------------------------------------------------------
  ee  = e/d
  ! calculate speed of sound.
  self%lnd => lnd
  self%ee  => ee
  call self%gas_pressure_sub (d, self%lnd, self%ee, pg)
  self%pg => pg
  cs  = sqrt(self%gamma*pg/d)
  Sij = cs**2
  !
  if (self%mhd) then                         
    ! Diagonal Reynolds stress
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    emf = up(B)                                                                 ! valid 3:n-3
    ! add magnetic pressure
    Ca  = sqrt(Sij + (Ex**2+Ey**2+Ez**2)/d)                                     ! alfven speed
  else
    emf = 0.0
    Ca  = sqrt(Sij)
  end if
  !
  if (self%Ca_max /= 0.0) then
    fudge = min(1.,self%Ca_max/Ca)
    Ca    = Ca * fudge
    fudge = fudge*fudge
  end if  
  !
  Uv(:,:,:,1) = self%nu(1)*Ca+self%nu(2)*abs(Uxc)
  Uv(:,:,:,2) = self%nu(1)*Ca+self%nu(2)*abs(Uyc)
  Uv(:,:,:,3) = self%nu(1)*Ca+self%nu(2)*abs(Uzc)
  
  scr1 = abs((Ca+U_center(:,:,:,1)))!/self%courant)
  scr2 = abs((Ca+U_center(:,:,:,2)))!/self%courant)
  scr3 = abs((Ca+U_center(:,:,:,3)))!/self%courant)
  scr4 = merge(scr1,scr2,mask=scr1>scr2)
  scr5 = merge(scr4,scr3,mask=scr4>scr2)
  cu = self%fmaxval(scr5,outer=.false.)
  !print *, 'u_max(cu)', cu, self%u_max
  self%u_max = max(self%u_max,cu)

  call self%difxyz_b3d(ds, Ex, Ey, Ez, Uxc, Uyc, Uzc, Um, self%nu(6))           ! valid 4:n-4
  
  Um = Um + Uv
  
  
  !call self%max_loc(Um(:,:,:,1),scr4)                                           ! valid 6:n-6
  !call self%max_loc(Um(:,:,:,2),scr5)
  !call self%max_loc(Um(:,:,:,3),scr6)
  
  !call self%smooth_loc(scr4,Um(:,:,:,1))                                        ! valid 6:n-6
  !call self%smooth_loc(scr5,Um(:,:,:,2))
  !call self%smooth_loc(scr6,Um(:,:,:,3))
  
  call self%cowboy_hat(Um(:,:,:,1),scr4)
  !call self%symmetric_xyz(scr4)
  Um(:,:,:,1)=scr4
  call self%cowboy_hat(Um(:,:,:,2),scr5)
  !call self%symmetric_xyz(scr5)
  Um(:,:,:,2)=scr5
  call self%cowboy_hat(Um(:,:,:,3),scr6)
  !call self%symmetric_xyz(scr6)
  Um(:,:,:,3)=scr6
  !
  call self%difxyz_u3d(ds, Ux, Uy, Uz, Uv, self%nu(3))                          ! valid 4:n-2
  
  !call self%max_loc(Uv(:,:,:,1),scr4)                                           ! valid 6:n-4
  !call self%max_loc(Uv(:,:,:,2),scr5)
  !call self%max_loc(Uv(:,:,:,3),scr6)
  
  !call self%smooth_loc(scr4,Uv(:,:,:,1))                                        ! valid 7:n-5
  !call self%smooth_loc(scr5,Uv(:,:,:,2))
  !call self%smooth_loc(scr6,Uv(:,:,:,3))
  
  call self%cowboy_hat(Uv(:,:,:,1),scr4)
  !call self%symmetric_xyz(scr4)
  Um(:,:,:,1)=scr4
  call self%cowboy_hat(Uv(:,:,:,2),scr5)
  !call self%symmetric_xyz(scr5)
  Um(:,:,:,2)=scr5
  call self%cowboy_hat(Uv(:,:,:,3),scr6)
  !call self%symmetric_xyz(scr6)
  Um(:,:,:,3)=scr6
  
  
  if (symmetric_e) then
    Uv(:,:,:,1)=max(Uv(:,:,:,1),Uv(:,:,:,2),Uv(:,:,:,3))
    Uv(:,:,:,2)=Uv(:,:,:,1)
    Uv(:,:,:,3)=Uv(:,:,:,1)

    Um(:,:,:,1)=max(Um(:,:,:,1),Um(:,:,:,2),Um(:,:,:,3))
    Um(:,:,:,2)=Um(:,:,:,1)
    Um(:,:,:,3)=Um(:,:,:,1)
  end if
  

  ! ----------------------------------------------------------------------------
  ! Mass conservation
  ! ----------------------------------------------------------------------------
  ! -- Density diffusion
  if ((self%nu(4) /= 0.0)) then
    if (loginterpolation) then
      fm = ddown1(ds,lnd)                                                       ! valid 4:n-2
    else
      fm = ddown1(ds,d)                                                          ! valid 4:n-2
    end if
    call self%quenchxyz(fm, fd)                                                 ! valid 6:n-4
    fd(:,:,:,1) = self%nu(4) * xdn1(Uv(:,:,:,1)) * fm(:,:,:,1) * fd(:,:,:,1) * ds(1)
    fd(:,:,:,2) = self%nu(4) * ydn1(Uv(:,:,:,2)) * fm(:,:,:,2) * fd(:,:,:,2) * ds(2)
    fd(:,:,:,3) = self%nu(4) * zdn1(Uv(:,:,:,3)) * fm(:,:,:,3) * fd(:,:,:,3) * ds(3)
    if (loginterpolation) then
      fd = fd*dd
    end if
    dddt = dddt + div (ds, fd)                                                  ! valid 8:n-7
    ! Energy diffusion as consequence of density diffusion
    !fm = down(log(e))                                                           ! valid 4:n-2
    !fm = fd*exp(fm)/dd                                                          ! valid 6:n-4
    dedt = dedt + div (ds, down1(ee)*fd)                                         ! valid 8:n-7
    ! Momentum diffusion as consequence of density diffusion in X
    fm = U_center*up1(fd)                                                       ! valid 8:n-7
    dpdtx = dpdtx + ddxdn(ds, fm(:,:,:,1)) + ddxdn(ds, fm(:,:,:,2)) + ddxdn(ds, fm(:,:,:,3))  ! valid 11:n-9
    dpdty = dpdty + ddydn(ds, fm(:,:,:,1)) + ddydn(ds, fm(:,:,:,2)) + ddydn(ds, fm(:,:,:,3))  ! valid 11:n-9
    dpdtz = dpdtz + ddzdn(ds, fm(:,:,:,1)) + ddzdn(ds, fm(:,:,:,2)) + ddzdn(ds, fm(:,:,:,3))  ! valid 11:n-9
  end if
  ! -- Advection of density --
  dddt = dddt - div(ds,p)                                                       ! valid 8:-7 if diff, 4:-2 if not
  ! print *,'umax(dddt)', u_max, self%u_max
  !
  ! ----------------------------------------------------------------------------
  ! Equation of motion
  ! ----------------------------------------------------------------------------  
  fm = dup1(ds,U)                                                               ! valid 3:-3
  call self%quenchxyz(fm,fd)                                                    ! valid 6:-4

  fd(:,:,:,1) = fd(:,:,:,1) * fm(:,:,:,1) * Uv(:,:,:,1) * ds(1) * d             !Txx  ! valid 7:-5
  fd(:,:,:,2) = fd(:,:,:,2) * fm(:,:,:,2) * Uv(:,:,:,2) * ds(2) * d             !Tyy  ! valid 7:-5
  fd(:,:,:,3) = fd(:,:,:,3) * fm(:,:,:,3) * Uv(:,:,:,3) * ds(3) * d             !Tzz  ! valid 7:-5
  dpdtx = dpdtx + ddxdn(ds, fd(:,:,:,1) - Uxc**2*d - pg)                      ! valid 11:-9
  dpdty = dpdty + ddydn(ds, fd(:,:,:,2) - Uyc**2*d - pg)                      ! valid 11:-9
  dpdtz = dpdtz + ddzdn(ds, fd(:,:,:,3) - Uzc**2*d - pg)                      ! valid 11:-9
  ! Contribution to internal energy gradient
  dedt  = dedt + (fd(:,:,:,1) - pg) * fm(:,:,:,1) &
               + (fd(:,:,:,2) - pg) * fm(:,:,:,2) &
               + (fd(:,:,:,3) - pg) * fm(:,:,:,3)                             ! valid 8:-7
  ! -- Cross terms --
  ! --    XY    --
  scr1 = ydn(dd(:,:,:,1))
  scr2 = ydn1(xdn1(Uv(:,:,:,2)))                                                ! valid 10:-7
  scr3 = ddydn1(ds,Ux)
  Sij  = scr3
  call self%quenchy(scr3,scr4)
  Tij(:,:,:) = scr2*scr3*scr4*ds(2)
  !
  scr2 = ydn1(xdn1(Uv(:,:,:,1)))
  scr3 = ddxdn1(ds,Uy)
  call self%quenchx(scr3,scr4)
  Sij(:,:,:) =  0.5*(scr3+Sij)         ! Sij         = 0.5*(ddx(Uy)+ddy(Ux))
  Tij(:,:,:) =  scr1*(Tij+scr2*scr3*scr4*ds(1))
  !
  scr2 = xdn(Uy)
  scr3 = ydn(Ux)
  scr4 = Tij - scr1*scr2*scr3
  !
  dpdtx = dpdtx + ddyup(ds,scr4)                                                ! valid 12:-9
  dpdty = dpdty + ddxup(ds,scr4)                                                ! valid 12:-9
  !
  scr3 = 2*Tij*Sij
  dedt = dedt+xup(yup(scr3))                                                    ! valid 12:-9
  ! --    XZ    --
  scr1 = zdn(dd(:,:,:,1))
  scr2 = zdn1(xdn1(Uv(:,:,:,3)))
  scr3 = ddzdn1(ds,Ux)
  Sij  = scr3
  call self%quenchz(scr3,scr4)
  Tij(:,:,:) = scr2*scr3*scr4*ds(3)
  !
  scr2 = zdn1(xdn1(Uv(:,:,:,1)))
  scr3 = ddxdn1(ds,Uz)
  call self%quenchx(scr3,scr4)
  Sij(:,:,:) =  0.5*(scr3+Sij)         ! Sij         = 0.5*(ddx(Uy)+ddy(Ux))
  Tij(:,:,:) =  scr1*(Tij+scr2*scr3*scr4*ds(1))
  !
  scr2 = xdn(Uz)
  scr3 = zdn(Ux)
  scr4 = Tij - scr1*scr2*scr3
  !
  dpdtx = dpdtx + ddzup(ds,scr4)
  dpdtz = dpdtz + ddxup(ds,scr4)
  !
  scr3 = 2*Tij*Sij
  dedt = dedt+xup(zup(scr3))
  ! --    YZ    --
  scr1 = zdn(dd(:,:,:,2))
  scr2 = zdn1(ydn1(Uv(:,:,:,3)))
  scr3 = ddzdn1(ds,Ux)
  Sij  = scr3
  call self%quenchz(scr3,scr4)
  Tij(:,:,:) = scr2*scr3*scr4*ds(3)
  !
  scr2 = zdn1(ydn1(Uv(:,:,:,2)))
  scr3 = ddydn1(ds,Uz)
  call self%quenchy(scr3,scr4)
  Sij(:,:,:) =  0.5*(scr3+Sij)         ! Sij         = 0.5*(ddx(Uy)+ddy(Ux))
  Tij(:,:,:) =  scr1*(Tij+scr2*scr3*scr4*ds(2))
  !
  scr2 = ydn(Uz)
  scr3 = zdn(Uy)
  scr4 = Tij - scr1*scr2*scr3
  !
  dpdty = dpdty + ddzup(ds,scr4)
  dpdtz = dpdtz + ddyup(ds,scr4)
  !
  scr3 = 2*Tij*Sij
  dedt = dedt+yup(zup(scr3))
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (allocated(self%force_per_unit_mass)) then
    dpdt = dpdt + dd*self%force_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (allocated(self%force_per_unit_volume)) then
    dpdt = dpdt + self%force_per_unit_volume
  end if
  ! ----------------------------------------------------------------------------
  ! Equation of energy
  ! ----------------------------------------------------------------------------
  ! -- Energy diffusion --
  if (self%nu(5) /= 0.0) then
    if (loginterpolation) then
      fm = ddown1(ds,log(ee))                                                   ! valid 2:n
    else
      fm = ddown1(ds,ee)                                                        ! valid 2:n
    end if
    call self%quenchxyz(fm,fd)                                                  ! valid 4:n-2

    fd(:,:,:,1) = dd(:,:,:,1)*self%nu(5)*fm(:,:,:,1)*fd(:,:,:,1)*xdn1(Uv(:,:,:,1))*ds(1)
    fd(:,:,:,2) = dd(:,:,:,2)*self%nu(5)*fm(:,:,:,2)*fd(:,:,:,2)*ydn1(Uv(:,:,:,2))*ds(2)
    fd(:,:,:,3) = dd(:,:,:,3)*self%nu(5)*fm(:,:,:,3)*fd(:,:,:,3)*zdn1(Uv(:,:,:,3))*ds(3)
    if (loginterpolation) then
      fm = exp(down(log(ee)))
      dedt = dedt + div(ds,fd*fm)
    else
      fm = down(ee)
      dedt = dedt + div(ds,fd)
    end if
  else
    if (loginterpolation) then
      fm = exp(down(log(ee)))
    else
      fm = down(ee)
    end if
  end if
  ! -- Energy advection --
  dedt = dedt - div(ds,fm*p)
  !-----------------------------------------------------------------------------
  ! External heating / cooling, if any; e.g. radiative or conductive
  !-----------------------------------------------------------------------------
  if (allocated(self%heating_per_unit_mass)) then
    dedt = dedt + d*self%heating_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (allocated(self%heating_per_unit_volume)) then
    dedt = dedt + self%heating_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  ! ----------------------------------------------------------------------------      
  if (self%mhd) then
    !J = curl_dn(ds,B)
    !if (self%Ca_max /= 0.0) J = fudge * J
    !emf = self%nu(6)*edge(Um)*J
    !
    J = curl_dn(ds,B)
    if (self%Ca_max /= 0.0) J = fudge * J
    ! -- Ex --
    scr1 = ddydn1(ds,Bz)
    scr2 = ddzdn1(ds,By)
    call self%quenchy(scr1, scr3)  !(Ix)
    call self%quenchz(scr2, scr4)  !(Ix)
    scr5 = zdn1(ydn1(Um(:,:,:,2)))                   !edge(Um_y)
    Ex   = ds(2)*scr5*scr3
    if (.not.symmetric_b) Ex = Ex*scr1
    !
    scr5 = zdn1(ydn1(Um(:,:,:,3)))                   !edge(Um_z)
    scr1 = ds(3)*scr5*scr4
    if (.not.symmetric_b) then
      Ex   = self%nu(6)*    (Ex - scr1*scr2)
    else
      Ex   = self%nu(6)*0.5*(Ex + scr1) * Jx
    end if
    ! -- Ey --
    scr1 = ddzdn1(ds,Bx)
    scr2 = ddxdn1(ds,Bz)
    call self%quenchz(scr1, scr3)  !(Iy)
    call self%quenchx(scr2, scr4)  !(Iy)
    scr5 = zdn1(xdn1(Um(:,:,:,3)))                   !edge(Um_z)
    Ey   = ds(3)*scr5*scr3
    if (.not.symmetric_b) Ey = Ey*scr1
    !
    scr5 = zdn1(xdn1(Um(:,:,:,1)))                   !edge(Um_x)
    scr1 = ds(1)*scr5*scr4
    if (.not.symmetric_b) then
      Ey   = self%nu(6)*    (Ey - scr1*scr2)
    else
      Ey   = self%nu(6)*0.5*(Ey + scr1) * Jy
    end if
    ! -- Ez --
    scr1 = ddxdn1(ds,By)
    scr2 = ddydn1(ds,Bx)
    call self%quenchx(scr1, scr3)  !(Iz)
    call self%quenchy(scr2, scr4)  !(Iz)
    scr5 = xdn(ydn(Um(:,:,:,1)))                   !edge(Um_x)
    Ez   = ds(1)*scr5*scr3
    if (.not.symmetric_b) Ez = Ez*scr1
    !
    scr5 = xdn1(ydn1(Um(:,:,:,2)))                   !edge(Um_y)
    scr1 = ds(2)*scr5*scr4
    if (.not.symmetric_b) then
      Ez   = self%nu(6)*    (Ez - scr1*scr2)
    else
      Ez   = self%nu(6)*0.5*(Ez + scr1) * Jz
    end if
    ! -- Qjoule --- yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))
    scr1 = Jx*Ex
    scr2 = Jy*Ey
    scr3 = Jz*Ez
    scr5 = yup(zup(scr1)+xup(scr3))
    Q    = merge(scr5,0.0,scr5>0.0)
    scr5 = zup(xup(scr2))
    Q    = Q + merge(scr5,0.0,scr5>0.0)
    !---------------------------------------------------------------------------
    ! Lorenz force = I x B, and electric field - u x B
    !---------------------------------------------------------------------------
    scr1 = ydn(Bz)
    scr2 = zdn(By)
    scr3 = zdn(Bx)
    scr4 = xdn(Bz)
    scr5 = xdn(By)
    scr6 = ydn(Bx)
    !
    Ex = Ex + ydn(Uz)*scr2 - zdn(Uy)*scr1                             
    Ey = Ey + zdn(Ux)*scr4 - xdn(Uz)*scr3                          
    Ez = Ez + xdn(Uy)*scr6 - ydn(Ux)*scr5                          
    !
    dpdtx = dpdtx + zup(Jy*scr4) - yup(Jz*scr5)
    dpdty = dpdty + xup(Jz*scr6) - zup(Jx*scr1)
    dpdtz = dpdtz + yup(Jx*scr2) - xup(Jy*scr3)
    !
    call non_ideal%update (self%gn, dsmin, d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dedt = dedt + Q
    !
    dBdt = -curl_up(ds,emf)
  end if
  !print *, 'dpdtx', minval(dpdtx), maxval(dpdtx), minloc(dpdtx), maxloc(dpdtx)
  !print *, 'ux', minval(ux), maxval(ux), minloc(ux), maxloc(ux)
  !print *, 'dddt', minval(dddt), maxval(dddt), minloc(dddt), maxloc(dddt)
  !print *, 'dedt', minval(dedt), maxval(dedt), minloc(dedt), maxloc(dedt)
  !print *, 'Uv', minval(uv), maxval(uv), minloc(uv), maxloc(uv)
  !print *, 'Um', minval(um), maxval(um), minloc(um), maxloc(um)
  !do k=self%mesh(1)%lb, self%mesh(1)%ub
  !  print *, k, 'uv', uv(k,self%mesh(2)%li,self%mesh(3)%li,1), 'dddt', dddt(k,self%mesh(2)%li,self%mesh(3)%li),'dpdtx', dpdtx(k,self%mesh(2)%li,self%mesh(3)%li)
  !end do
  !-----------------------------------------------------------------------------
  ! Calculate u_max
  !-----------------------------------------------------------------------------
  u_max = self%cdtd*dsmax*self%fmaxval(abs(dddt/d))                             ! ops: 0 0 1
  !print *, 'u_max(dddt)', u_max, self%u_max
  self%u_max = max(self%u_max,u_max)
  u_max = self%cdtd*dsmax*self%fmaxval(abs((dedt-ee*dddt)/d))                   ! ops: 1 1 1
  !print *, 'u_max(dedt)', u_max, self%u_max
  self%u_max = max(self%u_max,u_max)
  !
  sawtooth_factor = 6.2
  safety = self%qmax*sawtooth_factor*0.7/2.
  scr1 = max(abs(Uv(:,:,:,1)),abs(Uv(:,:,:,2)),abs(Uv(:,:,:,3)))
  u_max= self%fmaxval(scr1(:,:,:) *&
         safety)
  !print *,'umax(Uv)', u_max, self%u_max
  self%u_max = max(self%u_max,u_max)
  !
  scr1 = max(abs(Um(:,:,:,1)),abs(Um(:,:,:,2)),abs(Um(:,:,:,3)))
  u_max = self%cdtd*self%fmaxval(scr1(:,:,:) * &
           safety)
  !print *,'umax(Um)', u_max, self%u_max
  self%u_max = max(self%u_max,u_max)
  !do k=1,3
  !      u_max=self%cdtd*self%fmaxval(Uv(:,:,:,k)/ds(k),outer=.true.)
  !  if (minval(Uv) <= 0.0) call io%abort('uv leq 0')
  !  if (minval(Um) <= 0.0) call io%abort('um leq 0')
  !  print *,'umax(Uv)',k, u_max, self%u_max
  !  self%u_max = max(self%u_max,u_max)
  !  u_max = safety*self%fmaxval(Um(:,:,:,k)/ds(k),outer=.true.)
  !  print *,'umax(Um)',k, u_max, self%u_max
  !  self%u_max = max(self%u_max,u_max)
  !end do
  !-----------------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------------
  call deallocate_vectors_a (ld, dd, Uv, Um, U_center, fm, fd)
  call deallocate_scalars_a (cs, ca, lnd, Sij, fudge, ee, Q, Tij)
  call deallocate_scalars_a (scr1, scr2, scr3, scr4, scr5, scr6, scr7)
  if (self%mhd) call deallocate_vectors_a (J, emf)
  !  
  call deallocate_vectors_a (U)
  call deallocate_scalars (pg)
  
  !
  call trace%end (itimer)
  if (mpi%master .and. first_time) then
    print *, stagger%count, ' stagger calls', stagger%flops/product(real(self%n)),' flops per interior point'
    stagger%count=0
    stagger%flops=0
    first_time = .false.
  end if
END SUBROUTINE pde

!===============================================================================
!> This procedure is called by solver_t%update, before it calls mhd_t%update
!===============================================================================
SUBROUTINE pre_update (self)
  class(mhd_t), target:: self
  real, dimension(:,:,:), pointer:: d, e, lnd, ee
  integer:: iz
  !-----------------------------------------------------------------------------
  ! Compute initial temperature
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pre_update')
  if (do_temperature .and. .not. associated(self%tt)) then
    allocate (real:: self%tt(self%gn(1),self%gn(2),self%gn(3)))
    allocate (           lnd(self%gn(1),self%gn(2),self%gn(3)))
    allocate (            ee(self%gn(1),self%gn(2),self%gn(3)))
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    e => self%mem(:,:,:,self%idx%e,self%it,1)
    lnd = alog(d)
    ee = e/d
    call eos%lookup_table (shape(self%tt), lnd=lnd, ee=ee, tt=self%tt)
    deallocate (lnd, ee)
    if (verbose > 0) &
      write(io%output,*) self%id, 'mhd_t%pre_update:  tt =', &
        minval(self%tt), maxval(self%tt)
    call self%aux%register ('tt', self%tt)
  end if
  !-----------------------------------------------------------------------------
  ! Now we can call the extras_t%pre_update, which may need the temperature
  !-----------------------------------------------------------------------------
  call self%extras_t%pre_update
  call trace%end
END SUBROUTINE pre_update

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod, where boundary conditions may be inserted
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Output the state after guard zone downloads, selfgravity, and boundary
  ! conditions. If another call is placed at an earlier point in the task handling
  ! the call here will have no effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  if (self%gamma==1d0) then
    self%mem(:,:,:,self%idx%e,self%it,1) = &
    self%mem(:,:,:,self%idx%d,self%it,1)*self%csound**2
  end if
  !if (self%mhd) call self%divB_clean
  call self%pde
  !-----------------------------------------------------------------------------
  ! The call to courant_condition may optionally be interscepted in extras_t
  !-----------------------------------------------------------------------------
  call self%courant_condition
  if (verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  call self%check_density(' after timestep')
  call self%counter_update
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Perform a single iteration of a div(B) clean.  The limit of stability is
!> eps=1./6.  This is intended to prevent build-up of layers of large div(B)
!> between patches. An exact solution would be to solve an internal Laplace
!> equation, with boundary conditiosn = the boundary value of div(B), and then
!> subtract the gradient of the resulting potential.
!>
!> The problem can also be greatly reduced by using a higher order time
!> interpoation scheme, thus making the electric field computed from the
!> guard zone values of B and U a better match to the neighbors E-fields.
!===============================================================================
SUBROUTINE divb_clean1(self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: b
  real, dimension(:,:,:)  , pointer:: bx, by, bz, divb
  real, dimension(:,:,:,:)  , allocatable:: dndivb
  real, parameter:: eps=0.11
  !-----------------------------------------------------------------------------
  b  => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
  bx => b(:,:,:,1)
  by => b(:,:,:,2)
  bz => b(:,:,:,3)
  call allocate_scalars (self%gn, divb)
  call allocate_vectors_a (self%gn, dndivb)
  divb = ddup(bx,1) + ddup(by,2) + ddup(bz,3)
  bx = bx + eps*dddn(divb,1)
  by = by + eps*dddn(divb,2)
  bz = bz + eps*dddn(divb,3)
  print *, 'div_b clean min/max', minval(divb), maxval(divb)
  
  call deallocate_vectors_a (dndivb)
  call deallocate_scalars(divb)
    return
contains
  !-----------------------------------------------------------------------------
  function ddup(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = cshift(f,1,i) - f
  end function
  !-----------------------------------------------------------------------------
  function dddn(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = f - cshift(f,-1,i)
  end function
END SUBROUTINE divb_clean1

!===============================================================================
!> Enforce div(B)=0 in the guard zones, by succesively imposing the condition,
!> layer-by-layer, on the component perpendicular to each of the patch faces.
!===============================================================================
SUBROUTINE divb_clean2 (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: divb
  real, parameter:: eps=0.11
  integer:: l(3), u(3), ix, iy, iz, i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean2', itimer=itimer)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  l = self%mesh%lb
  u = self%mesh%ub-1
  i1 = merge(1, 0, self%n(1) > 1)
  j1 = merge(1, 0, self%n(2) > 1)
  k1 = merge(1, 0, self%n(3) > 1)
  !--- x-direction -------------------------------------------------------------
  do ix=self%mesh(1)%lo,self%mesh(1)%lb,-1
    bx(ix   ,l(2):u(2),l(3):u(3)) = &
    bx(ix+i1,l(2):u(2),l(3):u(3)) + ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  do ix=self%mesh(1)%ui,self%mesh(1)%ub-1
    bx(ix+i1,l(2):u(2),l(3):u(3)) = &
    bx(ix   ,l(2):u(2),l(3):u(3)) - ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  !--- y-direction -------------------------------------------------------------
  do iy=self%mesh(2)%lo,self%mesh(2)%lb,-1
    by(l(1):u(1),iy   ,l(3):u(3)) = &
    by(l(1):u(1),iy+j1,l(3):u(3)) + ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  do iy=self%mesh(2)%ui,self%mesh(2)%ub-1
    by(l(1):u(1),iy+j1,l(3):u(3)) = &
    by(l(1):u(1),iy   ,l(3):u(3)) - ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  !--- z-direction -------------------------------------------------------------
  do iz=self%mesh(3)%lo,self%mesh(3)%lb,-1
    bz(l(1):u(1),l(2):u(2),iz  ) = &
    bz(l(1):u(1),l(2):u(2),iz+k1) + ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  do iz=self%mesh(3)%ui,self%mesh(3)%ub-1
    bz(l(1):u(1),l(2):u(2),iz+k1) = &
    bz(l(1):u(1),l(2):u(2),iz   ) - ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  call trace%end (itimer)
END SUBROUTINE divb_clean2

!===============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!===============================================================================
SUBROUTINE output (self)
  class(mhd_t):: self
  character(len=64):: filename
  real(8):: out_next
  logical:: exists
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%output')
  out_next = self%out_next
  call self%extras_t%output
  if (self%time >= out_next .and. io%do_output .and. io%do_legacy) then
    !$omp critical (out_cr)
    write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
      self%id, self%iout
    inquire(file=filename, exist=exists)
    if (exists) then
      open (io%data_unit, file=trim(filename), form='formatted', status='old', &
        position='append')
      write (io%data_unit,'(8f10.4)') self%nu
      close (io%data_unit)
    end if
  !$omp end critical (out_cr)
  end if
  call trace%end
END SUBROUTINE output

!===============================================================================
SUBROUTINE gas_pressure_sub (self, d, lnd, ee, pg)
  class(mhd_t):: self
  real, dimension(:,:,:):: pg
  real, dimension(:,:,:), pointer:: d, lnd, ee
  integer, save:: nprint=1
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%gas_pressure_sub')
  if (trim(self%eos)=='ideal') then
    if (self%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = d*ee*(self%gamma-1)
    end if
  else
    if (do_temperature) then
      call eos%lookup_table (shape(lnd), lnd=lnd, ee=ee, pg=pg, tt=self%tt)
    else
      call eos%lookup_table (shape(lnd), lnd=lnd, ee=ee, pg=pg)
    end if
    if (verbose > 1 .or. nprint > 0) then
      nprint = nprint-1
      write(io%output,*) self%id, 'lookup_table: lnd =', minval(lnd), maxval(lnd)
      write(io%output,*) self%id, 'lookup_table:  ee =', minval(ee) , maxval(ee)
      write(io%output,*) self%id, 'lookup_table:  pg =', minval(pg) , maxval(pg)
      if (do_temperature) &
        write(io%output,*) self%id, 'lookup_table:  tt =', minval(self%tt), &
          maxval(self%tt)
    end if
  end if
  call trace%end()
END SUBROUTINE gas_pressure_sub

!===============================================================================
FUNCTION gas_pressure (self) result(pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  real, dimension(:,:,:), pointer:: d, lnd, e, ee
  integer, save:: nprint=1
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%gas_pressure_sub')
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  call allocate_scalars(self%gn,ee)
  ee = e/d
  if (trim(self%eos)=='ideal') then
    if (self%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = d*ee*(self%gamma-1)
    end if
  else
    if (do_temperature) then
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg, tt=self%tt)
    else
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg)
    end if
    if (verbose > 1 .or. nprint > 0) then
      nprint = nprint-1
      write(io%output,*) self%id, 'lookup_table: lnd =', minval(lnd), maxval(lnd)
      write(io%output,*) self%id, 'lookup_table:  ee =', minval(ee) , maxval(ee)
      write(io%output,*) self%id, 'lookup_table:  pg =', minval(pg) , maxval(pg)
      if (do_temperature) &
        write(io%output,*) self%id, 'lookup_table:  tt =', minval(self%tt), &
          maxval(self%tt)
    end if
  end if
  call deallocate_scalars(ee)
  call trace%end()
END FUNCTION gas_pressure

!> =============================================================================
!> =============================================================================
SUBROUTINE max_loc(self, f, g)
  class(mhd_t) :: self
  real, dimension (:,:,:) :: f, g
  ! ----------------------------------------------------------------------------
  if (max_r==3) then
    !g=stagger%max3(f)
    g=stagger%zup(f)
    !call self%symmetric_z(g)
    f=stagger%ymax3(g)
    !call self%symmetric_y(g)
    g=stagger%xmax3(f)
    !call self%symmetric_x(g)
  else if (max_r == 5) then
    !g=stagger%max5(f)
    g=stagger%zmax5(f)
    !call self%symmetric_z(g)
    f=stagger%ymax5(g)
    !call self%symmetric_y(g)
    g=stagger%xmax5(f)
    !call self%symmetric_x(g)
  else
    g = f
  end if
END SUBROUTINE max_loc

!> =============================================================================
!> =============================================================================
SUBROUTINE smooth_loc(self, f, g)
  class(mhd_t) :: self
  real, dimension (:,:,:) :: f, g
  ! ----------------------------------------------------------------------------
  if (smooth_r==3) then
    !g=stagger%xyzsm(f)
    g=stagger%zsm(f)
    !call self%symmetric_z(g)
    f=stagger%ysm(g)
    !call self%symmetric_y(g)
    g=stagger%xsm(f)
    !call self%symmetric_x(g)
  else if (smooth_r == 5) then
    !g=stagger%xyzsm(f)
    g=stagger%zsm(f)
    !call self%symmetric_z(g)
    f=stagger%ysm(g)
    !call self%symmetric_y(g)
    g=stagger%xsm(f)
    !call self%symmetric_x(g)
  else
    g=f
  end if
END SUBROUTINE smooth_loc

!> =============================================================================
!> The name of this subroutine is deliberately kept the same as it is in Bifrost
!> mhd sover, for clarity. 
!> Finite difference one zone left to one zones right along vector field bx,by,bz
!> =============================================================================
SUBROUTINE difxyz_b3d(self, ds, bx, by, bz, ux, uy, uz, g, cc)
  class(mhd_t) :: self
  real(8):: ds(3)
  real, dimension(:,:,:)   :: bx, by, bz
  real, dimension(:,:,:)   :: ux, uy, uz
  real, dimension(:,:,:,:) :: g
  real :: cc
  integer:: n(4), ix, iy, iz
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), allocatable :: udotb, bdotb
  ! ............................................................................
  !
  n=shape(g)
  associate (gn => self%mesh%gn)
  call allocate_scalars_a (gn, udotb, bdotb)
  udotb = ux*bx+uy*by+uz*bz                                                     ! valid 3:n-3
  bdotb = bx**2+by**2+bz**2+tiny(1.0)
  udotb = udotb/bdotb
  bx = ux - bx*udotb
  by = uy - by*udotb
  bz = uz - bz*udotb
  call deallocate_scalars_a (udotb, bdotb)
  !
  if (n(1) > 3) then
    g(   1,:,:,1) = 0.0
    g(n(1),:,:,1) = 0.0
    do iz=1,n(3)
      do iy=1,n(2)
        do ix=2,n(1)-1
          g(ix,iy,iz,1) = cc*0.5*max(0.,bx(ix-1,iy,iz) - bx(ix+1,iy,iz))
        end do
      end do
    end do
  else
    g(:,:,:,1) = 0.0
  endif
  !
  if (n(2) > 3) then
    g(:,   1,:,2) = 0.0
    g(:,n(2),:,2) = 0.0
    do iz=1,n(3)
      do iy=2,n(2)-1
        do ix=1,n(1)
          g(ix,iy,iz,2) = cc*0.5*max(0.,by(ix,iy-1,iz) - by(ix,iy+1,iz))
        end do
      end do
    end do
  else
    g(:,:,:,2) = 0.0
  endif
  !
  if (n(3) > 3) then
    g(:,:,   1,3) = 0.0
    g(:,:,n(3),3) = 0.0
    do iz=2,n(3)-1
      do iy=1,n(2)
        do ix=1,n(1)
          g(ix,iy,iz,3) = cc*0.5*max(0.,bz(ix,iy,iz-1) - bz(ix,iy,iz+1))
        end do
      end do
    end do
  else
    g(:,:,:,3) = 0.0
  endif
  end associate
END SUBROUTINE difxyz_b3d

!> =============================================================================
!> The name of this subroutine is deliberately kept the same as it is in Bifrost
!> mhd sover, for clarity. 
!> Finite difference one zone left to one zones right along vector field bx,by,bz
!> =============================================================================
SUBROUTINE difxyz_u3d(self, ds, ux, uy, uz, g, cc)
  class(mhd_t) :: self
  real(8) :: ds(3)
  real, dimension(:,:,:)   :: ux, uy, uz
  real, dimension(:,:,:,:) :: g
  integer:: ix, iy, iz, n(4)
  real :: cc
  ! ----------------------------------------------------------------------------
  n=shape(g)
  
  if (n(1) > 1) then
    do iz=1,n(3)
      do iy=1,n(2)
        do ix=1,n(1)-1
           g(ix,iy,iz,1) = g(ix,iy,iz,1) + cc*max(0.,(ux(ix,iy,iz) - ux(ix+1,iy,iz)))
        end do
      end do
    end do
  endif
  
  if (n(2) > 1) then
    do iz=1,n(3)
      do iy=1,n(2)-1
        do ix=1,n(1)
           g(ix,iy,iz,2) = g(ix,iy,iz,2) + cc*max(0.,(uy(ix,iy,iz) - uy(ix,iy+1,iz)))
        end do
      end do
    end do
  endif
  
  if (n(3) > 1) then
    do iz=1,n(3)-1
      do iy=1,n(2)
        do ix=1,n(1)
          g(ix,iy,iz,3) = g(ix,iy,iz,3) + cc*max(0.,(uz(ix,iy,iz) - uz(ix,iy,iz+1)))
        end do
      end do
    end do
  endif
END SUBROUTINE difxyz_u3d

!===============================================================================
!===============================================================================
subroutine quenchxyz (self, f, q)
  class (mhd_t) :: self
  real, intent(in), dimension(:,:,:,:) :: f
  real, intent(out), dimension(:,:,:,:) :: q
  ! ----------------------------------------------------------------------------
  call self%quenchx(f(:,:,:,1),q(:,:,:,1))
  call self%quenchy(f(:,:,:,2),q(:,:,:,2))
  call self%quenchz(f(:,:,:,3),q(:,:,:,3))
end subroutine quenchxyz

!===============================================================================
!===============================================================================
subroutine quenchx (self, f, qx)
  class (mhd_t) :: self
  real, intent(in), dimension(:,:,:) :: f
  real, intent(out), dimension(:,:,:) :: qx
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  real, dimension(:,:,:), allocatable :: d1,d2,fa,qq
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  real:: aqmax
  real(8)::ds(3)
  ! ............................................................................
  m => self%mesh
  aqmax=1./self%qmax
  ds = 1.
  if (self%qmax == 1.) then
    qx=1.
    return
  end if
  
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  
  allocate(d1(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(d2(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(fa(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(qq(m(1)%gn,m(2)%gn,m(3)%gn))
  
  if (m(1)%n>1) then 
    fa=abs(f)
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        do i=lb(1), ub(1)-1
          d1(i,j,k) = f(i+1,j,k) - f(i,j,k)                                     ! valid 4:n-3
        enddo
      end do
    end do
    
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        do i=2, ub(1)-1
          d2(i,j,k)=abs(d1(i,j,k) - d1(i-1,j,k))                                ! valid 5:n-3
        end do
      end do
    end do
    d2(lb(1)   ,:,:) = 0.0
    d2(ub(1)-1:,:,:) = 0.0    !
    !
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        do i=lb(1), ub(1)
          qq(i,j,k)=d2(i,j,k)/(fa(i,j,k)+aqmax*d2(i,j,k)+tiny(1.))
        enddo
        !
        qx(:,j,k) = cowboy_hat_1d(qq(:,j,k))
        !
      enddo
    enddo
  else
    qx(:,:,:)=0.0
  endif
  !
  deallocate(d2,fa,qq,d1)
 end subroutine quenchx
 
  !===============================================================================
!===============================================================================
subroutine quenchy (self, f, qy)
  class (mhd_t) :: self
  real, intent(in), dimension(:,:,:) :: f
  real, intent(out), dimension(:,:,:) :: qy
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  real, dimension(:,:,:), allocatable :: d1, d2,fa,qq
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  real:: aqmax
  real(8)::ds(3)
  ! ............................................................................
  m => self%mesh
  aqmax=1./self%qmax
  ds = 1.
  if (self%qmax == 1.) then
    qy=1.
    return
  end if
  
  allocate(d1(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(d2(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(fa(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(qq(m(1)%gn,m(2)%gn,m(3)%gn))
  
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  
  if (m(2)%n>1) then 
    fa=abs(f)
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)-1
        do i=lb(1), ub(1)
          d1(i,j,k) = f(i,j+1,k) - f(i,j,k)                                     ! valid 4:n-3
        enddo
      end do
    end do
    
    do k=lb(3), ub(3)
      do j=lb(2)+1, ub(2)-1
        do i=lb(1), ub(1)
          d2(i,j,k)=abs(d1(i,j,k) - d1(i,j-1,k))                                ! valid 5:n-3
        end do
      end do
    end do
    d2(:,lb(2)   ,:) = 0.0
    d2(:,ub(2)-1:,:) = 0.0    !
    
    do k=lb(3), ub(3)
      do j=lb(2),ub(2)
        do i=lb(1),ub(1)
          qq(i,j,k)=d2(i,j,k)/(fa(i,j,k)+aqmax*d2(i,j,k)+tiny(1.))
        end do
      end do
    end do
    !
    do k = lb(3), ub(3)
      do i=lb(1), ub(1)
        qy(i,:,k) = cowboy_hat_1d(qq(i,:,k))
      end do
    end do
  else
    qy(:,:,:)=0.0
  endif
  deallocate(d1,d2,fa,qq)
 end subroutine quenchy
 
  !===============================================================================
!===============================================================================
subroutine quenchz (self, f, qz)
  class (mhd_t) :: self
  real, intent(in), dimension(:,:,:) :: f
  real, intent(out), dimension(:,:,:) :: qz
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  real, dimension(:,:,:), allocatable :: d1,d2,fa,qq
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  real:: aqmax
  real(8)::ds(3)
  ! ............................................................................
  m => self%mesh
  aqmax=1./self%qmax
  ds=1.
  if (self%qmax == 1.) then
    qz=1.
    return
  end if
  
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  
  allocate(d1(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(d2(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(fa(m(1)%gn,m(2)%gn,m(3)%gn))
  allocate(qq(m(1)%gn,m(2)%gn,m(3)%gn))
  
  if (m(3)%n>1) then   
    fa=abs(f)
    do k=lb(3), ub(3)-1
      do j=lb(2), ub(2)
        do i=lb(1), ub(1)
          d1(i,j,k) = f(i,j,k+1) - f(i,j,k)                                     ! valid 4:n-3
        enddo
      end do
    end do
    
    do k=lb(3)+1, ub(3)-1
      do j=lb(2), ub(2)
        do i=lb(1), ub(1)
          d2(i,j,k)=abs(d1(i,j,k) - d1(i,j,k-1))                                ! valid 5:n-3
        end do
      end do
    end do
    d2(:,:,lb(3)   ) = 0.0
    d2(:,:,ub(3)-1:) = 0.0    !
    !
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        do i=lb(1),ub(1)
          qq(i,j,k)=d2(i,j,k)/(fa(i,j,k)+aqmax*d2(i,j,k)+tiny(1.))
        end do
      end do
    end do
   !

    do j=lb(2), ub(2)
      do i=lb(1), ub(1)
        qz(i,j,:) = cowboy_hat_1d(qq(i,j,:))
      end do
    end do

  else
    qz(:,:,:)=0.0
  endif
  !call self%symmetric_z(f)
  deallocate(d2,fa,qq)
 end subroutine quenchz
 
subroutine symmetric_xyz(self,f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  call self%symmetric_x(f)
  call self%symmetric_y(f)
  call self%symmetric_z(f)
end subroutine symmetric_xyz
 
!===============================================================================
!===============================================================================
subroutine symmetric_x (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  ! ............................................................................
  m => self%mesh
 
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  
  if (m(1)%n>1) then 
    !
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        ! lower symmetric guard zone extrapolation
        !do i=1,self%ng(1)
        !  f(li(1)-i,j,k) =f(li(1)+i-1,j,k)
        !  f(ui(1)+i,j,k) =f(ui(1)-i+1,j,k)
        !end do
        do i=lb(1),lo(1)
          f(i,j,k) =f(li(1)+(lo(1)-i),j,k)
        end do
        ! upper symmetric guard zone extrapolation
        do i=uo(1),ub(1)
          f(i,j,k) = f(ui(1)-(i-uo(1)),j,k)
        end do        
      enddo
    enddo
  endif
 end subroutine symmetric_x
 
  !===============================================================================
!===============================================================================
subroutine symmetric_y (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  ! ............................................................................
  m => self%mesh
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do

  if (m(2)%n>1) then 
    do k = lb(3), ub(3)
      ! lower symmetric guard zone extrapolation
      do j=lb(2),lo(2)
        do i=lb(1), ub(1)
          f(i,j,k) = f(i,li(2)+(lo(2)-j),k)
        end do
      end do
      ! upper symmetric guard zone extrapolation
      do j=uo(2),ub(2)
        do i=lb(1), ub(1)
          f(i,j,k) = f(i,ui(2)-(j-uo(2)),k)
        end do
      end do
     !do j=1,self%ng(2)
     !do i=lb(1), ub(1)
     !     f(i,li(2)-j,k) =f(i,li(2)+j-1,k)
     !     f(i,ui(2)+j,k) =f(i,ui(2)-j+1,k)
     ! end do    
    !end do
    end do
  endif
end subroutine symmetric_y
 
!===============================================================================
!===============================================================================
subroutine symmetric_z (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3)
  ! ............................................................................
  m => self%mesh

  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do

  if (m(3)%n>1) then   
    ! lower symmetric guard zone extrapolation
    do k = lb(3), lo(3)
      do j=lb(2),ub(2)
        do i=lb(1), ub(1)
          f(i,j,k) = f(i,j,li(3)+(lo(3)-k))
        end do
      end do
    end do
    ! upper symmetric guard zone extrapolation
    do k = uo(3), ub(3)
      do j=lb(2),ub(2)
        do i=lb(1), ub(1)
          f(i,j,k) = f(i,j,ui(3)-(k-uo(3)))
        end do
      end do    
    end do
    !do k=1,self%ng(3)
    !do j=lb(2), ub(2)
    !do i=lb(1), ub(1)
    !      f(i,j,li(3)-k) =f(i,j,li(3)+k-1)
    !      f(i,j,ui(3)+k) =f(i,j,ui(3)-k+1)
    !end do
    !end do
    !end do
  endif
 end subroutine symmetric_z

!===============================================================================
!===============================================================================
subroutine periodic_x (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3), ng
  ! ............................................................................
  m => self%mesh
 
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  ng=self%ng(1)
  if (m(1)%n>1) then 
    !
    do k=lb(3), ub(3)
      do j=lb(2), ub(2)
        do i=1,self%ng(1)
          f(lb(1)+i-1,j,k) =f(ui(1)-ng+i,j,k)
          f(ui(1)+i,  j,k) =f(lo(1)   +i,j,k)
        end do
      enddo
    enddo
  endif
 end subroutine periodic_x
 
!===============================================================================
!===============================================================================
subroutine periodic_y (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3), ng
  ! ............................................................................
  m => self%mesh
 
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  ng=self%ng(2)
  if (m(2)%n>1) then 
    !
    do k=lb(3), ub(3)
      do j=1, self%ng(2)
        do i=lb(1), ub(1)
          f(i,lb(2)+j-1,k) =f(i,ui(2)-ng+j,k)
          f(i,ui(2)+j,  k) =f(i,lo(2)   +j,k)
        end do
      enddo
    enddo
  endif
 end subroutine periodic_y

!===============================================================================
!===============================================================================
subroutine periodic_z (self, f)
  class (mhd_t) :: self
  real, dimension(:,:,:) :: f
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  integer :: i,j,k
  integer:: lb(3), li(3), ui(3), ub(3), lo(3), uo(3), ng
  ! ............................................................................
  m => self%mesh
 
  do i=1,3
    lb(i) = m(i)%lb; ub(i) = m(i)%ub
    lo(i) = m(i)%lo; uo(i) = m(i)%uo
    li(i) = m(i)%li; ui(i) = m(i)%ui
  end do
  ng=self%ng(3)
  if (m(3)%n>1) then 
    !
    do k=1, self%ng(3)
      do j=lb(2), ub(2)
        do i=lb(1), ub(1)
          f(i,j,lb(3)+k-1) =f(i,j,ui(3)-ng+k)
          f(i,j,ui(3)+k  ) =f(i,j,lo(3)   +k)
        end do
      enddo
    enddo
  endif
 end subroutine periodic_z
 
!*******************************************************************************

SUBROUTINE cowboy_hat(self,f,g)
  class(mhd_t):: self
  real, dimension(:,:,:) :: f,g
  real, dimension(:,:), allocatable :: f2
  real, dimension(:), allocatable :: f1
  ! ----------------------------------------------------------------------------
  class (mesh_t), dimension(:), pointer :: m
  ! ............................................................................
  integer:: itimer=0
  call trace%begin ('stagger_mod::cowboy_hat', itimer=itimer)
  m => self%mesh
  !
  if (all(self%n>2)) then
    g = cowboy_hat_3d(f,self%ds)
  else if ((self%n(1)>2).and.(self%n(2)>2)) then
    allocate (f2(size(g,1),size(g,2)))
    f2 = cowboy_hat_2d(f(:,:,m(3)%li),[m(1)%d,m(2)%d])
    g  = spread(f2,3,size(g,3))
    deallocate (f2)
  else if ((self%n(1)>2).and.(self%n(3)>2)) then
    allocate (f2(size(g,1),size(g,3)))
    f2 = cowboy_hat_2d(f(:,m(2)%li,:),[m(1)%d,m(3)%d])
    g  = spread(f2,2,size(g,2))
    deallocate (f2)
  else if ((self%n(2)>2).and.(self%n(3)>2)) then
    allocate (f2(size(g,2),size(g,3)))
    f2 = cowboy_hat_2d(f(m(1)%li,:,:),[m(2)%d,m(3)%d])
    g  = spread(f2,1,size(g,1))
    deallocate (f2)
  else if (self%n(1)>2) then
    allocate (f1(size(g,1)))
    allocate (f2(size(g,1),size(g,2)))
    f1 = cowboy_hat_1d(f(:,m(2)%li,m(3)%li))
    f2 = spread(f1,2,size(g,2))
    g  = spread(f2,3,size(g,3))
    deallocate(f1,f2)
  else if (self%n(2)>2) then
    allocate (f1(size(g,2)))
    allocate (f2(size(g,1),size(g,3)))
    f1 = cowboy_hat_1d(f(m(1)%li,:,m(3)%li))
    f2 = spread(f1,1,size(g,1))
    g  = spread(f2,3,size(g,3))
    deallocate(f1,f2)
  else if (self%n(3)>2) then
    allocate (f1(size(g,3)))
    allocate (f2(size(g,2),size(g,3)))
    f1 = cowboy_hat_1d(f(m(1)%li,m(2)%li,:))
    f2 = spread(f1,2,size(g,2))
    g  = spread(f2,1,size(g,1))
    deallocate(f1,f2)
  else
    call io%abort("mhd_mod::cowboy_hat failed finding the right dimensions")
  end if
  call trace%end (itimer)
END SUBROUTINE cowboy_hat


FUNCTION cowboy_hat_1d(a) RESULT (b)
  real, dimension(:), intent(in):: a
  real, dimension(size(a)):: b
  integer:: i, n
  real:: c1, c2
  integer:: itimer=0
  call trace%begin ('stagger_mod::cowboy_hat_1d', itimer=itimer)
  n = size(a)
  c1=0.2
  c2=0.4

  do i=2,n-1
    b(i) = a(i-1)*c2 + a(i)*c1 + a(i+1)*c2
  end do
  b(1)=b(2)
  b(n)=b(n)-1
  call trace%end (itimer)
END FUNCTION cowboy_hat_1d
 
FUNCTION cowboy_hat_2d(a,ds) RESULT (b)
  real, dimension(:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2)):: b
  real(8):: ds(2), nds(2)
  integer:: ix, iy, mx, my, i, n(2)
  real(8), dimension(3,3):: m
  integer:: itimer=0
  call trace%begin ('stagger_mod::cowboy_hat_2d', itimer=itimer)
  n = shape(a)
  nds=ds/minval(ds)
  m(:,1)=[0.4/sqrt(nds(1)**2+nds(2)**2),0.4/nds(2),0.4/sqrt(nds(1)**2+nds(2)**2)]
  m(:,2)=[0.4/nds(1),0.2_8,0.4/nds(1)]
  m(:,1)=[0.4/sqrt(nds(1)**2+nds(2)**2),0.4/nds(2),0.4/sqrt(nds(1)**2+nds(2)**2)]
  m=m/sum(m)
  !
  b=0.0
  do iy=2,n(2)-1
    do ix=2,n(1)-1
      do my=1,3
        do mx=1,3
          b(ix, iy) = b(ix, iy) + a (ix-2+mx,iy-2+my)*m(mx,my)
        end do
      end do
    end do
  end do
  b(:,1)   =a(:,2)
  b(:,n(2))=a(:,n(2)-1)
  b(1,   :)=a(2,   :)
  b(n(1),:)=a(n(1)-1,:)
  
  call trace%end (itimer)
END FUNCTION cowboy_hat_2d
 
FUNCTION cowboy_hat_3d(a,ds) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), nds(3)
  integer:: ix, iy, iz, mx, my, mz, n(3)
  real:: c1, c2, c3
  real(8), dimension(3,3,3):: m
  integer:: itimer=0
  call trace%begin ('stagger_mod::cowboy_hat_3d', itimer=itimer)
  n = shape(a)
  nds=ds/minval(ds)
  
  do mz = -1,1,1
    do my = -1,1,1
      do mx = -1,1,1
        m(mx+2,my+2,mz+2)= 0.4/(sqrt(nds(1)**2*mx**2+nds(2)**2*my**2+ &
                                     nds(3)**2*mx**2)+tiny(1.0_8))
      end do
    end do
  end do
  m(2,2,2) = 0.2
  m = m/sum(m)
  !
  b=0.0
  do iz=2,n(3)-1
    do iy=2,n(2)-1
      do ix=2,n(1)-1
        do mz=1,3
          do my=1,3
            do mx=1,3
              b(ix,iy,iz) = b(ix,iy,iz) + a (ix-2+mx,iy-2+my,iz-2+mz)*m(mx,my,mz)
            end do
          end do
        end do
      end do
    end do
  end do
  
  b(:,:,1   )=b(:,:,2   )
  b(:,:,n(3))=b(:,:,n(3)-1)
  b(:,1   ,:)=b(:,2   ,:)
  b(:,n(2),:)=b(:,n(2)-1,:)
  b(1,   :,:)=b(2,   :,:)
  b(n(1),:,:)=b(n(1)-1,:,:)
  
  call trace%end (itimer)
END FUNCTION cowboy_hat_3d
 
 
 
END MODULE mhd_mod
