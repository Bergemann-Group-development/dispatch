!===============================================================================
!> Staggered mesh solver using entropy per unit volume as the energy variable.
!>
!> Flop counting is enable via annotation at the end of lines. To count flops,
!> use a Linux pipeline such as
!>
!> % grep ops: ../../../solvers/stagger2/mhd_mod.f90 | sed 's/.*ops://' | sumcol 2
!>
!> where 'sumcol' is any script that sums up values in column 2
!===============================================================================
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
  USE divb_clean_mod
  USE task_mod
  USE timer_mod
  USE scratch_mod
  USE eos_mod
  implicit none
  private
  type, extends(extras_t), public:: mhd_t
    real, dimension(:,:,:),   pointer    :: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer    :: p, dpdt, B, dBdt
    real, dimension(:,:,:),   pointer    :: lnd, ee, pg=>null(), tt=>null()
    real, dimension(:,:,:),   allocatable:: qrad
    real:: nu(6)
    real:: cdtd
    real:: maxJ=0.0
    logical:: mhd=.false.
    logical:: first_time=.true.
  contains
    procedure:: pre_init
    procedure:: test
    procedure:: init
    procedure:: post_init
    procedure:: pde
    procedure:: pre_update
    procedure:: update
    procedure:: output
    procedure :: to_total_e
    procedure:: gas_pressure
    procedure:: gas_pressure_sub
    procedure :: nbor_relations
  end type
  logical, save:: unsigned=.true.
  logical, save:: do_test=.false.
  logical, save:: do_smooth=.false.
  logical, save:: do_maxwell = .true.
  logical, save:: do_temperature = .false.
  integer, save:: verbose=0
CONTAINS

!===============================================================================
!> Read the namelist info and set the default parameters.
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  type(index_t):: idx
  integer:: i, iv, iostat
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], csound=1.0, cdtd=0.5, courant=0.2
  logical, save:: first_time=.true.
  logical, save:: mhd=.true.
  namelist /stagger_params/ nu, courant, csound, cdtd, hardwire, &
    do_smooth, unsigned, do_maxwell, do_temperature, verbose, divb_cleaner, &
    do_test, mhd
  !----------------------------------------------------------------------------
  call trace%begin ('mhd_t%pre_init')
  call self%extras_t%pre_init
  !-----------------------------------------------------------------------------
  ! Now pass on these wishes to the next layers, which may modify them
  !-----------------------------------------------------------------------------
  call io%lock%set ('mhd_t%pre_init')
  if (first_time) then
    first_time = .false.
    mhd = self%mhd
    rewind (io%input)
    read (io%input, stagger_params, iostat=iostat)
    call io%namelist_warning ('stagger_params', iostat)
    write (io%namelist, stagger_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  call io%lock%unset ('mhd_t%pre_init')
  !-----------------------------------------------------------------------------
  ! Note that the experiment must set self%mhd, possibly depending on IC
  !-----------------------------------------------------------------------------
  self%mhd = mhd
  if (self%mhd) then
    self%kind = 'stagger2e_mhd_patch'
    self%nv = 8
  else
    self%kind = 'stagger2e_hd_patch'
    self%nv = 5
  end if
  self%nw = 2
  self%nu = nu
  self%csound = csound
  self%courant = courant
  self%cdtd = cdtd
  if (stagger%order==6) then
    self%ng = 5
  else if (stagger%order==4) then
    self%ng = 3
  else
    self%ng = 2
  end if
  call self%idx%init (energy=5, mhd=self%mhd)
  where (self%n<2) self%ng=0
  self%type = 'mhd_t'
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer :: i
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !-----------------------------------------------------------------------------
  ! Allocate memory and mesh, etc
  !-----------------------------------------------------------------------------
  call io%assert (self%nw == 2, &
    'mhd_t%pde: memory slot for time derivatives missing')
  call self%extras_t%init
  !-----------------------------------------------------------------------------
  ! Set the staggering, unsigned and pervolume flags before calling the rest
  ! of the cascade if %init(), where %gpatch_t%init may need them, in IC calls
  !-----------------------------------------------------------------------------
  do i=1,3
    self%mesh(i)%h(self%idx%px+i-1) = -0.5
    if (self%mhd) self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  if (verbose>1 .and. self%track) then
    do i=1,3
      print '("h:",8f5.2)', self%mesh(i)%h
    end do
  end if
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  !-----------------------------------------------------------------------------
  ! Initialize everything else, and do a sanity check
  !-----------------------------------------------------------------------------
  self%type = 'mhd_t'
  call non_ideal%init
  call divb%init
  !-----------------------------------------------------------------------------
  ! Self tests, only on one thread
  !-----------------------------------------------------------------------------
  if (do_test) then
    !$omp single
    call stagger_test
    call self%test
    !$omp end single
  end if
  call trace%end()
END SUBROUTINE init

!> =============================================================================
!> Post init, whatever is required
!> =============================================================================
SUBROUTINE post_init (self)
  class(mhd_t):: self
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! Call extras post-init
  ! ----------------------------------------------------------------------------
  call self%extras_t%post_init
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  self%boundaries%mhd = self%mhd
 call trace%end()
END SUBROUTINE post_init

!===============================================================================
SUBROUTINE test (self)
  class(mhd_t):: self
  integer      :: ix, iy, iz, m(3), n(3), n0(3), i1
  real         :: fx, fy, fz, eps
  logical      :: ok, allok
  class(mhd_t), pointer:: tmp
  !-----------------------------------------------------------------------------
  ! Allocate mem to a temporary patch
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%test')
  allocate (tmp)
  tmp%nv = self%nv
  tmp%mhd = self%mhd
  call tmp%idx%init (energy=5, mhd=tmp%mhd)
  call tmp%patch_t%init
  n0 = 2*(tmp%ncell/2)
  if (any(tmp%n /= n0)) then
    print *,io%hs
    print *, 'mhd_t%test not possible with n, no_mans_land =', tmp%n, tmp%no_mans_land
    print *,io%hs
    call tmp%dealloc
    call trace%end()
    return
  end if
  !-----------------------------------------------------------------------------
  ! Give it initial values (linear compression), also in the guard zones
  !-----------------------------------------------------------------------------
  do iz=tmp%mesh(3)%lb,tmp%mesh(3)%ub
    fz = sin(4.*cgs%pi*tmp%mesh(3)%r(iz)/tmp%size(3))
    do iy=tmp%mesh(2)%lb,tmp%mesh(2)%ub
      fy = sin(4.*cgs%pi*tmp%mesh(2)%r(iy)/tmp%size(2))
      do ix=tmp%mesh(1)%lb,tmp%mesh(1)%ub
        fx = sin(4.*cgs%pi*tmp%mesh(1)%r(ix)/tmp%size(1))
        tmp%mem(ix,iy,iz,tmp%idx%d,tmp%it,1) = 1.0+0.1*(fx+fy+fz)
        tmp%mem(ix,iy,iz,tmp%idx%s,tmp%it,1) = 0.1*(fx+fy+fz)
        tmp%mem(ix,iy,iz,tmp%idx%px,tmp%it,1) = fx
        tmp%mem(ix,iy,iz,tmp%idx%py,tmp%it,1) = fy
        tmp%mem(ix,iy,iz,tmp%idx%pz,tmp%it,1) = fz
        if (tmp%mhd) then
          tmp%mem(ix,iy,iz,tmp%idx%bx,tmp%it,1) = fy
          tmp%mem(ix,iy,iz,tmp%idx%by,tmp%it,1) = fz
          tmp%mem(ix,iy,iz,tmp%idx%bz,tmp%it,1) = fx
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Evaluate the PDE and print the time derivatives (which are placed in the 2nd
  ! time slot and the 2nd variable slot)
  !-----------------------------------------------------------------------------
  call tmp%pde
  m = (1+tmp%gn)/2
  allok = .true.
  eps = 3e-6/minval(tmp%ds)
  print *,io%hl
  print *,'minval(ds), eps =', minval(tmp%ds), eps
  n = nint(tmp%size*0.5_8/tmp%ds)
  do ix=tmp%mesh(1)%li,tmp%mesh(1)%ui
    i1 = mod(ix-tmp%mesh(1)%li+n(1),2*n(1)) + tmp%mesh(1)%li
    ok = all(abs(tmp%mem(ix,m(2),m(3),1:5,tmp%it,2)- &
                 tmp%mem(i1,m(2),m(3),1:5,tmp%it,2)) < eps)
    allok = allok .and. ok
    if (io%verbose>0 .or. io_unit%do_validate .or. .not.ok) &
      print 1,ix,tmp%mem(ix,m(2),m(3),1:5,tmp%it,2), &
                 tmp%mem(i1,m(2),m(3),1:5,tmp%it,2)
    1 format("mhd_t%test:",i4,1p,2(2x,5e14.5))
  end do
  do iy=tmp%mesh(2)%li,tmp%mesh(2)%ui
    i1 = mod(iy-tmp%mesh(2)%li+n(2),2*n(2)) + tmp%mesh(2)%li
    ok = all(abs(tmp%mem(m(1),iy,m(3),1:5,tmp%it,2)- &
                 tmp%mem(m(1),i1,m(3),1:5,tmp%it,2)) < eps)
    allok = allok .and. ok
    if (io%verbose>0 .or. io_unit%do_validate .or. .not.ok) &
      print 1,iy,tmp%mem(m(1),iy,m(3),1:5,tmp%it,2), &
                 tmp%mem(m(1),i1,m(3),1:5,tmp%it,2)
  end do
  do iz=tmp%mesh(3)%li,tmp%mesh(3)%ui
    i1 = mod(iz-tmp%mesh(3)%li+n(3),2*n(3)) + tmp%mesh(3)%li
    ok = all(abs(tmp%mem(m(1),m(2),iz,1:5,tmp%it,2)- &
                 tmp%mem(m(1),m(2),i1,1:5,tmp%it,2)) < eps)
    allok = allok .and. ok
    if (io%verbose>0 .or. io_unit%do_validate .or. .not.ok) &
      print 1,iz,tmp%mem(m(1),m(2),iz,1:5,tmp%it,2), &
                 tmp%mem(m(1),m(2),i1,1:5,tmp%it,2)
  end do
  if (allok) then
    print *,'mhd_t%test passed'
    print *,io%hl
  else
    print *,'mhd_t%test failed'
    print *,io%hl
    call io%abort('mhd_t%test')
  end if
  call tmp%dealloc
  deallocate (tmp)
  call trace%end()
END SUBROUTINE test

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:),   pointer:: d, e, dddt, dedt
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi, qr
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, pg
  real, dimension(:,:,:,:), pointer:: U, p, B, dpdt, dBdt
  real, dimension(:,:,:),   pointer:: Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   allocatable:: du, ss, cs, pa, pp, pb, Q, uu, divU
  real, dimension(:,:,:),   allocatable, target:: ee, lnd
  real, dimension(:,:,:),   allocatable:: Txx, Txy
  real, dimension(:,:,:),   allocatable::      Tyy, Tyz
  real, dimension(:,:,:),   allocatable:: Tzx,      Tzz
  real, dimension(:,:,:,:), allocatable, target:: emf, J
  real, dimension(:,:,:,:), allocatable, target:: fx, fy, fz, fd, ld, dd, fm
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true.
  !.............................................................................
  real(8):: ds(3), dsmax
  real:: u_max
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  associate (i => self%idx)
  self%d    => self%mem(:,:,:,i%d      ,self%it,1)
  self%e    => self%mem(:,:,:,i%e      ,self%it,1)
  self%p    => self%mem(:,:,:,i%px:i%pz,self%it,1)
  self%dddt => self%mem(:,:,:,i%d      ,self%it,2)
  self%dedt => self%mem(:,:,:,i%e      ,self%it,2)
  self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  d     => self%d
  e     => self%e
  p     => self%p
  dddt  => self%dddt
  dedt  => self%dedt
  dpdt  => self%dpdt
  px    => p(:,:,:,1)
  py    => p(:,:,:,2)
  pz    => p(:,:,:,3)
  dpdtx => dpdt(:,:,:,1)
  dpdty => dpdt(:,:,:,2)
  dpdtz => dpdt(:,:,:,3)
  if (self%mhd) then
    self%B    => self%mem(:,:,:,i%bx:i%bz,self%it,1)
    self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
    B     => self%B
    dBdt  => self%dBdt
    Bx    => B(:,:,:,1)
    By    => B(:,:,:,2)
    Bz    => B(:,:,:,3)
    dBdtx => dBdt(:,:,:,1)
    dBdty => dBdt(:,:,:,2)
    dBdtz => dBdt(:,:,:,3);
  end if
  end associate
  call allocate_vectors_a (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars_a (self%gn, lnd, pa, du, ee, cs, pp, pb, Q, uu, divU)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
   ! self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux => U(:,:,:,1)
  Uy => U(:,:,:,2)
  Uz => U(:,:,:,3)
  Ex => emf(:,:,:,1)
  Ey => emf(:,:,:,2)
  Ez => emf(:,:,:,3)
  Jx => J(:,:,:,1)
  Jy => J(:,:,:,2)
  Jz => J(:,:,:,3)
  fxx => fx(:,:,:,1); fxy => fx(:,:,:,2)
                    fyy => fy(:,:,:,2); fyz => fy(:,:,:,3)
  fzx => fz(:,:,:,1);                   fzz => fz(:,:,:,3)
  fdx => fd(:,:,:,1); fdy => fd(:,:,:,2); fdz => fd(:,:,:,3)
  ldx => ld(:,:,:,1); ldy => ld(:,:,:,2); ldz => ld(:,:,:,3)
  ddx => dd(:,:,:,1); ddy => dd(:,:,:,2); ddz => dd(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)                                                                  ! ops: 0 0 0 0 3
  ld = down(lnd)
  dd = exp(ld)                                                                  ! ops: 0 0 0 0 3
  U = p/dd                                                                      ! ops: 0 0 3
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
  !-----------------------------------------------------------------------------
  ee = e/d                                                                      ! ops: 0 0 1
  self%lnd => lnd
  self%ee  => ee
  call self%gas_pressure_sub (d, self%lnd, self%ee, pg)
  self%pg => pg
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  ! The diagonal fxx is valid in (2:gn-1)
  !-----------------------------------------------------------------------------
  call allocate_scalars_a (self%gn, Txx, Tyy, Tzz)
  ds = self%ds
  Txx = ddxup(ds,Ux)
  Tyy = ddyup(ds,Uy)
  Tzz = ddzup(ds,Uz)
  divU = Txx+Tyy+Tzz                                    ! convergence rate      ! ops: 0 2 0
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  dsmax = maxval(ds,self%n > 2)
  du = -dsmax*min(divU,0.0)                             ! positive part
  if (do_smooth) then
    pa = self%nu(2)*d*xyzsm(du**2)
  else
    pa = self%nu(2)*d*du**2                                                     ! ops: 0 3 0
  end if
  !-----------------------------------------------------------------------------
  ! Diagonal Reynolds stress
  !-----------------------------------------------------------------------------
  fxx = xup(Ux)**2
  fyy = yup(Uy)**2
  fzz = zup(Uz)**2
  uu = sqrt(fxx+fyy+fzz)
  self%u_max = 0.0
  if (self%mhd) then
    ! Diagonal Reynolds stress
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = xup(Bx)**2                                                             ! mhd: 0 1 0
    Ey = yup(By)**2                                                             ! mhd: 0 1 0
    Ez = zup(Bz)**2                                                             ! mhd: 0 1 0
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
    cs = sqrt((eos%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
    if (verbose > 3) then
      u_max = self%fmaxval(cs,outer=.true.)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(ca)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
    end if
    ! Isotropic pressure part
    pp = pg + pa + pb                                                           ! mhd: 1 0 0
    ! Total diagonal stress
    du = (dsmax*self%nu(1))*(cs + 2.*uu)                                        ! mhd: 0 2 0
    !either use Maxwell stress tensor or Lorentz force
    if (do_maxwell) then
      fxx = pp + d*(fxx - du*Txx) - Ex                                          ! mhd: 1 0 0
      fyy = pp + d*(fyy - du*Tyy) - Ey                                          ! mhd: 1 0 0
      fzz = pp + d*(fzz - du*Tzz) - Ez                                          ! mhd: 1 0 0
    else
      fxx = pp + d*(fxx - du*Txx)                                               ! mhd: 2 2 0
      fyy = pp + d*(fyy - du*Tyy)                                               ! mhd: 2 2 0
      fzz = pp + d*(fzz - du*Tzz)                                               ! mhd: 2 2 0
    end if
  else
    ! Sound speed + velocity
    cs = sqrt((eos%gamma*pg+pa)/d)                                             ! ops: 1 1 1 1 0
    if (verbose > 3) then
      u_max = self%fmaxval(cs,outer=.true.)
      self%u_max = max(self%u_max,u_max)
      print 1, 'u_max(cs)  ', u_max, self%u_max
    end if
    pp = pg + pa                                                                ! ops: 1 0 0
    du = (dsmax*self%nu(1))*(cs + 2.*uu)                                        ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! Total diagonal stress, with density factor
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = pp + d*(fxx - du*Txx)                                                 ! ops: 2 2 0
    fyy = pp + d*(fyy - du*Tyy)                                                 ! ops: 2 2 0
    fzz = pp + d*(fzz - du*Tzz)                                                 ! ops: 2 2 0
  end if
  self%u_phase = cs+uu
  u_max = self%fmaxval(cs+uu,outer=.true.)
  self%u_max = max(self%u_max,u_max)
  if (verbose > 2) print 1, 'u_max(cs+u)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)                                           ! ops: 2 5 0
  !if (verbose==1) print *, 'average(Q):', self%faver(Q), &
    !self%faver(e-d**eos%gamma/(eos%gamma-1.0)), self%faver(fm(:,:,:,1))
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux   ! ops: 0 9 0
  fm = p+fd                                             ! mass flux             ! ops: 3 0 0
  dddt = -div(ds,fm)                                    ! density derivative
  u_max = self%cdtd*dsmax*self%fmaxval(abs(dddt/d),outer=.true.)                ! ops: 0 0 1
  self%u_max = max(self%u_max,u_max)
  if (verbose > 2) print 1, 'u_max(dddt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !-----------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)                                                         ! ops: 0 1 0
  fyz = zdn(Uy)*ydn(Uz)                                                         ! ops: 0 1 0
  fzx = xdn(Uz)*zdn(Ux)                                                         ! ops: 0 1 0
  !-----------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !-----------------------------------------------------------------------------
  call allocate_scalars_a (self%gn, Txy, Tyz, Tzx)
  Txy = ddxdn(ds,Uy)+ddydn(ds,Ux)                                               ! ops: 0 1 0
  Tyz = ddydn(ds,Uz)+ddzdn(ds,Uy)                                               ! ops: 0 1 0
  Tzx = ddzdn(ds,Ux)+ddxdn(ds,Uz)                                               ! ops: 0 1 0
  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy                                 ! ops: 1 3 0
  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz                                 ! ops: 1 3 0
  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx                                 ! ops: 1 3 0
  Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &                             ! ops: 2 5 0
                               yup(zup(Tyz))**2 + &                             ! ops: 1 1 0
                               zup(xup(Tzx))**2)                                ! ops: 1 0 0
  call deallocate_scalars_a (Txy, Tyz, Tzx)
  !-----------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !-----------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy                                                       ! ops: 0 1 0 0 1
  fyz = exp(zdn(ldy))*fyz                                                       ! ops: 0 1 0 0 1
  fzx = exp(xdn(ldz))*fzx                                                       ! ops: 0 1 0 0 1
  !-----------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !-----------------------------------------------------------------------------
  if (self%mhd.and.do_maxwell) then
    fxy = fxy - ydn(Bx)*xdn(By)                                                 ! mhd: 1 1 0
    fyz = fyz - zdn(By)*ydn(Bz)                                                 ! mhd: 1 1 0
    fzx = fzx - xdn(Bz)*zdn(Bx)                                                 ! mhd: 1 1 0
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)                                                        ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
    fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)                                        ! ops: 1 2 0
    fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)                                        ! ops: 1 2 0
    fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)                                        ! ops: 1 2 0
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))                           ! ops: 2 2 0
    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))                           ! ops: 2 2 0
    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))                           ! ops: 2 2 0
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)                    ! ops: 2 0 0
  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)                    ! ops: 2 0 0
  dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)                    ! ops: 2 0 0
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(scratch%force_per_unit_mass)) then
    dpdt = dpdt + dd*scratch%force_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (associated(scratch%force_per_unit_volume)) then
    dpdt = dpdt + scratch%force_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of energy; reuse the mass flux for the energy flux
  !-----------------------------------------------------------------------------
  if (eos%gamma == 1d0) then                                                   ! unless isothermal
    dedt = dddt*ee                                                              ! ops: 0 1 0
  else
    block
    real:: afm
    afm = self%faver(fm(:,:,:,1))
    fm = fm*down(ee) - (self%nu(5)*dd*down(du))*ddown(ds,ee)                    ! ops: 1 4 0
    dedt = -div(ds,fm) - pg*divU + Q
    if (verbose == 1) print '(a,10f12.6)', 'average(Q):', &
      self%time, self%time*self%faver(Q), &
      self%faver(e-d**eos%gamma/(eos%gamma-1.0)), &
      self%faver(- pg*divU), &
      afm, self%faver(fm(:,:,:,1))
    end block
  end if
  call deallocate_scalars_a (Txx, Tyy, Tzz)
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    self%maxJ=maxval(abs(J))
    emf = self%nu(6)*edge(du)*J                                                 ! mhd: 0 6 0
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                     ! mhd: 2 3 0
    if (eos%gamma/=1.0) &
      dedt = dedt + Q                                                           ! mhd: 1 0 0
    !---------------------------------------------------------------------------
    ! Lorentz force, used instead of Maxwell stress.
    !---------------------------------------------------------------------------
    if (.not.do_maxwell) then
      ! borrowing fx for scratch cell centering Bx -> fxx, By -> fyy, Bz -> fzz
      fxx = xup(Bx)
      fyy = yup(By)
      fzz = zup(Bz)
      dpdtx = dpdtx + zup(jy)*xdn(fzz) - yup(jz)*xdn(fyy)
      dpdty = dpdty + xup(jz)*ydn(fxx) - zup(jx)*ydn(fzz)
      dpdtz = dpdtz + yup(jx)*zdn(fyy) - xup(jy)*zdn(fxx)
    end if
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
    call non_ideal%update (self%gn, minval(ds), d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  ! External heating / cooling, if any; e.g. radiative or conductive
  !-----------------------------------------------------------------------------
  if (associated(scratch%heating_per_unit_mass)) then
    dedt = dedt + d*scratch%heating_per_unit_mass                                  ! ops: 2 2 0
  end if
  if (associated(scratch%heating_per_unit_volume)) then
    if (verbose > 1) then
      associate (l3=>self%mesh(3)%li)
      write (io_unit%log,'(i6," heating:",1p,10e11.2)') self%id, &
      sum(sum(scratch%heating_per_unit_volume(:,:,l3:l3+9),1),1)
      end associate
    end if
    dedt = dedt + scratch%heating_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Time step limitation for rate of change dE/dt = (de/dt-(e/d)*dddt)/d
  !-----------------------------------------------------------------------------
  if (eos%gamma /= 1d0) then
     u_max = self%cdtd*dsmax*self%fmaxval(abs((dedt-ee*dddt)/d),outer=.true.)   ! ops: 1 1 1
    self%u_max = max(self%u_max,u_max)
    if (verbose > 2) print 1, 'u_max(dedt)', u_max, self%u_max
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors_a (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars_a (lnd, pa, du, ee, cs, pp, pb, Q, uu, divU)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
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
  real, dimension(:,:,:,:), pointer:: mem
  integer:: iz
  !-----------------------------------------------------------------------------
  ! Compute initial temperature
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pre_update')
  !-----------------------------------------------------------------------------
  ! Set boundary conditions on primary variables
  !-----------------------------------------------------------------------------
  mem => self%mem(:,:,:,:,self%it,1)
  call self%boundaries%primary (mem)
  !
  if (do_temperature .and. .not. associated(self%tt)) then
    allocate (real:: self%tt(self%gn(1),self%gn(2),self%gn(3)))
    allocate (           lnd(self%gn(1),self%gn(2),self%gn(3)))
    allocate (            ee(self%gn(1),self%gn(2),self%gn(3)))
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    e => self%mem(:,:,:,self%idx%e,self%it,1)
    lnd = alog(d)
    ee = e/d
    call eos%lookup_table (shape(self%tt), lnd=lnd, ee=ee, tt=self%tt, label='mhd')
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
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod, where boundary conditions may be inserted
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer:: mem
  real(kind=KindScalarVar), pointer:: B(:,:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Output the state after guard zone downloads, selfgravity, and boundary
  ! conditions. If another call is placed at an earlier point in the task handling
  ! the call here will have no effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  if (eos%gamma==1d0) then
    self%mem(:,:,:,self%idx%e,self%it,1) = &
    self%mem(:,:,:,self%idx%d,self%it,1)*self%csound**2
  end if
  if (self%mhd) then
    B => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
    call divB%clean(self%mesh,B)
  end if
  call self%pde
  !-----------------------------------------------------------------------------
  ! The call to courant_condition may optionally be interscepted in extras_t
  !-----------------------------------------------------------------------------
  call self%courant_condition
  if (verbose > 2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  call self%check_density(' after timestep')
  call self%counter_update
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

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
    call io%lock%set ('mhd_t%output')
    write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
      self%id, self%iout
    inquire(file=filename, exist=exists)
    if (exists) then
      open (io%data_unit, file=trim(filename), form='formatted', status='old', &
        position='append')
      write (io%data_unit,'(8f10.4)') self%nu
      close (io%data_unit)
    end if
  call io%lock%unset ('mhd_t%output')
  end if
  call trace%end()
END SUBROUTINE output

!===============================================================================
SUBROUTINE gas_pressure_sub (self, d, lnd, ee, pg)
  class(mhd_t):: self
  real, dimension(:,:,:):: pg
  real, dimension(:,:,:), pointer:: d, lnd, e, ee
  integer, save:: nprint=1
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%gas_pressure_sub')
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  call allocate_scalars(self%gn,ee)
  ee = e/d
  if (trim(self%eos)=='ideal') then
    if (eos%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = d*ee*(eos%gamma-1)
    end if
  else
    if (do_temperature) then
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg, tt=self%tt, label='mhd_gp_sub')
    else
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg, label='mhd_gp_sub2')
    end if
    if (verbose > 1 .or. nprint > 0) then
      nprint = nprint-1
      write(io%output,*) self%id, 'lookup_table: d =', minval(d), maxval(d)
      write(io%output,*) self%id, 'lookup_table:  ee =', minval(ee) , maxval(ee)
      write(io%output,*) self%id, 'lookup_table:  pg =', minval(pg) , maxval(pg)
      if (do_temperature) &
        write(io%output,*) self%id, 'lookup_table:  tt =', minval(self%tt), &
          maxval(self%tt)
    end if
  end if
  call deallocate_scalars(ee)
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
    if (eos%gamma==1d0) then
      pg = d*self%csound**2
    else
      pg = d*ee*(eos%gamma-1)
    end if
  else
    if (do_temperature) then
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg, tt=self%tt, label='mhd_gp')
    else
      call eos%lookup_table (shape(d), d=d, ee=ee, pg=pg, label='mhd_gp2')
    end if
    if (verbose > 1 .or. nprint > 0) then
      nprint = nprint-1
      write(io%output,*) self%id, 'lookup_table: lnd =', minval(d), maxval(d)
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

!===============================================================================
!> Set up nbor relations between MHD tasks
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(mhd_t):: self
  class(task_t), pointer:: another, task
  logical:: needed, needs_me, download, overlaps
  ! ----------------------------------------------------------------------------
  integer, save:: itimer=0
  class(mhd_t),  pointer:: nb_mhd
  ! ----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('mhd_t%nbor_relations', 2, itimer=itimer)
  task => another
  call self%task_t%nbor_relations(task,needed,needs_me,download)
  select type (another)
  class is (mhd_t)
    ! do nothing, all is done in the base task_t layer
  class default
    needed=.false.; needs_me=.false.; download=.false.
  end select

  if (timer%detailed) &
    call trace%end(itimer=itimer)

END SUBROUTINE nbor_relations

!===============================================================================
!> For Riemann solvers we need to convert internal energy to total energy
!===============================================================================
SUBROUTINE to_total_e(self, it)
  class(mhd_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n, gn
  integer:: ix, iy, iz, iz1, i, it
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz, emag
  real, dimension(:,:,:), pointer:: b1, b2, b3
  ! ----------------------------------------------------------------------------
  call io%abort('mhd_t: to_total_e not set up')

END subroutine to_total_e

END MODULE mhd_mod
