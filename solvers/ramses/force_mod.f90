!===============================================================================
!> Simple forcing module for tests
!===============================================================================
MODULE force_mod
  USE io_mod
  USE mpi_global_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE scaling_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    integer:: id
    integer :: seed
    integer :: verbose = 0
    real :: t_turb = 0.0
    real :: t_damp = 0.0
    real(8):: t_save=0d0
    real:: phase(3), a(3)
    real, pointer, contiguous, dimension(:,:,:,:,:)   :: fran => null()
    procedure(single_solenoidal), pointer:: selected=>null()
    class(mesh_t), dimension(:), pointer  :: mesh
    character(len=64):: solver
    type(random_t):: random
    real(8), dimension(4):: momentum=0d0
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: single_solenoidal
    procedure:: zero_force
  end type
  real :: ampl_turb = 0.1
  real :: a_helmh = 1.0
  real :: k1 = 1.0
  real :: k2 = 1.5
  real :: pk = 7./6.
  real :: t_turn = 0.1
  real :: t_turb = 0.0
  real :: t_damp = 0.01
  real :: pi=3.14159365
  logical :: do_force = .true.
  logical :: do_helmh = .true.
  integer :: verbose = 0
  integer :: seed = -3
  type(mpi_global_t):: global_momentum
CONTAINS

!===============================================================================
!> Choose type of initial condition, and initialize the random number generator.
!===============================================================================
SUBROUTINE init (self, solver, id, mesh)
  class(force_t):: self
  character(len=64):: solver
  integer:: id
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  integer:: iostat
  character(len=32), save:: type='none'
  logical, save:: first_time=.true.
  namelist /force_params/ type, seed
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, force_params, iostat=iostat)
    call io%namelist_warning ('force_params', iostat)
    write (io%namelist, force_params)
  end if
  call io%lock%unset
  select case (trim(type))
  case ('zero')
    self%selected => zero_force
  case ('single_solenoidal')
    self%selected => single_solenoidal
  case ('single_compressive')
    self%selected => single_compressive
  case ('turbulence')
    self%selected => turbulence
    call init_turbulence (self, solver, id, mesh)
  case default
    nullify(self%selected)
  end select
  self%solver = trim(solver)
  self%id = id
  call self%random%init (seed)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Empty procedure; nothing to deallocate
!===============================================================================
SUBROUTINE dealloc (self)
  class(force_t):: self
END SUBROUTINE dealloc

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!> Note that random%ran3() returns 3 numbers in the range (-1.0,1.0)
!===============================================================================
FUNCTION single_solenoidal (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  !.............................................................................
  integer        :: ix, iy, iz, iostat
  real(8)        :: kx, ky, kz
  real, save     :: t_turn=0.3, a0(3)=0.1, k(3)=1.0
  integer, save  :: itimer=0
  logical, save  :: first_time=.true.
  namelist /force_solenoidal_params/ k, a0, t_turn
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%single_solenoidal', itimer=itimer)
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, force_solenoidal_params, iostat=iostat)
    call io%namelist_warning ('force_solenoidal_params', iostat)
    write (io%namelist, force_solenoidal_params)
  end if
  call io%lock%unset
  if (time >= self%t_save) then
    self%t_save = self%t_save + t_turn
    self%a = a0*2./3.*(1.0+0.1*self%random%ran3())/t_turn
    self%phase = pi2*0.5*self%random%ran3()
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz = pi2*k(3)*(r3(iz)+m(3)%p)/m(3)%b + self%phase(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*k(2)*(r2(iy)+m(2)%p)/m(2)%b + self%phase(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*k(1)*(r1(ix)+m(1)%p)/m(1)%b + self%phase(1)
        ff(ix,iy,iz,1) = self%a(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,2) = self%a(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,3) = self%a(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END FUNCTION single_solenoidal

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!> Note that random%ran3() returns 3 numbers in the range (-1.0,1.0)
!===============================================================================
FUNCTION single_compressive (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  !.............................................................................
  integer        :: ix, iy, iz, iostat
  real(8)        :: kx, ky, kz
  logical, save  :: first_time=.true.
  real, save     :: t_turn=0.3, a(3)=0.1, k(3)=1.0, t_turb=1.0
  integer, save  :: itimer=0
  namelist /force_compressive_params/ k, a, t_turn, t_turb
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%single_compressive', itimer=itimer)
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, force_compressive_params, iostat=iostat)
    call io%namelist_warning ('force_compressive_params', iostat)
    write (io%namelist, force_compressive_params)
  end if
  call io%lock%unset
  self%t_turb = t_turb
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz = pi2*k(3)*(r3(iz)+m(3)%p)/m(3)%b + self%phase(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*k(2)*(r2(iy)+m(2)%p)/m(2)%b + self%phase(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*k(1)*(r1(ix)+m(1)%p)/m(1)%b + self%phase(1)
        ff(ix,iy,iz,1) = a(1)*sin(kx)
        ff(ix,iy,iz,2) = a(2)*sin(ky)
        ff(ix,iy,iz,3) = a(3)*sin(kz)
      end do
    end do
  end do
  end associate
  call trace_end (itimer)
END FUNCTION single_compressive

!===============================================================================
!> Initialize the force module
!===============================================================================
SUBROUTINE init_turbulence (self, solver, id, mesh)
  class(force_t):: self
  integer:: id
  character(len=64):: solver
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  integer :: mx, my, mz, iostat
  real    :: rhoav, xav, yav, zav
  real(8) :: momentum(4)
  namelist /turbulence_params/ do_force, do_helmh, a_helmh, k1, k2, pk, &
                               ampl_turb, t_turn, t_turb, seed, t_damp, verbose
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('force_t%init_turbulence')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, turbulence_params, iostat=iostat)
    call io%namelist_warning ('turbulence_params', iostat)
    write (io%namelist, turbulence_params)
    call global_momentum%init (4)
  end if
  if (t_turb == 0.) t_turb = -1.5*t_turn
  self%id     = id
  self%seed   = seed
  self%t_turb = t_turb
  self%t_damp = t_damp
  self%verbose = verbose
  self%mesh   => mesh
  mx = self%mesh(1)%gn
  my = self%mesh(2)%gn
  mz = self%mesh(3)%gn
  if (.not.associated(self%fran)) then
    allocate (self%fran(mx,my,mz,3,2))
    call io%bits_mem (storage_size(self%fran), product(shape(self%fran)), 'fran')
  end if
  self%fran = 0.0
  self%selected => turbulence
  call self%random%init (seed)
  call io%lock%unset
  call trace%end()
END SUBROUTINE init_turbulence

!===============================================================================
!> Calculate a random force, confined to a shell in k-space, at regular intervals
!> in time.  We assume rho, P, and the sound speed to be near unity.  The velocity
!> amplitude should be of the order ampl_turb.  The size (scale) of the drivings
!> motions is of the order of 1.0/k1.  Thus the turnover time is given by
!> t_turb*ampl_turb = 1.0/k1; and the acceleration is of the order ampl_turb/t_turb.
!>
!> Note:  In order to restart properly, one needs to retrace the force updates.
!>
!>  06-aug-93/aake:  bug corrected; seed was not negative at start
!>  09-aug-93/aake:  changed to allow restart w same random coeffs
!>  21-feb-95/paolo: changed in order to use segments of constant
!>                   time der. of acceleration (so continuous accel.)
!>  23-jan-96/aake:  bug corrected; ran1() -> 2.*ran1()-1. to span [-1,1]
!>
!>            Restart must be done AFTER t=t_turn.
!>            ------------------------------------
!===============================================================================
FUNCTION turbulence (self, time, d, p, Ux, Uy, Uz, mesh) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: mesh
  real, dimension(mesh(1)%gn,mesh(2)%gn,mesh(3)%gn,3)    :: ff
  !.............................................................................
  integer                :: i, j, k, jx, jy, jz, mx, my, mz, nrand, kmax
  complex                :: fxx, fyy, fzz, kx, ky, kz, corr
  complex                :: expikr
  real                   :: x, y, z, w, fact, fpow, accel
  real                   :: dx, dy, dz, eps, fk
  class(mesh_t), pointer :: m
  real(8)                :: momentum(4)
  integer                :: l(3), u(3)
  integer, save          :: itimer=0, nprint=5
  !-----------------------------------------------------------------------------
  if (.not.do_force) then
    ff(:,:,:,:) = 0.0
    return
  end if
  call trace_begin ('force_t%turbulence', itimer=itimer)
  eps = 1e-5
  mx = mesh(1)%gn
  my = mesh(2)%gn
  mz = mesh(3)%gn
  !-----------------------------------------------------------------------------
  !  Loop until self%t_turb > t
  !-----------------------------------------------------------------------------
  do while (self%t_turb <= time)
    self%t_turb = self%t_turb + t_turn
    !---------------------------------------------------------------------------
    ! Maintain two force snapshots
    !---------------------------------------------------------------------------
    self%fran(:,:,:,:,1) = self%fran(:,:,:,:,2)
    self%fran(:,:,:,:,2) = 0.0
    !---------------------------------------------------------------------------
    ! Count the number of wavenumbers
    !---------------------------------------------------------------------------
    nrand = 0
    kmax=ifix(k2+1.)
    do jx=-kmax,kmax
      do jy=-kmax,kmax
        do jz=-kmax,kmax
          fk = sqrt(float(jx**2+jy**2+jz**2))
          if (fk >= k1-eps .and. fk <= k2+eps) nrand=nrand+1
        enddo
      enddo
    enddo
    !---------------------------------------------------------------------------
    ! Normalized acceleration factor
    !---------------------------------------------------------------------------
    accel  = ampl_turb/t_turn/sqrt(float(nrand)/8.)    ! rms=1./8. per comp.
    !---------------------------------------------------------------------------
    ! To obtain a Kolmogorov slope of the driving alone, one should have amplitudes
    ! a(k) that drop with k^(-11./6.), to have a(k)^2*k^2 = k^(-5./3.).  This ASSUMES
    ! that the amplitudes are proportional to the driving.  But the energy drain may
    ! be rather inversely proportional to the turnover time, which goes as k^(2./3.)
    !---------------------------------------------------------------------------
    fpow = 0.
    do jx=-kmax,kmax
      kx = cmplx (0., jx*2.*pi/mesh(1)%b)
      do jy=-kmax,kmax
        ky = cmplx (0., jy*2.*pi/mesh(2)%b)
        do jz=-kmax,kmax
          kz = cmplx (0., jz*2.*pi/mesh(3)%b)
          fk = sqrt(float(jx**2+jy**2+jz**2))
          if (fk >= k1-eps .and. fk <= k2+eps) then
            fxx = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            fyy = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            fzz = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            !-------------------------------------------------------------------
            ! After refinement, fast-forward to the last two force calculations
            !-------------------------------------------------------------------
            if (self%t_turb <= time - 2.0*t_turn) &
              cycle
            !------------------
            ! solenoidal field:
            !------------------
            if (do_helmh) then
               corr=(kx*fxx+ky*fyy+kz*fzz)/(kx*kx+ky*ky+kz*kz+1e-20) 
               if (jx /= 0) fxx = fxx - a_helmh*corr*kx
               if (jy /= 0) fyy = fyy - a_helmh*corr*ky
               if (jz /= 0) fzz = fzz - a_helmh*corr*kz
            endif
            !------------------
            fact=1.
            if (jx /= 0) fact=fact*0.5
            if (jy /= 0) fact=fact*0.5
            if (jz /= 0) fact=fact*0.5
            fpow = fpow+fact*(cabs(fxx)**2+cabs(fyy)**2+cabs(fzz)**2)
            dx = mesh(1)%d
            dy = mesh(2)%d
            dz = mesh(3)%d
            do k=1,mz
              m => mesh(3)
              z = m%p + m%r(k)
              do j=1,my
                m => mesh(2)
                y = m%p + m%r(j)
                m => mesh(1)
                do i=1,mx
                  x = m%p + m%r(i)
                  expikr  = accel*cexp(kx*x+ky*y+kz*z)
                  self%fran(i,j,k,1,2) = self%fran(i,j,k,1,2) + real(fxx*expikr)
                  self%fran(i,j,k,2,2) = self%fran(i,j,k,2,2) + real(fyy*expikr)
                  self%fran(i,j,k,3,2) = self%fran(i,j,k,3,2) + real(fzz*expikr)
                end do
              end do
            end do
          endif
        enddo
      enddo
    enddo
    if (nprint>0 .and. self%id==1) then
      !$omp atomic
      nprint = nprint-1
      print '(a,i6,f10.4,i4,f10.4)', 'turbulence: id, t_turb, nrand, fpow =', &
        self%id, self%t_turb, nrand, fpow
    end if
  enddo
  !-----------------------------------------------------------------------------
  ! Update the total momentum vector, by adding the new mean and subtracting
  ! the previous mean
  !-----------------------------------------------------------------------------
  l = self%mesh%li
  u = self%mesh%ui
  block
    real(8):: new(4), c1, c2
    integer, save:: nprint=10
    !---------------------------------------------------------------------------
    ! Compute average momentum and density in this patch
    !---------------------------------------------------------------------------
    c1 = 1d0/product(u-l+1)
    new = 0d0
    do jz=l(3),u(3)
    do jy=l(2),u(2)
      do i=1,3
        new(i) = new(i) + sum(p(l(1):u(1),jy,jz,i))*c1
      end do
      new(4) = new(4) + sum(d(l(1):u(1),jy,jz))*c1
    end do
    end do
    !---------------------------------------------------------------------------
    ! Add and subtract the local contribution, weighted with the siz/box volume
    !---------------------------------------------------------------------------
    c2 = product(mesh%s/mesh%b)
    momentum = global_momentum%update (c2*(new - self%momentum))
    if (verbose > 0) then
      write (stdlog,'(i6,2x,a,3(2x,1p,4e12.3))') self%id, &
        'momentum = ', self%momentum, new, momentum
    else if (nprint > 0) then
      nprint = nprint-1
      write (stdlog,*) 'force: verbose =', verbose
    end if
    self%momentum = new
    if (any(isnan(new))) then
      print '(a,1p,5e12.3,2(2x,3i4))', 'NaN error:', new, c1, l, u
      !call io%abort ('NaN in force')
    end if
  end block
  !-----------------------------------------------------------------------------
  ! Time interpolation of the force
  !-----------------------------------------------------------------------------
  w = (t_turn-(self%t_turb-time))/t_turn
  w = 0.5*(1.-cos(w*pi))
  l = self%mesh%lb
  u = self%mesh%ub
  do i=1,3
    !---------------------------------------------------------------------------
    ! Friction force per unit mass, with t_damp in units of t_turn, which
    ! damps out bulk net momentum (momentum(4) is the average density)
    !---------------------------------------------------------------------------
    if (t_damp > 0.0 .and. time > 0d0) then
      accel = -momentum(i)/momentum(4)/(t_damp*t_turn)
    else
      accel = 0.0
    end if
    do jz=l(3),u(3)
    do jy=l(2),u(2)
    do jx=l(1),u(1)
      ff(jx,jy,jz,i) = accel + self%fran(jx,jy,jz,i,1) + &
        (self%fran(jx,jy,jz,i,2)-self%fran(jx,jy,jz,i,1))*w
    end do
    end do
    end do
  end do
  call trace_end (itimer)
END FUNCTION turbulence

!===============================================================================
!> Zero force
!===============================================================================
FUNCTION zero_force (self, time, d, p, Ux, Uy, Uz, mesh) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: mesh
  real, dimension(mesh(1)%gn,mesh(2)%gn,mesh(3)%gn,3)    :: ff
  !-----------------------------------------------------------------------------
  ff(:,:,:,:) = 0.0
END FUNCTION zero_force

END MODULE
