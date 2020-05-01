!===============================================================================
!> Generate an initial state as per case A of Politano et al (1995), with
!> additional random components added to the velocity and magnetic field.
!>
!> The module uses as little as possible from the general DISPATCH set of
!> modules, and has the option to use random numbers read from a file
!> random.dat.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: id=0
    integer:: seed
    logical:: mhd
    real(8):: gamma=1d0
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type
  integer, save:: iran=0
  integer, save:: seed=11
  integer, save:: seed1
  integer, save:: id=0
  real, save:: alpha=0.8
  real, save:: beta=1.0
  real, save:: epsilon=0.1
  real, allocatable:: ran(:)
  logical, save:: first_time=.true.
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  logical, save:: read_random=.false.
  type(index_t):: idx
CONTAINS

!===============================================================================
!> Initialize the initial condition module via namelist parameters
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  !............................................................................
  logical, save:: mhd=.true.
  character(len=32), save:: type='politano+95'
  namelist /initial_params/ alpha, beta, epsilon, seed, read_random
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if
  self%solver = trim(solver)
  call io%lock%unset
  self%mhd = .true.
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Ancient Numerical Recipes code, for transportability
!===============================================================================
FUNCTION ran1()
  real(8):: ran1
  integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
  real(8):: AM,EPS,RNMX
  parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer:: k,iy
  !.............................................................................
  idum = seed1
  k =idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum<0) idum = idum+IM
  iy = idum
  ran1 = min(AM*iy,RNMX)
  seed1 = idum
END FUNCTION

!===============================================================================
!> Generate and save (or read) random numbers
!===============================================================================
FUNCTION rand()
  real(8):: rand
  iran = iran+1
  if (read_random) then
    read (15,*) rand
  else
    rand = ran1()
    ran(iran) = rand
  end if
END FUNCTION rand

!===============================================================================
!> Generate 3D random vector in [0,1]
!===============================================================================
FUNCTION rand3()
  real(8):: rand3(3)
  integer:: j
  !.............................................................................
  do j=1,3
    rand3(j) = rand()
  end do
END FUNCTION rand3

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase.  The
!> mesh_t data type contains box sizes m%b(1:3), 3D position m%p(1:3), and cell
!> coordinates m%r(i), where the index limits are m%lb and m%ub.
!===============================================================================
SUBROUTINE condition (self, m, f)
  class(initial_t)       :: self
  class(mesh_t), pointer :: m(:)
  real, pointer          :: f(:,:,:,:)
  !.............................................................................
  real(8), allocatable   :: ff(:,:,:,:)
  integer                :: kx, ky, kz, ix, iy, iz, j
  real(8)                :: a, k, phi(3), phi0(3), aa(3), bb(3), kk(3)
  logical, save          :: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Initialize in an OMP critical region, to allow writing/reading random numbers
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  !$omp critical (initial_cr)
  id = id+1
  self%id = id
  seed1 = seed
  print *, 'initializing patch', self%id, omp%thread
  allocate (ff(size(f,1),size(f,2),size(f,3),6))
  !-----------------------------------------------------------------------------
  ! Politano et al (1995), case A
  !-----------------------------------------------------------------------------
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    phi(3)   =  pi2*(r3(iz)+m(3)%p)/m(3)%b
    do iy=m(2)%lb,m(2)%ub
      phi(2)   =  pi2*(r2(iy)+m(2)%p)/m(2)%b
      do ix=m(1)%lb,m(1)%ub
        phi(1)   =  pi2*(r1(ix)+m(1)%p)/m(1)%b
        f(ix,iy,iz,idx%d ) = beta
        f(ix,iy,iz,idx%e ) = beta
        ff(ix,iy,iz,1) = -2.0*sin(phi(2))                         ! ux
        ff(ix,iy,iz,2) = +2.0*sin(phi(1))                         ! uy
        ff(ix,iy,iz,3) = 0.0                                      ! uz
        ff(ix,iy,iz,4) = alpha*(-2.*sin(2.*phi(2))+sin(phi(3)))   ! bx
        ff(ix,iy,iz,5) = alpha*(+2.*sin(   phi(1))+sin(phi(3)))   ! by
        ff(ix,iy,iz,6) = alpha*(    sin(   phi(1))+sin(phi(2)))   ! bz
      end do
    end do
  end do
  f(:,:,:,idx%px:idx%pz) = ff(:,:,:,1:3)
  f(:,:,:,idx%bx:idx%bz) = ff(:,:,:,4:6)
  if (epsilon > 0.0) then
    !---------------------------------------------------------------------------
    ! Count random numbers needed
    !---------------------------------------------------------------------------
    iran = 0
    do kz=-3,3
    do ky=-3,3
    do kx=-3,3
      kk = [kx,ky,kz]
      if (all(kk==0)) cycle                           ! skip k==0 component
      k = sqrt(sum(kk**2))                            ! wavenumber magnitude
      if (k < 3D0+tiny(1.0)) then                     ! shell in k-space
        iran = iran+9                                 ! three 3D random numbers
      end if
    end do
    end do
    end do
    allocate (ran(iran))
    !---------------------------------------------------------------------------
    ! Generate random components
    !---------------------------------------------------------------------------
    open (15,file='random.dat',form='formatted',status='unknown')
    if (read_random) then
      read (15,*)
    end if
    iran = 0
    do kz=-3,3
    do ky=-3,3
    do kx=-3,3
      kk = [kx,ky,kz]
      if (all(kk==0)) cycle                               ! skip k==0 component
      k = sqrt(sum(kk**2))                                ! wavenumber magnitude
      if (k < 3D0+tiny(1.0)) then                         ! shell in k-space
        aa = rand3()/k**2                                 ! random 3D vector ~ k**(-2)
        bb = rand3()/k**2                                 ! random 3D vector ~ k**(-2)
        aa = aa - sum(kk*aa)/sqrt(sum(kk*kk)*sum(aa*aa))  ! incompressible U
        bb = bb - sum(kk*bb)/sqrt(sum(kk*kk)*sum(bb*bb))  ! divergence-free B
        aa = aa*epsilon
        bb = bb*epsilon
        phi0 = pi2*rand3()                                ! random phase offsets
        do iz=m(3)%lb,m(3)%ub
          phi(3)   =  pi2*kz*(r3(iz)+m(3)%p)/m(3)%b
          do iy=m(2)%lb,m(2)%ub
            phi(2)   =  pi2*ky*(r2(iy)+m(2)%p)/m(2)%b
            do ix=m(1)%lb,m(1)%ub
              phi(1)  = pi2*kx*(r1(ix)+m(1)%p)/m(1)%b
              ff(ix,iy,iz,1:3) = ff(ix,iy,iz,1:3) + aa*cos(phi+phi0)
              ff(ix,iy,iz,4:6) = ff(ix,iy,iz,4:6) + bb*cos(phi+phi0)
            end do
          end do
        end do
      end if
    end do
    end do
    end do
    f(:,:,:,idx%px:idx%pz) = ff(:,:,:,1:3)
    f(:,:,:,idx%bx:idx%bz) = ff(:,:,:,4:6)
    !---------------------------------------------------------------------------
    ! Write a file with the random numbers needed
    !---------------------------------------------------------------------------
    if (first_time) then
      first_time = .false.
      if (read_random) then
        print *,'number of random numbers read:', iran
      else
        write (15,*) iran
        do j=1,iran
          write (15,'(f19.16)') ran(j)
        end do
        print *,'number of random numbers written:', iran
      end if
    end if
    close (15)
    deallocate (ran)
  end if
  end associate
  !$omp end critical (initial_cr)
  call trace_end
END SUBROUTINE condition

END MODULE
