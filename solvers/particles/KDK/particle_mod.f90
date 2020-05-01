!===============================================================================
!> Particle data type, extends a dll_node data type, so it can become part of
!> a particles_t data type
!===============================================================================
MODULE particle_mod
  USE iso_fortran_env, only: int8
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  USE patch_forces_mod
  USE lagrange_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: particle_t
    integer:: it=1, nt=0
    integer(8):: nstep=1
    real:: ds=1.0
    real:: mass=1.0
    real(8):: time=0d0, dtime=0d0
    integer(kind=int8), allocatable:: iit(:)
    real, allocatable:: v(:,:)
    real(8), allocatable:: r(:,:), t(:), m(:)
    real(8):: f(3)
    type(patch_forces_t), pointer:: patch_forces => null()
  contains
    procedure:: init
    procedure:: predict
    procedure:: timeslots
    procedure:: interpolate
  end type
  integer, save:: id=0
  integer, save:: order=2
  integer, save:: verbose=0
  integer, save:: nt=4
  type(particle_t), public:: particle
  !-----------------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Read parameters, and initialize a particle with a unique ID
!===============================================================================
SUBROUTINE init (self)
  class (particle_t):: self
  logical, save:: first_time=.true.
  integer:: iostat, it
  namelist /particle_params/ verbose, nt, order
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_t%init')
  if (first_time) then
    !$omp critical (particle_input_cr)
    if (first_time) then
      rewind (io_unit%input)
      read (io_unit%input, particle_params,iostat=iostat)
      write (io_unit%output, particle_params)
      patch_forces%verbose = verbose
      first_time = .false.
    end if
    !$omp end critical (particle_input_cr)
  end if
  !$omp atomic
  id = id+1
  if (verbose > 1) &
    print *, 'particle_t%init: id =', id
  self%id = id
  if (allocated(self%r)) then
    write (stderr,*) 'WARNING, particle_t%init: attempt to allocate twice'
  else
    allocate (self%r(3,nt), self%v(3,nt), self%m(nt), self%t(nt), self%iit(nt))
  end if
  self%r  = 0.0_8
  self%v  = 0.0_8
  self%f  = 0.0_8
  self%t  = 0.0_8
  self%it = 1
  self%nt = nt
  allocate (self%patch_forces)
  call self%patch_forces%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Estimate the position and velocity half a step forward in time
!===============================================================================
SUBROUTINE predict (self)
  class (particle_t):: self
  integer:: prv, now, new, i
  real(8):: x0, x(self%nt)
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_t%predict')
  !$omp atomic read
  now = self%it
  new = mod   (now  ,nt) + 1
  prv = modulo(now-2,nt) + 1
  self%r(:,new) = self%r(:,now)*1.5_8 - self%r(:,prv)*0.5_8
  self%v(:,new) = self%v(:,now)*1.5_8 - self%v(:,prv)*0.5_8
  self%t(  new) = self%t(  now)*1.5_8 - self%t(  prv)*0.5_8
!  do i=1,self%nt
!    x(i) = i
!  end do
!  x0 = self%nt + 0.5d0
!  self%r(:,new) = self%interpolate (x0, x, self%r)
  call trace%end()
END SUBROUTINE predict

!===============================================================================
!> Interpolate a particle position to a given point in time
!===============================================================================
FUNCTION interpolate (self, t0, yin) RESULT(f)
  class(particle_t):: self
  real(8):: yin(:,:), f(3)
  !.............................................................................
  integer:: i, j, i0, i1, it, jt, iit(self%nt), order
  real(8):: t0, t(self%nt), y(3,self%nt), w(self%nt)
  !-----------------------------------------------------------------------------
  ! Get the timestlot indices
  !-----------------------------------------------------------------------------
  iit = self%timeslots()
!print *,'iit',iit
  !-----------------------------------------------------------------------------
  ! Pick up time and position values with atomic read
  !-----------------------------------------------------------------------------
  do it=1,self%nt
    jt = iit(it)
    !$omp atomic read
    t(it) = self%t(jt)
    do j=1,3
      !$omp atomic read
      y(j,it) = yin(j,jt)
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Use Lagrange inter- and extrapolation
  !-----------------------------------------------------------------------------
  if (self%nstep == 1) then
    f = y(:,1)
  else 
    order = min(self%nstep-1,3)
    call lagrange%sequence_weights (t0, t, i0, i1, w, order=order)
    if (verbose > 2) &
      write (stdout,'(a,i3,2i4,f8.4,2x,9f8.4)') &
        'particle_t%interpolate: order, i0, i1, sum, w =', &
        order, i0, i1, sum(w(1:1+i1-i0)), w(1:1+i1-i0)
    f = 0d0
    do i=i0,i1
      do j=1,3
        f(j) = f(j) + w(i-i0+1)*y(j,i)
      end do
    end do
  end if
  if (verbose > 1) &
    write (stdout,*) 'particle_t%interpolate: id, f =', self%id, f
END FUNCTION interpolate

!===============================================================================
!> Return index array, based on only one atomic read
!===============================================================================
FUNCTION timeslots(self) RESULT (iit)
  class(particle_t):: self
  integer:: iit(self%nt)
  integer:: i, it, n
  !-----------------------------------------------------------------------------
  !$omp atomic read
  it = self%it
  do i=1,self%nt
    iit(i) = mod(i-1+it,self%nt) + 1
  end do
  !do i=1,self%nt-self%nstep
  !  iit(i) = 1
  !end do
END FUNCTION timeslots

END MODULE particle_mod
