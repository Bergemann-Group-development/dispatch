!===============================================================================
!> Handle sink particle in RAMSES format
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE omp_timer_mod
  USE task_mod
  USE scaling_mod
  USE units_mod
  USE dll_mod
  implicit none
  private
  integer, parameter :: p8=kind(1.0D0)  ! double precision
  integer, parameter :: dp=kind(1.0D0)  ! double precision
  !-----------------------------------------------------------------------------
  ! Particle cache, for vectorized operations
  !-----------------------------------------------------------------------------
  type:: cache_t
    real, allocatable:: p(:,:)          ! positions
    real, allocatable:: v(:,:)          ! velocities
    real, allocatable:: m(:)            ! masses
  end type
  !-----------------------------------------------------------------------------
  ! Type holding data for all sinks
  !-----------------------------------------------------------------------------
  type, public :: sinks_t
    type(dll_t):: dll           ! doubly linked list
    integer:: version           ! I/O version
    integer:: center_star       ! index of center star
    integer:: nstar             ! number of sinks
    integer:: szstar            ! ?
    integer:: sizeofstar        ! ?
    integer:: i                 ! ?
    integer:: nstars            ! ?
    integer:: nmet              ! number of heavy elements
    integer:: dp1               ! ?
    logical:: do_translate      ! logical switch
    real(4):: t                 ! single precision time
    real(8):: t8                ! double precisions time
    real(8):: beta_create       ! plasma beta at creation
    real(8):: dtold             ! previous dt
    real(8):: tflush            ! flush time
    real(8):: r_translate(3)    ! translation position
    real(8):: v_translate(3)    ! translation velocity
    real(8):: x_refine          ! refinement center
    real(8):: y_refine          ! refinement center
    real(8):: z_refine          ! refinement center
    type(cache_t):: cache       ! cache for vectorized operations
  contains
    procedure:: init
    procedure:: update
    procedure, private:: init_nbors
  end type
  !-----------------------------------------------------------------------------
  ! Type holding data for one RAMSES sink particle
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: sink_t
    !SEQUENCE
    real(p8):: position(3)              ! positions
    real(dp):: momentum(3)              ! momenta
    real(dp):: velocity(3)              ! velocities
    real(dp):: dmom(3)                  ! momenta
    real(dp):: mass                     ! mass
    real(dp):: dmass                    ! mass
    real(dp):: potential                ! gravitational potential
    real(dp):: density                  ! neigborhood density
    real(dp):: t_create                 ! creation time
    real(dp):: t_explode                ! explosion time
    real(dp), allocatable:: metal(:)    ! metal density
    real(dp), allocatable:: dmetal(:)   ! metal density
    real(dp):: beta_create              ! Beta at creation
    real(dp):: dtold                    ! previous dt
    integer:: gid                       ! pointer to grid cell
    integer:: father                    ! pointer to father cell
    integer:: ilevel                    ! level of AMR hierachy where star lives
    integer:: old_ilevel                ! level of hierachy star used to have, if changed
    integer:: owner                     ! owner thread
    integer:: has_moved                 ! indicator
    type(cache_t):: cache               ! for vectorized operations
  end type
  logical :: on=.false.                 ! on/off
  integer :: verbose=0                  ! verbosity
  real    :: a_min = 10                 ! min acceleration
CONTAINS

!===============================================================================
!> Initialize by reading in a snapshot of sink particles
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  !.............................................................................
  integer:: sink_unit, nout, iostat, i, unit
  character(len=80):: fname, path
  class(sink_t), pointer:: sink, star
  class(dll_node_t), pointer:: dll_node
  namelist /sinks_params/ on, verbose, path, nout, a_min
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
  call io%header ('sinks_t%init:')
  rewind (io%input)
  read (io%input, sinks_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, sinks_params)
  end if
  call io%namelist_warning ('sinks_params', iostat)
  !-----------------------------------------------------------------------------
  ! If on, read sinks
  !-----------------------------------------------------------------------------
  if (on) then
    call self%dll%init
    !---------------------------------------------------------------------------
    ! Read header part
    !---------------------------------------------------------------------------
    if (verbose > 1) then
      unit = stdout
    else
      unit = stdlog
    end if
    write (fname,'(a,i5.5,a)') trim(path)//'/output_', nout, '/stars_output.dat'
    write (unit,*) fname
    open  (newunit=sink_unit, file=fname, form='unformatted', status='unknown')
    read  (sink_unit) self%version, self%nstars, self%t, self%szstar, self%sizeofstar
    write (unit,*) 'versions =', self%version
    read  (sink_unit) self%dp1, self%nmet, self%t8
    write (unit,*) 'dp, nmetal =', self%dp1, self%nmet
    write (unit,*) 't8 =', self%t8
    write (unit,*) 'nstars =', self%nstars
    read  (sink_unit) self%do_translate, self%center_star, self%r_translate, self%v_translate,  &
                     self%x_refine, self%y_refine, self%z_refine
    write (unit,*) 'center_star', self%center_star
    !---------------------------------------------------------------------------
    ! Read sinks
    !---------------------------------------------------------------------------
    nullify (star)
    do i=1,self%nstars
      allocate (sink)
      allocate (sink%metal(self%nmet))
      allocate (sink%dmetal(self%nmet))
      dll_node => sink
      read (sink_unit) sink%position, sink%momentum , &
                       sink%dmom, sink%mass, &
                       sink%dmass, sink%potential, sink%density, sink%t_create, &
                       sink%t_explode, sink%metal(1:self%nmet), sink%dmetal(1:self%nmet), &
                       sink%gid, sink%father, sink%ilevel, sink%old_ilevel, &
                       sink%owner, sink%has_moved, sink%beta_create, sink%dtold
      if (verbose > 1) &
        write (*,'(3f10.6,1pe11.2,i12)') sink%position, sink%mass, sink%ilevel
      if (sink%mass > 0.0) then
        call self%dll%append (dll_node)
        if (i == self%center_star) then
          star => sink
          write (stderr,*) 'found center star'
        end if
      else
        deallocate (sink%metal, sink%dmetal)
        deallocate (sink)
      end if
    enddo
    close(sink_unit)
    !...........................................................................
    if (verbose >= 0) &
      write (stdout,'(a,2x,a,i6)') &
        ' file: '//trim(fname),'sinks:',self%dll%n
    !---------------------------------------------------------------------------
    ! Construct initial nbor tables
    !---------------------------------------------------------------------------
    call self%init_nbors
  else
    if (verbose >= 0) &
      write (unit,'(a,i6)')  ' sink_t%init: sinks off'
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Construct tables with relevant neighbor sinks
!===============================================================================
SUBROUTINE init_nbors (self)
  class(sinks_t):: self
  class(dll_node_t), pointer:: sink, other, nbor
  real(8):: start
  real:: r(3), ri, a
  integer:: i, n, j, m, m_tot
  integer, allocatable:: jj(:)
  !---------------------------------------------------------------------------
  call trace%begin ('sinks_t%init_nbors')
  call timer%tic (start)
  n = self%dll%n
  allocate (self%cache%p(n,3))
  allocate (self%cache%v(n,3))
  allocate (self%cache%m(n))
  allocate (jj(n))
  !-----------------------------------------------------------------------------
  ! Cache positions and mass
  !-----------------------------------------------------------------------------
  i = 0
  sink => self%dll%head
  do while (associated(sink))
    select type (sink)
    class is (sink_t)
      i = i+1
      self%cache%m(i) = sink%mass
      self%cache%p(i,:) = sink%position
      self%cache%v(i,:) = sink%velocity
    end select
    sink => sink%next
  end do
  !-----------------------------------------------------------------------------
  ! Find the number of nbors and their indices
  !-----------------------------------------------------------------------------
  m_tot = 0
  sink => self%dll%head
  do i=1,n
    m = 0
    do j=1,n
      if (j==i) cycle
      r = self%cache%p(i,:) - self%cache%p(j,:)
      r = modulo (r+0.5d0,1d0) - 0.5d0
      ri = 1./sqrt(r(1)**2 + r(2)**2 + r(3)**2)
      a = scaling%grav*self%cache%m(j)*ri**2
      if (a > a_min) then
        m = m+1
        jj(m) = j
      end if
    end do
    m_tot = m_tot + m
    !---------------------------------------------------------------------------
    ! Allocate cache for nbors
    !---------------------------------------------------------------------------
    select type (sink)
    class is (sink_t)
      allocate (sink%cache%m(m))
      allocate (sink%cache%p(m,3))
      do j=1,m
        sink%cache%m(j) = self%cache%m(jj(j))
        sink%cache%p(j,:) = self%cache%p(jj(j),:)
      end do
    end select
    sink => sink%next
  end do
  !.............................................................................
  if (verbose > 0) then
    write (stdout,*) m_tot/self%dll%n, 'nbors per sink'
    call timer%toc ('sinks_t%init_nbors', 1, start)
  end if
  call trace%end ()
END SUBROUTINE init_nbors

!===============================================================================
!> Update all positions -- for now simple 1st order, constant dt
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  class(dll_node_t), pointer:: sink, other, nbor
  real(8):: start, mass, m
  real:: f(3), r(3), ri, a(3), dt
  integer:: j
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  call timer%tic (start)
  dt = 1e-10
  sink => self%dll%head
  do while (associated(sink))
    select type (sink)
    class is (sink_t)
      !-------------------------------------------------------------------------
      ! Run through all nbors
      !-------------------------------------------------------------------------
      do j=1,size(sink%cache%m)
        r = sink%position - sink%cache%p(j,:)
        r = modulo (r+0.5d0,1d0) - 0.5d0
        ri = 1./sqrt(r(1)**2 + r(2)**2 + r(3)**2)
        f = r*scaling%grav*sink%cache%m(j)*ri**3
        self%cache%v(j,1) = self%cache%v(j,1) + a(1)*dt
        self%cache%v(j,3) = self%cache%v(j,2) + a(2)*dt
        self%cache%v(j,3) = self%cache%v(j,3) + a(3)*dt
        self%cache%p(j,1) = self%cache%p(j,1) + self%cache%v(j,1)*dt
        self%cache%p(j,2) = self%cache%p(j,2) + self%cache%v(j,2)*dt
        self%cache%p(j,3) = self%cache%p(j,3) + self%cache%v(j,3)*dt
      end do
    end select
    sink => sink%next
  end do
  call timer%toc ('sinks_t%update', self%dll%n, start)
  call trace%end ()
END SUBROUTINE update

END MODULE sinks_mod
