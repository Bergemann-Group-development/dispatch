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
  type, public :: sinks_t
    type(dll_t):: dll
    integer:: center_star, version, nstar, szstar, sizeofstar, i, nstars, nmet, dp1
    logical:: do_translate
    real(4):: t, beta_create, dtold
    real(8):: t8, tflush, r_translate(3), v_translate(3), x_refine, y_refine, z_refine
    type(sink_t), pointer:: sink(:)
    real(8), allocatable:: sink_p(:,:)
    real(8), allocatable:: sink_v(:,:)
    real(8), allocatable:: sink_m(:)
  contains
    procedure:: init
    procedure:: init_nbors
    procedure:: update
  end type
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
    integer:: has_moved
    real, allocatable:: nbor_p(:,:)
    real, allocatable:: nbor_m(:)
    type(dll_t):: nbors
  end type
  integer :: verbose=0                  ! verbosity
  logical :: on=.false.                 ! on/off
  real    :: fraction=0.01              ! fraction of max force
  real    :: a_min = 1                  ! min acceleration
CONTAINS

!===============================================================================
!> Facilitate call to %read_snapshot().  All heavy lifting is done there
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  !.............................................................................
  integer:: sink_unit, nout, iostat, i, unit
  character(len=80):: fname, path
  class(sink_t), pointer:: sink, star
  class(dll_node_t), pointer:: dll_node
  namelist /sinks_params/ on, verbose, path, nout, fraction, a_min
  !.............................................................................
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
      call sink%nbors%init
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
    if (associated(star)) then
      write (unit,*) 'type(sink_t) :: star => center_star'
      write (unit,'(a,"["2(f14.9,","),f14.9,"]",2x,a,1i3.1)') &
        ' pos = ',star%position,', level = ',star%ilevel
    end if
    if (verbose >= 0) &
      write (stdout,'(a,2x,a,i6)') &
        ' file: '//trim(fname),'sinks:',self%dll%n
    call self%init_nbors
  else
    if (verbose >= 0) &
      write (unit,'(a,i6)')  ' sink_t%init: sinks off'
  end if
  call trace%end ()
END SUBROUTINE init

SUBROUTINE init_nbors (self)
  class(sinks_t):: self
  class(dll_node_t), pointer:: sink, other, nbor
  real(8):: start
  real:: r(3), ri, a
  integer:: i, n, j, m, m_tot
  integer, allocatable:: jj(:)
  !---------------------------------------------------------------------------
  ! Add nbor lists with most relevant other sinks
  !---------------------------------------------------------------------------
  call timer%tic (start)
  n = self%dll%n
  allocate (self%sink_p(n,3))
  allocate (self%sink_v(n,3))
  allocate (self%sink_m(n))
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
      self%sink_m(i) = sink%mass
      self%sink_p(i,:) = sink%position
      self%sink_v(i,:) = sink%velocity
    end select
    sink => sink%next
  end do
  !-----------------------------------------------------------------------------
  ! Find the number of nbors
  !-----------------------------------------------------------------------------
  m_tot = 0
  sink => self%dll%head
  do i=1,n
    m = 0
    do j=1,n
      if (j==i) cycle
      r = self%sink_p(i,:) - self%sink_p(j,:)
      r = modulo (r+0.5d0,1d0) - 0.5d0
      ri = 1./sqrt(r(1)**2 + r(2)**2 + r(3)**2)
      a = scaling%grav*self%sink_m(j)*ri**2
      if (a > a_min) then
        m = m+1
        jj(m) = j
      end if
    end do
    m_tot = m_tot + m
    select type (sink)
    class is (sink_t)
      allocate (sink%nbor_m(m))
      allocate (sink%nbor_p(m,3))
      do j=1,m
        sink%nbor_m(j) = self%sink_m(jj(j))
        sink%nbor_p(j,:) = self%sink_p(jj(j),:)
      end do
    end select
    sink => sink%next
  end do
  print *,m_tot
  if (verbose > 0) then
    write (stdout,*) m_tot/self%dll%n, 'nbors per sink'
    call timer%toc ('sinks_t%init_nbors', 1, start)
  end if
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
      do j=1,size(sink%nbor_m)
        r = sink%position - sink%nbor_p(j,:)
        r = modulo (r+0.5d0,1d0) - 0.5d0
        ri = 1./sqrt(r(1)**2 + r(2)**2 + r(3)**2)
        f = r*scaling%grav*sink%nbor_m(j)*ri**3
        self%sink_v(j,1) = self%sink_v(j,1) + a(1)*dt
        self%sink_v(j,3) = self%sink_v(j,2) + a(2)*dt
        self%sink_v(j,3) = self%sink_v(j,3) + a(3)*dt
        self%sink_p(j,1) = self%sink_p(j,1) + self%sink_v(j,1)*dt
        self%sink_p(j,2) = self%sink_p(j,2) + self%sink_v(j,2)*dt
        self%sink_p(j,3) = self%sink_p(j,3) + self%sink_v(j,3)*dt
      end do
    end select
    sink => sink%next
  end do
  call timer%toc ('sinks_t%update', self%dll%n, start)
  call trace%end ()
END SUBROUTINE update

END MODULE sinks_mod
