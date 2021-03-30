!===============================================================================
!> Handle sink particle in RAMSES format
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
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
  contains
    procedure:: init
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
  end type
  integer, save:: verbose=0             ! verbosity
  logical, save:: on=.false.            ! on/off
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
  namelist /sinks_params/ on, verbose, path, nout
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
    if (verbose > 0) then
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
      call self%dll%append (dll_node)
      if (i == self%center_star) then
        star => sink
        write (stderr,*) 'found center star'
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
  else
    if (verbose >= 0) &
      write (unit,'(a,i6)')  ' sink_t%init: sinks off'
  end if
  call trace%end ()
END SUBROUTINE init

END MODULE sinks_mod
