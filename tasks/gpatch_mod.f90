!===============================================================================
!> The gpath_t layer now essentially only handles restarts
!===============================================================================
MODULE gpatch_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE omp_lock_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE mesh_mod
  USE units_mod
  USE timer_mod
  USE omp_timer_mod
  USE bits_mod
  ! USE eos_mod
  USE data_io_mod
  USE validate_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE scratch_mod
  USE boundaries_mod
  ! USE stagger_mod
  USE connect_mod
  USE shared_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(connect_t):: connect
    real(8):: d0=0d0
    logical:: pre_updated=.false.
    class(list_t),       pointer:: task_list => null()
    class(boundaries_t), pointer:: boundaries => null()
  contains
    procedure, nopass:: cast2gpatch
    procedure:: pre_init
    procedure:: init
    procedure:: init_bdries
    procedure:: init_task_list
    procedure:: add_to_task_list
    procedure:: pre_update
    procedure:: post_update
    procedure:: input
    procedure:: output
    procedure:: counter_update
    procedure:: dnload
    procedure:: courant_condition
    procedure:: to_thermal
    ! procedure:: to_total_e
    procedure:: clone_mem_accounting
    procedure:: rotate
    procedure:: info                    ! print info to stdout
    procedure, private :: fminval8sph
    procedure, private :: fminval4sph
    procedure, private :: fmaxval8sph
    procedure, private :: fmaxval4sph
    generic, public :: fmaxvalsph => fmaxval4sph, fmaxval8sph
    generic, public :: fminvalsph => fminval4sph, fminval8sph
    procedure:: output_metadata
    procedure:: input_metadata
    procedure:: reset
    procedure:: initial_load_balance
    procedure:: debug_dump
  end type
  real(8), save :: restart_time = 0d0
  logical, save :: detailed_timer=.false.
  logical, save :: extrapolate_gz=.true.            ! extrapolate guard zones?
  integer, save :: verbose=0
  integer, save :: order=1
  type(gpatch_t), public:: gpatch
CONTAINS

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2gpatch (task) RESULT(gpatch)
  class(task_t), pointer:: task
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  select type (task)
  class is (gpatch_t)
  gpatch => task
  class default
  nullify(gpatch)
  call io%abort ('gpatch_t%cast: failed to cast a task to patch_t')
  end select
END FUNCTION cast2gpatch

!===============================================================================
!> Initialize boundary bits
!===============================================================================
SUBROUTINE pre_init (self)
  class(gpatch_t), target:: self
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%pre_init')
  !call eos%init()
  call self%patch_t%pre_init
  if (.not. associated (self%boundaries)) then
    allocate (self%boundaries)
  end if
  call self%boundaries%pre_init (self)
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize boundary bits
!===============================================================================
SUBROUTINE init_bdries (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%init_bdries')
  if (.not.associated(self%boundaries)) then
    write (stdout,'(a)') &
      'gpatch_t%init_bdries ERROR: boundaries not associated. Check pre_init sequences!'
  end if
  call self%boundaries%init ()
  call trace%end()
END SUBROUTINE init_bdries

!===============================================================================
!> Make a copy of the task list pointer
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(gpatch_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('gpatch_t%init_task_list')
  self%task_list => task_list
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Void subroutine, to be intercepted by extras_mod, if extras underlings 
!> need to be added to the task list
!===============================================================================
SUBROUTINE add_to_task_list (self, task_list)
  class(gpatch_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('gpatch_t%init_task_list')                                     ! sinkp
  call trace%end ()
END SUBROUTINE add_to_task_list

!===============================================================================
!> Post-initialize. For all solvers, copy the %it slot to the %new slot.
!===============================================================================
SUBROUTINE pre_update (self)
  class(gpatch_t), target:: self
  integer:: i1, i2, i3, i4
  type(lock_t), pointer:: lock
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%pre_updated) &
    return
  call trace%begin ('gpatch_t%pre_update', itimer=itimer)
  call self%patch_t%pre_update ()
  !-----------------------------------------------------------------------------
  ! Initialize scratch arrays to zero
  !-----------------------------------------------------------------------------
  call scratch%pre_update
  !-----------------------------------------------------------------------------
  ! Mark task pre_updated (cleared in post_update)
  !-----------------------------------------------------------------------------
  self%pre_updated = .true.
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
SUBROUTINE post_update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  self%pre_updated = .false.
END SUBROUTINE post_update

!===============================================================================
!> Interface routine to the data_io input procedure, signalling if the patch
!> could be read, or not
!===============================================================================
SUBROUTINE input (self, ok)
  class(gpatch_t):: self
  logical:: ok
  !-----------------------------------------------------------------------------
  ok = .false.
  call data_io%input (self, ok)
END SUBROUTINE input

!===============================================================================
!> Interface routine to the data_io output procedure
!===============================================================================
SUBROUTINE output (self)
  class(gpatch_t):: self
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%output')
  call data_io%output (self)
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Initialize patch_t and data_io
!===============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  character(len=129):: id= &
  '$Id$'
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%frozen)) &
    return
  call trace%begin ('gpatch_t%init', itimer=itimer)
  call io_unit%print_id (id)
  call self%patch_t%init
  call self%init_bdries ()
  call data_io%init (self)
  call trace%end(itimer)
END SUBROUTINE init

!===============================================================================
!> Increment the cell udpate counter.  This may be overloaded in more complex
!> solvers, which prefer a different way of counting.
!===============================================================================
SUBROUTINE counter_update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%counter_update')
  !$omp atomic
  timer%n_update = timer%n_update + product(self%n)
  call trace%end()
END SUBROUTINE counter_update

!===============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!===============================================================================
SUBROUTINE dnload (self, only)
  class(gpatch_t):: self
  integer, optional:: only
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%dnload', itimer=itimer)
  if (self%is_clear (bits%frozen)) then
    write (stderr,*) self%id, 'called gpatch_t%dnload, time=', self%time
    call io%warning ('gpatch_t%download is deprecated')
    call self%link%update
  end if
  call trace%end (itimer)
END SUBROUTINE dnload

!===============================================================================
!> Intercept to make validate%check call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(gpatch_t):: self
  logical, optional:: detailed_timer
  type(link_t), pointer:: nbor
  real(8):: dtime, nbdtime
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%courant_condition')
  call validate%check (self%link, self%u_max, 'courant')
  call self%patch_t%courant_condition (detailed_timer)
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Compute or assign thermal energy
!===============================================================================
SUBROUTINE to_thermal (self, it, eth, kind)
  class(gpatch_t):: self
  integer:: it
  real, contiguous:: eth(:,:,:)
  character(len=*), optional:: kind
  !.............................................................................
  integer:: i1, i2, i3
  real, pointer, contiguous, dimension(:,:,:):: d, e, px, py, pz, bx, by, bz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%to_thermal', itimer=itimer)
  e  => self%mem(:,:,:,self%idx%e ,it,1)
  if (present(kind) .and. kind(1:6)=='ramses') then
    d  => self%mem(:,:,:,self%idx%d ,it,1)
    px => self%mem(:,:,:,self%idx%px,it,1)
    py => self%mem(:,:,:,self%idx%py,it,1)
    pz => self%mem(:,:,:,self%idx%pz,it,1)
    eth = e - 0.5*(px**2+py**2+pz**2)/d
    if (self%nv > 5) then
      bx => self%mem(:,:,:,self%idx%bx,it,1)
      by => self%mem(:,:,:,self%idx%by,it,1)
      bz => self%mem(:,:,:,self%idx%bz,it,1)
      do i3=1,self%gn(3)-1
      do i2=1,self%gn(2)-1
      do i1=1,self%gn(1)-1
        eth(i1,i2,i3) = &
        eth(i1,i2,i3) - 0.125*((bx(i1,i2,i3)+bx(i1+1,i2,i3))**2 + &
                               (by(i1,i2,i3)+by(i1,i2+1,i3))**2 + &
                               (bz(i1,i2,i3)+bz(i1,i2,i3+1))**2)
      end do
      end do
      end do
    end if
  else
    eth = e
  end if
  call trace%end (itimer)
END SUBROUTINE to_thermal

! !===============================================================================
! !> Convert to total energy
! !===============================================================================
! SUBROUTINE to_total_e (self, it)
!   class(gpatch_t):: self
!   integer:: it
!   real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by ,bz
!   !-----------------------------------------------------------------------------
!   if (self%kind(1:6)=='ramses') then
!     d  => self%mem(:,:,:,self%idx%d ,it,1)
!     e  => self%mem(:,:,:,self%idx%e ,it,1)
!     px => self%mem(:,:,:,self%idx%px,it,1)
!     py => self%mem(:,:,:,self%idx%py,it,1)
!     pz => self%mem(:,:,:,self%idx%pz,it,1)
!     e = e + 0.5*(px**2+py**2+pz**2)/d
!     if (self%idx%bx > 0) then
!       bx => self%mem(:,:,:,self%idx%bx,it,1)
!       by => self%mem(:,:,:,self%idx%by,it,1)
!       bz => self%mem(:,:,:,self%idx%bz,it,1)
!       e = e + 0.5*(xup(bx)**2+yup(by)**2+zup(bz)**2)
!     end if
!   end if
! END SUBROUTINE to_total_e

!===============================================================================
!> Stub -- left for compatibility only
!===============================================================================
SUBROUTINE clone_mem_accounting (self)
  class(gpatch_t):: self
END SUBROUTINE clone_mem_accounting

!===============================================================================
FUNCTION fmaxval4sph (self, f, outer) RESULT (fmaxval)
  real:: fmaxval
  class(gpatch_t):: self
  real, dimension(:,:,:), intent(in):: f
  logical, optional:: outer
  logical:: outer_l
  !-----------------------------------------------------------------------------
  if (present(outer)) then
    outer_l = outer
  else
    outer_l = .false.
  end if
  if (self%boundaries%is_set (bits%spherical)) then
    fmaxval = maxval(f, mask=self%boundaries%mask(self%mesh))
  else if (outer_l) then
    fmaxval = maxval(f(self%mesh(1)%lo:self%mesh(1)%uo, &
                       self%mesh(2)%lo:self%mesh(2)%uo, &
                       self%mesh(3)%lo:self%mesh(3)%uo))
  else
    fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fmaxval4sph
!===============================================================================
FUNCTION fminval4sph (self, f) RESULT(fminval)
  real:: fminval
  class(gpatch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%spherical)) then
    fminval = minval(f, mask=self%boundaries%mask(self%mesh))
  else
    fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fminval4sph
!===============================================================================
FUNCTION fmaxval8sph (self, f, outer) RESULT (fmaxval)
  real(8):: fmaxval
  class(gpatch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  logical, optional:: outer
  logical:: outer_l
  !-----------------------------------------------------------------------------
  if (present(outer)) then
    outer_l = outer
  else
    outer_l = .false.
  end if
  if (self%boundaries%is_set (bits%spherical)) then
    fmaxval = maxval(f, mask=self%boundaries%mask(self%mesh))
  else if (outer_l) then
    fmaxval = maxval(f(self%mesh(1)%lo:self%mesh(1)%uo, &
                       self%mesh(2)%lo:self%mesh(2)%uo, &
                       self%mesh(3)%lo:self%mesh(3)%uo))
  else
    fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fmaxval8sph
!===============================================================================
FUNCTION fminval8sph (self, f) RESULT(fminval)
  real(8):: fminval
  class(gpatch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%spherical)) then
    fminval = minval(f, mask=self%boundaries%mask(self%mesh))
  else
    fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fminval8sph

!===============================================================================
!> Print info to stdout.  If
!>   verbose==0:        print only id=1 to stdout and io_unit%log
!>   verbose==1:        print only id=1 to stdout but all to io_unit%log
!>   verbose==2:        print iv=1 to stdout and io_unit%log for all id
!>   verbose==3:        print iv=1,nv to stdout and io_unit%log for all id
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(gpatch_t)    :: self
  integer            :: nq, ntask
  character(len=64)  :: experiment_name
  optional           :: nq, ntask, experiment_name
  !.............................................................................
  type(lock_t), save :: lock
  integer, parameter :: max_lines=50
  integer, save      :: itimer=0
  integer, save      :: counts(6)=0, id_prv=0
  integer            :: counts_l(6), i
  integer            :: it0, it1, it, iv, nv, iprint_l
  real, save         :: d_max=0.0, d_min=huge(1.)
  real, save         :: dt_min=huge(1.), dt_max=0.0
  real, save         :: courant_info=0.0
  real               :: fmin, fmax, dmin, dmax, dtime
  real, pointer      :: df(:,:,:)
  real(8)            :: wc, time, print_next
  logical            :: print_it
  character(len=120) :: fmt
  logical, save      :: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%info', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Double-if to initialize the shared lock since this is entered often
  !-----------------------------------------------------------------------------
  if (first_time) then
    call io%lock%set ('gpatch_t%info')
    if (first_time) then
      call lock%init
      first_time = .false.
    end if
    call io%lock%unset ('gpatch_t%info')
  end if
  !-----------------------------------------------------------------------------
  ! Keep track of time and timer level
  !-----------------------------------------------------------------------------
  if (self%level >= timer%time_level) then
    !$omp atomic update
    timer%time = max(timer%time,self%time)
    if (self%level > timer%time_level) then
      !$omp atomic write
      timer%time_level = self%level
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Low noise case for log-file based validation
  !-----------------------------------------------------------------------------
  if (io_unit%do_validate .and. associated(self%mem)) then
    it0 = max(self%nt-2, 1)
    it1 = max(self%nt-1, 1)
    it0 = self%iit(it0)
    it1 = self%iit(it1)
    if (self%nw==1) then
      allocate (df(size(self%mem,1),size(self%mem,2),size(self%mem,3)))
    end if
    do iv=1,self%nv
      if (self%nw==1) then
        df = self%mem(:,:,:,iv,it1,1) &
           - self%mem(:,:,:,iv,it0,1)
        fmin = self%fminval(df)
        fmax = self%fmaxval(df)
      else
        fmin = self%fminval(self%mem(:,:,:,iv,it0,self%nw))
        fmax = self%fmaxval(self%mem(:,:,:,iv,it0,self%nw))
      end if
      write(io%output,2) &
      'upd',self%id, self%level, iv, self%istep, self%time-self%dtime, &
      self%dtime, courant_info, fmin, fmax, &
      self%fminval(self%mem(:,:,:,iv,it0,1)), &
      self%fmaxval(self%mem(:,:,:,iv,it0,1)), &
      real(timer%n_update), ntask, nq, self%is_set(bits%boundary), self%is_set(bits%virtual)
    2 format(a,i8,2i4,i3,1p,e14.5,e12.4,0p,f7.2,2x,1p,2e12.4,1x,2e12.4,2x,e10.2,2i5,l3,l1,o12)
    end do
    if (self%nw==1) then
      deallocate (df)
    end if
    call trace%end (itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Normal case conditional output
  !-----------------------------------------------------------------------------
  if ((io%verbose >= 0) .and. self%is_clear(bits%no_density)) then
    !---------------------------------------------------------------------------
    ! After coming back to the saved ID, reset it, so we stop searching
    !---------------------------------------------------------------------------
    if (self%id == id_prv) then
      !$omp atomic write
      id_prv = 0
    end if
    !---------------------------------------------------------------------------
    ! If print_every==0, then check if print_time controls cadence, using
    ! an if-lock-if sequence to ensure only one thread prints and updates
    !---------------------------------------------------------------------------
    print_it = .false.
    if (io%print_seconds > 0.0) then
      wc = wallclock()
      if (wc >= io%print_next) then
        call lock%set ('print')
        if (wc >= io%print_next) then
          print_it = .true.
          print_next = (nint(wc/io%print_seconds)+1)*io%print_seconds
          io%print_next = print_next
        end if
        call lock%unset ('print')
      end if
    else if (io%print_every > 0) then
      !$omp atomic read
      iprint_l = io%iprint
      if (iprint_l >= io%print_every .and. io%print_every > 0) then
        !$omp atomic write
        io%iprint = 1
        print_it = .true.
      else
        !$omp atomic update
        io%iprint = io%iprint+1
        print_it = .false.
      end if
    !---------------------------------------------------------------------------
    ! If print_every==0, then check if print_time controls cadence, using
    ! an if-lock-if sequence to ensure only one thread prints and updates
    !---------------------------------------------------------------------------
    else if (io%print_time > 0.0) then
      if (io%print_every == 0 .and. self%time >= io%print_next) then
        call lock%set ('print')
        if (self%time >= io%print_next) then
          print_it = .true.
          print_next = (nint(self%time/io%print_time)+1)*io%print_time
          io%print_next = print_next
        end if
        call lock%unset ('print')
      end if
    !---------------------------------------------------------------------------
    ! Otherwise, fall back on printing the task update with the smallest dt
    !---------------------------------------------------------------------------
    else if (self%id==io%id) then
      print_it = .true.
    end if
    !---------------------------------------------------------------------------
    ! Accumulate extrema
    !---------------------------------------------------------------------------
    !$omp atomic update
    courant_info = max (self%courant_info, courant_info)
    dmin = self%fminval(self%idx%d)
    dmax = self%fmaxval(self%idx%d)
    dtime = self%dtime
    !$omp atomic update
    d_min = min (dmin, d_min)
    !$omp atomic update
    d_max = max (dmax, d_max)
    !$omp atomic update
    dt_min = min(dt_min,dtime)
    !$omp atomic update
    dt_max = max(dt_max,dtime)
    !---------------------------------------------------------------------------
    ! Print info for selected task, or tracked task
    !---------------------------------------------------------------------------
    print_it = print_it .or. self%track
    if (print_it) then
      nv = 1
      if (io%verbose>2) nv = self%nv
      !-------------------------------------------------------------------------
      counts_l = &
        [mpi_mesg%n_check, mpi_mesg%n_ready, mpi_mesg%n_update, &
         mpi_mesg%n_send , mpi_mesg%n_recv , mpi_mesg%n_unpk  ]
      do i=1,size(counts)
        !$omp atomic update
        counts(i) = counts(i) + counts_l(i)
      end do
      !$omp atomic update
      timer%n_lines = timer%n_lines-1
      if (timer%n_lines<=0) then
        !$omp critical (info_cr2)
        !$omp atomic write
        timer%n_lines = max_lines
        !$omp end atomic
        write(io%output,'(/a)') &
        '         id lvl thr        time       dtime   cinf      ' &
        //'dt_min    dt_max      rho_min    rho_max     cellupd  ntsk   nq ' &
        //' chk  rdy   upd  snd  rcv  unp  snq BV'
        !$omp end critical (info_cr2)
      end if
      !-----------------------------------------------------------------------
      ! time-dtime is the task time when update started; it is guaranteed to
      ! be increasing monotomically, since the ready queue is time sorted.
      ! To allow using the column for plotting, we show this time when all
      ! times are printed.
      !-----------------------------------------------------------------------
      time = self%time - self%dtime
      do iv=1,nv
        call write_info (io%output)
      end do
      flush (io%output)
      !$omp atomic write
      id_prv = io%id
      counts(:) = 0
      !-----------------------------------------------------------------------
      ! Reset the message counters after every timestep, if relevant
      !-----------------------------------------------------------------------
      mpi_mesg%n_check=0; mpi_mesg%n_ready=0; mpi_mesg%n_update=0
      mpi_mesg%n_send=0 ; mpi_mesg%n_recv=0 ; mpi_mesg%n_unpk=0
      mpi_mesg%n_delay=0
    end if
    !---------------------------------------------------------------------------
    ! Keep track of the shortest times step between updates of the io%id local
    ! task.  Use this task in patch%info, rather than the one with ip==1.  As
    ! long as shorter dtime are being found, io%id will not match the next self%id.
    ! But if the same io%id survices one update cycle, it will match.
    !---------------------------------------------------------------------------
    if (self%dtime < io%dtime .and. id_prv /= 0) then
      !$omp atomic write
      io%id = self%id
      !$omp atomic write
      io%dtime = self%dtime
    end if
  end if
  call trace%end (itimer)
contains
  !=============================================================================
  subroutine write_info (unit)
    integer:: unit
    character(64):: fmt_end
    !---------------------------------------------------------------------------
    fmt_end = '0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e9.2,2i6,2i5,i6,4i5,l2,l1)'
    if (self%time < 1d4 .and. self%time > 1e-4) then
      fmt = '(a,i8,2i4,f12.6,1p,e12.4,'//fmt_end
    else
      fmt = '(a,i8,2i4,1p,2e12.5,     '//fmt_end
    end if
    write(unit,fmt) &
      'upd',self%id, self%level, omp%thread, time, &
      self%dtime, courant_info, &
      dt_min, dt_max, d_min, d_max, &
      real(timer%n_update), ntask, nq, counts, mpi_mesg%nq_send,  &
      self%is_set(bits%boundary), self%is_set(bits%virtual)
    !---------------------------------------------------------------------------
    ! Reset extrema info
    !---------------------------------------------------------------------------
    !$omp atomic write
    courant_info = 0.0
    !$omp atomic write
    d_min = huge(1.)
    !$omp atomic write
    d_max = 0.0
    !$omp atomic write
    dt_min = huge(1.)
    !$omp atomic write
    dt_max = 0.0
  end subroutine write_info
END SUBROUTINE info

!===============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!===============================================================================
SUBROUTINE rotate (self)
  class(gpatch_t):: self
  integer:: iv, nv, new
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation
  !-----------------------------------------------------------------------------
  if (self%rotated) return
  call trace%begin ('gpatch_t%rotate', itimer=itimer)
  if (self%is_periodic()) then
    call self%periodic_grid                     ! periodic root grid
  else if (extrapolate_gz) then
    call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  !-----------------------------------------------------------------------------
  ! One line summary
  !-----------------------------------------------------------------------------
  if (io%verbose>1) then
    write (io_unit%log,'(a,i6,2i3,2e13.6,2x,a,1p,6e10.2)') 'update: id,j,k,l,n,dt,t', &
      self%id, self%it, self%new, self%dtime, self%time, 'minmax:', &
      minval(self%mem(:,:,:,1,self%it ,self%nw)), &
      maxval(self%mem(:,:,:,1,self%it ,self%nw)), &
      minval(self%mem(:,:,:,1,self%new,1)), &
      maxval(self%mem(:,:,:,1,self%new,1))
  end if
  !-----------------------------------------------------------------------------
  ! Update patch position for non-zero patch velocities;
  ! This happens *before* `rotate` because `dtime` stores the time step just taken.
  !-----------------------------------------------------------------------------
  call self%update_position
  !-----------------------------------------------------------------------------
  ! Task time slot rotation
  !-----------------------------------------------------------------------------
  call self%task_t%rotate
  call trace%end (itimer)
END SUBROUTINE rotate

!===============================================================================
!> Class interface to output of metadata
!===============================================================================
SUBROUTINE output_metadata (self, unit, first_time)
  class(gpatch_t):: self
  integer:: unit
  logical, optional:: first_time
  !-----------------------------------------------------------------------------
  call data_io%output_nml (self)
END SUBROUTINE output_metadata

!===============================================================================
!===============================================================================
SUBROUTINE input_metadata (self, unit)
  class(gpatch_t):: self
  integer:: unit
  !-----------------------------------------------------------------------------
  call data_io%input_nml (self, unit)
END SUBROUTINE input_metadata

!===============================================================================
!> Make sure pointers are nullified
!===============================================================================
SUBROUTINE reset (self)
  class (gpatch_t):: self
  !.............................................................................
  call trace%begin ('gpatch_t%reset')
  call self%patch_t%reset ()
  nullify (self%task_list)
  nullify (self%boundaries)
  self%d0 = 0.0
  self%pre_updated = .false.
  call trace%end ()
END SUBROUTINE reset

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE initial_load_balance (self, task_list)
  class(gpatch_t):: self
  class(list_t), pointer:: task_list
  !----------------------------------------------------------------------------
  call trace_begin('gpatch_t%initial_load_balance')
  call trace_end
END SUBROUTINE initial_load_balance

!===============================================================================
!> Write a debug_dump file with binary Fortran data, for either
!>       default: only the patch in the argument list
!>   all==.true.: in addition, all patches in the task list
!> nbors==.true.: in addition, all patches in the nbor list
!===============================================================================
SUBROUTINE debug_dump (self, label, nbors, all)
  class(gpatch_t):: self
  character(*):: label
  logical, optional:: nbors, all
  !.............................................................................
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  integer:: sentinel
  integer, save:: u=0
  !-----------------------------------------------------------------------------
  ! Dump patch data to unformatted Fortran file
  !-----------------------------------------------------------------------------
  call io%lock%set ('debug_dump')
  if (u /= 0) return
  write (stderr,*) mpi%rank, omp%thread, ' dumping patch', self%id 
  flush (stderr)
  sentinel = 0
  if (shared%option(nbors)) then
    sentinel = 1
    !---------------------------------------------------------------------------
    ! Dump nbor patches as well
    !---------------------------------------------------------------------------
    write (stderr,*) 'dumping all nbor patches as well'
    flush (stderr)
    call dump_patch (self)
    nbor => self%link%nbor
    do while (associated(nbor))
      task => nbor%task
      select type (task)
      class is (patch_t)
        call dump_patch (task)
      end select
      nbor => nbor%next
    end do
  else if (shared%option(all)) then
    !---------------------------------------------------------------------------
    ! Dump all other patches as well
    !---------------------------------------------------------------------------
    sentinel = 2
    write (stderr,*) 'dumping all other patches as well'
    flush (stderr)
    call dump_patch (self)
    link => self%task_list%head
    do while (associated(link))
      task => link%task
      if (task%id /= self%id) then
        select type (task)
        class is (patch_t)
          call dump_patch (task)
        end select
      end if
      link => link%next
    end do
  else
    call dump_patch (self)    
  end if
  write (stderr,*) 'dump done'
  flush (stderr)
  call io%lock%unset ('debug_dump')
contains
  !=============================================================================
  !> Dump a singel patch
  !=============================================================================
  subroutine dump_patch (p)
    class(patch_t):: p
    character(len=128):: f
    !---------------------------------------------------------------------------
    if (u==0) then
      write (f,'(a,"/debug_data.dat")') trim(io%run_dir)
      write (stderr,*) 'file: '//trim(f)
      open (newunit=u, file=trim(f), form='unformatted', status='unknown')
    end if
    write (u) sentinel
    write (u) p%id, p%it, p%gn, p%n, p%nv, p%nt, p%nw, mpi%rank, omp%thread
    write (u) p%time, p%dtime, p%t, p%dt, p%iit
    write (u) p%mem
    flush (u)
  end subroutine
END SUBROUTINE debug_dump

END MODULE gpatch_mod
