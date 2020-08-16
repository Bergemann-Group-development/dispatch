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
  USE initial_mod
  USE kinds_mod
  USE mesh_mod
  USE units_mod
  USE timer_mod
  USE omp_timer_mod
  USE bits_mod
  USE data_io_mod
  USE guard_zones_mod
  USE download_mod
  USE data_hub_mod
  USE validate_mod
  USE list_mod
  USE task_mod
  USE scratch_mod
  USE boundaries_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t):: initial
    type(data_hub_t):: data_hub
    class(list_t), pointer:: task_list => null()
    type(boundaries_t):: boundaries
    real(8):: d0
    character(len=8):: etype='thermal' ! or 'total', 'entropy'
  contains
    procedure, nopass:: cast2gpatch
    procedure:: init_bdries
    procedure:: init_task_list
    procedure:: pre_update
    procedure:: update
    procedure:: input
    procedure:: output
    procedure:: init
    procedure:: counter_update
    procedure:: dnload
    procedure:: courant_condition
    procedure:: to_thermal
    procedure:: clone_mem_accounting
    procedure:: rotate
    procedure:: info                    ! print info to stdout
    procedure, private :: fminval8sph
    procedure, private :: fminval4sph
    procedure, private :: fmaxval8sph
    procedure, private :: fmaxval4sph
    generic, public :: fmaxvalsph => fmaxval4sph, fmaxval8sph
    generic, public :: fminvalsph => fminval4sph, fminval8sph
  end type
  real(8), save :: restart_time = 0d0
  logical, save :: detailed_timer=.false.
  logical, save :: extrapolate_gz=.false.            ! extrapolate guard zones?
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
SUBROUTINE init_bdries (self)
  class(gpatch_t):: self
  real(8):: position(3), limit(3)
  !-----------------------------------------------------------------------------
  ! If at edge, set boundary bits
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%init_bdries')
  position = self%position - self%origin - 0.5d0 * self%box
  limit = 0.5d0 * (self%box - self%size - self%ds)
  if (.not.self%periodic(1)) then
    if (position(1)  <= -limit(1)) then
      call self%boundaries%set (bits%xl)
      self%mesh(1)%lower_boundary = .true.
    end if
    if (position(1) >= +limit(1)) then
      call self%boundaries%set (bits%xu)
      self%mesh(1)%upper_boundary = .true.
    end if
  end if
  if (.not.self%periodic(2)) then
    if (position(2)  <= -limit(2)) then
      call self%boundaries%set (bits%yl)
      self%mesh(2)%lower_boundary = .true.
    end if
    if (position(2) >= +limit(2)) then
      call self%boundaries%set (bits%yu)
      self%mesh(2)%upper_boundary = .true.
    end if
  end if
  if (.not.self%periodic(3)) then
    if (position(3)  <= -limit(3)) then
      call self%boundaries%set (bits%zl)
      self%mesh(3)%lower_boundary = .true.
    end if
    if (position(3) >= +limit(3)) then
      call self%boundaries%set (bits%zu)
      self%mesh(3)%upper_boundary = .true.
    end if
  end if
  call self%boundaries%init (self)
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
!> Post-initialize
!===============================================================================
SUBROUTINE pre_update (self)
  class(gpatch_t), target:: self
  if (associated (scratch%force_per_unit_mass)) &
    scratch%force_per_unit_mass = 0.0
  if (associated (scratch%force_per_unit_volume)) &
    scratch%force_per_unit_volume = 0.0
  if (associated (scratch%heating_per_unit_mass)) &
    scratch%heating_per_unit_mass = 0.0
  if (associated (scratch%heating_per_unit_volume)) &
    scratch%heating_per_unit_volume = 0.0
END SUBROUTINE pre_update

!===============================================================================
!> Interface allowing intercept
!===============================================================================
SUBROUTINE update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call self%patch_t%update
END SUBROUTINE update

!===============================================================================
!> Interface routine to the data_io input procedure, signalling if the patch
!> could be read, or not
!===============================================================================
SUBROUTINE input (self, ok)
  class(gpatch_t):: self
  logical:: ok
  !-----------------------------------------------------------------------------
  ok = .false.
  !call data_io%input (self, ok)
END SUBROUTINE input

!===============================================================================
!> Interface routine to the data_io output procedure
!===============================================================================
SUBROUTINE output (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call data_io%output (self)
END SUBROUTINE output

!===============================================================================
!> Add restart functionality to patch initialization.  This %init routine is
!> called for both internal and virtual patches, which both can be initialized
!> fully by the self%ininitial%condition() procedure. When instead calling the
!> self%input() procedure, we need to make a decision:  Either we must make sure
!> that the data_io%input() routine can handle off-rank reading, or we need to
!> nevertheless return ok=.true., and instead rely on a negative time to force
!> any task that relies on the patch for guard zone values to wait until data
!> has been recieved with MPI from the owner rank.
!===============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%frozen)) &
    return
  call trace%begin ('gpatch_t%init')
  call self%initial%init (self%kind, real(self%gamma))
  call data_io%init (self)
  call guard_zones%init
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
    self%istep = 0
  else 
    !---------------------------------------------------------------------------
    ! Make sure the default time for patches with IC values is zero, to ensure 
    ! boundary tasks are not considered ready for updating until their virtual
    ! task nbors have been updated, in evolved cases.
    !---------------------------------------------------------------------------
    restart_time = 0.0
    self%mem = 0.0
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff, self%idx)
    !---------------------------------------------------------------------------
    ! patch_t%setup sets time = -1, to prevent downloading until nbors have been
    ! initialized, and their time has been set to 0.0 here
    !---------------------------------------------------------------------------
    self%time = 0.0
    self%t = 0.0
  end if
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    call self%lock%set ('gpatch_t')
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = io%out_time*(1+nint(self%time/io%out_time))
    self%iout = self%restart+1
    if (io%do_output) then
      !$omp master
      write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
      if (io%iodir/=self%iout) call os%mkdir (trim(filename))
      io%iodir = self%iout
      !$omp end master
    end if
    call self%lock%unset ('gpatch_t')
  end if
  call trace%end()
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
    if (self%use_data_hub) then
      call self%data_hub%update (self%link, only=only)
    else
      call download%download_link (self%link, only=only)
    end if
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

!===============================================================================
SUBROUTINE clone_mem_accounting (self)
  class(gpatch_t):: self
  if (associated(scratch%force_per_unit_mass)) &
    call io%bits_mem (storage_size(scratch%force_per_unit_mass), &
                     product(shape(scratch%force_per_unit_mass)),'fpm')
  if (associated(scratch%force_per_unit_volume)) &
    call io%bits_mem (storage_size(scratch%force_per_unit_volume), &
                     product(shape(scratch%force_per_unit_volume)),'fpv')
  if (associated(scratch%heating_per_unit_mass)) &
    call io%bits_mem (storage_size(scratch%heating_per_unit_mass), &
                     product(shape(scratch%heating_per_unit_mass)),'qpm')
  if (associated(scratch%heating_per_unit_volume)) &
    call io%bits_mem (storage_size(scratch%heating_per_unit_volume), &
                     product(shape(scratch%heating_per_unit_volume)),'qpv')
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
  class(gpatch_t)     :: self
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
  real, save         :: d_min=0.0, d_max=huge(1.)
  real               :: fmin, fmax
  real, pointer      :: df(:,:,:)
  real(8)            :: time, print_next
  logical            :: print_it
  character(len=120) :: fmt
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%info', itimer=itimer)
  !$omp atomic write
  timer%time = self%time
  if (io_unit%do_validate) then
    it0 = merge (self%iit(self%nt-2), 1, self%nt>2)
    it1 = merge (self%iit(self%nt-1), 1, self%nt>1)
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
      self%dtime, self%u_max, fmin, fmax, &
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
  if ((io%verbose >= 0) .and. self%is_clear(bits%no_density)) then
    !call self%mem_lock(self%it)%set ('info')
    self%d_min = self%fminval(self%idx%d)
    self%d_max = self%fmaxval(self%idx%d)
    !call self%mem_lock(self%it)%unset ('info')
    !$omp atomic
    io%dmin = min(io%dmin, self%d_min)
    !$omp atomic
    io%dmax = max(io%dmax, self%d_max)
    !---------------------------------------------------------------------------
    ! After coming back to the saved ID, reset it, so we stop searching
    !---------------------------------------------------------------------------
    if (self%id == id_prv) then
      !$omp atomic write
      id_prv = 0
    end if
    !---------------------------------------------------------------------------
    ! Allow print every print_every update, if print_every > 0
    !---------------------------------------------------------------------------
    if (io%print_every > 0) then
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
    !$omp atomic update
    d_min = min(d_min, self%d_min)
    !$omp atomic update
    d_max = min(d_max, self%d_max)
    !---------------------------------------------------------------------------
    ! Write info to rank log for every time step if io%verbose > 0
    !---------------------------------------------------------------------------
    if (io%verbose > 0) then
      time = self%time
      iv = 1
      call write_info (io_unit%mpi)
    end if
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
      if (timer%n_lines==0) then
        !$omp critical (info_cr2)
        !$omp atomic write
        timer%n_lines = max_lines
        !$omp end atomic
        write(io%output,'(a)') &
        '0....+....1....+....2....+....3....+....4....+....5....+....' &
        //'6....+....7....+....8....+....9....+....0....+....1  ntsk   nq ' &
        //' chk  rdy   upd  snd  rcv  unp  snq rcq  BV'
        !$omp end critical (info_cr2)
      end if
      !-----------------------------------------------------------------------
      ! time-dtime is the task time when update started; it is guaranteed to
      ! be increasing monotomically, since the ready queue is time sorted.
      ! To allow using the column for plotting, we show this time when all
      ! times are printed.
      !-----------------------------------------------------------------------
      if (io%print_time == 0.0 .and. io%verbose == 0) then
        time = self%time - self%dtime
      else
        time = self%time
      end if
      do iv=1,nv
        call write_info (io%output)
      end do
      flush (io%output)
      !$omp atomic write
      io%dmin = huge(1.)
      !$omp atomic write
      io%dmax = tiny(1.)
      !$omp atomic write
      io%dtime = huge(1d0)
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
    real:: dmin, dmax
    !$omp atomic read
    dmin = io%dmin
    !$omp atomic read
    dmax = io%dmax
    if (self%time < 1d4 .and. self%time > 1e-4) then
      fmt = '(a,i8,2i4,i3,f12.6,1p,e12.4,0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,i6,3i5,i6,4i5,i4,l3,l1,o12)'
    else
      fmt = '(a,i8,2i4,i3,1p,2e12.5,     0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,i6,3i5,i6,4i5,i4,l3,l1,o12)'
    end if
    write(unit,fmt) &
      'upd',self%id, self%level, iv, omp_mythread, time, &
      self%dtime, log10(max(self%u_max,1e-30)), &
      dmin, dmax, &
      self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
      self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
      real(timer%n_update), ntask, nq, counts, &
      mpi_mesg%nq_send, mpi_mesg%nq_recv, &
      self%is_set(bits%boundary), self%is_set(bits%virtual)
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
  call trace%begin ('patch_t%rotate', itimer=itimer)
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

END MODULE gpatch_mod
