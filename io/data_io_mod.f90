!===============================================================================
!> Interface from gpatch_mod to a choice of binary data I/O methods, controlled
!> by the io%method text string -- optionally (but deprecated) by io%do_legacy
!> or io%do_direct.
!>
!> The output chain of routines should be called just before a solver updates,
!> at a point where both guard zone loading and boundary conditions have been
!> applied.  That ensures a consistent state, that includes the guard zones,
!> and any action that has taken place in the extras_t%pre_update() procedure.
!>
!> Calling sequence:
!>   task_list%t%update
!>     experiment_t%dnload
!>     experiment_t%bdries
!>     experiment_t%update
!>       solver_t%update
!>         solver_t%pre_update
!>           mhd_t%pre_update       ! [optional]
!>           extras_t%pre_update
!>       mhd_t%update
!>         self%output
!>           gpatch_t%output
!>             date_io_t%output
!>
!> At this point timeslot %it has been copied to timeslot%new, in preparation
!> for the update, and which timeslot data_io_t%output uses is unimportant,
!> since only the current thread is allowed to change the timeslot info.
!>
!> The correct time is held in self%time = self%t(self%it), where self%it =
!> self%iit(self%nt-1).
!>
!> To force the output to occur exactly at self%out_next, the time step in the
!> previous update, where the courant condition routine decided on a self%dtime
!> large enough to reach self%out_next, should be reduced accordingly.
!>
!> Concerning thread safety: Simultaneous calls from many threads are typical,
!> and it is the responsibility of each sub-module here (e.g. parallel_io_mod)
!> to introduce critical regions to protect agains simultaneus calls to non-
!> threadsafe procedures, or if/when accessing modified shared data.
!>
!> To prevent other threads from getting into trouble while collective MPI calls
!> are being made, the io%halt flag is set while this is happening, and since
!> the data_io%open_and_close() procedure is called in every task_list%update
!> iteration, trap can be activated there, causing all threads except the one
!> doing collective MPI I/O calls to wait.
!===============================================================================
MODULE data_io_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE timer_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE amr_io_mod
  USE legacy_io_mod
  USE background_io_mod
  USE parallel_io_mod
  USE pdf_io_mod
  USE counters_mod
  USE index_mod
  USE dll_mod
  USE time_slices_mod
  USE aux_mod
  USE eos_mod
  USE bits_mod
  implicit none
  private
  type, public:: io_kind_t
    character(len=64):: filename=''
    class(patch_t), pointer:: patch=>null()
    integer:: size=0, nv=0
  end type
  type, public:: data_io_t
    logical:: initialized=.false.
    integer:: iostat=0
    type(dll_t):: kind_list
    type(lock_t):: lock
    type(counters_t):: counters
  contains
    procedure:: init
    procedure:: update_counters
    procedure:: register_kind
    procedure:: output
    procedure, nopass:: output_txt
    procedure, nopass:: output_nml
    procedure:: input
    procedure:: input_nml
    procedure, nopass:: close
    procedure:: open_and_close
    procedure:: eof
  end type
  type(data_io_t), public:: data_io
CONTAINS

!===============================================================================
!> Initialize the I/O method, if not already done
!===============================================================================
SUBROUTINE init (self, task)
  class(data_io_t):: self
  class(patch_t):: task
  !.............................................................................
  logical, save:: first_time=.true.
  character(len=120):: id = &
  '$Id$'
  !-----------------------------------------------------------------------------
  call trace%begin('data_io_t%init')
  call trace%print_id (id)
  if (first_time) then
    !$omp critical (data_io_cr)
    if (first_time) then
      if (io%do_output) then
        call time_slices%init (task%nt)
        select type (task)
        class is (patch_t)
          select case (trim(io%method))
          case ('amr')
            call amr_io%init        (task)
          case ('parallel')
            call parallel_io%init   (task)
          case ('background')
          case ('legacy')
          case default
            call io%warning ("'amr', 'background', 'parallel' & 'legacy' supported")
            call mpi%abort ('unsupported io method = '//io%in%method)
          end select
        end select
      end if
      call aux%init
      call self%lock%init ('data_io')
      call self%kind_list%init
      self%initialized = .true.
      first_time = .false.
    end if
    !$omp end critical (data_io_cr)
  end if
  !-----------------------------------------------------------------------------
  call self%counters%init
  select type (task)
  class is (patch_t)
  call pdf_io%init (task)
  end select
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update I/O counters, when tasks are added or removed by refine_mod.
!===============================================================================
SUBROUTINE update_counters (self, patch, incr)
  class(data_io_t):: self
  class(patch_t):: patch
  integer:: incr, count1, count2, iout
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%update_counters')
  call self%lock%set ('update_counters')
  !-----------------------------------------------------------------------------
  ! Update the active counters, which may or may not already have been altered.
  ! The relevant counter index is most reliably determined by the very value
  ! that gets incremented in pdf_io_t%update(); i.e., patch%pdf_next.
  !-----------------------------------------------------------------------------
  if (pdf_io%on) then
    iout = nint(patch%pdf_next/pdf_io%out_time)
    call pdf_io%counters%update (      iout+1, io%nwrite, incr, count1)
  end if
  call   self%counters%update (patch%iout+1, io%nwrite, incr, count2)
  !-----------------------------------------------------------------------------
  ! When tasks are added or removed, the number of tasks to reset the counters
  ! to should be correspondingly updated.
  !-----------------------------------------------------------------------------
  io%ntask = io%ntask+incr
  call self%lock%unset ('update_counters')
  call trace%end ()
END SUBROUTINE update_counters

!===============================================================================
!> Register one type of task, for output into separate file
!===============================================================================
SUBROUTINE register_kind (self, patch)
  class(data_io_t):: self
  class(patch_t), pointer:: patch
  class(dll_node_t), pointer:: item
  type(io_kind_t), pointer:: io_kind
  !-----------------------------------------------------------------------------
  allocate (item, io_kind)
  item%car => io_kind
  io_kind%patch => patch
  io_kind%filename = 'snapshot_'//trim(patch%kind)//'.dat'
  io_kind%size = product(shape(patch%mem))
  io_kind%nv = patch%mv
  if (io%master) then
    print '(1x,a,4x,a,i8,4x,a,i3)', &
     'data_io_t%register_kind: filename = '//trim(io_kind%filename), &
     'array size (words) =', io_kind%size, 'nv =', io_kind%nv
  end if
  call self%kind_list%append (item)
END SUBROUTINE register_kind

!===============================================================================
!> Write results to disk
!===============================================================================
SUBROUTINE output (self, patch, experiment_name)
  class(data_io_t):: self
  class(patch_t):: patch
  !.............................................................................
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer, save:: created=-1
  integer, save:: itimer=0
  logical:: append, use_lock
  logical, save:: io_active=.false.
  integer:: it, count, np, prv, io_count
  real(8):: t
  !-----------------------------------------------------------------------------
  call pdf_io%update (patch)
  if (.not.io%do_output) return
  call trace%begin('data_io_t%output', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If io%needs_check is set, this rank needs to check if all other ranks are
  ! ready to do I/O, and if so, call the relevant I/O routine
  !-----------------------------------------------------------------------------
  call self%lock%set ('output1')
  if (io%needs_check) then
    call self%lock%unset ('output1')
  else
    call self%lock%unset ('output1')
    !---------------------------------------------------------------------------
    ! io%nwrite is the number of task writes this process is doing per snapshot
    ! it is given a default value here, which may be modified elsewhere in
    ! data_io_mod, if more than one write per task is being performed.
    !---------------------------------------------------------------------------
    if (io%ntask<=0) then
      call mpi%abort ('io%ntask, which should be the number of active tasks, is zero')
    end if
    if (io%verbose > 1) &
      write (io_unit%log,*) 'io%nwrite =', io%nwrite
    !---------------------------------------------------------------------------
    use_lock = .true.
    if (io%method=='legacy'     .or. &
        io%method=='background' .or. &
        time_slices%order <= 0) then
      np = 1
      use_lock = .false.
    else
      np = (time_slices%order+1)/2
    end if
    !$omp atomic read
    it = patch%iit(patch%nt-np)
    !$omp atomic read
    t = patch%t(it)
    if (t >= patch%out_next) then
      if (io%verbose > 0) &
        write (stdout,'(a,i6,3i4,2f12.6)') &
        'data_io_t%output: ', patch%id, np, it, patch%nt, t, patch%time
      call io%lock%set ('output2')
      if (io%verbose > 1) &
        write (io%output,*) wallclock(),' thread',omp%thread,' waitfor output'
      !-------------------------------------------------------------------------
      ! The data io refers to task_t%iit(:) and task_t%t(:), and either all of
      ! these places need to be protected with omp atomic read, or (simpler) the
      ! occasional visit here to do I/O need to block other accesses (mainly in
      ! download_mod) vi the task lock.  Certain types of output (but not for
      ! example legacy_io) may need to be protected by critical regions or a lock
      ! on data_io), but (FIXME) for now this is neglected.
      !-------------------------------------------------------------------------
      if (io%verbose > 1) &
        write (io%output,*) wallclock(),' thread',omp%thread,' locked output'
      do while (patch%t(it) >= patch%out_next)
        if (created < patch%iout) then
          created = patch%iout
          if (io%do_output) call os%mkdir (trim(io%outputname))
        end if
        if (io%do_output) call os%mkdir_next_output (patch%iout)
        !-----------------------------------------------------------------------
        count = self%counters%decrement (patch%iout+1, io%nwrite)
        if (io%verbose > 0) &
          write (io_unit%log,*) 'count, nwrite =', count, io%nwrite
        append = (count < io%nwrite-1)
        !-----------------------------------------------------------------------
        select case (io%method)
        case('amr')
          call amr_io%output        (patch, count)
        case('legacy')
          call legacy_io%output     (patch)
        case('background')
          call background_io%output (patch)
        case('parallel')
          if (io%format == 0) &
            io%format = 6
          call parallel_io%output   (patch, count)
        case('snapshot')
          if (io%format == 0) &
            io%format = 10
          call parallel_io%output   (patch, count)
        case default
          call mpi%abort ('unknown output method = '//io%method)
        end select
        !-----------------------------------------------------------------------
        if (count==0 .and. .not.io_unit%do_validate) then
        !  io_count = mod(io%counter%update(1),mpi%size)
        !  io_active = .true.
          !---------------------------------------------------------------------
          ! Write a line in the log with snapshot time info
          !---------------------------------------------------------------------
          if (io%out_time > 0d0) then
            prv = modulo(it-2,patch%nt) + 1
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,3g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout, &
              'size(GB):', io%gb_out, &
              'out_next,time:', patch%out_next, patch%t(it)
          else
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout-1, &
              'size(GB):', io%gb_out, 'time:', patch%time
          end if
        end if
        !-----------------------------------------------------------------------
        ! Output aux file, and namelist file with metadata
        !-----------------------------------------------------------------------
        call patch%aux%output (patch%iout, patch%id)
        call output_nml (patch, append)
        !-----------------------------------------------------------------------
        patch%iout = patch%iout+1
        if (io%out_time==0.0) exit
        !-----------------------------------------------------------------------
        if (io%verbose > 1) &
          write (io%output,'(a,i6,3f12.6,3i4)') &
            'io_t%output: id, iout, time, out_next, count =', &
            patch%id, patch%time, patch%time-patch%dtime, patch%out_next, &
            patch%iout, io%nwrite, count
        !-----------------------------------------------------------------------
        patch%out_next = patch%out_next + io%out_time
        if (io%do_output) call os%mkdir_next_output (patch%iout)
      end do
      if (io%verbose > 1) &
        write (io%output,*) wallclock(),' thread',omp%thread,' unlocked output'
      call io%lock%unset ('output2')
    end if
    !---------------------------------------------------------------------------
    ! Detect when an I/O operation is complete on all ranks
    !---------------------------------------------------------------------------
    if (io_active .and. use_lock) then
      call io%lock%set ('output3')
      if (io_active) then
        io_count = io%counter%update(0)
        write (io_unit%mpi,*) &
          wallclock(), 'data_io%output: io_count =', io_count
        if (mod(io_count,mpi%size) == 0) then
          !---------------------------------------------------------------------
          ! Write a line in the log with snapshot time info
          !---------------------------------------------------------------------
          if (io%out_time > 0d0) then
            prv = modulo(it-2,patch%nt) + 1
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,3g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout-1, &
              'size(GB):', io%gb_out, &
              'time:', patch%t(prv), patch%out_next-io%out_time, patch%t(it)
          else
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout-1, &
              'size(GB):', io%gb_out, 'time:', patch%time
          end if
          call io%counter%reset (mpi%size)
          io_active = .false.
        end if
      end if
      call io%lock%unset ('output3')
    end if
  end if    ! (io%needs_check)
  call trace%end (itimer)
END SUBROUTINE output

!===============================================================================
SUBROUTINE open_and_close (self)
  class(data_io_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%open_and_close')
  if (io%method == 'snapshot') then
    call parallel_io%open_and_close ()
  else if (io%method == 'amr') then
    if (io%halt) then
      call self%lock%set ('halt')
      call self%lock%unset ('halt')
    end if
  end if
  call trace%end ()
END SUBROUTINE open_and_close

!===============================================================================
!> Read snapshot from disk, using methods established by metadata
!===============================================================================
SUBROUTINE input (self, patch, ok)
  class(data_io_t):: self
  class(patch_t):: patch
  logical:: ok
  integer:: count
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input')
  ok = .false.
  if (io%restart >= 0) then
    call io%lock%set ('data_io_t%input')
    select case (trim(io%in%method))
    case ('amr')
      call amr_io%input        (patch, ok)
    case ('background')
      call background_io%input (patch, ok)
    case ('parallel')
      call parallel_io%input   (patch, ok)
    case ('legacy')
      call legacy_io%input     (patch, ok)
    case default
      call io%warning ("'amr', 'background', 'parallel' & 'legacy' supported")
      call mpi%abort ('unsupported io input method = '//io%in%method)
    end select
    if (ok) then
      patch%iout = patch%restart+1
      patch%t(patch%it) = patch%time
      call self%update_counters (patch, 0)
    end if
    call io%lock%unset ('data_io_t%input')
  end if
  call trace%end()
END SUBROUTINE input

!===============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!===============================================================================
SUBROUTINE output_txt (patch, append)
  class(patch_t):: patch
  logical:: append
  !.............................................................................
  real(8):: time
  integer:: l(3), u(3)
  character(len=128) filename
  integer:: nlines=6
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_txt')
  !$omp critical (output_cr)
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  if (io%format > 2) then
     write (filename,'(a,i5.5,"/rank_",i5.5,".txt")') trim(io%outputname), patch%iout, mpi%rank
  else
     write (filename,'(a,i5.5,"/",i5.5,".info")') trim(io%outputname), patch%iout, patch%id
  end if
  if (io%debug(2)) &
    write(io%output,'(a,3x,"time:",2f13.6)') &
      ' data_io_t%output_txt: '//trim(filename), patch%out_next, patch%time
  !-----------------------------------------------------------------------------
  ! Open the file, appending or not
  !-----------------------------------------------------------------------------
  open  (io_unit%data, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
  write (io_unit%data,'(i2.2)') io%format
  if (io%format==3 .or. io%format==4) write (io_unit%data,*) patch%id, nlines
  !-----------------------------------------------------------------------------
  ! Write the text info
  !-----------------------------------------------------------------------------
  write (io_unit%data,'(a)') trim(patch%kind)
  write (io_unit%data,'(a)') trim(patch%eos)
  write (io_unit%data,'(a)') trim(patch%opacity)
  write (io_unit%data,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') patch%out_next, &
        patch%dtime, patch%position, patch%ds, patch%velocity, patch%llc_nat, &
        patch%llc_cart, patch%mesh%centre_nat, patch%quality, eos%gamma
  write (io_unit%data,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(patch%box/patch%ds+0.5), &
        u-l+1, patch%li, patch%ui, patch%n, patch%gn, patch%ng, patch%mv, patch%level
  write (io_unit%data,'(3(1i8.1,1x),i2.1)') patch%id, patch%iout, patch%istep, patch%mesh_type
  close (io_unit%data)
  !$omp end critical (output_cr)
  call trace%end()
END SUBROUTINE output_txt

!===============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!===============================================================================
SUBROUTINE output_nml (patch, append, unit, first_time)
  class(patch_t):: patch
  logical :: append
  integer, optional:: unit
  logical, optional:: first_time
  !.............................................................................
  call io%lock%set ('output_nml')
  if (io%nml_version == 1) then
    call output_nml_v1(patch, append, unit, first_time)
  else if (io%nml_version == 2) then
    call output_nml_v2(patch, append, unit, first_time)
  else if (io%nml_version == 3) then
    call output_nml_v3(patch, append, unit, first_time)
  else if (io%nml_version == 4) then
    call output_nml_v4(patch, append, unit, first_time)
  else
    print*,io%nml_version
    call mpi%abort ('unknown nml_version!')
  end if
  call io%lock%unset ('output_nml')
CONTAINS

!===============================================================================
SUBROUTINE get_ids (nbor_ids, parent_id)
  integer, allocatable:: nbor_ids(:)
  integer:: parent_id, n
  class(nbor_t), pointer:: nbor
  !.............................................................................
  n = 0
  nbor => patch%link%nbor
  do while (associated(nbor))
    n = n+1
    nbor => nbor%next
  end do
  allocate (nbor_ids(n))
  n = 0
  nbor => patch%link%nbor
  do while (associated(nbor))
    n = n+1
    nbor_ids(n) = nbor%task%id
    nbor => nbor%next
  end do
  if (associated(patch%parent)) then
    parent_id  = patch%parent%id
  else
    parent_id  = 0
  end if
END SUBROUTINE get_ids

!===============================================================================
SUBROUTINE output_nml_v1 (patch, append, unit, snapshot)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  logical :: append
  integer, optional:: unit
  logical, optional:: snapshot
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format, nml_version
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs
  integer              :: omp_nthreads, mpi_size, mpi_dims(3)
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=16)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, mpi_dims
  namelist /patch_nml/ id, time, position, size, level, dtime
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml')
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
  nml_version = io%nml_version
  ntotal      = io%ntotal
  out_time    = io%out_time
  guard_zones = io%guard_zones
  time_derivs = io%time_derivs
  method      = io%method
  iout        = patch%iout
  time        = merge (patch%time, patch%out_next, io%out_time==0)
  dtime       = patch%dtime
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  ncell        = u-l+1
  ! Special consideration for ZEUS solvers because the include one extra cell
  ! for staggered quantities (all other variables should be zero in the extra
  ! cell).
  if (patch%kind(1:4) == 'zeus') then
    where (patch%n > 1)
      ncell = ncell + 1
    end where
  endif
  id           = patch%id
  istep        = patch%istep
  mesh_type    = patch%mesh_type
  position     = patch%position
  size         = patch%size
  ds           = patch%ds
  box          = patch%box
  velocity     = patch%velocity
  level        = patch%level
  quality      = patch%quality
  gamma        = eos_t%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%mv
  nv           = merge (io%nv, nv, io%nv>0 .and. io%nv<nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  mpi_dims     = mpi_coords%dims
  if (mpi%master .and. patch%iout > last_snapshot) then
    last_snapshot = patch%iout
    mpi_size = mpi%size
    omp_nthreads = omp%nthreads
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), patch%iout
    if (filename /= filename1) then
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opening ', trim(filename)
      flush (io_unit%log)
      open  (io_unit%nml1, file=trim(filename), form='formatted', status='unknown')
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opened  ', trim(filename)
      flush (io_unit%log)
      filename1 = filename
    end if
    write (io_unit%nml1, io_nml)
    write (io_unit%nml1, snapshot_nml)
    call patch%idx%output (io_unit%nml1)
    !call file_append (io_unit%nml, io_unit%nml1)
    flush (io_unit%nml1)
    if (io%debug(2)) write(io%output,'(a,3x,"time:",2f13.6)') &
      ' data_io_t%output_nml: '//trim(filename), patch%out_next, patch%time
  end if
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml", appending or not
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
         trim(io%outputname), patch%iout, mpi%rank
  if (filename /= filename2) then
    if (filename2 /= '') close (io_unit%nml2)
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opening ', trim(filename)
    flush (io_unit%log)
    open  (io_unit%nml2, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opened  ', trim(filename)
    flush (io_unit%log)
    filename2 = filename
  end if
  !-----------------------------------------------------------------------------
  ! Write the namelist info, including a leading idx_nml
  !-----------------------------------------------------------------------------
  if (mpi%master .and. .not.append) then
    call patch%idx%output (io_unit%nml2)
  end if
  write (io_unit%nml2, patch_nml)
  flush (io_unit%nml2)
  call trace%end()
END SUBROUTINE output_nml_v1

SUBROUTINE output_nml_v2 (patch, append, unit, snapshot)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  logical :: append
  integer, optional:: unit
  logical, optional:: snapshot
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format, nml_version
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs, refine_ratio
  integer              :: omp_nthreads, mpi_size, mpi_dims(3)
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=32)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, gn, n, &
    nv, nt, gamma, eos, opacity, periodic, guard_zones, time_derivs, no_mans_land, &
    omp_nthreads, mpi_size, mesh_type, mpi_dims
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, nw, &
    velocity, quality, mesh_type, kind
  integer, allocatable :: nbor_ids(:)
  integer              :: parent_id
  namelist /nbor_nml/ id, parent_id, nbor_ids
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml')
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
  nml_version = io%nml_version
  ntotal      = io%ntotal
  out_time    = io%out_time
  guard_zones = io%guard_zones
  time_derivs = io%time_derivs
  method      = io%method
  iout        = patch%iout
  time        = merge (patch%time, patch%out_next, io%out_time==0)
  dtime       = patch%dtime
  istep       = patch%istep
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  ncell        = u-l+1
  ! Special consideration for ZEUS solvers because the include one extra cell
  ! for staggered quantities (all other variables should be zero in the extra
  ! cell).
  if (patch%kind(1:4) == 'zeus') then
    where (patch%n > 1)
      ncell = ncell + 1
    end where
  endif
  id           = patch%id
  istep        = patch%istep
  mesh_type    = patch%mesh_type
  position     = patch%position
  size         = patch%size
  ds           = patch%ds
  box          = patch%box
  velocity     = patch%velocity
  level        = patch%level
  quality      = patch%quality
  gamma        = eos_t%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%mv
  nv           = merge (io%nv, nv, io%nv>0 .and. io%nv<nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  mpi_dims     = mpi_coords%dims
  if (mpi%master .and. patch%iout > last_snapshot) then
    mpi_size = mpi%size
    omp_nthreads = omp%nthreads
    last_snapshot = patch%iout
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), patch%iout
    if (filename /= filename1) then
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opening ', trim(filename)
      flush (io_unit%log)
      open  (io_unit%nml1, file=trim(filename), form='formatted', status='unknown')
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opened  ', trim(filename)
      flush (io_unit%log)
      filename1 = filename
    end if
    write (io_unit%nml1, io_nml)
    write (io_unit%nml1, snapshot_nml)
    call patch%idx%output (io_unit%nml1)
    !call file_append (io_unit%nml, io_unit%nml1)
    flush (io_unit%nml1)
    if (io%debug(2)) write(io%output,'(a,3x,"time:",2f13.6)') &
      ' data_io_t%output_nml: '//trim(filename), patch%out_next, patch%time
  end if
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml", appending or not
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
         trim(io%outputname), patch%iout, mpi%rank
  if (filename /= filename2) then
    if (filename2 /= '') close (io_unit%nml2)
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opening ', trim(filename)
    flush (io_unit%log)
    open  (io_unit%nml2, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opened  ', trim(filename)
    flush (io_unit%log)
    filename2 = filename
  end if
  !-----------------------------------------------------------------------------
  ! Write the namelist info, including a leading idx_nml
  !-----------------------------------------------------------------------------
  if (mpi%master .and. .not.append) then
    call patch%idx%output (io_unit%nml2)
  end if
  write (io_unit%nml2, patch_nml)
  call get_ids (nbor_ids, parent_id)
  write (io_unit%nml2, nbor_nml)
  flush (io_unit%nml2)
  deallocate (nbor_ids)
  call trace%end()
END SUBROUTINE output_nml_v2

SUBROUTINE output_nml_v3 (patch, append, unit, snapshot)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  logical :: append
  integer, optional:: unit
  logical, optional:: snapshot
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format, nml_version
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time, ms
  real(8)              :: times(2), dtimes(2)
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs, now, prv
  integer              :: omp_nthreads, mpi_size, mpi_dims(3)
  integer(8)           :: task_offset
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=16)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, mpi_dims
  namelist /patch_nml/ id, time, position, size, level, dtime, times, dtimes, &
    task_offset
  integer, allocatable :: nbor_ids(:)
  integer              :: parent_id
  namelist /nbor_nml/ id, parent_id, nbor_ids
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml_v3')
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
  nml_version = io%nml_version
  ntotal      = io%ntotal
  out_time    = io%out_time
  guard_zones = io%guard_zones
  time_derivs = io%time_derivs
  method      = io%method
  iout        = patch%iout
  time        = merge (patch%time, patch%out_next, io%out_time==0)
  dtime       = patch%dtime
  now         = patch%it
  task_offset = patch%amr_offset
  if (patch%istep > 1) then
    prv       = modulo(now-2,patch%nt) + 1
  else
    prv       = now
  end if
  times(1)    = patch%t (prv)
  times(2)    = patch%t (now)
  dtimes(1)   = patch%dt(prv)
  dtimes(2)   = patch%dt(now)
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  ncell        = u-l+1
  ! Special consideration for ZEUS solvers because the include one extra cell
  ! for staggered quantities (all other variables should be zero in the extra
  ! cell).
  if (patch%kind(1:4) == 'zeus') then
    where (patch%n > 1)
      ncell = ncell + 1
    end where
  endif
  id           = patch%id
  istep        = patch%istep
  mesh_type    = patch%mesh_type
  position     = patch%position
  size         = patch%size
  ds           = patch%ds
  box          = patch%box
  velocity     = patch%velocity
  level        = patch%level
  quality      = patch%quality
  gamma        = eos_t%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%mv
  nv           = merge (io%nv, nv, io%nv>0 .and. io%nv<nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  mpi_dims     = mpi_coords%dims
  if (mpi%master .and. patch%iout > last_snapshot) then
    mpi_size = mpi%size
    omp_nthreads = omp%nthreads
    last_snapshot = patch%iout
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), patch%iout
    if (filename /= filename1) then
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opening ', trim(filename)
      flush (io_unit%log)
      open  (io_unit%nml1, file=trim(filename), form='formatted', status='unknown')
      write(io_unit%log,'(a,i5,2x,2a)') &
       ' thread', omp%thread, 'opened  ', trim(filename)
      flush (io_unit%log)
      filename1 = filename
    end if
    write (io_unit%nml1, io_nml)
    write (io_unit%nml1, snapshot_nml)
    call patch%idx%output (io_unit%nml1)
    !call file_append (io_unit%nml, io_unit%nml1)
    flush (io_unit%nml1)
    if (io%debug(2)) write(io%output,'(a,3x,"time:",2f13.6)') &
      ' data_io_t%output_nml: '//trim(filename), patch%out_next, patch%time
  end if
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml", appending or not
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
         trim(io%outputname), patch%iout, mpi%rank
  if (filename /= filename2) then
    if (filename2 /= '') close (io_unit%nml2)
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opening ', trim(filename)
    flush (io_unit%log)
    open  (io_unit%nml2, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
    write(io_unit%log,'(a,i5,2x,2a)') &
     ' thread', omp%thread, 'opened  ', trim(filename)
    flush (io_unit%log)
    filename2 = filename
  end if
  !-----------------------------------------------------------------------------
  ! Write the namelist info, including a leading idx_nml
  !-----------------------------------------------------------------------------
  if (mpi%master .and. .not.append) then
    call patch%idx%output (io_unit%nml2)
  end if
  write (io_unit%nml2, patch_nml)
  call get_ids (nbor_ids, parent_id)
  write (io_unit%nml2, nbor_nml)
  flush (io_unit%nml2)
  deallocate (nbor_ids)
  call trace%end()
END SUBROUTINE output_nml_v3

SUBROUTINE output_nml_v4 (patch, append, unit, first_time)
  class(patch_t):: patch
  logical :: append
  integer, optional:: unit
  logical, optional:: first_time
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format, nml_version
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: quality, cost
  real(8)              :: gamma, time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs, refine_ratio, record, rank
  integer              :: omp_nthreads, mpi_size, mpi_dims(3)
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=16)    :: kind, eos_name, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, &
    gn, n, nv, nt, gamma, eos_name, opacity, periodic, guard_zones, &
    time_derivs, no_mans_land, omp_nthreads, mpi_size, mesh_type, mpi_dims, &
    refine_ratio
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, &
    nw, velocity, quality, mesh_type, kind, record, rank, cost
  integer, allocatable :: nbor_ids(:)
  integer              :: parent_id
  namelist /nbor_nml/ id, parent_id, nbor_ids
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml_v4')
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
  nml_version = io%nml_version
  ntotal      = io%ntotal
  out_time    = io%out_time
  guard_zones = io%guard_zones
  time_derivs = io%time_derivs
  method      = io%method
  iout        = patch%iout
  time        = merge (patch%time, patch%out_next, io%out_time==0)
  dtime       = patch%dtime
  istep       = patch%istep
  record      = patch%record
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  ncell        = u-l+1
  ! Special consideration for ZEUS solvers because the include one extra cell
  ! for staggered quantities (all other variables should be zero in the extra
  ! cell).
  if (patch%kind(1:4) == 'zeus') then
    where (patch%n > 1)
      ncell = ncell + 1
    end where
  endif
  id           = patch%id
  istep        = patch%istep
  mesh_type    = patch%mesh_type
  position     = patch%position
  size         = patch%size
  ds           = patch%ds
  box          = patch%box
  velocity     = patch%velocity
  level        = patch%level
  quality      = patch%quality
  gamma        = eos%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%mv
  nv           = merge (io%nv, nv, io%nv>0 .and. io%nv<nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos_name     = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  refine_ratio = patch%refine_ratio
  rank         = patch%rank
  cost         = patch%cost
  mpi_dims     = mpi_coords%dims
  if (present(unit)) then
    if (present(first_time)) then
      if (first_time) then
        write (unit, snapshot_nml)
      end if
    end if
    write (unit, patch_nml)
  else
    if (mpi%master .and. patch%iout > last_snapshot) then
      mpi_size = mpi%size
      omp_nthreads = omp%nthreads
      last_snapshot = patch%iout
      write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), patch%iout
      if (filename /= filename1) then
        write(io_unit%log,'(a,i5,2x,2a)') &
        ' thread', omp%thread, 'opening ', trim(filename)
        flush (io_unit%log)
        open  (io_unit%nml1, file=trim(filename), form='formatted', status='unknown')
        write(io_unit%log,'(a,i5,2x,2a)') &
        ' thread', omp%thread, 'opened  ', trim(filename)
        flush (io_unit%log)
        filename1 = filename
      end if
      write (io_unit%nml1, io_nml)
      write (io_unit%nml1, snapshot_nml)
      call patch%idx%output (io_unit%nml1)
      !call file_append (io_unit%nml, io_unit%nml1)
      flush (io_unit%nml1)
      if (io%debug(2)) write(io%output,'(a,3x,"time:",2f13.6)') &
        ' data_io_t%output_nml: '//trim(filename), patch%out_next, patch%time
    end if
    !---------------------------------------------------------------------------
    ! Open the file "data/run/rank_rrrrr_patches.nml", appending or not
    !---------------------------------------------------------------------------
    write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
          trim(io%outputname), patch%iout, mpi%rank
    if (filename /= filename2) then
      if (filename2 /= '') close (io_unit%nml2)
      write(io_unit%log,'(a,i5,2x,2a)') &
      ' thread', omp%thread, 'opening ', trim(filename)
      flush (io_unit%log)
      open  (io_unit%nml2, file=trim(filename), form='formatted', status='unknown', &
          access=trim(merge('append    ','sequential',append)))
      write (io_unit%nml2, snapshot_nml)
      write(io_unit%log,'(a,i5,2x,2a)') &
      ' thread', omp%thread, 'opened  ', trim(filename)
      flush (io_unit%log)
      filename2 = filename
    end if
    !---------------------------------------------------------------------------
    ! Write the namelist info, including a leading idx_nml
    !---------------------------------------------------------------------------
    if (mpi%master .and. .not.append) then
      call patch%idx%output (io_unit%nml2)
    end if
    write (io_unit%nml2, patch_nml)
    call get_ids (nbor_ids, parent_id)
    write (io_unit%nml2, nbor_nml)
    flush (io_unit%nml2)
    deallocate (nbor_ids)
  end if
  call trace%end()
END SUBROUTINE output_nml_v4

END SUBROUTINE output_nml

!===============================================================================
!> Read in the basic, common meta-data from the snapshot.nml file.  If the
!> patch%time is already set, assume that the metadata has been read already.
!===============================================================================
SUBROUTINE input_nml (self, patch, unit)
  class(data_io_t):: self
  class(patch_t):: patch
  integer:: unit
  !.............................................................................
  call trace%begin ('data_io_t%input_nml')
  if (patch%time <= 0d0) then
    !---------------------------------------------------------------------------
    ! Establish the nml_version
    !---------------------------------------------------------------------------
    call input_io (unit)
    select case (io%in%nml_version)
    case (1)
      call input_nml_v1(patch, unit)
    case (2)
      call input_nml_v2(patch, unit)
    case (3)
      call input_nml_v3(patch, unit)
    case (4)
      call input_nml_v4(self, patch, unit)
    case default
      print*,io%in%nml_version
      call mpi%abort ('unknown nml_version!')
    end select
    !---------------------------------------------------------------------------
    ! Make sure these tasks have the output bit and the out_next time set
    !---------------------------------------------------------------------------
    call patch%set (bits%has_output)
    patch%out_next = nint(patch%time/io%out_time + 1d0)*io%out_time
  end if
  call trace%end()
CONTAINS

!===============================================================================
!> Common part of metadata namelist reading
!===============================================================================
SUBROUTINE input_io (unit)
  implicit none
  integer:: unit
  !.............................................................................
  integer            :: format, ntotal, time_derivs, nml_version
  logical            :: guard_zones
  real(8)            :: out_time
  character(len=16)  :: method
  character(len=128) :: filename
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  !-----------------------------------------------------------------------------
  ! Use open unit, or open one
  !-----------------------------------------------------------------------------
  if (unit == 0) then
    unit = io_unit%nml1
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%inputdir), io%restart
print*,trim(filename)
    open  (unit=unit, file=filename, form='formatted', status='old')
  end if
  rewind (unit)
  read (unit, io_nml)
  io%in%format      = format
  io%in%ntotal      = ntotal
  io%in%out_time    = out_time
  io%in%guard_zones = guard_zones
  io%in%time_derivs = time_derivs
  io%in%method      = method
  io%in%nml_version = nml_version
END SUBROUTINE input_io

!===============================================================================
!> Read snapshot namelists, version 1.
!===============================================================================
SUBROUTINE input_nml_v1 (patch, unit)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  integer:: unit
  !.............................................................................
  integer              :: iostat, nml_version
  integer, save        :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format, mpi_size
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real   , save        :: gamma, quality
  real(8), save        :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer, save        :: time_derivs
  logical, save        :: guard_zones, no_mans_land, periodic(3)
  character(len=16)    :: kind, eos, opacity, method
  integer              :: omp_nthreads, mpize, mpi_dims(3)
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, mpi_dims
  namelist /patch_nml/ id, time, position, size, level, dtime
  !-----------------------------------------------------------------------------
  real(8), save        :: origin(3)
  integer, save        :: dims(3)
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml_v1')
  call io%lock%set ('input_nml_v1')
  patch%iout = io%restart
  !-----------------------------------------------------------------------------
  ! First time, rewind and read general snapshot metadata, then specific
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    read (unit, snapshot_nml)
  ! namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
  !   position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
  !   ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
  !   omp_nthreads, mpi_size, no_mans_land, mpi_dims
    io%in%format       = ioformat
    io%in%ntotal       = ntotal
    io%in%guard_zones  = guard_zones
    io%in%time_derivs  = time_derivs
    io%in%no_mans_land = no_mans_land
    io%in%mpi_dims     = mpi_dims
  end if
  patch%iout         = iout
  patch%time         = time
  patch%istep        = istep
  patch%mesh_type    = mesh_type
  patch%position     = position
  patch%size         = size
  patch%ds           = ds
  patch%box          = box
  patch%velocity     = velocity
  patch%level        = level
  patch%quality      = quality
  eos_t%gamma        = gamma
  patch%ncell        = ncell
  patch%li           = li
  patch%ui           = ui
  patch%n            = n
  patch%ng           = ng
  patch%gn           = gn
  patch%nv           = nv
  patch%nt           = nt
  patch%nw           = nw
  patch%kind         = kind
  patch%eos          = eos
  patch%opacity      = opacity
  patch%periodic     = periodic
  !-----------------------------------------------------------------------------
  ! Read patch-specific values
  !-----------------------------------------------------------------------------
  read (unit, patch_nml, iostat=iostat)
  ! namelist /patch_nml/ id, time, position, size, level, dtime
  if (iostat /= 0) then
    self%iostat = unit
  end if
  patch%id           = id
  patch%time         = time
  patch%position     = position
  patch%size         = size
  patch%level        = level
  patch%dtime        = dtime
  !-----------------------------------------------------------------------------
  call io%lock%unset ('input_nml_v1')
  call trace%end()
END SUBROUTINE input_nml_v1

!===============================================================================
!> Read snapshot namelists, version 2.  Relative to version 1, this has
!> individual values for patch attributes that differ between AMR levels
!===============================================================================
SUBROUTINE input_nml_v2 (patch, unit)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  integer:: unit
  !.............................................................................
  integer              :: iostat, unit_l
  integer,            save    :: ioformat, id, iout, istep, mesh_type, &
                                 level, nv, nt, nw, ntotal, format
  integer, dimension(3), save :: ncell, li, ui, n, ng, gn, l, u
  real   ,               save :: gamma, quality
  real(8),               save :: time, dtime, out_next, out_time, ms
  real(8), dimension(3), save :: size, position, ds, box, velocity, &
                                 llc_nat, llc_cart, centre_nat
  integer,               save :: time_derivs, refine_ratio
  logical,               save :: guard_zones, no_mans_land, periodic(3)
  character(len=32),     save :: kind, eos, opacity, method
  integer,               save :: dims(3), mpi_dims(3), mpi_size
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, &
    gn, n, nv, nt, gamma, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land, mpi_dims
  !-----------------------------------------------------------------------------
  real(8), save        :: origin(3)
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, nw, &
    velocity, quality, mesh_type, kind
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml_v2')
  call io%lock%set ('input_nml_v2')
  patch%iout = io%restart
  !-----------------------------------------------------------------------------
  ! First time, rewind and read general snapshot metadata, then specific
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    read (unit, snapshot_nml)
    io%in%format       = ioformat
    io%in%ntotal       = ntotal
    io%in%guard_zones  = guard_zones
    io%in%time_derivs  = time_derivs
    io%in%no_mans_land = no_mans_land
    io%in%mpi_dims     = mpi_dims
  end if
  ! namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, &
  !   gn, n, nv, nt, gamma, eos, opacity, periodic, guard_zones, time_derivs, &
  !   no_mans_land, mpi_dims
  patch%iout         = iout
  patch%time         = time
  patch%box          = box
  patch%li           = li
  patch%ui           = ui
  patch%ng           = ng
  patch%gn           = gn
  patch%nv           = nv
  patch%nt           = nt
  eos_t%gamma        = gamma
  patch%eos          = eos
  patch%opacity      = opacity
  patch%periodic     = periodic
  !-----------------------------------------------------------------------------
  ! Read patch-specific values
  !-----------------------------------------------------------------------------
  read (unit, patch_nml, iostat=iostat)
  if (iostat /= 0) then
    self%iostat = unit
  end if
  ! namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, nw, &
  !   velocity, quality, mesh_type, kind
  patch%id           = id
  patch%position     = position
  patch%size         = size
  patch%level        = level
  patch%dtime        = dtime
  patch%istep        = istep
  patch%ds           = ds
  patch%ncell        = ncell
  patch%n            = n
  patch%nw           = nw
  patch%velocity     = velocity
  patch%quality      = quality
  patch%mesh_type    = mesh_type
  patch%kind         = kind
  !-----------------------------------------------------------------------------
  call io%lock%unset ('input_nml_v2')
  call trace%end()
END SUBROUTINE input_nml_v2

!===============================================================================
!> Read snapshot namelists, version 3.  Relative to version 2, this has
!> less information in the patch-specific part (to reduce Python reading time)
!===============================================================================
SUBROUTINE input_nml_v3 (patch, unit)
  USE eos_mod, only: eos_t => eos
  class(patch_t):: patch
  integer:: unit
  !.............................................................................
  integer              :: iostat, omp_nthreads
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt
  integer              :: nw, ntotal, format, time_derivs
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u, mpi_size
  integer(8)           :: task_offset
  real                 :: quality
  real(8)              :: gamma, time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity
  real(8), dimension(3):: llc_nat, llc_cart, centre_nat
  logical              :: guard_zones, no_mans_land, periodic(3)
  character(len=32)    :: kind, eos, opacity, method
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, mpi_dims
  save                 :: ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, mpi_dims

  !-----------------------------------------------------------------------------
  real(8)              :: times(2), dtimes(2)
  real(8)              :: origin(3)
  integer              :: dims(3), mpi_dims(3)
  namelist /patch_nml/ id, time, position, size, level, dtime, times, dtimes, &
    task_offset
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml_v2')
  call io%lock%set ('input_nml_v2')
  patch%iout = io%restart
  !-----------------------------------------------------------------------------
  ! First time, rewind and read general snapshot metadata, then specific
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    read (unit, snapshot_nml)
    io%in%format       = ioformat
    io%in%ntotal       = ntotal
    io%in%guard_zones  = guard_zones
    io%in%time_derivs  = time_derivs
    io%in%no_mans_land = no_mans_land
    io%in%mpi_dims     = mpi_dims
  end if
  ! namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
  !   position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
  !   ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
  !   omp_nthreads, mpi_size, no_mans_land, mpi_dims
  patch%iout         = iout
  patch%box          = box
  patch%velocity     = velocity            ! FIXME
  patch%istep        = istep               ! FIXME
  patch%mesh_type    = mesh_type
  patch%ncell        = ncell
  patch%li           = li
  patch%ui           = ui
  patch%n            = n
  patch%ng           = ng
  patch%gn           = gn
  patch%nv           = nv
  patch%nt           = nt
  patch%nw           = nw
  eos_t%gamma        = gamma
  patch%eos          = eos
  patch%opacity      = opacity
  patch%periodic     = periodic
  !-----------------------------------------------------------------------------
  ! Read patch-specific values
  !-----------------------------------------------------------------------------
  read (unit, patch_nml, iostat=iostat)
  if (iostat /= 0) then
    self%iostat = unit
  end if
  ! namelist /patch_nml/ id, time, position, size, level, dtime, times, dtimes, &
  !   task_offset
  patch%id           = id
  patch%time         = time
  patch%position     = position
  patch%size         = size
  patch%level        = level
  patch%dtime        = dtime
  patch%quality      = quality
  patch%kind         = kind
  !-----------------------------------------------------------------------------
  ! Derived values
  !-----------------------------------------------------------------------------
  patch%ds           = patch%size/patch%n
  !-----------------------------------------------------------------------------
  call io%lock%unset ('input_nml_v3')
  call trace%end()
END SUBROUTINE input_nml_v3

!===============================================================================
!> Read snapshot namelists, version 4.
!===============================================================================
SUBROUTINE input_nml_v4 (self, patch, unit)
  class(data_io_t):: self
  class(patch_t):: patch
  integer:: unit
  !.............................................................................
  integer, save               :: unit_l=0, rank_in=0
  integer                     :: iostat
  integer                     :: ioformat, id, iout, istep, mesh_type, record
  integer                     :: level, nv, nt, nw, ntotal, format, rank
  integer                     :: time_derivs, refine_ratio, omp_nthreads
  integer, dimension(3)       :: ncell, li, ui, n, ng, gn, l, u
  real                        :: gamma, quality, cost
  real(8)                     :: time, dtime, out_next, out_time, ms
  real(8), dimension(3)       :: size, position, ds, box, velocity, &
                                 llc_nat, llc_cart, centre_nat
  logical                     :: guard_zones, no_mans_land, periodic(3)
  character(len=16)           :: kind, eos_name, opacity, method
  character(len=128)          :: filename
  integer                     :: dims(3), mpi_dims(3), mpi_size
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, &
    gn, n, nv, nt, gamma, eos_name, opacity, periodic, guard_zones, &
    time_derivs, no_mans_land, omp_nthreads, mpi_size, mesh_type, mpi_dims, &
    refine_ratio
  save                        :: ioformat, iout, time, ntotal, box, li, ui, ng, &
    gn, n, nv, nt, gamma, eos_name, opacity, periodic, guard_zones, &
    time_derivs, no_mans_land, omp_nthreads, mpi_size, mesh_type, mpi_dims, &
    refine_ratio
  !-----------------------------------------------------------------------------
  real(8), save        :: origin(3)
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, &
    nw, velocity, quality, mesh_type, kind, record, rank, cost
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml_v4')
  call io%lock%set ('input_nml_v4')
  patch%iout = io%restart
  !-----------------------------------------------------------------------------
  ! First time, rewind and read general snapshot metadata, then specific
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    read (unit, snapshot_nml)
  ! namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, &
  !   gn, n, nv, nt, gamma, eos_name, opacity, periodic, guard_zones, &
  !   time_derivs, no_mans_land, omp_nthreads, mpi_size, mesh_type, mpi_dims, &
  !   refine_ratio
    io%in%format       = ioformat
    io%in%ntotal       = ntotal
    io%in%guard_zones  = guard_zones
    io%in%time_derivs  = time_derivs
    io%in%no_mans_land = no_mans_land
    io%in%mpi_size     = mpi_size
    io%in%mpi_dims     = mpi_dims
  end if
  patch%time         = time
  patch%box          = box
  patch%li           = li
  patch%ui           = ui
  patch%ng           = ng
  patch%gn           = gn
  patch%n            = n
  io%nv              = nv
  patch%nt           = nt
  eos%gamma          = gamma
  patch%eos          = eos_name
  patch%opacity      = opacity
  patch%periodic     = periodic
  patch%guard_zones  = guard_zones
  patch%time_derivs  = time_derivs
  patch%no_mans_land = no_mans_land
  patch%mesh_type    = mesh_type
  patch%refine_ratio = refine_ratio
  !-----------------------------------------------------------------------------
  ! Read patch-specific values
  !-----------------------------------------------------------------------------
  if (unit_l == 0) then
    write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
      trim(io%inputdir), io%restart, rank_in
    open (file=filename, newunit=unit_l, form='formatted', status='unknown')
    rewind (unit_l)
  end if
  read (unit_l, patch_nml, iostat=iostat)
  if (iostat /= 0 .and. rank_in < mpi%size) then
    self%iostat = unit_l
    unit_l = 0
    rank_in = rank_in + 1
  end if
  if (rank_in == mpi%size) then
    self%iostat = unit
  else
! namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, &
!   nw, velocity, quality, mesh_type, kind, record, rank, cost
    patch%id        = id
    patch%position  = position
    patch%size      = size
    patch%level     = level
    patch%dtime     = dtime
    patch%istep     = istep
    patch%ds        = ds
    patch%ncell     = ncell
    patch%n         = n
    patch%nw        = nw
    patch%velocity  = velocity
    patch%quality   = quality
    patch%mesh_type = mesh_type
    patch%kind      = kind
    patch%record    = record
    patch%rank      = rank
    patch%cost      = cost
  end if
  !-----------------------------------------------------------------------------
  call io%lock%unset ('input_nml_v4')
  call trace%end()
END SUBROUTINE input_nml_v4

END SUBROUTINE input_nml

!===============================================================================
!> Copy one text file to another
!===============================================================================
SUBROUTINE file_append (unit1, unit2)
  integer:: unit1, unit2
  character(len=128):: line
  integer:: iostat
  !-----------------------------------------------------------------------------
  rewind (unit1)
  do while (.true.)
    read (unit1,'(a)',iostat=iostat) line
    if (iostat /= 0) exit
    write (unit2,'(a)') trim(line)
  end do
  backspace (unit1)
END SUBROUTINE file_append

!===============================================================================
!> Close the data I/O
!===============================================================================
SUBROUTINE close()
  !.............................................................................
  if (.not.io%do_output) return
  call trace%begin('patch_t%close')
  if      (io%method=='legacy') then
    !call legacy_io%close()
  else if (io%method=='pan') then
    call parallel_io%close()
  else if (io%method=='parallel') then
    call parallel_io%close()
  else if (io%method=='snapshot') then
    call parallel_io%close()
  end if
  call trace%end()
END SUBROUTINE close

!===============================================================================
!> Return true when a unit has recently reached end-of-file
!===============================================================================
LOGICAL FUNCTION eof (self, unit)
  class(data_io_t):: self
  integer:: unit
  !-----------------------------------------------------------------------------
  if (self%iostat == unit) then
    self%iostat = 0
    eof = .true.
  else
    eof = .false.
  end if
END FUNCTION eof

END MODULE data_io_mod
