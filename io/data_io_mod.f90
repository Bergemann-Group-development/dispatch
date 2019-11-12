!===============================================================================
!> Interface from gpatch_mod to a choice of binary data I/O methods, controlled
!> by the io%method text string -- optionally (but deprecated) by io%do_legacy
!> or io%do_direct.
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
  !USE aux_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE amr_io_mod
  USE direct_io_mod
  USE legacy_io_mod
  USE parallel_io_mod
  USE buffered_io_mod
  USE pdf_io_mod
  USE counters_mod
  USE index_mod
  USE dll_mod
  USE time_slices_mod
  implicit none
  private
  type, public:: io_kind_t
    character(len=64):: filename=''
    class(patch_t), pointer:: patch=>null()
    integer:: size=0, nv=0
  end type
  type, public:: data_io_t
    logical:: initialized=.false.
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
    procedure, nopass:: close
    procedure:: open_and_close
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
  '$Id$ io/data_io_mod.f90'
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
          if      (io%method=='direct') then
            call direct_io%init (task)
          else if (io%method=='amr') then
            call amr_io%init (task)
          else if (io%method=='pan') then
            call parallel_io%init (task)
          else if (io%method=='parallel') then
            call parallel_io%init (task)
          else if (io%method=='snapshot') then
            call parallel_io%init (task)
          end if
        end select
      end if
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
  io%nwrite = io%ntask
!write(io_unit%log,'(a,f12.6,9i7)') &
!'update_counters: clk, thread, id, incr, io%nwrite, counts, iouts =', &
!wallclock(), omp%thread, patch%id, incr, io%nwrite, count1, count2, patch%iout, iout
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
  io_kind%nv = patch%nv
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
  logical:: append
  integer:: it, count, np, prv
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
    if (io%method=='amr') then
      call amr_io%check
    end if
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
    if (io%ntotal<=0) then
      call mpi%abort ('io%ntotal, which should be the total number of tasks, is zero')
    end if
    if (io%verbose > 1) &
      write (io_unit%log,*) 'io%nwrite =', io%nwrite
    !---------------------------------------------------------------------------
    if (io%method=='legacy') then
      np = 1
    else
      np = (time_slices%order+1)/2
    end if
    !$omp atomic read
    it = patch%iit(patch%nt-np)
    !$omp atomic read
    t = patch%t(it)
    if (t >= patch%out_next) then
      call self%lock%set ('output2')
      if (io%verbose > 1) &
        write (io%output,*) wallclock(),' thread',omp%thread,' waitfor output'
      !-------------------------------------------------------------------------
      ! The data io refers to task_t%iit(:) and task_t%t(:), and either all of
      ! these places need to be protected wuth omp atomic read, or (simpler) the
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
          call os%mkdir (trim(io%outputname))
        end if
        call os%mkdir_no (patch%iout)
        !-----------------------------------------------------------------------
        count = self%counters%decrement (patch%iout+1, io%nwrite)
        if (io%verbose > 0) &
          write (io_unit%log,*) 'count, nwrite =', count, io%nwrite
        append = (count < io%nwrite-1)
        !-----------------------------------------------------------------------
!write (io_unit%log,'(a,f12.6,2i6,f10.6,2i6,l4)') &
!'io_data_t%output: clk, thread, id, time, count, nwrite, append =', &
!wallclock(), omp%thread, patch%id, patch%time, count, io%nwrite, append
        if      (io%method=='direct') then
          call direct_io%output (patch)
        else if (io%method=='amr') then
          call amr_io%output (patch, count)
        else if (io%method=='legacy') then
          call legacy_io%output (patch)
          !call output_txt (patch, append)
        else if (io%method=='buffered') then
          call buffered_io%output (patch)
        else if (io%method=='pan') then
          if (io%format == 0) &
            io%format = 6
          call parallel_io%output (patch, count)
        else if (io%method=='parallel') then
          if (io%format == 0) &
            io%format = 6
          call parallel_io%output (patch, count)
        else if (io%method=='snapshot') then
          if (io%format == 0) &
            io%format = 10
          call parallel_io%output (patch, count)
        else
          call mpi%abort ('unknown io%method = '//io%method)
        end if
        !-----------------------------------------------------------------------
        if (count==0 .and. .not.io_unit%do_validate) then
          if (io%out_time > 0d0) then
            prv = modulo(patch%it-2,patch%nt) + 1
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,3g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout, &
              'size(GB):', io%gb_out, &
              'time:', patch%t(prv), patch%out_next, patch%time
          else
            write(io%output,'(1x,a,2x,a,i5.5,6x,a,f9.3,6x,a,1p,g14.5)') &
              trim(io%method)//'_io:output', &
              'snapshot '//trim(io%outputname), patch%iout, &
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
        patch%out_next = (floor(patch%time/io%out_time+1.0+1e-6))*io%out_time
        !-----------------------------------------------------------------------
        if (io%verbose > 1) & 
          write (io%output,'(a,i6,3f12.6,3i4)') &
            'io_t%output: id, iout, time, out_next, count =', &
            patch%id, patch%time, patch%time-patch%dtime, patch%out_next, &
            patch%iout, io%nwrite, count
        !-----------------------------------------------------------------------
        call os%mkdir_no (patch%iout)
      end do
      if (io%verbose > 1) &
        write (io%output,*) wallclock(),' thread',omp%thread,' unlocked output'
      call self%lock%unset ('output2')
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
!> Read snapshot from disk
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
    !$omp critical (data_io_input_cr)
    if      (io%method=='direct') then
      call direct_io%input (patch, ok)
    else if (io%method=='legacy') then
      call input_nml (patch)
      call legacy_io%input (patch, ok)
    else if (io%method=='pan') then
      call input_nml (patch)
      call parallel_io%input (patch, ok)
    else if (io%method=='parallel') then
      call input_nml (patch)
      call parallel_io%input (patch, ok)
    else if (io%method=='snapshot') then
      call input_nml (patch)
      call parallel_io%input (patch, ok)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
    if (ok) then
      count = self%counters%decrement (io%restart+1, io%nwrite)
      if (count==0 .and. io%master) then
        write(io%output,'(1x,a,2x,a,i5.5,5x,"time:",1p,g14.5)') &
          trim(io%method)//'_io:input', &
          'snapshot '//trim(io%outputname), io%restart, patch%time
      end if
      patch%t(patch%it) = patch%time
    end if
    !$omp end critical (data_io_input_cr)
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
        patch%llc_cart, patch%mesh%centre_nat, patch%quality, patch%gamma
  write (io_unit%data,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(patch%box/patch%ds+0.5), &
        u-l+1, patch%li, patch%ui, patch%n, patch%gn, patch%ng, patch%nv, patch%level
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
SUBROUTINE output_nml (patch, append)
  class(patch_t):: patch
  logical :: append
  !.............................................................................
  if (io%nml_version == 1) then
    call output_nml_v1(patch, append)
  else if (io%nml_version == 2) then
    call output_nml_v2(patch, append)
  else if (io%nml_version == 3) then
    call output_nml_v3(patch, append)
  else
    print*,io%nml_version
    call mpi%abort ('unknown nml_version!')
  end if
CONTAINS

SUBROUTINE output_nml_v1 (patch, append)
  class(patch_t):: patch
  logical :: append
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=16)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land
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
  gamma        = patch%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%nv
  nv           = merge (io%nv, patch%nv, io%nv>0 .and. io%nv<patch%nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  if (mpi%master .and. patch%iout > last_snapshot) then
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
    call file_append (io_unit%nml, io_unit%nml1)
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

SUBROUTINE output_nml_v2 (patch, append)
  class(patch_t):: patch
  logical :: append
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
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=32)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, nml_version
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, gn, &
    nv, nt, gamma, eos, opacity, periodic, guard_zones, time_derivs, no_mans_land
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, nw, &
    velocity, quality, mesh_type, kind
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
  gamma        = patch%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%nv
  nv           = merge (io%nv, patch%nv, io%nv>0 .and. io%nv<patch%nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  if (mpi%master .and. patch%iout > last_snapshot) then
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
    call file_append (io_unit%nml, io_unit%nml1)
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
END SUBROUTINE output_nml_v2

SUBROUTINE output_nml_v3 (patch, append)
  class(patch_t):: patch
  logical :: append
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time, ms
  real(8)              :: times(2), dtimes(2)
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer              :: time_derivs, now, prv
  integer(8)           :: task_offset
  logical              :: guard_zones, no_mans_land, periodic(3)
  integer, save        :: last_snapshot=-1
  character(len=16)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land
  namelist /patch_nml/ id, time, position, size, level, dtime, times, dtimes, &
    task_offset
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml')
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
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
  gamma        = patch%gamma
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%nv
  nv           = merge (io%nv, patch%nv, io%nv>0 .and. io%nv<patch%nv)
  nv           = nv*(io%time_derivs+1)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  if (mpi%master .and. patch%iout > last_snapshot) then
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
    call file_append (io_unit%nml, io_unit%nml1)
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
END SUBROUTINE output_nml_v3

END SUBROUTINE output_nml

!===============================================================================
!> Read in the basic, common meta-data from the snapshot.nml file.
!===============================================================================
SUBROUTINE input_nml (patch)
  class(patch_t):: patch
  !.............................................................................
  if (io%nml_version == 1) then
    call input_nml_v1(patch)
  else if (io%nml_version == 2) then
    call input_nml_v2(patch)
  else
    print*,io%nml_version
    call mpi%abort ('unknown nml_version!')
  end if

CONTAINS
!===============================================================================
!> Read snapshot namelists, version 1.
!===============================================================================
SUBROUTINE input_nml_v1 (patch)
  class(patch_t):: patch
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer, save        :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real   , save        :: gamma, quality
  real(8), save        :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer, save        :: time_derivs
  logical, save        :: guard_zones, no_mans_land, periodic(3)
  character(len=16)    :: kind, eos, opacity, method
  integer              :: iostat
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land
  !-----------------------------------------------------------------------------
  real(8), save        :: origin(3)
  integer, save        :: dims(3), mpi_dims(3), per_rank(3)
  logical, save        :: face_nbors, fast_nbors, omp_init
  namelist /cartesian_params/ size, dims, mpi_dims, per_rank, origin, face_nbors, &
    fast_nbors, omp_init
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml')
  !$omp critical (output_cr)
  if (first_time) then
    first_time = .false.
    patch%iout = io%restart
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%inputdir), patch%iout
    open  (unit=io_unit%nml1, file=filename, form='formatted', status='old')
    read  (io_unit%nml1, snapshot_nml, iostat=iostat)
    read  (io_unit%nml1, cartesian_params, iostat=iostat)
    close (io_unit%nml1)
    io%mpi_odims = mpi_dims
  end if
  patch%t     = time
  patch%dt    = 0d0
  patch%time  = time
  patch%istep = 0
  patch%out_next = (floor(patch%time/io%out_time+1.0+1e-6))*io%out_time
  patch%guard_zones = guard_zones
  patch%time_derivs = time_derivs
  !$omp end critical (output_cr)
  call trace%end()
END SUBROUTINE input_nml_v1

END SUBROUTINE input_nml

!===============================================================================
!> Read snapshot namelists, version 2.
!===============================================================================
SUBROUTINE input_nml_v2 (patch)
  class(patch_t):: patch
  !.............................................................................
  character(len=128)   :: filename
  character(len=128), save:: filename1='', filename2=''
  integer, save        :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real   , save        :: gamma, quality
  real(8), save        :: time, dtime, out_next, out_time, ms
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  integer, save        :: time_derivs
  logical, save        :: guard_zones, no_mans_land, periodic(3)
  character(len=32)    :: kind, eos, opacity, method
  integer              :: iostat
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, box, li, ui, ng, gn, &
    nv, nt, gamma, eos, opacity, periodic, guard_zones, time_derivs, no_mans_land
  !-----------------------------------------------------------------------------
  real(8), save        :: origin(3)
  integer, save        :: dims(3), mpi_dims(3), per_rank(3)
  logical, save        :: face_nbors, fast_nbors, omp_init
  namelist /cartesian_params/ size, dims, mpi_dims, per_rank, origin, face_nbors, &
    fast_nbors, omp_init
  namelist /patch_nml/ id, position, size, level, dtime, istep, ds, ncell, n, nw, &
    velocity, quality, mesh_type, kind
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Read info from namelist
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input_nml')
  !$omp critical (output_cr)
  patch%iout = io%restart
  if (first_time) then
    first_time = .false.
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%inputdir), patch%iout
    open  (unit=io_unit%nml1, file=filename, form='formatted', status='old')
    read  (io_unit%nml1, snapshot_nml, iostat=iostat)
    read  (io_unit%nml1, cartesian_params, iostat=iostat)
    close (io_unit%nml1)
    io%mpi_odims = mpi_dims
    io%format = ioformat
    io%ntotal = ntotal
    io%time_derivs = time_derivs
    io%guard_zones = guard_zones
  end if
  patch%time  = time
  patch%t     = time
  patch%guard_zones = guard_zones
  patch%time_derivs = time_derivs
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml" for reading
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
         trim(io%inputdir), patch%iout, mpi%rank
  if (filename /= filename2) then
    if (filename2 /= '') close (io_unit%nml2)
    open (io_unit%nml2, file=trim(filename), form='formatted', status='old')
    filename2 = filename
  end if
  rewind (io_unit%nml2)
  do
    read (io_unit%nml2, patch_nml, iostat=iostat)
    if (id == patch%id) exit
    if (iostat < 0) then
      print*,'ID = ',patch%id
      call mpi%abort('input_nml: patch_nml not found! Abort!')
    end if
  end do
  patch%dt    = dtime
  patch%dtime = dtime
  patch%istep = istep
  patch%out_next = (floor(patch%time/io%out_time+1.0+1e-6))*io%out_time
  patch%velocity = velocity
  !$omp end critical (output_cr)
  call trace%end()
END SUBROUTINE input_nml_v2

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
  else if (io%method=='direct') then
    call direct_io%close()
  else if (io%method=='buffered') then
    call buffered_io%close()
  else if (io%method=='pan') then
    call parallel_io%close()
  else if (io%method=='parallel') then
    call parallel_io%close()
  else if (io%method=='snapshot') then
    call parallel_io%close()
  end if
  call trace%end()
END SUBROUTINE close

END MODULE data_io_mod
