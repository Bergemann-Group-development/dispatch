!> =============================================================================
!> Interface from gpatch_mod to a choice of binary data I/O methods, controlled
!> by the io%method text string -- optionally (but deprecated) by io%do_legacy
!> or io%do_direct.
!>
!> Concerning thread safety: Simultaneous calls from many threads are typical,
!> and it is the responsibility of each sub-module here (e.g. parallel_io_mod)
!> to introduce critical regions to protect agains simultaneus calls to non-
!> threadsafe procedures, or if/when accessing modified shared data.
!> =============================================================================
MODULE data_io_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE direct_io_mod
  USE legacy_io_mod
  USE parallel_io_mod
  USE buffered_io_mod
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
  contains
    procedure:: init
    procedure:: register_kind
    procedure:: output
    procedure, nopass:: output_txt
    procedure, nopass:: output_nml
    procedure:: input
    procedure, nopass:: close
    procedure, nopass:: open_and_close
  end type
  type(data_io_t), public:: data_io
CONTAINS

!> =============================================================================
!> Initialize the I/O method, if not already done
!> =============================================================================
SUBROUTINE init (self, task)
  class(data_io_t):: self
  class(patch_t):: task
  !.............................................................................
  if (self%initialized) return
  !$omp critical (data_io_cr)
  if (.not.self%initialized) then
    call trace%begin('data_io_t%init')
    if (io%do_output) then
      select type (task)
      class is (patch_t)
        if      (io%method=='direct') then
          call direct_io%init (task)
        else if (io%method=='pan') then
          call parallel_io%init (task)
        else if (io%method=='parallel') then
          call parallel_io%init (task)
        else if (io%method=='snapshot') then
          call parallel_io%init (task)
        end if
      end select
    end if
    call self%kind_list%init
    self%initialized = .true.
    if (io%master) print *,'data_io_t%init: io%nwrite =', io%nwrite
    call trace%end
  end if
  !$omp end critical (data_io_cr)
  END SUBROUTINE init

  !> =============================================================================
  !> Register one type of task, for output into separate file
  !> =============================================================================
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

!> =============================================================================
!> Write results to disk
!> =============================================================================
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
  integer:: it, count, np
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin('data_io_t%output', itimer=itimer)
  call self%init (patch)
  np = (time_slices%order+1)/2
  it = patch%iit(patch%nt-np)
  if (patch%t(it) >= patch%out_next) then
    !$omp critical (data_io_cr)
    do while (patch%t(it) >= patch%out_next)
      if (created < patch%iout) then
        created = patch%iout
        call os%mkdir (trim(io%outputname))
      end if
      call os%mkdir_no (patch%iout)
      !---------------------------------------------------------------------------
      count = counters%decrement (patch%iout+1, io%nwrite)
      if (io%verbose > 0) &
        write (io_unit%log,*) 'count, nwrite =', count, io%nwrite
      append = (count < io%nwrite-1)
      !---------------------------------------------------------------------------
      if      (io%method=='direct') then
        call direct_io%output (patch)
      else if (io%method=='legacy') then
        call legacy_io%output (patch)
        call output_txt (patch, append)
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
      !---------------------------------------------------------------------------
      if (count==0 .and. .not.io_unit%do_validate) then
        if (io%out_time > 0d0) then
          write(io%output,'(1x,a,2x,a,i5.5,5x,"time:",1p,2g14.5)') &
            trim(io%method)//'_io:output', &
            'snapshot '//trim(io%outputname), patch%iout, patch%out_next, patch%time
        else
          write(io%output,'(1x,a,2x,a,i5.5,5x,"time:",1p,g14.5)') &
            trim(io%method)//'_io:output', &
            'snapshot '//trim(io%outputname), patch%iout, patch%time
        end if
      end if
      !---------------------------------------------------------------------------
      call output_nml (patch, append)
      patch%iout = patch%iout+1
      if (io%out_time==0.0) exit
      patch%out_next = (nint(patch%out_next/io%out_time)+1)*io%out_time
      call os%mkdir_no (patch%iout)
    end do
    !$omp end critical (data_io_cr)
  end if
  call trace%end (itimer)
END SUBROUTINE output

!> =============================================================================
SUBROUTINE open_and_close ()
  !-----------------------------------------------------------------------------
  if (io%method == 'snapshot') then
    call parallel_io%open_and_close ()
  end if
END SUBROUTINE open_and_close 

!> =============================================================================
!> Read snapshot from disk
!> =============================================================================
SUBROUTINE input (self, patch, ok)
  class(data_io_t):: self
  class(patch_t):: patch
  logical:: ok
  integer:: count
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input')
  call self%init (patch)
  !$omp critical (data_io_cr)
  ok = .false.
  if (io%restart > 0) then
    if      (io%method=='direct') then
      call direct_io%input (patch, ok)
    else if (io%method=='legacy') then
      call legacy_io%input (patch, ok)
    else if (io%method=='pan') then
      call parallel_io%input (patch, ok)
      call input_nml (patch)
    else if (io%method=='parallel') then
      call parallel_io%input (patch, ok)
      call input_nml (patch)
    else if (io%method=='snapshot') then
      call parallel_io%input (patch, ok)
      call input_nml (patch)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
  end if
  if (ok) then
    count = counters%decrement (io%restart+1, io%nwrite)
    if (count==0 .and. io%master) then
      write(io%output,'(1x,a,2x,a,i5.5,5x,"time:",1p,g14.5)') &
        trim(io%method)//'_io:input', &
        'snapshot '//trim(io%outputname), io%restart, patch%time
    end if
    patch%t(patch%it) = patch%time
  end if
  !$omp end critical (data_io_cr)
  call trace%end
END SUBROUTINE input

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
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
     write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%iout, patch%id
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
  call trace%end
END SUBROUTINE output_txt

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
SUBROUTINE output_nml (patch, append)
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
  character(len=16)    :: kind, eos, opacity, method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land
  namelist /patch_nml/ id, position, dtime
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml')
  !$omp critical (output_cr)
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
  if (mpi%master .and. .not.append) then
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
  !$omp end critical (output_cr)
  call trace%end
END SUBROUTINE output_nml

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
SUBROUTINE input_nml (patch)
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
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    no_mans_land
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
    read  (io_unit%nml1, snapshot_nml)
    close (io_unit%nml1)
  end if
  patch%t     = time
  patch%dt    = 0d0
  patch%time  = time
  patch%istep = 0
  patch%out_next = (int(patch%time/io%out_time)+1)*io%out_time
  patch%guard_zones = guard_zones
  patch%time_derivs = time_derivs
  !$omp end critical (output_cr)
  call trace%end
END SUBROUTINE input_nml

!> =============================================================================
!> Copy one text file to another
!> =============================================================================
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

!> =============================================================================
!> Close the data I/O
!> =============================================================================
SUBROUTINE close
  !.............................................................................
  if (.not.io%do_output) return
  call trace%begin('patch_t%close')
  if      (io%method=='legacy') then
    !call legacy_io%close
  else if (io%method=='direct') then
    call direct_io%close
  else if (io%method=='buffered') then
    call buffered_io%close
  else if (io%method=='pan') then
    call parallel_io%close
  else if (io%method=='parallel') then
    call parallel_io%close
  else if (io%method=='snapshot') then
    call parallel_io%close
  end if
  call trace%end
END SUBROUTINE close

END MODULE data_io_mod
