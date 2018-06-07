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
  USE mpi_mod
  USE omp_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE direct_io_mod
  USE legacy_io_mod
  USE parallel_io_mod
  USE buffered_io_mod
  USE snapshot_list_mod
  implicit none
  private
  type, public:: data_io_t
    logical:: initialized=.false.
  contains
    procedure:: init
    procedure:: output
    procedure, nopass:: output_txt
    procedure:: input
    procedure, nopass:: close
  end type
  type(data_io_t), public:: data_io
CONTAINS

!> =============================================================================
!> Initialize the I/O method
!> =============================================================================
SUBROUTINE init (self, task)
  class(data_io_t):: self
  class(patch_t):: task
  !.............................................................................
  call trace%begin('data_io_t%init')
  if (io%do_output) then
    select type (task)
    class is (patch_t)
      if      (io%method=='direct') then
        call direct_io%init (task)
      else if (io%method=='parallel') then
        call parallel_io%init (task)
      end if
    end select
  end if
  self%initialized = .true.
  call trace%end
END SUBROUTINE init

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
  type(snapshot_t), pointer:: snapshot
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin('data_io_t%output', itimer=itimer)
  do while (patch%time >= patch%out_next)
    call io%mkdir (trim(io%outputname))
    call io%mkdir_no (patch%iout)
    call snapshot_list%count (patch, snapshot)
    if      (io%method=='direct') then
      call direct_io%output (patch)
    else if (io%method=='legacy') then
      call legacy_io%output (patch)
      call output_txt (patch, snapshot)
    else if (io%method=='buffered') then
      call buffered_io%output (patch)
    else if (io%method=='parallel') then
      io%format = 5
      call parallel_io%output (patch)
      call output_txt (patch, snapshot)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
    if (io%out_time==0.0) exit
    patch%iout = patch%iout+1
    patch%out_next = (nint(patch%out_next/io%out_time)+1)*io%out_time
    call io%mkdir_no (patch%iout)
  end do
  call trace%end (itimer)
END SUBROUTINE output

!> =============================================================================
!> Read snapshot from disk
!> =============================================================================
SUBROUTINE input (self, patch, ok)
  class(data_io_t):: self
  class(patch_t):: patch
  logical:: ok
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input')
  if (.not.self%initialized) then
    !$omp critical (data_io_cr)
    if (.not.self%initialized) then
      call self%init (patch)
    end if
    !$omp end critical (data_io_cr)
  end if
  ok = .false.
  if (patch%restart > 0) then
    if      (io%method=='direct') then
      call direct_io%input (patch, ok)
    else if (io%method=='legacy') then
      call legacy_io%input (patch, ok)
    else if (io%method=='parallel') then
      call parallel_io%input (patch, ok)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
  end if
  if (ok .and. io%verbose>2) print *, patch%id, ' LEGACY INPUT'
  call trace%end
END SUBROUTINE input

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
SUBROUTINE output_txt (patch, snapshot)
  class(patch_t):: patch
  real(8):: time
  integer:: l(3), u(3)
  type(snapshot_t), pointer:: snapshot
  !.............................................................................
  character(len=128) filename
  integer:: nlines=6
  logical:: append
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_txt')
  !$omp critical (output_cr)
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  append = .false.
  if (io%format > 2) then
     write (filename,'(a,i5.5,"/",i5.5,".rank")') trim(io%outputname), patch%iout, mpi%rank
     if (associated(snapshot)) then
       if (snapshot%ibuf < snapshot%nbuf-1) append=.true.
     else
       append = .true.
     end if
  else
     write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%iout, patch%id
  end if
  if (io%debug(2)) &
    print'(a,3x,"time:",2f13.6)', ' data_io_t%output_txt: '//trim(filename), &
    patch%out_next, patch%time
  !-----------------------------------------------------------------------------
  ! Open the file, appending or not
  !-----------------------------------------------------------------------------
  if (append) then
    open  (io_unit%data, file=trim(filename)//'.txt', form='formatted', &
           status='unknown', access='append')
  else
    open  (io_unit%data, file=trim(filename)//'.txt', form='formatted', &
           status='unknown')
    write (io_unit%data,'(i2.2)') io%format
  end if
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
!> Close the data I/O
!> =============================================================================
SUBROUTINE close
  !.............................................................................
  if (.not.io%do_output) return
  call trace%begin('patch_t%close')
  if      (io%method=='legacy') then
    !call legacy_io%close
  else if (io%method=='direct') then
    !call direct_io%close
  else if (io%method=='buffered') then
    call buffered_io%close
  else if (io%method=='parallel') then
    call parallel_io%close
  end if
  call trace%end
END SUBROUTINE close

END MODULE data_io_mod
