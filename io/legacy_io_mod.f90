MODULE legacy_io_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE kinds_mod
  implicit none
  private
  type, public:: legacy_io_t
  contains
    procedure, nopass:: output => output
    procedure, nopass:: input => input
    procedure, nopass:: check_error
  end type
  type(legacy_io_t), public:: legacy_io
CONTAINS

!===============================================================================
!> Write results to disk, placing the files belonging to each snapshot in a
!> separate directory, to avoid the slowdown sometimes associated with having
!> many files in the same directory.
!===============================================================================
SUBROUTINE output (patch, experiment_name)
  class(patch_t):: patch
  class(task_t), pointer:: task
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3), data_unit
  real:: pt, qt
  real(4), dimension(:,:,:,:), pointer:: buf
  real(8):: time, t1, t2
 ! integer, save:: iodir=-1
  !-----------------------------------------------------------------------------
  call trace%begin('legacy_io%output')
  !omp critical (output_cr)
  !-----------------------------------------------------------------------------
  ! Check number of patches written and issue a message when all done
  !-----------------------------------------------------------------------------
  ! = patch%max_files-1
  if (patch%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is >0
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    l = patch%mesh%lb
    u = patch%mesh%ub
  else
    l = patch%mesh%li
    u = patch%mesh%ui
  end if
  if (io%time_derivs>0) then
    no = min(2,patch%nw)
  else
    no = 1
  end if
  if (io%out_time==0.0) then
    !$omp atomic read
    time = patch%time
    !$omp end atomic
  else if (io%time_derivs>0) then
    !$omp atomic read
    time = patch%t(patch%nt-2)
    !$omp end atomic
  else
    !$omp atomic read
    time = patch%out_next
    !$omp end atomic
  end if
  !-----------------------------------------------------------------------------
  ! Set io%format
  !-----------------------------------------------------------------------------
  io%format = 1
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    call patch%stats (patch%mem(:,:,:,patch%idx%d,patch%it,1), 'legacy')
  else if (io%verbose > 2) then
    call patch%stats ('legacy')
  end if
  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  if (patch%id < 100000) then
    write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%iout, patch%id
  else if (patch%id < 1000000) then
    write (filename,'(a,i5.5,"/",i6.6)') trim(io%outputname), patch%iout, patch%id
  else
    write (filename,'(a,i5.5,"/",i7.7)') trim(io%outputname), patch%iout, patch%id
  end if
  n = u-l+1
  !-----------------------------------------------------------------------------
  ! If we write out time derivatives we want a pair of consistent values
  ! and time derivatives, which we have at iit(nt-2); time will select that.
  !-----------------------------------------------------------------------------
  !$omp atomic read
  it1 = patch%iit(patch%nt-2)
  !$omp end atomic
  !$omp atomic read
  it2 = patch%iit(patch%nt-1)
  !$omp end atomic
  !-----------------------------------------------------------------------------
  l = merge (patch%mesh%lb, patch%mesh%li, io%guard_zones)
  u = merge (patch%mesh%ub, patch%mesh%ui, io%guard_zones)
  allocate (buf(l(1):u(1),l(2):u(2),l(3):u(3),patch%nv))
  task => patch%link%task
  if (io%do_generic) then
    call task%allocate_generic()
    call task%to_generic ()
    buf = patch%generic
    call task%deallocate_generic()
  else
    buf = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),1:patch%mv,patch%it,1)
  end if
  open (newunit=data_unit, file=trim(filename)//'.dat', form='unformatted', &
        access='direct', status='unknown', recl=kind(buf)*product(shape(buf)))
  write (data_unit,rec=1) buf
  if (io%verbose>1 .and. iv==patch%idx%d) &
    print *,'rho minmax',minval(buf), maxval(buf)
  close (data_unit)
  deallocate (buf)
  !.............................................................................
  if (io%verbose > 0) then
    call patch%stats (buf, 'background_io')
  end if
  !.............................................................................
  !omp end critical (output_cr)
  call trace%end()
END SUBROUTINE output

!===============================================================================
!> Read snapshot from disk
!===============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  class(task_t), pointer:: task
  logical:: ok
  character(len=64) filename, line, fmt
  integer:: iv, l(3), u(3), ioformat, recl
  real(8):: time
  real, pointer:: buffer(:,:,:,:)
  logical, save:: do_print=.true.
  !-----------------------------------------------------------------------------
  call trace%begin('legacy_io_t%input')
  ok = .false.
  if (patch%id==io%id_debug) &
    print *,'MK input', patch%id, patch%restart
  if (patch%restart < 0) then
    call trace%end ()
    return
  end if
  !$omp critical (output_cr)
  !---------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.
  !---------------------------------------------------------------------------
  l = merge (patch%mesh%lb, patch%mesh%li, io%guard_zones)
  u = merge (patch%mesh%ub, patch%mesh%ui, io%guard_zones)
  !---------------------------------------------------------------------------
  ! Read the binary file with the data cube.
  !---------------------------------------------------------------------------
  if (patch%id < 100000) then
    write (filename,'(a,i5.5,"/",i5.5)') trim(io%inputdir), patch%restart, patch%id
  else if (patch%id < 1000000) then
    write (filename,'(a,i5.5,"/",i6.6)') trim(io%inputdir), patch%restart, patch%id
  else
    write (filename,'(a,i5.5,"/",i7.7)') trim(io%inputdir), patch%restart, patch%id
  end if
  inquire (file=trim(filename)//'.dat', exist=ok)
  if (.not.ok) then
    call io%abort('MK input_legacy: no .dat file'//trim(filename)//'.dat')
  end if
  if (ok) then
    if (io%verbose > 0) &
      write (stdout,*) 'legacy_io_t%input: id, restart, it', &
        patch%id, patch%restart, patch%it
    !...........................................................................
    buffer => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),1:patch%mv,patch%it,1)
    recl = kind(buffer)*product(shape(buffer))
    open (io_unit%data, file=trim(filename)//'.dat', form='unformatted', &
      access='direct', status='unknown', recl=recl)    
    read (io_unit%data,rec=1) buffer
    close (io_unit%data)
    if (io%in%do_generic) then
      task => patch%link%task
      call task%allocate_generic()
      patch%generic = buffer
      call task%from_generic ()
      call patch%deallocate_generic()
    end if
    !...........................................................................
    if (io%master .and. io%verbose>0 .or. patch%track) &
      print*, 'legacy_io%input_legacy: ', trim(filename), time
    if (patch%id==io%id_debug) &
      print *,'MK input_legacy', patch%id, ok, patch%mem(18,18,36,1,patch%it,1)
    !---------------------------------------------------------------------------
    ! Update the out_next for continued experiment
    !---------------------------------------------------------------------------
    patch%out_next = (nint(patch%time/io%out_time)+1)*io%out_time
    !-----------------------------------------------------------------------------
  end if
  !$omp end critical (output_cr)
  call trace%end ()
END SUBROUTINE input

!===============================================================================
!> Check for y-redundancy
!===============================================================================
SUBROUTINE check_error (patch, a, label, error)
  class(patch_t):: patch
  character(len=*), optional:: label
  integer, optional:: error(3)
  integer:: loc(3), it, l(3), u(3), n(3), iy
  real(4), dimension(:,:,:), pointer:: a, b
  integer, save:: nerror=4
  !-----------------------------------------------------------------------------
  if (.not.io%do_debug) return
  if (io%guard_zones) then
    l = patch%mesh%lb
    u = patch%mesh%ub
  else
    l = patch%mesh%li
    u = patch%mesh%ui
  end if
  n = u-l+1
  allocate (b(n(1),n(2),n(3)))
  do iy=l(2),u(2)
    b(:,iy,:)=a(:,iy,:)-a(:,1,:)
  end do
  if (any(b /= 0.0)) then
    loc = maxloc(abs(b)) + l-1
    if (present(label)) print*,'legacy_io::error ', label
    print*,'legacy_io::error loc        :', loc
    print*,'legacy_io::error it,new,iit :', patch%it, patch%new, patch%iit
    do iy=l(2),u(2)
      print*,'legacy_io::error iy, a:', iy, a(loc(1),iy,loc(3))
    end do
    if (present(error)) then
      error=loc
      return
    else
      if (nerror<=0) stop
      nerror=nerror-1
    end if
  else
    if (present(label)) then
      print*,'legacy_io::check_error maxval ', label, maxval(abs(a))
    else
      print*,'legacy_io::check_error maxval', maxval(abs(a))
    end if
  end if
  deallocate(b)
END SUBROUTINE check_error

END MODULE legacy_io_mod
