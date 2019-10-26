MODULE legacy_io_mod
  USE io_mod
  USE trace_mod
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
SUBROUTINE output (self, experiment_name)
  class(patch_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real(4), dimension(:,:,:), pointer:: buf
  real(8):: time, t1, t2
 ! integer, save:: iodir=-1
  !-----------------------------------------------------------------------------
  call trace%begin('legacy_io%output')
  !$omp critical (output_cr)
  !-----------------------------------------------------------------------------
  ! Check number of patches written and issue a message when all done
  !-----------------------------------------------------------------------------
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is >0
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (io%time_derivs>0) then
    no = min(2,self%nw)
  else
    no = 1
  end if
  if (io%out_time==0.0) then
    !$omp atomic read
    time = self%time
    !$omp end atomic
  else if (io%time_derivs>0) then
    !$omp atomic read
    time = self%t(self%nt-2)
    !$omp end atomic
  else
    !$omp atomic read
    time = self%out_next
    !$omp end atomic
  end if
  !-----------------------------------------------------------------------------
  ! Set io%format
  !-----------------------------------------------------------------------------
  if (self%no_mans_land) then
    io%format = 1
  else
    io%format = 2
  end if
  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  if (self%id < 100000) then
    write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), self%iout, self%id
  else if (self%id < 1000000) then
    write (filename,'(a,i5.5,"/",i6.6)') trim(io%outputname), self%iout, self%id
  else
    write (filename,'(a,i5.5,"/",i7.7)') trim(io%outputname), self%iout, self%id
  end if
  n = u-l+1
  open (io_unit%data, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl=kind(buf)*product(n))
  do jo=1,no
    do iv=1,self%nv
      allocate (buf(n(1),n(2),n(3)))
      !-------------------------------------------------------------------------
      ! If we write out time derivatives we want a pair of consistent values
      ! and time derivatives, which we have at iit(nt-2); time will select that.
      !-------------------------------------------------------------------------
      !$omp atomic read
      it1 = self%iit(self%nt-2)
      !$omp end atomic
      !$omp atomic read
      it2 = self%iit(self%nt-1)
      !$omp end atomic
      if (io%out_time == 0.0) then
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,jo)
      !-------------------------------------------------------------------------
      ! Normally, we write out only variable values, which can then be interpolated
      ! to the exact output time, self%out_next, the update of which must come
      ! after this point
      !-------------------------------------------------------------------------
      else
        !$omp atomic read
        t1 = self%t(it1)
        !$omp end atomic
        !$omp atomic read
        t2 = self%t(it2)
        !$omp end atomic
        pt = (time-t1)/max(t2-t1,1d-30)
        if (t2==t1) pt=0.0
        qt = 1.0-pt
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)*qt &
            + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
      end if
      call check_error(self, buf, "output 2")
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it1,jo,minval(buf),maxval(buf),pt
      write (io_unit%data,rec=iv+(jo-1)*self%nv) buf
      if (io%verbose>1 .and. iv==self%idx%d) &
        print *,'rho minmax',minval(buf), maxval(buf)
      deallocate (buf)
    end do
  end do
  close (io_unit%data)
  !$omp end critical (output_cr)
  call trace%end()
END SUBROUTINE output

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

!===============================================================================
!> Read snapshot from disk
!===============================================================================
SUBROUTINE input (self, ok)
  class(patch_t):: self
  logical:: ok
  character(len=64) filename, line, fmt
  integer:: iv, l(3), u(3), ioformat
  real(8):: time
  logical, save:: do_print=.true.
  !-----------------------------------------------------------------------------
  ok = .false.
  if (self%id==io%id_debug) &
    print *,'MK input', self%id, self%restart
  if (self%restart < 0) return
  call trace_begin('legacy_io_t%input')
  !$omp critical (output_cr)
  !---------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.
  !---------------------------------------------------------------------------
  l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
  u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
  !---------------------------------------------------------------------------
  ! Read the binary file with the data cube.
  !---------------------------------------------------------------------------
  if (self%id < 100000) then
    write (filename,'(a,i5.5,"/",i5.5)') trim(io%inputdir), self%restart, self%id
  else if (self%id < 1000000) then
    write (filename,'(a,i5.5,"/",i6.6)') trim(io%inputdir), self%restart, self%id
  else
    write (filename,'(a,i5.5,"/",i7.7)') trim(io%inputdir), self%restart, self%id
  end if
  inquire (file=trim(filename)//'.dat', exist=ok)
  if (.not.ok) then
    call io%abort('MK input_legacy: no .dat file'//trim(filename)//'.dat')
  end if
  if (ok) then
    open (io_unit%data, file=trim(filename)//'.dat', form='unformatted', &
      access='direct', status='unknown', recl=kind(self%mem)*product(u-l+1))
    do iv=1,self%nv
      read (io_unit%data,rec=iv) self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)
    end do
    if (io%master .and. io%verbose>0 .or. self%track) &
      print*, 'legacy_io%input_legacy: ', trim(filename), time
    close (io_unit%data)
    if (self%id==io%id_debug) &
      print *,'MK input_legacy', self%id, ok, self%mem(18,18,36,1,self%it,1)
    !---------------------------------------------------------------------------
    ! Update the iout and out_next for continued experiment
    !---------------------------------------------------------------------------
    self%out_next = (int(self%time/io%out_time)+1)*io%out_time
  end if
  !$omp end critical (output_cr)
  call trace_end
END SUBROUTINE input

END MODULE legacy_io_mod
