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

!> =============================================================================
!> Write results to disk, placing the files belonging to each snapshot in a
!> separate directory, to avoid the slowdown sometimes associated with having
!> many files in the same directory.
!> =============================================================================
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
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: buf
  real(8):: time
 ! integer, save:: iodir=-1
  !-----------------------------------------------------------------------------
  call trace_begin('patch_t%output')
  !$omp critical (output_cr)
  !-----------------------------------------------------------------------------
  ! Check number of patches written and issue a message when all done
  !-----------------------------------------------------------------------------
  if (io%nwrite<=0) then
    io%nwrite = io%ntask
  end if
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (io%time_derivs) then
    no = min(2,self%nw)
  else
    no = 1
  end if
  if (io%out_time==0.0) then
    time = self%time
  else if (io%time_derivs) then
    time = self%t(self%nt-2)
  else
    time = self%out_next
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
  write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout+1
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
      it1 = self%iit(self%nt-2)
      it2 = self%iit(self%nt-1)
      if (io%out_time == 0.0) then
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,jo)
      !-------------------------------------------------------------------------
      ! Normally, we write out only variable values, which can then be interpolated
      ! to the exact output time, self%out_next, the update of which must come
      ! after this point
      !-------------------------------------------------------------------------
      else
        pt = (time-self%t(it1))/max(self%t(it2)-self%t(it1),1d-30)
        if (self%t(it2)==self%t(it1)) pt=0.0
        qt = 1.0-pt
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)*qt &
            + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
      end if
      call check_error(self, buf, "output 2")
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it1,jo,minval(buf),maxval(buf),pt
      write (io_unit%data,rec=iv+(jo-1)*self%nv) buf
      if (io%verbose>1 .and. iv==self%idx%d) print *,'rho minmax',minval(buf), maxval(buf)
      deallocate (buf)
    end do
  end do
  close (io_unit%data)
  !$omp end critical (output_cr)
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> Check for y-redundancy
!> =============================================================================
SUBROUTINE check_error (patch, a, label, error)
  class(patch_t):: patch
  character(len=*), optional:: label
  integer, optional:: error(3)
  integer:: loc(3), it, l(3), u(3), n(3), iy
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: a, b
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
    if (present(label)) print*,'patch_t::error ', label
    print*,'patch_t::error loc        :', loc
    print*,'patch_t::error it,new,iit :', patch%it, patch%new, patch%iit
    do iy=l(2),u(2)
      print*,'patch_t::error iy, a:', iy, a(loc(1),iy,loc(3))
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
      print*,'patch_t::check_error maxval ', label, maxval(abs(a))
    else
      print*,'patch_t::check_error maxval', maxval(abs(a))
    end if
  end if
  deallocate(b)
END SUBROUTINE check_error

!> =============================================================================
!> Read snapshot from disk
!> =============================================================================
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
  !-----------------------------------------------------------------------------
  ! Read the patch time from the .txt file
  !-----------------------------------------------------------------------------
  write (filename,'(a,"/",i4.4,"_",i4.4)') trim(io%inputdir), self%id, self%restart
  inquire (file=trim(filename)//'.txt', exist=ok)
  if (.not.ok) then
    write (filename,'(a,"/",i5.5,"/",i5.5)') trim(io%inputdir), self%restart, self%id
    inquire (file=trim(filename)//'.txt', exist=ok)
    if (.not.ok) then
      if (self%id==io%id_debug) &
        print *,'MK input_legacy: no .txt file', trim(filename)//'.txt'
    end if
  end if
  if (ok) then
    do_print = do_print .and. io%master
    open (io_unit%data, file=trim(filename)//'.txt', form='formatted', status='unknown')
    read (io_unit%data,'(a)') line
    if (line(1:2)=='  ') then
      read (line,*) time
      ioformat = 0
      if (do_print) print '(a,i12)', 'snapshot ioformat:', ioformat
    else
      read (line,*) ioformat
      if (do_print) print '(a,i12)', 'snapshot ioformat:', ioformat
      if (ioformat==1 .or. ioformat==2) then
        fmt = '(a)'
        read (io_unit%data,fmt) line
        if (trim(line) /= trim(self%kind)) error stop 'wrong solver'
        read (io_unit%data,fmt) line
        if (do_print) print '(a,a12)','snapshot      EOS:', trim(line)
        read (io_unit%data,fmt) line
        if (do_print) print '(a,a12)','snapshot       RT:', trim(line)
        read (io_unit%data,*) time
      else
        error stop 'unsupported ioformat'
      end if
    end if
    do_print = .false.
    close (io_unit%data)
    self%time = time
    self%t(self%it) = time
    !---------------------------------------------------------------------------
    ! If guard_zones is true, include the guard zones.
    !---------------------------------------------------------------------------
    l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
    u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
    !---------------------------------------------------------------------------
    ! Read the binary file with the data cube.
    !---------------------------------------------------------------------------
    inquire (file=trim(filename)//'.dat', exist=ok)
    if (.not.ok) then
      if (self%id==io%id_debug) &
        print *,'MK input_legacy: no .dat file', trim(filename)//'.dat'
    end if
    if (ok) then
      open (io_unit%data, file=trim(filename)//'.dat', form='unformatted', &
        access='direct', status='unknown', recl=kind(self%mem)*product(u-l+1))
      do iv=1,self%nv
        read (io_unit%data,rec=iv) self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)
      end do
      if (io%master .and. io%verbose>0 .or. self%track) &
        print*, 'patch_t%input_legacy: ', trim(filename), time
      close (io_unit%data)
      if (self%id==io%id_debug) &
        print *,'MK input_legacy', self%id, ok, self%mem(18,18,36,1,self%it,1)
      !---------------------------------------------------------------------------
      ! Update the iout and out_next for continued experiment
      !---------------------------------------------------------------------------
      self%out_next = (int(self%time/io%out_time)+1)*io%out_time
    end if
  end if
  !$omp end critical (output_cr)
  call trace_end
END SUBROUTINE input

END MODULE legacy_io_mod
