!> =============================================================================
!> Wrapper module, which seamlessly incorporates RT into the stagger2/ solver.
!> =============================================================================
MODULE mhd_rt_mod
  USE io_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  USE rt_data_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_data_t), pointer:: rt_data => null()
    integer:: prev_iout=0
  contains
    procedure :: init => init_rt
    procedure :: update
    procedure :: output_legacy => rt_output_legacy
    procedure :: output => rt_output
  end type
  ! ----------------------------------------------------------------------------
CONTAINS

!> =============================================================================
!> Initialize the underlying MHD patch & allocate the extra memory needed for RT
!> =============================================================================
SUBROUTINE init_rt (self)
  class (mhd_rt_t) :: self
  integer :: nt = 5
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%init')
  self%type = 'mhd_rt_t'
  ! --- initialise heating tables ---
  call self%mhd_t%init
  call trace%end
END SUBROUTINE init_rt

!> =============================================================================
!> Update the MHD state and add the RT head exchange
!> =============================================================================
SUBROUTINE update (self)
  class (mhd_rt_t) :: self
  real, dimension(:,:,:), pointer:: q
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%update')
  ! --- Update MHD, whatever it means ---
  call self%mhd_t%update
  ! --- Update dsdt by adding Q ---
  q => self%rt_data%q(:,:,:,self%it)  ! FIXME: should be interpolation in time
  call self%mhd_t%apply_heating (q)
  call self%output
  call trace%end
END SUBROUTINE update


!> =====================================================================
!> Derived from `patch_t%output_legacy`; extended to output RT quantities
!> not present in `patch_t%mem`.
!> =====================================================================
SUBROUTINE rt_output_legacy (self, experiment_name)
  USE iso_fortran_env, only: real32
  class(mhd_rt_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_rt_mod::output_legacy')
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (self%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (self%time_derivs) then
    no = 2
  else
    no = 1
  end if

  !-----------------------------------------------------------------------------
  ! Write a log file, with files written by this rank
  !-----------------------------------------------------------------------------
  write (logname,'(a)') trim(io%outputname)//'log'
  if (first_time) then
    first_time = .false.
    if (io%master) print*, 'output data file log written to '//trim(logname)
    open (io%data_unit,file=trim(logname),status='unknown')
  else
    open (io%data_unit,file=trim(logname),status='unknown',position='append')
  end if
  write (filename,'(a,i4.4,"_",i4.4)') trim(io%outputname), self%id, self%iout-1
  write (io%data_unit,'(a)') trim(filename)
  ds = self%box/2.0**io%levelmax
  do i=1,3
    if (self%size(i) > 0.0) then
      li(i) = (self%position(i)+0.5*(self%box(i)-self%size(i)))/ds(i)
      ui(i) = (self%position(i)+0.5*(self%box(i)+self%size(i)))/ds(i) - 1.0
    end if
  end do
  write (io%data_unit,'(1p,e14.6,3(2x,3e14.6),1x,i3.1)') self%time, self%position, self%ds, &
        self%velocity, self%level
  write (io%data_unit,'(4(3i7,2x))') int(self%box/self%ds+0.5), u-l+1, li, ui
  close (io%data_unit)

  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  if (io%master .and. io%verbose>=0) &
    print'(a,3x,"time:",f9.5)', 'mhd_rt_mod::output: file= '//trim(filename), self%out_next
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  write (io%data_unit,'(1p,e18.10,6(2x,3e18.10),1x,i3.1)') self%time, self%position, self%ds, &
        self%velocity, self%llc_nat, self%llc_cart, self%mesh%centre_nat, self%level
  write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1)') int(self%box/self%ds+0.5), u-l+1, &
        li, ui, self%n, self%gn, self%ng, self%nv+1 ! +1 for mean intensity
  write (io%data_unit,'(3(1i8.1,1x),i2.1,1x,a)') self%id, self%iout, self%istep, &
        self%mesh_type, trim(self%kind)
  close (io%data_unit)

  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  n = u-l+1
  open (io%data_unit, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl=real32*product(n))
  do jo=1,no
    do iv=1,self%nv
      allocate (buf(n(1),n(2),n(3)))
      !-------------------------------------------------------------------------
      ! If we write out time derivatives we want a pair of consistent values
      ! and time derivatives, which we have at iit(nt-2).
      !-------------------------------------------------------------------------
      it1 = self%iit(self%nt-2)
      if (self%time_derivs) then
        it = it1
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)
      !-------------------------------------------------------------------------
      ! Normally, we write out only variable values, which can then be interpolated
      ! to the exact output time, self%out_next, the update of which must come
      ! after this point
      !-------------------------------------------------------------------------
      else
        it2 = self%iit(self%nt-1)
        it = it2
        pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1d-30)
        if (self%t(it2)==self%t(it1)) pt=0.0
        qt = 1.0-pt
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)*qt &
            + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
      end if
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*self%nv) buf
      deallocate (buf)
    end do
    !--------------------------------------------------------------------------
    ! integrated rt quantities
    allocate (buf(n(1),n(2),n(3)))
    it1 = self%iit(self%nt-2)
    it2 = self%iit(self%nt-1)
    pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1d-30)
    if (self%t(it2)==self%t(it1)) pt=0.0
    qt = 1.0-pt
    buf = self%rt_data%intensity(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it)
    if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*self%nv) buf
    deallocate (buf)
  end do
  close (io%data_unit)

  call trace_end
END SUBROUTINE rt_output_legacy

!> =============================================================================
SUBROUTINE rt_output (self, experiment_name)
  class(mhd_rt_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  !$omp critical (out_cr)
  call trace_begin('mhd_rt_mod::output')
  if (self%prev_iout < self%iout) then
    call self%output_legacy
  end if
  self%prev_iout = self%iout
  call trace_end
  !$omp end critical (out_cr)
END SUBROUTINE rt_output


END MODULE mhd_rt_mod
