!===============================================================================
!> Extend the nbor_tt data type for guard zones with the AMR level
!===============================================================================
MODULE guard_mod
  USE trace_mod
  USE io_mod
  USE kinds_mod
  USE lagrange_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  !USE face_mod
  implicit none
  private
  type, public, extends(nbor_t):: guard_t
    integer:: l(3)                                ! lower corner
    integer:: u(3)                                ! upper corner
    integer:: o(3)                                ! offset source to target
    integer:: nt, nv
    real, pointer, dimension(:,:,:):: den => null()
    real, pointer, dimension(:,:,:,:):: mem => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: range
    procedure:: same
    procedure:: differ
  end type
  integer:: order=1
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Register the source and target, and initialize the range
!===============================================================================
SUBROUTINE init (self, target)
  class(guard_t):: self
  class(task_t), pointer:: target, task
  integer:: iostat
  namelist /guard_params/ order, verbose
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_t%init')
  call io%lock%set ('guard_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, guard_params, iostat=iostat)
    call io%namelist_warning ('guard_params', iostat)
    write (io%namelist, guard_params)
    first_time = .false.
  end if
  call io%lock%unset ('guard_t%init')
  !-----------------------------------------------------------------------------
  ! Register a pointer to the target task
  !-----------------------------------------------------------------------------
  self%target => target                       ! cache source pointer
  self%nt = target%nt
  !-----------------------------------------------------------------------------
  ! Register, in the case where the nbor is a patch, a pointer to the patch
  !-----------------------------------------------------------------------------
  task => self%task
  select type (task)
  class is (patch_t)
    self%source => task
    self%nv = task%nv
  end select
  !-----------------------------------------------------------------------------
  ! Initialize the index ranges and offsets
  !-----------------------------------------------------------------------------
  call self%range
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Decide which type of guard zone update to use
!===============================================================================
SUBROUTINE update (self)
  class(guard_t):: self
  !-----------------------------------------------------------------------------
  if (self%source%level == self%target%level) then
    call self%same
  else
    call self%differ
  end if
END SUBROUTINE update

!===============================================================================
!> Setup the index range in the source, and the offset to the target
!===============================================================================
SUBROUTINE range (self)
  class(guard_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_t%range')
  call trace%end ()
END SUBROUTINE range

!===============================================================================
!> Update the guard zones covered by this guard data type 
!===============================================================================
SUBROUTINE same (self)
  class(guard_t):: self
  integer:: i, i1, i2, i3, i4, o(3), l(3), u(3), it, jt, jt1, jt2, nw, id
  real(var_k):: mem(self%l(1):self%u(1))
  real(8):: w(self%nt), times(self%nt)
  class(task_t), pointer:: task
  class(patch_t), pointer:: source, target
  integer:: slots(self%nt), or
  !-----------------------------------------------------------------------------
  task => self%source
  select type (task)
  class is (patch_t)
  source => task
  end select
  task => self%target
  select type (task)
  class is (patch_t)
  target => task
  end select
!  select type (source)
!  class is (patch_t)
  call source%timeslots (slots, times)
  call lagrange%interval (target%time, times, jt1, jt2, or, order)
  w(1:or+1) = lagrange%weights (target%time, times(i1:i2))
  nw = or+1
  allocate (self%mem(self%l(1):self%u(1), &
                     self%l(2):self%u(2), &
                     self%l(3):self%u(3),self%nv))
  allocate (self%den(self%l(1):self%u(1), &
                     self%l(2):self%u(2), &
                     self%l(3):self%u(3)))
  l = self%l
  u = self%u
  o = self%o
  !-----------------------------------------------------------------------------
  ! First time slice
  !-----------------------------------------------------------------------------
  jt = jt1
  it = slots(jt)
  call source%mem_lock(it)%set ('guard%update')
  do i4=1,target%nv
    if (source%unsigned(i4)) then
      call add (i4, w(1), set=.true., lg=.true.)
    else if (source%pervolume(i4)) then
      call add (i4, w(1), set=.true., den=.true.)
    else
      call add (i4, w(1), set=.true.)
    end if
  end do
  call source%mem_lock(it)%unset ('guard%update')
  !---------------------------------------------------------------------------
  ! Other time slices
  !---------------------------------------------------------------------------
  do i=2,nw
    jt = jt+1
    it = slots(jt)
    call source%mem_lock(it)%set ('guard%update')
    do i4=1,target%nv
      if (source%unsigned(i4)) then
        call add (i4, w(1), lg=.true.)
      else if (source%pervolume(i4)) then
        call add (i4, w(1), den=.true.)
      else
        call add (i4, w(1))
      end if
    end do
    call source%mem_lock(it)%unset ('guard%update')
  end do
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i4=1,self%nv
    if (source%unsigned(i4)) then
      self%mem(:,i2,i3,i4) = exp(self%mem(:,i2,i3,i4))
    else if (source%pervolume(i4)) then
      self%mem(:,i2,i3,i4) = &
      self%mem(:,i2,i3,i4) * self%den(:,i2,i3)
    end if
    !---------------------------------------------------------------------------
    ! Copy results to target, with offset
    !---------------------------------------------------------------------------
    do i1=l(1),u(1)
      target%mem (i1+self%o(1), i2+self%o(2), i3+self%o(3), i4, it, 1) = &
        self%mem(i1,i2,i3,i4)
    end do
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Wrap up
  !-----------------------------------------------------------------------------
  deallocate (self%mem, self%den)
!  end select
  call trace%end()
  return
contains
!===============================================================================
!> Add or set array valuesø
!===============================================================================
subroutine add (iv, w, lg, set, den)
  integer:: iv
  real(8):: w
  logical, optional:: lg, set, den
  integer:: i1, i2, i3, i4
  !-----------------------------------------------------------------------------
  if (present(set)) then
    if (present(lg)) then
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = w*log(source%mem(i1,i2,i3,i4,it,1))
      end do; end do; end do
    else if (present(den)) then
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = w*source%mem(i1,i2,i3,i4,it,1) &
                                 /source%mem(i1,i2,i3,id,it,1) 
      end do; end do; end do
    else
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = w*source%mem(i1,i2,i3,i4,it,1)
      end do; end do; end do
    end if
    if (present(den) .and. i4==id) then
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%den(i1,i2,i3) = w*source%mem(i1,i2,i3,i4,it,1)
      end do; end do; end do
    end if
  else
    if (present(lg)) then
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = &
        self%mem(i1,i2,i3,i4) + w*log(source%mem(i1,i2,i3,i4,it,1))
      end do; end do; end do
    else if (present(den)) then
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = &
        self%mem(i1,i2,i3,i4) + w*source%mem(i1,i2,i3,i4,it,1) &
                                 /source%mem(i1,i2,i3,id,it,1) 
      end do; end do; end do
    else
      do i3=l(3),u(3); do i2=l(2),u(2); do i1=l(1),u(1)
        self%mem(i1,i2,i3,i4) = &
        self%mem(i1,i2,i3,i4) + w*source%mem(i1,i2,i3,i4,it,1)
      end do; end do; end do
    end if
  end if
end subroutine add
END SUBROUTINE same

!===============================================================================
!> Guard zone updates requiring interpolation
!===============================================================================
SUBROUTINE differ (self)
  class(guard_t):: self
END SUBROUTINE differ

END MODULE guard_mod
