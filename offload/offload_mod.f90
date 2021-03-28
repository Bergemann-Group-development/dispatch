!===============================================================================
!> Module handling task updating via offloading controlled by OpenMP directives
!===============================================================================
MODULE offload_mod
  USE trace_mod
  USE io_mod
  USE task_mod
  USE patch_mod
  USE list_mod
  USE link_mod
  USE bits_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE device_mod
  implicit none
  private
  type, public:: offload_t
    logical:: on=.false.
    real, pointer, contiguous:: bunch(:,:,:,:,:)
    integer:: bunch_size=0
    real(8):: time=0d0
    real(8):: next_time=0d0
    class(list_t), pointer:: task_list => null()
    type(list_t):: bunch_list
    type(lock_t):: lock
  contains
    procedure:: init
    procedure:: update
    procedure, private:: copy_to_bunch
    procedure, private:: update_bunch
    procedure, private:: copy_from_bunch
  end type
  real(8):: interval=0.1
  integer:: bunch_size=300
  integer:: verbose=0
  logical:: on=.false.
  type(offload_t), public:: offload
CONTAINS

!===============================================================================
!> Initialize offload handler
!===============================================================================
SUBROUTINE init(self)
  class(offload_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /offload_params/ on, verbose, bunch_size, interval
  !-----------------------------------------------------------------------------
  call io%lock%set ('offload_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, offload_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, offload_params)
    end if
    call io%namelist_warning ('offload_params', iostat)
    write (io%namelist, offload_params)
    !...........................................................................
    call self%lock%init ('offload')
    first_time = .false.
  end if
  self%on = on
  call self%bunch_list%init
  call io%lock%unset ('offload_t%init')
END SUBROUTINE init

!===============================================================================
!> Add a task to the bunch and updare the bunch when full or next_time reached
!===============================================================================
SUBROUTINE update (self, task)
  class(offload_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------  
  call self%lock%set ('offload_t%update')
  call self%bunch_list%append (task, keep=.true.)
  call self%copy_to_bunch (task)
  self%time = wallclock()
  if (self%bunch_size >= bunch_size .or. self%time > self%next_time) then
    if (verbose > 1) then
      write(stderr,'(a,2i5,2f12.6)') &
        'offload_t%update: size,bunch_size,time,next_time=', &
        self%bunch_size, bunch_size, self%time, self%next_time
    end if
    call self%update_bunch 
    call self%copy_from_bunch
    self%next_time = self%time + interval
  end if
  call self%lock%unset ('offload_t%update')
END SUBROUTINE update

!===============================================================================
!> Update all tasks in the bunch via offloading
!===============================================================================
SUBROUTINE update_bunch (self)
  class(offload_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  if (device%type == 'CPU') then
    link => self%bunch_list%head
    do while (associated(link))
      task => link%task
      call task%update ()
      link => link%next
    end do
  else
    call device%update (self%bunch)
  end if
END SUBROUTINE update_bunch

!===============================================================================
!> Copy a time slot to the bunch, adding to the self%bunch_size
!===============================================================================
SUBROUTINE copy_to_bunch (self, task)
  class(offload_t):: self
  class(task_t), pointer:: task
  integer:: m(6)
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
    if (.not.associated(self%bunch)) then
      m = shape (task%mem)
      allocate (self%bunch(m(1),m(2),m(3),m(4),bunch_size))
    end if
    self%bunch_size = self%bunch_size+1
    self%bunch(:,:,:,:,self%bunch_size) = task%mem(:,:,:,:,task%new,1)
  end select
END SUBROUTINE copy_to_bunch

!===============================================================================
!> Copy the updated time slices back into place
!===============================================================================
SUBROUTINE copy_from_bunch (self)
  class(offload_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: i
  !-----------------------------------------------------------------------------
  ! Copy the updated time slots into place and clear the tasks' busy bit
  !-----------------------------------------------------------------------------
  i = 1
  link => self%bunch_list%head
  do while (associated(link))
    task => link%task
    if (.not. device%type=='CPU') then
      select type (task)
      class is (patch_t)
        task%mem(:,:,:,:,task%new,1) = self%bunch(:,:,:,:,i)
      end select
    end if
    call task%clear (bits%busy)
    call self%task_list%check_nbors (link)
    i = i+1
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Reset the task list and bunch storage
  !-----------------------------------------------------------------------------
  call self%bunch_list%reset (keep=.true.)
  self%bunch_size = 0
END SUBROUTINE copy_from_bunch

END MODULE offload_mod
