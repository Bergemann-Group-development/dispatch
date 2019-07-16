!===============================================================================
!> Module with which one can register any number of pointers to real or integer
!> arrays, and then output the contents of the links later
!>
!> USE aux_mod
!> real, pointer:: a1(10), a2(11,12), a3(13,14,15)
!> ...
!> call aux%register ('a1', a1)
!> call aux%register ('a2', a2)
!> call aux%register ('a3', a3)
!> ...
!> call aux%output ('test.aux')
!===============================================================================
MODULE aux_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  implicit none
  private
  !
  integer, save:: max_select=30
  type, public, extends(dll_t):: aux_t
    integer:: unit
    integer:: iout_prv=-1
  contains
    procedure:: prepare
    procedure:: new_item
    procedure, private:: register1
    procedure, private:: register2
    procedure, private:: register3
    procedure, private:: register4
    procedure, private:: register1r
    procedure, private:: register2r
    procedure, private:: register3r
    procedure, private:: register4r
    generic:: register => register1r, register2r, register3r, register4r
    procedure:: output => output_aux
    procedure:: test
  end type
  !
  type, extends(dll_node_t):: item_t
    class(*), pointer:: c1(:)=>null(), c2(:,:)=>null(), c3(:,:,:)=>null(), &
      c4(:,:,:,:)=>null()
    character(len=32):: name
    integer:: rank
  contains
    procedure:: output => output_item
  end type
  integer:: format=1, verbose=0, uniq=200
  logical:: first_time=.true., second_time=.true.
  character(len=32), allocatable:: select(:)
  !type(aux_t), public:: aux
CONTAINS

!===============================================================================
FUNCTION new_item (self, name)
  class(aux_t)     :: self
  character(len=*) :: name
  logical          :: new_item
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  new_item = .true.
  item => self%head
  do while (associated(item))
    select type (item)
    class is (item_t)
    if (trim(item%name) == trim(name)) then
      new_item = .false.
      exit
    end if
    end select
    item => item%next
  end do
  if (verbose > 0) &
    write (stdout,*) 'new aux item:', trim(name)
END FUNCTION

!===============================================================================
SUBROUTINE register1r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:)
  class(*), pointer  :: c(:)
  c => r
  call register1 (self, name, c)
END SUBROUTINE
!===============================================================================
SUBROUTINE register2r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:)
  class(*), pointer  :: c(:,:)
  c => r
  call register2 (self, name, c)
END SUBROUTINE
!===============================================================================
SUBROUTINE register3r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:,:)
  class(*), pointer  :: c(:,:,:)
  c => r
  call register3 (self, name, c)
END SUBROUTINE
!===============================================================================
SUBROUTINE register4r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:,:,:)
  class(*), pointer  :: c(:,:,:,:)
  c => r
  call register4 (self, name, c)
END SUBROUTINE

!===============================================================================
SUBROUTINE register1 (self, name, c)
  class(aux_t)       :: self
  character(len=*)   :: name
  class(*), pointer  :: c(:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = size(shape(c))
      item%name = name
      item%c1 => c
    end select
  call self%append (item)
  end if
END SUBROUTINE register1

!===============================================================================
SUBROUTINE register2 (self, name, c)
  class(aux_t)       :: self
  character(len=*)   :: name
  class(*), pointer  :: c(:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  allocate (item_t:: item)
  select type (item)
  class is (item_t)
    item%rank = size(shape(c))
    item%name = name
    item%c2 => c
  end select
  call self%append (item)
END SUBROUTINE register2

!===============================================================================
SUBROUTINE register3 (self, name, c)
  class(aux_t)       :: self
  character(len=*)   :: name
  class(*), pointer  :: c(:,:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%register3')
  allocate (item_t:: item)
  select type (item)
  class is (item_t)
    item%rank = size(shape(c))
    item%name = name
    item%c3 => c
  end select
  call self%append (item)
  call trace%end
END SUBROUTINE register3
!===============================================================================
SUBROUTINE register4 (self, name, c)
  class(aux_t)       :: self
  character(len=*)   :: name
  class(*), pointer  :: c(:,:,:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%register4')
  allocate (item_t:: item)
  select type (item)
  class is (item_t)
    item%rank = size(shape(c))
    item%name = name
    item%c4 => c
  end select
  call self%append (item)
  call trace%end
END SUBROUTINE register4

!===============================================================================
SUBROUTINE output_aux (self, iout, id, testfile)
  class(aux_t):: self
  integer:: iout, id
  character(len=*), optional:: testfile
  !.............................................................................
  character(len=64):: filename
  class(dll_node_t), pointer:: item
  logical:: ok
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%output')
  call self%prepare 
  if (present(testfile)) then
    filename = testfile
  else
    write (filename,'(a,i5.5"/",i5.5,".aux")') trim(io%outputname), iout, id
  end if
  !$omp critical (input_cr)
  open (self%unit, file=filename, form='unformatted', access='sequential', &
    status='unknown')
  if (io%verbose > 0) &
    write (io_unit%log,*) filename
  write (self%unit) format, id
  item => self%head
  do while (associated(item))
    select type (item)
    class is (item_t)
    ok = .false.
    do i=1,max_select
      if (trim(select(i)) == trim(item%name)) then
        ok = .true.
        exit
      end if
    end do
    if (ok) then
      if (verbose > 0) &
        write (stdout,*) iout, id, 'write '//trim(item%name)
      call item%output (self%unit) 
    end if
    end select
    item => item%next
  end do
  close (self%unit)
  !$omp end critical (input_cr)
END SUBROUTINE

!===============================================================================
SUBROUTINE prepare (self)
  class(aux_t):: self
  !.............................................................................
  namelist /aux_params/ select, verbose
  integer:: iostat, i
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%open')
  if (first_time) then
    !---------------------------------------------------------------------------
    ! The first time we acquire a unique unit number, and allocate select
    !---------------------------------------------------------------------------
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      uniq = uniq+1  
      self%unit = uniq
      allocate (select(max_select))
      select(:) = ''
      rewind (stdin)
      read (stdin, aux_params, iostat=iostat)
      do i=1,max_select
        if (trim(select(i)) /= '') &
          write (stdout,*) 'aux_%select: ', trim(select(i))
      end do
     end if
    !$omp end critical (input_cr)
  else
    !---------------------------------------------------------------------------
    ! Should add logics here to only do this once per new snapshot
    !---------------------------------------------------------------------------
    !$omp critical (input_cr)
    rewind (stdin)
    read (stdin, aux_params, iostat=iostat)
    !$omp end critical (input_cr)
  end if
  call trace%end
END SUBROUTINE prepare

!===============================================================================
SUBROUTINE output_item (self, unit)
  class(item_t):: self
  integer:: unit, rank
  !-----------------------------------------------------------------------------
  write (unit) self%name
  write (unit) self%rank
  rank = self%rank
  select case (rank)
  case(1)
    call write1 (self%c1)
  case(2)
    call write2 (self%c2)
  case(3)
    call write3 (self%c3)
  case(4)
    call write4 (self%c4)
  end select
  flush (unit)
!===============================================================================
contains
subroutine write1 (c)
  class(*), pointer:: c(:)
  write (unit) shape(c)
  select type (c)
    type is (real)
      write (unit) 'r1'
      write (unit) c
    type is (integer)
      write (unit) 'i1'
      write (unit) c
    end select
end subroutine write1
subroutine write2 (c)
  class(*), pointer:: c(:,:)
  write (unit) shape(c)
  select type (c)
    type is (real)
      write (unit) 'r2'
      write (unit) c
    type is (integer)
      write (unit) 'i2'
      write (unit) c
    end select
end subroutine write2
subroutine write3 (c)
  class(*), pointer:: c(:,:,:)
  write (unit) shape(c)
  select type (c)
    type is (real)
      write (unit) 'r3'
      write (unit) c
    type is (integer)
      write (unit) 'i3'
      write (unit) c
    end select
end subroutine write3
subroutine write4 (c)
  class(*), pointer:: c(:,:,:,:)
  write (unit) shape(c)
  select type (c)
    type is (real)
      write (unit) 'r4'
      write (unit) c
    type is (integer)
      write (unit) 'i4'
      write (unit) c
    end select
end subroutine write4
END SUBROUTINE 

!===============================================================================
SUBROUTINE test (self)
  class(aux_t):: self
  class(item_t), pointer:: item
  real, pointer:: r1(:), r2(:,:), r3(:,:,:)
  type(aux_t):: aux
  !-----------------------------------------------------------------------------
  allocate (real:: r1(10), r2(2,3), r3(3,4,5))
  call aux%register ('r1', r1)
  call aux%register ('r2', r2)
  call aux%register ('r3', r3)
  call aux%output (3, 1, 'test.aux')
END SUBROUTINE

END MODULE aux_mod
