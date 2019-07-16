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
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  implicit none
  private
  !
  integer, save:: max_select=30
  type, public, extends(dll_t):: aux_t
    integer:: unit
  contains
    procedure:: init
    procedure:: open
    procedure:: close
    procedure, private:: register1
    procedure, private:: register2
    procedure, private:: register3
    generic:: register => register1, register2, register3
    procedure:: output => output_aux
    procedure:: test
  end type
  !
  type, extends(dll_node_t):: item_t
    class(*), pointer:: c1(:)=>null(), c2(:,:)=>null(), c3(:,:,:)=>null()
    character(len=32):: name
    integer:: rank
    character(2):: type
  contains
    procedure:: output => output_item
  end type
  character(len=32), allocatable:: select(:)
  logical, save:: first_time=.true.
  integer:: uniq=200
  type(aux_t), public:: aux
CONTAINS

!===============================================================================
SUBROUTINE init (self, name)
  class(aux_t):: self
  character(len=*), optional:: name
  !.............................................................................
  namelist /aux_params/ select
  integer:: iostat
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      allocate (select(max_select))
      rewind (stdin)
      read (stdin, aux_params, iostat=iostat)
      write (stdout, aux_params)
    end if
    !$omp end critical (input_cr)
  end if
END SUBROUTINE

!===============================================================================
SUBROUTINE register1 (self, name, c)
  class(aux_t)       :: self
  character(len=*)   :: name
  class(*), pointer  :: c(:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  allocate (item_t:: item)
  select type (item)
  class is (item_t)
    item%rank = size(shape(c))
    item%name = name
    item%c1 => c
 end select
  call self%append (item)
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
SUBROUTINE output_aux (self, file)
  class(aux_t):: self
  character(len=*):: file
  class(dll_node_t), pointer:: item
  logical:: ok
  integer:: i
  !-----------------------------------------------------------------------------
  call self%init
  call self%open (file)
  item => aux%head
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
    if (ok) &
      call item%output (self%unit) 
    end select
    item => item%next
  end do
  !$omp atomic write
  first_time = .true.
  call self%close
END SUBROUTINE

!===============================================================================
SUBROUTINE open (self, file)
  class(aux_t):: self
  character(len=*):: file
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%open')
  !$omp atomic
  uniq = uniq+1  
  self%unit = uniq
  open (self%unit, file=file, form='unformatted', access='sequential', status='unknown')
  print *,'open '//trim(file)
  call trace%end
END SUBROUTINE

!===============================================================================
SUBROUTINE close (self)
  class(aux_t):: self
  close (self%unit)
END SUBROUTINE

!===============================================================================
SUBROUTINE output_item (self, unit)
  class(item_t):: self
  integer:: unit, rank
  !-----------------------------------------------------------------------------
  print *,'write '//trim(self%name)
  write (unit) self%type, self%rank, self%name
  rank = self%rank
  select case (rank)
  case(1)
    call write1 (self%c1)
  case(2)
    call write2 (self%c2)
  case(3)
    call write3 (self%c3)
  end select
  flush (unit)
!===============================================================================
contains
subroutine write1 (c1)
  class(*), pointer:: c1(:)
  write (unit) shape(c1)
  select type (c1)
    type is (real)
      write (unit) c1
    type is (integer)
      write (unit) c1
    end select
end subroutine write1
subroutine write2 (c2)
  class(*), pointer:: c2(:,:)
  write (unit) shape(c2)
  select type (c2)
    type is (real)
      write (unit) c2
    type is (integer)
      write (unit) c2
    end select
end subroutine write2
subroutine write3 (c3)
  class(*), pointer:: c3(:,:,:)
  write (unit) shape(c3)
  select type (c3)
    type is (real)
      write (unit) c3
    type is (integer)
      write (unit) c3
    end select
end subroutine write3
END SUBROUTINE 

!===============================================================================
SUBROUTINE test (self)
  class(aux_t):: self
  class(item_t), pointer:: item
  class(*), pointer:: r1(:), r2(:,:), r3(:,:,:)
  !-----------------------------------------------------------------------------
  allocate (real:: r1(10), r2(2,3), r3(3,4,5))
  call aux%register ('r1', r1)
  call aux%register ('r2', r2)
  call aux%register ('r3', r3)
  call aux%output  ('test.aux')
END SUBROUTINE

END MODULE aux_mod
