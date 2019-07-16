!===============================================================================
!> Module with which one can register any number of pointers to real or integer
!> arrays, and then output the contents of the links later
!>
!> USE aux_mod
!> real, pointer:: a1(10), a2(11,12), a3(13,14,15)
!> ...
!> call aux%open ('test.aux')
!> call aux%register ('a1', a1)
!> call aux%register ('a2', a2)
!> call aux%register ('a3', a3)
!> ...
!> call aux%output
!> call aux%close
!===============================================================================
MODULE aux_mod
  USE dll_mod
  implicit none
  private
  !
  type, public, extends(dll_t):: aux_t
    integer:: unit
  contains
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
  integer:: uniq=200
  type(aux_t), public:: aux
CONTAINS

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
  allocate (item_t:: item)
  select type (item)
  class is (item_t)
    item%rank = size(shape(c))
    item%name = name
    item%c3 => c
  end select
  call self%append (item)
END SUBROUTINE register3

!===============================================================================
SUBROUTINE output_aux (self)
  class(aux_t):: self
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  item => aux%head
  do while (associated(item))
    select type (item)
    class is (item_t)
    call item%output (self%unit)
    end select
    item => item%next
  end do
END SUBROUTINE

!===============================================================================
SUBROUTINE open (self, file)
  class(aux_t):: self
  character(len=*):: file
  !$omp atomic
  uniq = uniq+1  
  self%unit = uniq
  open (self%unit, file=file, form='unformatted', access='sequential', status='unknown')
  print *,'open '//trim(file)
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
  call aux%open ('test.aux')
  allocate (real:: r1(10), r2(2,3), r3(3,4,5))
  call aux%register ('r1', r1)
  call aux%register ('r2', r2)
  call aux%register ('r3', r3)
  call aux%output
  call aux%close
END SUBROUTINE

END MODULE aux_mod
