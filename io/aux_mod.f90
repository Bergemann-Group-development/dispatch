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
    integer:: unit=-1
  contains
    procedure:: prepare
    procedure:: new_item
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
    real, pointer:: r1(:)=>null(), r2(:,:)=>null(), r3(:,:,:)=>null(), &
      r4(:,:,:,:)=>null()
    character(len=32):: name
    integer:: rank
  contains
    procedure:: output => output_item
  end type
  integer:: version=1, verbose=0, uniq=200, iout_prv=-1
  logical:: first_time=.true., second_time=.true., on=.false., do_aux=.true.
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
      return
    end if
    end select
    item => item%next
  end do
  if (verbose > 2) &
    write (stdout,*) 'new aux item: ', trim(name)
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
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%register3')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 3
      item%name = name
      item%r3 => r
    end select
    call self%append (item)
    if (verbose > 2) &
      write (stdout,*) self%n, ' items registered'
  end if
  call trace%end()
END SUBROUTINE register3r

!===============================================================================
SUBROUTINE register4r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:,:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%register4')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 4
      item%name = name
      item%r4 => r
    end select
    call self%append (item)
    if (verbose > 2) &
      write (stdout,*) self%n, ' items registered'
  end if
  call trace%end()
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
      item%rank = 1
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
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 2
      item%name = name
      item%c2 => c
    end select
  call self%append (item)
  end if
END SUBROUTINE register2

!===============================================================================
SUBROUTINE output_aux (self, iout, id, singlefile)
  class(aux_t):: self
  integer:: iout, id
  character(len=*), optional:: singlefile
  !.............................................................................
  character(len=64):: filename
  class(dll_node_t), pointer:: item
  logical:: ok
  integer:: i
  !-----------------------------------------------------------------------------
  ! Make sure only one output runs at a time, so there is no interference of
  ! reading stdin and deciding of to write messages only once
  !-----------------------------------------------------------------------------
  if (.not. do_aux) &
    return
  call trace%begin ('aux_t%output')
  !$omp critical (aux_cr)
  call self%prepare (iout)
  if (present(singlefile)) then
    filename = singlefile
    open (self%unit, file=trim(singlefile), form='unformatted', &
      access='sequential', status='unknown')
  else
    write (filename,'(a,i5.5"/",i5.5,".aux")') trim(io%outputname), iout, id
    open (self%unit, file=trim(filename), form='unformatted', &
      access='sequential', status='unknown')
  end if
  if (verbose > 2) &
    write (io_unit%log,*) filename
  write (self%unit) version, id
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
      if (verbose > 0 .or. iout > iout_prv) &
        write (stdout,'(a,i5,2x,a)') 'aux_t%output: snapshot =', &
          iout, trim(item%name)
      call item%output (self%unit) 
    end if
    end select
    item => item%next
  end do
  iout_prv = iout
  close (self%unit)
  do_aux=on
  !$omp end critical (aux_cr)
  call trace%end()
END SUBROUTINE output_aux

!===============================================================================
SUBROUTINE prepare (self, iout)
  class(aux_t):: self
  integer:: iout
  !.............................................................................
  namelist /aux_params/ on, select, verbose
  integer:: iostat, i
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%open')
  if (self%unit < 0) then
    !---------------------------------------------------------------------------
    ! The first time we acquire a unique unit number, and allocate select
    !---------------------------------------------------------------------------
    uniq = uniq+1  
    self%unit = uniq
    if (.not.allocated(select)) then
      allocate (select(max_select))
      select(:) = ''
      rewind (stdin)
      read (stdin, aux_params, iostat=iostat)
      write (stdout,'(a,/,1x,a)') '&AUX_PARAMS', 'SELECT='
      do i=1,max_select
        if (trim(select(i)) /= '') then
          write (stdout,*) "'"//trim(select(i))//"',"
        end if
      end do
      write (stdout,*) '/'
    end if
  else
    !---------------------------------------------------------------------------
    ! Only do this once per new snapshot
    !---------------------------------------------------------------------------
    if (iout > iout_prv) then
      rewind (stdin)
      read (stdin, aux_params, iostat=iostat)
    end if
  end if
  call trace%end()
END SUBROUTINE prepare

!===============================================================================
SUBROUTINE output_item (self, unit)
  class(item_t):: self
  integer:: unit, rank
  !-----------------------------------------------------------------------------
  call trace%begin ('item_t%output_item')
  write (unit) self%name                                      ! write: name
  write (unit) self%rank                                      ! write: rank
  rank = self%rank
  select case (rank)
  case(1)
    call write1 (self%r1)
  case(2)
    call write2 (self%r2)
  case(3)
    if (associated (self%r3)) then
      if (verbose > 1) then
        write (stdout,*) trim(self%name), ' shape:', shape(self%r3)
        flush (stdout)
      end if
      call write3 (unit, self%r3)                             ! write: value
    else
      write (stderr,*) 'self%r3 pointer not associated for ', trim(self%name)
    end if
  case(4)
    if (associated (self%r4)) then
      if (verbose > 1) then
        write (stdout,*) trim(self%name), ' shape:', shape(self%r4)
        flush (stdout)
      end if
      call write4 (unit, self%r4)
    else
      write (stderr,*) 'self%r4 pointer not associated for ', trim(self%name)
    end if
  end select
  flush (unit)
  call trace%end()
!===============================================================================
contains
subroutine write1 (c)
  real, pointer:: c(:)
  write (unit) shape(c)
  write (unit) 'r1'
  write (unit) c
end subroutine write1
subroutine write2 (c)
  real, pointer:: c(:,:)
  write (unit) shape(c)
  write (unit) 'r2'
  write (unit) c
end subroutine write2
END SUBROUTINE 

!===============================================================================
subroutine write3 (unit, r)
  integer:: unit
  real, pointer:: r(:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write3')
  write (unit) shape(r)
  if (verbose > 2) write (stdout,*) 'r3', shape(r)
  write (unit) 'r3'
  write (unit) r
  call trace%end()
end subroutine write3

subroutine write4 (unit, r)
  integer:: unit
  real, pointer:: r(:,:,:,:)
  call trace%begin ('write4')
  write (unit) shape(r)
  if (verbose > 2) write (stdout,*) 'r4', shape(r)
  write (unit) 'r4'
  write (unit) r
  call trace%end()
end subroutine write4

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
