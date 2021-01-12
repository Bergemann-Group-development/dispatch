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
!>
!> From data_io%output, the call hierarchy is (for each task):
!>
!> aux_t%output_sequence
!>   aux_t%output => output_aux
!>     item_t%output => output_item
!>
!===============================================================================
MODULE aux_mod
  USE os_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  USE omp_timer_mod
  USE omp_mod
!  USE omp_lock_mod
  implicit none
  private
  !
  integer, save:: max_select=30
  type, public, extends(dll_t):: aux_t
    integer:: iout = 0
    real(8):: aux_next = 0.0
!    type(lock_t), pointer:: lock => null()
  contains
    procedure:: init
    procedure:: new_item
    procedure, private:: register0i
    procedure, private:: register1i
    procedure, private:: register2i
    procedure, private:: register0r
    procedure, private:: register1r
    procedure, private:: register2r
    procedure, private:: register3r
    procedure, private:: register4r
    procedure, private:: register2d
    generic:: register => &
      register0r, register1r, register2r, register3r, register4r, &
      register0i, register1i, register2i, register2d
    procedure:: list
    procedure:: output => output_aux
    procedure:: output_sequence
    procedure:: test
  end type
  !
  type, extends(dll_node_t):: item_t
    class(*), pointer:: &
      c0          => null(), &
      c1(:)       => null(), &
      c2(:,:)     => null(), &
      c3(:,:,:)   => null(), &
      c4(:,:,:,:) => null()
    integer, pointer:: &
      i0          => null(), &
      i1(:)       => null(), &
      i2(:,:)     => null()
    real, pointer:: &
      r0          => null(), &
      r1(:)       => null(), &
      r2(:,:)     => null(), &
      r3(:,:,:)   => null(), &
      r4(:,:,:,:) => null()
    real(8), pointer:: &
      d0          => null(), &
      d1(:)       => null(), &
      d2(:,:)     => null(), &
      d3(:,:,:)   => null(), &
      d4(:,:,:,:) => null()
    character(len=32):: name
    integer:: rank
  contains
    procedure:: output => output_item
  end type
  integer:: version=1, verbose=0, uniq=200, iout_prv=-1
  character(len=32), pointer:: select(:)
  type(aux_t), public:: aux
CONTAINS

!===============================================================================
!> This routine is called from data_io%init, already inside an OMP lock region
!===============================================================================
SUBROUTINE init (self, name)
  class(aux_t):: self
  character(len=*), optional:: name
  !.............................................................................
  integer:: iout
  real(8):: aux_time = 1.0
  integer:: iostat, i
  namelist /aux_params/ verbose, select
  logical, save:: first_time=.true.
  character(len=120):: id= &
  '$Id$'
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%init')
  call trace%print_id (id)
  call io%lock%set ('aux_t%init')
  if (first_time) then
    if (.not.associated(select)) then
      allocate (select(max_select))
      select(:) = ''
    end if
    rewind (stdin)
    read (stdin, aux_params, iostat=iostat)
    call io%namelist_warning ('aux_params', iostat)
    call write_namelist (io%namelist)
    first_time = .false.
  end if
  self%aux_next = io%start_time
  call self%dll_t%init
  call io%lock%unset ('aux_t%init')
  !-----------------------------------------------------------------------------
  ! Initialize the underlying data type, with its lock
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

!===============================================================================
FUNCTION new_item (self, name)
  class(aux_t)     :: self
  character(len=*) :: name
  logical          :: new_item
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (verbose > 0) call trace%begin ('aux_t%new_item')
  call self%lock%set ('new_item')
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
  if (verbose > 1) &
    write (io_unit%log,*) 'aux_t%new_item: ', trim(name), new_item
  call self%lock%unset ('new_item')
  if (verbose > 0) call trace%end ()
  END FUNCTION new_item

!===============================================================================
! Scalar output: Allow the argument to be non-pointer, and allocate a copy
!===============================================================================
SUBROUTINE register0i (self, name, i)
  class(aux_t)       :: self
  character(len=*)   :: name
  integer            :: i
  integer , pointer  :: ip
  class(*), pointer  :: c
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register0r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 0
      item%name = name
      allocate (ip)
      ip = i
      item%i0 => ip
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          deallocate (item%i0)
          allocate (ip)
          ip = i
          item%i0 => ip
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call trace%end()
END SUBROUTINE register0i

!===============================================================================
! Scalar output: Allow the argument to be non-pointer, and allocate a copy
!===============================================================================
SUBROUTINE register0r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real               :: r
  real    , pointer  :: rp
  class(*), pointer  :: c
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register0r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 0
      item%name = name
      allocate (rp)
      rp = r
      item%r0 => rp
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          deallocate (item%r0)
          allocate (rp)
          rp = r
          item%r0 => rp
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call trace%end()
  call self%list
END SUBROUTINE register0r

!===============================================================================
SUBROUTINE register1i (self, name, i)
  class(aux_t)       :: self
  character(len=*)   :: name
  integer , pointer  :: i(:)
  !class(*), pointer  :: c(:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register1i')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 1
      item%name = name
      item%i1 => i
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%i1 => i
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(1i) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register1i

!===============================================================================
SUBROUTINE register1r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:)
  class(*), pointer  :: c(:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register1r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 1
      item%name = name
      item%r1 => r
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%r1 => r
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(1r) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register1r

!===============================================================================
SUBROUTINE register2i (self, name, i)
  class(aux_t)       :: self
  character(len=*)   :: name
  integer , pointer  :: i(:,:)
  !class(*), pointer  :: c(:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register2i')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 2
      item%name = name
      item%i2 => i
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%i2 => i
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(2i) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register2i

!===============================================================================
SUBROUTINE register2r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register2r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 2
      item%name = name
      item%r2 => r
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%r2 => r
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(2r) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register2r

!===============================================================================
SUBROUTINE register2d (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real(8), pointer   :: r(:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register2r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 2
      item%name = name
      item%d2 => r
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%d2 => r
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(2r) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register2d

!===============================================================================
SUBROUTINE register3r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register3r')
  if (self%new_item (name)) then
    allocate (item_t:: item)
    select type (item)
    class is (item_t)
      item%rank = 3
      item%name = name
      item%r3 => r
    end select
    call self%append (item)
    if (verbose > 0) &
      write (io_unit%log,*) self%n, ' items registered'
  else
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        if (trim(item%name) == trim(name)) then
          item%r3 => r
          if (verbose > 1) &
            write (io_unit%log,*) 'aux_io_t: re-registered(3r) '//name
          exit
        end if
      end select
      item => item%next
    end do
  end if
  call self%list
  call trace%end()
END SUBROUTINE register3r

!===============================================================================
SUBROUTINE register4r (self, name, r)
  class(aux_t)       :: self
  character(len=*)   :: name
  real    , pointer  :: r(:,:,:,:)
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('aux_t%register4r')
  if (verbose > 1) &
    write (io_unit%log,*) 'aux_t%register4r: '//name
  if (associated(r)) then
    if (self%new_item (name)) then
      allocate (item_t:: item)
      select type (item)
      class is (item_t)
        item%rank = 4
        item%name = name
        item%r4 => r
        if (verbose > 1) &
          write (io_unit%log,*) 'aux_t%register4r: added '//name, &
            associated(r), associated(item%r4)
      end select
      call self%append (item)
    else
      item => self%head
      do while (associated(item))
        select type (item)
        class is (item_t)
          if (verbose > 2) &
            write (io_unit%log,*) 'aux_t%register4r: item =', trim(item%name)
          if (trim(item%name) == trim(name)) then
            item%r4 => r
            if (verbose > 1) &
              write (io_unit%log,*) 'aux_t%register4r: re-registered '//name
            exit
          end if
        end select
        item => item%next
      end do
    end if
    if (verbose > 0) &
      write (io_unit%log,*) 'aux_t%register4r:', self%n, ' items registered'
  else
    call io%warning ('aux_t%register4r called with pointer not associated')
  end if
  call self%list
  call trace%end()
END SUBROUTINE register4r

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
  call self%list
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
  call self%list
END SUBROUTINE register2

!===============================================================================
!===============================================================================
SUBROUTINE list (self)
  class(aux_t):: self
  class(dll_node_t), pointer:: item
  !-----------------------------------------------------------------------------
  if (verbose > 2) then
    item => self%head
    do while (associated(item))
      select type (item)
      class is (item_t)
        write (io_unit%log,*) 'aux_io_t: '//trim(item%name)
      end select
      item => item%next
    end do
  end if
END SUBROUTINE list

!===============================================================================
SUBROUTINE output_sequence (self, id, time, singlefile)
  class(aux_t):: self
  integer:: id
  character(len=*), optional:: singlefile
  real(8), optional:: time
  !-----------------------------------------------------------------------------
  call trace%begin ('aux_t%output_sequence')
  if (io%do_output) then
    call io%lock%set ('aux_output_seq')
    if (verbose > 1) then
      write (io_unit%log,'(a,i6,2f12.6)') &
        'aux_t%output_sequence: id, time, aux_next =', &
        id, time, self%aux_next
    end if
    !...........................................................................
    do while (time >= self%aux_next)
      call os%mkdir_next_output (self%iout+1)
      call output_aux (self, self%iout, id, singlefile)
      self%aux_next = self%aux_next + io%out_time
      self%iout = self%iout + 1
    end do
    !...........................................................................
    if (verbose > 2) then
      write (io_unit%log,'(a,i6,2f12.6)') &
        'aux_t%output_sequence: id, time, aux_next =', &
        id, time, self%aux_next
    end if
    call io%lock%unset ('aux_output_seq')
  end if
  call trace%end ()
END SUBROUTINE output_sequence

!===============================================================================
! This procedure is called from data_io_t%output, inside a locked region, and
! only when io%do_output is true, so it needs no lock and no test on do_output.
! When called from other places, it should also be from within locked regions.
!===============================================================================
SUBROUTINE output_aux (self, iout, id, singlefile)
  class(aux_t):: self
  integer:: iout, id
  character(len=*), optional:: singlefile
  !.............................................................................
  character(len=64):: filename
  class(dll_node_t), pointer:: item
  logical:: ok
  integer:: i, iostat
  namelist /aux_params/ verbose, select
  !-----------------------------------------------------------------------------
  ! Make sure only one output runs at a time, so there is no interference of
  ! reading stdin and deciding if to write messages only once
  !-----------------------------------------------------------------------------
  if (self%n == 0) &
    return
  call trace%begin ('aux_t%output')
  call io%lock%set ('output_aux')
  call self%list
  !-----------------------------------------------------------------------------
  ! Once per snapshot, re-read the input namelist, to allow edits
  !-----------------------------------------------------------------------------
  if (iout > iout_prv) then
    rewind (io%input)
    read (io%input, aux_params, iostat=iostat)
    iout_prv = iout
  end if
  !-----------------------------------------------------------------------------
  ! Open one file per task, or use a single file
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    write (io_unit%mpi,'(f12.6,i4,1x,a,i6)') wallclock(), omp%thread, &
    'aux_t%ouput_aux: opened unit', io_unit%aux
  if (present(singlefile)) then
    filename = singlefile
    open (io_unit%aux, file=trim(singlefile), form='unformatted', &
      access='sequential', status='unknown')
  else
    write (filename,'(a,i5.5"/",i5.5,".aux")') trim(io%outputname), iout, id
    open (io_unit%aux, file=trim(filename), form='unformatted', &
      access='sequential', status='unknown')
  end if
  if (verbose > 0) &
    write (io_unit%log,*) filename
  write (io_unit%aux) version, id
  !-----------------------------------------------------------------------------
  ! Run through the items, and write them out
  !-----------------------------------------------------------------------------
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
    if (verbose > 1 .or. iout > iout_prv) &
      write (io_unit%log,'(a,a,2x,a,l3)') &
        'aux_t%output: ', trim(filename), trim(item%name), ok
    if (ok) then
      call item%output (io_unit%aux) 
    end if
    end select
    item => item%next
  end do
  close (io_unit%aux)
  if (verbose > 0) &
    write (io_unit%mpi,'(f12.6,i4,1x,a,i6)') wallclock(), omp%thread, &
    'aux_t%ouput_aux: closed unit', io_unit%aux
  call io%lock%unset ('output_aux')
  call trace%end()
END SUBROUTINE output_aux

!===============================================================================
! Emulate namelist output; ifort prints out text items without quotation
! marks, which causes problems in the Python reader
!===============================================================================
SUBROUTINE write_namelist (unit)
  integer:: unit, i
  write (unit,'(a)') ' &AUX_PARAMS'
  write (unit,'(a,i10,",")') ' VERBOSE=', verbose
  write (unit,'(a,$)') ' SELECT='
  do i=1,max_select
    if (trim(select(i)) /= '') then
      write (unit,'(a,$)') "'"//trim(select(i))//"', "
    end if
  end do
  write (unit,'(/a)') '/'
  flush (unit)
END SUBROUTINE write_namelist

!===============================================================================
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
    if (associated (self%r1)) then
      call write1r (unit, self%r1)                             ! write: value
    else if (associated (self%i1)) then
      call write1i (unit, self%i1)                             ! write: value
    end if
  case(2)
    if (associated (self%r2)) then
      call write2r (unit, self%r2)                             ! write: value
    else if (associated (self%i2)) then
      call write2i (unit, self%i2)                             ! write: value
    else if (associated (self%d2)) then
      call write2d (unit, self%d2)                             ! write: value
    end if
  case(3)
    if (associated (self%r3)) then
      if (verbose > 1) then
        write (io_unit%log,*) trim(self%name), ' shape:', shape(self%r3)
        flush (io_unit%log)
      end if
      call write3 (unit, self%r3)                             ! write: value
    else
      write (stderr,*) 'self%r3 pointer not associated for ', trim(self%name)
      write (io_unit%log,*) 'self%r3 pointer not associated for ', trim(self%name)
    end if
  case(4)
    if (associated (self%r4)) then
      if (verbose > 1) then
        write (io_unit%log,*) trim(self%name), ' shape:', shape(self%r4)
        flush (io_unit%log)
      end if
      call write4 (unit, self%r4)
    else
      write (stderr,*) 'self%r4 pointer not associated for ', trim(self%name)
      write (io_unit%log,*) 'self%r4 pointer not associated for ', trim(self%name)
    end if
  end select
  flush (unit)
  call trace%end()
!===============================================================================
contains
END SUBROUTINE 

!===============================================================================
subroutine write1i (unit, r)
  integer:: unit
  integer, pointer:: r(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write1i')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'i1', shape(r)
  write (unit) 'i1'
  write (unit) r
  call trace%end()
end subroutine write1i

!===============================================================================
subroutine write1r (unit, r)
  integer:: unit
  real, pointer:: r(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write1r')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'r1', shape(r)
  write (unit) 'r1'
  write (unit) r
  call trace%end()
end subroutine write1r

!===============================================================================
subroutine write2r (unit, r)
  integer:: unit
  real, pointer:: r(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write2r')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'r2', shape(r)
  write (unit) 'r2'
  write (unit) r
  call trace%end()
end subroutine write2r

!===============================================================================
subroutine write2d (unit, d)
  integer:: unit
  real(8), pointer:: d(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write2d')
  write (unit) shape(d)
  if (verbose > 2) write (io_unit%log,*) 'd2', shape(d)
  write (unit) 'd2'
  write (unit) d
  call trace%end()
end subroutine write2d

!===============================================================================
subroutine write2i (unit, r)
  integer:: unit
  integer, pointer:: r(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write2i')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'i2', shape(r)
  write (unit) 'i2'
  write (unit) r
  call trace%end()
end subroutine write2i

!===============================================================================
subroutine write3 (unit, r)
  integer:: unit
  real, pointer:: r(:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('write3')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'r3', shape(r)
  write (unit) 'r3'
  write (unit) r
  call trace%end()
end subroutine write3

subroutine write4 (unit, r)
  integer:: unit
  real, pointer:: r(:,:,:,:)
  call trace%begin ('write4')
  write (unit) shape(r)
  if (verbose > 2) write (io_unit%log,*) 'r4', shape(r)
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
