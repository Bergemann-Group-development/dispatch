!===============================================================================
!> Hash table module for the use inside DISPATCH
!> - KEY: A tuple (ip,id) acts as hash key
!> - VALUE: pointer
!> - HASH FUNCTION: Simple hash function based on multiplication with constants
!> - COLLISIONS: A linked list is used to deal with collisions
!> Author: Troels Haugboelle
!===============================================================================
MODULE hash_table_mod
  USE iso_c_binding, only: c_loc, c_ptr
  USE trace_mod
  USE io_mod
  USE timer_mod
  USE task_mod
  USE link_mod
  USE shared_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! General module parameters
  !-----------------------------------------------------------------------------
  integer(kind=8), dimension(1:4), parameter :: constants = (/5, -1640531527, 97, 1003313/)
  !-----------------------------------------------------------------------------
  ! Define a bucket as a derived type (sequence statement!) for better
  ! cache efficiency.
  !-----------------------------------------------------------------------------
  type, public:: bucket_t
    integer, allocatable, dimension(:) :: key
    integer :: next_ibucket
    class(*), pointer :: value
  end type
  !-----------------------------------------------------------------------------
  ! The actual hash table is an array of buckets
  !-----------------------------------------------------------------------------
  type, public:: hash_table_t
    type(bucket_t), allocatable, dimension(:)  :: data
    integer         :: id=0
    integer         :: ndim=4
    integer         :: key_length
    integer         :: total_size, head_free, nfree_chain, nfree
    integer(8)      :: n_add = 0_8
    integer(kind=8) :: size
    integer(kind=8) :: bitmask
    integer, allocatable, dimension(:) :: next_free
  contains
    procedure:: init
    procedure:: reset_entire_hash
    procedure:: dealloc
    procedure:: set_int4s
    procedure:: set_int4
    procedure:: set_int8
    procedure:: set_real4
    procedure:: set_task
    procedure:: set_link
    procedure:: set_ptr
    generic:: set => set_int4s, set_int4, set_int8, set_real4, set_task, set_link
    generic:: add => set_int4s, set_int4, set_int8, set_real4, set_task, set_link
    procedure:: get_int4
    procedure:: get_task
    procedure:: get_link
    procedure:: get_ptr
    generic:: get => get_int4, get_task, get_link
    procedure:: existsn
    procedure:: exists1
    generic:: exists => existsn, exists1
    procedure:: free
    procedure:: stats
    procedure:: same_keys
    procedure:: reset_bucket
    procedure, nopass:: test
    procedure:: remove => free
  end type
  type(hash_table_t), public:: hash_table
CONTAINS

!===============================================================================
!> Allocate all hash table arrays and variables.
!> Chose size (excluding the chaining space) as the smallest
!> power of two >= the required_size.
!===============================================================================
SUBROUTINE init (self, size, ndim)
  class(hash_table_t), intent(inout) :: self
  integer, optional  , intent(in)    :: size
  integer, optional  , intent(in)    :: ndim
  !.............................................................................
  integer                            :: req_size
  integer, save                      :: itimer=0, id=0
  !-----------------------------------------------------------------------------
  call trace%begin ('hash_table_t%init', itimer=itimer)
  req_size  = shared%option (size, 2**15-1)
  self%ndim = shared%option (ndim, 4)
  !$omp atomic update
  id = id + 1
  self%id = id
  self%size = 2
  self%key_length = self%ndim*8
  do while (self%size < req_size*4)
     self%size = self%size * 2
  end do
  call reset_entire_hash (self, .false.)
  call trace%end (itimer)
END SUBROUTINE init

!===============================================================================
!> Hash function
!===============================================================================
PURE FUNCTION hash_func(key)
  integer, dimension(:), intent(in) :: key
  integer(kind=8)                   :: hash_func
  hash_func = dot_product(key(:), constants(1:size(key)))
END FUNCTION hash_func

!===============================================================================
!> Subroutine to reset the entire hash table
!> IMPORTANT: The new size of the hash table is adapted based on the
!> load factor before resetting the hash table.
!===============================================================================
SUBROUTINE reset_entire_hash (htable, resize)
  class(hash_table_t), intent(inout) :: htable
  logical, intent(in)                :: resize
  integer                            :: i
  real                               :: load_factor
  integer, save                      :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('hash_table_t%reset', itimer=itimer)
  if (resize) then
     load_factor = real(htable%size - htable%nfree,kind=4) / real(htable%size,kind=4)
     if (load_factor > 0.6) then
        htable%size = htable%size * 2
        call io%bits_mem (-storage_size(htable%data), product(shape(htable%data)))
        call io%bits_mem (-storage_size(htable%next_free), &
          product(shape(htable%next_free)), 'hash')
        deallocate(htable%data, htable%next_free)
     else if (load_factor < 0.2 .and. htable%size > 2)then
        htable%size = htable%size / 2
        call io%bits_mem (-storage_size(htable%next_free), &
          product(shape(htable%next_free)), 'hash')
        deallocate(htable%data, htable%next_free)
     end if
  end if
  !-----------------------------------------------------------------------------
  ! Compute sizes and allocate arrays
  !-----------------------------------------------------------------------------
  htable%total_size = int(htable%size/4,kind=4) + int(htable%size,kind=4)
  htable%nfree = int(htable%size,kind=4)
  htable%nfree_chain = htable%total_size - int(htable%size,kind=4)
  htable%head_free = int(htable%size,kind=4) + 1
  htable%bitmask = htable%size - 1
  if (.not. allocated(htable%data))then
     allocate(htable%data(1: htable%total_size))
     allocate(htable%next_free (htable%size + 1: htable%total_size))
     call io%bits_mem (storage_size(htable%data), &
       product(shape(htable%data)), 'hash')
     call io%bits_mem (storage_size(htable%next_free), &
       product(shape(htable%next_free)), 'hash')
  end if
  !-----------------------------------------------------------------------------
  ! Initialize data
  !-----------------------------------------------------------------------------
  do i = 1, int(htable%total_size,kind=4)
     call htable%reset_bucket (htable%data(i))
  end do
  do i = int(htable%size,kind=4) + 1, htable%total_size - 1
     htable%next_free(i) = i + 1
  end do
  htable%next_free(htable%total_size) = 0
  call trace%end (itimer)
END SUBROUTINE reset_entire_hash

!===============================================================================
!> Deallocate the hash table
!===============================================================================
SUBROUTINE dealloc (htable)
  class(hash_table_t), intent(inout) :: htable
  integer :: i
  real :: load_factor
  !-----------------------------------------------------------------------------
  call io%bits_mem (-storage_size(htable%data), product(shape(htable%data)))
  call io%bits_mem (-storage_size(htable%next_free), &
    product(shape(htable%next_free)), 'hash')
  deallocate (htable%data, htable%next_free)
END SUBROUTINE dealloc

!===============================================================================
!> Reset the content of a bucket
!===============================================================================
SUBROUTINE reset_bucket (self, buck)
  class(hash_table_t):: self
  type(bucket_t), intent(inout) :: buck
  !-----------------------------------------------------------------------------
  if (allocated(buck%key)) then
    deallocate(buck%key)
  end if
  allocate (buck%key(self%ndim))
  buck%next_ibucket = -1
  buck%key = 0
END SUBROUTINE reset_bucket

SUBROUTINE set_int4s (self, key1, value)
  class(hash_table_t),   intent(inout) :: self
  integer,               intent(in)    :: key1
  integer, target,       intent(in)    :: value
  !.............................................................................
  class(*),              pointer       :: p
  integer, dimension(1)                :: key
  integer, save                        :: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) call trace%begin ('hash_table_t%int4s', itimer=itimer)
  key(1) = key1
  call set_int4 (self, key, value)
  if (timer%detailed) call trace%end (itimer)
END SUBROUTINE set_int4s

SUBROUTINE set_int4 (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  integer, target,       intent(in)    :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  allocate (integer:: p)
  select type (p)
  type is (integer)
  p = value
  end select
  call set_ptr (self, key, p)
END SUBROUTINE set_int4

SUBROUTINE set_int8 (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  integer(8), target,    intent(in)    :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  allocate (integer(8):: p)
  select type (p)
  type is (integer(8))
  p = value
  end select
  call set_ptr (self, key, p)
END SUBROUTINE set_int8

SUBROUTINE set_real4 (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  real(4),  target,      intent(in)    :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  allocate (real:: p)
  select type (p)
  type is (real)
  p = value
  end select
  call set_ptr (self, key, p)
END SUBROUTINE set_real4

SUBROUTINE set_task (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  class(task_t),         pointer       :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  if (.not. self%exists(key)) then
    p => value
    call set_ptr (self, key, p)
  end if
END SUBROUTINE set_task

SUBROUTINE set_link (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  class(link_t),         pointer       :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  if (.not. self%exists(key)) then
    p => value
    call set_ptr (self, key, p)
  end if
END SUBROUTINE set_link

!===============================================================================
!> Add a key/value pair to the hash table. If there is already a key/value
!> pair stored for this key, return an error message.
!===============================================================================
SUBROUTINE set_ptr (htable, key, val)
  class(hash_table_t),         intent(inout) :: htable
  integer, dimension(:)     ,  intent(in)    :: key
  class(*), pointer,           intent(in)    :: val
  integer(kind=8) :: full_hash
  integer(kind=8) :: ibucket
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute ibucket
  !-----------------------------------------------------------------------------
  !call trace%begin ('hash_table_t%set', itimer=itimer)
  associate (htable_n_add => htable%n_add)
  !$omp atomic update
  htable_n_add = htable_n_add + 1
  end associate
  full_hash = hash_func(key)
  ibucket = IAND(full_hash, htable%bitmask) + 1
  if (htable%data(ibucket)%next_ibucket < 0) then          
     ! Bucket is empty, simply insert value       
     htable%data(ibucket)%next_ibucket = 0
     htable%data(ibucket)%value => val
     htable%data(ibucket)%key = key
     htable%nfree = htable%nfree - 1
  !-----------------------------------------------------------------------------
  ! Bucket is not empty, walk through linked list
  !-----------------------------------------------------------------------------
  else if (htable%nfree_chain>0)then
     do while (htable%data(ibucket)%next_ibucket .ne. 0)
        ! Check if key already exists - abort if so
        if (htable%same_keys(htable%data(ibucket)%key,key)) then
          write(*,*) "hash_taable: trying to insert already existing key: ", key
          call io%abort ('hash_table: double insert')
        end if
        ibucket = htable%data(ibucket)%next_ibucket
     end do
     !--------------------------------------------------------------------------
     ! Check again (at the end of linked list)
     !--------------------------------------------------------------------------
     if (htable%same_keys(htable%data(ibucket)%key,key))then
        write(*,*) "trying to insert already existing key: ",key
        call io%abort ('hash_table: double insert')
     end if
     !--------------------------------------------------------------------------
     ! Have reached end of chain, val not present yet -> add
     !--------------------------------------------------------------------------
     htable%data(ibucket)%next_ibucket = htable%head_free
     ibucket = htable%head_free
     htable%data(ibucket)%next_ibucket = 0
     htable%data(ibucket)%value => val
     htable%data(ibucket)%key = key
     !--------------------------------------------------------------------------
     ! remove bucket from head of free linked list
     !--------------------------------------------------------------------------
     htable%head_free   = htable%next_free(htable%head_free)
     htable%nfree_chain = htable%nfree_chain - 1
  else
     write(*,*)  'hash chaining space full'
     write(*,*)  'size, n_add =', htable%size, htable%n_add
     stop
  end if
  !call trace%end (itimer)
END SUBROUTINE set_ptr

SUBROUTINE get_int4 (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  integer, target,       intent(out)   :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  call get_ptr (self, key, p)
  if (associated(p)) then
    select type (p)
    type is (integer)
    value = p
    end select
  else
    value = -999
  endif
END SUBROUTINE get_int4

SUBROUTINE get_task (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  class(task_t),         pointer       :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  call get_ptr (self, key, p)
  if (associated(p)) then
    select type (p)
    class is (task_t)
    value => p
    end select
  else
    nullify(value)
  endif
END SUBROUTINE get_task

SUBROUTINE get_link (self, key, value)
  class(hash_table_t),   intent(inout) :: self
  integer, dimension(:), intent(in)    :: key
  class(link_t),         pointer       :: value
  class(*),              pointer       :: p
  !-----------------------------------------------------------------------------
  call get_ptr (self, key, p)
  if (associated(p)) then
    select type (p)
    class is (link_t)
    value => p
    end select
  else
    nullify(value)
  endif
END SUBROUTINE get_link

!===============================================================================
!> Function (not subroutine, could also be changed...? ) which retrieves the 
!> hash table value for a given key. If no entry exists, return 0
!===============================================================================
SUBROUTINE get_ptr (htable, key, value)
  class(hash_table_t),            intent(in) :: htable
  integer    , dimension(:)     , intent(in) :: key
  class(*), pointer                          :: value
  integer(kind=8) :: ibucket, full_hash
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace%begin ('hash_table_t%get', itimer=itimer)
  full_hash = hash_func(key)
  ibucket = IAND(full_hash, htable%bitmask) + 1
  if (htable%same_keys(htable%data(ibucket)%key, key))then
     value => htable%data(ibucket)%value
     !call trace%end (itimer)
     return
  end if
  !-----------------------------------------------------------------------------
  ! Walk linked list until key is found or to the end is reached
  !-----------------------------------------------------------------------------
  do while( htable%data(ibucket)%next_ibucket > 0)
     ibucket = htable%data(ibucket)%next_ibucket
     if (htable%same_keys(htable%data(ibucket)%key, key))then
        value => htable%data(ibucket)%value
        !call trace%end (itimer)
        return
     end if
  end do
  !-----------------------------------------------------------------------------
  ! Nothing found...
  !-----------------------------------------------------------------------------
  nullify(value)
  !call trace%end (itimer)
END SUBROUTINE get_ptr

!===============================================================================
!> Determine if a key entry exists
!===============================================================================
LOGICAL FUNCTION existsn (self, key)
  class(hash_table_t)               :: self
  integer, dimension(:), intent(in) :: key
  class(*), pointer                 :: value
  !-----------------------------------------------------------------------------
  call self%get_ptr (key, value)
  existsn = associated(value)
END FUNCTION existsn

LOGICAL FUNCTION exists1 (self, key1)
  class(hash_table_t) :: self
  integer, intent(in) :: key1
  !.............................................................................
  class(*), pointer   :: value
  integer             :: key(1)
  !-----------------------------------------------------------------------------
  key(1) = key1
  call self%get_ptr (key, value)
  exists1 = associated(value)
END FUNCTION exists1

!===============================================================================
!> Remove the hash table entry for a given key, after checking
!===============================================================================
SUBROUTINE free (self, key)
  class(hash_table_t),         intent(inout) :: self
  integer , dimension(:)     , intent(in)    :: key
  !-----------------------------------------------------------------------------
  if (self%exists(key)) then
    call real_free (self, key)
  end if
END SUBROUTINE free

!===============================================================================
!> Remove the hash table entry for a given key 
!===============================================================================
SUBROUTINE real_free (htable, key)
  class(hash_table_t),         intent(inout) :: htable
  integer , dimension(:)     , intent(in)    :: key
  integer(kind=8) :: ibucket, previous_ibucket=0, full_hash
  full_hash = hash_func(key)
  ibucket = IAND(full_hash, htable%bitmask) + 1
  ! No collision case
  if (htable%data(ibucket)%next_ibucket == 0) then     
     htable%data(ibucket)%next_ibucket = -1
     htable%data(ibucket)%key = 0
     htable%nfree = htable%nfree + 1
  else
     ! Collision case
     do while (.not. htable%same_keys(htable%data(ibucket)%key, key))
        previous_ibucket=ibucket
        ibucket=htable%data(ibucket)%next_ibucket
     end do
     if (ibucket <= htable%size) then           
        ! It's the first element we need to erase: Move first element from chaning 
        ! space into bucket and do as if the value to remove had been in the chaning space
        htable%data(ibucket)%value => htable%data(htable%data(ibucket)%next_ibucket)%value
        htable%data(ibucket)%key = htable%data(htable%data(ibucket)%next_ibucket)%key
        previous_ibucket = ibucket
        ibucket = htable%data(ibucket)%next_ibucket
     end if
     ! fill the hole and reconnect linked list
     htable%data(previous_ibucket)%next_ibucket = htable%data(ibucket)%next_ibucket
     htable%next_free(ibucket) = htable%head_free
     htable%head_free = int(ibucket,kind=4)
     htable%nfree_chain = htable%nfree_chain + 1
  end if
END SUBROUTINE real_free

!===============================================================================
!> Check if keys are equivalent
!===============================================================================
PURE FUNCTION same_keys(self, key1, key2)
  class(hash_table_t)       , intent(in) :: self
  logical :: same_keys
  integer, dimension(:)     , intent(in) :: key1, key2
  logical, dimension(1:self%ndim)        :: ok
  integer                                :: i
  !-----------------------------------------------------------------------------
  do i = 1, self%ndim
     ok(i) = (key1(i)==key2(i))
  end do
  same_keys = ALL(ok)
END FUNCTION same_keys

!===============================================================================
!> Statistics
!===============================================================================
SUBROUTINE stats (htable)
  class(hash_table_t)::htable
  !-----------------------------------------------------------------------------
  if (io%master) then 
    write(*,*)"Total values stored in hash table: "&
         ,htable%total_size - htable%nfree - htable%nfree_chain
    write(*,*)"Size of hash table (without chaning space): "&
         ,htable%size
    write(*,*)"Load factor: "&
         ,(htable%size - htable%nfree) * 1.D0 / (htable%size + tiny(0.D0))
    write(*,*)"Total collisions in hash table: "&
         ,htable%total_size - htable%size - htable%nfree_chain
    write(*,*)"Collision fraction: "&
         ,real(htable%total_size - htable%size - htable%nfree_chain,kind=4)&
         /(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Perfect collision fraction (assuming perfect randomness): "&
         ,(htable%total_size - htable%nfree - htable%nfree_chain - &
         htable%size * (1.d0 - ((htable%size - 1.d0)/(htable%size)) &
         **(htable%total_size - htable%nfree - htable%nfree_chain))) & 
         *1./(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Fraction of collision space used: "&
         ,(htable%total_size - htable%size - htable%nfree_chain)&
         * 1.D0 / (htable%total_size - htable%size + tiny(0.D0))
  endif
END SUBROUTINE stats

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test ()
  type(hash_table_t), allocatable:: htable
  integer:: size, ndim, i
  integer, allocatable:: key(:)
  class(*), pointer:: value
  logical:: ok  
  !-----------------------------------------------------------------------------
  size = 30
  ndim = 4
  allocate (key(ndim))
  !-----------------------------------------------------------------------------
  ! Initialize and set 1-D values
  !-----------------------------------------------------------------------------
  allocate (htable)
  call htable%init (size, ndim)
  do i=1,size
    allocate (real:: value)
    select type (value)
    type is (real)
      value = i
    end select
    key = i*13
    call htable%set_ptr (key, value)
  end do
  !-----------------------------------------------------------------------------
  ! Retrieve 1-D values and deallocate
  !-----------------------------------------------------------------------------
  ok = .true.
  do i=1,size
    key = i*13
    call htable%get_ptr (key, value)
    select type (value)
    type is (real)
      if (nint(value) /= i) ok = .false.
    end select
    deallocate (value)
  end do
  call htable%dealloc ()
  if (ok) then
    write (stdout,*) 'hash_table%test passed'
  else
    write (stdout,*) 'hash_table%test failed'
  end if
END SUBROUTINE test

END MODULE hash_table_mod
