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
  USE mpi_mod
  USE trace_mod
  USE io_mod
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
    integer         :: ndim=2
    integer         :: key_length
    integer         :: total_size, head_free, nfree_chain, nfree
    integer(kind=8) :: size
    integer(kind=8) :: bitmask
    integer, allocatable, dimension(:) :: next_free
  contains
    procedure:: init
    procedure:: reset_entire_hash
    procedure:: dealloc
    procedure:: set
    procedure:: get
    procedure:: free
    procedure:: stats
    procedure:: same_keys
    procedure:: reset_bucket
  end type
  type(hash_table_t), public:: hash_table
CONTAINS

!===============================================================================
!> Allocate all hash table arrays and variables.
!> Chose size (excluding the chaining space) as the smallest
!> power of two >= the required_size.
!===============================================================================
SUBROUTINE init (self, req_size, ndim)
  class(hash_table_t), intent(inout) :: self
  integer            , intent(in)    :: req_size
  integer, optional  , intent(in)    :: ndim
  !-----------------------------------------------------------------------------
  self%size = 2
  if (present(ndim)) then
    self%ndim = ndim
  else
    self%ndim = 2
  end if
  self%key_length = self%ndim*8
  do while (self%size < req_size*2)
     self%size = self%size * 2
  end do
  call reset_entire_hash (self, .false.)
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
  integer :: i
  real :: load_factor
  !-----------------------------------------------------------------------------
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
END SUBROUTINE reset_entire_hash

!===============================================================================
!> Subroutine to reset the entire hash table
!> IMPORTANT: The new size of the hash table is adapted based on the
!> load factor before resetting the hash table.
!===============================================================================
SUBROUTINE dealloc (htable, resize)
  class(hash_table_t), intent(inout) :: htable
  logical, intent(in)                :: resize
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

!===============================================================================
!> Add a key/value pair to the hash table. If there is already a key/value
!> pair stored for this key, return an error message.
!===============================================================================
SUBROUTINE set (htable, key, val)
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
          call mpi%abort ('hash_table: double insert')
        end if
        ibucket = htable%data(ibucket)%next_ibucket
     end do
     !--------------------------------------------------------------------------
     ! Check again (at the end of linked list)
     !--------------------------------------------------------------------------
     if (htable%same_keys(htable%data(ibucket)%key,key))then
        write(*,*) "trying to insert already existing key: ",key
        stop
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
     write(*,*)"hash chaining space full "
     stop
  end if
  !call trace%end (itimer)
END SUBROUTINE set

!===============================================================================
!> Function (not subroutine, could also be changed...? ) which retrieves the 
!> hash table value for a given key. If no entry exists, return 0
!===============================================================================
SUBROUTINE get (htable, key, value)
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
END SUBROUTINE get

!===============================================================================
!> Remove the hash table entry for a given key 
!===============================================================================
SUBROUTINE free (htable, key)
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
END SUBROUTINE free

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
  if (mpi%master) then 
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

END MODULE hash_table_mod
