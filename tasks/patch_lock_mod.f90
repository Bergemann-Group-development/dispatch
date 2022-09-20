!===============================================================================
!> Module that provides an atomic lock service for keys that consist of 
!> [rank,key], where key is the hash key of a proposed patch on the rank.
!>
!> The method is to get (and lock) access to an integer(n,5) array with slots
!> for proposed patches, where the 5 integers are [rank,hash_key].  If a
!> corresponding slot exists, such a patch is being created.  If not, the
!> values are added and the array is released.
!>
!> The slot is only cleared when the send buffer with the new patch has been
!> received.  This is done when the buffer is removed from the sent list.
!===============================================================================
MODULE patch_lock_mod
    USE io_mod
    USE trace_mod
    USE mpi_mod
    USE mpi_global_mod
    implicit none
    private
    !-----------------------------------------------------------------------------
    type, public:: patch_lock_t
      integer:: n                                   ! number of slots
      type(mpi_global_t), pointer:: global=>null()  ! global integer(n,5) array
    contains
      procedure:: init
      procedure:: set
      procedure:: unset
    end type
    type(patch_lock_t), public:: patch_lock
  CONTAINS
  
  !===============================================================================
  !> Initialize a global, atomic integer(n,5) array
  !===============================================================================
  SUBROUTINE init (self, n_per_rank)
    class(patch_lock_t):: self
    integer:: n_per_rank
    !.............................................................................
    integer, pointer, contiguous:: locks(:,:)
    integer:: n
    !-----------------------------------------------------------------------------
    call trace%begin ('patch_lock_t%init')
    allocate (self%global)
    self%n = n_per_rank
    n = mpi%size*self%n
    allocate (locks(n,5))
    locks(:,:) = -1
    call self%global%init (locks)
    deallocate (locks)
    call trace%end ()
  END SUBROUTINE init
  
  !===============================================================================
  !> Set a patch creation lock, by adding in an empty slot, the rank and hash key
  !> of the new patch
  !===============================================================================
  SUBROUTINE set (self, rank, key, ok)
    class(patch_lock_t):: self
    logical:: ok, exists
    integer:: rank, key(4), i
    integer:: lock(5)
    integer:: locks(self%n,5)
    type(mpi_global_t), pointer:: global
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    ! Acquire a copy of the global array, and keep the acess locked
    !-----------------------------------------------------------------------------
    call trace%begin ('patch_lock_t%set')
    global => self%global
    call global%get (locks)
    !-----------------------------------------------------------------------------
    ! Check if the key already exists
    !-----------------------------------------------------------------------------
    exists = .false.
    lock(1) = rank
    lock(2:5) = key
    do i=1,self%n
      if (all(locks(i,:)==lock)) then
         exists = .true.
         exit
      end if
    end do
    !-----------------------------------------------------------------------------
    ! If it does not, seach for an empty slot and add the proposed key
    !-----------------------------------------------------------------------------
    if (exists) then
      ok = .false.
    else
      ok = .false.
      do i=1,self%n
        if (locks(i,1) == -1) then
          locks(i,:) = lock
          ok = .true.
          exit
        end if
      end do
    end if
    !-----------------------------------------------------------------------------
    ! Put back (and release) the global array
    !-----------------------------------------------------------------------------
    call global%put (locks)
    call trace%end ()
  END SUBROUTINE set
  
  !===============================================================================
  !> Unset a patch creation lock
  !===============================================================================
  SUBROUTINE unset (self, rank, key, ok)
    class(patch_lock_t):: self
    integer:: rank, key(4), i
    logical:: ok
    !.............................................................................
    integer:: locks(self%n,5)
    integer:: lock(5)
    type(mpi_global_t), pointer:: global
    !-----------------------------------------------------------------------------
    ! Search for the lock and zap it, if found
    !-----------------------------------------------------------------------------
    call trace%begin ('patch_lock_t%unset')
    global => self%global
    call global%get (locks)
    lock(1) = rank
    lock(2:5) = key
    ok = .false.
    do i=1,self%n
      if (all(locks(i,:)==lock)) then
         locks(i,:) = -1
         ok = .true.
         exit
      end if
    end do
    global => self%global
    call global%put (locks)
    call trace%end ()
  END SUBROUTINE unset
  
  END MODULE patch_lock_mod