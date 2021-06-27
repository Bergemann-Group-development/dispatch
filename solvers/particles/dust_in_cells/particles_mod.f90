!===============================================================================
!> Particles represented by arrays.  A patch task, which is handled by one
!> thread, has an array of these data types, one for each cell, where before
!> each update the cell address particles_t%q(1:4) is the same as the index
!> of the data type instance in the array of instances.  After an update, this
!> is no longer true.  The next thread that starts an update must therefore,
!> as a first step, allocate a new array of instances, knowing the number of
!> particles in each cell because it was counted during the previous update.
!> After copying each particle (= its attributes) from the previous array of
!> particles_t instances, that the instances may be deallocated and when they
!> all are the array of instances can be deallocated as well.  In fact, one
!> can delay allocating new cell instances until the 1st time they are needed,
!> and one can deallocate each previous one as soon as the particle attributes
!> have been copied over.
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  implicit none
  private
  integer, parameter:: n_words = 2*3 + 2 + 1 + 1
  type, public:: particles_t
    real      , pointer, contiguous :: r(:,:)  => null()    ! 12 bytes
    real      , pointer, contiguous :: v(:,:)  => null()    ! 12 bytes
    real      , pointer, contiguous :: w(:)    => null()    !  4 bytes
    integer(2), pointer, contiguous :: q(:,:)  => null()    !  8 bytes
    integer(2), pointer, contiguous :: flag(:) => null()    !  2 bytes
    !.......................................................! 38 bytes/particle
    type(lock_t), pointer           :: lock    => null()
    integer                         :: verbose = 0
    integer                         :: id(5) = 0
    integer                         :: n = 0          ! actual count
    integer                         :: np_new = 0     ! new particle count
    integer                         :: lp = 0         ! array length
    integer                         :: mp = 0         ! max count
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: pack
    procedure:: validate
    procedure:: resize
  end type
CONTAINS

!===============================================================================
!> Initialize particle arrays
!===============================================================================
SUBROUTINE init (self, mp, verbose)
  class(particles_t):: self
  integer:: mp, verbose
  !-----------------------------------------------------------------------------
  allocate (self%r(mp,3), self%v(mp,3), self%q(mp,3), self%w(mp), &
            self%flag(mp))
  call io%bits_mem (storage_size(self%w), mp*n_words, 'p')
  allocate (self%lock)
  call self%lock%init ('particles_t%init')
  self%n  = 0
  self%lp = 0
  self%mp = mp
  self%verbose = verbose
END SUBROUTINE init

!===============================================================================
!> Initialize particle arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(particles_t):: self
  !-----------------------------------------------------------------------------
  call io%bits_mem (storage_size(self%w), -self%mp**n_words, 'p')
  deallocate (self%r, self%v, self%q, self%w, self%flag)
  call self%lock%destroy
  deallocate (self%lock)
END SUBROUTINE dealloc

!===============================================================================
!> Resize
!===============================================================================
SUBROUTINE resize (self, mp)
  class(particles_t):: self
  integer:: mp, mo
  real, pointer:: r2(:,:), r1(:)
  integer(kind=2), pointer:: i2(:,:)
  !-----------------------------------------------------------------------------
  mo = self%mp
  call io%bits_mem (storage_size(self%w), -mo*n_words, 'p')
  call io%bits_mem (storage_size(self%w), +mp*n_words, 'p')
  r2 => self%r; allocate (self%r(mp,3)); self%r(1:mo,:) = r2(:,:)
  r2 => self%v; allocate (self%v(mp,3)); self%v(1:mo,:) = r2(:,:)
  i2 => self%q; allocate (self%q(mp,3)); self%q(1:mo,:) = i2(:,:)
  r1 => self%w; allocate (self%w(mp  )); self%w(1:mo  ) = r1(:)
  self%mp = mp
END SUBROUTINE resize

!===============================================================================
!> Remove holes in particle arrays
!===============================================================================
SUBROUTINE pack (self, n)
  class(particles_t):: self
  integer:: n
  integer:: j1, j2, i, lp1, np0, lp0
  !-----------------------------------------------------------------------------
  if (self%lp > self%n) then
    if (self%verbose > 0) then
      if (count(self%w(1:self%lp)==0.0) > self%lp-self%n) then
        write (stdlog,'(a,i5,4i4,3i5)') &
          'pack:', self%id, self%n, self%lp, count(self%w(1:self%lp) == 0.0)
        call io%abort ('too many zero-w in pack')
      end if
    end if
    !---------------------------------------------------------------------------
    ! Search for additional holes, if any
    !---------------------------------------------------------------------------
    j2 = self%lp
    lp1 = self%lp
    j1 = 1
    n = self%lp-self%n
    np0 = self%n
    lp0 = self%lp
    do while (self%lp > self%n)
      if (self%w(j1) == 0.0)   then
        !-----------------------------------------------------------------------
        ! Shrink from the end, looking for a valid entry beyond j1
        !-----------------------------------------------------------------------
        do while (self%w(j2) == 0.0 .and. j2 > 0)
          if (j2 < self%n) then
            write (stdlog,'(a,i5,4i4,3i5)') &
            'pack:', self%id, self%n, self%lp, count(self%w(1:self%lp) == 0.0)
            call io%abort ('too many zero-w in pack(2)')
          end if
          j2 = j2-1
        end do
        !-----------------------------------------------------------------------
        ! Unless j2 has shrunk below j1, copy the valid entry to the hole
        !-----------------------------------------------------------------------
        if (j2 > j1) then
          do i=1,3
            self%q(j1,i) = self%q(j2,i)
            self%r(j1,i) = self%r(j2,i)
            self%v(j1,i) = self%v(j2,i)
          end do
          self%w(j1) = self%w(j2)
          j2 = j2-1                                   ! previous end slot
        end if
        self%lp = j2
      end if
      j1 = j1+1
    end do
  end if
  if (self%n /= self%lp) then
    write (stdlog,'(a,i6,3(2x,2i4))') &
    'pack ERROR: id,j1,j2,np,lp=', self%id, j1, j2, self%n, self%lp, lp0, np0
    call io%abort ('particles_t%pack')
  end if
END SUBROUTINE pack

!===============================================================================
!> Remove holes in particle arrays
!===============================================================================
SUBROUTINE validate (self, label)
  class(particles_t):: self
  character(len=*):: label
  integer:: j, n, c
  !-----------------------------------------------------------------------------
  n = self%lp - self%n
  c = count(self%w(1:self%lp)==0.0)
  if (c /= n) then
    write (stderr,'(a,i4,i6,7i4)') &
      'validate ERROR: zero-w particle remains at '//label, &
      j, self%id, self%n, self%lp, c
    write (stdlog,'(a,i4,i6,7i4)') &
      'validate ERROR: zero-w particle remains at '//label, &
      j, self%id, self%n, self%lp, c
    call io%abort ('particles_t%validate')
  end if
END SUBROUTINE validate

END MODULE particles_mod
