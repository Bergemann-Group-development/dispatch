!===============================================================================
!> Particles represented by arrays
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  implicit none
  private
  type, public:: particles_t
    real(8), pointer                    :: t(:)
    real           , pointer, contiguous:: r(:,:), v(:,:), w(:)
    integer(kind=2), pointer, contiguous:: q(:,:), flag(:)
    type(lock_t), pointer               :: lock => null()
    integer                             :: verbose = 0
    integer                             :: id(5)
    integer                             :: np = 0         ! actual count
    integer                             :: lp = 0         ! array length
    integer                             :: mp = 0         ! max count
  contains
    procedure:: init
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
  allocate (self%t(mp), self%r(mp,3), self%v(mp,3), self%q(mp,3), self%w(mp), &
            self%flag(mp))
  call io%bits_mem (storage_size(self%w), mp*(2+2+3*3), 'p')
  self%mp = mp
  self%np = 0
  self%lp = 0
  self%verbose = verbose
END SUBROUTINE init

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
  call io%bits_mem (storage_size(self%w), -mo*(2+2+3*3), 'p')
  call io%bits_mem (storage_size(self%w), +mp*(2+2+3*3), 'p')
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
  integer:: j1, j2, i, lp1
  !-----------------------------------------------------------------------------
  if (self%lp > self%np) then
    if (self%verbose > 0) then
      if (count(self%w(1:self%lp)==0.0) > self%lp-self%np) then
        write (io_unit%log,'(a,i5,4i4,3i5)') &
          'pack:', self%id, self%np, self%lp, count(self%w(1:self%lp) == 0.0)
        call io%abort ('too many zero-w in pack')
      end if
    end if
    !---------------------------------------------------------------------------
    ! Search for additional holes, if any
    !---------------------------------------------------------------------------
    j2 = self%lp
    lp1 = self%lp
    j1 = 1
    n = self%lp-self%np 
    do while (self%lp > self%np)
      if (self%w(j1) == 0.0)   then
        !-----------------------------------------------------------------------
        ! Shrink from the end, looking for a valid entry beyond j1
        !-----------------------------------------------------------------------
        do while (self%w(j2) == 0.0 .and. j2 > 0)
          if (j2 < self%np) then
            write (io_unit%log,'(a,i5,4i4,3i5)') &
            'pack:', self%id, self%np, self%lp, count(self%w(1:self%lp) == 0.0)
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
          self%t(j1) = self%t(j2)
          j2 = j2-1                                   ! previous end slot
        end if
        self%lp = j2
      end if
      j1 = j1+1
    end do
  end if
  if (self%np /= self%lp) then
    write (io_unit%log,'(a,i6,8i4)') &
    'pack ERROR: id,j1,j2,np,lp=', self%id, j1, j2, self%np, self%lp
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
  n = self%lp - self%np
  c = count(self%w(1:self%lp)==0.0)
  if (c /= n) then
    write (stderr,'(a,i4,i6,7i4)') &
      'validate ERROR: zero-w particle remains at '//label, &
      j, self%id, self%np, self%lp, c
    write (io_unit%log,'(a,i4,i6,7i4)') &
      'validate ERROR: zero-w particle remains at '//label, &
      j, self%id, self%np, self%lp, c
    call io%abort ('particles_t%validate')
  end if
END SUBROUTINE validate
END MODULE particles_mod
