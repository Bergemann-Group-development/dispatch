!===============================================================================
!> Particles represented by arrays.  A patch task, which is handled by one
!> thread, has an array of these data types, one for each cell, where before
!> each update the cell address particles_t%q(1:4) is the same as the index
!> of the data type instance in the array of instances.  The attribute %m is
!> used to indicate that a particle has moved to a new cell; this speeds up
!> testing a lot in the %trade() procedure.
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  USE dll_mod
  implicit none
  private
  integer, parameter:: n_words = 2*3 + 2 + 1 + 1
  type, public:: particles_t
    real      , pointer, contiguous :: r(:,:)  => null()    ! 12 bytes
    real      , pointer, contiguous :: v(:,:)  => null()    ! 12 bytes
    real      , pointer, contiguous :: w(:)    => null()    !  4 bytes
    integer(2), pointer, contiguous :: q(:,:)  => null()    !  6 bytes
    integer(2), pointer, contiguous :: m(:)    => null()    !  2 bytes
    real(8)   , pointer, contiguous :: t(:)    => null()    !  8 bytes
    !.......................................................! 44 bytes/particle
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
  integer, public:: resize_calls = 0
  !$omp threadprivate (resize_calls)
  !-----------------------------------------------------------------------------
  ! Since a patch may need to export particles multiple times before the nbors
  ! have picked up all particles from previous exports, we use a list of
  ! exported particles to handle this.  Once a list of particles have been
  ! added to the list of exports, it will not be touched / modified by the
  ! particle solver owning it until it is ready to be deleted = when it contains
  ! no more particles.  At that time no neighbor task has any reason to access
  ! it either, except for reading the %n attribute, which contains zero.
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: export_list_t
  contains
    procedure:: prune
  end type
  !-----------------------------------------------------------------------------
  ! To be able to add a list of export particle lists to the solver, we let it
  ! be an extension of the basic doubly-linked list (dll_t) data type.  Utility
  ! procedures are usd to maintain the list of exports, and the particle
  ! lists themselves.
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: export_item_t
    integer:: n
    class(particles_t), pointer:: exports(:) => null()
  contains
    procedure:: init => init_exports
    procedure:: dealloc => dealloc_exports
    procedure:: is_empty
  end type
CONTAINS

!===============================================================================
!> Initialize particle arrays
!===============================================================================
SUBROUTINE init (self, mp, verbose)
  class(particles_t):: self
  integer:: mp, verbose
  !-----------------------------------------------------------------------------
  allocate (self%r(mp,3), self%v(mp,3), self%q(mp,3), self%w(mp), self%t(mp), &
            self%m(mp))
  call io%bits_mem (storage_size(self%w), mp*n_words, 'p')
  allocate (self%lock)
  call self%lock%init ('particles_t%init')
  self%n  = 0
  self%lp = 0
  self%mp = mp
  self%verbose = verbose
END SUBROUTINE init

!===============================================================================
!> Deallocate particle arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(particles_t):: self
  !-----------------------------------------------------------------------------
  call io%bits_mem (storage_size(self%w), -self%mp*n_words, 'p')
  deallocate (self%r, self%v, self%q, self%w, self%t, self%m)
  call self%lock%destroy
  deallocate (self%lock)
END SUBROUTINE dealloc

!===============================================================================
!> Resize particle arrays
!===============================================================================
SUBROUTINE resize (self, mp, timeit)
  class(particles_t):: self
  integer, optional:: mp
  logical, optional:: timeit
  integer:: mo, mn
  real           , pointer, contiguous:: r1(:), r2(:,:)
  real(8)        , pointer, contiguous:: d1(:)
  integer(kind=2), pointer, contiguous:: i1(:), i2(:,:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (present(timeit)) call trace%begin ('particles_t%resize', itimer=itimer)
  resize_calls = resize_calls+1
  mo = self%mp
  if (present(mp)) then
    mn = mp
  else
    mn = 1.2*mo
  end if
  call io%bits_mem (storage_size(self%w), -mo*n_words, 'p')
  call io%bits_mem (storage_size(self%w), +mn*n_words, 'p')
  allocate (r2(mn,3)); r2(1:mo,:) = self%r; deallocate(self%r); self%r => r2
  allocate (r2(mn,3)); r2(1:mo,:) = self%v; deallocate(self%v); self%v => r2
  allocate (i2(mn,3)); i2(1:mo,:) = self%q; deallocate(self%q); self%q => i2
  allocate (r1(mn  )); r1(1:mo)   = self%w; deallocate(self%w); self%w => r1
  allocate (i1(mn  )); i1(1:mo)   = self%m; deallocate(self%m); self%m => i1
  allocate (d1(mn  )); d1(1:mo)   = self%t; deallocate(self%t); self%t => d1
  self%mp = mn
  if (present(timeit)) call trace%end (itimer)
END SUBROUTINE resize

!===============================================================================
!> Remove holes in particle arrays
!===============================================================================
SUBROUTINE pack (self, n, timeit)
  class(particles_t):: self
  integer:: n
  logical, optional:: timeit
  integer:: j1, j2, i, lp1, np0, lp0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (present(timeit)) call trace%begin ('particles_t%pack', itimer=itimer)
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
          self%t(j1) = self%t(j2)
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
  if (present(timeit)) call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Validate particle arrays
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

!===============================================================================
!> Prune empty members from a list of exports
!===============================================================================
SUBROUTINE prune (self)
  class(export_list_t):: self
  class(export_item_t), pointer:: item
  class(dll_node_t),pointer:: this, next
  type(particles_t),pointer:: pe, exports(:)
  integer:: i, ib, n_bin
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%prune', itimer=itimer)
  i = 0
  this => self%head
  do while (associated(this))
    next => this%next
    i = i+1
    select type (this)
    class is (export_item_t)
      item => this
    end select
    if (item%is_empty()) then
      write (stdlog,*) 'particles_t%prune: removing list', i
      call item%dealloc ()
      call self%remove (this)
      deallocate (this)
    end if
    this => next
  end do
  call trace%end (itimer)
END SUBROUTINE prune

!===============================================================================
!> Deallocate an item
!===============================================================================
SUBROUTINE init_exports (self, n)
  class(export_item_t) :: self
  integer       :: n
  !-----------------------------------------------------------------------------
  self%n = n
  allocate (self%exports(n))
END SUBROUTINE init_exports

!===============================================================================
!> Check if an item is empty6
!===============================================================================
FUNCTION is_empty (self) RESULT (yes)
  class(export_item_t):: self
  logical:: yes
  integer:: ib
  !-----------------------------------------------------------------------------
  yes = .true.
  if (associated(self%exports)) then
    do ib=1,self%n
      if (self%exports(ib)%n > 0) then
        yes = .false.
        exit
      end if
    end do
  end if
END FUNCTION is_empty

!===============================================================================
!> Deallocate an item
!===============================================================================
SUBROUTINE dealloc_exports (self)
  class(export_item_t):: self
  integer:: ib
  !-----------------------------------------------------------------------------
  if (associated(self%exports)) then
    do ib=1,self%n
      call self%exports(ib)%dealloc
    end do
    deallocate (self%exports)
  end if
END SUBROUTINE dealloc_exports

END MODULE particles_mod
