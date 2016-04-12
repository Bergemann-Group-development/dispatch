!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> A patch extension with pebbles.  There are so many patch properties that the
!> pebbles need, that it is much simpler to use an extension of an mhd_t, rather
!> than to pass (or/and recompute) those properties.  The key flag is do_pebbles,
!> which modifies slightly the variable allocations in mhd_mod.f90, so that
!> gas pressure and gas velocities are permanently allocated, and are available
!> to the pebble procedures.  To simplify neighbor access, the patch also has
!> an anonymous pointer to its primary list node, which has the nbor head of the
!> neighbor list.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebbles_mod
  USE io_mod
  USE omp_lock_mod
  USE trace_mod
  USE bits_mod
  USE pebble_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! A linked list of pebbles
  !-----------------------------------------------------------------------------
  type, public:: pebbles_t
    integer:: id = -1
    integer:: n = 0                             ! number of active pebbles
    integer:: n_old = 0                         ! have stepped to many times
    integer:: n_accrete = 0                     ! have accreted
    integer:: n_outside(3) = 0                  ! have left the patch
    integer:: n_behind = 0                      ! are behind self%time
    integer:: n_export = 0                      ! number of exported pebbles
    integer:: n_remove = 0                      ! number of removed pebbles
    integer:: n_frozen = 0                      ! number of frozen pebbles
    integer:: it_max                            ! max life time
    integer:: id_dbg = 0                        ! follow this pebble
    integer(8):: n_update = 0                   ! speed counter
    real:: ppc                                  ! particles per cell
    real:: mass                                 ! representative pebble mass
    real:: power                                ! random distribution
    real:: size                                 ! representative pebble size
    real:: density                              ! mass density
    real:: drag_coeff                           ! drag coefficient
    real:: t_min                                ! oldest t
    real:: dt_min, dt_max                       ! dt extrema
    real:: z_min                                ! smallest positive z
    real:: gamma                                ! gamma for sound speed
    real:: courant=0.1                          ! courant number for motion
    real:: cell_probability
    real:: position(3)                          ! relative to star
    class(pebble_t), pointer:: head => null()    ! first & last pebble
    class(pebble_t), pointer:: tail => null()    ! first & last pebble
    logical:: initialized=.false.               ! initialize flag
    type(lock_t):: lock
  contains
    procedure:: init                            ! initialize pebbles
    procedure:: append                          ! append a pebble
    procedure:: remove                          ! remove a pebble
    procedure:: freeze                          ! freeze a pebble
    procedure:: export                          ! export to other patches (or remove)
    procedure:: delete                          ! delete pebbles
    procedure:: check_list                      ! check list consistency
    procedure:: print_list_n                    ! print top of list
  end type
  integer, save:: id=0
  integer, save:: id_dbg=0
  logical, save:: debug=.false.
CONTAINS

!> =============================================================================
!> Append a pebble to the list
!> =============================================================================
SUBROUTINE init (self, id)
  class(pebbles_t):: self
  integer:: id
  !.............................................................................
  call self%lock%init ('pebb')
  if (debug) print*,'dbg: lock%id', id, self%lock%id
  self%id = id
  if (self%initialized) then
    print *,  self%id, 'ERROR: second init attempted'
    stop
  else
    if (io%verbose>0) &
      print*,'pebble_patch::init', self%id
  end if
  self%initialized = .true.
END SUBROUTINE init

!> =============================================================================
!> Append a pebble to the list
!> =============================================================================
SUBROUTINE append (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble
  !.............................................................................
  if (debug) print *,self%id,'dbg: append pebble',pebble%id,self%lock%id,self%n
  call self%lock%set ('append')
  !$omp flush (pebble,self)
  !-----------------------------------------------------------------------------
  ! Place the pebble at the tail, which might be empty, in which case the prev
  ! pointer becomes null, as it should be.  If not, it points to the currently
  ! last item
  !-----------------------------------------------------------------------------
  pebble%prev => self%tail  
  !-----------------------------------------------------------------------------
  ! If the head is null, then this is the first pebble and head should point to it.
  !-----------------------------------------------------------------------------
  if (.not.associated(self%head)) self%head => pebble
  !-----------------------------------------------------------------------------
  ! If the tail exists, the new item should be the next item after it
  !-----------------------------------------------------------------------------
  if (associated(self%tail)) self%tail%next => pebble 
  nullify (pebble%next)
  self%tail => pebble
  self%n = self%n + 1
  !$omp flush (pebble,self)
  call self%lock%unset ('append')
END SUBROUTINE append

!> =============================================================================
!> Remove a pebble from the list
!> =============================================================================
SUBROUTINE remove (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble, p
  !.............................................................................
  if (debug) print *,self%id,'dbg: remove pebble',pebble%id,self%lock%id,self%n
  call self%lock%set ('remove')
  !$omp flush (self)
  !-----------------------------------------------------------------------------
  ! If the pebble is the head item, its %prev should be null, and head should
  ! be pointed to the next item, which could be null(), if removing last item.
  !-----------------------------------------------------------------------------
  if (associated(pebble,self%head)) then
    if (associated(pebble,self%tail)) then
      nullify (self%head, self%tail, pebble%prev, pebble%next)
      if (io%verbose>1) print*,'removing last item'
      self%n = self%n - 1
      if (self%n /= 0) then
        print*,'ERROR: pebble list has non-zero size', self%n
        self%n = 0
      end if
      call self%lock%unset ('remove')
      return
    end if
    if (io%verbose>1) print*,'removing head'
    self%head => pebble%next
    if (associated(self%head)) nullify(self%head%prev)
  end if
  if (associated(pebble%prev)) then
    pebble%prev%next => pebble%next
  end if
!  print *,'head =', associated(self%head)
!  if (associated(self%head)) then
!    print *,'head%id =', self%head%id
!    print *,'head%next =', associated(self%head%next)
!    if (associated(self%head%next)) &
!      print *,'head%next%id =', self%head%next%id
!  end if
  !-----------------------------------------------------------------------------
  ! If the pebble is the tail item, its next should be null, and tail should
  ! be pointed to the previous item, which could be null(), if removing last item.
  !-----------------------------------------------------------------------------
  if (associated(pebble,self%tail)) then
    if (io%verbose>1) print*,'removing tail'
    self%tail => pebble%prev
    if (associated(self%tail)) nullify(self%tail%next)
  end if
  if (associated(pebble%next)) then
    pebble%next%prev => pebble%prev
  end if
  nullify (pebble%next,pebble%prev)
!  print *,'tail =', associated(self%tail)
!  if (associated(self%tail)) print *,'tail%id =', self%tail%id
  self%n = self%n - 1
  !$omp flush (self,pebble)
!  call self%print_list_n (3, 'after remove')
  call self%check_list ('after remove')
  call self%lock%unset ('remove')
END SUBROUTINE remove

!> =============================================================================
!> Freeze a pebble in place, to maintain indices in output files
!> =============================================================================
SUBROUTINE freeze (self, pebble)
  class(pebbles_t):: self
  class(pebble_t), pointer:: pebble, p
  !.............................................................................
  p%status = ior(p%status, bits%frozen)
END SUBROUTINE freeze

!> =============================================================================
SUBROUTINE export (self, pebble, nbor_list)
  class(pebbles_t), target:: self
  class(pebble_t), pointer:: pebble
  class(pebbles_t), pointer:: nbor_list
  integer, save:: nprint=20
  !.............................................................................
  if (debug) print*,'dbg: export', pebble%id, self%id, nbor_list%id, nbor_list%initialized
  call self%remove (pebble)
  if (associated(nbor_list)) then
    if (associated(nbor_list,self)) then
      print*,'ERROR: attempting to export to self', self%id, pebble%id
      stop
    end if
    call nbor_list%append (pebble)
    if (debug) print*,'dbg: export done', pebble%id, self%id, nbor_list%id, nbor_list%initialized
  else
    if (nprint > 0) then
      nprint = nprint-1
      print *,'exporting pebble to a black hole ...'
    end if
  end if
  call self%check_list ('self after export')
  call nbor_list%check_list ('nbor after export')
END SUBROUTINE export

!> =============================================================================
!> =============================================================================
SUBROUTINE delete (self)
  class(pebbles_t):: self
  !.............................................................................
  class(pebble_t), pointer:: pebble, old
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_t%delete', 2)
  pebble => self%head
  do while (associated(pebble))
    old => pebble
    pebble => pebble%next
    deallocate (old)
  end do
  self%n = 0
  call trace_end
END SUBROUTINE delete

!> =============================================================================
!> Check list consistency
!> =============================================================================
SUBROUTINE check_list (self, label)
  class(pebbles_t):: self
  character(len=*), optional:: label
  class(pebble_t), pointer:: p
  !.............................................................................
  integer:: n, id
  !-----------------------------------------------------------------------------
  if (io%verbose<=0) return
  if (self%n == 0 .and. associated(self%head)) then
    print *,'ERROR: check_list found empty list with head'
    return
  end if
  n=0
  id = -1
  p => self%head
  do while (associated(p))
    if (p%id == id) then
      print*,'ERROR: pebble repeated', p%id
    end if
    id = p%id
    p => p%next
    n=n+1
  end do
  if (n /= self%n) then
    if (present(label)) then
      print*, self%id,  'ERROR: inconsistent pebble list ', label, n, self%n
    else
      print*, self%id,  'ERROR: inconsistent pebble list ', n, self%n
    end if
    call self%print_list_n (max(n,self%n))
    self%n = n
  else if (io%verbose>3) then
    if (present(label)) then
      print*, 'check_list: '//label, self%n
    else
      print*, 'check_list:', self%n
    end if
  end if
END SUBROUTINE check_list

!> =============================================================================
!> Print top of list
!> =============================================================================
SUBROUTINE print_list_n (self, n, label)
  class(pebbles_t):: self
  integer:: n
  character(len=*), optional:: label
  !.............................................................................
  class(pebble_t), pointer:: p
  integer:: m
  !-----------------------------------------------------------------------------
  m = 0
  p => self%head
  if (present(label)) print*,'print_list_n: '//label
  do while (associated(p).and.m<n)
    if (associated(p%prev)) then
      if (associated(p%next)) then
        print*, m, p%prev%id, p%id, p%next%id, p%t
      else
        print*, m, p%prev%id, p%id
      end if
    else
      if (associated(p%next)) then
        print*, m, -1, p%id, p%next%id
      else
        print*, m, -1, p%id
      end if
    end if
    p => p%next
    m = m+1
  end do
END SUBROUTINE print_list_n

END MODULE pebbles_mod
