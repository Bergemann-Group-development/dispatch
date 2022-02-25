!===============================================================================
!> A patch will typically contain about 20-30 cells per side.  This means that
!> a typical ray packet for face diagonal rays (of which there are 6) have a
!> contain 20-30 rays, while the axis rays contain many more, and the space
!> diagonals (of which there are 4) contain fewer. 
!>
!> An important question is then: How much time does it take to find the ICs
!> for a ray packet, if these are indexed by some integer?  There will typically
!> be n_packet ~ 20-30 to search among.  If linear seach in a list is too slow
!> one can allocate a table, and use binary search in the table to find the index.
!===============================================================================
MODULE test_mod
  implicit none
  public
  type link_t
    type(link_t), pointer:: next => null()
    integer:: id
    real, dimension(:), pointer:: bndry => null()
  end type
CONTAINS

FUNCTION search (head, id) RESULT (link)
  type(link_t), pointer:: head, link
  integer:: id
  link => head
  do while (associated(link))
    if (link%id == id) return
    link => link%next
  end do
END FUNCTION
END MODULE

PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  USE test_mod
  USE random_mod
  implicit none
  type(link_t), pointer:: link, head, tail
  integer :: i, n, it, nit, id
  !-----------------------------------------------------------------------------
  call omp%init
  !-----------------------------------------------------------------------------
  ! Create a short linked list, with pointers to the bndry data arrays
  !-----------------------------------------------------------------------------
  allocate (tail)
  allocate (tail%bndry(n))
  head => tail
  head%id = 0
  n = 30
  do i=1,n
    allocate (link)
    allocate (link%bndry(n))
    link%id = i
    tail%next => link
    tail => link
  end do
  !-----------------------------------------------------------------------------
  ! Search repeatedly for a random id in the list
  !-----------------------------------------------------------------------------
  nit = 10000000
  call tic
  do it=1,nit
    id = n*random%ran1()
  end do
  call timer%toc ('  void', nit*n)
  call tic
  do it=1,nit
    id = n*random%ran1()
    link => search (head, id)
  end do
  call timer%toc ('search', nit*n)

END PROGRAM

