!> =============================================================================
!> Module with simple forward linked list
!> =============================================================================
MODULE test_mod
  USE omp_mod
  implicit none
  private
  type, public:: link_t
    integer:: i=0
    type (link_t), pointer:: next => null()
  end type
CONTAINS
END MODULE

!> =============================================================================
!> Test demonstrating how to loop in parallel over the list
!> =============================================================================
PROGRAM test
  USE test_mod
  USE omp_mod
  implicit none
  integer:: i
  type(link_t), pointer:: head, link, prev, mylink
  !
  print *,'allocate loop'
  allocate (link)
  link%i = 1
  head => link
  prev => link
  do i=2,50
    allocate (link)
    link%i = i
    prev%next => link
    prev => link
    print*,i
  end do
  !
  print *,'deallocate loop'
  call omp%init
  link => head
  !$omp parallel private(mylink)
  do while (associated(link))
    !$omp critical
    mylink => link
    link => link%next
    !$omp end critical
    print *, omp%thread, mylink%i
    deallocate (mylink)
  end do
  !$omp end parallel
END PROGRAM
