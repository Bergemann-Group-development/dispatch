!> =============================================================================
!> =============================================================================
PROGRAM test_mpi_counter
  USE io_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  implicit none
  !.............................................................................
  type(mpi_counter_t):: ticket_counter, access_counter
  integer:: ticket, access, i, try
  !-----------------------------------------------------------------------------
  call omp%init ()
  call mpi%init ()
  call io%init ()
  !-----------------------------------------------------------------------------
  call ticket_counter%init (0)
  call access_counter%init (0)
  if (mpi%rank == 0) then
    access = access_counter%update (1)
  end if
  ticket = ticket_counter%update (1)
  access = access_counter%update (0)
  call mpi%barrier ()
  !-----------------------------------------------------------------------------
  ! At this point all ranks have access==1, while one has ticket=1 and the other
  ! ones has larger tickets
  !-----------------------------------------------------------------------------
  do while (access < 31)
    !---------------------------------------------------------------------------
    ! The rank with ticket==1 gets to do I/O, and increases its ticket to 3,
    ! while also increasing the access to 2, which gives access to the 2nd rank.
    ! The rank with access tries 10 times before actually doing the (fake) I/O.
    !---------------------------------------------------------------------------
    try = try+1
    if (access == ticket) then
      if (try > 10) then
        PRINT *, wallclock(), mpi%rank, access, ticket, 'I/O'
        ticket = ticket_counter%update (1)
        access = access_counter%update (1)
        try = 0
      end if
    end if
    access = access_counter%update (0)
  end do
  call mpi%barrier ()
  call mpi%end ()
END PROGRAM
