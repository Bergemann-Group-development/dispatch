!> =============================================================================
!> Keep track of neighbor ranks and their loads, by sending and receiving short
!> messages, storing the info in a linked list
!> =============================================================================
MODULE rank_info_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  implicit none
  private
  type, public:: rank_info_t
    logical:: ok
    integer:: rank
    integer:: nq, n_swap
    real:: cadence, patch_cost
    real(8):: cost=0.0
    real(8):: time, wall, dtime
    real(8):: otime, owall
    type(rank_info_t), pointer:: next => null()
  contains
    procedure:: imbalance
    procedure:: measure_load
  end type
  type, public:: io_rank_info_t
    sequence
    logical:: ok
    integer:: rank
    integer:: nq, n_swap
    real(8):: cost=0.0
    real(8):: time, wall, dtime
    real(8):: otime, owall
  end type
  type, public:: rank_info_list_t
    type(rank_info_t), pointer:: head => null()
    type(rank_info_t), pointer:: tail => null()
    integer:: n=0
  contains
    procedure:: init
    procedure:: pack
    procedure:: unpack
    procedure:: add
    procedure:: append
    procedure:: remove
    procedure:: send
    procedure:: find
    procedure:: print
  end type
  type(rank_info_t), public:: rank_info
  type(rank_info_list_t), public:: rank_info_list
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(rank_info_list_t):: self
  !.............................................................................
  if (io%verbose<1) return
  !call self%test
END SUBROUTINE init

!> =============================================================================
!> Pack the rank_info into an mpi_mesg_t data type
!> =============================================================================
SUBROUTINE pack (self, info, mesg)
  class(rank_info_list_t):: self
  type(rank_info_t):: info
  class(mesg_t), pointer:: mesg
  type(io_rank_info_t):: io_rank_info
  integer:: n
  !.............................................................................
  call trace_begin ('rank_info_mod::pack')
  n = sizeof(io_rank_info)/4
  allocate (mesg)
  allocate (mesg%buffer(n))
  allocate (mesg%reqs(self%n))
  mesg%id = 1
  mesg%nbuf = n
  io_rank_info%ok    = info%ok
  io_rank_info%rank  = info%rank
  io_rank_info%cost  = info%cost
  io_rank_info%nq    = info%nq
  io_rank_info%n_swap= info%n_swap
  io_rank_info%time  = info%time
  io_rank_info%otime = info%otime
  io_rank_info%dtime = info%dtime
  io_rank_info%wall  = info%wall
  io_rank_info%owall = info%owall
  call anonymous_copy (n, io_rank_info, mesg%buffer)
  if (io%verbose>0) then
    write (io_unit%log,'(a,i6,1p,g12.3,2i6,g16.6)') 'rank_info_mod::pack rank,cost,nq,n,time =', &
       info%rank, info%cost, info%nq, n, info%time
    flush (io_unit%log)
  end if
  call trace_end
END SUBROUTINE pack

!> =============================================================================
!> Unpack an mpi_mesg into a rank_info data type
!> =============================================================================
SUBROUTINE unpack (self, buffer)
  class(rank_info_list_t):: self
  integer, dimension(:), pointer:: buffer
  class(rank_info_t), pointer:: nbor_info
  type(io_rank_info_t):: io_rank_info
  integer:: n
  !.............................................................................
  call trace_begin ('rank_info_mod::pack')
  n = sizeof(io_rank_info)/4
  call anonymous_copy (n, buffer, io_rank_info)
  if (io_rank_info%rank==mpi%rank) return
  !-----------------------------------------------------------------------------
  ! Search for the relevant nbor_info
  !-----------------------------------------------------------------------------
  !$omp critical (rank_info_list_cr)
  nbor_info => self%head
  do while (associated(nbor_info))
    if (nbor_info%rank == io_rank_info%rank) exit
    nbor_info => nbor_info%next
  end do
  !$omp end critical (rank_info_list_cr)
  !-----------------------------------------------------------------------------
  ! If no nbor_info was found, allocate a new one and append
  !-----------------------------------------------------------------------------
  if (.not.associated(nbor_info)) then
    allocate (nbor_info)
    call self%append (nbor_info)
  end if
  !-----------------------------------------------------------------------------
  ! Copy over info
  !-----------------------------------------------------------------------------
  nbor_info%ok    = io_rank_info%ok
  nbor_info%rank  = io_rank_info%rank
  nbor_info%cost  = io_rank_info%cost
  nbor_info%nq    = io_rank_info%nq
  nbor_info%n_swap= io_rank_info%n_swap
  nbor_info%time  = io_rank_info%time
  nbor_info%otime = io_rank_info%otime
  nbor_info%dtime = io_rank_info%dtime
  nbor_info%wall  = io_rank_info%wall
  nbor_info%owall = io_rank_info%owall
  if (io%verbose>1) then
    write (io_unit%log,'(a,i6,1p,g12.3,i6,g16.6)') 'rank_info_mod::unpack, rank,load,nq,n =', &
      io_rank_info%rank, io_rank_info%cost, io_rank_info%nq, io_rank_info%time
  end if
  call trace_end
END SUBROUTINE unpack

!> =============================================================================
!> Append to the rank_info_list -- called from init_nbors to initialize list
!> =============================================================================
SUBROUTINE add (self, rank)
  class(rank_info_list_t):: self
  class(rank_info_t), pointer:: nbor_info
  integer:: rank
  logical:: found
  !.............................................................................
  if (rank==mpi%rank) then
    print *,mpi%rank,'WARNING: trying to add same rank'
    return
  end if
  call trace_begin ('rank_info_mod::add')
  found = .false.
  !$omp critical (rank_info_list_cr)
  nbor_info => rank_info_list%head
  do while (associated(nbor_info))
    if (nbor_info%rank == rank) then
      found = .true.
      exit
    end if
    nbor_info => nbor_info%next
  end do
  !$omp end critical (rank_info_list_cr)
  if (.not. found) then
    allocate (nbor_info)
    nbor_info%rank = rank
    call rank_info_list%append (nbor_info)
    if (io%verbose>1) then
      write (io_unit%log,*) 'rank_info_mod::add', rank, self%n
    end if
  end if
  call trace_end
END SUBROUTINE add

!> =============================================================================
!> Append to the rank_info_list
!> =============================================================================
SUBROUTINE append (self, rank_info)
  class(rank_info_list_t):: self
  class(rank_info_t), pointer:: rank_info
  !.............................................................................
  call trace_begin ('rank_info_mod::append')
  !$omp critical (rank_info_list_cr)
  if (associated(self%head)) then
    self%tail%next => rank_info
  else 
    self%head => rank_info
  end if
  self%tail => rank_info
  self%n = self%n+1
  !$omp end critical (rank_info_list_cr)
  call trace_end
END SUBROUTINE append

!> =============================================================================
!> Remove an item from the rank_info_list.
!> =============================================================================
SUBROUTINE remove (self, this)
  class(rank_info_list_t):: self
  class(rank_info_t), pointer:: this, rank_info, prev
  !.............................................................................
  call trace_begin ('rank_info_mod::remove')
  !$omp critical (rank_info_list_cr)
  nullify (prev)
  rank_info => self%head
  do while (associated(rank_info))
    if (associated(rank_info,this)) then
      !-------------------------------------------------------------------------
      ! This item is not head, so jump over it from prev
      !-------------------------------------------------------------------------
      if (associated(prev)) then
        prev%next => this%next
      !-------------------------------------------------------------------------
      ! This item is head, and if this is the last item, head becomes null
      !-------------------------------------------------------------------------
      else
        self%head => this%next
      endif
      !-------------------------------------------------------------------------
      ! This item is tail, and if this is the last item, tail becomes null
      !-------------------------------------------------------------------------
      if (associated(this,self%tail)) then
        self%tail => prev
      end if
      self%n = self%n-1
      exit
    end if
    prev => rank_info
    rank_info => rank_info%next
  end do
  !-----------------------------------------------------------------------------
  ! If this item was not in the list, it is nevertheless quietly deallocated
  !-----------------------------------------------------------------------------
  deallocate (this)
  !$omp end critical (rank_info_list_cr)
  call trace_end
END SUBROUTINE remove

!> =============================================================================
!> =============================================================================
SUBROUTINE send (self)
  class(rank_info_list_t):: self
  class(rank_info_t), pointer:: nbor_info
  class(mesg_t), pointer:: mesg
  !.............................................................................
  if (.not. associated(self%head)) return
  call trace_begin('rank_info_mod::send')
  !-----------------------------------------------------------------------------
  ! Send to all ranks in the rank_info list (which is sorted by rank)
  !-----------------------------------------------------------------------------
  call self%pack (rank_info, mesg)                              ! pack into mesg
  !$omp critical (rank_info_list_cr)
  nbor_info => self%head                                        ! first rank_info
  do while (associated(nbor_info))                              ! until end
    if (io%verbose>1) then
      write (io_unit%log,'(f12.6,2x,a,2i5,1p,e12.3)') &
        wallclock(), 'rank_info_mod::send', nbor_info%rank, rank_info%rank, rank_info%cost
      flush (io_unit%log)
    end if
    call mesg%send (nbor_info%rank)                             ! send it
    nbor_info => nbor_info%next                                 ! next rank_info
  end do
  !$omp end critical (rank_info_list_cr)
  call mpi_mesg%sent (mesg)                                     ! add to sent_list
  call trace_end
END SUBROUTINE send

!> =============================================================================
!> Interpolate the local code time to the wall clock time nbor%wall (the latest
!> known for nbor), and return the difference between that and nbor%time. A
!> positive value means that nbor is ahead of the local rank.
!> =============================================================================
FUNCTION imbalance (self, nbor) RESULT (diff)
  class(rank_info_t):: self
  class(rank_info_t), pointer:: nbor
  real(8):: p, time, diff
  !.............................................................................
  !p = (nbor%wall-self%owall)/(self%wall-self%owall)
  !time = self%otime + p*(self%time-self%otime)
  !diff = (nbor%time-time)/self%dtime
  diff = (self%nq-nbor%nq)
END FUNCTION imbalance

!> =============================================================================
!> Compute a measure of the task load, which for patches is the number of cells
!> divided by the time step.
!> =============================================================================
SUBROUTINE measure_load (self, head)
  class (rank_info_t):: self
  class (link_t), pointer:: head
  class (link_t), pointer:: link
  class (task_t), pointer:: task
  real:: load, sum, sum_cost, cells, sum_cells, ready
  real(8):: dt, sum_dt, wc
  integer:: n, n1
  integer, save:: itimer
  logical:: ok, active
  !-----------------------------------------------------------------------------
  call trace_begin ('rank_info_mod::measure_load', itimer=itimer)
  !
  ok = .true.
  link => head
  do while (associated(link))
    ok = merge(.false., ok, link%task%rank==mpi%rank .and. link%task%dtime==0d0)
    link => link%next
  end do
  rank_info%ok = ok                             ! active patches have non-zero dtime
  if (.not. ok) then
    rank_info%dtime = 1.0
    call trace_end (itimer)
    return
  end if
  !
  sum_dt = 0.0
  sum_cost = 0.0
  sum_cells = 0.0
  link => head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      if (task%rank==mpi%rank) then
        cells = product(task%mesh%n)
        sum_cells = sum_cells + cells
        sum_cost = sum_cost + cells/task%dtime
      end if
    end select
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! cost = number of cell updates per code time unit
  !-----------------------------------------------------------------------------
  rank_info%cost = sum_cost                     ! total rank update cost
  rank_info%rank = mpi%rank                     ! local rank
  rank_info%dtime = sum_cells/sum_cost          ! cell-weighted dtime
  wc = wallclock()
  if (wc > rank_info%wall) then
    rank_info%owall = rank_info%wall            ! previous wall time
    rank_info%wall  = wc                        ! wall time
    rank_info%otime = rank_info%time            ! previsou queue time
    rank_info%time  = 0.9*rank_info%time  &     ! smoothed ..
                    + 0.1*head%task%time        ! .. queue time
  end if
  if (io%verbose > 1) &
    write (io_unit%log,'(a,3i5,1p,3e12.3)') 'measure_load: n,n1,nq,load,cost,dt =', &
      n, n1, rank_info%nq, rank_info%cost, rank_info%dtime
  call trace_end (itimer)
END SUBROUTINE measure_load

!> =============================================================================
!> Look up the load for a neighbor rank
!> =============================================================================
FUNCTION find (self, rank, debug) RESULT (nbor_info)
  class(rank_info_list_t):: self
  integer:: rank
  logical, optional:: debug
  class(rank_info_t), pointer:: nbor_info
  !.............................................................................
  nbor_info => self%head
  if (present(debug)) then
    write (io_unit%log,*) 'debug: associated =', associated(nbor_info)
    do while (associated(nbor_info))
      write (io_unit%log,*) 'debug: rank =', nbor_info%rank
      if (nbor_info%rank == rank) then
        return
      end if
      nbor_info => nbor_info%next
    end do
    flush (io_unit%log)
  else
    do while (associated(nbor_info))
      if (nbor_info%rank == rank) then
        return
      end if
      nbor_info => nbor_info%next
    end do
  end if
END FUNCTION find

!> =============================================================================
!> Print a rank_info list
!> =============================================================================
SUBROUTINE print (self, time)
  class(rank_info_list_t):: self
  real(8):: time
  class(rank_info_t), pointer:: nbor_info
  real:: load_diff, imbalance
  !.............................................................................
  if (.not.rank_info%ok) return
  !$omp critical (rank_info_list_cr)
  nbor_info => self%head
  do while (associated(nbor_info))
    load_diff = (rank_info%cost-nbor_info%cost)/rank_info%patch_cost
   !imbalance = (nbor_info%time-rank_info%time)/rank_info%dtime
    imbalance = rank_info%imbalance (nbor_info)
    write (io_unit%log,'(2f12.6,2x,a,i6,1p,4g12.3,3i6)') &
      wallclock(), &
      time, &
      'rank_info_list: rnk,load,time,cost[12],wall,nq =', &
      nbor_info%rank,load_diff,imbalance, &
      nbor_info%cost, rank_info%cost, &
      nbor_info%nq, rank_info%nq, nbor_info%n_swap
    nbor_info => nbor_info%next
  end do
  !$omp end critical (rank_info_list_cr)
END SUBROUTINE print

END MODULE rank_info_mod
