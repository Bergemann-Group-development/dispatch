!> =============================================================================
!> $Id$
!> MPI calls related to coordinate systems
!> =============================================================================
MODULE mpi_coords
  USE mpi_basic
  implicit none
PUBLIC
  logical:: cart_created=.false.! Set to true if cartesian mpi coordinates are set up
  logical:: mpi_reorder         ! Reordering of ranks allowed, or not
  logical:: mpi_periodic(3)     ! Periodic MPI dimensions, or not
  integer:: mpi_dims(3)         ! Cartesian MPI dimensions
  integer:: mpi_coord(3)        ! Our cordinates in Cartesian MPI space
  integer:: mpi_plane(3)        ! Communicator btw ranks in planes
  integer:: mpi_beam(3)         ! Communicator btw ranks along axis directions
  integer:: mpi_dn(3)           ! Ranks downwards (circular)
  integer:: mpi_up(3)           ! Ranks upwards (circular)
  integer:: mpi_comm_cart       ! Communicator for cartesian ordering
CONTAINS

!> =============================================================================
!> Create a Cartesian arrangement of MPI ranks
!> =============================================================================
SUBROUTINE cart_create_mpi (n_mpi, gn)
  implicit none
  integer, dimension(3):: n_mpi, gn, n
  integer:: i, color, size
  character(len=mch), save:: fmt = '(13x,a,3i5)'
  character(len=mch), save:: id= &
    'mpi_coords.f90 $Id: 5b66abdd71074a2181f902f590868f350cef1c30 $'
!-------------------------------------------------------------------------------
  mpi_dims = n_mpi
  mpi_reorder = .true.
  mpi_periodic = .true.
  if (master) print'(a,i7,2x,3i5)',' mpi_size, mpi_dims =', mpi_size, mpi_dims
  !-----------------------------------------------------------------------------
  ! Create a communicator (cart) with Cartesian ordering
  !-----------------------------------------------------------------------------
  if (mpi_size /= product(mpi_dims)) then
     if (product(gn)/mpi_size*mpi_size == product(gn)) then
       do i=1,16
         mpi_dims = (/i,mpi_size,i/); n = gn/mpi_dims
         if (master) print fmt,'trying mpi_dims =', mpi_dims
         if (all(n*mpi_dims == gn) .and. all(mpi_dims==1 .or. n>=5)) goto 1
       end do
       if (master) print *,'ERROR: could not find a valid MPI distribution'
       call end_mpi
     else
       if (master) print *,'ERROR: mpi_size must be a divisor in =', product(gn)
       call end_mpi
     end if
  end if
1 if (master) print fmt,' Using mpi_dims =', mpi_dims
  call MPI_Cart_create (mpi_comm_world, 3, mpi_dims, mpi_periodic, mpi_reorder, &
                        mpi_comm_cart, mpi_err)
  call assert_mpi ('cart_creat_mpi: cart', mpi_err)
  !-----------------------------------------------------------------------------
  ! Compute mpi_coord; our coordinates in Cartesian MPI space
  !-----------------------------------------------------------------------------
  call MPI_Cart_coords (mpi_comm_cart, mpi_rank, 3, mpi_coord, mpi_err)                  ! mpi_coord(i) = mpi_[xyz]
  call assert_mpi ('cart_creat_mpi: mpi_coord', mpi_err)
  do i=1,3
    !---------------------------------------------------------------------------
    ! Compute mpi_up and mpi_dn; the ranks of nearest neighbors
    !---------------------------------------------------------------------------
    call MPI_Cart_shift (mpi_comm_cart, i-1, 1, mpi_dn(i), mpi_up(i), mpi_err)           ! mpi_{up,dn} ranks
    call assert_mpi ('cart_creat_mpi: mpi_dn/up', mpi_err)
    !---------------------------------------------------------------------------
    ! Compute mpi_plane; the communicators for planes in MPI-space
    !---------------------------------------------------------------------------
    call MPI_Comm_split (mpi_comm_world, mpi_coord(i), mpi_rank, mpi_plane(i),& ! mpi_plane(
                         mpi_err)
    call assert_mpi ('cart_creat_mpi: mpi_plane', mpi_err)
    !---------------------------------------------------------------------------
    ! Compute mpi_beam; the communicators along MPI coordinat axes
    !---------------------------------------------------------------------------
    color = merge(0,mpi_coord(1),i==1) &
          + merge(0,mpi_coord(2),i==2)*mpi_dims(1) &
          + merge(0,mpi_coord(3),i==3)*mpi_dims(1)*mpi_dims(2)
    call MPI_Comm_split (mpi_comm_world, color, mpi_coord(i), mpi_beam(i), mpi_err)
    call assert_mpi ('cart_creat_mpi: mpi_beam', mpi_err)
  end do
  cart_created = .true.
END SUBROUTINE cart_create_mpi

END MODULE mpi_coords
!
SUBROUTINE finalize_cart
  USE mpi_coords
  implicit none
  integer i
  !
  if (cart_created) then
    call MPI_Comm_Free(mpi_comm_cart,mpi_err)
    do i=1,3
      call MPI_Comm_Free(mpi_plane(i),mpi_err)
      call assert_mpi ('finalize_cart: mpi_plane',mpi_err)
      call MPI_Comm_Free(mpi_beam(i),mpi_err)
      call assert_mpi ('finalize_cart: mpi_beam',mpi_err)
    enddo
  endif
END SUBROUTINE finalize_cart
