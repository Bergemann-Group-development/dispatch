!> =============================================================================
!> $Id$
!> MPI calls related to coordinate systems
!> =============================================================================
MODULE mpi_coords
  USE mpi_basic
  implicit none
  logical:: cart_created=.false.! Set to true if cartesian mpi coordinates are set up
  logical:: mpi_reorder         ! Reordering of ranks allowed, or not
  logical:: mpi_periodic(3)     ! Periodic MPI dimensions, or not
  integer:: mpi_dim3(3)         ! Cartesian MPI dimensions
  integer:: mpi_coord(3)        ! Our cordinates in Cartesian MPI space
  integer:: mpi_plane(3)        ! Communicator btw ranks in planes
  integer:: mpi_beam(3)         ! Communicator btw ranks along axis directions
  integer:: mpi_dn(3)           ! Ranks downwards (circular)
  integer:: mpi_up(3)           ! Ranks upwards (circular)
  integer:: mpi_comm_cart       ! Communicator for cartesian ordering
  type, extends(mpi_t):: mpi_coord_t
    integer:: coord(3)
    integer:: dims(3)
    integer:: plane(3)
    integer:: beam(3)
    integer:: dn(3)
    integer:: up(3)
  contains
    procedure:: init
    procedure:: print
  end type
  logical:: ok
PUBLIC mpi_coord_t
CONTAINS

!> =============================================================================
!> Extent the basic MPI object with MPI coordinate info
!> =============================================================================
SUBROUTINE init (mpi, mpi_dims, dims)
  class(mpi_coord_t):: mpi
  integer, optional:: mpi_dims(3), dims(3)
  interface
    subroutine mpi_create_mpi (mpi_dims, dims)
      integer, dimension(3), optional:: mpi_dims, dims
    end subroutine
  end interface
  !.............................................................................
  call mpi%mpi_t%init
  call cart_create_mpi (mpi_dims, dims=dims)
  mpi%coord = mpi_coord
  mpi%dims  = mpi_dim3
  mpi%plane = mpi_plane
  mpi%beam  = mpi_beam
  mpi%dn    = mpi_dn
  mpi%up    = mpi_up
  mpi%ok    = ok
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE print (mpi)
  class(mpi_coord_t):: mpi
  integer:: rank
  !.............................................................................
  if (mpi%ok) then
    call mpi%mpi_t%print
    do rank=0,mpi%size-1
      if (rank==mpi%rank) then
        print *,'...............................................'
        print *, 'mpi%rank  =', mpi%rank
        print *, 'mpi%dims  =', mpi%dims
        print *, 'mpi%coord =', mpi%coord
        print *, 'mpi%beam  =', mpi%beam
        print *, 'mpi%plane =', mpi%plane
        print *, 'mpi%dn    =', mpi%dn
        print *, 'mpi%up    =', mpi%up
      end if
      call mpi%barrier (delay=1)
    end do
  else
    print *, "WARNING: no MPI cartesian coordinates"
  end if
END SUBROUTINE

!> =============================================================================
!> Create a Cartesian arrangement of MPI ranks, with dimension n_mpi(3).  If
!> dims(3) is non-zero, these are taken as information about the global number
!> of mesh points, and an attempt is made to distribute the MPI-coordinate
!> dimension to get product(dims) = mpi_size.  If mpi_dims is present, the
!> dims(3) information is ignored, even if present.
!> =============================================================================
SUBROUTINE cart_create_mpi (mpi_dims, dims)
  implicit none
  integer, dimension(3), optional:: mpi_dims, dims
  integer:: i, m, color, size, n(3)
  character(len=mch), save:: id= &
    'mpi_coords.f90 $Id$'
!-------------------------------------------------------------------------------
  mpi_reorder = .true.
  mpi_periodic = .true.
  !-----------------------------------------------------------------------------
  ! Create a communicator (cart) with Cartesian ordering
  !-----------------------------------------------------------------------------
  ok = .false.
  if (present(mpi_dims)) then
    mpi_dim3 = mpi_dims
    if (master) print'(a,i7,2x,3i5)',' mpi_size, mpi_dims =', mpi_size, mpi_dims
  else
    if (present(dims)) then
      if (product(dims)/mpi_size*mpi_size == product(dims)) then
        do while (.true.)
          m = mpi_size**0.333334
          if (m**3 == mpi_size) then
            mpi_dim3 = [m,m,m]; if (test_is_ok()) exit
          end if
          m = mpi_size**0.500001
          if (m**2 == mpi_size) then
            mpi_dim3 = [m,m, 1]; if (test_is_ok()) exit
            mpi_dim3 = [m, 1,m]; if (test_is_ok()) exit
            mpi_dim3 = [ 1,m,m]; if (test_is_ok()) exit
          end if
          m = mpi_size
            mpi_dim3 = [m,1,1]; if (test_is_ok()) exit
            mpi_dim3 = [1,m,1]; if (test_is_ok()) exit
            mpi_dim3 = [1,1,m]; if (test_is_ok()) exit
          call fail ('ERROR: none of the simple choices worked')
        end do
      else
        call fail ('ERROR: mpi_size must be a divisor in =', product(dims))
      end if
    else
      call fail ('ERROR: either mpi_dims(3) or dims(3) must be present')
    end if
  end if
  if (mpi_size /= product(mpi_dim3)) then
    call fail ('ERROR: mpi_size must be equal to product(mpi_dims) =', &
      product(mpi_dim3))
  end if
  if (.not. ok) return
  call MPI_Cart_create (mpi_comm_world, 3, mpi_dim3, mpi_periodic, mpi_reorder, &
                        mpi_comm_cart, mpi_err)
  if (master) print '(13x,a,3i5)', ' Using mpi_dims =', mpi_dim3
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
          + merge(0,mpi_coord(2),i==2)*mpi_dim3(1) &
          + merge(0,mpi_coord(3),i==3)*mpi_dim3(1)*mpi_dim3(2)
    call MPI_Comm_split (mpi_comm_world, color, mpi_coord(i), mpi_beam(i), mpi_err)
    call assert_mpi ('cart_creat_mpi: mpi_beam', mpi_err)
  end do
  cart_created = .true.
  CONTAINS
    logical function test_is_ok ()
      n = dims/mpi_dim3
      if (master) print '(13x,a,3i5)', 'trying mpi_dims =', mpi_dim3
      test_is_ok = all(n*mpi_dim3 == dims) .and. product(mpi_dim3)==mpi_size
      ok = test_is_ok
    end function
    subroutine fail (label, i)
      character(len=*):: label
      integer, optional:: i
      if (master) then
        if (present(i)) then
          print *, label, i
        else
          print *, label
        end if
      end if
    end subroutine
END SUBROUTINE cart_create_mpi
END MODULE mpi_coords
!
!> =============================================================================
!> =============================================================================
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
