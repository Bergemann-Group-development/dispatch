!> =============================================================================
!> $Id$
!> MPI calls related to Cartesian MPI coordinate systems
!>
!> To make use of this module, make sure it is compiled together with your code,
!> using either the Makefile and the 'make' command, or compiling manually, with
!>
!> mpi = ../../../../../mpi             # or wherever
!> mpifort -c $mpi/mpi_mod.f90          # compiler mpi_mod.f90
!> mpifort -c $mpi/mpi_coords.f90       # add any other module your need
!> mpifort -c your_code.f90             # compile your code
!> mpifort *.o -o your_code.x           # link together into your_code.x
!>
!> In your_code.f90 you add lines such as:
!>
!> USE mpi                              ! makes the mpi% object available
!> USE mpi_coords                       ! makes the mpi_coords% object available
!> implicit none
!> integer:: dims(3)                    ! MPI cartesian dimension
!> ...
!> call mpi%init                        ! initializes the mpi%object
!> call mpi_coords%init (dims=dims)     ! initializes the mpi_coords%object
!> call mpi_coords%print                ! prints an overview of the geometry
!> ...
!> ...
!> call mpi%end                         ! closes MPI
!>
!> Too see which variables and procedures are inside mpi%, just look below!
!>
!> =============================================================================
MODULE mpi_coords_mod
  USE mpi_mod
  implicit none
  private
  include 'mpif.h'
  logical:: cart_created=.false.! Set to true if cartesian mpi coordinates are set up
  logical:: mpi_reorder         ! Reordering of ranks allowed, or not
  logical:: mpi_periodic(3)     ! Periodic MPI dimensions, or not
  integer:: mpi_dim3(3)         ! Cartesian MPI dimensions
  integer:: mpi_cord3(3)        ! Our cordinates in Cartesian MPI space
  integer:: mpi_plane(3)        ! Communicator btw ranks in planes
  integer:: mpi_beam(3)         ! Communicator btw ranks along axis directions
  integer:: mpi_dn(3)           ! Ranks downwards (circular)
  integer:: mpi_up(3)           ! Ranks upwards (circular)
  integer:: mpi_comm_cart       ! Communicator for cartesian ordering
  logical:: ok
  integer:: mpi_err
  type, public, extends(mpi_t):: mpi_coords_t
    integer:: coords(3)
    integer:: dims(3)
    integer:: plane(3)
    integer:: beam(3)
    integer:: dn(3)
    integer:: up(3)
    integer:: cart_comm
  contains
    procedure:: init
    procedure:: print
    procedure:: coords_to_rank
    procedure:: rank_to_coords
    procedure:: finalize
  end type
  type(mpi_coords_t), public:: mpi_coords
CONTAINS

!> =============================================================================
!> Set the coordinates from the rank
!> =============================================================================
FUNCTION rank_to_coords (self, rank) result (coords)
  class(mpi_coords_t):: self
  integer:: rank, coords(3)
  !.............................................................................
  if (self%size == 1) then
    coords = 0
  else
    call MPI_CART_COORDS (self%cart_comm, rank, 3, coords, mpi_err)
    call mpi%assert ('MPI_CART_COORDS', mpi_err)
  end if
END FUNCTION rank_to_coords

!> =============================================================================
!> Set the rank from the coordinates
!> =============================================================================
FUNCTION coords_to_rank (self, coords) result (rank)
  class(mpi_coords_t):: self
  integer:: coords(3), rank
  !.............................................................................
  call MPI_CART_RANK (self%cart_comm, coords, rank, mpi_err)
  call mpi%assert ('MPI_CART_RANK', mpi_err)
END FUNCTION coords_to_rank

!> =============================================================================
!> Extent the basic MPI object with MPI coordinate info
!> =============================================================================
SUBROUTINE init (mp, mpi_dims, dims)
  class(mpi_coords_t):: mp
  integer, optional:: mpi_dims(3), dims(3)
  interface
    subroutine mpi_create_mpi (mpi_dims, dims)
      integer, dimension(3), optional:: mpi_dims, dims
    end subroutine
  end interface
  !.............................................................................
  call mp%mpi_t%init
  call cart_create_mpi (mp, mpi_dims, dims=dims)
  mp%coords = mpi_cord3
  mp%dims   = mpi_dim3
  mp%plane  = mpi_plane
  mp%beam   = mpi_beam
  mp%dn     = mpi_dn
  mp%up     = mpi_up
  mp%ok     = ok
  mpi_coords%cart_comm = mp%cart_comm
  mpi_coords%coords = mp%coords
  mpi_coords%dims   = mp%dims
  mpi_coords%plane  = mp%plane
  mpi_coords%beam   = mp%beam
  mpi_coords%dn     = mp%dn
  mpi_coords%up     = mp%up
  mpi_coords%ok     = mp%ok
  if (mpi%size /= product(mpi_dim3)) then
    print*,'inconstent mpi%size:', mpi%size, mpi_dim3
    stop
  end if
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE print (mp)
  class(mpi_coords_t):: mp
  integer:: rank
  !.............................................................................
  if (mp%ok) then
    call mp%mpi_t%print
    do rank=0,mp%size-1
      if (rank==mp%rank) then
        print *,'...............................................'
        print *, 'mpi%rank   =', mp%rank
        print *, 'mpi%dims   =', mp%dims
        print *, 'mpi%coords =', mp%coords
        print *, 'mpi%beam   =', mp%beam
        print *, 'mpi%plane  =', mp%plane
        print *, 'mpi%dn     =', mp%dn
        print *, 'mpi%up     =', mp%up
      end if
      call mpi%barrier (delay=1)
    end do
  else if (mp%master) then
    print *, "WARNING: no MPI cartesian coordinates"
  end if
END SUBROUTINE

!> =============================================================================
!> Create a Cartesian arrangement of MPI ranks, with dimension n_mpi(3).  If
!> dims(3) is non-zero, these are taken as information about the global number
!> of mesh points, and an attempt is made to distribute the MPI-coordinate
!> dimension to get product(dims) = mpi%size.  If mpi_dims is present, the
!> dims(3) information is ignored, even if present.
!> =============================================================================
SUBROUTINE cart_create_mpi (mp, mpi_dims, dims)
  implicit none
  class(mpi_coords_t):: mp
  integer, dimension(3), optional:: mpi_dims, dims
  integer:: i, m, color, size, n(3)
  integer:: mpi_err
  character(len=120), save:: id= &
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
    if (present(dims)) then
      ok = test_is_ok()
    else
      ok = product(mpi_dims)==mp%size
    end if
    if (mp%rank==0) print'(a,i7,2x,3i5)',' mpi%size, mpi_dims =', mp%size, mpi_dims
  else
    if (present(dims)) then
      if (product(dims)/mp%size*mp%size == product(dims)) then
        do while (.true.)
          m = mp%size**0.333334
          if (m**3 == mp%size) then
            mpi_dim3 = [m,m,m]; if (test_is_ok()) exit
          end if
          m = mp%size**0.500001
          if (m**2 == mp%size) then
            mpi_dim3 = [m,m, 1]; if (test_is_ok()) exit
            mpi_dim3 = [m, 1,m]; if (test_is_ok()) exit
            mpi_dim3 = [ 1,m,m]; if (test_is_ok()) exit
          end if
          m = mp%size
            mpi_dim3 = [m,1,1]; if (test_is_ok()) exit
            mpi_dim3 = [1,m,1]; if (test_is_ok()) exit
            mpi_dim3 = [1,1,m]; if (test_is_ok()) exit
          call fail ('WARNING: none of the simple choices worked')
        end do
      else
        call fail ('WARNING: mpi%size must be a divisor in =', product(dims))
      end if
    else
      call fail ('WARNING: either mpi_dims(3) or dims(3) must be present')
    end if
  end if
  if (mp%size /= product(mpi_dim3)) then
    call fail ('WARNING: mpi%size must be equal to product(mpi_dims) =', &
      product(mpi_dim3))
  end if
  if (.not. ok) return
  call MPI_Cart_create (mpi_comm_world, 3, mpi_dim3, mpi_periodic, mpi_reorder, &
                        mpi_comm_cart, mpi_err)
  mp%cart_comm = mpi_comm_cart
  if (mp%rank==0) print '(13x,a,3i5)', ' Using mpi_dims =', mpi_dim3
  call mpi%assert ('cart_creat_mpi: cart', mpi_err)
  !-----------------------------------------------------------------------------
  ! Compute mpi_cord3; our coordinates in Cartesian MPI space
  !-----------------------------------------------------------------------------
  call MPI_Cart_coords (mpi_comm_cart, mp%rank, 3, mpi_cord3, mpi_err)          ! mpi_cord3(i) = mpi_[xyz]
  call mpi%assert ('cart_creat_mpi: mpi_cord3', mpi_err)
  do i=1,3
    !---------------------------------------------------------------------------
    ! Compute mpi_up and mpi_dn; the ranks of nearest neighbors
    !---------------------------------------------------------------------------
    call MPI_Cart_shift (mpi_comm_cart, i-1, 1, mpi_dn(i), mpi_up(i), mpi_err)  ! mpi_{up,dn} ranks
    call mpi%assert ('cart_creat_mpi: mpi_dn/up', mpi_err)
    !---------------------------------------------------------------------------
    ! Compute mpi_plane; the communicators for planes in MPI-space
    !---------------------------------------------------------------------------
    call MPI_Comm_split (mpi_comm_world, mpi_cord3(i), mp%rank, mpi_plane(i),&  ! mpi_plane(
                         mpi_err)
    call mpi%assert ('cart_creat_mpi: mpi_plane', mpi_err)
    !---------------------------------------------------------------------------
    ! Compute mpi_beam; the communicators along MPI coordinat axes
    !---------------------------------------------------------------------------
    color = merge(0,mpi_cord3(1),i==1) &
          + merge(0,mpi_cord3(2),i==2)*mpi_dim3(1) &
          + merge(0,mpi_cord3(3),i==3)*mpi_dim3(1)*mpi_dim3(2)
    call MPI_Comm_split (mpi_comm_world, color, mpi_cord3(i), mpi_beam(i), mpi_err)
    call mpi%assert ('cart_creat_mpi: mpi_beam', mpi_err)
  end do
  cart_created = .true.
  CONTAINS
    logical function test_is_ok ()
      n = dims/mpi_dim3
      test_is_ok = all(n*mpi_dim3 == dims) .and. product(mpi_dim3)==mp%size
      ok = test_is_ok
      if (mp%rank==0) print '(13x,a,3i5,l5)', 'trying mpi_dims =', mpi_dim3, ok
    end function
    subroutine fail (label, i)
      character(len=*):: label
      integer, optional:: i
      if (mp%rank==0) then
        if (present(i)) then
          print *, label, i
        else
          print *, label
        end if
      end if
    end subroutine
END SUBROUTINE cart_create_mpi
!
!> =============================================================================
!> =============================================================================
SUBROUTINE finalize (mp)
  implicit none
  class(mpi_coords_t):: mp
  integer:: mpi_err, i
  !
  if (mp%ok) then
    call MPI_Comm_Free(mpi_comm_cart,mpi_err)
    do i=1,3
      call MPI_Comm_Free(mpi_plane(i),mpi_err)
      call mpi%assert ('finalize_cart: mpi_plane',mpi_err)
      call MPI_Comm_Free(mpi_beam(i),mpi_err)
      call mpi%assert ('finalize_cart: mpi_beam',mpi_err)
    enddo
  endif
END SUBROUTINE finalize
END MODULE mpi_coords_mod
