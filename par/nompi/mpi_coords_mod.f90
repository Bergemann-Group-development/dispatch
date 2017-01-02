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
  USE timer_mod
  implicit none
  private
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
  coords = 0
END FUNCTION rank_to_coords

!> =============================================================================
!> Set the rank from the coordinates
!> =============================================================================
FUNCTION coords_to_rank (self, coords) result (rank)
  class(mpi_coords_t):: self
  integer:: coords(3), rank
  !.............................................................................

END FUNCTION coords_to_rank

!> =============================================================================
!> Extent the basic MPI object with MPI coordinate info
!> =============================================================================
SUBROUTINE init (self, mpi_dims, dims)
  class(mpi_coords_t):: self
  integer, optional:: mpi_dims(3), dims(3)
  interface
    subroutine mpi_create_mpi (mpi_dims, dims)
      integer, dimension(3), optional:: mpi_dims, dims
    end subroutine
  end interface
  !.............................................................................
  call self%mpi_t%init
  call cart_create_mpi (self, mpi_dims, dims=dims)
  self%coords = mpi_cord3
  self%dims   = mpi_dim3
  self%plane  = mpi_plane
  self%beam   = mpi_beam
  self%dn     = mpi_dn
  self%up     = mpi_up
  self%ok     = ok
  mpi_coords%cart_comm = self%cart_comm
  mpi_coords%coords = self%coords
  mpi_coords%dims   = self%dims
  mpi_coords%plane  = self%plane
  mpi_coords%beam   = self%beam
  mpi_coords%dn     = self%dn
  mpi_coords%up     = self%up
  mpi_coords%ok     = self%ok
  if (mpi%size /= product(mpi_dim3)) then
    print*,'inconstent mpi%size:', mpi%size, mpi_dim3
    stop
  end if
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE print (self)
  class(mpi_coords_t):: self
  integer:: rank
  !.............................................................................
  if (self%ok) then
    call self%mpi_t%print
    do rank=0,self%size-1
      if (rank==self%rank) then
        print *,'...............................................'
        print *, 'mpi%rank   =', self%rank
        print *, 'mpi%dims   =', self%dims
        print *, 'mpi%coords =', self%coords
        print *, 'mpi%beam   =', self%beam
        print *, 'mpi%plane  =', self%plane
        print *, 'mpi%dn     =', self%dn
        print *, 'mpi%up     =', self%up
      end if
      call mpi%barrier (delay=0.25)
    end do
  else if (self%master) then
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
SUBROUTINE cart_create_mpi (self, mpi_dims, dims)
  implicit none
  class(mpi_coords_t):: self
  integer, dimension(3), optional:: mpi_dims, dims
  integer:: i, m, color, size, n(3)
  integer:: mpi_err
  character(len=120), save:: id= &
    'mpi_coords.f90 $Id$'
  integer, save:: itimer=0
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
      ok = product(mpi_dims)==self%size
    end if
    if (self%rank==0) print'(a,i7,2x,3i5)',' mpi%size, mpi_dims =', self%size, mpi_dims
  else
    if (present(dims)) then
      if (product(dims)/self%size*self%size == product(dims)) then
        do while (.true.)
          m = self%size**0.333334
          if (m**3 == self%size) then
            mpi_dim3 = [m,m,m]; if (test_is_ok()) exit
          end if
          m = self%size**0.500001
          if (m**2 == self%size) then
            mpi_dim3 = [m,m, 1]; if (test_is_ok()) exit
            mpi_dim3 = [m, 1,m]; if (test_is_ok()) exit
            mpi_dim3 = [ 1,m,m]; if (test_is_ok()) exit
          end if
          m = self%size
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
  if (self%size /= product(mpi_dim3)) then
    call fail ('WARNING: mpi%size must be equal to product(mpi_dims) =', &
      product(mpi_dim3))
  end if
  if (.not. ok) return

  self%cart_comm = mpi_comm_cart
  if (self%rank==0) print '(13x,a,3i5)', ' Using mpi_dims =', mpi_dim3

  !-----------------------------------------------------------------------------
  ! Compute mpi_cord3; our coordinates in Cartesian MPI space
  !-----------------------------------------------------------------------------

  do i=1,3
    !---------------------------------------------------------------------------
    ! Compute mpi_up and mpi_dn; the ranks of nearest neighbors
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute mpi_plane; the communicators for planes in MPI-space
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Compute mpi_beam; the communicators along MPI coordinat axes
    !---------------------------------------------------------------------------
    color = merge(0,mpi_cord3(1),i==1) &
          + merge(0,mpi_cord3(2),i==2)*mpi_dim3(1) &
          + merge(0,mpi_cord3(3),i==3)*mpi_dim3(1)*mpi_dim3(2)

  end do
  cart_created = .true.
  CONTAINS
    logical function test_is_ok ()
      n = dims/mpi_dim3
      test_is_ok = all(n*mpi_dim3 == dims) .and. product(mpi_dim3)==self%size
      ok = test_is_ok
      if (self%rank==0) print '(13x,a,3i5,l5)', 'trying mpi_dims =', mpi_dim3, ok
    end function
    subroutine fail (label, i)
      character(len=*):: label
      integer, optional:: i
      if (self%rank==0) then
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
SUBROUTINE finalize (self)
  implicit none
  class(mpi_coords_t):: self
  integer:: mpi_err, i
  !

END SUBROUTINE finalize
END MODULE mpi_coords_mod
