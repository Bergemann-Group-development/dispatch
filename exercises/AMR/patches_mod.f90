!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Generic module for MHD patch.  To be overloaded by an MHD layer and a timestep
!> layer.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patches_mod
  USE trace_mod
  USE io_mod
  USE bits_mod
  USE central_mod
  USE upwind_mod
  USE derivatives_mod
  implicit none
  integer, parameter:: ntime=7
  private
  type, public:: patch_t
    integer:: id = 0
    integer:: nt, nv
    integer:: n(3), ng(3)=2, gn(3), li(3), ui(3)
    integer:: it=1, new=2
    integer:: i, j, k
    integer:: status
    integer, pointer:: iit(:)
    real, dimension(3):: position, size, gsize, box, ds, origin
    real, pointer:: t(:)
    real, pointer:: dt(:)
    real:: time=0.0, dtime=0.1, end_time=1.0, grace=0.0
    real:: out_time=0.1, out_next=0.0           ! controls output
    real:: courant=0.2                          ! courant constant
    integer:: istep=0                           ! time step counter
    integer:: iout=0                            ! output counter
    real, dimension(:,:,:,:,:,:), pointer:: mem
    integer:: verbose = 0
    logical:: is_leaf = .true.
    logical:: is_allocated = .false.
    class(central_t), pointer:: derivs
  contains
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: distance        ! return the distance btw two patches
    procedure:: overlaps        ! return true if they overlap (incl. ghost zones)
    procedure:: index_range     ! returns the index ranges btw two points
    procedure:: indices         ! computes indices and fractional indices
    procedure:: index           ! compute the index for a single point
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: ahead_of
  end type
  integer:: id = 0
CONTAINS
!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self, position, size)
  class(patch_t):: self
  integer:: dims(5) = 4
  real:: position(3), size(3)
  logical:: first_time = .true.
  namelist /patch_params/ dims, size
  !.............................................................................
  !$omp atomic
  id = id + 1                                   ! increment global
  self%id = id                                  ! new patch id
  self%position = position                      ! as per call
  call trace_begin ('patch_t%init', 3)
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
  end if
  self%size = size                              ! patch size
  self%n  = dims(1:3)                           ! patch dimensions
  where (self%n==1)                             ! for collapsed dimensions
    self%ng = 0                                 ! no guard zones
    self%size = 0.0                             ! zero size
  end where
  self%gn = self%n + 2*self%ng                  ! including guard cells
  self%li = 1       + self%ng                   ! lower interior index
  self%ui = self%gn - self%ng                   ! upper interior index
  self%nt = dims(4)                             ! number of time slices
  self%nv = dims(5)                             ! number of variables
  allocate (self%derivs)                        ! thread safe derivatives
  call self%derivs%init (dims(1:3), size/self%n)! pass on patch info to derivs
  if (self%is_allocated) then
    call trace_end
    return
  end if
  self%is_allocated = .true.                    ! mark
  allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
  self%t = 0.0
  self%dt = 1.0
  self%iit = 1
  associate (nx=>self%gn(1), ny=>self%gn(2), nz=>self%gn(3))
    print*,'mem:',nx,ny,nz,self%nv,self%nt
    allocate (self%mem(nx,ny,nz,self%nv,self%nt,2))
  end associate
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE update (self)
  class(patch_t):: self
  !.............................................................................
  call trace_begin ('patch_t%update', 3)
  if (self%time >= self%out_time) then
    call self%output
    self%out_next = self%out_next + self%out_time
  end if
  self%istep = self%istep + 1
  self%dt(self%it) = self%dtime
  self%t(self%new) = self%time
  self%it = self%new
  self%new = mod(self%new,ntime)+1
  self%t(self%new) = self%time
  call trace_end
END SUBROUTINE update
!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self)
  class(patch_t):: self
  !.............................................................................
  if (io%master) print*,'patch_t%output: should do output at t =', self%time
END SUBROUTINE output
!> =============================================================================
FUNCTION index_range (self, position)
  class(patch_t):: self
  real:: position(3)
  integer:: index_range(3,2)
  integer:: dir
  !.............................................................................
  do dir=1,3
    index_range(:,dir) = int(modulo(position(dir)-self%origin, self%box/self%ds)) + 1
    index_range(:,dir) = max(1, min(self%gn, index_range(:,dir)))
  end do
END FUNCTION index_range
!> =============================================================================
SUBROUTINE indices (self, np, p, ii, f)
  class(patch_t):: self
  real:: p(np), f(np,3)
  integer:: np, dir, ii(np,3)
  !.............................................................................
  do dir=1,3
    f(:,dir) = modulo(p-self%origin(dir), self%box(dir))/self%ds(dir)
    ii(:,dir) = max(1, min(self%gn(dir), int(f(:,dir))))
  end do
  f = f-ii
END SUBROUTINE indices
!> =============================================================================
FUNCTION index (self, position)
  class(patch_t):: self
  real, dimension(3):: position, f
  integer:: index(3)
  !.............................................................................
  f = modulo(position-self%origin, self%box)/self%ds + 1.0
  index = max(1, min(self%gn, int(f)))
  f = f - index
END FUNCTION index
!> =============================================================================
FUNCTION distance (self, patch)
  class(patch_t):: self, patch
  real:: distance(3)
  distance = self%position-patch%position+0.5*self%box
  distance = abs(modulo(distance+0.5*self%box,self%box)-0.5*self%box)
END FUNCTION distance
!> =============================================================================
FUNCTION overlaps (self, patch)
  class(patch_t):: self, patch
  logical:: overlaps
  !.............................................................................
  overlaps = all(distance(self,patch) <= 0.5000000001_8*(self%gsize+patch%size))
END FUNCTION

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%set')
  self%status = ior(self%status, bit)
!  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%clear')
  self%status = iand(self%status, not(bit))
!  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%is_set')
  is_set = (iand(self%status, bit) /= 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION ahead_of (self, patch)
  class (patch_t):: self, patch
  !.............................................................................
  ahead_of = self%time > patch%time + patch%grace
!  print*,'ahead_of: ', self%time, patch%time, patch%grace, ahead_of
END FUNCTION
END MODULE
