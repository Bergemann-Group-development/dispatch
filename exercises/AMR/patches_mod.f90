!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Generic module for MHD patch.  To be overloaded by an MHD layer and a timestep
!> layer.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patches_mod
  USE trace_mod
  USE io_mod
  USE bits_mod
  USE central_mod
  implicit none
  integer, parameter:: ntime=7
  private
  type, public:: patch_t
    integer:: id = 0
    integer:: nt = ntime, nv = 5
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
    logical:: mem_allocated = .false.
    class(central_t), pointer:: derivs
  contains
    procedure:: init
    procedure:: update => update_patch
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
    procedure:: is_relevant_for
    procedure:: get_guardcells_for
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
  dims(4) = self%nt
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
  end if
  self%size = size                              ! patch size
  self%n  = dims(1:3)                           ! patch dimensions
  self%ds = real(self%size,kind=8) / self%n     ! ds in patch
  where (self%n==1)                             ! for collapsed dimensions
    self%ng = 0                                 ! no guard zones
    self%size = 0.0                             ! zero size
    self%ds = 0.0                               ! first zero
  end where
  where (self%n==1)                             ! for collapsed dimensions
    self%ds = maxval(self%ds)                   ! then set to reasonable val
  end where
  self%gn = self%n + 2*self%ng                  ! including guard cells
  self%gsize = self%size*real(self%gn,kind=8)/self%n ! size incl guard cells
  self%li = 1       + self%ng                   ! lower interior index
  self%ui = self%gn - self%ng                   ! upper interior index
  self%nt = dims(4)                             ! number of time slices
  self%nv = dims(5)                             ! number of variables
  allocate (self%derivs)                        ! thread safe derivatives
  call self%derivs%init (dims(1:3), size/self%n)! pass on patch info to derivs
  if (self%mem_allocated) then
    call trace_end
    return
  end if
  self%mem_allocated = .true.                   ! mark
  print*,'MK1', associated(self%t) , associated(self%iit), &
                associated(self%dt), associated(self%mem)
  allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
  self%t = 0.0
  self%dt = 1.0
  self%iit = 1
  associate (nx=>self%gn(1), ny=>self%gn(2), nz=>self%gn(3))
    if (io%verbose>2) &
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
SUBROUTINE update_patch (self)
  class(patch_t):: self
  !.............................................................................
  call trace_begin ('patch_t%update', 3)
  if (self%time >= self%out_time) then
    call self%output
    self%out_next = self%out_next + self%out_time
  end if
  self%dt(self%it) = self%dtime                         ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = self%time                          ! initial time
  self%it = self%new                                    ! update time slot
  self%new = mod(self%new,ntime)+1                      ! increment / rotate
  self%iit(1:self%nt-1) = self%iit(2:self%nt)           ! rotate left
  self%iit(self%nt) = self%new                          ! new right-most slot
  self%istep = self%istep + 1
  if (io%verbose>2) &
   print*,'patch id:',self%id,'  it:',self%it,'  new:',self%new
  if (io%verbose>3) &
   print*,'patch iit:', self%iit
  if (io%verbose>3) &
   print*,'patch t(iit):', self%t(self%iit)
  call trace_end
END SUBROUTINE update_patch

!> =============================================================================
SUBROUTINE get_guardcells_for (self, patch)
  class(patch_t):: self
  class(patch_t), pointer:: patch
  !
  integer :: ls(3), us(3), lp(3), up(3), dir
  real(kind=8) :: rlb(3), rub(3), hs
  !
  ! Check if the two patches have same size
  if (maxval(abs(patch%ds-self%ds)/max(patch%ds,self%ds)) < 1e-5) then  ! we are on the same level
    ! Get real lower and upper boundary with possible wrap-around
    rlb = self%position-(self%n/2+1+1e-5_8)*self%ds ! go slightly below lower boundary
    rub = self%position+(self%n/2+1+1e-5_8)*self%ds ! go slightly above upper boundary
    rlb = modulo(rlb + self%box, self%box)          ! wrap around
    rub = modulo(rub + self%box, self%box)
    ! check if we can get guard cells in a given direction
    do dir=1,3
      hs = 0.5_8 * patch%size(dir)
      ! check lower boundary
      if (abs(patch%position(dir) - rub(dir)) < hs) then
        ! We are overlapping in this direction, exchange data
        ls = 1; us=self%gn
        ls(dir) = self%ui(dir) - patch%ng(dir) + 1
        us(dir) = self%ui(dir)

        lp = 1; up=patch%gn
        lp(dir) = 1
        up(dir) = patch%li(dir) - 1
        ! copy data
        patch%mem(lp(1):up(1),lp(2):up(2),lp(3):up(3),:,patch%it,:) = &
         self%mem(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,self%it,:)
      endif
      ! check upper boundary
      if (abs(rlb(dir)-patch%position(dir)) < hs) then
        ! We are overlapping in this direction, exchange data
        ls = 1; us=self%gn
        ls(dir) = 1 
        us(dir) = patch%li(dir) - 1
        
        lp = 1; up=patch%gn
        lp(dir) = patch%ui(dir) - patch%ng(dir) + 1
        up(dir) = patch%ui(dir) - 1
        ! copy data
        patch%mem(lp(1):up(1),lp(2):up(2),lp(3):up(3),:,patch%it,:) = &
         self%mem(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,self%it,:)
      endif
    enddo
  endif
END SUBROUTINE get_guardcells_for

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
  print '(1x,a,2(3x,3f8.3))',   'positions:', position
  do dir=1,3
    index_range(:,dir) = int(modulo(position(dir)-self%origin, self%box/self%ds)) + 1
    index_range(:,dir) = max(1, min(self%gn, index_range(:,dir)))
    print '(1x,a,2(3x,i3,a,i3))', 'index range:', &
      index_range(1,dir), ' - ', index_range(2,dir)
  end do
END FUNCTION index_range
!> =============================================================================
SUBROUTINE indices (self, np, p, ii, f)
  class(patch_t):: self
  real:: p(np), f(np,3)
  integer:: np, dir, ii(np,3), i
  !.............................................................................
  do dir=1,3
    f(:,dir) = modulo(p-self%origin(dir), self%box(dir))/self%ds(dir)
    ii(:,dir) = max(1, min(self%gn(dir), int(f(:,dir))))
  end do
  f = f-ii
  do i=1,np
    print'(i4,3x,a,3i4,3x,a,3f8.3)', i, 'index:', ii(i,:), 'fraction:', f(i,:)
  end do
END SUBROUTINE indices
!> =============================================================================
FUNCTION index (self, position)
  class(patch_t):: self
  real, dimension(3):: position, f
  integer:: index(3)
  !.............................................................................
  f = modulo(position-self%origin, self%box)/self%ds + 1.0
  index = max(1, min(self%gn-1, int(f)))
  f = f - index
  print*,'index:', index
  print*,'fraction:', f
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
  integer:: idir
  !.............................................................................
  if (minval(self%ds) <= minval(patch%ds)) then                                 ! make sure smallest guard zone is used
    overlaps = all(distance(self,patch) <= 0.5000000001_8*(self%gsize+patch%size)) ! prologation needs at least an additional zone
  else
    overlaps = all(distance(self,patch) <= 0.5000000001_8*(self%size+patch%gsize))
  endif
END FUNCTION
!> =============================================================================
!> This is the relevance criterion if we don't need corner nbors
!> =============================================================================
LOGICAL FUNCTION is_relevant_for (self, patch)
  class (patch_t):: self
  class (patch_t), pointer:: patch
  !............................................................................
  is_relevant_for = any(distance(self,patch) < 0.01*self%size)
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
