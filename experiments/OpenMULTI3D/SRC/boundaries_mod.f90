!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE patch_mod
  USE shared_mod
  implicit none
  private
  !===============================================================================
  !> Boundary flags of a patch: the faces on a non-periodic domain edge and its
  !> relation to a spherical boundary.
  !===============================================================================
  type, extends(has_bits_t), public:: boundaries_t
    integer:: id = 0
    logical:: unset1 = .true., unset2 = .true.
    class(mesh_t), pointer:: mesh(:)
    class(patch_t), pointer:: patch
    real(8):: position(3)=0d0, radius=0d0
    logical:: set_bc_bits=.true.
  contains
    procedure:: pre_init
    procedure:: init_bc_bits
    procedure:: init
    procedure:: mask
  end type
CONTAINS

  !===============================================================================
  !> Attaches the patch and sets its boundary bits from its position, before the
  !> patch meshes exist.
  !===============================================================================
  SUBROUTINE pre_init (self, patch)
    class(boundaries_t) :: self
    class(patch_t), target :: patch !< patch the boundaries belong to
    !-----------------------------------------------------------------------------
    call trace%begin ('standard_boundaries_t%init')
    self%patch => patch
    self%unset1 = .false.
    call self%init_bc_bits ()
    call trace%end ()
  END SUBROUTINE pre_init

  !===============================================================================
  !> Setting BC bits requires (only) that the patch position, origin, and box
  !> are set (at the pre_init stage we cannot assume that patch%ds is correct)
  !===============================================================================
  SUBROUTINE init_bc_bits (self)
    class(boundaries_t) :: self
    class(patch_t), pointer:: patch
    real(8)                      :: position(3), lower(3), upper(3)
    !-----------------------------------------------------------------------------
    ! Some non-cartesian setups might already set the boundary bits by some other
    ! recipe, so the bits should not be touched
    !-----------------------------------------------------------------------------
    if (self%set_bc_bits) then
      patch => self%patch
      position = patch%position
      lower = patch%origin + 0.501d0 * patch%size
      upper = patch%origin - 0.501d0 * patch%size + patch%box
      if (.not.patch%periodic(1)) then
        if (position(1)  <= lower(1)) then
          call self%set (bits%xl)
          call patch%set (bits%has_BC)
        end if
        if (position(1) >= upper(1)) then
          call self%set (bits%xu)
          call patch%set (bits%has_BC)
        end if
      end if
      if (.not.patch%periodic(2)) then
        if (position(2)  <= lower(2)) then
          call self%set (bits%yl)
          call patch%set (bits%has_BC)
        end if
        if (position(2) >= upper(2)) then
          call self%set (bits%yu)
          call patch%set (bits%has_BC)
        end if
      end if
      if (.not.patch%periodic(3)) then
        if (position(3)  <= lower(3)) then
          call self%set (bits%zl)
          call patch%set (bits%has_BC)
        end if
        if (position(3) >= upper(3)) then
          call self%set (bits%zu)
          call patch%set (bits%has_BC)
        end if
      end if
    end if
  END SUBROUTINE init_bc_bits

  !===============================================================================
  !> Initialize the boundary data type, while testing also if the patch contains
  !> a spherical boundary, or possibly is completely immersed in such a boundary.
  !===============================================================================
  SUBROUTINE init (self)
    class(boundaries_t) :: self
    !.............................................................................
    class(patch_t), pointer      :: patch
    real(8)                      :: position(3)
    real(8)                      :: size(3), p(3)
    integer                      :: i, j, k
    logical                      :: one_inside, all_inside
    !-----------------------------------------------------------------------------
    call trace%begin ('standard_boundaries_t%init')
    if (self%unset1) then
      call io%abort ('standard_boundaries_t%pre_init has not been called')
    end if
    patch => self%patch
    self%id = patch%id
    self%mesh => patch%mesh
    self%position = patch%position
    !-----------------------------------------------------------------------------
    ! If at edge, set boundary bits
    !-----------------------------------------------------------------------------
    call self%init_bc_bits
    !-----------------------------------------------------------------------------
    ! In case boundary bits are set elsewhere, make sure mesh flags are consistent
    ! (mesh flags can only be set after meshes have been allocated, while bits may
    ! be set already in the pre_init step).
    !-----------------------------------------------------------------------------
    if (associated(self%mesh)) then
      self%mesh(1)%lower_boundary = self%is_set (bits%xl)
      self%mesh(1)%upper_boundary = self%is_set (bits%xu)
      self%mesh(2)%lower_boundary = self%is_set (bits%yl)
      self%mesh(2)%upper_boundary = self%is_set (bits%yu)
      self%mesh(3)%lower_boundary = self%is_set (bits%zl)
      self%mesh(3)%upper_boundary = self%is_set (bits%zu)
    else
      call io%warning ('standard_boundaries_t%init: self%mesh not associated')
    end if
    !-----------------------------------------------------------------------------
    ! Spherical boundary flag
    !-----------------------------------------------------------------------------
    if (self%set_bc_bits) then
      position = patch%position
      one_inside = .false.
      all_inside = .true.
      do i=1,3
        size(i) = patch%size(i)
      end do
      do k=-1,1; do j=-1,1; do i=-1,1
            p(:) = position(:) + 0.5*size(:)*[i,j,k]
            if (sum(p**2) < self%radius**2) then
              one_inside = .true.
            else
              all_inside = .false.
            end if
          end do; end do; end do
      if (one_inside) call self%set (bits%spherical)
      if (all_inside) call self%set (bits%immersed)
    end if
    !-----------------------------------------------------------------------------
    self%unset2 = .false.
    call trace%end()
  END SUBROUTINE init

  !===============================================================================
  !> Mask of the interior points li:ui of the meshes m, false in the guard zones.
  !> On a spherical patch the points within the boundary radius of position are
  !> false too.
  !===============================================================================
  FUNCTION mask (self, m) RESULT (out)
    class(boundaries_t)      :: self
    class(mesh_t), pointer   :: m(:) !< meshes of the patch, one per axis
    logical                  :: out(m(1)%gn,m(2)%gn,m(3)%gn)
    !............................................................................
    integer                  :: ix, iy, iz
    real(8)                  :: x, y, z, r
    !-----------------------------------------------------------------------------
    call trace%begin('standard_boundaries_t%mask')
    out(:,:,:) = .false.
    out(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui) = .true.
    if (self%is_set (bits%spherical)) then
      associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
        do iz=m(3)%li,m(3)%ui
          z = m(3)%p + r3(iz) - self%position(3)
          do iy=m(2)%li,m(2)%ui
            y = m(2)%p + r2(iy) - self%position(2)
            do ix=m(1)%li,m(1)%ui
              x = m(1)%p + r1(ix) - self%position(1)
              r = sqrt(x**2+y**2+z**2)
              out(ix,iy,iz) = merge (out(ix,iy,iz), .false., r>self%radius)
            end do
          end do
        end do
      end associate
    end if
    call trace%end ()
  END FUNCTION mask

END MODULE boundaries_mod
