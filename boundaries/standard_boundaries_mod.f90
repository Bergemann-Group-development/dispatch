!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE standard_boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE patch_mod
  USE boundary_wavekill_mod
  USE shared_mod
  implicit none
  private
  type, public:: standard_boundaries_t
    integer:: id = 0
    integer:: bits = 0
    logical:: unset1 = .true., unset2 = .true.
    logical :: mhd = .false.
    class(mesh_t), pointer:: mesh(:)
    class(patch_t), pointer:: patch
    real(8):: position(3)=0d0, radius=0d0
    real:: vmax
    logical(kind=2), pointer:: filled(:,:,:,:)
    logical:: was_filled(3,2)=.false.
    logical:: check_filled=.false.
    logical:: set_bc_bits=.true.
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: set_attribute
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: mask
    procedure:: condition
    procedure:: symmetric
    procedure:: antisymmetric
    procedure:: extrapolate
    procedure:: extrapolate_constrained
    procedure:: log_extrapolate
    procedure:: constant
    procedure:: zero_gradient
    procedure:: spherical0
    procedure:: spherical3
    procedure:: spherical_keep
    procedure:: spherical_p_ns
    procedure:: spherical_ns
    procedure:: spherical_jagged_d
    procedure:: spherical_jagged_s0
    procedure:: spherical_jagged_s
    procedure:: spherical_jagged_p
    generic:: spherical => spherical0, spherical3
    procedure:: spherical_extrapolate
    procedure:: spherical_constant_entropy
    procedure, nopass:: sponge => wave_killing
    procedure:: set_filled
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(standard_boundaries_t) :: self
  class(patch_t), target       :: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_boundaries_t%init')
  self%patch => patch
  self%unset1 = .false.
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(standard_boundaries_t) :: self
  !.............................................................................
  class(patch_t), pointer      :: patch
  integer                      :: id
  real(8)                      :: position(3)
  real(8)                      :: size(3), p(3), limit(3)
  integer                      :: i, j, k
  logical                      :: one_inside, all_inside
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_boundaries_t%init')
  if (self%unset1) then
    call io%abort ('standard_boundaries_t%pre_init has not been called')
  end if
  patch => self%patch
  self%id = patch%id
  self%patch => patch
  self%mesh => patch%mesh
  self%position = patch%position
  !-----------------------------------------------------------------------------
  ! If at edge, set boundary bits
  ! FIXME: Consider position and origin and box definitions
  !-----------------------------------------------------------------------------
  position = patch%position
  limit = 0.5d0 * (patch%box - patch%size - patch%ds)
  ! ----------------------------------------------------------------------------
  ! some non-cartesian setups might already set the boundary bits by some other
  ! recipe, so the bits should not touched.
  ! ----------------------------------------------------------------------------
  if (self%set_bc_bits) then
    if (.not.patch%periodic(1)) then
      if (position(1)  <= -limit(1)) then
        call self%set (bits%xl)
      end if
      if (position(1) >= +limit(1)) then
        call self%set (bits%xu)
      end if
    end if
    if (.not.patch%periodic(2)) then
      if (position(2)  <= -limit(2)) then
        call self%set (bits%yl)
      end if
      if (position(2) >= +limit(2)) then
        call self%set (bits%yu)
      end if
    end if
    if (.not.patch%periodic(3)) then
      if (position(3)  <= -limit(3)) then
        call self%set (bits%zl)
      end if
      if (position(3) >= +limit(3)) then
        call self%set (bits%zu)
      end if
    end if
  end if
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
!> Set optional attributes
!===============================================================================
SUBROUTINE set_attribute (self, radius)
  class(standard_boundaries_t) :: self
  real(8), optional            :: radius
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_boundaries_t%set_attribute')
  if (present(radius)) &
    self%radius = radius
  call trace%end ()
END SUBROUTINE set_attribute

!===============================================================================
!> Set bits in the status word
!===============================================================================
SUBROUTINE set (self, bits)
  class(standard_boundaries_t) :: self
  integer                      :: bits
  !-----------------------------------------------------------------------------
  self%bits = ior (self%bits, bits)
END SUBROUTINE set

!===============================================================================
!> Set bits in the status word
!===============================================================================
SUBROUTINE clear (self, bits)
  class(standard_boundaries_t) :: self
  integer                      :: bits
  !-----------------------------------------------------------------------------
  self%bits = iand (self%bits, not(bits))
END SUBROUTINE clear

!===============================================================================
!> Check bits in the status word
!===============================================================================
FUNCTION is_set (self, bits) RESULT (out)
  class(standard_boundaries_t) :: self
  integer                      :: bits
  logical                      :: out
  !-----------------------------------------------------------------------------
  out = (iand (self%bits, bits) /=0)
END FUNCTION is_set

!===============================================================================
!> Check bits in the status word
!===============================================================================
FUNCTION is_clear (self, bits) RESULT (out)
  class(standard_boundaries_t) :: self
  integer                      :: bits
  logical                      :: out
  !-----------------------------------------------------------------------------
  out = (iand (self%bits, bits) ==0)
END FUNCTION is_clear

!===============================================================================
!> Check bits in the status word
!===============================================================================
FUNCTION mask (self, m) RESULT (out)
  class(standard_boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
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

!===============================================================================
!> Apply a boundary condition, using the procedure pointer "operator"
!===============================================================================
SUBROUTINE condition (self, f, iv, operator, value, select)
  class(standard_boundaries_t)       :: self
  real(var_k)                        :: f(:,:,:)
  integer                            :: iv
  real, optional                     :: value
  procedure(make_symmetric), pointer :: operator
  integer, optional                  :: select
  !.............................................................................
  class(mesh_t), pointer:: m(:), m1, m2
  integer:: i, saved_status
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_boundaries_t%condition')
  if (self%unset2) then
    call io%abort ('standard_boundaries_t%init has not been called')
  end if
  m => self%mesh
  if (present(select)) then
    saved_status = self%bits
    self%bits = iand (self%bits, select)
    if (self%is_set(bits%xl+bits%xu)) &
      call x_condition
    if (self%is_set(bits%yl+bits%yu)) &
      call y_condition
    if (self%is_set(bits%zl+bits%zu)) &
      call z_condition
    self%bits = saved_status
  else
    if (self%is_set(bits%xl+bits%xu)) &
      call x_condition
    if (self%is_set(bits%yl+bits%yu)) &
      call y_condition
    if (self%is_set(bits%zl+bits%zu)) &
      call z_condition
  end if
  call trace%end()
  return
contains
  subroutine x_condition
    logical:: lower, upper
    if (self%check_filled) call trace%begin ('x_condition')
    if (self%id==io%id_debug) &
      print *, 'MK BOUNDARY CONDITION: x-direction', self%bits, present(value)
    m1 => self%mesh(2)
    m2 => self%mesh(1)
    if (m2%n<2) print *, 'unconstrained boundaries:: more cells per X direction needed'
    lower = self%is_set(bits%xl)
    upper = self%is_set(bits%xu)
    do i=m(3)%lb,m(3)%ub
      call operator (self, f(:,:,i),m1,m2,iv,lower=lower, &
        upper=upper, transpose=.true., value=value)
    end do
    if (self%check_filled) then
      if (self%is_set(bits%xl)) then
        self%filled(self%mesh(1)%lb:self%mesh(1)%lo,:,:,iv) = self%patch%level
      end if
      if (self%is_set(bits%xu)) then
        self%filled(self%mesh(1)%uo:self%mesh(1)%ub,:,:,iv) = self%patch%level
      end if
    end if
    if (self%is_set(bits%xl)) self%was_filled(1,1) = .true.
    if (self%is_set(bits%xu)) self%was_filled(1,2) = .true.
    if (self%check_filled) call trace%end ()
  end subroutine x_condition
  subroutine y_condition
    if (self%check_filled) call trace%begin ('y_condition')
    if (self%id==io%id_debug) &
      print *, 'MK BOUNDARY CONDITION: y-direction', self%bits, present(value)
    m1 => self%mesh(1)
    m2 => self%mesh(2)
    if (m2%n<2) print *, 'unconstrained boundaries:: more cells per Y direction needed'
    do i=m(3)%lb,m(3)%ub
      call operator (self, f(:,:,i),m1,m2,iv,lower=self%is_set(bits%yl), &
        upper=self%is_set(bits%yu), transpose=.false., value=value)
    end do
    if (self%check_filled) then
      if (self%is_set(bits%yl)) then
        self%filled(:,self%mesh(2)%lb:self%mesh(2)%lo,:,iv) = .true.
      end if
      if (self%is_set(bits%yu)) then
        self%filled(:,self%mesh(2)%uo:self%mesh(2)%ub,:,iv) = .true.
      end if
    end if
    if (self%is_set(bits%yl)) self%was_filled(2,1) = .true.
    if (self%is_set(bits%yu)) self%was_filled(2,2) = .true.
    if (self%check_filled) call trace%end ()
  end subroutine y_condition
  subroutine z_condition
    if (self%check_filled) call trace%begin ('z_condition')
    if (self%id==io%id_debug) &
      print *, 'MK BOUNDARY CONDITION: z-direction', self%bits, present(value)
    m1 => self%mesh(1)
    m2 => self%mesh(3)
    if (m2%n<2) print *, 'unconstrained boundaries:: more cells per Z direction needed'
    do i=m(2)%lb,m(2)%ub
      call operator (self, f(:,i,:),m1,m2,iv,lower=self%is_set(bits%zl), &
        upper=self%is_set(bits%zu), transpose=.false., value=value)
    end do
    if (self%check_filled) then
      if (self%is_set(bits%zl)) then
        self%filled(:,:,self%mesh(3)%lb:self%mesh(3)%lo,iv) = .true.
      end if
      if (self%is_set(bits%zu)) then
        self%filled(:,:,self%mesh(3)%uo:self%mesh(3)%ub,iv) = .true.
      end if
    end if
    if (self%is_set(bits%zl)) self%was_filled(3,1) = .true.
    if (self%is_set(bits%zu)) self%was_filled(3,2) = .true.
    if (self%check_filled) call trace%end ()
  end subroutine z_condition
END SUBROUTINE condition

!===============================================================================
!> Apply symmetric boundary conditions in the directions indicated by status bits
!===============================================================================
SUBROUTINE symmetric (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_symmetric
  call trace%begin ('standard_boundaries_t%symmetric')
  if (self%id==io%id_debug) &
    print *, self%id, 'SYMM', iv, self%is_set(bits%yl), self%is_set(bits%yu)
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE symmetric

!===============================================================================
!> Apply antisymmetric boundary conditions
!===============================================================================
SUBROUTINE antisymmetric (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_antisymmetric
  call trace%begin ('standard_boundaries_t%antisymmetric')
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE antisymmetric

!===============================================================================
!> Apply extrapolating boundary conditions
!===============================================================================
SUBROUTINE extrapolate (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_extrapolated
  call trace%begin ('standard_boundaries_t%extrapolate')
  if (self%id==io%id_debug) &
    print *, self%id, 'EXTRAP', iv, self%is_set(bits%yl), self%is_set(bits%yu)
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE extrapolate

!===============================================================================
!> Apply extrapolating boundary conditions
!===============================================================================
SUBROUTINE extrapolate_constrained (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_extrapolated_constrained
  call trace%begin ('standard_boundaries_t%extrapolate_constrained')
  if (self%id==io%id_debug) &
    print *, self%id, 'EXTRAP', iv, self%is_set(bits%yl), self%is_set(bits%yu)
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE extrapolate_constrained

!===============================================================================
!> Apply extrapolating boundary conditions
!===============================================================================
SUBROUTINE log_extrapolate (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_log_extrapolated
  call trace%begin ('standard_boundaries_t%extrapolate')
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE log_extrapolate

!===============================================================================
!> Apply constant boundary conditions
!===============================================================================
SUBROUTINE constant (self, mem, iv, value, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  real                               :: value
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_constant
  call trace%begin ('standard_boundaries_t%constant')
  call self%condition (mem, iv, operator, value=value, select=select)
  call trace%end()
END SUBROUTINE constant

!===============================================================================
!> Apply constant/zero gradient boundary conditions
!===============================================================================
SUBROUTINE zero_gradient (self, mem, iv, select)
  class(standard_boundaries_t)       :: self
  real(var_k),               pointer :: mem(:,:,:)
  integer                            :: iv
  integer, optional                  :: select
  procedure(make_symmetric), pointer :: operator
  !-----------------------------------------------------------------------------
  operator => make_zero_grad
  call trace%begin ('standard_boundaries_t%constant')
  call self%condition (mem, iv, operator, select=select)
  call trace%end()
END SUBROUTINE zero_gradient

!===============================================================================
!> Set indices that define the location of the physical boundary, by determining
!> how symmtry conditions are applied across boundaries.  lo and uo are the points
!> in the boundary zones closest to the interior, while li and ui are the neartest
!> points in the interior, which are being used to construct the boundary zone
!> values.
!===============================================================================
SUBROUTINE bndry_indices (m, iv, lo, li, ui, uo)
  class(mesh_t), pointer :: m
  integer                :: iv, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  ! no-mans-land mesh, with staggered points at patch boundaries
  !                             |                     |
  !  centered variable:   x---x-|-o---o-- ... --o---o-|-x---x
  !                          lo | li               ui | uo
  ! staggered variable: x---x---o---o-- ... --o---o---x---x
  !                            lo  li            ui  uo
  !-----------------------------------------------------------------------------
  ! --- down-staggered variables ---
  if (m%h(iv) < 0.0) then
    lo = m%li
    li = lo+1
    uo = m%uo
    ui = uo-1
  ! --- up-staggered variables ---
  else if (m%h(iv) > 0.0) then
    lo = m%lo
    li = lo+1
    uo = m%ui
    ui = uo-1
  ! --- centered variables ---
  else
    lo = m%lo
    li = lo+1
    uo = m%uo
    ui = uo-1
  end if
END SUBROUTINE bndry_indices

!===============================================================================
!> Apply symmtric boundary condition to the 2nd index, looping over the 1st
!===============================================================================
SUBROUTINE make_symmetric (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = f(li+(lo-j),i)
        end do
      end do
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = f(i,li+(lo-j))
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = f(ui-(j-uo),i)
        end do
      end do
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = f(i,ui-(j-uo))
        end do
      end do
    end if
  end if
END SUBROUTINE make_symmetric

!===============================================================================
!> Make the solution antisymmetric across the boundary
!===============================================================================
SUBROUTINE make_antisymmetric (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = -f(li+(lo-j),i)
        end do
      end do
      if (li-lo == 2) then
        do i=m1%lb,m1%ub
          f(li-1,i) = 0.0
        end do
      end if
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = -f(i,li+(lo-j))
        end do
      end do
      if (li-lo == 2) then
        do i=m1%lb,m1%ub
          f(i,li-1) = 0.0
        end do
      end if
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = -f(ui-(j-uo),i)
        end do
      end do
      if (uo-ui == 2) then
        do i=m1%lb,m1%ub
          f(ui+1,i) = 0.0
        end do
      end if
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = -f(i,ui-(j-uo))
        end do
      end do
      if (uo-ui == 2) then
        do i=m1%lb,m1%ub
          f(i,ui+1) = 0.0
        end do
      end if
    end if
  end if
END SUBROUTINE make_antisymmetric

!===============================================================================
!> Make an extrapolated boundary condition, possibly transposing for vector
!> speedup
!===============================================================================
SUBROUTINE make_extrapolated (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = 2.*f(lo+1,i) - f(lo+2+(lo-j),i)
        end do
      end do
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = 2.*f(i,lo+1) - f(i,lo+2+(lo-j))
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = 2.*f(uo-1,i) - f(uo-2-(j-uo),i)
        end do
      end do
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = 2.*f(i,uo-1) - f(i,uo-2-(j-uo))
        end do
      end do
    end if
  end if
END SUBROUTINE make_extrapolated

!===============================================================================
!> Make an extrapolated boundary condition, possibly transposing for vector
!> speedup
!===============================================================================
SUBROUTINE make_extrapolated_constrained (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  real                         :: p
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        p = 1.-((lo+1-j)/(lo+1.-3.))**2
        do i=m1%lb,m1%ub
          f(j,i) = p*(f(lo+1,i)-f(lo+2+(lo-j),i)) + f(lo+1,i)
        end do
      end do
    else
      do j=1,lo
        p = 1.-((lo+1-j)/(lo+1.-3.))**2
        do i=m1%lb,m1%ub
          f(i,j) = p*(f(i,lo+1) - f(i,lo+2+(lo-j))) + f(i,lo+1)
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        p = 1.-((uo-1.-j)/real(m2%ub-uo-1.-2.))**2
        do i=m1%lb,m1%ub
          f(j,i) = p*(f(uo-1,i) - f(uo-2-(j-uo),i)) + f(uo-1,i)
        end do
      end do
    else
      do j=uo,m2%ub
        p = 1.-((uo-1.-j)/real(m2%ub-uo-1.-2.))**2
        do i=m1%lb,m1%ub
          f(i,j) = p*(f(i,uo-1) - f(i,uo-2-(j-uo))) + f(i,uo-1)
        end do
      end do
    end if
  end if
END SUBROUTINE make_extrapolated_constrained

!===============================================================================
!> Make an extrapolated boundary condition, possibly transposing for vector
!> speedup
!===============================================================================
SUBROUTINE make_log_extrapolated (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                                  :: iv
  logical                                  :: lower, upper, transpose
  real, optional                           :: value
  integer                                  :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = exp(2.*log(f(lo+1,i)) - log(f(lo+2+(lo-j),i)))
        end do
      end do
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = exp(2.*log(f(i,lo+1)) - log(f(i,lo+2+(lo-j))))
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = exp(2.*log(f(uo-1,i)) - log(f(uo-2-(j-uo),i)))
        end do
      end do
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = exp(2.*log(f(i,uo-1)) - log(f(i,uo-2-(j-uo))))
        end do
      end do
    end if
  end if
END SUBROUTINE make_log_extrapolated

!===============================================================================
!> Make an constant boundary condition, possibly transposing for vector
!> speedup
!===============================================================================
SUBROUTINE make_constant (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = value
        end do
      end do
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = value
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = value
        end do
      end do
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = value
        end do
      end do
    end if
  end if
END SUBROUTINE make_constant

!===============================================================================
!> Make an constant boundary condition, possibly transposing for vector
!> speedup
!===============================================================================
SUBROUTINE make_zero_grad (self, f, m1, m2, iv, lower, upper, transpose, value)
  class(standard_boundaries_t) :: self
  real(var_k), dimension(:,:)  :: f
  class(mesh_t), pointer       :: m1, m2
  integer                      :: iv
  logical                      :: lower, upper, transpose
  real, optional               :: value
  integer                      :: i, j, lo, li, ui, uo
  !-----------------------------------------------------------------------------
  call bndry_indices (m2, iv, lo, li, ui, uo)
  if (lower) then
    if (transpose) then
      do j=1,lo
        do i=m1%lb,m1%ub
          f(j,i) = f(li,i)
        end do
      end do
    else
      do j=1,lo
        do i=m1%lb,m1%ub
          f(i,j) = f(i,li)
        end do
      end do
    end if
  end if
  if (upper) then
    if (transpose) then
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = f(ui,i)
        end do
      end do
    else
      do j=uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = f(i,ui)
        end do
      end do
    end if
  end if
END SUBROUTINE make_zero_grad

!===============================================================================
!> Compute stagger offsets, taking care with buggy compilers, which may need
!> help from specific mesh information.
!===============================================================================
SUBROUTINE set_stagger (m, h, iv)
  class(mesh_t), pointer :: m(:)
  real(8)                :: h(3)
  integer                :: iv, i
  !-----------------------------------------------------------------------------
  select type (m)
  type is (Cartesian_mesh)
    do i=1,3
      h(i) = m(i)%h(iv)
    end do
  type is (cylindrical_mesh)
    do i=1,3
      h(i) = m(i)%h(iv)
    end do
  type is (spherical_mesh)
    do i=1,3
      h(i) = m(i)%h(iv)
    end do
  class default
    call io%abort ('boundaries_mod::set_stagger: staggering h undefined')
  end select
END SUBROUTINE set_stagger

!===============================================================================
!> Impose fixed boundary condition; scalar values
!===============================================================================
SUBROUTINE spherical0 (self, mem, v, m, iv, maxfill)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: mem
  real                          :: v
  class(mesh_t), pointer        :: m(:)
  integer                       :: iv
  logical, optional             :: maxfill
  !............................................................................
  integer                                    :: ix, iy, iz
  real(8)                                    :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical0')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius+m(1)%d*0.5) then
          mem(ix,iy,iz) = v
        end if
      end do
    end do
  end do
  end associate
  call trace%end ()
END SUBROUTINE spherical0

!===============================================================================
!> Copy the values inside a given radius from the previous time slice
!===============================================================================
SUBROUTINE spherical_keep (self, mem, iv, it, new, m, position, radius)
  class(standard_boundaries_t) :: self
  real(var_k), pointer         :: mem(:,:,:,:,:)
  integer                      :: iv, it, new
  class(mesh_t), pointer       :: m(:)
  real(8)                      :: position(3), radius
  !............................................................................
  integer                      :: ix, iy, iz, i
  real                         :: h(3), r(3), rad2
  class(mesh_t), pointer       :: mp
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_p_ns')
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do i=1,3
    mp => m(i)
    h(i) = mp%h(iv)*mp%d
  end do
  rad2 = (radius-m(1)%d)**2
  do iz=m(3)%lb,m(3)%ub
    r(3) = m(3)%p - position(3) + r3(iz) + h(3)
    do iy=m(2)%lb,m(2)%ub
      r(2) = m(2)%p - position(2) + r2(iy) + h(2)
      do ix=m(1)%lb,m(1)%ub
        r(1) = m(1)%p - position(1) + r1(ix) + h(1)
        if (sum(r**2) < rad2) then
          mem(ix,iy,iz,iv,new) = mem(ix,iy,iz,iv,it)
        end if
      end do
    end do
  end do
  end associate
  call trace%end ()
END SUBROUTINE spherical_keep

!===============================================================================
!> Impose vanishing mass flux on any face that is inside the radius; this
!> automatically makes density and entropy well behaved there also
!===============================================================================
SUBROUTINE spherical_p_ns (self, p, m, position, radius)
  class(standard_boundaries_t) :: self
  real(var_k), pointer         :: p(:,:,:,:)
  class(mesh_t), pointer       :: m(:)
  real(8)                      :: position(3), radius
  !............................................................................
  integer                      :: ix, iy, iz, iv, i
  real                         :: h(3), r(3)
  class(mesh_t), pointer       :: mm
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_p_ns')
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iv=1,3
    do i=1,3
      mm => m(i)
      h(i) = mm%h(idx%px+iv-1)*mm%d
    end do
    do iz=m(3)%lb,m(3)%ub
      r(3) = m(3)%p - position(3) + r3(iz) + h(3)
      do iy=m(2)%lb,m(2)%ub
        r(2) = m(2)%p - position(2) + r2(iy) + h(2)
        do ix=m(1)%lb,m(1)%ub
          r(1) = m(1)%p - position(1) + r1(ix) + h(1)
          if (sum(r**2) < radius**2) p(ix,iy,iz,iv) = 0.0
        end do
      end do
    end do
  end do
  end associate
  call trace%end ()
END SUBROUTINE spherical_p_ns

!===============================================================================
!> Impose no-slip boundary condition; scalar values
!===============================================================================
SUBROUTINE spherical_ns (self, mem, v, m, iv, maxfill)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: mem
  real                          :: v
  class(mesh_t), pointer        :: m(:)
  integer                       :: iv
  logical, optional             :: maxfill
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x, y, z, r, h(3)
  real(8)      :: d1, dw, x1, rp, p1, p2
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_ns')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (iv == idx%px) then                                          ! mass flux in X-direction
          if (r <= self%radius) then
            p1 = sqrt(y**2+z**2+(m(1)%p + r1(ix-1) + h(1)*m(1)%d - self%position(1))**2)
            p2 = sqrt(y**2+z**2+(m(1)%p + r1(ix+1) + h(1)*m(1)%d - self%position(1))**2)
            if ((p1 >= self%radius+m(1)%d*0.5).and.(p2 <= self%radius+m(1)%d*0.5)) then
              !entering the boundary
              rp = sqrt(abs(self%radius**2 - y**2 - z**2))              ! boundary location at X axis
              if (x<0) rp = -rp
              x1 = m(1)%p + r1(ix-1) + h(1)*m(1)%d - self%position(1)   ! last point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(1)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix-1,iy,iz) = mem(ix-2,iy,iz) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix  ,iy,iz) = &
                          extrapolate_antisym(x,x1,rp,mem(ix-1,iy,iz))  ! antisymmetric extrapolation from the last point before the bnd
            else if ((p1 <= self%radius+m(1)%d*0.5).and.(p2 >= self%radius+m(1)%d*0.5)) then
            ! leaving the boundary
              rp = sqrt(abs(self%radius**2 - y**2 - z**2))              ! boundary location at X axis
              if (x<0) rp = -rp
              x1 = m(1)%p + r1(ix+1) + h(1)*m(1)%d - self%position(1)   ! first point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and the point outside it
              dw = d1/(d1+m(1)%d)                                       ! interpolation weight for the first point after the bnd
              mem(ix+1,iy,iz) = mem(ix+2,iy,iz)*dw                      ! linear interpolation for the first point after the bnd
              mem(ix  ,iy,iz) = &
                          extrapolate_antisym(x,x1,rp,mem(ix+1,iy,iz))    ! antisymmetric extrapolation from the first point after the bnd
            else if ((p1 <= self%radius+m(1)%d*0.5).and.(p2 <= self%radius+m(1)%d*0.5)) then
              mem(ix,iy,iz) = 0.0                                       ! cell is well inside the boundary
            else
              rp = sqrt(abs(self%radius**2 - y**2 - z**2))              ! boundary location at X axis
              if (x<0) rp = -rp
              x1 = m(1)%p + r1(ix-1) + h(1)*m(1)%d - self%position(1)   ! last point INSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(1)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix-1,iy,iz) = mem(ix-2,iy,iz) * dw                    ! linear interpolation for the last point before the bnd
              x1 = m(1)%p + r1(ix+1) + h(1)*m(1)%d - self%position(1)   ! first point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(1)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix+1,iy,iz) = mem(ix+2,iy,iz) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix,iy,iz) = 0.0                                       ! zero for the bnd cell
            end if
          end if
        else if (iv == idx%py) then                                     ! mass flux in Y-direction
          if (r <= self%radius) then
            p1 = sqrt(x**2+(m(2)%p + r2(iy-1) + h(2)*m(2)%d - self%position(2))**2+z**2)
            p2 = sqrt(x**2+(m(2)%p + r2(iy+1) + h(2)*m(2)%d - self%position(2))**2+z**2)
            if ((p1 > self%radius+m(2)%d*0.5).and.(p2 < self%radius+m(2)%d*0.5)) then
              !entering the boundary
              rp = sqrt(abs(self%radius**2 - x**2 - z**2))              ! boundary location at Y axis
              if (y<0) rp = -rp
              x1 = m(2)%p + r2(iy-1) + h(2)*m(2)%d - self%position(2)   ! last point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(2)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy-1,iz) = mem(ix,iy-2,iz) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix,iy  ,iz) = &
                          extrapolate_antisym(y,x1,rp,mem(ix,iy-1,iz))  ! antisymmetric extrapolation from the last point before the bnd
            else if ((p1 < self%radius+m(2)%d*0.5).and.(p2 > self%radius+m(2)%d*0.5)) then
            ! leaving the boundary
              rp = sqrt(abs(self%radius**2 - x**2 - z**2))              ! boundary location at Y axis
              if (y<0) rp = -rp
              x1 = m(2)%p + r2(iy+1) + h(2)*m(2)%d - self%position(2)   ! first point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and the point outside it
              dw = d1/(d1+m(2)%d)                                       ! interpolation weight for the first point after the bnd
              mem(ix,iy+1,iz) = mem(ix,iy+2,iz)*dw                      ! linear interpolation for the first point after the bnd
              mem(ix,iy,iz) = &
                          extrapolate_antisym(y,x1,rp,mem(ix,iy+1,iz))  ! antisymmetric extrapolation from the first point after the bnd
            else if ((p1 < self%radius+m(2)%d*0.5).and.(p2 < self%radius+m(2)%d*0.5)) then
              mem(ix,iy,iz) = 0.0                                       ! cell is well inside the boundary
            else                                                        ! special case, where a single cell gets a bc
              rp = sqrt(abs(self%radius**2 - x**2 - z**2))              ! boundary location at Y axis
              if (y<0) rp = -rp
              x1 = m(2)%p + r2(iy-1) + h(2)*m(2)%d - self%position(2)   ! last point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(2)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy-1,iz) = mem(ix,iy-2,iz) * dw                    ! linear interpolation for the last point before the bnd
              x1 = m(2)%p + r2(iy+1) + h(2)*m(2)%d - self%position(2)   ! first point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(2)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy+1,iz) = mem(ix,iy+2,iz) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix,iy,iz) = 0.0                                       ! zero for the bnd cell
            end if
          end if
        else if (iv == idx%pz) then                                     ! mass flux in Z-direction
          if (r <= self%radius) then
            p1 = sqrt(x**2+y**2+(m(3)%p + r3(iz-1) + h(3)*m(3)%d - self%position(3))**2)
            p2 = sqrt(x**2+y**2+(m(3)%p + r3(iz+1) + h(3)*m(3)%d - self%position(3))**2)
            if ((p1 > self%radius+m(3)%d*0.5).and.(p2 < self%radius+m(3)%d*0.5)) then
              !entering the boundary
              rp = sqrt(abs(self%radius**2 - x**2 - y**2))              ! boundary location at Z axis
              if (z<0) rp = -rp
              x1 = m(3)%p + r3(iz-1) + h(3)*m(3)%d - self%position(3)   ! last point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(3)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy,iz-1) = mem(ix,iy,iz-2) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix,iy,iz  ) = &
                          extrapolate_antisym(z,x1,rp,mem(ix,iy,iz-1))  ! antisymmetric extrapolation from the last point before the bnd
            else if ((p1 < self%radius+m(3)%d*0.5).and.(p2 > self%radius+m(3)%d*0.5)) then
            ! leaving the boundary
              rp = sqrt(abs(self%radius**2 - x**2 - y**2))              ! boundary location at Z axis
              if (y<0) rp = -rp
              x1 = m(3)%p + r3(iz+1) + h(3)*m(3)%d - self%position(3)   ! last point INSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and the point outside it
              dw = d1/(d1+m(3)%d)                                       ! interpolation weight for the first point after the bnd
              mem(ix,iy,iz+1) = mem(ix,iy,iz+2)*dw                      ! linear interpolation for the first point after the bnd
              mem(ix,iy,iz) = &
                          extrapolate_antisym(z,x1,rp,mem(ix,iy,iz+1))    ! antisymmetric extrapolation from the first point after the bnd
            else if ((p1 < self%radius+m(3)%d*0.5).and.(p2 < self%radius+m(3)%d*0.5)) then
              mem(ix,iy,iz) = 0.0                                       ! cell is well inside the boundary
            else
              rp = sqrt(abs(self%radius**2 - x**2 - y**2))              ! boundary location at Z axis
              if (z<0) rp = -rp
              x1 = m(3)%p + r3(iz-1) + h(3)*m(3)%d - self%position(3)   ! last point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(3)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy,iz-1) = mem(ix,iy,iz-2) * dw                    ! linear interpolation for the last point before the bnd
              x1 = m(3)%p + r3(iz+1) + h(3)*m(3)%d - self%position(3)   ! first point OUTSIDE the bnd
              d1 = abs(rp-x1)                                           ! distance between bnd and last point before it
              dw = d1/(d1+m(3)%d)                                       ! interpolation weight for the last point before the bnd
              mem(ix,iy,iz+1) = mem(ix,iy,iz+2) * dw                    ! linear interpolation for the last point before the bnd
              mem(ix,iy,iz) = 0.0                                       ! zero for the bnd cell
            end if
          end if
        end if
      end do
    end do
  end do
 ! error stop "end"
  end associate
  call trace%end ()
END SUBROUTINE spherical_ns

!===============================================================================
!> Impose spherical density boundary condition, assuming that the sphere is jagged.
!===============================================================================
SUBROUTINE spherical_jagged_d (self, mem, m, iv, it, new)
  class(standard_boundaries_t)    :: self
  real(var_k), dimension(:,:,:,:) :: mem
  class(mesh_t), pointer          :: m(:), m1, m2, m3
  integer                         :: iv, it, new
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r(3)
  logical      :: inside(8)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_jagged_d')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
  do iy=m(2)%lb,m(2)%ub
  do ix=m(1)%lb,m(1)%ub
    r(1) = m(1)%p - self%position(1) + r1(ix)
    r(2) = m(2)%p - self%position(2) + r2(iy)
    r(3) = m(3)%p - self%position(3) + r3(iz)
    if (sum(r**2) < self%radius**2) mem(ix,iy,iz,new) = mem(ix,iy,iz,it)
  end do
  end do
  end do
  end associate
  nullify (m1, m2, m3)
  call trace%end ()
END SUBROUTINE spherical_jagged_d

!===============================================================================
!> Impose spherical boundary condition, assuming that the sphere is jagged
!===============================================================================
SUBROUTINE spherical_jagged_s (self, new, old, m)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: old, new
  class(mesh_t), pointer        :: m(:), m1, m2, m3
  !............................................................................
  type(index_t):: idx
  integer      :: ix, iy, iz
  real(8)      :: r(3)
  logical      :: inside(8)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_jagged_s')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
  do iy=m(2)%lb,m(2)%ub
  do ix=m(1)%lb,m(1)%ub
    r(1) = m(1)%p - self%position(1) + r1(ix)
    r(2) = m(2)%p - self%position(2) + r2(iy)
    r(3) = m(3)%p - self%position(3) + r3(iz)
    if (sum(r**2) < self%radius**2) new(ix,iy,iz) = old(ix,iy,iz)
  end do
  end do
  end do
  end associate
  nullify (m1, m2, m3)
  call trace%end ()
END SUBROUTINE spherical_jagged_s

!===============================================================================
!> Impose spherical density boundary condition, assuming that the sphere is jagged.
!===============================================================================
SUBROUTINE spherical_jagged_s0 (self, mem, m, iv, v)
  class(standard_boundaries_t)                        :: self
  real(var_k), dimension(:,:,:) :: mem
  class(mesh_t), pointer        :: m(:), m1, m2, m3
  integer                       :: iv, it
  real                          :: v
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r2
  logical      :: inside(8)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_jagged_s0')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  r2 = self%radius**2
  if (iv == idx%s) then
    do iz = m3%lb, m3%ub
      do iy = m2%lb, m2%ub
        do ix = m1%lb, m1%ub
          x1 = m1%p - self%position(1) + m1%r(ix) + 0.5*m1%d
          x2 = m1%p - self%position(1) + m1%r(ix) - 0.5*m1%d
          y1 = m2%p - self%position(2) + m2%r(iy) + 0.5*m2%d
          y2 = m2%p - self%position(2) + m2%r(iy) - 0.5*m2%d
          z1 = m3%p - self%position(3) + m3%r(iz) + 0.5*m3%d
          z2 = m3%p - self%position(3) + m3%r(iz) - 0.5*m3%d
          inside(1) = x1**2+y1**2+z1**2 < r2
          inside(2) = x2**2+y1**2+z1**2 < r2
          inside(3) = x1**2+y2**2+z1**2 < r2
          inside(4) = x2**2+y2**2+z1**2 < r2
          inside(5) = x1**2+y1**2+z2**2 < r2
          inside(6) = x2**2+y1**2+z2**2 < r2
          inside(7) = x1**2+y2**2+z2**2 < r2
          inside(8) = x2**2+y2**2+z2**2 < r2
          if (all(inside)) mem(ix,iy,iz) = v
        end do
      end do
    end do
  else
    print *, 'iv passed',iv, 'iv required', idx%s
    call io%abort ("boundaries_mod::spherical_jagged_d:: invalid memory index iv")
  end if
  nullify (m1, m2, m3)
  call trace%end ()
END SUBROUTINE spherical_jagged_s0

!===============================================================================
!> Same as above, but for momentum
!===============================================================================
SUBROUTINE spherical_jagged_p (self, mem, m, iv)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: mem
  class(mesh_t), pointer        :: m(:), m1, m2, m3
  integer                       :: iv
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r2, h(3)
  logical      :: inside(4)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_jagged_p')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  call set_stagger (m, h, iv)
  r2 = self%radius**2
  if (iv == idx%px) then
    do iz = m3%lb, m3%ub
      do iy = m2%lb, m2%ub
        do ix = m1%lb, m1%ub
          y1 = m2%p - self%position(2) + m2%r(iy) + 0.5*m2%d
          y2 = m2%p - self%position(2) + m2%r(iy) - 0.5*m2%d
          z1 = m3%p - self%position(3) + m3%r(iz) + 0.5*m3%d
          z2 = m3%p - self%position(3) + m3%r(iz) - 0.5*m3%d
          x1 = m1%p - self%position(1) + m1%r(ix) + h(1)*m1%d
          inside(1) = y1**2+z1**2 + x1**2 < r2
          inside(2) = y2**2+z1**2 + x1**2 < r2
          inside(3) = y1**2+z2**2 + x1**2 < r2
          inside(4) = y2**2+z2**2 + x1**2 < r2
          if (any(inside)) mem(ix,iy,iz) = 0.0
        end do
      end do
    end do
  else if (iv == idx%py) then
    do iz = m3%lb, m3%ub
      do iy = m2%lb, m2%ub
        do ix = m1%lb, m1%ub
          x1 = m1%p - self%position(1) + m1%r(ix) + 0.5*m1%d
          x2 = m1%p - self%position(1) + m1%r(ix) - 0.5*m1%d
          z1 = m3%p - self%position(3) + m3%r(iz) + 0.5*m3%d
          z2 = m3%p - self%position(3) + m3%r(iz) - 0.5*m3%d
          y1 = m2%p - self%position(2) + m2%r(iy) + h(2)*m2%d
          inside(1) = x1**2+z1**2 + y1**2 < r2
          inside(2) = x2**2+z1**2 + y1**2 < r2
          inside(3) = x1**2+z2**2 + y1**2 < r2
          inside(4) = x2**2+z2**2 + y1**2 < r2
          if (any(inside)) mem(ix,iy,iz) = 0.0
        end do
      end do
    end do
  else if (iv == idx%pz) then
    do iz = m3%lb, m3%ub
      do iy = m2%lb, m2%ub
        do ix = m1%lb, m1%ub
          x1 = m1%p - self%position(1) + m1%r(ix) + 0.5*m1%d
          x2 = m1%p - self%position(1) + m1%r(ix) - 0.5*m1%d
          y1 = m2%p - self%position(2) + m2%r(iy) + 0.5*m2%d
          y2 = m2%p - self%position(2) + m2%r(iy) - 0.5*m2%d
          z1 = m3%p - self%position(3) + m3%r(iz) + h(3)*m3%d
          inside(1) = x1**2+y1**2 + z1**2 < r2
          inside(2) = x2**2+y1**2 + z1**2 < r2
          inside(3) = x1**2+y2**2 + z1**2 < r2
          inside(4) = x2**2+y2**2 + z1**2 < r2
          if (any(inside)) mem(ix,iy,iz) = 0.0
        end do
      end do
    end do
  else
    print *, 'iv passed',iv, 'iv required', idx%px, idx%py, idx%pz
    call io%abort ("boundaries_mod::spherical_jagged_p:: invalid memory index iv")
  end if
  nullify (m1, m2, m3)
  call trace%end ()
END SUBROUTINE spherical_jagged_p

!=======================================================================
!> Linear antisymmetric extrapolation
!=======================================================================
FUNCTION extrapolate_antisym(ip, op, bp, yop) result(out)
  real(8) :: ip, op, bp
  real(var_k):: yop, out
  !.....................................................................
  out = yop+(ip-op)/(bp-op)*(-yop)
END FUNCTION extrapolate_antisym

!===============================================================================
!> Impose fixed boundary condition; array values
!===============================================================================
SUBROUTINE spherical3 (self, mem, v, m, iv, maxfill)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: mem, v
  class(mesh_t), pointer        :: m(:)
  integer                       :: iv
  logical, optional             :: maxfill
  !............................................................................
  integer                                    :: ix, iy, iz, n1, n2
  real(8)                                    :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical3')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  self%vmax = 0.0
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius+m(1)%d*0.5) then
          mem(ix,iy,iz) = v(ix,iy,iz)
        else if (present(maxfill)) then
          self%vmax = max(self%vmax,mem(ix,iy,iz))
        end if
      end do
    end do
  end do
  if (present(maxfill)) then
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius-m(1)%d*0.1) then
          mem(ix,iy,iz) = max(mem(ix,iy,iz),self%vmax)
        end if
      end do
    end do
  end do
  end if
  end associate
  call trace%end ()
END SUBROUTINE spherical3

!===============================================================================
!> Extrapolate values across a spherical boundary, sampling values 4 mesh points
!> from the boundary, and at the boundary, to establish a stable slope
!===============================================================================
SUBROUTINE spherical_extrapolate (self, mem, m, iv)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: mem
  class(mesh_t), pointer        :: m(:)
  integer                       :: iv
  !............................................................................
  integer                       :: ix, iy, iz, n1, n2
  real(8)                       :: x, y, z, r, h(3), a1, a2, ds, p
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%spherical_extrapolate')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  ds = m(1)%d
  a1 = 0d0; n1 = 0
  a2 = 0d0; n2 = 0
  !-----------------------------------------------------------------------------
  ! a1 and a2 sample the values in annuli around R-3*ds and R, respectively
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (abs(r-(self%radius+2.*ds)) <= ds) then
          n1 = n1+1
          a1 = a1 + mem(ix,iy,iz)
        else if (abs(r-self%radius) <= ds) then
          n2 = n2+1
          a2 = a2 + mem(ix,iy,iz)
        end if
      end do
    end do
  end do
  if (n1>0) a1 = log(a1/n1)
  if (n2>0) a2 = log(a2/n2)
  a2 = max(a2,a1)
  !-----------------------------------------------------------------------------
  ! For radii inside R-0.5*ds, set the value to be log-interpolated from a1, a2
  ! The value at a2 would be chosen for r=R
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        r = max(r,self%radius-4.0*ds)
        if (r < self%radius) then
          if (a1 /= 0d0 .and. a2 /= 0d0) then
            p = (self%radius-r)/(2.0*ds) + 1.0
            mem(ix,iy,iz) = exp(a1*(1d0-p) + a2*p)
          end if
        end if
      end do
    end do
  end do
  end associate
  call trace%end ()
END SUBROUTINE spherical_extrapolate

!===============================================================================
!> Set constant entropy inside the boundary radius (included)
!===============================================================================
SUBROUTINE spherical_constant_entropy (self, d, e, m, s, gamma)
  class(standard_boundaries_t)  :: self
  real(var_k), dimension(:,:,:) :: d, e
  class(mesh_t), pointer        :: m(:)
  real                          :: s, gamma
  !............................................................................
  integer                       :: ix, iy, iz
  real(8)                       :: x, y, z, r, ds, c
  !-----------------------------------------------------------------------------
  call trace%begin('standard_boundaries_t%sperical_entropy')
  !-----------------------------------------------------------------------------
  ! a1 and a2 sample the values in annuli around R-3*ds and R, respectively
  !-----------------------------------------------------------------------------
  c = exp(s*(gamma-1.0))/(gamma-1.0)
  ds = m(1)%d
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix)- self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius) then
          e(ix,iy,iz) = c*d(ix,iy,iz)**gamma
          if (iy==18.and.iz==18) print *,'BC: s', m(1)%id, ix, s
        end if
      end do
    end do
  end do
  end associate
  call trace%end ()
END SUBROUTINE spherical_constant_entropy

!===============================================================================
!> Set self%patch%filled(:,:,:) if associated
!===============================================================================
SUBROUTINE set_filled (self)
  class(standard_boundaries_t):: self
  integer:: i, j, l(3), u(3)
  logical:: boundary(2)
  character(5):: dir
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_boudaries_t%set_filled')
  if (associated(self%patch%filled)) then
    do i=1,3
      boundary(1) = self%mesh(i)%lower_boundary
      boundary(2) = self%mesh(i)%upper_boundary
      do j=1,2
        if (boundary(j)) then
          l = self%patch%mesh%lb
          u = self%patch%mesh%ub
          if (j==1) then
            dir = 'lower'
            u(i) = self%patch%mesh(i)%lo
          else
            dir = 'upper'
            l(i) = self%patch%mesh(1)%uo
          end if
          self%patch%filled(l(1):u(1),l(2):u(2),l(3):u(3)) = self%patch%level
          if (.not.self%was_filled(i,j)) then
            write (stderr,'(a,i2,a,i6)') &
              'WARNING: unfilled '//dir//' boundary for direction', i, &
              ', patch', self%patch%id
            self%was_filled(i,j) = .true.
          end if
        end if
      end do
    end do
  end if
  call trace%end ()
END SUBROUTINE set_filled

END MODULE standard_boundaries_mod
