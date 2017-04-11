!> =============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!> =============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE bits_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundaries_t
    integer:: status = 0
    logical:: unset = .true.
    class(mesh_t), pointer:: m(:)
    real(8):: position(3)=0d0, radius=0d0
    real:: vmax
  contains
    procedure:: init
    procedure:: set
    procedure:: is_set
    procedure:: is_clear
    procedure:: condition
    procedure:: symmetric
    procedure:: antisymmetric
    procedure:: extrapolate
    procedure:: spherical_extrapolate
    procedure:: spherical0
    procedure:: spherical3
    generic:: spherical => spherical0, spherical3
  end type
CONTAINS

!> =============================================================================
!> Initialize the boundary data type
!> =============================================================================
SUBROUTINE init (self, mesh, position, radius)
  class(boundaries_t):: self
  class(mesh_t), pointer, optional:: mesh(:)
  real(8), optional:: position(3), radius
  !-----------------------------------------------------------------------------
  call trace%begin ('boundaries_t%init')
  if (present(mesh))     self%m        => mesh
  if (present(position)) self%position =  position
  if (present(radius))   self%radius   =  radius
  self%unset = .false.
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Set bits in the status word
!> =============================================================================
SUBROUTINE set (self, bits)
  class(boundaries_t):: self
  integer:: bits
  !-----------------------------------------------------------------------------
  self%status = ior (self%status, bits)
END SUBROUTINE set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
FUNCTION is_set (self, bits) RESULT (out)
  class(boundaries_t):: self
  integer:: bits
  logical:: out
  !-----------------------------------------------------------------------------
  out = (iand (self%status, bits) /=0)
END FUNCTION is_set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
FUNCTION is_clear (self, bits) RESULT (out)
  class(boundaries_t):: self
  integer:: bits
  logical:: out
  !-----------------------------------------------------------------------------
  out = (iand (self%status, bits) ==0)
END FUNCTION is_clear

!> =============================================================================
!> Apply a boundary condition, using the procedure pointer "operator"
!> =============================================================================
SUBROUTINE condition (self, f, operator)
  class(boundaries_t):: self
  real, dimension(:,:,:):: f
  procedure(make_symmetric), pointer:: operator
  !.............................................................................
  class(mesh_t), pointer:: m(:), m1, m2
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('boundaries_t%condition')
  if (self%unset) then
    call mpi%abort ('boundaries_t%init has not been called')
  end if
  m => self%m
  if (self%is_set(bits%xl+bits%xu)) then
    if (io%verbose>1) print *, 'BOUNDARY CONDITION: x-direction'
    m1 => self%m(2)
    m2 => self%m(1)
    do i=m(3)%lb,m(3)%ub
      call operator (self, f(:,:,i),m1,m2,lower=self%is_set(bits%xl),upper=self%is_set(bits%xu),transpose=.true.)
    end do
  end if
  if (self%is_set(bits%yl+bits%yu)) then
    if (io%verbose>1) print *, 'BOUNDARY CONDITION: y-direction'
    m1 => self%m(1)
    m2 => self%m(2)
    do i=m(3)%lb,m(3)%ub
      call operator (self, f(:,:,i),m1,m2,lower=self%is_set(bits%yl),upper=self%is_set(bits%yu),transpose=.false.)
    end do
  end if
  if (self%is_set(bits%zl+bits%zu)) then
    if (io%verbose>1) print *, 'BOUNDARY CONDITION: z-direction'
    m1 => self%m(1)
    m2 => self%m(3)
    do i=m(2)%lb,m(2)%ub
      call operator (self, f(:,i,:),m1,m2,lower=self%is_set(bits%zl),upper=self%is_set(bits%zu),transpose=.false.)
    end do
  end if
  call trace%end
END SUBROUTINE condition

!> =============================================================================
!> Apply symmetric boundary conditions in the directions indicated by status bits
!> =============================================================================
SUBROUTINE symmetric (self, mem)
  class(boundaries_t):: self
  real, dimension(:,:,:), pointer:: mem
  procedure(make_symmetric), pointer:: operator => make_symmetric
  !-----------------------------------------------------------------------------
  call self%condition (mem, operator)
END SUBROUTINE symmetric

!> =============================================================================
!> Apply antisymmetric boundary conditions
!> =============================================================================
SUBROUTINE antisymmetric (self, mem)
  class(boundaries_t):: self
  real, dimension(:,:,:), pointer:: mem
  procedure(make_symmetric), pointer:: operator => make_antisymmetric
  !-----------------------------------------------------------------------------
  call self%condition (mem, operator)
END SUBROUTINE antisymmetric

!> =============================================================================
!> Apply extrapolating boundary conditions
!> =============================================================================
SUBROUTINE extrapolate (self, mem)
  class(boundaries_t):: self
  real, dimension(:,:,:), pointer:: mem
  procedure(make_symmetric), pointer:: operator => make_extrapolated
  !-----------------------------------------------------------------------------
  call self%condition (mem, operator)
END SUBROUTINE extrapolate

!> =============================================================================
!> Apply symmtric boundary condition to the 2nd index, looping over the 1st
!> =============================================================================
SUBROUTINE make_symmetric (self, f, m1, m2, lower, upper, transpose)
  class(boundaries_t):: self
  real, dimension(:,:):: f
  class(mesh_t), pointer:: m1, m2
  logical:: lower, upper, transpose
  integer:: i, j
  !-----------------------------------------------------------------------------
  if (lower) then
    !print *, 'lower boundary'
    if (transpose) then
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(j,i) = f(m2%li+(m2%lo-j),i)
        end do
      end do
    else
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(i,j) = f(i,m2%li+(m2%lo-j))
        end do
      end do
    end if
  end if
  if (upper) then
    !print *, 'upper boundary'
    if (transpose) then
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = f(m2%ui-(j-m2%uo),i)
        end do
      end do
    else
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = f(i,m2%ui-(j-m2%uo))
        end do
      end do
    end if
  end if
END SUBROUTINE make_symmetric

!> =============================================================================
!> Make the solution antisymmetric across the boundary
!> =============================================================================
SUBROUTINE make_antisymmetric (self, f, m1, m2, lower, upper, transpose)
  class(boundaries_t):: self
  real, dimension(:,:):: f
  class(mesh_t), pointer:: m1, m2
  logical:: lower, upper, transpose
  integer:: i, j
  !-----------------------------------------------------------------------------
  if (lower) then
    !print *, 'lower boundary'
    if (transpose) then
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(j,i) = -f(m2%li+(m2%lo-j),i)
        end do
      end do
    else
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(i,j) = -f(i,m2%li+(m2%lo-j))
        end do
      end do
    end if
  end if
  if (upper) then
    !print *, 'upper boundary'
    if (transpose) then
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = -f(m2%ui-(j-m2%uo),i)
        end do
      end do
    else
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = -f(i,m2%ui-(j-m2%uo))
        end do
      end do
    end if
  end if
END SUBROUTINE make_antisymmetric

!> =============================================================================
!> Make an extrapolated boundary condition, possibly transposing for vector
!> speedup
!> =============================================================================
SUBROUTINE make_extrapolated (self, f, m1, m2, lower, upper, transpose)
  class(boundaries_t):: self
  real, dimension(:,:):: f
  class(mesh_t), pointer:: m1, m2
  logical:: lower, upper, transpose
  integer:: i, j
  !-----------------------------------------------------------------------------
  if (lower) then
    !print *, 'lower boundary'
    if (transpose) then
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(j,i) = exp(2.*log(f(m2%li,i)) - log(f(m2%li+1+(m2%lo-j),i)))
        end do
      end do
    else
      do j=1,m2%lo
        do i=m1%lb,m1%ub
          f(i,j) = exp(2.*log(f(i,m2%li)) - log(f(i,m2%li+1+(m2%lo-j))))
        end do
      end do
    end if
  end if
  if (upper) then
    !print *, 'upper boundary'
    if (transpose) then
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(j,i) = exp(2.*log(f(m2%ui,i)) - log(f(m2%ui-1-(j-m2%uo),i)))
        end do
      end do
    else
      do j=m2%uo,m2%ub
        do i=m1%lb,m1%ub
          f(i,j) = exp(2.*log(f(i,m2%ui)) - log(f(i,m2%ui-1-(j-m2%uo))))
        end do
      end do
    end if
  end if
END SUBROUTINE make_extrapolated

!> =============================================================================
!> Compute stagger offsets, taking care with buggy compilers, which may need
!> help from specific mesh information.
!> =============================================================================
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
  class default
    call mpi%abort ('boundaries_mod::set_stagger: staggering h undefined')
  end select
END SUBROUTINE set_stagger

!> =============================================================================
!> Impose fixed boundary condition; scalar values
!> =============================================================================
SUBROUTINE spherical0 (self, mem, v, m, iv, maxfill)
  class(boundaries_t)      :: self
  real, dimension(:,:,:) :: mem
  real                   :: v
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  logical, optional      :: maxfill
  !............................................................................
  integer                :: ix, iy, iz
  real(8)                :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace_begin('boundaries_t%fixed0')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius-m(1)%d*0.1) then
          mem(ix,iy,iz) = v
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE spherical0

!> =============================================================================
!> Impose fixed boundary condition; array values
!> =============================================================================
SUBROUTINE spherical3 (self, mem, v, m, iv, maxfill)
  class(boundaries_t)      :: self
  real, dimension(:,:,:) :: mem, v
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  logical, optional      :: maxfill
  !............................................................................
  integer                :: ix, iy, iz, n1, n2
  real(8)                :: x, y, z, r, h(3)
  !-----------------------------------------------------------------------------
  call trace_begin('boundaries_t%fixed3')
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
        if (r < self%radius-m(1)%d*0.1) then
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
  call trace_end
END SUBROUTINE spherical3

!> =============================================================================
!> Extrapolate values across a spherical boundary, sampling values 4 mesh points
!> from the boundary, and at the boundary, to establish a stable slope
!> =============================================================================
SUBROUTINE spherical_extrapolate (self, mem, v, m, iv)
  class(boundaries_t)      :: self
  real, dimension(:,:,:) :: mem, v
  class(mesh_t), pointer :: m(:)
  integer                :: iv
  !............................................................................
  integer                :: ix, iy, iz, n1, n2
  real(8)                :: x, y, z, r, h(3), a1, a2, ds, p
  !-----------------------------------------------------------------------------
  call trace_begin('boundaries_t%fixed3')
  call set_stagger (m, h, iv)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  ds = m(1)%d
  a1 = 0d0; n1 = 0
  a2 = 0d0; n2 = 0
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (abs(r-(self%radius+4.*ds)) <= ds) then
          n1 = n1+1
          a1 = a1 + mem(ix,iy,iz)
        else if (abs(r-(self%radius+ds)) <= ds) then
          n2 = n2+1
          a2 = a2 + mem(ix,iy,iz)
        end if
      end do
    end do
  end do
  if (n1>0) a1 = log(a1/n1)
  if (n2>0) a2 = log(a2/n2)
  a2 = max(a2,a1)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz) + h(3)*m(3)%d - self%position(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy) + h(2)*m(2)%d - self%position(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix) + h(1)*m(1)%d - self%position(1)
        r = sqrt(x**2+y**2+z**2)
        if (r < self%radius-0.1*ds) then
          p = (self%radius+ds-r)/(3.*ds) + 1.0
          if (a1 /= 0d0 .and. a2 /= 0d0) then
            mem(ix,iy,iz) = exp(a1*(1d0-p) + a2*p)
          else
            mem(ix,iy,iz) = v(ix,iy,iz)
          end if
        end if
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE spherical_extrapolate

END MODULE boundaries_mod
