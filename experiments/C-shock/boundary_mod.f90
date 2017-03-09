!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    logical:: unset = .true.
  contains
    procedure:: init
    procedure:: set
    procedure condition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Initialize the boundary data type
!> =============================================================================
SUBROUTINE init (self, mesh, position, radius)
  class(boundary_t):: self
  class(mesh_t), pointer, optional:: mesh(:)
  real(8), optional:: position(3), radius
  !-----------------------------------------------------------------------------
  call trace%begin ('boundary_t%init')

  self%unset = .false.
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Set bits in the status word
!> =============================================================================
SUBROUTINE set (self, bits)
  class(boundary_t):: self
  integer:: bits
  !-----------------------------------------------------------------------------
  self%status = ior (self%status, bits)
END SUBROUTINE set

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, mem, m, h)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  !real, dimension(size(mem,1),size(mem,2),size(mem,3)):: infl
  !integer, dim1          :: size(mem,1)
  !integer, dim2          :: size(mem,2)
  !integer, dim3          :: size(mem,3)
  real,save,dimension(64,32,1):: inrho, ins, inpx, inpy, inpz, inbx, inby, inbz
  class(mesh_t), pointer, dimension(:) :: m
  integer,save           :: ini=1
  logical,save           :: first_time=.true.
  !............................................................................
  integer      :: ix, iy, iz, iv
  real(8)      :: x, y, z, r
  real(8)      :: h(3)

  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')
  if (first_time) then
    if (ini .eq. 1) then
      inrho = mem
    else if (ini .eq. 2) then
      ins   = mem
    else if (ini .eq. 3) then
      inpx  = mem
    else if (ini .eq. 4) then
      inpy  = mem
    else if (ini .eq. 5) then
      inpz  = mem
    else if (ini .eq. 6) then
      inbx  = mem
    else if (ini .eq. 7) then
      inby  = mem
    else if (ini .eq. 8) then
      print*,'entering first_time'
      inbz  = mem
      first_time = .false.
    end if
  end if
  !print*, 'dims of inrho, lb and ub',size(inrho,1),size(inrho,2),size(inrho,3),m(1)%lb,m(2)%lb,m(3)%lb,m(1)%ub,m(2)%ub,m(3)%ub
  do ix=m(1)%lb,m(1)%ub
   do iy=m(2)%lb,m(2)%ub
     do iz=m(3)%lb,m(3)%ub
        if (iy .gt. (m(2)%ub - 2) ) then
          mem(ix,iy,iz) = mem(ix,m(2)%ub - 2,iz)
        end if
        if (iy .lt. (m(2)%lb + 2) ) then
          mem(ix,iy,iz) = mem(ix,m(2)%lb + 2,iz)
        endif
        if (ix .gt. (m(1)%ub - 2) ) then
          mem(ix,iy,iz) = mem(m(1)%ub - 2,m(2)%ub/2,iz)
        end if
        if (ix .lt. (m(1)%lb + 2) ) then
          !mem(ix,iy,iz) = mem(m(1)%lb + 2,m(2)%ub/2,iz)
          if (ini .eq. 1) mem(ix,iy,iz) = inrho(ix,m(2)%ub/2,iz)
          if (ini .eq. 2) mem(ix,iy,iz) = ins(ix,m(2)%ub/2,iz)
          if (ini .eq. 3) mem(ix,iy,iz) = inpx(ix,m(2)%ub/2,iz)
          if (ini .eq. 4) mem(ix,iy,iz) = inpy(ix,m(2)%ub/2,iz)
          if (ini .eq. 5) mem(ix,iy,iz) = inpz(ix,m(2)%ub/2,iz)
          if (ini .eq. 6) mem(ix,iy,iz) = inbx(ix,m(2)%ub/2,iz)
          if (ini .eq. 7) mem(ix,iy,iz) = inby(ix,m(2)%ub/2,iz)
          if (ini .eq. 8) mem(ix,iy,iz) = inbz(ix,m(2)%ub/2,iz)
        end if
      end do
    end do
  end do
  if (ini .eq. 8) ini = 0
  ini = ini+1
  call trace_end
END SUBROUTINE condition

END MODULE boundary_mod
