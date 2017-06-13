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
SUBROUTINE condition (self, mem, m, lower)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  real, save, dimension(:,:,:), pointer:: inrho, ins, inpx, inpy, inpz, inbx, inby, inbz
  !integer, dim1          :: size(mem,1)
  !integer, dim2          :: size(mem,2)
  !integer, dim3          :: size(mem,3)
  !real,save,dimension(64,32,1):: inrho, ins, inpx, inpy, inpz, inbx, inby, inbz
  class(mesh_t), pointer, dimension(:) :: m
  integer,save           :: ini=1
  logical,save           :: first_time=.true.
  !............................................................................
  integer      :: ix, iy, iz, iv
  real(8)      :: x, y, z, r
  logical, optional:: lower

  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%condition')
  if (first_time) then
    if (ini .eq. 1) then
      allocate (inrho(size(mem,1),size(mem,2),size(mem,3)))
      inrho = mem
    else if (ini .eq. 2) then
      allocate (ins(size(mem,1),size(mem,2),size(mem,3)))
      ins   = mem
    else if (ini .eq. 3) then
      allocate (inpx(size(mem,1),size(mem,2),size(mem,3)))
      inpx  = mem
    else if (ini .eq. 4) then
      allocate (inpy(size(mem,1),size(mem,2),size(mem,3)))
      inpy  = mem
    else if (ini .eq. 5) then
      allocate (inpz(size(mem,1),size(mem,2),size(mem,3)))
      inpz  = mem
    else if (ini .eq. 6) then
      allocate (inbx(size(mem,1),size(mem,2),size(mem,3)))
      inbx  = mem
    else if (ini .eq. 7) then
      allocate (inby(size(mem,1),size(mem,2),size(mem,3)))
      inby  = mem
    else if (ini .eq. 8) then
      print*,'entering first_time'
      allocate (inbz(size(mem,1),size(mem,2),size(mem,3)))
      inbz  = mem
      first_time = .false.
    end if
  end if
  !print*, 'dims of inrho, lb and ub',size(inrho,1),size(inrho,2),size(inrho,3),m(1)%lb,m(2)%lb,m(3)%lb,m(1)%ub,m(2)%ub,m(3)%ub
  do ix=m(1)%lb,m(1)%ub
   do iy=m(2)%lb,m(2)%ub
     do iz=m(3)%lb,m(3)%ub
        if (m(2)%n .gt. 1 .and. iy .gt. (m(2)%ub - 2) ) then
          mem(ix,iy,iz) = mem(ix,m(2)%ub - 2,iz)
        end if
        if (m(2)%n .gt. 1 .and. iy .lt. (m(2)%lb + 2) ) then
          mem(ix,iy,iz) = mem(ix,m(2)%lb + 2,iz)
        endif
        if (lower) then
          !mem(ix,iy,iz) = mem(m(1)%lb + 2,m(2)%ub/2,iz)
          if (ini .eq. 1) mem(m(1)%lb,iy,iz) = inrho(m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 2) mem(m(1)%lb,iy,iz) = ins  (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 3) mem(m(1)%lb,iy,iz) = inpx (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 4) mem(m(1)%lb,iy,iz) = inpy (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 5) mem(m(1)%lb,iy,iz) = inpz (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 6) mem(m(1)%lb,iy,iz) = inbx (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 7) mem(m(1)%lb,iy,iz) = inby (m(1)%lb,(m(2)%ub+1)/2,iz)
          if (ini .eq. 8) mem(m(1)%lb,iy,iz) = inbz (m(1)%lb,(m(2)%ub+1)/2,iz)
        else
          mem(m(1)%ub,iy,iz) = mem(m(1)%ub - 1,(m(2)%ub+1)/2,iz)
        end if
      end do
    end do
  end do
  if (ini .eq. 8) ini = 0
  ini = ini+1
  call trace_end
END SUBROUTINE condition

END MODULE boundary_mod