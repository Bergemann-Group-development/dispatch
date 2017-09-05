!> =============================================================================
!> =============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE variable_mod
  USE mhd_mod
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
    procedure:: gas_pressure
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
    procedure:: grav_potential
    procedure:: gas_velocity
  end type
  type(solver_t):: solver
  
  interface norm
    module procedure:: Norm_ArrayToArray, Norm_VectorToArray
  end interface
  interface down
    module procedure:: Down_ScalarToArray
  end interface
  interface MomToVel
    module procedure:: MomToVel_VectorToArray
  end interface
CONTAINS

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%d%q
  else
    pg = self%e%q * (self%gamma - 1.0)
  end if
END FUNCTION gas_pressure

!> =============================================================================
FUNCTION gas_velocity (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  v = MomToVel(self%s,self%d)
END FUNCTION gas_velocity

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%d%q)
END SUBROUTINE log_density

!> =============================================================================
SUBROUTINE log_pressure (self, v)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(gas_pressure(self))
END SUBROUTINE log_pressure

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(MomToVel(self%s,self%d))
END SUBROUTINE velocity_magnitude

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%B)
END SUBROUTINE magnetic_field_magnitude

!> =============================================================================
SUBROUTINE grav_potential (self, v)
  class (solver_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = self%gp%q
END SUBROUTINE grav_potential

!> =============================================================================
FUNCTION Norm_VectorToArray (in) RESULT (out)
  type(vector_t), intent(in) :: in
  real(kind=KindScalarVar), dimension(size(in%x%q,1),size(in%x%q,2),size(in%x%q,3)) :: out
  !-----------------------------------------------------------------------------
  out = sqrt(in%x%q(:,:,:)**2 + in%y%q(:,:,:)**2 + in%z%q(:,:,:)**2)
END FUNCTION Norm_VectorToArray

!> =============================================================================
FUNCTION Norm_ArrayToArray (in) RESULT (out)
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in) :: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3)) :: out
  !-----------------------------------------------------------------------------
  out = sqrt(in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)**2)
END FUNCTION Norm_ArrayToArray

!> =============================================================================
FUNCTION MomToVel_VectorToArray (s,d) RESULT (out)
  type(vector_t), intent(in):: s
  type(scalar_t), intent(in):: d
  real(kind=KindScalarVar), dimension(size(d%q,1),size(d%q,2),size(d%q,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = s%x%q / xdn_scalar(d%q)
  out(:,:,:,2) = s%y%q / ydn_scalar(d%q)
  out(:,:,:,3) = s%z%q / zdn_scalar(d%q)
END FUNCTION MomToVel_VectorToArray

!> =============================================================================
FUNCTION Down_ScalarToArray (in) RESULT (out)
  type(scalar_t), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in%q,1),size(in%q,2),size(in%q,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xdn_scalar (in%q)
  out(:,:,:,2) = ydn_scalar (in%q)
  out(:,:,:,3) = zdn_scalar (in%q)
END FUNCTION Down_ScalarToArray

!> =============================================================================
FUNCTION xdn_scalar (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, ione, jone, kone, n(3)
  !.............................................................................
  c = 0.5
  n = ubound(a)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)
  if (n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(1,:,:) = a(1,:,:)
    b(n(1),:,:) = a(n(1)-1,:,:)

    do iz=1,n(3)-kone
      do iy=1,n(2)-jone
        do ix=2,n(1)-ione
          b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
        end do
      end do
    end do
  else
    b = a
  end if
  if (n(2) > 1) b(:,n(2),:) = a(:,n(2)-1,:)
  if (n(3) > 1) b(:,:,n(3)) = a(:,:,n(3)-1)
  b(n(1),n(2),n(3)) = a(n(1)-ione,n(2)-jone,n(3)-kone)
END FUNCTION xdn_scalar

!> =============================================================================
FUNCTION ydn_scalar (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, ione, jone, kone, n(3)
  !.............................................................................
  c = 0.5
  n = ubound(a)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,1,:) = a(:,1,:)
    b(:,n(2),:) = a(:,n(2)-1,:)
    do iz=1,n(3)-kone
      do iy=2,n(2)-jone
        do ix=1,n(1)-ione
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
        end do
      end do
    end do
  else
    b = a
  end if
  if (n(3) > 1) b(:,:,n(3)) = a(:,:,n(3)-1)
  if (n(1) > 1) b(n(1),:,:) = a(n(1)-1,:,:)
  b(n(1),n(2),n(3)) = a(n(1)-ione,n(2)-jone,n(3)-kone)
END FUNCTION ydn_scalar

!> =============================================================================
FUNCTION zdn_scalar (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, ione, jone, kone, n(3)
!...............................................................................
  c = 0.5
  n = ubound(a)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,1) = a(:,:,1)
    b(:,:,n(3)) = a(:,:,n(3)-1)
    do iz=2,n(3)-kone
      do iy=1,n(2)-jone
        do ix=1,n(1)-ione
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
        end do
      end do
    end do
  else
    b = a
  end if
  if (n(1) > 1) b(n(1),:,:) = a(n(1)-1,:,:)
  if (n(2) > 1) b(:,n(2),:) = a(:,n(2)-1,:)
  b(n(1),n(2),n(3)) = a(n(1)-ione,n(2)-jone,n(3)-kone)
END FUNCTION zdn_scalar

!> =============================================================================
FUNCTION xdn_1vector (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, jone, kone
  !.............................................................................
  c = 0.5
  jone = merge(1, 0, size(a,2) > 1)
  kone = merge(1, 0, size(a,3) > 1)
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(1,:,:) = a(1,:,:)
    do iz=1,size(a,3)-kone
      do iy=1,size(a,2)-jone
        do ix=2,size(a,1)
          b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
        end do
      end do
    end do
  else
    b = a
  end if
END FUNCTION xdn_1vector

!> =============================================================================
FUNCTION ydn_2vector (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, kone, ione
  !.............................................................................
  c = 0.5
  kone = merge(1, 0, size(a,3) > 1)
  ione = merge(1, 0, size(a,1) > 1)
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,1,:) = a(:,1,:)
    do iz=1,size(a,3)-kone
      do iy=2,size(a,2)
        do ix=1,size(a,1)-ione
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
        end do
      end do
    end do
  else
    b = a
  end if
END FUNCTION ydn_2vector

!> =============================================================================
FUNCTION zdn_3vector (a) RESULT (b)
  real(kind=KindScalarVar), dimension(:,:,:):: a
  real(kind=KindScalarVar), dimension(size(a,1),size(a,2),size(a,3)):: b
  real(kind=KindScalarVar) :: c
  integer :: ix, iy, iz, ione, jone
  !.............................................................................
  c = 0.5
  ione = merge(1, 0, size(a,1) > 1)
  jone = merge(1, 0, size(a,2) > 1)
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,1) = a(:,:,1)
    do iz=2,size(a,3)
      do iy=1,size(a,2)-jone
        do ix=1,size(a,1)-ione
          b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
        end do
      end do
    end do
  else
    b = a
  end if
END FUNCTION zdn_3vector

END MODULE solver_mod
