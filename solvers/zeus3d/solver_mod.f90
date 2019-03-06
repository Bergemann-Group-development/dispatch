!> =============================================================================
!> =============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE variable_mod
  USE mhd_mod
  USE mesh_mod, only: mesh_types
  USE validate_mod
  implicit none
  private
  type, public, extends(mhd_t):: solver_t
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
    procedure:: grav_potential
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
    generic:: gas_velocity => gas_velocity_vector, gas_velocity_scalar
    procedure:: p2U
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

!> =====================================================================
!> Organize calls to the extras and the mhd solver init
!> =====================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.....................................................................
  call self%mhd_t%init
  call self%extras_t%init
  call validate%init
END SUBROUTINE init

!> =====================================================================
!> Organize calls to the extras and the mhd solver
!> =====================================================================
SUBROUTINE update (self)
  class(solver_t) :: self
  associate (d=>self%mem(:,:,:,self%idx%d,self%it,1))
  !.....................................................................
  call self%extras_t%pre_update
  call validate%check (self, d, 'before update')
  call self%mhd_t%update
  call validate%check (self, d, ' after update')
  call self%extras_t%post_update
  end associate
END SUBROUTINE update

!> =============================================================================
SUBROUTINE p2U (self, v, it)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: v
  integer:: it
  !---------------------------------------------------------------------
  v = self%gas_velocity_vector()
END SUBROUTINE p2U

!> =============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  v = MomToVel(self)
END FUNCTION gas_velocity_vector

!> =============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  class (solver_t):: self
  integer:: idir
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,self%idx%d,self%it,1))
  select case (idir)
  case (1)
    v(:,:,:) = self%mem(:,:,:,self%idx%px,self%it,1) / xdn_scalar(d)
  case (2)
    v(:,:,:) = self%mem(:,:,:,self%idx%py,self%it,1) / ydn_scalar(d)
  case (3)
    v(:,:,:) = self%mem(:,:,:,self%idx%pz,self%it,1) / zdn_scalar(d)
  case default
    error stop "solver_mod::gas_velocity_scalar:: invalid value of idir"
  end select
  end associate
END FUNCTION gas_velocity_scalar

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%mem(:,:,:,self%idx%d,self%it,1))
END SUBROUTINE log_density

!> =============================================================================
SUBROUTINE log_pressure (self, v)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%gas_pressure())
END SUBROUTINE log_pressure

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (solver_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(MomToVel(self))
END SUBROUTINE velocity_magnitude

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (solver_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1))
END SUBROUTINE magnetic_field_magnitude

!> =============================================================================
SUBROUTINE grav_potential (self, v)
  class (solver_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = self%mem(:,:,:,self%idx%phi,self%it,1)
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
FUNCTION MomToVel_VectorToArray (self) RESULT (out)
  class(solver_t):: self
  real(kind=KindScalarVar), dimension(self%gn(1),self%gn(2),self%gn(3),3):: out
  integer:: i1, i2, i3
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,self%idx%d,self%it,1))
  out(:,:,:,1) = self%mem(:,:,:,self%idx%px,self%it,1) / xdn_scalar(d)
  out(:,:,:,2) = self%mem(:,:,:,self%idx%py,self%it,1) / ydn_scalar(d)
  out(:,:,:,3) = self%mem(:,:,:,self%idx%pz,self%it,1) / zdn_scalar(d)
  ! in curvilinear coords., it is the angular momentum that is stored
  if (self%mesh_type /= mesh_types%Cartesian) then
    do i3=self%mesh(3)%lb,self%mesh(3)%ub
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
        do i1=self%mesh(1)%lb,self%mesh(1)%ub
          !out(i1,i2,i3,1) = out(i1,i2,i3)
          out(i1,i2,i3,2) = out(i1,i2,i3,2) / self%mesh(1)%h2c(i1)
          out(i1,i2,i3,3) = out(i1,i2,i3,3) / (self%mesh(1)%h31c(i1) * self%mesh(2)%h32c(i2))
        end do
      end do
    end do
  end if

  end associate
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
