!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE io_mod
  USE trace_mod
  USE timestep_mod
  implicit none
  private
  type:: refine_t
    type(refine_t), pointer:: next
    character(len=16):: name
    integer:: ivar = 0
    real:: max_diff = 1.0
  contains
    procedure:: init
    procedure:: check
    procedure:: check_max_diff
    procedure:: check_log_density
    procedure:: check_log_pressure
    procedure:: check_velocity
    procedure:: check_magnetic_field
  end type
  type(refine_t), public:: refine
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  type(refine_t):: self
  namelist /refine_params/ max_change, refine_on
  !.............................................................................
  rewind (io%input)
  read (io%input, refine_params)
  if (io%master) write (*, refine_params)
END SUBROUTINE init
!> =============================================================================
FUNCTION check (self, task)
  class(refine_t):: self
  logical:: check
  !.............................................................................
  check = check_log_density    (self, task); if (check) return
  check = check_log_pressure   (self, task); if (check) return
  check = check_velocity       (self, task); if (check) return
  check = check_magnetic_field (self, task)
END FUNCTION
!> =============================================================================
FUNCTION check_log_density (self, task)
  class(refine_t):: self
  class(timestep_t):: task
  real, dimension(:,:,:), pointer:: v, d
  logical:: check_log_density
  !.............................................................................
  allocate (v)
  call task%log_density (v)
  call check_max_diff (self, task, v)
  deallocate (p)
END FUNCTION
!> =============================================================================
FUNCTION check_log_pressure (self, task)
  class(refine_t):: self
  class(timestep_t):: task
  type(mesh_t):: m
  real, dimension(:,:,:), pointer:: p, v
  logical:: check_log_pressure
  !.............................................................................
  allocate (v)
  call task%log_pressure (v)
  call check_max_diff (self, task, v)
  deallocate (p)
END FUNCTION
!> =============================================================================
FUNCTION check_velocity (self, task)
  class(refine_t):: self
  class(timestep_t):: task
  real, dimension(:,:,:), pointer:: v
  logical:: check_velocity
  !.............................................................................
  allocate (v)
  call task%velocity_magnitude (v)
  call check_max_diff (self, task, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_magnetic_field (self, task)
  class(refine_t):: self
  class(timestep_t):: task
  real, dimension(:,:,:), pointer:: v
  logical:: check_magnetic_field
  !.............................................................................
  allocate (v)
  call task%magnetic_field_magnitude (v)
  call check_max_diff (self, task, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_max_diff (self, task, v)
  class(refine_t):: self
  real, dimension(:,:,:), pointer:: v
  integer:: n(3), i1, i2, i3
  real:: max_diff
  logical:: check_max_diff
  !.............................................................................
  check_max_diff = .false.
  assoiate (li=>task%mesh%li, ui=>task&mesh&ui, n=>task%mesh%n)
    do i1=-1,1
    do i2=-1,1
    do i3=-1,1
      if (i1==0.and.i2==0.and.i3==0) cycle
      if (maxval(abs(v(li(1)+1+i1:ui(1)-1+i1, &
                       li(1)+1+i2:ui(2)-1+i2, &
                       li(1)+1+i3:ui(3)-1+i3) &
                   - v(li(1)+1   :ui(1)-1,   &
                       li(2)+1   :ui(2)-1,   &
                       li(3)+1   :ui(3)-1))) &
          > self%max_diff) then
        check_max_diff = .true.
        return
      end if
    end do
    end do
    end do
  end associate
END FUNCTION
!> =============================================================================
END MODULE refine_mod
