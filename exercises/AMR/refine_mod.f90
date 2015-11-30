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
    procedure:: needed
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
  class(refine_t):: self
  real:: max_change
  namelist /refine_params/ max_change
  !.............................................................................
  call trace_begin ('refine_t%init')
  rewind (io%input)
  read (io%input, refine_params)
  if (io%master) write (*, refine_params)
  call trace_end
END SUBROUTINE init
!> =============================================================================
FUNCTION needed (self, task)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  logical:: needed
  !.............................................................................
  call trace_begin ('refine_t%needed', 2)
  needed = check_log_density    (self, task); if (needed) return
  needed = check_log_pressure   (self, task); if (needed) return
  needed = check_velocity       (self, task); if (needed) return
  needed = check_magnetic_field (self, task)
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_density (self, task) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  real, dimension(:,:,:), pointer:: v, d
  logical:: out
  !.............................................................................
  call trace_begin ('refine_t%check_log_densisty', 2)
  allocate (v(task%li(1):task%ui(1),task%li(2):task%ui(2),task%li(3):task%ui(3)))
  call task%log_density (v)
  out = check_max_diff (self, task, v)
  deallocate (v)
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_pressure (self, task) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  real, dimension(:,:,:), pointer:: p, v
  logical:: out
  !.............................................................................
  allocate (v(task%li(1):task%ui(1),task%li(2):task%ui(2),task%li(3):task%ui(3)))
  call task%log_pressure (v)
  out = check_max_diff (self, task, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_velocity (self, task) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  real, dimension(:,:,:), pointer:: v
  logical:: out
  !.............................................................................
  allocate (v(task%li(1):task%ui(1),task%li(2):task%ui(2),task%li(3):task%ui(3)))
  call task%velocity_magnitude (v)
  out = check_max_diff (self, task, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_magnetic_field (self, task) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  real, dimension(:,:,:), pointer:: v
  logical:: out
  !.............................................................................
  allocate (v(task%li(1):task%ui(1),task%li(2):task%ui(2),task%li(3):task%ui(3)))
  call task%magnetic_field_magnitude (v)
  out = check_max_diff (self, task, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_max_diff (self, task, v) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: task
  real, dimension(:,:,:), pointer:: v
  integer:: n(3), i1, i2, i3
  real:: max_diff
  logical:: out
  !.............................................................................
  out = .false.
  associate (li=>task%li, ui=>task%ui, n=>task%n)
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
        out = .true.
        return
      end if
    end do
    end do
    end do
  end associate
END FUNCTION
!> =============================================================================
END MODULE refine_mod
