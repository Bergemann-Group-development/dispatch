!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE timestep_mod
  implicit none
  private
  type:: refine_t
    type(refine_t), pointer:: next
    character(len=16):: name
    integer:: ivar = 0
    real:: max_diff = 1.0
    real:: refine_from=999.
    real:: refine_until=0.0
    integer:: levelmax=100
  contains
    procedure:: init
    procedure:: needed
    procedure:: check_max_diff
    procedure:: check_log_density
    procedure:: check_log_pressure
    procedure:: check_velocity
    procedure:: check_magnetic_field
  end type
  integer, public:: msplit=2                   ! shared with patch_list_mod
  real, public:: max_change=100
  type(refine_t), public:: refine
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  real:: refine_from = 0.0, refine_until = 1e9
  integer:: levelmax = 100
  logical:: first_time = .true.
  namelist /refine_params/ max_change, msplit, refine_from, refine_until, levelmax
  !.............................................................................
  call trace_begin ('refine_t%init')
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, refine_params)
    if (io%master) write (*, refine_params)
  end if
  self%max_diff = max_change
  self%refine_from = refine_from
  self%refine_until = refine_until
  self%levelmax = levelmax
  call trace_end
END SUBROUTINE init
!> =============================================================================
FUNCTION needed (self, patch)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  logical:: needed
  !.............................................................................
   needed = .false.
  if (patch%is_set(bits%not_leaf) .or. patch%level >= self%levelmax) return
  if (patch%time < self%refine_from .or. patch%time > self%refine_until) return
  call trace_begin ('refine_t%needed', 2)
  if (io%verbose>2) print*,'refine%needed? patch:', patch%id, patch%it, patch%new
  needed = check_log_density    (self, patch)
  if (needed) then
    call trace_end
    return
  end if
  !needed = check_log_pressure   (self, patch); if (needed) return
  !needed = check_velocity       (self, patch); if (needed) return
  !needed = check_magnetic_field (self, patch)
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_density (self, patch) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v, d
  logical:: out
  !.............................................................................
  call trace_begin ('refine_t%check_log_densisty', 2)
  if (io%verbose>3) then
    print*, 'check_log', patch%id
    print*, 'check_log', patch%li, patch%ui
  end if
  !allocate (v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3)))
  allocate (v(1:patch%n(1),1:patch%n(2),1:patch%n(3)))
  call patch%log_density (v)
  !print*, 'check_log', minval(v), maxval(v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_pressure (self, patch) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  real, dimension(:,:,:), pointer:: p, v
  logical:: out
  !.............................................................................
  allocate (v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3)))
  call patch%log_pressure (v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_velocity (self, patch) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v
  logical:: out
  !.............................................................................
  allocate (v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3)))
  call patch%velocity_magnitude (v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_magnetic_field (self, patch) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v
  logical:: out
  !.............................................................................
  allocate (v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3)))
  call patch%magnetic_field_magnitude (v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_max_diff (self, patch, v) RESULT(out)
  class(refine_t):: self
  class(timestep_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v
  !real:: v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3))
  integer:: n(3), i1, i2, i3, m(3)=1
  real:: max_diff, diff
  logical:: out
  !.............................................................................
  out = .false.
  where (patch%n==1)
    m = 0
  end where
  associate (li=>patch%li, ui=>patch%ui, n=>patch%n)
    do i1=-m(1),m(1)
    do i2=-m(2),m(2)
    do i3=-m(3),m(3)
      if (i1==0.and.i2==0.and.i3==0) cycle
      diff = maxval(abs(v(1+m(1)+i1:n(1)-m(1)+i1, &
                          1+m(2)+i2:n(2)-m(2)+i2, &
                          1+m(3)+i3:n(3)-m(3)+i3) &
                      - v(1+m(1)   :n(1)-m(1),   &
                          1+m(2)   :n(2)-m(2),   &
                          1+m(3)   :n(3)-m(3))))
      if (diff >= self%max_diff) then
        if (io%verbose > 1) &
          print *,'refine_t%chck_max_diff: id, index, diff, max_diff', &
          patch%id, diff, self%max_diff
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
