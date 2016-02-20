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
  USE patches_mod
  USE experiment_mod
  implicit none
  private
  type:: refine_t
    type(refine_t), pointer:: next
    character(len=16):: name
    integer:: ivar = 0
    real:: max_diff
    real:: refine_from
    real:: refine_until
    real:: min_jeans
    real:: max_dln_rho
    integer:: levelmax
  contains
    procedure:: init
    procedure:: needed
    procedure:: check_max_diff
    procedure:: check_log_density
    procedure:: check_log_pressure
    procedure:: check_velocity
    procedure:: check_potential
    procedure:: check_jeans
    procedure:: check_magnetic_field
  end type
  real, public:: max_change=100
  type(refine_t), public:: refine
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  real:: refine_from=0.0, refine_until=1e9, min_jeans=0.0,max_dln_rho=2.3
  integer, save:: levelmax
  logical:: first_time = .true.
  namelist /refine_params/ msplit, refine_from, refine_until, levelmax, &
    min_jeans, max_dln_rho
  !.............................................................................
  call trace_begin ('refine_t%init')
  if (first_time) then
    first_time=.false.
    levelmax = io%levelmax
    rewind (io%input)
    read (io%input, refine_params)
    if (io%master) write (*, refine_params)
  end if
  self%max_diff = max_change
  self%refine_from = refine_from
  self%refine_until = refine_until
  self%levelmax = levelmax
  self%min_jeans = min_jeans
  self%max_dln_rho = max_dln_rho
  io%levelmax = min(io%levelmax, self%levelmax)
  call trace_end
END SUBROUTINE init
!> =============================================================================
FUNCTION needed (self, patch)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
  logical:: needed
  !.............................................................................
  needed = .false.
  if (io%verbose>1) &
    print*, 'needed:', patch%is_set(bits%not_leaf), patch%level >= self%levelmax, &
            patch%time < self%refine_from, patch%time > self%refine_until, &
            patch%time, self%refine_until
  if (patch%is_set(bits%not_leaf) .or. patch%level >= self%levelmax) return
  if (patch%time < self%refine_from .or. patch%time > self%refine_until) return
  call trace_begin ('refine_t%needed', 2)
  if (io%verbose>2) print*,'refine%needed? patch:', patch%id, patch%it, patch%new
  !
  !needed = check_log_density    (self, patch)
  !needed = self%check_potential (patch)
  if (self%min_jeans   > 0.0) needed = needed.or.self%check_Jeans (patch)
!  if (self%max_dln_rho > 0.0) needed = needed.or.self%check_max_diff (patch)
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
FUNCTION check_jeans (self, patch) RESULT(out)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v, d
  logical:: out
  real:: c_sound, jeans_length, ratio
  !.............................................................................
  call trace_begin ('refine_t%check_check_jeans', 2)
  c_sound = 1.0     ! FIXME
  jeans_length = c_sound/sqrt(maxval(patch%mem(:,:,:,patch%i_d,patch%it,1)*patch%G))
  ratio = jeans_length/minval(patch%ds)
  out = ratio < self%min_jeans
  if (io%verbose > 2) &
    print*,'Jeans:', patch%id, ratio, self%min_jeans, out
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_density (self, patch) RESULT(out)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
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
FUNCTION check_potential (self, patch) RESULT(out)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v, d
  logical:: out
  !.............................................................................
  call trace_begin ('refine_t%potential', 2)
  if (io%verbose>3) then
    print*, 'check_pot', patch%id
    print*, 'check_pot', patch%li, patch%ui
  end if
  allocate (v(1:patch%n(1),1:patch%n(2),1:patch%n(3)))
  call patch%grav_potential (v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION check_log_pressure (self, patch) RESULT(out)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
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
  class(experiment_t), pointer:: patch
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
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v
  logical:: out
  !.............................................................................
  allocate (v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3)))
  call patch%magnetic_field_magnitude (v)
  out = check_max_diff (self, patch, v)
  deallocate (v)
END FUNCTION
!> =============================================================================
FUNCTION check_max_diff (self, patch,v ) RESULT(out)
  class(refine_t):: self
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), pointer:: v
  !real:: v(patch%li(1):patch%ui(1),patch%li(2):patch%ui(2),patch%li(3):patch%ui(3))
  integer:: n(3), i1, i2, i3, m(3)=1
  real:: max_diff, diff
  logical:: out
  !.............................................................................
associate (n=>patch%n)
  out = .false.
  where (n==1)
    m = 0
  end where
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
