!===============================================================================
!> Experiment initialization and update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE data_io_mod
  USE solver_mod
  !USE particle_patch_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: input
    procedure:: update
    procedure:: dnload
    procedure:: task_log
  end type
  !type (particle_patch_t):: particle_patch
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%solver_t%init
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Setup a uniform initial state
!===============================================================================
SUBROUTINE input (self, ok)
  class(experiment_t):: self
  logical:: ok
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%input')
  call data_io%input (self, ok)
  call trace_end
END SUBROUTINE input

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real:: dmin
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if (self%verbose > 0) then
    dmin = minval(self%mem(:,:,:,self%idx%d,self%it,1))
    if (dmin <= 0.0) then
      print *, self%id, 'density is non-positive before update', &
        dmin, minloc(self%mem(:,:,:,self%idx%d,self%it,1))
    end if
  end if
  call self%solver_t%update
  if (self%verbose > 0) then
    dmin = self%fminval(self%mem(:,:,:,self%idx%d,self%new,1))
    if (dmin <= 0.0) then
      print *, self%id, 'density is non-positive after update', &
        dmin, self%fminloc(self%mem(:,:,:,self%idx%d,self%new,1))
    end if
  end if
  call self%task_log ('aft upd')
  call trace_end
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%dnload')
  call self%task_log ('bef dnl')
  call self%solver_t%dnload (only)
  call self%task_log ('aft dnl')
  call trace_end
END SUBROUTINE dnload

!===============================================================================
!===============================================================================
SUBROUTINE task_log (self, l)
  class(experiment_t):: self
  character(len=*):: l
  real:: vmax(self%nv)
  real:: vmin(self%nv)
  integer:: iv
  !----------------------------------------------------------------------------
  if (io%task_logging > 4) then
    do iv=1,self%nv
        vmax(iv) = self%fmaxval(self%mem(:,:,:,iv,self%it,1))
        vmin(iv) = self%fminval(self%mem(:,:,:,iv,self%it,1))
    end do
    !$omp critical (tasklog_cr)
    write(io_unit%task,'(f12.6,i4,i6,i4,f12.6,1p,e11.2,2(2x,5e10.2),2x,a)') &
      wallclock(), omp%thread, self%id, self%level, self%time, self%dtime, &
      vmin, vmax, l
    flush (io_unit%task)
    !$omp end critical (tasklog_cr)
  end if
END SUBROUTINE task_log

END MODULE
