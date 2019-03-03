!> =============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object solver_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE initial_mod
  USE timer_mod
  USE link_mod
  USE index_mod
  USE task_mod
  USE trace_particles_mod
  USE validate_mod
  !USE particle_patch_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(trace_particles_t):: trace_particles
    real, pointer:: v(:,:,:,:) => null()
  contains
    procedure:: init
    procedure:: update
  end type
  !type (particle_patch_t):: particle_patch
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, pointer:: d(:,:,:), p(:,:,:,:)
  integer:: i
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%solver_t%init
  self%d0 = self%initial%d0
  call self%set (bits%static)
  if (self%initial%time > 0d0) then
    self%time = self%initial%time
    self%t(self%it) = self%initial%time
  end if
  !----------------------------------------------------------------------------
  ! Startup support for trace particles
  !----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  allocate (self%v(self%gn(1),self%gn(2),self%gn(3),3))
  do i=1,3
    self%v(:,:,:,i) = self%mem(:,:,:,self%idx%px+i-1,self%it,1) &
                    / self%mem(:,:,:,self%idx%d     ,self%it,1)
  end do
  call self%trace_particles%init_nt (self%nt)
  call self%trace_particles%add (self%mesh, self%time, self%dtime, d, p)
  !----------------------------------------------------------------------------
  ! Support (disabled for now) for particle patches (carrying sink particles)
  !----------------------------------------------------------------------------
  !call particle_patch%init
  !call particle_patch%test
  call validate%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/solver_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the solver_t or experiment_t
!> level, by simply adding a courant_condition procedure there).
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  real, pointer:: d(:,:,:)
  integer:: i
  type(link_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  class(trace_particles_t), pointer:: nbpart
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if (first_time) then
    !$omp critical (cleanup_cr)
    if (first_time) then
      first_time = .false.
      call self%initial%cleanup
    end if
    !$omp end critical (cleanup_cr)
  end if
  associate (d => self%mem(:,:,:,self%idx%d,self%it,1))
  call validate%dump (self, d, 'before update')
  call self%solver_t%update
  call validate%dump (self, d, ' after update')
  end associate
  !----------------------------------------------------------------------------
  ! Startup support for trace particles; this could be improved, by reusing
  ! the velocity from the primitive variables in the solver.  In any case,
  ! this part is solver specific, so the support should ultimately be moved
  ! to the solver subdirectories.
  !----------------------------------------------------------------------------
  do i=1,3
    self%v(:,:,:,i) = self%mem(:,:,:,self%idx%px+i-1,self%it,1) &
                    / self%mem(:,:,:,self%idx%d     ,self%it,1)
  end do
  call self%trace_particles%update (self%mesh, self%time, self%dtime, &
    self%mem(:,:,:,self%idx%d,self%it,1), self%v)
  !----------------------------------------------------------------------------
  ! If particles have been acciumlated in the "export" particle list, loop over
  ! nbors and transfer particles to their "import" particle lists.  This is to
  ! minimize the need for locking or critical regions.  It should be enough to
  ! protect the import lists with OMP locks, while they are being updated, and
  ! while they are being read.
  !----------------------------------------------------------------------------
  if (self%trace_particles%export%n > 0) then
    if (self%trace_particles%verbose > 0) &
      write (io%output,*) &
        self%id, self%trace_particles%export%n, ' particles to export'
    nbor => self%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type (nbtask)
      type is (experiment_t)
      nbpart => nbtask%trace_particles
      nbpart%offset = nint((mod(self%position-nbtask%position+1.5*self%box,self%box)-0.5*self%box)/self%ds)
      call self%trace_particles%trade (nbpart)
      end select
      nbor => nbor%next
    end do
  end if
  call trace_end
END SUBROUTINE update

END MODULE
