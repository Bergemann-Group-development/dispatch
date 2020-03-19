!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : attributes for a single particle
!>   particle_list_t    : attributes for lists of particles
!>   particles_t        : attributes for particle solvers
!===============================================================================
MODULE particle_solver_mod
  USE io_mod
  USE trace_mod
  !USE timer_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE particles_mod
  USE link_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Everything needed by a particle solver
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: particle_solver_t
    type(particles_t)       :: particles
    type(particles_t)       :: exports
    class(patch_t), pointer :: patch => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: import_particles
    procedure:: is_outside
  end type
  integer:: verbose = 0
  real:: courant = 1.0
CONTAINS

!===============================================================================
!> Initialize particle solver
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t):: self
  integer:: iostat
  namelist /particle_solver_params/ verbose, courant
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%output, particle_solver_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  !.............................................................................
  real, dimension(3)         :: r, dr, v
  integer, dimension(3)      :: q, ir
  real   , dimension(3)      :: ods, t, p1
  class(particle_t), pointer :: pa, next
  real(8)                    :: tmin, tmax, time(2)
  logical                    :: outside
  integer                    :: ip, iter, np, it
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%update')
  call self%import_particles
  !-----------------------------------------------------------------------------
  ods = 1.0/self%patch%ds
  np = self%particles%n
  time(1) = wallclock()
  !-----------------------------------------------------------------------------
  !> Non-vectorized Vay pusher
  !-----------------------------------------------------------------------------
  ip = 0
  iter = 0
  pa => self%particles%head
  do while (associated(pa))
    !---------------------------------------------------------------------------
    ! Get pointer to next, in case the particle is removed
    !---------------------------------------------------------------------------
    next => pa%next
    !---------------------------------------------------------------------------
    ! Move the particle until it is ahead of the patch time
    !---------------------------------------------------------------------------
    it = pa%it
    do while (pa%time < self%patch%time)
      pa%r(:,it) = pa%r(:,it) + pa%dt*pa%v(:)
      pa%time = pa%time + pa%dt
    end do
    q = floor(pa%r(:,it))
    pa%r(:,it) = pa%r(:,it) - q
    pa%q(:,it) = pa%q(:,it) + q
    !---------------------------------------------------------------------------
    ! Move the particle to the export list if it falls outside
    !---------------------------------------------------------------------------
    outside = self%is_outside (pa)
    if (outside) then
      call self%exports%append (pa)
      call self%particles%remove (pa)
    end if
    pa => next
    ip = ip+1
  end do
  if (verbose > 0) then
    time(2) = wallclock()-time(1)
    write (stdout,'(i6,a,9p,f8.1,a)') &
      ip, ' updates', time(2)/iter, ' ns/update non-vectorized'
  end if
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Import particles from nbors export lists
!===============================================================================
SUBROUTINE import_particles (self)
  class(particle_solver_t)     :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: task, nbtask
  class(patch_t)   , pointer :: nbpatch
  class(particle_t), pointer :: pa, next
  real(8)                    :: p(3)
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%import_particles')
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (particle_solver_t)
      nbtask => nbor%task
      select type (nbtask)
      class is (particle_solver_t)
        nbpatch => nbtask%patch
        !-----------------------------------------------------------------------
        ! Loop over particles in export list and grab those that belong here
        !-----------------------------------------------------------------------
        pa => nbtask%exports%head
        do while (associated(pa))
          next => pa%next
          select type (pa)
          class is (particle_t)
            p = nbpatch%llc_cart + (pa%r(:,pa%it) + pa%q(:,pa%it))*nbpatch%ds
          end select
          if (self%patch%contains(p)) then
            call self%particles%append (pa)
            call task%exports%remove (pa)
          end if
          pa => next
        end do
      end select
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE import_particles

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_outside (self, particle)
  class(particle_solver_t):: self
  class(particle_t):: particle
  !.............................................................................
  is_outside = any(abs((particle%r + particle%q) - 0.5) > 0.5)
END FUNCTION is_outside

END MODULE particle_solver_mod
