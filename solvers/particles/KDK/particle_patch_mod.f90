!===============================================================================
!> Particle patches are collections (lists) of particles histories, which may
!> be used for example to compute the forces at other particles at a particular
!> time.  They may also be evolved, by performing kick-drift-kick on each of its
!> members.
!>
!> Each particle patch has a pointer to list of nbor, whose particles are also
!> treated explicitly, and also a pointer to a list with all the other particle
!> patches, which are added as a collective force.  Procedures:
!>
!>  call particle_patch%init (n)                          : initialize n particles
!>  call particle_patch%force_on (p)                      : total force on particle
!>  call particle_patch%add_particle_forces_to (owner, p) : forces from all particles
!>  call particle_patch%add_patch_force_to (owner, p)     : collective force on particle
!>  call particle_patch%move_past (time)                  : move all particles past time
!===============================================================================
MODULE particle_patch_mod
  USE io_mod
  USE trace_mod
  USE particle_history_mod
  USE dll_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: particle_patch_t
    integer      :: id = 0
    integer      :: n = 0
    real         :: moves = 0
    integer(8)   :: forces = 0
    real         :: mass = 0.0
    type (dll_t) :: particles
    type (dll_t) :: nbor
    type (dll_t) :: other
    real(8), dimension(3) :: position = 0.0_8, size = 1.0_8
  contains
    procedure :: init
    procedure :: create
    procedure :: force_on
    procedure :: add_particle_forces_to
    procedure :: add_patch_force_to
    procedure :: move_past
    procedure :: randomize
    procedure :: info
    procedure :: test
    procedure :: delete
  end type
  real :: test_time  = 10.0
  integer :: n_star  = 100
  integer :: n_nbor  = 10
  integer :: n_other = 100
  integer :: verbose = 0
  logical, save:: do_test = .false.
  logical, save:: first_time = .true.
  !type (particle_patch_t), pointer:: particle_patches(:)
CONTAINS

!===============================================================================
!> Initialize a list of particle histories
!===============================================================================
SUBROUTINE init (self)
  class (particle_patch_t)           :: self
  !.............................................................................
  real, save    :: size = 1.0
  integer, save :: id = 0
  integer       :: iostat
  namelist /particle_patch_params/ do_test, n_star, n_nbor, n_other, test_time, &
                                   size, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%init')
  !$omp critical (particle_patch_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, particle_patch_params, iostat=iostat)
    if (io%master) write (io%output, particle_patch_params)
  end if
  id = id+1
  !$omp end critical (particle_patch_cr)
  self%id = id
  self%size = size
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize a list of particle histories
!===============================================================================
SUBROUTINE create (self, n)
  class (particle_patch_t)           :: self
  integer                            :: n
  !.............................................................................
  class (dll_node_t), pointer        :: item
  type (particle_history_t), pointer :: history
  integer                            :: i
  namelist /particle_patch_params/ do_test
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%create')
  self%mass = 0.0
  do i=1,n
    allocate (history, item)
    call history%init
    item%car => history
    call self%particles%append (item)
    self%mass = self%mass + history%mass
  end do
  if (verbose>0) &
    print *, 'ppatch created: id, n, mass =', history%id, self%particles%n, self%mass
  call trace%end()
END SUBROUTINE create

!===============================================================================
!> Add the force on a particle from particles in the local particle patch, and
!> those in the list of nbor particle patches
!===============================================================================
SUBROUTINE force_on (self, particle)
  class (particle_patch_t)            :: self
  class (particle_t), pointer         :: particle
  !.............................................................................
  class (particle_history_t), pointer :: history
  class (dll_node_t), pointer         :: item
  class (*), pointer                  :: particle_patch
  real(8)                             :: mass
  !-----------------------------------------------------------------------------
  ! Add forces from particles in the current particle patch
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%force_on',2)
  particle%a = 0.0
  call self%add_particle_forces_to (self, particle)
  !-----------------------------------------------------------------------------
  ! Add forces from particles in nbor patches.
  !-----------------------------------------------------------------------------
  item => self%nbor%head
  do while (associated(item))
    particle_patch => item%car
    select type (particle_patch)
    class is (particle_patch_t)
      call particle_patch%add_particle_forces_to (self, particle)
    end select
    item => item%next
  end do
  !-----------------------------------------------------------------------------
  ! Finally, add forces from other patches, further away
  !-----------------------------------------------------------------------------
  item => self%other%head
  do while (associated(item))
    particle_patch => item%car
    select type (particle_patch)
    class is (particle_patch_t)
      call particle_patch%add_patch_force_to (self, particle)
    end select
    item => item%next
  end do
  call trace%end()
END SUBROUTINE force_on

!===============================================================================
!> Add the force on a particle from particles in the given particle patch
!===============================================================================
SUBROUTINE add_particle_forces_to (self, owner, particle)
  class (particle_patch_t)    :: self
  class (particle_patch_t)    :: owner
  class (particle_t), pointer :: particle
  !.............................................................................
  class (dll_node_t), pointer :: item
  class (*), pointer          :: history
  !-----------------------------------------------------------------------------
  ! Add forces from particles in the current patch
  !-----------------------------------------------------------------------------
  item => self%particles%head
  do while (associated(item))
    history => item%car
    select type (history)
    class is (particle_history_t)
      history%forces = 0
      call history%add_particle_force_to (particle)
      owner%forces = owner%forces + history%forces
    end select
    item => item%next
  end do
END SUBROUTINE add_particle_forces_to

!===============================================================================
!> Add the force from a particle_patch to a particle in the owner patch
!===============================================================================
SUBROUTINE add_patch_force_to (self, owner, particle)
  class (particle_patch_t)    :: self
  class (particle_patch_t)    :: owner
  class (particle_t), pointer :: particle
  real                        :: r(3), ri3
  !-----------------------------------------------------------------------------
  r = particle%p + (owner%position - self%position)
  ri3 = (1./sqrt(sum(r**2)))**3
  particle%a = particle%a - (self%mass*ri3)*r
  !$omp atomic
  owner%forces = owner%forces + 1
END SUBROUTINE add_patch_force_to

!===============================================================================
!> Update the positions of all particles, maintaining their histories
!===============================================================================
SUBROUTINE move_past (self, time)
  class (particle_patch_t)    :: self
  real                        :: time
  !.............................................................................
  integer                     :: n, nh, it
  class (dll_node_t), pointer :: item
  class (*), pointer          :: history
  class (particle_t), pointer  :: particle
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%move_past')
  n = self%n
  item => self%particles%head
  do while (associated(item))
    history => item%car
    select type (history)
    class is (particle_history_t)
      it = history%it
      particle => history%particle(it)
      do while (particle%t < time)
        if (verbose>1) &
          print *, 'move_past', particle%id, particle%t, particle%dt, history%moves
        call self%force_on (particle)
        call particle%kick
        call particle%drift
        call self%force_on (particle)
        call particle%kick
        call history%rotate
        history%moves = history%moves+1
      end do
      self%moves = self%moves + history%moves
    end select
    item => item%next
  end do
  call trace%end()
END SUBROUTINE move_past

!===============================================================================
!> Randomize the particle velocities
!===============================================================================
SUBROUTINE randomize (self)
  class (particle_patch_t)    :: self
  !.............................................................................
  class (dll_node_t), pointer :: item
  class (*), pointer          :: history
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%randomize')
  item => self%particles%head
  do while (associated(item))
    history => item%car
    select type (history)
    class is (particle_history_t)
      call history%randomize
    class default
      error stop 'randomize'
    end select
    item => item%next
  end do
  call trace%end()
END SUBROUTINE randomize

!===============================================================================
!> Print information about the particles in the particle patch
!===============================================================================
SUBROUTINE info (self)
  class (particle_patch_t)    :: self
  !.............................................................................
  class (dll_node_t), pointer :: item
  class (*), pointer          :: history
  integer(8)                  :: moves
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_patch_t%info')
  moves = 0
  if (verbose>0) &
    print *, 'info:',self%id, self%particles%n, associated(self%particles%head)
  item => self%particles%head
  do while (associated(item))
    history => item%car
    select type (history)
    class is (particle_history_t)
      call history%info
      moves = moves+history%moves
    end select
    item => item%next
  end do
  if (verbose>0) &
    print *,'total moves:', real(moves)
  self%moves = moves
  call trace%end()
END SUBROUTINE info

!===============================================================================
!> Create n_nbor neighbors and n_other other particle patches, all with n_star
!> particles in each.  Then move all particles until time>test_time
!===============================================================================
SUBROUTINE test (self)
  class (particle_patch_t)          :: self
  !.............................................................................
  integer                           :: i
  class (particle_patch_t), pointer :: nbor
  class (dll_node_t), pointer       :: item
  real(8)                           :: wc
  !-----------------------------------------------------------------------------
  if (.not.do_test) return
  !$omp critical (test_cr)
  do_test = .false.
  print *, '===================== TEST: ========================================'
  print *, '--------------------- Generating central ppatches ------------------'
  call trace%begin ('particle_patch_t%test')
  call self%create (n_star)
  call self%randomize
  !-----------------------------------------------------------------------------
  ! nbors are placed at distances from 1 to 2.5 size units
  !-----------------------------------------------------------------------------
  print *, '--------------------- Generating nbor ppatches ---------------------'
  call self%nbor%init (' nbor')
  do i=1,n_nbor
    allocate (nbor, item)
    item%car => nbor
    call self%nbor%append (item)
    call nbor%init
    call nbor%create (n_star)
    nbor%position = 0d0
    do while (sqrt(sum(nbor%position**2)) < 1d0)
      call random_number (nbor%position)
      nbor%position = 5d0*(nbor%position-0.5d0)
    end do
    nbor%position = nbor%position*nbor%size
    if (verbose>0) &
      print *, 'nbor ', nbor%id, sqrt(sum(nbor%position**2)), nbor%mass
  end do
  !-----------------------------------------------------------------------------
  ! others are placed at distances from 2.5 to 10 size units
  !-----------------------------------------------------------------------------
  print *, '--------------------- Generating other ppatches --------------------'
  call self%other%init ('other')
  do i=1,n_other
    allocate (nbor, item)
    item%car => nbor
    call self%other%append (item)
    call nbor%init
    call nbor%create (n_star)
    nbor%position = 0d0
    do while (sqrt(sum(nbor%position**2)) < 2.5d0)
      call random_number (nbor%position)
      nbor%position = 10d0*(nbor%position-0.5d0)
    end do
    nbor%position = nbor%position*nbor%size
    if (verbose>0) &
      print *, 'other', nbor%id, sqrt(sum(nbor%position**2)), nbor%mass
  end do
  print *, '--------------------- Moving particles -----------------------------'
  print "(' direct forces from',i5,' particles, indirect from',i6)", &
    n_star*(n_nbor+1)-1, n_star*n_other
  self%moves = 0
  wc = wallclock()
  call self%move_past (test_time)
  wc = wallclock()-wc
  print "(' moves, forces/move:',f7.0,f8.0,f8.2,' ms/move',f8.2,' mus/force')", &
    self%moves, self%forces/self%moves, 1e3*wc/self%moves, 1e6*wc/self%forces
  print *, '--------------------- Cleanup --------------------------------------'
  call delete_ppatch_list (self%nbor)
  call delete_ppatch_list (self%other)
  call self%delete
  print *, '===================== TEST: Done ==================================='
  call trace%end()
  !$omp end critical (test_cr)
END SUBROUTINE test

!===============================================================================
SUBROUTINE delete_ppatch_list (dll)
  class (dll_t)              :: dll
  !.............................................................................
  type (dll_node_t), pointer :: item
  class (*)        , pointer :: ppatch
  !-----------------------------------------------------------------------------
  print *,'deleting ppatch list: ', trim(dll%name), '   with', dll%n, ' ppatches'
  item => dll%head
  do while (associated(item))
    if (associated(item%car)) then
      ppatch => item%car
      select type (ppatch)
      class is (particle_patch_t)
        call ppatch%delete
      end select
      deallocate (item%car)
    end if
    item => item%next
  end do
  call dll%init
END SUBROUTINE delete_ppatch_list

!===============================================================================
!> Delete all particle history instances in the linked list self%particles
!===============================================================================
SUBROUTINE delete (self)
  class (particle_patch_t)   :: self
  !.............................................................................
  type (dll_node_t), pointer :: item
  class (*)        , pointer :: history
  !-----------------------------------------------------------------------------
  if (verbose>0) &
    print *,'deleting particles in ppatch', self%id
  item => self%particles%head
  do while (associated(item))
    if (associated(item%car)) then
      history => item%car
      select type (history)
      class is (particle_history_t)
        call history%delete
      end select
      deallocate (item%car)
    end if
    item => item%next
  end do
  call self%particles%init
END SUBROUTINE delete

END MODULE particle_patch_mod
