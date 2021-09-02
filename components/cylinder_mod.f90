!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Cylinder-organized Cartesian patches.  Patches are assumed to be cubical,
!> With radii defined at the outer edges, in such a way that there are exactly
!> n_phi patches covering 2\pi, and the size of patches is thus exactly equal to
!> radius*dphi*cos(dphi/2), whith radius counted at the outer edge. The ratio
!> between subsequent radii is thus (1-dphi*cos(dhpi/2)) when decreasing, and
!> the inverse of that when increasing, and hence we can define dlnr to be
!> dlnr = -ln(1-dphi*cos(dphi/2)), which is close to, but slightly smaller than,
!> dlnr = -ln(1-dphi).
!>
!> Since we may want to specify dlnr rather than dphi, we can compute a first
!> approximation of dphi as dphi = 1-exp(-dlnr), which is then marginally too
!> large.  Correcting dphi to be the next smaller value of 2*pi/n_phi, we can
!> then recompute an exact value of dlnr = -ln(1-dphi*cos(dphi/2))
!>
!> The inner edge overlap is approximately n*dphi cells, which is assumed to be
!> less than one cell, and hence dphi must be less than 1/n, where n is the
!> number of internal cells in a patch.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE cylinder_mod
  USE omp_lib
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE timer_mod
  USE mpi_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE halo_mod
  USE experiment_mod
  USE region_mod
  USE scaling_mod
  USE math_mod
  USE task_list_mod
  USE hash_table_mod
  USE bits_mod
  USE units_mod
  USE disk_model_mod
  USE shared_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Data type representing a single cylinder shell
  !-----------------------------------------------------------------------------
  type, public:: cylinder_t
    type(cylinder_t), pointer:: next => null()  ! pointer to the next cylinder
    type(cylinder_t), pointer:: prev => null()  ! pointer to the previous cylinder
    integer:: n_phi                             ! azimuthal number of patches
    integer:: n_z                               ! vertical patch layers
    real(8):: radius                            ! cylinder radius
    !real(8):: height                            ! cylinder height
    real(8):: omega                             ! angular rotation speed
    real(8):: dphi                              ! azimuth step
    real(8):: m_star                            ! stellar mass
    type(task_list_t), pointer:: task_list => null()  ! demoted task list
    real(8), pointer, contiguous, dimension(:,:,:,:)   :: r => null()
    real(8), pointer, contiguous, dimension(:,:,:,:,:) :: r_hat, f_grav, r_inv
    real(8), pointer, contiguous, dimension(:,:,:,:)   :: phi_hat, v_patch
  contains
    procedure:: init => init_cylinder
  end type
  !-----------------------------------------------------------------------------
  ! Data type representing a series of cylinder shells
  !-----------------------------------------------------------------------------
  type, public, extends(task_list_t):: cylinders_t
    integer                    :: n_r               ! number of cylinders
    real(8)                    :: lnr_min, dlnr     ! radial spaching
    type(cylinder_t),  pointer :: cylinder(:) => null()
    class(task_list_t), pointer:: task_list   => null()
  contains
    procedure:: init => init_cylinders
    procedure:: init_nbors
  end type
  !-----------------------------------------------------------------------------
  ! Shared module parameters
  !-----------------------------------------------------------------------------
  integer, save:: n_task(2)=0, n_z=0, n_nbors=0
  real,    save:: n_h=1.0
  integer, save:: n_r=10, n_phi=0, verbose=0
  real(8), save:: dlnr=0.1, r_min=1.0, dphi
CONTAINS

!===============================================================================
!> Initialize a number of cylinder shells
!===============================================================================
SUBROUTINE init_cylinders (self, name)
  class(cylinders_t), target:: self
  character(len=*), optional:: name
  namelist /cylinder_params/ verbose, n_r, n_h, dlnr, r_min
  !.............................................................................
  class(task_list_t), pointer:: task_list
  class(task_t), pointer:: task
  integer:: i_r, iostat, n, gn
  real(8):: r_in, r_out, size, radius, cos2, sin2, v_phi
  class(cylinder_t), pointer:: cylinder, prev, next
  !-----------------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinders_t%init')
  call io%lock%set ('cylinders_t%init')
  rewind (io%input)
  read (io%input, cylinder_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, cylinder_params)
  end if
  call io%namelist_warning ('cylinder_params', iostat)
  write (io%namelist, cylinder_params)
  !call self%hash_table%init ()
  call disk_model%init ()
  call io%lock%unset ('cylinders_t%init')
  self%n_r   = n_r
  self%lnr_min = log(r_min)
  !-----------------------------------------------------------------------------
  ! Compute exact values of dphi and dlnr
  !-----------------------------------------------------------------------------
  dphi = 1d0-exp(-dlnr)
  n_phi = 2d0*math%pi/dphi + 1
  dphi = 2d0*math%pi/n_phi
  cos2 = cos(dphi/2d0)
  sin2 = sin(dphi/2d0)
  dlnr = log(cos2/(1d0-2d0*sin2/cos2))
  !-----------------------------------------------------------------------------
  ! Increase patch size by 1%, to avoid corner cases with missing guard cell
  !-----------------------------------------------------------------------------
  dlnr = dlnr*1.01
  self%dlnr = dlnr
  !-----------------------------------------------------------------------------
  ! Create the cylinders
  !-----------------------------------------------------------------------------
  allocate (self%cylinder(n_r))
  !-----------------------------------------------------------------------------
  ! Loop over radii, creating cylinder shells
  !-----------------------------------------------------------------------------
  call self%task_list_t%init
  nullify (prev)
  task_list => self
  self%task_list => task_list
  do i_r=1,n_r
    cylinder => self%cylinder(i_r)
    cylinder%dphi = dphi
    cylinder%n_phi = n_phi
    cylinder%task_list => task_list
    r_in  = exp(self%lnr_min + (i_r-1)*self%dlnr)
    r_out = exp(self%lnr_min + i_r*self%dlnr)
    size  = (r_out-r_in)/cos2
    call cylinder%init (i_r, r_out, size, self%hash_table)
    !---------------------------------------------------------------------------
    ! Set up doubly linked list relations
    !---------------------------------------------------------------------------
    if (associated(prev)) &
      prev%next => cylinder
      cylinder%prev => prev
    prev => cylinder
  end do
  !-----------------------------------------------------------------------------
  ! The angle opening needs to be smaller than 1 cell, so we need n*dphi < 1
  !-----------------------------------------------------------------------------
  task => self%task_list%head%task
  select type (task)
  class is (patch_t)
    n  = task%n(1)
    gn = task%gn(1)
    v_phi = task%velocity(2)
  end select
  write (stdout,'(1x,a,3i5,2i9)') 'n_r,n_phi,n_z,n_task=', n_r,n_phi,n_z,n_task
  write (stdout,'(4(1x,a,f6.3,3x))') &
    'dlnr:', dlnr, 'dphi:', dphi, 'n*dphi:', n*dphi, 'gn*dphi:', gn*dphi
  write (stdout,'(1x,a,f5.1,4x,a,f6.2)') &
    'v_phi:', v_phi, 'dv_r:', v_phi*dlnr/2.0
  !-----------------------------------------------------------------------------
  ! Initialize  nbor lists
  !-----------------------------------------------------------------------------
  !$omp parallel default(none) shared(self)
  call init_all_nbors (self)
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Remove external tasks
  !-----------------------------------------------------------------------------
  !FIXME: call  remove_external (self)
  !-----------------------------------------------------------------------------
  ! Initialize the nbor lists again, after tasks have been removed
  !-----------------------------------------------------------------------------
  call self%task_list%reset_status
  !$omp parallel default(none) shared(self)
  call init_all_nbors (self)
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Reset status bits and finalize, allocating task memory
  !-----------------------------------------------------------------------------
  call self%task_list%reset_status
  !$omp parallel default(none) shared(self)
  call finalize_tasks (self)
  !$omp end parallel
  !-----------------------------------------------------------------------------
  call self%task_list%count_status
  call self%task_list%set_io_counts
  call self%task_list%info
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init_cylinders

!===============================================================================
!> Remove external tasks, and add remaining tasks to a new task list
!===============================================================================
SUBROUTINE remove_external (self)
  class(cylinders_t):: self
  !-----------------------------------------------------------------------------
  class(task_list_t), pointer:: task_list
  class(link_t), pointer :: link, next, halo
  class(task_t), pointer:: task
  integer:: n
  !-----------------------------------------------------------------------------
  ! Identify real, virtual, and external tasks
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinders_t%remove_external')
  call self%task_list%reset_status
  allocate (task_list)
  call task_list%init ()
  link => self%task_list%head
  do while (associated(link))
    next => link%next
    task => link%task
    if (task%is_set (bits%external)) then
      if (verbose > 0) write (io_unit%mpi,*) 'remove', task%id
      call self%hash_table%free (task%key)
      call self%task_list%remove (link)
      call task%dealloc
      deallocate (task)
      nullify (task)
    else
      !-------------------------------------------------------------------------
      ! This is a work-around, since the task_list%remove does not seem to
      ! work as intended.   Calling tmp_list%append (task) does not make the
      ! link a class halo_t data type.
      !-------------------------------------------------------------------------
      allocate (halo_t:: halo)
      halo%task => task
      select type (task)
      class is (patch_t)
      task%link => halo
      end select
      call task_list%append (halo)
      if (verbose > 1) write (io_unit%mpi,*) '  keep', task%id
    end if
    deallocate (link)
    link => next
  end do
  self%task_list => task_list
  call trace%end ()
END SUBROUTINE remove_external

!===============================================================================
!> Finalize the task setup, allocating memory in remaining tasks.
!> This procedure is called in a parallel region, in order to share the work
!> and to make sure memory allocation is spread out over core-affinity.
!===============================================================================
SUBROUTINE finalize_tasks (self)
  class(cylinders_t) :: self
  !.............................................................................
  class(link_t), pointer :: link
  class(task_t), pointer:: task
  integer:: n, l(3), u(3), key(4)
  real(8):: time, centre_nat(3), phi
  !-----------------------------------------------------------------------------
  ! Allocate tasks that remain, making sure to spread allocation over threads
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinders_t%finalize_tasks')
  call timer%tic (time)
  n = 0
  link => self%task_list%head
  do while (associated(link))
    if (mod(n,omp%nthreads) == omp%thread) then
      task => link%task
      select type (task)
      class is (experiment_t)
        key = task%key              ! FIXME
        centre_nat = task%centre_nat
        call task%init ()
        task%centre_nat = centre_nat
        task%key = key              ! FIXME
        task%mesh%d = task%ds
        task%mesh%p = task%position
        task%mesh%llc_cart = task%llc_cart
        !-----------------------------------------------------------------------
        ! Unit vector matrix and inverse
        !-----------------------------------------------------------------------
        phi = centre_nat(2)
        task%erot(:,1) = [ cos(phi), sin(phi), 0d0]
        task%erot(:,2) = [-sin(phi), cos(phi), 0d0]
        task%erot(:,3) = [0d0, 0d0, 1d0]
        !-----------------------------------------------------------------------
        ! Boundaries and initial conditions
        !-----------------------------------------------------------------------
        call task%init_bdries ()
        l = 1
        u = task%gn
        call task%initial_condition (l, u)
      end select
    end if
    n = n+1
    link => link%next
  end do
  call timer%toc ("time to initialize tasks", n, time)
  call trace%end ()
END SUBROUTINE finalize_tasks

!===============================================================================
!> Initialize a number of Cartesian patches in a cylinder arrangement
!===============================================================================
SUBROUTINE init_cylinder (self, i_r, r_out, size, hash_table)
  class(cylinder_t) :: self
  integer           :: i_r
  real(8)           :: r_out, size
  type(hash_table_t):: hash_table
  !.............................................................................
  class(experiment_t), pointer :: cpatch
  class(halo_t)      , pointer :: halo
  class(*)           , pointer :: ptr
  real(8)                      :: dphi, phi, z
  integer                      :: i_z, i_phi, unit, l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinder_t%init')
  if (verbose > 1) then
    unit = stdout
  else
    unit = stdlog
  end if
  !-----------------------------------------------------------------------------
  ! Patches align at outer radial edges in azimuthal dimension, i.e. azimuthal
  ! size at outer edges is equal to radial size
  !-----------------------------------------------------------------------------
  self%radius = r_out-size/2d0
  self%omega = disk_model%v_kepler(self%radius)/self%radius
  n_z = calc_n_z (self%radius, size, n_h)
  self%n_z = n_z
  do i_z=1,n_z
    !$omp parallel do default(none) private(cpatch,halo,phi,z,l,u,ptr) &
    !$omp shared(i_r,i_z,n_r,n_phi,n_z,math,size,self,verbose,unit,bits,mpi) &
    !$omp shared(n_task,hash_table)
    do i_phi=1,n_phi
      write (stdlog,*) 'init_cylinder: i_phi,thread=', i_phi, omp_get_thread_num()
      allocate (cpatch, halo)
      !-------------------------------------------------------------------------
      ! Local coordinate system of cpatch.  The box size is set arbitrarily,
      ! but large, to avoid any unintended wrapping of coordinates, and to set
      ! a basis for AMR level assignment
      !-------------------------------------------------------------------------
      cpatch%box = 4.0*math%pi
      cpatch%size = size
      phi = (i_phi-1)*self%dphi
      z = (i_z - (1+n_z)*0.5d0) * cpatch%size(3)
      !-------------------------------------------------------------------------
      ! Unit vectors and inverse of unit vector matrix
      !-------------------------------------------------------------------------
      allocate (cpatch%erot(3,3))
      cpatch%erot(:,1) = [ cos(phi), sin(phi), 0d0]
      cpatch%erot(:,2) = [-sin(phi), cos(phi), 0d0]
      cpatch%erot(:,3) = [0d0, 0d0, 1d0]
      cpatch%omega(3) = self%omega
      !-----------------------------------------------------------------------
      ! Global Cartesian coordinate system
      !-----------------------------------------------------------------------
      cpatch%box = self%radius*2d0*math%pi
      cpatch%radius = self%radius
      cpatch%position = [cpatch%radius*cos(phi), cpatch%radius*sin(phi), z]
      call cpatch%pre_init ()
      call cpatch%init_unique (same=.true.)
      if (i_r==1 .and. i_phi==1 .and. i_z==1) &
        call check_n_phi(cpatch%n(1), cpatch%ng(1), self%n_phi)
      !-----------------------------------------------------------------------
      ! Overwrite incorrectly initialised local and global coordinates
      !-----------------------------------------------------------------------
      cpatch%position = [cpatch%radius*cos(phi), cpatch%radius*sin(phi), z]
      cpatch%velocity = [0d0, self%omega*cpatch%radius, 0d0]
      cpatch%periodic = [.false., .true., .false.]
      cpatch%size = size
      cpatch%ds = size/cpatch%n
      cpatch%centre_nat = [cpatch%radius, phi, z]
      cpatch%llc_nat = cpatch%centre_nat-0.5d0*cpatch%size
      cpatch%llc_nat(2) = (cpatch%llc_nat(1)*phi - 0.5d0*cpatch%size(2)) / &
        cpatch%llc_nat(1)
      cpatch%centre_cart = cpatch%position
      cpatch%llc_cart = &
        [cpatch%llc_nat(1)*cos(cpatch%llc_nat(2)), &
        cpatch%llc_nat(1)*sin(cpatch%llc_nat(2)), cpatch%llc_nat(3)]
      !-------------------------------------------------------------------------
      ! Add static array pointers
      !-------------------------------------------------------------------------
      call static_arrays (self, cpatch)
      !-------------------------------------------------------------------------
      ! Add task link to task list and cylinder hash_table
      !-------------------------------------------------------------------------
      cpatch%link => halo
      cpatch%link%task => cpatch
      ptr => cpatch
      cpatch%key = [cpatch%level, i_r, i_phi-1, i_z]
      call self%task_list%append (cpatch%link)
      !if (verbose > 0) &
        !write(unit,'(a,i6,2x,4i3,l3)') &
         !'cylinder_t%init: id, key', cpatch%id, cpatch%key, associated(ptr)
      !-------------------------------------------------------------------------
      ! Set boundary bits
      !-------------------------------------------------------------------------
      if (i_r == 1) then
        call cpatch%boundaries%set (bits%xl)
      endif
      if (i_r == n_r) then
        call cpatch%boundaries%set (bits%xu)
      end if
      if (i_z == 1) then
        call cpatch%boundaries%set (bits%zl)
      endif
      if (i_z == n_z) then
        call cpatch%boundaries%set (bits%zu)
      end if
      !-------------------------------------------------------------------------
      ! Assign ranks
      !-------------------------------------------------------------------------
      cpatch%rank = (mpi%size-0.0001)*(i_phi-1.0)/(n_phi-1.0)
      !$omp atomic update
      n_task(1) = n_task(1) + 1
      if (cpatch%rank == mpi%rank) then
       !$omp atomic update
        n_task(2) = n_task(2) + 1
      end if
    end do
  end do
  call trace%end()
END SUBROUTINE init_cylinder

!===============================================================================
!> Create nbor lists, using the hash table to speed up
!===============================================================================
SUBROUTINE init_all_nbors (self)
  class(cylinders_t), target:: self
  class(link_t), pointer :: link
  integer:: n
  real(8):: time
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinders_t%init_all_nbors')
  call timer%tic (time)
  n = 0
  link => self%task_list%head
  do while (associated(link))
    if (mod(n,omp%nthreads) == omp%thread) then
      call self%init_nbors (link)
    end if
    n = n+1
    link => link%next
  end do
  call timer%toc ("time to initialize nbors", n, time)
  call trace_end
END SUBROUTINE init_all_nbors

!===============================================================================
!> Prepare an nbor-list for a cylinder task link.  Since the basic data type
!> (list_t) has an init_nbors(), this one will overload it, for task lists
!> of the cylinders_t type.
!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(cylinders_t)         :: self
  class(link_t),     pointer :: link
  class(cylinder_t), pointer :: my_cylinder, cylinder
  class(link_t),     pointer :: nblink
  class(task_t),     pointer :: task
  class(task_t),     pointer :: nbtask
  type(hash_table_t)         :: hash_table
  logical                    :: added
  logical,              save :: use_hash = .false.
  integer                    :: i_r, i1, i2, i3, j2, key(4), unit, tid
  integer,              save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinders_t%init_nbors', itimer=itimer)
  if (verbose > 4) then
    unit = stdout
  else if (verbose > 1) then
    unit = io_unit%mpi
  else
    unit = stdlog
  end if
  call link%remove_nbor_list (link%nbor)
  task => link%task
  tid = task%id
  select type (task)
  class is (patch_t)
  i_r = task%key(2)
  if (verbose > 1) &
    write (unit,'(a,i6,":",i2.2,3i4)') 'cylinder_t%init_nbors: target=', &
      tid, task%key
  !-----------------------------------------------------------------------------
  ! Initialize a hash table to protect agains double add
  !-----------------------------------------------------------------------------
  if (use_hash) call hash_table%init (100)
  n_nbors = 0
  my_cylinder => self%cylinder(i_r)
  do i1=-1,1
    select case (i1)
    case (-1)
      cylinder => my_cylinder%prev
    case ( 0)
      cylinder => my_cylinder
    case (+1)
      cylinder => my_cylinder%next
    end select
    if (associated(cylinder)) then
      if (shared%debug) &
        print *, 'checking cylinder:', i1, i1+key(2)
      do i2=-2,2
        !-----------------------------------------------------------------------
        ! Estimate how much further the nbor cylinder has rotated
        !-----------------------------------------------------------------------
        j2 = i2 + &
          nint((cylinder%omega - my_cylinder%omega)*task%time/dphi)
          ! (cylinder%omega - my_cylinder%omega)*task%time/dphi
        do i3=-2,2
          if (i1==0 .and. j2==0 .and. i3==0) cycle
          ! print *, '     checking height:', i3, i3+key(4)
          key = task%key + [0,i1,j2,i3]
          key(3) = modulo(key(3),n_phi)           ! wrap around
          call self%hash_table%get (key, nblink)
          added = .false.
          if (associated(nblink)) then
            !-------------------------------------------------------------------
            ! Check if the task is already added
            !-------------------------------------------------------------------
            if (use_hash) then
              call hash_table%get (key, nbtask)
            else
              nullify (nbtask)
            end if
            if (.not. associated(nbtask)) then
              nbtask => nblink%task
              !-----------------------------------------------------------------
              ! If not, add it, and add it to local hash table
              !-----------------------------------------------------------------
              if (use_hash) call hash_table%set (key, nbtask)
              select type (nbtask)
                class is (patch_t)
                if (nbtask%id > 0) then
                  if (overlaps_l (nbtask, task)) then
                    nblink => nbtask%link
                    if (verbose > 1) write (unit,'(i6,":",i2.2,3i4,a)') &
                      nbtask%id, key,' added as nbor'
                    call link%add_nbor (nblink)
                    !$omp atomic update
                    n_nbors = n_nbors + 1
                  else
                    added = .false.
                    if (verbose > 3) write (unit,'(i6,":",i2.2,3i4,a)') &
                      nbtask%id, key,' not added'
                  end if
                else
                  write (stderr,*) 'nbtask with non-positive ID:', nbtask%id
                end if
              end select
            else
              print *, 'skipping existing', nbtask%id
            end if
            if (shared%debug) &
              print *, '  azimuth:', i2, j2, key(3), nbtask%id, added, &
                merge('  ','no',added)//' overlap'
          end if
        end do
      end do
    end if
  end do
  if (use_hash) call hash_table%dealloc ()
  end select
  if (shared%debug) &
    print *, task%id, ' found', n_nbors
  call trace%end (itimer)
contains

!===============================================================================
!> Local definition, to ensure full control; the distance to the source must
!> be computed at the target%time, and hence cannot be used with distance_to()
!===============================================================================
logical function overlaps_l (source, target)
  class(patch_t):: source, target
  type(region_t):: authority
  real(8):: d(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) call trace%begin ('cylinder_t%overlaps_l', itimer=itimer)
  select type (target)
  class is (experiment_t)
    call target%get_position (source, d)
  end select
  overlaps_l = target%interest%overlaps (source%authority, d)
  if (timer%detailed) call trace%end (itimer)
end function overlaps_l

END SUBROUTINE init_nbors

!===============================================================================
! Domain covers n_h gas scale heights above and below mid-plane
!===============================================================================
FUNCTION calc_n_z (radius, size, n_h) RESULT (n_z)
  real(8), intent(in) :: radius, size
  real,    intent(in) :: n_h
  integer             :: n_z
  real(8)             :: m, v_k, r_au, t, csound, h
  !-----------------------------------------------------------------------------
  m = disk_model%m_star                         ! stellar mass
  v_k = sqrt (scaling%grav*m/radius)            ! Kepler speed
  r_au = radius*scaling%l/cgs%AU                ! radius in au
  t = disk_model%temperature (r_au)             ! temperature
  csound = disk_model%c_sound (r_au)            ! sound speed
  h = disk_model%h_scale (r_au)                 ! scale height
  n_z = ceiling (2d0*n_h*h/size)                ! number of patches vertically
END FUNCTION calc_n_z

!===============================================================================
!> Check if first and last patch in azimuthal dimension overlap at outer radial
!> edges, if so decrease dphi and dlnr to remove overlap
!===============================================================================
SUBROUTINE check_dphi (dlnr, dphi, n_phi)
  real(8), intent(inout) :: dlnr, dphi
  integer, intent(in)    :: n_phi
  !.............................................................................
  if (n_phi*dphi /= 2d0*math%pi) then
    dphi = 2d0*math%pi/n_phi
    dlnr = -log(-dphi+1d0)
    write (io%output,*) "To ensure that first and last patch in azimuthal &
    &dimension do not overlap at outer radial edges, dlnr is decreased to:", &
    dlnr
  end if
END SUBROUTINE check_dphi

!===============================================================================
!> Check if number of patches in azimuthal dimension is sufficiently large for
!> patches to overlap by at most one guard cell at inner radial edges
!===============================================================================
SUBROUTINE check_n_phi (n, ng, n_phi)
  integer, intent(in) :: n, ng, n_phi
  real(8)             :: n_phimin, d_phimax, dlnr_max
  !.............................................................................
  n_phimin = ceiling(2d0*math%pi/asin(1d0/(n+2d0*ng)))
  if (n_phi < n_phimin) then
    d_phimax = 2d0*math%pi/n_phimin
    dlnr_max = -log(1-d_phimax)
    write (io%output,*) "Patches overlapping by more than one guard &
    &cell at inner radial edges. Decrease dlnr to at most:", &
    dlnr_max
    call io%abort("cylinder%init")
  end if
END SUBROUTINE check_n_phi

!===============================================================================
!> Cache a few often used scalar and vector field, shared by all patches in the
!> same cylinder shell.  The first four are different for patches at different
!> vertical positions while the last two do not depend on vertical position.
!>         r: the exact radius
!>     r_hat: the unit vector in the exact radial direction
!>     r_inv: the exact inverse radius, times the radial unit vector
!>    f_grav: the force of gravity in the local coordinate system
!>   phi_hat: the unit vector in the phi direction
!>   v_patch: patch velocity relative to the global coordinate system
!===============================================================================
SUBROUTINE static_arrays (self, cpatch)
  class(cylinder_t):: self
  class(experiment_t):: cpatch
  integer:: m, i1, i2, i3, iz, j
  real(8):: x(3), phi_hat(3), p, r, rc, ri, g_times_mass, o, omega
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinder_t%static_arrays')
  call io%lock%set ('static_arrays')
  if (.not.associated(self%r)) then
    m  = cpatch%gn(1)
    g_times_mass = scaling%grav*disk_model%m_star
    allocate (self%r      (m,m,m,  self%n_z))
    allocate (self%r_hat  (m,m,m,3,self%n_z))
    allocate (self%r_inv  (m,m,m,3,self%n_z))
    allocate (self%f_grav (m,m,m,3,self%n_z))
    allocate (self%phi_hat(m,m,m,3))
    allocate (self%v_patch(m,m,m,3))
    o = (1+m)/2.0
    do i3=1,m
    do i2=1,m
    do i1=1,m
      x(1) = (i1-o)*cpatch%ds(1) + self%radius
      x(2) = (i2-o)*cpatch%ds(2)
      rc = sqrt(x(1)**2 + x(2)**2)
      phi_hat(1) = -x(2)/rc
      phi_hat(2) = +x(1)/rc
      phi_hat(3) = 0d0
      !-------------------------------------------------------------------------
      ! Compute static arrays in thin-disk approximation
      ! v_patch = the speed -- in the Eulerioan frame -- of the patch cells
      !-------------------------------------------------------------------------
      do j=1,3
        self%phi_hat(i1,i2,i3,j) = phi_hat(j)
        self%v_patch(i1,i2,i3,j) = phi_hat(j)*self%omega*rc
      end do
      do iz=1,self%n_z
        x(3) =  (i3-o)*cpatch%ds(3) + (iz-(1+self%n_z)*0.5d0)*cpatch%size(3)
        r = sqrt(rc**2 + x(3)**2)
        do j=1,3
          self%r    (i1,i2,i3,  iz) = r
          self%r_hat(i1,i2,i3,j,iz) = x(j)/r
          ri = 1d0/r
          self%f_grav(i1,i2,i3,j,iz) = -g_times_mass &
                                     * self%r_hat(i1,i2,i3,j,iz)*ri**2
          self%r_inv (i1,i2,i3,j,iz) = self%r_hat(i1,i2,i3,j,iz)*ri
        end do
      end do
    end do
    end do
    end do
  end if
  cpatch%r       => self%r
  cpatch%r_hat   => self%r_hat
  cpatch%phi_hat => self%phi_hat
  cpatch%f_grav  => self%f_grav
  cpatch%r_inv   => self%r_inv
  cpatch%v_patch => self%v_patch
  call io%lock%unset ('static_arrays')
  call trace%end ()
END SUBROUTINE static_arrays

END MODULE cylinder_mod
