!===============================================================================
!> Module to set up angles and bins for an RT patch with only absorption.
!> This means that the only extra memory with time slot that is needed is
!> one slot for the net radiation heating, \Sum_bins \Sum_angles Q_bin_angle
!>
!> The data structure consists of self, an extension of the gpatch_t data
!> type, which also contains n_omega clones of itself, with the specific info
!> for each ray direction.  Since the data type is an extension of patch_t, we
!> can use the standard guard zone loading procedures in download_mod or
!> data_hub_mod to fill the RT guard zones, and since the sub-tasks and main
!> task have the same data type, the sub-tasks can be cloned from the main task.
!>
!> Since directional information becomes available at different times, the
!> guard zone filling must take place for each angle, just before the solver
!> is called.
!>
!> The self%time, self%rt_timestep, and self%rt_next control the scheduling of
!> the tasks.  The MHD tasks require that the times of their nbors are "ahead",
!> allowing for a grace fraction of the nbor time step.  The RT sub-tasks
!> advance their time by self%rt_timestep (effectively their self%dtime) when
!> done, hence they must update when their MHD master task passes their own time,
!> which conforms with standard nbor relations.   To clarify, with the times
!>
!> MHD time: mhd%time
!>  RT time: mhd%rt%time
!> ray time: mhd%rt%omega%time
!>
!> The RT_time is a "rabbit time", which gives the time until which the MHD task
!> should still be happy with the RT solution.  When the parent MHD task passes
!> that time, the RT sub-tasks become out-of-date, and they perform solutions.
!> Afterwards, they advance their "rabbit time" again, and when all of them have
!> done so, the main RT task becomes ready to update, and when it runs it sums
!> up the sub-task contributions, and advances its own "rabbit time", which then
!> allows the MHD task to update.   Clearly, the EOS call that gives values to
!> the RT tasks at the new MHD time has to be made after the MHD task has passed
!> by, but before the first sub-task executes.  This could be detected either in
!> the MHD task, or in the RT sub-tasks, where the first detection should lead
!> to locking the EOS while executing it, followed by another test (so a double
!> critical region type of thing).
!>
!> Initially, we wish that only the RT states are ready to update, and we set
!> things up to compute the EOS, solve the RT, sum up Q, and only then let the
!> MHD update.  This means that initially, the EOS_time should be negative, the
!> RT times should be RT_timestep/2 negative, and the MHD_time should be zero.
!===============================================================================
MODULE rt_mod
  USE rt_boundaries_mod
  USE rt_integral_mod
  USE radau_mod
  USE boundaries_mod
  USE eos_mod
  USE list_mod
  USE gpatch_mod
  USE link_mod
  USE task_mod
  USE index_mod
  USE mesh_mod
  USE io_mod
  USE aux_mod
  USE timer_mod
  USE scaling_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE math_mod
  USE bits_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from gpatch, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: rt_t
    class(gpatch_t), pointer :: patch                 ! pointer to MHD patch
    type(rt_boundaries_t)    :: rt_boundaries         ! RT boundary conditions
    integer                  :: n_warmup=5            ! number of initial solutions
    integer                  :: n_bin=4               ! number of frequency bins
    integer                  :: n_omega=0             ! number of space-angles
    integer                  :: i_omega=0             ! space-angle index
    integer                  :: axis=0                ! neares axis
    real, pointer            :: w_bin(:)              ! frequency weight
    type(rt_t), pointer      :: omega(:)  => null()   ! different RT directions
    type(rt_t), pointer      :: main => null()        ! point to main RT task
    real                     :: cdtd = 1.0
    real                     :: mu, rt_phi, w         ! mu, phi, weights
    real                     :: rt_hat(3)             ! RT direction
    real                     :: epsilon =1.0          ! scattering coefficient
    real(8)                  :: eos_time              ! laste EOS time
    real                     :: rt_grace=0.01         ! deal with roundoff
    real(8)                  :: rt_llc(3)=0d0         ! lower left corner
    real(8)                  :: rt_urc(3)=0d0         ! upper right corner
    real(8)                  :: rt_timestep=0.01_8    ! constant time step
    real(8)                  :: warmup_end_time=0d0   ! time when RT time should start varying
    integer                  :: n_mu=3, n_phi=4       ! no. of angles
    logical                  :: on=.false.            ! allows local turn-off
    logical                  :: vertical_propagation=.false. ! see code
    logical                  :: warmup_done=.false.
    real, pointer            :: q(:,:,:,:)            ! net radiative heating
    real, pointer            :: tt(:,:,:)             ! EOS temperature
    real, pointer            :: rk(:,:,:,:)           ! EOS rho*kappa
    real, pointer            :: src(:,:,:,:)          ! EOS source function
    real, pointer            :: qr(:,:,:)             ! net radiative heating
  contains
    procedure:: pre_init                              ! read parameters
    procedure, nopass:: init_task                     ! init RT task
    procedure:: add_to_task_list                      ! add to task list
    procedure:: init_omega                            ! init omega task
    procedure:: needs                                 ! check upstream
    procedure:: is_upstream_of                        ! check upstream
    procedure:: pre_update                            ! pre MHD
    procedure:: post_update                           ! post MHD
    procedure:: update                                ! RT task update
    procedure:: rt_eos                                ! EOS for RT
    procedure:: solve                                 ! RT solver
    procedure:: setup_boundaries                      ! BCs
    procedure:: task_info                             ! formatted info print
    procedure:: nbor_info                             ! formatted info print
    procedure:: courant_condition                     ! amend Courant condition
    procedure:: output                                ! intercept output
  end type
  integer, save:: verbose=0                           ! verbosity
  type(rt_t), public:: rt                             ! public instance
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, link)
  class (rt_t):: self
  class (link_t), pointer :: link
  !-----------------------------------------------------------------------------
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n
  integer, save:: n_mu=3, n_phi=4, n_bin=4, nt=3, n_warmup=3, ng=1
  real, save   :: courant=0.5, cdtd=1.0
  real,    save:: rt_grace=0.01
  real(8), save:: rt_llc(3)=-huge(1d0)
  real(8), save:: rt_urc(3)=+huge(1d0)
  real(8), save:: rt_timestep=0.01_8
  logical, save:: vertical_propagation=.false.
  logical, save:: first_time=.true.
  logical, save:: on=.true.
  !-----------------------------------------------------------------------------
  namelist /sc_rt_params/ on, verbose, nt, ng, rt_llc, rt_urc, rt_timestep, &
    n_mu, n_phi, n_bin, rt_grace, vertical_propagation, n_warmup, courant, cdtd
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%pre_init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, sc_rt_params)
    if (io%master) write (io%output, sc_rt_params)
  end if
  !$omp end critical (input_cr)
  self%patch => gpatch%cast2gpatch(link%task)
  self%on = on
  self%nt = nt
  self%ng = ng
  self%n_warmup = n_warmup
  self%n_bin = n_bin
  self%rt_grace = rt_grace
  self%rt_llc = rt_llc
  self%rt_urc = rt_urc
  self%rt_timestep = rt_timestep
  self%vertical_propagation = vertical_propagation
  self%n_mu  = n_mu
  self%n_phi = n_phi
  self%courant = courant
  self%cdtd    = cdtd
  call rt_integral%init
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the rt_t data types.  To setup for n_warmup initial solutions, we
!> set the time rt_t times = -n_warmup * rt_t%timestep, making sure to call
!> rt_t%rt_eos for the initial MHD state.
!===============================================================================
SUBROUTINE init_task (self)
  class (rt_t), pointer:: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
  class (rt_t),     pointer :: omega
  class (link_t),   pointer :: link
  class (*),        pointer :: dummy
  integer:: n_omega, i_omega, dir, i_mu, i_phi, iostat, n, l(3), u(3), id_mpi
  real, dimension(:), allocatable:: mu, w_mu
  ! ----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_task')
  p => self%patch
  !-----------------------------------------------------------------------------
  ! interconnect itself
  !-----------------------------------------------------------------------------
  nullify(self%link)
  allocate(link)
  link%task => self
  self%link => link
  self%track   = .false.
  self%type    = 'rt_t'
  self%kind    = 'rt_sc'
  self%box     = p%box
  !-----------------------------------------------------------------------------
  ! Number of ray directions
  !-----------------------------------------------------------------------------
  n_omega = merge (1+self%n_phi*(self%n_mu-1), 0, self%n_mu > 0)
  n_omega = 2*n_omega
  self%n_omega = n_omega
  !-----------------------------------------------------------------------------
  ! Set the RT task ID in a predictable way, spacing by (1+n_omega) behind the
  ! largest MPI-wide ID used so far, to make room for the omega tasks
  !-----------------------------------------------------------------------------
  id_mpi = mpi%id%update(0)
  self%id = id_mpi + (p%id-1)*(1+n_omega) + 1
  if (verbose == -1) &
    write (io_unit%log,*) 'rt_t%init: id, mhd%id, mpi%id =', self%id, p%id, id_mpi
  call self%patch_t%setup (position=p%position, size=p%size, n=p%n, nt=self%nt,&
                           ng=self%ng, nv=self%n_bin, nw=1)
  !-----------------------------------------------------------------------------
  ! Synchronize with MHD patch
  ! ----------------------------------------------------------------------------
  self%rank       = p%rank
  self%level      = p%level
  self%size       = p%size
  self%origin     = p%origin
  self%llc_cart   = p%llc_cart
  self%llc_nat    = p%llc_nat
  self%position   = p%position
  self%centre_nat = p%position
  self%restart    = p%restart
  self%warmup_end_time = p%time
  self%iit      = 1
  self%mu       = 1.0
  self%rt_phi   = 0.0
  self%t(self%it) = self%time
  ! -- add an `eps` to avoid round-off errors
  self%out_next = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time
  self%dt_fixed = self%rt_timestep
  self%dtime    = self%rt_timestep
  self%dt       = self%rt_timestep
  self%grace    = self%rt_grace
  self%time     = p%time - (self%n_warmup-self%grace)*self%dtime
  self%eos_time = self%time - self%rt_timestep
  self%iout     = self%restart+1
  self%periodic = p%periodic
  self%verbose  = verbose
  call self%set (bits%no_density)
  !-----------------------------------------------------------------------------
  ! Check if this MHD patch needs RT
  !-----------------------------------------------------------------------------
  call self%setup_boundaries()
  if (self%on) then
    call self%rt_boundaries%init (self%mesh)
    !---------------------------------------------------------------------------
    ! Allocate heating for MHD
    !---------------------------------------------------------------------------
    if (.not.allocated(self%patch%heating_per_unit_volume)) then
      allocate(self%patch%heating_per_unit_volume(p%gn(1), p%gn(2), p%gn(3)))
      self%patch%heating_per_unit_volume = 0.0
    end if
    !---------------------------------------------------------------------------
    ! Integration weights (w_bin is baked into the source function in this case)
    !---------------------------------------------------------------------------
    allocate (self%w_bin (self%n_bin))
    self%w_bin = 1.0
    !---------------------------------------------------------------------------
    ! Assign angles and integration weight and allocate the required RT solvers
    ! for all angles
    !---------------------------------------------------------------------------
    if (self%n_mu > 0) then
      allocate (mu(self%n_mu),w_mu(self%n_mu))
      call radau%init (self%n_mu, mu, w_mu)
    end if
    !---------------------------------------------------------------------------
    ! Allocate data, belonging to the main task
    !---------------------------------------------------------------------------
    allocate (self%rk (self%gn(1),self%gn(2),self%gn(3),self%n_bin))
    allocate (self%tt (self%gn(1),self%gn(2),self%gn(3)))
    allocate (self%qr (self%gn(1),self%gn(2),self%gn(3)))
    allocate (self%src(self%gn(1),self%gn(2),self%gn(3),self%n_bin))
    call io%bits_mem(storage_size(self%tt  ),product(shape(self%tt )),'tt')
    call io%bits_mem(storage_size(self%rk  ),product(shape(self%rk )),'rk')
    call io%bits_mem(storage_size(self%qr  ),product(shape(self%qr )),'qr')
    call io%bits_mem(storage_size(self%src ),product(shape(self%src)),'src')
    !---------------------------------------------------------------------------
    ! Allocate n_omega copies of self, cloning self
    !---------------------------------------------------------------------------
    if (n_omega > 0) allocate (self%omega(n_omega), source=self)
    self%i_omega = 0
    self%main => self%omega(1)%main
    !---------------------------------------------------------------------------
    ! Assign angles and initialize clones
    !----------------------------------------------------------------------------
    i_omega = 0
    do i_mu=1,self%n_mu
      n = self%n_phi
      if (i_mu == self%n_mu) n = 1
      do i_phi=1,n
        do dir=-1,1,2
          i_omega = i_omega+1
          omega => self%omega(i_omega)
          omega%mu  = dir*mu(i_mu)
          omega%w   = w_mu(i_mu)*math%pi4/n
          omega%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + &
            (dir+1)/2.0, 2.0) - 1.0)
          if (i_mu == self%n_mu) omega%rt_phi = 0.0
          call self%init_omega(i_omega)
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    ! Check task angles
    !---------------------------------------------------------------------------
    if (verbose > 0 .and. self%id == 1) then
      do i_omega=1,n_omega
        omega => self%omega(i_omega)
        if (verbose > 0) &
          write(io%output,'(a,3i7,f8.3)') &
            'rt_t%init: i_omega, elevation, phi, mu =', i_omega, &
            nint(asin(omega%mu)*180./math%pi), nint(omega%rt_phi*180./math%pi), &
            omega%mu
      end do
    end if
    !---------------------------------------------------------------------------
    ! Point mem to a patch mem section -- download_mod should be able handle a
    ! target with fewer guard cells than the source
    !---------------------------------------------------------------------------
    l = p%mesh%li - self%mesh%ng
    u = p%mesh%ui + self%mesh%ng
    self%mem  => p%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,:,:)
  else
    !---------------------------------------------------------------------------
    ! Deallocate the self%mem for this task, but keep the task, for simplicity
    !---------------------------------------------------------------------------
    call self%dealloc()
    self%n_omega = 0
  end if
  call trace%end()
END SUBROUTINE init_task

!===============================================================================
!> Add the (1+n_omega) RT tasks to the task list
!===============================================================================
SUBROUTINE add_to_task_list (self)
  class(rt_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  integer:: i_omega
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%add_tasks_to')
  !-----------------------------------------------------------------------------
  task_list => self%patch%task_list
  if (self%is_clear(bits%no_rt)) then
    call task_list%prepend_link (self%link)
    do i_omega=1,self%n_omega
      call task_list%prepend_link (self%omega(i_omega)%link)
    end do
  end if
  call trace%end()
END SUBROUTINE add_to_task_list

!===============================================================================
!> Initialize rt sub-tasks; note that all self values are initially cloned
!===============================================================================
SUBROUTINE init_omega (self, i_omega)
  class (rt_t), target:: self
  integer:: i_omega
  ! ............................................................................
  class (rt_t), pointer:: omega
  real:: mu(3)
  integer:: loc(1)
  ! ----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_omega')
  omega => self%omega(i_omega)
  omega%main => self
  ! -- assign unique and globally valid ID
  omega%id = self%id + i_omega
  call omega%task_t%init
  if (verbose == -1) &
    write (io_unit%log,*) 'rt_t%init_omega: id, self%id, mhd%id =', omega%id, &
      self%id, self%patch%id
  ! -- set values that should be different from self
  omega%i_omega = i_omega
  omega%it      = max(omega%nt-1,1)
  omega%new     = omega%nt
  omega%nw      = 1
  omega%nv      = self%n_bin
  ! -- allocate unsigned and pervolume
  allocate (omega%unsigned(omega%nv), omega%pervolume(omega%nv))
  omega%unsigned  = .false.
  omega%pervolume = .false.
  ! -- allocate time slots
  allocate (omega%t(omega%nt), omega%dt(omega%nt), omega%iit(omega%nt))
  !-----------------------------------------------------------------------------
  ! Setting omega%time beyond self%time to trigger first call to self%update(),
  ! which is needed to compute the EOS data
  !-----------------------------------------------------------------------------
  omega%time    = self%time + 2.0*self%grace*self%dtime
  omega%t       = omega%time
  omega%dt      = omega%dtime
  omega%iit     = 1
  allocate (omega%mem(self%gn(1),self%gn(2),self%gn(3),self%n_bin,self%nt,1))
  call io%bits_mem(storage_size(omega%mem)   ,product(shape(omega%mem))   ,'q')
  omega%q => omega%mem(:,:,:,:,omega%it,1)
  omega%mem = 0.0
  !-----------------------------------------------------------------------------
  ! Determine closest axis direction
  !-----------------------------------------------------------------------------
  omega%rt_hat(1) = sin(acos(omega%mu))*cos(omega%rt_phi)
  omega%rt_hat(2) = sin(acos(omega%mu))*sin(omega%rt_phi)
  omega%rt_hat(3) = omega%mu
  loc = maxloc(abs(omega%rt_hat))
  omega%axis = loc(1)
  !-----------------------------------------------------------------------------
  ! Allocate task link and set link%task 
  !-----------------------------------------------------------------------------
  allocate (omega%link)
  omega%link%task => omega
  call trace%end()
END SUBROUTINE init_omega

!===============================================================================
!> Determine if a task is needed (is upstream).  We are being presented with
!> the nbors of the MHD task, so there is no task above the top-most task, but
!> there may well be tasks below the rt_bottom, which are not needed as RT nbors
!===============================================================================
FUNCTION needs (self, omega, task)
  logical:: needs
  class(rt_t):: self, omega
  class(rt_t):: task
  !.............................................................................
  real:: ds1(3), ds2(3), sinth
  real, parameter:: eps=1e-4
  integer:: is1(3), is2(3), i, lbit(3), ubit(3)
  logical:: lower, upper, upstream
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and omega points inward, or the task has no RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  ds1 = self%position-task%position
  where (self%periodic)
    ds1 = modulo(ds1 + 0.5_8*self%mesh%b, self%mesh%b) - 0.5_8*self%mesh%b
  end where
  ds1 = ds1/self%box
  sinth = merge(0.0, sqrt(1.0-omega%mu**2), abs(omega%mu)==1.0)
  ds2 = [sinth*cos(omega%rt_phi),sinth*sin(omega%rt_phi),omega%mu]
  is1 = 0
  is2 = 0
  do i=1,3
    if (ds1(i) > +eps) is1(i)=1
    if (ds2(i) > +eps) is2(i)=1
    if (ds1(i) < -eps) is1(i)=-1
    if (ds2(i) < -eps) is2(i)=-1
  end do
  needs = .true.
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%no_rt)) then
    needs = .false.
    return
  end if
  !-----------------------------------------------------------------------------
  ! If the nbor lies outside the boundaries it is not needed
  !-----------------------------------------------------------------------------
  lbit = [bits%xl, bits%yl, bits%zl]
  ubit = [bits%xu, bits%yu, bits%zu]
  do i=1,3
    lower = self%boundaries%is_set(lbit(i))
    upper = self%boundaries%is_set(ubit(i))
    if (upper .and. ds1(i)==1) then
      needs = .false.
    end if
    if (lower .and. ds1(i)==-1) then
      needs = .false.
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Vertical rays depend only on exact alignment
  !-----------------------------------------------------------------------------
  if (abs(omega%mu)==1.0) then
    needs = all(is1==is2)
  !-----------------------------------------------------------------------------
  ! For inclined rays, define upstream as "not downstream in any direction"
  !-----------------------------------------------------------------------------
  else if (any(is1*is2 < 0)) then
    needs = .false.
  !-----------------------------------------------------------------------------
  ! If a horizontal ray component is axis aligned, exclude non-aligned nbors
  !-----------------------------------------------------------------------------
  else if (is2(1)==0 .and. is1(1)/=0) then
    needs = .false.
  else if (is2(2)==0 .and. is1(2)/=0) then
    needs = .false.
  end if
  upstream = task%is_upstream_of(omega)
if (needs .neqv. upstream) then
  print '(a,2i4,2f8.3,l4,2(3f8.3,2x),2i4)', &
    'needs and is_upstream of disagree on', omega%id, task%id, omega%mu, &
    omega%rt_phi, needs, real(omega%task_t%distance(task),kind=4), &
    omega%mesh%s, omega%i_omega, task%i_omega
  !needs = upstream
end if
END FUNCTION needs

!===============================================================================
!> Determine if the self task is "upstream" of the other task; i.e., if its
!> guard cells that lie between the two patches overlap with the internal
!> region of the other patch
!>
!> In the most general case one needs to first identify the "faces" of the task
!> that are at play.  Those correspond to the components of the ray direction
!> that are non-zero.  For each upstream candidate, one should run through the
!> candidate faces, identified by squares (or more generally flat volumes), and
!> check if they have overlap in all three directions with the internal region
!> of the candidate upstream task.
!===============================================================================
LOGICAL FUNCTION is_upstream_of (self, other)
  class(rt_t):: self, other
  real(8):: p(3), d(3), dh(3)
  integer:: i
  logical:: debug
  !-----------------------------------------------------------------------------
  ! Check if the direction vector d is in the same octant; to be correct for
  ! periodic meshes this must use the patch_t%distance() function!
  !-----------------------------------------------------------------------------
  d = other%distance (self)                       ! direction towards other
  dh = d*self%rt_hat                              ! projection onto rt_hat
  is_upstream_of = all(dh > -self%ds)             ! all components "not negative"
debug=.false.
!debug = self%position(2)==other%position(2) !&
  !.and. abs(self%position(2)-0.5d0) < 0.01
if (debug) then
print *, self%id, other%id
print 1, 'sp:',  self%position      , 'sm:',  self%mesh%p
print 1, 'op:', other%position      , 'om:', other%mesh%p
print 1, ' d:', d
print 1, 'dh:', dh                  , 'rt:', self%rt_hat, is_upstream_of
1 format(2(a,3f8.3,2x),l3)
end if
  !-----------------------------------------------------------------------------
  ! If we have a possible upstream task
  !-----------------------------------------------------------------------------
  if (is_upstream_of) then
    p = self%position
    do i=1,3
      associate (r => self%mesh(i)%r)
      !-------------------------------------------------------------------------
      ! If the direction vector projected on the ray is non-zero (w roundoff),
      ! then move the point to the guard zone in the other task direction
      !-------------------------------------------------------------------------
      if (abs(dh(i)) > self%ds(i)) then
        if      (d(i) > +self%ds(i)) then
          p(i) = p(i) + r(self%mesh(i)%uo)
        else if (d(i) < -self%ds(i)) then
          p(i) = p(i) + r(self%mesh(i)%lo)
        end if
      end if
      end associate
    end do
    !---------------------------------------------------------------------------
    ! If that doesn't overlap with the other task, self is not upstream
    !---------------------------------------------------------------------------
    if (.not. other%overlaps_point (p)) then
      is_upstream_of = .false.
    end if
  end if
if (debug) &
print '(3i4,l4,5(2x,a,3f8.3))', self%id, other%id, self%axis, is_upstream_of, &
'sp:',self%mesh%p, 'op:',other%mesh%p, 'd:',abs(d), &
's2:',other%mesh%s/2d0
END FUNCTION is_upstream_of

!===============================================================================
!> Check if patch belongs to an RT box or is at its boundary.  The mesh flags
!> lower_boundary and upper_boundary are set for the tasks that contain bndries,
!> and self%on is set false for tasks that are entirely outside.
!===============================================================================
SUBROUTINE setup_boundaries(self)
  class(rt_t):: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3)
  integer :: lbit(3), ubit(3), i
  ! ............................................................................
  call self%init_bdries
  llc = self%position-0.5_8*self%size
  urc = self%position+0.5_8*self%size
  !
  lbit = [bits%xl, bits%yl, bits%zl]        ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]        ! upper boundary bits -- convenience
  !
  do i = 1,3
    if ((urc(i) < self%rt_llc(i))) then
      call self%patch%set(bits%no_rt)
      self%on = .false.
    else if (llc(i) <= self%rt_llc(i)) then
      call self%boundaries%set (lbit(i))
      self%mesh(i)%lower_boundary = .true.
    end if
    !
    if ((llc(i) > self%rt_urc(i))) then
      call self%patch%set(bits%no_rt)
      self%on = .false.
    else if (urc(i) >= self%rt_urc(i)) then
      call self%boundaries%set (ubit(i))
      self%mesh(i)%upper_boundary = .true.
    end if
  end do
END SUBROUTINE setup_boundaries

!===============================================================================
!> This is called only for MHD tasks, just before they update, with self = main,
!> so this is suitable place to add the resulting net radiative heating to the
!> heating per unit volume.
!===============================================================================
SUBROUTINE pre_update (self)
  class(rt_t):: self
  !.............................................................................
  class(rt_t), pointer:: omega
  integer:: l(3), u(3), i_omega, i_bin
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) &
    return
  call trace%begin ('rt_t%pre_update', itimer=itimer)
  !---------------------------------------------------------------------------
  ! -- Sum Q over angles and bin
  self%qr = 0.0
  ! -- multiply with the stored net intensities
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    do i_bin=1,self%n_bin
      self%qr = self%qr + self%rk(:,:,:,i_bin)*self%w_bin(i_bin)*omega%w* &
                omega%mem(:,:,:,i_bin,omega%it,1)
    end do
  end do
  call self%patch%aux%register ('qr', self%qr)
  !---------------------------------------------------------------------------
  l = self%patch%mesh%li - self%mesh%ng
  u = self%patch%mesh%ui + self%mesh%ng
  self%patch%heating_per_unit_volume(l(1):u(1),l(2):u(2),l(3):u(3)) = &
  self%patch%heating_per_unit_volume(l(1):u(1),l(2):u(2),l(3):u(3)) + self%qr
  !---------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> This is called from MHD tasks, after updates, with self = main RT task.
!> Sets the time step of the RT tasks equal to the MHD timestep, which was
!> recently computed, and collects u_max and modifies the mhd%u_max
!===============================================================================
SUBROUTINE post_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(gpatch_t), pointer:: mhd
  class(rt_t), pointer:: omega
  integer:: i
  !-----------------------------------------------------------------------------
  mhd => self%patch
  if (verbose == -1) then
    !$omp atomic
    timer%n_mhd = timer%n_mhd+1
    write (io_unit%queue,*) 'mhd', timer%n_mhd, wallclock(), &
      self%patch%id, self%id, mhd%time+mhd%dtime, self%on, self%rotated
    if (timer%n_mhd > 1000) then
      !$omp atomic write
      timer%n_mhd = 0
      !$omp atomic write
      timer%n_solve = 0
    end if
  end if
  if (.not.self%on) return
  call trace%begin ('rt_t%post_update')
  !-----------------------------------------------------------------------------
  ! We now know the timestep size and the time, and we can set these for all the
  ! RT solvers, in such a way that they update in the desired order.
  !-----------------------------------------------------------------------------
  self%dtime = mhd%dtime
  self%time  = mhd%time + self%grace*self%dtime
  do i=1,self%n_omega
    omega => self%omega(i)
    omega%time = mhd%time
    !-------------------------------------------------------------------------
    ! Inclined tasks may advance their time a factor mu slower, and thus forced
    ! to update more often, to propagate the solution by about the same amount
    ! vertically.  NOTE (FIXME):  Purely horizontal rays need special treatment
    !-------------------------------------------------------------------------
    if (self%vertical_propagation) then
      omega%dtime = mhd%dtime*abs(omega%mu)
    else
      omega%dtime = mhd%dtime + 2*self%grace*self%dtime
    end if
  end do
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE post_update

!===============================================================================
!> The rt_t%update procedure gets called directly by the task dispatcher, when
!> it has concluded that the task is ready to update.  An update of the main RT
!> task consists of summing up the results from the RT sub-tasks, while the
!> sub-task updates perform the actual RT solutions along rays.
!===============================================================================
SUBROUTINE update (self)
  class(rt_t):: self
  !.............................................................................
  class(rt_t), pointer :: omega
  integer:: i_omega, i_bin
  real, dimension(:,:,:), pointer:: q
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) &
    return
  call trace%begin ('rt_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! The omega0 task is responsible for calling rt_t%rt_eos to compute the RT EOS
  ! data. The nbor list has triggered download of nbor MHD data, before this call
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) then
    call self%rt_eos (self%time+self%dtime)
  !-----------------------------------------------------------------------------
  ! RT sub-task, solving for the radiation field in a given direction
  !-----------------------------------------------------------------------------
  else
    ! -- Call RT solver, storing result in %new slot, at time = t(%it)
    self%q => self%mem(:,:,:,:,self%new,1)
    call self%solve ()
    if (verbose == -1) then
      !$omp atomic
      timer%n_solve = timer%n_solve+1
      write (io_unit%queue,*) 'ome', timer%n_solve, wallclock(), &
        self%patch%id, self%id, self%time+self%dtime, self%i_omega
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Compute the RT EOS quantities for the patch, using self%eos_time and an OMP
!> critical region to guard against duplicate calculations when multi-threading.
!===============================================================================
SUBROUTINE rt_eos (self, time)
  USE units_mod
  class (rt_t):: self
  real(8):: time
  !.............................................................................
  class (gpatch_t), pointer :: patch
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e, pgas, fmax
  !-----------------------------------------------------------------------------
  if (time > self%eos_time) then
    !$omp critical (eos_cr)
    if (time > self%eos_time) then
      call omp_eos
      self%eos_time = time
      if (verbose == -1) &
        write (io_unit%queue,*) 'eos', timer%n_mhd, wallclock(), &
          self%patch%id, self%id, time
    else if (verbose == -1) then
        write (io_unit%queue,*) 'eo0', timer%n_mhd, wallclock(), &
          self%patch%id, self%id, time, self%eos_time
    end if
    !$omp end critical (eos_cr)
  end if
  return
contains
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine omp_eos
  integer:: m(4), ib
  real:: u_max, stefan
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  patch => self%patch
  m = shape(self%rk)
  allocate (pgas(m(1),m(2),m(3)),fmax(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,patch%idx%d,patch%it,1)
  e => self%mem(:,:,:,patch%idx%e,patch%it,1)
  if (verbose > 2) then
    call self%stats (d, 'd')
    call self%stats (e, 'e')
  end if
  call eos%lookup_table (m, d=d, e=e, rk=self%rk, src=self%src, tt=self%tt, pg=pgas)
  ! ----------------------------------------------------------------------------
  ! Compute radiative diffusion speed, and set default omega%u_max
  ! ----------------------------------------------------------------------------
  self%u_max = 0.0
  stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
  do ib=1, eos%n_lambda
    fmax = math%pi2*self%src(:,:,:,ib)/pgas*(self%rk(:,:,:,ib)*minval(self%ds))
    !fmax = (pgas/d)**4/pgas*(self%rk(:,:,:,ib)*minval(self%ds))
    fmax = fmax/(1.+(self%rk(:,:,:,ib)*minval(self%ds))**2)
    u_max = 16./3.*self%fmaxval(fmax)*self%courant/self%cdtd*self%epsilon
    self%u_max = max(self%u_max,u_max)
  end do
  self%omega(:)%u_max =  self%u_max
  if (verbose > 2) then
    call self%stats (self%tt  , 'tt')
    call self%stats (pgas, 'pgas')
    call self%stats (self%rk(:,:,:,1), 'rkap')
    call self%stats (fmax, 'fmax')
  end if
  ! ----------------------------------------------------------------------------
  deallocate (fmax, pgas)
  if (.not.self%warmup_done) then
    if (self%time > self%warmup_end_time) then
      self%warmup_done=.true.
      self%dt_fixed = -1.0
    end if
  end if
  call trace%end (itimer)
end subroutine omp_eos
END SUBROUTINE rt_eos

!===============================================================================
!> Call integral short characteristics RT solver.
!===============================================================================
SUBROUTINE solve (self)
  class(rt_t):: self
  !.............................................................................
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) return
  call trace%begin ('rt_t%solve', itimer=itimer)
  !-----------------------------------------------------------------------------
  call self%rt_boundaries%condition (mem=self%q, src=self%src, rk=self%rk, &
                                     tt=self%tt, mu=self%rt_hat)
  call rt_integral%solve (self%rk, self%src, self%q, self%mesh, self%axis, &
                          self%rt_hat)
  call self%patch%aux%register ('rkap', self%rk )
  call self%patch%aux%register ('src',  self%src)
  call trace%end (itimer)
END SUBROUTINE solve

!===============================================================================
!> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
!===============================================================================
SUBROUTINE task_info (self)
  class(rt_t):: self
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(1x,a,12x,i6,2i4,7x,2l1,2x,3f9.3,i4,2f7.3)') &
    'rt_t%task_info: id, rank, level, BV, pos, i_omega, mu, phi =', &
    self%id, self%rank, self%level, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), self%position, &
    self%i_omega, self%mu, self%rt_phi
END SUBROUTINE task_info

!===============================================================================
!> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
!===============================================================================
SUBROUTINE nbor_info (self, nbor)
  class(rt_t):: self
  class(link_t):: nbor
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(3x,a,i6,2i4,2x,3l1,2x,2l1,2x,3f9.3,i4,2f7.3)') &
    'rt_t%nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
    self%id, self%rank, self%level, nbor%needed, nbor%needs_me, nbor%download, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), &
    self%position, self%i_omega, self%mu, self%rt_phi
END SUBROUTINE nbor_info

!===============================================================================
!> Intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(rt_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%courant_condition')
  self%patch%u_max = max (self%patch%u_max, self%u_max)
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Auxiliary output
!===============================================================================
SUBROUTINE output (self)
  class(rt_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE output

END MODULE rt_mod
