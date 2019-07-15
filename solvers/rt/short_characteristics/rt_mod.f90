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
  USE boundaries_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE index_mod
  USE gpatch_mod
  USE radau_mod
  USE link_mod
  USE eos_mod
  USE bits_mod
  USE math_mod
  USE task_mod
  USE rt_integral_mod
  USE rt_boundaries_mod
  USE scaling_mod
  USE mesh_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from gpatch, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: rt_t
    class(gpatch_t), pointer   :: patch                                         ! pointer to MHD patch
    type(rt_boundaries_t)      :: rt_boundaries
    integer                    :: n_lambda = 4                                  ! number of frequency bins
    integer                    :: n_omega=0                                     ! number of rt space-angles
    integer                    :: i_omega=0                                     ! space-angle ID for particular task
    real, pointer              :: w_bin(:)                                      ! frequency bin statistical weight
    type(rt_t), pointer        :: omega(:)                                      ! different RT directions
    real                       :: mu, rt_phi, w
    real(8)                    :: eos_time
    real   :: rt_grace=0.0
    real(8):: rt_llc(3)    = 0d0
    real(8):: rt_extent(3) = 0d0
    real(8):: rt_timestep=0.01_8
    integer:: n_mu=3, n_phi=4
    logical:: on = .false.
    logical:: vertical_propagation=.false.
    logical:: constant_q = .true.
    real, pointer              :: q(:,:,:,:)
    real, pointer              :: tt(:,:,:)
    real, pointer              :: rk(:,:,:,:)
    real, pointer              :: src(:,:,:,:)
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_omega
    procedure:: init_nbors
    procedure:: needs
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
    procedure:: rt_eos
    procedure:: solve
    procedure:: setup_boundaries
    procedure:: warm_up
    procedure:: nbor_relations
    procedure:: task_info
    procedure:: nbor_info
  end type
  integer, save:: verbose=0
  real(8), save:: dbg_position(2) = -1e30
  type(rt_t), public:: rt
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, link)
  class (rt_t):: self
  class (link_t), pointer :: link
  !-----------------------------------------------------------------------------
  integer, save:: n_mu=3, n_phi=4, n_lambda=4, nt=3
  real,    save:: rt_grace=0.0
  real(8), save:: rt_llc(3)    = 0d0
  real(8), save:: rt_extent(3) = 0d0
  real(8), save:: rt_timestep=0.01_8
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n
  logical, save:: vertical_propagation=.false.
  logical, save:: constant_q = .true.
  logical, save:: first_time=.true.
  logical, save:: on =.false.
  !-----------------------------------------------------------------------------
  namelist /sc_rt_params/ on, nt, rt_llc, rt_extent, rt_timestep, n_mu, n_lambda, n_phi, verbose, &
                       rt_grace, vertical_propagation, constant_q
                       
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
  if (on) then
    if (.not. io%master) verbose=0
    self%on = on
    self%nt = nt
    self%n_lambda = n_lambda
    self%rt_grace = rt_grace
    self%rt_llc = rt_llc
    self%rt_extent = rt_extent
    self%rt_timestep = rt_timestep
    self%vertical_propagation = vertical_propagation
    self%constant_q = constant_q
    self%n_mu  = n_mu
    self%n_phi = n_phi
    !---------------------------------------------------------------------------
    ! Request additional variable slot indices
    !---------------------------------------------------------------------------
    if (self%patch%idx%qr <1) then
      call self%patch%idx%request_index (self%patch%idx%qr, self%patch%nv)
    end if
  end if
  if (self%patch%idx%tt <1) then
    call self%patch%idx%request_index (self%patch%idx%tt, self%patch%nv)
  end if
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the rt_t data types.  To setup for n_warmup initial solutions, we
!> set the time rt_t times = -n_warmup * rt_t%timestep, making sure to call
!> rt_t%rt_eos for the initial MHD state.
!===============================================================================
SUBROUTINE init (self)
  USE link_mod
  class (rt_t):: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
  class (rt_t),     pointer :: omega
  class (link_t),   pointer :: link
  class (*), pointer :: dummy
  integer :: n_omega, i_omega
  integer:: dir, i_mu, i_phi, iostat, n
  real, dimension(:), allocatable:: mu, w_mu
  ! ----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%init')
  p => self%patch
  !-----------------------------------------------------------------------------
  ! interconnect itself
  !-----------------------------------------------------------------------------
  nullify(self%link)
  allocate(link)
  dummy => p%connect%rt
  select type(dummy)
  class is (task_t)
    link%task => dummy
  end select
  self%link => link
  !-----------------------------------------------------------------------------
  ! Initialize the patch
  !-----------------------------------------------------------------------------
  call self%init_unique                                                         ! get an unique ID
  self%track   = .false.
  self%kind    = 'rt_t'                                                         ! to be used in data_io_mod
  call self%patch_t%setup(position=p%position, size=p%size, n=p%n, ng=p%ng, &
                          nt=self%nt, nv = self%n_lambda, nw = 1)
  !-----------------------------------------------------------------------------
  ! Synchronize with MHD patch
  ! ---------------------------------------------------------------------------- 
  self%box      = p%box
  self%size     = p%size
  self%origin   = p%origin
  self%llc_cart = p%llc_cart
  self%llc_nat  = p%llc_nat
  self%restart  = p%restart
  self%time     = p%time
  self%t(self%it) = self%time
  self%out_next = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time             ! add an `eps` to avoid round-off errors
  self%dt_fixed = self%rt_timestep                                              ! Fixed timestep for RT
  self%iout = self%restart+1
  self%grace    = self%rt_grace
  self%periodic = p%periodic
  call self%set (bits%no_density)
  !-----------------------------------------------------------------------------
  ! Check if this MHD patch needs RT
  !-----------------------------------------------------------------------------
  call self%setup_boundaries()
  if(self%on) then
    !---------------------------------------------------------------------------
    ! Allocate heating for MHD
    !---------------------------------------------------------------------------
    if (.not.allocated(self%patch%heating_per_unit_volume)) then
      allocate(self%patch%heating_per_unit_volume(self%gn(1), self%gn(2), self%gn(3)))
      self%patch%heating_per_unit_volume = 0.0
    end if
    !---------------------------------------------------------------------------
    ! Integration weights (w_bin is baked into the source function in this case)
    !---------------------------------------------------------------------------
    allocate (self%w_bin (self%n_lambda))
    self%w_bin = 1.0                                                            ! Frequency bin weights are baked into EOS
    n_omega      = merge (1+self%n_phi*(self%n_mu-1), 0, self%n_mu > 0)
    n_omega      = 2*n_omega
    self%n_omega = n_omega
    self%dt_fixed = self%rt_timestep
    !-----------------------------------------------------------------------------
    ! Assign angles and integration weight and allocate the required RT solvers
    ! for all angles
    !-----------------------------------------------------------------------------
    if (self%n_mu > 0) then
      allocate (mu(self%n_mu),w_mu(self%n_mu))
      call radau%init (self%n_mu, mu, w_mu)
    end if
    !---------------------------------------------------------------------------
    ! Allocate data, belonging to the main task
    !---------------------------------------------------------------------------
    allocate (self%rk (self%gn(1),self%gn(2),self%gn(3),self%n_lambda))
    allocate (self%tt (self%gn(1),self%gn(2),self%gn(3)))
    allocate (self%src(self%gn(1),self%gn(2),self%gn(3),self%n_lambda))
    call io%bits_mem(storage_size(self%tt  ),product(shape(self%tt )),'tt')
    call io%bits_mem(storage_size(self%rk  ),product(shape(self%rk )),'rk')
    call io%bits_mem(storage_size(self%src ),product(shape(self%src)),'src')
    !---------------------------------------------------------------------------
    ! Allocate n_omega copies of self, cloning self
    !---------------------------------------------------------------------------
    if (n_omega > 0) allocate (self%omega(n_omega), source=self)
    self%i_omega = 0
    !---------------------------------------------------------------------------
    ! Initialize the clones of self
    !---------------------------------------------------------------------------
    do i_omega=1,self%n_omega
      call self%init_omega(i_omega)
    end do
    !-----------------------------------------------------------------------------
    ! Assign angles -- Should change to 
    !-----------------------------------------------------------------------------
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
          omega%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + (dir+1)/2.0, 2.0) - 1.0)
          if (i_mu == self%n_mu) omega%rt_phi = 0.0
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
    ! Finalize the main RT patch. For now, it needs only one nv entry for q with
    ! nt time slots.
    !---------------------------------------------------------------------------
    self%nv = 1
    if (associated(self%mem)) deallocate (self%mem)
    nullify (self%mem)
    allocate (self%mem(self%mesh(1)%gn, &
                       self%mesh(2)%gn, &
                       self%mesh(3)%gn, &
                       self%nv,self%nt,self%nw))
    self%mem_allocated = .true.               ! mark
    call io%bits_mem (storage_size(self%mem), product(shape(self%mem)), 'mem')
  else
    call self%dealloc()
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize rt sub-tasks; note that all self values are initially cloned
!===============================================================================
SUBROUTINE init_omega (self, i_omega)
  class (rt_t):: self
  integer:: i_omega
  ! ............................................................................
  class (rt_t), pointer:: omega
  ! ----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%init')
  !
  omega => self%omega(i_omega)
  ! -- assign MPI-unique ID
  omega%rank = self%rank
  ! -- assign unique and globally valid ID
  call omega%init_unique
  omega%track   = .false.
  omega%kind    = 'rt_t'   ! to be used in data_io_mod
  omega%nt      = self%nt
  allocate (omega%t(self%nt), omega%dt(self%nt), omega%iit(self%nt))
  omega%t       = self%time
  omega%dt      = self%dt_fixed
  omega%iit     = 1
  omega%nv      = self%n_lambda
  allocate (omega%unsigned(self%n_lambda), omega%pervolume(self%n_lambda))
  omega%unsigned  = .false.
  omega%pervolume = .false.
  omega%it      = max(omega%nt-1,1)
  omega%new     = omega%nt
  omega%nw      = 1
  omega%i_omega = i_omega
  omega%time    = self%time
  omega%status  = self%status
  omega%verbose = verbose
  call omega%set (bits%no_density)
  omega%grace   = self%grace
  omega%tt      => self%tt
  omega%rk      => self%rk
  omega%src     => self%src
  omega%rt_grace    = self%rt_grace
  omega%rt_llc      = self%rt_llc
  omega%rt_extent   = self%rt_extent
  omega%rt_timestep = self%rt_timestep
  omega%vertical_propagation = self%vertical_propagation
  omega%constant_q  = self%constant_q
  nullify (omega%link)
  allocate (omega%link)
  omega%link%task => omega
  allocate (omega%mem   (self%gn(1),self%gn(2),self%gn(3),self%n_lambda,self%nt,1))
  call io%bits_mem(storage_size(omega%mem)   ,product(shape(omega%mem))   ,'q')
  omega%q => omega%mem(:,:,:,:,omega%it,1)
  omega%q = 0.0
  omega%patch => self%patch
  call omega%setup_boundaries()
  call omega%rt_boundaries%init(omega%mesh)
  !if (associated(omega%mesh)) call MeshRecycler (omega%mesh)
  !omega%mesh=> self%mesh
  !
  call trace%end()
END SUBROUTINE init_omega

!===============================================================================
!> Determine if a task is needed (is upstream).  We are being presented with
!> the nbors of the MHD task, so there is no task above the top-most task, but
!> there may well be tasks below the rt_bottom, which are not needed as RT nbors
!===============================================================================
FUNCTION needs (self, omega, task)
  logical:: needs
  class(rt_t) :: self
  class(rt_t):: omega
  class(task_t), pointer:: task
  !class(gpatch_t), pointer :: patch
  real:: ds1(3), ds2(3), sinth
  real, parameter:: eps=1e-4
  integer:: is1(3), is2(3), i, lbit(3), ubit(3)
  logical:: lower, upper
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
END FUNCTION needs

!===============================================================================
!> Check if patch belongs to an RT box or is at its boundary.  The mesh flags
!> lower_boundary and upper_boundary are set for the tasks that contain bndries,
!> and self%on is set false for tasks that are entirely outside.
!===============================================================================
SUBROUTINE setup_boundaries(self)
  class (rt_t) :: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3), rt_urc(3)
  integer :: lbit(3), ubit(3), i
  ! ............................................................................
  call self%init_bdries
  llc = self%position-0.5_8*self%size
  urc = self%position+0.5_8*self%size
  rt_urc = self%rt_llc+self%rt_extent
  !
  lbit = [bits%xl, bits%yl, bits%zl]                                            ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]                                            ! upper boundary bits -- convenience
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
    if ((llc(i) > rt_urc(i))) then
      call self%patch%set(bits%no_rt)
      self%on = .false.
    else if (urc(i) >= rt_urc(i)) then
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
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('rt_t%pre_update', itimer=itimer)
  if (self%constant_q) then
    q => self%mem(:,:,:,1,self%it,1)
    self%patch%heating_per_unit_volume = &
    self%patch%heating_per_unit_volume + q
  else
    call io%abort("time-interpolated q in rt_t%pre_update is not implemented yet")
  end if
  call trace%end(itimer)
END SUBROUTINE pre_update

!===============================================================================
!> This is called from MHD tasks, after updates, with self = main, and we use it
!> to copy the net heating and temperature arrays over to the MHD memory, for
!> output.
!===============================================================================
SUBROUTINE post_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('rt_t%post_update', itimer=itimer)
  if (self%constant_q) then
    q => self%mem(:,:,:,1,self%it,1)
    self%patch%mem(:,:,:,self%patch%idx%qr,self%patch%new,1) = q
    self%patch%mem(:,:,:,self%patch%idx%tt,self%patch%new,1) = self%tt
  else
  
  end if
  call trace%end(itimer)
END SUBROUTINE post_update

!===============================================================================
!> The rt_t%update procedure gets called directly by the task dispatcher, when
!> it has concluded that the task is ready to update.  An update of the main RT
!> task consists of summing up the results from the RT sub-tasks, while the
!> sub-task updates perform the actual RT solutions along rays.
!===============================================================================
SUBROUTINE update (self)
  class (rt_t):: self
  !.............................................................................
  class(rt_t), pointer :: omega
  integer:: i_omega, i_bin
  real, dimension(:,:,:), pointer:: q
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) then
  call io%abort ('some rts are off - check if that is expected and remove this abort')
  return
  end if
  call trace%begin ('rt_t%update', itimer=itimer)
  !
  if (self%id == io%id_debug) &
    write (io%output,*) 'DBG rt_t%update: id, i_omega, BV, time =', &
      self%id, self%i_omega, self%is_set(bits%boundary), self%is_set(bits%virtual), &
      self%time
  !-----------------------------------------------------------------------------
  ! Main RT task, summing up the contributions.  Since all RT sub-tasks then
  ! must have finished, we can use the opportunity to call rt_eos at this time.
  !-----------------------------------------------------------------------------
  if (self%i_omega == 0) then
    ! -- Sum Q over angles and bin
    q       => self%mem(:,:,:,1,self%new,1)
    q = 0.0
    ! -- multiply with the stored net intensities
    do i_omega=1,self%n_omega
      omega => self%omega(i_omega)
      do i_bin=1,self%n_lambda
        q = q + self%rk(:,:,:,i_bin)*self%w_bin(i_bin)*omega%w* &
                omega%mem(:,:,:,i_bin,omega%it,1)
      end do
    end do
    ! -- Update EOS
    call self%rt_eos(self%time)
  !-----------------------------------------------------------------------------
  ! RT sub-task, solving for the radiation field in a given direction
  !-----------------------------------------------------------------------------
  else
    ! -- Call RT solver, storing result in %it slot, at time = t(%it)
    self%q => self%mem(:,:,:,:,self%it,1)
    call self%solve ()
  end if
  !
  if ((verbose > 1) .or. (verbose > 0 .and. self%i_omega == 0)) &
    write(io%output,'(i6,2f7.3,2x,a,f8.4)') &
      self%id, self%position(3), self%mu, &
      'rt_t%update:  RT task updated to t =', self%time
  !-----------------------------------------------------------------------------
  ! By advancing the time a factor mu slower for inclined tasks they are forced
  ! to update more often, thus propagation the solution by about the same amount
  ! vertically.  NOTE (FIXME):  Purely horizontal rays need special treatment
  !-----------------------------------------------------------------------------
  if ((self%vertical_propagation).and.(self%i_omega>0)) then
    self%dtime = self%dt_fixed*abs(self%mu)
  else
    self%dtime = self%dt_fixed
  end if
  call self%rotate()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Compute the RT EOS quantities for the patch, using self%eos_time and an OMP
!> critical region to guard against duplicate calculations when multi-threading.
!===============================================================================
SUBROUTINE rt_eos (self, time)
  class (rt_t):: self
  real(8):: time
  !.............................................................................
  class (gpatch_t), pointer :: patch
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e
  integer:: m(4), ib, i3, jt(2)
  real:: z, pt(2)
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  patch => self%patch
  !if (self%eos_time < time) then
    !$omp critical (eos_cr)
  !  if (self%eos_time < time) then
      ! -- interpolate to the self%eos_time
      call patch%time_interval (time, jt, pt, all=.true.)
      m = shape(self%rk)
      allocate (d(m(1),m(2),m(3)),e(m(1),m(2),m(3)))
      d = patch%mem(:,:,:,patch%idx%d,jt(1),1)*pt(1) &
        + patch%mem(:,:,:,patch%idx%d,jt(2),1)*pt(2)
      e = patch%mem(:,:,:,patch%idx%e,jt(1),1)*pt(1) &
        + patch%mem(:,:,:,patch%idx%e,jt(2),1)*pt(2)
      if (verbose > 1) then
        call patch%stats (d, 'rt_eos')
        call patch%stats (e, 'rt_eos')
      end if
      call eos%lookup_table (m, d=d, e=e, rk=self%rk, src=self%src, tt=self%tt)
      if (verbose > 3) then
        write(io%output,*) patch%id,'rt_t%rt_eos:  d',minval(d),maxval(d)
        write(io%output,*) patch%id,'rt_t%rt_eos:  e',minval(e),maxval(e)
        write(io%output,*) patch%id,'rt_t%rt_eos: rk',minval(self%rk),maxval(self%rk)
        write(io%output,*) patch%id,'rt_t%rt_eos: tt',minval(self%tt),maxval(self%tt),scaling%temp
      end if
      deallocate (d, e)
      self%eos_time = time
    !else
    !  if (debug(patch,0)) then
    !    write(io%output,*) patch%id, 'DBG rt_t%update: EOS tim', self%eos_time
    !    write(io%output,*) patch%id, 'DBG rt_t%update: EOS REJ', time
    !  end if
    !end if
    !$omp end critical (eos_cr)
  !else
  !  write(io%output,*) patch%id, 'rt_t%rt_eos: WARNING, called w/o need', self%eos_time-time
  !end if
  call trace%end (itimer)
END SUBROUTINE rt_eos

!===============================================================================
!> Call integral short characteristics RT solver.
!===============================================================================
SUBROUTINE solve (self)
  class(rt_t):: self
  integer:: i3, ib
  integer, save:: itimer=0
  !.............................................................................
  real, pointer :: mem(:,:,:,:)
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) return
  call trace%begin ('rt_t%solve', itimer=itimer)
  !-----------------------------------------------------------------------------
  rt_integral%verbose = verbose
  !if (self%mesh(3)%lower_boundary) print *, 'should be lower', self%mesh(3)%lower_boundary, self%patch%mesh(3)%lower_boundary
  !if (self%boundaries%is_set(bits%zl+bits%zu)) print *, 'bits in bnd set properly' --VORKS
  call self%rt_boundaries%condition (mem=self%q, src=self%src, rk=self%rk, tt=self%tt, mu=self%mu, w=self%w)
  call rt_integral%solve (self%rk, self%src, self%q, self%mesh, self%mu, self%rt_phi)
  call trace%end (itimer)
END SUBROUTINE solve

!===============================================================================
!===============================================================================
SUBROUTINE warm_up(self)
  class(rt_t):: self
  !.............................................................................
  class(rt_t), pointer :: omega
  integer:: i_omega, i_bin, itimer=0
  real, dimension(:,:,:), pointer:: q
  real(8):: time
  !-----------------------------------------------------------------------------
  if ((.not.self%on)) return
  call trace%begin ('rt_t%warm_up', itimer=itimer)
  if (self%i_omega == 0) then
    ! -- Sum Q over angles and bin
    q       => self%mem(:,:,:,1,self%new,1)
    q = 0.0
    ! -- multiply with the stored net intensities
    do i_omega=1,self%n_omega
      omega => self%omega(i_omega)
      do i_bin=1,self%n_lambda
        q = q + self%rk(:,:,:,i_bin)*self%w_bin(i_bin)*omega%w* &
                omega%mem(:,:,:,i_bin,omega%it,1)
      end do
    end do
    ! -- Update EOS
    call self%rt_eos(self%time)
  else
    ! -- Call RT solver, storing result in %it slot, at time = t(%it)
    self%q => self%mem(:,:,:,:,self%it,1)
    call self%solve ()
  end if
  call trace%end (itimer)
END SUBROUTINE warm_up

!===============================================================================
!> NOTE: This procedure is currently NOT called.  It should be activated only
!> after checking carefully that it gives the same result as when nbor relations
!> are specified explicitly.  The procedure would only be relevant inside a loop
!> over the RT task list.
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(rt_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !.............................................................................
  class(rt_t), pointer:: rt_another
  class(task_t), pointer:: other
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%nbor_relations')
  !
  ! ----------------------------------------------------------------------------
  ! First check the main relations, based on location and level of refinement
  ! ----------------------------------------------------------------------------
  call self%task_t%nbor_relations(another, needed, needs_me, download)
  ! ----------------------------------------------------------------------------
  ! Now, if another is RT, add it appropriately, if it is MHD, it has to overlap
  ! and have different kind of relations
  ! ----------------------------------------------------------------------------
  other => another
  select type (another)
  class is (rt_t)
  rt_another => another
  ! ----------------------------------------------------------------------------
  ! for no chains solver nothing has to be done here - all logicals are true, 
  ! unless are made false by other reasons.
  ! ----------------------------------------------------------------------------
  if (self%needs(self,other)) then
    needs_me = .false.
  end if
  if (rt_another%needs(rt_another,self%link%task)) then
    needed   = .false.
    download = .false.
  end if
  if (self%vertical_propagation) then
    if (any(self%position(1:2) /= rt_another%position(1:2))) then
      needed = .false.
      needs_me = .false.
    end if
  else if (self%position(3) == rt_another%position(3)) then
    needed = .false.
    needs_me = .false.
  end if
  !
  if (self%i_omega /= rt_another%i_omega) download = .false.
  !
  class default
  other => another
  ! ----------------------------------------------------------------------------
  ! Check if this is the "master" MHD patch. It might not be associated yet, i.e.
  ! after moving to another MPI rank or changing level of refinement. in that 
  ! case, associate it.
  ! ----------------------------------------------------------------------------
  if (self%overlaps(other)) then
    download = .false.  ! don't download 
    if (.not.(associated(self%patch))) self%patch => gpatch%cast2gpatch(other)
  end if
  end select
  !
  call trace%end()
END SUBROUTINE nbor_relations

!===============================================================================
!> Return true if the patch is the horizontally closest to a give position, or
!> if the patch ID equals io%id_debug, unless level is less than verbose.
!===============================================================================
FUNCTION debug (self, level) RESULT (yes)
  class(gpatch_t):: self
  integer:: level
  real(8):: position(3)
  logical:: yes
  integer:: id
  !-----------------------------------------------------------------------------
  yes = maxval(abs(self%position(1:2)-dbg_position)/self%size(1:2)) < 0.55
  if (associated(self%parent)) then
    id = self%parent%id
  else
    id = self%id
  end if
  yes = yes .or. id == io%id_debug
  yes = yes .and. (verbose >= level)
END FUNCTION debug

!===============================================================================
SUBROUTINE task_info (self)
  class(rt_t):: self
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(1x,a,i6,2x,2l1,2x,3f9.3,i4,2f7.3)') &
    'rt_t%task_info: id, BV, pos, i_omega, mu, phi =', self%id, self%is_set(bits%boundary), &
    self%is_set(bits%virtual), self%position, self%i_omega, self%mu, self%rt_phi
END SUBROUTINE task_info

!===============================================================================
SUBROUTINE nbor_info (self, nbor)
  class(rt_t):: self
  class(link_t):: nbor
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(3x,a,i6,i4,2x,3l1,2x,2l1,2x,3f9.3,i4)') &
    'rt_t%nbor_info: id, rank, needed, needs_me, download, BV, pos =', &
    self%id, self%rank, nbor%needed, nbor%needs_me, nbor%download, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), &
    self%position, self%i_omega
END SUBROUTINE nbor_info

END MODULE rt_mod
