!> =============================================================================
!> Module to set up angles and bins for an RT patch with only absorption.
!> This means that the only extra memory with time slot that is needed is
!> one slot for the net radiation heating, \Sum_bins \Sum_angles Q_bin_angle
!>
!> The data structure consists of self, an extension of the MHD solver_t data
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
!> The self%time, self%rt_time, and self%rt_next control the scheduling of the
!> tasks.  The MHD tasks require that the times of their nbors are "ahead",
!> allowing for a grace fraction of the nbor time step.  The RT sub-tasks
!> advance their time by self%rt_time (effectively their self%dtime) when done,
!> hence they must update when their MHD master task passes their own time,
!> which conforms with standard nbor relations.
!>
!> Can this also be consistent with the "always update the oldest taske"
!> strategy of dispather0?  It will update the MHD task, as long as it does
!> not pass beyond the rt_next time (which generally should be the same while
!> the RT sub-task are idle).  Once it passes rt_next (which should also be
!> the self%time of all RT sub-tasks), they will be "older" than the MHD task,
!> and will all be updated (subject to their inter-dependencies).
!>
!> Initially, we wish that only the MHD states are ready to update, and we set
!> things up so the only thing they do is compute the EOS for the RT solvers.
!> So, we need self%time < omega%time, and self%eos_time < self%rt_time < 0.
!> =============================================================================
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE index_mod
  USE math_mod
  USE radau_mod
  USE solver_mod
  USE task_mod
  USE link_mod
  USE eos_mod
  USE download_mod
  USE scalar_mod
  USE bits_mod
  USE scaling_mod
  USE rt_integral_mod
  USE timestep_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from MHD, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(solver_t):: rt_t
    logical                    :: done=.false.
    logical                    :: first=.true.
    logical                    :: needs_nbors=.true.
    type(rt_t), pointer        :: omega(:)
    integer                    :: n_omega=0
    integer                    :: i_omega=0
    integer                    :: n_bin
    integer                    :: offset_o(2)
    real                       :: offset_p(2)
    real                       :: mu, phi, w
    real(8), pointer           :: rt_tm(:)
    real(8)                    :: mhd_time=0d0
    real(8)                    :: rt_time=0d0
    real(8)                    :: rt_next=0d0
    real(8)                    :: eos_time=-1d0
    real(8)                    :: rt_grace_time=0.0
    real, pointer              :: w_bin(:)
    real, pointer              :: q(:,:,:,:)
    real, pointer              :: tt(:,:,:)
    real, pointer              :: rk(:,:,:,:)
    real, pointer              :: src(:,:,:,:)
    real, pointer              :: layers(:,:,:,:,:)
  contains
    procedure:: init
    procedure:: update
    procedure:: init_nbors
    procedure:: rt_eos
    procedure:: boundary_conditions
    procedure:: solve
    procedure:: needs
    procedure, nopass:: task2rt
    procedure:: is_ahead_of
  end type
  real, save:: dtau_min=1e-3, dtau_max=30., dtau_dif=10.
  ! -- tolerate that RT is up to one rt_time behind MHD
  real, save:: grace=0.0
  real(8):: rt_bottom=0.45
  real(8):: dbg_position(2) = -1e30
  integer(8):: n_min=0, n_max=0, n_dif, n_ful=0
  integer, save:: verbose=0, order=1
  logical, save:: vertical_propagation=.false.
  type(rt_t), public:: rt
CONTAINS

!> =============================================================================
!> Read parameters
!> =============================================================================
SUBROUTINE init (self)
  class (rt_t):: self
  real(8):: rt_time=0.01_8
  integer, save:: n_mu=3, n_phi=4, n_bin=4, nt=3
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n
  real, dimension(:), allocatable:: mu, w_mu
  type(rt_t), pointer:: omega
  class(task_t), pointer:: task
  class(rt_t),   pointer:: rt_task
  logical, save:: first_time=.true.
  namelist /rt_params/ n_mu, n_phi, n_bin, verbose, rt_time, dtau_min, dtau_max, &
                       dtau_dif, vertical_propagation, dbg_position, grace, nt, &
                       order, rt_bottom
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_params, iostat=iostat)
    if (io%master) write (io%output, rt_params)
  end if
  !$omp end critical (input_cr)
  if (.not. io%master) verbose=0
  n_omega      = merge (1+n_phi*(n_mu-1), 0, n_mu > 0)
  n_omega      = 2*n_omega
  self%n_bin   = n_bin
  self%n_omega = n_omega
  self%rt_time = rt_time
  !-----------------------------------------------------------------------------
  ! Punt to the next level, initializing patch memory and meshes
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  self%nv  = 10
  self%nt  = 5
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Initialize variable indices, using slot 9 for temperature and slot 10 for Q
  !-----------------------------------------------------------------------------
  call self%idx%init (5, self%mhd)
  self%idx%tt = 10
  self%idx%qr = 9
  call self%idx%copy
  !-----------------------------------------------------------------------------
  ! Integration weights (w_bin is baked into the source function in this case)
  !-----------------------------------------------------------------------------
  allocate (self%w_bin (self%n_bin))
  self%w_bin = 1.0
  !-----------------------------------------------------------------------------
  ! Assign angles and integration weight and allocate the required RT solvers
  ! for all angles
  !-----------------------------------------------------------------------------
  if (n_mu > 0) then
    allocate (mu(n_mu),w_mu(n_mu))
    call radau%init (n_mu, mu, w_mu)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate data common to all sub-tasks
  !-----------------------------------------------------------------------------
  allocate (self%rt_tm(self%nt))
  allocate (self%rk (self%gn(1),self%gn(2),self%gn(3),n_bin))
  allocate (self%src(self%gn(1),self%gn(2),self%gn(3),n_bin))
  call io%bits_mem(storage_size(self%rk  ),product(shape(self%rk )),'rk')
  call io%bits_mem(storage_size(self%src ),product(shape(self%src)),'src')
  !-----------------------------------------------------------------------------
  ! Patches with top boundaries below rt_bottom have no RT and are not listed
  !-----------------------------------------------------------------------------
  if (self%position(3)-0.5_8*self%size(3) > rt_bottom) then
    call self%set (bits%no_rt)
    if (debug(self,1)) &
      write(io%output,*) self%id, 'does not need to do RT'
  !-----------------------------------------------------------------------------
  ! Patches with bottom boundaries below rt_bottom are considered RT bdry
  !-----------------------------------------------------------------------------
  else if (self%position(3)+0.5_8*self%size(3) > rt_bottom) then
    call self%boundaries%set (bits%zu)
    if (debug(self,1)) &
      write(io%output,*) self%id, 'is at the bottom RT boundary'
  end if
  !-----------------------------------------------------------------------------
  ! Allocate n_omega copies of self, cloning self via self
  !-----------------------------------------------------------------------------
  if (n_omega > 0) &
    allocate (self%omega(n_omega), source=self)
  !-----------------------------------------------------------------------------
  ! Initialize the directional solver connection links and data. 
  !-----------------------------------------------------------------------------
  self%rt_next  = self%time+self%rt_time
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    omega = self
    ! -- assign MPI-unique ID
    omega%rank = self%rank
    call omega%task_t%init_unique
    omega%track   = .false.
    omega%kind    = 'rt_t'   ! to be used in data_io_mod
    omega%nt      = nt
    allocate (omega%t(nt), omega%dt(nt), omega%iit(nt))
    omega%t       = self%time
    omega%dt      = self%rt_time
    omega%iit     = 1
    omega%nv      = n_bin
    allocate (omega%unsigned(n_bin), omega%pervolume(n_bin))
    omega%unsigned  = .false.
    omega%pervolume = .false.
    omega%it      = max(omega%nt-1,1)
    omega%new     = omega%nt
    omega%nw      = 1
    omega%i_omega = i_omega
    omega%grav    = self%grav
    omega%time    = self%rt_next
    omega%status  = self%status
    call omega%set (bits%no_density)
    omega%parent => self%link%task
    !! -- force 3 solutions before Q is valid in MHD update
    !omega%time    = -3.0_8*self%rt_time
    !! -- force starting at boundaries
    !if (self%mesh(3)%lower_boundary .or. self%mesh(3)%upper_boundary) &
      !omega%time    = -4.0_8*self%rt_time
    omega%grace   = grace
    omega%tt      => self%tt
    omega%rk      => self%rk
    omega%src     => self%src
    omega%pg      => self%pg
    nullify (omega%link)
    allocate (omega%link)
    omega%link%task => omega
    ! -- only allocate RT mem on patches above rt_bottom
    if (self%is_clear(bits%no_rt)) then
      allocate (omega%layers(self%gn(1),self%gn(2),           n_bin,3,3))
      allocate (omega%mem   (self%gn(1),self%gn(2),self%gn(3),n_bin,nt,1))
      call io%bits_mem(storage_size(omega%layers),product(shape(omega%layers)),'l')
      call io%bits_mem(storage_size(omega%mem)   ,product(shape(omega%mem))   ,'q')
      omega%q => omega%mem(:,:,:,:,self%it,1)
      omega%q = 0.0
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Assign angles
  !-----------------------------------------------------------------------------
  i_omega = 0
  do i_mu=1,n_mu
    n = n_phi
    if (i_mu == n_mu) n = 1
    do i_phi=1,n
      do dir=-1,1,2
        i_omega = i_omega+1
        omega => self%omega(i_omega)
        omega%mu  = dir*mu(i_mu)
        omega%w   = w_mu(i_mu)*math%pi4/n
        omega%phi = math%pi*(modulo(real(2*(i_phi-1))/n_phi + (dir+1)/2.0, 2.0) - 1.0)
        if (i_mu == n_mu) omega%phi = 0.0
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Check task angles
  !-----------------------------------------------------------------------------
  if (verbose > 0 .and. self%id == 1) then
    do i_omega=1,n_omega
      omega => self%omega(i_omega)
      if (verbose > 0) &
        write(io%output,'(a,3i7,f8.3)') &
          'rt_t%init: i_omega, elevation, phi, mu =', i_omega, &
          nint(asin(omega%mu)*180./math%pi), nint(omega%phi*180./math%pi), &
          omega%mu
    end do
  end if
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Run through all angular solvers and attempt to update them.  Solvers that
!> lack upstream info (based on internal criteria) will do nothing, while
!> solvers that have enough information will be updated.  The angular solvers
!> are so lightweight (taking of the order 100-200 ns) that it would be too
!> much overhead to consider them separate tasks.
!>
!> It would be easy to spawn independent OMP tasks for the ones that are
!> ready, and then wait on them, but for now we avoid complications.
!>
!> The update step computes the current RT patch time, based on the current
!> time in the attached RT solvers for each angle
!> =============================================================================
SUBROUTINE update (self)
  class (rt_t):: self
  integer:: i_omega, i_bin, n_solve, count, id, i3, l(3), u(3), jt(2), loc(3), iv
  real:: n_tot, pt(2)
  real, dimension(:,:,:), pointer:: q, d, e, rk
  class(rt_t), pointer:: omega, mhd
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%update', itimer=itimer)
  if (self%id == io%id_debug) &
    write (io%output,*) 'DBG rt_t%udpate: time =', self%time, &
    self%is_set(bits%boundary), self%is_set(bits%virtual)
  !-----------------------------------------------------------------------------
  ! The MHD tasks are called repeatedly, until their time has surpassed (one of)
  ! the RT tasks that it depends upon.  These should all have formal self%time
  ! set to the previous EOS time + rt_time.  Therefore, as soon as the MHD task
  ! time is actually set to the new time, all RT tasks should be behind it, and
  ! at least the boundary RT tasks should be able to execute.  When they do,
  ! they should be using EOS values from the time the MHD task has reached.  
  ! These must therefore be computed before setting mhd%time, to ensure that the
  ! EOS values are actually available when the RT tasks are released.
  !
  ! The next phase should then consist of all RT sub-tasks updating, in the order
  ! needed to ensure upstream tasks preceded downstream tasks.  A dead-lock 
  ! situation to avoid is to have an upstream RT task that cannot run, because
  ! it is actually too far ahead of other MHD tasks.   But this should not be
  ! able to happen, since by construction, after landing on the right side of
  ! the rabbit time, it did compute the new EOS values, so the RT tasks then
  ! only depend on being able to extrapolate nbor tasks Q values ahead in time,
  ! if necessary.  This will indeed be necessary, if the situations is the
  ! following:
  !
  ! ------------- BC --------------------------------------
  ! ------------- MHD1 ---------------|--+
  ! --------------RT1  ---------------|--+
  ! ------------- MHD2 ---------------|---+
  ! --------------RT2  ---------------|---+
  !
  ! Now RT1 can run, because it relies on the upstream BC.  But RT2 cannot run,
  ! unless a grace period is allow, because RT1 is at an earlier time.
  !
  ! Hence it is only possible to enforce strict time order btw the RT tasks if
  ! they all are solved at the same time.  So, either we maintain a strict
  ! cadence of RT times, or we must allow a grace time, at least as large as
  ! the largest difference btw MHD times.
  !
  ! A strict RT cadence has two (related) advantages:  1) no time interpolation
  ! will be needed when communicating RT data btw patches.  2) for the same
  ! reason, it won't be necessary to keep a number of time slices of the angular
  ! and bin dependent radiation field.
  !
  ! A slight disadvantage is the need to interpolate backwards in MHD time to
  ! the common EOS/RT time, plus the fact that the RT field is not quite as far
  ! advanced in time as the MHD variables are.
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) then
    if (verbose > 1) &
      write(io%output,*) self%id, 'rt_t%update: time, dtime, it, i_omega', &
        self%time, self%dtime, self%it, self%i_omega
    !---------------------------------------------------------------------------
    ! If the main task is called with self%time > self%rt_next, then the
    ! scheduler has concluded that the pre-requisites are in order, and hence
    ! we can update self%q, which is later used (and reused) as the next
    ! radiative heating.
    !---------------------------------------------------------------------------
    if (self%time >= self%rt_next) then
      ! -- Sum Q over angles and bin
      if (debug(self,0)) &
        write(io%output,*) self%id, 'DBG rt_t%update: compute Q at', self%rt_next
      self%tt => self%mem(:,:,:,self%idx%tt,self%it,1)
      q => self%mem(:,:,:,self%idx%qr,self%new,1)
      q = 0.0
      self%mem(:,:,:,self%idx%qr,self%it,2) = 0.0
      if (self%is_clear (bits%no_rt)) then
        ! -- multiply with the stored net intensities
        do i_omega=1,self%n_omega
          omega => self%omega(i_omega)
          do i_bin=1,self%n_bin
            q = q + self%rk(:,:,:,i_bin)*self%w_bin(i_bin)*omega%w* &
                    omega%mem(:,:,:,i_bin,omega%it,1)
          end do
          if (debug(self,1)) &
            write(io%output,'(i7,3x,2i3,2f7.3,3x,f8.4,3x,a,1p,2e12.4)') &
              self%id, jt, pt, self%time, 'q', self%fminval(q), self%fmaxval(q)
        end do
        if (debug(self,1)) then
          l = self%li
          u = self%ui
          loc = maxloc(abs(q(l(1):u(1),l(2):u(2),l(3):u(3)))) + l-1
          write(io%output,*) 'min/max loc, value =', loc, q(loc(1),loc(2),loc(3))
          do i_bin=1,self%n_bin
            write(io%output,*) 'i_bin, q1, q2 =', i_bin, &
              omega%mem(loc(1),loc(2),loc(3),i_bin,omega%it,1), &
              self%rk(loc(1),loc(2),loc(3),i_bin)
          end do
          do iv=1,self%nv
            write(io%output,*) 'iv, v', iv, self%mem(loc(1),loc(2),loc(3),iv,self%it,1)
          end do
          ! -- outgoing vertical ray
          omega => self%omega(self%n_omega)
          do i3=l(3),u(3)
            write(io%output,'(a,i4,1p,10e10.2)') &
              'rt_t%update: i3, q, q_cont, qmin, qmax =', i3, &
              q(loc(1),loc(2),i3), &
              omega%mem(loc(1),loc(2),i3,1,omega%it,1), &
              minval(q(l(1):u(1),l(2):u(2),i3)), &
              maxval(q(l(1):u(1),l(2):u(2),i3))
          end do
        end if
      end if
    else
      self%mem(:,:,:,self%idx%qr,self%new,1) = self%mem(:,:,:,self%idx%qr,self%it,1)
    end if
    !---------------------------------------------------------------------------
    ! Update the MHD, limiting the suggested time step to at most rt_time,
    ! to avoid synchronization problems
    !---------------------------------------------------------------------------
    call mhd_update (self)
    !---------------------------------------------------------------------------
    ! If self%time, after update to happen in self%rotate, will be larger than
    ! self%rt_next, then call the EOS, and increment self%rt_next to the next 
    ! anticipated RT time.
    !---------------------------------------------------------------------------
    if (self%is_clear(bits%no_rt)) then
      time = self%time + self%dtime
      if (debug(self,0)) then
        write(io%output,*) self%id, 'rt_t%update: time   ', time
        write(io%output,*) self%id, 'rt_t%update: rt_next', self%rt_next
      end if
      if (time >= self%rt_next) then
        ! -- make sure the new time slot is up-to-date
        !$omp atomic write
        self%t(self%new) = time
        ! -- make a new EOS calculation at self%time
        call self%rt_eos(self%rt_next)
        ! -- increase the "rabbit" time one RT time
        self%rt_next = self%rt_next + self%rt_time
      end if
    end if
    !---------------------------------------------------------------------------
    ! To avoid RT sub-tasks starting before the EOS data is updated, delay the
    ! call to rotate until after the rt_eos call
    !---------------------------------------------------------------------------
    call self%rotate
!    write(io_unit%mpi,'(f12.6,i4,i7,2x,a,f12.6)') &
!      wallclock(), omp%thread, self%id, 'rt_t%update: MHD time =', self%time
    if (debug(self,0)) &
      write(io%output,*) self%id, 'DBG rt_t%update: MHD at ', self%time
  !-----------------------------------------------------------------------------
  ! Call the RT solver for RT sub-tasks when self is not the main task, for as
  ! long as self%time is less than the "rabbit time" rt_next:  For the RT solver,
  ! self%time is the time of the RT solution, which generally lags behind MHD
  !-----------------------------------------------------------------------------
  else
    !---------------------------------------------------------------------------
    ! Get rt_t qualified mhd task pointer
    !---------------------------------------------------------------------------
    mhd => task2rt (self%parent)
    !---------------------------------------------------------------------------
    ! Local equivalent of rotate, where we do not use %new, but make sure that
    ! the Q result is stored at the time
    !---------------------------------------------------------------------------
    if (debug(self,1)) write(io%output,*) mhd%id, &
      'rt_t%update: update  RT, id, t = ', self%id, self%time, mhd%time
    self%it = mod(self%it,self%nt)+1
    self%iit(1:self%nt-1) = self%iit(2:self%nt)
    self%iit(self%nt) = self%it
    if (mhd%id==io%id_debug) &
      write(io%output,*) self%id, 'DBG rt_t%update: RT time =', self%time
    self%t  (self%it) = self%time
    self%dt (self%it) = self%rt_time
    if (debug(self,0)) write(io%output,*) mhd%id, &
      'rt_t%update: RT  at ', self%time, self%id
    if (mhd%eos_time /= self%time) then
      write(io%output,'(i7,2x,a,i7,2f12.6)') &
        self%id, 'WARNING: EOS out of sync: mhd%id, eos%time, rt%time =', &
        mhd%id, mhd%eos_time, self%time
!      write(io_unit%mpi,'(f12.6,i4,i7,2x,a,f12.6)') &
!        wallclock(), omp%thread, self%id, 'rt_t%update: WARNING  =', self%time
    end if
    !---------------------------------------------------------------------------
    ! -- Call RT solver, storing result in %it slot, at time = t(%it)
    self%q => self%mem(:,:,:,:,self%it,1)
    call self%solve ()
    !---------------------------------------------------------------------------
    self%rotated = .true.
    if (vertical_propagation) then
      self%dtime = self%rt_time*abs(self%mu)
      self%time = self%time + self%dtime
    else
      ! -- choose an RT time step that is at least the parent MHD time step
      self%dtime = self%rt_time
      self%time = self%time + self%dtime
    end if
    !---------------------------------------------------------------------------
!    write(io_unit%mpi,'(f12.6,i4,i7,2x,a,f12.6)') &
!      wallclock(), omp%thread, self%id, 'rt_t%update: RT  time =', self%time
    if (verbose > 2) &
      write(io%output,'(i6,2f7.3,2x,a,f8.4)') &
        self%id, self%position(3), self%mu, &
        'rt_t%update:  RT patch updated to t =', self%time
  end if
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> This replaces mhd_t%update, and is identical, except for the extra limit on
!> the timestep
!> =============================================================================
SUBROUTINE mhd_update (self)
  class(rt_t):: self
  integer, save:: itimer=0
  integer:: it, iv, nt
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !call self%gravity
  call self%output
  call self%pde
  call self%courant_condition
  if (self%is_clear(bits%no_rt)) then
    self%dtime = min(self%dtime, 0.999999d0*self%rt_time)
  end if
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  call self%gpatch_t%update
  call trace%end (itimer)
END SUBROUTINE mhd_update

!> =============================================================================
!> Add upstream and downstream tasks to an nbor list.  The upstream nbors should
!> have the settings needed=T, needs_me=F, download=T, while the downstream 
!> tasks should have needed=F, needs_me=T, download=F.
!>
!> Each omega subtask of an MHD task has, for each MHD nbor, either
!> 1) no connection, or
!> 2) an upstream connection, or
!> 3) a downstream connection
!>
!> The nbor's sub-task in the omega direction should not be added to then MHD
!> task's nbor list, and should, in the local sub-tasks nbor list, either
!> 1) not be added
!> 2) be added with needed=T, needs_me=F, download=T, or
!> 3) be added with needed=F, needs_me=t, download=F
!>
!> =============================================================================
SUBROUTINE init_nbors (self)
  class (rt_t):: self
  !.............................................................................
  integer:: i_omega
  type(rt_t), pointer:: omega
  class(link_t), pointer:: nbor, new
  class(task_t), pointer:: nbtask, raytask
  integer:: n_nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbors')
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    n_nbor = 0
    nbor => self%link%nbor
    do while (associated(nbor))
      ! if (self%id == io%id_debug) &
        ! write(io_unit%mpi,'(1x,a,i6,3f7.3,2x,3f7.3,l3)') &
        ! 'DBG rt_t%init_nbors: nbor, relpos, mu, is needed', nbor%task%id, &
        ! (nbor%task%position-self%position)/self%size, &
        ! sin(acos(omega%mu))*[cos(omega%phi),sin(omega%phi)],omega%mu, &
        ! self%needs(omega,nbor%task)
      if (self%needs(omega,nbor%task)) then
        !-----------------------------------------------------------------------
        ! Add a new nbor, a new nbor%link, and pointers from link to task and
        ! from task to link.  The omega-task on nbor already has its own link.
        !-----------------------------------------------------------------------
        nbtask => nbor%task
        select type (nbtask)
        class is (rt_t)
          raytask => nbtask%omega(i_omega)
          if (raytask%is_clear(bits%no_rt)) then
            select type (raytask)
            class is (rt_t)
            !-------------------------------------------------------------------
            ! Create an "upstrean type" nbor node pointing to the nbor sub-task
            !-------------------------------------------------------------------
            allocate (new)
            new%link => raytask%link
            new%task  => raytask
            ! -- set "upstream attributes"
            new%needed = .true.
            new%needs_me = .false.
            new%download = .true.
            ! -- add to the sub-tasks nbor list, in MPI rank order
            call omega%link%add_nbor_by_rank (new)
            !-------------------------------------------------------------------
            ! Create a "downstrean type" nbor node pointing back to omega
            !-------------------------------------------------------------------
            allocate (new)
            new%link => omega%link
            new%task  => omega%link%task
            new%needed = .false.
            new%needs_me = .true.
            new%download = .false.
            !! -- add it to the nbor's nbor chain, in rank order
            call raytask%link%add_nbor_by_rank (new)
            end select
          end if
        end select
        !-----------------------------------------------------------------------
        ! Add the new nbor to the sub-tasks nbor list, sorted by MPI rank
        !-----------------------------------------------------------------------
        n_nbor = n_nbor+1
      end if
      nbor => nbor%next
    end do
    if (verbose > 2) &
      write(io%output,*) omega%id, 'rt_t%init_nbors: n_nbor=', n_nbor
  end do
  self%needs_nbors = .false.
  call trace%end
END SUBROUTINE init_nbors

!> =============================================================================
!> Compute the RT EOS quantities for the patch, using self%eos_time to guard
!> againts duplicate calculations when multi-threading.  Check both before and
!> after entering a critical region, so most threads never hit the critical region.
!> =============================================================================
SUBROUTINE rt_eos (self, time)
  class (rt_t):: self
  real(8):: time
  !.............................................................................
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e, tt
  integer:: m(4), ib, i3, jt(2)
  real:: z, pt(2)
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  if (self%eos_time < time) then
    !$omp critical (eos_cr)
    if (self%eos_time < time) then
      if (debug(self,0)) then
        write(io%output,*) self%id, 'DBG rt_t%update: EOS at ', time
      end if
      ! -- interpolate to the self%eos_time
      call self%time_interval (time, jt, pt, all=.true.)
      m = shape(self%rk)
      allocate (d(m(1),m(2),m(3)),e(m(1),m(2),m(3)))
      d = self%mem(:,:,:,self%idx%d,jt(1),1)*pt(1) &
        + self%mem(:,:,:,self%idx%d,jt(2),1)*pt(2)
      e = self%mem(:,:,:,self%idx%e,jt(1),1)*pt(1) &
        + self%mem(:,:,:,self%idx%e,jt(2),1)*pt(2)
      if (verbose > 1) then
        call self%stats (d, 'rt_eos')
        call self%stats (e, 'rt_eos')
      end if
      tt => self%mem(:,:,:,self%idx%tt,self%new,1)
      call eos%lookup_table (m, d=d, e=e, rk=self%rk, src=self%src, tt=tt)
      if (verbose > 3) then
        write(io%output,*) self%id,'rt_t%rt_eos:  d',minval(d),maxval(d)
        write(io%output,*) self%id,'rt_t%rt_eos:  e',minval(e),maxval(e)
        write(io%output,*) self%id,'rt_t%rt_eos: rk',minval(self%rk),maxval(self%rk)
        write(io%output,*) self%id,'rt_t%rt_eos: tt',minval(tt),maxval(tt),scaling%temp
      end if
      deallocate (d, e)
      self%eos_time = time
!      write(io_unit%mpi,'(f12.6,i4,i7,2x,a,f12.6)') &
!        wallclock(), omp%thread, self%id, 'rt_t%rt_eos: eos_time =', self%eos_time
    else
      if (debug(self,0)) then
        write(io%output,*) self%id, 'DBG rt_t%update: EOS tim', self%eos_time
        write(io%output,*) self%id, 'DBG rt_t%update: EOS REJ', time
      end if
    end if
    !$omp end critical (eos_cr)
  else
    write(io%output,*) self%id, 'rt_t%rt_eos: WARNING, called w/o need', self%eos_time-time
  end if
  call trace%end (itimer)
END SUBROUTINE rt_eos

!> =============================================================================
!> Call integral short characteristics RT solver.
!> =============================================================================
SUBROUTINE solve (self)
  class(rt_t):: self
  integer:: i3, ib
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%solve', itimer=itimer)
  !-----------------------------------------------------------------------------
  rt_integral%verbose = verbose
  call self%boundary_conditions
  call rt_integral%solve (self%rk, self%src, self%q, self%mesh, self%mu, self%phi)
  call trace%end (itimer)
END SUBROUTINE solve

!> =============================================================================
!> Vertical boundary conditions; here we explicitly assume the vertical direction
!> to be downwards, so the lower boundary index-wise is the top physical boundary,
!> where I=0 => Q=-S.
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class (rt_t):: self
  real:: h_scale, tau
  integer:: i1, i2, i3, ib
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%boundary_conditions')
  if (self%mesh(3)%lower_boundary .and. self%mu > 0) then
    i3 = self%mesh(3)%lo
    do ib=1,self%n_bin
      do i2=self%mesh(2)%lo,self%mesh(2)%uo
      do i1=self%mesh(1)%lo,self%mesh(1)%uo
        ! -- Estimate the boundary Q from an estimate of the total optical depth
        !h_scale = self%tt(i1,i2,i3)/(self%grav*scaling%temp)
        h_scale = 0.15
        tau = h_scale*self%rk(i1,i2,i3+1,ib)/self%mu
        self%mem(i1,i2,i3,ib,self%it,1) = -self%src(i1,i2,i3,ib)*exp(-tau)
      end do
      end do
    end do
    if (verbose > 2) &
      write(io%output,*) self%id, 'rt_t%lower_BC', scaling%temp
  else if (self%mesh(3)%upper_boundary .and. self%mu < 0) then
    self%mem(:,:,self%mesh(3)%uo,:,self%it,1) = 0.0
    if (verbose > 2) &
      write(io%output,*) self%id, 'rt_t%upper_BC'
  end if
  call trace%end
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Determine if a task is needed (is upstream).  We are being presented with
!> the nbors of the MHD task, so there is no task above the top-most task, but
!> there may well be tasks below the rt_bottom, which are not needed as RT nbors
!> =============================================================================
FUNCTION needs (self, omega, task)
  logical:: needs
  class(rt_t):: self, omega
  class(task_t), pointer:: task
  class(rt_t), pointer:: rt_task
  real:: ds1(3), ds2(3), sinth
  real, parameter:: eps=1e-4
  integer:: i, is1(3), is2(3)
  logical:: lower, upper
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  select type (task)
  class is (rt_t)
    rt_task => task
  end select
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and omega points inward, or the task has no RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  ds1 = self%position-task%position
  where (self%periodic)
    ds1 = modulo(ds1 + 0.5*self%box, self%box) - 0.5*self%box
  end where
  ds1 = ds1/self%box
  sinth = merge(0.0, sqrt(1.0-omega%mu**2), abs(omega%mu)==1.0)
  ds2 = [sinth*cos(omega%phi),sinth*sin(omega%phi),omega%mu]
  is1 = 0
  is2 = 0
  do i=1,3
    if (ds1(i) > +eps) is1(i)=1
    if (ds2(i) > +eps) is2(i)=1
    if (ds1(i) < -eps) is1(i)=-1
    if (ds2(i) < -eps) is2(i)=-1
  end do
  lower = self%boundaries%is_set(bits%zl)
  upper = self%boundaries%is_set(bits%zu)
  needs = .true.
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%no_rt)) then
    needs = .false.
  !-----------------------------------------------------------------------------
  ! If the nbor lies outside the boundaries it is not needed
  !-----------------------------------------------------------------------------
  else if (upper .and. ds1(3)==1) then
    needs = .false.
  else if (lower .and. ds1(3)==-1) then
    needs = .false.
  !-----------------------------------------------------------------------------
  ! Vertical rays depend only on exact alignment
  !-----------------------------------------------------------------------------
  else if (abs(omega%mu)==1.0) then
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
  i = omega%i_omega
  if (needs) &
    write(io_unit%mpi, &
    '("DBG rt_t",i7,4x,"needs task",i7,4x,"rel pos:",3f8.3,4x,"ray:",3f8.3,6i3,l3)') &
    omega%id, rt_task%omega(i)%id, ds1, ds2, is1, is2, needs
END FUNCTION needs

!> =============================================================================
!> self is the source task
!> =============================================================================
FUNCTION is_ahead_of (self, target)
  logical:: is_ahead_of
  class(rt_t):: self
  class(task_t):: target
  !-----------------------------------------------------------------------------
  ! If target is an MHD task, require that RT tasks are at least at the same
  ! time, while accepting a grace period for other MHD tasks
  !-----------------------------------------------------------------------------
  select type (target)
  class is (rt_t)
  call trace%begin ('rt_t%is_ahead_of',2)
  if (target%i_omega==0) then
    if (self%i_omega>0) then
      is_ahead_of = self%time >=  target%time
    else
      is_ahead_of = self%time + self%grace*self%dtime >= target%time
    end if
  !-----------------------------------------------------------------------------
  ! If target is an RT task, require that other RT tasks are strictly ahead,
  ! while accepting also equal MHD task times
  !-----------------------------------------------------------------------------
  else
    if (self%i_omega>0) then
      is_ahead_of = self%time >  target%time
    else
      is_ahead_of = self%time >= target%time
    end if
  end if
  end select
  call trace%end
END FUNCTION is_ahead_of

!> =============================================================================
!> Return true if the patch is the horizontally closest to a give position, or
!> if the patch ID equals io%id_debug, unlesslevel is less than verbose.
!> =============================================================================
FUNCTION debug (self, level) RESULT (yes)
  class(rt_t):: self
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

!> =============================================================================
!> Upgrade task pointer to experiment level
!> =============================================================================
FUNCTION task2rt (task) RESULT (rt)
  class(task_t), pointer:: task
  class(rt_t), pointer:: rt
  !...........................................................................
  select type (task)
  class is (rt_t)
  rt => task
  end select
END FUNCTION task2rt

END MODULE rt_mod
