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
  end type
  real, save:: dtau_min=1e-3, dtau_max=30., dtau_dif=10.
  ! -- tolerate that RT is up to one rt_time behind MHD
  real, save:: grace=0.0
  real(8):: rt_bottom=0.45
  real(8):: dbg_position(2) = -1e30
  integer(8):: n_min=0, n_max=0, n_dif, n_ful=0
  integer, save:: verbose=0, order=1
  logical, save:: vertical_propagation=.false.
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
    if (verbose > 0) &
      print *, self%id, 'does not need to do RT'
  !-----------------------------------------------------------------------------
  ! Patches with bottom boundaries below rt_bottom are considered RT bdry
  !-----------------------------------------------------------------------------
  else if (self%position(3)+0.5_8*self%size(3) > rt_bottom) then
    call self%boundaries%set (bits%zu)
    if (verbose > 0) &
      print *, self%id, 'is at the bottom RT boundary'
  end if
  !-----------------------------------------------------------------------------
  ! Allocate n_omega copies of self, cloning self via self
  !-----------------------------------------------------------------------------
  if (n_omega > 0) &
    allocate (self%omega(n_omega), source=self)
  !-----------------------------------------------------------------------------
  ! Initialize the directional solver connection links and data.  For now, we
  ! give the task IDs that make them easy to recognize in debugging.  Later,
  ! in rt_setup, we optionally give them unique task IDs across MPI ranks.
  !-----------------------------------------------------------------------------
  self%eos_time = self%time
  self%rt_next  = self%time+self%rt_time
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    omega = self
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
    omega%id      = self%id*1000 + i_omega
    omega%i_omega = i_omega
    omega%status  = self%status
    omega%grav    = self%grav
    call omega%set (bits%no_density)
    omega%time    = self%rt_next
    !! -- force 3 solutions before Q is valid in MHD update
    !omega%time    = -3.0_8*self%rt_time
    !! -- force starting at boundaries
    !if (self%mesh(3)%lower_boundary .or. self%mesh(3)%upper_boundary) &
      !omega%time    = -4.0_8*self%rt_time
    omega%grace   = grace
    omega%rt_next = self%rt_next
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
        print '(a,3i7,f8.3)', 'rt_t%init: i_omega, elevation, phi, mu =', &
        i_omega, nint(asin(omega%mu)*180./math%pi), nint(omega%phi*180./math%pi), &
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
  integer:: i_omega, i_bin, n_solve, count, id, i3, l(3), u(3), jt(2)
  real:: n_tot, pt(2)
  real, dimension(:,:,:), pointer:: q
  class(rt_t), pointer:: omega, mhd
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Call (only) the MHD solver, until it has passed the rt_next time, then
  ! work repeatedly with the RT solver, until that is ahead again
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Call the MHD solver update when self is the main task
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) then
    !---------------------------------------------------------------------------
    ! Allow the MHD solution to advance at most rt_time ahead of the RT, then
    ! update the rt_next "rabbit time", which is the time RT should advance to
    !---------------------------------------------------------------------------
    if (verbose > 1) &
      print *, self%id, 'rt_t%update: time, dtime, it, i_omega', &
        self%time, self%dtime, self%it, self%i_omega
    ! -- Sum Q over angles and bin
    if (verbose > 1) &
      print *, self%id, 'rt_t%update: compute Q at', self%rt_next
    self%tt => self%mem(:,:,:,self%idx%tt,self%it,1)
    q => self%mem(:,:,:,self%idx%qr,self%it,1)
    q = 0.0
    self%mem(:,:,:,self%idx%qr,self%it,2) = 0.0
    if (self%is_clear (bits%no_rt)) then
      do i_omega=1,self%n_omega
        omega => self%omega(i_omega)
        call omega%time_interval (self%time, jt, pt, all=.true.)
        if (order<1) then
          pt(2) = min(pt(2),1.0)
          pt(1) = 1.0-pt(2)
        end if
        ! -- these should ALL be extrapolations forward in time
        do i_bin=1,self%n_bin
          q = q + self%rk(:,:,:,i_bin)*self%w_bin(i_bin)*omega%w* &
                  (omega%mem(:,:,:,i_bin,jt(1),1)*pt(1) + &
                   omega%mem(:,:,:,i_bin,jt(2),1)*pt(2))
        end do
        if (verbose > 0) &
          print '(i7,3x,2i3,2f7.3,3x,f8.4,3x,a,1p,2e12.4)', &
            self%id, jt, pt, self%time, 'q', self%fminval(q), self%fmaxval(q)
      end do
      if (self%id == io%id_debug) then
        l = self%li
        u = self%ui
        do i3=l(3),u(3)
        print '(a,i4,1p,10e10.2)', 'rt_t%update: i3, q =', i3, &
          minval(q(l(1):l(1),l(2):u(2),i3)), &
          maxval(q(l(1):l(1),l(2):u(2),i3))
        end do
      end if
      if (verbose > 1) &
        call rt_integral%diagnostics ('at MHD update')
    end if
    ! -- Update the MHD state
    if (verbose > 1) &
      print *, self%id, 'rt_t%update: update MHD at', self%time
    !---------------------------------------------------------------------------
    ! Update the MHD, and also rotate the time slots (which updates self%time)
    !---------------------------------------------------------------------------
    call self%solver_t%update
    call self%rotate
    if (verbose > 0) &
      print *, self%id, 'rt_t%update: advanced time to', self%time
    if (self%id == io%id_debug) then
      l = self%li
      u = self%ui
      do i3=l(3),u(3)
      print '(a,i4,1p,10e10.2)', 'rt_t%update: i3, q =', i3, &
        q(l(1):l(1)+9,(l(2)+u(2))/2,i3)
      end do
    end if
    if (self%is_clear(bits%no_rt)) then
      if (self%time > self%rt_next) then
        ! -- maake a new EOS calculation, interpolated to self%rt_next
        call self%rt_eos
        ! -- increase the "rabbit" time one RT time
        self%rt_next = self%rt_next + self%rt_time
        do i_omega=1,self%n_omega
          self%omega(i_omega)%rt_next = &
          self%omega(i_omega)%rt_next + self%rt_time
        end do
        if (verbose > 0) &
          print *, self%id, 'rt_t%update: (1) increased rt_next to', self%rt_next
      end if
    end if
  !-----------------------------------------------------------------------------
  ! Call the RT solver for RT sub-tasks when self is not the main task, for as
  ! long as self%time is less than the "rabbit time" rt_next:  For the RT solver,
  ! self%time is the time of the RT solution, which generally lags behind MHD
  !-----------------------------------------------------------------------------
  else
    !---------------------------------------------------------------------------
    ! If the RT task has updated before, but is still behind MHD, then
    ! increment its "rabbit time" rt_next
    !---------------------------------------------------------------------------
    if (self%time==self%rt_next .and. self%time<self%parent%time) then
      self%rt_next = self%rt_next + self%rt_time
      if (verbose > 0) &
        print *, self%id, 'rt_t%update: (3) increased rt_next to', &
          self%rt_next, self%time, self%parent%time
    end if
    ! -- Prepare a new time slot for Q
    self%it = mod(self%it,self%nt)+1
    self%iit(1:self%nt-1) = self%iit(2:self%nt)
    self%iit(self%nt) = self%it
    self%t  (self%it) = self%time
    self%dt (self%it) = self%dtime
    ! -- Call RT solver, storing result in %it slot, at time = t(%it)
    self%q => self%mem(:,:,:,:,self%it,1)
    call self%solve ()
!    print *, self%id, 'min, max', &
!      self%fminval(self%q(:,:,:,1)), &
!      self%fmaxval(self%q(:,:,:,1))
    if (verbose > 2) &
      call rt_integral%diagnostics ('after omega solve')
    !---------------------------------------------------------------------------
    ! Local equivalent of rotate, where we do not use %new, but make sure that
    ! the Q result is stored at the time 
    !---------------------------------------------------------------------------
    if (vertical_propagation) then
      self%dtime = self%rt_time*abs(self%mu)
      self%time = self%time + self%dtime
      if (abs(self%time-self%rt_next) < 0.001*self%rt_time) then
        self%time = self%rt_next
      end if
    else
      self%dtime = self%rt_time
      self%time = self%rt_next
    end if
    if (verbose > 0) &
      print '(i6,2x,a,i3,4f8.4)', self%id, &
        'rt_t%update: it, t(it), time, dtime, rt_next, t =', &
        self%it, self%t(self%it), self%time, self%dtime, self%rt_next
    self%rotated = .true.
    !---------------------------------------------------------------------------
    if (self%time > self%rt_next) then
      print *, self%id, 'WARNING: RT time > rt_next', self%time, self%rt_next
    end if
    if (debug1(self)) &
      print '(i6,2f7.3,2x,a,f8.4)', self%id, self%position(3), self%mu, &
        'rt_t%update:  RT patch updated to t =', self%time
  end if
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Add upstream tasks to an nbor list
!> =============================================================================
SUBROUTINE init_nbors (self)
  class (rt_t):: self
  !.............................................................................
  integer:: i_omega
  type(rt_t), pointer:: omega
  type(link_t), pointer:: nbor, new, prev
  class(task_t), pointer:: nbtask, raytask
  integer:: n_nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbors')
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    n_nbor = 0
    nbor => self%link%nbor
    do while (associated(nbor))
      if (self%needs(omega,nbor%task)) then
        !-----------------------------------------------------------------------
        ! Add a new nbor, a new nbor%link, and pointers from link to task and
        ! from task to link.  The omega-task on nbor already has its own link.
        !-----------------------------------------------------------------------
        nbtask => nbor%task
        select type (nbtask)
        class is (rt_t)
          raytask => nbtask%omega(i_omega)
          select type (raytask)
          class is (rt_t)
            allocate (new)
            new%link => raytask%link
            new%task  => raytask
            if (verbose > 2) &
              print *, 'init_nbors: om=', omega%id, omega%link%task%id, 'nb=', &
            nbtask%id, nbtask%link%task%id, 'ray=', raytask%id, raytask%link%task%id
          end select
        end select
        !-----------------------------------------------------------------------
        ! Point to the new nbor,from the task%link, or from the previous nbor
        !-----------------------------------------------------------------------
        if (.not.associated(omega%link%nbor)) then
          omega%link%nbor => new
          if (verbose > 3) print *, omega%id, 'nbor => new'
        else
          if (verbose > 3) print *, omega%id, 'prev%next => new'
          prev%next => new
        end if
        if (verbose > 3) print *, omega%id, 'prev => new'
        prev => new
        n_nbor = n_nbor+1
      end if
      nbor => nbor%next
    end do
    if (verbose > 2) &
      print *, omega%id, 'rt_t%init_nbors: n_nbor=', n_nbor
  end do
  self%needs_nbors = .false.
  call trace%end
END SUBROUTINE init_nbors

!> =============================================================================
!> Compute the RT EOS quantities for the patch.  Delay checking in critical
!> region to here, so most threads never hit the critical region.
!> =============================================================================
SUBROUTINE rt_eos (self)
  class (rt_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e, tt
  integer:: m(4), ib, i3, jt(2)
  real:: pt(2), z
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  !$omp critical (eos_cr)
  ! -- interpolate to time = self%rt_next
  call self%time_interval (self%rt_next, jt, pt)
  if (verbose > 0) &
    print '(i7,3x,2i3,2f7.3,3x,f8.4,3x,5f8.4,3x,a)', &
      self%id, jt, pt, self%rt_next, self%t, 'eos'
  m = shape(self%rk)
  call allocate_scalars (m, d, e)
  d = self%mem(:,:,:,self%idx%d,jt(1),1)*pt(1) &
    + self%mem(:,:,:,self%idx%d,jt(2),1)*pt(2)
  e = self%mem(:,:,:,self%idx%e,jt(1),1)*pt(1) &
    + self%mem(:,:,:,self%idx%e,jt(2),1)*pt(2)
  if (verbose > 1) then
    call self%stats (d, 'rt_eos')
    call self%stats (e, 'rt_eos')
  end if
  tt => self%mem(:,:,:,self%idx%tt,self%it,1)
  call eos%lookup_table (m, d=d, e=e, rk=self%rk, src=self%src, tt=tt)
  if (verbose > 3) then
    print *,self%id,'rt_t%rt_eos:  d',jt,minval(d),maxval(d)
    print *,self%id,'rt_t%rt_eos:  e',jt,minval(e),maxval(e)
    print *,self%id,'rt_t%rt_eos: rk',jt,minval(self%rk),maxval(self%rk)
    print *,self%id,'rt_t%rt_eos: tt',jt,minval(tt),maxval(tt),scaling%temp
  end if
  call deallocate_scalars (d, e)
  self%eos_time = self%rt_next
  !$omp end critical (eos_cr)
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
      print *, self%id, 'rt_t%lower_BC', scaling%temp
  else if (self%mesh(3)%upper_boundary .and. self%mu < 0) then
    self%mem(:,:,self%mesh(3)%uo,:,self%it,1) = 0.0
    if (verbose > 2) &
      print *, self%id, 'rt_t%upper_BC'
  end if
  call trace%end
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Determine if a task is needed (is upstream)
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
  ! A task at a boundary with infalling intensity is always considered OK, even
  ! if the angle also requires wall intensities, because in a horizontally 
  ! periodic system nothing would ever start, otherwise.  The intensity coming
  ! in trough side walls is going to be slightly delayed in time, in some cases,
  ! depending on the order of task execution.
  !-----------------------------------------------------------------------------
  select type (task)
  class is (rt_t)
    rt_task => task
  end select
  lower = self%boundaries%is_set(bits%zl) .and. omega%mu > 0.0
  upper = self%boundaries%is_set(bits%zu) .and. omega%mu < 0.0
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and omega points inward, or the task has not RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  if (lower .or. upper .or. task%is_set(bits%no_rt)) then
    needs = .false.
  else
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
    !---------------------------------------------------------------------------
    ! Vertical rays depend only on exact alignment
    !---------------------------------------------------------------------------
    if (abs(omega%mu)==1.0) then
      needs = all(is1==is2)
    !---------------------------------------------------------------------------
    ! Define upstream as "not downstream in any direction"
    !---------------------------------------------------------------------------
    else
      needs = .not.any(is1*is2 < 0)
    end if
    !---------------------------------------------------------------------------
    ! If a horizontal ray component is axis aligned, exclude non-aligned nbors
    !---------------------------------------------------------------------------
    if (is2(1)==0 .and. is1(1)/=0) needs = .false.
    if (is2(2)==0 .and. is1(2)/=0) needs = .false.
    if (verbose > 2) then
      i = omega%i_omega
      select type(task)
      class is (rt_t)
      if (needs .or. verbose > 3) print &
        '("RT task",i7,4x,"needs task",i7,4x,"rel pos:",3f8.3,4x,"ray:",3f8.3,6i3,l3)', &
        omega%id, task%omega(i)%id, ds1, ds2, is1, is2, needs
      end select
    end if
  end if
END FUNCTION needs

!> =============================================================================
!> Return true if the patch is the horizontally closest to the give position
!> =============================================================================
FUNCTION debug1 (self) RESULT (close)
  class(rt_t):: self
  real(8):: position(3)
  logical:: close
  !-----------------------------------------------------------------------------
  close = maxval(abs(self%position(1:2)-dbg_position)/self%size(1:2)) < 0.55
  close = close .and. verbose>0
END FUNCTION debug1

END MODULE rt_mod
