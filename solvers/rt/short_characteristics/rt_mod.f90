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
!> =============================================================================
MODULE rt_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE index_mod
  USE radau_mod
  USE solver_mod
  USE task_mod
  USE link_mod
  USE eos_mod
  USE download_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from MHD, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(solver_t):: rt_t
    logical                    :: done=.true.
    logical                    :: needs_nbors=.true.
    type(rt_t), pointer        :: omega(:)
    integer                    :: n_omega=0
    integer                    :: i_omega=0
    integer                    :: n_bin
    integer                    :: offset_o(2)
    real                       :: offset_p(2)
    real                       :: mu, phi, w
    real(8), pointer           :: rt_t(:)
    real(8)                    :: rt_time=0d0
    real(8)                    :: rt_next=0d0
    real(8)                    :: eos_time=-1d0
    real(8)                    :: rt_grace_time=0.0
    real, pointer              :: w_bin(:)
    real, pointer              :: q(:,:,:,:)
    real, pointer              :: rk(:,:,:,:)
    real, pointer              :: src(:,:,:,:)
    real, pointer              :: layers(:,:,:,:,:)
  contains
    procedure:: init
    procedure:: update
    procedure:: test_d
    procedure:: init_nbors
    procedure:: upstream_is_ok
    procedure:: rt_eos
    procedure:: solve
    procedure:: integral_solver
    procedure:: prepare_layers
    procedure:: calculate_offsets
    procedure:: bilinear
    procedure:: needs
  end type
  real, parameter:: pi=acos(-1.0)
  real, save:: dtau_min=1e-3, dtau_max=30., dtau_dif=10.
  integer(8):: n_min=0, n_max=0, n_dif, n_ful=0
  integer, save:: verbose=1
  logical, save:: vertical_propagation=.false.
CONTAINS

!> =============================================================================
!> Read parameters
!> =============================================================================
SUBROUTINE init (self)
  class (rt_t):: self
  real(8):: rt_time=0.01_8
  integer, save:: n_mu=3, n_phi=4, n_bin=4
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat
  real, dimension(:), allocatable:: mu, w_mu
  type(rt_t), pointer:: omega
  class(task_t), pointer:: task
  class(rt_t),   pointer:: rt_task
  logical, save:: first_time=.true.
  namelist /rt_params/ n_mu, n_phi, n_bin, verbose, rt_time, dtau_min, dtau_max, &
                       dtau_dif, vertical_propagation
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
  n_omega      = 1 + n_phi*(n_mu-1)
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
  self%idx%tt = 9
  self%idx%qr = 10
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
  allocate (mu(n_mu),w_mu(n_mu))
  call radau%init (n_mu, mu, w_mu)
  !-----------------------------------------------------------------------------
  ! Allocate data common to all sub-tasks
  !-----------------------------------------------------------------------------
  allocate (self%rt_t(self%nt))
  allocate (self%rk (self%gn(1),self%gn(2),self%gn(3),n_bin))
  allocate (self%src(self%gn(1),self%gn(2),self%gn(3),n_bin))
  call io%bits_mem(storage_size(self%rk  ),product(shape(self%rk )),'rk')
  call io%bits_mem(storage_size(self%src ),product(shape(self%src)),'src')
  !-----------------------------------------------------------------------------
  ! Allocate n_omega copies of self, cloning self via self
  !-----------------------------------------------------------------------------
  allocate (self%omega(n_omega), source=self)
  !-----------------------------------------------------------------------------
  ! Initialize the directional solver connection links and data
  !-----------------------------------------------------------------------------
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    omega = self
    allocate (omega%link)
    omega%link%task => omega
    omega%nt      = 1
    omega%id      = self%id
    omega%i_omega = i_omega
    omega%nv      = n_bin
    omega%time    = -self%rt_time
    allocate (omega%mem(self%gn(1),self%gn(2),self%gn(3),n_bin,1,1))
    call io%bits_mem(storage_size(omega%mem),product(shape(omega%mem)),'q')
    omega%q => omega%mem(:,:,:,:,1,1)
    omega%q = 0.0
  end do
  !-----------------------------------------------------------------------------
  ! Assign angles
  !-----------------------------------------------------------------------------
  i_omega = 0
  do i_mu=1,n_mu-1
    do i_phi=1,n_phi
      do dir=-1,1,2
        i_omega = i_omega+1
        omega => self%omega(i_omega)
        omega%mu  = dir*mu(i_mu)
        omega%w   = w_mu(i_mu)
        omega%phi = pi*(modulo(real(2*(i_phi-1))/n_phi + (dir+1)/2.0, 2.0) - 1.0)
      end do
    end do
  end do
  do dir=-1,1,2
    i_omega = i_omega+1
    omega => self%omega(i_omega)
    omega%mu  = dir*mu(n_mu)
    omega%w   = w_mu(n_mu)
    omega%phi = 0.0
  end do
  !-----------------------------------------------------------------------------
  ! Check task angles
  !-----------------------------------------------------------------------------
  if (verbose > 0 .and. self%id == 1) then
    do i_omega=1,n_omega
      omega => self%omega(i_omega)
      print *, 'rt_t%init: i_omega, elevation, phi, mu =', &
        i_omega, nint(asin(omega%mu)*180./pi), nint(omega%phi*180./pi), &
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
  integer:: i_omega, i_bin, n_solve
  real:: n_tot
  real, dimension(:,:,:), pointer:: q
  class(rt_t), pointer:: omega
  !-----------------------------------------------------------------------------
  ! Call (only) the MHD solver, until it has passed the rt_next time, then
  ! work repeatedly with the RT solver, until that is ahead again
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%update')
                                                            call self%test_d (1)
  !-----------------------------------------------------------------------------
  ! Guard against round-off, and catch the MHD task at exactly self%rt_next
  !-----------------------------------------------------------------------------
  if (abs(self%time-self%rt_next) < 1e-6*self%dtime) then
    self%time = self%rt_next
    if (self%dt_fixed > 0.0) then
      if (verbose > 0) &
        print *, self%id, 'rt_t%update: arrived at =', self%rt_next
      self%dt_fixed = -1.0
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If self%rt_next is within reach, then set a constant time step to reach it,
  ! and if not, take a normal MHD step.
  !-----------------------------------------------------------------------------
  if (self%time < self%rt_next) then
    if (self%time + self%dtime > self%rt_next) then
      self%dt_fixed = self%rt_next-self%time
      if (verbose > 0) &
        print *, self%id, 'rt_t%update: setting dt_fixed =', self%dt_fixed
    end if
    call self%solver_t%update
  else
    !---------------------------------------------------------------------------
    ! Cancel the fixed time step, and prevent time slot rotation and time update, 
    ! as long as only RT is done
    !---------------------------------------------------------------------------
    self%dt_fixed = -1.0
    self%rotated = .true.
    !---------------------------------------------------------------------------
    ! Compute the opacity and source function for all bins
    !---------------------------------------------------------------------------
    if (self%eos_time < self%time) then
      call self%rt_eos
      self%eos_time = self%time
    end if
    self%done = .true.
    n_solve = 0
    !-------------------------------------------------------------------------
    ! Loop over directions, make sure nbor list exist, check it, and solve if OK
    !-------------------------------------------------------------------------
    do i_omega=1,self%n_omega
      omega => self%omega(i_omega)
      !-----------------------------------------------------------------------
      ! Check if time has advanced past the next RT time
      !-----------------------------------------------------------------------
      if (omega%time < self%rt_next) then
        if (self%upstream_is_ok(i_omega)) then
          call self%solve (omega)
          if (vertical_propagation) then
            omega%time = omega%time + self%rt_time*abs(omega%mu)
            omega%time = min(omega%time, self%rt_next)
          else
            omega%time = self%rt_next
          end if
          n_solve = n_solve+1
          if (verbose > 0) &
            print *, self%id, i_omega, 'rt_t%update: time =', &
              omega%time, self%rt_next
        else if (verbose > 1) then
            print *, self%id, i_omega, 'rt_t%update: upstream_not_OK' 
        end if
      end if
      self%done = self%done .and. (omega%time == self%rt_next)
      if (verbose > 1) &
        print *, self%id, i_omega, 'rt_t%update: time, rt_next =', omega%time, self%rt_next 
    end do
    !---------------------------------------------------------------------------
    ! If (and only if) all omega are done, sum up the net heating rates and
    ! update the time for the next RT solution
    !---------------------------------------------------------------------------
    if (self%done) then
      self%rt_next = self%rt_next + self%rt_time
      q => self%mem(:,:,:,self%idx%qr,self%it,1)
      q = 0.0
      do i_omega=1,self%n_omega
        omega => self%omega(i_omega)
        do i_bin=1,self%n_bin
          q = q + omega%mem(:,:,:,i_bin,1,1)*self%w_bin(i_bin)*omega%w
        end do
      end do
      if (self%id==1 .and. verbose > 0) then
        n_tot = max(n_min+n_max+n_dif+n_ful,1)
        print '(1x,a,4f8.3)', 'rt_t%update: n_min, n_max, n_dif, n_ful =', &
          n_min/n_tot, n_max/n_tot, n_dif/n_tot, n_ful/n_tot
      end if
    end if
    if (verbose > 1) &
      print *, self%id, 'rt_t%update: n_solve, rt_next, done =', &
        n_solve, self%rt_next, self%done
  end if
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
SUBROUTINE test_d (self, label)
  class(rt_t):: self
  real:: dmin, dmax
  integer:: label
  dmin = minval(self%mem(:,:,:,1,self%it,1))
  dmax = maxval(self%mem(:,:,:,1,self%it,1))
  if (dmin <= 0.0) then
    print *,'error at ',label
    stop
  else
    !print *, self%id, 'dmin,max', label, dmin, dmax
  end if
END SUBROUTINE test_d

!> =============================================================================
!> Check upstream relations. Upstream tasks are added to an nbor list
!> =============================================================================
SUBROUTINE init_nbors (self, i_omega)
  class (rt_t):: self
  integer:: i_omega
  !.............................................................................
  type(rt_t), pointer:: omega
  type(link_t), pointer:: nbor, new, prev
  class(task_t), pointer:: nbtask, raytask
  integer:: n_nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbors')
  nbor => self%link%nbor
  omega => self%omega(i_omega)
  n_nbor = 0
  do while (associated(nbor))
    if (self%needs(omega,nbor%task)) then
      !-------------------------------------------------------------------------
      ! Add a new nbor, a new nbor%link, and pointers from link to task and from
      ! task to link.  The omega-task on nbor already has its own link.
      !-------------------------------------------------------------------------
      nbtask => nbor%task
      select type (nbtask)
      class is (rt_t)
        raytask => nbtask%omega(i_omega)
        select type (raytask)
        class is (rt_t)
          allocate (new)
          new%link => raytask%link
          new%task  => raytask
        end select
      end select
      !-------------------------------------------------------------------------
      ! Point to the new nbor,from the task%link, or from the previous nbor
      !-------------------------------------------------------------------------
      if (.not.associated(omega%link%nbor)) then
        omega%link%nbor => new
        if (verbose>2) print *, self%id, omega%i_omega, 'nbor => new'
      else
        if (verbose>2) print *, self%id, omega%i_omega, 'prev%next => new'
        prev%next => new
      end if
      if (verbose>2) print *, self%id, omega%i_omega, 'prev => new'
      prev => new
      n_nbor = n_nbor+1
    end if
    nbor => nbor%next
  end do
  omega%needs_nbors = .false.
  if (verbose > 0) &
    print *, self%id, omega%i_omega, 'rt_t%init_nbors: n_nbor=', n_nbor
  call trace%end
END SUBROUTINE init_nbors

!> =============================================================================
!> Check upstream times, returning OK if they are ahead.  "self" is the i_omega
!> task in the local patch, and we compare its time with the time of the upstream
!> task.
!> =============================================================================
FUNCTION upstream_is_ok (self, i_omega)
  logical:: upstream_is_ok
  class(rt_t):: self
  integer:: i_omega
  !.............................................................................
  class(rt_t), pointer:: omega
  type(link_t), pointer:: link, nbor, next
  class(task_t), pointer:: rt, nbtask
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%upstream_is_ok')
  omega => self%omega(i_omega)
  if (omega%needs_nbors) &
    call self%init_nbors (i_omega)
  upstream_is_ok = .true.
  nbor => omega%link%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
      type is (rt_t)
      !-------------------------------------------------------------------------
      ! The upstream time must be larger than the local omega task time.
      ! Equal is not OK, to force initial start from the boundaries and in.
      !-------------------------------------------------------------------------
      if (nbtask%time <= omega%time) then
        !-----------------------------------------------------------------------
        ! Ignore purely horizontal dependency, since it would never be resolved
        !-----------------------------------------------------------------------
        if (nbtask%position(3) /= self%position(3)) then
          upstream_is_ok = .false.
          if (verbose > 1) &
            print *, self%id, i_omega, &
              'rt_t%upstream_is_ok: fail on id, mu, pos =', &
              nbtask%id, nbtask%mu, nbtask%time, omega%time
          exit
        end if
      end if
    end select
    nbor => nbor%next
  end do
  if (verbose > 1) &
    print *, self%id, i_omega, 'rt_t%upstream_is_ok:', upstream_is_ok
  if (upstream_is_ok) &
    call download%download_link (omega%link)
  call trace%end
END FUNCTION upstream_is_ok

!> =============================================================================
!> Compute the RT EOS quantities for the patch
!> =============================================================================
SUBROUTINE rt_eos (self)
  class (rt_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e
  integer:: m(4), ib, i3
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  d   => self%mem(:,:,:,self%idx%d,self%it,1)
  e   => self%mem(:,:,:,self%idx%e,self%it,1)
  m = shape(self%rk)
  call eos%lookup_table (m, d=d, e=e, rk=self%rk, src=self%src)
  if (verbose > 2) then
    print *,minval(d),maxval(d)
    print *,minval(e),maxval(e)
    print *,minval(self%rk),maxval(self%rk)
  end if
  call trace%end (itimer)
END SUBROUTINE rt_eos

!> =============================================================================
!> Integral short characteristics RT solver.  Flips between 1st and 2nd solution,
!> with the order of forward and revers unknown.
!> =============================================================================
SUBROUTINE solve (self, omega)
  class(rt_t):: self, omega
  integer:: i3, ib, m(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%solve', itimer=itimer)
  if (omega%done) then
    !---------------------------------------------------------------------------
    ! Allocate scratch space, and compute three depth layers of everything
    !---------------------------------------------------------------------------
    m = self%mesh%gn
    allocate (omega%layers(m(1),m(2),m(3),3,3))
    call io%bits_mem (storage_size(omega%layers), product(shape(omega%layers)),'layers')
  end if
  !-----------------------------------------------------------------------------
  ! Forward-index solution
  !-----------------------------------------------------------------------------
  if (omega%mu > 0) then
    do ib=1,self%n_bin
      call self%prepare_layers (omega, ib)
      do i3=self%mesh(3)%li,self%mesh(3)%ui
        call self%integral_solver (omega%layers(:,:,:,:,1)   , &   ! k = opacity
                                   omega%layers(:,:,:,:,2)   , &   ! s = source
                                   omega%layers(:,:,:,:,3)   , &   ! q = input
                                   omega%mem   (:,:,:,ib,1,1), &   ! q = output
                                   i3, 1, 3)
      end do
    end do
  else
  !-----------------------------------------------------------------------------
  ! Reverse-index solution
  !-----------------------------------------------------------------------------
    do ib=1,self%n_bin
      call self%prepare_layers (omega, ib)
      do i3=self%mesh(3)%ui,self%mesh(3)%li,-1
        call self%integral_solver (omega%layers(:,:,:,:,1)   , &   ! k = opacity
                                   omega%layers(:,:,:,:,2)   , &   ! s = source
                                   omega%layers(:,:,:,:,3)   , &   ! q = input
                                   omega%mem   (:,:,:,ib,1,1), &   ! q = output
                                   i3, 3, 1)
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Check if both directions done, and deallocate layers if so
  !-----------------------------------------------------------------------------
  if (omega%done) then
    omega%done = .false.
  else
    deallocate (omega%layers)
    nullify (omega%layers)
    call io%bits_mem (-storage_size(omega%layers), product(shape(omega%layers)),'layers')
    omega%done = .true.
  end if
  if (verbose > 1) then
    print *, self%id, omega%i_omega, 'rt_t%solve: mu, time, done', &
      omega%mu, omega%time, omega%done
  end if
  call trace%end (itimer)
END SUBROUTINE solve

!> =============================================================================
!> Integral solver for one layer.  The exact integral solution Q=I-S for a 2nd
!> order polynomial has coefficient for the upstream Q, and for the 1st and 2nd
!> derivative of S.  However, if one is doing a solution both forwards and 
!> backward, the term in the 1st derivative cancels.
!> =============================================================================
SUBROUTINE integral_solver (self, k, s, qin, qout, i3, l1, l3)
  class (rt_t):: self
  real, dimension(:,:,:,:) :: k, s, qin
  real, dimension(:,:,:)   :: qout
  integer:: i3, l1, l3, iv
  !.............................................................................
  integer:: l(3), u(3), i1, i2
  real:: dz1, dz2, dsdtau1, dsdtau2, d2sdtau2, ex0, ex2, dtaumin, dtaumax
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  real, dimension(self%mesh(1)%gn):: dtau1, dtau2
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%integral_solver', itimer=itimer)
  dz2 = self%mesh(3)%d/2.0
  l = self%mesh(1:3)%li
  u = self%mesh(1:3)%ui
  z => self%mesh(3)%r
  dz1      = abs(z(l1)-z(2))*0.5
  dz2      = abs(z(2)-z(l3))*0.5
  do i2=l(2),u(2)
    dtau1(l(1):u(1)) = dz1*(k(l(1):u(1),i2,i3,2)+k(l(1):u(1),i2,i3,l1))
    dtau2(l(1):u(1)) = dz2*(k(l(1):u(1),i2,i3,2)+k(l(1):u(1),i2,i3,l3))
    dtaumin = minval(dtau1(l(1):u(1)))
    dtaumax = maxval(dtau1(l(1):u(1)))
    if (dtaumin > dtau_max) then
      do i1=l(1),u(1)
        qout(i1,i2,i3) = 0.0
      end do
      !$omp atomic
      n_max = n_max+1
    else if (dtaumin > dtau_dif) then
      do i1=l(1),u(1)
        dsdtau1  =     (s(i1,i2,i3, 2) - s(i1,i2,i3,l1))/dtau1(i1)
        dsdtau2  =     (s(i1,i2,i3,l3) - s(i1,i2,i3, 2))/dtau2(i1)
        d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
        qout(i1,i2,i3) = d2sdtau2
      end do
      !$omp atomic
      n_dif = n_dif+1
    else if (dtaumax < dtau_min) then
      do i1=l(1),u(1)
        ex0      = 1.0-dtau1(i1)
        qout(i1,i2,i3) = qin(i1,i2,i3,l1)*ex0
      end do
      !$omp atomic
      n_min = n_min+1
    else 
      do i1=l(1),u(1)
        dsdtau1  =     (s(i1,i2,i3, 2) - s(i1,i2,i3,l1))/dtau1(i1)
        dsdtau2  =     (s(i1,i2,i3,l3) - s(i1,i2,i3, 2))/dtau2(i1)
        d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
        ex0      = exp(-dtau1(i1))
        ex2      = (1.0-ex0) - dtau1(i1)*ex0
        qout(i1,i2,i3) = qin(i1,i2,i3,l1)*ex0 &
                          + d2sdtau2     *ex2
      end do
      !$omp atomic
      n_ful = n_ful+1
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE integral_solver

!> =============================================================================
!> Prepare three layers; one centered on the mesh points, one above and below
!> The indices in the layers are (i1,i2,i3,il,iv,1:3), where i1-i3 are 3D coords,
!> ib is bin index, iv is variable index (k,s,q), and the last is layer index
!> =============================================================================
SUBROUTINE prepare_layers (self, omega, ib)
  class (rt_t):: self, omega
  integer:: ib
  integer:: iv, l(3), u(3), i3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%prepare_layers', itimer=itimer)
  l = self%mesh(1:3)%li
  u = self%mesh(1:3)%ui
  !-----------------------------------------------------------------------------
  ! Vertical boundary conditions; here we explicitly assume the vertical direction
  ! to be downwards, so the lower boundar index-wise os the top physical booundary,
  ! where I=0 => Q=-S, and the "upper" boundary is the physical "bottom" boundary,
  ! where I=S => Q=0.
  !-----------------------------------------------------------------------------
  if     (self%mesh(3)%lower_boundary .and. omega%mu > 0.0) then
    omega%mem(:,:,l(3)-1,:,1,1) = -self%src(:,:,l(3)-1,:)
    if (verbose > 1) &
      print *, self%id, omega%i_omega, 'rt_t%prepare_layers: using lower BC'
  else if(self%mesh(3)%upper_boundary .and. self%mu < 0.0) then
    omega%mem(:,:,u(3)+1,:,1,1) = 0.0
    if (verbose > 1) &
      print *, self%id, omega%i_omega, 'rt_t%prepare_layers: using upper BC'
  end if
  !-----------------------------------------------------------------------------
  do i3=2,self%mesh(3)%gn-1
    call self%calculate_offsets (i3, i3-1)
    call self%bilinear (self%rk (:,:,i3,ib)     , omega%layers(:,:,i3,1,1))
    call self%bilinear (self%src(:,:,i3,ib)     , omega%layers(:,:,i3,1,2))
    call self%bilinear (omega%mem(:,:,i3,ib,1,1), omega%layers(:,:,i3,1,3))
    !---------------------------------------------------------------------------
    omega%layers(:,:,i3,2,1) = self%rk (:,:,i3,ib)
    omega%layers(:,:,i3,2,2) = self%src(:,:,i3,ib)
    omega%layers(:,:,i3,2,3) = omega%mem(:,:,i3,ib,1,1)
    !---------------------------------------------------------------------------
    call self%calculate_offsets (i3, i3+1)
    call self%bilinear (self%rk (:,:,i3,ib)     , omega%layers(:,:,i3,3,1))
    call self%bilinear (self%src(:,:,i3,ib)     , omega%layers(:,:,i3,3,2))
    call self%bilinear (omega%mem(:,:,i3,ib,1,1), omega%layers(:,:,i3,3,3))
  end do
  call trace%end (itimer)
END SUBROUTINE prepare_layers

!> =============================================================================
!> Calculate offsets from one layer to the next
!> =============================================================================
SUBROUTINE calculate_offsets (self, i3a, i3b)
  class (rt_t):: self
  integer:: i3a, i3b
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  real:: dxdz, dydz, dz, sinmu
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%calculate_offsets', 2, itimer=itimer)
  z => self%mesh(3)%r
  dz = z(i3b)-z(i3a)
  sinmu = sqrt(1.0-self%mu**2)
  dxdz = sinmu*self%mesh(1)%d/dz
  dydz = sinmu*self%mesh(2)%d/dz
  self%offset_p = [dxdz, dydz]
  self%offset_o = floor(self%offset_p)
  self%offset_p = self%offset_p - self%offset_o
  call trace%end (itimer)
END SUBROUTINE calculate_offsets

!> =============================================================================
!> Bilinear interpolation in a plane
!> =============================================================================
SUBROUTINE bilinear (self, in, out)
  class (rt_t):: self
  real, dimension(:,:), intent(in)  :: in
  real, dimension(:,:), intent(out) :: out
  real:: p(2), q(2)
  integer:: i1, i2, i3, l(3), u(3), o(2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%bilinear', 2, itimer=itimer)
  l = self%mesh%li
  u = self%mesh%ui
  o = self%offset_o
  p = self%offset_p
  q = 1.0-p
  !-----------------------------------------------------------------------------
  ! Purely vertical case
  !-----------------------------------------------------------------------------
  if (all(abs(o+p) < 1e-5)) then
    out = in
  !-----------------------------------------------------------------------------
  ! Interpolation only in y
  !-----------------------------------------------------------------------------
  else if (abs(o(1)+p(1)) < 1e-5) then
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2) = &
        q(2)*in(i1,i2+o(2)) + p(2)*in(i1,i2+o(2)+1)
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Interpolation only in x
  !-----------------------------------------------------------------------------
  else if (abs(o(2)+p(2)) < 1e-5) then
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2) = &
        q(1)*in(i1+o(1),i2) + p(1)*in(i1+o(1)+1,i2)
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Interpolation in x and y
  !-----------------------------------------------------------------------------
  else
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    out(i1,i2) = &
        q(1)*q(2)*in(i1+o(1)  ,i2+o(2)  ) + &
        p(1)*q(2)*in(i1+o(1)+1,i2+o(2)  ) + &
        q(1)*p(2)*in(i1+o(1)+1,i2+o(2)+1) + &
        p(1)*p(2)*in(i1+o(1)+1,i2+o(2)+1)
  end do
  end do
  end if
  call trace%end (itimer)
END SUBROUTINE bilinear

!> =============================================================================
!> Determine if a task is needed (is upstream)
!> =============================================================================
FUNCTION needs (self, omega, task)
  logical:: needs
  class(rt_t):: self, omega
  class(task_t), pointer:: task
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
  lower = self%mesh(3)%lower_boundary .and. omega%mu > 0.0
  upper = self%mesh(3)%upper_boundary .and. omega%mu < 0.0
  if (lower .or. upper) then
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
    ! I a horizontal ray component is axis aligned, exclude non-aligned nbors
    !---------------------------------------------------------------------------
    if (is2(1)==0 .and. is1(1)/=0) needs = .false.
    if (is2(2)==0 .and. is1(2)/=0) needs = .false.
    if (verbose > 0) then
      if (needs .or. verbose > 1) print &
      '("RT task",2i5,4x,"needs patch",i5,4x,"rel pos:",3f8.3,4x,"ray:",3f8.3,6i3,l3)', &
      self%id, omega%i_omega, task%id, ds1, ds2, is1, is2, needs
     end if
  end if
END FUNCTION needs

END MODULE rt_mod
