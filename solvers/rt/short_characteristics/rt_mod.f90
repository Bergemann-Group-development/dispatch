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
  integer:: verbose=1
CONTAINS

!> =============================================================================
!> Read parameters
!> =============================================================================
SUBROUTINE init (self)
  class (rt_t):: self
  real(8):: rt_time=0.01_8
  integer, save:: n_mu=3, n_phi=4, n_bin=4
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat
  real, dimension(:), allocatable:: mu, w_mu, phi
  type(rt_t), pointer:: omega
  class(task_t), pointer:: task
  class(rt_t),   pointer:: rt_task
  logical, save:: first_time=.true.
  namelist /rt_params/ n_mu, n_phi, n_bin, verbose, rt_time
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
  allocate (mu(n_mu),w_mu(n_mu), phi(n_phi))
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
    omega%nt = 1
    omega%id = i_omega
    omega%nv = n_bin
    omega%time = -1d0
    allocate (omega%mem(self%gn(1),self%gn(2),self%gn(3),n_bin,1,1))
    call io%bits_mem(storage_size(omega%mem),product(shape(omega%mem)),'q')
    omega%q => omega%mem(:,:,:,:,1,1)
    omega%q = 0.0
  end do
  !-----------------------------------------------------------------------------
  ! Assign angles
  !-----------------------------------------------------------------------------
  i_omega = 0
  do i_phi=1,n_phi
    do i_mu=1,n_mu-1
      do dir=-1,1,2
        i_omega = i_omega+1
        omega => self%omega(i_omega)
        omega%mu  = dir*mu(i_mu)
        omega%w   = w_mu(i_mu)
        omega%phi = phi(i_phi)
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
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Run through all angular solvers, and attempt to update them.  Solvers that
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
  integer:: i_omega, n_solve, iter
  real, dimension(:,:,:), pointer:: q
  class(rt_t), pointer:: omega
  integer, parameter:: max_iter=3
  !-----------------------------------------------------------------------------
  ! Call (only) the MHD solver, until it has passed the rt_next time, then
  ! work repeatedly with the RT solver, until that is ahead again
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%update')
call self%test_d (1)
  if (self%time < self%rt_next) then
    call self%solver_t%update
  else
    !---------------------------------------------------------------------------
    ! Compute the opacity and source function for all bins
    !---------------------------------------------------------------------------
    if (self%eos_time < self%time) then
      call self%rt_eos
      self%eos_time = self%time
call self%test_d (2)
    end if
    do iter=1,max_iter
      !-------------------------------------------------------------------------
      ! Loop over directions, make sure nbor list exist, check it, and solve if OK
      !-------------------------------------------------------------------------
      self%done = .true.
      n_solve = 0
      do i_omega=1,self%n_omega
        omega => self%omega(i_omega)
        !-----------------------------------------------------------------------
        ! Check if time has advanced past the next RT time
        !-----------------------------------------------------------------------
        if (omega%time < self%rt_next) then
          if (self%upstream_is_ok(i_omega)) then
            call self%solve (omega)
  call self%test_d (2+i_omega)
            n_solve = n_solve+1
          else if (verbose > 0) then
              print *, self%id, omega%id, 'rt_t%update: upstream_not_OK' 
          end if
        end if
        self%done = self%done .and. (omega%time == self%rt_next)
      end do
      if (verbose > 0) &
        print *, self%id, 'rt_t%update: iter, n_solve, done =', iter, n_solve, self%done
    end do
    !---------------------------------------------------------------------------
    ! If (and only if) all omega are done, update the time for the next RT solution
    !---------------------------------------------------------------------------
    if (self%done) then
      self%rt_next = self%rt_next + self%rt_time
      if (verbose > 0) &
        print *, self%id, omega%id, 'rt_t%update: rt_next =', self%rt_next, self%rt_time
    end if
  end if
  call trace%end
END SUBROUTINE update

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
        if (verbose>2) print *, self%id, omega%id, 'nbor => new'
      else
        if (verbose>2) print *, self%id, omega%id, 'prev%next => new'
        prev%next => new
      end if
      if (verbose>2) print *, self%id, omega%id, 'prev => new'
      prev => new
      n_nbor = n_nbor+1
    end if
    nbor => nbor%next
  end do
  omega%needs_nbors = .false.
  if (verbose > 0) &
    print *, self%id, omega%id, 'rt_t%init_nbors: n_nbor=', n_nbor
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
  logical:: lower, upper
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%upstream_is_ok')
  omega => self%omega(i_omega)
  if (omega%needs_nbors) &
    call self%init_nbors (i_omega)
  upstream_is_ok = .true.
  !-----------------------------------------------------------------------------
  ! A boundary with infalling intensity is always considered OK
  !-----------------------------------------------------------------------------
  lower = self%mesh(3)%lower_boundary .and. omega%mu > 0.0
  upper = self%mesh(3)%upper_boundary .and. omega%mu < 0.0
  if (.not.(lower .or. upper)) then
    nbor => omega%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type (nbtask)
        type is (rt_t)
        !---------------------------------------------------------------------
        ! The upstream time must be larger than the local omega task time.
        ! Equal is not OK, to force initial start from the boundaries and in.
        !---------------------------------------------------------------------
        if (nbtask%time <= omega%time) then
          upstream_is_ok = .false.
          if (verbose > 1) &
            print *, self%id, omega%id, 'rt_t%upstream_is_ok: fail on id, mu, pos =', &
              nbtask%id, nbtask%mu, nbtask%position
          exit
        end if
      end select
      nbor => nbor%next
    end do
  end if
  if (verbose > 1) &
    print *, self%id, omega%id, 'rt_t%upstream_is_ok:', upstream_is_ok
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
    omega%time = self%rt_next
    omega%done = .true.
  end if
  if (verbose > 1) then
    print *, self%id, omega%id, 'rt_t%solve: mu, time, done', &
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
  real:: dz1, dz2, dtau1, dtau2, dsdtau1, dsdtau2, d2sdtau2, ex0, ex2
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%integral_solver', itimer=itimer)
  dz2 = self%mesh(3)%d/2.0
  l = self%mesh(1:3)%li
  u = self%mesh(1:3)%ui
  z => self%mesh(3)%r
  dz1      = abs(z(l1)-z(2))*0.5
  dz2      = abs(z(2)-z(l3))*0.5
  do i2=l(2),u(2)
    do i1=l(1),u(2)
      dtau1    = dz1*(k(i1,i2,i3, 2) + k(i1,i2,i3,l1))
      dtau2    = dz2*(k(i1,i2,i3,l3) + k(i1,i2,i3, 2))
      dsdtau1  =     (s(i1,i2,i3, 2) - s(i1,i2,i3,l1))/dtau1
      dsdtau2  =     (s(i1,i2,i3,l3) - s(i1,i2,i3, 2))/dtau2
      d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1+dtau2)
      ex0      = exp(-dtau1)
      ex2      = (1.0-ex0) - dtau1*ex0
      qout(i1,i2,i3) = qin(i1,i2,i3,l1)*ex0 &
                        + d2sdtau2     *ex2
    end do
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
  ! Vertical boundary conditions
  !-----------------------------------------------------------------------------
  if     (self%mesh(3)%lower_boundary .and. omega%mu > 0.0) then
    self%mem(:,:,l(3)-1,:,1,1) = -self%src(:,:,l(3)-1,:)
    if (verbose > 1) &
      print *, self%id, omega%id, 'rt_t%prepare_layers: using lower BC'
  else if(self%mesh(3)%upper_boundary .and. self%mu < 0.0) then
    self%mem(:,:,u(3)+1,:,1,1) = 0.0
    if (verbose > 1) &
      print *, self%id, omega%id, 'rt_t%prepare_layers: using upper BC'
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
  real:: ds1(3), ds2(3), sinmu
  real, parameter:: eps=1e-4
  integer:: i, is1(3), is2(3)
  !-----------------------------------------------------------------------------
  needs = .true.
  ds1 = (task%position-omega%mesh%p)/omega%mesh%b
  sinmu = sqrt(1d0-omega%mu**2)
  ds2 = [sinmu*cos(omega%phi),sinmu*sin(omega%phi),omega%mu]
  is1 = 0
  is2 = 0
  do i=1,3
    if (ds1(i) > +eps) is1(i)=1
    if (ds2(i) > +eps) is2(i)=1
    if (ds1(i) < -eps) is1(i)=-1
    if (ds2(i) < -eps) is2(i)=-1
  end do
  !-----------------------------------------------------------------------------
  ! Define upstream as "not downstream in any direction"
  !-----------------------------------------------------------------------------
  if (abs(omega%mu)==1.0) then
    needs = all(is1==is2)
  else
    needs = .not.any(is1*is2 < 0)
  end if
  if (verbose>0) then
    if (needs) print &
    '("RT task",2i5,4x,"needs patch",i5,4x,"rel pos:",3f8.3,4x,"ray:",3f8.3,6i3)', &
    self%id, omega%id, task%id, ds1, ds2, is1, is2
   end if
END FUNCTION needs

END MODULE rt_mod
