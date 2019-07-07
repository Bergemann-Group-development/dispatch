!===============================================================================
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
    ! Assign the MHD patch and request for two more mem indices
    !---------------------------------------------------------------------------
    self%patch => gpatch%cast2gpatch(link%task)
    if (self%patch%idx%qr <1) then
      call self%patch%idx%request_index(self%patch%idx%qr, self%patch%nv)
      call self%patch%idx%request_index(self%patch%idx%tt, self%patch%nv)
    end if
  end if
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the object
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
    ! Allocate n_omega copies of self, cloning self via self
    !---------------------------------------------------------------------------
    if (n_omega > 0) allocate (self%omega(n_omega), source=self)
    self%i_omega = 0
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
    ! Initialize the copies of self
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
    !-----------------------------------------------------------------------------
    ! Check task angles
    !-----------------------------------------------------------------------------
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
!> Initialize rt sub-tasks
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
!===============================================================================
SUBROUTINE init_nbors (self)
  class (rt_t):: self
  !.............................................................................
  integer:: i_omega
  type(rt_t), pointer:: omega
  class(*), pointer:: nb_rt_task
  class(link_t), pointer:: nbor, new
  class(task_t), pointer:: nbtask, raytask
  class(gpatch_t), pointer:: patch, nb_patch
  integer:: n_nbor
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%init_nbors')
  patch => self%patch
  call self%clear(bits%init_nbors)
  do i_omega=1,self%n_omega
    omega => self%omega(i_omega)
    n_nbor = 0
    nbor => patch%link%nbor
    do while (associated(nbor))
      if (self%needs(omega,nbor%task)) then
        !-----------------------------------------------------------------------
        ! Add a new nbor, a new nbor%link, and pointers from link to task and
        ! from task to link.  The omega-task on nbor already has its own link.
        !-----------------------------------------------------------------------
        nbtask => nbor%task
        select type (nbtask)
        class is (gpatch_t)
          nb_patch   => nbtask
          nb_rt_task => nbtask%connect%rt
          select type(nb_rt_task)
          class is (rt_t)
            raytask => nb_rt_task%omega(i_omega)
            select type (raytask)
            class is (rt_t)
              if (raytask%on) then
                !-------------------------------------------------------------------
                ! Create an "upstrean type" nbor node pointing to the nbor sub-task
                !-------------------------------------------------------------------
                allocate (new)
                new%link => raytask%link
                new%task  => raytask
                ! -- add to the sub-tasks nbor list, in MPI rank order
                !call omega%link%add_nbor_by_rank (omega%link%nbors_by_level, new)
                call omega%link%add_nbor_by_rank (omega%link%nbor, new)
                !-------------------------------------------------------------------
                ! Create a "downstrean type" nbor node pointing back to omega
                !-------------------------------------------------------------------
                allocate (new)
                new%link => omega%link
                new%task  => omega%link%task
                !! -- add it to the nbor's nbor chain, in rank order
                !call raytask%link%add_nbor_by_rank (raytask%link%nbors_by_level, new)
                call raytask%link%add_nbor_by_rank (raytask%link%nbor, new)
                n_nbor = n_nbor+1
              end if
            end select
           end select
        end select
      end if
      nbor => nbor%next
    end do
    if (verbose > 2) &
      write(io%output,*) omega%id, 'rt_t%init_nbors: n_nbor=', n_nbor
  end do
  call trace%end
END SUBROUTINE init_nbors

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
!> Check if patch belongs to an RT box or is at its boundary.
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
    if (self%rt_extent(i) > 0.) then
      if ((urc(i) < self%rt_llc(i))) then
        call self%patch%set(bits%no_rt)
        self%on = .false.
      else if (llc(i) <= self%rt_llc(i)) then
        call self%boundaries%set (lbit(i))
      end if
      !
      if ((llc(i) > rt_urc(i))) then
        call self%patch%set(bits%no_rt)
        self%on = .false.
      else if (urc(i) >= rt_urc(i)) then
        call self%boundaries%set (ubit(i))
      end if
    end if
  end do
END SUBROUTINE setup_boundaries

!===============================================================================
!>
!===============================================================================
SUBROUTINE pre_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  !---------------------------------------------------------------------------
  ! If the main task is called with self%time > self%rt_next, then the
  ! scheduler has concluded that the pre-requisites are in order, and hence
  ! we can update self%q, which is later used (and reused) as the next
  ! radiative heating.
  !---------------------------------------------------------------------------
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
!>
!===============================================================================
SUBROUTINE post_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  !---------------------------------------------------------------------------
  ! If the main task is called with self%time > self%rt_next, then the
  ! scheduler has concluded that the pre-requisites are in order, and hence
  ! we can update self%q, which is later used (and reused) as the next
  ! radiative heating.
  !---------------------------------------------------------------------------
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
    write (io%output,*) 'DBG rt_t%udpate: time =', self%time, &
    self%is_set(bits%boundary), self%is_set(bits%virtual)
  !
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
  !
  if (verbose > 2) &
    write(io%output,'(i6,2f7.3,2x,a,f8.4)') &
      self%id, self%position(3), self%mu, &
      'rt_t%update:  RT patch updated to t =', self%time
  !
  if ((self%vertical_propagation).and.(self%i_omega>0)) then
    self%dtime = self%dt_fixed*abs(self%mu)
  else
    ! -- choose an RT time step that is at least the parent MHD time step
    self%dtime = self%dt_fixed
  end if
  call self%rotate()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Compute the RT EOS quantities for the patch, using self%eos_time to guard
!> againts duplicate calculations when multi-threading.  Check both before and
!> after entering a critical region, so most threads never hit the critical region.
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


SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(rt_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !.............................................................................
  class(rt_t), pointer:: rt_another
  class(task_t), pointer:: other
  !-----------------------------------------------------------------------------
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
  if(self%vertical_propagation) then
    if (any(self%position(1:2) /= rt_another%position(1:2))) &
      needed = .false.
  end if
  
  if (self%i_omega /= rt_another%i_omega) download = .false.
  
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

END MODULE rt_mod
