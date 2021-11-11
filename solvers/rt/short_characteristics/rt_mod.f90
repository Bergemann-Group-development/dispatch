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

  USE rt_integral_mod
  USE radau_mod
  USE boundaries_mod
  USE eos_mod
  USE list_mod
  USE halo_mod
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
  USE units_mod
  USE scratch_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from gpatch, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: rt_t
    class(gpatch_t), pointer :: patch                                           ! pointer to MHD patch
    integer                  :: n_warmup=5                                      ! number of initial solutions
    integer                  :: n_bin=4                                         ! number of frequency bins
    integer                  :: n_omega=0                                       ! number of space-angles
    integer                  :: i_omega=0                                       ! space-angle index
    integer                  :: axis=0                                          ! neares axis
    real, pointer            :: w_bin(:)                                        ! frequency weight
    class(rt_t), pointer     :: omeg  => null()                                 ! different RT directions
    class(rt_t), pointer     :: main  => null()                                 ! point to main RT task
    class(rt_t), pointer     :: next  => null()                                 ! point to next omega task
    class(rt_t), pointer     :: prev  => null()                                 ! point to prev omega task
    class(rt_t), pointer     :: tail  => null()                                 ! point to prev omega task
    real                     :: cdtd  = 1.0
    real                     :: mu, rt_phi, w                                   ! mu, phi, weights
    real                     :: rt_hat(3)                                       ! RT direction
    real                     :: epsilon =1.0                                    ! scattering coefficient
    real                     :: rt_freq
    real(8)                  :: rt_llc(3)=0d0                                   ! lower left corner
    real(8)                  :: rt_urc(3)=0d0                                   ! upper right corner
    real(8)                  :: rt_timestep=0.01_8                              ! constant time step
    real(8)                  :: qr_time                                         ! time when qr was last collected
    real(8)                  :: warmup_end_time=0d0                             ! time when RT time should start varying
    integer                  :: n_mu=3, n_phi=4                                 ! no. of angles
    integer                  :: irk(2), isrc(2)                                 ! lower/upper indices of rk and src in mem
    integer                  :: itt, iqr                                        ! temperature and qr memory indices
!
    integer:: master_id
!
    logical                  :: on=.false.                                      ! allows local turn-off
    logical                  :: vertical_propagation=.false. ! see code
    logical                  :: warmup_done=.false.
    logical                  :: pre_initialized_rt=.false.
    logical                  :: initialized=.false.
    logical                  :: first_eos  =.true.
    logical                  :: in_task_list = .false.
    real, pointer            :: q(:,:,:,:)                                      ! pointer to net radiative heating
    real, pointer            :: tt(:,:,:)                                       ! EOS temperature
    real, pointer            :: rk(:,:,:,:)                                     ! EOS rho*kappa
    real, pointer            :: src(:,:,:,:)                                    ! EOS source function
    real, pointer            :: qr(:,:,:)                                       ! net radiative heating
  contains
    procedure:: dealloc                               ! clean deallocate RT
    procedure:: pre_init                              ! read parameters
    procedure:: init                                  ! init RT task
    procedure:: init_task_list                        ! add to task list
    procedure:: needs                                 ! check upstream
    procedure:: is_upstream_of                        ! check upstream
    procedure:: pre_update                            ! pre MHD
    procedure:: post_update                           ! post MHD
    procedure:: update                                ! RT task update
    procedure:: rt_eos                                ! EOS for RT
    procedure:: solve                                 ! RT solver
    procedure:: check_boundaries                      ! BCs
    procedure:: init_bdries                           ! properlyy initialize boundaries
    procedure:: task_info                             ! formatted info print
    procedure:: nbor_info                             ! formatted info print
    procedure:: nbor_relations                        ! dependency relationships between different tasks from RT point of view
    procedure:: output                                ! intercept output
  end type
  integer, save:: verbose=0                           ! verbosity
  real,    save:: rt_freq=0                           ! how often RT should be updated WRT MHD -- 0=unlimited
  logical,save:: dnload_mhd=.true.
  type(rt_t), public:: rt                             ! public instance
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
recursive SUBROUTINE pre_init (self)
  class (rt_t), target:: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  class(link_t), pointer::link
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n
  integer, save:: n_mu=3, n_phi=4, n_bin=4, nt=1, n_warmup=3, ng(3)=1, loc(1)
  real, save   :: courant=0.5, cdtd=1.0
  real, dimension(:), allocatable:: mu, w_mu
  real,    save:: rt_grace=0.01
  real(8), save:: rt_llc(3)=-huge(1d0)
  real(8), save:: rt_urc(3)=+huge(1d0)
  real(8), save:: rt_timestep=0.01_8
  logical, save:: vertical_propagation=.false.
  logical, save:: first_time=.true.
  logical, save:: on=.true.
  !-----------------------------------------------------------------------------
  namelist /sc_rt_params/ on, verbose, nt, ng, rt_llc, rt_urc, rt_timestep, &
    n_mu, n_phi, n_bin, rt_grace, vertical_propagation, n_warmup, courant, cdtd, rt_freq, dnload_mhd
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%pre_init')
  call io%lock%set ('rt_t%pre_init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, sc_rt_params, iostat=iostat)
    call io%namelist_warning ('sc_rt_params', iostat)
    write (io%namelist, sc_rt_params)
  end if
  if (self%pre_initialized_rt) call io%abort ('double rt pre_init')
  self%pre_initialized_rt = .true.
  if (self%pre_initialized) print *, 'pre_init issued somewhere'
  call io%lock%unset ('rt_t%pre_init')
  ! ----------------------------------------------------------------------------
  ! Mute a warning about bad nw, nv, and ng - these are set later.
  ! nw and ng - just below,
  ! nv - in the init procedure, where we will know, how many RT directions there
  ! will be. 
  ! ----------------------------------------------------------------------------
  self%nw = self%patch%nw
  self%nv = self%patch%nv
  self%ng = self%patch%ng
  call self%gpatch_t%pre_init()
  self%on = on
  self%nt = nt
  self%ng = ng
  self%nw = 1                                                                   ! no time derivatives
  self%n_warmup = n_warmup
  self%n_bin = n_bin
  self%grace = rt_grace
  self%rt_llc = rt_llc
  self%rt_urc = rt_urc
  self%rt_timestep = rt_timestep
  self%vertical_propagation = vertical_propagation
  self%n_mu  = n_mu
  self%n_phi = n_phi
  self%courant = courant
  self%cdtd    = cdtd
  self%rt_freq = rt_freq
  self%in_task_list = .false.
  ! ----------------------------------------------------------------------------
  ! relevant synchronization items for cartesian setup
  ! ----------------------------------------------------------------------------
  self%ipos       = self%patch%ipos
  self%rank       = self%patch%rank
  self%box        = self%patch%box
  self%size       = self%patch%size
  self%position   = self%patch%position
  self%ds         = self%patch%ds
  self%origin     = self%patch%origin
  self%llc_cart   = self%patch%llc_cart
  self%centre_nat = self%patch%centre_nat
  self%llc_nat    = self%patch%llc_nat
  self%status     = self%patch%status
  self%periodic   = self%patch%periodic
  self%on = on
  !-----------------------------------------------------------------------------
  ! Check if this MHD patch needs RT
  !-----------------------------------------------------------------------------
  call self%check_boundaries
  if (self%on) then 
    call self%task_t%init_unique (same=.true.)                                  ! obtain an unique ID
    !---------------------------------------------------------------------------
    self%type    = 'rt_t'
    if (self%i_omega==0) then
      self%kind    = 'rt_sc0'
    else
      self%kind    = 'rt_sc'
    end if
    !
    if (self%i_omega==0) then
      self%main => self         ! convenience
      call rt_integral%init
      call self%patch%aux%init
      !-------------------------------------------------------------------------
    ! Number of ray directions
      !-------------------------------------------------------------------------
      n_omega = merge (1+self%n_phi*(self%n_mu-1), 0, self%n_mu > 0)
      n_omega = 2*n_omega
      self%n_omega = n_omega
      !-------------------------------------------------------------------------
      ! Assign angles and integration weight and allocate the required RT 
      ! solvers for all angles
      !-------------------------------------------------------------------------
      if (self%n_mu > 0) then
        allocate (mu(self%n_mu),w_mu(self%n_mu))
        call radau%init (self%n_mu, mu, w_mu)
        !-----------------------------------------------------------------------
        ! Assign angles and initialize clones
        !-----------------------------------------------------------------------
        i_omega = 0
        do i_mu=1,self%n_mu
          n = self%n_phi
          if (i_mu == self%n_mu) n = 1
          do i_phi=1,n
            do dir=-1,1,2
              i_omega = i_omega+1
              allocate(omega)
              omega%main => self
              omega%master_id = self%id
              if (associated(self%omeg)) then
                omega%prev => self%tail
                omega%prev%next => omega
                self%tail => omega
              else
                self%omeg => omega
                self%tail => omega
              end if
              omega%mu      = dir*mu(i_mu)
              omega%w       = w_mu(i_mu)*math%pi4/n/2.
              omega%patch   => self%patch
              omega%i_omega = i_omega
              omega%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + &
                (dir+1)/2.0, 2.0) - 1.0)
              if (i_mu == self%n_mu) omega%rt_phi = 0.0
              call omega%pre_init
            end do
          end do
        end do
      end if
    else
      ! ------------------------------------------------------------------------
      ! Determine closest axis direction for omega/=0
      ! ------------------------------------------------------------------------
      self%rt_hat(1) = sin(acos(self%mu))*cos(self%rt_phi)                      ! sin(theta)*cos(phi), dx
      self%rt_hat(2) = sin(acos(self%mu))*sin(self%rt_phi)                      ! sin(theta)*sin(phi), dy
      self%rt_hat(3) = self%mu                                                  ! cos(theta), dz
      loc = maxloc(abs(self%rt_hat))
      self%axis = loc(1)                                                        ! direction of max propagation
    end if
  end if
  !call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the rt_t data types.  To setup for n_warmup initial solutions, we
!> set the time rt_t times = -n_warmup * rt_t%timestep, making sure to call
!> rt_t%rt_eos for the initial MHD state. This subroutine is recursive, as
!> omega0 initializes all the omegaX. 
!===============================================================================
SUBROUTINE init (self)
  class (rt_t):: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
  class (rt_t),     pointer :: omega
  class (link_t),   pointer :: link
  class (*),        pointer :: dummy
  integer:: n_omega, i_omega, dir, i_mu, i_phi, iostat, n
  ! ----------------------------------------------------------------------------
  if (self%initialized) then
    print *, self%id, self%i_omega
    call io%abort ('double rt init')
  end if 
  call trace%begin ('rt_t%init')
  self%initialized = .true.
  p => self%patch
  !
  ! ----------------------------------------------------------------------------
  ! omega0 stores EOS data, plus the collected heating rates. omegaX stores just
  ! the directional heating rates.
  ! ----------------------------------------------------------------------------
  if (self%i_omega==0) then
    self%nv = self%n_bin*2+2                                                    ! nbin*[rk,src]+tt+qr
    self%irk(1)  = 1
    self%irk(2)  = self%irk(1) +self%n_bin-1
    self%isrc(1) = self%irk(2) +1
    self%isrc(2) = self%isrc(1)+self%n_bin-1
    self%itt     = self%isrc(2)+1
    self%iqr     = self%itt+1
  else
    self%nv = self%n_bin
  end if
  self%mv = self%nv
  !
  write (io_unit%mpi,'(3i4,3f8.3,i4,a)') p%id, self%id,self%i_omega,p%position,p%rank,'pID, sID, omega, pos, rank'
  !
  if (verbose == -1) &
    write (io_unit%log,*) 'rt_t%init: id, mhd%id i_omega =', self%id, p%id, self%i_omega
  !
  call self%gpatch_t%init
  self%mem = 0.0
  if (self%i_omega==0) then
    self%mem(:,:,:,self%irk(1):self%irk(2),:,:)=tiny(1.)
    self%mem(:,:,:,self%isrc(1):self%isrc(2),:,:)=tiny(1.)
    self%mem(:,:,:,self%itt,:,:)=tiny(1.)
  end if
  write (io_unit%log,*) 'call self%patch_t%setup'
  !-----------------------------------------------------------------------------
  ! Synchronize with MHD patch
  ! ----------------------------------------------------------------------------
  self%level      = p%level
  self%restart    = p%restart
  self%warmup_end_time = p%time
  self%t(self%it) = self%time
  ! -- add an `eps` to avoid round-off errors
  self%out_next = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time
  self%dt_fixed = self%rt_timestep
  self%dtime    = self%rt_timestep
  self%dt       = self%rt_timestep
  self%time     = p%time - (self%n_warmup-self%grace)*self%dtime
  self%t        = self%time
  self%qr_time  = self%time - 2.0*self%grace*self%dtime
  ! ----------------------------------------------------------------------------
  ! Set omegaX time to be marginally after omega0, so EOS would be done before
  ! omegaX
  ! ----------------------------------------------------------------------------
  if (self%i_omega/=0) then
    self%time     = self%time + 2.0*self%grace*self%dtime
  end if
  self%iout     = self%restart+1
  self%verbose  = verbose
  ! ----------------------------------------------------------------------------
  ! Do not check density or validate with it. This is rather important to be set
  ! ----------------------------------------------------------------------------
  call self%set (bits%no_density)
  ! ----------------------------------------------------------------------------
  ! Allocate and set unsigned and pervolume. qr, src and tt are definite
  ! positive, whereas Qr can be positive or negative. omegaX stores just Qr.
  ! Although quantities are generally per volume, but density is not stored
  ! locally, thus this flag must not be set.
  ! ----------------------------------------------------------------------------
  allocate (self%unsigned(self%nv), self%pervolume(self%nv))
  self%pervolume = .false.
  if (self%i_omega==0) then
    self%unsigned  = .false.
    self%unsigned(self%nv) = .false.
    self%pervolume(self%irk(1):self%irk(2))   = .true.
    self%pervolume(self%isrc(1):self%isrc(2)) = .true.
    self%pervolume(self%itt)   = .true.
  else
    self%unsigned  = .false.
  end if
  ! ----------------------------------------------------------------------------
  ! Init boundary conditions
  ! ----------------------------------------------------------------------------  
  call self%init_bdries
  !
  if (self%on) then
    !---------------------------------------------------------------------------
    ! Integration weights (w_bin is baked into the source function in this case)
    !---------------------------------------------------------------------------
    allocate (self%w_bin (self%n_bin))
    self%w_bin = 1.0
  else
    !---------------------------------------------------------------------------
    ! Deallocate the self%mem for this task, but keep the task, for simplicity
    !---------------------------------------------------------------------------
    call self%dealloc()
    self%n_omega = 0
  end if
  call trace%end()
END SUBROUTINE init


!===============================================================================
!> The rt_t%update procedure gets called directly by the task dispatcher, when
!> it has concluded that the task is ready to update.  An update of the main RT
!> task consists of looking up the EOS, while the
!> sub-task updates perform the actual RT solutions along rays.
!===============================================================================
SUBROUTINE update (self)
  class(rt_t):: self
  !.............................................................................
  class(rt_t), pointer :: omega
  integer:: i_omega, i_bin
  real, dimension(:,:,:), pointer:: q, mem3, mem4
  real, dimension(:,:,:,:), pointer:: mem1, mem2
  real(8):: time
  real::pt(2)
  integer:: jt(2), ir(2),is(2), itt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) &
    return
  !
  call trace%begin ('rt_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! The omega0 task is responsible for calling rt_t%rt_eos to compute the RT EOS
  ! data. 
  !-----------------------------------------------------------------------------
  if (self%i_omega==0) then
    call self%rt_eos()
    call self%patch_t%courant_condition()
  !-----------------------------------------------------------------------------
  ! RT sub-task, solving for the radiation field in a given direction
  !-----------------------------------------------------------------------------
  else
    ! --------------------------------------------------------------------------
    ! Allocate RT scratches
    ! --------------------------------------------------------------------------
    allocate(self%rk (self%gn(1),self%gn(2),self%gn(3),self%n_bin))
    allocate(self%src(self%gn(1),self%gn(2),self%gn(3),self%n_bin))
    allocate(self%tt (self%gn(1),self%gn(2),self%gn(3)))
    ! --------------------------------------------------------------------------
    ! get the right time slices for rk, src and temperature
    ! --------------------------------------------------------------------------
    call self%main%time_interval (self%time, jt, pt, extrap=.true.)
    ! --------------------------------------------------------------------------
    ! get the right time slices for rk, src and temperature
    ! --------------------------------------------------------------------------
    ir  = self%main%irk
    is  = self%main%isrc
    itt = self%main%itt
    call self%main%mem_lock(jt(1))%set ('rt_eos_collect')
    if (jt(2) /= jt(1)) &
      call self%main%mem_lock(jt(2))%set ('rt_eos_collect')
    mem1 => self%main%mem(:,:,:,ir(1):ir(2),jt(1),1)
    mem2 => self%main%mem(:,:,:,ir(1):ir(2),jt(2),1)
    self%rk  = exp(log(mem1)*pt(1) + log(mem2)*pt(2))
    mem1 => self%main%mem(:,:,:,is(1):is(2),jt(1),1)
    mem2 => self%main%mem(:,:,:,is(1):is(2),jt(2),1)
    self%src = exp(log(mem1)*pt(1) + log(mem2)*pt(2))
    mem3 => self%main%mem(:,:,:,itt,jt(1),1)
    mem4 => self%main%mem(:,:,:,itt,jt(2),1)
    self%tt  = exp(log(mem3)*pt(1) + log(mem4)*pt(2))
    call self%main%mem_lock(jt(1))%unset ('rt_eos_collect')
    if (jt(2) /= jt(1)) &
      call self%main%mem_lock(jt(2))%unset ('rt_eos_collect')
    ! --------------------------------------------------------------------------
    ! Store the solution in a new time slot
    ! -------------------------------------------------------------------------- 
    self%q => self%mem(:,:,:,:,self%new,1)
    self%q = self%mem(:,:,:,:,self%it,1)
    call self%solve ()
    ! --------------------------------------------------------------------------
    ! Estimate dtime
    ! --------------------------------------------------------------------------
    self%u_max     = self%main%u_max
    self%dt_fixed = self%main%dtime
    !print *, 'minmax q', minval(self%q), maxval(self%q), self%time, self%main%time
    !print *, 'minmax rk', minval(self%rk), maxval(self%rk)
    !print *, 'minmax src', minval(self%src), maxval(self%src)
    ! --------------------------------------------------------------------------
    ! If there is a limit, how much RT can advance without MHD update
    ! --------------------------------------------------------------------------

    call self%patch_t%courant_condition
write (io_unit%log,'(a,2i6,1p,3e11.2,l4,e11.2)') &                          !DEBUG
'upd: thread,id,time,dtime,rotated =', &                                    !DEBUG
omp%thread, self%id, self%time, self%dtime, self%dt_fixed, self%rotated, &  !DEBUG
self%u_max                                                                  !DEBUG
     ! --------------------------------------------------------------------------
    ! cleanup
    ! --------------------------------------------------------------------------
    deallocate (self%rk, self%src, self%tt)
  end if
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> This is called only for MHD tasks, just before they update, with self = main,
!> so this is suitable place to add the resulting net radiative heating to the
!> heating per unit volume.
!===============================================================================
SUBROUTINE pre_update (self)
  class(rt_t), target:: self
  !.............................................................................
  class(rt_t), pointer:: omega
  integer:: l(3), u(3), n(3), i_omega, i_bin, jt(2), x, y, z, gn(3)
  real   :: pt(2)
  real(8):: time
  integer, save:: itimer=0
  real, pointer:: qr(:,:,:), rk(:,:,:,:), aux(:,:,:)
  real, dimension(:,:,:), pointer:: qr1, qr2
  !-----------------------------------------------------------------------------
  if (.not. self%on) then
    gn = self%patch%mesh%gn
    call scratch%allocate (scratch%heating_per_unit_volume, gn)
    aux => scratch%heating_per_unit_volume
    call self%patch%aux%register ('qr', aux)
    return
  end if
  call trace%begin ('rt_t%pre_update', itimer=itimer)
  ! ----------------------------------------------------------------------------
  ! Check, if new qr has to be collected.
  ! ----------------------------------------------------------------------------
  l  = self%mesh%lb
  u  = self%mesh%ub
  n  = self%mesh%gn
  allocate(rk(n(1),n(2),n(3),self%n_bin))
  allocate(qr(n(1),n(2),n(3)))
  qr=0.0
  call self%time_interval(self%patch%time, jt, pt, extrap=.true.)
  call self%mem_lock(jt(1))%set ('rt_collect')
  if (jt(2) /= jt(1)) &
    call self%mem_lock(jt(2))%set ('rt_collect')
  rk = exp(log(self%mem(:,:,:,self%irk(1):self%irk(2),jt(1),1))*pt(1) + &
           log(self%mem(:,:,:,self%irk(1):self%irk(2),jt(2),1))*pt(2))
  call self%mem_lock(jt(1))%unset ('rt_collect')
  if (jt(2) /= jt(1)) &
    call self%mem_lock(jt(2))%unset ('rt_collect')
  omega => self%omeg
  do while(associated(omega))
  call omega%time_interval(self%patch%time,jt,pt, extrap=.true.)
  call omega%mem_lock(jt(1))%set ('rt_assembly')
  if (jt(2) /= jt(1)) &
      call omega%mem_lock(jt(2))%set ('rt_assembly')
    do i_bin=1,self%n_bin
      qr1 => omega%mem(:,:,:,i_bin,jt(1),1)
      qr2 => omega%mem(:,:,:,i_bin,jt(2),1)
      do z= l(3),  u(3)
        do y= l(2), u(2)
          do x= l(1), u(1)
            qr(x,y,z) = qr(x,y,z) + rk(x,y,z,i_bin)*self%w_bin(i_bin)*omega%w* &
                (qr1(x,y,z)*pt(1) + qr2(x,y,z)*pt(2))
          end do
        end do
      end do 
    end do
    call omega%mem_lock(jt(1))%unset ('rt_assembly')
    if (jt(2) /= jt(1)) &
      call omega%mem_lock(jt(2))%unset ('rt_assembly')
      omega => omega%next
    end do
  ! ----------------------------------------------------------------------------
  ! Interpolate qr to MHD patch time and apply to heating_per_unit_volume
  !-----------------------------------------------------------------------------
  l = self%patch%mesh%li - self%mesh%ng
  u = self%patch%mesh%ui + self%mesh%ng
  gn = self%patch%mesh%gn
  call scratch%allocate (scratch%heating_per_unit_volume, gn)
  n = shape(scratch%heating_per_unit_volume)
  scratch%heating_per_unit_volume(l(1):u(1),l(2):u(2),l(3):u(3)) = &
  scratch%heating_per_unit_volume(l(1):u(1),l(2):u(2),l(3):u(3)) + qr
  !print *, 'minmax collect qr', minval(qr), maxval(qr), self%time, self%dtime
  !print *, 'minmax collect rk', minval(rk), maxval(rk), self%time, self%dtime
  !-----------------------------------------------------------------------------
  aux => scratch%heating_per_unit_volume
  !print *, 'minmax q', minval(aux(l(1):u(1),l(2):u(2),l(3):u(3))), maxval(aux(l(1):u(1),l(2):u(2),l(3):u(3))), self%patch%id
  call self%patch%aux%register ('qr', aux)
  ! ----------------------------------------------------------------------------
  ! cleanup
  ! ----------------------------------------------------------------------------
  deallocate(rk, qr)
  ! ----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> 
!===============================================================================
SUBROUTINE post_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(gpatch_t), pointer:: mhd
  class(rt_t), pointer:: omega
  integer:: i
  !-----------------------------------------------------------------------------

END SUBROUTINE post_update

!===============================================================================
!> Compute the RT EOS quantities for the patch.
!===============================================================================
SUBROUTINE rt_eos (self)
  class (rt_t):: self
  !.............................................................................
  class (gpatch_t), pointer :: patch
  real, dimension(:,:,:), pointer:: lnd, pgas
  real, dimension(:,:,:), allocatable:: e1,e2, fmax, lnd1, lnd2
  integer:: n(3), ib, loc(3), pl(3), pu(3),l(3), u(3), it, jt(2), pgn(3)
  real:: u_max, stefan, fmaxv, pt(2), dsmin
  real, pointer:: eth(:,:,:), tt(:,:,:), rk(:,:,:,:), src(:,:,:,:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  patch => self%patch
  ! ----------------------------------------------------------------------------
  ! In the initial EOS lookup, nbors have not been downloaded, so everything has
  ! to be looked up locally
  ! ----------------------------------------------------------------------------
  it=self%new
  if (self%first_eos.or..not.self%warmup_done) then
    if (.not.dnload_mhd) self%first_eos=.false.
    pl = patch%mesh%li-self%mesh%ng
    pu = patch%mesh%ui+self%mesh%ng
    l  = self%mesh%li-self%mesh%ng
    u  = self%mesh%ui+self%mesh%ng
    n  = merge(self%n+2*self%mesh%ng,self%n,self%n>1)
    call patch%time_interval (self%time, jt, pt, extrap=.true.)
  else
    pl = patch%mesh%li
    pu = patch%mesh%ui
    l  = self%mesh%li
    u  = self%mesh%ui
    n  = self%n
    call patch%time_interval (self%time, jt, pt, extrap=.true.)
  end if
  pgn = patch%mesh%gn
  if (dnload_mhd) call self%patch%dnload()
  allocate (pgas(n(1),n(2),n(3)))
  allocate (fmax(n(1),n(2),n(3)))
  allocate (lnd(n(1),n(2),n(3)))
  allocate (lnd1(n(1),n(2),n(3)))
  allocate (lnd2(n(1),n(2),n(3)))
  allocate (eth(n(1),n(2),n(3)))
  allocate (e1(pgn(1),pgn(2),pgn(3)))
  allocate (e2(pgn(1),pgn(2),pgn(3)))
  !-----------------------------------------------------------------------------
  call patch%mem_lock(jt(1))%set ('rt_eos')
    if (jt(2) /= jt(1)) &
      call patch%mem_lock(jt(2))%set ('rt_eos')
  lnd1 = log(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(1),1))
  !print *, minval(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1)), &
  !         maxval(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1)), &
  !         minloc(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1)), &
  !         maxloc(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1))
  !print *, jt, pt, self%time, patch%time
  !print *, 'pl, pu', pl, pu
  lnd2 = log(patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1))
  lnd = lnd1*pt(1) + lnd2*pt(2)
  call patch%to_thermal(jt(1),eth=e1, kind=patch%kind)
  call patch%to_thermal(jt(2),eth=e2, kind=patch%kind)
  eth = exp(log(e1(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3)))*pt(1) + &
            log(e2(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3)))*pt(2))
  call patch%mem_lock(jt(1))%unset ('rt_eos')
    if (jt(2) /= jt(1)) &
      call patch%mem_lock(jt(2))%unset ('rt_eos')
  rk  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), &
         self%irk(1):self%irk(2),it,1)
  src => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), &
         self%isrc(1):self%isrc(2),it,1)
  tt  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%itt,it,1)
  if (verbose > 2) then
    call self%stats (lnd, 'lnd')
    call self%stats (eth, 'e')
  end if
  !call patch%to_thermal_begin (patch%it, kind=patch%kind) 
  call eos%lookup_table (n, lnd=lnd, e=eth, rk=rk, src=src, &
                         tt=tt, pg=pgas, label='rt', id=self%id)
  !call patch%to_thermal_end (kind=patch%kind)
  self%rk  => self%mem(:,:,:,self%irk(1):self%irk(2),it,1)
  self%src => self%mem(:,:,:,self%isrc(1):self%isrc(2),it,1)
  ! ----------------------------------------------------------------------------
  ! Compute radiative diffusion speed, and set default omega%u_max
  ! ----------------------------------------------------------------------------
  self%u_max = 0.0
  dsmin = minval(self%ds,self%n>1)
  stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
  do ib=1, eos%n_lambda
    fmax = math%pi2*src(:,:,:,ib)/pgas*(rk(:,:,:,ib)*dsmin)
    !fmax = (pgas/d)**4/pgas*(self%rk(:,:,:,ib)*minval(self%ds))
    fmax = fmax/(1.+(rk(:,:,:,ib)*dsmin)**2)
    fmaxv = maxval(fmax)
    u_max = 16./3.*fmaxv*self%courant/self%cdtd*self%epsilon
    self%u_max = max(self%u_max,u_max)
  end do
  ! ----------------------------------------------------------------------------
  ! If there is a limit, how much RT can advance without MHD update
  ! ----------------------------------------------------------------------------
  if (rt_freq>0) then
    self%u_max = max(self%u_max,self%patch%u_max/rt_freq)
  end if
  if (verbose > 2) then
    call self%stats (tt  , 'tt')
    call self%stats (pgas, 'pgas')
    call self%stats (rk(:,:,:,1), 'rkap')
    call self%stats (fmax, 'fmax')
  end if
  ! ----------------------------------------------------------------------------
  deallocate (fmax, pgas, lnd, eth, e1, e2, lnd1, lnd2)
  if (.not.self%warmup_done) then
    if (self%time > self%warmup_end_time) then
      self%warmup_done=.true.
      self%dt_fixed = -1.0
    end if
  end if
  call trace%end (itimer)
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
  call self%boundaries%radiation (mem=self%q, src=self%src, rk=self%rk, &
                                     tt=self%tt, dir=self%axis, mu=self%rt_hat)

  call rt_integral%solve (self%rk, self%src, self%q, self%mesh, self%axis, &
                          self%rt_hat)
  call trace%end (itimer)
END SUBROUTINE solve


!===============================================================================
!> Add the (1+n_omega) RT tasks to the task list
!===============================================================================
SUBROUTINE init_task_list(self, task_list)
  class(rt_t):: self
  class(task_t),pointer:: rt_task
  class(list_t), pointer:: task_list
  class(rt_t), pointer:: omega
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%add_tasks_to')
  !-----------------------------------------------------------------------------
  self%task_list => task_list
  if (self%on) then
    rt_task => self%main        ! point to self
    if (.not.associated(self%link)) &
      call self%task_list%prepare(rt_task)
    if (.not.self%in_task_list) &
      call self%task_list%append_link (self%link)
    self%in_task_list = .true.
    omega => self%omeg
    do while(associated(omega))
      omega%task_list => task_list
      rt_task => omega
      if (.not.associated(omega%link)) &
        call self%task_list%prepare(rt_task)
      if (.not.omega%in_task_list) &
        call self%task_list%append_link (omega%link)
      omega%in_task_list = .true.
      omega => omega%next
    end do
  end if
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Determine if a task is needed (is upstream).  
!===============================================================================
FUNCTION needs (self, omega)
  logical:: needs
  class(rt_t):: self, omega
  !.............................................................................
  real(8):: ds1(3), ds2(3)
  real:: sinth
  real, parameter:: eps=1e-4
  integer:: is1(3), is2(3), is3(3), is4(3), i, lbit(3), ubit(3)
  logical:: lower, upper, upstream
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and omega points inward, or the task has no RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  if (.not.self%on) then
    needs=.false.
    return
  end if
  !
  if (self%i_omega /= omega%i_omega) then
    needs=.false.
    return
  end if
  !
  call trace%begin ('rt_t%needs')
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.omega%on) then
    needs = .false.
    call trace%end ()
    return
  end if
  !
  needs=.true.
  ds1 = self%position-omega%position
  where (self%periodic)
    ds1 = modulo(ds1 + 0.5_8*self%box, self%box) - 0.5_8*self%box
  end where
  ds1 = ds1/self%box                                                            ! relative task offsets
  sinth = merge(0.0, sqrt(1.0-self%mu**2), abs(self%mu)==1.0)
  ds2 = [sinth*cos(self%rt_phi),sinth*sin(omega%rt_phi),self%mu]
  is1 = 0; is2=0
  do i=1,3
    if (ds1(i)         > +eps) is1(i)=1
    if (self%rt_hat(i) > +eps) is2(i)=1
    if (ds1(i)         < -eps) is1(i)=-1
    if (self%rt_hat(i) < -eps) is2(i)=-1
  end do
  if (all(is1*is2 <= 0)) needs = .false.
  ! ----------------------------------------------------------------------------
  ! If the nbor lies outside the boundaries it is not needed
  ! ----------------------------------------------------------------------------
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
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.omega%on) then
    needs = .false.
    call trace%end ()
    return
  end if
  call trace%end ()
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
SUBROUTINE check_boundaries(self)
  class(rt_t):: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3)
  integer :: lbit(3), ubit(3), i
  ! ............................................................................
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
    else if (llc(i) <= self%rt_llc(i).and. .not.self%periodic(i)) then
      call self%boundaries%set (lbit(i))
    end if
    !
    if ((llc(i) > self%rt_urc(i))) then
      call self%patch%set(bits%no_rt)
      self%on = .false.
    else if (urc(i) >= self%rt_urc(i).and. .not.self%periodic(i)) then
      call self%boundaries%set (ubit(i))
    end if
  end do
END SUBROUTINE check_boundaries

!===============================================================================
!> Override the init_bndries, residing in patch_t
!===============================================================================
SUBROUTINE init_bdries(self)
  class(rt_t):: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3)
  integer :: lbit(3), ubit(3), i
  ! ............................................................................
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
    else if (llc(i) <= self%rt_llc(i).and. .not.self%periodic(i)) then
      call self%boundaries%set (lbit(i))
      self%mesh(i)%lower_boundary = .true.
    end if
    !
    if ((llc(i) > self%rt_urc(i))) then
      call self%patch%set(bits%no_rt)
      self%on = .false.
    else if (urc(i) >= self%rt_urc(i).and. .not.self%periodic(i)) then
      call self%boundaries%set (ubit(i))
      self%mesh(i)%upper_boundary = .true.
    end if
  end do
  call self%boundaries%init ()
END SUBROUTINE init_bdries

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
  class(nbor_t):: nbor
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(3x,a,i6,2i4,2x,3l1,2x,2l1,2x,3f9.3,i4,2f7.3)') &
    'rt_t%nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
    self%id, self%rank, self%level, nbor%needed, nbor%needs_me, nbor%download, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), &
    self%position, self%i_omega, self%mu, self%rt_phi
END SUBROUTINE nbor_info

!===============================================================================
!> This procedure is called from rt_setup, at a point where nbor lists have
!> already been made (based on needs() btw RT tasks), so the only thing to do
!> is to set the three flags appropriately.
!===============================================================================
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
  ! If this RT task has no RT on, it should have no nbor relations
  ! ----------------------------------------------------------------------------
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if
  ! ----------------------------------------------------------------------------
  ! First check the main relations, based on location and level of refinement
  ! ----------------------------------------------------------------------------
  call self%task_t%nbor_relations(another, needed, needs_me, download)
  ! ----------------------------------------------------------------------------
  ! Now, if another is RT, add it appropriately, if it is MHD, it has to overlap
  ! and have different kind of relations
  ! ----------------------------------------------------------------------------
  select type (another)
  class is (rt_t)
    rt_another => another
    if (.not.(rt_another%on)) then
      ! ------------------------------------------------------------------------
      ! if nbor RT is not on, we don't need any relations to it.
      ! ------------------------------------------------------------------------
      needed = .false.; needs_me = .false.; download = .false.
      call trace%end()
      return
    end if
    if ((rt_another%level<self%level) .and. (rt_another%i_omega==self%i_omega)) then
      if (all(rt_another%size>(abs(self%position-rt_another%position)))) then
        print *, 'found RT parent', self%i_omega, self%level, rt_another%level
        self%parent => rt_another
        self%link%parent => rt_another%link
        self%parentid = rt_another%id
        return
      end if
    end if
    !
    if ((rt_another%i_omega == 0).and.(self%i_omega==0)) then
      !-------------------------------------------------------------------------
      ! omega0 need all the relations with another omega0, thus defaut setting
      ! ------------------------------------------------------------------------
      if (dnload_mhd) download = .false.
      call trace%end()
      return
    else if ((rt_another%i_omega == 0).and.self%i_omega/=0) then
      if (self%main%id == rt_another%id) then
        ! ----------------------------------------------------------------------
        ! non-omega0 needs the master omega0
        ! ----------------------------------------------------------------------
        needs_me = .false.
        download = .false. 
      else 
        ! ----------------------------------------------------------------------
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        !-----------------------------------------------------------------------
        needed = .false.; needs_me = .false.; download = .false.
      end if
    else if (self%i_omega==0) then
      if (self%id==rt_another%main%id) then
        ! ----------------------------------------------------------------------
        ! non-omega0 needs the master omega0, which is self
        ! ----------------------------------------------------------------------
        !needed = .false.
        download = .false. 
      else 
        ! ----------------------------------------------------------------------
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        !-----------------------------------------------------------------------
        needed = .false.; needs_me = .false.; download = .false.
      end if
    !
    else if (self%needs(rt_another)) then
      !-------------------------------------------------------------------
      ! The omega task needs the upstream task,
      ! but the upstream task does not need the downstream omega task.
      !-------------------------------------------------------------------
      if ((self%rank==mpi%rank).and.(rt_another%rank/=mpi%rank)) needed=.true.
      needs_me=.false.
    ! --------------------------------------------------------------------------
    else if (rt_another%needs(self)) then
      ! Downstream RT task will need me, but it is not needed and I should not
      ! download data from it.
      ! ------------------------------------------------------------------------
      if ((self%rank==mpi%rank).and.(rt_another%rank/=mpi%rank)) needs_me=.true.
      needed = .false.
      download = .false.
    else
      needed = .false.; needs_me = .false.; download = .false.
    end if
    
    if ((self%i_omega >0).and.(rt_another%i_omega>0)) then
      if (self%i_omega /= rt_another%i_omega) then
      download = .false.; needed = .false.; needs_me=.false.
      end if
    end if
    !
  class is (gpatch_t)
    other => another
    if (all(self%position==other%position)) then
  ! ----------------------------------------------------------------------------
  ! Check if this is the "master" MHD patch. It might not be associated yet, i.e.
  ! after moving to another MPI rank or changing level of refinement. in that 
  ! case, associate it.
  ! ----------------------------------------------------------------------------
      if (self%i_omega==0) then
        !needs_me = .false.
        download = .false.
      else
        needs_me=.false.
        needed = .false.
        download = .false.
      end if
      if (.not.(associated(self%patch))) self%patch => gpatch%cast2gpatch(other)
    else
      needed = .false.; needs_me = .false.; download = .false.
    end if 
  class default
    needed = .false.; needs_me = .false.; download = .false.
  end select

  !
  call trace%end()
END SUBROUTINE nbor_relations

!===============================================================================
!> Auxiliary output
!===============================================================================
SUBROUTINE output (self)
  class(rt_t):: self
  !-----------------------------------------------------------------------------
END SUBROUTINE output

!===============================================================================
!===============================================================================
recursive SUBROUTINE dealloc(self)
  class(rt_t):: self
  class(rt_t), pointer:: omega, next
  ! ----------------------------------------------------------------------------
  integer:: i
  ! ............................................................................
  if (self%i_omega==0) then
    if (associated(self%w_bin)) deallocate (self%w_bin)
    omega => self%omeg
    do while(associated(omega))
      next => omega%next
      call omega%dealloc()
      omega => next
    end do
  else
  !
  end if
  if (associated(self%link)) then
    call self%link%remove_nbor_list2 (self%link%nbor)
    if (associated(self%task_list)) &
      call self%task_list%remove (self%link)
  end if
  call self%patch_t%dealloc()
END SUBROUTINE dealloc

END MODULE rt_mod
