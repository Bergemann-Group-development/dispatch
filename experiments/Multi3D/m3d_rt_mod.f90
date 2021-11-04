!===============================================================================
!> Extension of short charadcteristis RT module with Multi3D specifics
!===============================================================================
MODULE m3d_rt_mod
  USE io_mod
  USE trace_mod
  USE rt_mod
  USE rt_integral_mod
  USE multi3d_mod
  USE radau_mod
  USE boundaries_mod
  USE math_mod
  USE bits_mod
  USE gpatch_mod
  USE scs_mod, only: dis_package_t
  implicit none
  private
  type, public, extends(rt_t):: m3d_rt_t
    !---------------------------------------------------------------------------
    ! Before adding things here, check of the variable already exists in rt_t,
    ! and consider if it should be shared, or Multi3D specific (new name)
    !---------------------------------------------------------------------------
    real:: task_variable
    !---------------------------------------------------------------------------
    ! Add arrays here to store population numbers are whatever else
    !---------------------------------------------------------------------------
    !real, pointer, contiguous:: Jint(:,:,:,:) => null()
    !real, pointer, contiguous:: Npop(:,:,:,:) => null()
    integer:: fbin
    type(m3d_rt_t), pointer:: sheep  => null()        ! different RT directions
    type(m3d_rt_t), pointer:: shepherd => null()      ! point to main RT task
    type(m3d_rt_t), pointer:: next_sh  => null()      ! point to next omega task
    type(m3d_rt_t), pointer:: prev_sh  => null()      ! point to prev omega task
    type(m3d_rt_t), pointer:: tail_sh  => null()      ! point to prev omega task
    type(dis_package_t) :: dispackage
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_bdries
    procedure:: update
    procedure:: output
  end type
  integer, save:: verbose=0                           ! verbosity
  real,    save:: rt_freq=0                           ! how often RT should be updated WRT MHD -- 0=unlimited
  logical,save:: dnload_mhd=.true.
CONTAINS

!===============================================================================
!> Pre-initialize: read parameters, if any.  Probably we can rely on the
!> existing parameters inside Multi3D and DISPATCH.
!===============================================================================
recursive SUBROUTINE pre_init (self)
  class(m3d_rt_t), target:: self
  !-----------------------------------------------------------------------------
  class(m3d_rt_t), pointer:: omega
  character(128), save:: folder='input_multi3d'
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n, fbin
  integer, save:: n_mu=3, n_phi=4, n_bin=4, nt=1, n_warmup=3, ng(3)=[2,2,1], loc(1)
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
  namelist /multi3d_params/ folder, on, verbose, nt, ng, rt_llc, rt_urc, rt_timestep, &
  n_mu, n_phi, n_bin, rt_grace, vertical_propagation, n_warmup, courant, cdtd, rt_freq, dnload_mhd
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%pre_init')

  if (first_time) then
    call io%lock%set ('multi3d_t')
    rewind (io%input)
    read (io%input, multi3d_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, multi3d_params)
    end if
    call io%namelist_warning ('multi3d_params', iostat)
    
    ! ----------------------------------------------------------------------------
    ! Starting up the Multi3D Module
    ! ----------------------------------------------------------------------------
    multi3d%dir = folder
    call multi3d%pre_init(self%dispackage)

    first_time = .false.
    call io%lock%unset ('multi3d_t')
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

  ! ----------------------------------------------------------------------------
  ! Fill namelist variables
  ! ----------------------------------------------------------------------------
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
    self%type    = 'm3d_rt_t'
    if (self%i_omega==0) then
      self%kind    = 'm3d_rt_shepherd'
    else
      self%kind    = 'm3d_rt_sheep'
    end if
    !
    if (self%i_omega==0) then
      self%shepherd => self         ! convenience
      call rt_integral%init
      call self%patch%aux%init
      !-------------------------------------------------------------------------
      ! Number of ray directions
      !-------------------------------------------------------------------------
      n_omega = merge (1+self%n_phi*(self%n_mu-1), 0, self%n_mu > 0)
      self%n_omega = 2*n_omega*n_bin
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
              do fbin=1, n_bin
                i_omega = i_omega+1
                allocate(omega)
                omega%shepherd => self
                omega%master_id = self%id
                if (associated(self%omeg)) then
                  omega%prev_sh => self%tail_sh
                  omega%prev_sh%next_sh => omega
                  self%tail_sh => omega
                else
                  self%sheep => omega
                  self%tail_sh => omega
                end if
                omega%mu      = dir*mu(i_mu)
                omega%w       = w_mu(i_mu)*math%pi4/n/2.
                omega%patch   => self%patch
                omega%i_omega = i_omega
                omega%fbin    = fbin
                omega%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + &
                  (dir+1)/2.0, 2.0) - 1.0)
                if (i_mu == self%n_mu) omega%rt_phi = 0.0
                call omega%pre_init
              end do
            end do
          end do
        end do
      end if
    else
      ! ------------------------------------------------------------------------
      ! Determine closest axis direction for omega/=0
      ! ------------------------------------------------------------------------
      self%rt_hat(1) = sin(acos(self%mu))*cos(self%rt_phi)   ! sin(theta)*cos(phi), dx
      self%rt_hat(2) = sin(acos(self%mu))*sin(self%rt_phi)   ! sin(theta)*sin(phi), dy
      self%rt_hat(3) = self%mu                               ! cos(theta), dz
      loc = maxloc(abs(self%rt_hat))
      self%axis = loc(1)                                          ! direction of max propagation

      self%dispackage%geometry%nx = self%n(1)
      self%dispackage%geometry%ny = self%n(2)
      self%dispackage%geometry%nz = self%n(3)
      self%dispackage%mux   = dble(self%rt_hat(1))
      self%dispackage%muy   = dble(self%rt_hat(2))
      self%dispackage%muz   = dble(self%rt_hat(3))
      self%dispackage%fbin  = self%fbin
      self%dispackage%nfint = self%n_bin
      self%dispackage%wm    = dble(self%w)
      self%dispackage%ds    = dble(self%ds)

      call multi3d%init(self%dispackage)
    end if
  end if
  
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the rt_t data types.  To setup for n_warmup initial solutions, we
!> set the time rt_t times = -n_warmup * rt_t%timestep, making sure to call
!> rt_t%rt_eos for the initial MHD state. This subroutine is recursive, as
!> omega0 initializes all the omegaX. 
!===============================================================================
SUBROUTINE init (self)
  class (m3d_rt_t):: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
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
    self%nv = self%n_bin*2+2 ! nbin*[rk,src]+tt+qr
    self%irk(1)  = 1
    self%irk(2)  = self%irk(1) +self%n_bin-1
    self%isrc(1) = self%irk(2) +1
    self%isrc(2) = self%isrc(1)+self%n_bin-1
    self%itt     = self%isrc(2)+1
    self%iqr     = self%itt+1
  else
    self%nv = 1
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
    !allocate (self%w_bin (self%n_bin))
    !self%w_bin = 1.0
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
!> Update the population numbers and other dynamic task data.  Consider spliting
!> this into three steps; put data to Multi3D in pre_update, update, and get
!> data back from Mulit3D in post_update.
!===============================================================================
SUBROUTINE update (self)
  class(m3d_rt_t) :: self
  real, dimension(:,:,:), pointer :: temp, rho, pe, vx, vy, vz
  integer :: it
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%update')

  it = 1

  temp => self%patch%mem(:,:,:,self%patch%idx%e, it, 1)
  rho  => self%patch%mem(:,:,:,self%patch%idx%d, it, 1)
  pe   => self%patch%mem(:,:,:,6               , it, 1)

  vx   => self%patch%vgas(:,:,:,1,it)
  vy   => self%patch%vgas(:,:,:,2,it)
  vz   => self%patch%vgas(:,:,:,3,it)

  call multi3d%pre_update(self%dispackage, temp, rho, pe, vx, vy, vz)

  !call self%rt_t%update
  !call multi3d%update
  !call multi3d%post_update ! (... argument list ...)
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Do Multi3D specific output
!===============================================================================
SUBROUTINE output (self)
  class(m3d_rt_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%output')
  !call multi3d%output
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Override the init_bndries, residing in patch_t
!===============================================================================
SUBROUTINE init_bdries(self)
  class(m3d_rt_t):: self
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

END MODULE m3d_rt_mod
