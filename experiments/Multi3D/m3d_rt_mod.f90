!===============================================================================
!> Extension of short charadcteristis RT module with Multi3D specifics
!===============================================================================
MODULE m3d_rt_mod
  USE io_mod
  use omp_lib
  USE trace_mod
  USE rt_mod
  USE rt_integral_mod
  USE multi3d_mod
  USE radau_mod
  use mpi_mod
  USE boundaries_mod
  USE math_mod
  USE bits_mod
  USE gpatch_mod
  USE scs_mod, only: dis_package_t
  use task_mod
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
    integer:: fbin                                     ! Frequency bin
    integer,dimension(:,:),allocatable :: nu_regions   ! Frequency index limits
    integer :: nu_s,nu_e
    type(dis_package_t)  :: dispackage    ! Dispatch <-> Multi communicator
    integer :: i_n,i_n_lte,i_ne           ! Index range of pops and ne in mem arr
    integer :: i_chi,i_snu,i_inu,n_nu     ! Index range of chi, snu and inu 
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_bdries
    procedure:: update
    procedure:: output
    procedure:: get_nu_regions
    procedure:: reduce_omegaX
    procedure:: link_dispackage
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
    call multi3d%pre_init()

    first_time = .false.
    call io%lock%unset ('multi3d_t')
  end if

  ! Multi3D pre init is ran only once per OMP thread
  ! This means that another dispackages of the same thread can not be filled 
  ! within the same routine. Fill them here with thread global variables
  self%dispackage%ggpgh         = multi3d%ggpgh
  self%dispackage%ref_atom     => multi3d%ref_atom
  self%dispackage%ref_spectrum => multi3d%ref_spectrum

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
      self%kind    = 'm3d_rt_sc0'
    else
      self%kind    = 'm3d_rt_sc'
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
      self%n_omega = 2*n_omega*n_bin

      call self%get_nu_regions(n_bin)            ! Get frequency domains for all tasks
      
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

                omega%fbin    = fbin
                omega%nu_s    = self%nu_regions(fbin,1)  ! Starting frequency index in spectrum freq. grid
                omega%nu_e    = self%nu_regions(fbin,2)  ! Ending frequency index in spectrum freq. grid

                omega%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + &
                  (dir+1)/2.0, 2.0) - 1.0)
                if (i_mu == self%n_mu) omega%rt_phi = 0.0
                call omega%pre_init
              end do
            end do
          end do
        end do
      end if
      self%dispackage%geometry%nx = self%n(1)
      self%dispackage%geometry%ny = self%n(2)
      self%dispackage%geometry%nz = self%n(3)
      self%dispackage%fbin  = self%fbin
      self%dispackage%nu_s  = 1
      self%dispackage%nu_e  = self%dispackage%ref_spectrum%nnu
      self%dispackage%nfint = 0
      self%dispackage%wm    = dble(self%w)
      self%dispackage%ds    = dble(self%ds)
      call multi3d%init(self%dispackage,skip_ray_in=.true.)
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
      self%dispackage%nu_s  = self%nu_s
      self%dispackage%nu_e  = self%nu_e
      self%dispackage%nfint = self%n_bin
      self%dispackage%wm    = dble(self%w)
      self%dispackage%ds    = dble(self%ds)
      call multi3d%init(self%dispackage)
    end if
  end if
  call self%init ()
  call trace%end ()
END SUBROUTINE pre_init


! PE: obtain frequency splitted regions
subroutine get_nu_regions(self,n_splits) 
    class(m3d_rt_t),intent(inout):: self
    integer,intent(in) :: n_splits
    integer :: i,nu_bin_size

    ! Simple approach, can be adjusted in the future if needed (based on active lines etc.)
    if (allocated(self%nu_regions)) deallocate(self%nu_regions)
    allocate(self%nu_regions(n_splits,2))

    nu_bin_size =  floor(self%dispackage%ref_spectrum%nnu / real(n_splits))
    self%nu_regions(1,1)        = 1                         ! Starpoint
    self%nu_regions(1,2)        = nu_bin_size               ! Starpoint
    do i = 2, n_splits
        self%nu_regions(i,1) = self%nu_regions(i-1,2) + 1
        self%nu_regions(i,2) = self%nu_regions(i,1) + nu_bin_size - 1
    enddo
    self%nu_regions(n_splits,2) = self%dispackage%ref_spectrum%nnu  ! Endpoint
end subroutine get_nu_regions

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
  integer :: idiff
  ! ----------------------------------------------------------------------------
  if (self%initialized) then
    print *, self%id, self%i_omega
    call io%abort ('double rt init')
  end if 
  call trace%begin ('rt_t%init')
  self%initialized = .true.
  p => self%patch
  
  !======================================================================================
  ! Omega0: > Collects the transition rates from omegaX tasks
  !         > Stores population numbers + rates 
  !
  ! OmegaX: > Computes radiative transfer
  !         > Stores chi, snu and inu in its frequency range
  !
  ! The storage that needs to be communicated between the tasks is chi,snu and inu !!!
  ! The angle tasks do not need to know the population numbers in the ghost zones
  ! if the opacity is passed!
  !=====================================================================================
  if (self%i_omega==0) then
    self%nv      = 2   ! population numbers (+ heating)
    self%i_n     = 1
    self%i_n_lte = 2
    self%nw      = self%dispackage%ref_atom%nlevel
    !self%itt    = self%i_n(2) + 1
    !self%i_ne   = self%itt + 1
    !self%iqr    = self%i_ne + 1

    !self%nv = self%n_bin*2+2 ! nbin*[rk,src]+tt+qr
    !self%irk(1)  = 1
    !self%irk(2)  = self%irk(1) +self%n_bin-1
    !self%isrc(1) = self%irk(2) +1
    !self%isrc(2) = self%isrc(1)+self%n_bin-1
    !self%itt     = self%isrc(2)+1
    !self%iqr     = self%itt+1
  else
    idiff = self%dispackage%nu_e-self%dispackage%nu_s + 1
    self%n_nu = idiff
    self%nv = 3
    self%nw = idiff  ! chi, snu and inu (for every nu in this patch)
    self%i_chi = 1
    self%i_snu = 2
    self%i_inu = 3

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
    !self%mem(:,:,:,self%irk(1):self%irk(2),:,:)=tiny(1.)
    !self%mem(:,:,:,self%isrc(1):self%isrc(2),:,:)=tiny(1.)
    !self%mem(:,:,:,self%itt,:,:)=tiny(1.)
    
    ! Initialize the populations to -1 so that their absense can be caught in m3d
    self%mem(:,:,:,self%i_n,:,:)     = -1   
    self%mem(:,:,:,self%i_n_lte,:,:) = -1   

  end if
  write (io_unit%log,*) 'call self%patch_t%setup'

  !print *, 'IOMEGA, FBIN, SHAPE(MEM) ',self%i_omega, self%fbin, shape(self%mem)
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
  ! Set omegaX time to be marginally before omega0, so omega0 would be done after
  ! omegaX
  ! ----------------------------------------------------------------------------
  if (self%i_omega==0) then
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
    !self%pervolume(self%irk(1):self%irk(2))   = .true.
    !self%pervolume(self%isrc(1):self%isrc(2)) = .true.
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

  if (self%i_omega .eq. 0) call self%rt_t%init_task_list (self%patch%task_list)

  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update the population numbers and other dynamic task data.  Consider spliting
!> this into three steps; put data to Multi3D in pre_update, update, and get
!> data back from Mulit3D in post_update.
!===============================================================================
SUBROUTINE update (self)
  class(m3d_rt_t) :: self
  real, dimension(:,:,:),contiguous, pointer :: temp, rho, pe, vx, vy, vz
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

  if (self%i_omega .gt. 0) then
    call self%link_dispackage(i_omega = self%i_omega, it = it)                                    ! Establish links to the dispackage
    call multi3d%pre_update (self%dispackage, temp, rho, pe, vx, vy, vz, skip_ray_in = .false.)
    call multi3d%update (self%dispackage, solver = 'rad')                                         ! call rad. transfer solver
                                                                                                  ! The solver will compute the radiative transfer for one angle.
                                                                                                  ! it will assume that the upwind tasks are already done and 
                                                                                                  ! that the mem object is downloaded! (So Chi,Inu and Snu are available in the ghost zones)
                                                                                                  ! The boundary contidion needs to be filled in the bottom ghost zone (chi, snu and inu need to be present!)
  else
    call self%link_dispackage(i_omega = self%i_omega, it = it)                                    ! Establish links to the dispackage
    call multi3d%pre_update (self%dispackage, temp, rho, pe, vx, vy, vz, skip_ray_in = .true.)
    call self%reduce_omegaX ()                                                                    ! Reduce the rates etc. obtained from angle tasks
    call multi3d%update (self%dispackage, solver = 'pop')                                         ! call stat. eq. solver
  endif

  !call multi3d%post_update ! (... argument list ...)


  call trace%end ()
END SUBROUTINE update

subroutine reduce_omegaX(self)
  class(m3d_rt_t),target :: self
  class(rt_t),pointer:: omeg
  ! Loop through all the angle tasks and save the rates etc.
  omeg => self%main%omeg
  do while (associated(omeg))
    if ((omeg%master_id.ne.self%master_id).or.(omeg%i_omega.eq.0)) cycle ! This is not a task for this master
    
    ! Add the neighbor quants to the quants of this task
    select type (o=>omeg)
    class is (m3d_rt_t)
      self%dispackage%atom_rad%gamma_ij   = self%dispackage%atom_rad%gamma_ij   + o%dispackage%atom_rad%gamma_ij
      self%dispackage%atom_rad%gamma_ji   = self%dispackage%atom_rad%gamma_ji   + o%dispackage%atom_rad%gamma_ji
      self%dispackage%atom_rad%gamma_ij_c = self%dispackage%atom_rad%gamma_ij_c + o%dispackage%atom_rad%gamma_ij_c
      self%dispackage%atom_rad%gamma_ji_c = self%dispackage%atom_rad%gamma_ji_c + o%dispackage%atom_rad%gamma_ji_c

      self%dispackage%atom_rad%Rij   = self%dispackage%atom_rad%Rij   + o%dispackage%atom_rad%Rij
      self%dispackage%atom_rad%Rji   = self%dispackage%atom_rad%Rji   + o%dispackage%atom_rad%Rji
      self%dispackage%atom_rad%Rij_c = self%dispackage%atom_rad%Rij_c + o%dispackage%atom_rad%Rij_c
      self%dispackage%atom_rad%Rji_c = self%dispackage%atom_rad%Rji_c + o%dispackage%atom_rad%Rji_c

      self%dispackage%atom_rad%qrad    = self%dispackage%atom_rad%qrad   + o%dispackage%atom_rad%qrad
      self%dispackage%atom_rad%qrad_c  = self%dispackage%atom_rad%qrad_c + o%dispackage%atom_rad%qrad_c 
    end select  
  enddo

end subroutine reduce_omegaX

subroutine link_dispackage(self,i_omega,it)
  class(m3d_rt_t),intent(inout) :: self
  integer, optional :: i_omega,it
  integer :: i_omega_l,it_l

  ! Defaults
  i_omega_l = 0
  it_l      = 1
  if (present(i_omega)) i_omega_l = i_omega
  if (present(it)) it_l = it

  ! Link
  if (i_omega_l.eq.0) then  
    self%dispackage%n     => self%mem(:,:,:,self%i_n,it_l,:)
    self%dispackage%nstar => self%mem(:,:,:,self%i_n_lte,it_l,:)
  else
    self%dispackage%chi_disp => self%mem(:,:,:,self%i_chi,it_l,:)
    self%dispackage%snu_disp => self%mem(:,:,:,self%i_snu,it_l,:)
    self%dispackage%inu_disp => self%mem(:,:,:,self%i_inu,it_l,:)
  endif

end subroutine link_dispackage

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

!===============================================================================
!> The nbor relations among the omega
!
! The omega0 task has a different role than in the default setup!
! It will be called AFTER the omegaX tasks have solved the radiative transfer,
! because it will collect the rates and sets + solves the stat. eq. equations.
!
!===============================================================================
subroutine nbor_relations(self, another, needed, needs_me, download)
  class(m3d_rt_t)       :: self            ! current M3d rt task
  class(task_t),pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  class(m3d_rt_t),pointer :: rt_another
  class(task_t),pointer   :: other
  

  call trace%begin ('rt_t%nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  ! First check the main relations, based on location and level of refinement
  call self%task_t%nbor_relations(another, needed, needs_me, download)

  ! Now, if another is RT, add it appropriately, if it is MHD, it has to overlap
  ! and have different kind of relations
  select type (another)
  class is (m3d_rt_t)

    rt_another => another
    if (.not.(rt_another%on)) then
      ! if nbor RT is not on, we don't need any relations to it.
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
    
    if ((rt_another%i_omega == 0).and.(self%i_omega==0)) then                   ! Both are omega0
      ! omega0 need all the relations with another omega0, thus defaut setting
      if (dnload_mhd) download = .false.
      call trace%end()
      return
    else if ((rt_another%i_omega == 0).and.self%i_omega/=0) then                ! Self: X, other: 0
      if (self%main%id == rt_another%id) then
        ! non-omega0 doesnt need the master omega0
        needed   = .false.
        needs_me = .true.
        download = .false. 
      else 
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        needed = .false.; needs_me = .false.; download = .false.
      end if
    else if (self%i_omega==0) then                                              ! Self: 0, other: X
      if (self%id==rt_another%main%id) then
        ! omega0 needs other omegax tasks
        needed   = .true.
        needs_me = .false.
        download = .false. 
      else 
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        needed = .false.; needs_me = .false.; download = .false.
      end if
  
    else if (self%needs(rt_another)) then                                       ! rt_another is upstream
      
      ! The omega task needs the upstream task,
      ! but the upstream task does not need the downstream omega task.
      if ((self%rank==mpi%rank).and.(rt_another%rank/=mpi%rank)) needed=.true.
      needs_me=.false.
    
    else if (rt_another%needs(self)) then
      ! Downstream RT task will need me, but it is not needed and I should not
      ! download data from it.
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
    
  class is (gpatch_t)
    other => another
    if (all(self%position==other%position)) then
      ! Check if this is the "master" MHD patch. It might not be associated yet, i.e.
      ! after moving to another MPI rank or changing level of refinement. in that 
      ! case, associate it.
      needed   = .true.
      needs_me = .false.
      download = .false.
      if (self%i_omega .eq. 0) needs_me = .true.   ! For returning the heating rate
      if (.not.(associated(self%patch))) self%patch => gpatch%cast2gpatch(other)
    else
      needed = .false.; needs_me = .false.; download = .false.
    end if 
  class default
    needed = .false.; needs_me = .false.; download = .false.
  end select

  call trace%end()
end subroutine nbor_relations


END MODULE m3d_rt_mod
