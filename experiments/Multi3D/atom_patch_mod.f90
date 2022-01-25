module atom_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE solver_mod
    USE link_mod
    USE task_mod
    USE omp_timer_mod
    !USE index_mod
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    !USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE opac_patch_mod
    USE rad_patch_mod
    USE rate_patch_mod
    USE atom_var, only: atom_type, read_atom
    USE trans_mod, only: line_type
    USE spectrum_mod, only:spectrum_type
    USE composition_mod, only: composition_type
    USE setup_mod
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: atom_patch_t
      class(gpatch_t), pointer :: patch            ! pointer to MHD patch
      type(atom_type) :: atom
      type(spectrum_type) :: sp
      class(composition_type), pointer :: composition
      class(rate_patch_t), pointer, dimension(:,:) :: rate
      class(opacity_t), pointer, dimension(:,:) :: opac
      real(8), dimension(:,:,:), pointer :: temp => null()
      real(8), dimension(:,:,:), pointer :: rho  => null()
      real(8), dimension(:,:,:), pointer :: ne   => null()
      real(8), dimension(:,:,:,:), pointer :: nh => null()
      real(8) :: t_offset=0.1d0
      integer :: iter, maxiter
    contains
      procedure :: pre_init => atom_task_pre_init
      procedure :: init => atom_task_init
      procedure :: update => atom_task_update
      procedure :: nbor_relations => atom_task_nbor_relations
      procedure pre_init_atom
      procedure init_atom
      procedure append_to_task_list
      procedure init_level
      procedure task_info
      procedure normalize_phi
      procedure sum_rates
      procedure update_populations
      procedure freeze_atom
    end type atom_patch_t

    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_task_pre_init (self)
    class(atom_patch_t), target:: self
    ! ----------------------------------------------------------------------------
    call self%gpatch_t%pre_init()
    call self%pre_init_atom()
    call self%task_t%init_unique (same=.true.) ! obtain an unique ID
    
    ! ----------------------------------------------------------------------------
    ! relevant synchronization items for cartesian setup
    ! ----------------------------------------------------------------------------
    self%n          = self%patch%n
    self%ng         = 0 !self%patch%ng
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
  END SUBROUTINE atom_task_pre_init
  
  !===============================================================================
  SUBROUTINE atom_task_init (self)
    class(atom_patch_t) :: self
    !real(8) :: grace
    integer :: i
    ! ----------------------------------------------------------------------------
    call self%gpatch_t%init()
    self%kind = trim(self%atom%id) // '_patch'
    self%nv = 0
    !-----------------------------------------------------------------------------
    ! Synchronize with MHD patch
    ! ----------------------------------------------------------------------------
    self%dt_fixed = self%patch%dtime
    self%dtime    = self%patch%dtime
    self%dt       = self%patch%dtime

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)
    call self%set (bits%internal)
  
    !grace = 0.1d0
  
    self%time = self%patch%time + self%t_offset*self%dtime
  END SUBROUTINE atom_task_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(atom_patch_t) :: self
    type(atom_type), save :: atom
    type(spectrum_type), save :: sp
    logical, save :: first_time=.true.
    character(len=128), save :: atom_file='./input_multi3d/atoms/atom.ba06'
    real(8), save :: qnorm=10.0
    integer, save :: maxiter=30
    integer :: iostat
    namelist /atom_params/ atom_file, qnorm, maxiter
    !-----------------------------------------------------------------------------
    if (first_time) then
      write (*,*) ' - Reading atom'
      rewind (io%input)
      read (io%input, atom_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, atom_params)
      end if
      call io%namelist_warning ('atom_params', iostat)

      call read_atom(atom, atom_file, qnorm)
      call sp%init(atom%line, atom%cont, atom%qnorm)
      atom%sp = sp

      first_time = .false.
    endif

    self%atom = atom
    self%sp   = sp
    self%maxiter = maxiter + 1 ! we start from zeroth iteration

    allocate(self%atom%C(atom%nlevel, atom%nlevel))
    allocate(self%atom%W(atom%nlevel, atom%nlevel))
    allocate(self%atom%e(atom%nlevel))
    !allocate(self%atom%pos(3,atom%nlevel))
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE init_atom (self)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    !-----------------------------------------------------------------------------
    atom => self%atom

    call atom%compute_line_damping(self%ne, self%temp, self%nh)
    call self%composition%get_totn(trim(atom%id), atom%abnd, atom%totn)
    call atom%compute_lte_populations(self%temp, self%ne)
  END SUBROUTINE init_atom

  !===============================================================================
  SUBROUTINE atom_task_update (self)
    class(atom_patch_t) :: self
    real(8), dimension(:,:,:,:), allocatable :: gamma_ij, gamma_ji
    character(len=64) :: text
    real(8), save :: time=0
    ! ----------------------------------------------------------------------------
    call trace_begin('atom_patch_t%update')

    if (self%iter == 0) then
      time = wallclock()
      call self%init_atom()
    elseif (self%iter == 1) then
      call self%normalize_phi()
    endif

    if (self%iter > 0) then
      allocate(gamma_ij(self%n(1), self%n(2), self%n(3), self%atom%nrad))
      allocate(gamma_ji(self%n(1), self%n(2), self%n(3), self%atom%nrad))

      call self%sum_rates(gamma_ij, gamma_ji)
      call self%update_populations(gamma_ij, gamma_ji)
    endif

    if (self%iter < self%maxiter) then
      write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter-1, '.'
      call io%header (text)

    elseif(self%iter == self%maxiter) then
      write(text,'(A,I2,A,A)') 'End of iteration ', self%iter-1, '.'
      call io%header (text)
      call self%freeze_atom()
      print *, 'Dispatch finished after ', wallclock() - time, ' seconds'
      stop
    endif

    write(*,'(A10, 1(I4), 3(ES18.7))'), trim(self%kind), self%patch%id, self%atom%n(1,1,1,3)
    write(*,'(A10, 1(I4), 3(ES18.7))'), trim(self%kind), self%patch%id, self%atom%n(1,1,self%n(3),3)
  
    self%iter = self%iter + 1
  
    call trace%end()
  END SUBROUTINE atom_task_update

  !===============================================================================
  SUBROUTINE update_populations (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    integer, dimension(3, self%atom%nlevel) :: pos3
    integer, dimension(self%atom%nlevel) :: maxerr
    real(8), pointer :: temp, ne, totn
    real(8), dimension(:), pointer :: nh, Gij, Gji
    real(4), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    integer :: i, j, k, ii
    integer :: cswitch=1
    ! ----------------------------------------------------------------------------
    atom => self%atom

    maxerr   = 0.0d0

    do k=1, self%n(3)
      do j=1, self%n(2)
        do i=1, self%n(1)
          temp  => self%temp (i,j,k)
          ne    => self%ne   (i,j,k)
          totn  => atom%totn (i,j,k)
          nh    => self%nh   (i,j,k,:)
          n     => atom%n    (i,j,k,:)
          nstar => atom%nstar(i,j,k,:)
          Gij   => gamma_ij  (i,j,k,:)
          Gji   => gamma_ji  (i,j,k,:)

          call atom%solve_stat_equilibrium(temp, ne, totn, n, nstar, nh, Gij, Gji, cswitch)

          ! keep track of maximum population change
          do ii=1,atom%nlevel
            if (abs(atom%e(ii)) > maxerr(ii)) then
                maxerr(ii) = atom%e(ii)
                pos3(1,ii) = i  
                pos3(2,ii) = j 
                pos3(3,ii) = k
            end if
          end do
        end do
      end do
    end do
  END SUBROUTINE update_populations

  !===============================================================================
  SUBROUTINE normalize_phi (self)
    class(atom_patch_t), target :: self
    class(line_type), pointer :: line
    class(opacity_t), pointer :: opac
    integer :: i, i_omega, i_nu
    ! ----------------------------------------------------------------------------
    do i=1, self%atom%nline
      line => self%atom%line(i)
      allocate(line%wphi(self%n(1), self%n(2), self%n(3)))
      line%wphi = 0

      do i_nu=1, size(self%opac(1:,1:), 2)
        do i_omega=1, size(self%opac(1:,1:), 1)
          opac => self%opac(i_omega, i_nu)
          line%wphi = line%wphi + opac%wmu * opac%wphi(:,:,:,i)
        end do
      end do

      line%wphi = math%pi4 / line%wphi
    end do
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE sum_rates (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(rate_patch_t), pointer :: rate
    integer :: i, i_omega, i_nu
    real(8), dimension(:,:,:,:) :: gamma_ij, gamma_ji
    ! ----------------------------------------------------------------------------
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0

    do i_nu=1, size(self%rate, 2)
      do i_omega=1, size(self%rate, 1)
        rate => self%rate(i_omega, i_nu)
        gamma_ij = gamma_ij + rate%wmu * rate%mem(:,:,:,:,rate%it,rate%i_Gij)
        gamma_ji = gamma_ji + rate%wmu * rate%mem(:,:,:,:,rate%it,rate%i_Gji)
      end do
    end do

    do i=1, self%atom%nline
      gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * self%atom%line(i)%wphi
      gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * self%atom%line(i)%wphi
    end do
  END SUBROUTINE
  
  !===============================================================================
  SUBROUTINE atom_task_nbor_relations(self, another, needed, needs_me, download)
    class(atom_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    
    call trace%begin ('atom_patch_t%nbor_relations')

    needed = .false.; needs_me = .false.; download = .false.
  
    if (another%id == self%patch%id) then
      needed = .true.; needs_me = .false.; download = .false.
    else
      select type (another)
      class is (opacity_t)
        if (associated(another%atom)) then
          if (another%atom%id == self%atom%id) then
            needed = .false.; needs_me = .true.; download = .false.
          endif
        endif
      class is (rate_patch_t)
        if (associated(another%atom)) then
          if (another%atom%id == self%atom%id) then
            needed = .true.; needs_me = .false.; download = .false.
          endif
        endif
      end select
    endif
  
    call trace%end()
  END SUBROUTINE atom_task_nbor_relations

  !===============================================================================
  SUBROUTINE append_to_task_list(self, task_list)
    class (atom_patch_t), target :: self
    class (list_t) :: task_list
    class(link_t),  pointer :: link
    class (task_t), pointer :: task
    !-----------------------------------------------------------------------------
    allocate (halo_m3d_t:: link)
    link%task => self
    self%link => link
    call task_list%append (link)
  END SUBROUTINE append_to_task_list

  !===============================================================================
  !> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
  !===============================================================================
  SUBROUTINE task_info (self)
    class(atom_patch_t):: self
    !-----------------------------------------------------------------------------
    write(io_unit%mpi,'(1x,a,2x,a10,i6,2i4,7x,2l1,2x,3f9.3,i4,i4,i4)') &
      'atom_patch_t%task_info: kind, id, rank, level, BV, pos =', &
      self%kind, self%id, self%rank, self%level, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info

  !===============================================================================
  SUBROUTINE init_level (self)
    class(atom_patch_t):: self
    !-----------------------------------------------------------------------------
    if (any(self%n > 1)) then
      ! The level does not depend on the patch size in non-regular grids
      self%level = 1
    end if
  END SUBROUTINE init_level

  !===============================================================================
  SUBROUTINE freeze_atom (self)
    class(atom_patch_t):: self
    integer :: i_nu, i_omega
    !-----------------------------------------------------------------------------
    call self%set(bits%frozen)
    self%time = 0 !self%patch%time + self%t_offset

    do i_nu=1, size(self%rate, 2)
      call self%opac(0, i_nu)%set(bits%frozen)
      self%opac(0, i_nu)%time = 0 !self%patch%time + self%opac(0, i_nu)%t_offset

      do i_omega=1, size(self%rate, 1)
        call self%opac(i_omega, i_nu)%set(bits%frozen)
        call self%rate(i_omega, i_nu)%set(bits%frozen)
        self%opac(i_omega, i_nu)%time = 0 !self%patch%time + self%opac(i_omega, i_nu)%t_offset
        self%rate(i_omega, i_nu)%time = 0 !self%patch%time + self%rate(i_omega, i_nu)%t_offset
      end do
    end do
  END SUBROUTINE freeze_atom
  
  end module