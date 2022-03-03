module atom_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
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
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, ilambd
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      type(atom_type), pointer :: atom
      class(composition_type), pointer :: composition
      class(rate_patch_t), pointer, dimension(:,:) :: rate
      class(opacity_t), pointer, dimension(:,:) :: opac
      real(8) :: maxerr=huge(1.0)
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: init => atom_patch_init
      procedure :: update => atom_patch_update
      procedure :: post_update => atom_patch_post_update
      procedure :: finish => atom_patch_finish
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure pre_init_atom
      procedure add_to_task_list
      procedure init_atom
      procedure normalize_phi
      procedure sum_rates
      procedure update_populations
    end type atom_patch_t

    real(8), public:: global_maxerr=0
    real(8), public:: convlim=1d-3

    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_patch_pre_init (self)
    class(atom_patch_t), target:: self
    integer :: m, n
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_pre_init')
    call self%m3d_patch_t%pre_init()
    call self%set_nv_and_nw(0, 0)
    call self%pre_init_atom()

    allocate(self%opac(0:self%n_mu, self%n_nu))
    allocate(self%rate(self%n_omega, self%n_nu))

    do n=1, self%n_nu
      self%opac(0,n)%kind = 'cont_opac'
      self%opac(0,n)%atom => self%atom
      call self%opac(0,n)%pre_init_and_link(self)
    end do

    do m=1, self%n_mu
      do n=1, self%n_nu
        self%opac(m,n)%kind = 'line_opac'
        self%opac(m,n)%atom => self%atom
        call self%opac(m,n)%pre_init_and_link(self, i_mu=m, i_nu=n)
        
        if (m > self%n_omega) cycle
        self%rate(m,n)%kind = 'rad_rate'
        self%rate(m,n)%atom => self%atom
        self%rate(m,n)%cont_opac => self%opac(0,n)
        self%rate(m,n)%line_opac => self%opac(m,n)
        call self%rate(m,n)%pre_init_and_link(self, i_mu=m, i_nu=n)
      end do
    end do

    call trace%end()
  END SUBROUTINE atom_patch_pre_init
  
  !===============================================================================
  SUBROUTINE atom_patch_init (self)
    class(atom_patch_t) :: self
    integer :: m, n
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_pre_init')
    call self%m3d_patch_t%init()

    do n=1, self%n_nu
      call self%opac(0,n)%init()
    end do

    do m=1, self%n_mu
      do n=1, self%n_nu
        call self%opac(m,n)%init()
        if (m > self%n_omega) cycle
        call self%rate(m,n)%init()
      end do
    end do
    call trace%end()
  END SUBROUTINE atom_patch_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(atom_patch_t) :: self
    type(atom_type), save :: atom
    type(spectrum_type), save, target :: sp
    logical, save :: first_time=.true.
    character(len=128), save :: atom_file='./input_multi3d/atoms/atom.ba06'
    real(8), save :: qnorm=10.0
    integer :: iostat
    namelist /atom_params/ atom_file, qnorm, convlim
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_atom')
    if (first_time) then
      write (*,*) ' - Reading atom'
      rewind (io%input)
      read (io%input, atom_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, atom_params)
      end if
      call io%namelist_warning ('atom_params', iostat)
      write (io%namelist, atom_params)

      call read_atom(atom, atom_file, qnorm)
      call sp%init(atom%line, atom%cont, atom%qnorm)
      atom%sp = sp

      first_time = .false.
    endif

    allocate(self%atom)
    self%atom = atom
    self%sp => sp

    allocate(self%atom%C(atom%nlevel, atom%nlevel))
    allocate(self%atom%W(atom%nlevel, atom%nlevel))
    allocate(self%atom%e(atom%nlevel))
    call trace%end()
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE init_atom (self)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    real(4), dimension(:,:,:), pointer :: ne, temp
    real(4), dimension(:,:,:,:), pointer :: nh
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_atom')
    atom => self%atom

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)

    call atom%compute_line_damping(ne, temp, nh)
    call self%composition%get_totn(trim(atom%id), atom%abnd, atom%totn)
    call atom%compute_lte_populations(temp, ne)
    call trace%end()
  END SUBROUTINE init_atom

  !===============================================================================
  SUBROUTINE atom_patch_update (self)
    class(atom_patch_t) :: self
    real(8), dimension(:,:,:,:), allocatable :: gamma_ij, gamma_ji
    character(len=64) :: text
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    if (self%finished) return
    call trace%begin('atom_patch_t%atom_patch_update', itimer=itimer)

    if (self%iter == 0) then
      write(text,'(A,I2,A,I2,A)') 'LTE Iteration .'
      call io%header (text)

      call self%init_atom()
      self%dtime = self%dtime + self%dtime*ilambd
      self%rate%time = self%rate%time + self%rate%dtime*ilambd
      self%opac(:,:)%dtime = self%dtime
    else
      write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter, '.'
      call io%header (text)

      if (self%iter == 1) then
        self%dtime = self%dt_fixed
        self%opac(:,:)%dtime = self%opac%dt_fixed
        call self%normalize_phi()
      endif

      allocate(gamma_ij(self%n(1), self%n(2), self%n(3), self%atom%nrad))
      allocate(gamma_ji(self%n(1), self%n(2), self%n(3), self%atom%nrad))

      call self%sum_rates(gamma_ij, gamma_ji)
      call self%update_populations(gamma_ij, gamma_ji)
    endif
  
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE atom_patch_update

  !=============================================================================
  SUBROUTINE atom_patch_post_update (self)
    class(atom_patch_t) :: self
    ! --------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_post_update')
    if (self%verbose > 1) then
      call io%lock%set
       write(*,'(F5.1, L2, A12, 1(I4), 3(ES15.7))') &
         self%time, self%finished, trim(self%kind), self%patch%id, &
         self%atom%n(1,1,1,3), self%atom%n(1,1,self%n(3),3)
      write(*,'(A19, ES19.7, A10, 3ES15.7)') 'Maxerr =', self%maxerr, &
        ', global =', global_maxerr
      write(*,'(1X, A80)') &
      '................................................................................'
      call io%lock%unset
    end if

    if (self%iter == self%maxiter) then
      call self%finish()
    endif

    if ((global_maxerr > 0) .and. global_maxerr < convlim) then
      print *, 'Converged!'
      call self%finish()
    endif

    self%iter = self%iter + 1
    call trace%end()
  END SUBROUTINE atom_patch_post_update

  !===============================================================================
  SUBROUTINE atom_patch_finish (self, rotated)
    class(atom_patch_t) :: self
    logical, optional :: rotated
    integer :: i, j
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_finish')
    self%finished = .true.
    self%time = io%end_time + self%t_offset*self%dtime
    !call self%set(bits%frozen)

    do j=1, self%n_nu
      !call self%opac(0,j)%finish()

      do i=1, self%n_omega
        call self%opac(i,j)%finish()
        call self%rate(i,j)%finish()
      end do
    end do
    call trace%end()
  END SUBROUTINE atom_patch_finish

  !===============================================================================
  SUBROUTINE update_populations (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    integer, dimension(3, self%atom%nlevel) :: pos3
    real, dimension(self%atom%nlevel) :: maxerr
    real(4), pointer :: temp, ne, totn
    real(4), dimension(:), pointer :: nh
    real(8), dimension(:), pointer :: Gij, Gji
    real(4), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    integer :: i, j, k, ii
    integer :: cswitch=1
    ! ----------------------------------------------------------------------------
    atom => self%atom

    maxerr = 0.0d0

    do k=1, self%n(3)
      do j=1, self%n(2)
        do i=1, self%n(1)
          ne    => self%m3d_parent%mem(i,j,k,self%i_ne,            self%it,1)
          temp  => self%m3d_parent%mem(i,j,k,self%i_temp,          self%it,1)
          nh    => self%m3d_parent%mem(i,j,k,self%i_nh1:self%i_nh6,self%it,1)
          totn  => atom%totn  (i,j,k)
          n     => atom%n     (i,j,k,:)
          nstar => atom%nstar (i,j,k,:)
          Gij   => gamma_ij   (i,j,k,:)
          Gji   => gamma_ji   (i,j,k,:)

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

    if (global_maxerr == self%maxerr) then
      self%maxerr = maxval(abs(maxerr))
      global_maxerr = self%maxerr
    else
      self%maxerr = maxval(abs(maxerr))
      global_maxerr = max(global_maxerr, self%maxerr)
    endif
  END SUBROUTINE update_populations

  !===============================================================================
  SUBROUTINE normalize_phi (self)
    class(atom_patch_t), target :: self
    class(line_type), pointer :: line
    class(opacity_t), pointer :: opac
    integer :: i, i_mu, i_nu
    ! ----------------------------------------------------------------------------
    do i=1, self%atom%nline
      line => self%atom%line(i)
      allocate(line%wphi(self%n(1), self%n(2), self%n(3)))
      line%wphi = 0

      do i_nu=1, size(self%opac(1:self%n_omega,1:), 2)
        do i_mu=1, size(self%opac(1:self%n_omega,1:), 1)
          opac => self%opac(i_mu, i_nu)
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
    integer :: i, i_mu, i_nu, m(4)
    real(8), dimension(:,:,:,:), allocatable:: gamma_ij, gamma_ji
    ! ----------------------------------------------------------------------------
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0

    do i_nu=1, size(self%rate, 2)
      do i_mu=1, size(self%rate, 1)
        rate => self%rate(i_mu, i_nu)
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
  SUBROUTINE atom_patch_nbor_relations(self, another, needed, needs_me, download)
    class(atom_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    
    call trace%begin ('atom_patch_t%nbor_relations')

    needed = .false.; needs_me = .false.; download = .false.
  
    if (another%id == self%m3d_parent%id) then
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
  END SUBROUTINE atom_patch_nbor_relations

  !===============================================================================
  SUBROUTINE add_to_task_list(self, task_list)
    class (atom_patch_t) :: self
    class(list_t), pointer :: task_list
    class(link_t),  pointer :: link
    class (task_t), pointer :: task
    integer :: m,n
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%add_to_task_list')
    call self%append_to_task_list(task_list)

    do n=1, self%n_nu
      call self%opac(0,n)%append_to_task_list(task_list)
    end do

    do m=1, self%n_mu
      do n=1, self%n_nu
        call self%opac(m,n)%append_to_task_list(task_list)
        if (m > self%n_omega) cycle
        call self%rate(m,n)%append_to_task_list(task_list)
      end do
    end do
    call trace%end()
  END SUBROUTINE add_to_task_list

  end module
