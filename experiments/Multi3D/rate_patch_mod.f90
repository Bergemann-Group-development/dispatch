module rate_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    !USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    !USE index_mod
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    USE os_mod
    !USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE m3d_patch_mod, only: conditioning
    USE ray_patch_mod
    USE opac_patch_mod
    USE rad_patch_mod
    USE trans_mod, only: transition_type
    USE atom_mod, only: atom_type
    USE lines_patch_mod, only: lines_patch_t
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(ray_patch_t) :: rate_patch_t
      class(atom_type), pointer :: atom => null()
      class(radiation_t), pointer :: rad => null()
      class(lines_patch_t), pointer :: lines => null()
      integer, dimension(:), pointer, contiguous :: idx_map, idx_remap
      contains
        procedure :: pre_init => rate_patch_pre_init
        procedure :: update => rate_patch_update
        procedure :: post_update => rate_patch_post_update
        procedure :: nbor_relations => rate_patch_nbor_relations
        procedure :: integrate_rates
    end type rate_patch_t
  
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE rate_patch_pre_init (self)
    class(rate_patch_t), target:: self
    class(transition_type), pointer :: trans
    integer :: i, k, nrad, fs, fe
    logical :: skip
    ! ----------------------------------------------------------------------------
    call trace%begin('rate_patch_t%rate_patch_pre_init')
    call self%ray_patch_t%pre_init()
    call self%append_to_task_list()

    call os%allocate(self%idx_map, [self%atom%nrad], zero=.true.)
    nrad = 0
    do i=1, size(self%atom%line)
      trans => self%atom%line(i)
      call trans%check_range(self%fs,self%fe,fs,fe,skip)
      if (skip) cycle
      nrad = nrad+1
      self%idx_map(i) = nrad
    end do
    do i=1, size(self%atom%cont)
      trans => self%atom%cont(i)
      call trans%check_range(self%fs,self%fe,fs,fe,skip)
      if (skip) cycle
      nrad = nrad+1
      self%idx_map(i+self%atom%nline) = nrad
    end do
    call os%allocate(self%idx_remap, [nrad])
    do i=1, self%atom%nrad
      k = self%idx_map(i)
      if (k > 0) self%idx_remap(k) = i
    end do
    self%linecheck = self%idx_map(self%linecheck)

    call self%set_nv_and_nw(nrad, 2)

    call self%post_pre_init()
    call trace%end()
  END SUBROUTINE rate_patch_pre_init

  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji, Inu, Psi, Ieff
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) then
      return
    endif
    call trace%begin('rate_patch_t%rate_patch_update', itimer=itimer)

    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    if (.not. self%rad%stiff) then
      self%static_mem = .false.
    endif

    if (.not. self%static_mem) then
      Psi  => self%rad%Psi
      Ieff => self%rad%Ieff
      Inu  => self%rad%authority_mem(:,:,:,:,self%rad%it,self%rad%i_inu)
      temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)

      ! Zero the rates
      gamma_ij = 0.0d0
      gamma_ji = 0.0d0

      if (conditioning == 'full') then
        call self%atom%full_preconditioning(self%fs, temp, Ieff, Psi, gamma_ij, gamma_ji, self%idx_map, self%lines%line, self%lines%idx_map)
      elseif (conditioning == 'self') then
        if (self%i_mu==0) then
          ! Ieff is actually Jeff here
          call self%atom%self_cont_preconditioning(self%fs, temp, Ieff, Psi, gamma_ij, gamma_ji, self%idx_map)
        else
          call self%atom%self_line_preconditioning(self%fs, temp, Ieff, Psi, gamma_ij, gamma_ji, self%idx_map, self%lines%line)
        endif
      else
        print *, conditioning, ' is an unknown type of rate equation conditioning. Stopping...'
        stop
      endif
    endif

    call self%integrate_rates(gamma_ij, gamma_ji)

    call trace%end(itimer)
    call self%post_update()
  END SUBROUTINE rate_patch_update

  !===============================================================================
  SUBROUTINE rate_patch_post_update (self)
    class(rate_patch_t) :: self
    integer :: k, m(6), p(3)
    !-----------------------------------------------------------------------------
    if (self%verbose > 1) then
      if (self%static_mem) then
        if (self%check_in_bounds .or. self%verbose > 2) then
          write(stdout,'(A, A)') self%m3d_patch_info(), 'update skipped'
        endif
      elseif (self%check_in_bounds) then
        p = self%pcheck
        if (self%linecheck>0) then
          write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), &
          self%mem(p(1),p(2),p(3),self%linecheck,1,:)
        else
          write(stdout,'(A, A)') self%m3d_patch_info(), &
          'transition not active here'
        endif
      elseif (self%verbose > 2) then
        write(stdout,'(A,A)') self%m3d_patch_info()
      endif
    endif
    if (self%rad%stiff) then
      self%static_mem = .true.
      call self%rad%compress_mem(self)
    endif
    call self%m3d_patch_t%post_update()
    if (.not. self%static_mem) then
      call self%deallocate_mem(self)
    endif
  END SUBROUTINE rate_patch_post_update

  !===============================================================================
  SUBROUTINE integrate_rates (self, gamma_ij, gamma_ji)
    class(rate_patch_t) :: self
    integer :: m(4), i, k
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji
    !-----------------------------------------------------------------------------
    call self%m3d_parent%mem_lock(1)%set('integrate_rates')
    if (.not. associated(self%atom%gamma_ij)) then
      m(1:3) = self%n
      m(4) = self%atom%nrad
      call os%allocate(self%atom%gamma_ij, m, zero=.true.)
      call os%allocate(self%atom%gamma_ji, m, zero=.true.)
      self%atom%stiff=.true.
    endif
    if (self%i_mu>0) then
      do i=1, self%nv
        k = self%idx_remap(i)
        self%atom%gamma_ij(:,:,:,k) = self%atom%gamma_ij(:,:,:,k) + self%wmu * gamma_ij(:,:,:,i)
        self%atom%gamma_ji(:,:,:,k) = self%atom%gamma_ji(:,:,:,k) + self%wmu * gamma_ji(:,:,:,i)
      end do
    else
      do i=1, self%nv
        k = self%idx_remap(i)
        self%atom%gamma_ij(:,:,:,k) = self%atom%gamma_ij(:,:,:,k) + gamma_ij(:,:,:,i)
        self%atom%gamma_ji(:,:,:,k) = self%atom%gamma_ji(:,:,:,k) + gamma_ji(:,:,:,i)
      end do
    endif
    if (.not. self%static_mem) then
      self%atom%stiff=.false.
    endif
    call self%m3d_parent%mem_lock(1)%unset('integrate_rates')
  END SUBROUTINE integrate_rates
  
  !===============================================================================
  SUBROUTINE rate_patch_nbor_relations(self, another, needed, needs_me, download)
    class(rate_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    
    call trace%begin ('radiation_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type (another)
      class is (radiation_t)
        if (self%i_mu == another%i_mu) then
          needs_me=.false.
          self%rad => another
        else
          needed = .false.; needs_me = .false.
        endif
      class is (lines_patch_t)
        if (self%i_mu > 0) then
          needs_me=.false.
          self%lines => another
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        if (another%kind == 'atom_patch') then
          ! is needed and needs me
        else
          needed = .false.; needs_me = .false.
        endif
      end select
    endif
    call trace%end()
  END SUBROUTINE rate_patch_nbor_relations
  
  end module
