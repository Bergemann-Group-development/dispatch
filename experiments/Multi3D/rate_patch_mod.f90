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
    USE ray_patch_mod
    USE opac_patch_mod
    USE rad_patch_mod
    USE atom_var, only: atom_type
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(ray_patch_t) :: rate_patch_t
      class(atom_type), pointer :: atom => null()
      class(radiation_t), pointer :: rad => null()
      class(opacity_t), pointer :: cont_opac => null()
      class(opacity_t), pointer :: line_opac => null()
      contains
        procedure :: pre_init => rate_patch_pre_init
        procedure :: update => rate_patch_update
        procedure :: nbor_relations => rate_patch_nbor_relations
    end type rate_patch_t
  
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE rate_patch_pre_init (self)
    class(rate_patch_t), target:: self
    ! ----------------------------------------------------------------------------
    call trace%begin('rate_patch_t%rate_patch_pre_init')
    call self%ray_patch_t%pre_init()
    call self%append_to_task_list()
    call self%set_nv_and_nw(self%atom%nrad, 2)
    call trace%end()
  END SUBROUTINE rate_patch_pre_init
  
  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    use mconst, only:cc
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji
    real(4), dimension(:,:,:), pointer :: Inu, chi_l, chi_c, eta_l, eta_c, Psi
    real(4), dimension(:,:,:), pointer, contiguous :: temp, zero
    integer :: f, f2, ff
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) then
      return
    endif
    
    call trace%begin('rate_patch_t%rate_patch_update', itimer=itimer)

    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    ! Zero the rates etc
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0
    
    if (self%iter == 0) then
      call self%comp_vproj(self%vproj)
      self%vproj = (1 + self%vproj/cc) / self%atom%doppler_factor
    endif

    call os%allocate(zero, shape(temp))
    zero = 0.0

    if (.not. associated(self%cont_opac)) then
      chi_c => zero
      eta_c => zero
    endif

    if (.not. associated(self%line_opac)) then
      chi_l => zero
      eta_l => zero
    endif

    do f=1, self%rad%nv
      Inu => self%rad%authority_mem(:,:,:,f,self%rad%it,self%rad%i_inu)
      Psi => self%rad%Psi(:,:,:,f)

      ff = self%fs + f - 1

      if (associated(self%cont_opac)) then
        if (ff >= self%cont_opac%fs .and. ff <= self%cont_opac%fe) then
          f2 = ff - self%cont_opac%fs + 1
          chi_c => self%cont_opac%mem(:,:,:,f2,self%it,self%cont_opac%i_chi)
          eta_c => self%cont_opac%mem(:,:,:,f2,self%it,self%cont_opac%i_eta)
        else
          chi_c => zero
          eta_c => zero
        endif
      endif
  
      if (associated(self%line_opac)) then
        if (ff >= self%line_opac%fs .and. ff <= self%line_opac%fe) then
          f2 = ff - self%line_opac%fs + 1
          chi_l => self%line_opac%mem(:,:,:,f2,self%it,self%line_opac%i_chi)
          eta_l => self%line_opac%mem(:,:,:,f2,self%it,self%line_opac%i_eta)
        else
          chi_l => zero
          eta_l => zero
        endif
      endif

      call self%atom%add_all_gammas(ff, temp, self%vproj, Inu, Psi, &
           chi_l, chi_c, eta_l, eta_c, gamma_ij, gamma_ji)
    end do

    call os%deallocate(zero)
    call trace%end(itimer)
    if (.not. self%finished) call self%post_update()
  END SUBROUTINE rate_patch_update
  
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
        if (another%kind == 'angle_rad') then
          needs_me=.false.
          self%rad => another
        else
          needed = .false.; needs_me = .false.
        endif
      class is (opacity_t)
        select case (another%kind)
        case ('cont_opac')
          needs_me=.false.
          self%cont_opac => another
        case ('line_opac')
          needs_me=.false.
          self%line_opac => another
        case default
          needed = .false.; needs_me = .false.
        end select
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
