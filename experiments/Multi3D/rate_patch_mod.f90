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
    call self%set_nv_and_nw(self%atom%nrad, 2)
    call trace%end()
  END SUBROUTINE rate_patch_pre_init
  
  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: Inu, gamma_ij, gamma_ji, chi_l, chi_c, eta_l, eta_c, Psi
    real(4), dimension(:,:,:,:), contiguous, pointer :: Ieff
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    integer :: fs, fe, i,j,k,f
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) then
      return
    endif
    
    call trace%begin('rate_patch_t%rate_patch_update', itimer=itimer)

    fs = self%fs
    fe = self%fe

    Inu => self%rad%authority_mem(:,:,:,:,self%rad%it,self%rad%i_inu)
    Psi => self%rad%Psi

    chi_c => self%cont_opac%mem(:,:,:,:,self%it,self%cont_opac%i_chi)
    eta_c => self%cont_opac%mem(:,:,:,:,self%it,self%cont_opac%i_eta)

    chi_l => self%line_opac%mem(:,:,:,:,self%it,self%line_opac%i_chi)
    eta_l => self%line_opac%mem(:,:,:,:,self%it,self%line_opac%i_eta)

    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)
    
    call os%allocate(Ieff, shape(Psi))
    
    do f=1, self%rad%nv
    do k=1, self%n(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      Ieff(i,j,k,f) = Inu(i,j,k,f) - Psi(i,j,k,f) * (eta_c(i,j,k,f) + eta_l(i,j,k,f))
    end do
    end do
    end do
    end do

    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    if (self%iter == 0) then
      call self%comp_vproj(self%vproj)
    endif

    call self%atom%add_all_gammas(fs, fe, temp, self%vproj, Ieff, Psi, &
                                  chi_l, chi_c, eta_l, eta_c, gamma_ij, gamma_ji)

    call os%deallocate(Ieff)
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
      select case (another%kind)
      case ('angle_rad', 'cont_opac', 'line_opac')
        needs_me=.false.
        
        select type(another)
        class is (radiation_t)
        self%rad => another
        end select
      case ('atom_patch')
        ! is needed and needs me
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE rate_patch_nbor_relations
  
  end module
