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
    USE atom_mod, only: atom_type
    USE lines_patch_mod, only: lines_patch_t
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(ray_patch_t) :: rate_patch_t
      class(atom_type), pointer :: atom => null()
      class(radiation_t), pointer :: rad => null()
      class(opacity_t), pointer :: cont_opac => null()
      class(opacity_t), pointer :: line_opac => null()
      class(lines_patch_t), pointer :: lines => null()
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
    call self%post_pre_init()
    call trace%end()
  END SUBROUTINE rate_patch_pre_init

  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: gamma_ij, gamma_ji, Inu, Psi
    real(4), dimension(:,:,:), pointer, contiguous :: temp
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) then
      return
    endif
    call trace%begin('rate_patch_t%rate_patch_update', itimer=itimer)

    Psi => self%rad%Psi
    Inu => self%rad%authority_mem(:,:,:,:,self%rad%it,self%rad%i_inu)
    temp => self%m3d_parent%m3d_parent%mem(:,:,:,self%i_temp, self%it,1)

    gamma_ij => self%mem(:,:,:,:, self%it, self%i_Gij)
    gamma_ji => self%mem(:,:,:,:, self%it, self%i_Gji)

    ! Zero the rates etc
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0

    if (conditioning == 'full') then
      call self%atom%full_preconditioning(self%fs, temp, Inu, Psi, gamma_ij, gamma_ji, self%lines%line, self%lines%idx_map)
    elseif (conditioning == 'self') then
      call self%atom%self_preconditioning(self%fs, temp, Inu, Psi, gamma_ij, gamma_ji, self%lines%line)
    else
      print *, conditioning, ' is an unknown type of rate equation conditioning. Stopping...'
      stop
    endif

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
      class is (lines_patch_t)
        needs_me=.false.
        self%lines => another
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
