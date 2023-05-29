MODULE back_line_patch_mod
    USE os_mod
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
    USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE mconst
    USE handyfuncs_mod
    USE m3d_patch_mod, only: m3d_patch_t
    USE ray_patch_mod, only: ray_patch_t
    USE composition_mod, only: composition_type
    USE element_mod, only: element_type, isotope_type
    USE ion_mod, only: ion_type
    USE molecule_mod, only: molec_type
    USE bbtrans_mod, only: bbtrans_type
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(ray_patch_t) :: back_line_patch_t
      !class(composition_type), pointer :: composition
      integer :: nline, ired, iblue
    contains
      procedure :: pre_init => back_line_patch_pre_init
      procedure :: update => back_line_patch_update
      procedure :: nbor_relations => back_line_patch_nbor_relations
      procedure :: init_nbors_hash
    end type back_line_patch_t

    !-----------------------------------------------------------------------------
  
  contains

  !===============================================================================
  SUBROUTINE back_line_patch_pre_init (self)
    class(back_line_patch_t), target :: self
    integer :: nv, nw
    integer, save :: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin('back_line_t%pre_init', itimer=itimer)
    self%kind = 'back_line'
    call self%ray_patch_t%pre_init()

    nw = 2
    nv = self%fe-self%fs+1
    call self%set_nv_and_nw (nv,nw)
    self%ng = 0
    !call self%mask_linelist()
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE back_line_patch_pre_init

  !===============================================================================
  SUBROUTINE back_line_patch_update (self)
    class(back_line_patch_t) :: self
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
    real(4), dimension(:,:,:,:), pointer, contiguous :: eta, chi
    real(8), dimension(:), pointer, contiguous :: nu
    integer, save:: itimer=0, itimer2=0
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('back_line_patch_t%back_line_patch_update', itimer=itimer)

    eta => self%mem(:,:,:,:,1,self%i_eta)
    chi => self%mem(:,:,:,:,1,self%i_chi)
    nu  => self%sp%nu(self%fs:self%fe)
    
    call os%allocate(vproj, self%n)
    call self%comp_vproj(vproj)
    call self%composition%add_linelist_opacities(nu, eta, chi, vproj)
    call os%deallocate(vproj)
    call trace%end(itimer)
  END SUBROUTINE back_line_patch_update

  !===============================================================================
  SUBROUTINE back_line_patch_nbor_relations(self, another, needed, needs_me, download)
    class(back_line_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                            ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_line_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type(another)
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          needed = .false.
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE back_line_patch_nbor_relations
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(back_line_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: i, i_nu, i_mu, key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('back_line_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)
  
    key(1:3) = self%key(1:3)
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE
