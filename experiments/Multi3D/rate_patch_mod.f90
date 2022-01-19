module rate_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    !USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE solver_mod
    USE link_mod
    USE task_mod
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
    USE m3d_patch_mod
    USE opac_patch_mod
    USE rad_patch_mod
    USE atom_var, only: atom_type
    implicit none

    !-----------------------------------------------------------------------------

    type, public, extends(m3d_t) :: rate_patch_t
      integer :: i_Gij=1, i_Gji=2
      class(atom_type), pointer :: atom => null()
      real(8), dimension(:,:,:), allocatable :: vproj
      class(radiation_t), pointer :: rad
      class(opacity_t), pointer :: cont_opac, line_opac
      contains
        procedure :: pre_init => rate_patch_pre_init
        procedure :: init => rate_patch_init
        procedure :: update => rate_patch_update
        procedure :: nbor_relations => rate_patch_nbor_relations
    end type rate_patch_t
  
    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE rate_patch_pre_init (self)
    class(rate_patch_t), target:: self
    ! ----------------------------------------------------------------------------
    call self%m3d_t%pre_init()

    self%nv = 2

    self%mv = self%nv
    self%nw = self%atom%nrad
    self%nt = 1
    self%it = 1
    self%new = self%it
    
    ! ----------------------------------------------------------------------------
    ! Allocate and set unsigned and pervolume. Rates are all positive definite
    ! Although quantities are generally per volume, but density is not stored
    ! locally, thus this flag must not be set.
    ! ----------------------------------------------------------------------------
    allocate (self%unsigned(self%nv), self%pervolume(self%nv))
    self%pervolume = .false.
    self%unsigned  = .true.
  END SUBROUTINE rate_patch_pre_init
  
  !===============================================================================
  SUBROUTINE rate_patch_init (self)
    class(rate_patch_t) :: self
    integer :: i
    ! ----------------------------------------------------------------------------
    call self%m3d_t%init()
    
    self%t_offset = 0.9d0
  
    self%time = self%patch%time + self%t_offset*self%dtime
  END SUBROUTINE rate_patch_init
  
  !===============================================================================
  SUBROUTINE rate_patch_update (self)
    class(rate_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: Inu, gamma_ij, gamma_ji
    real(4), dimension(:,:,:,:), allocatable :: Ieff, eta
    integer :: fs, fe
    ! ----------------------------------------------------------------------------
    !call timer%tic ()

    fs = self%fs
    fe = self%fe

    Inu => self%rad%authority_mem(:,:,:,self%rad%i_inu,:)

    eta = self%cont_opac%mem(:,:,:,self%cont_opac%i_eta,self%it,:)
    eta = self%line_opac%mem(:,:,:,self%line_opac%i_eta,self%it,:) + eta

    Ieff = Inu - self%rad%Psi * eta

    gamma_ij => self%mem(:,:,:, self%i_Gij, self%it, :)
    gamma_ji => self%mem(:,:,:, self%i_Gji, self%it, :)

    call self%atom%add_all_gammas4(fs, fe, self%temp, self%vproj, &
    Ieff, self%rad%Psi, gamma_ij, gamma_ji)
  
    self%iter = self%iter + 1

    !select case (self%dir)
    !case(1)
    !  write(*,'(A10, 2(I4), 3(ES18.7))'), trim(self%kind), self%i_omega, self%i_nu, gamma_ij(1,1,1,4), gamma_ji(1,1,1,4)
    !case(-1)
    !  write(*,'(A10, 2(I4), 3(ES18.7))'), trim(self%kind), self%i_omega, self%i_nu, gamma_ij(1,1,self%n(3),4), gamma_ji(1,1,self%n(3),4)
    !end select
  
    !call timer%toc ('update took ', 1)
  END SUBROUTINE rate_patch_update
  
  !===============================================================================
  SUBROUTINE rate_patch_nbor_relations(self, another, needed, needs_me, download)
    class(rate_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    
    call trace%begin ('radiation_t%nbor_relations')
  
    ! If this RT task has no RT on, it should have no nbor relations
    if (.not.(self%on)) then
      needed = .false.; needs_me = .false.; download = .false.
      call trace%end()
      return
    end if
  
    ! First check the main relations, based on location and level of refinement
    call self%m3d_t%nbor_relations(another, needed, needs_me, download)

    if (another%id == self%rad%id) then
      needed = .true.; needs_me = .false.; download = .false.
    endif

    if (another%id == self%cont_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.
    endif

    if (another%id == self%line_opac%id) then
      needed = .true.; needs_me = .false.; download = .false.
    endif
  
    call trace%end()
  END SUBROUTINE rate_patch_nbor_relations
  
  end module