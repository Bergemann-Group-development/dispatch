MODULE out_rad_patch_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE io_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    USE math_mod
    USE bits_mod
    USE handyfuncs_mod
    USE m3d_patch_mod, only: m3d_patch_t, long_vec, long_wmu
    USE ray_patch_mod, only: ray_patch_t
    USE mconst, only: cc
    USE os_mod
    USE omp_mod
    USE omp_lock_mod
    implicit none
  
    type, public, extends(ray_patch_t) :: out_rad_patch_t
    contains
      procedure :: pre_init => out_rad_pre_init
      procedure :: update => out_rad_update
      procedure :: nbor_relations
      procedure :: save_quadrature
      procedure :: print_post_update_info
      procedure :: init_nbors_hash
    end type out_rad_patch_t
  
  contains
    
  !===============================================================================
  SUBROUTINE out_rad_pre_init (self)
    class(out_rad_patch_t), target:: self
    logical, save :: first_time=.true.
    integer :: iostat, nv, nw
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('out_rad_patch_t%radiation_pre_init', itimer=itimer)
    self%kind = 'out_rad'
    call self%ray_patch_t%pre_init()
  
    self%n(3) = 1
    nw = 3
    call self%set_nv_and_nw (self%sp%nnu, nw)
    call self%save_quadrature()
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE out_rad_pre_init
  
  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(out_rad_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    ! ----------------------------------------------------------------------------
    call trace%begin ('out_rad_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
    if (needed .or. needs_me) then
      if (another%kind=='angle_rad') then
          needs_me=.false.
      else
          needed=.false. ; needs_me=.false.
      endif
    endif
    call trace%end()
  END SUBROUTINE nbor_relations
  
  !===============================================================================
  SUBROUTINE out_rad_update(self)
    class(out_rad_patch_t) :: self
    class(nbor_t), pointer:: nbor
    character(len=128) :: filename
    real(8), allocatable :: lam(:)
    real(4), dimension(:,:,:), pointer :: Inu
    real(4), dimension(:,:,:,:), pointer :: mem
    real(4), dimension(:,:,:,:,:), pointer, contiguous :: II
    integer :: lunit
    integer :: i, iw
    real(4) :: norm
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('out_rad_patch_t%save_radiation_output')
    write(stdout, '(A, I3)') 'Writing output for angle ', abs(self%i_mu)
  
    nbor => self%link%nbor
    do while (associated(nbor))
      select type(nb => nbor%task)
      class is (ray_patch_t)
        call self%check_access(nb)
        do iw=1, self%nw
          Inu => self%mem(:,:,1, nb%fs:nb%fe, self%it, iw)
          do i=1, nb%nv
            Inu(:,:,i) = nb%authority_mem(:,:,nb%n(3), i, nb%it, iw)
          end do
        end do
      end select
      nbor => nbor%next
    end do

    call os%allocate(II, [1,1,1,self%nv,self%nw])
    norm = 1.0 / product(self%n(1:2))
    II(1,1,1,:,:) = sum(sum(self%mem(:,:,1, :, self%it, :), 1),1) * norm
    call self%save(II, 'ie')
    call os%deallocate(II)

    call trace%end()
  END SUBROUTINE out_rad_update
  
  !===============================================================================
  SUBROUTINE save_quadrature(self)
    class(out_rad_patch_t),intent(inout):: self
    integer :: meta_unit, iostat, i, nmu
    character(len=128) :: meta_file
    ! ----------------------------------------------------------------------------
    meta_file = trim(io%run_dir) // 'angle_quadrature.txt'
    call io%lock%set('out_rad')
    if (.not. os%exists(meta_file)) then
      open(newunit=meta_unit, file=meta_file, form='formatted')
      nmu = size(long_wmu)
      write(meta_unit, '(A3, 4(A12))') '#', 'weight', 'mux', 'muy', 'muz'
      do i=1, nmu
        write(meta_unit, '(I3, F12.8, SP, 3(F12.8))') i, long_wmu(i), long_vec(:,i)
      end do
      close(meta_unit)
    endif
    call io%lock%unset('out_rad')
  END SUBROUTINE
  
  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(out_rad_patch_t) :: self
    character(48) :: text
    integer :: p(3), f
    real(4) :: norm
    !-----------------------------------------------------------------------------
    f = self%fcheck
    norm = product(self%n(1:2))
    call self%print_to_stdout(arr=[sum(self%mem(:,:,1,f,1,1)) / norm])
  END SUBROUTINE print_post_update_info
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(out_rad_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: i_mu, i_nu, key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('out_rad_patch_t%init_nbors_hash', itimer=itimer)
    key(1) = self%key(1)
    key(3) = self%key(3)
    do i_nu=-self%n_nu, -1
      key(2) = i_nu
      key(4) = hash_string('angle_rad')
      call self%check_nbor(key, nbors)
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE out_rad_patch_mod
  