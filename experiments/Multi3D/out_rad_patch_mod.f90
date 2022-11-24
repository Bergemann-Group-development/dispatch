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
    USE m3d_patch_mod, only: m3d_patch_t, hash_string
    USE ray_patch_mod, only: ray_patch_t
    USE angle_rad_patch_mod, only: rad_patch_t
    USE mconst, only: cc
    USE os_mod
    USE omp_mod
    USE omp_lock_mod
    implicit none
  
    type, public, extends(ray_patch_t) :: out_rad_patch_t
      character(len=128) :: meta_file, ray_name
    contains
      procedure :: pre_init => out_rad_pre_init
      procedure :: update => out_rad_update
      procedure :: nbor_relations
      procedure :: setup_metafile
      procedure :: save_meta
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
    call self%setup_metafile()
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
    integer :: lunit
    integer :: i, iw
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('out_rad_patch_t%save_radiation_output')
    write(stdout, '(A64)') 'Writing output for angle ' // self%ray_name
  
    if (self%i_out == 1) then
      lam = 1d8 * dble(cc) / self%sp%nu
      filename = trim(io%run_dir) // 'save/' // 'lam.bin'
      open(newunit=lunit, file=filename, access='stream')
      write(lunit) lam
      close(lunit)
    endif
  
    nbor => self%link%nbor
    do while (associated(nbor))
      select type(nb => nbor%task)
      class is (rad_patch_t)
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
  
    mem => self%mem(:,:,1,:, self%it, :)
    filename = trim(self%ray_name) // '_inu.bin'
    call self%save_meta(filename, storage_size(mem), shape(mem))
    filename = trim(io%run_dir) // 'save/' // trim(self%ray_name) // '_inu.bin'
    open(newunit=lunit, file=filename, access='stream')
    write(lunit) mem
    close(lunit)
    call trace%end()
    call self%post_update()
  END SUBROUTINE out_rad_update
  
  !===============================================================================
  SUBROUTINE setup_metafile(self)
    class(out_rad_patch_t),intent(inout):: self
    integer :: meta_unit, iostat
    ! ----------------------------------------------------------------------------
    write(self%ray_name,'(A3, I3.3)') 'ray', self%i_out
  
    self%meta_file = trim(io%run_dir) // 'save/' // trim(self%ray_name) // '_meta.txt'
    open(newunit=meta_unit, file=self%meta_file, form='formatted')
  
    write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'vec', self%vec
    write(meta_unit, '(A8, 1X, SP, F12.8)') 'weight', self%wmu
    close(meta_unit)
  END SUBROUTINE
  
  !===============================================================================
  SUBROUTINE save_meta(self, filename, store, shpe)
    class(out_rad_patch_t):: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: store
    integer, dimension(:), intent(in) :: shpe
    integer :: meta_unit
    ! ----------------------------------------------------------------------------
    open(newunit=meta_unit, file=self%meta_file, form='formatted', position='append')
    write(meta_unit, '(A64, I4, 8(I10))') filename, store, shpe
    close(meta_unit)
  END SUBROUTINE save_meta
  
  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(out_rad_patch_t) :: self
    character(48) :: text
    integer :: p(3), f
    !-----------------------------------------------------------------------------
    f = self%fcheck
    call self%print_to_stdout(arr=[sum(self%mem(:,:,1,f,1,1)) / product(self%n(1:2))])
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
    do i_nu=self%n_nu, 2*self%n_nu
    key(2) = i_nu
    key(4) = hash_string('angle_rad')
    call self%check_nbor(key, nbors)
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE out_rad_patch_mod
  