MODULE proj_patch_mod
    USE os_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE bits_mod
    USE os_mod
    USE handyfuncs_mod
    USE m3d_patch_mod, only:maxiter, m3d_patch_t, use_vel
    USE ray_patch_mod, only:ray_patch_t
    USE net_opac_patch_mod, only:net_opac_t
    USE op_math_mod, only:interp1d_linear
    USE mconst, only: CC
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: proj_patch_t
      class(ray_patch_t), pointer :: net_opac
    contains
      procedure :: pre_init
      procedure :: update
      procedure :: nbor_relations
      procedure :: init_nbors_hash
    end type proj_patch_t

    public :: interpolate_to_ray
      
  contains
  
  !===============================================================================
  SUBROUTINE pre_init (self)
    class(proj_patch_t), target:: self
    integer :: nw, fs, fe
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('proj_patch_t%pre_init', itimer=itimer)
    self%kind = 'proj_patch'
    call self%ray_patch_t%pre_init()
    call self%set_nv_and_nw (self%fe-self%fs+1, 2)
    call self%set_guard_zones()
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE pre_init
  
  !===============================================================================
  SUBROUTINE update (self)
    class(proj_patch_t) :: self
    real(4), dimension(:,:,:,:), pointer :: yy, new_yy
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
    real(8), dimension(:), pointer, contiguous :: nu_in, nu_out
    integer, save:: itimer=0
    integer :: iw, it, f
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('proj_patch_t%update', itimer=itimer)
    call self%check_access(self%net_opac)
    
    if (use_vel) then
      if (.not. self%sp%sorted) call io%abort('Comoving frame for unsorted spectrum!')
      call os%allocate(vproj, self%n)
      call self%comp_vproj(vproj)
      vproj = 1 / (1.0 + vproj)
      nu_in => self%sp%nu(self%net_opac%fs:self%net_opac%fe)
      nu_out => self%sp%nu(self%fs:self%fe)

      do iw=1, self%nw
      do it=1, self%nt
        yy => self%net_opac%authority_mem(:,:,:,:,it,iw)
        new_yy => self%authority_mem(:,:,:,:,it,iw)
        call interpolate_to_ray(vproj, nu_in, nu_out, yy, new_yy)
      end do
      end do
      call os%deallocate(vproj)
    else
      do iw=1, self%nw
      do it=1, self%nt
      do f=1,  self%nv
        self%authority_mem(:,:,:,f,it,iw) = self%net_opac%authority_mem(:,:,:,f,it,iw)
      end do
      end do
      end do
    endif

    call self%boundaries%mirror()
    call self%periodic_grid()
    call trace%end(itimer)
  END SUBROUTINE update
  
  !===============================================================================
  SUBROUTINE interpolate_to_ray(proj, nu_in, nu_out, yy, new_yy, weight)
    real(4), intent(in), dimension(:,:,:), contiguous :: proj
    real(4), intent(in), dimension(:,:,:,:) :: yy
    real(8), intent(in), dimension(:) :: nu_in, nu_out
    real(4), intent(in), optional :: weight
    real(4), intent(inout), dimension(:,:,:,:) :: new_yy
    integer, dimension(:,:,:), pointer, contiguous :: ff
    integer :: n(3), f, nf, i, x, y, z, nf_in
    real(4) :: yy1, yy2, wmu, res
    real(8) :: new_nu, dnu
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin('proj_patch_t%interpolate_to_ray', itimer=itimer)
    n = shape(proj)
    nf = size(nu_out)
    nf_in = size(nu_in)
    call os%allocate(ff, n, zero=.true.)

    if (present(weight)) then
      wmu = weight
    else
      wmu = 0
    endif

    do f=1, nf
      do z=1, n(3)
      do y=1, n(2)
      do x=1, n(1)
        new_nu = nu_out(f) * dble(proj(x,y,z))
        if (new_nu <= nu_in(1)) then
          res = yy(x,y,z,1)
        elseif (new_nu >= nu_in(nf_in)) then
          res = yy(x,y,z,nf_in)
        else
          i = ff(x,y,z)
          do while (i < nf_in .and. new_nu >= nu_in(i+1))
            i = i + 1
          end do
          yy1 = yy(x,y,z,i)
          yy2 = yy(x,y,z,i+1)
          dnu = (new_nu - nu_in(i)) / (nu_in(i+1) - nu_in(i))
          res = (yy2 - yy1) * dnu + yy1
          ff(x,y,z) = i
        endif
        if (wmu == 0) then
          new_yy(x,y,z,f) = res
        else
          new_yy(x,y,z,f) = new_yy(x,y,z,f) + wmu * res
        endif
      end do
      end do
      end do
    end do
    call os%deallocate(ff)
    call trace%end(itimer)
  END SUBROUTINE interpolate_to_ray
  
  !===============================================================================
  SUBROUTINE nbor_relations(self, another, needed, needs_me, download)
    class(proj_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('proj_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type(another)
      class is (ray_patch_t)
        if (another%kind == 'net_opac') then
            needs_me = .false.
            self%net_opac => another
        elseif (another%kind == 'angle_rad') then
            needed = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif

    call trace%end()
  END SUBROUTINE nbor_relations
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(proj_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('proj_patch_t%init_nbors_hash', itimer=itimer)
    call self%ray_patch_t%init_nbors_hash(nbors)
  
    key(1:3) = self%key(1:3)
    key(4) = hash_string('angle_rad')
    call self%check_nbor(key, nbors)

    key(3) = 0
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE proj_patch_mod
  