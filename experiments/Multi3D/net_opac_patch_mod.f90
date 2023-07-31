MODULE net_opac_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    !USE list_mod
    USE os_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    !USE index_mod
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    !USE scaling_mod
    USE trace_mod
    !USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    USE os_mod
    !USE units_mod
    !USE scratch_mod
    USE mmath, only: cc_kernel
    USE handyfuncs_mod
    USE ray_patch_mod, only:ray_patch_t
    USE back_cont_patch_mod, only:back_cont_t
    USE atom_cont_patch_mod, only:atom_cont_t
    USE atom_line_patch_mod, only:atom_line_t
    USE m3d_halo_mod, only: lock_mem, unlock_mem
    USE m3d_patch_mod, only:maxiter, m3d_patch_t, comoving_frame, linear_long, &
                            interp_log, ray_offset, pivot_z, dlam
    USE ray_mod, only: one_dimensional
    USE m3d_debug_mod
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: net_opac_t
    contains
      procedure :: pre_init => net_opac_pre_init
      procedure :: execute
      procedure :: rotate_cube
      procedure :: nbor_relations => net_opac_nbor_relations
      procedure :: init_nbors_hash
    end type net_opac_t
  
    !-----------------------------------------------------------------------------
    integer, save:: verbose=0                           ! verbosity
      
  contains
  
  !===============================================================================
  SUBROUTINE net_opac_pre_init (self)
    class(net_opac_t), target:: self
    integer :: nw, fs, fe
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('net_opac_t%opacity_pre_init', itimer=itimer)
    self%kind = 'net_opac'
    call self%ray_patch_t%pre_init()
  
    if (self%out_ray) then
      nw = 6
    else
      nw = 2
    endif
    call self%set_nv_and_nw (self%fe-self%fs+1,nw)

    if (comoving_frame .and. .not. self%out_ray) then
      self%ng = 0
    else
      call self%set_guard_zones()
    endif
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE net_opac_pre_init
  
  !===============================================================================
  SUBROUTINE execute (self)
    class(net_opac_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    real(4), dimension(:,:,:,:,:,:), pointer :: mem
    real(4), dimension(:,:,:,:), pointer :: chi, eta, Snu
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
    real(8), dimension(:), pointer, contiguous :: nu
    integer :: iw, it, f, fs, fe, fn, x, y, z, iiw, nw, p(3)
    integer, save :: itimer=0, itimer2=0
    !-----------------------------------------------------------------------------
    call trace%begin (trim(self%kind)//'%update', itimer=itimer)
    mem => self%authority_mem
    p = self%pcheck
    memcheck = 0

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. .not. nbor%download) then
        nbtask => nbor%task
        select type (nbtask)
        class is (back_cont_t)
          call lock_mem( nbtask, 'net_opac_update', read=.true.)
          call self%check_access(nbtask)
          fs = nbtask%fs - self%fs + 1
          fe = fs + nbtask%nv - 1
          nw = self%nw
          do iw=1, nw
            iiw = mod(iw-1,2)+1
          do it=1, self%nt
          do f=fs,fe
            fn = f - fs + 1
            do z=1, self%n(3)
            do y=1, self%n(2)
            do x=1, self%n(1)
              mem(x,y,z,f,it,iw) = nbtask%mem(x,y,z,fn,it,iiw)
            end do
            end do
            end do
          end do
          end do
          end do
          call unlock_mem( nbtask, 'net_opac_update', read=.true.)
        end select
      endif
      nbor => nbor%next
    end do

    nu  => self%sp%nu(self%fs:self%fe)
    if (self%out_ray .and. self%composition%ll_initialised) then
      chi => self%authority_mem(:,:,:,:,1,1)
      eta => self%authority_mem(:,:,:,:,1,2)
      call os%allocate(vproj, self%n)
      call self%comp_vproj(vproj)
      call trace%begin ('add_linelist_opacities', itimer=itimer2)
      call self%composition%add_linelist_opacities(nu, eta, chi, vproj, dlam)
      call trace%end (itimer2)
      call os%deallocate(vproj)
      do f=1, self%nv
        self%authority_mem(:,:,:,f,1,3) = chi(:,:,:,f)
        self%authority_mem(:,:,:,f,1,4) = eta(:,:,:,f)
      end do
    endif

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. .not. nbor%download) then
        nbtask => nbor%task
        select type (nbtask)
        class is (back_cont_t)
          ! already added
        class is (atom_line_t)
          it = self%it
          if (self%out_ray) then
            chi => self%authority_mem(:,:,:,:,it,1)
            eta => self%authority_mem(:,:,:,:,it,2)
            call nbtask%add_line_opacity(nu, chi, eta, LTE=.false.)
            chi => self%authority_mem(:,:,:,:,it,3)
            eta => self%authority_mem(:,:,:,:,it,4)
            call nbtask%add_line_opacity(nu, chi, eta, LTE=.true.)
          else
            chi => self%authority_mem(:,:,:,:,it,1)
            eta => self%authority_mem(:,:,:,:,it,2)
            call nbtask%interpolate_lines(chi, eta, LTE=.false.)
          endif
        class is (ray_patch_t)
          call lock_mem( nbtask, 'net_opac_update', read=.true.)
          call self%check_access(nbtask)
          fs = nbtask%fs - self%fs + 1
          fe = fs + nbtask%nv - 1
          nw = self%nw
          do iw=1, nw
            iiw = mod(iw-1,2)+1
          do it=1, self%nt
          do f=fs,fe
            fn = f - fs + 1
            do z=1, self%n(3)
            do y=1, self%n(2)
            do x=1, self%n(1)
              mem(x,y,z,f,it,iw) = mem(x,y,z,f,it,iw) + nbtask%mem(x,y,z,fn,it,iiw)
            end do
            end do
            end do
          end do
          end do
          end do
          call unlock_mem( nbtask, 'net_opac_update', read=.true.)
        end select
      endif
      nbor => nbor%next
    end do
  
    !$omp simd
    do iw=1, self%nw, 2
      do it=1, self%nt
        do f=1, self%nv
          !Snu = eta / chi
          mem(:,:,:,f,it,iw+1) = mem(:,:,:,f,it,iw+1) / mem(:,:,:,f,it,iw)
        end do
      end do
    end do
    
    call self%bndries()
    call self%periodic_grid()

    if (self%out_ray .and. .not. one_dimensional) then
      call self%rotate_cube()
    endif
    call trace%end(itimer)
  END SUBROUTINE execute

  !===============================================================================
  SUBROUTINE rotate_cube (self)
    class(net_opac_t) :: self
    real(4), dimension(:,:,:), contiguous, pointer :: kern
    real(4), dimension(:,:), contiguous, pointer :: shifted, plane
    integer, dimension(:,:), contiguous, pointer :: ipart
    real(4), dimension(2) :: delta, deltak, fpart
    integer, dimension(:,:,:), contiguous, pointer :: ix, iy
    real(4) :: distance2pivot, pmin, pmax, pup, p
    integer :: i, j, k, f, it, iw, up1, up2, nx, ny
    integer :: li, ui, ii, jj, x1, x2, y1, y2, kn
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('net_opac_t%rotate_cube', itimer=itimer)
    call os%allocate(shifted, self%n(1:2))
    call os%allocate(plane, self%gn(1:2))

    delta = self%vec(:2) / (self%vec(3) * self%ds(:2)) ! shift in units of grid cells
    distance2pivot = self%ds(3)*pivot_z - (self%position(3) - 0.5*self%size(3))

    nx = self%n(1)
    ny = self%n(2)
    li = self%li(3)
    ui = self%ui(3)

    if (linear_long) then
      kn = 2
    else
      kn = 4
    endif

    call os%allocate(ix, [kn, nx, ui])
    call os%allocate(iy, [kn, ny, ui])
    call os%allocate(kern,[kn, kn, ui])

    do k=li, ui
      deltak = (distance2pivot - (k-1)*self%ds(3)) * delta + ray_offset
      up1 = floor(deltak(1)) !integer part of shift
      up2 = floor(deltak(2)) !integer part of shift
      fpart = deltak - [up1, up2] !fractional part
      call set_interpolation_kernel(kern(:,:,k), fpart(1), fpart(2))

      do ii=1,kn
        ix(ii,:,k) = [(modulo(i+up1-kn/2+ii-1, nx)+1, i=1, nx)]
        iy(ii,:,k) = [(modulo(i+up2-kn/2+ii-1, ny)+1, i=1, ny)]
      end do
    end do

    do iw=1, self%nw
    do it=1, self%nt
    do f=1,  self%nv
      do k=li, ui
        plane = self%mem(:,:,k,f,it,iw)

        if (interp_log) plane = log(plane)

        if (linear_long) then
          do j=1, ny
            y1 = iy(1,j,k)
            y2 = iy(2,j,k)
        
            do i=1, nx
              x1 = ix(1,i,k)
              x2 = ix(2,i,k)

              ! indices might cross periodic boundary, hence no array sum possible
              shifted(i,j) = kern(1,1,k)*plane(x1,y1) + kern(2,1,k)*plane(x2,y1) &
                           + kern(1,2,k)*plane(x1,y2) + kern(2,2,k)*plane(x2,y2)
            end do
          end do
        else
          do j=1, ny
            do i=1, nx
              pmin = huge(1.0)
              pmax = tiny(1.0)
              pup = 0
              do jj=1,4
                do ii=1,4
                  p = plane(ix(ii,i,k), iy(jj,j,k))
                  if (kern(ii,jj,k) > 0) then
                    pup = pup + kern(ii,jj,k) * p
                    if (ii>1 .and. ii<4 .and. jj>1 .and. jj<4) then
                      pmin = min(pmin, p)
                      pmax = max(pmax, p)
                    endif
                  endif
                end do
              end do
              shifted(i,j) = min(max(pup, pmin), pmax)
            end do
          end do
        endif

        if (interp_log) then
          self%authority_mem(:,:,k-li+1,f,it,iw) = exp(shifted)
        else
          self%authority_mem(:,:,k-li+1,f,it,iw) = shifted
        endif
      end do
    end do
    end do
    end do

    call os%deallocate(shifted)
    call os%deallocate(plane)
    call os%deallocate(kern)
    call os%deallocate(ix)
    call os%deallocate(iy)
    call trace%end(itimer)
  END SUBROUTINE rotate_cube
  
  !===============================================================================
  SUBROUTINE net_opac_nbor_relations(self, another, needed, needs_me, download)
    class(net_opac_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('net_opac_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type(another)
      class is (net_opac_t)
        ! do nothing
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.; needed = .false.
        else
          select case(another%kind)
          case ('atom_line', 'back_cont', 'atom_cont', 'atom_patch')
            needs_me = .false.
          case ('angle_rad')
            if (self%out_ray .or. .not. comoving_frame) then
              needed = .false.
            else
              needed = .false.; needs_me = .false.
            endif
          case ('proj_patch', 'mean_rad')
            needed = .false.
          case default
            needed = .false.; needs_me = .false.
          end select
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif

    call trace%end()
  END SUBROUTINE net_opac_nbor_relations
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(net_opac_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: key(4), i, mu_start, mu_stop
    !-----------------------------------------------------------------------------
    call trace%begin('net_opac_t%init_nbors_hash', itimer=itimer)
    if (comoving_frame .and. .not. self%out_ray) then
      call self%m3d_patch_t%init_nbors_hash(nbors)
    else
      call self%ray_patch_t%init_nbors_hash(nbors)
    endif

    key(1:3) = self%key(1:3)
    key(4) = hash_string('atom_line')
    call self%check_nbor(key, nbors)

    mu_start = merge(-self%n_out, 1, self%out_ray)
    mu_stop  = merge(-1, self%n_mu, self%out_ray)

    if (comoving_frame .and. .not. self%out_ray) then
      key(4) = hash_string('mean_rad')
      call self%check_nbor(key, nbors)

      do i=mu_start, mu_stop
        key(3) = i
        key(4) = hash_string('proj_patch')
        call self%check_nbor(key, nbors)
      end do
    else
      key(4) = hash_string('angle_rad')
      call self%check_nbor(key, nbors)
    endif
    
    key(3) = 0
    key(4) = hash_string('back_cont')
    call self%check_nbor(key, nbors)

    key(4) = hash_string('atom_cont')
    call self%check_nbor(key, nbors)

    ! because of atom_lines
    key(2) = 0
    key(4) = hash_string('atom_patch')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE set_interpolation_kernel(kern,h1,h2)
    real(4),intent(in) :: h1,h2
    real(4),dimension(:,:),intent(out) :: kern
    real(8),dimension(4) :: xkern,ykern
    integer :: i,j
    !-----------------------------------------------------------------------------
    if (linear_long) then
      ! bilinear interpolation to avoid overshooting
      kern(1,1) = (1.0-h1)*(1.0-h2)
      kern(2,1) =      h1 *(1.0-h2)
      kern(1,2) = (1.0-h1)*     h2 
      kern(2,2) =      h1 *     h2
    else
      call cc_kernel(dble(h1),xkern)
      call cc_kernel(dble(h2),ykern)
      forall (i=1:4,j=1:4)
        kern(i,j) = xkern(i)*ykern(j)
      endforall
    endif
  END SUBROUTINE set_interpolation_kernel
  
  END MODULE net_opac_patch_mod
  