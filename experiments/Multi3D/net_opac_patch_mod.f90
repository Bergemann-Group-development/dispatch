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
    USE m3d_patch_mod, only:maxiter, m3d_patch_t, hash_string, comoving_frame
    USE ray_patch_mod, only:ray_patch_t
    USE back_cont_patch_mod, only:back_cont_patch_t
    USE atom_cont_patch_mod, only:atom_cont_patch_t
    USE atom_line_patch_mod, only:atom_line_t
    USE back_line_patch_mod, only:back_line_patch_t
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: net_opac_t
    contains
      procedure :: pre_init => net_opac_pre_init
      procedure :: update => net_opac_update
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
      if (self%out_ray) then
        self%ng = [0,0,1]
      else
        if (self%ray%long) then
          self%ng = [3,3,1]
        else
          self%ng = [2,2,1]
        endif
      endif
    endif
    call self%post_pre_init()
    call trace%end(itimer)
  END SUBROUTINE net_opac_pre_init
  
  !===============================================================================
  SUBROUTINE net_opac_update (self)
    class(net_opac_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    real(4), dimension(:,:,:,:,:,:), pointer :: mem
    real(4), dimension(:,:,:,:,:), pointer :: mem1
    real(4), dimension(:,:,:,:), pointer :: chi, eta, Snu
    integer, save:: itimer=0
    integer :: iw,it, f, fs, fe, fn, x, y, z, iiw, nw
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    !-----------------------------------------------------------------------------
    ! Boundary conditions
    !-----------------------------------------------------------------------------
    call self%bndries()
    call self%periodic_grid(only=merge(1, self%gn(3), self%dir>0))
    call self%decompress_mem(self)

    call trace%begin('net_opac_t%net_opac_update', itimer=itimer)
    mem => self%authority_mem

    do iw=1, self%nw
      do it=1, self%nt
        do f=1, self%nv
          mem(:,:,:,f,it,iw) = 0
        end do
      end do
    end do

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. .not. nbor%download) then
        nbtask => nbor%task
        call self%check_access(nbtask)
        select type (nbtask)
        class is (atom_line_t)
          it = self%it
          chi => self%authority_mem(:,:,:,:,it,1)
          eta => self%authority_mem(:,:,:,:,it,2)
          call nbtask%interpolate_lines(chi, eta, LTE=.false.)
          if (self%out_ray) then
            chi => self%authority_mem(:,:,:,:,it,3)
            eta => self%authority_mem(:,:,:,:,it,4)
            call nbtask%interpolate_lines(chi, eta, LTE=.true.)
          endif
        class is (ray_patch_t)
          fs = nbtask%fs - self%fs + 1
          fe = fs + nbtask%nv - 1
          if (self%out_ray) then
            nw = merge(4, self%nw, nbtask%kind=='back_line')
          else
            nw = self%nw
          endif
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
        end select
      endif
      nbor => nbor%next
    end do
  
    !$omp simd
    do iw=1, self%nw, 2
      do it=1, self%it
        do f=1, self%nv
            !Snu(:,:,:,f) = eta(:,:,:,f) / chi(:,:,:,f)
            mem(:,:,:,f,it,iw+1) = mem(:,:,:,f,it,iw+1) / mem(:,:,:,f,it,iw)
        end do
      end do
    end do
    
    mem1 => self%mem(:,:,:,:,self%it,:)
    if (self%ng(3) > 0) then ! not called when comoving frame
      call self%boundaries%mirror(mem1)
    endif
    if (self%is_periodic()) then
      call self%periodic_grid()
    end if
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE net_opac_update
  
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
          case ('atom_line', 'back_cont', 'back_line', 'atom_cont', 'atom_patch')
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

    mu_start = merge(self%n_omega+1, 1, self%out_ray)
    mu_stop  = merge(self%n_mu, self%n_omega, self%out_ray)

    if (self%out_ray) then
      key(4) = hash_string('back_line')
      call self%check_nbor(key, nbors)
    endif

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

    ! because of atom_linees
    key(2) = 0
    key(4) = hash_string('atom_patch')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash
  
  END MODULE net_opac_patch_mod
  