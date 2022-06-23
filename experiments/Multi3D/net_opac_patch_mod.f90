module net_opac_patch_mod
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
    USE m3d_patch_mod, only:maxiter, m3d_patch_t
    USE ray_patch_mod, only:ray_patch_t
    USE opac_patch_mod, only:opacity_t
    USE lines_patch_mod, only:lines_patch_t
    USE back_line_patch_mod, only:back_line_patch_t
    implicit none
    private
  
    !-----------------------------------------------------------------------------
    type, public, extends(ray_patch_t) :: net_opac_t
    contains
      procedure :: pre_init => net_opac_pre_init
      procedure :: update => net_opac_update
      procedure :: nbor_relations => net_opac_nbor_relations
    end type net_opac_t
  
    !-----------------------------------------------------------------------------
    integer, save:: verbose=0                           ! verbosity
      
  contains
  
  !===============================================================================
  SUBROUTINE net_opac_pre_init (self)
    class(net_opac_t), target:: self
    integer :: nw, fs, fe
    ! ----------------------------------------------------------------------------
    call trace%begin('net_opac_t%opacity_pre_init')
    call self%ray_patch_t%pre_init()
    call self%append_to_task_list()
  
    nw = 2
    call self%set_nv_and_nw (self%fe-self%fs+1,nw)
  
    if (self%out_ray) then
      self%ng = [0,0,1]
    else
      if (self%ray%long) then
        self%ng = [3,3,1]
      else
        self%ng = [2,2,1]
      endif
    endif
    call self%post_pre_init()
    call trace%end()
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
    integer :: iw,it, f, fs, fe, fn, x, y, z
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) then
      return
    endif

    call trace%begin('net_opac_t%net_opac_update', itimer=itimer)
    it = 1
    mem => self%authority_mem
    chi => self%authority_mem(:,:,:,:,it,self%i_chi)
    eta => self%authority_mem(:,:,:,:,it,self%i_eta)
    Snu => self%authority_mem(:,:,:,:,it,self%i_Snu)
    mem = 0
  
    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. .not. nbor%download) then
        nbtask => nbor%task
        select type (nbtask)
        class is (lines_patch_t)
          call nbtask%interpolate_lines(chi, eta)
        class is (opacity_t)
          fs = nbtask%fs - self%fs + 1
          fe = fs + nbtask%nv - 1
          do iw=1, self%nw
          do it=1, self%nt
          do f=fs,fe
            fn = f - fs + 1
            do z=1, self%n(3)
            do y=1, self%n(2)
            do x=1, self%n(1)
              mem(x,y,z,f,it,iw) = mem(x,y,z,f,it,iw) + nbtask%mem(x,y,z,fn,it,iw)
            end do
            end do
            end do
          end do
          end do
          end do
        class is (back_line_patch_t)
          fs = nbtask%fs - self%fs + 1
          fe = fs + nbtask%nv - 1
          do iw=1, self%nw
          do it=1, self%nt
          do f=fs,fe
            fn = f - fs + 1
            do z=1, self%n(3)
            do y=1, self%n(2)
            do x=1, self%n(1)
              mem(x,y,z,f,it,iw) = mem(x,y,z,f,it,iw) + nbtask%mem(x,y,z,fn,it,iw)
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
    do f=1, self%nv
      do z=1, self%n(3)
      do y=1, self%n(2)
      do x=1, self%n(1)
        Snu(x,y,z,f) = eta(x,y,z,f) / chi(x,y,z,f)
      end do
      end do
      end do
    end do
    
    mem1 => self%mem(:,:,:,:,self%it,:)
    call self%boundaries%mirror(mem1)
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
      
      class is (lines_patch_t)
        needs_me=.false.

      class is (back_line_patch_t)
        needs_me=.false.
        
      class is(opacity_t)
        if (another%kind/='net_opac') then
          needs_me=.false.
        endif
        
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.; needed = .false.
        else
          if (another%kind=='angle_rad') then
            needed = .false.
          else
            needed = .false.; needs_me = .false.
          endif
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif

    call trace%end()
  END SUBROUTINE net_opac_nbor_relations
  
  end module net_opac_patch_mod
  