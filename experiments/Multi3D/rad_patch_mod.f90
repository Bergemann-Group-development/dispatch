MODULE rad_patch_mod
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
  !USE mpi_mod
  !USE omp_mod
  !USE omp_timer_mod
  USE math_mod
  USE bits_mod
  USE units_mod, only: cgs
  USE mconst, only: over_4pi
  !USE scratch_mod
  USE m3d_patch_mod, only: m3d_patch_t, hash_string, comoving_frame
  USE ray_patch_mod, only: ray_patch_t, formal_solver
  USE proj_patch_mod, only:interpolate_to_ray
  USE os_mod
  USE omp_mod
  USE omp_lock_mod
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(ray_patch_t) :: rad_patch_t
    class(ray_patch_t), pointer :: opac => null()
    class(ray_patch_t), pointer :: back_opac => null()
    class(rad_patch_t), pointer :: mean_rad => null()
    real(4), dimension(:,:,:,:), contiguous, pointer :: Inu, chi, Snu
    real(4), dimension(:,:,:,:), contiguous, pointer :: Psi, Ieff
    real(4), dimension(:,:,:,:), contiguous, pointer :: Inu_old
    integer :: stiff_count=1
    logical :: stiff=.false.
  contains
    procedure :: pre_init => radiation_pre_init
    procedure :: init => radiation_init
    procedure :: update => radiation_update
    procedure :: nbor_relations => radiation_nbor_relations
    !procedure :: long_ray_update
    procedure :: long_char_update
    procedure :: short_ray_update
    procedure :: Jnu_update
    procedure :: bndries
    procedure :: periodic_plane
    procedure :: interp_long
    procedure :: integrate_intensity
    !procedure :: interp_xy_up
    !procedure :: interp_xy_up_plane
    !procedure :: solve_long_plane
    procedure :: remove_lower_layers_mem
    procedure :: deallocate_mem => rad_patch_deallocate_mem
    procedure :: allocate_mem => rad_patch_allocate_mem
    procedure :: print_post_update_info
    procedure :: init_nbors_hash
  end type rad_patch_t

  !-----------------------------------------------------------------------------

  real(4), save :: stiff_lim=1e-8
  logical, public, save :: lock_Inu=.false.

contains
  
!===============================================================================
SUBROUTINE radiation_pre_init (self)
  class(rad_patch_t), target:: self
  logical, save :: first_time=.true.
  integer :: iostat, nv, nw
  integer, save :: itimer=0
  namelist /radiation_params/ stiff_lim, lock_Inu
  ! ----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%radiation_pre_init', itimer=itimer)
  self%kind = merge('angle_rad', 'mean_rad ', self%i_mu>0)
  if (first_time) then
    rewind (io%input)
    read (io%input, radiation_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, radiation_params)
    end if
    call io%namelist_warning ('radiation_params', iostat)
    write (io%namelist, radiation_params)
    first_time = .false.
  endif

  call self%ray_patch_t%pre_init()

  if (self%kind == 'angle_rad') then
    call self%ray%init(self%vec, self%wmu)
    call self%ray%set_stencil(self%ds * scaling%l)

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

  nv = self%fe-self%fs+1
  call self%set_nv_and_nw (nv, 1)
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE radiation_pre_init

!===============================================================================
SUBROUTINE radiation_init (self)
  class(rad_patch_t) :: self
  integer :: i, m(6)
  ! ----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%radiation_pre_init')
  call self%ray_patch_t%init()

  if (.not. self%out_ray) then
    m = shape(self%authority_mem)
    if (self%kind=='angle_rad') then
      call self%allocate_mem(self)
    end if
    if (lock_Inu) then
      m = shape(self%authority_mem)
      m(3) = 2
      call os%allocate(self%Inu_old, m(1:4), label=trim(self%kind)//'%Inu_old')
      self%Inu_old = tiny(1.0)
    endif
  endif

  call trace%end()
END SUBROUTINE radiation_init

!===============================================================================
SUBROUTINE radiation_update (self)
  class(rad_patch_t) :: self
  real(4), dimension(:,:,:,:), pointer :: Inu_new
  integer :: k, f, p(3)
  ! ----------------------------------------------------------------------------
  if (self%is_finished()) return
  !-----------------------------------------------------------------------------
  ! Boundary conditions
  !-----------------------------------------------------------------------------
  call self%bndries()
  call self%periodic_grid(only=merge(1, self%gn(3), self%dir>0))

  if (self%kind=='angle_rad') then
    call self%check_access(self%opac)
    if (self%out_ray) then
      call self%long_char_update()
    else
      call self%decompress_mem(self)
      call self%check_access(self%back_opac)
      if (self%ray%long) then
        print *, 'long ray detected! Stopping...'
        stop
        !call self%long_ray_update()
      else
        call self%short_ray_update()
        if (comoving_frame) then
          call self%integrate_intensity()
        endif
      endif
    endif
  else
    call self%Jnu_update()
  endif

  if (lock_Inu .and. .not. self%out_ray) then
    k = self%n(3)
    Inu_new => self%authority_mem(:,:,1:k:k-1,:,self%it,self%i_Inu)
    if (all(abs(1 - self%Inu_old/Inu_new) < stiff_lim)) then
      self%stiff_count = self%stiff_count*2
      self%stiff = .true.
      if (.not. self%opac%finished) then
        call self%opac%compress_mem(self)
        call self%opac%skip_iterations(self%stiff_count-1)
      endif

    elseif(self%stiff) then
      !print *, self%m3d_patch_info(), 'Inu has changed, but used to be stiff'
      self%stiff_count = 1
      self%stiff = .false.
    endif

    self%Inu_old = Inu_new
  endif

  call self%post_update()

  if (self%stiff .and. .not. self%is_finished()) then
    call self%skip_iterations(self%stiff_count-1)
  endif

  if (self%kind == 'angle_rad' .and. self%out_ray .and. self%boundaries%upper) then
    call self%remove_lower_layers_mem()
  endif
END SUBROUTINE radiation_update

!===============================================================================
SUBROUTINE Jnu_update (self)
  class(rad_patch_t) :: self
  class(nbor_t), pointer:: nbor
  integer :: i,j,k,f
  real(4), dimension(:,:,:,:), pointer :: Jnu, Psi, Ieff
  real(4), dimension(:,:,:,:), pointer :: Inu
  real(4), dimension(:,:,:,:,:), allocatable :: opac
  integer, save:: itimer=0
  real(4) :: wmu
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%Jnu_update', itimer=itimer)

  Jnu => self%mem(:,:,:,:, self%it, self%i_inu)
  Psi => self%Psi
  Ieff => self%Ieff
  Jnu = 0
  Psi = 0
  Ieff = 0

  nbor => self%link%nbor
  do while (associated(nbor))
    select type(nb => nbor%task)
    class is (rad_patch_t)
      Inu => nb%authority_mem(:,:,:,:, nb%it, self%i_inu)
      wmu = nb%wmu * over_4pi
      do f=1, self%nv
      do k=1, self%n(3)
      do j=1, self%n(2)
      do i=1, self%n(1)
        Jnu(i,j,k,f)  = Jnu(i,j,k,f)  + wmu * Inu(i,j,k,f)
        Psi(i,j,k,f)  = Psi(i,j,k,f)  + wmu * nb%Psi(i,j,k,f)
        Ieff(i,j,k,f) = Ieff(i,j,k,f) + wmu * nb%Ieff(i,j,k,f)
      end do
      end do
      end do
      end do
    end select
    nbor => nbor%next
  end do

  call trace%end(itimer)
END SUBROUTINE Jnu_update

!===============================================================================
SUBROUTINE short_ray_update (self)
  class(rad_patch_t), target :: self
  real(4), dimension(:,:), pointer, contiguous :: kernel
  real(4), dimension(:,:,:,:), pointer, contiguous :: eta_back
  real(4) :: dt, u0, ex, Inuup, xxup, ssup, xxmin, ssmin, Inumin, dsup05, Psi, xx, ss, Inu
  integer :: f, k, i, j, dx,dy,i1,j1,x1,x2,y1,y2, kup, k1, k2, dir, x, y, z, ii, jj
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%short_ray_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)
  eta_back => self%back_opac%mem(:,:,:,:, self%back_opac%it, self%back_opac%i_eta)

  call os%allocate (kernel, shape(self%ray%upkern))
  kernel = real(self%ray%upkern)
  dsup05 = 0.5*real(self%ray%dsup)

  dx = size(kernel, 1)
  dy = size(kernel, 2)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)

  do f=1, self%nv
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1
      kup = k-self%dir
      
      do j=self%li(2),self%ui(2)
        y = j - self%li(2) + 1
        y1 = j + self%ray%upo2 - 1
    
        !$omp simd
        do i=self%li(1),self%ui(1)
          x = i - self%li(1) + 1
          x1 = i + self%ray%upo1 - 1

          xxmin = huge(1.0)
          ssmin = huge(1.0)
          xxup = 0
          ssup = 0
          do jj=1,dy
            do ii=1,dx
              xx = self%chi(ii+x1, jj+y1, kup,f)
              ss = self%Snu(ii+x1, jj+y1, kup,f)
              xxup = xxup + kernel(ii,jj) * xx
              ssup = ssup + kernel(ii,jj) * ss
              if (ii>1 .and. ii<dx .and. jj>1 .and. jj<dy) then
                xxmin = min(xxmin, xx)
                ssmin = min(ssmin, ss)
              endif
            end do
          end do
          xxup = max(xxup, xxmin)
          ssup = max(ssup, ssmin)

          xx = self%chi(i,j,k,f)
          ss = self%Snu(i,j,k,f)
          dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
          if (dt < 80) then
            ex = exp(-dt) 
            if ( dt > 1.0e-3 ) then
              u0 = 1.0 - ex
              Psi = 1 - u0 / dt
            else
              Psi = dt * (0.5 - 0.16666667*dt)
              u0 = dt * (1 - Psi)
            endif

            Inuup = 0
            Inumin = huge(1.0)
            do jj=1,dy
              do ii=1,dx
                Inu = self%Inu(ii+x1, jj+y1, kup,f)
                Inuup = Inuup + kernel(ii,jj) * Inu
                if (ii>1 .and. ii <dx .and. jj>1 .and. jj<dy) then
                  Inumin = min(Inumin, Inu)
                endif
              end do
            end do
            Inuup = max(Inuup, Inumin)

            self%Inu(i,j,k,f) = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
          else
            !ex = 0
            u0 = 1 / dt
            Psi = 1 - u0
            self%Inu(i,j,k,f) = u0*ssup + Psi*ss
          endif

          if (.not. comoving_frame) then
            self%Psi(x,y,z,f) = Psi / xx
            self%Ieff(x,y,z,f) = self%Inu(i,j,k,f) - self%Psi(x,y,z,f) * (ss*xx - eta_back(x,y,z,f))
          endif
        end do
      end do

      call self%periodic_plane(self%Inu(:,:,k,f))
    end do ! k-loop
  end do ! f-loop

  call os%deallocate(kernel)
  call trace%end(itimer)
END SUBROUTINE short_ray_update

!===============================================================================
SUBROUTINE long_char_update (self)
  class(rad_patch_t), target :: self
  real(4), dimension(:,:,:), pointer, contiguous :: xx, ss
  real(4), dimension(:,:), contiguous, pointer :: kernel
  integer, dimension(:,:), contiguous, pointer :: ipart
  real(4), dimension(2) :: delta, deltak, fpart
  real(4) :: dsup, ex, dt, u0, Psi
  integer :: f, i, j, k
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%long_char_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)
  
  call self%bndries() ! reinforce boundary condition

  call os%allocate(kernel, [4,self%ui(3)])
  call os%allocate(ipart, [2,self%ui(3)])
  call os%allocate(xx, self%gn)
  call os%allocate(ss, self%gn)

  delta = self%vec(:2) / (self%vec(3) * self%ds(:2)) ! shift in units of grid cells
  dsup = real(self%ray%dsup)

  do k=1, self%ui(3)
    deltak = (self%box(3) - (k-1)*self%ds(3) - self%origin(3)) * delta
    ipart(:,k) = int(deltak) !integer part of shift
    fpart = deltak - ipart(:,k) !fractional part
    where (deltak < 0) 
      fpart = 1.0 + fpart
      ipart(:,k) = ipart(:,k) - 1.0
    end where

    ! bilinear interpolation to avoid overshooting
    kernel(1,k) = (1.0-fpart(1))*(1.0-fpart(2))
    kernel(2,k) =      fpart(1) *(1.0-fpart(2))
    kernel(3,k) = (1.0-fpart(1))*     fpart(2) 
    kernel(4,k) =      fpart(1) *     fpart(2)
  end do

  do f=1, self%nv
    ! ----------------------------------------------------------------------------
    ! frequency vectorisation is slower, because k varies quicker in mem than f.
    ! ----------------------------------------------------------------------------
    do k=1, self%ui(3)
      call self%interp_long(self%chi(:,:,k,f), xx(:,:,k), ipart(1,k), ipart(2,k), kernel(:,k))
      call self%interp_long(self%Snu(:,:,k,f), ss(:,:,k), ipart(1,k), ipart(2,k), kernel(:,k))
    end do

    do k=1, self%n(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      dt = 0.5 * dsup * ( xx(i,j,k) + xx(i,j,k+1) ) ! upwind optical pathlength
      if (dt < 80) then
        ex = exp(-dt) 
      else
        ex = 0
      endif
      ! Compute the U functions, use Taylor expansion for small dt
      ! to avoid accuracy issues with 1-(number slightly smaller than one)
      if ( dt > 1.0e-3 ) then
              u0 = 1.0 - ex
              Psi = 1 - u0 / dt
      else
              Psi = dt * (0.5 - 0.16666667*dt)
              u0 = dt * (1 - Psi)
      endif

      ! intensity at center point
      self%Inu(i,j,k+1,f) = self%Inu(i,j,k,f)*ex + (u0 - Psi) * ss(i,j,k) + Psi*ss(i,j,k+1)
      !Psi = Psi / xx(i,j,k+1)
    end do
    end do
    end do
  end do ! f-loop

  call os%deallocate(kernel)
  call os%deallocate(ipart)
  call os%deallocate(xx)
  call os%deallocate(ss)

  call trace%end(itimer)
END SUBROUTINE long_char_update

!===============================================================================
SUBROUTINE integrate_intensity(self)
  class(rad_patch_t), intent(in) :: self
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: nbtask
  class(rad_patch_t), pointer :: mean
  real(4), dimension(:,:,:), pointer :: vproj
  real(4), dimension(:,:,:,:), pointer :: Jnu, Inu
  real(8), dimension(:), pointer, contiguous :: nu_in, nu_out
  real(4) :: wmu
  integer :: f
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%integrate_intensity', itimer=itimer)
  nbor => self%link%nbor
  do while (associated(nbor))
    if (nbor%download) then
      nbtask => nbor%task
      select type (nbtask)
      class is (rad_patch_t)
        call nbtask%mem_lock(1)%set ('free_unused_memory')
        call nbtask%deallocate_mem(self, print=.true.)
        call nbtask%mem_lock(1)%unset ('free_unused_memory')
      end select
    endif
    nbor => nbor%next
  end do
  mean => self%mean_rad
  call mean%opac%mem_lock(1)%set('integrate_intensity')
  if (.not. associated(mean%mem)) then
    call mean%check_mem()
    do f=1, self%nv
      mean%mem(:,:,:,f,mean%it,1) = 0
    end do
    do f=1, self%nv
      mean%Psi(:,:,:,f) = 0
    end do
  endif
  call mean%mem_lock(1)%set('integrate_intensity')
  call mean%opac%mem_lock(1)%unset('integrate_intensity')

  call os%allocate(vproj, self%n)
  call self%comp_vproj(vproj)
  vproj = 1.0 + vproj
  wmu = self%wmu * over_4pi
  nu_in => mean%sp%nu(self%fs:self%fe)
  nu_out => mean%sp%nu(mean%fs:mean%fe)
  Jnu => mean%mem(:,:,:,:,mean%it,1)
  Inu => self%authority_mem(:,:,:,:,self%it,1)

  call interpolate_to_ray(vproj, nu_in, nu_out, Inu, Jnu, wmu)
  call interpolate_to_ray(vproj, nu_in, nu_out, self%Psi, mean%Psi, wmu)
  call mean%mem_lock(1)%unset('integrate_intensity')
  call trace%end(itimer)
END SUBROUTINE integrate_intensity

!===============================================================================
SUBROUTINE interp_long(self, Varin, Varup, up1, up2, kernel)
  class(rad_patch_t), intent(in) :: self
  integer :: x1,x2,y1,y2,i,j
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:), intent(in), contiguous :: VarIn
  real(4), dimension(:,:), intent(out), contiguous :: VarUp
  real(4), dimension(:), intent(in), contiguous :: kernel
  !-----------------------------------------------------------------------------
  do j=1,self%n(2)
    y1=modulo(j+up2-1,self%n(2))+1
    y2=modulo(j+up2  ,self%n(2))+1

    do i=1,self%n(1)
      x1=modulo(i+up1-1,self%n(1))+1
      x2=modulo(i+up1  ,self%n(1))+1

      ! indices might cross periodic boundary, hence no array sum possible
      VarUp(i,j) = kernel(1)*VarIn(x1,y1) + kernel(2)*Varin(x2,y1) &
                 + kernel(3)*VarIn(x1,y2) + kernel(4)*Varin(x2,y2)
    end do
  end do
END SUBROUTINE interp_long

!===============================================================================
SUBROUTINE periodic_plane (self, arr)
  class(rad_patch_t), intent(in) :: self
  real(4), dimension(:,:), intent(inout) :: arr
  integer, dimension(3) :: li, ui, lb, ub
  integer :: dir, is, ie, off
  !-----------------------------------------------------------------------------
  li = self%mesh%li; ui = self%mesh%ui
  lb = self%mesh%lb; ub = self%mesh%ub

  do dir=-1,1,2
    is = dir*max(dir*li(1), dir*ui(1))
    ie = dir*min(dir*li(1), dir*ui(1))

    do off=dir, dir*self%ng(1), dir
      arr(is+off,li(2):ui(2)) = arr(ie-dir+off,li(2):ui(2))
    end do
  end do

  do dir=-1,1,2
    is = dir*max(dir*li(2), dir*ui(2))
    ie = dir*min(dir*li(2), dir*ui(2))

    do off=dir, dir*self%ng(2), dir
      arr(:,is+off) = arr(:,ie-dir+off)
    end do
  end do
END SUBROUTINE periodic_plane

!===============================================================================
SUBROUTINE radiation_nbor_relations(self, another, needed, needs_me, download)
  class(rad_patch_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  ! ----------------------------------------------------------------------------
  call trace%begin ('rad_patch_t%nbor_relations')
  call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    select type(another)
    class is (ray_patch_t)
      if (self%kind==another%kind) then
        ! do nothing
      elseif(another%kind=='rate_patch' .and. self%i_mu==another%i_mu &
        .and. .not. comoving_frame) then
        needed=.false.
      elseif (self%kind=='angle_rad') then
          if (another%kind=='net_opac' .and. (self%out_ray .or. .not. comoving_frame)) then
            needs_me=.false.
            self%opac => another
          elseif (another%kind=='mean_rad'  .and. .not. self%out_ray) then
            needed=.false.
          elseif (another%kind=='back_opac' .and. .not. self%out_ray) then
            needs_me=.false.
            self%back_opac => another
          elseif (another%kind=='proj_patch') then
            needs_me=.false.
            self%opac => another
          elseif (another%kind=='out_rad') then
            needed=.false.
          else
            needed=.false. ; needs_me=.false.
          endif
      elseif (self%kind=='mean_rad') then
        if (another%kind=='angle_rad' .and. .not. another%out_ray) then
          needs_me=.false.
        elseif(another%kind=='back_opac') then
          needed=.false.
        else
          needed=.false. ; needs_me=.false.
        endif
      else
        needed=.false. ; needs_me=.false.
      endif
    class default
      needed=.false. ; needs_me=.false.
    end select
  endif
  call trace%end()
END SUBROUTINE radiation_nbor_relations

!===============================================================================
SUBROUTINE bndries (self)
  class(rad_patch_t):: self
  integer :: li, ui, f
  !----------------------------------------------------------------------------
  if (self%kind == 'angle_rad') then
    call self%check_access(self%opac)
    self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
    self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
    self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

    if (self%boundaries%upper) then
      ui = ubound(self%mem,3)
      do f=1, self%nv
        self%Inu(:,:,ui,f) = 0
      end do
    endif
  
    if (self%boundaries%lower) then
      li = lbound(self%mem,3)
      do f=1, self%nv
        self%Inu(:,:,li,f) = self%Snu(:,:,li,f)
      end do
    endif
  else
    call self%ray_patch_t%bndries()
  endif
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE remove_lower_layers_mem (self)
  class(rad_patch_t):: self
  real(4), dimension(:,:,:,:,:,:), contiguous, pointer :: top_slice 
  integer :: m(6)
  !----------------------------------------------------------------------------
  m = shape(self%authority_mem)
  m(3) = 1
  call os%allocate(top_slice, m, label='top_slice')
  top_slice = self%authority_mem(:,:,self%n(3):self%n(3),:,:,:)
  call os%deallocate(self%mem, label=trim(self%kind)//'%mem')
  self%mem => top_slice
  self%authority_mem => top_slice
  self%n(3)=1
END SUBROUTINE remove_lower_layers_mem

!===============================================================================
SUBROUTINE rad_patch_deallocate_mem (self, deallocator, print)
  class(rad_patch_t):: self
  class(m3d_patch_t), intent(in), target :: deallocator
  logical, optional :: print
  integer :: m(6), li(3), ui(3)
  !-----------------------------------------------------------------------------
  call self%ray_patch_t%deallocate_mem(deallocator, print)
  if (associated(self%Psi) .and. self%is_clear(bits%busy)) then
    call os%deallocate(self%Psi, label=trim(self%kind)//'%Psi')
    if (associated(self%Ieff)) then
      call os%deallocate(self%Ieff,label=trim(self%kind)//'%Ieff')
    endif
  endif
END SUBROUTINE rad_patch_deallocate_mem

!===============================================================================
SUBROUTINE rad_patch_allocate_mem (self, allocator, m_in)
  class(rad_patch_t):: self
  class(m3d_patch_t), intent(in), target :: allocator
  integer, intent(in), optional :: m_in(6)
  integer :: m(4)
  !-----------------------------------------------------------------------------
  call self%ray_patch_t%allocate_mem(allocator, m_in)
  if (.not. self%out_ray) then
    if (.not. (associated(self%Psi) .or. self%is_set(bits%virtual))) then
      m(1:3) = self%n
      m(4) = self%nv
      call os%allocate(self%Psi, m,label=trim(self%kind)//'%Psi')
      self%Psi  = tiny(1.0)
      if (.not. comoving_frame .or. self%i_mu == 0) then
        call os%allocate(self%Ieff,m,label=trim(self%kind)//'%Ieff')
        self%Ieff = tiny(1.0)
      endif
    endif
  endif
END SUBROUTINE rad_patch_allocate_mem

!===============================================================================
SUBROUTINE print_post_update_info (self)
  class(rad_patch_t) :: self
  character(48) :: text
  integer :: p(3), p2(3), f, w
  real(4) :: arr(3)
  !-----------------------------------------------------------------------------
  f = self%fcheck
  p = self%pcheck
  if (self%stiff) then
    write(text, '(A,I4)') 'intensity has not changed, skipping', self%stiff_count-1
    call self%print_to_stdout(text=text)
  else
    w = 1
    arr = 0
    arr(1) = self%mem(p(1),p(2),p(3),f,1,1)
    if (associated(self%Psi)) then
      w = w + 1
      p2 = p - self%ng
      if (all(p2>0)) then
        arr(2) = self%Psi(p2(1),p2(2),p2(3),f)
        if (associated(self%Ieff)) then
          w = w + 1
          arr(3) = self%Ieff(p2(1),p2(2),p2(3),f)
        endif
      endif
    endif
    call self%print_to_stdout(arr=arr(1:w))
  endif
END SUBROUTINE print_post_update_info

!===============================================================================
!> Fast method to generate nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(rad_patch_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  integer :: i_mu, i_nu, key(4)
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%init_nbors_hash', itimer=itimer)

  key(1) = self%key(1)
  select case (self%kind)
  case ('angle_rad')
    call self%ray_patch_t%init_nbors_hash(nbors)
    key(2:3) = self%key(2:3)
    if (comoving_frame .and. .not. self%out_ray) then
      key(4) = hash_string('proj_patch')
    else
      key(4) = hash_string('net_opac')
    endif
    call self%check_nbor(key, nbors)
    if (.not. self%out_ray) then
      key(4) = hash_string('rate_patch')
      call self%check_nbor(key, nbors)
      key(3) = 0
      key(4) = hash_string('back_opac')
      call self%check_nbor(key, nbors)
      key(4) = hash_string('mean_rad')
      call self%check_nbor(key, nbors)
    else
      key(2) = 0
      key(4) = hash_string('out_rad')
      call self%check_nbor(key, nbors)
    endif
  case ('mean_rad')
    key(2) = self%key(2)
    do i_mu=1, self%n_omega
      key(3) = i_mu
      key(4) = hash_string('angle_rad')
      call self%check_nbor(key, nbors)
    end do
  end select
  call trace%end (itimer)
END SUBROUTINE init_nbors_hash

!===============================================================================
!SUBROUTINE long_ray_update (self)
!  class(rad_patch_t), target :: self
!  !real(4), dimension(:,:,:,:), allocatable :: xxup, ssup
!  real(4), dimension(:,:,:,:,:), pointer :: xx_ss
!  real(4), dimension(:,:,:,:,:), allocatable :: xx_ss_up
!  real(4), dimension(:,:), allocatable :: Inuup, iil, kernel
!  real(4), dimension(:,:,:,:), pointer :: Inu, chi, Snu
!  integer :: k1, k2
!  real(4) :: dsup
!  integer :: f, k, z
!  integer, save:: itimer=0
!  !-----------------------------------------------------------------------------
!  call trace%begin('rad_patch_t%long_ray_update', itimer=itimer)
!
!  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
!  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
!  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)
!
!  Inu => self%authority_mem(:,:,:,:,self%it,self%i_inu)
!  chi => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)
!  Snu => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_Snu)
!
!  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
!  k2 = merge(self%li(3), self%ui(3), self%dir < 0)
!
!  !allocate(xxup, ssup, mold=chi)
!  allocate(Inuup(self%n(1), self%n(2)))
!  allocate(iil, mold=Inuup)
!  allocate(xx_ss_up(self%n(1), self%n(2), self%n(3), self%nv, self%opac%nw))
!  xx_ss => self%opac%mem(:,:,:,:,self%opac%it,:)
!  
!  call self%interp_xy_up(xx_ss, xx_ss_up, self%ray%longo1, self%ray%longo2, real(self%ray%longkern), .false.)
!
!  dsup = real(self%ray%dsup)
!  kernel = real(self%ray%upkern(2:3,2:3))
!  do f=1, self%nv
!  ! ----------------------------------------------------------------------------
!  ! frequency vectorisation is slower, because k varies quicker in mem than f.
!  ! ----------------------------------------------------------------------------
!    do k=k1,k2,self%dir
!      z = k - self%li(3) + 1
!
!      call self%interp_xy_up_plane(self%Inu(:,:,k-self%dir,f), Inuup, self%ray%upo1, self%ray%upo2, kernel, .false.)
!      call self%solve_long_plane(k, f, xx_ss_up(:,:,z,f,1), xx_ss_up(:,:,z,f,2), Inuup)
!      call linear_formsolv(dsup, xx_ss_up(:,:,z,f,1), chi(:,:,z,f), xx_ss_up(:,:,z,f,2), &
!                            Snu(:,:,z,f), Inuup, self%Psi(:,:,z,f), Inu(:,:,z,f))
!      call self%periodic_plane(self%Inu(:,:,k,f))
!    end do ! k-loop
!  end do ! f-loop
!
!  !self%Psi = self%Psi / chi
!
!  call trace%end(itimer)
!END SUBROUTINE long_ray_update

!===============================================================================
!pure subroutine interp_xy_up(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
!  class(rad_patch_t), intent(in) :: self
!  integer :: z1,k,f,v,k1
!  integer, intent(in) :: up1, up2
!  real(4), dimension(:,:,:,:,:), intent(in) :: VarIn
!  real(4), dimension(:,:,:,:,:), intent(out) :: VarUp
!  real(4), dimension(:,:), intent(in) :: kernel
!  logical, intent(in) :: check_undershoot
!  !-----------------------------------------------------------------------------
!  k1 = self%li(3) - 1 - self%dir
!
!  do v=1,2
!    do f=1,self%nv
!      do k=1,self%n(3)
!        z1 = k + k1
!        call self%interp_xy_up_plane(VarIn(:,:,z1,f,v), VarUp(:,:,k,f,v), up1, up2, kernel, check_undershoot)
!      end do
!    end do
!  end do
!end subroutine interp_xy_up

!===============================================================================
!pure subroutine interp_xy_up_plane(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
!  class(rad_patch_t), intent(in) :: self
!  integer :: x1,x2,y1,y2,i,j, dx,dy, i1,j1
!  integer, intent(in) :: up1, up2
!  real(4), dimension(:,:), intent(in) :: VarIn, kernel
!  real(4), dimension(:,:), intent(out) :: VarUp
!  logical, intent(in) :: check_undershoot
!  real(4) :: VarMin
!  !-----------------------------------------------------------------------------
!  dx = size(kernel, 1) - 1
!  dy = size(kernel, 2) - 1
!  i1 = self%li(1) - 1 + up1
!  j1 = self%li(2) - 1 + up2
!
!  do j=1,self%n(2)
!    y1 = j + j1
!    y2 = y1 + dy
!
!  !$omp simd
!  do i=1,self%n(1)
!      x1 = i + i1
!      x2 = x1 + dx
!
!      VarUp(i,j) = sum(kernel * VarIn(x1:x2, y1:y2))
!
!      if (check_undershoot) then
!        VarMin = minval(VarIn(x1+1:x1+2, y1+1:y1+2))
!        if (VarUp(i,j) < VarMin) then
!          VarUp(i,j) = VarMin
!        endif
!      endif
!  end do
!  end do
!end subroutine interp_xy_up_plane

!===============================================================================
!subroutine solve_long_plane(self, k, f, xxl, ssl, Inuup)
!  !---------------------------------------------------------------------
!  !
!  ! Perform a formal solution with linear interpolation of the source
!  ! function from the upwind XY crossing to the last cell intersection
!  ! before the center point
!  !
!  !---------------------------------------------------------------------
!  class(rad_patch_t), intent(in) :: self
!  real(4), dimension(:,:), intent(inout) :: xxl, ssl, Inuup
!  real(4), dimension(:,:), allocatable :: xx, ss
!  real(4), dimension(:,:), allocatable :: psi_tmp 
!  integer, intent(in) :: k, f
!  integer ::  x1,x2,y1,y2,z1,z2, i,j,m, i1,j1,k1
!  ! ----------------------------------------------------------------------------
!  allocate(xx, ss, psi_tmp, mold=Inuup)
!  ! linear interpolation of snu and chi at each xz or yz intersection
!  do m=1,self%ray%Nintersect
!    xx = 0.0
!    ss = 0.0
!    do j=1,2
!      do i=1,2 
!        x1 = self%ray%int(m)%iof(i,j) + self%li(1)
!        x2 = self%ray%int(m)%iof(i,j) + self%ui(1)
!        y1 = self%ray%int(m)%jof(i,j) + self%li(2)
!        y2 = self%ray%int(m)%jof(i,j) + self%ui(2)
!        z1 = k + self%ray%int(m)%kof(i,j)
!    
!        xx(:,:) = xx(:,:) + self%ray%int(m)%wght(i,j) * self%chi( x1:x2, y1:y2, z1, f)
!        ss(:,:) = ss(:,:) + self%ray%int(m)%wght(i,j) * self%Snu( x1:x2, y1:y2, z1, f)
!      end do
!    end do
!    call linear_formsolv(real(self%ray%int(m)%ds), xxl, xx, ssl, ss, Inuup, psi_tmp, Inuup)
!    xxl = xx
!    ssl = ss
!  end do
!end subroutine solve_long_plane

!====================================================================================
! Does the formal solution between the upwind and central point
!====================================================================================
!dir$ ATTRIBUTES FORCEINLINE :: linear_formsolv
!pure elemental subroutine linear_formsolv(ds, xx1, xx2, ss1, ss2, Inuup, Psi, Inu)
!  !!$omp declare simd(linear_formsolv) uniform(ds) 
!  real(4), intent(in) :: ds
!  real(4), intent(in) :: xx1, xx2, ss1, ss2, Inuup
!  real(4), intent(out) :: Psi
!  real(4), intent(out) :: Inu
!  real(4) :: dt, ex, u0
!  ! ----------------------------------------------------------------------------
!  dt = 0.5 * ds * ( xx1 + xx2 ) ! upwind optical pathlength
!  ex = exp(-dt) 
!
!  ! Compute the U functions, use Taylor expansion for small dt
!  ! to avoid accuracy issues with 1-(number slightly smaller than one)
!  if ( dt > 1.0e-3 ) then
!          u0 = 1.0 - ex
!          Psi = 1 - u0 / dt
!  else
!          Psi = dt * (0.5 - 0.16666667*dt)
!          u0 = dt * (1 - Psi)
!  endif
!
!  ! intensity at center point
!  Inu = Inuup*ex + (u0 - Psi)*ss1 + Psi*ss2
!end subroutine linear_formsolv

END MODULE
