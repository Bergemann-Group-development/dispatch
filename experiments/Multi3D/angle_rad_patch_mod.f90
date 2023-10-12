MODULE angle_rad_patch_mod
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
  USE mconst, only: one_third, one_sixth
  !USE scratch_mod
  USE handyfuncs_mod
  USE m3d_halo_mod, only:lock_mem, unlock_mem
  USE m3d_patch_mod, only: m3d_patch_t, comoving_frame, ilambd, decouple_continuum, &
  line_blanketing, use_vel, preconditioning, save_Qrad, binned_opac, make_opac_table
  USE ray_patch_mod, only: ray_patch_t, formal_solver
  USE proj_patch_mod, only: interpolate_to_ray
  USE ray_mod, only: one_dimensional
  USE os_mod
  USE omp_mod
  USE omp_lock_mod
  USE m3d_debug_mod
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(ray_patch_t) :: rad_patch_t
    class(ray_patch_t), pointer :: opac => null()
    class(ray_patch_t), pointer :: back_cont => null()
    class(rad_patch_t), pointer :: mean_rad => null()
    real(4), dimension(:,:,:,:), contiguous, pointer :: Inu=>null(), chi=>null(), Snu=>null()
    real(4), dimension(:,:,:,:), contiguous, pointer :: Psi=>null(), Ieff=>null()
    real(4), dimension(:,:,:), pointer, contiguous :: Qrad=>null()
    integer :: stiff_count=1
    logical :: stiff=.false.
    real(8) :: dsup=0
  contains
    procedure :: pre_init => radiation_pre_init
    procedure :: init => radiation_init
    procedure :: nbor_relations => radiation_nbor_relations
    !procedure :: long_ray_update
    procedure :: execute
    procedure :: long_char_update
    procedure :: short_ray_update
    procedure :: bndries
    procedure :: integrate_intensity
    !procedure :: interp_xy_up
    !procedure :: interp_xy_up_plane
    !procedure :: solve_long_plane
    procedure :: linear_solver
    procedure :: hermite_solver
    procedure :: remove_lower_layers_mem
    procedure :: deallocate_mem
    procedure :: allocate_mem
    procedure :: print_post_update_info
    procedure :: init_nbors_hash
    procedure :: save
  end type rad_patch_t

  !-----------------------------------------------------------------------------

  real(4), save :: stiff_lim=1e-8
  logical, public, save :: lock_Inu=.false., hermite=.false.

contains
  
!===============================================================================
SUBROUTINE radiation_pre_init (self)
  class(rad_patch_t), target:: self
  logical, save :: first_time=.true.
  integer :: iostat, nv, nw
  integer, save :: itimer=0
  namelist /radiation_params/ stiff_lim, lock_Inu, hermite
  ! ----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%radiation_pre_init', itimer=itimer)
  self%kind = 'angle_rad'
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
  call self%ray%init(self%vec, self%wmu)
  call self%ray%set_stencil(self%ds * scaling%l)

  nw = merge(3, 1, self%out_ray)
  nv = self%fe-self%fs+1
  call self%set_nv_and_nw (nv, nw)
  call self%set_guard_zones()
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
    call self%allocate_mem(self)
  endif

  call trace%end()
END SUBROUTINE radiation_init

!===============================================================================
SUBROUTINE execute (self)
  class(rad_patch_t) :: self
  real(4), dimension(:,:,:,:), pointer :: Inu_new
  integer :: k, f, p(3)
  integer, save :: itimer=0
  if (make_opac_table) return
  !-----------------------------------------------------------------------------
  call trace%begin (trim(self%kind)//'%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Boundary conditions
  !-----------------------------------------------------------------------------
  call lock_mem (self%opac, 'angle_rad_update', read=.true.)
  call self%check_access(self%opac)

  call self%bndries()
  call self%periodic_grid(only=merge(1, self%gn(3), self%dir>0))

  if (self%out_ray) then
    call self%long_char_update()
  else
    call self%decompress_mem(self)
    call lock_mem (self%back_cont, 'angle_rad_update', read=.true.)
    call self%check_access(self%back_cont)
    if (self%ray%long) then
      print *, 'long ray detected! Stopping...'
      stop
      !call self%long_ray_update()
    elseif (one_dimensional) then
      call self%long_char_update()
    else
      call self%short_ray_update()
    endif
    call unlock_mem (self%back_cont, 'angle_rad_update', read=.true.)
    if (comoving_frame .or. self%iter < ilambd .or. decouple_continuum .or. line_blanketing) then
      call self%integrate_intensity()
    endif
  endif
  call self%opac%compress_mem(self, static=.false., locked=.true.)
  call unlock_mem (self%opac, 'angle_rad_update', read=.true.)

  if (self%print_update_info) then
    f = self%fcheck
    p = self%pcheck
    memcheck(1) = self%mem(p(1),p(2),p(3),f,1,1)
  endif
  
  if (.not. self%out_ray .or. .not. self%save_output) then
    call self%compress_mem(self, static=.false., locked=.true.)
  endif
  
  call trace%end(itimer)
END SUBROUTINE execute

!===============================================================================
SUBROUTINE short_ray_update (self)
  class(rad_patch_t), target :: self
  real(4), dimension(:,:), pointer, contiguous :: kernel
  real(4), dimension(:,:,:,:), pointer, contiguous :: eta_back
  real(4) :: Inuup, xxup, ssup, xxmin, ssmin, Inumin, dsup05, xx, ss, Inu, u0, Psi, Ieff
  real(8) :: dt, ex ! important
  integer :: f, k, i, j, dx,dy,i1,j1,x1,x2,y1,y2, kup, k1, k2, dir, x, y, z, ii, jj, t, li(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%short_ray_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)
  eta_back => self%back_cont%mem(:,:,:,:, self%back_cont%it, self%back_cont%i_eta)

  !call self%init_periodic_idx()

  call os%allocate (kernel, shape(self%ray%upkern))
  kernel = real(self%ray%upkern)
  dsup05 = 0.5*real(self%ray%dsup)

  dx = size(kernel, 1)
  dy = size(kernel, 2)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)
  li = self%ng + 1

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
            if ( dt > 1d-3 ) then
              u0 = 1 - ex
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
            Ieff = Inuup*ex + (u0 - Psi)*ssup
          else
            !ex = 0
            u0 = 1 / dt
            Psi = 1 - u0
            Ieff = u0*ssup
          endif

          self%Inu(i,j,k,f) = Ieff + Psi*ss
          !self%Inu(i,j,k,f) = ss
          !Psi = 1.0
          !self%Ieff(x,y,z,f) = 0.

          if (comoving_frame) then
            self%Psi(x,y,z,f) = Psi
            !Ieff is computed in mean_rad%update after Psi has been integrated
          else
            self%Psi(x,y,z,f) = Psi / xx
            if (preconditioning) then
              self%Ieff(x,y,z,f) = Ieff + self%Psi(x,y,z,f) * eta_back(x,y,z,f)
            else
              self%Ieff(x,y,z,f) = self%Inu(i,j,k,f)
            endif
          endif
        end do
      end do

      do t=1, size(self%periodic_idx, 2)
        i = self%periodic_idx(1,t)
        j = self%periodic_idx(2,t)
        self%Inu(i,j,k,f) = self%Inu(self%periodic_i(i),self%periodic_j(j),k,f)
      end do
    end do ! k-loop
  end do ! f-loop

  call os%deallocate(kernel)
  call trace%end(itimer)
END SUBROUTINE short_ray_update

!===============================================================================
SUBROUTINE long_char_update (self)
  class(rad_patch_t), target :: self
  real(4), dimension(:,:,:,:), pointer, contiguous :: eta_back
  integer :: iw, iiw
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%long_char_update', itimer=itimer)

  do iw=1, self%opac%nw, 2
    self%chi => self%opac%mem(:,:,:,:, self%opac%it, iw)
    self%Snu => self%opac%mem(:,:,:,:, self%opac%it, iw+1)
    iiw = (iw+1)/2
    self%Inu => self%mem(:,:,:,:, self%it, iiw)
    if (.not. self%out_ray) then
      eta_back => self%back_cont%mem(:,:,:,:, self%back_cont%it, self%back_cont%i_eta)
    endif

    !if (hermite) then
    !  call self%hermite_solver(self%chi, self%Snu, self%Inu)
    !else
    call self%linear_solver(self%chi, self%Snu, self%Inu, self%Psi, self%Ieff, eta_back)
    !endif
  end do

  call trace%end(itimer)
END SUBROUTINE long_char_update

!===============================================================================
SUBROUTINE linear_solver(self, xx, ss, Inu, Psi, Ieff, eta_back)
  class(rad_patch_t), intent(in) :: self
  real(4), intent(in), dimension(:,:,:,:), contiguous :: xx, ss
  real(4), intent(inout), dimension(:,:,:,:), contiguous :: Inu
  real(4), intent(inout), dimension(:,:,:,:), pointer, contiguous, optional :: Psi, Ieff, eta_back
  real(4) :: dsup05, u0, Psi0, ss1, ss2, xx1, xx2
  real(8) :: ex, dt
  integer :: i,j,k,f, k1, k2, kup, z
  logical :: save_Psi
  !-----------------------------------------------------------------------------
  dsup05 = 0.5*real(self%ray%dsup)
  save_Psi = .not. self%out_ray .and. present(Psi) .and. present(Ieff) .and. present(eta_back)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)

  do f=1, self%nv
    do k=k1, k2, self%dir
      kup = k - self%dir
      z = k - self%ng(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      ss1 = ss(i,j,kup,f)
      ss2 = ss(i,j,k,f)
      xx1 = xx(i,j,kup,f)
      xx2 = xx(i,j,k,f)

      dt = dsup05 * ( xx1 + xx2 ) ! upwind optical pathlength
      if (dt < 80) then
        ex = exp(-dt)
      else
        ex = 0
      endif
      ! Compute the U functions, use Taylor expansion for small dt
      ! to avoid accuracy issues with 1-(number slightly smaller than one)
      if ( dt > 1d-3 ) then
        u0 = 1 - ex
        Psi0 = 1 - u0 / dt
      else
        Psi0 = dt * (0.5 - 0.16666667*dt)
        u0 = dt * (1 - Psi0)
      endif

      if (save_Psi) then
        Ieff(i,j,z,f) = Inu(i,j,kup,f)*ex + (u0 - Psi0)*ss1
        Inu (i,j,k,f) = Ieff(i,j,z,f) + Psi0*ss2
        if (comoving_frame) then
          Psi (i,j,z,f) = Psi0
          !Ieff is computed in mean_rad%update after Psi has been integrated
        else
          Psi (i,j,z,f) = Psi0 / xx2
          Ieff(i,j,z,f) = Ieff(i,j,z,f) + Psi (i,j,z,f) * eta_back(i,j,z,f)
        endif
      else
        ! intensity at center point
        Inu(i,j,k,f) = Inu(i,j,kup,f)*ex + (u0 - Psi0)*ss(i,j,kup,f) + Psi0*ss(i,j,k,f)
      endif
    end do
    end do
    end do
  end do ! f-loop
END SUBROUTINE linear_solver


!===============================================================================
SUBROUTINE hermite_solver(self, xx, ss, Inu)
  class(rad_patch_t), intent(in) :: self
  real(4), intent(in), dimension(:,:,:,:), pointer, contiguous :: xx, ss
  real(4), intent(inout), dimension(:,:,:,:), pointer, contiguous :: Inu
  real(4) :: dsup, u0, Psi, rdsup, dsup1, dsup2
  real(4) :: p0, pup, pdn
  real(4) :: alpha, beta, alphaprim, betaprim
  real(8) :: ex, dt, dtdn, t1, t2, t3, dt2, dt3, a
  integer :: i,j,k,f
  !-----------------------------------------------------------------------------
  dsup = real(self%ray%dsup)
  rdsup = 1./ dsup
  dsup1 = 0.5*dsup
  dsup2 = 1./12.*dsup*dsup
  do f=1, self%nv
    do k=1, self%n(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      pup = xx(i,j,k,f)
      p0  = xx(i,j,k+1,f)

      if (k < self%n(3)) then
        pdn = xx(i,j,k+2,f)

        ! one sided dchi/ds at upwind point
        t1 = (p0 - pup) * rdsup
        ! one sided dchi/ds at dnwind point
        t3 = (pdn - p0) * rdsup
        ! harmonic mean  dchi/ds at central point
        t2 = t1*t3
        if (t2 > 0) then
          t2 = t2 / (0.5d0*(t1 + t3))
        else
          t2 = 0
        endif

        ! upwind optical pathlength
        dt   = dsup1*(pup + p0) + dsup2*(t1-t2)
        dtdn = dsup1*(pdn + p0) + dsup2*(t2-t3)
      else
        dt = dsup1*(pup + p0)
      endif

      if (dt > 0.05) then
        ex=exp(-dt)
        dt2 = dt*dt
        dt3 = dt*dt2
        ! compute interpolation params
        alpha     =  ( 6.0*(dt-2.0)       + (-dt3+ 6.0*(dt+2.0))*ex  ) / dt3
        beta      =  ( (dt3-6.0*(dt-2.0)) - 6.0*(dt+2.0)*ex          ) / dt3
        alphaprim =  ( (2.0*dt-6.0)       + (dt2+4.0*dt+6.0)*ex      ) / dt2
        betaprim  =  ( (-dt2+4.0*dt-6.0)  + (2.0*dt+6.0)*ex          ) / dt2
      else
        dt2 = 0.5*dt*dt
        dt3 = one_third*dt*dt2
        ex = 1. - dt + dt2 - dt3
        alpha     =   0.8*dt3  - 0.7*dt2 + 0.5*dt
        beta      =   0.2*dt3  - 0.3*dt2 + 0.5*dt
        alphaprim = - 0.3*dt3  + one_sixth*dt2
        betaprim  =   0.2*dt3  - one_sixth*dt2
      endif

      pup = ss(i,j,k,f)
      p0  = ss(i,j,k+1,f)

      if (k < self%n(3)) then
        pdn = ss(i,j,k+2,f)

        ! one sided dS/dt at upwind point
        t1 = (p0 - pup) / dt
        t3 = (pdn - p0) / dtdn
        ! harmonic mean of dS/dt at central point
        t2 = t1*t3
        if (t2 > 0) then
          a  = one_third * (1 + dtdn/dt)
          t2 = t2 / ((1-a)*t1 + a*t3)
        else
          t2 = 0
        endif
      else
        t1 = (p0 - pup) / dt
        t2 = t1
      endif

      ! intensity at center point
      Inu(i,j,k+1,f) = Inu(i,j,k,f)*ex + alpha*pup + beta*p0 + alphaprim*t1 + betaprim*t2
    end do
    end do
    end do
  end do ! f-loop
END SUBROUTINE hermite_solver

!===============================================================================
SUBROUTINE integrate_intensity(self)
  class(rad_patch_t), intent(in) :: self
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: nbtask
  class(rad_patch_t), pointer :: mean
  real(4), dimension(:,:,:), pointer :: vproj
  real(4), dimension(:,:,:,:), pointer :: Jnu, Inu, Snu, chi
  real(8), dimension(:), pointer, contiguous :: nu_in, nu_out
  real(4) :: wmu, wmunu
  integer :: f, i, j, k, w
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('rad_patch_t%integrate_intensity', itimer=itimer)
  mean => self%mean_rad
  call lock_mem(mean, 'integrate_intensity')
  if (.not. associated(mean%mem)) then
    call mean%check_mem(locked=.true.)
    do w=1, mean%nw
    do f=1, mean%nv
      mean%mem(:,:,:,f,mean%it,w) = 0
    end do
    end do
  endif

  wmu = self%wmu
  Jnu => mean%mem(:,:,:,:,mean%it,1)
  Inu => self%authority_mem(:,:,:,:,self%it,1)

  if (comoving_frame .and. use_vel) then
    if (.not. mean%sp%sorted) call io%abort('Comoving frame for unsorted spectrum!')
    call os%allocate(vproj, self%n)
    call self%comp_vproj(vproj)
    vproj = 1.0 + vproj
    nu_in => mean%sp%nu(self%fs:self%fe)
    nu_out => mean%sp%nu(mean%fs:mean%fe)
    call interpolate_to_ray(vproj, nu_in, nu_out, Inu, Jnu, wmu)
    call interpolate_to_ray(vproj, nu_in, nu_out, self%Psi, mean%Psi, wmu)
    call interpolate_to_ray(vproj, nu_in, nu_out, self%Ieff, mean%Ieff, wmu)
    call os%deallocate(vproj)
  else
    do f=1, self%nv
    do k=1, self%n(3)
    do j=1, self%n(2)
    do i=1, self%n(1)
      Jnu(i,j,k,f)       = Jnu(i,j,k,f)       + wmu * Inu(i,j,k,f)
      mean%Psi(i,j,k,f)  = mean%Psi(i,j,k,f)  + wmu * self%Psi(i,j,k,f)
      mean%Ieff(i,j,k,f) = mean%Ieff(i,j,k,f) + wmu * self%Ieff(i,j,k,f)
    end do
    end do
    end do
    end do
  endif

  if (save_Qrad) then
    chi => self%opac%authority_mem(:,:,:,:,mean%it,self%i_chi)
    Snu => self%opac%authority_mem(:,:,:,:,mean%it,self%i_Snu)
    if (binned_opac) then
      wmunu = wmu
    endif
    do f=1, self%nv
      if (.not. binned_opac) then
        wmunu = wmu * mean%sp%wnu(mean%fs-f+1)
      endif
      do k=1, self%n(3)
      do j=1, self%n(2)
      do i=1, self%n(1)
        mean%Qrad(i,j,k) = mean%Qrad(i,j,k) + wmunu * chi(i,j,k,f) * (Inu(i,j,k,f) - Snu(i,j,k,f))
      end do
      end do
      end do
    end do
  endif

  call unlock_mem(mean, 'integrate_intensity')
  call trace%end(itimer)
END SUBROUTINE integrate_intensity

!===============================================================================
SUBROUTINE save (self, mem, label)
  class(rad_patch_t) :: self
  real(4), dimension(:,:,:,:,:) :: mem
  real(4), dimension(:,:,:,:,:), pointer :: mem2
  character(len=*) :: label
  integer, dimension(4) :: m
  integer :: fs, fe
  !-----------------------------------------------------------------------------
  call self%ray_patch_t%save(mem, label)
  if (self%padded) then
    fs = self%pre_pad_fs - self%fs + 1
    fe = fs + self%pre_pad_fe - self%pre_pad_fs
  else
    fs = 1
    fe = self%nv
  endif
  if (associated(self%Ieff)) then
    m = shape(self%Ieff)
    mem2(1:m(1),1:m(2),1:m(3),fs:fe,1:1) => self%Ieff(:,:,:,fs:fe)
    self%saved = .false.
    call self%ray_patch_t%save(mem2, label='Ieff')
  endif
  if (associated(self%Psi)) then
    m = shape(self%Psi)
    mem2(1:m(1),1:m(2),1:m(3),fs:fe,1:1) => self%Psi(:,:,:,fs:fe)
    self%saved = .false.
    call self%ray_patch_t%save(mem2, label='Psi')
  endif
  if (associated(self%Qrad)) then
    m(1:3) = shape(self%Qrad)
    mem2(1:m(1),1:m(2),1:m(3),1:1,1:1) => self%Qrad
    self%saved = .false.
    call self%ray_patch_t%save(mem2, label='Qrad')
  endif
END SUBROUTINE save

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
      else
        if (another%kind=='net_opac' .and. (self%out_ray .or. .not. comoving_frame)) then
          needs_me=.false.
          self%opac => another
        elseif (another%kind=='mean_rad'  .and. .not. self%out_ray) then
          needed=.false.
        elseif (another%kind=='back_cont' .and. .not. self%out_ray) then
          needs_me=.false.
          self%back_cont => another
        elseif (another%kind=='proj_patch') then
          needs_me=.false.
          self%opac => another
        elseif(another%kind=='atom_rate' .and. .not. comoving_frame) then
          needed=.false.
        elseif (another%kind=='out_rad') then
          needed=.false.
        else
          needed=.false. ; needs_me=.false.
        endif
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
  integer :: lb, ub, f, iw, iiw, li(3), ui(3), o1(2), o2(2)
  !----------------------------------------------------------------------------
  if (self%out_ray) then
    li = self%opac%li
    ui = self%opac%ui
  else
    li = self%opac%mesh%lb
    ui = self%opac%mesh%ub
  endif

  do iw=1, self%opac%nw, 2
    self%chi => self%opac%mem(:,:,:,:, self%opac%it, iw)
    self%Snu => self%opac%mem(:,:,:,:, self%opac%it, iw+1)
    iiw = (iw+1)/2
    self%Inu => self%mem(:,:,:,:, self%it, iiw)

    if (self%boundaries%is_set(bits%zu)) then
      ub = ubound(self%mem,3)
      do f=1, self%nv
        self%Inu(:,:,ub,f) = 0
      end do
    endif

    if (self%boundaries%is_set(bits%zl)) then
      lb = lbound(self%mem,3)
      do f=1, self%nv
        self%Inu(:,:,lb,f) = self%Snu(li(1):ui(1),li(2):ui(2),lb,f)
      end do
    endif
  end do
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE remove_lower_layers_mem (self)
  class(rad_patch_t):: self
  real(4), dimension(:,:,:,:,:,:), contiguous, pointer :: top_slice 
  integer :: m(6)
  !----------------------------------------------------------------------------
  m = shape(self%authority_mem)
  m(3) = 1
  call os%allocate(top_slice, m, label=trim(self%kind)//'%mem')
  top_slice = self%authority_mem(:,:,self%n(3):self%n(3),:,:,:)
  call os%deallocate(self%mem, label=trim(self%kind)//'%mem')
  self%mem => top_slice
  self%authority_mem => top_slice
  self%n(3)=1
END SUBROUTINE remove_lower_layers_mem

!===============================================================================
SUBROUTINE deallocate_mem (self, deallocator, print)
  class(rad_patch_t):: self
  class(m3d_patch_t), intent(in), target :: deallocator
  logical, optional :: print
  integer :: m(6), li(3), ui(3)
  !-----------------------------------------------------------------------------
  call self%ray_patch_t%deallocate_mem(deallocator, print)
  if (self%is_clear(bits%busy)) then
    if (associated(self%Psi)) then
      call os%deallocate(self%Psi, label=trim(self%kind)//'%Psi')
    endif
    if (associated(self%Ieff)) then
      call os%deallocate(self%Ieff,label=trim(self%kind)//'%Ieff')
    endif
  endif
END SUBROUTINE deallocate_mem

!===============================================================================
SUBROUTINE allocate_mem (self, allocator, m_in, label)
  class(rad_patch_t):: self
  class(m3d_patch_t), intent(in), target :: allocator
  integer, intent(in), optional :: m_in(6)
  character(len=*), optional :: label
  integer :: m(4)
  !-----------------------------------------------------------------------------
  call self%ray_patch_t%allocate_mem(allocator, m_in)
  if (.not. self%out_ray .and. .not. self%compressed) then
    if (.not. (associated(self%Psi) .or. self%is_set(bits%virtual))) then
      m(1:3) = self%n
      m(4) = self%nv
      call os%allocate(self%Psi, m,label=trim(self%kind)//'%Psi')
      self%Psi  = tiny(1.0)
      !if (.not. comoving_frame .or. self%i_mu == 0) then
        call os%allocate(self%Ieff,m,label=trim(self%kind)//'%Ieff')
        self%Ieff = tiny(1.0)
    endif
  endif
END SUBROUTINE allocate_mem

!===============================================================================
SUBROUTINE print_post_update_info (self)
  class(rad_patch_t) :: self
  character(48) :: text
  integer :: p(3), p2(3), f, w
  real(4) :: arr(3)
  !-----------------------------------------------------------------------------
  f = self%fcheck
  p = self%pcheck
  w = 1
  arr = 0
  arr(1) = memcheck(1)
  p2 = p - self%ng
  if (all(p2>0) .and. all(p2<=self%n)) then
    if (associated(self%Psi)) then
      w = w + 1
      arr(2) = self%Psi(p2(1),p2(2),p2(3),f)
    endif
    if (associated(self%Ieff)) then
      w = w + 1
      arr(3) = self%Ieff(p2(1),p2(2),p2(3),f)
    endif
  endif
  call self%print_to_stdout(arr=arr(1:w))
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
  call self%ray_patch_t%init_nbors_hash(nbors)
  key(1:3) = self%key(1:3)
  if (comoving_frame .and. .not. self%out_ray) then
    key(4) = hash_string('proj_patch')
  else
    key(4) = hash_string('net_opac')
  endif
  call self%check_nbor(key, nbors)
  if (.not. self%out_ray) then
    key(4) = hash_string('atom_rate')
    call self%check_nbor(key, nbors)
    key(3) = 0
    key(4) = hash_string('back_cont')
    call self%check_nbor(key, nbors)
    key(4) = hash_string('mean_rad')
    call self%check_nbor(key, nbors)
  else
    key(2) = 0
    key(4) = hash_string('out_rad')
    call self%check_nbor(key, nbors)
  endif
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
