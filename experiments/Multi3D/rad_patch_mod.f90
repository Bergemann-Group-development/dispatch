module rad_patch_mod
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
  USE math_mod
  USE bits_mod
  USE units_mod, only: cgs
  !USE scratch_mod
  USE m3d_patch_mod
  USE opac_patch_mod
  USE ray_mod, only:ray_type
  USE os_mod
  USE omp_mod
  USE omp_lock_mod
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_t) :: radiation_t
    integer :: i_Inu=1
    class(opacity_t), pointer :: opac => null()
    class(radiation_t), pointer :: sheep(:) => null()
    class(radiation_t), pointer :: mean_rad => null()
    class(radiation_t), pointer :: out_rad => null()
    real(4), dimension(:,:,:,:), contiguous, pointer :: Inu, chi, Snu
    real(4), dimension(:,:,:,:), contiguous, pointer :: Psi
    character(len=64) :: meta_file, ray_name
  contains
    procedure :: pre_init => radiation_pre_init
    procedure :: init => radiation_init
    procedure :: update => radiation_update
    procedure :: nbor_relations => radiation_nbor_relations
    procedure :: long_ray_update
    procedure :: long_char_update
    procedure :: short_ray_update
    procedure :: Jnu_update
    procedure :: bndries
    procedure :: periodic_index
    procedure :: periodic_plane
    procedure :: interp_long
    procedure :: interp_xy_up
    procedure :: interp_xy_up2
    procedure :: interp_xy_up3
    procedure :: interp_xy_up_plane
    procedure :: interp_xy_up_plane2
    procedure :: solve_long_plane
    procedure :: setup_metafile
    procedure :: save_meta4
    procedure :: save_radiation_output
  end type radiation_t

  !-----------------------------------------------------------------------------

  logical :: save_output=.false.
  character(len=64), save:: out_folder = './out/'

contains
  
!===============================================================================
SUBROUTINE radiation_pre_init (self)
  class(radiation_t), target:: self
  logical, save :: first_time=.true.
  integer :: iostat
  namelist /radiation_params/ save_output
  ! ----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, radiation_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, radiation_params)
    end if
    call io%namelist_warning ('radiation_params', iostat)
    first_time = .false.
  endif

  call self%m3d_t%pre_init()

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

  elseif (self%kind == 'out_rad') then
    self%n(3) = 1
    self%nv = self%sp%nnu
  endif

  self%nw = 1
  self%mv = self%nv
  !self%nv is assigned in experiment_mod : init_sub_patches
  self%nt = 1
  self%it = 1
  self%new = self%it
  
  ! ----------------------------------------------------------------------------
  ! Allocate and set unsigned and pervolume. Opacity, source function,
  ! intensity and populations are all positive definite
  ! Although quantities are generally per volume, but density is not stored
  ! locally, thus this flag must not be set.
  ! ----------------------------------------------------------------------------
  allocate (self%unsigned(self%nv), self%pervolume(self%nv))
  self%pervolume = .false.
  self%unsigned  = .true.

  if (save_output .and. self%kind=='out_rad') call self%setup_metafile()
END SUBROUTINE radiation_pre_init

!===============================================================================
SUBROUTINE radiation_init (self)
  class(radiation_t) :: self
  integer :: i
  ! ----------------------------------------------------------------------------
  call self%m3d_t%init()
  call self%lock%init ()

  if (self%kind=='angle_rad') then
    call os%allocate(self%Psi, [self%n(1),  self%n(2),  self%n(3),  self%nv])
  end if
END SUBROUTINE radiation_init

!===============================================================================
SUBROUTINE radiation_update (self)
  class(radiation_t) :: self
  ! ----------------------------------------------------------------------------
  call self%pre_update()
  if (self%finished) then
    return
  endif

  select case (self%kind)
  case ('angle_rad')
    if (self%opac%finished .and. self%opac%iter < self%opac%maxiter) then
      call self%finish(rotated=.true.)
      return
    endif

    if (self%out_ray) then
      call self%long_char_update()
    else
      if (self%ray%long) then
        call self%long_ray_update()
      else
        call self%short_ray_update()
      endif
    endif
  case ('mean_rad')
    call self%Jnu_update()
  case ('out_rad')
    if (save_output) call self%save_radiation_output()
  end select

  call self%post_update()

  if (self%kind == 'angle_rad') then
    call self%task_t%rotate
    self%rotated=.true.
  endif
END SUBROUTINE radiation_update

!===============================================================================
SUBROUTINE Jnu_update (self)
  class(radiation_t) :: self
  class(radiation_t), pointer :: sheep
  integer :: i
  real(4), dimension(:,:,:,:), pointer :: Jnu
  real(4), dimension(:,:,:,:,:), allocatable :: opac
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%Jnu_update', itimer=itimer)

  Jnu => self%mem(:,:,:,:, self%it, self%i_inu)
  Jnu = 0

  do i=1, size(self%sheep)
    Jnu = Jnu + self%sheep(i)%wmu * self%sheep(i)%authority_mem(:,:,:,:, self%sheep(i)%it, self%i_inu)
  end do

  Jnu = Jnu / math%pi4

  call trace%end(itimer)
END SUBROUTINE Jnu_update

!===============================================================================
SUBROUTINE long_ray_update (self)
  class(radiation_t), target :: self
  !real(4), dimension(:,:,:,:), allocatable :: xxup, ssup
  real(4), dimension(:,:,:,:,:), pointer :: xx_ss
  real(4), dimension(:,:,:,:,:), allocatable :: xx_ss_up
  real(4), dimension(:,:), allocatable :: Inuup, iil, kernel
  real(4), dimension(:,:,:,:), pointer :: Inu, chi, Snu
  integer :: k1, k2
  real(4) :: dsup
  integer :: f, k, z
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%long_ray_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

  Inu => self%authority_mem(:,:,:,:,self%it,self%i_inu)
  chi => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)
  Snu => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_Snu)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)

  !allocate(xxup, ssup, mold=chi)
  allocate(Inuup(self%n(1), self%n(2)))
  allocate(iil, mold=Inuup)
  allocate(xx_ss_up(self%n(1), self%n(2), self%n(3), self%nv, self%opac%nw))
  xx_ss => self%opac%mem(:,:,:,:,self%opac%it,:)
  
  call self%interp_xy_up(xx_ss, xx_ss_up, self%ray%longo1, self%ray%longo2, real(self%ray%longkern), .false.)

  dsup = real(self%ray%dsup)
  kernel = real(self%ray%upkern(2:3,2:3))
  do f=1, self%nv
  ! ----------------------------------------------------------------------------
  ! frequency vectorisation is slower, because k varies quicker in mem than f.
  ! ----------------------------------------------------------------------------
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1

      call self%interp_xy_up_plane(self%Inu(:,:,k-self%dir,f), Inuup, self%ray%upo1, self%ray%upo2, kernel, .false.)
      call self%solve_long_plane(k, f, xx_ss_up(:,:,z,f,1), xx_ss_up(:,:,z,f,2), Inuup)
      call linear_formsolv(dsup, xx_ss_up(:,:,z,f,1), chi(:,:,z,f), xx_ss_up(:,:,z,f,2), &
                            Snu(:,:,z,f), Inuup, self%Psi(:,:,z,f), Inu(:,:,z,f))
      call self%periodic_plane(self%Inu(:,:,k,f))
    end do ! k-loop
  end do ! f-loop

  self%Psi = self%Psi / chi

  call trace%end(itimer)
END SUBROUTINE long_ray_update

!===============================================================================
SUBROUTINE short_ray_update (self)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:),   pointer, contiguous :: ex, source, xxup, ssup
  real(4), dimension(:,:),       pointer, contiguous :: Inuup, kernel
  real(4), dimension(:,:),       pointer, contiguous :: II
  real(4), dimension(:,:,:,:),   pointer             :: Inu, xx,ss
  integer :: k1, k2, m(4), i2, i3, i4
  real(4) :: Inumin, dsup
  integer :: f, k, z, i, j, dx,dy,i1,j1,x1,x2,y1,y2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%short_ray_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

  Inu => self%authority_mem(:,:,:,:,self%it,self%i_inu)
  xx => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)
  ss => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_Snu)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)

  m = shape(ss)
  call os%allocate (Inuup, self%n(1:2))
  call os%allocate (kernel, shape(self%ray%upkern))
  call os%allocate (ex, m)
  call os%allocate (xxup, m)
  call os%allocate (ssup, m)
  call os%allocate (source, m)

  kernel = real(self%ray%upkern)
  dsup   = real(self%ray%dsup)
  call self%interp_xy_up3(self%chi, xxup, self%ray%upo1, self%ray%upo2, kernel, .true.)
  call self%interp_xy_up3(self%Snu, ssup, self%ray%upo1, self%ray%upo2, kernel, .true.)
  call extinction_and_source(dsup, xx, xxup, ss, ssup, self%Psi, ex, source)

  dx = size(kernel, 1) - 1
  dy = size(kernel, 2) - 1
  i1 = self%li(1) - 1 + self%ray%upo1
  j1 = self%li(2) - 1 + self%ray%upo2

  do f=1, self%nv
  ! ----------------------------------------------------------------------------
  ! frequency vectorisation is slower, because k varies quicker in mem than f.
  ! ----------------------------------------------------------------------------
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1
      II => self%Inu(:,:, k-self%dir,f)
      
      do j=1,self%n(2)
        y1 = j + j1
        y2 = y1 + dy
    
        !$omp simd
        do i=1,self%n(1)

          if (ex(i,j,z,f) == 0) then
            Inu(i,j,z,f) = source(i,j,z,f)
          else
            x1 = i + i1
            x2 = x1 + dx
            Inuup(i,j) = sum(kernel * II(x1:x2, y1:y2))

            !supress undershoot
            Inumin = minval(II(x1+1:x1+2, y1+1:y1+2))
            if (Inuup(i,j) < Inumin) then
              Inuup(i,j) = Inumin
            endif

            Inu(i,j,z,f) = Inuup(i,j)*ex(i,j,z,f) + source(i,j,z,f)
          endif
        end do
      end do

      call self%periodic_plane(self%Inu(:,:,k,f))
      !if (self%verbose>1 .and. z==self%pcheck(3) .and. f==self%fcheck) then
      !  write(*,'(F5.1, A12, 4(I4), 6(ES15.7))') self%time, 'Inuup', self%patch%id, self%i_nu, self%dir, self%i_mu, &
      !  Inuup(self%pcheck(1),self%pcheck(2))
      !endif
    end do ! k-loop
  end do ! f-loop

  !if (self%verbose>1) then
  !  write(*,'(F5.1, A12, 4(I4), 6(ES15.7))') self%time, 'xxup ssup', self%patch%id, self%i_nu, self%dir, self%i_mu, &
  !  xxup(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  ssup(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  ex(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  source(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck)
  !endif
  call os%deallocate(Inuup)
  call os%deallocate(kernel)
  call os%deallocate(ex)
  call os%deallocate(xxup)
  call os%deallocate(ssup)
  call os%deallocate(source)

  call trace%end(itimer)
END SUBROUTINE short_ray_update

!===============================================================================
SUBROUTINE long_char_update (self)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:), pointer, contiguous :: xx, ss, Psi, ex, source
  real(4), dimension(:,:,:,:), pointer :: xxup, ssup, xxl, ssl
  real(4), dimension(2) :: delta, deltak, fpart
  integer, dimension(2) :: ipart
  real(4), dimension(4) :: kernel
  integer :: f, k
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%long_char_update', itimer=itimer)

  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

  !xx => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)
  !ss => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_Snu)
  !xxup => self%chi(:,:,:self%n(3),:)
  !ssup => self%Snu(:,:,:self%n(3),:)
  call os%allocate(xx,     shape(self%mem))
  call os%allocate(ss,     shape(self%mem))
  call os%allocate(Psi,    shape(self%authority_mem))
  call os%allocate(ex,     shape(self%authority_mem))
  call os%allocate(source, shape(self%authority_mem))

  delta = self%vec(:2) / (self%vec(3) * self%ds(:2)) ! shift in units of grid cells

  do f=1, self%nv
    ! ----------------------------------------------------------------------------
    ! frequency vectorisation is slower, because k varies quicker in mem than f.
    ! ----------------------------------------------------------------------------
    do k=1, self%ui(3)
      deltak = (self%box(3) - (k-1)*self%ds(3) - self%origin(3)) * delta
      ipart = int(deltak) !integer part of shift
      fpart = deltak - ipart !fractional part
      where (deltak < 0) 
        fpart = 1.0 + fpart
        ipart = ipart - 1.0
      end where

      ! bilinear interpolation to avoid overshooting
      kernel(1) = (1.0-fpart(1))*(1.0-fpart(2))
      kernel(2) =      fpart(1) *(1.0-fpart(2))
      kernel(3) = (1.0-fpart(1))*     fpart(2) 
      kernel(4) =      fpart(1) *     fpart(2) 

      call self%interp_long(self%chi(:,:,k,f), xx(:,:,k,f), ipart(1), ipart(2), kernel)
      call self%interp_long(self%Snu(:,:,k,f), ss(:,:,k,f), ipart(1), ipart(2), kernel)
    end do
  end do ! f-loop

  xxl => xx(:,:,self%li(3):self%ui(3),:) 
  ssl => ss(:,:,self%li(3):self%ui(3),:) 
  xxup => xx(:,:,1:self%n(3),:) 
  ssup => ss(:,:,1:self%n(3),:) 
  
  call self%bndries() ! reinforce boundary condition

  call extinction_and_source(real(self%ray%dsup), xxl, xxup, ssl, ssup, Psi, ex, source)

  do f=1, self%nv
    do k=1, self%n(3)
      self%Inu(:,:,k+1,f) = self%Inu(:,:,k,f)*ex(:,:,k,f) + source(:,:,k,f)
      !if (self%verbose>1 .and. k==self%pcheck(3) .and. f==self%fcheck) then
      !  write(*,'(F5.1, A12, 4(I4), 6(ES15.7))') self%time, 'Inuup', self%patch%id, self%i_nu, self%dir, self%i_mu, &
      !  self%Inu(self%pcheck(1),self%pcheck(2), k, f)
      !endif
    end do
  end do ! f-loop

  !if (self%verbose>1) then
  !  write(*,'(F5.1, A12, 4(I4), 8(ES15.7))') self%time, 'xxup ssup', self%patch%id, self%i_nu, self%dir, self%i_mu, &
  !  self%chi(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  self%Snu(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  xxup(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  ssup(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  xxl(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), &
  !  ssl(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck), self%origin(3)
  !endif

  call os%deallocate(xx)
  call os%deallocate(ss)
  call os%deallocate(Psi)
  call os%deallocate(ex)
  call os%deallocate(source)

  call trace%end(itimer)
END SUBROUTINE long_char_update

!===============================================================================
pure subroutine interp_xy_up(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
  class(radiation_t), intent(in) :: self
  integer :: z1,k,f,v,k1
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:,:,:,:), intent(in) :: VarIn
  real(4), dimension(:,:,:,:,:), intent(out) :: VarUp
  real(4), dimension(:,:), intent(in) :: kernel
  logical, intent(in) :: check_undershoot
  !-----------------------------------------------------------------------------
  k1 = self%li(3) - 1 - self%dir

  do v=1,2
    do f=1,self%nv
      do k=1,self%n(3)
        z1 = k + k1
        call self%interp_xy_up_plane(VarIn(:,:,z1,f,v), VarUp(:,:,k,f,v), up1, up2, kernel, check_undershoot)
      end do
    end do
  end do
end subroutine interp_xy_up

!===============================================================================
subroutine interp_xy_up2(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
  class(radiation_t), intent(in) :: self
  integer :: z1,k,f,v,k1
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:,:,:), intent(in) :: VarIn
  real(4), dimension(:,:,:,:), intent(out) :: VarUp
  real(4), dimension(:,:), intent(in) :: kernel
  logical, intent(in) :: check_undershoot
  !-----------------------------------------------------------------------------
  k1 = self%li(3) - 1 - self%dir

  do f=1,self%nv
    do k=1,self%n(3)
      z1 = k + k1
      call self%interp_xy_up_plane(VarIn(:,:,z1,f), VarUp(:,:,k,f), up1, up2, kernel, check_undershoot)
    end do
  end do
end subroutine interp_xy_up2

!===============================================================================
subroutine interp_xy_up3(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
  class(radiation_t):: self
  integer :: x1,x2,y1,y2,z1,i,j,k,f,v,i1,j1,k1,dx,dy
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:,:,:) :: VarIn
  real(4), dimension(:,:,:,:) :: VarUp
  real(4), dimension(:,:) :: kernel
  logical, intent(in) :: check_undershoot
  real(4) :: VarMin
  !-----------------------------------------------------------------------------
  dx = size(kernel, 1) - 1
  dy = size(kernel, 2) - 1
  i1 = self%li(1) - 1 + up1
  j1 = self%li(2) - 1 + up2
  k1 = self%li(3) - 1 - self%dir

  do f=1,self%nv
    do k=1,self%n(3)
      z1 = k + k1

      do j=1,self%n(2)
        y1 = j + j1
        y2 = y1 + dy

        !$omp simd
        do i=1,self%n(1)
            x1 = i + i1
            x2 = x1 + dx

            VarUp(i,j,k,f) = sum(kernel * VarIn(x1:x2, y1:y2, z1, f))

            if (check_undershoot) then
              VarMin = minval(VarIn(x1+1:x1+2, y1+1:y1+2, z1, f))
              if (VarUp(i,j,k,f) < VarMin) then
                VarUp(i,j,k,f) = VarMin
              endif
            endif
        end do
      end do
    end do
  end do
end subroutine interp_xy_up3

!===============================================================================
pure subroutine interp_xy_up_plane(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,i,j, dx,dy, i1,j1
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:), intent(in) :: VarIn, kernel
  real(4), dimension(:,:), intent(out) :: VarUp
  logical, intent(in) :: check_undershoot
  real(4) :: VarMin
  !-----------------------------------------------------------------------------
  dx = size(kernel, 1) - 1
  dy = size(kernel, 2) - 1
  i1 = self%li(1) - 1 + up1
  j1 = self%li(2) - 1 + up2

  do j=1,self%n(2)
    y1 = j + j1
    y2 = y1 + dy

  !$omp simd
  do i=1,self%n(1)
      x1 = i + i1
      x2 = x1 + dx

      VarUp(i,j) = sum(kernel * VarIn(x1:x2, y1:y2))

      if (check_undershoot) then
        VarMin = minval(VarIn(x1+1:x1+2, y1+1:y1+2))
        if (VarUp(i,j) < VarMin) then
          VarUp(i,j) = VarMin
        endif
      endif
  end do
  end do
end subroutine interp_xy_up_plane

!===============================================================================
pure subroutine interp_xy_up_plane2(self, VarIn, VarUp, up1, up2, kernel, check_undershoot)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,i,j, kx,ky
  integer, intent(in) :: up1, up2
  real(4), dimension(:,:), intent(in) :: VarIn, kernel
  real(4), dimension(:,:), intent(out) :: VarUp
  logical, intent(in) :: check_undershoot
  real(4), dimension(:,:), allocatable :: VarMin
  !-----------------------------------------------------------------------------
  kx = size(kernel, 1)
  ky = size(kernel, 2)

  VarUp = 0.0

  do i=1,kx
    do j=1,ky

       x1 = self%li(1) + up1 + i - 1
       x2 = self%ui(1) + up1 + i - 1              
       y1 = self%li(2) + up2 + j - 1
       y2 = self%ui(2) + up2 + j - 1

       VarUp = VarUp + kernel(i,j) * VarIn(x1:x2, y1:y2)

    end do
  end do

  if (check_undershoot) then
    ! suppress undershoot
    allocate(Varmin, mold=VarUp)
    Varmin = huge(1.0e0)

    do i=1,2
        do j=1,2

          x1 = self%li(1) + up1 + i
          x2 = self%ui(1) + up1 + i               
          y1 = self%li(2) + up2 + j
          y2 = self%ui(2) + up2 + j 

          where ( VarMin > VarIn(x1:x2, y1:y2) ) VarMin = VarIn(x1:x2, y1:y2)
        end do
    end do

    where ( VarUp < VarMin ) VarUp = VarMin
  endif
end subroutine interp_xy_up_plane2

!===============================================================================
subroutine solve_long_plane(self, k, f, xxl, ssl, Inuup)
  !---------------------------------------------------------------------
  !
  ! Perform a formal solution with linear interpolation of the source
  ! function from the upwind XY crossing to the last cell intersection
  ! before the center point
  !
  !---------------------------------------------------------------------
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:), intent(inout) :: xxl, ssl, Inuup
  real(4), dimension(:,:), allocatable :: xx, ss
  real(4), dimension(:,:), allocatable :: psi_tmp 
  integer, intent(in) :: k, f
  integer ::  x1,x2,y1,y2,z1,z2, i,j,m, i1,j1,k1
  ! ----------------------------------------------------------------------------
  allocate(xx, ss, psi_tmp, mold=Inuup)
  ! linear interpolation of snu and chi at each xz or yz intersection
  do m=1,self%ray%Nintersect
    xx = 0.0
    ss = 0.0
    do j=1,2
      do i=1,2 
        x1 = self%ray%int(m)%iof(i,j) + self%li(1)
        x2 = self%ray%int(m)%iof(i,j) + self%ui(1)
        y1 = self%ray%int(m)%jof(i,j) + self%li(2)
        y2 = self%ray%int(m)%jof(i,j) + self%ui(2)
        z1 = k + self%ray%int(m)%kof(i,j)
    
        xx(:,:) = xx(:,:) + self%ray%int(m)%wght(i,j) * self%chi( x1:x2, y1:y2, z1, f)
        ss(:,:) = ss(:,:) + self%ray%int(m)%wght(i,j) * self%Snu( x1:x2, y1:y2, z1, f)
      end do
    end do
    call linear_formsolv(real(self%ray%int(m)%ds), xxl, xx, ssl, ss, Inuup, psi_tmp, Inuup)
    xxl = xx
    ssl = ss
  end do
end subroutine solve_long_plane

!===============================================================================
SUBROUTINE interp_long(self, Varin, Varup, up1, up2, kernel)
  class(radiation_t), intent(in) :: self
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

!====================================================================================
! Does the formal solution between the upwind and central point
!====================================================================================
pure elemental subroutine linear_formsolv(ds, xx1, xx2, ss1, ss2, Inuup, Psi, Inu)
  !!$omp declare simd(linear_formsolv) uniform(ds) 
  real(4), intent(in) :: ds
  real(4), intent(in) :: xx1, xx2, ss1, ss2, Inuup
  real(4), intent(out) :: Psi
  real(4), intent(out) :: Inu
  real(4) :: dt, ex, u0
  ! ----------------------------------------------------------------------------

  dt = 0.5 * ds * ( xx1 + xx2 ) ! upwind optical pathlength
  ex = exp(-dt) 

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
  Inu = Inuup*ex + (u0 - Psi)*ss1 + Psi*ss2
end subroutine linear_formsolv

!===============================================================================
ELEMENTAL SUBROUTINE extinction_and_source(dsup, xx, xxup, ss, ssup, Psi, ex, source)
  real(4), intent(in) :: dsup
  real(4), intent(in) :: xxup, ssup, xx, ss
  real(4), intent(out) :: Psi, ex, source
  real(4) :: u0, dt
  ! ----------------------------------------------------------------------------
  dt = 0.5 * dsup * ( xxup + xx ) ! upwind optical pathlength
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

  source = (u0 - Psi) * ssup + Psi*ss
  Psi = Psi / xx
END SUBROUTINE extinction_and_source

!===============================================================================
SUBROUTINE periodic_plane (self, arr)
  class(radiation_t), intent(in) :: self
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
  class(radiation_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  ! ----------------------------------------------------------------------------
  
  call trace%begin ('radiation_t%nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  ! First check the main relations, based on location and level of refinement
  needed = .false.; needs_me = .false.; download = .false.

  select case (self%kind)
  case ('angle_rad')
    if (another%id == self%opac%id) then
      needed = .true.; needs_me = .false.; download = .false.
    elseif(another%id == self%mean_rad%id) then
      needed = .false.; needs_me = .true.; download = .false.
    elseif (associated(self%out_rad)) then
      if (another%id == self%out_rad%id) then
        needed = .false.; needs_me = .true.; download = .false.
      endif
    endif
    select type (another)
    class is (radiation_t)
      if (another%i_nu == self%i_nu) then
        if (another%i_mu == self%i_mu) then
          if (another%kind == self%kind) then
            if (self%is_upstream(another)) then
              needed = .true.; needs_me = .false.; download = .true.
            endif
          endif
        endif
      endif
    end select
  case ('mean_rad')
    if (any(self%sheep%id == another%id)) then
      needed = .true.; needs_me = .false.; download = .false.
    endif
  case ('out_rad')
    if (any(self%sheep%id == another%id)) then
      needed = .true.; needs_me = .false.; download = .false.
    endif
  end select

  call trace%end()
END SUBROUTINE radiation_nbor_relations

!===============================================================================
SUBROUTINE bndries (self)
  class(radiation_t):: self
  integer :: li(3), ui(3)
  !----------------------------------------------------------------------------
  if (self%kind == 'angle_rad') then
    self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
    self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
    self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

    if (self%boundaries%topbot(1)) then
      self%Inu(:,:,self%gn(3),:) = 0
    endif
  
    if (self%boundaries%topbot(2)) then
      self%Inu(:,:,1,:) = self%Snu(:,:,1,:)
    endif
  else
    call self%m3d_t%bndries()
  endif
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE save_radiation_output(self)
  USE mconst, only: cc
  class(radiation_t),intent(inout):: self
  class(radiation_t), pointer :: sheep
  character(len=64) :: filename
  real(4), allocatable :: lam(:)
  real(4), dimension(:,:,:), pointer :: Inu
  integer :: lunit
  integer :: i
  !-----------------------------------------------------------------------------
  call trace_begin('radiation_t%save_radiation_output')
  write(*,'(A64)') 'Writing output for angle ' // self%ray_name

  do i=1, size(self%sheep)
    sheep => self%sheep(i)
    Inu => self%mem(:,:,1, sheep%fs:sheep%fe, self%it, self%i_inu)
    Inu = sheep%authority_mem(:,:,sheep%n(3), :, sheep%it, sheep%i_inu)
  end do

  lam = 1e8 * cc / self%sp%nu
  Inu => self%mem(:,:,1,:, self%it, self%i_inu)

  if (self%i_out == 1) then
    filename = trim(out_folder) // 'lam.bin'
    open(newunit=lunit, file=filename, access='stream')
    write(lunit) lam
    !call self%save_meta4(filename, pack(lam, .true.), shape(lam))
    close(lunit)
  endif

  filename = trim(self%ray_name) // '_inu.bin'
  call self%save_meta4(filename, pack(Inu, .true.), shape(Inu))
  filename = trim(out_folder) // trim(self%ray_name) // '_inu.bin'
  open(newunit=lunit, file=filename, access='stream')
  write(lunit) Inu
  close(lunit)
  
  call trace%end()
END SUBROUTINE save_radiation_output

!===============================================================================
SUBROUTINE setup_metafile(self)
  class(radiation_t),intent(inout):: self
  integer :: meta_unit, iostat
  ! ----------------------------------------------------------------------------

  write(self%ray_name,'(A3, I3.3)') 'ray', self%i_out

  self%meta_file = trim(out_folder) // trim(self%ray_name) // '_meta.txt'
  open(newunit=meta_unit, file=self%meta_file, form='formatted')

  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'vec', self%vec
  write(meta_unit, '(A8, 1X, SP, F12.8)') 'weight', self%wmu

  close(meta_unit)
END SUBROUTINE

!===============================================================================
SUBROUTINE save_meta4(self, filename, arr, shpe)
  class(radiation_t):: self
  character(len=*), intent(in) :: filename
  real(4), dimension(:), intent(in) :: arr
  integer, dimension(:), intent(in) :: shpe
  integer :: meta_unit
  ! ----------------------------------------------------------------------------

  open(newunit=meta_unit, file=self%meta_file, form='formatted', position='append')
  write(meta_unit, '(A64, I4, 8(I5))') filename, storage_size(arr), shpe
  close(meta_unit)
END SUBROUTINE save_meta4

!===============================================================================
integer pure FUNCTION periodic_index(self, i, axis)
  class(radiation_t), intent(in) :: self
  integer, intent(in) :: i, axis
  ! ----------------------------------------------------------------------------
  if (self%periodic_xy) then
    periodic_index = mod(i, self%n(axis))
    if (periodic_index == 0) periodic_index = self%n(axis)
  else
    periodic_index = i
  endif
END FUNCTION periodic_index

end module
