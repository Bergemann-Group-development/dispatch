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
  USE geometry_var, only: geometry_type, ray_type
  USE mgeometry, only: set_stencil
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_t) :: radiation_t
    integer :: i_Inu=1
    class(opacity_t), pointer :: opac => null()
    class(radiation_t), pointer :: sheep(:) => null()
    class(radiation_t), pointer :: mean_rad => null()
    class(radiation_t), pointer :: out_rad => null()
    type(geometry_type) :: geometry
    type(ray_type) :: rays(1)
    type(ray_type) :: ray
    real(4), dimension(:,:,:,:), pointer :: Inu, chi, Snu
    real(4), dimension(:,:,:,:), allocatable :: Psi
    character(len=64) :: meta_file, ray_name
  contains
    procedure :: pre_init => radiation_pre_init
    procedure :: init => radiation_init
    procedure :: update => radiation_update
    procedure :: nbor_relations => radiation_nbor_relations
    procedure :: Inu_update
    procedure :: Jnu_update
    procedure :: bndries
    procedure :: set_ray_ang
    procedure :: interp_z_up_plane
    procedure :: interp_xy_up_plane
    procedure :: interp_xx_ss_up
    procedure :: interp_Inuup_plane
    procedure :: interp_xy_up_long_plane
    procedure :: solveLong_plane
    procedure :: linear_formsolv
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
    self%ng = [4,4,1]
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
  use parallel_mod, only: xs,ys,zs,xe,ye,ze
  class(radiation_t) :: self
  integer :: i
  ! ----------------------------------------------------------------------------
  call self%m3d_t%init()

  select case (self%kind)
  case ('angle_rad')
    allocate(self%Psi(self%n(1),  self%n(2),  self%n(3),  self%nv))
    allocate(self%Inu(self%gn(1), self%gn(2), self%gn(3), self%nv))
    allocate(self%chi, self%Snu, mold=self%Inu)

    self%geometry%n_angles = 1

    self%geometry%x = [( self%ds(1)*(i-0.5d0), i=1,self%gn(1) )] * scaling%l
    self%geometry%y = [( self%ds(2)*(i-0.5d0), i=1,self%gn(2) )] * scaling%l
    self%geometry%z = [( self%ds(3)*(i-0.5d0), i=1,self%gn(3) )] * scaling%l

    !allocate(self%ray(self%li(3):self%ui(3)))

    ! and fill angle info into the dispackage%ray struct
    self%ray%mux = self%vec(1)
    self%ray%muy = self%vec(2)
    self%ray%muz = self%vec(3)
    self%ray%wmu = self%wmu

    call self%set_ray_ang(self%ray)

    self%rays = [self%ray]
    ! compute short characteristic ray lengths and interpolation weights
    xs = self%li(1); ys = self%li(2); zs = self%li(3)
    xe = self%ui(1); ye = self%ui(2); ze = self%ui(3)
    call set_stencil(self%geometry, self%rays, self%dir)
    self%ray = self%rays(1)

    self%t_offset = 0.6d0
  case ('mean_rad', 'out_rad')
    self%t_offset = 0.8d0
  end select

  self%time = self%patch%time + self%t_offset*self%dtime

  if (self%kind == 'out_rad') then
    self%time = self%time + self%maxiter * self%dtime
  endif
END SUBROUTINE radiation_init

!===============================================================================
SUBROUTINE radiation_update (self)
  class(radiation_t) :: self
  ! ----------------------------------------------------------------------------
  select case (self%kind)
  case ('angle_rad')
    call self%Inu_update()
  case ('mean_rad')
    call self%Jnu_update()
  case ('out_rad')
    if (save_output) call self%save_radiation_output()
  end select

  self%iter = self%iter + 1
END SUBROUTINE radiation_update

!===============================================================================
SUBROUTINE Inu_update (self)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:), allocatable :: xxup, ssup
  real(4), dimension(:,:), allocatable :: Inuup
  real(4), dimension(:,:,:,:), pointer :: Inu, chi
  integer :: k1, k2, f, k, z, i, j
  ! ----------------------------------------------------------------------------
  call trace_begin('radiation_t%Inu_update')

  Inu => self%authority_mem(:,:,:,:,self%it,self%i_inu)
  chi => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)

  allocate(xxup, ssup, mold=chi)
  allocate(Inuup, mold=chi(:,:,1,1))

  if (self%dir == 1) then                       ! top to bottom of atmosphere, bottom to top in index space
    k1=self%li(3) ; k2=self%ui(3)
  else 
    k1=self%ui(3) ; k2=self%li(3)
  endif

  select case(self%ray%upplane)
  case(1,2) ! yz-plane or xz-plane
    print *, 'SolveLong_plane needs to be updated. Stopping...'
    stop
    !call self%interp_z_up_plane(k, xx, ss, self%ray%upplane)
    !call self%solveLong_plane(k, xx, ss, Inuup)
  case(3) ! xy-plane
    call self%interp_xx_ss_up(xxup,ssup)      
  end select

  do f=1, self%nv
  ! ----------------------------------------------------------------------------
  ! frequency vectorisation is slower, because k varies quicker in mem than f.
  ! ----------------------------------------------------------------------------
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1

      select case(self%ray%upplane)
          case(1,2) ! yz-plane or xz-plane
            !call self%interp_z_up_plane(k, xxup, ssup, self%ray%upplane)
            !call self%solveLong_plane(k, xxup, ssup, Inuup)
          case(3) ! xy-plane
            call self%interp_Inuup_plane(k, f, Inuup)         
      end select
      
      call self%linear_formsolv(xxup(:,:,z,f), self%chi(:,:,k,f), ssup(:,:,z,f), self%Snu(:,:,k,f), &
                                Inuup, self%Psi(:,:,z,f), Inu(:,:,z,f))

      ! HOT-FIX must be updated later
      ! Otherwise Inu remains zero in the guard zones
      ! Only valid for plane parallel atmospheres 
      self%Inu(:,:,k,f) = Inu(1,1,z,f)
    end do ! k-loop
  end do ! f-loop

  self%Psi = self%Psi / chi

  if (self%verbose > 1) then
    if (self%dir > 0) then
      write(*,'(A10, 3(I4), 3(ES18.7))'), trim(self%kind), self%patch%id, self%i_omega, self%i_nu, &
      self%mem(5,5,1,1,1,1), self%mem(5,5,self%ui(3),1,1,1)
    else
      write(*,'(A10, 3(I4), 3(ES18.7))'), trim(self%kind), self%patch%id, self%i_omega, self%i_nu, &
      self%mem(5,5,self%li(3),1,1,1), self%mem(5,5,self%gn(3),1,1,1)
    endif
  endif
  call trace%end()
END SUBROUTINE Inu_update

!===============================================================================
SUBROUTINE Jnu_update (self)
  class(radiation_t) :: self
  class(radiation_t), pointer :: sheep
  integer :: i
  real(4), dimension(:,:,:,:), pointer :: Jnu
  real(4), dimension(:,:,:,:,:), allocatable :: opac
  ! ----------------------------------------------------------------------------
  call trace_begin('radiation_t%Jnu_update')

  Jnu => self%mem(:,:,:,:, self%it, self%i_inu)
  Jnu = 0

  do i=1, size(self%sheep)
    Jnu = Jnu + self%sheep(i)%wmu * self%sheep(i)%authority_mem(:,:,:,:, self%sheep(i)%it, self%i_inu)
  end do

  Jnu = Jnu / math%pi4

  if (self%verbose > 1) then
    write(*,'(A10, 3(I4), 3(ES18.7))'), trim(self%kind), self%patch%id, self%i_omega, self%i_nu, Jnu(1,1,1,1), Jnu(1,1,self%n(3),1)
  endif
  call trace%end()
END SUBROUTINE Jnu_update

!===============================================================================
elemental subroutine set_ray_ang(self, ray)
  class(radiation_t), intent(in) :: self
  class(ray_type), intent(inout) :: ray
  ! ----------------------------------------------------------------------------

  if (ray%mux.gt.0) then
    ! x sweep in direction of increasing x index
    ray%iup = -1
    ray%idn =  1
    ray%x1 = self%li(1)
    ray%x2 = self%ui(1)
    ray%x3 = 1
  else
    ! x sweep in direction of decreasing x index
    ray%iup =  1
    ray%idn = -1
    ray%x1 = self%ui(1)
    ray%x2 = self%li(1)
    ray%x3 = -1
  end if

  if (ray%muy.gt.0) then
    ! y sweep in direction of increasing y index
    ray%jup = -1
    ray%jdn =  1
    ray%y1 = self%li(2)
    ray%y2 = self%ui(2)
    ray%y3 = 1
  else
    ! y sweep in direction of decreasing y index
    ray%jup = 1
    ray%jdn = -1
    ray%y1 = self%ui(2)
    ray%y2 = self%li(2)
    ray%y3 = -1
  end if

  if (ray%muz.gt.0) then
    ! z sweep in direction of increasing z index
    ray%kup = -1
    ray%kdn =  1
    ray%z1 = self%li(3)
    ray%z2 = self%ui(3)
    ray%z3 = 1
  else
    ! z sweep in direction of decreasing z index
    ray%kup =  1
    ray%kdn = -1
    ray%z1 = self%ui(3)
    ray%z2 = self%li(3)
    ray%z3 = -1
  end if
end subroutine set_ray_ang

!===============================================================================
pure subroutine interp_z_up_plane(self, k, xx, ss, axis)
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: xx,ss
  integer, intent(in) :: k, axis
  integer :: x1,x2,y1,y2,z1,i,j
  ! ----------------------------------------------------------------------------

  if (axis == 1) then
    x1 = self%li(1) + self%ray%iup
    x2 = self%ui(1) + self%ray%iup
  elseif (axis == 2) then
    y1 = self%li(2) + self%ray%jup
    y2 = self%ui(2) + self%ray%jup
  endif

  do i=2,3
  do j=2,3
      if (axis /= 1) then
        x1 = self%li(1) + self%ray%upo1 + i - 1
        x2 = self%ui(1) + self%ray%upo1 + i - 1     
      elseif (axis /= 2) then  
        y1 = self%li(2) + self%ray%upo1 + i - 1
        y2 = self%ui(2) + self%ray%upo1 + i - 1
      endif

      if (self%dir<0) then
        z1 = k + (j-2)
        !z1=4-j ! i.e. 2 then 1
      else
        z1 = k + (j-3)
        !z1=j-1 ! ie. 1 then 2
      endif

      xx(:,:,1,:) = xx(:,:,1,:) + self%ray%upkern(i,j) * self%chi( x1:x2, y1:y2, z1, :)
      ss(:,:,1,:) = ss(:,:,1,:) + self%ray%upkern(i,j) * self%Snu( x1:x2, y1:y2, z1, :)
  end do
  end do
end subroutine interp_z_up_plane

!===============================================================================
pure subroutine interp_xx_ss_up(self, xxup, ssup)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,z1,z2,i,j,k,f
  real(4),dimension(:,:,:,:),intent(out) :: xxup,ssup
  real(4),dimension(:,:,:,:), allocatable :: xxmin,ssmin
  ! ----------------------------------------------------------------------------
  xxup = 0.0
  ssup = 0.0

  do f=1,self%nv
  do k=1,self%n(3)
  do j=1,self%n(2)
  do i=1,self%n(1)

      x1 = i + self%li(1) - 1 + self%ray%upo1
      y1 = j + self%li(2) - 1 + self%ray%upo2
      z1 = k + self%li(3) - 1 - self%dir

      xxup(i,j,k,f) = xxup(i,j,k,f) + sum(self%ray%upkern * self%chi ( x1:x1+3, y1:y1+3, z1, f))
      ssup(i,j,k,f) = ssup(i,j,k,f) + sum(self%ray%upkern * self%Snu ( x1:x1+3, y1:y1+3, z1, f))

  end do
  end do
  end do
  end do

  ! suppress undershoot

  allocate(xxmin, ssmin, mold=xxup)
  
  z1 = self%li(3)-self%dir
  z2 = self%ui(3)-self%dir

  xxmin = huge(1.0e0)
  ssmin = huge(1.0e0)

  do i=1,2
  do j=1,2

      x1 = self%li(1) + self%ray%upo1 + i
      x2 = self%ui(1) + self%ray%upo1 + i               
      y1 = self%li(2) + self%ray%upo2 + j
      y2 = self%ui(2) + self%ray%upo2 + j 

      where ( xxmin > self%chi( x1:x2, y1:y2, z1:z2, :) )    xxmin = self%chi( x1:x2, y1:y2, z1:z2, :)
      where ( ssmin > self%Snu( x1:x2, y1:y2, z1:z2, :) )    ssmin = self%Snu( x1:x2, y1:y2, z1:z2, :)

  end do
  end do

  where ( xxup < xxmin ) xxup = xxmin
  where ( ssup < ssmin ) ssup = ssmin
end subroutine interp_xx_ss_up

!===============================================================================
pure subroutine interp_Inuup_plane(self, k, f, Inuup)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,z1,i,j
  integer, intent(in) :: k, f
  real(4),dimension(:,:), allocatable :: Inuupmin
  real(4),dimension(:,:), intent(out) :: Inuup
  ! ----------------------------------------------------------------------------
  z1 = k-self%dir ! 2 + ry%kup 

  Inuup = 0.0

  do j=1,self%n(2)
  do i=1,self%n(1)

      x1 = i + self%li(1) - 1 + self%ray%upo1
      y1 = j + self%li(2) - 1 + self%ray%upo2

      Inuup(i,j) = Inuup(i,j) + sum(self%ray%upkern * self%Inu ( x1:x1+3, y1:y1+3, z1, f))

  end do
  end do

  ! suppress undershoot

  allocate(Inuupmin, mold=Inuup)

  z1 = k-self%dir ! 2 + ry%kup 

  Inuupmin = huge(1.0e0)

  do i=1,2
  do j=1,2

      x1 = self%li(1) + self%ray%upo1 + i
      x2 = self%ui(1) + self%ray%upo1 + i               
      y1 = self%li(2) + self%ray%upo2 + j
      y2 = self%ui(2) + self%ray%upo2 + j 

      where ( Inuupmin > self%Inu( x1:x2, y1:y2, z1, f) ) Inuupmin = self%Inu( x1:x2, y1:y2, z1, f)

  end do
  end do

  where ( Inuup < Inuupmin ) Inuup = Inuupmin
end subroutine interp_Inuup_plane

!===============================================================================
pure subroutine interp_xy_up_plane(self, k, xx, ss, Inuup)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,z1,i,j
  integer, intent(in) :: k
  real(4),dimension(:,:,:,:),intent(inout) :: xx,ss
  real(4),dimension(:,:,:), allocatable :: xxmin,ssmin,Inuupmin
  real(4),dimension(:,:,:), intent(out) :: Inuup
  ! ----------------------------------------------------------------------------

  z1 = k-self%dir ! 2 + ry%kup 

  ! convolve with cubic convolution kernel
  do i=1,4
  do j=1,4

      x1 = self%li(1) + self%ray%upo1 + i - 1
      x2 = self%ui(1) + self%ray%upo1 + i - 1              
      y1 = self%li(2) + self%ray%upo2 + j - 1
      y2 = self%ui(2) + self%ray%upo2 + j - 1

      xx(:,:,1,:) = xx(:,:,1,:) + self%ray%upkern(i,j) * self%chi ( x1:x2, y1:y2, z1, :)
      ss(:,:,1,:) = ss(:,:,1,:) + self%ray%upkern(i,j) * self%Snu ( x1:x2, y1:y2, z1, :)
      Inuup       = Inuup       + self%ray%upkern(i,j) * self%Inu ( x1:x2, y1:y2, z1, :)

  end do
  end do

  ! suppress undershoot

  x1 = self%li(1) + self%ray%upo1 + 1
  x2 = self%ui(1) + self%ray%upo1 + 1               
  y1 = self%li(2) + self%ray%upo2 + 1
  y2 = self%ui(2) + self%ray%upo2 + 1 

  allocate(xxmin,ssmin,Inuupmin, mold=xx(:,:,1,:))

  xxmin    = huge(1.0e0)
  ssmin    = huge(1.0e0)
  Inuupmin = huge(1.0e0)

  do i=1,2
  do j=1,2

      x1 = self%li(1) + self%ray%upo1 + i
      x2 = self%ui(1) + self%ray%upo1 + i               
      y1 = self%li(2) + self%ray%upo2 + j
      y2 = self%ui(2) + self%ray%upo2 + j 

      where (    xxmin > self%chi( x1:x2, y1:y2, z1, :) )    xxmin = self%chi( x1:x2, y1:y2, z1, :)
      where (    ssmin > self%Snu( x1:x2, y1:y2, z1, :) )    ssmin = self%Snu( x1:x2, y1:y2, z1, :)
      where ( Inuupmin > self%Inu( x1:x2, y1:y2, z1, :) ) Inuupmin = self%Inu( x1:x2, y1:y2, z1, :)

  end do
  end do

  where ( xx(:,:,1,:) < xxmin    ) xx(:,:,1,:) = xxmin
  where ( ss(:,:,1,:) < ssmin    ) ss(:,:,1,:) = ssmin
  where ( Inuup       < Inuupmin ) Inuup       = Inuupmin
end subroutine interp_xy_up_plane

!===============================================================================
pure subroutine interp_xy_up_long_plane(self, k, xxl, ssl, iil)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,z1,i,j
  integer, intent(in) :: k
  real(4), dimension(:,:,:,:), intent(inout) :: xxl,ssl
  real(4), dimension(:,:,:,:), intent(inout) :: iil
  real(4), dimension(:,:,:), allocatable :: xxmin,ssmin,iimin
  ! ----------------------------------------------------------------------------

  xxl = 0.0
  ssl = 0.0
  iil = 0.0

  z1 = k-self%dir

  ! convolve with cubic convolution kernel
  do i=1,4
  do j=1,4

      x1 = self%li(1) + self%ray%longo1 + i - 1
      x2 = self%ui(1) + self%ray%longo1 + i - 1              
      y1 = self%li(2) + self%ray%longo2 + j - 1
      y2 = self%ui(2) + self%ray%longo2 + j - 1

      xxl(:,:,1,:) = xxl(:,:,1,:) + self%ray%longkern(i,j) * self%chi ( x1:x2, y1:y2, z1, :)
      ssl(:,:,1,:) = ssl(:,:,1,:) + self%ray%longkern(i,j) * self%Snu ( x1:x2, y1:y2, z1, :)
      iil(:,:,1,:) = iil(:,:,1,:) + self%ray%longkern(i,j) * self%Inu( x1:x2, y1:y2, z1, :)

  end do
  end do

  allocate(xxmin,ssmin,iimin, mold=xxl(:,:,1,:))

  ! suppress undershoot
  xxmin = huge(1.0e0)
  ssmin = huge(1.0e0)
  iimin = huge(1.0e0)

  do i=1,2
  do j=1,2

      x1 = self%li(1) + self%ray%longo1 + i
      x2 = self%ui(1) + self%ray%longo1 + i               
      y1 = self%li(2) + self%ray%longo2 + j
      y2 = self%ui(2) + self%ray%longo2 + j 

      where ( xxmin > self%chi ( x1:x2, y1:y2, z1, :) ) xxmin = self%chi ( x1:x2, y1:y2, z1, :)
      where ( ssmin > self%Snu ( x1:x2, y1:y2, z1, :) ) ssmin = self%Snu ( x1:x2, y1:y2, z1, :)
      where ( iimin > self%Inu( x1:x2, y1:y2, z1, :) ) iimin = self%Inu( x1:x2, y1:y2, z1, :)

  end do
  end do

  where ( xxl(:,:,1,:) < xxmin    ) xxl(:,:,1,:) = xxmin
  where ( ssl(:,:,1,:) < ssmin    ) ssl(:,:,1,:) = ssmin
  where ( iil(:,:,1,:) < iimin    ) iil(:,:,1,:) = iimin
end subroutine interp_xy_up_long_plane

!===============================================================================
subroutine solveLong_plane(self, k, xxl, ssl, Inuup)
  !---------------------------------------------------------------------
  !
  ! Perform a fomral solution with linear interpolation of the source
  ! function from the upwind XY crossing to the last cell intersection
  ! before the center point
  !
  !---------------------------------------------------------------------
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: xxl,ssl
  real(4), dimension(:,:,:,:), allocatable :: iil
  real(4), dimension(:,:,:), intent(inout) :: Inuup
  real(4), dimension(:,:,:), allocatable :: psi_tmp
  integer, intent(in) :: k
  integer ::  x1,x2,y1,y2,z1,i,j,m
  ! ----------------------------------------------------------------------------

  allocate(iil, mold=xxl)

  ! X,S and incoming intensity  at the upwind XY plane
  call self%interp_xy_up_long_plane(k, xxl,ssl,iil)

  ! linear interpolation of snu and chi at each xz or yz intersection
  do m=1,self%ray%Nintersect

      xxl(:,:,2,:) = 0.0
      ssl(:,:,2,:) = 0.0

    do i=1,2
        do j=1,2 

          x1=self%ray%int(m)%iof(i,j) + self%li(1)
          x2=self%ray%int(m)%iof(i,j) + self%ui(1)
          y1=self%ray%int(m)%jof(i,j) + self%li(2)
          y2=self%ray%int(m)%jof(i,j) + self%ui(2)

        if (self%dir<0) then
          z1=k-abs(self%ray%int(m)%kof(i,j))
          !z1=2-abs(self%ray%int(m)%kof(i,j))
        else
          z1=k-abs(self%ray%int(m)%kof(i,j))
          !z1=2-abs(self%ray%int(m)%kof(i,j))
        !  z1=1+abs(ry%int(m)%kof(i,j))                              
        endif
          xxl(:,:,2,:) = xxl(:,:,2,:) + self%ray%int(m)%wght(i,j) * self%chi( x1:x2, y1:y2, z1 ,:)
          ssl(:,:,2,:) = ssl(:,:,2,:) + self%ray%int(m)%wght(i,j) * self%Snu( x1:x2, y1:y2, z1 ,:)
        end do
    end do

    allocate(psi_tmp, mold=xxl(:,:,1,:))

    call self%linear_formsolv(xxl(:,:,1,:), xxl(:,:,2,:), ssl(:,:,1,:), ssl(:,:,2,:), &
                              iil(:,:,1,:), psi_tmp, iil(:,:,2,:))

    ! set new start point to old end point
    xxl(:,:,1,:) = xxl(:,:,2,:)
    ssl(:,:,1,:) = ssl(:,:,2,:)
    iil(:,:,1,:) = iil(:,:,2,:)
  end do

  Inuup = iil(:,:,2,:)
end subroutine solvelong_plane

!====================================================================================
! Does the formal solution between the upwind and central point
!====================================================================================
elemental subroutine linear_formsolv(self, xx1, xx2, ss1, ss2, Inuup, Psi, Inu)
  use mconst
  class(radiation_t), intent(in) :: self
  real(4), intent(in) :: xx1, xx2, ss1, ss2, Inuup
  real(4), intent(out) :: Psi
  real(4), intent(out) :: Inu
  real(8) :: dt, ex, u0
  ! ----------------------------------------------------------------------------

  dt = 0.5d0 * self%ray%dsup * ( xx1 + xx2 ) ! upwind optical pathlength
  ex = exp(-dt) 

  ! Compute the U functions, use Taylor expansion for small dt
  ! to avoid accuracy issues with 1-(number slightly smaller than one)
  if ( dt > 1.0d-3 ) then
          u0 = 1.0e0 - ex
          Psi = 1 - u0 / dt
  else
          Psi = dt * (half - one_sixth*dt)
          u0 = dt * (1 - Psi)
  endif

  ! intensity at center point
  Inu = Inuup*ex + (u0 - Psi)*ss1 + Psi*ss2
end subroutine linear_formsolv

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
        if (another%i_omega == self%i_omega) then
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

  filename = trim(out_folder) // trim(self%ray_name) // '_lam.bin'
  open(newunit=lunit, file=filename, access='stream')
  write(lunit) lam
  call self%save_meta4(filename, lam)
  close(lunit)

  filename = trim(out_folder) // trim(self%ray_name) // '_inu.bin'
  open(newunit=lunit, file=filename, access='stream')
  write(lunit) Inu
  call self%save_meta4(filename, Inu)
  close(lunit)
  call trace%end()
END SUBROUTINE save_radiation_output

!===============================================================================
SUBROUTINE setup_metafile(self)
  class(radiation_t),intent(inout):: self
  integer :: meta_unit, iostat
  ! ----------------------------------------------------------------------------

  write(self%ray_name,'(A1, I3.3, A4, I3.3)') 'p', self%patch%id, '_ray', self%i_omega/2

  self%meta_file = trim(out_folder) // trim(self%ray_name) // '_meta.txt'
  open(newunit=meta_unit, file=self%meta_file, form='formatted')

  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'position', self%patch%position
  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'vec', self%vec

  close(meta_unit)
END SUBROUTINE

!===============================================================================
SUBROUTINE save_meta4(self, filename, arr)
  class(radiation_t):: self
  character(len=*), intent(in) :: filename
  real(4), intent(in) :: arr(..)
  integer :: meta_unit
  ! ----------------------------------------------------------------------------

  open(newunit=meta_unit, file=self%meta_file, form='formatted', access='append')
  write(meta_unit, '(A64, I4, 8(I5))') filename, storage_size(arr), shape(arr)
  close(meta_unit)
END SUBROUTINE save_meta4

end module