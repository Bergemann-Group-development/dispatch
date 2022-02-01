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
  USE geometry_mod, only:ray_type
  implicit none

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_t) :: radiation_t
    integer :: i_Inu=1
    class(opacity_t), pointer :: opac => null()
    class(radiation_t), pointer :: sheep(:) => null()
    class(radiation_t), pointer :: mean_rad => null()
    class(radiation_t), pointer :: out_rad => null()
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
    procedure :: long_ray_update
    procedure :: short_ray_update
    procedure :: Jnu_update
    procedure :: bndries
    procedure :: periodic_plane
    procedure :: interp_xy_up_plane
    procedure :: interp_xx_ss_up
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

    if (self%ray%long) then
      self%ng = [3,3,1]
    else
      self%ng = [2,2,1]
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

  select case (self%kind)
  case ('angle_rad')
    allocate(self%Psi(self%n(1),  self%n(2),  self%n(3),  self%nv))
    allocate(self%Inu(self%gn(1), self%gn(2), self%gn(3), self%nv))
    allocate(self%chi, self%Snu, mold=self%Inu)

    self%t_offset = 0.6d0
  case ('mean_rad')
    self%t_offset = 0.8d0
  case ('out_rad')
    self%iter = self%maxiter
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

  call self%post_update()
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
SUBROUTINE Inu_update (self)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:), pointer :: Inu, chi, Snu
  integer :: k1, k2
  !-----------------------------------------------------------------------------
  self%Inu => self%mem(:,:,:,:, self%it, self%i_inu)
  self%chi => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_chi)
  self%Snu => self%opac%mem(:,:,:,:, self%opac%it, self%opac%i_snu)

  Inu => self%authority_mem(:,:,:,:,self%it,self%i_inu)
  chi => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_chi)
  Snu => self%opac%authority_mem(:,:,:,:,self%opac%it,self%opac%i_Snu)

  k1 = merge(self%li(3), self%ui(3), self%dir > 0)
  k2 = merge(self%li(3), self%ui(3), self%dir < 0)

  if (self%ray%long) then
    call self%long_ray_update(k1,k2, chi, Snu, Inu)
  else
    call self%short_ray_update(k1, k2, chi, Snu, Inu)
  endif

  self%Psi = self%Psi / chi
END SUBROUTINE Inu_update

!===============================================================================
SUBROUTINE long_ray_update (self, k1, k2, chi, Snu, Inu)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:), allocatable :: xxup, ssup
  real(4), dimension(:,:), allocatable :: Inuup
  real(4), dimension(:,:,:,:) :: Inu, chi, Snu
  integer, intent(in) :: k1, k2
  integer :: f, k, z
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%long_ray_update', itimer=itimer)

  allocate(xxup, ssup, mold=chi)
  allocate(Inuup, mold=chi(:,:,1,1))

  call self%interp_xx_ss_up(xxup,ssup)

  do f=1, self%nv
  ! ----------------------------------------------------------------------------
  ! frequency vectorisation is slower, because k varies quicker in mem than f.
  ! ----------------------------------------------------------------------------
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1

      call self%solve_long_plane(k, f, Inuup)
      call linear_formsolv(self%ray%dsup, xxup(:,:,z,f), chi(:,:,z,f), ssup(:,:,z,f), Snu(:,:,z,f), Inuup, self%Psi(:,:,z,f), Inu(:,:,z,f))
      call self%periodic_plane(self%Inu(:,:,k,f))
    end do ! k-loop
  end do ! f-loop

  call trace%end(itimer)
END SUBROUTINE long_ray_update

!===============================================================================
SUBROUTINE short_ray_update (self, k1, k2, chi, Snu, Inu)
  class(radiation_t), target :: self
  real(4), dimension(:,:,:,:), allocatable :: xxup, ssup
  real(4), dimension(:,:), allocatable :: Inuup
  real(4), dimension(:,:,:,:) :: Inu, chi, Snu
  integer, intent(in) :: k1, k2
  integer :: f, k, z
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('radiation_t%short_ray_update', itimer=itimer)

  allocate(xxup, ssup, mold=chi)
  allocate(Inuup, mold=chi(:,:,1,1))

  call self%interp_xx_ss_up(xxup,ssup)

  do f=1, self%nv
  ! ----------------------------------------------------------------------------
  ! frequency vectorisation is slower, because k varies quicker in mem than f.
  ! ----------------------------------------------------------------------------
    do k=k1,k2,self%dir
      z = k - self%li(3) + 1

      call self%interp_xy_up_plane(self%Inu(:,:,k-self%dir,f), Inuup, self%ray%upo1, self%ray%upo2, self%ray%upkern)
      call linear_formsolv(self%ray%dsup, xxup(:,:,z,f), chi(:,:,z,f), ssup(:,:,z,f), Snu(:,:,z,f), Inuup, self%Psi(:,:,z,f), Inu(:,:,z,f))
      call self%periodic_plane(self%Inu(:,:,k,f))
    end do ! k-loop
  end do ! f-loop

  call trace%end(itimer)
END SUBROUTINE short_ray_update

!===============================================================================
pure subroutine interp_xx_ss_up(self, xxup, ssup)
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:,:,:), intent(inout) :: xxup,ssup
  real(8), dimension(:,:), allocatable :: kernel
  real :: xxmin, ssmin
  integer :: x1,x2,y1,y2,z1,z2,i,j,k,f, i1,j1,k1, d1,d2
  ! ----------------------------------------------------------------------------
  xxup = 0.0
  ssup = 0.0

  if (self%ray%upplane /= 3) then
    kernel = self%ray%upkern(2:3,2:3)
  else
    kernel = self%ray%upkern
    i1 = self%li(1) - 1 + self%ray%upo1
    j1 = self%li(2) - 1 + self%ray%upo2
  endif

  k1 = self%li(3) - 1 - self%dir

  d1 = size(kernel, 1) - 1
  d2 = size(kernel, 2) - 1

  if (self%ray%upplane == 1) then
    i1 = self%li(1) - 1 + self%ray%iup
    j1 = self%li(2) - 1 + self%ray%upo1

    do f=1,self%nv
    do k=1,self%n(3)
    do j=1,self%n(2)
    do i=1,self%n(1)
      x1 = i + i1
      y1 = j + j1
      y2 = y1 + d1
      z1 = k + k1
      z2 = z1 + d2
      xxup(i,j,k,f) = xxup(i,j,k,f) + sum(kernel * self%chi ( x1, y1:y2, z1:z2:self%dir, f))
      ssup(i,j,k,f) = ssup(i,j,k,f) + sum(kernel * self%Snu ( x1, y1:y2, z1:z2:self%dir, f))
    end do
    end do
    end do
    end do

  elseif(self%ray%upplane == 2) then
    i1 = self%li(1) - 1 + self%ray%upo1
    j1 = self%li(2) - 1 + self%ray%jup

    do f=1,self%nv
    do k=1,self%n(3)
    do j=1,self%n(2)
    do i=1,self%n(1)
      x1 = i + i1
      x2 = x1 + d1
      y1 = j + j1
      z1 = k + k1
      z2 = z1 + d2
      xxup(i,j,k,f) = xxup(i,j,k,f) + sum(kernel * self%chi ( x1:x2, y1, z1:z2:self%dir, f))
      ssup(i,j,k,f) = ssup(i,j,k,f) + sum(kernel * self%Snu ( x1:x2, y1, z1:z2:self%dir, f))
    end do
    end do
    end do
    end do

  elseif(self%ray%upplane == 3) then
    do f=1,self%nv
    do k=1,self%n(3)
    do j=1,self%n(2)
    do i=1,self%n(1)
      x1 = i + i1
      x2 = x1 + d1
      y1 = j + j1
      y2 = y1 + d2
      z1 = k + k1

      xxup(i,j,k,f) = xxup(i,j,k,f) + sum(self%ray%upkern * self%chi ( x1:x2, y1:y2, z1, f))
      ssup(i,j,k,f) = ssup(i,j,k,f) + sum(self%ray%upkern * self%Snu ( x1:x2, y1:y2, z1, f))

      ! suppress undershoot
      xxmin = minval(self%chi ( x1+1:x2-1, y1+1:y2-1, z1, f))
      ssmin = minval(self%Snu ( x1+1:x2-1, y1+1:y2-1, z1, f))
      if (xxup(i,j,k,f) < xxmin) xxup(i,j,k,f) = xxmin
      if (ssup(i,j,k,f) < ssmin) ssup(i,j,k,f) = xxmin
    end do
    end do
    end do
    end do
  endif
end subroutine interp_xx_ss_up

!===============================================================================
pure subroutine interp_xy_up_plane(self, VarIn, VarUp, up1, up2, kernel)
  class(radiation_t), intent(in) :: self
  integer :: x1,x2,y1,y2,i,j, dx,dy, i1,j1
  integer, intent(in) :: up1, up2
  real(4) :: VarMin
  real(4), dimension(:,:), intent(in) :: VarIn
  real(4), dimension(:,:), intent(out) :: VarUp
  real(8), dimension(:,:), intent(in) :: kernel
  !-----------------------------------------------------------------------------
  dx = size(kernel, 1) - 1
  dy = size(kernel, 2) - 1
  i1 = self%li(1) - 1 + up1
  j1 = self%li(2) - 1 + up2

  VarUp = 0.0

  do j=1,self%n(2)
    y1 = j + j1
    y2 = y1 + dy

  do i=1,self%n(1)

      x1 = i + i1
      x2 = x1 + dx

      VarUp(i,j) = VarUp(i,j) + sum(kernel * VarIn(x1:x2, y1:y2))

      ! suppress undershoot
      VarMin = minval(VarIn(x1+1:x2-1, y1+1:y2-1))
      if (VarUp(i,j) < VarMin) VarUp(i,j) = VarMin
  end do
  end do
end subroutine interp_xy_up_plane

!===============================================================================
subroutine solve_long_plane(self, k, f, Inuup)
  !---------------------------------------------------------------------
  !
  ! Perform a formal solution with linear interpolation of the source
  ! function from the upwind XY crossing to the last cell intersection
  ! before the center point
  !
  !---------------------------------------------------------------------
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:), allocatable :: xxup,ssup, xx, ss, iiup, ii
  real(4), dimension(:,:), intent(out) :: Inuup
  real(4), dimension(:,:), allocatable :: psi_tmp
  integer, intent(in) :: k, f
  integer ::  x1,x2,y1,y2,z1,i,j,m
  ! ----------------------------------------------------------------------------
  allocate(xxup, ssup, iiup, xx, ss, psi_tmp, mold=Inuup)
  ! X,S and incoming intensity  at the upwind XY plane
  call self%interp_xy_up_plane(self%Chi(:,:,k-self%dir,f), xxup, self%ray%longo1, self%ray%longo2, self%ray%longkern)
  call self%interp_xy_up_plane(self%Snu(:,:,k-self%dir,f), ssup, self%ray%longo1, self%ray%longo2, self%ray%longkern)
  call self%interp_xy_up_plane(self%Inu(:,:,k-self%dir,f), iiup, self%ray%longo1, self%ray%longo2, self%ray%longkern)

  ! linear interpolation of snu and chi at each xz or yz intersection
  do m=1,self%ray%Nintersect
    xx = 0.0
    ss = 0.0

    do i=1,2
      do j=1,2 
        x1 = self%ray%int(m)%iof(i,j) + self%li(1)
        x2 = self%ray%int(m)%iof(i,j) + self%ui(1)
        y1 = self%ray%int(m)%jof(i,j) + self%li(2)
        y2 = self%ray%int(m)%jof(i,j) + self%ui(2)
        z1 = k - self%dir*abs(self%ray%int(m)%kof(i,j))

        xx(:,:) = xx(:,:) + self%ray%int(m)%wght(i,j) * self%chi( x1:x2, y1:y2, z1, f)
        ss(:,:) = ss(:,:) + self%ray%int(m)%wght(i,j) * self%Snu( x1:x2, y1:y2, z1, f)
      end do
    end do

    call linear_formsolv(self%ray%int(m)%ds, xxup, xx, ssup, ss, iiup, psi_tmp, Inuup)

    ! set new start point to old end point
    xxup = xx
    ssup = ss
    iiup = Inuup
  end do
end subroutine solve_long_plane

!====================================================================================
! Does the formal solution between the upwind and central point
!====================================================================================
elemental subroutine linear_formsolv(ds, xx1, xx2, ss1, ss2, Inuup, Psi, Inu)
  use mconst
  real(8), intent(in) :: ds
  real(4), intent(in) :: xx1, xx2, ss1, ss2, Inuup
  real(4), intent(out) :: Psi
  real(4), intent(out) :: Inu
  real(8) :: dt, ex, u0
  ! ----------------------------------------------------------------------------

  dt = 0.5d0 * ds * ( xx1 + xx2 ) ! upwind optical pathlength
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
SUBROUTINE periodic_plane (self, arr)
  class(radiation_t), intent(in) :: self
  real(4), dimension(:,:), intent(inout) :: arr
  integer, dimension(3) :: li, ui, lb, ub
  integer :: dir, is, ie, off
  ! ----------------------------------------------------------------------------
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

  filename = trim(out_folder) // trim(self%ray_name) // '_lam.bin'
  open(newunit=lunit, file=filename, access='stream')
  write(lunit) lam
  call self%save_meta4(filename, pack(lam, .true.), shape(lam))
  close(lunit)

  filename = trim(out_folder) // trim(self%ray_name) // '_inu.bin'
  open(newunit=lunit, file=filename, access='stream')
  write(lunit) Inu
  call self%save_meta4(filename, pack(Inu, .true.), shape(Inu))
  close(lunit)
  
  call trace%end()
END SUBROUTINE save_radiation_output

!===============================================================================
SUBROUTINE setup_metafile(self)
  class(radiation_t),intent(inout):: self
  integer :: meta_unit, iostat
  ! ----------------------------------------------------------------------------

  write(self%ray_name,'(A1, I3.3, A4, I3.3)') 'p', self%patch%id, '_ray', self%i_omega

  self%meta_file = trim(out_folder) // trim(self%ray_name) // '_meta.txt'
  open(newunit=meta_unit, file=self%meta_file, form='formatted')

  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'position', self%patch%position
  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'vec', self%vec

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

end module
