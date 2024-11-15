MODULE read_atmos_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE handyfuncs_mod
  USE const_mod
  USE op_math_mod, only: linspace, linterpolate_3d, pchip_mono4, interp1d_sorted
  USE m3d_patch_mod, only: pivot_z, ray_offset, amr, composition, use_ne, use_rho, make_eos
  USE ray_mod, only: one_dimensional
  USE composition_mod, only: composition_type
  implicit none
  private

  type, public :: input_atmos_t
    character(64), allocatable :: kind
    integer :: n(3), nx, ny, nz
    integer :: gz=0 ! ghost zones
    logical :: is_3d=.false.
    real(4), dimension(:,:,:), allocatable :: ltemp, lpe, lrho
    real(4), dimension(:,:,:), allocatable :: vx, vy, vz, vmic
    real(4), dimension(:),     allocatable :: xx, yy, zz
    real(8), dimension(:),     allocatable :: sub_pos, sub_size
  contains
    procedure :: allocate   => allocate_input_atmos
    procedure :: deallocate => deallocate_input_atmos
    procedure :: get_box
    procedure :: thick2thin
    procedure :: setup_amr_mesh
  end type input_atmos_t

  ! stagger coefficients
  real(4), parameter :: uul = 1.0e+8      ! length
  real(4), parameter :: uur = 1.0e-7      ! density
  real(4), parameter :: luul = log(uul)
  real(4), parameter :: luur = log(uur)
  real(4), parameter :: uut = 1.0e2       ! time
  real(4), parameter :: uuv = uul / uut   ! velocity
  real(4), parameter :: uum = uuv * uur   ! momentum
  real(4), parameter :: uue = uuv**2      ! internal energy per unit mass
  real(4), parameter :: uup = uur * uue   ! pressure
  real(4), parameter :: uua = 1./uul      ! extinction coefficient (alpha)
  !character(len=128), public :: stagger_file='/lustre/astro/rhoppe/multi-3d/input/atmos/t5777g44m00/t5777g44m0005'
  character(len=256), save :: mesh_file, snap_file, data_file, aux_file
  real(4), save :: maxTemp=huge(1.0)
  real(4), save :: ee_offset=3.2122
  logical, public, save :: legacy_stagger_interp=.false.
  type(input_atmos_t), public, save :: ia

  public :: read_stagger
  public :: read_stagger2
  public :: read_co5bold
  public :: read_marcs
  public :: read_multi_snapshot
  public :: read_1d_snapshot
  public :: interpolate_atmos
    
contains

!===============================================================================
SUBROUTINE read_co5bold(atmos_file)
  character(len=*), intent(in) :: atmos_file
  real(4), dimension(:,:,:), allocatable :: over_rho
  integer :: lunit, iostat, nsub
  real(4) :: model_time
  real(8) :: ds(3)
  integer :: n(3), nx, ny, nz
  !-----------------------------------------------------------------------------
  open(newunit=lunit, file=atmos_file, form='unformatted', access='stream')
  read(lunit) model_time
  read(lunit) nx, ny, nz

  call ia%allocate(nx, ny, nz)

  read(lunit) ia%xx, ia%yy, ia%zz
  read(lunit) ia%ltemp
  read(lunit) ia%lrho 
  read(lunit) ia%lpe  
  read(lunit) ia%vx ! this is px
  read(lunit) ia%vy ! this is py
  read(lunit) ia%vz ! this is pz

  allocate(over_rho, source=exp(-ia%lrho))
  ia%vx = ia%vx * over_rho
  ia%vy = ia%vy * over_rho
  ia%vz = ia%vz * over_rho
  close(lunit)
END SUBROUTINE read_co5bold

!===============================================================================
SUBROUTINE read_stagger(stagger_file, snap, use_ne)
  character(len=*), intent(in) :: stagger_file
  integer, intent(in) :: snap
  logical, intent(inout) :: use_ne
  real(8) :: ds(3)
  integer :: n(3), nx, ny, nz
  !-----------------------------------------------------------------------------
  mesh_file = trim(stagger_file) // '.msh'
  write(snap_file, '(A,A1,I5.5)') trim(stagger_file), '_', snap
  write(data_file, '(A,A)') trim(snap_file), '.dat'
  write(aux_file,  '(A,A)') trim(snap_file), '.aux'
  call read_stagger_mesh()
  call read_stagger_data(use_ne)
END SUBROUTINE read_stagger

!===============================================================================
SUBROUTINE read_stagger_dim(nx, ny, nz, lunit_out)
  integer, intent(out) :: nx, ny, nz
  integer :: lunit, iostat
  integer, optional :: lunit_out
  !-----------------------------------------------------------------------------
  call os%check_file(mesh_file)
  open(newunit=lunit, file=mesh_file, iostat=iostat, convert='big_endian', form='unformatted')

  ! read first entry in mesh file
  read(unit=lunit, iostat=iostat) nx, nz, ny ! y is z in stagger

  if (present(lunit_out)) then
    lunit_out = lunit
  else
    close(lunit)
  endif
END SUBROUTINE read_stagger_dim

!===============================================================================
!> Reading mesh file
!===============================================================================
SUBROUTINE read_stagger_mesh()
  real(4), allocatable, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
  real(4), allocatable, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
  real(4), allocatable, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
  real(4), allocatable, dimension(:) :: zm
  integer :: lunit, iostat, nx, ny, nz
  !-----------------------------------------------------------------------------
  call read_stagger_dim(nx, ny, nz, lunit_out=lunit)
  call ia%allocate(nx, ny, nz)
  ia%gz = 5 ! stagger has 5 ghost zones on each side

  allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=ia%xx)
  allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=ia%yy)
  allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=ia%zz)

  read (unit=lunit, iostat=iostat) dx, dxdn, ia%xx, xdn, dxidxup, dxidxdn
  read (unit=lunit, iostat=iostat) dz, dzdn, ia%zz, zdn, dzidzup, dzidzdn ! y is z in stagger
  read (unit=lunit, iostat=iostat) dy, dydn, ia%yy, ydn, dyidyup, dyidydn
  ! close mesh file
  close(unit=lunit)
  ia%xx = uul * ia%xx
  ia%yy = uul * ia%yy
  ia%zz = uul * ia%zz
END SUBROUTINE read_stagger_mesh

!===============================================================================
!> Reading data and aux file
!===============================================================================
SUBROUTINE read_stagger_data(use_ne)
  logical, intent(inout) :: use_ne
  real(4), dimension(:,:,:), allocatable :: ppx, ppy, ppz, lne, rho
  integer :: iostat, lunit, recl
  integer :: nx, ny, nz, m(3), i, j, k
  integer, save :: itimer
  real(4) :: const, over_rho
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_data', itimer=itimer)
  call read_stagger_dim(nx, ny, nz)

  m = [nx, ny, nz]
  allocate(ppx(nx,ny,nz))
  allocate(ppy(nx,ny,nz))
  allocate(ppz(nx,ny,nz))
  allocate(rho(nx,ny,nz))
  allocate(lne(nx,ny,nz))

  recl = kind(1.0)*product(m)
  use_ne = use_ne .and. os%exists(aux_file)

  !=============================================================================
  ! Reading data file
  !-----------------------------------------------------------------------------
  call os%check_file(data_file)
  open(newunit=lunit, file=data_file, iostat=iostat, convert='big_endian', &
       form='unformatted', status='old', access='direct', recl=recl)

    ! 1: density
    ! 2: momentum per unit volume in direction x
    ! 3: momentum per unit volume in direction y
    ! 4: momentum per unit volume in direction z
    ! 5: internal energy per unit volume
    ! 6: temperature

  read (unit=lunit, iostat=iostat, rec=1) rho
  read (unit=lunit, iostat=iostat, rec=2) ppx
  read (unit=lunit, iostat=iostat, rec=3) ppz ! y is z in stagger
  read (unit=lunit, iostat=iostat, rec=4) ppy ! z is y in stagger
  !read (unit=lunit, iostat=iostat, rec=5) ie
  if (.not. use_ne) then
    read (unit=lunit, iostat=iostat, rec=6) ia%ltemp ! this is temp
    ia%ltemp = log(ia%ltemp)
  endif
  ! close data file
  close(unit=lunit)

  if (use_ne) then
    !=============================================================================
    ! Reading aux file
    !-----------------------------------------------------------------------------
    call os%check_file(aux_file)
    open(newunit=lunit, file=aux_file, iostat=iostat, convert='big_endian', &
          form='unformatted', status='old', access='direct', recl=recl)

      ! 1: ln total pressure (gas + rad)
      ! 2: ln Rosseland extinction coeff.
      ! 3: ln Temperature
      ! 4: ln electron number density
      ! 5: ln Planck-averaged extinction coeff.
      ! 6: ln extinction coeff at 500 nm

    read (unit=lunit, iostat=iostat, rec=3) ia%ltemp
    read (unit=lunit, iostat=iostat, rec=4) lne
    ! close aux file
    close(unit=lunit)
    !=============================================================================
  else
    ia%lpe = 0
  endif

  ! "un"-staggering vectors relative to scalars
  call xup_set(ppx, ia%vx)

  ! transposing data (y to z)
  call transpose_stagger(rho)
  call transpose_stagger(ia%ltemp)
  if (use_ne) call transpose_stagger(lne)

  call transpose_stagger(ia%vx)
  call transpose_stagger(ppy)
  call transpose_stagger(ppz)

  call yup_set(ppy, ia%vy)
  call zup_set(ppz, ia%vz)

  const = -3.0 * luul + log(BK)
  do k=1, nz
  do j=1, ny
  do i=1, nx
    over_rho = uuv / rho(i,j,k)
    ia%lrho (i,j,k) = log(rho(i,j,k)) + luur
    ia%vx   (i,j,k) =   ia%vx(i,j,k) * over_rho
    ia%vy   (i,j,k) =   ia%vy(i,j,k) * over_rho
    ia%vz   (i,j,k) =   ia%vz(i,j,k) * over_rho
    if (use_ne) then
      ia%lpe(i,j,k) = lne(i,j,k) + const + ia%ltemp(i,j,k)
    endif
  end do
  end do
  end do

  call trace%end(itimer)
END SUBROUTINE read_stagger_data

!===============================================================================
SUBROUTINE read_multi_snapshot(atmos_file)
  integer :: lunit, nx, ny, nz
  character(len=*) :: atmos_file
  real(4), dimension(:,:,:), allocatable :: ne, temp, rho, vx, vy, vz
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  call os%check_file(trim(atmos_file) // '.mesh')
  open(file=trim(atmos_file) // '.mesh', newunit=lunit, form='FORMATTED')
  read(lunit, *) nx, ny, nz
  
  call ia%allocate(nx, ny, nz)
  allocate(temp, ne, vx, vy, vz, rho, mold=ia%lrho)

  read(lunit, *) ia%xx
  read(lunit, *) ia%yy
  read(lunit, *) ia%zz
  close(lunit)

  call os%check_file(trim(atmos_file) // '.atmos')
  open(file=trim(atmos_file) // '.atmos', newunit=lunit, form='UNFORMATTED', access='stream')
  read(lunit) ne
  read(lunit) temp
  read(lunit) ia%vx
  read(lunit) ia%vy
  read(lunit) ia%vz
  read(lunit) rho
  close(lunit)

  ia%lpe = log(ne * BK * temp)
  ia%ltemp = log(temp)
  ia%lrho  = log(rho)
END SUBROUTINE read_multi_snapshot

!===============================================================================
SUBROUTINE read_stagger2(atmos_file)
  integer :: lunit, nx, ny, nz, pos
  character(len=*) :: atmos_file
  real(4), dimension(:), allocatable :: xx, yy, zz
  real(4), dimension(:,:,:), allocatable :: ne, temp, rho, vx, vy, vz
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  pos = index(atmos_file, 'atm3d.')
  mesh_file = atmos_file(1:pos-1) // 'mesh.' // atmos_file(pos+6:)
  call os%check_file(mesh_file)
  open(file=mesh_file, newunit=lunit, form='FORMATTED')

  read(lunit, *) nx
  allocate(xx(nx))
  read(lunit, *) xx

  read(lunit, *) ny
  allocate(yy(ny))
  read(lunit, *) yy

  read(lunit, *) nz
  allocate(zz(nz))
  read(lunit, *) zz
  close(lunit)

  call ia%allocate(nx, ny, nz)
  ia%xx = xx
  ia%yy = yy
  ia%zz = zz
  
  allocate(temp, ne, vx, vy, vz, rho, mold=ia%lrho)
  call os%check_file(atmos_file)
  open(file=atmos_file, newunit=lunit, form='unformatted', access='stream')
  read(lunit) ne
  read(lunit) temp
  read(lunit) vx
  read(lunit) vy
  read(lunit) vz
  read(lunit) rho
  close(lunit)

  ia%vx = vx * 1e5
  ia%vy = vy * 1e5
  ia%vz = vz * 1e5
  ia%lpe = log(ne * BK * temp)
  ia%ltemp = log(temp)
  ia%lrho  = log(rho)
END SUBROUTINE read_stagger2

!===============================================================================
SUBROUTINE read_marcs (file, FeH, vmic, abund)
  character(len=*), intent(in) :: file
  real(8), intent(out) :: FeH
  real(4), intent(in) :: vmic
  real(4), intent(inout), dimension(:), allocatable :: abund
  integer :: lunit, nz, kk, k
  real(8) :: FeH_marcs
  real(4) :: teff, flux, vmic_Marcs, g, rmass, alpha, radius, rluminosity
  real(4) :: xx, yy, zz, palfa, pny, py, pbeta
  real(4) :: lgTauR, lgTau5, pg, prad, pturb, xkapr, emu, vconv, fconv, rhox
  real(4), dimension(:), allocatable :: temp, pe, rho
  character(512) :: fileID
  !-----------------------------------------------------------------------------
  allocate(abund(92))
  call os%check_file(file)
  open(newunit=lunit, file=file)
  read(lunit,*) fileID
  read(lunit,*) teff                  ! Teff [K].
  read(lunit,*) flux                  ! Flux [erg/cm2/s]
  read(lunit,*) g                     ! Surface gravity [cm/s2]
  read(lunit,*) vmic_Marcs            ! Microturbulence parameter [km/s]
  read(lunit,*) rmass                 ! No mass for plane-parallel models
  read(lunit,*) FeH_marcs, alpha      ! Metallicity [Fe/H] and [alpha/Fe]
  read(lunit,*) radius                ! 1 cm radius for plane-parallel models
  read(lunit,*) rluminosity           ! Luminosity [Lsun] FOR A RADIUS OF 1 cm!
  read(lunit,*) palfa, pny, py, pbeta ! Convection parameters: alpha, nu, y and beta
  read(lunit,*) xx, yy, zz            ! X, Y and Z, 12C/13C=89 (=solar)
  read(lunit,*) ! skip comment
  read(lunit,*) abund
  read(lunit,*) nz
  read(lunit,*) ! skip comment
  read(lunit,*) ! skip comment

  call ia%allocate(1,1,nz)
  allocate(temp(nz), pe(nz), rho(nz))

  if (FeH < -99) then
    FeH = FeH_marcs
    write(stdout, '(A, F6.2)') 'FeH was overwritten by Marcs model: [Fe/H] =', FeH
  endif

  if (vmic < 0) then
    ia%vmic (1,1,:) = vmic_Marcs
    write(stdout, '(A, F6.2)') 'Vmic was overwritten by Marcs model: Vmic =', vmic_Marcs
  else
    ia%vmic (1,1,:) = vmic
  endif
  
  do kk=1, nz
    read(lunit,*) k, lgTauR, lgTau5, ia%zz(k), temp(k), pe(k), pg, prad, pturb
    if(kk /= k) stop 'model and loop indices differ'
  end do
  
  read(lunit,*) ! skip comment
  
  do kk=1, nz
    read(lunit,*) k, lgTauR, xkapr, rho(k), emu, vconv, fconv, rhox
    if(kk /= k) stop 'model and loop indices differ'
  end do
  close(lunit)

  ia%ltemp(1,1,:) = log(temp)
  ia%lpe  (1,1,:) = log(pe)
  ia%lrho (1,1,:) = log(rho)
END SUBROUTINE read_marcs

!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
SUBROUTINE read_1d_snapshot(atmos_file, vmic)
  character(len=*) :: atmos_file
  real(4), intent(in) :: vmic
  integer :: nz, lunit, i
  real(4), dimension(:), allocatable :: temp, pe, rho
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  call os%check_file(atmos_file)
  open(file=atmos_file, newunit=lunit)
  call nextentry(lunit)
  read(lunit, *) ! skip name
  call nextentry(lunit)
  read(lunit, *) nz

  call ia%allocate(1,1,nz)
  allocate(temp(nz), pe(nz), rho(nz))

  do i=1, nz
    call nextentry(lunit)
    read(lunit, *) ia%zz(i), temp(i), pe(i), rho(i), ia%vmic(1,1,i)
  end do
  close(lunit)

  if (vmic >= 0) then
    ia%vmic = vmic
    write(stdout, '(A, F6.2)') 'Vmic was overwritten by namelist: Vmic =', vmic
  endif

  use_ne = use_ne .and. pe(1) /= 0
  use_rho = use_rho .and. rho(1) /= 0
  where (temp > 0) ia%ltemp(1,1,:) = log(temp)
  where (rho  > 0) ia%lrho (1,1,:) = log(rho)
  where (pe   > 0) ia%lpe  (1,1,:) = log(pe)
END SUBROUTINE read_1d_snapshot

!===============================================================================
SUBROUTINE xup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: j, k
  !--------------------------------------------------------------------------
  do k=1, size(f, 3)
    do j=1, size(f, 2)
      call destagger_vector(f(:,j,k), fp(:,j,k))
    end do
  end do
END SUBROUTINE xup_set
    
!===============================================================================
SUBROUTINE yup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: i, k
  !--------------------------------------------------------------------------
  do k=1, size(f, 3)
    do i=1, size(f, 1)
      call destagger_vector(f(i,:,k), fp(i,:,k))
    end do
  end do
END SUBROUTINE yup_set

!===============================================================================
SUBROUTINE zup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: i, j
  !--------------------------------------------------------------------------
  do j=1, size(f, 2)
    do i=1, size(f, 1)
      call destagger_vector(f(i,j,:), fp(i,j,:))
    end do
  end do
END SUBROUTINE zup_set

!===============================================================================
SUBROUTINE destagger_vector(f, fout)
  real, dimension(:), intent(in)  :: f
  real, dimension(:), intent(out) :: fout
  integer :: i, n, j, idx
  real(4), parameter :: c = 3./256.
  real(4), parameter :: b = -25./256.
  real(4), parameter :: a = 0.5-b-c
  real(4) :: wts(-2:3)=[c,b,a,a,b,c]
  ! ----------------------------------------------------------------------------
  n = size(f,1)
  do i=1, n
    do j=-2,3
      idx = i + j
      if (idx<1) then
        idx = n + idx
      elseif (idx>n) then
        idx = idx-n
      endif  
      fout(i) = fout(i) + wts(j) * f(idx)
    end do
  end do
END SUBROUTINE destagger_vector

!===============================================================================
SUBROUTINE transpose_stagger(f)
  real, dimension(:,:,:), intent(inout)  :: f
  integer :: i, itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('transpose_stagger', itimer=itimer)
  do i=1, size(f, 1)
    f(i, :, :) = transpose(f(i, :, :))
  end do
  call trace%end(itimer)
END SUBROUTINE transpose_stagger

!===============================================================================
SUBROUTINE allocate_input_atmos(self, nx, ny, nz)
  class(input_atmos_t), intent(inout) :: self
  integer, intent(in) :: nx, ny, nz
  !-----------------------------------------------------------------------------
  self%n  = [nx, ny, nz]
  self%nx = nx
  self%ny = ny
  self%nz = nz
  self%is_3d = nx*ny > 1
  allocate(self%xx(nx), source=0.)
  allocate(self%yy(ny), source=0.)
  allocate(self%zz(nz), source=0.)
  allocate(self%vx   (nx,ny,nz), source=0.)
  allocate(self%vy   (nx,ny,nz), source=0.)
  allocate(self%vz   (nx,ny,nz), source=0.)
  allocate(self%vmic (nx,ny,nz), source=0.)
  allocate(self%lrho (nx,ny,nz))
  allocate(self%ltemp(nx,ny,nz))
  allocate(self%lpe  (nx,ny,nz))
END SUBROUTINE allocate_input_atmos

!===============================================================================
SUBROUTINE deallocate_input_atmos(self)
  class(input_atmos_t), intent(inout) :: self
  !-----------------------------------------------------------------------------
  if (allocated(self%lrho))  deallocate(self%lrho)
  if (allocated(self%ltemp)) deallocate(self%ltemp)
  if (allocated(self%lpe))   deallocate(self%lpe)
  if (allocated(self%vx))    deallocate(self%vx)
  if (allocated(self%vy))    deallocate(self%vy)
  if (allocated(self%vz))    deallocate(self%vz)
  if (allocated(self%vmic))  deallocate(self%vmic)
END SUBROUTINE deallocate_input_atmos

!===============================================================================
SUBROUTINE get_box(self, box)
  class(input_atmos_t), intent(inout) :: self
  real(8), intent(out) :: box(3)
  real(8) :: ds(3)
  integer :: n(3)
  !-----------------------------------------------------------------------------
  n = [self%nx, self%ny, self%nz-2*self%gz]
  box(1) = self%xx(self%nx)
  box(2) = self%yy(self%ny)
  box(3) = self%zz(self%nz-self%gz)
  box(1:2) = box(1:2) * n(1:2) / max(n(1:2)-1, 1)
  if (self%nx == 1) box(1) = 1d8
  if (self%ny == 1) box(2) = 1d8
END SUBROUTINE get_box

!===============================================================================
! This makes sure that the atmosphere is flipped properly
!===============================================================================
SUBROUTINE thick2thin(self)
  class(input_atmos_t), intent(inout) :: self
  type(composition_type) :: comp
  real(8), dimension(:,:,:), allocatable :: temp8, ne8, rho8
  real(8), dimension(:), allocatable :: chi, tau
  real(8) :: dz
  integer :: i, j, nx, ny, nz
  !-----------------------------------------------------------------------------
  nx = self%nx
  ny = self%ny
  nz = self%nz
  if (self%lrho(1,1,1) < self%lrho(1,1,nz)) then
    do j=1, ny
    do i=1, nx
      self%ltemp(i,j,:) =  self%ltemp(i,j,nz:1:-1)
      self%lrho (i,j,:) =  self%lrho (i,j,nz:1:-1)
      self%lpe  (i,j,:) =  self%lpe  (i,j,nz:1:-1)
      self%vx   (i,j,:) =  self%vx   (i,j,nz:1:-1)
      self%vy   (i,j,:) =  self%vy   (i,j,nz:1:-1)
      self%vz   (i,j,:) = -self%vz   (i,j,nz:1:-1)
      self%vmic (i,j,:) =  self%vmic (i,j,nz:1:-1)
    end do
    end do
    self%zz = self%zz (nz:1:-1)
  endif

  if (nx*ny==1 .and. abs(ia%zz(1) - ia%zz(nz)) < 1e4) then
    ! log tau scale
    allocate(temp8(nx, ny, nz))
    allocate(rho8 (nx, ny, nz))
    allocate(ne8  (nx, ny, nz))
    allocate(chi(nz))
    allocate(tau(nz))
    temp8 = exp(dble(self%ltemp))
    rho8  = exp(dble(self%lrho))
    ne8   = exp(dble(self%lpe)) / (BK * temp8)
    comp = composition
    call comp%reassign_pointer()
    call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne)
    chi=comp%chi500(1,1,:)
    tau=10**(dble(ia%zz))

    ia%zz(nz) = 0
    do i=nz-1, 1, -1
      dz = (tau(i) - tau(i+1)) * 2 / (chi(i) + chi(i+1))
      ia%zz(i) = ia%zz(i+1) + dz
    end do
  endif

  if (self%zz(1) > self%zz(nz)) then
    self%zz = -self%zz
  endif

  self%xx = self%xx - self%xx(1)         ! let the origin be zero
  self%yy = self%yy - self%yy(1)         ! let the origin be zero
  self%zz = self%zz - self%zz(self%gz+1) ! let the origin be the layer above ghost zone
END SUBROUTINE thick2thin

!===============================================================================
SUBROUTINE setup_amr_mesh(self, nsub, use_ne)
  class(input_atmos_t), intent(inout), target :: self
  integer, intent(in) :: nsub
  logical, intent(in) :: use_ne
  integer :: i, j, k, nx, ny, nz
  real(8), dimension(nsub+1) :: var_edge
  real(4), dimension(nsub+1) :: sub_edge
  real(4), dimension(:), allocatable :: amr_ref, tmp
  real(4), dimension(:), pointer, contiguous :: zz
  real(4), dimension(:,:,:), pointer, contiguous :: lpe
  real(4) :: vmin, vmax
  !-----------------------------------------------------------------------------
  nx = self%nx
  ny = self%ny
  nz = self%nz - 2*self%gz
  allocate(amr_ref(nz), source=0.0)
  allocate(tmp(nz), source=0.0)
  allocate(self%sub_pos(nsub), self%sub_size(nsub))

  zz  => self%zz     (1+self%gz:nz+self%gz)
  if (use_ne) then
    lpe => self%lpe(:,:,1+self%gz:nz+self%gz)
  else
    lpe => self%ltemp(:,:,1+self%gz:nz+self%gz)
  endif

  do k=1, nz
  do j=1, ny
  do i=1, nx
    ! we use minus here, because interp1d_sorted expects increasing values
    amr_ref(k) = amr_ref(k) - lpe(i,j,k)
  end do
  end do
  end do
  amr_ref = amr_ref / (nx*ny)

  ! we want to sample finely where the gradient is large
  ! it is ok to ensure monotonicity here by taking the absolute value of the gradient
  tmp(1) = amr_ref(1)
  do k=2, nz
    tmp(k) = tmp(k-1) + abs(amr_ref(k) - amr_ref(k-1))
  end do
  amr_ref = tmp

  vmin = amr_ref(1)
  vmax = amr_ref(nz)
  if (.not. use_ne) then
    amr_ref = (amr_ref - vmin) / (vmax - vmin)
    amr_ref = 1 - (1-amr_ref)**(2.d0/3.d0)
    vmin = 0
    vmax = 1
  endif
  var_edge = linspace(vmin, vmax, nsub+1)
  call interp1d_sorted(amr_ref, zz, var_edge, sub_edge)
  !sub_edge = sub_edge * nz / (nz-1)
  self%sub_size = sub_edge(2:nsub+1) - sub_edge(1:nsub)
  self%sub_pos  = sub_edge(1:nsub) + 0.5*self%sub_size
END SUBROUTINE setup_amr_mesh

!===============================================================================
SUBROUTINE interpolate_atmos(n, sub, pos, rho, px, py, pz, temp, pe, vmic, vec, bot2pivot)
  integer, intent(in) :: n(3)
  real(8), intent(in) :: sub(3), pos(3)
  real(4), intent(out), dimension(:,:,:) :: rho, px, py, pz, temp, pe, vmic
  real(8), intent(in), optional :: vec(3), bot2pivot
  real(8), dimension(:), allocatable :: newx, newy, newz
  real(8) :: lower(3), ds(3)
  real(8) :: dx, dy, dz
  integer :: i, nx, ny, nz, itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('interpolate_atmos', itimer=itimer)
  nx = n(1); ny = n(2); nz = n(3)

  lower = pos - 0.5d0*sub
  ds = sub / n
  
  allocate(newz(nz))
  if (amr) then
    ! this way the distance to the patch below is ds of the current patch, otherwise it is ds of the patch below
    ! this is needed for the amr setup in order to avoid additional interpolations for the final spectrum.  
    newz = [ (ds(3)*i + lower(3), i=1, nz) ]
  else
    newz = [ (ds(3)*i + lower(3), i=0, nz-1) ]
  endif
  if (make_eos) then
    ! no interpolation
    temp = ia%ltemp
    rho  = ia%lrho 
    pe   = ia%lpe  
    vmic = ia%vmic 
    px = 0
    py = 0
    pz = 0
  elseif (one_dimensional) then
    call pchip_mono4(ia%zz, ia%ltemp(1,1,:), real(newz), temp(1,1,:), sorted=.true.)
    call pchip_mono4(ia%zz, ia%lrho (1,1,:), real(newz), rho (1,1,:), sorted=.true.)
    call pchip_mono4(ia%zz, ia%lpe  (1,1,:), real(newz), pe  (1,1,:), sorted=.true.)
    call pchip_mono4(ia%zz, ia%vmic (1,1,:), real(newz), vmic(1,1,:), sorted=.true.)
    px = 0
    py = 0
    pz = 0
  else
    allocate(newx(nx))
    allocate(newy(ny))
    newx = [ (ds(1)*i + lower(1), i=0, nx-1) ]
    newy = [ (ds(2)*i + lower(2), i=0, ny-1) ]
    
    dx = ia%xx(2)
    dy = ia%yy(2)

    if (present(vec)) then
      call linterpolate_3d_rotated(ia%ltemp, temp, dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      call linterpolate_3d_rotated(ia%lrho,  rho,  dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      call linterpolate_3d_rotated(ia%lpe,   pe,   dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      call linterpolate_3d_rotated(ia%vx,    px,   dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      call linterpolate_3d_rotated(ia%vy,    py,   dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      call linterpolate_3d_rotated(ia%vz,    pz,   dx, dy, ia%zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    else
      call linterpolate_3d(ia%ltemp, temp, dx, dy, ia%zz, newx, newy, newz)
      call linterpolate_3d(ia%lrho,  rho,  dx, dy, ia%zz, newx, newy, newz)
      call linterpolate_3d(ia%lpe,   pe,   dx, dy, ia%zz, newx, newy, newz)
      call linterpolate_3d(ia%vx,    px,   dx, dy, ia%zz, newx, newy, newz)
      call linterpolate_3d(ia%vy,    py,   dx, dy, ia%zz, newx, newy, newz)
      call linterpolate_3d(ia%vz,    pz,   dx, dy, ia%zz, newx, newy, newz)
    endif
  endif

  temp = exp(temp)
  rho  = exp(rho)
  pe   = exp(pe)
  px   = px * rho
  py   = py * rho
  pz   = pz * rho
  call trace%end(itimer)
END SUBROUTINE interpolate_atmos

END MODULE read_atmos_mod