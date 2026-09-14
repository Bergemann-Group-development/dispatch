MODULE read_atmos_mod
  USE, intrinsic :: ieee_arithmetic
  USE io_mod
  USE os_mod
  USE trace_mod
  USE handyfuncs_mod
  USE const_mod
  USE op_math_mod, only: linspace
  USE interpol_mod, only: linterpolate_3d, linterpolate_3d_rotated, pchip_mono4, interp1d_sorted
  USE m3d_patch_mod, only: composition, use_ne, use_rho, make_eos, prune_top, prune_tol
  USE rt_params_mod, only: pivot_z, ray_offset
  USE ray_mod, only: one_dimensional, rstar, rorigin
  USE composition_mod, only: composition_type
  USE atmos_reader_mod, only: atmos_data_t, read_marcs_file, read_text_1d_file, read_multi_1d_file
  USE units_mod, only: cgs
  implicit none
  private

  type, public :: input_atmos_t
    character(64) :: kind
    integer :: n(3), nx, ny, nz
    integer :: gz=0 ! ghost zones
    logical :: is_3d=.false.
    real(4), dimension(:,:,:), allocatable :: ltemp, lpe, lrho
    real(4), dimension(:,:,:), allocatable :: vx, vy, vz, vmic
    real(4), dimension(:),     allocatable :: xx, yy, zz
    real(8), dimension(:),     allocatable :: sub_pos, sub_size
    ! The snapshot is shared by every patch on this rank, and each interpolates
    ! from it exactly once.  Patches announce themselves with acquire(), and the
    ! last one to release() frees it.
    integer :: n_users=0
  contains
    procedure :: allocate   => allocate_input_atmos
    procedure :: deallocate => deallocate_input_atmos
    procedure :: acquire    => acquire_input_atmos
    procedure :: release    => release_input_atmos
    procedure :: scale_box
    procedure :: thick2thin
    procedure :: prune_top_planes
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
  real(4), public :: timestamp=-1.0
  logical, public, save :: log_interpolation=.true.
  logical, public, save :: pchip_interpolation=.true.
  type(input_atmos_t), public, save, target :: ia

  public :: read_stagger
  public :: read_stagger2
  public :: read_co5bold
  public :: read_marcs
  public :: read_multi_snapshot
  public :: read_text_1d
  public :: read_multi_1d
  public :: read_siu_bin
  public :: interpolate_atmos

contains

  !===============================================================================
  SUBROUTINE read_co5bold(atmos_file)
    character(len=*), intent(in) :: atmos_file
    real(4), dimension(:,:,:), allocatable :: over_rho
    integer :: lunit, iostat, nsub
    real(8) :: ds(3)
    integer :: n(3), nx, ny, nz
    !-----------------------------------------------------------------------------
    open(newunit=lunit, file=atmos_file, form='unformatted', access='stream')
    read(lunit) timestamp
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
          ia%vx   (i,j,k) = ia%vx(i,j,k) * over_rho
          ia%vy   (i,j,k) = ia%vy(i,j,k) * over_rho
          ia%vz   (i,j,k) = ia%vz(i,j,k) * over_rho
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
    if (endswith(trim(atmos_file), '.atmos')) then
      atmos_file = atmos_file(1:len_trim(atmos_file)-6)
    elseif (endswith(trim(atmos_file), '.mesh')) then
      atmos_file = atmos_file(1:len_trim(atmos_file)-5)
    endif
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
  !===============================================================================
  !> The 1D text formats are parsed by atmos_reader_mod (shared with hop_prune.x);
  !> the routines below apply the namelist overrides and copy the result into ia.
  !===============================================================================
  SUBROUTINE read_marcs (file, FeH, vmic, abund)
    character(len=*), intent(in) :: file
    real(8), intent(out) :: FeH
    real(4), intent(in) :: vmic
    real(4), intent(inout), dimension(:), allocatable :: abund
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_marcs_file, file, atmos)
    allocate(abund, source=atmos%abund)
    rstar = atmos%radius

    if (FeH < -99) then
      FeH = atmos%feh
      write(stdout, '(A, F6.2)') 'FeH was overwritten by Marcs model: [Fe/H] =', FeH
    endif

    if (vmic < 0) then
      write(stdout, '(A, F6.2)') 'Vmic was overwritten by Marcs model: Vmic =', atmos%vmic(1,1,1)
    else
      ia%vmic = vmic
    endif
  END SUBROUTINE read_marcs

!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
  SUBROUTINE read_text_1d(atmos_file, vmic, FeH)
    character(len=*) :: atmos_file
    real(4), intent(in) :: vmic
    real(8), intent(out) :: FeH
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_text_1d_file, atmos_file, atmos)
    if (atmos%has_feh .and. FeH < -99) then
      FeH = atmos%feh
      write(stdout, '(A, F6.2)') 'FeH was overwritten by model: [Fe/H] =', FeH
    endif

    if (vmic >= 0) then
      ia%vmic = vmic
      write(stdout, '(A, F6.2)') 'Vmic was overwritten by namelist: Vmic =', vmic
    endif

    use_ne = use_ne .and. atmos%has_pe
    use_rho = use_rho .and. atmos%has_rho
  END SUBROUTINE read_text_1d

!===============================================================================
!> Reads the actual text files and stores variables under ia
!===============================================================================
  SUBROUTINE read_multi_1d(atmos_file, vmic)
    character(len=*) :: atmos_file
    real(4), intent(in) :: vmic
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_multi_1d_file, atmos_file, atmos)
    if (vmic >= 0) then
      ia%vmic = vmic
      write(stdout, '(A, F6.2)') 'Vmic was overwritten by namelist: Vmic =', vmic
    endif

    use_ne = .true.
    use_rho = .false.
    ia%lrho = ieee_value( ia%lrho, ieee_signaling_nan )
  END SUBROUTINE read_multi_1d

  !===============================================================================
  !> Run one atmos_reader_mod parser and move its columns into ia, stopping on a
  !> parse failure the way the readers always have.
  !===============================================================================
  SUBROUTINE read_parsed(reader, atmos_file, atmos)
    interface
      subroutine reader(atmos_file, atmos, ierr, msg)
        import :: atmos_data_t
        character(len=*), intent(in) :: atmos_file
        type(atmos_data_t), intent(out) :: atmos
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg
      end subroutine reader
    end interface
    character(len=*), intent(in) :: atmos_file
    type(atmos_data_t), intent(out) :: atmos
    integer :: ierr
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call reader(atmos_file, atmos, ierr, msg)
    if (ierr /= 0) then
      write(stdout,*) 'Error: ', trim(msg)
      error stop 'atmosphere read failed'
    endif
    call ia%allocate(atmos%nx, atmos%ny, atmos%nz)
    ia%zz    = atmos%zz
    ia%vz    = atmos%vz
    ia%vmic  = atmos%vmic
    ia%ltemp = atmos%ltemp
    ia%lpe   = atmos%lpe
    ia%lrho  = atmos%lrho
  END SUBROUTINE read_parsed

  !===============================================================================
  SUBROUTINE read_siu_bin(atmos_file, vmic)
    character(len=*) :: atmos_file
    real(4), intent(in) :: vmic
    integer :: lunit, nz, nel, nmol, nbreak
    integer :: i, j, k, recl
    character(128) :: text
    real(4), parameter :: ln10 = log(10.0)
    integer, dimension(:), allocatable :: at_idx, at_nion, break
    real(4), dimension(:), allocatable :: kappa500, pgas
    !-----------------------------------------------------------------------------
    call os%check_file(atmos_file)
    open(file=atmos_file, newunit=lunit, form='unformatted', access='stream')

    read(lunit) recl ! skip the record length in stream access
    read(lunit) nz, nel, nmol, nbreak

    call ia%allocate(1,1,nz)
    allocate(kappa500(nz), pgas(nz))
    allocate(at_idx(nel), at_nion(nel))
    allocate(break(nbreak))

    read(lunit) at_idx   ! Atomic number of the elements
    read(lunit) at_nion  ! Number of ionisation stages considered for each element
    read(lunit) break    ! break points in tau scale for MARCS
    read(lunit) ia%zz    ! log10(tau500)
    read(lunit) ia%ltemp ! log10(temperature)
    read(lunit) kappa500 ! log10(kappa500)
    read(lunit) ia%lpe   ! log10(electron pressure)
    read(lunit) pgas     ! p gas MARCS
    read(lunit) pgas     ! p gas MARCS (not sure why this appears twice)
    read(lunit) ia%lrho  ! log10(density)
    ! the data goes on, but we have everything we need
    close(lunit)

    ! convert from log10() to ln
    ia%ltemp = ia%ltemp * ln10
    ia%lpe   = ia%lpe   * ln10
    ia%lrho  = ia%lrho  * ln10
    if (vmic >= 0) then
      ia%vmic  = vmic
    endif
  END SUBROUTINE read_siu_bin

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
  !> Announce one more patch that will interpolate from the snapshot.
  !===============================================================================
  SUBROUTINE acquire_input_atmos(self)
    class(input_atmos_t), intent(inout) :: self
    !-----------------------------------------------------------------------------
    !$omp atomic
    self%n_users = self%n_users + 1
  END SUBROUTINE acquire_input_atmos

  !===============================================================================
  !> Give up this patch's claim on the snapshot; the last one out frees it.  When
  !> no snapshot was read at all (atmos_format=='save'), deallocate is a no-op.
  !===============================================================================
  SUBROUTINE release_input_atmos(self)
    class(input_atmos_t), intent(inout) :: self
    logical :: last
    !-----------------------------------------------------------------------------
    !$omp critical (input_atmos_cr)
    self%n_users = self%n_users - 1
    last = self%n_users == 0
    if (last) call self%deallocate()
    !$omp end critical (input_atmos_cr)
  END SUBROUTINE release_input_atmos

  !===============================================================================
  SUBROUTINE scale_box(self, box, scale)
    class(input_atmos_t), intent(inout) :: self
    real(8), intent(out) :: box(3)
    real(8), intent(in) :: scale
    real(8) :: ds(3)
    integer :: n(3)
    !-----------------------------------------------------------------------------
    n = [self%nx, self%ny, self%nz-2*self%gz]
    self%xx = self%xx/scale
    self%yy = self%yy/scale
    self%zz = self%zz/scale
    box(1) = self%xx(self%nx)
    box(2) = self%yy(self%ny)
    box(3) = self%zz(self%nz-self%gz)
    box(1:2) = box(1:2) * n(1:2) / max(n(1:2)-1, 1)
    if (self%nx == 1) box(1) = 0.01
    if (self%ny == 1) box(2) = 0.01
  END SUBROUTINE scale_box

!===============================================================================
! This makes sure that the atmosphere is flipped properly
  !===============================================================================
  SUBROUTINE thick2thin(self)
    class(input_atmos_t), intent(inout) :: self
    type(composition_type) :: comp
    real(8), dimension(:,:,:), allocatable :: temp8, ne8, rho8
    real(8), dimension(:), allocatable :: chi, tau
    real(8) :: dz, v1, vn
    integer :: i, j, nx, ny, nz
    !-----------------------------------------------------------------------------
    nx = self%nx
    ny = self%ny
    nz = self%nz
    if (use_rho) then
      v1 = self%lrho(1,1,1)
      vn = self%lrho(1,1,nz)
    else
      v1 = self%lpe(1,1,1)
      vn = self%lpe(1,1,nz)
    endif

    if (v1 < vn) then
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

    if (one_dimensional .and. abs(ia%zz(1) - ia%zz(nz)) < 1e4) then
      ! log tau scale
      allocate(temp8(nx, ny, nz))
      allocate(rho8 (nx, ny, nz))
      allocate(ne8  (nx, ny, nz))
      allocate(chi(nz))
      allocate(tau(nz))
      temp8 = exp(dble(self%ltemp))
      if (use_rho) rho8  = exp(dble(self%lrho))
      if (use_ne)  ne8   = exp(dble(self%lpe)) / (BK * temp8)
      comp = composition
      call comp%reassign_pointer()
      call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne)
      chi=comp%chi500(1,1,:)
      tau=10**(dble(self%zz))

      self%zz(nz) = 0
      do i=nz-1, 1, -1
        dz = (tau(i) - tau(i+1)) * 2 / (chi(i) + chi(i+1))
        self%zz(i) = self%zz(i+1) + dz
      end do
    endif

    if (self%zz(1) > self%zz(nz)) then
      self%zz = -self%zz
    endif

    if (rstar > 2) then
      ! for spherical geometries
      rorigin = self%zz(self%gz+1) + rstar
      rstar = rstar + self%zz(nz)
    endif

    self%xx = self%xx - self%xx(1)         ! let the origin be zero
    self%yy = self%yy - self%yy(1)         ! let the origin be zero
    self%zz = self%zz - self%zz(self%gz+1) ! let the origin be the layer above ghost zone
  END SUBROUTINE thick2thin

  !===============================================================================
  SUBROUTINE setup_amr_mesh(self, amr, nsub, use_ne)
    class(input_atmos_t), intent(inout), target :: self
    character(len=*), intent(in) :: amr
    integer, intent(in) :: nsub
    logical, intent(in) :: use_ne
    integer :: i, j, k, nx, ny, nz
    real(4), dimension(nsub+1) :: sub_edge, var_edge
    real(4), dimension(:), allocatable :: amr_ref, tmp
    real(4), dimension(:), pointer, contiguous :: zz
    real(4), dimension(:,:,:), pointer, contiguous :: lpe
    real(4) :: vmin, vmax
    integer :: zs, ze
    !-----------------------------------------------------------------------------
    nx = self%nx
    ny = self%ny
    nz = self%nz - 2*self%gz
    allocate(self%sub_pos(nsub), self%sub_size(nsub))

    if (uppercase(amr) == 'Z') then
      var_edge(1) = 0
      do i=1, nsub
        call unequal_intervals(nz, i, nsub, zs, ze)
        var_edge(i+1) = ze
      end do
      call interp1d_sorted([(real(i-self%gz), i=1, self%nz)], self%zz, var_edge, sub_edge)
      self%sub_size = sub_edge(2:nsub+1) - sub_edge(1:nsub)
      self%sub_pos  = sub_edge(1:nsub) + 0.5*self%sub_size
      return
    endif

    allocate(amr_ref(nz), source=0.0)
    allocate(tmp(nz), source=0.0)
    zz => self%zz(1+self%gz:nz+self%gz)

    select case (uppercase(amr))
    case ('P')
      if (.not. use_ne) stop 'AMR based on pressure is not supported when use_ne is false'
      lpe => self%lpe(:,:,1+self%gz:nz+self%gz)
    case ('T')
      lpe => self%ltemp(:,:,1+self%gz:nz+self%gz)
    end select

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
    if (uppercase(amr) == 'T') then
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
  !> Lower the top of the box to the highest plane whose columns all carry an
  !> electron pressure above min_k <ln Pe>(k) -- the topmost plane's own
  !> horizontal mean. Above that height the spread of ln Pe across a plane grows
  !> much faster than the mean falls, so a large fraction of the columns sit
  !> orders of magnitude below it; those cells have no counterpart in any 1D
  !> model and their statistical equilibrium is radiatively controlled with
  !> essentially no collisional damping. prune_tol relaxes "all" to "all but this
  !> fraction".
  !>
  !> This lowers self%nz, and must therefore run after thick2thin (which decides
  !> which end is the top) but BEFORE scale_box, which takes the box height from
  !> zz(nz). The snapshot arrays are left at their full length; only the box
  !> shrinks, so the interpolation simply stops sampling the planes above it.
  !===============================================================================
  SUBROUTINE prune_top_planes (self)
    class(input_atmos_t), intent(inout), target :: self
    real(4), dimension(:,:,:), pointer, contiguous :: lpe
    real(4) :: thr
    integer :: k, ktop, nxy, nbad, nz
    !-----------------------------------------------------------------------------
    if (.not. prune_top) return
    if (self%nx*self%ny <= 1) return
    if (.not. use_ne) then
      write(stdout,*) 'WARNING: prune_top needs use_ne=T; no planes removed'
      return
    endif
    nz = self%nz - 2*self%gz
    lpe => self%lpe(:,:,1+self%gz:nz+self%gz)
    nxy = self%nx * self%ny
    ! thick2thin has already put the deepest plane at index 1, so the mean is
    ! smallest at the top and this is just <ln Pe> of the outermost plane
    thr = huge(1.0)
    do k=1, nz
      thr = min(thr, sum(lpe(:,:,k)) / nxy)
    end do

    ktop = nz
    do k=nz, 1, -1
      nbad = count(lpe(:,:,k) <= thr)
      if (nbad <= nint(prune_tol*nxy)) then
        ktop = k
        exit
      endif
    end do

    if (ktop < nz) then
      write(stdout,'(A,I0,A,I0,A,F6.2,A)') &
        ' prune_top: dropping ', nz-ktop, ' of ', nz, ' planes (', &
        100.0*(nz-ktop)/nz, '% of the box height)'
      write(stdout,'(A,ES10.3,A,ES10.3)') &
        '            min ln Pe at the new top ', minval(lpe(:,:,ktop)), &
        ' vs ', minval(lpe(:,:,nz))
      self%nz = ktop + 2*self%gz
    else
      write(stdout,*) 'prune_top: no planes removed'
    endif
  END SUBROUTINE prune_top_planes

  !===============================================================================
  !=============================================================================
  !> Interpolate the input snapshot onto a patch mesh, in the units and variables
  !> the composition expects: temperature, electron density, mass density, the
  !> three velocity components in units of c, and the microturbulence.  The
  !> electron pressure that the readers store in ia%lpe is a scratch quantity
  !> here -- it is interpolated, then inverted to n_e and dropped.
  !>
  !> The optional gradient arguments return the partial derivatives of ln(T)
  !> and of Pe (linear, not log -- the hermite solvers assume chi/Pe linear
  !> along the ray) evaluated on the PCHIP interpolant of the snapshot, per
  !> code length unit; for a rotated call (vec present) the z-derivative is the
  !> along-ray derivative parameterized by z, and the horizontal partials are
  !> unavailable.  With use_ne=.false. there is no Pe field, and the Pe
  !> gradient outputs are zeroed -- the caller is expected to rebuild them from
  !> Pe = n_e*k_B*T after the thermal equilibrium is solved.
  !=============================================================================
  SUBROUTINE interpolate_atmos(n, sub, pos, temp, ne, rho, vx, vy, vz, vmic, vec, bot2pivot, &
                               dlnt_dx, dlnt_dy, dlnt_dz, dpe_dx, dpe_dy, dpe_dz)
    integer, intent(in) :: n(3)
    real(8), intent(in) :: sub(3), pos(3)
    real(4), intent(out), dimension(:,:,:) :: temp, ne, rho, vx, vy, vz, vmic
    real(4), intent(out), dimension(:,:,:), optional :: dlnt_dx, dlnt_dy, dlnt_dz
    real(4), intent(out), dimension(:,:,:), optional :: dpe_dx, dpe_dy, dpe_dz
    real(4), dimension(:,:,:), pointer, contiguous :: rho_in, pe_in, temp_in
    real(4), dimension(:,:,:), allocatable :: pe
    real(8), intent(in), optional :: vec(3), bot2pivot
    real(4), dimension(:), allocatable :: newz
    real(8) :: lower(3), ds(3)
    real(8) :: dx, dy, dz
    integer :: i, ijk(3), nx, ny, nz, itimer=0
    logical :: do_pchip, get_grad_t, get_grad_pe
    !-----------------------------------------------------------------------------
    call trace%begin('interpolate_atmos', itimer=itimer)
    nx = n(1); ny = n(2); nz = n(3)
    allocate(pe(nx,ny,nz))
    get_grad_t  = present(dlnt_dz)
    get_grad_pe = present(dpe_dz) .and. use_ne
    if (present(dlnt_dx)) dlnt_dx = 0.
    if (present(dlnt_dy)) dlnt_dy = 0.
    if (present(dlnt_dz)) dlnt_dz = 0.
    if (present(dpe_dx))  dpe_dx  = 0.
    if (present(dpe_dy))  dpe_dy  = 0.
    if (present(dpe_dz))  dpe_dz  = 0.
    ! with use_ne=.false. pe is never interpolated, and we hand composition ne=0,
    ! which is its signal to build a charge-neutrality initial guess for
    ! iterate_thermal_eq rather than trusting the value
    if (.not. use_ne) pe = 0

    ds = sub / n
    lower = pos - 0.5d0*sub
    lower(1:2) = lower(1:2) + ray_offset*ds(1:2)
    do_pchip = pchip_interpolation

    if (log_interpolation) then
      temp_in => ia%ltemp
      if (use_rho) rho_in => ia%lrho
      if (use_ne)  pe_in => ia%lpe
    else
      allocate(temp_in, source=exp(ia%ltemp))
      if (use_rho) allocate(rho_in, source=exp(ia%lrho))
      if (use_ne)  allocate(pe_in, source=exp(ia%lpe))
    endif

    allocate(newz(nz))
    newz = [ (ds(3)*i + lower(3), i=1, nz) ]
    if (make_eos) then
      ! no interpolation
      temp = ia%ltemp
      rho  = ia%lrho
      pe   = ia%lpe
      vmic = ia%vmic
      vx = 0
      vy = 0
      vz = 0
    elseif (one_dimensional) then
      if (.false.) then ! set to true for debugging original atmosphere points
        temp(1,1,:) = temp_in(1,1,1:nz)
        if (use_rho) rho(1,1,:) = rho_in(1,1,1:nz)
        if (use_ne ) pe (1,1,:) = pe_in (1,1,1:nz)
        vmic(1,1,:) = ia%vmic(1,1,1:nz)
      else
        if (use_rho) call pchip_mono4(ia%zz, rho_in(1,1,:), newz, rho (1,1,:), sorted=.true.)
        if (get_grad_pe) then
          call pchip_mono4(ia%zz, pe_in (1,1,:), newz, pe  (1,1,:), sorted=.true., newdydx=dpe_dz(1,1,:))
        elseif (use_ne) then
          call pchip_mono4(ia%zz, pe_in (1,1,:), newz, pe  (1,1,:), sorted=.true.)
        endif
        if (get_grad_t) then
          call pchip_mono4(ia%zz, temp_in(1,1,:), newz, temp(1,1,:), sorted=.true., newdydx=dlnt_dz(1,1,:))
        else
          call pchip_mono4(ia%zz, temp_in(1,1,:), newz, temp(1,1,:), sorted=.true.)
        endif
        call pchip_mono4(ia%zz, ia%vmic(1,1,:), newz, vmic(1,1,:), sorted=.true.)
      endif
      vx = 0
      vy = 0
      vz = 0
    else
      dx = ia%xx(2)
      dy = ia%yy(2)

      if (present(vec)) then
        if (use_rho) call linterpolate_3d_rotated(rho_in,  rho,  dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip)
        if (use_ne ) call linterpolate_3d_rotated(pe_in,   pe,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip, dfdz=dpe_dz)
        call linterpolate_3d_rotated(temp_in, temp, dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip, dfdz=dlnt_dz)
        call linterpolate_3d_rotated(ia%vx,   vx,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip)
        call linterpolate_3d_rotated(ia%vy,   vy,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip)
        call linterpolate_3d_rotated(ia%vz,   vz,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip)
        call linterpolate_3d_rotated(ia%vmic, vmic, dx, dy, ia%zz, lower(1:2), ds(1:2), newz, vec, bot2pivot, pchip=do_pchip)
      else
        if (use_rho) call linterpolate_3d(rho_in,  rho,  dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip)
        if (use_ne ) call linterpolate_3d(pe_in,   pe,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip, &
                                          dfdx=dpe_dx, dfdy=dpe_dy, dfdz=dpe_dz)
        call linterpolate_3d(temp_in, temp, dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip, &
                             dfdx=dlnt_dx, dfdy=dlnt_dy, dfdz=dlnt_dz)
        call linterpolate_3d(ia%vx,   vx,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip)
        call linterpolate_3d(ia%vy,   vy,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip)
        call linterpolate_3d(ia%vz,   vz,   dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip)
        call linterpolate_3d(ia%vmic, vmic, dx, dy, ia%zz, lower(1:2), ds(1:2), newz, pchip=do_pchip)
      endif
    endif

    if (log_interpolation) then
      temp = exp(temp)
      if (use_rho) rho = exp(rho)
      if (use_ne ) pe  = exp(pe)
      ! the interpolant ran on the log fields: dlnT is already what the callers
      ! expect, but the Pe convention is linear -- convert dlnPe -> dPe
      if (get_grad_pe) then
        if (present(dpe_dx)) dpe_dx = dpe_dx * pe
        if (present(dpe_dy)) dpe_dy = dpe_dy * pe
        dpe_dz = dpe_dz * pe
      endif
    else
      ! the interpolant ran on the linear fields, so its Pe derivatives are
      ! already dPe/d; only the T derivatives need the log conversion
      if (get_grad_t) then
        if (present(dlnt_dx)) dlnt_dx = dlnt_dx / temp
        if (present(dlnt_dy)) dlnt_dy = dlnt_dy / temp
        dlnt_dz = dlnt_dz / temp
      endif
      deallocate(temp_in)
      if (use_rho) deallocate(rho_in)
      if (use_ne)  deallocate(pe_in)
    endif

    ! ---- hand over in the composition's own variables and units ---------------
    ne = pe / (cgs%k_b * temp)
    vx = vx / cgs%c
    vy = vy / cgs%c
    vz = vz / cgs%c

    call trace%end(itimer)
  END SUBROUTINE interpolate_atmos

END MODULE read_atmos_mod
