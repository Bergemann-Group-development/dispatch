!===============================================================================
!> The model atmosphere of the run: the snapshot ia read from the model file, and
!> its interpolation onto the patch meshes.
!>
!> The model formats are parsed by atmos_reader_mod (shared with hop_prune.x and
!> libhop); the read routines here apply the namelist overrides and move the result
!> into ia.
!===============================================================================
MODULE read_atmos_mod
  USE, intrinsic :: ieee_arithmetic
  USE io_mod
  USE os_mod
  USE trace_mod
  USE handyfuncs_mod
  USE const_mod
  USE op_math_mod, only: linspace
  USE interpol_mod, only: linterpolate_3d, linterpolate_3d_rotated, pchip_mono4, interp1d_sorted, &
                          heights_from_tau
  USE m3d_patch_mod, only: composition, use_ne, use_rho, make_eos, prune_top, prune_tol
  USE rt_params_mod, only: pivot_z, ray_offset
  USE ray_mod, only: one_dimensional
  USE composition_mod, only: composition_type
  USE atmos_reader_mod, only: atmos_data_t, read_marcs_file, read_text_1d_file, read_multi_1d_file, &
    read_stagger_file, read_stagger2_file, read_multi_snapshot_file, read_co5bold_file, &
    read_siu_bin_file
  USE console_mod
  implicit none
  private

  !===============================================================================
  !> The input snapshot as read, shared by all patches of a rank.
  !===============================================================================
  type, public :: input_atmos_t
    integer :: n(3)       !< number of points (nx, ny, nz)
    integer :: nx         !< number of points in x
    integer :: ny         !< number of points in y
    integer :: nz         !< number of points in z, ghost zones included
    integer :: gz=0       !< ghost zones
    logical :: is_3d=.false.
    !> The snapshot carries an electron pressure in lpe.  The mesh geometry (AMR
    !> slabs, top pruning) may use it even when use_ne=F takes ne from the EOS.
    logical :: has_pe=.false.
    real(4), dimension(:,:,:), allocatable :: ltemp    !< ln T [K]
    real(4), dimension(:,:,:), allocatable :: lpe      !< ln Pe [dyn cm^-2]
    real(4), dimension(:,:,:), allocatable :: lrho     !< ln rho [g cm^-3]
    real(4), dimension(:,:,:), allocatable :: vx       !< velocity in x [cm s^-1]
    real(4), dimension(:,:,:), allocatable :: vy       !< velocity in y [cm s^-1]
    real(4), dimension(:,:,:), allocatable :: vz       !< velocity in z, outward positive [cm s^-1]
    real(4), dimension(:,:,:), allocatable :: vmic     !< microturbulence [km s^-1]
    real(4), dimension(:),     allocatable :: xx       !< x axis
    real(4), dimension(:),     allocatable :: yy       !< y axis
    real(4), dimension(:),     allocatable :: zz       !< z axis, up after thick2thin
    real(8), dimension(:),     allocatable :: sub_pos  !< centre heights of the depth slabs
    real(8), dimension(:),     allocatable :: sub_size !< heights of the depth slabs
    !> The snapshot is shared by every patch on this rank, and each interpolates
    !> from it exactly once.  Patches announce themselves with acquire(), and the
    !> last one to release() frees it.
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

  real(4), public :: timestamp=-1.0                   !< snapshot time of a CO5BOLD model
  logical, public, save :: log_interpolation=.true.   !< interpolate ln T, ln Pe and ln rho
  logical, public, save :: pchip_interpolation=.true. !< monotone cubic interpolation, else linear
  type(input_atmos_t), public, save, target :: ia     !< the snapshot of this rank

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
  !> Reads a CO5BOLD snapshot into ia and keeps its time stamp.
  !===============================================================================
  SUBROUTINE read_co5bold(atmos_file)
    character(len=*), intent(in) :: atmos_file !< model file
    type(atmos_data_t) :: atmos
    integer :: ierr
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call read_co5bold_file(atmos_file, atmos, ierr, msg)
    call adopt(atmos, ierr, msg)
    timestamp = atmos%timestamp
  END SUBROUTINE read_co5bold

  !===============================================================================
  !> Reads snapshot snap of a Stagger model into ia. The .aux file is read whenever
  !> it exists, so that the mesh grading can use its electron pressure; use_ne
  !> leaves true only where it was requested and the snapshot has one.
  !===============================================================================
  SUBROUTINE read_stagger(stagger_file, snap, use_ne)
    character(len=*), intent(in) :: stagger_file !< model file
    integer, intent(in) :: snap                  !< snapshot number
    logical, intent(inout) :: use_ne             !< use the model's electron pressure
    type(atmos_data_t) :: atmos
    integer :: ierr, itimer=0
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call trace%begin('read_stagger_data', itimer=itimer)
    call read_stagger_file(stagger_file, snap, .true., atmos, ierr, msg)
    call adopt(atmos, ierr, msg)
    use_ne = use_ne .and. atmos%has_pe
    call trace%end(itimer)
  END SUBROUTINE read_stagger

  !===============================================================================
  !> Reads a MULTI3D snapshot into ia.
  !===============================================================================
  SUBROUTINE read_multi_snapshot(atmos_file)
    character(len=*) :: atmos_file !< model file
    type(atmos_data_t) :: atmos
    integer :: ierr
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call read_multi_snapshot_file(atmos_file, atmos, ierr, msg)
    call adopt(atmos, ierr, msg)
  END SUBROUTINE read_multi_snapshot

  !===============================================================================
  !> Reads a Stagger2 snapshot into ia.
  !===============================================================================
  SUBROUTINE read_stagger2(atmos_file)
    character(len=*) :: atmos_file !< model file
    type(atmos_data_t) :: atmos
    integer :: ierr
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call read_stagger2_file(atmos_file, atmos, ierr, msg)
    call adopt(atmos, ierr, msg)
  END SUBROUTINE read_stagger2

  !===============================================================================
  !> Reads a MARCS model into ia and returns its abundances; an FeH below -99 is
  !> replaced by the model's [Fe/H].
  !===============================================================================
  SUBROUTINE read_marcs (file, FeH, abund)
    character(len=*), intent(in) :: file                       !< model file
    real(8), intent(out) :: FeH                                !< [Fe/H], replaced if below -99
    real(4), intent(inout), dimension(:), allocatable :: abund !< log eps of Z = 1..92, allocated
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_marcs_file, file, atmos)
    allocate(abund, source=atmos%abund)

    if (FeH < -99) then
      FeH = atmos%feh
      if (console_on(ch_atmos, 1)) write(console_unit,'(A, F6.2)') 'FeH was overwritten by Marcs model: [Fe/H] =', FeH
    endif
  END SUBROUTINE read_marcs

  !===============================================================================
  !> Reads a 1D text model into ia. An FeH below -99 is replaced by the model's
  !> [Fe/H] if it states one, and use_ne and use_rho are cleared when the model has
  !> no electron pressure or density.
  !===============================================================================
  SUBROUTINE read_text_1d(atmos_file, FeH)
    character(len=*) :: atmos_file !< model file
    real(8), intent(out) :: FeH    !< [Fe/H], replaced if below -99
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_text_1d_file, atmos_file, atmos)
    if (atmos%has_feh .and. FeH < -99) then
      FeH = atmos%feh
      if (console_on(ch_atmos, 1)) write(console_unit,'(A, F6.2)') 'FeH was overwritten by model: [Fe/H] =', FeH
    endif

    use_ne = use_ne .and. atmos%has_pe
    use_rho = use_rho .and. atmos%has_rho
  END SUBROUTINE read_text_1d

  !===============================================================================
  !> Reads a MULTI 1D model into ia. Its electron densities are used and its
  !> densities are not, so ia%lrho is set to a signalling NaN.
  !===============================================================================
  SUBROUTINE read_multi_1d(atmos_file)
    character(len=*) :: atmos_file !< model file
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_multi_1d_file, atmos_file, atmos)

    use_ne = .true.
    use_rho = .false.
    ia%lrho = ieee_value( ia%lrho, ieee_signaling_nan )
  END SUBROUTINE read_multi_1d

  !===============================================================================
  !> Runs reader, one atmos_reader_mod parser of a single-file format, and moves the
  !> result into ia.
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
    character(len=*), intent(in) :: atmos_file !< model file
    type(atmos_data_t), intent(out) :: atmos   !< parsed model
    integer :: ierr
    character(len=512) :: msg
    !-----------------------------------------------------------------------------
    call reader(atmos_file, atmos, ierr, msg)
    call adopt(atmos, ierr, msg)
  END SUBROUTINE read_parsed

  !===============================================================================
  !> Stops on a parse failure, else copies the parsed columns and mesh into ia.
  !===============================================================================
  SUBROUTINE adopt(atmos, ierr, msg)
    type(atmos_data_t), intent(in) :: atmos !< parsed model
    integer, intent(in) :: ierr             !< parser status, nonzero on failure
    character(len=*), intent(in) :: msg     !< parser error message
    !-----------------------------------------------------------------------------
    if (ierr /= 0) then
      call console_error (trim(msg))
      error stop 'atmosphere read failed'
    endif
    call ia%allocate(atmos%nx, atmos%ny, atmos%nz)
    ia%gz    = atmos%gz
    ia%xx    = atmos%xx
    ia%yy    = atmos%yy
    ia%zz    = atmos%zz
    ia%vx    = atmos%vx
    ia%vy    = atmos%vy
    ia%vz    = atmos%vz
    ia%vmic  = atmos%vmic
    ia%ltemp = atmos%ltemp
    ia%lpe   = atmos%lpe
    ia%lrho  = atmos%lrho
    ia%has_pe = atmos%has_pe
  END SUBROUTINE adopt

  !===============================================================================
  !> Reads a SIU binary model into ia.
  !===============================================================================
  SUBROUTINE read_siu_bin(atmos_file)
    character(len=*) :: atmos_file !< model file
    type(atmos_data_t) :: atmos
    !-----------------------------------------------------------------------------
    call read_parsed(read_siu_bin_file, atmos_file, atmos)
  END SUBROUTINE read_siu_bin


  !===============================================================================
  !> Allocates the snapshot for nx*ny*nz points; axes and velocities start at zero.
  !===============================================================================
  SUBROUTINE allocate_input_atmos(self, nx, ny, nz)
    class(input_atmos_t), intent(inout) :: self
    integer, intent(in) :: nx !< number of points in x
    integer, intent(in) :: ny !< number of points in y
    integer, intent(in) :: nz !< number of points in z
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
  !> Frees the snapshot fields; the axes stay.
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
  !> Divides the axes by scale and returns the box size: the horizontal extent of a
  !> periodic grid, whose n points span n intervals, and the height above the lower
  !> ghost zone. A direction with a single point gets 0.01.
  !===============================================================================
  SUBROUTINE scale_box(self, box, scale)
    class(input_atmos_t), intent(inout) :: self
    real(8), intent(out) :: box(3) !< box size [scale]
    real(8), intent(in) :: scale   !< length unit [cm]
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
  !> Puts the snapshot bottom first with z pointing up, and the origin at the first
  !> horizontal point and at the lowest plane above the ghost zone. A 1D model on
  !> a log tau scale gets heights from its 500 nm opacity (heights_from_tau).
  !===============================================================================
  SUBROUTINE thick2thin(self)
    class(input_atmos_t), intent(inout) :: self
    type(composition_type) :: comp
    real(8), dimension(:,:,:), allocatable :: temp8, ne8, rho8
    real(8), dimension(:), allocatable :: chi, tau, z
    real(8) :: v1, vn
    integer :: i, j, nx, ny, nz
    logical :: z_down
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

    ! The mesh has z pointing up.  The vz of a 3D model is the component along the
    ! model's own z axis, which changes sign when that axis points into the star,
    ! whichever end of the arrays the model stores first (Stagger2 stores the top
    ! first with z up).  The 1D formats state vz positive outward whatever their
    ! depth scale, so it keeps its sign.
    z_down = (v1 < vn) .eqv. (self%zz(1) < self%zz(nz))
    if (z_down .and. self%is_3d) self%vz = -self%vz

    if (v1 < vn) then
      do j=1, ny
        do i=1, nx
          self%ltemp(i,j,:) =  self%ltemp(i,j,nz:1:-1)
          self%lrho (i,j,:) =  self%lrho (i,j,nz:1:-1)
          self%lpe  (i,j,:) =  self%lpe  (i,j,nz:1:-1)
          self%vx   (i,j,:) =  self%vx   (i,j,nz:1:-1)
          self%vy   (i,j,:) =  self%vy   (i,j,nz:1:-1)
          self%vz   (i,j,:) =  self%vz   (i,j,nz:1:-1)
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
      allocate(tau(nz), z(nz))
      temp8 = exp(dble(self%ltemp))
      if (use_rho) rho8  = exp(dble(self%lrho))
      if (use_ne)  ne8   = exp(dble(self%lpe)) / (BK * temp8)
      comp = composition
      call comp%reassign_pointer()
      call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne)
      chi=comp%chi500(1,1,:)
      ! the column is bottom first here; heights_from_tau takes it top first
      tau = 10**(dble(self%zz(nz:1:-1)))
      call heights_from_tau(tau, chi(nz:1:-1), z)
      self%zz = real(z(nz:1:-1), kind(self%zz))
    endif

    if (self%zz(1) > self%zz(nz)) then
      self%zz = -self%zz
    endif

    self%xx = self%xx - self%xx(1)         ! let the origin be zero
    self%yy = self%yy - self%yy(1)         ! let the origin be zero
    self%zz = self%zz - self%zz(self%gz+1) ! let the origin be the layer above ghost zone
  END SUBROUTINE thick2thin

  !===============================================================================
  !> Splits the depth into nsub slabs for the patches and sets sub_pos and sub_size:
  !> equal numbers of points (amr='Z'), or equal steps in the accumulated change of
  !> the horizontally averaged ln Pe (amr='P'), which makes slabs thin where Pe
  !> changes fast. 'P' falls back to 'Z' for a model without electron pressure.
  !===============================================================================
  SUBROUTINE setup_amr_mesh(self, amr, nsub)
    class(input_atmos_t), intent(inout), target :: self
    character(len=*), intent(in) :: amr !< 'P' or 'Z', in either case
    integer, intent(in) :: nsub         !< number of slabs
    integer :: i, j, k, nx, ny, nz
    real(4), dimension(nsub+1) :: sub_edge, var_edge
    real(4), dimension(:), allocatable :: amr_ref, tmp
    real(4), dimension(:), pointer, contiguous :: zz
    real(4), dimension(:,:,:), pointer, contiguous :: lpe
    real(4) :: vmin, vmax
    integer :: zs, ze
    character :: grading
    !-----------------------------------------------------------------------------
    nx = self%nx
    ny = self%ny
    nz = self%nz - 2*self%gz
    allocate(self%sub_pos(nsub), self%sub_size(nsub))

    select case (uppercase(amr))
    case ('P')
      if (.not. self%has_pe) then
        call console_warn (ch_atmos, "amr='P' needs a model with an electron pressure; using amr='Z'")
        grading = 'Z'
      else
        grading = 'P'
      endif
    case ('Z')
      grading = 'Z'
    case default
      call console_error ("amr='"//trim(amr)//"' is not one of 'P', 'Z' or ''")
      error stop 'invalid amr'
    end select

    if (grading == 'Z') then
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
    lpe => self%lpe(:,:,1+self%gz:nz+self%gz)

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
    var_edge = linspace(vmin, vmax, nsub+1)
    call interp1d_sorted(amr_ref, zz, var_edge, sub_edge)
    !sub_edge = sub_edge * nz / (nz-1)
    self%sub_size = sub_edge(2:nsub+1) - sub_edge(1:nsub)
    self%sub_pos  = sub_edge(1:nsub) + 0.5*self%sub_size
  END SUBROUTINE setup_amr_mesh

  !===============================================================================
  !> Lower the top of the box to the highest plane whose columns all carry an
  !> electron pressure above min_k `<ln Pe>`(k), the topmost plane's own
  !> horizontal mean. Above that height the spread of ln Pe across a plane grows
  !> much faster than the mean falls, so a large fraction of the columns sit
  !> orders of magnitude below it; those cells have no counterpart in any 1D
  !> model and their statistical equilibrium is radiatively controlled with
  !> essentially no collisional damping. prune_tol relaxes "all" to "all but this
  !> fraction".
  !>
  !> This lowers self%nz, and must therefore run after thick2thin (which decides
  !> which end is the top) but before scale_box, which takes the box height from
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
    if (.not. self%has_pe) then
      call console_warn (ch_atmos, 'prune_top needs a snapshot with an electron pressure; no planes removed')
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
      if (console_on(ch_atmos, 1)) then
        write(console_unit,'(A,I0,A,I0,A,F6.2,A)') &
          'prune_top: dropping ', nz-ktop, ' of ', nz, ' planes (', &
          100.0*(nz-ktop)/nz, '% of the box height)'
        write(console_unit,'(A,ES10.3,A,ES10.3)') &
          '           min ln Pe at the new top ', minval(lpe(:,:,ktop)), &
          ' vs ', minval(lpe(:,:,nz))
      end if
      self%nz = ktop + 2*self%gz
    else
      if (console_on(ch_atmos, 1)) write(console_unit,'(a)') 'prune_top: no planes removed'
    endif
  END SUBROUTINE prune_top_planes

  !===============================================================================
  !> Interpolate the input snapshot onto a patch mesh, in the units and variables
  !> the composition expects: temperature, electron density, mass density, the
  !> three velocity components in units of c, and the microturbulence.  The
  !> electron pressure that the readers store in ia%lpe is a scratch quantity
  !> here: it is interpolated, then inverted to n_e and dropped.
  !>
  !> The optional gradient arguments return the partial derivatives of ln(T)
  !> and of Pe (linear, not log, because the hermite solvers assume chi/Pe linear
  !> along the ray) evaluated on the PCHIP interpolant of the snapshot, per
  !> code length unit; for a rotated call (vec present) the z-derivative is the
  !> along-ray derivative parameterized by z, and the horizontal partials are
  !> unavailable.  With use_ne=.false. there is no Pe field, and the Pe
  !> gradient outputs are zeroed; the caller is expected to rebuild them from
  !> Pe = n_e*k_B*T after the thermal equilibrium is solved.
  !===============================================================================
  SUBROUTINE interpolate_atmos(n, sub, pos, temp, ne, rho, vx, vy, vz, vmic, vec, bot2pivot, &
                               dlnt_dx, dlnt_dy, dlnt_dz, dpe_dx, dpe_dy, dpe_dz)
    integer, intent(in) :: n(3)                                 !< patch mesh points
    real(8), intent(in) :: sub(3)                               !< patch size
    real(8), intent(in) :: pos(3)                               !< patch centre
    real(4), intent(out), dimension(:,:,:) :: temp              !< temperature [K]
    real(4), intent(out), dimension(:,:,:) :: ne                !< electron number density [cm^-3]
    real(4), intent(out), dimension(:,:,:) :: rho               !< mass density [g cm^-3]
    real(4), intent(out), dimension(:,:,:) :: vx                !< velocity in x [c]
    real(4), intent(out), dimension(:,:,:) :: vy                !< velocity in y [c]
    real(4), intent(out), dimension(:,:,:) :: vz                !< velocity in z [c]
    real(4), intent(out), dimension(:,:,:) :: vmic              !< microturbulence [km s^-1]
    real(8), intent(in), optional :: vec(3)                     !< ray direction of a rotated patch
    real(8), intent(in), optional :: bot2pivot                  !< height of the pivot plane above the patch bottom
    real(4), intent(out), dimension(:,:,:), optional :: dlnt_dx !< d ln T/dx
    real(4), intent(out), dimension(:,:,:), optional :: dlnt_dy !< d ln T/dy
    real(4), intent(out), dimension(:,:,:), optional :: dlnt_dz !< d ln T/dz
    real(4), intent(out), dimension(:,:,:), optional :: dpe_dx  !< dPe/dx [dyn cm^-2]
    real(4), intent(out), dimension(:,:,:), optional :: dpe_dy  !< dPe/dy [dyn cm^-2]
    real(4), intent(out), dimension(:,:,:), optional :: dpe_dz  !< dPe/dz [dyn cm^-2]
    real(4), dimension(:,:,:), pointer, contiguous :: rho_in, pe_in, temp_in
    real(4), dimension(:,:,:), allocatable :: pe
    real(4), dimension(:), allocatable :: newz
    real(8) :: lower(3), ds(3)
    real(8) :: dx, dy
    integer :: i, nx, ny, nz, itimer=0
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
      if (any(ia%vz /= 0)) then
        call pchip_mono4(ia%zz, ia%vz(1,1,:), newz, vz(1,1,:), sorted=.true.)
      else
        vz = 0
      endif
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
      ! expect, but the Pe convention is linear: convert dlnPe -> dPe
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
    ne = pe / (bk * temp)
    vx = vx / cc
    vy = vy / cc
    vz = vz / cc

    call trace%end(itimer)
  END SUBROUTINE interpolate_atmos

END MODULE read_atmos_mod
