module read_stagger_snapshot_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE eos_mod
  USE scaling_mod
  implicit none

  ! stagger ghost zones
  integer, parameter :: gz=5

  ! stagger coefficients
  real(4), parameter :: uul = 1.0e+8      ! length
  real(4), parameter :: uur = 1.0e-7      ! density
  real(4), parameter :: uut = 1.0e2       ! time
  real(4), parameter :: uuv = uul / uut   ! velocity
  real(4), parameter :: uum = uuv * uur   ! momentum
  real(4), parameter :: uue = uuv**2      ! internal energy per unit mass
  real(4), parameter :: uup = uur * uue   ! pressure
  real(4), parameter :: uua = 1./uul      ! extinction coefficient (alpha)
  real(4), parameter :: c = 3./256.
  real(4), parameter :: b = -25./256.
  real(4), parameter :: a = 0.5-b-c
    
contains

!===============================================================================
SUBROUTINE read_stagger_snapshot(ds, d, px, py, pz, ee, tt, ne)
  real(4), intent(out) :: ds(3)
  real(4), pointer, contiguous, dimension(:,:,:), intent(out) :: d, px, py, pz, ee
  real(4), pointer, contiguous, dimension(:,:,:), intent(out), optional :: tt, ne
  logical, save :: first_time=.true.
  character(len=128), save :: stagger_file='/lustre/astro/rhoppe/multi-3d/input/atmos/t5777g44m00/t5777g44m0005'
  character(len=128), save :: mesh_file, snap_file, data_file, aux_file
  integer, save :: snap=20, inx=5, iny=5, inz=230
  real(4), save :: zcut=1.5, ee_offset=3.2122 ! in Mm below the surface
  integer :: iostat, lunit, recl
  integer :: nx, ny, nz, nnz, m(3)
  real(4), pointer, contiguous, dimension(:) :: xx, yy, zm, zz
  real(4), pointer, contiguous, dimension(:) :: newx, newy, newz
  real(4), pointer, contiguous, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
  real(4), pointer, contiguous, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
  real(4), pointer, contiguous, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
  real(4), pointer, contiguous, dimension(:,:,:) :: irho, ipx, ipy, ipz, iie, itemp, ilne
  real(4), pointer, contiguous, dimension(:,:,:) :: ppx, ppy, ppz
  namelist /stagger_params/ stagger_file, snap, zcut, ee_offset, inx, iny, inz
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_snapshot_mod%read_stagger')

  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, stagger_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, stagger_params)
    end if
    call io%namelist_warning ('stagger_params', iostat)
    mesh_file = trim(stagger_file) // '.msh'
    write(snap_file, '(A,A1,I5.5)') trim(stagger_file), '_', snap
    write(data_file, *) trim(snap_file), '.dat'
    write(aux_file,  *) trim(snap_file), '.aux'
  end if

  !=============================================================================
  ! Reading mesh file
  !-----------------------------------------------------------------------------
  open(newunit=lunit, file=mesh_file, iostat=iostat, convert='big_endian', form='unformatted')

  ! read first entry in mesh file
  read(unit=lunit, iostat=iostat) nx, nnz, ny ! y is z in stagger

  call os%allocate(xx, [nx])
  call os%allocate(yy, [ny])
  call os%allocate(zm, [nnz])
  allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=xx)
  allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=yy)
  allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=zm)

  read (unit=lunit, iostat=iostat) dx, dxdn, xx, xdn, dxidxup, dxidxdn
  read (unit=lunit, iostat=iostat) dz, dzdn, zm, zdn, dzidzup, dzidzdn ! y is z in stagger
  read (unit=lunit, iostat=iostat) dy, dydn, yy, ydn, dyidyup, dyidydn
  ! close mesh file
  close(unit=lunit)
  !=============================================================================

  m = [nx, ny, nnz]
  call os%allocate(irho, m)
  call os%allocate(ipx, m)
  call os%allocate(ipy, m)
  call os%allocate(ipz, m)
  call os%allocate(iie, m)

  recl = kind(1.0)*nx*ny*nnz
  !=============================================================================
  ! Reading data file
  !-----------------------------------------------------------------------------
  open(newunit=lunit, file=data_file, iostat=iostat, convert='big_endian', &
       form='unformatted', status='old', access='direct', recl=recl)

    ! 1: density
    ! 2: momentum per unit volume in direction x
    ! 3: momentum per unit volume in direction y
    ! 4: momentum per unit volume in direction z
    ! 5: internal energy per unit volume
    ! 6: temperature

  read (unit=lunit, iostat=iostat, rec=1) irho
  read (unit=lunit, iostat=iostat, rec=2) ipx
  read (unit=lunit, iostat=iostat, rec=3) ipz ! y is z in stagger
  read (unit=lunit, iostat=iostat, rec=4) ipy ! z is y in stagger
  read (unit=lunit, iostat=iostat, rec=5) iie
  if (present(tt)) then
    call os%allocate(itemp, m)
    read (unit=lunit, iostat=iostat, rec=6) itemp
  endif
  ! close data file
  close(unit=lunit)

  !=============================================================================
  ! Reading aux file
  !-----------------------------------------------------------------------------
  if (present(ne)) then
    call os%allocate(ilne, m)
    open(newunit=lunit, file=aux_file, iostat=iostat, convert='big_endian', &
         form='unformatted', status='old', access='direct', recl=recl)

      ! 1: ln total pressure (gas + rad)
      ! 2: ln Rosseland extinction coeff.
      ! 3: ln Temperature
      ! 4: ln electron number density
      ! 5: ln Planck-averaged extinction coeff.
      ! 6: ln extinction coeff at 500 nm

    read (unit=lunit, iostat=iostat, rec=4) ilne
    ! close aux file
    close(unit=lunit)
  endif
  !=============================================================================

  ! "un"-staggering vectors relative to scalars
  call os%allocate(ppx, m)
  call os%allocate(ppy, m)
  call os%allocate(ppz, m)
  call xup_set(ipx, ppx)

  ! transposing data (y to z)
  call transpose_stagger(irho)
  call transpose_stagger(ppx)
  call transpose_stagger(ipy)
  call transpose_stagger(ipz)
  call transpose_stagger(iie)

  call yup_set(ipy, ppy)
  call zup_set(ipz, ppz)

  ! Index the ghost zones to be ommited 
  nz  = nnz - (2 * gz)

  call os%allocate(zz, [nz])

  zz = -zm(nnz-gz:gz+1:-1) ! flipping z-axis

  call os%allocate(newx, [inx])
  call os%allocate(newy, [iny])
  call os%allocate(newz, [inz])

  newx = linspace(xx(1), xx(nx), inx)
  newy = linspace(yy(1), yy(ny), iny)
  newz = linspace(zz(1), zz(nz), inz)

  m = [inx, iny, inz]
  call os%allocate(d, m)
  call os%allocate(px, m)
  call os%allocate(py, m)
  call os%allocate(pz, m)
  call os%allocate(ee, m)

  call linterpolate_3d(irho (:,:, nnz-gz:gz+1:-1), d,  xx, yy, zz, newx, newy, newz)
  call linterpolate_3d( ppx (:,:, nnz-gz:gz+1:-1), px, xx, yy, zz, newx, newy, newz)
  call linterpolate_3d( ppy (:,:, nnz-gz:gz+1:-1), py, xx, yy, zz, newx, newy, newz)
  call linterpolate_3d(-ppz (:,:, nnz-gz:gz+1:-1), pz, xx, yy, zz, newx, newy, newz)
  call linterpolate_3d( iie (:,:, nnz-gz:gz+1:-1), ee, xx, yy, zz, newx, newy, newz)

  ee = ee / d + ee_offset

  if (present(tt)) then
    call transpose_stagger(itemp)
    call os%allocate(tt, m)
    call linterpolate_3d(itemp(:,:, nnz-gz:gz+1:-1), tt, xx, yy, zz, newx, newy, newz)
  endif

  if (present(ne)) then
    call transpose_stagger(ilne)
    call os%allocate(ne, m)
    ilne = exp(ilne - 3*log(uul))
    call linterpolate_3d(ilne(:,:, nnz-gz:gz+1:-1), ne, xx, yy, zz, newx, newy, newz)
  endif

  call os%deallocate(irho)
  call os%deallocate(ppx)
  call os%deallocate(ppy)
  call os%deallocate(ppz)
  call os%deallocate(iie)
  call os%deallocate(itemp)
  call os%deallocate(ilne)

  ds(1) = (newx(inx) - newx(1)) / (inx - 1)
  ds(2) = (newy(iny) - newy(1)) / (iny - 1)
  ds(3) = (newz(inz) - newz(1)) / (inz - 1)

  call trace%end()
END SUBROUTINE read_stagger_snapshot

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
PURE SUBROUTINE transpose_stagger(f)
  real, dimension(:,:,:), intent(inout)  :: f
  integer :: i
  ! ----------------------------------------------------------------------------
  do i=1, size(f, 1)
    f(i, :, :) = transpose(f(i, :, :))
  end do
END SUBROUTINE transpose_stagger

!===============================================================================
FUNCTION interp1d_sorted (xx, yy, newx) result(newy)
  real(4), dimension(:), intent(in) :: xx, yy, newx
  real(4), dimension(:), allocatable :: newy
  real(4) :: m
  integer :: ii, i, nnew, nold
  !-----------------------------------------------------------------------------
  nnew = size(newx, 1)
  nold = size(xx, 1)
  allocate(newy(nnew))

  i = 1
  do ii=1, nnew
    do while (xx(i+1) < newx(ii)) 
        i = i + 1
        if (i > nold) then
          i = nold
          exit
        endif
    end do

    if (i == nold) then
        newy(ii) = yy(i)
    elseif (xx(1) > newx(ii)) then
        newy(ii) = yy(1)
    else
        m = (yy(i+1) - yy(i)) /  (xx(i+1) - xx(i))
        newy(ii) = m * (newx(ii) - xx(i)) + yy(i)
    endif
  end do
END FUNCTION interp1d_sorted

!===============================================================================
PURE FUNCTION linspace (x1, x2, n) result(xx)
  real(4), intent(in) :: x1, x2
  integer, intent(in) :: n
  real(4), dimension(n) :: xx
  real(4) :: dx
  integer :: i
  !-----------------------------------------------------------------------------
  dx = (x2 - x1) / real(n-1)
  xx = (/ (x1 + i*dx, i=0, n-1) /)
END FUNCTION linspace

!===============================================================================
SUBROUTINE linterpolate_3d(f, fout, oldx, oldy, oldz, newx, newy, newz)
  real(4), dimension(:,:,:), intent(in)  :: f
  real(4), dimension(:,:,:), intent(out) :: fout
  real(4), dimension(:), intent(in) :: oldx, oldy, oldz, newx, newy, newz
  real(4), dimension(:,:,:), allocatable :: fz, fy
  integer :: i,j,k, nx, ny, nz, ox,oy,oz
  !-----------------------------------------------------------------------------
  ox = size(oldx,1); nx = size(newx,1)
  oy = size(oldy,1); ny = size(newy,1)
  oz = size(oldz,1); nz = size(newz,1)
  
  allocate(fz(ox,oy,nz))
  allocate(fy(ox,ny,nz))

  ! z-dimension linear interpolation
  do j = 1, oy
    do i = 1, ox
       fz(i,j,:) = interp1d_sorted(oldz, f(i,j,:), newz)
    end do
  end do

  ! y-dimension linear interpolation
  do k = 1, nz
    do i = 1, ox
       fy(i,:,k) = interp1d_sorted(oldy, fz(i,:,k), newy)
    end do
  end do

  ! x-dimension linear interpolation
  do k = 1, nz
    do j = 1, ny
       fout(:,j,k) = interp1d_sorted(oldx, fy(:,j,k), newx)
    end do
  end do
END SUBROUTINE linterpolate_3d

end module read_stagger_snapshot_mod