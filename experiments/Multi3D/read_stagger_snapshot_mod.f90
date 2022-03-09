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
SUBROUTINE read_stagger
  logical, save :: first_time=.true.
  character(len=128), save :: stagger_file='/lustre/astro/rhoppe/multi-3d/input/atmos/t5777g44m00/t5777g44m0005'
  character(len=128), save :: mesh_file, snap_file, data_file, aux_file
  integer, save :: snap=20, ix=1, iy=1, iz=1
  real(4), save :: zcut=1.5 ! in Mm below the surface
  integer :: iostat, lunit, recl
  integer :: nx, ny, nz, nnz, m(3)
  real(4), pointer, contiguous, dimension(:) :: xm, xx, x
  real(4), pointer, contiguous, dimension(:) :: ym, yy, y
  real(4), pointer, contiguous, dimension(:) :: zm, zz, z
  real(4), pointer, contiguous, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
  real(4), pointer, contiguous, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
  real(4), pointer, contiguous, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
  real(4), pointer, contiguous, dimension(:,:,:) :: temp, rho, px, py, pz, ee, lne
  real(4), pointer, contiguous, dimension(:,:,:) :: eos_lne, eos_tt
  real(4), pointer, contiguous, dimension(:,:,:) :: irho, ipx, ipy, ipz, iie, itemp, ilne
  real(4), pointer, contiguous, dimension(:,:,:) :: ppx, ppy, ppz
  namelist /stagger_params/ stagger_file, snap, zcut, ix, iy, iz
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
  call os%allocate(zz, [nnz])
  allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=xx)
  allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=yy)
  allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=zz)

  read (unit=lunit, iostat=iostat) dx, dxdn, xx, xdn, dxidxup, dxidxdn
  read (unit=lunit, iostat=iostat) dz, dzdn, zz, zdn, dzidzup, dzidzdn ! y is z in stagger
  read (unit=lunit, iostat=iostat) dy, dydn, yy, ydn, dyidyup, dyidydn
  ! close mesh file
  close(unit=lunit)
  !=============================================================================

  m = [nx, ny, nnz]
  call os%allocate(itemp, m)
  call os%allocate(irho, m)
  call os%allocate(ipx, m)
  call os%allocate(ipy, m)
  call os%allocate(ipz, m)
  call os%allocate(iie, m)
  call os%allocate(ilne, m)

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
  read (unit=lunit, iostat=iostat, rec=6) itemp
  ! close data file
  close(unit=lunit)

  !=============================================================================
  ! Reading aux file
  !-----------------------------------------------------------------------------
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
  !=============================================================================

  ! "un"-staggering vectors relative to scalars
  call os%allocate(ppx, m)
  call os%allocate(ppy, m)
  call os%allocate(ppz, m)
  call xup_set(ipx, ppx, nx, ny, nnz)
  call yup_set(ipy, ppy, nx, ny, nnz)
  call zup_set(ipz, ppz, nx, ny, nnz)

  ! transposing data (y to z)
  call transpose_stagger(irho)
  call transpose_stagger(ppx)
  call transpose_stagger(ppy)
  call transpose_stagger(ppz)
  call transpose_stagger(iie)
  call transpose_stagger(itemp)
  call transpose_stagger(ilne)

  ! Index the ghost zones to be ommited 
  nz  = nnz - (2 * gz)
    
  m = [nx, ny, nz]
  call os%allocate(x, [nx])
  call os%allocate(y, [ny])
  call os%allocate(z, [nz])
  call os%allocate(temp, m)
  call os%allocate(rho, m)
  call os%allocate(px, m)
  call os%allocate(py, m)
  call os%allocate(pz, m)
  call os%allocate(ee, m)
  call os%allocate(lne, m)

  ! Converting to cgs from internal units
  ! and flipping z-axis
  x = xx ! * uul
  y = yy ! * uul
  z = zz(nnz-gz:gz+1:-1) ! * uul

  rho  = irho  (:,:, nnz-gz:gz+1:-1 )! * uur
  ee   = iie   (:,:, nnz-gz:gz+1:-1 ) / rho
  px   = ppx   (:,:, nnz-gz:gz+1:-1 )! * uum
  py   = ppy   (:,:, nnz-gz:gz+1:-1 )! * uum
  pz   = ppz   (:,:, nnz-gz:gz+1:-1 )! * uum
  temp = itemp (:,:, nnz-gz:gz+1:-1 )
  lne  = ilne  (:,:, nnz-gz:gz+1:-1 )
  !ne   = exp(ilne(:,:, nnz-gz:gz+1:-1 ) - 3.0*alog(uul)) ! 1 / length [1/cm^3]

  print *, 'temp min = ', minval(temp), 'temp max = ', maxval(temp)
  print *, 'rho  min = ', minval(rho ), 'rho  max = ', maxval(rho )
  print *, 'ee   min = ', minval(ee  ), 'ee   max = ', maxval(ee  )

  allocate(eos_lne, eos_tt, mold=temp(ix:ix, iy:iy, iz:iz))
  !call eos%lookup (d=rho(ix:ix, iy:iy, iz:iz), ee=ee(ix:ix, iy:iy, iz:iz), &
  !                 tt=eos_tt, lne=eos_lne, id=1, label='read_stagger')

  call eos%lookup (d=rho(ix:ix, iy:iy, iz:iz), ee=ee(ix:ix, iy:iy, iz:iz), &
  tt=eos_tt, id=1, label='read_stagger')

  print *, nx, ny, nz
  print *, temp(ix,iy,iz), eos_tt
  !print *, lne(ix,iy,iz), eos_lne
  stop

  call trace%end()
END SUBROUTINE read_stagger

!===============================================================================
SUBROUTINE xup_set (f, fp, nx, ny, nz)
  ! f is centered on (i-.5,j,k), fifth order stagger
  integer, intent(in) :: nx, ny, nz
  real(4), dimension(nx,ny,nz), intent(in) :: f
  real(4), dimension(nx,ny,nz), intent(out) :: fp
  integer :: i, j, k
  integer :: ix, iy, iz
  !--------------------------------------------------------------------------
  ix=1; iy=1; iz=1

  if (nz==1) then
    fp = f
  else
    do k=iz, nz
      do j=iy, ny
        fp(nx-2, j, k) = (                          &
            a * ( f(nx-2, j, k) + f(nx-1, j, k) ) + &
            b * ( f(nx-3, j, k) + f(nx  , j, k) ) + &
            c * ( f(nx-4, j, k) + f(ix  , j, k) )   )

        fp(nx-1, j, k) = (                          &
            a * ( f(nx-1, j, k) + f(nx  , j, k) ) + &
            b * ( f(nx-2, j, k) + f(ix  , j, k) ) + &
            c * ( f(nx-3, j, k) + f(ix+1, j, k) )   )

        fp(nx  , j, k) = (                          &
            a * ( f(nx  , j, k) + f(ix  , j, k) ) + &
            b * ( f(nx-1, j, k) + f(ix+1, j, k) ) + &
            c * ( f(nx-2, j, k) + f(ix+2, j, k) )   )

        fp(ix  , j, k) = (                          & 
            a * ( f(ix  , j, k) + f(ix+1, j, k) ) + &
            b * ( f(nx  , j, k) + f(ix+2, j, k) ) + &
            c * ( f(nx-1, j, k) + f(ix+3, j, k) )   )

        fp(ix+1, j, k) = (                          &
            a * ( f(ix+1, j, k) + f(ix+2, j, k) ) + &
            b * ( f(ix  , j, k) + f(ix+3, j, k) ) + &
            c * ( f(nx  , j, k) + f(ix+4, j, k) )   )
      end do

      do j = iy, ny
        do i = ix+2, nx-3
          fp(i  , j, k) = (                           &
                a * ( f(i  , j, k) + f(i+1, j, k) ) + &
                b * ( f(i-1, j, k) + f(i+2, j, k) ) + &
                c * ( f(i-2, j, k) + f(i+3, j, k) )   )
        end do
      end do
    end do
  end if
END SUBROUTINE xup_set
    
!===============================================================================
SUBROUTINE yup_set (f, fp, nx, ny, nz)
  ! f is centered on (i-.5,j,k), fifth order stagger
  integer, intent(in) :: nx, ny, nz
  real(4), dimension(nx,ny,nz), intent(in) :: f
  real(4), dimension(nx,ny,nz), intent(out) :: fp
  integer :: i, j, k
  integer :: ix, iy, iz
  !--------------------------------------------------------------------------
  ix=1; iy=1; iz=1

  if (ny==1) then
    fp = f
  else
    do k=iz, nz
      do i=ix, nx
        fp(i, ny-2, k) = (                          &
            a * ( f(i, ny-2, k) + f(i, ny-1, k) ) + &
            b * ( f(i, ny-3, k) + f(i, ny  , k) ) + &
            c * ( f(i, ny-4, k) + f(i, iy  , k) )   )

        fp(i, ny-1, k) = (                          &
            a * ( f(i, ny-1, k) + f(i, ny  , k) ) + &
            b * ( f(i, ny-2, k) + f(i, iy  , k) ) + &
            c * ( f(i, ny-3, k) + f(i, iy+1, k) )   )

        fp(i, ny  , k) = (                          &
            a * ( f(i, ny  , k) + f(i, iy  , k) ) + &
            b * ( f(i, ny-1, k) + f(i, iy+1, k) ) + &
            c * ( f(i, ny-2, k) + f(i, iy+2, k) )   )

        fp(i, iy  , k) = (                          &
            a * ( f(i, iy  , k) + f(i, iy+1, k) ) + &
            b * ( f(i, ny  , k) + f(i, iy+2, k) ) + &
            c * ( f(i, ny-1, k) + f(i, iy+3, k) )   )

        fp(i, iy+1, k) = (                          &
            a * ( f(i, iy+1, k) + f(i, iy+2, k) ) + &
            b * ( f(i, iy  , k) + f(i, iy+3, k) ) + &
            c * ( f(i, ny  , k) + f(i, iy+4, k) )   )
      end do

      do j = iy+2, ny-3
        do i = ix, nx
          fp(i, j  ,k) = (                            &
              a * ( f(i, j   , k) + f(i, j+1 , k) ) + &
              b * ( f(i, j-1 , k) + f(i, j+2 , k) ) + &
              c * ( f(i, j-2 , k) + f(i, j+3 , k) )   )
        end do
      end do
    end do
  end if
END SUBROUTINE yup_set

!===============================================================================
SUBROUTINE zup_set (f, fp, nx, ny, nz)
  ! f is centered on (i-.5,j,k), fifth order stagger
  integer, intent(in) :: nx, ny, nz
  real(4), dimension(nx,ny,nz), intent(in) :: f
  real(4), dimension(nx,ny,nz), intent(out) :: fp
  integer :: i, j, k
  integer :: ix, iy, iz
  !--------------------------------------------------------------------------
  ix=1; iy=1; iz=1

  if (nz==1) then
    fp = f
  else
    do j = iy, ny
      do i = ix, nx
        fp(i, j, nz-2) = (                          &
            a * ( f(i, j, nz-2) + f(i, j, nz-1) ) + &
            b * ( f(i, j, nz-3) + f(i, j, nz  ) ) + &
            c * ( f(i, j, nz-4) + f(i, j, iz  ) )   )

        fp(i, j, nz-1) = (                          &
            a * ( f(i, j, nz-1) + f(i, j, nz  ) ) + &
            b * ( f(i, j, nz-2) + f(i, j, iz  ) ) + &
            c * ( f(i, j, nz-3) + f(i, j, iz+1) )   )

        fp(i, j, nz  ) = (                          &
            a * ( f(i, j, nz  ) + f(i, j, iz  ) ) + &
            b * ( f(i, j, nz-1) + f(i, j, iz+1) ) + &
            c * ( f(i, j, nz-2) + f(i, j, iz+2) )   )

        fp(i, j, iz  ) = (                          &
            a * ( f(i, j, iz  ) + f(i, j, iz+1) ) + &
            b * ( f(i, j, nz  ) + f(i, j, iz+2) ) + &
            c * ( f(i, j, nz-1) + f(i, j, iz+3) )   )

        fp(i, j, iz+1) = (                          &
            a * ( f(i, j, iz+1) + f(i, j, iz+2) ) + &
            b * ( f(i, j, iz  ) + f(i, j, iz+3) ) + &
            c * ( f(i, j, nz  ) + f(i, j, iz+4) )   )
      end do
    end do

    do k = iz+2, nz-3
      do j = iy, ny
          do i = ix, nx
            fp(i, j, k) = (                             &
                  a * ( f(i, j, k  ) + f(i, j, k+1) ) + &
                  b * ( f(i, j, k-1) + f(i, j, k+2) ) + &
                  c * ( f(i, j, k-2) + f(i, j, k+3) )   )
          end do
      end do
    end do
  end if
END SUBROUTINE zup_set

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
PURE FUNCTION interp1d_sorted (xx, yy, newx) result(newy)
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

end module read_stagger_snapshot_mod