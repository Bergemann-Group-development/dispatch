!===============================================================================
!> Initializiation for the shocked cylinder experiment
!===============================================================================
MODULE ic_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE mhd_mod
  USE kinds_mod
  implicit none
  private

  type, public :: ic_t
    ! initial shock parameters
    real(4) :: x0, rhol, rhor, presl, presr, vxl, vxr, vyl, vyr, vzl, vzr
  contains
    procedure :: init => init_initial
    procedure :: init_mhd_vars
  end type
  type(ic_t), public:: ics

  real(4), save :: x0, rhol, rhor, presl, presr, vxl, vxr, vyl, vyr, vzl, vzr

CONTAINS
!===============================================================================
!===============================================================================
SUBROUTINE init_initial(self)
  class(ic_t) :: self
  logical, save :: first_time = .true.
  namelist / IC_params / x0, rhol, rhor, presl, presr, vxl, vxr, vyl, vyr, vzl, vzr
  !----------------------------------------------------------------------------
  call trace_begin('ic_t%init')

  if (first_time) then
    ! default shock front parameters;
    ! downwind conditions given by the shock-jump conditions for M = 1.5
    x0 = -1.5
    rhol = 2.6069
    rhor = 1.4
    presl = 2.4583
    presr = 1.0
    vxl = 0.6944
    vxr = 0.0
    vyl = 0.0
    vyr = 0.0
    vzl = 0.0
    vzr = 0.0

    first_time = .false.
    ! read initial discontinuity parameters
    rewind (io%input); read (io%input, initial_params)
    if (io%master) write (*, initial_params)
  end if

  self%x0 = x0
  self%rhol = rhol
  self%rhor = rhor
  self%presl = presl
  self%presr = presr
  self%vxl = vxl
  self%vxr = vxr
  self%vyl = vyl
  self%vyr = vyr
  self%vzl = vzl
  self%vzr = vzr

  call trace_end
END SUBROUTINE init_initial

!===============================================================================
!> Initialise the MHD variables in a patch using the absolute position of
!> a cell.
!===============================================================================
SUBROUTINE init_mhd_vars(self, patch)
  class(ic_t) :: self
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz
  real(4) :: x, y, z, r, phi
  real(4), dimension(3) :: vec, vcyll, vcylr
  real(4) :: zeta_offset, gamma
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: q
  !-----------------------------------------------------------------------------
  call trace_begin('ic_t%init_mhd_vars')

  gamma = patch%gamma
  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)
  q  => patch%mem(:,:,:,:,patch%it,1)

  if (patch%mesh(1)%type == 'Cartesian_mesh') then
    ! Cartesian initialisation
    do iz=mz%lb,mz%ub
      z = mz%centre_nat + patch%mesh(3)%r(iz)
      do iy=my%lb,my%ub
        y = my%centre_nat + patch%mesh(2)%r(iy)
        do ix=mx%lb,mx%ub
          x = mx%centre_nat + patch%mesh(1)%r(ix)
          if (x <= x0) then
            q(ix,iy,iz,patch%idx%d) = rhol
            q(ix,iy,iz,patch%idx%e) = presl / (patch%gamma - 1.0)
            q(ix,iy,iz,patch%idx%py) = vyl * rhol
            q(ix,iy,iz,patch%idx%pz) = vzl * rhol
          else
            q(ix,iy,iz,patch%idx%d) = rhor
            q(ix,iy,iz,patch%idx%e) = presr / (patch%gamma - 1.0)
            q(ix,iy,iz,patch%idx%px) = vxr * rhor
            q(ix,iy,iz,patch%idx%py) = vyr * rhor
            q(ix,iy,iz,patch%idx%pz) = vzr * rhor
          end if
          ! because the velocity is staggered, special treatment is required
          ! if a value is co-spatial with the initial discontinuity.
          if ( (x-0.5*mx%d) < x0 ) then
            q(ix,iy,iz,patch%idx%px) = vxl * rhol
          else if ( (x-0.5*mx%d) >= x0 .and. (x-1.5*mx%d) < x0) then
            q(ix,iy,iz,patch%idx%px) = 0.25 * (vxl + vxr) * (rhol + rhor)
          else if ( (x-0.5*mx%d) > x0 ) then
            q(ix,iy,iz,patch%idx%px) = vxr * rhor
          end if
          ! in the case that a cell face is co-spatial with the initial discontinuity.
          if ( abs(x - x0 - patch%s%h(1)) < 1.0e-10 ) &
            q(ix,iy,iz,patch%idx%px) = 0.5 * (vxl + vxr) * 0.5 * (rhor + rhol)
          q(ix,iy,iz,patch%idx%bx) = 0.0
          q(ix,iy,iz,patch%idx%by) = 0.0
          q(ix,iy,iz,patch%idx%bz) = 0.0
        end do
      end do
    end do
  else if (patch%mesh(1)%type == 'cylindrical_mesh') then
    select type (mz)
    type is (cylindrical_mesh)
      zeta_offset = mz%zeta
    end select
    ! cylindrical initialisation
    do iz=mz%lb,mz%ub
      phi = mz%centre_nat + patch%mesh(3)%r(iz)
      vcyll(1) =  vzl
      vcyll(2) =  vxl * cos(phi) + vyl * sin(phi)
      vcyll(3) = -vxl * sin(phi) + vyl * cos(phi)
      vcylr(1) =  vzr
      vcylr(2) =  vxr * cos(phi) + vyr * sin(phi)
      vcylr(3) = -vxr * sin(phi) + vyr * cos(phi)
      do iy=my%lb,my%ub
        r = my%centre_nat + patch%mesh(2)%r(iy)
        do ix=mx%lb,mx%ub
          z = mx%centre_nat + patch%mesh(1)%r(ix)
          ! determine cell position in Cartesian coords.
          vec = CylindricalToCartesian(z,r,phi+zeta_offset) + patch%position
          if (vec(1) <= x0) then
            q(ix,iy,iz,patch%idx%d) = rhol
            q(ix,iy,iz,patch%idx%e) = presl / (patch%gamma - 1.0)
            q(ix,iy,iz,patch%idx%px) = vcyll(1) * rhol
            q(ix,iy,iz,patch%idx%py) = vcyll(2) * rhol
            q(ix,iy,iz,patch%idx%pz) = vcyll(3) * rhol
          else
            q(ix,iy,iz,patch%idx%d) = rhor
            q(ix,iy,iz,patch%idx%e) = presr / (patch%gamma - 1.0)
            q(ix,iy,iz,patch%idx%px) = vcylr(1) * rhor
            q(ix,iy,iz,patch%idx%py) = vcylr(2) * rhor
            q(ix,iy,iz,patch%idx%pz) = vcylr(3) * rhor
          end if
          ! the ZEUS module stores the angular momentum rather than the linear momentum
          q(ix,iy,iz,patch%idx%py) = q(ix,iy,iz,patch%idx%py) * mx%h2c(ix)
          q(ix,iy,iz,patch%idx%pz) = q(ix,iy,iz,patch%idx%pz) * mx%h31c(ix) * my%h32c(iy)
          q(ix,iy,iz,patch%idx%bx) = 0.0
          q(ix,iy,iz,patch%idx%by) = 0.0
          q(ix,iy,iz,patch%idx%bz) = 0.0
        end do
      end do
    end do
  end if

  nullify(mx, my, mz)
  call trace_end
END SUBROUTINE init_mhd_vars

END MODULE ic_mod
