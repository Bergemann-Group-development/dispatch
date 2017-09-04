!> =============================================================================
!> Boundary conditions for a 2-D box
!> =============================================================================
MODULE pboundary_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE scaling_mod
  USE mesh_mod
  USE mhd_mod
  USE parameters_mod
  USE boundary_wavekill_mod
  implicit none
  private

  type, public:: pboundary_t
  contains
    procedure conditions
  end type
  type(pboundary_t), public:: pboundaries

  real(8), parameter :: pi = 3.14159265358979323846d0

CONTAINS

!> =============================================================================
!> Impose boundary conditions, taking staggering into account.
!> Only for use with Cartesian patches.
!> =============================================================================
SUBROUTINE conditions (self, patch, cit)
  class(pboundary_t) :: self
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer :: i1, i2, i3, jay3, jay2
  integer, intent(in):: cit
  real(8) :: r
  real(8) :: sigma0, p0, vkep, rmin, rmax, tau_in, tau_out, rbnd, rkill, tau, &
             dt, skep
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: q
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%conditions')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  q  => patch%mem(:,:,:,:,cit,1)

  sigma0 = 2.0e-3 * experiment_params%mstar / pi / experiment_params%aplanet**2
  rmin = m2%centre_nat + m2%r(m2%li)
  rmax = m2%centre_nat + m2%r(m2%ui)
  tau_in  = 2.0 * pi * sqrt(rmin**3/(scaling%grav * experiment_params%mtotal))
  tau_out = 2.0 * pi * sqrt(rmax**3/(scaling%grav * experiment_params%mtotal))
  dt = patch%dtime

  ! periodic in z (is this necessary?)

  ! radial boundaries; "wave killing" conditions (following §3.2 of de Val-Borro et al.)
  ! in the regions determined by `kill_rmin` and `kill_rmax`, plus reset to initial
  ! values in the ghost zones.
  if (patch%n(2) > 1) then
    do i3=m3%lb,m3%ub
      do i1=m1%lb,m1%ub
        ! "wave killing" for active zones
        do i2=m2%li,m2%ui
          r = m2%centre_nat + patch%mesh(2)%r(i2)
          vkep = sqrt(experiment_params%mstar * scaling%grav / r)
          skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
          p0 = (experiment_params%aspect_rat * vkep)**2 * sigma0 / (patch%gamma - 1.0)
          if ( r <= experiment_params%kill_rmin ) then
            rkill = experiment_params%kill_rmin
            rbnd  = rmin
            tau   = tau_in
          else if ( r >= experiment_params%kill_rmax ) then
            rkill = experiment_params%kill_rmax
            rbnd  = rmax
            tau   = tau_out
          else
            cycle
          end if
          q(i1,i2,i3,patch%i_d)  = q(i1,i2,i3,patch%i_d) &
                                 + dt * wave_killing(q(i1,i2,i3,patch%i_d), sigma0, tau, rkill, rbnd, r)
          q(i1,i2,i3,patch%i_e)  = q(i1,i2,i3,patch%i_e) &
                                 + dt * wave_killing(q(i1,i2,i3,patch%i_e), p0, tau, rkill, rbnd, r)
          q(i1,i2,i3,patch%i_px) = q(i1,i2,i3,patch%i_px) &
                                 + dt * wave_killing(q(i1,i2,i3,patch%i_px), 0.0d0, tau, rkill, rbnd, r)
          q(i1,i2,i3,patch%i_py) = q(i1,i2,i3,patch%i_py) &
                                 + dt * wave_killing(q(i1,i2,i3,patch%i_py), 0.0d0, tau, rkill, rbnd, r) 
          q(i1,i2,i3,patch%i_pz) = q(i1,i2,i3,patch%i_pz) &
                                 + dt * wave_killing(q(i1,i2,i3,patch%i_pz), skep, tau, rkill, rbnd, r)
        end do
        if (m2%llc_nat < experiment_params%roffset ) then 
          ! reset ghost zones to initial values, except for s2, which is anti-symmetric.
          do i2=m2%lb,m2%lo
            r = m2%centre_nat + patch%mesh(2)%r(i2)
            vkep = sqrt(experiment_params%mstar * scaling%grav / r)
            skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
            p0 = (experiment_params%aspect_rat * vkep)**2 * sigma0 / (patch%gamma - 1.0)
            q(i1,i2,i3,patch%i_d)  = sigma0
            q(i1,i2,i3,patch%i_e)  = p0
            q(i1,i2,i3,patch%i_px) = 0.0d0
            jay2 = 2*m2%lo + 1 - i2
            q(i1,i2,i3,patch%i_py) = -q(i1,jay2,i3,patch%i_py)
            q(i1,i2,i3,patch%i_pz) = skep
          end do
        end if
        if ((m2%llc_nat + m2%s) >= (experiment_params%roffset + patch%box(2))) then
          do i2=m2%uo,m2%ub
            r = m2%centre_nat + patch%mesh(2)%r(i2)
            vkep = sqrt(experiment_params%mstar * scaling%grav / r)
            skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
            p0 = (experiment_params%aspect_rat * vkep)**2 * sigma0 / (patch%gamma - 1.0)
            q(i1,i2,i3,patch%i_d)  = sigma0
            q(i1,i2,i3,patch%i_e)  = p0
            q(i1,i2,i3,patch%i_px) = 0.0d0
            jay2 = 2*m2%uo - 1 - i2
            q(i1,i2,i3,patch%i_py) = -q(i1,jay2,i3,patch%i_py)
            q(i1,i2,i3,patch%i_pz) = skep
          end do
        end if
      end do
    end do
  end if

  ! periodic in phi (is this necessary? yes!)
  if (patch%n(3) > 1 .and. patch%size(3) >= patch%box(3)) then
    ! inner phi boundary
    do i3=m3%lb,m3%lo
      jay3 = modulo(i3 - m3%li, m3%n) + m3%li
      do i2=m2%lb,m2%ub
        do i1=m1%lb,m1%ub
          q(i1,i2,i3,:) = q(i1,i2,jay3,:)
        end do
      end do
    end do
    ! outer phi boundary
    do i3=m3%uo,m3%ub
      jay3 = modulo(i3 - m3%li, m3%n) + m3%li
      do i2=m2%lb,m2%ub
        do i1=m1%lb,m1%ub
          q(i1,i2,i3,:) = q(i1,i2,jay3,:)
        end do
      end do
    end do
  end if

  nullify (m1, m2, m3, q)
  call trace_end
END SUBROUTINE conditions

END MODULE pboundary_mod
