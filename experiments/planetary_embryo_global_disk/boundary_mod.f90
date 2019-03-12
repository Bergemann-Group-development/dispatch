!===============================================================================
!> Boundary conditions for a 2-D box
!===============================================================================
MODULE boundary_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE mhd_mod
  USE parameters_mod
  implicit none
  private

  type, public:: boundary_t
  contains
    procedure conditions
  end type
  type(boundary_t), public:: boundary

  real(8), parameter :: pi = 3.14159265358979323846d0

CONTAINS

!===============================================================================
!> Impose boundary conditions, taking staggering into account.
!> Only for use with Cartesian patches.
!===============================================================================
SUBROUTINE conditions (self, patch)
  class(boundary_t) :: self
  class(mhd_t), target :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer :: i1, i2, i3, jay3, jay2
  real(8) :: sigma0, e0, vkep, rmin, rmax, tau_in, tau_out, rbnd, rkill, tau, &
             dt, skep, mtotal, r
  !-----------------------------------------------------------------------------
  call trace_begin('boundary_t%conditions')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  ! FIXME: HACK to prevent non-cylindrical meshes from invoking wave-killing boundaries.
  ! If the mesh centred on the star is ever switched to spherical, then this will
  ! need to be fixed.
  select type (m1)
  type is (Cartesian_mesh)
    print*, "WARNING: Boundary conditions have not been implemented for Cartesian meshes!"
    return
  type is (spherical_mesh)
    print*, "WARNING: Boundary conditions have not been implemented for spherical meshes!"
    return
  end select

  rmin = m2%centre_nat + m2%r(m2%li)
  rmax = m2%centre_nat + m2%r(m2%ui)
  mtotal = experiment_params%mstar + experiment_params%planetesimal%mass
  tau_in  = 2.0 * pi * sqrt(rmin**3/(scaling%grav * mtotal))
  tau_out = 2.0 * pi * sqrt(rmax**3/(scaling%grav * mtotal))
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
          vkep = Keplerian_speed(r)
          sigma0 = MMSN_surface_density(r)
          skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
          e0 = MMSN_temperature(r) * sigma0 / (patch%gamma - 1.0)

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
          patch%d%q  (i1,i2,i3) = patch%d%q(i1,i2,i3) &
                                + dt * wave_kill(patch%d%q(i1,i2,i3), sigma0, tau, rkill, rbnd, r)
          patch%e1%q (i1,i2,i3) = patch%e1%q(i1,i2,i3) &
                                + dt * wave_kill(patch%e1%q(i1,i2,i3), e0, tau, rkill, rbnd, r)
          patch%s%x%q(i1,i2,i3) = patch%s%x%q(i1,i2,i3) &
                                + dt * wave_kill(patch%s%x%q(i1,i2,i3), 0.0d0, tau, rkill, rbnd, r)
          patch%s%y%q(i1,i2,i3) = patch%s%y%q(i1,i2,i3) &
                                + dt * wave_kill(patch%s%y%q(i1,i2,i3), 0.0d0, tau, rkill, rbnd, r) 
          patch%s%z%q(i1,i2,i3) = patch%s%z%q(i1,i2,i3) &
                                + dt * wave_kill(patch%s%z%q(i1,i2,i3), skep, tau, rkill, rbnd, r)
        end do
        ! reset ghost zones to initial values, except for s2, which is anti-symmetric.
        if (patch%btypes%ijb == 2) then
          do i2=m2%lb,m2%lo
            r = m2%centre_nat + patch%mesh(2)%r(i2)
            vkep = Keplerian_speed(r)
            sigma0 = MMSN_surface_density(r)
            skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
            e0 = MMSN_temperature(r) * sigma0 / (patch%gamma - 1.0)
            patch%d%q  (i1,i2,i3) = sigma0
            patch%e1%q (i1,i2,i3) = e0
            patch%s%x%q(i1,i2,i3) = 0.0d0
            jay2 = 2*m2%lo + 1 - i2
            patch%s%y%q(i1,i2,i3) = -patch%s%y%q(i1,jay2,i3)
            patch%s%z%q(i1,i2,i3) = skep
          end do
        end if
        if (patch%btypes%ojb == 2) then
          do i2=m2%uo,m2%ub
            r = m2%centre_nat + patch%mesh(2)%r(i2)
            vkep = Keplerian_speed(r)
            sigma0 = MMSN_surface_density(r)
            skep = sigma0 * vkep * m1%h31c(i1) * m2%h32c(i2)
            e0 = MMSN_temperature(r) * sigma0 / (patch%gamma - 1.0)
            patch%d%q  (i1,i2,i3) = sigma0
            patch%e1%q (i1,i2,i3) = e0
            patch%s%x%q(i1,i2,i3) = 0.0d0
            jay2 = 2*m2%uo - 1 - i2
            patch%s%y%q(i1,i2,i3) = -patch%s%y%q(i1,jay2,i3)
            patch%s%z%q(i1,i2,i3) = skep
          end do
        end if
      end do
    end do
  end if

  ! periodic in phi (is this necessary? yes!)
  if (patch%n(3) > 1) then
    ! inner phi boundary
    do i3=m3%lb,m3%lo
      jay3 = modulo(i3 - m3%li, m3%n) + m3%li
      do i2=m2%lb,m2%ub
        do i1=m1%lb,m1%ub
          patch%mem(i1,i2,i3,:,patch%it,1) = patch%mem(i1,i2,jay3,:,patch%it,1)
        end do
      end do
    end do
    ! outer phi boundary
    do i3=m3%uo,m3%ub
      jay3 = modulo(i3 - m3%li, m3%n) + m3%li
      do i2=m2%lb,m2%ub
        do i1=m1%lb,m1%ub
          patch%mem(i1,i2,i3,:,patch%it,1) = patch%mem(i1,i2,jay3,:,patch%it,1)
        end do
      end do
    end do
  end if

  nullify (m1, m2, m3)
  call trace_end
END SUBROUTINE conditions

!===============================================================================
FUNCTION wave_kill (x, x0, tau, rkill, rbnd, r) RESULT (rate)
  real(8), intent(in) :: x, x0, tau, rkill, rbnd, r
  real(8) :: rate, quad
  !-----------------------------------------------------------------------------
  quad = ((r - rkill) / (rbnd - rkill))**2
  rate = -(x - x0) * quad / tau

END FUNCTION wave_kill

END MODULE boundary_mod
