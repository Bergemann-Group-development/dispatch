!===============================================================================
!> Formal solvers for rt_patch_t (short/long characteristics; linear, quadratic
!> and hermite), split out of rt_patch_mod as a submodule with full host
!> association.
!===============================================================================
SUBMODULE (rt_patch_mod) rt_solvers_smod
  implicit none
contains

  !===============================================================================
  module subroutine short_char (self, vec)
    class(rt_patch_t) :: self
    type(ray_type) :: ray, amr_ray
    real(8), intent(in) :: vec(3)
    real(4) :: Inuup, xxup, ssup
    real(4) :: dsup1, dsup2, dsup05, xx, ss, Inu, Psi, kx, ky, Irem, dt, E, W, Esc
    real(8) :: ds(3)
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:,:,:), allocatable :: row_red
    integer :: xup_t(4,self%nx), yup_t(4,self%ny)
    integer :: zs, ze, xup, yup, zup, i, ii, jj, upo1, upo2
    logical :: amr_boundary, rebuild_tables, get_rates, split_rates, off
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char', itimer=itimer)
    allocate(work_up(self%fs:self%fe,3))
    allocate(work_min(self%fs:self%fe,3))
    allocate(row_red(self%fs:self%fe, 3, self%nx, self%ny))
    get_rates = self%compute_psi
    split_rates = self%split_rates

    zs = merge(1, self%nz, self%dir > 0)
    ze = merge(self%nz, 1, self%dir > 0)

    call ray%set_stencil(self%ds * scaling%l, vec)
    xkernel1 = real(ray%xkern)
    ykernel1 = real(ray%ykern)
    dsup1 = 0.5*real(ray%dsup)

    amr_boundary = amr/='' .and. self%dir < 0
    if (amr_boundary) then
      ds = self%ds
      if (associated(self%upstream)) ds = self%upstream%ds
      call amr_ray%set_stencil(ds * scaling%l, vec)
      xkernel2 = real(amr_ray%xkern)
      ykernel2 = real(amr_ray%ykern)
      dsup2 = 0.5*real(amr_ray%dsup)
    endif

    rebuild_tables = .true.
    do z=zs, ze, self%dir
      zup = z - self%dir

      if (amr_boundary .and. z==zs) then
        xkern => xkernel2
        ykern => ykernel2
        dsup05 = dsup2
        upo1 = amr_ray%upo1 - 2
        upo2 = amr_ray%upo2 - 2
      else
        xkern => xkernel1
        ykern => ykernel1
        dsup05 = dsup1
        upo1 = ray%upo1 - 2
        upo2 = ray%upo2 - 2
      endif

      if (rebuild_tables) then
        do x=1,self%nx; do ii=1,4
          xup_t(ii,x) = modulo(ii + x + upo1, self%nx) + 1
        end do; end do
        do y=1,self%ny; do jj=1,4
          yup_t(jj,y) = modulo(jj + y + upo2, self%ny) + 1
        end do; end do
        rebuild_tables = amr_boundary .and. z==zs
      endif

      !--- ROW PASS: x-reduce the whole upstream plane ----------------------------
      row_red = 0.0
      do yup = 1, self%ny
        do x = 1, self%nx
          do ii = 1, 4
            xup = xup_t(ii,x)
            kx = xkern(ii)
            do i=1, 3
              do f = self%fs, self%fe
                row_red(f,i,x,yup) = row_red(f,i,x,yup) + kx*self%work(f,xup,yup,zup,i)
              end do
            end do
          end do
        end do
      end do

      do y=1, self%ny
        do x=1, self%nx
          work_up = 0.0
          work_min = huge(1.0)

          ! 4x4 kernel sum
          do jj=1,4
            yup = yup_t(jj,y)
            ky = ykern(jj)
            do i=1, 3
              do f=self%fs,self%fe
                work_up(f,i) = work_up(f,i) + ky*row_red(f,i,x,yup)
              end do
            end do
          end do

          ! 2x2 min over the inner cells
          do jj=2,3
            do ii=2,3
              xup = xup_t(ii,x)
              yup = yup_t(jj,y)
              do i=1, 3
                do f=self%fs,self%fe
                  work_min(f,i) = min(work_min(f,i), self%work(f,xup,yup,zup,i))
                end do
              end do
            end do
          end do

          work_up = max(work_up, work_min)

          ! One f-loop per rate formulation, so that neither carries the other's
          ! branch.  Both are branch-free apart from selects, and vectorise.
          if (split_rates) then
            do f=self%fs, self%fe
              xxup  = work_up(f,1)
              ssup  = work_up(f,2)
              Inuup = work_up(f,3)
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              dt = dsup05 * ( xxup + xx )     ! upwind optical pathlength
              dt = max(dt, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dt<0)
              call linear_weights(dt, E, W, Psi, Esc)
              Inu = Inuup*E + W*ssup + Psi*ss
              ! Hold the non-local part and the escape factor 1-Psi on their own
              ! as well: the preconditioned rates need them, and both Inu - Psi*ss
              ! and 1 - Psi would cancel in an optically thick interval.  Psi=0
              ! switches the preconditioning off at a point, so Irem then carries
              ! the whole (clamped) intensity.
              off = dt <= 0 .or. xx <= 0 .or. abs(Inu) > real(maser_Inu_max)
              Irem = Inuup*E + W*ssup
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              Psi = merge(0., Psi, off)
              self%Irem(f,x,y,z) = merge(Inu, Irem, off)
              self%escape(f,x,y,z) = merge(1., Esc, off)
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning there
                self%Psi(f,x,y,z) = merge(Psi / xx, 0., dt > 0 .and. xx > 0)
              endif
            end do ! f-loop
          else
            do f=self%fs, self%fe
              xxup  = work_up(f,1)
              ssup  = work_up(f,2)
              Inuup = work_up(f,3)
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              dt = dsup05 * ( xxup + xx )     ! upwind optical pathlength
              dt = max(dt, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dt<0)
              call linear_weights(dt, E, W, Psi, Esc)
              Inu = Inuup*E + W*ssup + Psi*ss
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning there
                self%Psi(f,x,y,z) = merge(Psi / xx, 0., dt > 0 .and. xx > 0)
              endif
            end do ! f-loop
          endif
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE short_char

  !===============================================================================
  module subroutine long_char_linear (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, Inuup, xxup, ssup, dsup05, dsup1, dsup2, xx, ss, Inu, Psi, Irem, Esc
    real(4) :: dt, E, W
    integer :: zs, ze, zup
    logical :: amr_boundary, get_rates, split_rates, off
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_linear', itimer=itimer)
    get_rates = self%compute_psi
    split_rates = self%split_rates

    ds = self%ds(3) * scaling%l
    dsup1 = 0.5*ds/abs(vec(3))

    amr_boundary = amr/='' .and. self%dir < 0
    if (amr_boundary) then
      if (associated(self%upstream)) ds = self%upstream%ds(3) * scaling%l
      dsup2 = 0.5*ds/abs(vec(3))
    endif

    zs = merge(1, self%nz, self%dir > 0)
    ze = merge(self%nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      zup = z - self%dir

      if (amr_boundary) then
        dsup05 = dsup2
        amr_boundary = .false.
      else
        dsup05 = dsup1
      endif

      do y=1, self%ny
        do x=1, self%nx
          ! One f-loop per rate formulation, so that neither carries the other's
          ! branch.  Both are branch-free apart from selects, and vectorise.
          if (split_rates) then
            do f=self%fs, self%fe
              xxup  = self%work(f,x,y,zup,1)
              ssup  = self%work(f,x,y,zup,2)
              Inuup = self%work(f,x,y,zup,3)
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              dt = dsup05 * ( xxup + xx )     ! upwind optical pathlength
              dt = max(dt, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dt<0)
              call linear_weights(dt, E, W, Psi, Esc)
              Inu = Inuup*E + W*ssup + Psi*ss
              ! Hold the non-local part and the escape factor 1-Psi on their own
              ! as well: the preconditioned rates need them, and both Inu - Psi*ss
              ! and 1 - Psi would cancel in an optically thick interval.  Psi=0
              ! switches the preconditioning off at a point, so Irem then carries
              ! the whole (clamped) intensity.
              off = dt <= 0 .or. xx <= 0 .or. abs(Inu) > real(maser_Inu_max)
              Irem = Inuup*E + W*ssup
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              Psi = merge(0., Psi, off)
              self%Irem(f,x,y,z) = merge(Inu, Irem, off)
              self%escape(f,x,y,z) = merge(1., Esc, off)
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning there
                self%Psi(f,x,y,z) = merge(Psi / xx, 0., dt > 0 .and. xx > 0)
              endif
            end do ! f-loop
          else
            do f=self%fs, self%fe
              xxup  = self%work(f,x,y,zup,1)
              ssup  = self%work(f,x,y,zup,2)
              Inuup = self%work(f,x,y,zup,3)
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              dt = dsup05 * ( xxup + xx )     ! upwind optical pathlength
              dt = max(dt, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dt<0)
              call linear_weights(dt, E, W, Psi, Esc)
              Inu = Inuup*E + W*ssup + Psi*ss
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning there
                self%Psi(f,x,y,z) = merge(Psi / xx, 0., dt > 0 .and. xx > 0)
              endif
            end do ! f-loop
          endif
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear

  !===============================================================================
  !> Weights of the linear step over an interval of optical depth t,
  !> Inu = E*Inu_up + W*S_up + Psi*S, in fp32 and branch-free:
  !>
  !>   E = exp(-t),  u0 = 1 - E,  Psi = 1 - u0/t,  W = u0 - Psi.
  !>
  !> The reduction t = n*ln2 + r (ln2 split in two, as in cephes expf, |r| <= ln2/2)
  !> gives E = s*(1-m) and u0 = (1-s) + s*m, with s = 2^-n spliced into the exponent
  !> field and m = 1-exp(-r) = r*q(r), so neither is formed as 1 - exp(-t).  Psi and W
  !> still cancel as t -> 0, so for |t| <= 1 Psi = t*P(t) and W = u0 - Psi; above
  !> that Psi = 1 - u0/t and W = u0/t - E, which keeps W ~ 1/t at large t from being
  !> formed as 1 - Psi.  Esc = 1 - Psi is returned for the same reason: it is u0/t
  !> above t = 1, and the preconditioned rates need it where Psi rounds to 1 in fp32.
  !> q and P are relative-error minimax fits of degree 5 on
  !> |r| <= ln2/2 and degree 7 on |t| <= 1; every output is within 4 fp32 ulp of the
  !> exact value over -2 <= t <= 90.  Only the exponential sees the clamp at t = 87,
  !> past which exp(-t) is below every quantity it multiplies.
  !>
  !> Both forms are evaluated at every point and one is kept by a bitwise select:
  !> merge() on the reals would let gfortran move each form into its own arm, and
  !> under its default -ftrapping-math it will not speculate floating-point arms,
  !> which keeps the calling loop scalar.  The work is split over exp_neg, lw_series
  !> and lw_select so that each piece, and then linear_weights, stays under
  !> gfortran's inlining limits even in short_char and short_char_hermite.
  !===============================================================================
  pure subroutine linear_weights (t, E, W, Psi, Esc)
    real(4), intent(in) :: t
    real(4), intent(out) :: E, W, Psi, Esc
    real(4) :: u0, ut, ps
    logical :: thin
    !-----------------------------------------------------------------------------
    call exp_neg(t, E, u0)
    ut = u0 / t
    ps = lw_series(t)
    thin = abs(t) <= 1.0
    Psi = lw_select(ps, 1.0 - ut, thin)
    Esc = lw_select(1.0 - ps, ut, thin)
    W = lw_select(u0 - ps, ut - E, thin)
  end subroutine linear_weights

  include '../HOP/exp_neg.inc'

  !===============================================================================
  !> Psi = t*P(t), the series form of linear_weights for |t| <= 1.
  pure real(4) function lw_series (t)
    real(4), intent(in) :: t
    ! LW_P is P, lowest order first
    real(4), parameter :: LW_P(0:7) = [ &
      +5.000000000e-01, -1.666666716e-01, +4.166666791e-02, -8.333330974e-03, &
      +1.388740959e-03, -1.983609400e-04, +2.522676368e-05, -2.834696488e-06 ]
    !-----------------------------------------------------------------------------
    lw_series = t*(LW_P(0) + t*(LW_P(1) + t*(LW_P(2) + t*(LW_P(3) + t*(LW_P(4) + t*(LW_P(5) + t*(LW_P(6) + t*LW_P(7))))))))
  end function lw_series

  !===============================================================================
  !> a where first, else b, by bit mask: both are already evaluated, so there is no
  !> floating-point arm for gfortran to guard (see linear_weights).
  pure real(4) function lw_select (a, b, first)
    real(4), intent(in) :: a, b
    logical, intent(in) :: first
    integer :: sel
    !-----------------------------------------------------------------------------
    sel = merge(-1, 0, first)
    lw_select = transfer(ior(iand(transfer(a, 1), sel), iand(transfer(b, 1), not(sel))), 1.0)
  end function lw_select

  !===============================================================================
  !> Parabolic transfer moments over an interval of optical depth t, in fp32 and
  !> branch-free:
  !>
  !>   E = exp(-t),  U0 = 1 - E,  U1 = t - U0,  q = U2/t with U2 = t^2 - 2*U1.
  !>
  !> E and U0 come from exp_neg's range reduction, so neither is formed from the
  !> other.  U1 and q cancel as t -> 0.  U1 = t*Psi with linear_weights' Psi, so it
  !> is t*lw_series(t) for |t| <= 1 and t - U0 above.  q = t^2*R(t) for |t| <= 2 and
  !> t - 2*U1/t above; its series reaches further than Psi's because t - 2*U1/t
  !> still loses ~10 ulp near t = -1.5.  R is a relative-error minimax fit of degree
  !> 8 on |t| <= 2; every output is within 4 fp32 ulp of the exact value over
  !> -2 <= t <= 87, and E is clamped there as in exp_neg.  The forms are combined by
  !> lw_select for the reason given at linear_weights.
  !===============================================================================
  pure module subroutine parab_moments (t, E, U0, U1, q)
    real(4), intent(in)  :: t
    real(4), intent(out) :: E, U0, U1, q
    ! PM_R is R = q/t^2, lowest order first
    real(4), parameter :: PM_R(0:8) = [ &
      +3.333333433e-01, -8.333329856e-02, +1.666668802e-02, -2.777901944e-03, &
      +3.968053206e-04, -4.948913192e-05, +5.514686109e-06, -5.894584660e-07, &
      +5.153443894e-08 ]
    real(4) :: qs
    !-----------------------------------------------------------------------------
    call exp_neg(t, E, U0)
    U1 = lw_select(t*lw_series(t), t - U0, abs(t) <= 1.0)
    qs = t*t*(PM_R(0) + t*(PM_R(1) + t*(PM_R(2) + t*(PM_R(3) + t*(PM_R(4) + t*(PM_R(5) &
       + t*(PM_R(6) + t*(PM_R(7) + t*PM_R(8)))))))))
    q = lw_select(qs, t - 2.0*(U1/t), abs(t) <= 2.0)
  end subroutine parab_moments

!===============================================================================
!> Quadratic formal solver, equivalent to the 'bezier' Sinterpol case of Multi3D's
!> shortchar_3D.f90.  The parabola over the integrated interval gives weights
!> (A, B, C) on (S_upwind, S_target, dS/dtau) from the moments U0 = 1-exp(-dt),
!> U1 = dt-U0, U2 = dt^2-2*U1, q = U2/dt.  The MC (monotonized-central) limiter on
!> dS/dtau keeps the parabola from overshooting and driving Inu or the Lambda
!> diagonal negative, without forming an explicit control point.
!>
!> The z==ze plane has no downwind neighbour, so its stencil is shifted one cell
!> upwind (interval z0 -> z, with dS/dtau at the upwind end), which swaps the roles
!> of the weights; it is peeled out of the z-loop.
!===============================================================================
  module subroutine long_char_quadratic (self, vec)
    class(rt_patch_t), target :: self
    real(8) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2
    real(4) :: Inu, ppsi
    real(4), parameter :: dtmin=1e-20
    ! hoist self%work / self%Psi and the loop bounds into locals: through a
    ! polymorphic dummy gfortran reloads the array descriptor in the f-loop and
    ! cannot prove unit stride, which blocks vectorisation
    real(4), dimension(:,:,:,:,:), contiguous, pointer :: w
    real(4), dimension(:,:,:,:), contiguous, pointer :: psi
    integer :: zs, ze, z0, zup, zdn, dir, nx, ny, fs, fe
    logical :: get_rates
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_quadratic', itimer=itimer)
    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./12.*dsup*dsup     ! Euler-Maclaurin correction to the trapezoid
    get_rates = self%compute_psi
    dir = self%dir
    nx = self%nx
    ny = self%ny
    fs = self%fs
    fe = self%fe
    zs = merge(1, self%nz, dir > 0)
    ze = merge(self%nz, 1, dir > 0)
    w => self%work
    psi => null()
    if (get_rates) psi => self%Psi

    ! interior points: dS/dtau sits at the target (downwind) end of the interval
    do z=zs, ze-dir, dir
      zup = z - dir
      zdn = z + dir
      do y=1, ny
        do x=1, nx
          do f=fs, fe
            call quad_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                            w(f,x,y,zup,3), .true., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    end do ! z-loop

    ! last point: the stencil is shifted one cell upwind, so the interval
    ! integrated is z0 -> z and dS/dtau sits at its upwind end, z0
    z = ze
    z0 = z - dir
    zup = z0 - dir
    do y=1, ny
      do x=1, nx
        do f=fs, fe
          call quad_step (w(f,x,y,zup,1), w(f,x,y,z0,1), w(f,x,y,z,1), &
                          w(f,x,y,zup,2), w(f,x,y,z0,2), w(f,x,y,z,2), &
                          w(f,x,y,z0,3), .false., Inu, ppsi)
          w(f,x,y,z,3) = Inu
          if (get_rates) psi(f,x,y,z) = ppsi
        end do ! f-loop
      end do
    end do
    call trace%end(itimer)
  contains
    !---------------------------------------------------------------------------
    !> One cell of the sweep, on the stencil (up, c, dn) centred at c.  With
    !> at_target the interval integrated is up -> c; otherwise it is c -> dn, with
    !> dS/dtau taken at c (then the upwind end).  Either way Inu comes out at the
    !> downwind end of the interval and plam is the Lambda diagonal there.
    !---------------------------------------------------------------------------
    pure subroutine quad_step (xxup, xxc, xxdn, ssup, ssc, ssdn, Inuup, at_target, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xxdn, ssup, ssc, ssdn, Inuup
      logical, intent(in) :: at_target
      real(4), intent(out) :: Inu, plam
      real(4) :: t1, t2, t3, dtup, dtdn, dsdtau
      real(4) :: dt, rdt, ex, U0, U1, q, A, B, C
      !-------------------------------------------------------------------------
      ! one sided dchi/ds either side of c, and the MC-limited central slope
      t1 = (xxc - xxup) * rdsup
      t3 = (xxdn - xxc) * rdsup
      t2 = 0.5*(t1 + t3)
      t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))
      ! maser interval (chi<0 -> dtau<0): keep the sign, hold |dtau| away from
      ! zero (it divides the slopes below) and cap the amplification
      dtup = dsup1*(xxup + xxc) + dsup2*(t1-t2)
      dtdn = dsup1*(xxdn + xxc) + dsup2*(t2-t3)
      dtup = max(sign(max(abs(dtup), dtmin), dtup), -real(maser_dtau_max))
      dtdn = max(sign(max(abs(dtdn), dtmin), dtdn), -real(maser_dtau_max))

      dt = merge(dtup, dtdn, at_target)
      rdt = 1.0 / dt
      call parab_moments(dt, ex, U0, U1, q)

      ! one sided dS/dtau either side of c, and the MC-limited central slope
      if (at_target) then
        t1 = (ssc - ssup) * rdt
        t3 = (ssdn - ssc) / dtdn
      else
        t1 = (ssc - ssup) / dtup
        t3 = (ssdn - ssc) * rdt
      endif
      t2 = 0.5*(t1 + t3)
      dsdtau = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

      ! for dt<0 (bounded maser) or chi<0 the rates are left unpreconditioned
      ! (plam=0), and Inu is capped, as in the other solvers
      plam = 0.0
      if (at_target) then
        B = (2.0*U1 - q) * rdt        ! on ssc  (the target)
        A = U0 - B                    ! on ssup (the upwind point)
        C = q - U1                    ! on dsdtau
        Inu = Inuup*ex + A*ssup + B*ssc + C*dsdtau
        if (get_rates .and. dt > 0 .and. xxc > 0) &
          plam = B / xxc            ! diagonal of the Lambda operator
      else
        A = q * rdt                   ! on ssdn (the target)
        B = U0 - A                    ! on ssc  (the upwind point)
        C = U1 - q                    ! on dsdtau
        Inu = Inuup*ex + A*ssdn + B*ssc + C*dsdtau
        if (get_rates .and. dt > 0 .and. xxdn > 0) &
          plam = A / xxdn           ! diagonal of the Lambda operator
      endif
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
    end subroutine quad_step
  END SUBROUTINE long_char_quadratic

!===============================================================================
!> Chain-rule cubic-Hermite formal solver ('hermite'), long-characteristics
!> version.  The source over each integrated interval [up, c] is a cubic Hermite
!> in optical depth with slopes at BOTH ends, obtained by product rules on the
!> RATIOS chi/Pe and Snu/Bnu rather than finite differences of chi and Snu:
!>
!>   dchi/ds   = (chi/Pe) * (dPe/ds)  + Pe  * d(chi/Pe)/ds
!>   dSnu/dtau = (Snu/Bnu) * (dBnu/dT) * (dT/ds) / chi + Bnu * d(Snu/Bnu)/ds / chi
!>
!> dPe/ds and dT/ds are the PCHIP-interpolant z-gradients of the snapshot
!> (comp%dpe_dz, comp%dlnt_dz; the ray runs along z for rotated/column runs);
!> dBnu/dT is analytic in the precomputed beta = exp(-h nu/kT).  The ratio slopes
!> d(chi/Pe)/ds and d(Snu/Bnu)/ds come from a PARABOLA through the upwind, central
!> and downwind samples (centred estimate at the target end, one-sided at the far
!> end), as in the quadratic solver.  dtau is the trapezoid with the
!> Euler-Maclaurin correction from the chi slopes at the two ends.
!>
!> The source cubic's two Bezier control points are clamped into
!> [min(Sup,Sc), max(Sup,Sc)] (monotone), the cubic analogue of the quadratic
!> solver's MC limiter; the weights are integrals of the Bernstein basis against
!> exp(-t), from series below dt = 3 and the moments Mn = int_0^dt t^n e^-t dt
!> above (bezier3_weights).
!>
!> The first plane integrates from the boundary ghost plane (chi/Snu/Inu from
!> set_upstream_boundary via halo exchange).  With a real upstream neighbour
!> (comp_up) or a deep BOTTOM boundary (ghost Pe/T extrapolated by ghost_thermo)
!> the full step (herm_plane) is used.  Only at the zeroed TOP boundary (chi==0, where the
!> analytic dBnu/dtau=.../chi is singular) does the upwind-end source slope
!> degrade to the interval secant (herm_first).
!===============================================================================
  module subroutine long_char_hermite (self, comp, vec, comp_up)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(composition_type), pointer :: comp_up
    real(8), intent(in) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2, gfac
    real(4) :: Inu, ppsi, aTc
    real(4) :: pec, pedn, rpec, rpedn, gpc
    real(4), parameter :: dtmin=1e-20
    real(4), dimension(:,:,:,:,:), contiguous, pointer :: w
    real(4), dimension(:,:,:,:), contiguous, pointer :: psi, bk
    real(4), dimension(:), allocatable, target :: hknu, bghost
    integer :: zs, ze, zup, zdn, zu, z1, z2, z3, dir, nx, ny, fs, fe
    logical :: get_rates, use_upstream, extrap_bottom
    real(4) :: xx3g(3)
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_hermite', itimer=itimer)
    if (.not. allocated(comp%dpe_dz)) &
      error stop 'long_char_hermite: composition gradients missing -- set_gradients was not called'
    if (self%nz < 2) &
      error stop 'long_char_hermite: needs nz >= 2'
    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./12.*dsup*dsup     ! Euler-Maclaurin correction to the trapezoid
    ! converts the stored d/dz (per code length unit) to d/ds along propagation
    gfac = real(self%dir * abs(vec(3)) / scaling%l, 4)
    get_rates = self%compute_psi
    dir = self%dir
    nx = self%nx
    ny = self%ny
    fs = self%fs
    fe = self%fe
    zs = merge(1, self%nz, dir > 0)
    ze = merge(self%nz, 1, dir > 0)
    w => self%work
    bk => self%beta
    psi => null()
    if (get_rates) psi => self%Psi
    allocate(hknu(fs:fe))
    hknu = real(h_k * self%nu(fs:fe), 4)

    ! first plane: integrate from the boundary ghost plane, whose chi/Snu are
    ! already set (a real neighbour, or set_upstream_boundary's bottom
    ! extrapolation / top zeroing).  With a neighbour or bottom boundary the ghost
    ! Pe/T are available (or extrapolated by ghost_thermo) and the full step runs
    ! in herm_plane; the zeroed top ghost (chi==0) keeps the secant fallback
    ! herm_first.
    use_upstream = associated(comp_up)
    ! ghost extrapolated (not zeroed) exactly when set_upstream_boundary
    ! extrapolates chi/Snu: the zl branch with zu not also set
    extrap_bottom = (.not. use_upstream) .and. &
                    self%boundaries%is_set(bits%zl) .and. &
                    (.not. self%boundaries%is_set(bits%zu))
    if (use_upstream) zu = merge(self%upstream%nz, 1, dir > 0)
    if (extrap_bottom) then
      z1 = zs                                   ! three innermost planes, nearest
      z2 = zs + dir                             ! the ghost first (distances 1,2,3
      z3 = zs + 2*dir                           ! * ds from the ghost plane)
      xx3g = [1.0, 2.0, 3.0]
    endif
    allocate(bghost(fs:fe))
    if (use_upstream .or. extrap_bottom) then
      call herm_plane (zs, zs+dir, .true.)
    else
      ! top boundary: the zeroed ghost has no thermodynamic state (chi==0)
      z = zs
      zup = z - dir
      zdn = z + dir
      do y=1, ny
        do x=1, nx
          pec   = comp%pe4(x,y,z)
          pedn  = comp%pe4(x,y,zdn)
          rpec  = 1.0 / pec
          rpedn = 1.0 / pedn
          gpc   = gfac * comp%dpe_dz(x,y,z)
          aTc   = gfac * comp%dlnt_dz(x,y,z) / real(comp%temp(x,y,z), 4)
          do f=fs, fe
            call herm_first (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                             bk(f,x,y,z), bk(f,x,y,zdn), hknu(f), &
                             w(f,x,y,zup,3), pec, rpec, pedn, rpedn, gpc, &
                             aTc, Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    endif

    ! interior planes with a downwind neighbour: the parabola is centred at the
    ! target z (third sample = zdn), so its centred slope lands at the downwind
    ! end of the integrated interval [zup, z] and its one-sided slope at zup
    do z=zs+dir, ze-dir, dir
      call herm_plane (z, z+dir, .true.)
    end do

    ! last plane z==ze: no downwind sample exists, so the parabola is centred at
    ! the upwind end zup (third sample two planes upwind) and its centred slope
    ! lands at zup, the one-sided slope at the target z -- the same stencil shift
    ! the quadratic solver uses at its final point
    call herm_plane (ze, ze-2*dir, .false.)
    call trace%end(itimer)
  contains
    !---------------------------------------------------------------------------
    !> One plane z, integrated over [zup, z] with zup = z - dir.  zt is the third
    !> parabola sample: downwind with center_target, two planes upwind without.
    !> Every plane that takes the full step goes through here, so herm_dtau,
    !> bezier3_weights and herm_source each have a single call site, and gfortran
    !> inlines them into the f-loop whatever their size; called from several
    !> places they exceed its inlining limits and the loop stays scalar.
    !>
    !> The upwind point is the boundary ghost on the first plane.  Its Pe/T come
    !> from the neighbour or from ghost_thermo, and its beta is rebuilt from T as
    !> lte_background_opacities does, since the ghost carries chi/Snu/Inu but not
    !> beta.  The third sample is the ghost only on the last plane of a patch with
    !> nz == 2: from the neighbour, or at the zeroed top ghost the target's
    !> beta/Pe, which need only be finite since chi/Pe and Snu/Bnu vanish there.
    !---------------------------------------------------------------------------
    subroutine herm_plane (z, zt, center_target)
      integer, intent(in) :: z, zt
      logical, intent(in) :: center_target
      real(4), dimension(:,:,:,:,:), contiguous, pointer :: wp
      real(4), dimension(:,:,:,:), contiguous, pointer :: bp, pp
      real(4), dimension(:), contiguous, pointer :: hp, bu, bt
      real(4) :: peu, pec, rpeu, rpec, rpet, gpu, gpc, aTu, aTc, rTup
      real(4) :: Inu, ppsi, dt, ex, w0, w1, w2, w3
      logical :: ghost_up, ghost_t, rates
      integer :: x, y, f, zup
      !-------------------------------------------------------------------------
      ! local copies of the host pointers and of get_rates: through the host
      ! frame, which the f-loop's stores may alias, gfortran reloads them every
      ! iteration, and then cannot vectorise the loads or hoist the psi branch
      wp => w
      bp => bk
      pp => psi
      hp => hknu
      rates = get_rates
      zup = z - dir
      ghost_up = z == zs
      ghost_t = zt < 1 .or. zt > self%nz
      do y=1, ny
        do x=1, nx
          pec  = comp%pe4(x,y,z)
          rpec = 1.0 / pec
          gpc  = gfac * comp%dpe_dz(x,y,z)
          aTc  = gfac * comp%dlnt_dz(x,y,z) / real(comp%temp(x,y,z), 4)
          if (ghost_up) then
            if (use_upstream) then
              peu  = comp_up%pe4(x,y,zu)
              rpeu = 1.0 / peu
              gpu  = gfac * comp_up%dpe_dz(x,y,zu)
              rTup = 1.0 / real(comp_up%temp(x,y,zu), 4)
              aTu  = gfac * comp_up%dlnt_dz(x,y,zu) * rTup
            else
              ! ghost Pe/T extrapolated from the three innermost planes,
              ! mirroring the chi/Snu extrapolation in set_upstream_boundary
              call ghost_thermo (x, y, peu, rpeu, gpu, rTup, aTu)
            endif
            bghost(fs:fe) = max(exp(-hp*rTup), tiny(1.0))
            bu(fs:) => bghost
          else
            peu  = comp%pe4(x,y,zup)
            rpeu = 1.0 / peu
            gpu  = gfac * comp%dpe_dz(x,y,zup)
            aTu  = gfac * comp%dlnt_dz(x,y,zup) / real(comp%temp(x,y,zup), 4)
            bu(fs:) => bp(fs:fe,x,y,zup)
          endif
          if (.not. ghost_t) then
            rpet = 1.0 / comp%pe4(x,y,zt)
            bt(fs:) => bp(fs:fe,x,y,zt)
          else if (use_upstream) then
            rpet = 1.0 / comp_up%pe4(x,y,zu)
            bghost(fs:fe) = max(exp(-hp/real(comp_up%temp(x,y,zu), 4)), tiny(1.0))
            bt(fs:) => bghost
          else
            rpet = rpec
            bt(fs:) => bp(fs:fe,x,y,z)
          endif
          do f=fs, fe
            dt = herm_dtau (wp(f,x,y,zup,1), wp(f,x,y,z,1), wp(f,x,y,zt,1), &
                            peu, rpeu, gpu, pec, rpec, gpc, rpet, center_target)
            call bezier3_weights (dt, ex, w0, w1, w2, w3)
            call herm_source (wp(f,x,y,zup,1), wp(f,x,y,z,1), wp(f,x,y,zup,2), wp(f,x,y,z,2), wp(f,x,y,zt,2), &
                              bu(f), bp(f,x,y,z), bt(f), hp(f), wp(f,x,y,zup,3), &
                              aTu, aTc, center_target, dt, ex, w0, w1, w2, w3, Inu, ppsi)
            wp(f,x,y,z,3) = Inu
            if (rates) pp(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    end subroutine herm_plane
    !---------------------------------------------------------------------------
    !> Ghost-plane thermodynamics at a deep (bottom) domain boundary, where there
    !> is no upstream neighbour.  set_upstream_boundary fills the ghost chi/Snu by
    !> monotone extrapolation of the three innermost planes; this does the same
    !> for Pe and T (and their d/dz gradients), returning the exact quantities the
    !> upstream branch reads from a real neighbour, so herm_plane uses it unchanged.
    !---------------------------------------------------------------------------
    subroutine ghost_thermo (x, y, pe_g, rpe_g, gp_g, rT_g, aT_g)
      integer, intent(in) :: x, y
      real(4), intent(out) :: pe_g, rpe_g, gp_g, rT_g, aT_g
      real(4) :: yv(3), gv(1), T_g, dlnt_g, dpe_g
      real(4), parameter :: xg(1) = [0.0]        ! ghost sits one plane below z1
      !-------------------------------------------------------------------------
      yv = [comp%pe4(x,y,z1), comp%pe4(x,y,z2), comp%pe4(x,y,z3)]
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      pe_g = max(gv(1), tiny(1.0))
      yv = [comp%dpe_dz(x,y,z1), comp%dpe_dz(x,y,z2), comp%dpe_dz(x,y,z3)]
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      dpe_g = gv(1)
      yv = real([comp%temp(x,y,z1), comp%temp(x,y,z2), comp%temp(x,y,z3)], 4)
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      T_g = max(gv(1), tiny(1.0))
      yv = [comp%dlnt_dz(x,y,z1), comp%dlnt_dz(x,y,z2), comp%dlnt_dz(x,y,z3)]
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      dlnt_g = gv(1)
      rpe_g = 1.0 / pe_g
      gp_g  = gfac * dpe_g
      rT_g  = 1.0 / T_g
      aT_g  = gfac * dlnt_g * rT_g
    end subroutine ghost_thermo
    !---------------------------------------------------------------------------
    !> Cubic-Bezier weights over an interval of optical depth t, in fp32 and
    !> branch-free: integrals of the Bernstein basis against exp(-tau), tau
    !> measured backward from the target.  w0 multiplies S_up, w1/w2 the upwind-/
    !> target-side control points, w3 S_target, and E = exp(-t) is returned with
    !> them.  With u = tau/t,
    !>
    !>   w_k = t * int_0^1 b_k(u) exp(-t*u) du,   b = u^3, 3u^2(1-u), 3u(1-u)^2, (1-u)^3,
    !>
    !> and reflecting u -> 1-u gives w0(t) = t*E*P3(-t) and w1(t) = t*E*P2(-t), with
    !> P_k = w_k/t.  For t <= 3 the weights are these series.  P2 and P3 are
    !> relative-error minimax fits of degree 10 and 9 on |t| <= 3 (b3_p2, b3_p3);
    !> the reflection evaluates w0 and w1, whose basis functions peak at the far
    !> end, where their series has terms of one sign.  Above t = 3 the weights come
    !> from the moments Mn = int_0^t tau^n exp(-tau) dtau, by the upward recurrence
    !> Mn = n*M(n-1) - t^n*E from M0 = 1 - E, with w3 = M0 - w0 - w1 - w2 so that
    !> the weights sum to M0 there.  The moment combinations cancel below t ~ 2.5
    !> and for every t < 0.  E and M0 come from exp_neg.  Every weight is within
    !> 8 fp32 ulp of the exact value over -3 <= t <= 87.  The moment arm is formed
    !> at max(t, 3) so that it stays finite where it is not selected, and the arms
    !> are combined by lw_select (see linear_weights).
    !---------------------------------------------------------------------------
    pure subroutine bezier3_weights (t, E, w0, w1, w2, w3)
      real(4), intent(in) :: t
      real(4), intent(out) :: E, w0, w1, w2, w3
      real(4), parameter :: tser = 3.0
      real(4) :: M0, M1, M2, M3, tb, tb2, tb3, r3, tE, b0, b1, b2
      logical :: ser
      !-------------------------------------------------------------------------
      call exp_neg(t, E, M0)
      tb  = max(t, tser)
      tb2 = tb*tb
      tb3 = tb2*tb
      r3  = 1.0 / tb3
      M1 = M0 - tb*E
      M2 = 2.0*M1 - tb2*E
      M3 = 3.0*M2 - tb3*E
      b0 = M3*r3
      b1 = 3.0*(M2*tb - M3)*r3
      b2 = 3.0*(M1*tb2 - 2.0*M2*tb + M3)*r3
      ser = t <= tser
      tE  = t*E
      w0 = lw_select(tE*b3_p3(-t), b0, ser)
      w1 = lw_select(tE*b3_p2(-t), b1, ser)
      w2 = lw_select(t*b3_p2(t), b2, ser)
      w3 = lw_select(t*b3_p3(t), M0 - b0 - b1 - b2, ser)
    end subroutine bezier3_weights
    !---------------------------------------------------------------------------
    !> P2 = w2/t of bezier3_weights, for |t| <= 3.
    pure real(4) function b3_p2 (t)
      real(4), intent(in) :: t
      ! B3_C2 is P2, lowest order first
      real(4), parameter :: B3_C2(0:10) = [ &
        +2.500000000e-01, -1.000000238e-01, +2.499995567e-02, -4.761847202e-03, &
        +7.440865738e-04, -9.924791084e-05, +1.156294547e-05, -1.190877015e-06, &
        +1.137804801e-07, -1.100420821e-08, +7.658396139e-10 ]
      !-------------------------------------------------------------------------
      b3_p2 = B3_C2(0) + t*(B3_C2(1) + t*(B3_C2(2) + t*(B3_C2(3) + t*(B3_C2(4) + t*(B3_C2(5) &
            + t*(B3_C2(6) + t*(B3_C2(7) + t*(B3_C2(8) + t*(B3_C2(9) + t*B3_C2(10))))))))))
    end function b3_p2
    !---------------------------------------------------------------------------
    !> P3 = w3/t of bezier3_weights, for |t| <= 3.
    pure real(4) function b3_p3 (t)
      real(4), intent(in) :: t
      ! B3_C3 is P3, lowest order first
      real(4), parameter :: B3_C3(0:9) = [ &
        +2.500000000e-01, -4.999998957e-02, +8.333289996e-03, -1.190488110e-03, &
        +1.488493726e-04, -1.653099389e-05, +1.640762775e-06, -1.503805294e-07, &
        +1.414234507e-08, -1.025371787e-09 ]
      !-------------------------------------------------------------------------
      b3_p3 = B3_C3(0) + t*(B3_C3(1) + t*(B3_C3(2) + t*(B3_C3(3) + t*(B3_C3(4) + t*(B3_C3(5) &
            + t*(B3_C3(6) + t*(B3_C3(7) + t*(B3_C3(8) + t*B3_C3(9)))))))))
    end function b3_p3
    !---------------------------------------------------------------------------
    !> One integrated interval [up, c] is herm_dtau, then bezier3_weights, then
    !> herm_source, each called from the f-loop.  chi/Pe and Snu/Bnu are each
    !> fitted by a PARABOLA through the upwind, central and third samples
    !> (xx3/ss3/...); their derivatives set the product-rule chi slopes (-> dtau)
    !> and source slopes (-> Bezier control points).  center_target picks which
    !> end of [up,c] gets the centred 3-point estimate:
    !>   .true.  -> parabola centred at the target c (interior sweep; third sample
    !>              downwind), one-sided slope at up.
    !>   .false. -> parabola centred at up (final plane z==ze; third sample two
    !>              planes upwind), one-sided slope at c -- the quadratic end shift.
    !> The step is split in three so that each piece stays under gfortran's
    !> inlining limits and the f-loop can vectorise; as one routine it is not
    !> inlined.
    !---------------------------------------------------------------------------
    !> dtau over [up, c]: dchi/ds at both ends by the product rule
    !> chi = (chi/Pe)*Pe, clamped to the sign of the interval secant and twice its
    !> magnitude (MC envelope), then the Euler-Maclaurin trapezoid.
    pure real(4) function herm_dtau (xxup, xxc, xx3, peu, rpeu, gpu, pec, rpec, gpc, rpe3, &
                                     center_target) result (dt)
      real(4), intent(in) :: xxup, xxc, xx3, peu, rpeu, gpu, pec, rpec, gpc, rpe3
      logical, intent(in) :: center_target
      real(4) :: t1, sgn, mup, mc, rxpu, rxpc, rxp3, drxpu, drxpc, hh
      !-------------------------------------------------------------------------
      hh   = 0.5 * rdsup                      ! 1/(2*ds), the parabola half-step
      t1   = (xxc - xxup) * rdsup             ! chi secant over [up,c] (limiter ref)
      sgn  = sign(1.0, t1)
      rxpu = xxup * rpeu                      ! chi/Pe at the upwind point
      rxpc = xxc  * rpec                      ! chi/Pe at the target
      rxp3 = xx3  * rpe3                      ! chi/Pe at the third parabola sample
      if (center_target) then                ! centre = c (third sample downwind)
        drxpc = (rxp3 - rxpu) * hh                     ! centred at c
        drxpu = (4.0*rxpc - 3.0*rxpu - rxp3) * hh       ! left one-sided at up
      else                                   ! centre = up (third sample 2 upwind)
        drxpu = (rxpc - rxp3) * hh                      ! centred at up
        drxpc = (rxp3 - 4.0*rxpu + 3.0*rxpc) * hh        ! right one-sided at c
      endif
      mup = rxpu*gpu + peu*drxpu
      mc  = rxpc*gpc + pec*drxpc
      mup = sgn * max(0.0, min(sgn*mup, 2.0*abs(t1)))
      mc  = sgn * max(0.0, min(sgn*mc,  2.0*abs(t1)))
      ! maser interval (chi<0 -> dtau<0): keep the sign, hold |dtau| away from
      ! zero and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(mup - mc)
      dt = max(sign(max(abs(dt), dtmin), dt), -real(maser_dtau_max))
    end function herm_dtau
    !---------------------------------------------------------------------------
    !> Inu and the Lambda diagonal at c from the Bezier weights over [up, c].
    !> dSnu/dtau at both ends by the product rule S = rho*Bnu (rho = S/Bnu):
    !>   dS/dtau = rho*(dBnu/dtau) + Bnu*(drho/dtau)
    !> the analytic dBnu/dtau = hk*B*r1mb*aT/chi always enters as rho*dBnu/dtau
    !> (= S*(dBnu/dtau)/B) or as (dBnu/dtau)/B in the Lambda term, so its B factor
    !> cancels rho's 1/B: gu/gc hold (dBnu/dtau)/Bnu = hk*r1mb*aT/chi.  The
    !> departure part is the parabola's d(rho)/ds divided by chi.
    !---------------------------------------------------------------------------
    pure subroutine herm_source (xxup, xxc, ssup, ssc, ss3, betau, betac, beta3, hk, &
                                 Inuup, aTu, aTc, center_target, dt, ex, w0, w1, w2, w3, &
                                 Inu, plam)
      real(4), intent(in) :: xxup, xxc, ssup, ssc, ss3, betau, betac, beta3, hk, Inuup
      real(4), intent(in) :: aTu, aTc, dt, ex, w0, w1, w2, w3
      logical, intent(in) :: center_target
      real(4), intent(out) :: Inu, plam
      real(4) :: r1mbu, r1mbc, r1mb3, Bu, Bc, B3, gu, gc, rxxu, rxxc, hh
      real(4) :: dt3, smin, smax, P1, P2, rhou, rhoc, rho3, drhu, drhc, d1, d2
      !-------------------------------------------------------------------------
      ! source departure ratios rho = Snu/Bnu at the three samples, and their
      ! 3-point parabola derivatives d(rho)/ds at both ends of [up,c]
      hh = 0.5 * rdsup
      r1mbu = 1.0 / (1.0 - betau)
      r1mbc = 1.0 / (1.0 - betac)
      r1mb3 = 1.0 / (1.0 - beta3)
      Bu = betau * r1mbu
      Bc = betac * r1mbc
      B3 = beta3 * r1mb3
      rhou = ssup / Bu
      rhoc = ssc  / Bc
      rho3 = ss3  / B3
      if (center_target) then
        drhc = (rho3 - rhou) * hh
        drhu = (4.0*rhoc - 3.0*rhou - rho3) * hh
      else
        drhu = (rhoc - rho3) * hh
        drhc = (rho3 - 4.0*rhou + 3.0*rhoc) * hh
      endif

      rxxu = 1.0 / xxup                      ! 1/chi at up and c
      rxxc = 1.0 / xxc
      gu = hk * r1mbu * aTu * rxxu           ! (dBnu/dtau)/Bnu at up
      gc = hk * r1mbc * aTc * rxxc           ! (dBnu/dtau)/Bnu at c
      dt3 = dt * real(one_third,4)
      ! Sup + (dtau/3)*(dS/dtau)|up and Sc - (dtau/3)*(dS/dtau)|c, with
      ! (drho/dtau) = (drho/ds)/chi
      P1 = ssup + dt3*(ssup*gu + Bu*drhu*rxxu)
      P2 = ssc  - dt3*(ssc*gc  + Bc*drhc*rxxc)
      smin = min(ssup, ssc)                  ! monotone clamp: overshoot-free
      smax = max(ssup, ssc)
      P1 = min(max(P1, smin), smax)
      P2 = min(max(P2, smin), smax)

      Inu = Inuup*ex + w0*ssup + w1*P1 + w2*P2 + w3*ssc
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      ! Lambda diagonal from dP1/dSc = Bu/(3*Bc), dP2/dSc = 2/3 - dt3*gc (cross-
      ! terms and clamp neglected), both bounded into [0,1] so it cannot exceed
      ! the attenuated total or go negative.  For dt<0/chi<0 plam=0.  It is
      ! formed everywhere and kept by lw_select, which leaves the calling f-loop
      ! free of branches (see linear_weights).
      d1 = min(max(Bu*real(one_third,4)/Bc, 0.0), 1.0)
      d2 = min(max(2.0*real(one_third,4) - dt3*gc, 0.0), 1.0)
      plam = lw_select((w1*d1 + w2*d2 + w3) * rxxc, 0.0, get_rates .and. dt > 0 .and. xxc > 0)
    end subroutine herm_source
    !---------------------------------------------------------------------------
    !> The boundary interval [ghost, zs]: state-space slopes from the adjacent
    !> interior interval [zs, zs+dir] (the ghost plane has no beta/lnPe), and
    !> the upwind-end source slope degrades to the interval secant.
    !---------------------------------------------------------------------------
    pure subroutine herm_first (xxup, xxc, xxdn, ssup, ssc, ssdn, betac, betadn, &
                                hk, Inuup, pec, rpec, pedn, rpedn, gpc, aTc, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xxdn, ssup, ssc, ssdn, betac, betadn, hk, Inuup
      real(4), intent(in) :: pec, rpec, pedn, rpedn, gpc, aTc
      real(4), intent(out) :: Inu, plam
      real(4) :: t1, sgn, mc, r1mbc, r1mbdn, Bc, Bdn, Gc, dSc, dSu, rxpc
      real(4) :: dt3, smin, smax, P1, P2, dss, rhoc, rhodn, taudn
      real(4) :: dt, ex, w0, w1, w2, w3
      !-------------------------------------------------------------------------
      ! target-end chi slope by the same chi = (chi/Pe)*Pe product rule as
      ! herm_dtau; the ghost plane has no Pe, so the one-sided d(chi/Pe)/ds
      ! comes from the adjacent DOWN interval
      t1 = (xxc - xxup) * rdsup
      sgn = sign(1.0, t1)
      rxpc = xxc*rpec                          ! chi/Pe at c (reused in the secant)
      mc = rxpc*gpc + pec*(xxdn*rpedn - rxpc)*rdsup
      mc = sgn * max(0.0, min(sgn*mc, 2.0*abs(t1)))
      ! maser interval: keep the sign of dtau, bound it away from zero (it
      ! divides the secant slope below) and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(t1 - mc)
      dt = max(sign(max(abs(dt), dtmin), dt), -real(maser_dtau_max))
      call bezier3_weights (dt, ex, w0, w1, w2, w3)

      r1mbc  = 1.0 / (1.0 - betac)
      r1mbdn = 1.0 / (1.0 - betadn)
      Bc  = betac  * r1mbc
      Bdn = betadn * r1mbdn
      ! target-end slope by the same product rule S = rho*Bnu as herm_source, the
      ! Bnu again cancelling rho's 1/Bnu (Gc = (dBnu/dtau)/Bnu); the departure
      ! ratio secant across the DOWN interval uses the leading-trapezoid dtau
      Gc = hk * r1mbc * aTc / xxc              ! (dBnu/dtau)/Bnu at c
      rhoc  = ssc  / Bc
      rhodn = ssdn / Bdn
      taudn = dsup1 * (xxc + xxdn)             ! leading trapezoid dtau, down side
      taudn = sign(max(abs(taudn), dtmin), taudn)
      dSc = ssc*Gc + Bc*(rhodn - rhoc)/taudn
      dss = ssc - ssup
      dSu = dss / dt                   ! interval secant slope wrt tau

      dt3 = dt * real(one_third,4)
      P1 = ssup + dt3*dSu
      P2 = ssc  - dt3*dSc
      smin = min(ssup, ssc)                    ! monotone clamp: overshoot-free
      smax = max(ssup, ssc)
      P1 = min(max(P1, smin), smax)
      P2 = min(max(P2, smin), smax)

      Inu = Inuup*ex + w0*ssup + w1*P1 + w2*P2 + w3*ssc
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      ! P1's secant slope contributes dP1/dSc = 1/3; P2's sensitivity is >= 1
      ! and is bounded to 1, as in herm_source; for dt<0 or chi<0 the rates are
      ! left unpreconditioned, as in the other solvers
      plam = lw_select((w1*real(one_third,4) + w2 + w3) / xxc, 0.0, &
                       get_rates .and. dt > 0 .and. xxc > 0)
    end subroutine herm_first
  END SUBROUTINE long_char_hermite

!===============================================================================
!> Hybrid formal solver ('hermite2'), long-characteristics.  Optical depth is
!> built exactly as long_char_hermite's chain-rule cubic-Hermite chi treatment:
!>
!>   dchi/ds = (chi/Pe) * (dPe/ds)  +  Pe * d(chi/Pe)/ds
!>
!> with the product-rule slopes from a parabola through the upwind/central/
!> downwind samples, clamped to the cubic-Bezier control-point bound
!> (|dchi/ds| <= 3*|secant|, the monotonicity limit matching the cubic-Hermite
!> chi that the Euler-Maclaurin trapezoid dtau integrates exactly).  The same
!> Euler-Maclaurin trapezoid dtau as herm_dtau.  The source term, however, is
!> NOT the cubic Bezier / Bnu-chain-rule reconstruction of long_char_hermite:
!> it is long_char_quadratic's plain parabolic Sinterpol -- a single MC-limited
!> dS/dtau at the target, built from the one-sided S-secants on the two
!> intervals adjacent to the target, combined with parab_moments' (A,B,C)
!> weights exactly as in quad_step's at_target branch.  The second interval's
!> dtau (needed only for that S-secant) is obtained the same chain-rule way,
!> using the backward-one-sided chi/Pe slope at the third sample from the very
!> same parabola used for the primary interval -- no extra grid point needed.
!>
!> At the outermost plane there is no sample beyond the target, but the two-
!> planes-upwind sample still is real domain data (nz>=3), so dS/dtau at the
!> target is recovered to second order by a backward one-sided difference on
!> that stencil -- keeping the EMERGENT surface source second order, as
!> quad_step's final-plane branch does.  Only where the third sample is a
!> boundary ghost (the zeroed-top first plane, or the nz==2 last plane) does
!> dS/dtau degrade to the plain interval secant, mirroring herm_first's
!> analogous fallback.  Because the source treatment never needs Bnu/beta/T,
!> this solver only requires the Pe gradients (comp%dpe_dz), not the ln(T) ones.
!===============================================================================
  module subroutine long_char_hermite2 (self, comp, vec, comp_up)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(composition_type), pointer :: comp_up
    real(8), intent(in) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2, gfac
    real(4) :: Inu, ppsi
    real(4) :: peu, pec, pedn, rpeu, rpec, rpedn, gpu, gpc, gpdn
    real(4), parameter :: dtmin=1e-20
    real(4), dimension(:,:,:,:,:), contiguous, pointer :: w
    real(4), dimension(:,:,:,:), contiguous, pointer :: psi
    integer :: zs, ze, zup, zdn, zff, zu, z1, z2, z3, dir, nx, ny, fs, fe
    logical :: get_rates, use_upstream, extrap_bottom
    real(4) :: xx3g(3)
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_hermite2', itimer=itimer)
    if (.not. allocated(comp%dpe_dz)) &
      error stop 'long_char_hermite2: composition gradients missing -- set_gradients was not called'
    if (self%nz < 2) &
      error stop 'long_char_hermite2: needs nz >= 2'
    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./12.*dsup*dsup     ! Euler-Maclaurin correction to the trapezoid
    ! converts the stored d/dz (per code length unit) to d/ds along propagation
    gfac = real(self%dir * abs(vec(3)) / scaling%l, 4)
    get_rates = self%compute_psi
    dir = self%dir
    nx = self%nx
    ny = self%ny
    fs = self%fs
    fe = self%fe
    zs = merge(1, self%nz, dir > 0)
    ze = merge(self%nz, 1, dir > 0)
    w => self%work
    psi => null()
    if (get_rates) psi => self%Psi

    ! first plane: as long_char_hermite, integrate from the boundary ghost
    ! plane.  With a neighbour or bottom boundary the ghost Pe is available (or
    ! extrapolated by ghost_pe) and herm2_step runs with a real downwind
    ! sample; the zeroed top ghost (chi==0) keeps the secant fallback herm2_first.
    use_upstream = associated(comp_up)
    z = zs
    zup = z - dir
    zdn = z + dir
    extrap_bottom = (.not. use_upstream) .and. &
                    self%boundaries%is_set(bits%zl) .and. &
                    (.not. self%boundaries%is_set(bits%zu))
    if (use_upstream) zu = merge(self%upstream%nz, 1, dir > 0)
    if (extrap_bottom) then
      z1 = zs
      z2 = zs + dir
      z3 = zs + 2*dir
      xx3g = [1.0, 2.0, 3.0]
    endif
    do y=1, ny
      do x=1, nx
        if (use_upstream) then
          peu  = comp_up%pe4(x,y,zu)
          rpeu = 1.0 / peu
          gpu  = gfac * comp_up%dpe_dz(x,y,zu)
        else if (extrap_bottom) then
          call ghost_pe (x, y, peu, rpeu, gpu)
        endif
        if (use_upstream .or. extrap_bottom) then
          pec   = comp%pe4(x,y,z)
          rpec  = 1.0 / pec
          gpc   = gfac * comp%dpe_dz(x,y,z)
          pedn  = comp%pe4(x,y,zdn)
          rpedn = 1.0 / pedn
          gpdn  = gfac * comp%dpe_dz(x,y,zdn)
          do f=fs, fe
            call herm2_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                             w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, &
                             pedn, rpedn, gpdn, .true., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        else
          ! top boundary: the zeroed ghost has no thermodynamic state (chi==0)
          pec   = comp%pe4(x,y,z)
          pedn  = comp%pe4(x,y,zdn)
          rpec  = 1.0 / pec
          rpedn = 1.0 / pedn
          gpc   = gfac * comp%dpe_dz(x,y,z)
          do f=fs, fe
            call herm2_first (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                              w(f,x,y,zup,2), w(f,x,y,z,2), &
                              w(f,x,y,zup,3), pec, rpec, rpedn, gpc, &
                              Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        endif
      end do
    end do

    ! interior planes with a downwind neighbour: the parabola is centred at
    ! the target z (third sample = zdn), matching long_char_hermite
    do z=zs+dir, ze-dir, dir
      zup = z - dir
      zdn = z + dir
      do y=1, ny
        do x=1, nx
          peu   = comp%pe4(x,y,zup)
          pec   = comp%pe4(x,y,z)
          pedn  = comp%pe4(x,y,zdn)
          rpeu  = 1.0 / peu
          rpec  = 1.0 / pec
          rpedn = 1.0 / pedn
          gpu   = gfac * comp%dpe_dz(x,y,zup)
          gpc   = gfac * comp%dpe_dz(x,y,z)
          gpdn  = gfac * comp%dpe_dz(x,y,zdn)
          do f=fs, fe
            call herm2_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                             w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, &
                             pedn, rpedn, gpdn, .true., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    end do ! z-loop

    ! last plane z==ze: no downwind sample exists, so the parabola is centred
    ! at the upwind end zup (third sample = zff, two planes upwind), as in
    ! long_char_hermite.  Snu has no data beyond the target, so herm2_step's
    ! secant fallback applies -- but rxp3=xx3*rpe3 still feeds the (unchanged,
    ! herm_dtau-identical) chi/Pe slope estimate at both ends even when
    ! center_target is false, so rpe3 must be the real Pe at the third sample,
    ! not a placeholder
    z   = ze
    zup = z - dir
    zff = z - 2*dir
    if (zff >= 1 .and. zff <= self%nz) then
      ! nz >= 3: the two-planes-upwind sample is in the domain
      do y=1, ny
        do x=1, nx
          peu  = comp%pe4(x,y,zup)
          pec  = comp%pe4(x,y,z)
          rpeu = 1.0 / peu
          rpec = 1.0 / pec
          gpu  = gfac * comp%dpe_dz(x,y,zup)
          gpc  = gfac * comp%dpe_dz(x,y,z)
          rpedn = 1.0 / comp%pe4(x,y,zff)   ! rpe3: real Pe at the third sample
          do f=fs, fe
            call herm2_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zff,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zff,2), &
                             w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, &
                             0.0, rpedn, 0.0, .false., Inu, ppsi, .true.)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    else
      ! nz == 2: the two-planes-upwind sample is the boundary ghost (zff ==
      ! zup of the first plane); its Pe comes from the neighbour when one
      ! exists, or (extrapolated-bottom / zeroed-top ghost) falls back to the
      ! target's own Pe -- as long_char_hermite's rpedn does in this branch
      do y=1, ny
        do x=1, nx
          peu  = comp%pe4(x,y,zup)
          pec  = comp%pe4(x,y,z)
          rpeu = 1.0 / peu
          rpec = 1.0 / pec
          gpu  = gfac * comp%dpe_dz(x,y,zup)
          gpc  = gfac * comp%dpe_dz(x,y,z)
          if (use_upstream) then
            rpedn = 1.0 / comp_up%pe4(x,y,zu)      ! ghost Pe from neighbour
          else
            rpedn = rpec                           ! finite fallback
          endif
          do f=fs, fe
            call herm2_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zff,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zff,2), &
                             w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, &
                             0.0, rpedn, 0.0, .false., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    endif
    call trace%end(itimer)
  contains
    !---------------------------------------------------------------------------
    !> Ghost-plane Pe (and its d/dz) at a deep (bottom) domain boundary, where
    !> there is no upstream neighbour: monotone extrapolation of the three
    !> innermost planes, mirroring long_char_hermite's ghost_thermo but
    !> without the T/beta part that herm2_step never needs.
    !---------------------------------------------------------------------------
    subroutine ghost_pe (x, y, pe_g, rpe_g, gp_g)
      integer, intent(in) :: x, y
      real(4), intent(out) :: pe_g, rpe_g, gp_g
      real(4) :: yv(3), gv(1), dpe_g
      real(4), parameter :: xg(1) = [0.0]        ! ghost sits one plane below z1
      !-------------------------------------------------------------------------
      yv = [comp%pe4(x,y,z1), comp%pe4(x,y,z2), comp%pe4(x,y,z3)]
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      pe_g = max(gv(1), tiny(1.0))
      yv = [comp%dpe_dz(x,y,z1), comp%dpe_dz(x,y,z2), comp%dpe_dz(x,y,z3)]
      call pchip_mono (xx3g, yv, xg, gv, sorted=.true.)
      dpe_g = gv(1)
      rpe_g = 1.0 / pe_g
      gp_g  = gfac * dpe_g
    end subroutine ghost_pe
    !---------------------------------------------------------------------------
    !> One integrated interval [up, c].  chi/Pe slopes as in herm_dtau: a
    !> parabola through (up, c, 3) gives the product-rule dchi/ds at up
    !> (one-sided) and at c (centred when center_target, one-sided otherwise),
    !> here clamped to the cubic-Bezier bound |dchi/ds| <= 3*|secant| (herm_dtau
    !> itself keeps the tighter 2* MC clamp).
    !>
    !> Snu is long_char_quadratic's parabolic Sinterpol: when center_target
    !> (interior sweep, 3 = real downwind sample), the SAME parabola also
    !> yields the backward-one-sided dchi/Pe/ds at 3, giving the [c,3]
    !> interval's chain-rule dtau; the one-sided S-secants on [up,c] and [c,3]
    !> are then MC-combined into a single dS/dtau at c, exactly as quad_step's
    !> at_target branch.  Otherwise (outermost plane: no data beyond the
    !> target) dS/dtau degrades to the plain [up,c] secant, and the
    !> third-sample Pe/gradient arguments (pe3, rpe3, gp3) are unused.
    !---------------------------------------------------------------------------
    pure subroutine herm2_step (xxup, xxc, xx3, ssup, ssc, ss3, Inuup, &
                                peu, rpeu, gpu, pec, rpec, gpc, pe3, rpe3, gp3, &
                                center_target, Inu, plam, src_onesided)
      real(4), intent(in) :: xxup, xxc, xx3, ssup, ssc, ss3, Inuup
      real(4), intent(in) :: peu, rpeu, gpu, pec, rpec, gpc, pe3, rpe3, gp3
      logical, intent(in) :: center_target
      real(4), intent(out) :: Inu, plam
      ! at the outermost plane the third sample (two planes upwind) is real
      ! domain data, so dS/dtau at the target can still be reconstructed to
      ! second order by a backward one-sided difference; set when that sample
      ! is NOT a boundary ghost (nz>=3 last plane).  Absent -> plain secant.
      logical, intent(in), optional :: src_onesided
      real(4) :: dchi, t1, t3, sgn, sgn3, mup, mc, m3, rxpu, rxpc, rxp3
      real(4) :: drxpu, drxpc, drxp3, hh, dss, dsdtau
      real(8) :: dt, dt3int, rdt, A, B, C
      real(4) :: ex, U0, U1, q
      !-------------------------------------------------------------------------
      hh   = 0.5 * rdsup                      ! 1/(2*ds), the parabola half-step
      dchi = xxc - xxup
      t1   = dchi * rdsup                     ! chi secant over [up,c] (limiter ref)
      sgn  = sign(1.0, t1)
      rxpu = xxup * rpeu                      ! chi/Pe at the upwind point
      rxpc = xxc  * rpec                      ! chi/Pe at the target
      rxp3 = xx3  * rpe3                      ! chi/Pe at the third parabola sample

      if (center_target) then                ! centre = c (third sample downwind)
        drxpc = (rxp3 - rxpu) * hh                     ! centred at c
        drxpu = (4.0*rxpc - 3.0*rxpu - rxp3) * hh       ! left one-sided at up
      else                                   ! centre = up (third sample 2 upwind)
        drxpu = (rxpc - rxp3) * hh                      ! centred at up
        drxpc = (rxp3 - 4.0*rxpu + 3.0*rxpc) * hh        ! right one-sided at c
      endif

      mup = rxpu*gpu + peu*drxpu
      mc  = rxpc*gpc + pec*drxpc
      ! chi = cubic Hermite in s -> Euler-Maclaurin dtau below is its exact
      ! integral; the matching monotonicity bound is the cubic-Bezier control-
      ! point-in-bracket limit |dchi/ds| <= 3*|secant| (not the 2* MC limit,
      ! which is the quadratic-Bezier bound).  3* keeps dtau > 0 -- worst case
      ! dt = h*(3*chi_up + chi_c)/4 > 0 -- while retaining more of the analytic
      ! Pe gradient in steep, under-resolved cells.
      mup = sgn * max(0.0, min(sgn*mup, 3.0*abs(t1)))
      mc  = sgn * max(0.0, min(sgn*mc,  3.0*abs(t1)))
      ! maser interval (chi<0 -> dtau<0): keep the sign, hold |dtau| away from
      ! zero and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(mup - mc)
      dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
      rdt = 1.0d0 / dt
      call parab_moments(real(dt,4), ex, U0, U1, q)

      dss = ssc - ssup
      if (center_target) then
        ! backward one-sided chi/Pe slope at the third sample (dn), from the
        ! same parabola, gives the [c,3] interval's chain-rule dtau the same
        ! way as [up,c]'s
        t3    = (xx3 - xxc) * rdsup
        sgn3  = sign(1.0, t3)
        drxp3 = (3.0*rxp3 - 4.0*rxpc + rxpu) * hh
        m3    = rxp3*gp3 + pe3*drxp3
        m3    = sgn3 * max(0.0, min(sgn3*m3, 3.0*abs(t3)))
        dt3int = dsup1*(xxc + xx3) + dsup2*(mc - m3)
        dt3int = max(sign(max(abs(dt3int), dble(dtmin)), dt3int), -maser_dtau_max)
        ! MC-limited central dS/dtau at the target, exactly as quad_step's
        ! at_target branch (one-sided secants on the two adjacent intervals)
        t1 = dss * real(rdt, 4)
        t3 = (ss3 - ssc) / real(dt3int, 4)
        dsdtau = 0.5*(t1 + t3)
        dsdtau = (sign(1.0,t1) + sign(1.0,t3)) * min(abs(t1), abs(t3), abs(0.5*dsdtau))
      else
        ! outermost plane: no data beyond the target.  When the third parabola
        ! sample is real domain data (src_onesided, nz>=3 last plane) recover a
        ! second-order dS/dtau at the target by a backward one-sided difference
        ! on the uniform-s stencil (3=two-up, up, c), /chi_c for per-tau, then
        ! MC-limited to the near-interval secant -- keeping the emergent surface
        ! source second order, as quad_step's final-plane branch does.  When the
        ! third sample is a boundary ghost (src_onesided absent) fall back to
        ! the plain secant, as herm_first does for its chain-rule slope.
        t1 = dss * real(rdt, 4)                       ! near-interval S secant
        dsdtau = t1
        if (present(src_onesided)) then
          if (src_onesided .and. xxc > 0.0) then
            dsdtau = (3.0*ssc - 4.0*ssup + ss3) * hh / xxc
            sgn = sign(1.0, t1)
            dsdtau = sgn * max(0.0, min(sgn*dsdtau, 2.0*abs(t1)))
          endif
        endif
      endif

      ! for dt<0 (bounded maser) or chi<0 the rates are left unpreconditioned
      ! (plam=0), and Inu is capped, as in the other solvers
      B = (2.0d0*U1 - q) * rdt      ! on ssc  (the target)
      A = U0 - B                    ! on ssup (the upwind point)
      C = q - U1                    ! on dsdtau
      Inu = Inuup*ex + real(A,4)*ssup + real(B,4)*ssc + real(C,4)*dsdtau
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      plam = 0.0
      if (get_rates .and. dt > 0 .and. xxc > 0) &
        plam = real(B,4) / xxc    ! diagonal of the Lambda operator
    end subroutine herm2_step
    !---------------------------------------------------------------------------
    !> The boundary interval [ghost, zs], zeroed-top case.  As herm_first, the
    !> target-end chi/Pe slope has no ghost Pe to centre on, so it falls back
    !> to the one-sided estimate from the adjacent DOWN interval; Snu is the
    !> plain interval secant, as herm2_step's outermost-plane branch.
    !---------------------------------------------------------------------------
    pure subroutine herm2_first (xxup, xxc, xxdn, ssup, ssc, Inuup, &
                                 pec, rpec, rpedn, gpc, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xxdn, ssup, ssc, Inuup
      real(4), intent(in) :: pec, rpec, rpedn, gpc
      real(4), intent(out) :: Inu, plam
      real(4) :: t1, sgn, mc, rxpc, dss, dsdtau
      real(8) :: dt, rdt, A, B, C
      real(4) :: ex, U0, U1, q
      !-------------------------------------------------------------------------
      ! target-end chi slope by the same chi = (chi/Pe)*Pe product rule as
      ! herm2_step; the ghost plane has no Pe, so the one-sided d(chi/Pe)/ds
      ! comes from the adjacent DOWN interval
      t1   = (xxc - xxup) * rdsup
      sgn  = sign(1.0, t1)
      rxpc = xxc*rpec
      mc = rxpc*gpc + pec*(xxdn*rpedn - rxpc)*rdsup
      mc = sgn * max(0.0, min(sgn*mc, 3.0*abs(t1)))
      ! maser interval: keep the sign of dtau, bound it away from zero (it
      ! divides the secant slope below) and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(t1 - mc)
      dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
      rdt = 1.0d0 / dt
      call parab_moments(real(dt,4), ex, U0, U1, q)

      dss = ssc - ssup
      dsdtau = dss * real(rdt, 4)      ! no data beyond the target: interval secant

      B = (2.0d0*U1 - q) * rdt
      A = U0 - B
      C = q - U1
      Inu = Inuup*ex + real(A,4)*ssup + real(B,4)*ssc + real(C,4)*dsdtau
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      plam = 0.0
      if (get_rates .and. dt > 0 .and. xxc > 0) &
        plam = real(B,4) / xxc
    end subroutine herm2_first
  END SUBROUTINE long_char_hermite2

!===============================================================================
!> Chain-rule hermite variant of the 3D short-characteristics solver
!> (rt_iteration_params solver='hermite' on the true-3D iteration path).  The 4x4 upwind
!> kernel interpolation and its 2x2 min limiter are exactly those of
!> short_char; the linear integration step is replaced by quad_step's weight
!> algebra with the product-rule slopes at the target cell, assuming the
!> ratios chi/Pe and Snu/Bnu vary linearly along the ray (same construction
!> as long_char_hermite):
!>
!>   dchi/ds   = (chi/Pe) * (dPe/ds)  + Pe  * d(chi/Pe)/ds
!>   dSnu/dtau = (Snu/Bnu) * (dBnu/dT) * (dT/ds) / chi + Bnu * d(Snu/Bnu)/dtau
!>
!> Because the upwind point is horizontally displaced, d/ds is the DIRECTIONAL
!> gradient vec.grad projected from the stored partials (dlnt_dx/dy/dz,
!> dpe_dx/dy/dz); the upwind-point ratios chi/Pe and Snu/Bnu are formed at the
!> upwind grid points and interpolated through the same kernels as chi and Snu,
!> so the kernels act on the smooth ratios, not on separately interpolated steep
!> fields.  The first plane integrates from the boundary ghost plane, which has
!> no beta/Pe, so it keeps short_char's plain linear step.
!===============================================================================
  module subroutine short_char_hermite (self, comp, vec)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    type(ray_type) :: ray, amr_ray
    real(8), intent(in) :: vec(3)
    real(4) :: Inuup, xxup, ssup, xx, ss, Inu, kx, ky, betc, bet, rpe
    real(4) :: dsup1, dsup2, dsup05, rpath, dsc2, Psi4, Irem, dtl, El, Wl, Esc4
    real(4) :: v1, v2, v3, gpe, aT, pec, rpec, rxpc, rxpup
    real(4) :: dchi, t1, sgn, mc, r1mbc, Bc, kf, dsdtau, dss
    real(4) :: rhou, rhoc, gc, t1s, sgns, rdt4, srem, Psi_loc
    real(4) :: ex, U0, U1, q
    real(8) :: ds(3), dt, rdt, A, B, C ! keep in double precision
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:,:,:), allocatable :: plane_r
    real(4), dimension(:,:,:,:), allocatable :: row_red
    real(4), dimension(:), allocatable :: hknu, rem, esc
    real(4), parameter :: dtmin=1e-20
    integer :: xup_t(4,self%nx), yup_t(4,self%ny)
    integer :: zs, ze, xup, yup, zup, i, ii, jj, upo1, upo2
    logical :: amr_boundary, rebuild_tables, get_rates, with_chain, split_rates, off
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char_hermite', itimer=itimer)
    if (.not. (allocated(comp%dpe_dz) .and. allocated(comp%dpe_dx))) &
      error stop 'short_char_hermite: directional composition gradients missing -- set_gradients was not called'
    allocate(work_up(self%fs:self%fe,5))
    allocate(work_min(self%fs:self%fe,5))
    allocate(row_red(self%fs:self%fe, 5, self%nx, self%ny))
    allocate(plane_r(self%fs:self%fe, self%nx, self%ny, 4:5))
    allocate(hknu(self%fs:self%fe))
    allocate(rem(self%fs:self%fe), esc(self%fs:self%fe))
    hknu = real(h_k * self%nu(self%fs:self%fe), 4)
    get_rates = self%compute_psi
    split_rates = self%split_rates

    zs = merge(1, self%nz, self%dir > 0)
    ze = merge(self%nz, 1, self%dir > 0)

    ! directional projection vec.grad -> d/ds along the propagation direction,
    ! per cm; the sign factor makes the projection follow the sweep whichever
    ! way the stored vector points
    v1 = real(vec(1) * self%dir * sign(1.d0, vec(3)) / scaling%l, 4)
    v2 = real(vec(2) * self%dir * sign(1.d0, vec(3)) / scaling%l, 4)
    v3 = real(vec(3) * self%dir * sign(1.d0, vec(3)) / scaling%l, 4)

    call ray%set_stencil(self%ds * scaling%l, vec)
    xkernel1 = real(ray%xkern)
    ykernel1 = real(ray%ykern)
    dsup1 = 0.5*real(ray%dsup)

    amr_boundary = amr/='' .and. self%dir < 0
    if (amr_boundary) then
      ds = self%ds
      if (associated(self%upstream)) ds = self%upstream%ds
      call amr_ray%set_stencil(ds * scaling%l, vec)
      xkernel2 = real(amr_ray%xkern)
      ykernel2 = real(amr_ray%ykern)
      dsup2 = 0.5*real(amr_ray%dsup)
    endif

    rebuild_tables = .true.
    do z=zs, ze, self%dir
      zup = z - self%dir
      with_chain = (z /= zs)  ! the boundary ghost plane has no beta / lnPe

      if (amr_boundary .and. z==zs) then
        xkern => xkernel2
        ykern => ykernel2
        dsup05 = dsup2
        upo1 = amr_ray%upo1 - 2
        upo2 = amr_ray%upo2 - 2
      else
        xkern => xkernel1
        ykern => ykernel1
        dsup05 = dsup1
        upo1 = ray%upo1 - 2
        upo2 = ray%upo2 - 2
      endif
      rpath = 0.5 / dsup05                ! 1 / path length of the step
      dsc2 = dsup05*dsup05 / 3.0          ! (path)^2/12, Euler-Maclaurin term

      if (rebuild_tables) then
        do x=1,self%nx; do ii=1,4
          xup_t(ii,x) = modulo(ii + x + upo1, self%nx) + 1
        end do; end do
        do y=1,self%ny; do jj=1,4
          yup_t(jj,y) = modulo(jj + y + upo2, self%ny) + 1
        end do; end do
        rebuild_tables = amr_boundary .and. z==zs
      endif

      ! upwind-plane smooth ratios: the model assumes chi/Pe and Snu/Bnu vary
      ! linearly, so the kernels interpolate the RATIOS formed at the grid
      ! points, not ratios reconstructed from separately interpolated steep
      ! fields
      if (with_chain) then
        do yup = 1, self%ny
          do xup = 1, self%nx
            rpe = 1.0 / comp%pe4(xup,yup,zup)
            do f = self%fs, self%fe
              bet = self%beta(f,xup,yup,zup)
              plane_r(f,xup,yup,4) = self%work(f,xup,yup,zup,1) * rpe               ! chi/Pe
              plane_r(f,xup,yup,5) = self%work(f,xup,yup,zup,2) * (1.0 - bet)/bet   ! Snu/Bnu
            end do
          end do
        end do
      endif

      !--- ROW PASS: x-reduce the whole upstream plane ----------------------------
      row_red = 0.0
      do yup = 1, self%ny
        do x = 1, self%nx
          do ii = 1, 4
            xup = xup_t(ii,x)
            kx = xkern(ii)
            do i=1, 3
              do f = self%fs, self%fe
                row_red(f,i,x,yup) = row_red(f,i,x,yup) + kx*self%work(f,xup,yup,zup,i)
              end do
            end do
            if (with_chain) then
              do i=4, 5
                do f = self%fs, self%fe
                  row_red(f,i,x,yup) = row_red(f,i,x,yup) + kx*plane_r(f,xup,yup,i)
                end do
              end do
            endif
          end do
        end do
      end do

      do y=1, self%ny
        do x=1, self%nx
          work_up = 0.0
          work_min = huge(1.0)

          ! 4x4 kernel sum
          do jj=1,4
            yup = yup_t(jj,y)
            ky = ykern(jj)
            do i=1, 5
              do f=self%fs,self%fe
                work_up(f,i) = work_up(f,i) + ky*row_red(f,i,x,yup)
              end do
            end do
          end do

          ! 2x2 min over the inner cells
          do jj=2,3
            do ii=2,3
              xup = xup_t(ii,x)
              yup = yup_t(jj,y)
              do i=1, 3
                do f=self%fs,self%fe
                  work_min(f,i) = min(work_min(f,i), self%work(f,xup,yup,zup,i))
                end do
              end do
              if (with_chain) then
                do i=4, 5
                  do f=self%fs,self%fe
                    work_min(f,i) = min(work_min(f,i), plane_r(f,xup,yup,i))
                  end do
                end do
              endif
            end do
          end do
          work_up = max(work_up, work_min)

          if (with_chain) then
            ! per-cell directional gradients at the target
            gpe = v1*comp%dpe_dx(x,y,z) + v2*comp%dpe_dy(x,y,z) + v3*comp%dpe_dz(x,y,z)
            aT  = (v1*comp%dlnt_dx(x,y,z) + v2*comp%dlnt_dy(x,y,z) + v3*comp%dlnt_dz(x,y,z)) &
                  / real(comp%temp(x,y,z), 4)
            pec  = comp%pe4(x,y,z)
            rpec = 1.0 / pec

            do f=self%fs, self%fe
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              xxup  = work_up(f,1)
              ssup  = work_up(f,2)
              Inuup = work_up(f,3)
              rxpup = work_up(f,4)      ! kernel-interpolated upwind chi/Pe
              rhou  = work_up(f,5)      ! kernel-interpolated upwind Snu/Bnu
              betc  = self%beta(f,x,y,z)

              ! target-end dchi/ds by the product rule chi = (chi/Pe)*Pe with
              ! chi/Pe linear along the ray (one-sided ratio difference), then
              ! the MC positivity envelope of the interval secant
              dchi = xx - xxup
              t1 = dchi * rpath
              sgn = sign(1.0, t1)
              rxpc = xx * rpec
              mc = rxpc*gpe + pec*(rxpc - rxpup)*rpath
              mc = sgn * max(0.0, min(sgn*mc, 2.0*abs(t1)))
              ! maser interval: keep the sign of dtau, bound it away from zero
              ! and cap the amplification
              dt = dsup05*(xxup + xx) + dsc2*(t1 - mc)
              dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
              rdt = 1.0d0 / dt
              rdt4 = real(rdt, 4)
              call parab_moments(real(dt,4), ex, U0, U1, q)

              ! target-end dSnu/dtau by the product rule S = rho*Bnu
              ! (rho = S/Bnu):  dS/dtau = rho*(dBnu/dtau) + Bnu*(drho/dtau).  The
              ! analytic dBnu/dtau = hk*B*r1mb*aT/chi enters as rho*dBnu/dtau, so
              ! its B factor cancels rho's 1/B: gc holds (dBnu/dtau)/Bnu.  The
              ! departure part is the one-sided ratio secant over the interval,
              ! clamped into the MC positivity envelope of the S secant
              r1mbc = 1.0 / (1.0 - betc)
              Bc = betc * r1mbc
              gc = hknu(f) * r1mbc * aT / xx         ! (dBnu/dtau)/Bnu at c
              rhoc = ss / Bc
              dss = ss - ssup
              t1s = dss * rdt4                       ! interval secant dS/dtau
              sgns = sign(1.0, t1s)
              dsdtau = ss*gc + Bc*(rhoc - rhou)*rdt4
              ! Resolve the MC positivity clamp by branch rather than by min/max,
              ! so that the slope in ss and the part that does not depend on it
              ! are both known exactly.  kf is d(dsdtau)/dss on the branch taken,
              ! srem is what is left of dsdtau once kf*ss is removed, written out
              ! directly: Bc*rhoc*rdt cancels ss*rdt analytically, so the pair is
              ! never formed as a difference of two nearly equal numbers.
              if (sgns*dsdtau <= 0.0) then
                dsdtau = 0.0                         ! clamped to zero slope
                kf   = 0.0
                srem = 0.0
              elseif (sgns*dsdtau >= 2.0*abs(t1s)) then
                dsdtau = 2.0*t1s                     ! clamped to twice the secant
                kf   = 2.0*rdt4
                srem = -2.0*rdt4*ssup
              else
                kf   = gc + rdt4                     ! product-rule slope, unclamped
                srem = -Bc*rhou*rdt4
              endif

              B = (2.0d0*U1 - q) * rdt
              A = U0 - B
              C = q - U1
              Inu = Inuup*ex + real(A,4)*ssup + real(B,4)*ss + real(C,4)*dsdtau
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                ! Lambda diagonal from the product-rule slope on the branch the
                ! clamp selected.  Outside [0,1] the local term is not a weight
                ! and the scheme has no operator to precondition with, so drop it
                ! there, as for dt<=0 and chi<=0.
                ! rem is the non-local part the split rates need; it is formed
                ! unconditionally and copied out after the loop, which keeps
                ! the rate_formulation branch out of this one.
                Psi_loc = real(B + C*dble(kf), 4)
                if (dt > 0 .and. xx > 0 .and. Psi_loc >= 0. .and. Psi_loc <= 1. &
                    .and. abs(Inu) < real(maser_Inu_max)) then
                  self%Psi(f,x,y,z) = Psi_loc / xx
                  rem(f) = Inuup*ex + real(A,4)*ssup + real(C,4)*srem
                  ! 1 - Psi_loc, formed from the real(8) weights before rounding
                  esc(f) = real(1d0 - (B + C*dble(kf)), 4)
                else
                  ! Psi=0 switches the preconditioning off at this point, so Irem
                  ! then has to carry the whole (clamped) intensity.
                  self%Psi(f,x,y,z) = 0
                  rem(f) = Inu
                  esc(f) = 1.
                endif
              endif
            end do ! f-loop
            if (get_rates .and. split_rates) then
              do f=self%fs, self%fe
                self%Irem(f,x,y,z) = rem(f)
                self%escape(f,x,y,z) = esc(f)
              end do
            endif
          else
            ! boundary plane: plain linear step, as in short_char.
            ! One f-loop per rate formulation, so that neither carries the other's
            ! branch.  Both are branch-free apart from selects, and vectorise.
            if (split_rates) then
              do f=self%fs, self%fe
                xxup  = work_up(f,1)
                ssup  = work_up(f,2)
                Inuup = work_up(f,3)
                xx = self%work(f,x,y,z,1)
                ss = self%work(f,x,y,z,2)
                dtl = dsup05 * ( xxup + xx )     ! upwind optical pathlength
                dtl = max(dtl, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dtl<0)
                call linear_weights(dtl, El, Wl, Psi4, Esc4)
                Inu = Inuup*El + Wl*ssup + Psi4*ss
                ! Hold the non-local part and the escape factor 1-Psi4 on their own
                ! as well: the preconditioned rates need them, and both Inu - Psi4*ss
                ! and 1 - Psi4 would cancel in an optically thick interval.  Psi4=0
                ! switches the preconditioning off at a point, so Irem then carries
                ! the whole (clamped) intensity.
                off = dtl <= 0 .or. xx <= 0 .or. abs(Inu) > real(maser_Inu_max)
                Irem = Inuup*El + Wl*ssup
                Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
                Psi4 = merge(0., Psi4, off)
                self%Irem(f,x,y,z) = merge(Inu, Irem, off)
                self%escape(f,x,y,z) = merge(1., Esc4, off)
                self%work(f,x,y,z,3) = Inu
                if (get_rates) then
                  ! 0<=Psi4<=1 fails for dtl<0/chi<0; no preconditioning there
                  self%Psi(f,x,y,z) = merge(Psi4 / xx, 0., dtl > 0 .and. xx > 0)
                endif
              end do ! f-loop
            else
              do f=self%fs, self%fe
                xxup  = work_up(f,1)
                ssup  = work_up(f,2)
                Inuup = work_up(f,3)
                xx = self%work(f,x,y,z,1)
                ss = self%work(f,x,y,z,2)
                dtl = dsup05 * ( xxup + xx )     ! upwind optical pathlength
                dtl = max(dtl, -real(maser_dtau_max,4))  ! cap maser amplification (chi<0 -> dtl<0)
                call linear_weights(dtl, El, Wl, Psi4, Esc4)
                Inu = Inuup*El + Wl*ssup + Psi4*ss
                Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
                self%work(f,x,y,z,3) = Inu
                if (get_rates) then
                  ! 0<=Psi4<=1 fails for dtl<0/chi<0; no preconditioning there
                  self%Psi(f,x,y,z) = merge(Psi4 / xx, 0., dtl > 0 .and. xx > 0)
                endif
              end do ! f-loop
            endif
          endif
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE short_char_hermite

END SUBMODULE rt_solvers_smod
