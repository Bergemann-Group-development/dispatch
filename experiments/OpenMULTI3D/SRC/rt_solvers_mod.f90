!===============================================================================
!> Formal solvers for rt_patch_t: split out of rt_patch_mod as a submodule so the
!> RT physics (short/long characteristics, linear/quadratic/hermite) lives apart
!> from patch bookkeeping, while keeping full host association with rt_patch_mod
!> (types, constants, USE'd symbols) and unchanged self%-based call syntax.
!===============================================================================
SUBMODULE (rt_patch_mod) rt_solvers_smod
  implicit none
contains

  !===============================================================================
  module subroutine short_char (self, vec)
    class(rt_patch_t) :: self
    type(ray_type) :: ray, amr_ray
    real(8), intent(in) :: vec(3)
    real(4) :: Inuup, xxup, ssup, xxmin, ssmin, Inumin
    real(4) :: dsup1, dsup2, dsup05, xx, ss, Inu, u0, Psi, kx, ky
    real(8) :: ds(3), dt, ex ! important
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:,:,:), allocatable :: row_red
    integer :: xup_t(4,self%nx), yup_t(4,self%ny)
    integer :: zs, ze, xup, yup, zup, x1, y1, i, ii, jj, upo1, upo2
    logical :: inner, amr_boundary, rebuild_tables, get_rates
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char', itimer=itimer)
    allocate(work_up(self%fs:self%fe,3))
    allocate(work_min(self%fs:self%fe,3))
    allocate(row_red(self%fs:self%fe, 3, self%nx, self%ny))
    get_rates = any(self%get_rates)

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

          do f=self%fs, self%fe
            xx = self%work(f,x,y,z,1)
            ss = self%work(f,x,y,z,2)
            xxup  = work_up(f,1)
            ssup  = work_up(f,2)
            Inuup = work_up(f,3)

            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            ! An inverted population gives chi<0, hence dt<0 and exp(+|dt|)
            ! amplification.  Cap it, so that a maser cannot overflow inside a
            ! single formal solution; the population feedback then gets a chance
            ! to quench the inversion between iterations.
            dt = max(dt, -maser_dtau_max)
            if (dt < 80) then
              ex = exp(-dt)
              ! the series expansion is only valid for |dt| << 1, and dt<0 must not
              ! fall into it
              if ( abs(dt) > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif
            Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))

            self%work(f,x,y,z,3) = Inu
            if (get_rates) then
              if (dt > 0 .and. xx > 0) then
                self%Psi(f,x,y,z) = Psi / xx
              else
                ! Rybicki & Hummer operator splitting assumes 0 <= Psi <= 1, which
                ! is what makes the preconditioned rate matrix well behaved.  For
                ! dt<0, Psi is negative and exponentially large, and for chi<0 the
                ! correction would change sign; leave the rates unpreconditioned.
                self%Psi(f,x,y,z) = 0
              endif
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE short_char

  !===============================================================================
  module subroutine long_char_linear_spherical (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, Inuup, xxup, ssup, dsup05, dsup1, dsup2, xx, ss, Inu, u0, Psi
    real(8) :: dt, ex ! important
    real(4) :: mu, p, lower, ds1, ds2, r1, r2, s1, s2
    integer :: dx, dy, zs, ze, zup
    logical :: amr_boundary
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_linear_spherical', itimer=itimer)

    p = rstar * sqrt(1.d0 - vec(3)**2) ! impact parameter
    lower = (self%position(3) - 0.5*self%size(3))*scaling%l + rorigin
    ds1 = self%ds(3) * scaling%l

    amr_boundary = amr/='' .and. self%dir < 0 .and. associated(self%upstream)
    if (amr_boundary) then
      ds2 = self%upstream%ds(3) * scaling%l
    endif

    zs = merge(1, self%nz, self%dir > 0)
    ze = merge(self%nz, 1, self%dir > 0)

    do z=zs, ze, self%dir
      zup = z - self%dir
      r2 = lower + z*ds1

      if (amr_boundary) then
        r1 = r2 - self%dir*ds2
        amr_boundary = .false.
      else
        r1 = r2 - self%dir*ds1
      endif

      if (r2 < p .or. r1 < p) then
        self%work(:,:,:,z,3) = 0
        cycle
      endif

      s1 = sqrt(r1**2 - p**2)
      s2 = sqrt(r2**2 - p**2)
      dsup05 = 0.5d0 * abs(s2-s1)

      do y=1, self%ny
        do x=1, self%nx
          do f=self%fs, self%fe
            xxup = self%work(f,x,y,zup,1)
            ssup = self%work(f,x,y,zup,2)
            xx = self%work(f,x,y,z,1)
            ss = self%work(f,x,y,z,2)
            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            ! An inverted population gives chi<0, hence dt<0 and exp(+|dt|)
            ! amplification; cap it (see short_char), rather than clamping to 0,
            ! which would pass the maser through untouched
            dt = max(dt, -maser_dtau_max)
            if (dt < 80) then
              Inuup = self%work(f,x,y,zup,3)
              ex = exp(-dt)
              ! the series expansion is only valid for |dt| << 1, and dt<0 must not
              ! fall into it
              if ( abs(dt) > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif
            Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))

            self%work(f,x,y,z,3) = Inu
            if (any(self%get_rates)) then
              if (dt > 0 .and. xx > 0) then
                self%Psi(f,x,y,z) = Psi / xx
              else
                ! 0 <= Psi <= 1 no longer holds; leave the rates unpreconditioned
                self%Psi(f,x,y,z) = 0
              endif
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear_spherical

  !===============================================================================
  module subroutine long_char_linear (self, vec)
    class(rt_patch_t) :: self
    real(8) :: vec(3)
    real(4) :: ds, Inuup, xxup, ssup, dsup05, dsup1, dsup2, xx, ss, Inu, u0, Psi
    real(8) :: dt, ex ! important
    integer :: dx, dy, zs, ze, zup
    logical :: amr_boundary
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_linear', itimer=itimer)

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
          do f=self%fs, self%fe
            xxup = self%work(f,x,y,zup,1)
            ssup = self%work(f,x,y,zup,2)
            xx = self%work(f,x,y,z,1)
            ss = self%work(f,x,y,z,2)
            dt = dsup05 * ( xxup + xx ) ! upwind optical pathlength
            ! An inverted population gives chi<0, hence dt<0 and exp(+|dt|)
            ! amplification; cap it (see short_char), rather than clamping to 0,
            ! which would pass the maser through untouched
            dt = max(dt, -maser_dtau_max)
            if (dt < 80) then
              Inuup = self%work(f,x,y,zup,3)
              ex = exp(-dt)
              ! the series expansion is only valid for |dt| << 1, and dt<0 must not
              ! fall into it
              if ( abs(dt) > 1d-3 ) then
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              !ex = 0
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif
            Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))

            self%work(f,x,y,z,3) = Inu
            if (any(self%get_rates)) then
              !Psi = Psi + Inuup*ex*exp(-dsup05*xx*ss)*dsup05*xx
              if (dt > 0 .and. xx > 0) then
                self%Psi(f,x,y,z) = Psi / xx
              else
                ! 0 <= Psi <= 1 no longer holds; leave the rates unpreconditioned
                self%Psi(f,x,y,z) = 0
              endif
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear

  !===============================================================================
  !> Parabolic transfer moments over an interval of optical depth dt (with
  !> ex=exp(-dt), rdt=1/dt): U0 = 1-exp(-dt), U1 = dt-U0, q = U2/dt with
  !> U2 = dt^2-2*U1.  Below dt=1e-3 they are taken from their Taylor series
  !> (U1 = dt^2/2 - dt^3/6, q = dt^2/3 - dt^3/12) to avoid cancellation in 1-ex.
  !> Shared by every parabolic-weight formal solver; kept a leaf pure routine so
  !> -Ofast inlines it into the hot f-loop just as the in-place copies were.
  !===============================================================================
  pure module subroutine parab_moments (dt, ex, rdt, U0, U1, q)
    real(8), intent(in)  :: dt, ex, rdt
    real(8), intent(out) :: U0, U1, q
    !-----------------------------------------------------------------------------
    ! the Taylor branch is only valid for |dt| << 1; a maser interval (dt<0,
    ! bounded by maser_dtau_max) must use the exact moments, which hold for
    ! either sign of dt
    if (abs(dt) > 1e-3) then
      U0 = 1.0d0 - ex
      U1 = dt - U0
      q  = dt - 2.0d0*U1*rdt        ! q = U2/dt
    else
      U1 = dt*dt*(0.5d0 - dt*one_sixth)
      U0 = dt - U1
      q  = dt*dt*(one_third - dt*one_twelfth)
    end if
  end subroutine parab_moments

!===============================================================================
!> Quadratic formal solver, equivalent to the 'bezier' Sinterpol case of Multi3D's
!> shortchar_3D.f90.  The parabola over the integrated interval is written with the
!> limited slope dS/dtau at one end, giving weights (A, B, C) on (S_upwind, S_target,
!> dS/dtau) in terms of the moments
!>
!>   U0 = 1-exp(-dt),  U1 = dt-U0,  U2 = dt^2-2*U1,  q = U2/dt
!>
!> Writing the same parabola in Bezier form, S(u) = Sup*(1-u)^2 + 2*Sc*u*(1-u) +
!> S0*u^2, the control point Sc = S0 -/+ dt/2*dS/dtau stays inside [Sup,S0] -- i.e.
!> the parabola cannot overshoot and drive Inu or the Lambda diagonal negative --
!> exactly when dS/dtau has the sign of the one-sided slope across the integrated
!> interval and is no larger than twice its magnitude.  The monotonized-central
!> limiter below satisfies both bounds unconditionally (it never exceeds
!> 2*min(|t1|,|t3|) and vanishes at an extremum), so no explicit control point needs
!> to be formed or range-checked.
!>
!> The z==ze plane has no downwind neighbour, so its stencil is shifted one cell
!> upwind: the interval integrated is z0 -> z, and dS/dtau then sits at the UPWIND
!> end of it rather than at the target end, which swaps the roles of the weights.
!> It is peeled out of the z-loop to keep the hot loop free of the branch.
!===============================================================================
  module subroutine long_char_quadratic (self, vec)
    class(rt_patch_t), target :: self
    real(8) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2
    real(4) :: Inu, ppsi
    real(4), parameter :: dtmin=1e-20
    ! reached through a polymorphic dummy, self%work / self%Psi force gfortran to
    ! reload the array descriptor inside the f-loop, and it then cannot prove the
    ! base address is affine ("evolution of base is not affine") -- which blocks
    ! vectorisation.  Hoisting them, and the loop bounds, into local pointers and
    ! plain integers restores unit stride.
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
    get_rates = any(self%get_rates)
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
    !> at_target the interval integrated is up -> c and dS/dtau is taken at c;
    !> otherwise it is c -> dn and dS/dtau is taken at c, which is then the upwind
    !> end.  Either way Inu comes out at the downwind end of the integrated
    !> interval, and plam is the Lambda diagonal there.  Inlined by the compiler,
    !> so at_target folds away and the f-loop stays free of the branch.
    !---------------------------------------------------------------------------
    pure subroutine quad_step (xxup, xxc, xxdn, ssup, ssc, ssdn, Inuup, at_target, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xxdn, ssup, ssc, ssdn, Inuup
      logical, intent(in) :: at_target
      real(4), intent(out) :: Inu, plam
      real(4) :: t1, t2, t3, dtup, dtdn, dsdtau
      real(8) :: dt, rdt, ex, U0, U1, q, A, B, C ! important -- keep in double precision
      !-------------------------------------------------------------------------
      ! one sided dchi/ds either side of c, and the limited central slope
      t1 = (xxc - xxup) * rdsup
      t3 = (xxdn - xxc) * rdsup
      t2 = 0.5*(t1 + t3)
      t2 = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))
      ! a maser interval (chi<0) gives dtau<0: keep the sign, but hold |dtau|
      ! away from zero (it divides the slopes below) and cap the amplification
      ! at maser_dtau_max, consistent with the other solvers
      dtup = dsup1*(xxup + xxc) + dsup2*(t1-t2)
      dtdn = dsup1*(xxdn + xxc) + dsup2*(t2-t3)
      dtup = max(sign(max(abs(dtup), dtmin), dtup), -real(maser_dtau_max))
      dtdn = max(sign(max(abs(dtdn), dtmin), dtdn), -real(maser_dtau_max))

      dt = merge(dtup, dtdn, at_target)
      rdt = 1.0d0 / dt
      ex = exp(-min(dt, 80.d0))
      call parab_moments(dt, ex, rdt, U0, U1, q)

      ! one sided dS/dtau either side of c, and the limited central slope
      if (at_target) then
        t1 = (ssc - ssup) * real(rdt,4)
        t3 = (ssdn - ssc) / dtdn
      else
        t1 = (ssc - ssup) / dtup
        t3 = (ssdn - ssc) * real(rdt,4)
      endif
      t2 = 0.5*(t1 + t3)
      dsdtau = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

      ! for dt<0 (bounded maser) or a negative target chi the Lambda diagonal
      ! leaves 0 <= Psi <= 1, so the rates are left unpreconditioned (plam=0),
      ! and Inu is capped, exactly as in the other solvers
      plam = 0.0
      if (at_target) then
        B = (2.0d0*U1 - q) * rdt      ! on ssc  (the target)
        A = U0 - B                    ! on ssup (the upwind point)
        C = q - U1                    ! on dsdtau
        Inu = Inuup*ex + real(A,4)*ssup + real(B,4)*ssc + real(C,4)*dsdtau
        if (get_rates .and. dt > 0 .and. xxc > 0) &
          plam = real(B,4) / xxc    ! diagonal of the Lambda operator
      else
        A = q * rdt                   ! on ssdn (the target)
        B = U0 - A                    ! on ssc  (the upwind point)
        C = U1 - q                    ! on dsdtau
        Inu = Inuup*ex + real(A,4)*ssdn + real(B,4)*ssc + real(C,4)*dsdtau
        if (get_rates .and. dt > 0 .and. xxdn > 0) &
          plam = real(A,4) / xxdn   ! diagonal of the Lambda operator
      endif
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
    end subroutine quad_step
  END SUBROUTINE long_char_quadratic

!===============================================================================
!> Chain-rule cubic-Hermite formal solver ('hermite'), long-characteristics
!> version.  The source function over each integrated interval [up, c] is a
!> cubic Hermite in optical depth with slopes at BOTH ends, and both the chi and
!> Snu slopes are obtained by the chain rule on smoother variables instead of
!> from limited finite differences of chi and Snu themselves:
!>
!>   dchi/ds  = (dchi/dlnPe) * (dlnPe/ds)
!>   dSnu/dtau = (dSnu/dBnu) * (dBnu/dT) * (dT/ds) / chi
!>
!> dlnPe/ds and dT/ds are the PCHIP-interpolant gradients of the input snapshot
!> (comp%dlnpe_dz, comp%dlnt_dz -- for rotated compositions the z-columns run
!> along the ray, and for 1D/column runs the characteristic stays in the
!> column, so the z-gradient is all that is needed); dBnu/dT is analytic in the
!> precomputed beta = exp(-h nu/kT); and dchi/dlnPe, dSnu/dBnu are one-sided
!> linear estimates across the integrated interval, which keeps the result
!> invariant under domain decomposition.  Because the pchip interpolant is
!> monotone, dlnPe/ds carries the sign of the interval's dlnPe, so the chain
!> slope automatically has the sign of the interval secant; the clamps below
!> only bound its magnitude.
!>
!> Writing the source cubic in Bezier form, its two control points are clamped
!> into [min(Sup,Sc), max(Sup,Sc)], which keeps Inu and the Lambda diagonal
!> positive (the cubic analogue of the quadratic solver's MC limiter) while
!> passing the chain-rule slopes through untouched in smooth regions.  The
!> weights are integrals of the Bernstein basis against exp(-t), expressed in
!> the moments Mn = int_0^dt t^n e^-t dt with a Taylor branch at small dt.
!> dtau itself is the trapezoid with the Euler-Maclaurin correction built from
!> the chain-rule chi slopes at the two ends.
!>
!> The first plane of the sweep integrates from the boundary ghost plane, which
!> has chi/Snu/Inu (set_upstream_boundary) via halo exchange.  When a genuine
!> upstream neighbour patch exists (comp_up, passed in from radiative_transfer),
!> its composition and back-ground opacity give the real thermodynamic state at
!> the plane coincident with the ghost, so the interval is folded into the same
!> herm_step used for interior planes.  Only at a true domain boundary (no
!> neighbour, comp_up unassociated) does the upwind-end source slope degrade to
!> the interval secant (herm_first), borrowing the adjacent downwind interval
!> for the state-space slopes since no thermodynamic data exists upstream at all.
!===============================================================================
  module subroutine long_char_hermite (self, comp, vec, comp_up)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(composition_type), pointer :: comp_up
    real(8), intent(in) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2, gfac
    real(4) :: Inu, ppsi, rdlp, cpeu, cpec, aTu, aTc, rTup, betup
    real(4), parameter :: dtmin=1e-20, eps_pe=1e-10, eps_B=1e-25
    real(4), dimension(:,:,:,:,:), contiguous, pointer :: w, bk
    real(4), dimension(:,:,:,:), contiguous, pointer :: psi
    real(4), dimension(:), allocatable :: hknu
    integer :: zs, ze, zup, zdn, zu, dir, nx, ny, fs, fe
    logical :: get_rates, use_upstream
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_hermite', itimer=itimer)
    if (.not. allocated(comp%dlnpe_dz)) &
      stop 'long_char_hermite: composition gradients missing -- set_gradients was not called'
    ds = self%ds(3) * scaling%l
    dsup = ds / abs(vec(3))
    rdsup = 1./ dsup
    dsup1 = 0.5*dsup
    dsup2 = 1./12.*dsup*dsup     ! Euler-Maclaurin correction to the trapezoid
    ! converts the stored d/dz (per code length unit) to d/ds along propagation
    gfac = real(self%dir * abs(vec(3)) / scaling%l, 4)
    get_rates = any(self%get_rates)
    dir = self%dir
    nx = self%nx
    ny = self%ny
    fs = self%fs
    fe = self%fe
    zs = merge(1, self%nz, dir > 0)
    ze = merge(self%nz, 1, dir > 0)
    w => self%work
    bk => self%back
    psi => null()
    if (get_rates) psi => self%Psi
    allocate(hknu(fs:fe))
    hknu = real(h_k * self%nu(fs:fe), 4)

    ! first plane: integrate from the boundary ghost plane
    use_upstream = associated(comp_up)
    z = zs
    zup = z - dir
    zdn = z + dir
    if (use_upstream) zu = merge(self%upstream%nz, 1, dir > 0)
    do y=1, ny
      do x=1, nx
        if (use_upstream) then
          rdlp = comp%lnpe(x,y,z) - comp_up%lnpe(x,y,zu)
          rdlp = 1.0 / (rdlp + sign(eps_pe, rdlp))
          cpeu = gfac * comp_up%dlnpe_dz(x,y,zu) * rdlp
          cpec = gfac * comp%dlnpe_dz(x,y,z)     * rdlp
          rTup = 1.0 / real(comp_up%temp(x,y,zu), 4)   ! upwind-plane beta from T
          aTu  = gfac * comp_up%dlnt_dz(x,y,zu) * rTup
          aTc  = gfac * comp%dlnt_dz(x,y,z)     / real(comp%temp(x,y,z), 4)
          do f=fs, fe
            ! the ghost plane carries chi/Snu/Inu (set_upstream_boundary) but not
            ! beta -- the upstream patch's own back-array is long deallocated by
            ! now, so beta at the upwind point is rebuilt from the neighbour's T,
            ! exactly as lte_background_opacities forms it (max(exp(-hnu/kT),tiny))
            betup = max(exp(-hknu(f)*rTup), tiny(1.0))
            call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), &
                            betup, bk(f,x,y,z,3), hknu(f), &
                            w(f,x,y,zup,3), cpeu, cpec, aTu, aTc, Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        else
          rdlp = comp%lnpe(x,y,zdn) - comp%lnpe(x,y,z)
          rdlp = 1.0 / (rdlp + sign(eps_pe, rdlp))
          cpec = gfac * comp%dlnpe_dz(x,y,z) * rdlp
          aTc  = gfac * comp%dlnt_dz(x,y,z) / real(comp%temp(x,y,z), 4)
          do f=fs, fe
            call herm_first (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                             w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                             bk(f,x,y,z,3), bk(f,x,y,zdn,3), hknu(f), &
                             w(f,x,y,zup,3), cpec, aTc, Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        endif
      end do
    end do

    ! interior planes: everything comes from the integrated interval [zup, z]
    do z=zs+dir, ze, dir
      zup = z - dir
      do y=1, ny
        do x=1, nx
          rdlp = comp%lnpe(x,y,z) - comp%lnpe(x,y,zup)
          rdlp = 1.0 / (rdlp + sign(eps_pe, rdlp))
          cpeu = gfac * comp%dlnpe_dz(x,y,zup) * rdlp
          cpec = gfac * comp%dlnpe_dz(x,y,z)   * rdlp
          aTu  = gfac * comp%dlnt_dz(x,y,zup) / real(comp%temp(x,y,zup), 4)
          aTc  = gfac * comp%dlnt_dz(x,y,z)   / real(comp%temp(x,y,z),   4)
          do f=fs, fe
            call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), &
                            bk(f,x,y,zup,3), bk(f,x,y,z,3), hknu(f), &
                            w(f,x,y,zup,3), cpeu, cpec, aTu, aTc, Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  contains
    !---------------------------------------------------------------------------
    !> Cubic-Bezier weights: integrals of the Bernstein basis against exp(-t)
    !> over the interval, t measured backward from the target.  w0 multiplies
    !> S_up, w1/w2 the upwind-/target-side control points, w3 S_target;
    !> w0+w1+w2+w3 = 1-exp(-dt) exactly (enforced by construction of w3).
    !---------------------------------------------------------------------------
    pure subroutine bezier3_weights (dt, ex, w0, w1, w2, w3)
      real(8), intent(in) :: dt, ex
      real(8), intent(out) :: w0, w1, w2, w3
      real(8) :: M0, M1, M2, M3, r3
      !-------------------------------------------------------------------------
      ! the Taylor branch is only valid for |dt| << 1; a maser interval (dt<0,
      ! bounded by maser_dtau_max) must use the exact moments, which hold for
      ! either sign of dt
      if (abs(dt) > 0.05d0) then
        M0 = 1.d0 - ex
        M1 = 1.d0 - ex*(1.d0 + dt)
        M2 = 2.d0 - ex*(2.d0 + dt*(2.d0 + dt))
        M3 = 6.d0 - ex*(6.d0 + dt*(6.d0 + dt*(3.d0 + dt)))
        r3 = 1.d0 / (dt*dt*dt)
        w0 = M3*r3
        w1 = 3.d0*(M2*dt - M3)*r3
        w2 = 3.d0*(M1*dt*dt - 2.d0*M2*dt + M3)*r3
        w3 = M0 - w0 - w1 - w2
      else
        ! cancellation in the moment combinations grows with the moment order,
        ! so below dt=0.05 the weights are evaluated from their Taylor series
        w0 = dt*(0.25d0 - dt*(0.20d0 - dt*(1.d0/12.d0  - dt/42.d0)))
        w1 = dt*(0.25d0 - dt*(0.15d0 - dt*(1.d0/20.d0  - dt/84.d0)))
        w2 = dt*(0.25d0 - dt*(0.10d0 - dt*(1.d0/40.d0  - dt/210.d0)))
        w3 = dt*(0.25d0 - dt*(0.05d0 - dt*(1.d0/120.d0 - dt/840.d0)))
      endif
    end subroutine bezier3_weights
    !---------------------------------------------------------------------------
    !> One interior interval [up, c].  Chain-rule chi slopes at both ends set
    !> dtau; chain-rule source slopes at both ends set the Bezier control
    !> points, which are then clamped into [min,max](Sup,Sc).
    !---------------------------------------------------------------------------
    pure subroutine herm_step (xxup, xxc, ssup, ssc, betau, betac, hk, Inuup, &
                               cpeu, cpec, aTu, aTc, Inu, plam)
      real(4), intent(in) :: xxup, xxc, ssup, ssc, betau, betac, hk, Inuup
      real(4), intent(in) :: cpeu, cpec, aTu, aTc
      real(4), intent(out) :: Inu, plam
      real(4) :: dchi, t1, sgn, mup, mc, r1mbu, r1mbc, Bu, Bc, rdB, gu, gc
      real(4) :: dt3, smin, smax, P1, P2, dss, d1, d2
      real(8) :: dt, ex, w0, w1, w2, w3
      !-------------------------------------------------------------------------
      ! chain-rule dchi/ds at both ends, clamped to the sign of the interval
      ! secant and twice its magnitude (the MC positivity envelope)
      dchi = xxc - xxup
      t1 = dchi * rdsup
      sgn = sign(1.0, t1)
      mup = sgn * max(0.0, min(sgn*(dchi*cpeu), 2.0*abs(t1)))
      mc  = sgn * max(0.0, min(sgn*(dchi*cpec), 2.0*abs(t1)))
      ! a maser interval (chi<0) gives dtau<0: keep the sign, but hold |dtau|
      ! away from zero and cap the amplification, consistent with the other
      ! solvers
      dt = dsup1*(xxup + xxc) + dsup2*(mup - mc)
      dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
      ex = exp(-min(dt, 80.d0))
      call bezier3_weights (dt, ex, w0, w1, w2, w3)

      ! chain-rule dSnu/dtau at both ends, written as dss*g so the Lambda
      ! diagonal can carry the exact sensitivity to the local S: dSnu/dBnu
      ! estimated across the interval, dBnu/dT analytic in beta, dT/ds from
      ! the snapshot interpolant
      r1mbu = 1.0 / (1.0 - betau)
      r1mbc = 1.0 / (1.0 - betac)
      Bu = betau * r1mbu
      Bc = betac * r1mbc
      dss = ssc - ssup
      rdB = Bc - Bu
      rdB = 1.0 / (rdB + sign(eps_B, rdB))
      gu = hk * Bu * r1mbu * aTu / xxup * rdB
      gc = hk * Bc * r1mbc * aTc / xxc  * rdB

      dt3 = real(dt,4) * one_third
      smin = min(ssup, ssc)
      smax = max(ssup, ssc)
      P1 = min(max(ssup + dt3*dss*gu, smin), smax)
      P2 = min(max(ssc  - dt3*dss*gc, smin), smax)

      Inu = Inuup*real(ex,4) + real(w0,4)*ssup + real(w1,4)*P1 &
                             + real(w2,4)*P2   + real(w3,4)*ssc
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      plam = 0.0
      if (get_rates .and. dt > 0 .and. xxc > 0) then
        ! Lambda diagonal consistent with the control points (clamp neglected):
        ! dP1/dSc = dt3*gu, dP2/dSc = 1 - dt3*gc, both bounded into [0,1] so
        ! the diagonal can never exceed the attenuated total or go negative.
        ! For dt<0 or chi<0 (bounded maser) 0 <= Psi <= 1 no longer holds, and
        ! the rates are left unpreconditioned (plam=0), as in the other solvers.
        d1 = min(max(dt3*gu, 0.0), 1.0)
        d2 = min(max(1.0 - dt3*gc, 0.0), 1.0)
        plam = real(w1*dble(d1) + w2*dble(d2) + w3, 4) / xxc
      endif
    end subroutine herm_step
    !---------------------------------------------------------------------------
    !> The boundary interval [ghost, zs]: state-space slopes from the adjacent
    !> interior interval [zs, zs+dir] (the ghost plane has no beta/lnPe), and
    !> the upwind-end source slope degrades to the interval secant.
    !---------------------------------------------------------------------------
    pure subroutine herm_first (xxup, xxc, xxdn, ssup, ssc, ssdn, betac, betadn, &
                                hk, Inuup, cpec, aTc, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xxdn, ssup, ssc, ssdn, betac, betadn, hk, Inuup
      real(4), intent(in) :: cpec, aTc
      real(4), intent(out) :: Inu, plam
      real(4) :: t1, sgn, mc, r1mbc, r1mbdn, Bc, Bdn, FS, dSc, dSu
      real(4) :: dt3, smin, smax, P1, P2, dss
      real(8) :: dt, ex, w0, w1, w2, w3
      !-------------------------------------------------------------------------
      t1 = (xxc - xxup) * rdsup
      sgn = sign(1.0, t1)
      mc = sgn * max(0.0, min(sgn*((xxdn - xxc)*cpec), 2.0*abs(t1)))
      ! maser interval: keep the sign of dtau, bound it away from zero (it
      ! divides the secant slope below) and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(t1 - mc)
      dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
      ex = exp(-min(dt, 80.d0))
      call bezier3_weights (dt, ex, w0, w1, w2, w3)

      r1mbc  = 1.0 / (1.0 - betac)
      r1mbdn = 1.0 / (1.0 - betadn)
      Bc  = betac  * r1mbc
      Bdn = betadn * r1mbdn
      FS = Bdn - Bc
      FS = (ssdn - ssc) / (FS + sign(eps_B, FS))
      dSc = FS * hk * Bc * r1mbc * aTc / xxc
      dss = ssc - ssup
      dSu = dss * real(1.d0/dt, 4)      ! interval secant slope wrt tau

      dt3 = real(dt,4) * one_third
      smin = min(ssup, ssc)
      smax = max(ssup, ssc)
      P1 = min(max(ssup + dt3*dSu, smin), smax)
      P2 = min(max(ssc  - dt3*dSc, smin), smax)

      Inu = Inuup*real(ex,4) + real(w0,4)*ssup + real(w1,4)*P1 &
                             + real(w2,4)*P2   + real(w3,4)*ssc
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      plam = 0.0
      ! P1's secant slope contributes dP1/dSc = 1/3; P2's sensitivity is >= 1
      ! and is bounded to 1, as in herm_step; for dt<0 or chi<0 the rates are
      ! left unpreconditioned, as in the other solvers
      if (get_rates .and. dt > 0 .and. xxc > 0) &
        plam = real(w1*(1.d0/3.d0) + w2 + w3, 4) / xxc
    end subroutine herm_first
  END SUBROUTINE long_char_hermite

!===============================================================================
!> Chain-rule hermite variant of the 3D short-characteristics solver
!> (short_solver='hermite' on the true-3D iteration path).  The 4x4 upwind
!> kernel interpolation and its 2x2 min limiter are exactly those of
!> short_char; the linear integration step is replaced by quad_step's weight
!> algebra with the chain-rule slope at the target cell.  Because the upwind
!> point is horizontally displaced, d/ds is the DIRECTIONAL gradient
!> vec.grad projected from the stored partials (dlnt_dx/dy/dz, dlnpe_dx/dy/dz
!> of the unrotated composition); lnPe and beta at the upwind point come
!> through the same interpolation kernels as chi and Snu.  The first plane of
!> the sweep integrates from the boundary ghost plane, which has no beta/lnPe:
!> it keeps short_char's plain linear step.
!===============================================================================
  module subroutine short_char_hermite (self, comp, vec)
    class(rt_patch_t) :: self
    class(composition_type) :: comp
    type(ray_type) :: ray, amr_ray
    real(8), intent(in) :: vec(3)
    real(4) :: Inuup, xxup, ssup, xx, ss, Inu, kx, ky, betup, betc
    real(4) :: dsup1, dsup2, dsup05, rpath, dsc2, u04, Psi4
    real(4) :: v1, v2, v3, gpe, aT, cpe, rdlp, peup
    real(4) :: dchi, t1, sgn, mc, r1mbu, r1mbc, Bu, Bc, dB, gf, kf, dsdtau, dss
    real(8) :: ds(3), dt, rdt, ex, U0, U1, q, A, B, C ! important
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:), allocatable :: row_pe
    real(4), dimension(:,:,:,:), allocatable :: row_red
    real(4), dimension(:), allocatable :: hknu
    real(4), parameter :: dtmin=1e-20, eps_pe=1e-10, eps_B=1e-25
    integer :: xup_t(4,self%nx), yup_t(4,self%ny)
    integer :: zs, ze, xup, yup, zup, i, ii, jj, upo1, upo2
    logical :: amr_boundary, rebuild_tables, get_rates, with_chain
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char_hermite', itimer=itimer)
    if (.not. (allocated(comp%dlnpe_dz) .and. allocated(comp%dlnpe_dx))) &
      stop 'short_char_hermite: directional composition gradients missing -- set_gradients was not called'
    allocate(work_up(self%fs:self%fe,4))
    allocate(work_min(self%fs:self%fe,4))
    allocate(row_red(self%fs:self%fe, 4, self%nx, self%ny))
    allocate(row_pe(self%nx, self%ny))
    allocate(hknu(self%fs:self%fe))
    hknu = real(h_k * self%nu(self%fs:self%fe), 4)
    get_rates = any(self%get_rates)

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

      !--- ROW PASS: x-reduce the whole upstream plane ----------------------------
      row_red = 0.0
      if (with_chain) row_pe = 0.0
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
              do f = self%fs, self%fe
                row_red(f,4,x,yup) = row_red(f,4,x,yup) + kx*self%back(f,xup,yup,zup,3)
              end do
              row_pe(x,yup) = row_pe(x,yup) + kx*comp%lnpe(xup,yup,zup)
            endif
          end do
        end do
      end do

      do y=1, self%ny
        do x=1, self%nx
          work_up = 0.0
          work_min = huge(1.0)
          peup = 0.0

          ! 4x4 kernel sum
          do jj=1,4
            yup = yup_t(jj,y)
            ky = ykern(jj)
            do i=1, 4
              do f=self%fs,self%fe
                work_up(f,i) = work_up(f,i) + ky*row_red(f,i,x,yup)
              end do
            end do
            if (with_chain) peup = peup + ky*row_pe(x,yup)
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
                do f=self%fs,self%fe
                  work_min(f,4) = min(work_min(f,4), self%back(f,xup,yup,zup,3))
                end do
              endif
            end do
          end do
          work_up = max(work_up, work_min)

          if (with_chain) then
            ! per-cell directional gradients at the target
            gpe = v1*comp%dlnpe_dx(x,y,z) + v2*comp%dlnpe_dy(x,y,z) + v3*comp%dlnpe_dz(x,y,z)
            aT  = (v1*comp%dlnt_dx(x,y,z) + v2*comp%dlnt_dy(x,y,z) + v3*comp%dlnt_dz(x,y,z)) &
                  / real(comp%temp(x,y,z), 4)
            rdlp = comp%lnpe(x,y,z) - peup
            cpe = gpe / (rdlp + sign(eps_pe, rdlp))

            do f=self%fs, self%fe
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              xxup  = work_up(f,1)
              ssup  = work_up(f,2)
              Inuup = work_up(f,3)
              betup = min(work_up(f,4), 1.0 - 1e-6)  ! kernel-overshoot guard
              betc  = self%back(f,x,y,z,3)

              dchi = xx - xxup
              t1 = dchi * rpath
              sgn = sign(1.0, t1)
              mc = sgn * max(0.0, min(sgn*(dchi*cpe), 2.0*abs(t1)))
              ! maser interval: keep the sign of dtau, bound it away from zero
              ! and cap the amplification, consistent with the other solvers
              dt = dsup05*(xxup + xx) + dsc2*(t1 - mc)
              dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
              rdt = 1.0d0 / dt
              ex = exp(-min(dt, 80.d0))
              call parab_moments(dt, ex, rdt, U0, U1, q)

              ! the clamped chain slope factors as dsdtau = dss*kf with
              ! kf in [0, 2/dt] (the MC positivity envelope): kf is the exact
              ! dS/dtau sensitivity to the local S, so the Lambda diagonal
              ! B + C*kf stays consistent with the scheme -- Psi = B alone
              ! overestimates it (C < 0) and destabilizes the ALI iteration
              r1mbu = 1.0 / (1.0 - betup)
              r1mbc = 1.0 / (1.0 - betc)
              Bu = betup * r1mbu
              Bc = betc * r1mbc
              dss = ss - ssup
              dB = Bc - Bu
              gf = (hknu(f) * Bc * r1mbc * aT / xx) / (dB + sign(eps_B, dB))
              ! the positivity envelope 2/dt only makes sense for dt>0; on a
              ! maser interval bound kf by 2/|dt| so dsdtau stays finite
              kf = min(max(gf, 0.0), 2.0*abs(real(rdt, 4)))
              dsdtau = dss * kf

              B = (2.0d0*U1 - q) * rdt
              A = U0 - B
              C = q - U1
              Inu = Inuup*real(ex,4) + real(A,4)*ssup + real(B,4)*ss + real(C,4)*dsdtau
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                if (dt > 0 .and. xx > 0) then
                  self%Psi(f,x,y,z) = real(B + C*dble(kf), 4) / xx
                else
                  ! 0 <= Psi <= 1 no longer holds; leave the rates unpreconditioned
                  self%Psi(f,x,y,z) = 0
                endif
              endif
            end do ! f-loop
          else
            ! boundary plane: plain linear step, as in short_char
            do f=self%fs, self%fe
              xx = self%work(f,x,y,z,1)
              ss = self%work(f,x,y,z,2)
              xxup  = work_up(f,1)
              ssup  = work_up(f,2)
              Inuup = work_up(f,3)

              dt = dsup05 * ( xxup + xx )
              ! cap the maser amplification, as in short_char
              dt = max(dt, -maser_dtau_max)
              if (dt < 80) then
                ex = exp(-dt)
                ! the series expansion is only valid for |dt| << 1, and dt<0
                ! must not fall into it
                if ( abs(dt) > 1d-3 ) then
                  u04 = real(1.0d0 - ex, 4)
                  Psi4 = real(1.0d0 - (1.0d0 - ex) / dt, 4)
                else
                  Psi4 = real(dt * (0.5d0 - 0.16666667d0*dt), 4)
                  u04 = real(dt, 4) * (1.0 - Psi4)
                endif
                Inu = Inuup*real(ex,4) + (u04 - Psi4)*ssup + Psi4*ss
              else
                u04 = real(1.0d0 / dt, 4)
                Psi4 = 1.0 - u04
                Inu = u04*ssup + Psi4*ss
              endif
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                if (dt > 0 .and. xx > 0) then
                  self%Psi(f,x,y,z) = Psi4 / xx
                else
                  ! 0 <= Psi <= 1 no longer holds; leave the rates unpreconditioned
                  self%Psi(f,x,y,z) = 0
                endif
              endif
            end do ! f-loop
          endif
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE short_char_hermite

END SUBMODULE rt_solvers_smod
