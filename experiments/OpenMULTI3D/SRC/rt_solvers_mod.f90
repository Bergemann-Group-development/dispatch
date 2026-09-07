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
    get_rates = self%compute_psi

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
            dt = max(dt, -maser_dtau_max)     ! cap maser amplification (chi<0 -> dt<0)
            if (dt < 80) then
              ex = exp(-dt)
              if ( abs(dt) > 1d-3 ) then       ! Taylor branch only for |dt| << 1
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
                self%Psi(f,x,y,z) = 0   ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning
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
            dt = max(dt, -maser_dtau_max)     ! cap maser amplification (chi<0 -> dt<0)
            if (dt < 80) then
              Inuup = self%work(f,x,y,zup,3)
              ex = exp(-dt)
              if ( abs(dt) > 1d-3 ) then       ! Taylor branch only for |dt| << 1
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif
            Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))

            self%work(f,x,y,z,3) = Inu
            if (self%compute_psi) then
              if (dt > 0 .and. xx > 0) then
                self%Psi(f,x,y,z) = Psi / xx
              else
                self%Psi(f,x,y,z) = 0   ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning
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
            dt = max(dt, -maser_dtau_max)     ! cap maser amplification (chi<0 -> dt<0)
            if (dt < 80) then
              Inuup = self%work(f,x,y,zup,3)
              ex = exp(-dt)
              if ( abs(dt) > 1d-3 ) then       ! Taylor branch only for |dt| << 1
                u0 = 1 - ex
                Psi = 1 - u0 / dt
              else
                Psi = dt * (0.5 - 0.16666667*dt)
                u0 = dt * (1 - Psi)
              endif

              Inu = Inuup*ex + (u0 - Psi)*ssup + Psi*ss
            else
              u0 = 1 / dt
              Psi = 1 - u0
              Inu = u0*ssup + Psi*ss
            endif
            Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))

            self%work(f,x,y,z,3) = Inu
            if (self%compute_psi) then
              if (dt > 0 .and. xx > 0) then
                self%Psi(f,x,y,z) = Psi / xx
              else
                self%Psi(f,x,y,z) = 0   ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning
              endif
            endif
          end do ! f-loop
        end do
      end do
    end do ! z-loop
    call trace%end(itimer)
  END SUBROUTINE long_char_linear

  !===============================================================================
  !> Parabolic transfer moments over an interval of optical depth dt (ex=exp(-dt),
  !> rdt=1/dt): U0 = 1-exp(-dt), U1 = dt-U0, q = U2/dt with U2 = dt^2-2*U1.  Below
  !> dt=1e-3 they use the Taylor series (U1 = dt^2/2 - dt^3/6, q = dt^2/3 -
  !> dt^3/12) to avoid cancellation in 1-ex.
  !===============================================================================
  pure module subroutine parab_moments (dt, ex, rdt, U0, U1, q)
    real(8), intent(in)  :: dt, ex, rdt
    real(8), intent(out) :: U0, U1, q
    !-----------------------------------------------------------------------------
    if (abs(dt) > 1e-3) then       ! exact moments hold for either sign of dt
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
      real(8) :: dt, rdt, ex, U0, U1, q, A, B, C ! keep in double precision
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
      rdt = 1.0d0 / dt
      ex = exp(-min(dt, 80.d0))
      call parab_moments(dt, ex, rdt, U0, U1, q)

      ! one sided dS/dtau either side of c, and the MC-limited central slope
      if (at_target) then
        t1 = (ssc - ssup) * real(rdt,4)
        t3 = (ssdn - ssc) / dtdn
      else
        t1 = (ssc - ssup) / dtup
        t3 = (ssdn - ssc) * real(rdt,4)
      endif
      t2 = 0.5*(t1 + t3)
      dsdtau = (sign(1.0, t1) + sign(1.0, t3)) * min(abs(t1), abs(t3), abs(0.5*t2))

      ! for dt<0 (bounded maser) or chi<0 the rates are left unpreconditioned
      ! (plam=0), and Inu is capped, as in the other solvers
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
!> exp(-t), from the moments Mn = int_0^dt t^n e^-t dt with a Taylor branch at
!> small dt.
!>
!> The first plane integrates from the boundary ghost plane (chi/Snu/Inu from
!> set_upstream_boundary via halo exchange).  With a real upstream neighbour
!> (comp_up) or a deep BOTTOM boundary (ghost Pe/T extrapolated by ghost_thermo)
!> the full herm_step is used.  Only at the zeroed TOP boundary (chi==0, where the
!> analytic dBnu/dtau=.../chi is singular) does the upwind-end source slope
!> degrade to the interval secant (herm_first).
!===============================================================================
  module subroutine long_char_hermite (self, comp, vec, comp_up)
    class(rt_patch_t), target :: self
    class(composition_type) :: comp
    class(composition_type), pointer :: comp_up
    real(8), intent(in) :: vec(3)
    real(4) :: ds, dsup, rdsup, dsup1, dsup2, gfac
    real(4) :: Inu, ppsi, aTu, aTc, rTup, betup
    real(4) :: peu, pec, pedn, rpeu, rpec, rpedn, gpu, gpc
    real(4), parameter :: dtmin=1e-20
    real(4), dimension(:,:,:,:,:), contiguous, pointer :: w
    real(4), dimension(:,:,:,:), contiguous, pointer :: psi, bk
    real(4), dimension(:), allocatable :: hknu
    integer :: zs, ze, zup, zdn, zff, zu, z1, z2, z3, dir, nx, ny, fs, fe
    logical :: get_rates, use_upstream, extrap_bottom
    real(4) :: xx3g(3)
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%long_char_hermite', itimer=itimer)
    if (.not. allocated(comp%dpe_dz)) &
      stop 'long_char_hermite: composition gradients missing -- set_gradients was not called'
    if (self%nz < 2) &
      stop 'long_char_hermite: needs nz >= 2'
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
    ! Pe/T are available (or extrapolated by ghost_thermo) and the full herm_step
    ! runs; the zeroed top ghost (chi==0) keeps the secant fallback herm_first.
    use_upstream = associated(comp_up)
    z = zs
    zup = z - dir
    zdn = z + dir
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
    do y=1, ny
      do x=1, nx
        if (use_upstream) then
          peu  = comp_up%pe4(x,y,zu)
          rpeu = 1.0 / peu
          gpu  = gfac * comp_up%dpe_dz(x,y,zu)
          rTup = 1.0 / real(comp_up%temp(x,y,zu), 4)   ! upwind-plane beta from T
          aTu  = gfac * comp_up%dlnt_dz(x,y,zu) * rTup
        else if (extrap_bottom) then
          ! ghost Pe/T extrapolated from the three innermost planes, mirroring
          ! the chi/Snu extrapolation in set_upstream_boundary
          call ghost_thermo (x, y, peu, rpeu, gpu, rTup, aTu)
        endif
        if (use_upstream .or. extrap_bottom) then
          pec   = comp%pe4(x,y,z)
          rpec  = 1.0 / pec
          gpc   = gfac * comp%dpe_dz(x,y,z)
          aTc   = gfac * comp%dlnt_dz(x,y,z) / real(comp%temp(x,y,z), 4)
          rpedn = 1.0 / comp%pe4(x,y,zdn)              ! downwind parabola sample
          do f=fs, fe
            ! the ghost carries chi/Snu/Inu but not beta; rebuild beta at the
            ! upwind point from the (neighbour or extrapolated) T, as
            ! lte_background_opacities does (zdn is the third parabola sample)
            betup = max(exp(-hknu(f)*rTup), tiny(1.0))
            call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                            betup, bk(f,x,y,z), bk(f,x,y,zdn), hknu(f), &
                            w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, rpedn, &
                            aTu, aTc, .true., Inu, ppsi)
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
        endif
      end do
    end do

    ! interior planes with a downwind neighbour: the parabola is centred at the
    ! target z (third sample = zdn), so its centred slope lands at the downwind
    ! end of the integrated interval [zup, z] and its one-sided slope at zup
    do z=zs+dir, ze-dir, dir
      zup = z - dir
      zdn = z + dir
      do y=1, ny
        do x=1, nx
          peu   = comp%pe4(x,y,zup)
          pec   = comp%pe4(x,y,z)
          rpeu  = 1.0 / peu
          rpec  = 1.0 / pec
          rpedn = 1.0 / comp%pe4(x,y,zdn)
          gpu   = gfac * comp%dpe_dz(x,y,zup)
          gpc   = gfac * comp%dpe_dz(x,y,z)
          aTu   = gfac * comp%dlnt_dz(x,y,zup) / real(comp%temp(x,y,zup), 4)
          aTc   = gfac * comp%dlnt_dz(x,y,z)   / real(comp%temp(x,y,z),   4)
          do f=fs, fe
            call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zdn,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zdn,2), &
                            bk(f,x,y,zup), bk(f,x,y,z), bk(f,x,y,zdn), hknu(f), &
                            w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, rpedn, &
                            aTu, aTc, .true., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    end do ! z-loop

    ! last plane z==ze: no downwind sample exists, so the parabola is centred at
    ! the upwind end zup (third sample = zff, two planes upwind) and its centred
    ! slope lands at zup, the one-sided slope at the target z -- the same stencil
    ! shift the quadratic solver uses at its final point
    z   = ze
    zup = z - dir
    zff = z - 2*dir
    if (zff >= 1 .and. zff <= self%nz) then
      ! nz >= 3: the two-planes-upwind sample is in the domain
      do y=1, ny
        do x=1, nx
          peu   = comp%pe4(x,y,zup)
          pec   = comp%pe4(x,y,z)
          rpeu  = 1.0 / peu
          rpec  = 1.0 / pec
          rpedn = 1.0 / comp%pe4(x,y,zff)
          gpu   = gfac * comp%dpe_dz(x,y,zup)
          gpc   = gfac * comp%dpe_dz(x,y,z)
          aTu   = gfac * comp%dlnt_dz(x,y,zup) / real(comp%temp(x,y,zup), 4)
          aTc   = gfac * comp%dlnt_dz(x,y,z)   / real(comp%temp(x,y,z),   4)
          do f=fs, fe
            call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zff,1), &
                            w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zff,2), &
                            bk(f,x,y,zup), bk(f,x,y,z), bk(f,x,y,zff), hknu(f), &
                            w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, rpedn, &
                            aTu, aTc, .false., Inu, ppsi)
            w(f,x,y,z,3) = Inu
            if (get_rates) psi(f,x,y,z) = ppsi
          end do ! f-loop
        end do
      end do
    else
      ! nz == 2: the two-planes-upwind sample is the boundary ghost (zff == zup of
      ! the first plane).  Its Pe/beta come from the neighbour, or (zeroed top
      ! ghost) only need to be finite since chi/Pe and Snu/Bnu vanish there
      do y=1, ny
        do x=1, nx
          peu   = comp%pe4(x,y,zup)
          pec   = comp%pe4(x,y,z)
          rpeu  = 1.0 / peu
          rpec  = 1.0 / pec
          gpu   = gfac * comp%dpe_dz(x,y,zup)
          gpc   = gfac * comp%dpe_dz(x,y,z)
          aTu   = gfac * comp%dlnt_dz(x,y,zup) / real(comp%temp(x,y,zup), 4)
          aTc   = gfac * comp%dlnt_dz(x,y,z)   / real(comp%temp(x,y,z),   4)
          if (use_upstream) then
            rTup  = 1.0 / real(comp_up%temp(x,y,zu), 4)   ! ghost beta from neighbour T
            rpedn = 1.0 / comp_up%pe4(x,y,zu)             ! ghost Pe from neighbour
            do f=fs, fe
              betup = max(exp(-hknu(f)*rTup), tiny(1.0))
              call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zff,1), &
                              w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zff,2), &
                              bk(f,x,y,zup), bk(f,x,y,z), betup, hknu(f), &
                              w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, rpedn, &
                              aTu, aTc, .false., Inu, ppsi)
              w(f,x,y,z,3) = Inu
              if (get_rates) psi(f,x,y,z) = ppsi
            end do ! f-loop
          else
            do f=fs, fe    ! zeroed top ghost: beta/Pe from the target plane (finite)
              call herm_step (w(f,x,y,zup,1), w(f,x,y,z,1), w(f,x,y,zff,1), &
                              w(f,x,y,zup,2), w(f,x,y,z,2), w(f,x,y,zff,2), &
                              bk(f,x,y,zup), bk(f,x,y,z), bk(f,x,y,z), hknu(f), &
                              w(f,x,y,zup,3), peu, rpeu, gpu, pec, rpec, gpc, rpec, &
                              aTu, aTc, .false., Inu, ppsi)
              w(f,x,y,z,3) = Inu
              if (get_rates) psi(f,x,y,z) = ppsi
            end do ! f-loop
          endif
        end do
      end do
    endif
    call trace%end(itimer)
  contains
    !---------------------------------------------------------------------------
    !> Ghost-plane thermodynamics at a deep (bottom) domain boundary, where there
    !> is no upstream neighbour.  set_upstream_boundary fills the ghost chi/Snu by
    !> monotone extrapolation of the three innermost planes; this does the same
    !> for Pe and T (and their d/dz gradients), returning the exact quantities the
    !> upstream branch reads from a real neighbour, so herm_step is used unchanged.
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
      ! exact moments hold for either sign of dt; Taylor branch only for |dt| << 1
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
        ! Taylor series below dt=0.05, where the moment combinations cancel
        w0 = dt*(0.25d0 - dt*(0.20d0 - dt*(1.d0/12.d0  - dt/42.d0)))
        w1 = dt*(0.25d0 - dt*(0.15d0 - dt*(1.d0/20.d0  - dt/84.d0)))
        w2 = dt*(0.25d0 - dt*(0.10d0 - dt*(1.d0/40.d0  - dt/210.d0)))
        w3 = dt*(0.25d0 - dt*(0.05d0 - dt*(1.d0/120.d0 - dt/840.d0)))
      endif
    end subroutine bezier3_weights
    !---------------------------------------------------------------------------
    !> One integrated interval [up, c].  chi/Pe and Snu/Bnu are each fitted by a
    !> PARABOLA through the upwind, central and downwind samples (xx3/ss3/... is
    !> the third sample); their derivatives set the product-rule chi slopes
    !> (-> dtau) and source slopes (-> Bezier control points).  center_target
    !> picks which end of [up,c] gets the centred 3-point estimate:
    !>   .true.  -> parabola centred at the target c (interior sweep; third sample
    !>              downwind), one-sided slope at up.
    !>   .false. -> parabola centred at up (final plane z==ze; third sample two
    !>              planes upwind), one-sided slope at c -- the quadratic end shift.
    !---------------------------------------------------------------------------
    pure subroutine herm_step (xxup, xxc, xx3, ssup, ssc, ss3, betau, betac, beta3, &
                               hk, Inuup, peu, rpeu, gpu, pec, rpec, gpc, rpe3, &
                               aTu, aTc, center_target, Inu, plam)
      real(4), intent(in) :: xxup, xxc, xx3, ssup, ssc, ss3, betau, betac, beta3, hk, Inuup
      real(4), intent(in) :: peu, rpeu, gpu, pec, rpec, gpc, rpe3, aTu, aTc
      logical, intent(in) :: center_target
      real(4), intent(out) :: Inu, plam
      real(4) :: dchi, t1, sgn, mup, mc, r1mbu, r1mbc, r1mb3, Bu, Bc, B3, gu, gc
      real(4) :: rxpu, rxpc, rxp3, drxpu, drxpc, hh, rxxu, rxxc
      real(4) :: dt3, smin, smax, P1, P2, rhou, rhoc, rho3, drhu, drhc, depu, depc, d1, d2
      real(8) :: dt, ex, w0, w1, w2, w3
      !-------------------------------------------------------------------------
      hh   = 0.5 * rdsup                      ! 1/(2*ds), the parabola half-step
      dchi = xxc - xxup
      t1   = dchi * rdsup                     ! chi secant over [up,c] (limiter ref)
      sgn  = sign(1.0, t1)
      rxpu = xxup * rpeu                      ! chi/Pe at the upwind point
      rxpc = xxc  * rpec                      ! chi/Pe at the target
      rxp3 = xx3  * rpe3                      ! chi/Pe at the third parabola sample

      ! source departure ratios rho = Snu/Bnu at the three samples
      r1mbu = 1.0 / (1.0 - betau)
      r1mbc = 1.0 / (1.0 - betac)
      r1mb3 = 1.0 / (1.0 - beta3)
      Bu = betau * r1mbu
      Bc = betac * r1mbc
      B3 = beta3 * r1mb3
      rhou = ssup / Bu
      rhoc = ssc  / Bc
      rho3 = ss3  / B3

      ! 3-point parabola derivatives d(chi/Pe)/ds and d(rho)/ds at both ends of
      ! [up,c]: the centred estimate at the parabola centre, the one-sided
      ! estimate (same parabola) at the far end
      if (center_target) then                ! centre = c (third sample downwind)
        drxpc = (rxp3 - rxpu) * hh                     ! centred at c
        drxpu = (4.0*rxpc - 3.0*rxpu - rxp3) * hh       ! left one-sided at up
        drhc  = (rho3 - rhou) * hh
        drhu  = (4.0*rhoc - 3.0*rhou - rho3) * hh
      else                                   ! centre = up (third sample 2 upwind)
        drxpu = (rxpc - rxp3) * hh                      ! centred at up
        drxpc = (rxp3 - 4.0*rxpu + 3.0*rxpc) * hh        ! right one-sided at c
        drhu  = (rhoc - rho3) * hh
        drhc  = (rho3 - 4.0*rhou + 3.0*rhoc) * hh
      endif

      ! dchi/ds at both ends by the product rule chi = (chi/Pe)*Pe, then clamped
      ! to the sign of the interval secant and twice its magnitude (MC envelope)
      mup = rxpu*gpu + peu*drxpu
      mc  = rxpc*gpc + pec*drxpc
      mup = sgn * max(0.0, min(sgn*mup, 2.0*abs(t1)))
      mc  = sgn * max(0.0, min(sgn*mc,  2.0*abs(t1)))
      ! maser interval (chi<0 -> dtau<0): keep the sign, hold |dtau| away from
      ! zero and cap the amplification
      dt = dsup1*(xxup + xxc) + dsup2*(mup - mc)
      dt = max(sign(max(abs(dt), dble(dtmin)), dt), -maser_dtau_max)
      ex = exp(-min(dt, 80.d0))
      call bezier3_weights (dt, ex, w0, w1, w2, w3)

      ! dSnu/dtau at both ends by the product rule S = rho*Bnu (rho = S/Bnu):
      !   dS/dtau = rho*(dBnu/dtau) + Bnu*(drho/dtau)
      ! the analytic dBnu/dtau = hk*B*r1mb*aT/chi always enters as rho*dBnu/dtau
      ! (= S*(dBnu/dtau)/B) or as (dBnu/dtau)/B in the Lambda term, so its B factor
      ! cancels rho's 1/B: gu/gc hold (dBnu/dtau)/Bnu = hk*r1mb*aT/chi.  The
      ! departure part is the parabola's d(rho)/ds divided by chi.
      rxxu = 1.0 / xxup                      ! 1/chi at up and c, reused below
      rxxc = 1.0 / xxc
      gu = hk * r1mbu * aTu * rxxu           ! (dBnu/dtau)/Bnu at up
      gc = hk * r1mbc * aTc * rxxc           ! (dBnu/dtau)/Bnu at c
      depu = drhu * rxxu                     ! (drho/dtau)|up = (drho/ds)/chi
      depc = drhc * rxxc                     ! (drho/dtau)|c

      dt3 = real(dt,4) * one_third
      P1 = ssup + dt3*(ssup*gu + Bu*depu)    ! Sup + (dtau/3)*(dS/dtau)|up
      P2 = ssc  - dt3*(ssc*gc  + Bc*depc)    ! Sc  - (dtau/3)*(dS/dtau)|c
      smin = min(ssup, ssc)                  ! monotone clamp: overshoot-free
      smax = max(ssup, ssc)
      P1 = min(max(P1, smin), smax)
      P2 = min(max(P2, smin), smax)

      Inu = Inuup*real(ex,4) + real(w0,4)*ssup + real(w1,4)*P1 &
                             + real(w2,4)*P2   + real(w3,4)*ssc
      Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
      plam = 0.0
      if (get_rates .and. dt > 0 .and. xxc > 0) then
        ! Lambda diagonal from dP1/dSc = Bu/(3*Bc), dP2/dSc = 2/3 - dt3*gc (cross-
        ! terms and clamp neglected), both bounded into [0,1] so it cannot exceed
        ! the attenuated total or go negative.  For dt<0/chi<0 leave plam=0.
        d1 = min(max(Bu*real(one_third,4)/Bc, 0.0), 1.0)
        d2 = min(max(2.0*real(one_third,4) - dt3*gc, 0.0), 1.0)
        plam = real(w1*dble(d1) + w2*dble(d2) + w3, 4) * rxxc
      endif
    end subroutine herm_step
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
      real(8) :: dt, ex, w0, w1, w2, w3
      !-------------------------------------------------------------------------
      ! target-end chi slope by the same chi = (chi/Pe)*Pe product rule as
      ! herm_step; the ghost plane has no Pe, so the one-sided d(chi/Pe)/ds
      ! comes from the adjacent DOWN interval
      t1 = (xxc - xxup) * rdsup
      sgn = sign(1.0, t1)
      rxpc = xxc*rpec                          ! chi/Pe at c (reused in the secant)
      mc = rxpc*gpc + pec*(xxdn*rpedn - rxpc)*rdsup
      mc = sgn * max(0.0, min(sgn*mc, 2.0*abs(t1)))
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
      ! target-end slope by the same product rule S = rho*Bnu as herm_step, the
      ! Bnu again cancelling rho's 1/Bnu (Gc = (dBnu/dtau)/Bnu); the departure
      ! ratio secant across the DOWN interval uses the leading-trapezoid dtau
      Gc = hk * r1mbc * aTc / xxc              ! (dBnu/dtau)/Bnu at c
      rhoc  = ssc  / Bc
      rhodn = ssdn / Bdn
      taudn = dsup1 * (xxc + xxdn)             ! leading trapezoid dtau, down side
      taudn = sign(max(abs(taudn), real(dtmin,4)), taudn)
      dSc = ssc*Gc + Bc*(rhodn - rhoc)/taudn
      dss = ssc - ssup
      dSu = dss * real(1.d0/dt, 4)      ! interval secant slope wrt tau

      dt3 = real(dt,4) * one_third
      P1 = ssup + dt3*dSu
      P2 = ssc  - dt3*dSc
      smin = min(ssup, ssc)                    ! monotone clamp: overshoot-free
      smax = max(ssup, ssc)
      P1 = min(max(P1, smin), smax)
      P2 = min(max(P2, smin), smax)

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
!> Hybrid formal solver ('hermite2'), long-characteristics.  Optical depth is
!> built exactly as long_char_hermite's chain-rule cubic-Hermite chi treatment:
!>
!>   dchi/ds = (chi/Pe) * (dPe/ds)  +  Pe * d(chi/Pe)/ds
!>
!> with the product-rule slopes from a parabola through the upwind/central/
!> downwind samples, clamped to the cubic-Bezier control-point bound
!> (|dchi/ds| <= 3*|secant|, the monotonicity limit matching the cubic-Hermite
!> chi that the Euler-Maclaurin trapezoid dtau integrates exactly).  The same
!> Euler-Maclaurin trapezoid dtau as herm_step.  The source term, however, is
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
      stop 'long_char_hermite2: composition gradients missing -- set_gradients was not called'
    if (self%nz < 2) &
      stop 'long_char_hermite2: needs nz >= 2'
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
    ! herm_step-identical) chi/Pe slope estimate at both ends even when
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
    !> One integrated interval [up, c].  chi/Pe slopes as in herm_step: a
    !> parabola through (up, c, 3) gives the product-rule dchi/ds at up
    !> (one-sided) and at c (centred when center_target, one-sided otherwise),
    !> here clamped to the cubic-Bezier bound |dchi/ds| <= 3*|secant| (herm_step
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
      real(8) :: dt, dt3int, rdt, ex, U0, U1, q, A, B, C
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
      ex = exp(-min(dt, 80.d0))
      rdt = 1.0d0 / dt
      call parab_moments(dt, ex, rdt, U0, U1, q)

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
      Inu = Inuup*real(ex,4) + real(A,4)*ssup + real(B,4)*ssc + real(C,4)*dsdtau
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
      real(8) :: dt, rdt, ex, U0, U1, q, A, B, C
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
      ex = exp(-min(dt, 80.d0))
      rdt = 1.0d0 / dt
      call parab_moments(dt, ex, rdt, U0, U1, q)

      dss = ssc - ssup
      dsdtau = dss * real(rdt, 4)      ! no data beyond the target: interval secant

      B = (2.0d0*U1 - q) * rdt
      A = U0 - B
      C = q - U1
      Inu = Inuup*real(ex,4) + real(A,4)*ssup + real(B,4)*ssc + real(C,4)*dsdtau
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
    real(4) :: dsup1, dsup2, dsup05, rpath, dsc2, u04, Psi4
    real(4) :: v1, v2, v3, gpe, aT, pec, rpec, rxpc, rxpup
    real(4) :: dchi, t1, sgn, mc, r1mbc, Bc, kf, dsdtau, dss
    real(4) :: rhou, rhoc, gc, t1s, sgns, rdt4
    real(8) :: ds(3), dt, rdt, ex, U0, U1, q, A, B, C ! keep in double precision
    real(4), dimension(4), target :: xkernel1, ykernel1, xkernel2, ykernel2
    real(4), dimension(:), pointer, contiguous :: xkern, ykern
    real(4), dimension(:,:), allocatable :: work_up, work_min
    real(4), dimension(:,:,:,:), allocatable :: plane_r
    real(4), dimension(:,:,:,:), allocatable :: row_red
    real(4), dimension(:), allocatable :: hknu
    real(4), parameter :: dtmin=1e-20
    integer :: xup_t(4,self%nx), yup_t(4,self%ny)
    integer :: zs, ze, xup, yup, zup, i, ii, jj, upo1, upo2
    logical :: amr_boundary, rebuild_tables, get_rates, with_chain
    integer, save :: itimer=0
    integer :: x, y, z, f
    !-----------------------------------------------------------------------------
    call trace%begin('rt_patch%short_char_hermite', itimer=itimer)
    if (.not. (allocated(comp%dpe_dz) .and. allocated(comp%dpe_dx))) &
      stop 'short_char_hermite: directional composition gradients missing -- set_gradients was not called'
    allocate(work_up(self%fs:self%fe,5))
    allocate(work_min(self%fs:self%fe,5))
    allocate(row_red(self%fs:self%fe, 5, self%nx, self%ny))
    allocate(plane_r(self%fs:self%fe, self%nx, self%ny, 4:5))
    allocate(hknu(self%fs:self%fe))
    hknu = real(h_k * self%nu(self%fs:self%fe), 4)
    get_rates = self%compute_psi

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
              ex = exp(-min(dt, 80.d0))
              call parab_moments(dt, ex, rdt, U0, U1, q)

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
              dsdtau = sgns * max(0.0, min(sgns*dsdtau, 2.0*abs(t1s)))

              B = (2.0d0*U1 - q) * rdt
              A = U0 - B
              C = q - U1
              Inu = Inuup*real(ex,4) + real(A,4)*ssup + real(B,4)*ss + real(C,4)*dsdtau
              Inu = max(-real(maser_Inu_max), min(real(maser_Inu_max), Inu))
              self%work(f,x,y,z,3) = Inu
              if (get_rates) then
                if (dt > 0 .and. xx > 0) then
                  ! Lambda diagonal from the product-rule slope (clamp neglected):
                  ! d(dsdtau)/dss = gc + rdt, bounded into [0, 2/dt] so B + C*kf
                  ! stays within the attenuated total
                  kf = min(max(gc + rdt4, 0.0), 2.0*abs(rdt4))
                  self%Psi(f,x,y,z) = real(B + C*dble(kf), 4) / xx
                else
                  self%Psi(f,x,y,z) = 0   ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning
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
              dt = max(dt, -maser_dtau_max)     ! cap maser amplification
              if (dt < 80) then
                ex = exp(-dt)
                if ( abs(dt) > 1d-3 ) then       ! Taylor branch only for |dt| << 1
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
                  self%Psi(f,x,y,z) = 0   ! 0<=Psi<=1 fails for dt<0/chi<0; no preconditioning
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
