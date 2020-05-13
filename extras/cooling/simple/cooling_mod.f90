!===============================================================================
!> Cooling / heating prescriptions, analogous to the treatment in RAMSES
!> (patch/mhd/troels/cooling.f90).
!===============================================================================
MODULE cooling_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE scaling_mod
  USE units_mod
  USE kinds_mod
  implicit none
  private
  integer, parameter:: dp=KindScalarVar
  type, public:: cooling_t
    real(dp):: dtime
    procedure(cooling_Gnedin ), nopass, pointer:: function => null()
    real(dp), dimension(:,:,:), pointer:: tt => null()
    real(dp), dimension(:,:,:), pointer:: qq => null()
    real(dp), dimension(:,:,:), pointer:: e  => null()
    class(gpatch_t), pointer:: patch => null()
    contains
   procedure:: pre_init
   procedure:: pre_update
    procedure:: courant_condition
    procedure, private:: subcycle
    procedure, private:: test
  end type
  real(dp):: G0 = 3e-25                            ! heating prop to nH1^2
  real(dp):: T_G0 = 2e4                            ! G0 term above T_G0
  real(dp):: G1 = 1.e-25                           ! heating prop to nH1^1
  real(dp):: G2 = 0.0                              ! heating prop to nH1^0 (disabled)
  real(dp):: G_MC = 1.e-24                         ! heating below T_MC
  real(dp):: T_MC = 10.                            ! stable low (MC) temperature
  real(dp):: T_lim = 1e8                           ! limiting temperature
  real(dp):: N_crit = 500.                         ! cut-off number density
  real(dp):: max_fraction = 0.2                    ! max fractional cooling per iteration step
  real(dp):: end_fraction = 0.01                   ! max fractional cooling per iteration step
  real(dp):: mu = 1.2                              ! Molecular weight used for translating from T/mu to temperature
  real(dp):: courant = 0.5                         ! Cooling Courant number 
  integer:: niter_max = 1000                   ! Iteration limit
  integer:: verbose=0
  logical:: on=.false., do_heat=.true., do_subcycle=.true.
  character(len=16):: type='Gnedin'
CONTAINS

!===============================================================================
!> Initialize parameters for cooling / heating
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(cooling_t):: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /cooling_params/ verbose, type, courant, &
    on, do_heat, do_subcycle, G0, G1, G2, G_MC, T_G0, T_MC, N_crit, &
    niter_max, max_fraction, end_fraction, mu
  !-----------------------------------------------------------------------------
  call trace%begin ('simple/cooling_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, cooling_params, iostat=iostat)
    call io%namelist_warning ('cooling_params', iostat)
    write (io%output, cooling_params)
    write (io%namelist, cooling_params)
    if (verbose > 0) then
      self%function => cooling_Gnedin
      call self%test
    end if
    first_time = .false.
  end if
  call io%lock%unset
  select case (trim(type))
  case ('Casiana')
    self%function => cooling_Casiana
  case ('Dalgarno')
    self%function => cooling_Dalgarno
  case ('Rosen')
    self%function => cooling_Rosen_etal
  case ('Shure')
    self%function => cooling_Shure
  case ('Gnedin')
    self%function => cooling_Gnedin
  case default
    self%function => cooling_Gnedin
  end select
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Implement actual cooling
!===============================================================================
SUBROUTINE pre_update (self)
  class (cooling_t):: self
  !.............................................................................
  class (gpatch_t), pointer:: patch
  integer:: i1, i2, i3, l(3), u(3)
  real(dp):: q1, q2, eth, eth0, T, qq, d, e, px, py, pz, v, ekin, eth_MC
  real(dp):: crate, hrate, eth_lim
  real(dp):: T_min, T_max, Q_min, Q_max, r_min, r_max
  integer, save:: itimer=0
  logical, save:: info=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('cooling_t%pre_update', itimer=itimer)
  if (on) then
    patch => self%patch
    if (.not.associated(patch%heating_per_unit_volume)) then
      u = patch%gn
      allocate (patch%heating_per_unit_volume(u(1),u(2),u(3)))
      patch%heating_per_unit_volume = 0.0
      if (info) then
        !$omp critical (info_cr)
        write (io_unit%output,*) &
          'cooling_t%pre_update INFO: allocated heating_per_unit_volume'
        info = .false.
        !$omp end critical (info_cr)
      end if
    end if
    if (do_subcycle) then
      self%dtime = patch%dtime
    else
      self%dtime = huge(1.0)
    end if
    eth_MC = T_MC/(scaling%temp*mu)
    eth_lim = T_lim/(scaling%temp*mu)
    T_min = huge(1.)
    T_max = 0.0
    Q_min = huge(1.)
    Q_max = 0.0
    r_min = huge(1.)
    r_max = 0.0
    l = patch%mesh%lb
    u = patch%mesh%ub
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      !-------------------------------------------------------------------------
      ! Subtract kinetic energy and convert to temperature
      !-------------------------------------------------------------------------
      d  = patch%mem(i1,i2,i3,patch%idx%d ,patch%new,1)
      e  = patch%mem(i1,i2,i3,patch%idx%e ,patch%new,1)
      px = patch%mem(i1,i2,i3,patch%idx%px,patch%new,1)
      py = patch%mem(i1,i2,i3,patch%idx%py,patch%new,1)
      pz = patch%mem(i1,i2,i3,patch%idx%pz,patch%new,1)
      v  = 1./d
      ekin = 0.5*(px**2+py**2+pz**2)*v
      eth = (e-ekin)*v
      eth = max(eth,eth_MC)
      hrate = heat (d, eth)
      crate = cool (self, d, eth)
      qq = hrate-crate
      !-------------------------------------------------------------------------
      ! If sub-cycling, define the heating rate by the change, otherwise set
      ! a max timestep in self%dtime
      !-------------------------------------------------------------------------
      if (do_subcycle) then
        if (abs(patch%dtime*qq)/eth > max_fraction) then
          eth0 = eth
          call self%subcycle (d, eth, patch%dtime)
          eth = min(eth,eth_lim)
          qq = (eth-eth0)/patch%dtime
        else
          eth0 = eth
          eth = eth + patch%dtime*qq
          eth = min(eth,eth_lim)
          qq = (eth-eth0)/patch%dtime          
        end if
      else if (patch%dt_fixed < 0d0) then
          self%dtime = min(self%dtime, courant*eth/max(crate,hrate))
      end if
      patch%heating_per_unit_volume(i1,i2,i3) = d*qq
      T = eth*scaling%temp*mu
      T_min = min(T_min,T)
      T_max = max(T_max,T)
      Q_min = min(Q_min,qq)
      Q_max = max(Q_max,qq)
      r_min = min(r_min,qq/eth)
      r_max = max(r_max,qq/eth)
    end do
    end do
    end do
    if (patch%istep < 2) then
      patch%heating_per_unit_volume = 0.0
    end if
    if (verbose > 0) then
      write (io_unit%log,'(1x,a,i6,f12.6,1p,3(2x,2e10.2))') &
        'cooling_t%pre_update: id, {T,Q,rate}_min,max =', &
        patch%id, patch%time, T_min, T_max, Q_min, Q_max, r_min, r_max
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> Limit on timestep from cooling
!===============================================================================
SUBROUTINE courant_condition (self)
  class(cooling_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  integer:: i1, i2, i3, l(3), u(3), unit
  real(dp), dimension(:,:,:), pointer:: e
  !----------------------------------------------------------------------------
  call trace%begin ('cooling_t%courant_condition')
  if (on) then
    patch => self%patch
    if (self%dtime < patch%dtime) then
      unit = merge (stdout, io_unit%log, verbose > 1)
      if (verbose > 0) then
        write (unit,'(a,i6,1p,2e12.3)') &
          'cooling_t%courant_condition: WARNING', &
          patch%id, self%dtime, patch%dtime
      end if
      !$omp atomic write
      patch%dtime = self%dtime
      patch%u_max = maxval (patch%ds, mask=patch%n > 1)/patch%dtime
    end if
  end if
  call trace%end ()
END SUBROUTINE courant_condition

!===============================================================================
! Cooling function, returning cooling per unit mass
!===============================================================================
FUNCTION cool (self, d, ee) RESULT (crate)
  class(cooling_t):: self
  real(kind=dp), intent(in):: d, ee
  real(kind=dp):: crate, nH1, temp
  !-----------------------------------------------------------------------------
  if (on) then
    nH1 = d*scaling%d/(mu*cgs%amu)
    temp = ee*scaling%temp*mu
    crate = nH1**2*self%function(temp)
    crate = crate/(d*scaling%d)         ! cooling per unit mass in CGS
    crate = crate*scaling%t/scaling%ee  ! cooling per unit mass in code units
  else
    crate = 0.0
  end if
END FUNCTION cool

!===============================================================================
!> Heating rate per unit mass by photoionization (by G0*n) and by
!> G1*n^0 whose origin could be cosmic rays or photoelectric heating
!> (Reynolds & Cox 1992).  The shielding term "onoff" is inspired by
!> Franco & Cox (1996).  To be consistent with MacLow et al the G1 term
!> is turned off over 2e4 K.
!===============================================================================
PURE FUNCTION heat (d, ee) RESULT (hrate)
  real(kind=dp), intent(in):: d, ee
  real(kind=dp):: hrate
  real(kind=dp):: temp, hrate2, nH1, onoff, qq, maxrate, dt_sub, t_sub
  !-----------------------------------------------------------------------------
  if (do_heat) then
    temp = ee*scaling%temp*mu
    nH1 = d*scaling%d/(mu*cgs%amu)
    hrate = merge(G0*nH1**2, G1*nH1, temp>T_G0)
    hrate2 = G_MC*nH1*max(1.-temp/T_MC,0.0_dp) 
    onoff = exp(-nH1/N_crit)
    hrate = onoff*hrate + hrate2
    hrate = hrate/(d*scaling%d)         ! heating per unit mass in CGS
    hrate = hrate*scaling%t/scaling%ee  ! heating per unit mass in code units
  else
    hrate = 0.0
  end if
END FUNCTION heat

!==============================================================================
!> Sub-cycling the cooling / heating balance over a timestep dt_cool
!==============================================================================
SUBROUTINE subcycle (self, d, ee, dt_cool)
  class(cooling_t):: self
  class(gpatch_t), pointer:: patch
  real(dp):: d, ee, n
  real(8):: dt_cool, t_sub, dt_sub
  real(dp):: ee_sub, dee, crate, hrate, hrate2, maxrate, nH1, onoff, qmaxrate, T0, T
  integer:: niter
  !-----------------------------------------------------------------------------
  niter = 0
  t_sub = 0d0
  T = ee*scaling%temp*mu
  T0 = T
  do while (niter < niter_max .and. t_sub < dt_cool)
    !---------------------------------------------------------------------------
    ! Heating and cooling rates are per unit mass
    !---------------------------------------------------------------------------
    hrate = heat (d, ee)
    crate = cool (self, d, ee)
    !---------------------------------------------------------------------------
    ! Limit fractional change and avoid overshoot
    !---------------------------------------------------------------------------
    maxrate = max(crate,hrate)
    dt_sub  = max_fraction*ee/maxrate
    dt_sub  = min(dt_sub,dt_cool-t_sub)
    dee     = dt_sub*(hrate-crate)
    ee      = ee + dee
    t_sub   = t_sub + dt_sub
    T       = ee*scaling%temp*mu
    niter   = niter + 1
    !---------------------------------------------------------------------------
    ! Exit if fractional change is small or temperature goes below T_MC
    !---------------------------------------------------------------------------
    if (abs(dee)/ee < end_fraction) then
      t_sub = dt_cool
      exit
    end if
    if (T < T_MC) then
      T = T_MC
      ee = T_MC/(scaling%temp*mu)
      exit
    end if
    if (verbose > 2) &
      write (io_unit%log,'(1x,a,i5,1p,5e12.4)') &
        'cooling_t%subcycle: niter, dt_sub, t_sub, hrate, crate, T =', &
        niter, dt_sub, t_sub, hrate, crate, T
  end do
  if (verbose > 1) then
    n = d*scaling%d/(cgs%amu*mu)
    write (io_unit%log,'(1x,a,i5,1p,4e12.4)') &
      'cooling_t%subcycle: niter, n, dt_cool, T0, T =', niter, n, dt_cool, T0, T
  end if
END SUBROUTINE subcycle

!===============================================================================
!> Print cooling function values for two different densities
!===============================================================================
SUBROUTINE test (self)
  class(cooling_t):: self
  real(dp):: T, crates(5), n, d, eth, hrate
  real(8):: dt
  integer:: i
  logical:: save(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('simple/cooling_t%test')
  save = [on, do_heat, do_subcycle]
  on = .true.; do_heat = .true.; do_subcycle = .true.
  n = 1e3
  write (io_unit%log,'(1x,a,1p,3e11.2)') 'simple/cooling_t%test: n, scale =', &
    n, scaling%t, scaling%ee
  do i=2,16
    T = 10.**(0.5*i)
    crates(1) = cooling_Casiana (T)
    crates(2) = cooling_Dalgarno (T)
    crates(3) = cooling_Rosen_etal (T)
    crates(4) = cooling_Shure (T)
    crates(5) = cooling_Gnedin (T)
    d = n*cgs%amu/scaling%d
    eth = T/(scaling%temp*mu)
    hrate = heat (d, eth)
    write (io_unit%log,'(1x,a,1p,2(e11.2,2x),5e11.2)') &
      'cooling_t%test: T, hrate, crates =', T, hrate, crates
  end do
  if (verbose > 1) then
    dt = cgs%myr/scaling%t
    do_subcycle = .true.
    !--- case 1 ----------------------------------------------------------------
    eth = 3e4/(scaling%temp*mu)
    d = 1e3*cgs%amu*mu/scaling%d
    call self%subcycle (d, eth, dt)
    !--- case 2 ----------------------------------------------------------------
    eth = 3e4/(scaling%temp*mu)
    d = 1e1*cgs%amu*mu/scaling%d
    call self%subcycle (d, eth, dt)
    !--- case 3 ----------------------------------------------------------------
    eth = 3e7/(scaling%temp*mu)
    d = 1e1*cgs%amu*mu/scaling%d
    call self%subcycle (d, eth, dt)
  end if
  on = save(1); do_heat = save(2); do_subcycle = save(3)
  call trace%end ()
END SUBROUTINE test

!===============================================================================
!> Include cooling function definitions
!===============================================================================
include "cooling_Casiana.f90"
include "cooling_Dalgarno.f90"
include "cooling_Gnedin.f90"
include "cooling_Rosen_etal.f90"
include "cooling_Shure.f90"

END MODULE cooling_mod
