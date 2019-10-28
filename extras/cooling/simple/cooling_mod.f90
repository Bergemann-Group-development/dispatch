!==============================================================================
!==============================================================================
MODULE cooling_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE scaling_mod
  implicit none
  integer, parameter:: dp=4
  real, parameter:: T_MC=10.0
  private
  type, public:: cooling_t
    real:: dtime
    procedure(cooling_Casiana), pointer:: function => null()
    real, dimension(:,:,:), pointer:: tt => null()
    real, dimension(:,:,:), pointer:: qq => null()
    real, dimension(:,:,:), pointer:: e  => null()
  contains
    procedure:: init
    procedure:: pre_update
    procedure:: courant_condition
  end type
  real:: courant=0.5
  logical:: on=.false.
  integer:: type=0, verbose=0
  logical, save:: first_time=.true.
  type(cooling_t), public:: cooling
CONTAINS

!==============================================================================
!==============================================================================
SUBROUTINE init (self)
 class (cooling_t):: self
 integer:: iostat
 namelist /cooling_params/ on, verbose, type, courant
 !------------------------------------------------------------------------------
 call trace%begin ('simple/cooling_t%init')
 !$omp critical (input_cr)
 if (first_time) then
   first_time = .false.
   rewind (io%input)
   read (io%input, cooling_params, iostat=iostat)
   write (io%output, cooling_params)
 end if
 !$omp end critical (input_cr)
 select case (type)
 case (1)
   self%function => cooling_Casiana
 case (2)
   self%function => cooling_Dalgarno
 case (3)
   self%function => cooling_Rosen_etal
 case (4)
   self%function => cooling_Shure
 case default
   self%function => cooling_Gnedin
 end select
 call trace%end
END SUBROUTINE init

!==============================================================================
!> Implement actual cooling
!==============================================================================
SUBROUTINE pre_update (self, patch)
  class (cooling_t):: self
  class (gpatch_t), pointer:: patch
  integer:: i1, i2, i3, l(3), u(3)
  real:: scaling_constant, q1, q2, eth, tt, qq, d, e, px, py, pz
  real:: T_min, T_max, Q_min, Q_max, E_min, E_max
  !-----------------------------------------------------------------------------
  call trace%begin ('simple/cooling_t%pre_update')
  if (on) then
    if (.not.allocated(patch%heating_per_unit_volume)) then
      u = patch%gn
      allocate (patch%heating_per_unit_volume(u(1),u(2),u(3)))
      patch%heating_per_unit_volume = 0.0
    end if
    self%dtime=huge(1.)
    scaling_constant = scaling%t/scaling%p
    T_min = huge(1.)
    T_max = 0.0
    Q_min = huge(1.)
    Q_max = 0.0
    E_min = huge(1.)
    E_max = 0.0
    l = patch%mesh%lb
    u = patch%mesh%ub
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      !-------------------------------------------------------------------------
      ! Subtract kinetic energy and convert to temperature
      !-------------------------------------------------------------------------
      d  = patch%mem(i1,i2,i3,patch%idx%d ,patch%it,1)
      e  = patch%mem(i1,i2,i3,patch%idx%e ,patch%it,1)
      px = patch%mem(i1,i2,i3,patch%idx%px,patch%it,1)
      py = patch%mem(i1,i2,i3,patch%idx%py,patch%it,1)
      pz = patch%mem(i1,i2,i3,patch%idx%pz,patch%it,1)
      eth = e - 0.5*(px**2+py**2+pz**2)/d
      tt = real(scaling%temp*(patch%gamma-1d0))*eth/d
      !-------------------------------------------------------------------------
      ! Cool above T_MC, heat below
      !-------------------------------------------------------------------------
      q1 = self%function(tt)*d*scaling_constant
      q2 = (T_MC/tt-1.0)*eth/patch%dtime*courant
      qq = merge(q1,q2,tt>T_MC)
      !patch%heating_per_unit_volume(i1,i2,i3) = &
      !patch%heating_per_unit_volume(i1,i2,i3) + qq
      self%dtime = min(self%dtime, courant*eth/qq)
      T_min = min(T_min,tt)
      T_max = max(T_max,tt)
      Q_min = min(Q_min,qq)
      Q_max = max(Q_max,qq)
      E_min = min(E_min,eth)
      E_max = max(E_max,eth)
    end do
    end do
    end do
    if (patch%istep < 2) then
      patch%heating_per_unit_volume = 0.0
    end if
    if (verbose > 0) then
      write(stdout,*) 'simple/cooling_t%pre_update: T_min,max =', T_min, T_max
      write(stdout,*) 'simple/cooling_t%pre_update: Q_min,max =', Q_min, Q_max
    end if
  end if
  call trace%end
END SUBROUTINE pre_update

!==============================================================================
!> Implement actual cooling
!==============================================================================
SUBROUTINE courant_condition (self, patch)
  class(cooling_t):: self
  class(gpatch_t), pointer:: patch
  integer:: i1, i2, i3, l(3), u(3)
  real, dimension(:,:,:), pointer:: e
  !----------------------------------------------------------------------------
  call trace%begin ('simple/cooling_t%courant_condition')
  if (on .and. self%dtime < patch%dtime) then
    if (verbose > 0) then
      write (stdout,'(a,i6,1p,2e12.3)') &
        'simple/cooling_t%courant_condition: WARNING', &
        patch%id, self%dtime, patch%dtime
    end if
    !!!$omp atomic write
    !!patch%dtime = self%dtime
  end if
  call trace%end
END SUBROUTINE courant_condition

!==============================================================================
!==============================================================================
include "cooling_Casiana.f90"
include "cooling_Dalgarno.f90"
include "cooling_Gnedin.f90"
include "cooling_Rosen_etal.f90"
include "cooling_Shure.f90"

END MODULE cooling_mod
