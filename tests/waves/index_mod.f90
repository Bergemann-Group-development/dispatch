!> =============================================================================
!> This index file has slot indices for all solver, all initially equal to zero
!> It is the responsibility of the solver to set all values it needs to use.
!> NOTE: The only reference to indices in tasks should be via task%idx, never to
!> a temporary instance of index_t.
!> =============================================================================
MODULE index_mod
  USE io_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=0, e=0, px=0, py=0, pz=0, bx=0, by=0, bz=0, phi=0, tt=0, s=0, et=0, dphi=0
  contains
    procedure:: init
    procedure:: write
  end type
  type(index_t), public:: idx
CONTAINS

!> =============================================================================
!> Initialize the indices to default values, possibly to be modified by solvers
!> =============================================================================
SUBROUTINE init (self, ie, mhd)
  class(index_t):: self
  integer:: ie
  logical:: mhd
  !-----------------------------------------------------------------------------
  if (ie==2) then
    self%d  = 1
    self%e  = 2
    self%s  = 2
    self%px = 3
    self%py = 4
    self%pz = 5
  else
    self%d  = 1
    self%px = 2
    self%py = 3
    self%pz = 4
    self%e  = 5
    self%s  = 5
  end if
  if (mhd) then
    self%bx = 6
    self%by = 7
    self%bz = 8
  else
    self%bx = -1
    self%by = -1
    self%bz = -1
  end if
  !-----------------------------------------------------------------------------
  ! Ugly fix, necessary because there are a large number of places that refer to
  ! a local instance of index_t, which is assumed to be authentative.  In cases
  ! where the default values have been modified by the server, this may lead to
  ! inconsistencies.  As a fall back, such places may instead refer to this one
  ! pubclic instance.
  !-----------------------------------------------------------------------------
  !$omp critical (index_cr)
  idx = self
  print *, 'IDX', idx%e
  !$omp end critical (index_cr)
END SUBROUTINE init

!> =============================================================================
!> Write out index value namelist, for python/dispatch/
!> =============================================================================
SUBROUTINE write (self)
  class(index_t):: self
  integer:: d, e, et, s, px, py, pz, bx, by, bz, phi, dphi, p1, p2, p3, b1, b2, b3
  integer:: dddt, dedt, dsdt, dpxdt, dpydt, dpzdt, dbxdt, dbydt, dbzdt, dphidt
  namelist /idx_nml/ d, e, s, px, py, pz, bx, by, bz, phi, dphi, p1, p2, p3, b1, b2, b3
  if (io%time_derivs) then
    call output_time_derivs (self)
  else
    d    = self%d   - 1
    e    = self%e   - 1
    s    = self%s   - 1
    et   = self%et  - 1
    px   = self%px  - 1
    py   = self%py  - 1
    pz   = self%pz  - 1
    bx   = self%bx  - 1
    by   = self%by  - 1
    bz   = self%bz  - 1
    p1   = self%px  - 1
    p2   = self%py  - 1
    p3   = self%pz  - 1
    b1   = self%bx  - 1
    b2   = self%by  - 1
    b3   = self%bz  - 1
    phi  = self%phi - 1
    dphi = self%dphi - 1
    write (io_unit%data, idx_nml)
  end if
END SUBROUTINE write

!> =============================================================================
!> Write out time derivative index value namelist, for python/dispatch/
!> =============================================================================
SUBROUTINE output_time_derivs (self)
  class(index_t):: self
  integer:: d, e, et, s, px, py, pz, bx, by, bz, phi, p1, p2, p3, b1, b2, b3, dphi
  integer:: dpxdt, dpydt, dpzdt, dbxdt, dbydt, dbzdt, dphidt, dddt, dedt, dsdt
  integer:: dp1dt, dp2dt, dp3dt, db1dt, db2dt, db3dt
  namelist /idx_nml/ d, e, s, phi, dphi, &
                     dddt, dedt, dsdt, dphidt, &
                     px, py, pz, bx, by, bz, &
                     p1, p2, p3, b1, b2, b3, &
                     dpxdt, dpydt, dpzdt, dbxdt, dbydt, dbzdt, &
                     dp1dt, dp2dt, dp3dt, db1dt, db2dt, db3dt
  d   = self%d   - 1
  e   = self%e   - 1
  s   = self%s   - 1
  et  = self%et  - 1
  px  = self%px  - 1
  py  = self%py  - 1
  pz  = self%pz  - 1
  bx  = self%bx  - 1
  by  = self%by  - 1
  bz  = self%bz  - 1
  p1  = self%px  - 1
  p2  = self%py  - 1
  p3  = self%pz  - 1
  b1  = self%bx  - 1
  b2  = self%by  - 1
  b3  = self%bz  - 1
  phi = self%phi - 1
  !-----------------------------------------------------------------------------
  ! These values are preliminary, and are reassigned in python/dispatch/, based
  ! on knowledge of the mhd and selfgravity swwitches (not available here)
  !-----------------------------------------------------------------------------
  dddt   = self%d   + 5
  dedt   = self%e   + 5
  dsdt   = self%s   + 5
  dpxdt  = self%px  + 5
  dpydt  = self%py  + 5
  dpzdt  = self%pz  + 5
  dbxdt  = self%bx  + 8
  dbydt  = self%by  + 8
  dbzdt  = self%bz  + 8
  dphidt = self%phi + 9
  write (io_unit%data, idx_nml)
END SUBROUTINE output_time_derivs

END MODULE index_mod
