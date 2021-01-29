!===============================================================================
!> This index file has slot indices for all solver, all initially equal to zero
!> It is the responsibility of the solver to set all values it needs to use.
!> NOTE: The only reference to indices in tasks should be via task%idx, never to
!> a temporary instance of index_t.
!===============================================================================
MODULE index_mod
  USE io_mod
  implicit none
  private
  type, public:: index_t
    integer:: d=0, e=0, px=0, py=0, pz=0, bx=0, by=0, bz=0, phi=0, tt=0, s=0, &
              et=0, dphi=0, rphi=0, qr=0, ex=0, ey=0, ez=0, jx=0, jy=0, jz=0
    integer:: dtr=0
  contains
    procedure:: init
    procedure:: copy
    procedure:: request_index
    procedure:: output
    procedure:: print
  end type
  type(index_t), public:: idx
CONTAINS

!===============================================================================
!> Initialize the indices to default values, possibly to be modified by solvers
!===============================================================================
SUBROUTINE init (self, ienergy, isentropy, mhd)
  class(index_t):: self
  integer, optional:: ienergy, isentropy
  logical, optional:: mhd
  integer:: ienergy_l, isentropy_l
  logical:: mhd_l
  !-----------------------------------------------------------------------------
  if (present(ienergy)) then
    ienergy_l = ienergy
  else
    ienergy_l = -1
  end if
  if (present(isentropy)) then
    isentropy_l = isentropy
  else
    isentropy_l = -1
  end if
  if (present(mhd)) then
    mhd_l = mhd
  else
    mhd_l = .false.
  end if
  if (ienergy_l==2 .or. isentropy_l==2) then
    self%d  = 1
    self%px = 3
    self%py = 4
    self%pz = 5
  else if (ienergy_l==5 .or. isentropy_l==5) then
    self%d  = 1
    self%px = 2
    self%py = 3
    self%pz = 4
  else
    call io%abort ('ienergy or isentropy != 2 or 5. You are required to supply one!')
  end if
  self%e  = ienergy_l
  self%s  = isentropy_l
  if (mhd_l) then
    self%bx = 6
    self%by = 7
    self%bz = 8
  else
    self%bx = -1
    self%by = -1
    self%bz = -1
  end if
  call self%copy
END SUBROUTINE init

!===============================================================================
!> Request a new, unique index
!> idx should be index itself, nv is patch%nv.
!===============================================================================
SUBROUTINE request_index (self, idx, nv, success)
  class(index_t):: self
  integer :: idx, nv
  logical, optional :: success
  !-----------------------------------------------------------------------------
  if (present(success)) success = .false.
  if (idx <= 0) then
    !$omp critical (index_cr)
    nv  = nv+1
    idx = nv
    if (present(success)) success = .true.
    !$omp end critical (index_cr)
    call self%copy
  end if
  io%nv = max(io%nv,nv)
END SUBROUTINE request_index

!===============================================================================
!> Initialize the indices to default values, possibly to be modified by solvers
!===============================================================================
SUBROUTINE copy (self)
  class(index_t):: self
  !-----------------------------------------------------------------------------
  ! Ugly fix, necessary because there are a large number of places that refer to
  ! a local instance of index_t, which is assumed to be authentative.  In cases
  ! where the default values have been modified by the server, this may lead to
  ! inconsistencies.  As a fall back, such places may instead refer to this one
  ! pubclic instance.
  !-----------------------------------------------------------------------------
  !$omp critical (index_cr)
  select type(self)
  type is (index_t)
    idx = self
  end select
  !$omp end critical (index_cr)
END SUBROUTINE copy

!===============================================================================
!> Write out index value namelist, for python/dispatch/
!===============================================================================
SUBROUTINE output (self, unit)
  class(index_t):: self
  integer:: unit
  !.............................................................................
  integer:: &
    d, e, et, s, px, py, pz, bx, by, bz, qr, tt, phi, p1, p2, p3, b1, b2, b3, &
    ex, ey, ez, jx, jy, jz, rphi
  namelist /idx_nml/ &
    d, e, et, s, px, py, pz, bx, by, bz, qr, tt, phi, p1, p2, p3, b1, b2, b3, &
    ex, ey, ez, jx, jy, jz, rphi
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    call output_time_derivs (self, unit)
  else
    d    = self%d   - 1
    e    = self%e   - 1
    s    = self%s   - 1
    et   = self%et  - 1
    qr   = self%qr  - 1
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
    qr   = self%qr  - 1
    tt   = self%tt  - 1
    phi  = self%phi - 1
    rphi = self%rphi- 1
    ex   = self%ex  - 1
    ey   = self%ey  - 1
    ez   = self%ez  - 1
    jx   = self%jx  - 1
    jy   = self%jy  - 1
    jz   = self%jz  - 1
    write (unit, idx_nml)
  end if
END SUBROUTINE output

!===============================================================================
!> Write out time derivative index value namelist, for python/dispatch/
!===============================================================================
SUBROUTINE output_time_derivs (self, unit)
  class(index_t):: self
  integer:: unit
  !.............................................................................
  integer:: &
    d, e, et, s, px, py, pz, bx, by, bz, qr, tt, phi, p1, p2, p3, b1, b2, b3, &
    dpxdt, dpydt, dpzdt, dbxdt, dbydt, dbzdt, dphidt, dddt, dedt, dsdt, &
    dp1dt, dp2dt, dp3dt, db1dt, db2dt, db3dt, ex, ey, ez, jx, jy, jz, &
    dexdt, deydt, dezdt, djxdt, djydt, djzdt, rphi
  namelist /idx_nml/ &
    d, e, et, s, px, py, pz, bx, by, bz, qr, tt, phi, p1, p2, p3, b1, b2, b3, &
    dpxdt, dpydt, dpzdt, dbxdt, dbydt, dbzdt, dphidt, dddt, dedt, dsdt, &
    dp1dt, dp2dt, dp3dt, db1dt, db2dt, db3dt, ex, ey, ez, jx, jy, jz, &
    dexdt, deydt, dezdt, djxdt, djydt, djzdt, rphi
  d   = self%d   - 1
  e   = self%e   - 1
  s   = self%s   - 1
  et  = self%et  - 1
  qr  = self%qr  - 1
  tt  = self%tt  - 1
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
  rphi = self%rphi- 1
  ex   = self%ex  - 1
  ey   = self%ey  - 1
  ez   = self%ez  - 1
  jx   = self%jx  - 1
  jy   = self%jy  - 1
  jz   = self%jz  - 1
  !-----------------------------------------------------------------------------
  ! These values are preliminary, and are reassigned in python/dispatch/, based
  ! on knowledge of the mhd and selfgravity swwitches (not available here)
  !-----------------------------------------------------------------------------
  dddt   = d   + io%nv
  dedt   = e   + io%nv
  dsdt   = s   + io%nv
  dpxdt  = px  + io%nv
  dpydt  = py  + io%nv
  dpzdt  = pz  + io%nv
  dphidt = phi + io%nv
  dp1dt  = dpxdt
  dp2dt  = dpydt
  dp3dt  = dpzdt
  if (bx > 0) then
    dbxdt  = bx + io%nv
    dbydt  = by + io%nv
    dbzdt  = bz + io%nv
  else
    dbxdt  = -1
    dbydt  = -1
    dbzdt  = -1
  end if
  if (ex > 0) then
    dexdt  = ex + io%nv
    deydt  = ey + io%nv
    dezdt  = ez + io%nv
  else
    dexdt  = -1
    deydt  = -1
    dezdt  = -1
  end if
  if (jx > 0) then
    djxdt  = jx + io%nv
    djydt  = jy + io%nv
    djzdt  = jz + io%nv
  else
    djxdt  = -1
    djydt  = -1
    djzdt  = -1
  end if
  db1dt  = dbxdt
  db2dt  = dbydt
  db3dt  = dbzdt
  write (unit, idx_nml)
END SUBROUTINE output_time_derivs

!===============================================================================
!> Write out index value namelist, for python/dispatch/
!===============================================================================
SUBROUTINE print (self, unit)
  class(index_t):: self
  integer:: unit
  !.............................................................................
  write (stdout,1) 'variable indices:   d' , self%d
  write (stdout,1) ' p1', self%px
  write (stdout,1) ' p2', self%py
  write (stdout,1) ' p3', self%pz
  if (self%e > 0) write (stdout,1) ' e', self%e
  if (self%s > 0) write (stdout,1) ' s', self%s
  if (idx%bx > 0) then
    write (stdout,1) ' b1', self%bx
    write (stdout,1) ' b2', self%by
    write (stdout,1) ' b3', self%bz
  end if
  if (self%phi  > 0) write (stdout,1) ' phi' , self%phi
  if (self%dphi > 0) write (stdout,1) ' dphi', self%dphi
  write (stdout,1) ''
  1 format(1x,a,':',i2,2x,$)
END SUBROUTINE print

END MODULE index_mod
