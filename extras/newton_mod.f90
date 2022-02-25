!===============================================================================
!> Add Newton cooling to any solver.  Equation for heating per unit volume:´
!>
!> $$ Q = (E(s) - E_Q - E_s s)/t_Q f/(1+f) $$
!>
!> where f = exp((s-s_0)/s_1) and f/(1+f) is thus a function that turns on from
!> 0 to 1 near s_0, over a distance s_1, and the heating/cooling rate is 1/t_Q
!> well beyond s_0.  The E_s term alloww cooling towards a slope.
!>
!> The name of the corresponding namelist parameters are:t
!>
!> s_0: position
!> s_1: scale
!> E_0: ee0
!> E_1: ee1
!> t_Q: time
!===============================================================================
MODULE newton_mod
  USE io_mod
  USE link_mod
  USE gpatch_mod
  USE mesh_mod
  USE io_unit_mod
  USE trace_mod
  USE bits_mod
  USE scratch_mod
  USE shared_mod
  implicit none
  private
  type, public:: newton_t
    real(8) :: ee0       = 0.0
    real(8) :: ee1       = 0.0
    real(8) :: time      = 1.0
    real(8) :: scale     = 1.0
    real(8) :: position  = 0.0
    logical :: on        = .false.
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: get
    procedure:: set
    procedure:: pre_update
  end type
  integer, save:: verbose=0, axis=3
  ! -- the default values are reasonable for the solar photosphere
  real, save:: position=0.5, scale=0.1, ee0=0.0, ee1=-0.5, time=0.1
CONTAINS

!===============================================================================
!> Initialize Newton cooling
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(newton_t):: self
  class(gpatch_t), target:: patch
  ! ............................................................................
  integer:: iostat, m(3)
  logical, save:: on=.false.
  real(8), save:: position=-99., scale=1, ee0=0.0, ee1=0.5, time=0.1
  logical, save:: first_time=.true.
  namelist /newton_params/ on, axis, position, scale, ee0, ee1, time, &
    verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('newton_t%pre_init')
  self%patch => patch
  call io%lock%set ('newton_t%pre_init')
  if (first_time) then
    rewind (io_unit%input)
    read (io%input, newton_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io%input, newton_params)
    end if
    call io%namelist_warning ('newton_params', iostat)
    write (io%namelist, newton_params)
    first_time = .false.
  end if
  call io%lock%unset ('newton_t%pre_init')
  self%on = on
  self%time = time
  self%scale = scale
  self%ee0 = ee0
  self%ee1 = ee1
  self%position = position
  call trace%end()
END SUBROUTINE pre_init

SUBROUTINE get (self, ee0, ee1, position, time)
  class(newton_t)   :: self
  real(8), optional :: ee0, ee1, position, time
  !-----------------------------------------------------------------------------
  if (present(time)) time = self%time
  if (present(ee0)) ee0 = self%ee0
  if (present(ee1)) ee1 = self%ee1
  if (present(position)) position = self%position
END SUBROUTINE get

SUBROUTINE set (self, ee0, ee1, position, time)
  class(newton_t):: self
  real(8), optional :: ee0, ee1, position, time
  !-----------------------------------------------------------------------------
  self%time = shared%option (time, self%time)
  self%ee0 = shared%option (ee0, self%ee0)
  self%ee1 = shared%option (ee1, self%ee1)
  self%position = shared%option (position, self%position)
END SUBROUTINE set

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE init (self)
  class(newton_t):: self
  ! ............................................................................
END SUBROUTINE init

!===============================================================================
!> Apply Newton cooling
!===============================================================================
SUBROUTINE pre_update (self)
  class(newton_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, newton
  real:: f, r, ee, s, ee0_loc, em, ek, eth, b1, b2, b3
  integer:: ix, iy, iz, lo(3), uo(3), mprint
  integer, save:: nprint=5
  ! ----------------------------------------------------------------------------
  if (self%on) then
    if (verbose > 0) then
      !$omp atomic read
      mprint = nprint
      if (mprint > 0 .and. mod(self%patch%id,5100)==1) then
        !$omp atomic capture
        nprint = nprint-1
        mprint = nprint
        !$omp end atomic
      end if
      if (mprint > 0 .and. mod(self%patch%id,5100)==1) then
        PRINT *, self%patch%id, trim(self%patch%etype)
      end if
    end if
    !...........................................................................
    m1 => self%patch%mesh(1)
    m2 => self%patch%mesh(2)
    m3 => self%patch%mesh(3)
    call scratch%allocate (scratch%heating_per_unit_volume, self%patch%gn)
    if (trim(self%patch%etype)=='thermal') then
      d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
      e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
      if (ee0 /= 0.0) then
        ee0_loc = ee0
      else
        ee0_loc = self%ee0
      end if
      lo=self%patch%mesh%lb
      uo=self%patch%mesh%ub
      if (axis==1) then
        do iz = lo(3), uo(3)
          do iy = lo(2), uo(2)
            do ix = lo(1), uo(1)
              s = self%position - (m1%p + m1%r(ix))
              f = exp(min(-s/self%scale,50d0))
              r = f/((1.0+f)*self%time)
              ee = ee0_loc - self%ee1*s
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-e(ix,iy,iz))*r
            end do
          end do
        end do
      else
        do iz = lo(3), uo(3)
          if (axis==3) then
            s = self%position - (m3%p + m3%r(iz))
            f = exp(min(-s/self%scale,50d0))
            r = f/((1.0+f)*self%time)
            ee = ee0_loc - self%ee1*s
            if (verbose > 0) then
              if (mprint > 0 .and. mod(self%patch%id,5100)==1) then
                write (stdout,'(a,i6,i4,2x,1p,6g12.3)') &
                  'newton: id,iz,pos,rad,s,f,r,ee=', &
                  self%patch%id, iz, self%position, m3%p+m3%r(iz), s, f, r, ee
              end if
            end if
          end if
          do iy = lo(2), uo(2)
            if (axis==2) then
              s = self%position - (m2%p + m2%r(iy))
              f = exp(min(-s/self%scale,50d0))
              r = f/((1.0+f)*self%time)
              ee = ee0_loc - self%ee1*s
            end if
            do ix = lo(1), uo(1)
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-e(ix,iy,iz))*r
            end do
          end do
        end do
      end if
    else if (trim(self%patch%etype)=='total') then
      lo=self%patch%mesh%lo-1
      uo=self%patch%mesh%uo+1
      d  => self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%it,1)
      e  => self%patch%mem(:,:,:,self%patch%idx%e ,self%patch%it,1)
      px => self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1)
      py => self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1)
      pz => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1)
      bx => self%patch%mem(:,:,:,self%patch%idx%bx,self%patch%it,1)
      by => self%patch%mem(:,:,:,self%patch%idx%by,self%patch%it,1)
      bz => self%patch%mem(:,:,:,self%patch%idx%bz,self%patch%it,1)
      select case(axis)
      case(1)
        do iz = lo(3), uo(3)
          do iy = lo(2), uo(2)
            do ix = lo(1), uo(1)
              s = self%position - (m1%p + m1%r(ix))
              f = exp(min(-s/self%scale,50d0))
              r = f/((1.0+f)*self%time)
              ee = ee0_loc - self%ee1*s
              b1 = 0.5*(bx(ix,iy,iz)+bx(ix+1,iy,iz))
              b2 = 0.5*(by(ix,iy,iz)+by(ix,iy+1,iz))
              b3 = 0.5*(bz(ix,iy,iz)+bz(ix,iy,iz+1))
              em = 0.5*(b1**2+b2**2+b3**2)
              ek = 0.5*(px(ix,iy,iz)**2+py(ix,iy,iz)**2+pz(ix,iy,iz)**2)/d(ix,iy,iz)
              eth = e(ix,iy,iz)-em-ek
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-eth)*r
            end do
          end do
        end do
      case(2)
        do iz = lo(3), uo(3)
          do iy = lo(2), uo(2)
              s = self%position - (m2%p + m2%r(iy))
              f = exp(min(-s/self%scale,50d0))
              r = f/((1.0+f)*self%time)
              ee = ee0_loc - self%ee1*s
            do ix = lo(1), uo(1)
              b1 = 0.5*(bx(ix,iy,iz)+bx(ix+1,iy,iz))
              b2 = 0.5*(by(ix,iy,iz)+by(ix,iy+1,iz))
              b3 = 0.5*(bz(ix,iy,iz)+bz(ix,iy,iz+1))
              em = 0.5*(b1**2+b2**2+b3**2)
              ek = 0.5*(px(ix,iy,iz)**2+py(ix,iy,iz)**2+pz(ix,iy,iz)**2)/d(ix,iy,iz)
              eth = e(ix,iy,iz)-em-ek
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-eth)*r
            end do
          end do
        end do
      case(3)
        do iz = lo(3), uo(3)
          s = self%position - (m3%p + m3%r(iz))
          f = exp(min(-s/self%scale,50d0))
          r = f/((1.0+f)*self%time)
          ee = ee0_loc - self%ee1*s
          do iy = lo(2), uo(2)
            do ix = lo(1), uo(1)
              b1 = 0.5*(bx(ix,iy,iz)+bx(ix+1,iy,iz))
              b2 = 0.5*(by(ix,iy,iz)+by(ix,iy+1,iz))
              b3 = 0.5*(bz(ix,iy,iz)+bz(ix,iy,iz+1))
              em = 0.5*(b1**2+b2**2+b3**2)
              ek = 0.5*(px(ix,iy,iz)**2+py(ix,iy,iz)**2+pz(ix,iy,iz)**2)/d(ix,iy,iz)
              eth = e(ix,iy,iz)-em-ek
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-eth)*r
            end do
          end do
        end do
      end select
    else
      call io%abort('newton_mod:: etype not known'//trim(self%patch%etype))
    end if
    newton => scratch%heating_per_unit_volume
    call self%patch%aux%register('newton',newton)
  end if
END SUBROUTINE pre_update

END MODULE newton_mod
