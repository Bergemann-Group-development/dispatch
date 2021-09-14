!===============================================================================
!> Add Newton cooling to any solver
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
  implicit none
  private
  type, public:: newton_t
    real:: ee_top = 0.0
    real:: position=0.5
    logical:: on=.false.
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
  end type
  integer, save:: verbose=0, axis=3
  ! -- the default values are reasonable for the solar photosphere
  real, save:: position=-99., scale=0.1, ee0=0.0, ee1=0.5, time=0.1
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
  logical, save:: first_time=.true.
  namelist /newton_params/ on, axis, position, scale, ee0, ee1, time
  !-----------------------------------------------------------------------------
  call trace%begin ('newton_t%pre_init')
  self%patch => patch
  call io%lock%set ('newton_t%pre_init')
  if (first_time) then
    rewind (io_unit%input)
    read (io%input, newton_params, iostat=iostat)
    call io%namelist_warning ('newton_params', iostat)
    write (io%output, newton_params)
    write (io%namelist, newton_params)
    first_time = .false.
  end if
  call io%lock%unset ('newton_t%pre_init')
  self%on = on
  self%position=position
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize Newton cooling, allocating heating array knowing patch%gn
!===============================================================================
SUBROUTINE init (self)
  class(newton_t):: self
  ! ............................................................................
  integer:: m(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('newton_t%init')
  ! nothing to do here...
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Apply Newton cooling
!===============================================================================
SUBROUTINE pre_update (self)
  class(newton_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz
  real:: f, r, ee, s, ee_top, em, ek, eth, b1, b2, b3
  integer:: ix, iy, iz, lo(3), uo(3)
  ! ----------------------------------------------------------------------------
  if (self%on) then
    if (self%patch%kind(1:6)/='ramses') then
      d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
      e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
      call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
      if (ee0/=0.0) then
        ee_top = ee0
      else
        ee_top = self%ee_top
      end if
      m1 => self%patch%mesh(1)
      m2 => self%patch%mesh(2)
      m3 => self%patch%mesh(3)
      if      (axis==1) then
        do iz = m3%lo, m3%uo
          do iy = m2%lo, m2%uo
            do ix = m1%lo, m1%uo
              s = position - (m1%p + m1%r(ix))
              f = exp(-s/scale)
              r = f/((1.0+f)*time)
              ee = ee_top + ee1*s
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-e(ix,iy,iz))*r
            end do
          end do
        end do
      else
        do iz = m3%lo, m3%uo
          if (axis==3) then
            s = position - (m3%p + m3%r(iz))
            f = exp(-s/scale)/time
            r = f/((1.0+f)*time)
            ee = ee_top + ee1*s
          end if
          do iy = m2%lo, m2%uo
            if (axis==2) then
              s = position - (m2%p + m2%r(iy))
              f = exp(-s/scale)
              r = f/((1.0+f)*time)
              ee = ee_top + ee1*s
            end if
            do ix = m1%lo, m1%uo
              scratch%heating_per_unit_volume(ix,iy,iz) = &
              scratch%heating_per_unit_volume(ix,iy,iz) + &
                (ee*d(ix,iy,iz)-e(ix,iy,iz))*r
            end do
          end do
        end do
      end if
    else
      lo=self%patch%mesh%lo-1
      uo=self%patch%mesh%uo+1
      d  => self%patch%mem(:,:,:,self%patch%idx%d ,self%patch%new,1)
      e  => self%patch%mem(:,:,:,self%patch%idx%e ,self%patch%new,1)
      px => self%patch%mem(:,:,:,self%patch%idx%px,self%patch%new,1)
      py => self%patch%mem(:,:,:,self%patch%idx%py,self%patch%new,1)
      pz => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%new,1)
      bx => self%patch%mem(:,:,:,self%patch%idx%bx,self%patch%new,1)
      by => self%patch%mem(:,:,:,self%patch%idx%by,self%patch%new,1)
      bz => self%patch%mem(:,:,:,self%patch%idx%bz,self%patch%new,1)
      call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
      if (ee0/=0.0) then
        ee_top = ee0
      else
        ee_top = self%ee_top
      end if
      m1 => self%patch%mesh(1)
      m2 => self%patch%mesh(2)
      m3 => self%patch%mesh(3)
      select case(axis)
      case(1)
        do iz = lo(3), uo(3)
          do iy = lo(2), uo(2)
            do ix = lo(1), uo(1)
              s = position - (m1%p + m1%r(ix))
              f = exp(-s/scale)
              r = f/((1.0+f)*time)
              ee = ee_top + ee1*s
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
              s = position - (m2%p + m2%r(iy))
              f = exp(-s/scale)
              r = f/((1.0+f)*time)
              ee = ee_top + ee1*s
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
            s = position - (m3%p + m3%r(iz))
            f = exp(-s/scale)/time
            r = f/((1.0+f)*time)
            ee = ee_top + ee1*s
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
    end if
  end if
END SUBROUTINE pre_update

END MODULE newton_mod
