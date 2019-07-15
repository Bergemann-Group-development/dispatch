!===============================================================================
!> Add Newton cooling to any solver
!===============================================================================
MODULE newton_mod
  USE link_mod
  USE gpatch_mod
  USE mesh_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
  type, public:: newton_t
    class(gpatch_t), pointer:: patch
  contains
    procedure:: init
    procedure:: pre_update
  end type
  logical:: on=.false.
  integer:: verbose=0, axis=3
  ! -- the default values are reasonable for the solar photosphere
  real:: position=-0.5, scale=0.1, ee0=5.3, ee1=0.5, time=0.1
CONTAINS

!===============================================================================
!> Initialize Newton cooling
!===============================================================================
SUBROUTINE init (self, link)
  class(newton_t):: self
  class(link_t), pointer :: link
  ! ............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /newton_params/ on, axis, position, scale, ee0, ee1, time
  !-----------------------------------------------------------------------------
  call trace%begin ('newton_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, newton_params, iostat=iostat)
    write (io_unit%output, newton_params)
  end if
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  associate (task => link%task)
  select type (task)
  class is (gpatch_t)
  self%patch => task
  if (.not.allocated(task%heating_per_unit_volume)) &
    allocate (task%heating_per_unit_volume(task%gn(1),task%gn(2),task%gn(3)))
  end select
  end associate
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Apply Newton cooling
!===============================================================================
SUBROUTINE pre_update (self)
  class(newton_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  real, dimension(:,:,:), pointer:: d, e
  real:: f, r, ee, s
  integer:: ix, iy, iz
  ! ----------------------------------------------------------------------------
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  if (on) then
    m1 => self%patch%mesh(1)
    m2 => self%patch%mesh(2)
    m3 => self%patch%mesh(3)
    if      (axis==1) then
      do iz = m3%li, m3%ui
        do iy = m2%li, m2%ui
          do ix = m1%li, m1%ui
            s = m1%p + m1%r(ix)-position
            f = exp(-s/scale)
            r = f/((1.0+f)*time)
            ee = ee0 + ee1*s
            self%patch%heating_per_unit_volume(ix,iy,iz) = &
            self%patch%heating_per_unit_volume(ix,iy,iz) + &
              (ee-e(ix,iy,iz)/d(ix,iy,iz))*d(ix,iy,iz)*r
          end do
        end do
      end do
    else
      do iz = m3%li, m3%ui
        if (axis==3) then
          s = m3%p + m3%r(iz)-position
          f = exp(-s/scale)/time
          r = f/((1.0+f)*time)
          ee = ee0 + ee1*s
        end if
        do iy = m2%li, m2%ui
          if (axis==2) then
            s = m2%p + m2%r(iy)-position
            f = exp(-s/scale)
            r = f/((1.0+f)*time)
            ee = ee0 + ee1*s
          end if
          do ix = m1%li, m1%ui
            self%patch%heating_per_unit_volume(ix,iy,iz) = &
            self%patch%heating_per_unit_volume(ix,iy,iz) + &
              (ee-e(ix,iy,iz)/d(ix,iy,iz))*d(ix,iy,iz)*r
          end do
        end do
      end do
    end if
  end if
END SUBROUTINE pre_update

END MODULE newton_mod
