!===============================================================================
!> Add friction to any solver.  FIXME: Should be relative to the mean.
!===============================================================================
MODULE friction_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE mesh_mod
  USE scratch_mod
  implicit none
  private
  type, public:: friction_t
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
  end type
  logical, save:: on=.false.
  integer, save:: verbose=0, axis=3
  real,    save:: time=0.1
CONTAINS

!===============================================================================
!> Initialize Newton cooling
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(friction_t):: self
  class(gpatch_t), target:: patch
  ! ............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /friction_params/ on, axis, time
  !-----------------------------------------------------------------------------
  call trace%begin ('friction_t%pre_init')
  self%patch => patch
  call io%lock%set ('friction_t%pre_init')
  if (first_time) then
    rewind (io_unit%input)
    read (io%input, friction_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io%input, friction_params)
    end if
    call io%namelist_warning ('friction_params', iostat)
    write (io%namelist, friction_params)
    first_time = .false.
  end if
  call io%lock%unset ('friction_t%pre_init')
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE init (self)
  class(friction_t):: self
  ! ............................................................................
END SUBROUTINE init

!===============================================================================
!> Friction force
!===============================================================================
SUBROUTINE pre_update (self)
  class(friction_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  real, dimension(:,:,:), pointer:: pz
  integer:: ix, iy, iz, lo(3), uo(3)
  ! ----------------------------------------------------------------------------
  if (on) then
    m1 => self%patch%mesh(1)
    m2 => self%patch%mesh(2)
    m3 => self%patch%mesh(3)
    call scratch%allocate (scratch%force_per_unit_volume, self%patch%gn)
    if (trim(self%patch%kind(1:4))=='stag') then
      pz => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1)
      lo=self%patch%mesh%lb
      uo=self%patch%mesh%ub
      if (axis==3) then
        do iz = lo(3), uo(3)
          do iy = lo(2), uo(2)
            do ix = lo(1), uo(1)
              scratch%force_per_unit_volume(ix,iy,iz,3) = &
              scratch%force_per_unit_volume(ix,iy,iz,3) - pz(ix,iy,iz)/time
            end do
          end do
        end do
      end if
    else
      call io%abort('friction_mod:: solver kind not known'//trim(self%patch%kind))
    end if
  end if
END SUBROUTINE pre_update

END MODULE friction_mod
