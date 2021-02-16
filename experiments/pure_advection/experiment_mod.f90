!===============================================================================
!> Initializiation and update of each ZEUS patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE kinds_mod
  USE mesh_mod
  USE mhd_mod
  USE index_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    integer:: ptype
    real:: pcentre(3), psize(3), pbufferfrac
    real:: pvx, pvy, pvz
    logical:: add_patch, pmoving
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type
  logical, save:: first_time = .true.
PUBLIC initial_conditions
CONTAINS

!===============================================================================
!> Initialise this experiment.
!===============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')

  self%nv = 8
  call self%mhd_t%init()
  call initial_conditions(self)

  call self%output

  call trace_end
END SUBROUTINE init_experiment

!===============================================================================
!> Update procedure for this experiment.
!===============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call self%mhd_t%update()

  call trace_end
END SUBROUTINE update_experiment

!===============================================================================
!> Initial conditions generator.
!===============================================================================
SUBROUTINE initial_conditions (self)
  class(experiment_t):: self
  class(mesh_t), pointer:: m1, m2, m3
  integer, save:: ptype=1
  real, save:: rho0=1.0, rho1=4.0, vx0=0.0, vy0=0.0, vz0=0.0, &
               vx1=0.0, vy1=0.0, vz1=0.0, p0=2.5, p1=2.5, bx0=0.0, by0=0.0, bz0=0.0, &
               bx1=0.0, by1=0.0, bz1=0.0, pbufferfrac=0.2
  real, dimension(3), save:: pcentre=0.0, psize=0.5
  logical, save:: add_patch=.true., pmoving=.true.
  integer:: i, j, k
  real(8):: x, y, z
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  namelist /experiment_params/ ptype, pcentre, psize, rho0, rho1, vx0, vy0, vz0, &
                               vx1, vy1, vz1, p0, p1, bx0, by0, bz0, bx1, by1, bz1, &
                               pbufferfrac, add_patch, pmoving
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::initial_conditions')

  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  call io%lock%unset

  self%ptype = ptype
  self%pcentre = pcentre
  self%psize = psize
  self%pbufferfrac = pbufferfrac
  self%add_patch = add_patch
  self%pmoving = pmoving
  self%pvx = vx1
  self%pvy = vy1
  self%pvz = vz1

  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  f => self%mem(:,:,:,:,self%it,1)

  select case (ptype)
  case (1) ! hydrostatic square pulse
    do k=m3%lb,m3%ub
      z = m3%p + m3%r(k)
      do j=m2%lb,m2%ub
        y = m2%p + m2%r(j)
        do i=m1%lb,m1%ub
          x = m1%p + m1%r(i)

          ! ambient
          f(i,j,k,self%idx%d ) = rho0
          f(i,j,k,self%idx%e ) = p0 / (eos%gamma - 1.0)
          f(i,j,k,self%idx%px) = rho0 * vx0
          f(i,j,k,self%idx%py) = rho0 * vy0
          f(i,j,k,self%idx%pz) = rho0 * vz0
          f(i,j,k,self%idx%bx) = bx0
          f(i,j,k,self%idx%by) = by0
          f(i,j,k,self%idx%bz) = bz0

          ! pulse
          if (x >= (pcentre(1)-psize(1)*0.5) .and. x <= (pcentre(1)+psize(1)*0.5) .and. &
              y >= (pcentre(2)-psize(2)*0.5) .and. y <= (pcentre(2)+psize(2)*0.5) .and. &
              z >= (pcentre(3)-psize(3)*0.5) .and. z <= (pcentre(3)+psize(3)*0.5) ) then
            f(i,j,k,self%idx%d ) = rho1
            f(i,j,k,self%idx%e ) = p1 / (eos%gamma - 1.0)
            f(i,j,k,self%idx%px) = rho1 * vx1
            f(i,j,k,self%idx%py) = rho1 * vy1
            f(i,j,k,self%idx%pz) = rho1 * vz1
            f(i,j,k,self%idx%bx) = bx1
            f(i,j,k,self%idx%by) = by1
            f(i,j,k,self%idx%bz) = bz1
          end if
          f(i,j,k,self%idx%et) = f(i,j,k,self%idx%e ) + 0.5 * rho0 * (vx0**2 + vy0**2 + vz0**2) &
                          + 0.5 * (f(i,j,k,self%idx%bx)**2 + f(i,j,k,self%idx%by)**2 + f(i,j,k,self%idx%bz)**2)
        end do
      end do
    end do

  case default
  end select

  self%periodic = .true.
  self%mhd = .false.

  nullify (m1, m2, m3, f)
  call trace_end
END SUBROUTINE initial_conditions

END MODULE experiment_mod
