!===============================================================================
!> Radiative transfer shadow test.
!> http://adsabs.harvard.edu/abs/2015A&A...574A..81R
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private

  type, public:: initial_t
    real:: clump_pos(3), clump_a, clump_b, d0, t0, clump_dfact, gamma
    character(len=64):: solver
    logical          :: mhd = .false.
  contains
    procedure :: init
    procedure :: condition
  end type

  logical, save:: first_time=.true.

CONTAINS
!===============================================================================
!> Read namelist
!===============================================================================
SUBROUTINE init (self, solver, gam)
  class(initial_t):: self
  character(len=64):: solver
  real, optional :: gam
  real, save:: clump_pos(3)=0.0, clump_a=0.1, clump_b=0.06, d0=1.0, t0=290.0, clump_dfact=10.0
  real, save:: gamma=1.4
  namelist /IC_params/ clump_pos, clump_a, clump_b, d0, t0, clump_dfact, gamma
  integer:: iostat
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, initial_params, iostat=iostat)
    if (io%master) write (*, initial_params)
    call scaling%init()
  end if
  call io%lock%unset
  self%solver = trim(solver)
  self%clump_pos = clump_pos/scaling%l
  self%clump_a = clump_a/scaling%l
  self%clump_b = clump_b/scaling%l
  self%d0 = d0 /scaling%d
  self%t0 = t0 /scaling%temp
  self%clump_dfact = clump_dfact
  self%gamma=gamma

  call trace_end
END SUBROUTINE init

SUBROUTINE condition (self, m, mem, i)
!===============================================================================
!> Set the initial conditions
!===============================================================================
  class(initial_t) :: self
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: mem
  type(index_t) :: i
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz
  real :: p0, q, x, y, z, d
  !----------------------------------------------------------------------------
  call trace%begin('initial_t::initial_conditions')

  mx => m(1)
  my => m(2)
  mz => m(3)

  ! note: I've assumed a unit system where k_B / mu * m_pro is 1.
  p0 = self%t0 * self%d0
  do iz=mz%lb,mz%ub
    z = mz%p + mz%r(iz)
    do iy=my%lb,my%ub
      y = my%p  + my%r(iy)
      do ix=mx%lb,mx%ub
        x = mx%p  + mx%r(ix)

        ! currently, the y and z coordinates are treated symmetrically; this
        ! leads to a cylindrically symmetric problem in 3-D.
        q = ((x - self%clump_pos(1)) / self%clump_a)**2 &
          + ((y - self%clump_pos(2)) / self%clump_b)**2 &
          + ((z - self%clump_pos(3)) / self%clump_b)**2
		q = 10.0 * (q - 1.0)
        ! the if-then-else is to prevent overflow when (x,y,z) is very from zero.
        if (q > 88.0) then
          mem(ix,iy,iz,i%d) = self%d0
        else
          mem(ix,iy,iz,i%d) = self%d0 + (self%clump_dfact * self%d0 - self%d0) &
                                  / (1.0 + exp(q))
        end if
        ! the clump is initially in pressure equilibrium with the ambient medium
        d = mem(ix,iy,iz,i%d)
        if (self%solver(1:7)=='ramses_') then
          mem(ix,iy,iz,i%e) = p0/(self%gamma-1.0)!*d
        else if (self%solver(1:9)=='stagger2e') then
          mem(ix,iy,iz,i%e) = p0/(self%gamma-1.0)!*d
        else if (self%solver(1:9)=='stagger2_') then
          mem(ix,iy,iz,i%s) = d*(log(p0)-self%gamma*log(d))/(self%gamma-1.0)
        else
          print *,trim(self%solver)
          error stop "initial_mod (1): solver type not known"
        end if

        ! momentum components
        mem(ix,iy,iz,i%px) = 0.0
        mem(ix,iy,iz,i%py) = 0.0
        mem(ix,iy,iz,i%pz) = 0.0
        ! magnetic field components
        if (self%mhd) then
          mem(ix,iy,iz,i%bx) = 0.0
          mem(ix,iy,iz,i%by) = 0.0
          mem(ix,iy,iz,i%bz) = 0.0
        end if
      end do
    end do
  end do
  nullify (mx, my, mz, mem)
  call trace%end()
END SUBROUTINE condition

END MODULE
