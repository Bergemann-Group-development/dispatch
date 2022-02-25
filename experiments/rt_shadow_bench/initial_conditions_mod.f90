!===============================================================================
!> Radiative transfer shadow test.
!> http://adsabs.harvard.edu/abs/2015A&A...574A..81R
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE mesh_mod
  USE index_mod
  USE gpatch_mod
  USE io_unit_mod
  USE math_mod
  USE kinds_mod
  USE eos_mod
  implicit none
  private

  type, public:: initial_conditions_t
    real:: clump_pos(3), clump_a, clump_b, d0, t0, clump_dfact
    character(len=64):: solver
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure :: pre_init
    procedure :: init
  end type

  logical, save:: first_time=.true.

CONTAINS
!===============================================================================
!> Read namelist
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  real, save:: clump_pos(3)=0.5, clump_a=0.1, clump_b=0.06, d0=1.0, t0=290.0, clump_dfact=10.0
  namelist /IC_params/ clump_pos, clump_a, clump_b, d0, t0, clump_dfact
  integer:: iostat
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set('initial_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params)
    if (io%master) write (*, IC_params)
    call scaling%init()
  end if
  call io%lock%unset('initial_t%init')
  self%patch => patch
  self%solver = trim(patch%kind)
  self%clump_pos = clump_pos/scaling%l
  self%clump_a = clump_a/scaling%l
  self%clump_b = clump_b/scaling%l
  self%d0 = d0 /scaling%d
  self%t0 = t0 /scaling%temp
  self%clump_dfact = clump_dfact

  call trace_end
END SUBROUTINE pre_init

SUBROUTINE init (self)
!===============================================================================
!> Set the initial conditions
!===============================================================================
  class(initial_conditions_t) :: self
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: mem
  real, dimension(:,:,:), pointer :: tt, src
  type(index_t), pointer :: idx
  class(mesh_t), pointer :: mx, my, mz
  integer :: ix, iy, iz
  real :: p0, q, x, y, z, d, stefan
  !----------------------------------------------------------------------------
  call trace%begin('initial_t::initial_conditions')

  mx  => self%patch%mesh(1)
  my  => self%patch%mesh(2)
  mz  => self%patch%mesh(3)
  idx => self%patch%idx
  mem => self%patch%mem(:,:,:,:,1,1)
!  stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
!allocate(tt(self%patch%gn(1),self%patch%gn(2),self%patch%gn(3)))
!allocate(src(self%patch%gn(1),self%patch%gn(2),self%patch%gn(3)))
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
        q = sqrt(((x - self%clump_pos(1)) / self%clump_a)**2 &
          + ((y - self%clump_pos(2)) / self%clump_b)**2 &
          + ((z - self%clump_pos(3)) / self%clump_b)**2)
        !print *, q, x, y, z
        !q = 10.0 * (q - 1.0)
        ! the if-then-else is to prevent overflow when (x,y,z) is very from zero.
        if (q > 1.0) then
          d = self%d0+self%d0*sin(x)*cos(y)*cgs%pi/2
        else
          d = self%d0 + (self%clump_dfact * self%d0 - self%d0) &
                                  / (1.0 + exp(10*(q-1.0)))
        end if
        ! the clump is initially in pressure equilibrium with the ambient medium
        mem(ix,iy,iz,idx%d) = d
        if (idx%e>0) then
          mem(ix,iy,iz,idx%e) = p0/(eos%gamma-1.0)!*d
      !    tt(ix,iy,iz)= mem(ix,iy,iz,idx%e)/d*(eos%gamma-1.0)*scaling%temp
      !    src(ix,iy,iz)= stefan*real((mem(ix,iy,iz,idx%e)/d*(eos%gamma-1.0)*scaling%temp),kind=8)**4
        else if (idx%s>0) then
          mem(ix,iy,iz,idx%s) = d*(log(p0)-eos%gamma*log(d))/(eos%gamma-1.0)
      !    tt(ix,iy,iz)= d**eos%gamma*exp(mem(ix,iy,iz,idx%s)/d*(eos%gamma-1.0))/d*scaling%temp
      !    src(ix,iy,iz)= stefan*real(d**eos%gamma*exp(mem(ix,iy,iz,idx%s)/d*(eos%gamma-1.0))/d*scaling%temp,kind=8)**4
        else
          print *,trim(self%solver)
          error stop "initial_mod (1): solver type not known"
        end if
        

        ! momentum components
        mem(ix,iy,iz,idx%px) = 0.0
        mem(ix,iy,iz,idx%py) = 0.0
        mem(ix,iy,iz,idx%pz) = 0.0
        ! magnetic field components
        if (idx%bx>0) then
          mem(ix,iy,iz,idx%bx) = 0.0
          mem(ix,iy,iz,idx%by) = 0.0
          mem(ix,iy,iz,idx%bz) = 0.0
        end if
      end do
    end do
  end do
  !print *, 'IC minmax d', minval(mem(:,:,:,idx%d)), maxval(mem(:,:,:,idx%d))
  !print *, 'IC minmax tt', minval(tt), maxval(tt)
  !print *, 'IC minmax src', minval(src), maxval(src)
  nullify (mx, my, mz, mem,idx)
  !error stop "end"
  call trace%end()
END SUBROUTINE init

END MODULE
