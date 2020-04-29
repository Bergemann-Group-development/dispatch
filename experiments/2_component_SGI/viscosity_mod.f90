!===============================================================================
!> This is a start on a module that could become a general viscosity module
!===============================================================================
MODULE viscosity_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE force_mod
  USE parameters_mod
  USE scaling_mod
  USE index_mod
  USE vector_ops
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: viscosity_t
    integer:: id
    integer:: i_gas, i_dust, nc, i=1
    character(len=64):: solver
    real(4):: nu
    real(8), dimension(:,:,:), pointer:: cs
    procedure(turbulent_viscosity), pointer:: selected=>null()
    procedure(dust_diffusivity), pointer:: selected2=>null()
    procedure(dust_dispersion), pointer:: selected3=>null()
  contains
    procedure init
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Initialize viscosity_mod
!===============================================================================
SUBROUTINE init (self, solver)
   class(viscosity_t):: self
   character(len=64):: solver
   character(len=32), save:: type='zero', type2='zero', type3='zero'
   integer, save:: i_gas=1, i_dust=2
   namelist /viscosity_params/ type, type2, type3, i_gas, i_dust
   !----------------------------------------------------------------------------
   call trace_begin('viscosity_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, viscosity_params)
     if (io%master) write (*, viscosity_params)
   end if
   call io%lock%unset
   self%solver = trim(solver)
   self%nu     = experiment_params%nu
   self%nc     = experiment_params%nc
   self%i_gas  = i_gas
   self%i_dust = i_dust
   select case (type)
   case ('turbulent_viscosity')
     self%selected => turbulent_viscosity
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN VISCOSITY TYPE'
   end select
   select case (type3)
   case ('dust_dispersion')
     self%selected3 => dust_dispersion
   case ('zero')
     nullify(self%selected3)
   case default
     print*,'UNKNOWN DISPERSION TYPE'
   end select
   select case (type2)
   case ('dust_diffusivity')
     self%selected2 => dust_diffusivity
   case ('zero')
     nullify(self%selected2)
   case default
     print*,'UNKNOWN DIFFUSIVITY TYPE'
   end select
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Apply turbulent viscosity parameter
!===============================================================================
FUNCTION turbulent_viscosity (self, time, mem, U, cs, m) RESULT (nu)
  class(viscosity_t)    :: self
  real(8)               :: time, r, omega
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: nu
  real, dimension(:,:,:):: cs
  real, dimension(:,:,:,:):: mem, U
  class(mesh_t), pointer :: m1, m2, m3
  integer:: i1, i2, i3
  !-----------------------------------------------------------------------------
  call trace_begin('viscosity_t:turbulent_viscosity')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  if (self%i==self%i_gas) then
    do i2 = m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      omega = sqrt(scaling%grav*experiment_params%mstar/r**3)
      do i3 = m3%lb,m3%ub
        do i1 = m1%lb,m1%ub
          nu(i1,i2,i3) = self%nu*cs(i1,i2,i3)**2 / omega
        end do
      end do
    end do
  else if (self%i==self%i_dust) then
    nu = 0.0
  end if
  nullify(m1,m2,m3)
  if (.not.associated(self%selected2) .and. .not.associated(self%selected3)) then
    self%i = self%i + 1
    if (self%i > self%nc) self%i=1
  end if
  call trace_end
END FUNCTION turbulent_viscosity


FUNCTION dust_diffusivity (self, time, mem, U, m, force) RESULT (dddt)
  class(viscosity_t)     :: self
  real(8)                :: time, r, omega, Dt
  class(mesh_t), pointer, dimension(:) :: m
  type(force_t):: force
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: dddt, t_stop, D, pg
  real, dimension(:,:,:), pointer:: d_gas, d_dust, cs, s
  real, dimension(:,:,:,:):: mem, U
  class(mesh_t), pointer :: m1, m2, m3
  integer:: i1, i2, i3, ic
  !-----------------------------------------------------------------------------
  call trace_begin('viscosity_t:dust_diffusivity')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  if (self%i==self%i_dust) then
    allocate(d_gas(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(d_dust(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(cs(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(s(m(1)%gn,m(2)%gn,m(3)%gn))
    ic     = 5*(self%i_gas-1)
    d_gas  = mem(:,:,:,idx%d+ic)
    s      =  mem(:,:,:,idx%s +ic        )
    if (force%gamma==1.0) then
      pg    = d_gas*force%csound**2
    else
      pg    = d_gas**force%gamma*exp(s/d_gas*(force%gamma-1.0))
    end if
    cs     = sqrt((force%gamma*pg)/d_gas)
    ic     = 5*(self%i_dust-1)
    d_dust = mem(:,:,:,idx%d+ic)
    t_stop = force%stopping_time(d_gas,cs,m)
    do i2 = m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      omega = sqrt(scaling%grav*experiment_params%mstar/r**3)
      do i3 = m3%lb,m3%ub
        do i1 = m1%lb,m1%ub
          Dt = (1.0 +  t_stop(i1,i2,i3)*omega + 4.0*(t_stop(i1,i2,i3)*omega)**2) &
                      / (1.0 + (t_stop(i1,i2,i3)*omega)**2)**2
          D(i1,i2,i3) = Dt * self%nu * cs(i1,i2,i3)**2 / omega
        end do
      end do
    end do
    dddt = D * laplacian(m,d_dust)
    deallocate(d_gas)
    deallocate(d_dust)
    deallocate(cs)
    deallocate(s)
  else if (self%i==self%i_gas) then
    dddt = 0.0
  end if
  nullify(m1,m2,m3)
  if (.not.associated(self%selected3)) then
    self%i = self%i + 1
    if (self%i > self%nc) self%i=1
  end if
  call trace_end
END FUNCTION dust_diffusivity

FUNCTION dust_dispersion (self, time, mem, U, m, force) RESULT (dpdt)
  class(viscosity_t)     :: self
  real(8)                :: time, r, omega, ct2
  class(mesh_t), pointer, dimension(:) :: m
  type(force_t):: force
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: dpdt
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: t_stop, cd2, pg
  real, dimension(:,:,:), pointer:: d_gas, d_dust, cs, s
  real, dimension(:,:,:,:):: mem, U
  class(mesh_t), pointer :: m1, m2, m3
  integer:: i1, i2, i3, ic
  !-----------------------------------------------------------------------------
  call trace_begin('viscosity_t:dust_dispersion')
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  if (self%i==self%i_dust) then
    allocate(d_gas(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(d_dust(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(cs(m(1)%gn,m(2)%gn,m(3)%gn))
    allocate(s(m(1)%gn,m(2)%gn,m(3)%gn))
    ic     = 5*(self%i_gas-1)
    d_gas  = mem(:,:,:,idx%d+ic)
    s      =  mem(:,:,:,idx%s +ic        )
    if (force%gamma==1.0) then
      pg    = d_gas*force%csound**2
    else
      pg    = d_gas**force%gamma*exp(s/d_gas*(force%gamma-1.0))
    end if
    cs     = sqrt((force%gamma*pg)/d_gas)
    ic     = 5*(self%i_dust-1)
    d_dust = mem(:,:,:,idx%d+ic)
    t_stop = force%stopping_time(d_gas,cs,m)
    do i2 = m2%lb,m2%ub
      r = m2%centre_nat + m2%r(i2)
      omega = sqrt(scaling%grav*experiment_params%mstar/r**3)
      do i3 = m3%lb,m3%ub
        do i1 = m1%lb,m1%ub
          ct2 = (1.0 +  2.0*(t_stop(i1,i2,i3)*omega)**2 + (5.0/4.0)*(t_stop(i1,i2,i3)*omega)**3) &
                        / (1.0 + (t_stop(i1,i2,i3)*omega)**2)**2
          cd2(i1,i2,i3) = ct2 * self%nu * cs(i1,i2,i3)**2
        end do
      end do
    end do
    dpdt = - down(cd2) * vector_scalar_grad(m,d_dust)
    deallocate(d_gas)
    deallocate(d_dust)
    deallocate(cs)
    deallocate(s)
  else if (self%i==self%i_gas) then
    dpdt = 0.0
  end if
  nullify(m1,m2,m3)
  self%i = self%i + 1
  if (self%i > self%nc) self%i=1
  call trace_end
END FUNCTION dust_dispersion
END MODULE
