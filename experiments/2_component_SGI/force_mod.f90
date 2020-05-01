!===============================================================================
!> This is a start on a module that could become a general coupling force
!> between dust and gas fluids when using the stagger2_multifluid solver.
!> This force also contains the direct gravity from a central star at origin,
!> as well as a combined coupling and gravity force.
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE index_mod
  USE stagger_mod
  USE vector_ops
  USE scalar_mod
  USE scaling_mod
  USE star_mod
  USE potential_mod
  USE parameters_mod
  USE scalar_mod
  implicit none
  private
  type, public:: force_t
    integer:: i_gas, i_dust, nc, i=1
    real::    K, d_int, a, da, gamma, csound, b, sigma, tmin=1.0e6
    real::    pi=4.0*atan(1.0)
    real, dimension(:,:,:,:), allocatable:: ff_drag, ff_grav
    logical:: explicit_star
    character(len=32):: drag
    character(len=64):: solver
    procedure(star_and_drag),     pointer:: selected=>null()
  contains
    procedure init
    procedure dust_gas_interaction
    procedure stopping_time
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup initial values and decide what kind of drag regime to be used. The
!> options are: linear (default), quadratic, Epstein (see Kwon (1975) or
!> Laibe & Price (2012b)), power law, third-order expansion, and mixed
!===============================================================================
SUBROUTINE init (self, solver)
   class(force_t):: self
   class(mesh_t), pointer, dimension(:):: mesh
   character(len=64):: solver
   integer, save:: i_gas=1, i_dust=2
   real, save:: K=0.5, sigma=1.589e-8, csound=1.0, gamma=1.0, &
                b=0.5
   logical, save:: first_time = .true., explicit_star = .false.
   character(len=32), save:: type='zero', drag='linear'
   namelist /force_params/ type, i_gas, i_dust, K, gamma, &
                           drag, csound, b, sigma, explicit_star
   !----------------------------------------------------------------------------
   call trace%begin('force_t%init')
   self%nc = experiment_params%nc
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params)
     if (io%master) write (*, force_params)
   end if
   call io%lock%unset
   self%solver  = trim(solver)
   self%drag    = trim(drag)
   self%i_gas   = i_gas
   self%i_dust  = i_dust
   self%K       = K
   self%sigma   = sigma
   self%b       = b
   self%csound  = csound
   self%gamma   = gamma
   self%explicit_star = explicit_star
   self%da      = (experiment_params%d_int/scaling%d)*(experiment_params%a/scaling%l)
   call potential%init()
   select case (type)
   case ('star_and_drag')
     self%selected => star_and_drag
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   call trace%end()
END SUBROUTINE init

!>==============================================================================
!> Stopping time for the Epstein drag regime
!>==============================================================================
FUNCTION stopping_time (self, d_gas, cs, m) RESULT (t_stop)
  class(force_t):: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn):: t_stop, v_therm
  real, dimension(:,:,:):: d_gas, cs
  !-----------------------------------------------------------------------------
  call trace%begin('force_t%stopping_time')
  v_therm = sqrt(8.0/(self%pi*self%gamma))*cs
  t_stop  = (self%da)/(d_gas*v_therm)
  call trace%end()
END FUNCTION stopping_time

!===============================================================================
!> Coupling force between dust and gas fluids
!===============================================================================
FUNCTION dust_gas_interaction (self, time, d, mem, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer, contiguous :: mem
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  real, dimension(:,:,:,:), pointer:: p, U_gas, U_dust, dd, ld, KK
  !.............................................................................
  real, dimension(:,:,:), pointer:: d_gas, d_dust, lnd, s, pg, t_stop, cs, &
                                    v_therm, K, dU2, mu, nu, R, cd
  real, dimension(:,:,:):: Ux, Uy, Uz, d
  integer:: ic, nv
  integer, dimension(3):: gn
  integer:: ix, iy, iz
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%dust_gas_interaction')
  !!!i = idx
  gn(1)=size(mem,1); gn(2)=size(mem,2); gn(3)=size(mem,3)
  nv = size(mem,4)
  !--------------------------------------------------------------------------
  ! Start calculation of drag force if the current fluid component is the
  ! first of the coupled pair of fluids
  !--------------------------------------------------------------------------
  if (self%i==min(self%i_gas,self%i_dust)) then
    !----------------------------------------------------------------------------
    ! Allocate vector and scalar pointers
    !----------------------------------------------------------------------------
    call allocate_scalars(gn, d_gas, d_dust, lnd, pg, s, cs, t_stop, v_therm, K)
    call allocate_scalars(gn, dU2, mu, nu, R, cd)
    call allocate_vectors(gn, p, ld, dd, U_gas, U_dust, KK)
    !--------------------------------------------------------------------------
    ! Compute velocity=momentum/density for the gas fluid component, leaving
    ! gas density for reuse
    !--------------------------------------------------------------------------
    ic     =  5*(self%i_gas-1)
    d_gas  =  mem(:,:,:,idx%d +ic        )
    p      =  hinv_vector(m,mem(:,:,:,idx%px+ic:idx%pz+ic))
    s      =  mem(:,:,:,idx%s +ic        )
    lnd    =  log(d_gas)
    ld     =  down(lnd)
    dd     =  exp(ld)
    U_gas  =  p/dd
    !--------------------------------------------------------------------------
    ! Compute velocity=momentum/density for the dust fluid component, leaving
    ! dust density for reuse
    !--------------------------------------------------------------------------
    ic     =  5*(self%i_dust-1)
    d_dust =  mem(:,:,:,idx%d +ic        )
    p      =  hinv_vector(m,mem(:,:,:,idx%px+ic:idx%pz+ic))
    lnd    =  log(d_dust)
    ld     =  down(lnd)
    dd     =  exp(ld)
    U_dust =  p/dd
    !--------------------------------------------------------------------------
    ! Compute absolute velocity difference squared |dv|^2
    !--------------------------------------------------------------------------
    dU2    = norm2(U_gas - U_dust)
    !--------------------------------------------------------------------------
    ! Compute drag coefficient K for a given drag regime
    !--------------------------------------------------------------------------
    if      (self%drag=='linear') then
      K = self%K
    else if (self%drag=='epstein') then
      if (self%gamma==1.0) then
        pg    = d_gas*self%csound**2
      else
        pg    = d_gas**self%gamma*exp(s/d_gas*(self%gamma-1.0))
      end if
      cs      = sqrt((self%gamma*pg)/d_gas)
      t_stop  = stopping_time(self, d_gas, cs, m)
      K       = (d_dust/t_stop) * sqrt( 1.0 + (9.0*self%pi/128.0)*(dU2/(cs**2)) )
    else if (self%drag=='stokes') then !STOKES REGIME MIGHT BE BROKEN
      if (self%gamma==1.0) then
        pg = d_gas*self%csound**2
      else
        pg = d_gas**self%gamma*exp(s/d_gas*(self%gamma-1.0))
      end if
      cs   = sqrt((self%gamma*pg)/d_gas)
      mu   = (5.0/65.0)*self%sigma*sqrt(self%pi/self%gamma)*cs
      nu   = mu/d_gas
      R    = 2.0*self%a*sqrt(dU2)/nu
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            if (R(ix,iy,iz)>1) then
              cd(ix,iy,iz) = 24.0 * (R(ix,iy,iz)**(-1.0))
            else if (R(ix,iy,iz)>800) then
              cd(ix,iy,iz) = 0.44
            else
              cd(ix,iy,iz) = 24.0 * (R(ix,iy,iz)**(-0.6))
            end if
          end do
        end do
      end do
      K    = 0.5*cd*self%pi*(self%a**2)*d_gas*sqrt(dU2)
    else if (self%drag=='quadratic') then
      K = self%K * sqrt(dU2)
    else if (self%drag=='power') then
      K = self%K * ( sqrt(dU2) )**self%b
    else if (self%drag=='third_order') then
      K = self%K * ( 1.0 + self%b*dU2 )
    else if (self%drag=='mixed') then
      K = self%K * sqrt( 1.0 + self%b*dU2 )
    else
      K = self%K
    end if
    self%tmin = minval((d_gas*d_dust)/(K*(d_gas+d_dust)),K>0.0)
    KK = down(K)
    !--------------------------------------------------------------------------
    ! Compute drag force on the first relevant fluid component
    !--------------------------------------------------------------------------
    ff = KK*( U_gas - U_dust )
    if (self%i==self%i_gas) ff = -ff
    !--------------------------------------------------------------------------
    ! Save drag force for use on the second component
    !--------------------------------------------------------------------------
    allocate(self%ff_drag(m(1)%gn,m(2)%gn,m(3)%gn,3))
    self%ff_drag(:,:,:,1) =       ff(:,:,:,1)
    self%ff_drag(:,:,:,2) = h2c(m,ff(:,:,:,2))
    self%ff_drag(:,:,:,3) = h3c(m,ff(:,:,:,3))
    lnd   =  log(d)
    ld    =  down(lnd)
    dd    =  exp(ld)
    ff    = self%ff_drag/dd
    !--------------------------------------------------------------------------
    ! Deallocate pointers
    !--------------------------------------------------------------------------
    call deallocate_scalars (d_gas, d_dust, lnd, pg, s, cs, t_stop, v_therm, K)
    call deallocate_scalars (dU2, mu, nu, R, cd)
    call deallocate_vectors (p, ld, dd, U_gas, U_dust, KK)
  !----------------------------------------------------------------------------
  ! Return reverse drag force if the current fluid component is the second one
  ! of the coupled pair of fluids
  !----------------------------------------------------------------------------
  else if (self%i==max(self%i_gas,self%i_dust)) then
    !----------------------------------------------------------------------------
    ! Allocate vector and scalar pointers
    !----------------------------------------------------------------------------
    call allocate_scalars (gn, lnd)
    call allocate_vectors (gn, ld, dd)
    lnd   =  log(d)
    ld    =  down(lnd)
    dd    =  exp(ld)
    ff    = -self%ff_drag/dd
    !--------------------------------------------------------------------------
    ! Deallocate pointers
    !--------------------------------------------------------------------------
    call deallocate_scalars (lnd)
    call deallocate_vectors (ld, dd)
  !----------------------------------------------------------------------------
  ! Return no drag force if the current fluid component is not part of the
  ! coupled pair. Then update counter
  !----------------------------------------------------------------------------
  else
    ff     = 0d0
  end if
  call trace_end
END FUNCTION dust_gas_interaction
!===============================================================================
!> Implement both gravity from central star and dust-gas interaction
!===============================================================================
FUNCTION star_and_drag (self, time, d, mem, Ux, Uy, Uz, m) RESULT (ff)
  USE scalar_mod
  class(force_t):: self
  real(8):: time, dt
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:):: mem
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3):: ff
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn):: pot
  real, dimension(:,:,:):: Ux, Uy, Uz, d
  integer:: nv
  !-----------------------------------------------------------------------------
  call trace_begin('force_mod::star_and_drag')
  !-----------------------------------------------------------------------------
  ff = 0.0
  if (self%explicit_star) then
    ! Add the gravitational potential from the star once per time step
    if (self%i==1) then
      pot = potential%selected(m)
      allocate(self%ff_grav(m(1)%gn,m(2)%gn,m(3)%gn,3))
      self%ff_grav = - vector_scalar_grad(m,pot)
    end if
    ff = ff + self%ff_grav
  else
    if (self%i==1) then
      nv = size(mem,4)
      mem(:,:,:,nv) = mem(:,:,:,nv) + potential%selected(m)
    end if
  end if
  if (self%nc > 1) then
    ! Calculate drag force if more than 1 fluid is included
    ff = ff + dust_gas_interaction(self, time, d, mem, Ux, Uy, Uz, m)
  end if
  !ff =      dust_gas_interaction(self, time, d, mem, Ux, Uy, Uz, m)
  !-----------------------------------------------------------------------------
  ! Update internal fluid counter
  !-----------------------------------------------------------------------------
  self%i = self%i + 1
  !-----------------------------------------------------------------------------
  ! Reset internal fluid component counter self%i if it has surpassed number of
  ! fluid components self%nc
  !-----------------------------------------------------------------------------
  if (self%i > self%nc) then
    self%i=1
    if (self%nc > 1)   deallocate(self%ff_drag)
    if (self%explicit_star) deallocate(self%ff_grav)
  end if
  call trace_end
END FUNCTION star_and_drag
END MODULE
