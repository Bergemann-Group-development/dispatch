!===============================================================================
!> This is a start on a module that could become a general coupling force
!> between dust and gas fluids when using the stagger2_multifluid solver
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE index_mod
  USE vector_ops
  USE scalar_mod
  implicit none
  private
  type, public:: force_t
    integer:: seed=-11
    integer:: i_gas, i_dust, nc, i=1
    real::    K, d_int, a, gamma, csound, b, sigma, tmin
    real::    pi=4.0*atan(1.0)
    real, dimension(:,:,:,:), allocatable:: ff_drag
    character(len=32):: drag
    character(len=64):: solver
    procedure(dust_gas_interaction), pointer:: selected=>null()
  contains
    procedure init
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
   character(len=64):: solver
   integer, save:: seed=-11, iostat
   integer, save:: i_gas=1, i_dust=2, nc=2
   real, save:: K=0.5, d_int=3.0, sigma=1.589e-8, a=0.3, csound=1.0, gamma=1.0, b=0.5
   logical, save:: first_time = .true.
   character(len=32), save:: type='zero', drag='linear'
   namelist /force_params/ type, i_gas, i_dust, K, d_int, a, nc, gamma, &
                           drag, csound, b, sigma
   !----------------------------------------------------------------------------
   call trace%begin('force_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params)
     if (io%master) write (*, force_params)
   end if
   call io%lock%unset
   self%solver = trim(solver)
   self%drag   = trim(drag)
   self%seed   = seed
   self%i_gas  = i_gas
   self%i_dust = i_dust
   self%K      = K
   self%d_int  = d_int
   self%sigma  = sigma
   self%a      = a
   self%b      = b
   self%csound = csound
   self%nc     = nc
   self%gamma  = gamma
   select case (type)
   case ('dust_gas_interaction')
     self%selected => dust_gas_interaction
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   call trace%end()
END SUBROUTINE init

!===============================================================================
!> Coupling force between dust and gas fluids
!===============================================================================
FUNCTION dust_gas_interaction (self, time, d, mem, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer, contiguous :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  real, dimension(:,:,:,:), pointer:: p, U_gas, U_dust, dd, ld, KK
  real, dimension(:,:,:), pointer:: d_gas, d_dust, lnd, s, pg, t_stop, cs, &
                                    v_therm, K, dU2, mu, nu, R, cd
  real, dimension(:,:,:):: Ux, Uy, Uz, d
  integer:: ic, nv
  integer, dimension(3):: gn
  !.............................................................................
  integer:: ix, iy, iz
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%dust_gas_interaction')

  gn(1)=size(mem(:,1,1,1)); gn(2)=size(mem(1,:,1,1)); gn(3)=size(mem(1,1,:,1))
  nv = size(mem(1,1,1,:))
  !-----------------------------------------------------------------------------
  ! Start calculation of drag force if the current fluid component is the
  ! first of the coupled pair of fluids
  !-----------------------------------------------------------------------------
  if (self%i==min(self%i_gas,self%i_dust)) then
    !---------------------------------------------------------------------------
    ! Allocate vector and scalar pointers
    !---------------------------------------------------------------------------
    call allocate_scalars(gn, d_gas, d_dust, lnd, pg, s, cs, t_stop, v_therm, K)
    call allocate_scalars(gn, dU2, mu, nu, R, cd)
    call allocate_vectors(gn, p, ld, dd, U_gas, U_dust, KK)
    !---------------------------------------------------------------------------
    ! Compute velocity=momentum/density for the gas fluid component, leaving
    ! gas density for reuse
    !---------------------------------------------------------------------------
    ic    =  5*(self%i_gas-1)
    d_gas =  mem(:,:,:,idx%d +ic        )
    p     =  mem(:,:,:,idx%px+ic:idx%pz+ic)
    s     =  mem(:,:,:,idx%s +ic        )
    lnd   =  log(d_gas)
    ld    =  down(lnd)
    dd    =  exp(ld)
    U_gas =  p/dd
    !--------------------------------------------------------------------------
    ! Compute velocity=momentum/density for the dust fluid component, leaving
    ! dust density for reuse
    !--------------------------------------------------------------------------
    ic     =  5*(self%i_dust-1)
    d_dust =  mem(:,:,:,idx%d +ic        )
    p      =  mem(:,:,:,idx%px+ic:idx%pz+ic)
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
        pg    = d_gas**self%gamma*exp(s/d_gas*(self%gamma-1))
      end if
      cs      = sqrt((self%gamma*pg)/d_gas)
      v_therm = sqrt(8.0/(self%pi*self%gamma))*cs
      t_stop  = (self%d_int*self%a)/(d_gas*v_therm)
      K       = (d_dust/t_stop) * sqrt( 1 + (9.0*self%pi/128.0)*(dU2/(cs**2)) )
    else if (self%drag=='stokes') then !STOKES REGIME MIGHT BE BROKEN
      if (self%gamma==1.0) then
        pg = d_gas*self%csound**2
      else
        pg = d_gas**self%gamma*exp(s/d_gas*(self%gamma-1))
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
    ! Save drag force for use on the second component, and update counter
    !--------------------------------------------------------------------------
    allocate(self%ff_drag(m(1)%gn,m(2)%gn,m(3)%gn,3))
    self%ff_drag = ff
    lnd   =  log(d)
    ld    =  down(lnd)
    dd    =  exp(ld)
    ff    = self%ff_drag/dd
    self%i= self%i + 1
    !--------------------------------------------------------------------------
    ! Deallocate pointers
    !--------------------------------------------------------------------------
    call deallocate_scalars (d_gas, d_dust, lnd, pg, s, cs, t_stop, v_therm, K)
    call deallocate_scalars (dU2, mu, nu, R, cd)
    call deallocate_vectors (p, ld, dd, U_gas, U_dust, KK)
  !----------------------------------------------------------------------------
  ! Return reverse drag force if the current fluid component is the second one
  ! of the coupled pair of fluids. Then update counter
  !----------------------------------------------------------------------------
  else if (self%i==max(self%i_gas,self%i_dust)) then
    !--------------------------------------------------------------------------
    ! Allocate vector and scalar pointers
    !--------------------------------------------------------------------------
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
    self%i       =  self%i + 1
  !----------------------------------------------------------------------------
  ! Return no drag force if the current fluid component is not part of the
  ! coupled pair. Then update counter
  !----------------------------------------------------------------------------
  else
    ff = 0d0
    self%i = self%i + 1
  end if
  !----------------------------------------------------------------------------
  ! Reset internal fluid component counter self%i if it has surpassed number of
  ! fluid components self%nc
  !----------------------------------------------------------------------------
  if (self%i > self%nc) then
    self%i=1
    deallocate(self%ff_drag)
  end if
  call trace_end
END FUNCTION dust_gas_interaction

END MODULE
