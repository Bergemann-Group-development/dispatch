!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  !USE const
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  USE index_mod
  USE io_unit_mod
  USE math_mod
  USE kinds_mod
  USE pboundary_mod
  implicit none
  private
  type, public:: initial_t
    integer:: seed
    real:: gamma=1.4
    real:: d0
    logical:: mhd=.true.
    real(8):: time=0d0
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  !---------------------------------------------------------------------
  call trace_begin('initial_t%init')
  self%solver = trim(solver)
  self%gamma=1.4
  self%mhd = .true.
  call trace%end()
END subroutine init


SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:), pointer :: ff
  ! ....................................................................
  integer:: i, j, k
  real:: x, y, z, ds2(3), xs, rd
  logical:: first_time = .true.
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: d, s
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: e1, et
  real(kind=KindScalarVar):: rhov2, pxl
  real(4), save:: P_l=10.33333  , P_r=1.0, &
                  rho_l=3.857143, &
                  vx_l=2.629369, vx_r=0.0, &
                  a_d=0.2, f_d=5.0, x_disc = -4.0
  namelist /initial_params/ P_l, P_r, rho_l, &
                           vx_l, vx_r, a_d, f_d, x_disc
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, initial_params)
    if (io%master) write (*,initial_params)
    first_time = .false.
  end if
  call io%lock%unset
  !---------------------------------------------------------------------
  ! Point amd allocate memory
  !---------------------------------------------------------------------
    d => ff(:,:,:,idx%d)
  if (trim(self%solver(1:7)) == 'stagger') s => ff(:,:,:,idx%s)
  if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
    e1 => ff(:,:,:,idx%e)
  end if
  if (self%solver(1:9) == 'stagger2_') then
    s => ff(:,:,:,idx%s)
  end if
  if (self%solver(1:14) == 'zeus_mhd_patch') then
    et => ff(:,:,:,idx%et)
  end if
  if (self%solver(1:6) == 'ramses') then
    et => ff(:,:,:,idx%e)
    allocate (e1(m(1)%gn,m(2)%gn,m(3)%gn))
  end if
  ff(:,:,:,idx%py) = 0.0
  ff(:,:,:,idx%pz) = 0.0
  if (self%solver(1:7)=='stagger') then
    ds2 = 0.5*m%d
  else if (self%solver(1:6)=='zeus3d') then
    ds2 = 0.5*m%d
  else
    ds2 = 0.0
  endif
  ! --------------------------------------------------------------------
  ! Set initial conditions
  ! --------------------------------------------------------------------
  do i=m(1)%lb,m(1)%ub
    x = m(1)%p + m(1)%r(i)
    xs = x - ds2(1)
    if (x < x_disc) then
      d (i,:,:)= rho_l
      if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
        e1(i,:,:) = P_l / (self%gamma - 1.0)
      else if (self%solver(1:9) == 'stagger2_') then
        s (i,:,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
      else if (self%solver(1:6) == 'ramses')  then
        e1(i,:,:) = P_l / (self%gamma - 1.0)
      end if
      ff(i,:,:,idx%px) = vx_l * rho_l
    else
      d (i,:,:) = 1.0*(1.+a_d*sin(f_d*x))
      if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
        e1(i,:,:) = P_r / (self%gamma - 1.0)
      else if (self%solver(1:9) == 'stagger2_') then
        s (i,:,:) = (log(P_r) - log(d(i,:,:)) * self%gamma) * d(i,:,:) / (self%gamma - 1.0)
      else if (self%solver(1:6) == 'ramses')  then
        e1(i,:,:) = P_r / (self%gamma - 1.0)
      end if
      ff(i,:,:,idx%px) = 0.0
    end if
        ! 
    ! correction for staggered momentum
    if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:7) == 'stagger') then
      if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
        rd = 1.0*(1.+a_d*sin(f_d*x))
        ff(i,:,:,idx%px) = (rho_l + rd) * (vx_l + 0.0) * 0.25
      end if
    end if
  end do
  
  
  ! set total energy; it needs to be done after because it uses the `i+1` values.
  ! only applies to ZEUS and RAMSES solvers.
  if (self%solver(1:7) /= 'stagger') then
  do k=m(3)%lb,m(3)%ub
  do j=m(2)%lb,m(2)%ub
    do i=m(1)%lb,m(1)%uo
      if (self%solver(1:14) == 'zeus_mhd_patch') then
        rhov2 = (0.25 * (ff(i,j,k,idx%px) + ff(i+1,j,k,idx%px))**2) / d(i,j,k) 
      else if (self%solver(1:6) == 'ramses') then
        rhov2 = ff(i,j,k,idx%px)**2 / d(i,j,k) 
      end if
      et(i,j,k) = e1(i,j,k) + 0.5 * rhov2
    end do
  end do
  end do
  end if
  
  if (idx%bx >0) then
  ff(:,:,:,idx%bx:idx%bz) = 0.0
  end if
  
  ! --------------------------------------------------------------------
  ! Set boundary values
  ! --------------------------------------------------------------------
  pboundary%rho_l = rho_l
  pboundary%e_l   = ff(m(1)%li,m(2)%li,m(3)%li,idx%e)
  pboundary%vx_l  = vx_l
  pboundary%P_l   = p_l
  !
  if (self%solver(1:6) == 'ramses') then
    deallocate (e1)
  else
    nullify (e1)
  end if
  nullify (d, et)
  call trace%end()
END SUBROUTINE condition

END MODULE initial_mod
