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
  real(8), optional:: gamma
  !---------------------------------------------------------------------
  call trace_begin('initial_t%init')
  self%solver = trim(solver)
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
  real:: x, y, z, ds2(3), xs, ys, zs, rd
  logical:: first_time = .true.
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: d, s
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: e1, et
  real(kind=KindScalarVar):: rhov2, pxl
  real(4), save:: P_l=10.33333  , P_r=1.0, &
                  rho_l=3.857143, &
                  vx_l=2.629369, vx_r=0.0, &
                  a_d=0.2, f_d=5.0, x_disc = -4.0, bx=0.0, gamma=1.4
  integer, save:: idir=1
  namelist /IC_params/ P_l, P_r, rho_l, &
                           vx_l, vx_r, a_d, f_d, x_disc, idir, bx, gamma
  call io%lock%set('initial_t%condition')
  if (first_time) then
    rewind (io%input)
    read (io%input, IC_params)
    if (io%master) write (*,IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_t%condition')
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
    if (idx%e>0) then
      et => ff(:,:,:,idx%e)
      allocate (e1(m(1)%gn,m(2)%gn,m(3)%gn))
    else
      s => ff(:,:,:,idx%s)
    end if
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
  select case(idir)
  case(1)
    do i=m(1)%lb,m(1)%ub
      x = m(1)%p + m(1)%r(i)
      xs = x - ds2(1)
      if (x < x_disc) then
        d (i,:,:)= rho_l
        if (idx%e>0) then
          e1(i,:,:) = P_l / (self%gamma - 1.0)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**self%gamma)  * rho_l / (self%gamma - 1.0)
          s (i,:,:) = (log(P_l)-log(rho_l)*self%gamma)  * rho_l / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = vx_l * rho_l
      else
        d (i,:,:) = 1.0*(1.+a_d*sin(f_d*x))
        if (idx%e>0) then
          e1(i,:,:) = P_r / (self%gamma - 1.0)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * self%gamma) * d(i,:,:) / (self%gamma - 1.0)
          s (i,:,:) = (log(P_r)-log(d(i,:,:))*self%gamma) * d(i,:,:) / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = 0.0
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          rd = 1.0*(1.+a_d*sin(f_d*x))
          ff(i,:,:,idx%px) = (rho_l + rd) * (vx_l + 0.0) * 0.25
        end if
      end if
    end do

  case(-1)
    do i=m(1)%lb,m(1)%ub
      x = m(1)%p + m(1)%r(i)
      xs = x - ds2(1)
      if (x > x_disc) then
        d (i,:,:)= rho_l
        if (idx%e>0) then
          e1(i,:,:) = P_l / (self%gamma - 1.0)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**self%gamma)  * rho_l / (self%gamma - 1.0)
          s (i,:,:) = (log(P_l)-log(rho_l)*self%gamma)  * rho_l / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = -vx_l * rho_l
      else
        d (i,:,:) = 1.0*(1.-a_d*sin(f_d*x))
        if (idx%e>0) then
          e1(i,:,:) = P_r / (self%gamma - 1.0)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * self%gamma) * d(i,:,:) / (self%gamma - 1.0)
          s (i,:,:) = (log(P_r)-log(d(i,:,:))*self%gamma) * d(i,:,:) / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = 0.0
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          rd = 1.0*(1.-a_d*sin(f_d*(x - m(1)%d)))
          ff(i,:,:,idx%px) = (rho_l + rd) * (-vx_l + 0.0) * 0.25
        end if
      end if
    end do

  case(2)
    do i=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(i)
      ys = y - ds2(2)
      if (y < x_disc) then
        d (:,i,:)= rho_l
        if (idx%e>0) then
          e1(:,i,:) = P_l / (self%gamma - 1.0)
        else if (idx%s>0) then
          s (:,i,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,i,:,idx%py) = vx_l * rho_l
      else
        d (:,i,:) = 1.0*(1.+a_d*sin(f_d*y))
        if (idx%e>0) then
          e1(:,i,:) = P_r / (self%gamma - 1.0)
        else if (idx%s>0) then
          s (:,i,:) = (log(P_r) - log(d(:,i,:)) * self%gamma) * d(:,i,:) / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,i,:,idx%py) = 0.0
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (y >= x_disc .and. (y - m(2)%d) < x_disc) then
          rd = 1.0*(1.+a_d*sin(f_d*y))
          ff(:,i,:,idx%py) = (rho_l + rd) * (vx_l + 0.0) * 0.25
        end if
      end if
    end do

  case(3)
    do i=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(i)
      zs = z - ds2(3)
      if (z < x_disc) then
        d (:,:,i)= rho_l
        if (idx%e>0) then
          e1(:,:,i) = P_l / (self%gamma - 1.0)
        else if (idx%s>0) then
          s (:,:,i) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,:,i,idx%pz) = vx_l * rho_l
      else
        d (:,:,i) = 1.0*(1.+a_d*sin(f_d*z))
        if (idx%e>0) then
          e1(:,:,i) = P_r / (self%gamma - 1.0)
        else if (idx%s>0) then
          s (:,:,i) = (log(P_r) - log(d(:,:,i)) * self%gamma) * d(:,:,i) / (self%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,:,i,idx%pz) = 0.0
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (z >= x_disc .and. (z - m(3)%d) < x_disc) then
          rd = 1.0*(1.+a_d*sin(f_d*z))
          ff(:,:,i,idx%pz) = (rho_l + rd) * (vx_l + 0.0) * 0.25
        end if
      end if
    end do
  end select
      
  
  ! set total energy; it needs to be done after because it uses the `i+1` values.
  ! only applies to ZEUS and RAMSES solvers.
  if (associated(et)) then
  do k=m(3)%lb,m(3)%ub
  do j=m(2)%lb,m(2)%ub
    do i=m(1)%lb,m(1)%ub-1
      if (self%solver(1:14) == 'zeus_mhd_patch') then
        rhov2 = (0.25 * (ff(i,j,k,idx%px) + ff(i+1,j,k,idx%px))**2) / d(i,j,k) 
      else if (self%solver(1:6) == 'ramses') then
        rhov2 = (ff(i,j,k,idx%px)**2 + ff(i,j,k,idx%py)**2 + ff(i,j,k,idx%pz)**2) / d(i,j,k) 
      end if
      et(i,j,k) = e1(i,j,k) + 0.5 * rhov2
    end do
     et(m(1)%ub,j,k)=et(m(1)%ub-1,j,k)
  end do
  end do
  end if


  
  if (idx%bx >0) then
  ff(:,:,:,idx%bx:idx%bz) = 0.0
  ff(:,:,:,idx%bx) = bx
  end if
  
  ! --------------------------------------------------------------------
  ! Set boundary values
  ! -------------------------------------------------------------------
  !
  if (associated(et)) then
    deallocate (e1)
  else
    nullify (e1)
  end if
  nullify (d, et)
  call trace%end()
END SUBROUTINE condition

END MODULE initial_mod
