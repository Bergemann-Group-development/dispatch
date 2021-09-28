!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_conditions_mod
  !USE const
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  USE index_mod
  USE patch_mod
  USE gpatch_mod
  USE io_unit_mod
  USE math_mod
  USE kinds_mod
  USE eos_mod
  implicit none
  private
  type, public:: initial_conditions_t
    integer:: seed
    real:: d0
    logical:: mhd=.true.
    real(8):: time=0d0
    character(len=64):: solver
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: init
  end type
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !---------------------------------------------------------------------
  call trace_begin('initial_t%init')
  self%patch => patch
  self%mhd = .true.
  call trace%end()
END subroutine pre_init

 
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  ! ............................................................................
  type(index_t), pointer:: idx
  class(mesh_t), pointer    :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  ! ....................................................................
  integer:: i, j, k
  real:: x, y, z, ds2(3), xs, ys, zs, rd, b2
  logical:: first_time = .true.
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: d, s
  real(kind=KindScalarVar), dimension(:,:,:),   pointer:: e1, et
  real(kind=KindScalarVar):: rhov2, pxl
  real(4), save:: P_l=1.0  , P_r=0.1, &
                  rho_l=1.0, rho_r=0.125,&
                  vx_l=0., vx_r=0.0, &
                  vy_l=0., vy_r=0.0, &
                  vz_l=0., vz_r=0.0, &
                  Bx_l=0.75, Bx_r=0.75, &
                  By_l=1., By_r=-1.0, &
                  Bz_l=0., Bz_r=0.0, &
                  e_l=1.0, e_r=0.1, &
                  x_disc = 0.5
  integer, save:: idir=1
  namelist /IC_params/ p_l, p_r, rho_l, rho_r, bx_l, bx_r, by_l, by_r,bz_l, bz_r, &
                           vx_l, vx_r,vy_l, vy_r, vz_l, vz_r,x_disc, idir
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
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !
    d => ff(:,:,:,idx%d)
    if (idx%s>0)   s => ff(:,:,:,idx%s)
    if (idx%e>0)  then
      if (self%patch%kind(1:6) == 'ramses') then
        et => ff(:,:,:,idx%e)
        allocate (e1(m(1)%gn,m(2)%gn,m(3)%gn))
      else
        e1 => ff(:,:,:,idx%e)
      end if
    end if
    if (idx%et>0) et => ff(:,:,:,idx%et)

  
  ff(:,:,:,idx%py) = 0.0
  ff(:,:,:,idx%pz) = 0.0
  if (self%patch%kind(1:7)=='stagger') then
    ds2 = 0.5*m%d
  else if (self%patch%kind(1:6)=='zeus3d') then
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
          e1(i,:,:) = p_l/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * eos%gamma) * rho_l / (eos%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**eos%gamma)  * rho_l / (eos%gamma - 1.0)
          !s (i,:,:) = (log(p_l)-log(rho_l)*eos%gamma)  * rho_l / (eos%gamma - 1.0)
          s (i,:,:) = (eos%s0+(log(P_l)-log(rho_l)*eos%gamma)) * rho_l / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = vx_l * rho_l
        ff(i,:,:,idx%py) = vy_l * rho_l
        ff(i,:,:,idx%pz) = vz_l * rho_l
        ff(i,:,:,idx%bx) = bx_l
        ff(i,:,:,idx%by) = by_l
        ff(i,:,:,idx%bz) = bz_l
      else
        d (i,:,:) = rho_r
        if (idx%e>0) then
          e1(i,:,:) = p_r/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * eos%gamma) * d(i,:,:) / (eos%gamma - 1.0)
          !s (i,:,:) = (log(p_r)-log(rho_r)*eos%gamma)  * rho_r / (eos%gamma - 1.0)
          s (i,:,:) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma)) * rho_r / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = vx_r * rho_r
        ff(i,:,:,idx%px) = vy_r * rho_r
        ff(i,:,:,idx%px) = vz_r * rho_r
        ff(i,:,:,idx%bx) = bx_r 
        ff(i,:,:,idx%by) = by_r
        ff(i,:,:,idx%bz) = bz_r
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          ff(i,:,:,idx%px) = (rho_l + rho_r) * (vx_l + vx_r) * 0.25
        end if
      end if
      ! correction for staggered magnetic field
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          ff(i,:,:,idx%bx) = (bx_l + bx_r) * 0.5
        end if
    end do

  case(-1)
    do i=m(1)%lb,m(1)%ub
      x = m(1)%p + m(1)%r(i)
      xs = x - ds2(1)
      if (x > x_disc) then
        d (i,:,:)= rho_l
        if (idx%e>0) then
          e1(i,:,:) = p_l/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * eos%gamma) * rho_l / (eos%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**eos%gamma)  * rho_l / (eos%gamma - 1.0)
          !s (i,:,:) = (log(p_l)-log(rho_l)*eos%gamma)  * rho_l / (eos%gamma - 1.0)
          s (i,:,:) = (eos%s0+(log(P_l)-log(rho_l)*eos%gamma)) * rho_l / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = vx_l * rho_l
        ff(i,:,:,idx%py) = vy_l * rho_l
        ff(i,:,:,idx%pz) = vz_l * rho_l
        ff(i,:,:,idx%bx) = bx_l
        ff(i,:,:,idx%by) = by_l
        ff(i,:,:,idx%bz) = bz_l
      else
        d (i,:,:) = rho_r
        if (idx%e>0) then
          e1(i,:,:) = p_r/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * eos%gamma) * d(i,:,:) / (eos%gamma - 1.0)
          !s (i,:,:) = (log(p_r)-log(rho_r)*eos%gamma)  * rho_r / (eos%gamma - 1.0)
          s (i,:,:) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma)) * rho_r / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(i,:,:,idx%px) = vx_r * rho_r
        ff(i,:,:,idx%px) = vy_r * rho_r
        ff(i,:,:,idx%px) = vz_r * rho_r
        ff(i,:,:,idx%bx) = bx_r 
        ff(i,:,:,idx%by) = by_r
        ff(i,:,:,idx%bz) = bz_r
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          ff(i,:,:,idx%px) = (rho_l + rho_r) * (vx_l + vx_r) * 0.25
        end if
      end if
      ! correction for staggered magnetic field
        if (x >= x_disc .and. (x - m(1)%d) < x_disc) then
          ff(i,:,:,idx%bx) = (bx_l + bx_r) * 0.5
        end if
    end do

  case(2)
    do i=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(i)
      ys = y - ds2(2)
      if (y < x_disc) then
        d (:,i,:)= rho_l
        if (idx%e>0) then
          e1(:,i,:) = p_l/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l /
          !(self%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**self%gamma)  * rho_l / (self%gamma - 1.0)
         ! s (:,i,:) = (log(p_l)-log(rho_l)*eos%gamma)  * rho_l / (eos%gamma -1.0)
          s (:,i,:) = (eos%s0+(log(P_l)-log(rho_l)*eos%gamma)) * rho_l / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,i,:,idx%px) = vx_l * rho_l
        ff(:,i,:,idx%py) = vy_l * rho_l
        ff(:,i,:,idx%pz) = vz_l * rho_l
        ff(:,i,:,idx%bx) = bx_l
        ff(:,i,:,idx%by) = by_l
        ff(:,i,:,idx%bz) = bz_l
      else
        d (:,i,:) = rho_r
        if (idx%e>0) then
          e1(:,i,:) = p_r/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * self%gamma) * d(i,:,:) /
          !(self%gamma - 1.0)
          !s (:,i,:) = (log(p_r)-log(rho_r)*eos%gamma)  * rho_r / (eos%gamma -1.0)
          s (:,i,:) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma)) * rho_r / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,i,:,idx%px) = vx_r * rho_r
        ff(:,i,:,idx%px) = vy_r * rho_r
        ff(:,i,:,idx%px) = vz_r * rho_r
        ff(:,i,:,idx%bx) = bx_r 
        ff(:,i,:,idx%by) = by_r
        ff(:,i,:,idx%bz) = bz_r
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (y >= x_disc .and. (y - m(2)%d) < x_disc) then
          ff(:,i,:,idx%py) = (rho_l + rho_r) * (vy_l + vy_r) * 0.25
        end if
      end if
      ! correction for staggered magnetic field
        if (y >= x_disc .and. (y - m(2)%d) < x_disc) then
          ff(:,i,:,idx%by) = (by_l + by_r) * 0.5
        end if
    end do

  case(3)
    do i=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(i)
      zs = z - ds2(3)
      if (z < x_disc) then
        d (:,:,i)= rho_l
        if (idx%e>0) then
          e1(:,:,i) = p_l/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_l) - log(rho_l) * self%gamma) * rho_l /
          !(self%gamma - 1.0)
          !s (i,:,:) = log(P_l/rho_l**self%gamma)  * rho_l / (self%gamma - 1.0)
          !s (:,:,i) = (log(p_l)-log(rho_l)*eos%gamma)  * rho_l / (eos%gamma -1.0)
          s (:,:,i) = (eos%s0+(log(P_l)-log(rho_l)*eos%gamma)) * rho_l / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,:,i,idx%px) = vx_l * rho_l
        ff(:,:,i,idx%py) = vy_l * rho_l
        ff(:,:,i,idx%pz) = vz_l * rho_l
        ff(:,:,i,idx%bx) = bx_l
        ff(:,:,i,idx%by) = by_l
        ff(:,:,i,idx%bz) = bz_l
      else
        d (:,:,i) = rho_r
        if (idx%e>0) then
          e1(:,:,i) = p_r/(eos%gamma-1)
        else if (idx%s>0) then
          !s (i,:,:) = (log(P_r) - log(d(i,:,:)) * self%gamma) * d(i,:,:) /
          !(self%gamma - 1.0)
          !s (:,:,i) = (log(p_r)-log(rho_r)*eos%gamma)  * rho_r / (eos%gamma -1.0)
          s (:,:,i) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma)) * rho_r / (eos%gamma - 1.0)
        else 
          error stop "inital_mod: neither energy nor entropy in mem indices"
        end if
        ff(:,:,i,idx%px) = vx_r * rho_r
        ff(:,:,i,idx%px) = vy_r * rho_r
        ff(:,:,i,idx%px) = vz_r * rho_r
        ff(:,:,i,idx%bx) = bx_r 
        ff(:,:,i,idx%by) = by_r
        ff(:,:,i,idx%bz) = bz_r
      end if
          ! 
      ! correction for staggered momentum
      if (any(ds2>0)) then
        if (z >= x_disc .and. (z - m(3)%d) < x_disc) then
          ff(:,:,i,idx%pz) = (rho_l + rho_r) * (vz_l + vz_r) * 0.25
        end if
      end if
      ! correction for staggered magnetic field
        if (z >= x_disc .and. (z - m(3)%d) < x_disc) then
          ff(:,:,i,idx%bz) = (bz_l + bz_r) * 0.5
        end if
    end do


  end select
      
  
  ! set total energy; it needs to be done after because it uses the `i+1` values.
  ! only applies to ZEUS and RAMSES solvers.
  if (associated(et)) then
  do k=m(3)%lb,m(3)%ub-1
  do j=m(2)%lb,m(2)%ub-1
    do i=m(1)%lb,m(1)%ub-1 
      if (self%patch%kind(1:14) == 'zeus_mhd_patch') then
        rhov2 = (0.25 * (ff(i,j,k,idx%px) + ff(i+1,j,k,idx%px))**2) / d(i,j,k) 
      else if (self%patch%kind(1:6) == 'ramses') then
        rhov2 = (ff(i,j,k,idx%px)**2 + ff(i,j,k,idx%py)**2 + ff(i,j,k,idx%pz)**2) / d(i,j,k) 
        if (idx%bx > 0) then
              B2 = 0.25*(ff(i,j,k,idx%bx) + ff(i+1,j,k,idx%bx))**2 + &
                  0.25*(ff(i,j,k,idx%by) + ff(i,j+1,k,idx%by))**2 + &
                  0.25*(ff(i,j,k,idx%bz) + ff(i,j,k+1,idx%bz))**2
        else
            B2 = 0.0
        end if
      end if
      et(i,j,k) = e1(i,j,k) + 0.5 * rhov2 + 0.5*B2
    end do
     et(m(1)%ub,j,k)=et(m(1)%ub-1,j,k)
  end do
  end do
  end if

  
  ! --------------------------------------------------------------------
  ! Set boundary values
  ! -------------------------------------------------------------------
  !
  if (associated(e1)) then
    deallocate (e1)
  else
    nullify (e1)
  end if
  nullify (d, et)
  call trace%end()
END SUBROUTINE init

END MODULE initial_conditions_mod
