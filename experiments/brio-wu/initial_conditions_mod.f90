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
  USE gpatch_mod
  USE kinds_mod
  USE eos_mod
  implicit none
  private
  type, public:: initial_conditions_t
    integer:: seed
    integer:: idir
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
  !.............................................................................
  type(index_t), pointer:: idx
  class(mesh_t), pointer:: m
  integer:: i, j(3), iostat, iv
  real:: x, ds2, xs, rd, b2
  logical:: first_time = .true.
  real(var_k), pointer:: ff(:,:,:,:), f(:,:,:)
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
  !-----------------------------------------------------------------------------
  ! Read experiment parameters
  !-----------------------------------------------------------------------------
  call io%lock%set('initial_t%condition')
  if (first_time) then
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, IC_params)
    end if
    write (io%namelist,IC_params)
    first_time = .false.
  end if
  self%idir=idir
  call io%lock%unset('initial_t%condition')
  !-----------------------------------------------------------------------------
  ! Shorthand and staggering choices
  !-----------------------------------------------------------------------------
  if (io%do_generic) then
    call self%patch%allocate_generic()
    ff => self%patch%generic
  else
    ff => self%patch%mem(:,:,:,:,1,1)
  end if
  idx => self%patch%idx
  if (self%patch%kind(1:7)=='stagger' .or. &
      self%patch%kind(1:6)=='zeus3d') then
    ds2 = 0.5*self%patch%mesh(idir)%d
  else
    ds2 = 0.0
  endif
  !-----------------------------------------------------------------------------
  ! Initial conditions for the three possible orientations
  !-----------------------------------------------------------------------------
  m => self%patch%mesh(idir)
  do i=m%lb,m%ub
    select case(idir)
    case(1)
      j = [0,1,2]
      f => ff(i,:,:,:)
      call values
    case(2)
      j = [1,2,0]
      f => ff(:,i,:,:)
      call values
    case(3)
      j = [2,0,1]
      f => ff(:,:,i,:)
      call values
    end select
  end do
  call trace%end()
contains
  !=============================================================================
  !> Left and right values of the generic variables
  !=============================================================================
  subroutine values
    x = m%p + m%r(i)
    xs = x - ds2
    if (x < x_disc) then
      f(:,:,idx%d) = rho_l
      if (io%do_generic) then
        if (idx%e>0) f(:,:,idx%e) = p_l/rho_l
        if (idx%s>0) f(:,:,idx%s) = p_l/rho_l
      else
        if (idx%e>0) f(:,:,idx%e) = p_l/(merge(eos%gamma-1d0,1d0,eos%gamma>1))
        if (idx%s>0) f(:,:,idx%s) = (eos%s0+(log(P_l)-log(rho_l)*eos%gamma))  &
                                  * rho_l / (merge(eos%gamma-1d0,1d0,eos%gamma>1))
      end if
      f(:,:,idx%px+j(1)) = vx_l * rho_l
      f(:,:,idx%px+j(2)) = vy_l * rho_l
      f(:,:,idx%px+j(3)) = vz_l * rho_l
      f(:,:,idx%bx+j(1)) = bx_l
      f(:,:,idx%bx+j(2)) = by_l
      f(:,:,idx%bx+j(3)) = bz_l
    else
      f(:,:,idx%d) = rho_r
      if (idx%e>0) f(:,:,idx%e) = p_r/rho_r
      if (idx%s>0) f(:,:,idx%s) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma))  &
                                * rho_r / (merge(eos%gamma-1d0,1d0,eos%gamma>1))
      if (io%do_generic) then
        if (idx%e>0) f(:,:,idx%e) = p_l/rho_l
        if (idx%s>0) f(:,:,idx%s) = p_l/rho_l
      else
        if (idx%e>0) f(:,:,idx%e) = p_r/(merge(eos%gamma-1d0,1d0,eos%gamma>1))
      if (idx%s>0) f(:,:,idx%s) = (eos%s0+(log(P_r)-log(rho_r)*eos%gamma))  &
                                * rho_r / (merge(eos%gamma-1d0,1d0,eos%gamma>1))
      end if
      f(:,:,idx%px+j(1)) = vx_r * rho_r
      f(:,:,idx%px+j(2)) = vy_r * rho_r
      f(:,:,idx%px+j(3)) = vz_r * rho_r
      f(:,:,idx%bx+j(1)) = bx_r
      f(:,:,idx%bx+j(2)) = by_r
      f(:,:,idx%bx+j(3)) = bz_r
    end if
    !--- correction for staggered momentum and magnetic field ---
    if (x >= x_disc .and. (x - m%d) < x_disc) then
      f(:,:,idx%bx+j(1)) = (bx_l + bx_r) * 0.5
      if (ds2>0) then
        f(:,:,idx%px+j(1)) = (rho_l + rho_r) * (vx_l + vx_r) * 0.25
      end if
    end if
  end subroutine values
END SUBROUTINE init

END MODULE initial_conditions_mod
