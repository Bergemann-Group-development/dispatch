!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the C-shock problem to test ambipolar diffusion
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  USE mesh_mod
  USE pboundary_mod
  USE index_mod
  USE solver_mod
  USE mhd_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type

  logical, save:: first_time = .true.

CONTAINS

!===============================================================================
!> Set up a patch as part of a 1-D MHD shock tube.
!> Default values are from Ryu & Jones (1995), Fig. 4a.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i, j, k
  real(8):: x, y, z
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, &
                                                        px, py, pz, e1, et, s
  real(4), save:: P_l=1.0  , P_r=0.1, &
                  rho_l=1.0, rho_r=0.2, &
                  vx_l=0.0, vx_r=0.0, vy_l=0.0, vy_r=0.0, vz_l=0.0, vz_r=0.0, &
                  Bx_l=1.0, Bx_r=1.0, By_l=1.0, By_r=0.0, Bz_l=0.0, Bz_r=0.0, &
                  xdisc=0.5
  real(4), parameter:: eps = 1.0e-10
  real(kind=KindScalarVar):: rhov2, B2
  namelist /experiment_params/ P_l, P_r, rho_l, rho_r, &
                               vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, &
                               Bx_l, Bx_r, By_l, By_r, Bz_l, Bz_r, &
                               xdisc
  !----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  ! intialise solver
  self%mhd = .true.
  call self%solver_t%init
  self%periodic = [.false., .true., .true.]

  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)

  ! fool-proofing
  if (Bx_l .ne. Bx_r) then
    call mpi%abort("Shock tubes in the x-direction must have a constant Bx!")
  end if

  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  if (self%kind(1:4) == 'zeus' .or. self%kind(1:10) == 'stagger2e_') then
    e1 => self%mem(:,:,:,self%idx%e ,self%it,1)
  end if
  if (self%kind(1:9) == 'stagger2_') then
    s => self%mem(:,:,:,self%idx%s ,self%it,1)
  end if
  if (self%kind(1:14) == 'zeus_mhd_patch') then
    et => self%mem(:,:,:,self%idx%et,self%it,1)
  end if
  if (self%kind(1:6) == 'ramses') then
    et => self%mem(:,:,:,self%idx%e,self%it,1)
    allocate (e1(self%gn(1),self%gn(2),self%gn(3)))
  end if
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  Bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  By => self%mem(:,:,:,self%idx%by,self%it,1)
  Bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)

  do k=m3%lb,m3%ub
    y = m3%p + m3%r(k)
    do j=m2%lb,m2%ub
      y = m2%p + m2%r(j)
      do i=m1%lb,m1%ub
        x = m1%p + m1%r(i)
        ! for the moment, only initial discontinuities in the x-direction are supported.
        if (x < xdisc) then
          d (i,j,k) = rho_l
          if (self%kind(1:4) == 'zeus' .or. self%kind(1:10) == 'stagger2e_') then
            e1(i,j,k) = P_l / (self%gamma - 1.0)
          else if (self%kind(1:9) == 'stagger2_') then
            s (i,j,k) = (log(P_l) - log(rho_l) * self%gamma) * rho_l / (self%gamma - 1.0)
          else if (self%kind(1:6) == 'ramses')  then
            e1(i,j,k) = P_l / (self%gamma - 1.0)
          end if
          px(i,j,k) = vx_l * rho_l
          py(i,j,k) = vy_l * rho_l
          pz(i,j,k) = vz_l * rho_l
          Bx(i,j,k) = Bx_l
          By(i,j,k) = By_l
          Bz(i,j,k) = Bz_l
        else
          d (i,j,k) = rho_r
          if (self%kind(1:4) == 'zeus' .or. self%kind(1:10) == 'stagger2e_') then
            e1(i,j,k) = P_r / (self%gamma - 1.0)
          else if (self%kind(1:9) == 'stagger2_') then
            s (i,j,k) = (log(P_r) - log(rho_r) * self%gamma) * rho_r / (self%gamma - 1.0)
          else if (self%kind(1:6) == 'ramses')  then
            e1(i,j,k) = P_r / (self%gamma - 1.0)
          end if
          px(i,j,k) = vx_r * rho_r
          py(i,j,k) = vy_r * rho_r
          pz(i,j,k) = vz_r * rho_r
          Bx(i,j,k) = Bx_r
          By(i,j,k) = By_r
          Bz(i,j,k) = Bz_r
        end if
        ! correction for staggered momentum
        if (self%kind(1:4) == 'zeus' .or. self%kind(1:7) == 'stagger') then
          if (x >= xdisc .and. (x - m1%d) < xdisc) then
            px(i,j,k) = (rho_l + rho_r) * (vx_l + vx_r) * 0.25
          end if
        end if
      end do
      ! set total energy; it needs to be done after because it uses the `i+1` values.
      ! only applies to ZEUS3D and RAMSES solvers.
      if (self%kind(1:7) /= 'stagger' .and. self%kind(1:16) /= 'zeus_tw_mhd_patch') then
        do i=m1%lb,m1%uo
          if (self%kind(1:14) == 'zeus_mhd_patch') then
            rhov2 = (0.25 * (px(i,j,k) + px(i+1,j,k))**2 + py(i,j,k)**2 + pz(i,j,k)**2) &
                  / d(i,j,k) 
          else if (self%kind(1:6) == 'ramses') then
            rhov2 = (px(i,j,k)**2 + py(i,j,k)**2 + pz(i,j,k)**2) / d(i,j,k) 
          end if
          B2 = 0.25 * (Bx(i,j,k) + Bx(i+1,j,k))**2 + By(i,j,k)**2 + Bz(i,j,k)**2
          et(i,j,k) = e1(i,j,k) + 0.5 * rhov2 + 0.5 * B2
        end do
      end if
    end do
  end do

  ! assign boundary conditions.
  call pboundary%condition(self%mem(:,:,:,:,self%it,1), self%mesh, .true.)
  call pboundary%condition(self%mem(:,:,:,:,self%it,1), self%mesh, .false.)

  if (self%kind(1:6) == 'ramses') then
    deallocate (e1)
  else
    nullify (e1)
  end if
  nullify (d, et, px, py, pz, Bx, By, Bz)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Experiment update procedure
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if ((self%llc_nat(1)-self%ds(1)) .lt. self%mesh(1)%origin) then
    call pboundary%condition (self%mem(:,:,:,:,self%it,1), self%mesh, .true.)
  endif
  if ((self%llc_nat(1)+self%size(1)+self%ds(1)) .gt. (self%mesh(1)%origin + self%box(1))) then
    call pboundary%condition (self%mem(:,:,:,:,self%it,1), self%mesh, .false.)
  endif
  call self%mhd_t%update
  call trace_end
END SUBROUTINE update

!===============================================================================
END MODULE experiment_mod
