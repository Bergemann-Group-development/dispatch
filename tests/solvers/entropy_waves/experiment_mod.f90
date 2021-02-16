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
  USE index_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
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
  integer:: i, j, k, iostat
  real(8):: x, y, z, g1, g2
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, &
                                                        px, py, pz, e1, et, s
  real(4), save:: P_l=1.0  , P_r=1.0, &
                  rho_l=0.9, rho_r=1.1, &
                  vx_l=0.01, vx_r=0.01, vy_l=0.0, vy_r=0.0, vz_l=0.0, vz_r=0.0, &
                  Bx_l=0.0, Bx_r=0.0, By_l=0.0, By_r=0.0, Bz_l=0.0, Bz_r=0.0, &
                  x_disc=0.0
  real(4), parameter:: eps = 1.0e-10
  real(kind=KindScalarVar):: rhov2, B2
  namelist /experiment_params/ P_l, P_r, rho_l, rho_r, &
                               vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, &
                               Bx_l, Bx_r, By_l, By_r, Bz_l, Bz_r, &
                               x_disc
  !-----------------------------------------------------------------------------
  ! Initialize solver
  !-----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  self%mhd = .true.
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Input namlist and fool-proofing
  !-----------------------------------------------------------------------------
  call io%lock%set ('experiment_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params,iostat=iostat)
    call io%namelist_warning ('experiment_t%init', iostat)
    write (io%namelist,experiment_params)
    first_time = .false.
  end if
  call io%lock%unset ('experiment_t%init')
  if (Bx_l /= Bx_r) then
    call mpi%abort("Shock tubes in the x-direction must have a constant Bx!")
  end if
  !-----------------------------------------------------------------------------
  ! Variables for initial conditions
  !-----------------------------------------------------------------------------
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
  if (self%idx%bx > 0) then
    Bx => self%mem(:,:,:,self%idx%bx,self%it,1)
    By => self%mem(:,:,:,self%idx%by,self%it,1)
    Bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  end if
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  !-----------------------------------------------------------------------------
  ! Initial conditions
  !-----------------------------------------------------------------------------
  g1 = eos%gamma-1d0
  if (g1==0d0) then
    g1 = 1d0
    g2 = 1d0
  else
    g2 = 1d0/g1
  end if
  do k=m3%lb,m3%ub
    y = m3%p + m3%r(k)
    do j=m2%lb,m2%ub
      y = m2%p + m2%r(j)
      do i=m1%lb,m1%ub
        x = m1%p + m1%r(i)
        ! for the moment, only initial discontinuities in the x-direction are supported.
        if (x < x_disc) then
          d (i,j,k) = rho_l
          select case (trim(self%etype))
          case ('thermal')
            e1(i,j,k) = P_l*g2
          case ('total')
            e1(i,j,k) = P_l*g2
          case ('entropy')
            s (i,j,k) = (log(P_l)-log(rho_l)*eos%gamma)*rho_l*g2
          case default
            write (stderr,'(a)') 'unknown self%etype: '//trim(self%etype)
          end select
          px(i,j,k) = vx_l * rho_l
          py(i,j,k) = vy_l * rho_l
          pz(i,j,k) = vz_l * rho_l
          if (self%idx%bx > 0) then
            Bx(i,j,k) = Bx_l
            By(i,j,k) = By_l
            Bz(i,j,k) = Bz_l
          end if
        else
          d (i,j,k) = rho_r
          select case (trim(self%etype))
          case ('thermal')
            e1(i,j,k) = P_r*g2
          case ('total')
            e1(i,j,k) = P_r*g2
          case ('entropy')
            s (i,j,k) = (log(P_r)-log(rho_r)*eos%gamma)*rho_r*g2
          case default
            write (stderr,'(a)') 'unknown self%etype: '//trim(self%etype)
          end select
          px(i,j,k) = vx_r * rho_r
          py(i,j,k) = vy_r * rho_r
          pz(i,j,k) = vz_r * rho_r
          if (self%idx%bx > 0) then
            Bx(i,j,k) = Bx_r
            By(i,j,k) = By_r
            Bz(i,j,k) = Bz_r
          end if
        end if
        ! correction for staggered momentum
        if (self%kind(1:4) == 'zeus' .or. self%kind(1:7) == 'stagger') then
          if (x >= x_disc .and. (x - m1%d) < x_disc) then
            px(i,j,k) = (rho_l + rho_r) * (vx_l + vx_r) * 0.25
          end if
        end if
      end do
      ! set total energy; it needs to be done after because it uses the `i+1` values.
      ! only applies to ZEUS3D and RAMSES solvers.
      if (trim(self%etype) == 'total') then
        do i=m1%lb,m1%uo
          if (self%kind(1:14) == 'zeus_mhd_patch') then
            rhov2 = (0.25*(px(i,j,k) + px(i+1,j,k))**2 + py(i,j,k)**2 + pz(i,j,k)**2) &
                  / d(i,j,k) 
          else
            rhov2 = (px(i,j,k)**2 + py(i,j,k)**2 + pz(i,j,k)**2)/d(i,j,k) 
          end if
          if (self%idx%bx > 0) then
            B2 = 0.25*(Bx(i,j,k) + Bx(i+1,j,k))**2 + By(i,j,k)**2 + Bz(i,j,k)**2
          else
            B2 = 0.0
          end if
          et(i,j,k) = e1(i,j,k) + 0.5*rhov2 + 0.5*B2
        end do
      end if
    end do
  end do
  if (self%kind(1:6) == 'ramses') then
    deallocate (e1)
  else
    nullify (e1)
  end if
  nullify (d, et, px, py, pz, Bx, By, Bz)
  call trace_end ()
END SUBROUTINE init

!===============================================================================
END MODULE experiment_mod
