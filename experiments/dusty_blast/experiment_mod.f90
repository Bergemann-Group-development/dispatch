!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup 1D dusty waves to test hydrodynamics and gas-dust coupling
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE timestep_mod
  USE index_mod
  USE solver_mod
  USE pboundary_mod
  USE scalar_mod
  USE eos_mod
  implicit none
  private
    real:: pi2 = 8.0*atan(1.0)
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup a patch with a periodic sinusoidal wave in the x-direction
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k, nvc, nc
  real::    x, y, z, N
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: pp , pp2
  real, dimension(:,:,:), pointer::   d  , px , py , pz , s , dtr , &
                                      d2 , px2, py2, pz2, s2, &!dtr2, &
                                      r2 , lnPg
  real,save:: P  =0.1, rho  =0.1 , vx  =0d0, vy  =0d0, vz  =0d0, rhotr  =0d0, &
              Pc =1.0, rhoc =1.0 , vxc =0d0, vyc =0d0, vzc =0d0, rhotrc =0d0, &
              P2 =0d0, rho2 =0.01, vx2 =0d0, vy2 =0d0, vz2 =0d0, &!rhotr2 =0d0, &
              P2c=0d0, rho2c=0.01, vx2c=0d0, vy2c=0d0, vz2c=0d0, &!rhotr2c=0d0, &
              r0 =0.125, center(3)=0.0, Etherm=0.0, wgauss=0.0, Pgsum, pgauss=2.0
  namelist /experiment_params/ P  , rho  , vx  , vy  , vz  , rhotr  , &
                               Pc , rhoc , vxc , vyc , vzc , rhotrc , &
                               P2 , rho2 , vx2 , vy2 , vz2 , &!rhotr2 , &
                               P2c, rho2c, vx2c, vy2c, vz2c, &!rhotr2c, &
                               r0 , center, Etherm, wgauss, pgauss

  call trace_begin ('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values.  self%nv must be the number of
  ! variables per component, times the number of components.  io%nv is the
  ! number of variables to write to the output file(s).
  !----------------------------------------------------------------------------
  self%nv  = 10
  io%nv = 9
  self%nc  = 2
  self%mhd = .false.
  call self%solver_t%init

  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time=.false.
  end if

  d   => self%mem(:,:,:,self%idx%d             ,self%it,1)
  s   => self%mem(:,:,:,self%idx%s             ,self%it,1)
  pp  => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  dtr => self%mem(:,:,:,self%idx%dtr           ,self%it,1)
  px  => pp(:,:,:,1); py => pp(:,:,:,2); pz => pp(:,:,:,3);

  self%no_pressure = 2
  d2   => self%mem(:,:,:,self%idx%d+5               ,self%it,1)
  s2   => self%mem(:,:,:,self%idx%s+5               ,self%it,1)
  pp2  => self%mem(:,:,:,self%idx%px+5:self%idx%pz+5,self%it,1)
  !dtr2 => self%mem(:,:,:,self%idx%dtr+5        ,self%it,1)
  px2  => pp2(:,:,:,1); py2 => pp2(:,:,:,2); pz2 => pp2(:,:,:,3);

  call allocate_scalars(self%gn, r2, lnPg)

  do k = self%mesh(3)%lb , self%mesh(3)%ub
    z  = self%mesh(3)%p  + self%mesh(3)%r(k)
    do j = self%mesh(2)%lb , self%mesh(2)%ub
      y  = self%mesh(2)%p  + self%mesh(2)%r(j)
      do i = self%mesh(1)%lb , self%mesh(1)%ub
        x  = self%mesh(1)%p  + self%mesh(1)%r(i)
        r2(i,j,k) = (x - center(1))**2 + (y - center(2))**2 + (z - center(3))**2
        ! set up Gaussian pressure distribution on otherwise uniform background
        if (wgauss>0.0) then
          d   (i,j,k) = rhoc
          px  (i,j,k) = rhoc  * vxc
          py  (i,j,k) = rhoc  * vyc
          pz  (i,j,k) = rhoc  * vzc
          dtr (i,j,k) = rhotrc
          d2  (i,j,k) = rho2c
          px2 (i,j,k) = rho2c * vx2c
          py2 (i,j,k) = rho2c * vy2c
          pz2 (i,j,k) = rho2c * vz2c
          !dtr2(i,j,k) = rhotr2
          lnPg(i,j,k) = log(Pc) + (-0.5*(r2(i,j,k)/wgauss**2)**(pgauss/2.0))
        ! set up central region.
        else if (r2(i,j,k) .le. real(r0,kind=8)**2) then
          d   (i,j,k) = rhoc
          px  (i,j,k) = rhoc  * vxc
          py  (i,j,k) = rhoc  * vyc
          pz  (i,j,k) = rhoc  * vzc
          dtr (i,j,k) = rhotrc
          d2  (i,j,k) = rho2c
          px2 (i,j,k) = rho2c * vx2c
          py2 (i,j,k) = rho2c * vy2c
          pz2 (i,j,k) = rho2c * vz2c
          !dtr2(i,j,k) = rhotr2c
          lnPg(i,j,k) = log(Pc)
        ! set up central region.
        else
          d   (i,j,k) = rho
          px  (i,j,k) = rho  * vx
          py  (i,j,k) = rho  * vy
          pz  (i,j,k) = rho  * vz
          dtr (i,j,k) = rhotr
          d2  (i,j,k) = rho2
          px2 (i,j,k) = rho2 * vx2
          py2 (i,j,k) = rho2 * vy2
          pz2 (i,j,k) = rho2 * vz2
          !dtr2(i,j,k) = rhotr2
          lnPg(i,j,k) = log(P)
        end if
      end do
    end do
  end do
  print*, 'Loop done!'
  print*, 'Min and max r:', minval(sqrt(r2)), maxval(sqrt(r2))
  print*, 'Min and max d:', minval(d), maxval(d)
  print*, 'Min and max ln(Pg):', minval(lnPg), maxval(lnPg)
  ! Calculate entropy
  ! If total internal energy is stated, rescale pressure appropriately
  if (Etherm>0.0) then
    ! Normalize Pg
    N = 1.0/( wgauss**3 * sqrt(self%pi2)**3 )
    lnPg(:,:,:) = lnPg(:,:,:) + log(N) - log(Pc)
    ! Multiply by Etherm to make Pg the energy density profile
    lnPg(:,:,:) = lnPg(:,:,:) + log(Etherm)
    ! Use ideal gas law to turn Pg back into pressure profile
    lnPg(:,:,:) = log(eos%gamma-1.0) + lnPg(:,:,:) !Pg is now again the thermal pressure
    print*, 'Rescaling done!'
  end if
  ! Reset ambient pressure values to the chosen minimum value
  do k = self%mesh(3)%lb , self%mesh(3)%ub
    do j = self%mesh(2)%lb , self%mesh(2)%ub
      do i = self%mesh(1)%lb , self%mesh(1)%ub
        lnPg(i,j,k) = max(lnPg(i,j,k),log(P))
      end do
    end do
  end do
  print*, 'Min and max ln(Pg):', minval(lnPg), maxval(lnPg)
  ! Calculate entropy
  if (eos%gamma .eq. 1.) then
    s  = -log(d)*d
    s2 = 0d0
  else
    s  = (lnPg-log(d)*eos%gamma)/(eos%gamma-1.0) * d
    s2 = 0d0
  end if

  call deallocate_scalars(lnPg)
  call self%clear(bits%root_grid)
  self%periodic = .false.
  call pboundaries%conditions(self, self%it)
  print*, 'done!'
  call trace_end
END SUBROUTINE init

!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call pboundaries%conditions(self, self%it)
  call self%solver_t%update()
  call pboundaries%conditions(self, self%new)

  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
