!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup 1D dusty waves to test hydrodynamics and gas-dust coupling
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE scalar_mod
  USE timestep_mod
  USE index_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    real:: pi2=8.0*atan(1.0)
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
  real::    x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: pp, pp2
  real, dimension(:,:,:), pointer::   d  , pg , ux , uy , uz  , &
                                      px , py , pz , s  , dtr , &
                                      d2 , pg2, ux2, uy2, uz2 , &
                                      px2, py2, pz2, s2! , dtr2
  real,save:: P  =1.0, rho  =1.0, vx  =0d0, vy  =0d0, vz  =0d0, rhotr  =0d0, &
              dP =0.1, drho =0.1, dvx =0.1, dvy =0d0, dvz =0d0, drhotr =0d0, &
              P2 =0d0, rho2 =1.0, vx2 =0d0, vy2 =0d0, vz2 =0d0, &!rhotr2 =0d0, &
              dP2=0d0, drho2=0.1, dvx2=0.1, dvy2=0d0, dvz2=0d0!, drhotr2=0d0
  namelist /experiment_params/ P  , rho  , vx  , vy  , vz  , rhotr  , &
                               dP , drho , dvx , dvy , dvz , drhotr , &
                               P2 , rho2 , vx2 , vy2 , vz2 , &!rhotr2 , &
                               dP2, drho2, dvx2, dvy2, dvz2!, drhotr2

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  self%nv  = 10
  self%nc  = 2
  self%mhd = .false.
  call self%mhd_t%init

  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if

  call allocate_scalars(self%gn, pg, pg2, ux, uy, uz, ux2, uy2, uz2)

  d   => self%mem(:,:,:,self%idx%d        ,self%it,1)
  s   => self%mem(:,:,:,self%idx%s        ,self%it,1)
  pp  => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  dtr => self%mem(:,:,:,self%idx%dtr      ,self%it,1)
  px  => pp(:,:,:,1); py => pp(:,:,:,2); pz => pp(:,:,:,3);

  self%no_pressure = 2
  d2   => self%mem(:,:,:,self%idx%d+5          ,self%it,1)
  s2   => self%mem(:,:,:,self%idx%s+5          ,self%it,1)
  pp2  => self%mem(:,:,:,self%idx%px+5:self%idx%pz+5,self%it,1)
  !dtr2 => self%mem(:,:,:,idx%dtr+5        ,self%it,1)
  px2  => pp2(:,:,:,1); py2 => pp2(:,:,:,2); pz2 => pp2(:,:,:,3);

  print*, 'gn is:', self%gn(1), self%gn(2), self%gn(3)
  print*, 'P and rho and gamma:'  , P , dP , rho , drho , eos%gamma
  print*, 'P2 and rho2 and gamma:', P2, dP2, rho2!, drho2, eos%gamma
  print*, 'box size', self%box(1)

  do k = self%mesh(3)%lb , self%mesh(3)%ub
    z  = self%mesh(3)%p  + self%mesh(3)%r(k)
    do j = self%mesh(2)%lb , self%mesh(2)%ub
      y  = self%mesh(2)%p  + self%mesh(2)%r(j)
      do i = self%mesh(1)%lb , self%mesh(1)%ub
        x  = self%mesh(1)%p  + self%mesh(1)%r(i)
        !r2 = x**2 + y**2 + z**2

        d   (i,j,k) = rho    + drho    * sin( self%pi2*x/self%box(1) )
        dtr (i,j,k) = rhotr  + drhotr  * sin( self%pi2*x/self%box(1) )
        d2  (i,j,k) = rho2   + drho2   * sin( self%pi2*x/self%box(1) )
        !dtr2(i,j,k) = rhotr2 + drhotr2 * sin( self%pi2*x/self%box(1) )

        pg  (i,j,k) = P      + dP      * sin( self%pi2*x/self%box(1) )
        pg2 (i,j,k) = P2     + dP2     * sin( self%pi2*x/self%box(1) )

        ux  (i,j,k) = vx     + dvx     * sin( self%pi2*x/self%box(1) )
        uy  (i,j,k) = vy     + dvy     * sin( self%pi2*x/self%box(1) )
        uz  (i,j,k) = vz     + dvz     * sin( self%pi2*x/self%box(1) )
        ux2 (i,j,k) = vx2    + dvx2    * sin( self%pi2*x/self%box(1) )
        uy2 (i,j,k) = vy2    + dvy2    * sin( self%pi2*x/self%box(1) )
        uz2 (i,j,k) = vz2    + dvz2    * sin( self%pi2*x/self%box(1) )

        if (eos%gamma .eq. 1.) then
          s (i,j,k) = -log(d(i,j,k))*d(i,j,k)
          s2(i,j,k) = 0d0
        else
          s (i,j,k) = (log(pg(i,j,k))-log(d(i,j,k))*eos%gamma)/(eos%gamma-1.0) * d(i,j,k)
          s2(i,j,k) = 0d0
        end if

        px (i,j,k) = ux (i,j,k) * d (i,j,k)
        py (i,j,k) = uy (i,j,k) * d (i,j,k)
        pz (i,j,k) = uz (i,j,k) * d (i,j,k)
        px2(i,j,k) = ux2(i,j,k) * d2(i,j,k)
        py2(i,j,k) = uy2(i,j,k) * d2(i,j,k)
        pz2(i,j,k) = uz2(i,j,k) * d2(i,j,k)
      end do
    end do
  end do
  call deallocate_scalars(pg, pg2, ux, uy, uz, ux2, uy2, uz2)
  call self%clear(bits%root_grid)
  self%periodic = .true.
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
  call self%mhd_t%update
  !call self%mhd_t%nonideal_MHD_courant_condition(self%dt_Ohm,self%dt_AD)

  call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density('after timestep')

  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
