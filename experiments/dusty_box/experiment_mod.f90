!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup the dusty box problem of a gas fluid colliding with a dust fluid
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
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup a patch with a uniform gas and dust backgrounds
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k, nc!, nvc
  real::    x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: pp, pp2
  real, dimension(:,:,:), pointer::   d  , pg , ux , uy , uz  , &
                                      px , py , pz , s  , dtr , &
                                      d2 , pg2, ux2, uy2, uz2 , &
                                      px2, py2, pz2, s2! , dtr2
  real,save:: P =1.0, rho =1.0, vx = 0.5, vy =0d0, vz =0d0, rhotr =0d0, &
              P2=0d0, rho2=1.0, vx2=-0.5, vy2=0d0, vz2=0d0!, rhotr2=0d0
  namelist /experiment_params/ P , rho , vx , vy , vz , rhotr , &
                               P2, rho2, vx2, vy2, vz2!, rhotr2

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  nc = 2!; nvc = 6;
  self%nv  = 10!nvc*nc                        ! MUST BE NÙMBER OF VARIABLES PER COMPONENT TIMES NUMBER OF COMPONENTS
  self%nc  = nc
  self%mhd = .false.
  call self%mhd_t%init

  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) then
      write (*,experiment_params)
      write (io_unit%nml, experiment_params)
    end if
    first_time = .false.
  end if

  call allocate_scalars(self%gn, pg, pg2, ux, uy, uz, ux2, uy2, uz2)

  d   => self%mem(:,:,:,self%idx%d             ,self%it,1)
  s   => self%mem(:,:,:,self%idx%s             ,self%it,1)
  pp  => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  dtr => self%mem(:,:,:,self%idx%dtr           ,self%it,1)
  px  => pp(:,:,:,1); py => pp(:,:,:,2); pz => pp(:,:,:,3);

  self%no_pressure = 2
  d2   => self%mem(:,:,:,self%idx%d +5              ,self%it,1)
  s2   => self%mem(:,:,:,self%idx%s +5              ,self%it,1)
  pp2  => self%mem(:,:,:,self%idx%px+5:self%idx%pz+5,self%it,1)
  !dtr2 => self%mem(:,:,:,self%idx%dtr+self%nv-5          ,self%it,1)
  px2  => pp2(:,:,:,1); py2 => pp2(:,:,:,2); pz2 => pp2(:,:,:,3);

  print*, 'P and rho and gamma:'  , P , rho , eos%gamma
  print*, 'P2 and rho2 and gamma:', P2, rho2, eos%gamma

  do k = self%mesh(3)%lb , self%mesh(3)%ub
    do j = self%mesh(2)%lb , self%mesh(2)%ub
      do i = self%mesh(1)%lb , self%mesh(1)%ub

        d   (i,j,k) = rho
        dtr (i,j,k) = rhotr
        d2  (i,j,k) = rho2
        !dtr2(i,j,k) = rhotr2

        pg  (i,j,k) = P
        pg2 (i,j,k) = P2

        ux  (i,j,k) = vx
        uy  (i,j,k) = vy
        uz  (i,j,k) = vz
        ux2 (i,j,k) = vx2
        uy2 (i,j,k) = vy2
        uz2 (i,j,k) = vz2

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

  !call self%check_density('before timestep')
  !call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density('after timestep')

  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
