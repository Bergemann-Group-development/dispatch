!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Shock-tube problem to test hydrodynamics
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE boundary_mod
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
!> Setup a patch with a sharp shock front as present in a C-shock
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k, nvc, nc                                              ! ADDED INTEGER FOR NUMBER OF VARIABLES PER COMPONENT AND NUMBER OF COMPONENTS
  real:: x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: p, p2                               !
  real, dimension(:,:,:), pointer:: d, vx, vy, vz, px, py, pz, s, dtr, &  ! ADDED TRACER DENSITY
                                    d2, vx2, vy2, vz2, px2, py2, pz2, s2!, dtr2 ! ADDED SECOND COMPONENT
  real,save:: P_l, P_r, rho_i, rho_l, vx_l, vy_l, vz_l=0., &
              rho_r, vx_r, vy_r, vz_r=0., rhotr_l, rhotr_r, &             ! ADDED TRACER DENSITY
              P_l2, P_r2, rho_i2, rho_l2, vx_l2, vy_l2, vz_l2=0., &       ! ADDED SECOND COMPONENT
              rho_r2, vx_r2, vy_r2, vz_r2=0.!, rhotr_l2, rhotr_r2
  namelist /experiment_params/ P_l, P_r, rho_l, vx_l, vy_l, rho_r, &
                               vx_r, vy_r, rhotr_l, rhotr_r, &            ! ADDED TRACER DENSITY TO NAMELIST
                               P_l2, P_r2, rho_l2, vx_l2, vy_l2, rho_r2, &! ADDED SECOND COMPONENTS TO NAMELIST
                               vx_r2, vy_r2!, rhotr_l2, rhotr_r2

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  self%nv  = 10                                                        ! MUST BE NÙMBER OF VARIABLES PER COMPONENT TIMES NUMBER OF COMPONENTS
  self%nc  = 2                                                         ! ADDED NUMBER OF COMPONENTS
  self%mhd = .false.
  call self%mhd_t%init

  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if

  d   => self%mem(:,:,:,self%idx%d        ,self%it,1)
  s   => self%mem(:,:,:,self%idx%s        ,self%it,1)
  p   => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  dtr => self%mem(:,:,:,self%idx%dtr      ,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);                       ! ADDED TRACER DENSITY

  self%no_pressure = 2
  d2   => self%mem(:,:,:,self%idx%d +5         ,self%it,1)                 ! ADDED SECOND COMPONENT VARIABLES
  s2   => self%mem(:,:,:,self%idx%s +5         ,self%it,1)
  p2   => self%mem(:,:,:,self%idx%px+5:self%idx%pz+5,self%it,1)
  !dtr2 => self%mem(:,:,:,self%idx%dtr+nvc         ,self%it,1)
  px2=>p2(:,:,:,1); py2=>p2(:,:,:,2); pz2=>p2(:,:,:,3);

  print*, 'gn is:', self%gn(1)
  print*, 'P and rho and gamma:', P_l, P_r, rho_l, rho_r, eos%gamma
  print*, 'P2 and rho2 and gamma:', P_l2, P_r2, rho_l2, rho_r2, eos%gamma

  do i=self%mesh(1)%lb,self%mesh(1)%ub
    x = self%mesh(1)%p + self%mesh(1)%r(i)
    if (x < self%box(1)/2.0) then
      d   (i,:,:) = rho_l
      dtr (i,:,:) = rhotr_l                                              ! ADDED TRACER DENSITY

      d2  (i,:,:) = rho_l2                                               ! DEFINES SECOND COMPONENT DENSITY AND TRACER DENSITY
      !dtr2(i,:,:) = rhotr_l2
      if (eos%gamma .eq. 1.) then
        s (i,:,:) = -log(rho_l)*rho_l
        s2(i,:,:) = 0d0
      else
        s (i,:,:) = (log(P_l)-log(rho_l)*eos%gamma)/(eos%gamma-1.0) * rho_l
        s2(i,:,:) = 0d0
      end if
      px(i,:,:) = vx_l*rho_l
      py(i,:,:) = vy_l*rho_l
      pz(i,:,:) = 0d0*rho_l

      px2(i,:,:) = vx_l2*rho_l2                                          ! DEFINES SECOND COMPONENT MOMENTUM
      py2(i,:,:) = vy_l2*rho_l2
      pz2(i,:,:) = 0d0*rho_l2
    else
      d  (i,:,:) = rho_r
      dtr(i,:,:) = rhotr_r

      d2  (i,:,:) = rho_r2                                               ! DEFINES SECOND COMPONENT DENSITY AND TRACER DENSITY
      !dtr2(i,:,:) = rhotr_r2                                             ! ADDED TRACER DENSITY
      if (eos%gamma .eq. 1.) then
        s (i,:,:) = -log(rho_r)*rho_r
        s2(i,:,:) = 0d0
      else
        s (i,:,:) = (log(P_r)-log(rho_r)*eos%gamma)/(eos%gamma-1.0) * rho_r
        s2(i,:,:) = 0d0
      end if
      px(i,:,:) = vx_r*rho_r
      py(i,:,:) = vy_r*rho_r
      pz(i,:,:) = 0d0*rho_r

      px2(i,:,:) = vx_r2*rho_r2                                          ! DEFINES SECOND COMPONENT MOMENTUM
      py2(i,:,:) = vy_r2*rho_r
      pz2(i,:,:) = 0d0*rho_r2
    endif
  end do

  call self%clear(bits%root_grid)
  self%periodic = .false.
  print*, 'done!'

  ! assign constant boundary conditions.
  boundary%P_l      = P_l
  boundary%rho_l    = rho_l
  boundary%P_l2     = P_l2
  boundary%rho_l2   = rho_l2
  if (eos%gamma .eq. 1.) then
    boundary%s_l    = -log(rho_l)*rho_l
    boundary%s_l2   = 0d0
  else
    boundary%s_l    = (log(P_l)-log(rho_l)*eos%gamma)/(eos%gamma-1.0) * rho_l
    boundary%s_l2   = 0d0
  end if
  boundary%vx_l     = vx_l
  boundary%vy_l     = vy_l
  boundary%vz_l     = vz_l
  boundary%vx_l2    = vx_l2
  boundary%vy_l2    = vy_l2
  boundary%vz_l2    = vz_l2
  boundary%rhotr_l  = rhotr_l                                              ! ADDED TRACER DENSITY
  !boundary%rhotr_l2 = rhotr_l2

  boundary%P_r      = P_r
  boundary%rho_r    = rho_r
  boundary%P_r2     = P_r2
  boundary%rho_r2   = rho_r2
  if (eos%gamma .eq. 1.) then
    boundary%s_r    = -log(rho_r)*rho_r
    boundary%s_r2   = 0d0
  else
    boundary%s_r    = (log(P_r)-log(rho_r)*eos%gamma)/(eos%gamma-1.0) * rho_r
    boundary%s_r2   = 0d0
  end if
  boundary%vx_r     = vx_r
  boundary%vy_r     = vy_r
  boundary%vz_r     = vz_r
  boundary%vx_r2    = vx_r2
  boundary%vy_r2    = vy_r2
  boundary%vz_r2    = vz_r2
  boundary%rhotr_r  = rhotr_r                                              ! ADDED TRACER DENSITY
  !boundary%rhotr_r2 = rhotr_r2

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

  if (self%mesh(1)%p .lt. self%size(1)) then
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.)
    end do
  endif
  if (self%mesh(1)%p .gt. self%box(1) - self%size(1)) then
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.)
    end do
  endif

  call self%mhd_t%update

  if (self%mesh(1)%p .lt. self%size(1)) then
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%new,1), self%mesh, .true.)
    end do
  endif
  if (self%mesh(1)%p .gt. self%box(1) - self%size(1)) then
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%new,1), self%mesh, .false.)
    end do
  endif

  !call self%mhd_t%nonideal_MHD_courant_condition(self%dt_Ohm,self%dt_AD)
  !call self%check_density('before timestep')

  !call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  !call self%check_density('after timestep')
  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
