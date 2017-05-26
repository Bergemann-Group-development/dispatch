!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the C-shock problem to test ambipolar diffusion
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE scalar_mod
  USE boundary_mod
  USE timestep_mod
  USE mhd_mod
  USE initial_mod
  USE boundary_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    type(boundary_t) :: boundary
    real:: cs, gamma_exp, M, M_A, rho_l, vx_l, vy_l, vz_l, Bx_l, By_l, Bz_l, P_l
    real:: rho_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r, P_r
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Setup a patch with a sharp shock front as present in a C-shock
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  integer, save:: levelmax=0
  real:: x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: b, p
  real, dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, px, py, pz, s
  real,save:: cs, gamma_exp, P_l, P_r, rho_i, rho_l, vx_l, vy_l, vz_l, &
   Bx_l, By_l, Bz_l=0., & 
    rho_r, vx_r, vy_r, vz_r, Bx_r, By_r, Bz_r=0. 
  namelist /experiment_params/ gamma_exp, P_l, P_r, rho_l, vx_l, vy_l, Bx_l, By_l, rho_r, vx_r, vy_r, Bx_r, By_r
  !namelist /MHD_params/ gamma
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call io%init
   call trace_begin ('experiment_t%init')
  self%nv = 8
  call self%mhd_t%init
 
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)

    self%d => self%mem(:,:,:,  1,self%it,1)
    self%s => self%mem(:,:,:,  2,self%it,1)
    self%p => self%mem(:,:,:,3:5,self%it,1)
    self%B => self%mem(:,:,:,6:8,self%it,1)
    d=>self%d; s=>self%s
    p=>self%p; B=>self%B
    px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    print*, 'gn is:', self%gn(1)
    print*, 'P and rho and gamma:', P_l, P_r, rho_l, rho_r, gamma_exp
    do i=1,self%gn(1)
      !x = (i-(self%gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      !y = (j-(self%gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      !z = (k-(self%gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      if (i < self%gn(1)/2+1) then
        d (i,:,:)  = rho_l
        if (gamma_exp .eq. 1.) then
          s (i,:,:) = -log(rho_l)
        else
          s (i,:,:) = (log(P_l)-log(rho_l)*gamma_exp)/(gamma_exp-1.0)
          !s (i,:,:) = -(1.381e-16)/(gamma-1)*log(rho_l)
        end if
        px(i,:,:) = vx_l*rho_l
        py(i,:,:) = vy_l*rho_l
        pz(i,:,:) = 0d0*rho_l
        Bx(i,:,:) = Bx_l
        By(i,:,:) = By_l
        Bz(i,:,:) = Bz_l
      else
        d (i,:,:) = rho_r
        if (gamma_exp .eq. 1.) then
          s (i,:,:) = -log(rho_r)
        else
          s (i,:,:) = (log(P_r)-log(rho_r)*gamma_exp)/(gamma_exp-1.0)
          !s (i,:,:) = -(1.381e-16)/(gamma-1)*log(rho_r)
        end if
        px(i,:,:) = vx_r*rho_r 
        py(i,:,:) = vy_r*rho_r
        pz(i,:,:) = 0d0*rho_r 
        Bx(i,:,:) = Bx_r 
        By(i,:,:) = By_r
        Bz(i,:,:) = Bz_r 
      endif
    end do
  end if
  call boundary%condition (self%mem(:,:,:,1,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,2,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,6,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,7,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,8,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  first_time = .false.
  self%mhd = .true.
  print*, 'done!'
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

 !!  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  !!!FIXME: Estimate adequate timestep for AD by taking the minimum of Courant number
  !!!and t_AD = 0.1*min(gamma_AD*rho_i*rho/B^2*Delta_x^2)!!!
  call self%mhd_t%nonideal_MHD_courant_condition(self%dt_Ohm,self%dt_AD)
  call self%check_density('before timestep')
  call boundary%condition (self%mem(:,:,:,1,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,2,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,6,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,7,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call boundary%condition (self%mem(:,:,:,8,self%it,1), self%mesh, [0d0,-0.5d0,-0.5d0])
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density('after timestep')
  call trace_end
END SUBROUTINE update
!!!> =============================================================================
END MODULE experiment_mod
