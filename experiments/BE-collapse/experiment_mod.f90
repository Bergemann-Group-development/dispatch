!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE gmhd_mod
  USE boundary_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(gmhd_t):: experiment_t
    real:: w_BE
    real:: p_BE
    real:: a_BE
    real:: r_BE
    real:: amplitude
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  integer, save:: levelmax=0
  real:: x, y, z
  logical:: first_time = .true.
  real, save:: w_BE=5.0, p_BE=1.0, a_BE=0.10, r_BE=0.45, amplitude=10., r
  namelist /experiment_params/ w_BE, p_BE, a_BE, r_BE, amplitude, levelmax
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call io%init
  call trace_begin ('experiment_t%init')
  self%nv = 8
  call self%solver_t%init

  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
  end if

  self%d => self%mem(:,:,:,1,self%it,1)
  self%s => self%mem(:,:,:,  2,self%it,1)
  self%p => self%mem(:,:,:,3:5,self%it,1)
  self%B => self%mem(:,:,:,6:8,self%it,1)
  !  d=>self%d; s=>self%s
  !  p=>self%p; B=>self%B
  !  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  !  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);

  self%w_BE      = w_BE
  self%p_BE      = p_BE
  self%a_BE      = a_BE
  self%r_BE      = r_BE
  !self%mhd       = .false.
  boundary%r_BE  = r_BE
  self%amplitude = amplitude

  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Approximate a Bonnor-Ebert sphere
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=self%mesh(3)%lb,self%mesh(3)%ub
    do j=self%mesh(2)%lb,self%mesh(2)%ub
    do i=self%mesh(1)%lb,self%mesh(1)%ub
      x = self%mesh(1)%r(i) + (self%position(1) - self%box(1)*0.5d0)
      y = self%mesh(2)%r(j) + (self%position(2) - self%box(2)*0.5d0)
      z = self%mesh(3)%r(k) + (self%position(3) - self%box(3)*0.5d0)
      r = sqrt(x**2+y**2+z**2)/self%w_BE
      r = min(r,r_BE/w_BE)
      self%d(i,j,k)  = self%amplitude*exp(-r**self%p_BE/(1.+r**self%p_BE/log(r**2/self%a_BE)))
      self%p(i,j,k,1) = 0. ; self%p(i,j,k,2) = 0. ; self%p(i,j,k,3) = 0.
      self%B(i,j,k,1) = 0. ; self%B(i,j,k,2) = 0. ; self%B(i,j,k,3) = 0.
    end do
    end do
    end do
  end if
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Enforce constant values outside a radius r=r_BE (not equal to the BE-radius)
!===============================================================================
SUBROUTINE update(self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  integer:: i, j, k, ii(3), iv
  real:: x, y, z, r, ntot, nout
  !.............................................................................
  call trace_begin('experiment_t%update')
associate (li=>self%li,ui=>self%ui,gn=>self%gn,ds=>self%ds,pos=>self%position)
!  if (self%it .ne. 1) then
!    self%mem(:,:,:,:,self%it,:) = self%mem(:,:,:,:,self%it-1,:)
!    self%d => self%mem(:,:,:,1,self%it,1)
!  endif

!  print*,'before mhd self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  call self%mhd_t%update
!  print*,'after mhd self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  !call self%mhd_t%nonideal_MHD_courant_condition(self%dt_Ohm,self%dt_AD)
  call self%check_density('before timestep')
  self%d => self%mem(:,:,:,  1,self%it,1)
  self%s => self%mem(:,:,:,  2,self%it,1)
  self%p => self%mem(:,:,:,3:5,self%it,1)
  self%B => self%mem(:,:,:,6:8,self%it,1)
  do k=1,gn(3)
  do j=1,gn(2)
  do i=1,gn(1)
    x = self%mesh(1)%r(i) + (self%position(1) - self%box(1)*0.5d0)
    y = self%mesh(2)%r(j) + (self%position(2) - self%box(2)*0.5d0)
    z = self%mesh(3)%r(k) + (self%position(3) - self%box(3)*0.5d0)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%d(i,j,k) =   max(self%d(i,j,k),1e-7)
      self%s(i,j,k) =   0.
      self%p(i,j,k,:) = 0.
      self%B(i,j,k,:) = 0.
    end if
  end do
  end do
  end do
  call self%patch_t%courant_condition()
!  print*,'before timestep self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  call self%solver_t%update
  call self%check_density('after timestep')

end associate
  !
  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
