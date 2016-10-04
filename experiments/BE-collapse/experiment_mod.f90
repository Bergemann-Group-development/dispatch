!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mhd_mod
  USE boundary_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    real:: w_BE
    real:: p_BE
    real:: a_BE
    real:: r_BE
    real:: amplitude
    real:: phi
    real:: i_phi
    real:: i_s
    real:: i_px
    real:: i_py
    real:: i_pz
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  integer, save:: levelmax=0
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_BE=5.0, p_BE=1.0, a_BE=0.10, r_BE=0.45, amplitude=10., r
  namelist /experiment_params/ w_BE, p_BE, a_BE, r_BE, amplitude, levelmax
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  print*, 'start init'
  call io%init
  call trace_begin ('experiment_t%init')
  self%nv = 8
  call self%mhd_t%init

  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)

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
  end if

  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Approximate a Bonnor-Ebert sphere
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,self%gn(3)
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      y = (j-(self%gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      z = (k-(self%gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      r = sqrt(x**2+y**2+z**2)/self%w_BE
      r = min(r,r_BE/w_BE)
      self%d(i,j,k)  = self%amplitude*exp(-r**self%p_BE/(1.+r**self%p_BE/log(r**2/self%a_BE)))
      self%p(i,j,k,1) = 0. ; self%p(i,j,k,2) = 0. ; self%p(i,j,k,3) = 0.
      self%B(i,j,k,1) = 0. ; self%B(i,j,k,2) = 0. ; self%B(i,j,k,3) = 0.      
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Initial condition for phi
  !-----------------------------------------------------------------------------
  self%phi = 0.0
  !
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Enforce constant values outside a radius r=r_BE (not equal to the BE-radius)
!> =============================================================================
SUBROUTINE update(self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  integer:: i, j, k, ii(3), iv
  real:: x, y, z, r, ntot, nout
  !.............................................................................
  call trace_begin('experiment_t%update')
associate (li=>self%li,ui=>self%ui,gn=>self%gn,ds=>self%ds,pos=>self%position)
  !-----------------------------------------------------------------------------
  ! Set phi=0 at and outside the boundary
  !-----------------------------------------------------------------------------
!  if (self%it .ne. 1) then
!    self%mem(:,:,:,:,self%it,:) = self%mem(:,:,:,:,self%it-1,:)
!    self%d => self%mem(:,:,:,1,self%it,1)
!  endif
  !-----------------------------------------------------------------------------
  ! This is where d(f)/dt is computed, so zap to zero outside the boundary,
  ! treating staggered components separately
  !-----------------------------------------------------------------------------
!  print*,'before mhd self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  call self%mhd_t%update
!  print*,'after mhd self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  !call self%mhd_t%nonideal_MHD_courant_condition(self%dt_Ohm,self%dt_AD)
  call self%patch_t%courant_condition()
  call self%check_density('before timestep')
  do k=1,gn(3)
  do j=1,gn(2)
  do i=1,gn(1)
    x = (i-(gn(1)+1)*0.5)*ds(1) + pos(1)
    y = (j-(gn(2)+1)*0.5)*ds(2) + pos(2)
    z = (k-(gn(3)+1)*0.5)*ds(3) + pos(3)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%d(i,j,k) = 2e-11
    end if
  end do
  end do
  end do
!  print*,'before timestep self%dt and minval(self%mem(:,:,:,1,self%it,1))',self%dt, minval(self%mem(:,:,:,1,self%it,1))
  call timestep%update(self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density('after timestep')

end associate
  !
  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE experiment_mod
