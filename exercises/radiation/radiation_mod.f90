!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Radiative Transfer problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE scaling_mod
  USE microphysics_mod
  USE rt_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_tr=2
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, tt
    real:: w_gauss
    real:: p_gauss
    real:: amplitude
    real(kind=8):: tolerance
    type(rt_t)  :: rt
    type(microphysics_t) :: microphysics
  contains
    procedure:: init
    procedure:: solve
  end type
  character(len=12):: solver='integral'
CONTAINS

!> =============================================================================
!> Setup a patch with density and source function for the radiative transfer solver
!> =============================================================================
SUBROUTINE init (self)
  class(pde_t):: self
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_gauss=5.0, p_gauss=1.0, amplitude=1.0
  namelist /rt_params/ w_gauss, p_gauss, amplitude, solver
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  experimentname = 'rt'
  call io%init (experimentname)
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    call scaling%init
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_params)
    if (io%master) write (*,rt_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%w_gauss   = w_gauss
  self%p_gauss   = p_gauss
  self%amplitude = amplitude
  self%nv        = 2
  call self%patch_t%init
  ! Setup microphysics for RT. One microphysics object per wavelength is needed
  ! All new calls will get default value for wavelength, and this has to be updated
  ! "externally" if multi-frequency RT is desired.
  call self%microphysics%init
  !-----------------------------------------------------------------------------
  ! Memory mapping og density and temperature
  !-----------------------------------------------------------------------------
  call   self%d%init (self%gn, self%size/self%n, 'd')
  call  self%tt%init (self%gn, self%size/self%n, 'tt')
  self%d%f  => self%mem(:,:,:,i_d ,self%it,1)
  self%tt%f => self%mem(:,:,:,i_tr,self%it,1)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: gaussian top hat
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,self%gn(3)
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      y = (j-(self%gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      z = (k-(self%gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      self%d%f(i,j,k) = &
        exp(self%amplitude*exp(-(sqrt(x**2+y**2+z**2)/self%w_gauss)**self%p_gauss))
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Initial and boundary conditions for temperature [Kelvin]
  !-----------------------------------------------------------------------------
  self%tt%f = 10.
  !-----------------------------------------------------------------------------
  ! Set up the right dimensions from current patch, etc for radiative transfer
  !-----------------------------------------------------------------------------
  call self%microphysics%init
  call self%rt%init (solver)
  !
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Compute the solution to the radiative transfer equation
!> =============================================================================
SUBROUTINE solve (self)
  class(pde_t):: self
  real, pointer, dimension(:,:,:) :: qq
  real, pointer, dimension(:,:) :: dtau, source, qcool
  real, pointer, dimension(:) :: ifw0, irv0, ifw1, irv1
  integer:: iz
  !
  !.............................................................................
  call trace_begin ('radiation_mod:pde_t%update', 2)
associate(gn=>self%gn, n=>self%n, l=>self%li, u=>self%ui)
  allocate (ifw0(n(1)), irv0(n(1)))
  allocate (ifw1(n(1)), irv1(n(1)))
  allocate (dtau  (n(1),n(2)))
  allocate (source(n(1),n(2)))
  allocate (qcool (n(1),n(2)))
  allocate (qq    (n(1),n(2),n(3)))
  !.............................................................................
  call tic
  !dtau = self%dtau (direction, ...)
  do iz=l(3),u(3)
    ifw0 = 0.0; irv1 = 0.0
    call self%rt%solver (n(1), n(2), dtau, source, ifw0, ifw1, qcool, irv0, irv1)
    qq(:,:,iz) = qcool
  end do
  deallocate (dtau, source, qcool, ifw0, irv0, ifw1, irv1)
  deallocate (qq)
  call toc ('radiation', product(self%n))
end associate
  call trace_end
END SUBROUTINE solve
!> =============================================================================
END MODULE pde_mod
