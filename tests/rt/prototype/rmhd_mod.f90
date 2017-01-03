!> =====================================================================
!> This module combines the MHD and RT solver solutions, to produce
!> time derivatives that include the net radiative cooling
!>
!> Radiative transfer module structure:
!>
!> experiment_mod                       ! generic time advance modul
!>   rmhd_mod                           ! wrapper for MHD and RT
!>     mhd_mod                          ! MHD solver
!>     eos_mod                          ! equation-of-state module
!>     rt_mod                           ! RT organizer
!>       rt_bundle_mod                  ! handles bundles of rays
!>         rt_solver_mod                ! solver the RT equation
!> =====================================================================
MODULE rmhd_mod
  USE io_mod
  USE trace_mod
  USE eos_mod
  USE mhd_mod
  USE rt_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: rmhd_t
    type(rt_t)                        :: rt
    real, dimension(:,:,:),   pointer :: lnd, lne
    real, dimension(:,:,:,:), pointer :: rk, src
  contains
    procedure :: init
    procedure :: update
    procedure :: solve
    procedure :: log_interpolate
    procedure :: time_interval
  end type
  real, parameter:: pi4 = 16.*atan(1.)
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialization of radiative transfer, preparing for two main tasks:
!> 1) most of the time, just inter-/extrapolate in the exising RT solutions
!> 2) on the rt%dtime intervals, solve the RT and add the new net 
!>    heating / cooling in a time slot
!> =====================================================================
SUBROUTINE init (self)
  class (rmhd_t) :: self
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t_%init')
  !---------------------------------------------------------------------
  ! Initialize the MHD data type, and allocate arrays for RT
  !---------------------------------------------------------------------
  call self%mhd_t%init
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  allocate (self%lnd (mx,my,mz), &                                      ! log density
            self%lne (mx,my,mz), &                                      ! log internal energy
            self%rk  (mx,my,mz,eos%n_lambda), &                         ! opacity
            self%src (mx,my,mz,eos%n_lambda))                           ! source function
  end associate
  !---------------------------------------------------------------------
  ! Initialize the rt_mod object (one could create more than one)
  !---------------------------------------------------------------------
  self%rt%dim      = self%mesh%gn                                       ! RT dimensions
  self%rt%n_lambda = eos%n_lambda                                       ! number of wavelengths
  self%rt%n_bundle = 1                          ! FIXME                 ! one direction
  self%rt%ds       => self%mesh(2)%dsup                                 ! mesh spacing
  call self%rt%init                                                     ! allocate weights
  self%rt%w_lambda = 1.0/self%rt%n_lambda       ! FIXME                 ! wavelength weights
  self%rt%w_bundle = pi4/2./self%rt%n_bundle    ! FIXME                 ! space angle weights
  !
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> Update the PDE time derivatives, by calling the MHD module, and
!> adjusting by time-interpolation of the RT memory values
!> =====================================================================
SUBROUTINE update (self)
  class (rmhd_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%update')
  call self%mhd_t%update                                                ! get df/dt from MHD
  if (.not.self%rt%is_ahead_of (self)) then                             ! if RT is not up-to-date
    call self%solve                                                     ! solve the RT for a new slot
  end if
  call self%rt%interpolate (self%time)                                  ! interpolate in time and space
  print *, shape(self%rt%dedt)
  self%dedt = self%dedt + self%rt%dedt                                  ! add the RT contribution to de/dt
  call trace%end
END SUBROUTINE update

!> =====================================================================
!> Solve the RT problem, by first obtaining the thermodynamics at the
!> designated RT time, and using them to compute the opacity and source
!> function values.  Then solve for the net heating rate, storing it in
!> a slot that rotates with time.
!> =====================================================================
SUBROUTINE solve (self)
  class (rmhd_t) :: self
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%solve')
  self%rt%time = self%rt%time + self%rt%dtime                           ! increment RT solver time
  call self%log_interpolate (self%rt%time, self%i_d, self%lnd)          ! interpolate in time and space
  call self%log_interpolate (self%rt%time, self%i_e, self%lne)          ! interpolate in time and space
  self%lne = self%lne - self%lnd                                        ! internal energy per unit mass
  call eos%lookup (shape(self%lne), lne=self%lne, lnd=self%lnd, &       ! call the table lookup, ...
                   rk=self%rt%rk, src=self%rt%src)                      ! ... results go into rt_t
  call self%rt%solve                                                    ! all the RT solver to get Q
  call self%rotate                                                      ! rotate the time slots
  call trace%end
END SUBROUTINE solve

!> =====================================================================
!> Interpolate -- only in time for now -- from the MHD time slots save in
!> this module to the time requested from the solve procedure
!> =====================================================================
SUBROUTINE log_interpolate (self, time, i, a)
  class (rmhd_t)                  :: self
  real(8)                         :: time
  integer                         :: i
  real, dimension(:,:,:), pointer :: a
  !.....................................................................
  integer                         :: jt(2)
  real                            :: pt(2)
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%log_interpolate')
  call self%time_interval (time, jt, pt)
  a = log(self%mem(:,:,:,i,jt(1),1))*pt(1) &
    + log(self%mem(:,:,:,i,jt(2),1))*pt(2)
  call trace%end
END SUBROUTINE log_interpolate

!> =====================================================================
!> Search for a time interval to interpolate in, and return the time
!> slot indices and weight factors
!> =====================================================================
SUBROUTINE time_interval (self, time, jt, pt)
  class (rmhd_t)                  :: self
  real(8)                         :: time
  integer                         :: jt(2)
  real                            :: pt(2)
  !.....................................................................
  integer                         :: i, j
  real(8)                         :: t(self%nt), dt
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%time_interval')
  t = self%t(self%iit)
  j = 1
  do i=1,self%nt-2
    if (t(i) < time) j=i
  end do
  jt = self%iit(j:j+1)
  dt = t(j+1)-t(j)
  if (dt==0.0) then
    pt(2) = 1.0
  else
    pt(2) = (time-t(j))/dt
  end if
  pt(1) = 1.-pt(2)
  print *, 'time_interval: times  =', t(j), time, t(j+1)
  print *, 'time_interval: times  =', self%t(jt(1)), time, self%t(jt(2))
  print *, 'time_interval: jt, pt =', jt, pt
  call trace%end
END SUBROUTINE time_interval

END MODULE rmhd_mod
