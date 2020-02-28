!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE math_mod
  implicit none
  private
  ! ----------------------------------------------------------------------------
  ! Unit conversion factors between MHD and PIC
  ! ----------------------------------------------------------------------------
  type, public, extends(units_t):: convert_t
    real(kind=8):: d = 1.0, e = 1.0, ee = 1.0
    real(kind=8):: p = 1.0, u = 1.0, temp = 1.0 
    real(kind=8):: grav = 1.0, mu = 1.0, kr = 1.0
    real(kind=8):: b = 1.0, emf = 1.0, q = 1.0, momentum = 1.0, j = 1.0
  contains
    procedure:: init => init_convert
  end type
  ! ----------------------------------------------------------------------------
  ! Scaling for PIC
  ! ----------------------------------------------------------------------------
  type, public, extends(const_t):: pic_scaling_t
    real(kind=8):: l_scale = 1.0                  ! length scaling
    real(kind=8):: q_scale = 1.0                  ! electron charge scaling
    real(kind=8):: e_mass_scale = 1.0             ! electron mass  scaling
    real(kind=8):: eps_scale = 1.0                ! vacuum permittivity scaling
    real(kind=8):: mu0_scale = 1.0                ! vacuum permeability scaling
    real(kind=8):: l, d, l_debye,  mu, temp, v_e, n_p, t, u
    real(kind=8):: m,  B, wce, wci, wpe, wpi, delta, ek, v_a, p
    real(kind=8):: v_p, l_int, b0, cl, J, ve_drift, r_e, h, emf, ee, q_e
    real(kind=8):: ke, kb, kf
    real(kind=8):: q
    type(convert_t):: convert                     ! conversion factors MHD to PIC
  contains
    procedure:: init_pic => init_pic_cgs
  end type
  !-----------------------------------------------------------------------------
  ! Extend the units_t data type, which contains only the three fundamental
  ! units: length, mass, and time (cf. CGS and MKS)
  !-----------------------------------------------------------------------------
  type, public, extends(units_t):: scaling_t
    real(kind=8):: d, e, ee, p, u, temp, grav, mu, kr, b, emf, q, j
  contains
    procedure:: init
  end type
  type(scaling_t), public:: scaling
CONTAINS
!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace%begin ('scaling_t%init')
  if (io%master) print *, &
 '-------------------------------- scaling_mod -------------------------'
  first_time = .false.
  call cgs%init
  !--------------------------- defining quantities ---------------------
  self%system = 'cgs'                             ! pass on to Python
  self%l  = 1.0!1e5!1.7e7                         ! length = 1Mm = 1e8 cm
  self%d  = 1.0!d-9                               ! photospheric density
  self%t  = 1.0!0.1!24.89!1d0                     ! time = 100 sec
  !--------------------------- derived quantities ----------------------
  self%u  = self%l/self%t                         ! speed
  self%m  = self%d*self%l**3                      ! mass
  self%p  = self%d*self%u**2                      ! pressure [dyne / cm2]
  self%kr = 1.0/(self%d*self%l)                   ! rosseland opacity [cm2/g]
  self%ee = self%u**2                             ! specific energy
  self%e  = self%m*self%u**2                      ! energy density
  self%b = self%u*sqrt(self%d)
  !self%b  = self%l**(-0.5)*self%m**0.5*self%t**(-1)! magnetic flux density
  self%emf  = self%m**0.5/self%l**0.5/self%t      ! electric field
  self%q  = self%m**(1/2)*self%l**(3/2)/self%t    ! charge
  self%j  = self%m**0.5*self%l**0.5/self%t**2     ! current density
  self%mu = 0.5
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e)  ! temperature for given mu
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2      ! constant of gravity
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'               length: ',self%l
    print 1,'                 mass: ',self%m
    print 1,'                 time: ',self%t
    print 1,'             velocity: ',self%u
    print 1,'              density: ',self%d
    print 1,'             pressure: ',self%p
    print 1,'magnetic flux density: ',self%b
    print 1,'               energy: ',self%e
    print 1,'              gravity: ',self%grav
    print 1,'          temperature: ',self%temp
    print 1,'       electric field: ',self%emf
    print 1,'               charge: ',self%q
  1 format(1x,a,1p,e11.3)
  end if
  call self%output
  call trace%end()
END SUBROUTINE init
!===============================================================================
!> Initialize code units for PIC. 
!===============================================================================
SUBROUTINE init_pic_SI (self, ds)
  class(pic_scaling_t):: self
  real(8) :: ds, tmp
  ! ----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call si%init
  call trace%begin ('pic_scaling_t%init')
  if (io%master) print *, &
 '--------------------------------pic scaling_mod -------------------------'
  first_time = .false.
  self%name = 'SI'                                ! PIC scaled units | SI system
  ! --------------------------------------------------
  ! Defining quantities
  ! --------------------------------------------------
  self%m_e   = si%m_e    * self%e_mass_scale      ! electron mass
  self%m_u   = si%m_u                             ! proton mass
  self%q     = si%e      * self%q_scale           ! electron charge
  self%eps_0 = si%eps_0  * self%eps_scale         ! vacuum permittivity
  self%mu_0  = si%mu_0   * self%mu0_scale         ! vacuum permeability
  self%k_b     = si%k_b                           ! Boltzmann's constant
  self%mu      = 0.5                              ! mean molecular weight
  self%l       = scaling%l*1d-2                   ! length
  self%d       = scaling%d/1d3                    ! mass density
  self%t       = scaling%t                        ! time
  self%m       = self%l**3*self%d                 ! mass
  self%u       = self%l/self%t                    ! speed
  self%B       = self%l/self%t*sqrt(self%d)
  !self%B       = self%m/self%t/self%q
  self%e       = self%m*self%u**2                 ! energy
  self%ee      = self%u**2                        ! energy density
  self%p       = self%d*self%u**2                 ! pressure = energy density [dyne / cm2]
  self%q_e     = self%e/(self%B*self%l**2/self%t)
  self%temp    = 1.0!self%mu*(self%m_u/self%m)/(self%k_b/self%e)
  self%grav    = si%grav*(self%m/self%l**3)*self%t**2 ! constant of gravity
  ! ----------------------------------------------------------------------------
  ! We first have to scale down the physical constants for pic, then convert
  ! them to code units. We also need to rescale the MHD defining quantities
  ! ----------------------------------------------------------------------------
  self%c       = 1/sqrt(si%mu_0*si%eps_0)                     ! speed of light in SI
  self%c       = self%c*1/sqrt(self%mu0_scale*self%eps_scale) ! speed of light, with modified constants
  self%v_e     = sqrt(self%k_b*self%temp/self%m_e)            ! electron thermal velocity
  self%v_p     = sqrt(self%k_b*self%temp/self%m_u)            ! proton thermal velocity
  !self%l_debye = 0.5_8*ds*scaling%l/1d2                      ! Debye length, conversion CGS to SI
  self%n_p     = self%d / self%mu / si%m_u                    ! number density
  !
  !self%l_debye = 0.5_8*ds*scaling%l/1d2                      ! Debye length, conversion CGS to SI
  !self%q_e       = self%k_b*self%temp/(self%m_e*sqrt(self%n_p/self%m_e/self%eps_0))
  self%l_int   = 1./(self%n_p**(1./3.))                       ! interparticle distance
  !self%p       = 2.*self%n_p*self%k_b*self%temp              ! pressure
  self%v_a     = self%b/sqrt(self%mu_0*self%d)                ! alfven speed
  self%wce     = self%q*self%B/self%m_e                       ! electron gyro frequency
  self%wci     = self%q*self%B/self%m_u                       ! proton gyro frequency
  self%wpe     = self%q*sqrt(self%n_p/self%m_e/self%eps_0)    ! electron plasma frequency
  self%wpi     = self%q*sqrt(self%n_p/self%m_u/self%eps_0)    ! ion plasma frequency
  self%l_debye = self%v_e/self%wpe                            ! Debye length
  self%delta   = self%c/self%wpe                              ! electron skin depth
  self%h       = si%h_p  / self%m/self%u**2*self%t            ! Planck's constant
    if (io%master) then
    print 1,' scaled SI UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
  end if
  ! ----------------------------------------------------------------------------
  ! And now we scale the constants to code units
  ! ----------------------------------------------------------------------------
  self%m_e     = self%m_e     / self%m
  self%m_u     = self%m_u     / self%m                        ! proton mass
  self%n_p     = self%n_p     / self%l**3
  self%v_a     = self%l       / self%t
  self%u       = self%l       / self%t
  self%wce     = self%wce     / self%t
  self%wci     = self%wci     / self%t
  self%wpe     = self%wpe     / self%t
  self%wpi     = self%wpi     / self%t
  self%delta   = self%delta   / self%l
  self%v_e     = self%v_e     / self%u
  self%v_p     = self%v_p     / self%u
  self%c       = self%c       / self%u
  self%l_debye = self%l_debye / self%l
  self%l_int   = self%l_int   / self%l
  self%k_b     = self%k_b     / self%e
  self%q       = self%q       / self%q_e                      ! q scaling
  self%emf     = self%m*self%l/self%t**2/self%q
  self%J       = self%q/self%l**2/self%t
  self%eps_0 = self%eps_0 / (self%q_e*self%t/(self%l**3*self%B)) ! check!!
  !self%mu_0  = self%mu_0   / (self%m*self%l/self%q**2)          ! check!!
  self%mu_0  = self%mu_0  / (self%B*self%l*self%t/self%q_e)      ! check!!

  self%ke      = 1./(math%pi4*self%eps_0)
  self%kb      = self%mu_0/(math%pi4)
  self%kf      = 1.

  if (io%master) then
    print 1,' CODE UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,' scld electron_charge: ',self%q_e
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
  1 format(1x,a,1p,e11.3)
  end if
  call trace%end()
END SUBROUTINE init_pic_SI

!===============================================================================
!> Initialize code units for PIC. 
!===============================================================================
SUBROUTINE init_pic_SI2 (self, ds)
  class(pic_scaling_t):: self
  real(8) :: ds, tmp
  ! ----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call si%init
  call trace%begin ('pic_scaling_t%init')
  if (io%master) print *, &
 '--------------------------------pic scaling_mod -------------------------'
  first_time = .false.
  self%name = 'SI'                               ! PIC scaled units | SI system
  ! --------------------------------------------------
  ! Defining quantities
  ! --------------------------------------------------
  self%m_e     = si%m_e    * self%e_mass_scale   ! electron mass
  self%m_u     = si%m_u                          ! proton mass
  self%q       = si%e      * self%q_scale        ! electron charge
  self%eps_0   = si%eps_0  * self%eps_scale      ! vacuum permittivity
  self%mu_0    = si%mu_0   * self%mu0_scale      ! vacuum permeability
  self%m       = scaling%m/1d3                   ! mass
  self%l       = scaling%l*1d-2                  ! length
  self%t       = scaling%t                       ! time
  !!!!!!!!!!
  self%d       = self%m/self%l**2                ! density
  self%u       = self%l/self%t                   ! speed
  self%B       = self%m/self%t/self%q            ! magnetic induction
  self%e       = self%m*self%l**2/self%t**2      ! energy
  self%ee      = self%m/self%l/self%t**2         ! energy density
  self%p       = self%m/self%l/self%t**2         ! pressure = energy density
  
  self%k_b     = si%k_b                          ! Boltzmann's constant
  self%mu      = 0.5                             ! mean molecular weight
  self%temp    = self%mu*(self%m_u/self%m)/(self%k_b/self%e)  
  self%grav    = si%grav*(self%m/self%l**3)*self%t**2  ! constant of gravity
  ! ----------------------------------------------------------------------------
  ! We first have to scale down the physical constants for pic, then convert
  ! them to code units. We also need to rescale the MHD defining quantities
  ! ----------------------------------------------------------------------------
  self%c       = 1/sqrt(si%mu_0*si%eps_0)          ! speed of light in SI
  self%c       = self%c*1/sqrt(self%mu0_scale*self%eps_scale) ! modified constants
  self%v_e     = sqrt(self%k_b*self%temp/self%m_e) ! electron thermal velocity
  self%v_p     = sqrt(self%k_b*self%temp/self%m_u) ! proton thermal velocity
  !self%l_debye = 0.5_8*ds*scaling%l/1d2           ! Debye length, conversion CGS to SI
  self%n_p     = self%d / self%mu / self%m_u       ! number density
  !
  !self%l_debye = 0.5_8*ds*scaling%l/1d2           ! Debye length, conversion CGS to SI
  !self%q_e       = self%k_b*self%temp/(self%m_e*sqrt(self%n_p/self%m_e/self%eps_0))
  self%l_int   = 1./(self%n_p**(1./3.))            ! interparticle distance
  !self%p       = 2.*self%n_p*self%k_b*self%temp   ! pressure
  self%v_a     = self%b/sqrt(self%mu_0*self%d)     ! alfven speed
  self%wce     = self%q*self%B/self%m_e            ! electron gyro frequency
  self%wci     = self%q*self%B/self%m_u            ! proton gyro frequency
  self%wpe     = self%q*sqrt(self%n_p/self%m_e/self%eps_0) ! electron plasma frequency
  self%wpi     = self%q*sqrt(self%n_p/self%m_u/self%eps_0) ! ion plasma frequency
  self%l_debye = self%v_e/self%wpe                 ! Debye length
  self%delta   = self%c/self%wpe                   ! electron skin depth
  self%h       = si%h_p  / self%m/self%u**2*self%t ! Planck's constant
    if (io%master) then
    print 1,' scaled SI UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
  end if
  ! ----------------------------------------------------------------------------
  ! And now we scale the constants to code units
  ! ----------------------------------------------------------------------------
  self%m_e     = self%m_e     / self%m
  self%m_u     = self%m_u     / self%m                                          ! proton mass
  self%n_p     = self%n_p     / self%l**3
  self%v_a     = self%l       / self%t
  self%wce     = self%wce     / self%t
  self%wci     = self%wci     / self%t
  self%wpe     = self%wpe     / self%t
  self%wpi     = self%wpi     / self%t
  self%delta   = self%delta   / self%l
  self%v_e     = self%v_e     / self%u
  self%v_p     = self%v_p     / self%u
  self%c       = self%c       / self%u
  self%l_debye = self%l_debye / self%l
  self%l_int   = self%l_int   / self%l
  self%k_b     = self%k_b     / self%e
  self%emf     = self%m*self%l/self%t**2/self%q
  self%J       = self%q/self%l**2/self%t
  self%eps_0 = self%eps_0 / (self%q*self%t/(self%l**3*self%B))                  ! check!!
  !self%mu_0  = self%mu_0   / (self%m*self%l/self%q**2)                           ! check!!
  self%mu_0  = self%mu_0  / (self%B*self%l*self%t/self%q)                           ! check!!

  self%ke      = 1./(math%pi4*self%eps_0)
  self%kb      = self%mu_0/(math%pi4)
  self%kf      = 1.

  if (io%master) then
    print 1,' CODE UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,' scld electron_charge: ',self%q_e
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
  1 format(1x,a,1p,e11.3)
  end if
  call trace%end()
END SUBROUTINE init_pic_SI2

!===============================================================================
!> Initialize code units for PIC. 
!===============================================================================
SUBROUTINE init_pic_cgs (self, ds)
  class(pic_scaling_t):: self
  real(8) :: ds
  ! ----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace%begin ('pic_scaling_t%init')
  if (io%master) print *, &
 '--------------------------------pic scaling_mod -------------------------'
  first_time = .false.
  self%name = 'cgs'                              ! PIC scaled units | Gaussian
  ! --------------------------------------------------
  ! Defining quantities
  ! --------------------------------------------------
  self%m_e   = cgs%m_e    * self%e_mass_scale    ! electron mass
  self%m_u   = cgs%m_u                           ! proton mass
  self%q_e   = cgs%e      * self%q_scale         ! electron charge
  self%eps_0 = 1.                                ! vacuum permittivity
  self%mu_0  = 1.                                ! vacuum permeability
  self%mu      = 0.5                             ! mean molecular weight
  self%l       = scaling%l                       ! length
  self%d       = scaling%d                       ! mass density
  self%m       = self%l**3*self%d                ! mass
  self%temp    = scaling%temp
  self%B       = scaling%B
  self%t       = scaling%t                       ! time
  self%u        = self%l/self%t                  ! speed
  self%p        = self%d*self%u**2               ! pressure = energy density [dyne / cm2]
  self%e        = self%m*self%u**2               ! energy
  self%ee       = self%u**2                      ! energy density
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2  ! constant of gravity
  self%q       = self%m**0.5*self%l**(3./2.)/self%t
  self%j  = self%m**0.5*self%l**0.5/self%t**2        ! current density
  ! ----------------------------------------------------------------------------
  ! We first have to scale down the physical constants for pic, then convert
  ! them to code units. We also need to rescale the MHD defining quantities
  ! ----------------------------------------------------------------------------
  self%c       = 100*self%u                          ! speed of light
  self%k_b     = cgs%k_b                             ! Boltzmann's constant
  self%v_e     = sqrt(self%k_b*self%temp/self%m_e)   ! electron thermal velocity
  self%v_p     = sqrt(self%k_b*self%temp/self%m_u)   ! proton thermal velocity
  self%n_p     = self%d / self%mu / cgs%m_u          ! number density
  self%l_int    = 1/(self%n_p**(1/3))                ! interparticle distance
  self%p       = 2.*self%n_p*self%k_b*self%temp
  self%v_a     = self%b/sqrt(math%pi4*self%d)        ! alfven speed
  self%wce     = self%q*self%B/self%m_e              ! electron gyro frequency
  self%wci     = self%q*self%B/self%m_u              ! proton gyro frequency
  self%wpe     = self%q*sqrt(self%n_p/self%m_e/self%eps_0)  ! electron plasma frequency
  self%wpi     = self%q*sqrt(self%n_p/self%m_u/self%eps_0)  ! ion plasma frequency
  self%l_debye = self%v_e/self%wpe                          ! Debye length
  self%delta   = self%c/self%wpe                            ! electron skin depth
  self%h       = si%h_p  / self%m/self%u**2*self%wce        ! Planck's constant
    if (io%master) then
    print 1,' scaled SI UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,' der  electron_charge: ',self%q_e
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
  end if
  ! ----------------------------------------------------------------------------
  ! And now we scale to code units
  ! ----------------------------------------------------------------------------
  self%temp=scaling%temp
  self%m_e     = self%m_e     / self%m
  self%m_u     = self%m_u     / self%m                  ! proton mass
  self%n_p     = self%n_p     / self%l**3
  self%B       = self%u*sqrt(cgs%pi*4.0*self%d)
  !self%B       = self%m/self%t/self%q
  self%e       = self%m       * self%u**2
  self%v_a     = self%l       / self%t
  self%u       = self%l       / self%t
  self%p       = self%m       *self%v_a**2
  self%wce     = self%wce     / self%t
  self%wci     = self%wci     / self%t
  self%wpe     = self%wpe     / self%t
  self%wpi     = self%wpi     / self%t
  self%delta   = self%delta   / self%l
  self%v_e     = self%v_e     / self%u
  self%v_p     = self%v_p     / self%u
  self%c       = self%c       / self%u
  self%l_debye = self%l_debye / self%l
  self%l_int   = self%l_int   / self%l
  self%k_b     = self%k_b     / self%e
  self%emf     = self%m**0.5/self%l**0.5/self%t         ! electric field
  self%q_e     = self%q_e     / self%q                  ! q scaling
  self%q     = self%q_e                                 ! q scaling
  !self%emf     = self%m*self%l**3*self%B/self%t**3/self%e ! electric field
  self%eps_0 = 1.0
  self%mu_0  = 1.0

  self%ke      = 1.
  self%kb      = 1./self%c
  self%kf      = 1./self%c

  
  if (io%master) then
    print 1,' CODE UNITS:   (PIC):'
    print 1,'       number density: ',self%n_p
    print 1,'               length: ',self%l
    print 1,'      electron_charge: ',self%q
    print 1,' der  electron_charge: ',self%q_e
    print 1,'        electron_mass: ',self%m_e
    print 1,'vacuum permeability  : ',self%mu_0
    print 1,'vacuum permittivity  : ',self%eps_0
    print 1,'mean molecular weight: ',self%mu
    print 1,'temperature          : ',self%temp
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'energy               : ',self%e
    print 1,'atomic mass unit     : ',self%m_u
    print 1,'Boltzmanns constant  : ',self%k_b
    print 1,'speed of light       : ',self%c
    print 1,'magnetic field streng: ',self%B
    print 1,'electron gyro freq   : ',self%wce
    print 1,'proton gyro freq     : ',self%wci
    print 1,'electron plasma freq : ',self%wpe
    print 1,'ion plasma freq      : ',self%wpi
    print 1,'electron skin depth  : ',self%delta
    print 1,'electron thermal vel.: ',self%v_e
    print 1,'kinetic energy       : ',self%ek
    print 1,'alfven speed         : ',self%v_a
    print 1,'pressure             : ',self%p
    print 1,'proton thermal vel   : ',self%v_p
    print 1,'Debye length         : ',self%l_debye
    print 1,'Debye length ref     : ',ds
    print 1,'interparticle dist   : ',self%l_int
    print 1,'b0 in Rutherford scat: ',self%b0
    print 1,'coulomb log          : ',self%cl
    print 1,'electric curent      : ',self%J
    print 1,'electron drift speed : ',self%ve_drift
    print 1,'electron gyro radius : ',self%r_e
    print 1,'Planck constant      : ',self%h
    print 1,'EMF                  : ',self%emf
  1 format(1x,a,1p,e11.3)
  end if
  call trace%end()
END SUBROUTINE init_pic_cgs

SUBROUTINE init_convert (self, pic)
  class(convert_t):: self
  class(pic_scaling_t):: pic
  ! ----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace%begin ('convert_t%init')
  if (io%master) print *, &
 '--------------------------------unit conversion factors -------------------------'
  first_time = .false.
  ! ----------------------------------------------------------------------------
  ! select the conversion factor case. More conversions and systems will be
  ! added on demand or per-need basis.
  ! ----------------------------------------------------------------------------
  select case(trim(pic%name))
  case('SI')
    self%l   = scaling%l   / pic%l   * 1d-2         ! length  CGS to SI - based code units
    self%d   = scaling%d   / pic%d   * 1d3          ! density CGS to SI - based code units
    self%t   = scaling%t   / pic%t                  ! time unit, if it is different
    self%u   = scaling%u   / pic%u   * 1d-2         ! velocity
    self%m   = scaling%m   / pic%m   * 1d-3         ! mass
    self%p   = scaling%p   / pic%p   * 1d-1         ! pressure
    self%ee  = scaling%ee  / pic%ee  * 1d-1         ! energy density
    self%e   = scaling%e   / pic%e   * 1d-7         ! energy
    self%b   = scaling%b   / pic%b   * 1d-4         ! magnetic induction
    self%emf = scaling%emf / pic%emf / (1./3.*1d-4) ! electric field
    self%j   = scaling%j   / pic%j   / 3d5          ! current density
    self%temp= scaling%temp/ pic%temp               ! temperature for given mu and scaled Boltzmann constant
    self%momentum = scaling%m*scaling%u / pic%m / pic%u * 1d-5
    self%grav= scaling%grav/ pic%grav               ! constant of gravity
  case('cgs')
    self%l   = 1.0!scaling%l   / pic%l              ! length  CGS to SI - based code units
    self%d   = 1.0!scaling%d   / pic%d              ! density CGS to SI - based code units
    self%t   = 1.0!scaling%t   / pic%t              ! time unit, if it is different
    self%u   = 1.0!scaling%u   / pic%u              ! velocity
    self%m   = 1.0!scaling%m   / pic%m              ! mass
    self%p   = 1.0!scaling%p   / pic%p              ! pressure
    self%ee  = 1.0!scaling%ee  / pic%ee             ! energy density
    self%e   = 1.0!scaling%e   / pic%e              ! energy
    self%b   = 1.0!scaling%b   / pic%b              ! magnetic induction
    self%emf = 1.0!scaling%emf / pic%emf            ! electric field
    self%j   = 1.0!scaling%j   / pic%j              ! current density
    self%temp= 1.0!scaling%temp/ pic%temp           ! temperature for given mu and scaled Boltzmann constant
    self%momentum = scaling%m*scaling%u / pic%m / pic%u 
    self%grav= scaling%grav/ pic%grav      
  case default
    call io%abort ("scaling_mod%init_convert:: PIC system not understood")
  end select

  if (io%master) then
    print 1,' MHD-PIC conversion:'
    print 1,'               length: ',self%l
    print 1,'mass density         : ',self%d
    print 1,'time                 : ',self%t
    print 1,'speed                : ',self%u
    print 1,'mass                 : ',self%m
    print 1,'pressure             : ',self%p
    print 1,'energy density       : ',self%ee
    print 1,'energy               : ',self%e
    print 1,'magnetic induction   : ',self%b
    print 1,'electric field       : ',self%emf
    print 1,'current density      : ',self%j
    print 1,'temperature          : ',self%temp
    print 1,'momentum             : ',self%momentum
    print 1,'gravity              : ',self%grav
  1 format(1x,a,1p,e11.3)
  end if
  call trace%end()
END SUBROUTINE init_convert

END MODULE scaling_mod

