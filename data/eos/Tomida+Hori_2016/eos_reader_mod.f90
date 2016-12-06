!#-----------------------------------------
!# RHO-T table based on ideal EOS.
!#    log RHO = -22 to 1 , dlogRHO = 0.05
!#      log T = 0.4 to 8 ,   dlogT = 0.01
!# All the data were computed under the assumption of ideal gas.
!# We do not make use of e = Cv*T.
!# Ionization of He is not taken into account. 
!# 
!# by Kengo Tomida and Yasunori Hori,  Jul. 1st ,2016.
!#    X = 0.7, Y = 0.3 : H2 ortho : para = equilibrium
!#
!# 1st column : log RHO (g/cc) -- Density
!# 2nd column : log T (K) -- Temperature
!# 3rd column : log P (dyn/cm2) -- Pressure
!# 4th column : log U (erg/g) -- Specific Internal energy
!# 5th column : log cT (cm/s) -- isothermal sound speed
!# 6th column : log cS (cm/s) -- adiabatic sound speed
!# 7th column : gamma  -- (dlogP/dlogRHO)s
!# 8th column : log Cv (erg/K cc) -- heat capacity
!#----------------------------------------- 
!# 761  41 801
!> =============================================================================
!> =============================================================================
MODULE eos_reader_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
  type, public:: eos_reader_t
    integer:: nd=461, nT=761, nv=3
    real(8), dimension(:), pointer:: T, d
    real(8), dimension(:,:,:), pointer:: table
  contains
    procedure:: init
  end type
  type (eos_reader_t), public:: eos_reader
CONTAINS

!> ============================================================================
!> Skip lines
!> ============================================================================
SUBROUTINE skip_lines (n)
  integer:: n, i
  character(len=4):: line
  do i=1,n
    read (io_unit%tmp,'(a)') line
  end do
END SUBROUTINE skip_lines

!> =============================================================================
!> Read in the data from tomida.dat
!> =============================================================================
SUBROUTINE init (self)
  class (eos_reader_t):: self
  real:: row(8)
  integer:: id, iT
  associate (nT=>self%nT, nd=>self%nd, nv=>self%nv)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_reader_t%init')
  open (io_unit%tmp, form='formatted', status='old', &
    file=trim(io_unit%top)//'/data/eos/Tomida+Hori_2016/tomida.dat')
  !-----------------------------------------------------------------------------
  allocate (self%table(nd,nT,3), self%T(nT), self%d(nd))
  call skip_lines (19)
  do id=1,nd
    call skip_lines (2)
    do iT=1,nT
      ! row =            log_d, log_T, log_P, log_U, log_cT, log_cS, gamma, log_Cv
      read (io_unit%tmp,*) row
      self%table(iT,id,1) = 10.0*row(3)         ! pressure       P(T,rho)
      self%table(iT,id,2) = 10.0*row(4)         ! temperature    U(T,rho)
      self%table(iT,id,3) = row(7)              ! gas gamma   amma(T,rho)
      self%T(iT) = 10.**row(2)                  ! Temperature
    end do
    self%d(id) = 10.**row(1)                    ! density
  end do
  close (io_unit%tmp)
  print *,'log(T):', log10(eos_reader%T(1)), log10(eos_reader%T(nT))
  print *,'log(d):', log10(eos_reader%d(1)), log10(eos_reader%d(nd))
  call trace%end
  end associate
END SUBROUTINE init

END MODULE eos_reader_mod
