!===============================================================================
!> Opacity module for the piece-wise analytic opacities of Bell & Lin
!> (1994, ApJ, 427, 987).
!===============================================================================
MODULE opacities_mod
  USE trace_mod
  USE mpi_mod
  implicit none
  private

  type, public::  opacity_t
    real(4):: grain_rho, grain_radius
    real(4), dimension(:), allocatable:: wavelength, kappa_lamb
  contains
    procedure:: init => init_draine_lee
    procedure:: solve => void
    procedure:: deallocate => deallocate_kappa
  end type

  real(8):: pi = asin(1.0) * 2.0

CONTAINS
!===============================================================================
SUBROUTINE void (self)
  class(opacity_t):: self
  call mpi%abort('opacity_t%solve: NotYetImplemented Error!')
END SUBROUTINE void

!===============================================================================
!> Read the provided opacity data file.
!> Reference: Draine & Lee (1984), astronomical silicates with grain radius
!> 0.12 microns and dust bulk density of 3.6 g cm**-3.
!===============================================================================
SUBROUTINE init_draine_lee (self)
  class(opacity_t):: self
  logical:: lexists
  integer:: nbin_act, nheader, ios, i
  character(len=9):: opacfile
  character(len=1):: dummy
  character(len=64):: error_msg
  real(4), dimension(:), allocatable:: ext_cs, scat_cs
  real(4):: q1
  ! --------------------------------------------------------------------
  call trace_begin('opacity_t%init')

  nheader = 9
  opacfile = 'optSi.dat'
  self%grain_rho = 3.6 ! g cm**-3

  inquire (file=opacfile, iostat=ios, exist=lexists)
  if (.not. lexists) then
    write (error_msg,'(a,a,a)') 'Opacity file ', trim(opacfile), ' not found! Abort!'
    call mpi%abort(trim(error_msg))
  end if
  open (unit=99, file=opacfile, status='old')

  ! skip past header
  do i=1,nheader
    read (99)
  end do
  ! read the grain radius, the number of wavelength bins and then allocate arrays
  read (99, fmt='(a,e16.8)') dummy, self%grain_radius ! microns
  self%grain_radius = self%grain_radius * 1.0e-4 ! convert to cm
  read (99, fmt='(a,i5.1)') dummy, nbin_act
  allocate (self%wavelength(nbin_act), self%kappa_lamb(nbin_act))
  allocate (scat_cs(nbin_act), ext_cs(nbin_act))

  ! read in the wavelength-dependent opacity data
  i = 1
  do
    read (99, fmt='(e16.8,e16.8,e16.8)', iostat=ios) self%wavelength(i), scat_cs(i), ext_cs(i)
    if (ios > 0) then
      call mpi%abort('There was a problem reading the opacity file. Abort!')
    else if (ios < 0) then
      exit ! end of file reached
    end if
    i = i + 1
  end do

  ! Calculate the absorption opacity at each wavelength. This requires
  ! knowledge of the grain radius and bulk density.
  q1 = 0.75 / pi / self%grain_rho / self%grain_radius**3
  do i=1,nbin_act
    ! the factor of 100 converts the cross-sections from m**2 to cm**2.
    ! the wavelength is also converted from microns to cm.
    ext_cs(i) = ext_cs(i) * 100.0**2
    scat_cs(i) = scat_cs(i) * 100.0**2
    self%wavelength(i) = self%wavelength(i) * 1.0e-4

    ! NOTE!: No correction for the dust-to-gas ratio is made here!
    self%kappa_lamb(i) = (ext_cs(i) - scat_cs(i)) * q1
  end do

  close (99)
  deallocate (scat_cs, ext_cs)
  call trace_end
END SUBROUTINE init_draine_lee

!===============================================================================
!> Deallocate memory used for the opacity.
!===============================================================================
SUBROUTINE deallocate_kappa (self)
  class(opacity_t):: self
  deallocate(self%kappa_lamb)
END SUBROUTINE deallocate_kappa

END MODULE opacities_mod
