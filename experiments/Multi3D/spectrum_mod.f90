!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE mconst
    USE op_math_mod, only: air2vac
    implicit none

    !-------------------------------------------------------------------------------
    type :: spectrum_type
        ! total number of frequencies
        integer :: nnu=0
        ! frequency [Hz]
        real(8),dimension(:),allocatable :: nu
    contains
      procedure :: init
      procedure :: set_weights
    end type spectrum_type
    !---------------------------------------------------------------------

  contains
    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), save:: daa=0, aa_red=0, aa_blue=0
      integer :: iostat, i, nnu
      namelist /spectrum_params/ daa, aa_red, aa_blue
      !-----------------------------------------------------------------------------
      rewind (io%input)
      read (io%input, spectrum_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, spectrum_params)
      end if
      call io%namelist_warning ('spectrum_params', iostat)
      write (io%namelist, spectrum_params)

      if (daa==0) return

      aa_red = air2vac(aa_red)
      aa_blue = air2vac(aa_blue)

      nnu = ceiling((aa_red - aa_blue) / daa)
      self%nnu = nnu
      allocate(self%nu(nnu))
      do i=0,nnu-1
        self%nu(nnu-i) = CC / (aa_blue + i*daa) * 1e8
      end do
    END SUBROUTINE init

    !===============================================================================
    SUBROUTINE set_weights (self, weight, nu, start, end)
      class(spectrum_type) :: self
      real(8), dimension(:), allocatable :: weight
      real(8), dimension(:), intent(in) :: nu
      integer, intent(in) :: start, end
      integer :: j, n
      !-----------------------------------------------------------------------------
      n = size(nu, 1)
      if (allocated(weight)) deallocate(weight)
      allocate( weight(start:end) )

      weight(start) = 0.5 * ( nu(2) - nu(1) )

      do  j=2,n-1
          weight(start - 1 + j) = 0.5* ( nu(j+1) - nu(j-1) )
      end do

      weight(end) = 0.5 * ( nu(n) - nu(n-1) )
    END SUBROUTINE set_weights

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  