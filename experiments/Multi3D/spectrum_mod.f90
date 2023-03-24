!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE mconst
    USE op_math_mod, only: air2vac, get_ivl=>get_interval
    USE trans_mod, only: transition_type
    implicit none

    !-------------------------------------------------------------------------------
    type :: spectrum_type
        ! total number of frequencies
        integer :: nnu=0
        ! frequency [Hz]
        real(8), dimension(:), allocatable :: nu
        logical :: sorted=.false.
        character(len=32) :: name=''
    contains
      procedure :: init
      procedure :: get_interval
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
    SUBROUTINE get_interval (self, trans, ired, iblue)
      class(spectrum_type), intent(in):: self
      class(transition_type), intent(in) :: trans
      integer, intent(out) :: ired, iblue
      !-----------------------------------------------------------------------------
      if (self%sorted) then
        call get_ivl(self%nu, trans%nu_min, trans%nu_max, ired, iblue, limit=.false.)
      else
        ired = trans%ired
        iblue = trans%iblue
      endif
    END SUBROUTINE get_interval

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  