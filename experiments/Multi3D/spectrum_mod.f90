!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE mconst
    USE handyfuncs_mod, only:air2vac
    USE op_math_mod, only: get_ivl=>get_interval
    USE trans_mod, only: transition_type
    implicit none

    !-------------------------------------------------------------------------------
    type :: spectrum_type
        ! total number of frequencies
        integer :: nnu=0
        ! frequency [Hz]
        real(8), dimension(:), allocatable :: nu, wnu
        logical :: sorted=.false., saved=.false.
        character(len=32) :: name=''
    contains
      procedure :: init
      procedure :: get_interval
      procedure :: save
    end type spectrum_type
    !---------------------------------------------------------------------

    public :: set_weights
    real(8), public :: nu_tolerance = 1d6

  contains
    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), save:: daa=0, aa_red=0, aa_blue=0
      real(8), dimension(:), allocatable :: lam
      integer :: iostat, i, nnu=0
      logical :: give_lam=.false.
      namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam
      namelist /sp_wavelengths/ lam
      !-----------------------------------------------------------------------------
      rewind (io%input)
      read (io%input, spectrum_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, spectrum_params)
      end if
      call io%namelist_warning ('spectrum_params', iostat)
      write (io%namelist, spectrum_params)

      if (daa==0 .and. nnu==0) return
      if (give_lam .and. nnu>0) then
        allocate(lam(nnu))
        rewind (io%input)
        read (io%input, sp_wavelengths, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, sp_wavelengths)
        end if
        call io%namelist_warning ('spectrum_params', iostat)
        write (io%namelist, sp_wavelengths)
        self%nnu = nnu
        allocate(self%nu(nnu))
        lam = air2vac(lam)
        do i=0,nnu-1
          self%nu(nnu-i) = CC / (lam(i+1)) * 1e8
        end do
      elseif (aa_red /= 0 .and. aa_blue /= 0) then
        aa_red = air2vac(aa_red)
        aa_blue = air2vac(aa_blue)
        if (daa /= 0) then
          nnu = ceiling((aa_red - aa_blue) / daa)
        else
          daa = (aa_red - aa_blue) / real(nnu-1)
        endif
        self%nnu = nnu
        allocate(self%nu(nnu))
        do i=0,nnu-1
          self%nu(nnu-i) = CC / (aa_blue + i*daa) * 1e8
        end do
      else
        stop 'Give either wavelength set or aa_red and aa_blue in spectrum_params'
      endif
      allocate(self%wnu(nnu))
      call set_weights(self%wnu, self%nu)
      self%sorted=.true.
    END SUBROUTINE init

    !===============================================================================
    SUBROUTINE set_weights (weight, nu)
      real(8), dimension(:), intent(out) :: weight
      real(8), dimension(:), intent(in) :: nu
      integer :: j, n
      !-----------------------------------------------------------------------------
      n = size(nu, 1)
      if (n==1) then
        print *, 'weight for single frequency wanted!'
        weight(1) = 1.0
        return
      elseif (size(weight) /= n) then
        stop 'weight and nu do not have the same number of elements!'
      endif
  
      weight(1) = 0.5 * ( nu(2) - nu(1) )
      do j=2,n-1
          weight(j) = 0.5* ( nu(j+1) - nu(j-1) )
      end do
      weight(n) = 0.5 * ( nu(n) - nu(n-1) )
    END SUBROUTINE set_weights

    !===============================================================================
    SUBROUTINE get_interval (self, trans, ired, iblue)
      class(spectrum_type), intent(in):: self
      class(transition_type), intent(in) :: trans
      integer, intent(out) :: ired, iblue
      !-----------------------------------------------------------------------------
      if (self%sorted) then
        call get_ivl(self%nu, trans%nu_min, trans%nu_max, ired, iblue, limit=.false.)
        if (ired > 1) then
          if((trans%nu_min - self%nu(ired-1)) < nu_tolerance) then
            ired = ired - 1
          endif
        endif
        if (iblue < self%nnu) then
          if((self%nu(iblue+1) - trans%nu_max) < nu_tolerance) then
            iblue = iblue + 1
          endif
        endif
      else
        ired = trans%ired
        iblue = trans%iblue
      endif
    END SUBROUTINE get_interval

    !===============================================================================
    SUBROUTINE save (self, label)
      class(spectrum_type), intent(in):: self
      character(len=*), intent(in) :: label
      character(len=128) :: filename
      integer :: lunit
      !-----------------------------------------------------------------------------
      if (.not. self%saved) then
        filename = trim(io%run_dir) // 'save/' // trim(label) // '_lam.bin'
        open(newunit=lunit, file=filename, access='stream')
        write(lunit) 1d8 * dble(cc) / self%nu
        close(lunit)
      endif
    END SUBROUTINE save

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  