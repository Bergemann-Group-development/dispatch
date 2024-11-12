!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE const_mod
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
        logical :: sorted=.false.
        character(len=32) :: name=''
    contains
      procedure :: init
      procedure :: get_interval
      procedure :: save
    end type spectrum_type
    !---------------------------------------------------------------------

    public :: set_weights
    real(8), public :: nu_tolerance = 1d6
    logical :: give_lam=.false.
    logical :: in_log=.false.
    real(8) :: daa=0, aa_red=0, aa_blue=0, R=0
    integer :: nnu=0
    character(len=256) :: lam_file=''

  contains
    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), dimension(:), allocatable :: lam
      real(8) :: last, fac
      integer :: iostat, i, lunit
      namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam, in_log, R, lam_file
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

      if (daa==0 .and. nnu==0 .and. R==0 .and. lam_file=='') return
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
      elseif (lam_file /= '') then
        nnu = 0
        open(newunit=lunit, file=trim(lam_file), form='FORMATTED')
        read(lunit,'(A)', iostat=iostat)
        do while (iostat == 0)
          nnu = nnu + 1
          read(lunit,'(A)', iostat=iostat)
        end do
        rewind(lunit)
        allocate(lam(nnu))
        allocate(self%nu(nnu))
        read(lunit,*) lam
        close(lunit)
        self%nnu = nnu
        lam = air2vac(lam)
        do i=0,nnu-1
          self%nu(nnu-i) = CC / (lam(i+1)) * 1e8
        end do
      elseif (aa_red /= 0 .and. aa_blue /= 0) then
        if (in_log) then
          aa_red = log(air2vac(exp(aa_red)))
          aa_blue = log(air2vac(exp(aa_blue)))
        else
          aa_red = air2vac(aa_red)
          aa_blue = air2vac(aa_blue)
        endif
        if (R /= 0) then
          fac = 1d0 + 1d0/R
          nnu = ceiling(log(aa_red/aa_blue) / log(fac))
          allocate(self%nu(nnu))
          last = aa_blue
          do i=0,nnu-1
            last = last * fac
            self%nu(nnu-i) = last
          end do
        else
          if (daa /= 0) then
            nnu = ceiling((aa_red - aa_blue) / daa)
          else
            daa = (aa_red - aa_blue) / real(nnu-1)
          endif
          allocate(self%nu(nnu))
          do i=0,nnu-1
            if (in_log) then
              self%nu(nnu-i) = exp(aa_blue + i*daa)
            else
              self%nu(nnu-i) = aa_blue + i*daa
            endif
          end do
        endif
        self%nu = CC / self%nu * 1e8
        self%nnu = nnu
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
      if (n==0) then
        return
      elseif (n==1) then
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
      class(spectrum_type), intent(inout):: self
      character(len=*), intent(in) :: label
      character(len=256) :: filename
      integer :: lunit
      logical :: exists
      !-----------------------------------------------------------------------------
      filename = trim(io%run_dir) // 'save/' // trim(label) // '_lam.bin'
      inquire(file=filename, exist=exists)
      if (.not. exists) then
        open(newunit=lunit, file=filename, access='stream')
        write(lunit) 1d8 * dble(cc) / self%nu
        close(lunit)
      endif
    END SUBROUTINE save

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  