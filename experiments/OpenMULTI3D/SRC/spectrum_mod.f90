!---------------------------------------------------------------------
MODULE spectrum_mod 
    USE io_mod
    USE os_mod
    USE timer_mod
    USE trace_mod
    USE const_mod
    USE handyfuncs_mod, only: air2vac, remove_close_entries
    USE op_math_mod, only: get_ivl=>get_interval
    USE trans_mod, only: transition_type
    USE mmath_mod, only: refsor
    implicit none

    !-------------------------------------------------------------------------------
    type :: spectrum_type
        ! total number of wavelengths
        integer :: nnu=0
        real(8), dimension(:), allocatable :: nu, wnu, lam, llam
        logical :: sorted=.false.
        character(len=32) :: name=''
    contains
      procedure :: init
      procedure :: reset
      procedure :: get_interval
      procedure :: sort
      procedure :: save
    end type spectrum_type
    !---------------------------------------------------------------------

    public :: set_weights
    real(8), public :: lam_tolerance = 0.001
    logical :: give_lam=.false., in_air=.true.
    logical :: in_log=.false.
    real(8) :: daa=0, aa_red=0, aa_blue=0, R=0
    integer :: nnu=0
    character(len=256) :: lam_file=''

  contains
    !===============================================================================
    SUBROUTINE init(self)
      class(spectrum_type) :: self
      real(8), dimension(:), allocatable :: lam
      real(8) :: fac
      integer :: iostat, i, lunit
      namelist /spectrum_params/ daa, aa_red, aa_blue, nnu, give_lam, in_log, R, &
      lam_file, in_air
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
      elseif (lam_file /= '') then
        nnu = 0
        call os%check_file(lam_file)
        open(newunit=lunit, file=trim(lam_file), form='FORMATTED')
        read(lunit,'(A)', iostat=iostat)
        do while (iostat == 0)
          nnu = nnu + 1
          read(lunit,'(A)', iostat=iostat)
        end do
        rewind(lunit)
        allocate(lam(nnu))
        read(lunit,*) lam
        close(lunit)
      elseif (aa_red /= 0 .and. aa_blue /= 0) then
        if (R /= 0) then
          fac = 1d0 + 1d0/R
          nnu = ceiling(log(aa_red/aa_blue) / log(fac))
          allocate(lam(nnu))
          lam(1) = aa_blue
          do i=2, nnu
            lam(i) = lam(i-1) * fac
          end do
        else
          if (daa /= 0) then
            nnu = nint((aa_red - aa_blue) / daa) + 1
          else
            daa = (aa_red - aa_blue) / real(nnu-1)
          endif
          allocate(lam(nnu))
          lam = [(aa_blue + i*daa, i=0,nnu-1)]
        endif
      else
        stop 'Give either wavelength set or aa_red and aa_blue in spectrum_params'
      endif
      if (in_log) lam = exp(lam)
      if (in_air) lam = air2vac(lam)
      call self%reset(lam=lam)
      allocate(self%wnu(nnu))
      call set_weights(self%wnu, self%nu)
      self%sorted=.true.
      deallocate(lam)
    END SUBROUTINE init

    !===============================================================================
    SUBROUTINE reset(self, lam, nu)
      class(spectrum_type) :: self
      real(8), intent(in), dimension(:), optional :: lam, nu
      real(8), allocatable, dimension(:) :: tmp_lam, tmp_nu
      !-----------------------------------------------------------------------------
      ! 1. Capture inputs into local temporaries immediately
      if (present(lam)) allocate(tmp_lam, source=lam)
      if (present(nu)) allocate(tmp_nu, source=nu)

      if (allocated(self%lam)) deallocate(self%lam)
      if (allocated(self%nu)) deallocate(self%nu)
      
      if (present(lam)) then
        allocate(self%lam, source=tmp_lam)
        allocate(self%nu, source=HZ_TO_AA / self%lam)
      endif
      if (present(nu)) then
        allocate(self%nu, source=tmp_nu)
        allocate(self%lam, source=HZ_TO_AA / self%nu)
      endif
      self%nnu = size(self%nu)
    END SUBROUTINE reset

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
          weight(j) = 0.5 * ( nu(j+1) - nu(j-1) )
      end do
      weight(n) = 0.5 * ( nu(n) - nu(n-1) )
      weight = abs(weight)
    END SUBROUTINE set_weights

    !===============================================================================
    SUBROUTINE get_interval (self, trans, ired, iblue)
      class(spectrum_type), intent(in):: self
      class(transition_type), intent(in) :: trans
      integer, intent(out) :: ired, iblue
      !-----------------------------------------------------------------------------
      if (self%sorted) then
        call get_ivl(self%lam, trans%lam_min, trans%lam_max, iblue, ired, limit=.false.)
        if (iblue > 1) then
          if((trans%lam_min - self%lam(iblue-1)) < lam_tolerance) then
            iblue = iblue - 1
          endif
        endif
        if (ired < self%nnu) then
          if((self%lam(ired+1) - trans%lam_max) < lam_tolerance) then
            ired = ired + 1
          endif
        endif
      else
        ired = trans%ired
        iblue = trans%iblue
      endif
    END SUBROUTINE get_interval

    !===============================================================================
    SUBROUTINE sort (self)
      class(spectrum_type), intent(inout):: self
      !-----------------------------------------------------------------------------
      self%sorted = .true.
      call refsor(self%lam)
      ! discard duplicates, nnu contains number of unique elements on exit
      call remove_close_entries(self%lam, lam_tolerance)
      call self%reset(lam=self%lam)
    END SUBROUTINE

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
        write(lunit) self%lam
        close(lunit)
      endif
    END SUBROUTINE save

  END MODULE spectrum_mod
  !---------------------------------------------------------------------
  