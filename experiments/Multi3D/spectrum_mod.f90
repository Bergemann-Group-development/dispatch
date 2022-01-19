!---------------------------------------------------------------------
module spectrum_mod 
    implicit none
    !---------------------------------------------------------------------
    !
    ! initialises the spectrum structure
    !
    !---------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    type spectrum_type
        ! total number of frequencies
        integer :: nnu
        ! maximum number of active lines or continua at any given frequency
        integer :: maxal,maxac
        ! frequency [Hz] and frequency integration weights
        real(8),dimension(:),allocatable :: nu,wnu
        ! which continua (ac) and lines (al) are active at each frequency
        integer,dimension(:,:),  allocatable :: ac,al
        ! how many continua (nac) and lines (nal) are active at each frequency
        integer,dimension(:),  allocatable :: nac,nal
        ! has this frequency an active prd line
        logical,dimension(:),allocatable :: hasPRDline
    contains
        procedure :: init => initialise_spectrum
    end type spectrum_type
    !---------------------------------------------------------------------

    type(spectrum_type),target :: sp

  contains
  
    !---------------------------------------------------------------------
    subroutine initialise_spectrum(sp, line, cont, qnorm)
      use m3io_mod,only:stop
      use mmath
      use trans_mod
      implicit none

      class(spectrum_type) :: sp
      class(line_type), intent(inout), dimension(:), target :: line
      class(cont_type), intent(inout), dimension(:), target :: cont
      real(8), intent(in) :: qnorm
      type(line_type),pointer :: lin
      type(cont_type),pointer :: cnt
      integer :: i,j,nnu,ired,iblue,js,je
      character(len=13),parameter ::rname='initialise_spectrum'
      real(8),dimension(:),allocatable :: nu
  
      ! work space for tautsp
      integer :: iflag,nnutmp,k,impi_error, nline, ncont
      !real(8),dimension(:,:),allocatable :: work,coef
      !real(8),dimension(:)  ,allocatable :: break
      real(8),dimension(:)  ,allocatable :: dalpha_dnu,alpha_fromfile,nu_fromfile
  
      nline = size(line, 1)
      ncont = size(cont, 1)

      ! compute total number of frequencies
      nnu=0
      do i=1,nline
         nnu=nnu+line(i)%nnu
      end do
      do i=1,ncont
         nnu=nnu+cont(i)%nnu
      end do
  
      ! save all frequencies
      allocate(nu(nnu))
      je=0
      do i=1,nline
         js=je+1
         je=js+line(i)%nnu-1
         nu(js:je)=line(i)%nu
      end do
      do i=1,ncont
         js=je+1
         je=js+cont(i)%nnu-1
         nu(js:je)=cont(i)%nu
      end do
  
      ! sort frequencies into ascending order 
      call refsor(nu)
   
     ! discard duplicates, nnu contains number of unique elements on exit
      call unista(nu,nnu)
   
     ! now set overall frequency grid
      sp%nnu=nnu
      allocate(sp%nu(nnu)) 
      sp%nu=nu(1:nnu)

      ! and compute weights
      call set_weights(sp%wnu, sp%nu, 1, nnu)
     
      ! allocate list of active continua and lines
      allocate(sp%ac(ncont,sp%nnu)) ; sp%ac=0
      allocate(sp%nac(sp%nnu)) ; sp%nac=0
  
      allocate(sp%al(nline,sp%nnu)) ; sp%al=0
      allocate(sp%nal(sp%nnu)) ; sp%nal=0
  
      !allocate(sp%hasPRDline(sp%nnu)) ; sp%hasPRDline=.false.
      
      do i=1,ncont
  
         cnt => cont(i)
  
         ! find minimum and maximum index of this continuum in the sp%nu grid
         call locate(sp%nnu,sp%nu,cnt%nu(1),ired)
         call locate(sp%nnu,sp%nu,cnt%nu(cnt%nnu),iblue)
  
         ! total number of frequencies where this continuum is active
         nnu=iblue-ired+1
  
         ! min and max frequency where this continuum is active
         cnt%ired  = ired
         cnt%iblue = iblue
       
         ! add to number of active continua
         do j=ired,iblue
            sp%nac(j) = sp%nac(j) + 1 ! number of active continua at frequency index j
            sp%ac( sp%nac(j), j ) = i ! the index of this active continuum
         end do
         
         ! recompute the bf cross sections
         select case(cnt%bf_type)
  
         case('hydrogenic')
            ! use nu**-3 dependency
            cnt%nnu   = nnu
  
            deallocate( cnt%nu      )
            allocate(   cnt%nu(ired:iblue) )
            cnt%nu=sp%nu(ired:iblue)
  
            deallocate( cnt%alpha      )
            allocate(   cnt%alpha(ired:iblue) )
            call freqc(cnt)
  
         case('detail')
            
            ! interpolate from given values
            if (allocated(dalpha_dnu)) deallocate(dalpha_dnu)
            allocate( dalpha_dnu(cnt%nnu) )
            call hermite_derivative(dalpha_dnu, cnt%nu, cnt%alpha)
            
            if (allocated(alpha_fromfile)) deallocate(alpha_fromfile)
            allocate( alpha_fromfile(cnt%nnu) )
            alpha_fromfile = cnt%alpha
  
            if (allocated(nu_fromfile)) deallocate(nu_fromfile)
            allocate( nu_fromfile(cnt%nnu) )
            nu_fromfile = cnt%nu
            
            ! set new frequency grid        
            cnt%nnu   = nnu
            deallocate( cnt%nu      )
            allocate(   cnt%nu(ired:iblue) )
            cnt%nu=sp%nu(ired:iblue)
         
            ! allocate space for new cross sections
            deallocate( cnt%alpha             )
            allocate(   cnt%alpha(ired:iblue) )
  
            ! fill new cross sections 
            call hermite_interpolation(cnt%alpha, cnt%nu, nu_fromfile, alpha_fromfile, dalpha_dnu)
  
         case default
  
            call stop('initialise_spectrum : cont(i)%bf_type' &
                 // cnt%bf_type // 'not recognized')
  
         end select

         call set_weights(cnt%wnu, cnt%nu, ired, iblue)
      end do
  
      ired  =  sp%nnu/2
      iblue =  sp%nnu/2
  
  
      do i=1,nline
  
         lin => line(i)
  
         ! find minimum and maximum index of this line in the sp%nu grid
         call locate(sp%nnu,sp%nu,lin%nu(1),ired)
         call locate(sp%nnu,sp%nu,lin%nu(lin%nnu),iblue)
  
         ! total number of frequencies where this line is active
         nnu=iblue-ired+1
  
         ! min and max frequency where this line is active
         lin%ired  = ired
         lin%iblue = iblue
       
         ! add to number of active lines
         do j=ired,iblue
            sp%nal(j) = sp%nal(j) + 1
            sp%al( sp%nal(j), j ) = i
         end do
  
         ! set new frequency grid for the line      
         lin%nnu   = nnu
         deallocate( lin%nu      )
         allocate(   lin%nu(ired:iblue) )
         lin%nu=sp%nu(ired:iblue) ! [Hz]
         allocate(   lin%q(ired:iblue) )
         lin%q = (lin%nu/lin%nu0-1.0d0) * cc / (qnorm * KM_TO_CM) ! [Doppler units]
  
         ! fill frequency integration weights for the line
         call set_weights(lin%wnu, lin%nu, ired, iblue)
         call set_weights(lin%wq , lin%q , ired, iblue)
  
         ! flag if PRD line at this frequency
         !if (lin%profile_type == 'prd') then
         !   do j=ired,iblue
         !      sp%hasPRDline(j) = .true.
         !   end do
         !end if
  
      end do
  
      sp%maxac=maxval(sp%nac)
      sp%maxal=maxval(sp%nal)

      write(*,'(a,I0,a)') 'initialise_spectrum: frequency grid has ',sp%nnu,' points.'
  
    end subroutine initialise_spectrum
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    subroutine set_weights (weight, nu, start, end)
        real(8), dimension(:), allocatable :: weight
        real(8), dimension(:), intent(in) :: nu
        integer, intent(in) :: start, end
        integer :: j, n

        n = size(nu, 1)

        if (allocated(weight)) deallocate(weight)
        allocate( weight(start:end) )

        weight(start) = 0.5 * ( nu(2) - nu(1) )

        do  j=2,n-1
            weight(start - 1 + j) = 0.5* ( nu(j+1) - nu(j-1) )
        end do

        weight(end) = 0.5 * ( nu(n) - nu(n-1) )
    end subroutine set_weights
    !---------------------------------------------------------------------
  
  
  end module spectrum_mod
  !---------------------------------------------------------------------
  