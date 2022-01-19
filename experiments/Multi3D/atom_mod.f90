!---------------------------------------------------------------------
module atom_mod

    use m3io_mod
    use mconst
    use atom_var, only: atom_type, line_type, cont_type
    implicit none

contains

    !---------------------------------------------------------------------
    subroutine read_atom (atom, atom_file, qnorm, abnd)
        !---------------------------------------------------------------------
        !
        !  reads an atomic model and calculates einstein a,b coefficients.
        !
        !---------------------------------------------------------------------
        use mclu
        use coll_mod, only: get_ncoll, read_collisions
        !use mem_mod
        type(atom_type), intent(out), target :: atom
        character(len=*), intent(in) :: atom_file
        real(8), intent(in) :: qnorm
        real(8), intent(in), optional :: abnd

        ! Declaration of the local variables
        integer             :: i, j, io, latom, strl, kr, krs, iterm, nterm
        integer             :: ncoll, max_n_temp
        character(len = 80) :: text, mulfmt
        character(len = 9), parameter :: rname='read_atom'
        character(len = :), allocatable :: myfmt ! temporary format string
        type(line_type), pointer :: lin => null()
        type(cont_type), pointer :: cnt => null()

        open(newunit=latom, file=atom_file, form='formatted')

        call nextentry(latom)
        read(latom,fmt='(a)') atom%id
        atom%id=lowercase(adjustl(atom%id),len_trim(atom%id))

        call nextentry(latom)
        read(latom,*) atom%abnd, atom%awgt
        atom%awgt=atom%awgt*uu ! atomic weight [g]
        ! AJG: quickly change abundance without having to edit atom
        if (present(abnd)) atom%abnd = abnd

        call nextentry(latom)
        read(latom,*) atom%nlevel,atom%nline,atom%ncont,atom%nrfix
        atom%nrad = atom%nline + atom%ncont

        atom%qnorm = qnorm

        allocate(  atom%line  ( atom%nline ) ) 
        allocate(  atom%cont  ( atom%ncont ) )

        allocate(  atom%ev    ( atom%nlevel) ) ; atom%ev = 0.0d0
        allocate(  atom%g     ( atom%nlevel) ) ; atom%g = 0
        allocate(  atom%label ( atom%nlevel) ) ; atom%label = ' '
        allocate(  atom%ion   ( atom%nlevel) ) ; atom%ion = 0

        allocate(  atom%ilin  ( atom%nlevel,atom%nlevel) ) ; atom%ilin=0
        allocate(  atom%icon  ( atom%nlevel,atom%nlevel) ) ; atom%icon=0

        do i=1, atom%nlevel
            call nextentry(latom)
            read(latom,*) atom%ev(i),atom%g(i), atom%label(i), atom%ion(i)
            atom%ev(i)=atom%ev(i)*hc_e  ! [eV]  
        end do

        ! sanity check on level energies
        do i = 2, atom%nlevel
            if(atom%ev(i) < atom%ev(i - 1)) then
            myfmt = "('atom: levels are: ', a, 2x, a, ' energy not monotonically &
                &increasing.')"
            write(*, myfmt) atom%label(i), atom%label(i - 1)
            call stop(' ')
            endif
        end do

        do i=1,atom%nline
            lin => atom%line(i)

            call nextentry(latom) 

            ! io WAS present for historical reasons, in MULTI io determines
            ! whether a line is wide or not. In the RH formulation all lines are wide

            read(latom,*) lin%j   , lin%i , lin%f , &
                lin%nnu, lin%qmax, lin%q0, lin%io, &
                lin%ga , lin%gw  , lin%gq, lin%profile_type

            ! io is now modified so that multicomponent lines can be read
            ! io = 0: Normal line
            ! io = 1: For now is the same as io = 0
            ! io > 1: Number of components in term
            ! Like in Multi, thereafter follows one extra input line per component:

            !   determ  wterm  fterm  gaterm  gwterm  gqterm
            ! determ    energy difference for component i of transition kr. unit cm-1.
            ! wterm     weight of component kterm in weighted average
            ! fterm     oscillator strength
            ! gaterm    damping constants
            ! gwterm
            ! gqterm
            ! adterm    adamp for component, calculated in damp

            nterm = lin%io

            if (nterm .ge. 2) then

                mulfmt = ' '
                write(mulfmt, fmt = '(a)') uppercase(rname, len_trim(rname)) &
                    // ' found ' // string2(nterm) // &
                    'components in transition KR = ' // string3(i)
                call messg(mulfmt)

                allocate(lin%determ (nterm) )
                allocate(lin%wterm  (nterm) )
                allocate(lin%fterm  (nterm) )
                allocate(lin%gaterm (nterm) )
                allocate(lin%gwterm (nterm) )
                allocate(lin%gqterm (nterm) )
                !allocate(lin%adterm, xs, xe, ys, ye, zs1, ze1, 1, nterm, rname)

                do iterm = 1, nterm
                
                read(latom,*) lin%determ(iterm)          , lin%wterm(iterm) , &
                        lin%fterm(iterm) , lin%gaterm(iterm), lin%gwterm(iterm), &
                        lin%gqterm(iterm)

                end do

            end if

            lin%lambda0 = hc_e / ( atom%ev(lin%j) - atom%ev(lin%i) ) ! [cm]
            lin%nu0   = cc / lin%lambda0    ! [Hz]
            lin%Grat  = atom%g(lin%i) / atom%g(lin%j)
            lin%Aji   = f_to_Aji * lin%f * lin%Grat / lin%lambda0**2
            lin%Bji   = c2_2h / lin%nu0**3 * lin%Aji 
            lin%Bij   = lin%Bji / lin%Grat
            lin%profile_type = &
                adjustl(lowercase(lin%profile_type,len_trim(lin%profile_type)))

            lin%ired=1
            lin%iblue=lin%nnu

            lin%ntrans=i
            atom%ilin(lin%j,lin%i) = i
            atom%ilin(lin%i,lin%j) = i

            if ( (trim(lin%profile_type) .ne. 'voigt'  ) .and. &
                (trim(lin%profile_type) .ne. 'prd'    ) .and. &
                (trim(lin%profile_type) .ne. 'doppler') .and. &
                (trim(lin%profile_type) .ne. 'hlinop' )         ) then 

                if ( atom%id .eq. 'h' ) then

                write (*,fmt = '(a, x, a, a)') 'Profile not recognised:', &
                        uppercase(lin%profile_type, len_trim(lin%profile_type)), &
                        ' Setting it to HLINOP'

                lin%profile_type = adjustl(lowercase('hlinop', len_trim('hlinop')))
                
                else

                write (*, fmt = '(a, x, a, a)') 'Profile not recognised:', &
                        uppercase(lin%profile_type, len_trim(lin%profile_type)), &
                        ' Setting it to VOIGT'

                lin%profile_type = adjustl(lowercase('voigt', len_trim('voigt')))
                
                end if

            end if

            ! make sure we always have odd number of frequency points for a line
            if (mod(lin%nnu,2).eq.0) lin%nnu=lin%nnu+1
            allocate( lin%nu(  lin%nnu ) )

        end do

        !  bound-free transitions in detail
        !  if numax.lt.0.0 frequency points in angstrom (starting with
        !  threshold and decreasing) and crossections in cm2 are read.
        do i=1,atom%ncont
            cnt => atom%cont(i)
            call nextentry(latom) 

            read(latom,*) cnt%j,cnt%i,cnt%alpha0,cnt%nnu,cnt%numax
        
            cnt%nu0   = ( atom%ev(cnt%j) - atom%ev(cnt%i) ) / h_e ! [Hz]
            cnt%numax = cc/(cnt%numax*AA_TO_CM) ! [from AA to Hz]
            cnt%ntrans=atom%nline+i ! transition identifier for continua

            cnt%ired=1
            cnt%iblue=cnt%nnu
            atom%icon( cnt%j, cnt%i ) = i
            atom%icon( cnt%i, cnt%j ) = i

            allocate( cnt%nu(1:cnt%nnu)    )
            allocate( cnt%alpha(1:cnt%nnu) )

            if ( cnt%numax .lt. 0.0 ) then

                ! cross sections given in detail
                cnt%bf_type='detail'
                do j=1,cnt%nnu
                call nextentry(latom) 
                read(latom,*)  cnt%nu(j), cnt%alpha(j)
                end do
                ! sanity check
                do j=2,cnt%nnu
                if ( cnt%nu(j) .gt. cnt%nu(j-1) ) then
                    write(text,*) 'read_atom: continuum transition nr', i, ' wavelengths not decreasing'
                    call stop(text)
                endif
                end do

                ! and convert to frequency [Hz]
                cnt%nu = cc / ( cnt%nu * AA_TO_CM )
                cnt%numax = cnt%nu(cnt%nnu)

            else

                ! hydrogenic transitions
                cnt%bf_type='hydrogenic'

            endif

        end do

        ! read choice of collisional subroutine
        ! GENCOL and HCOL are currently the only collisions subroutines.
        call nextentry(latom)
        read(latom,fmt='(a)') atom%crout

        ! sanity check
        atom%crout = adjustl(atom%crout)
        strl=len_trim(atom%crout)
        atom%crout = lowercase(atom%crout,strl)

        ! read collisional data
        call get_ncoll(latom, atom%crout, ncoll, max_n_temp)
        allocate( atom%coll(ncoll) )
        call read_collisions(latom, atom%coll, max_n_temp, atom%nlevel)

        close(latom)

        ! frequency grid in lines
        do i=1,atom%nline
            call freql(atom%line(i), qnorm)
        end do

        !  frequency grid in continua
        do i=1,atom%ncont
            call freqc(atom%cont(i))
        end do
    end subroutine read_atom
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    subroutine freql(line, qnorm)
        !***********************************************************************
    
        !  gives  frequency quadrature points and corresponding weights.
        !  one- or two-sided cases.  the points and weights are based on
        !  the trapezoidal rule applied to a mapping function  that maps
        !  equidistant  points in  x to  equidistant points in  q  in the
        !  core (q .lt. q0), and equidistant points in log(q) in the wings
        !  (q .gt. q0).
    
        !  nq     : total number of quadrature points                (in)
        !  qmax   : maximum value of frequency (doppler units)       (in)
        !  ind    : =1 for one-sided and =2 for two-sided cases      (in)
        !  q0     : transition from core to wing (doppler units)     (in)
        !  q      : quadrature points                               (out)
        !  wq     : quadrature weights                              (out)
    
        !  coded by: ake nordlund (oct-1981).
    
        !  revised by mats carlsson (jan-1984): new options
        !  qmax.le.q0   q linear in frequency
        !  qmax or q0.lt.0  q given and only weights are calculated
    
        !:
        !: freql  88-05-03  modifications: (mats carlsson)
        !:        added options gave incorrect weights for ind=1
        !:        this affected cooling rates for lines with frequency
        !:        quadrature given with qmax=q0 or qmax.lt.0 or q0.lt.0
        !:        the corrected weights are twice the old ones
        !:
        !***********************************************************************
    
        class(line_type), intent(inout)  :: line
        real(8) :: qnorm
        real(8),dimension(:),allocatable :: q,wq
    
        real(8)  :: dq, x, x10, al10, a, half, xmax, dx,qmax,q0
        integer :: nnu, j
    
    
        q0=line%q0
        qmax=line%qmax
        nnu=line%nnu
        allocate(q(nnu))
        allocate(wq(nnu))
    
        if (qmax.le.q0) then
        ! linear spacing
        dq = 2 * qmax / (nnu-1)
        q(1)=-qmax
        do j=2,nnu
            q(j) = q(j-1) + dq
        end do
    
        else if(qmax.ge.0.0 .and. q0.ge.0.0) then
    
        ! constants
        al10=log(10.)
        a=10.**(q0+.5)
        half=0.5
        xmax=log10(a*max(half,qmax-q0-half))
    
        dx=2.*xmax/(nnu-1)
        do j=1,nnu
            x=-xmax+(j-1)*dx
            x10=10.**x
            q(j)=x+(x10-1./x10)/a
        end do
    
        end if
    
        ! convert to [Hz]
        line%nu = line%nu0*(1.0+q*qnorm*KM_TO_CM/cc)
    
    end subroutine freql
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    subroutine freqc(cnt)
        !---------------------------------------------------------------------   
        !
        !  gives quadrature points and weights for continua
        !  the points are equidistant in frequency, the weights are
        !  trapezoidal
        !  calculates crossections from nn**-3 dependence if not given
        !  explicitly
        !
        !  nnu      total number of quadrature points    
        !  nu     quadrature points (Hz)                
        !  wnu      quadrature weights                  
        !
        !---------------------------------------------------------------------

        class(cont_type), intent(inout) :: cnt
        
        integer :: j,nnu,ired,iblue
        real(8) :: dnu

        select case(cnt%bf_type)
        case('hydrogenic')
            nnu   = cnt%nnu
            ired  = cnt%ired
            iblue = cnt%iblue

            ! equidistant grid in nu, between nu0 and numax 
            cnt%nu(ired)  = cnt%nu0
            cnt%nu(iblue) = cnt%numax

            dnu = (cnt%numax-cnt%nu0) / (nnu-1)

            do j=ired,iblue
                cnt%nu(j)    = cnt%nu0 + (j-ired) * dnu
                cnt%alpha(j) = cnt%alpha0 * ( cnt%nu(j) / cnt%nu0 )**(-3)
            end do
        case('detail')
            ! do nothing, frequencies were computed in read_atom
        case default
            call stop('freqc: bf_type ' // cnt%bf_type // ' not recognized')
        end select
    end subroutine freqc
    !---------------------------------------------------------------------



end module atom_mod
!---------------------------------------------------------------------

