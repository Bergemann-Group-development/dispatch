module composition_mod
    use prec_mod
    use op_math_mod
    use element_mod
    use ion_mod
    use level_mod
    use transition_mod
    use omp_lib
    implicit none

    type :: composition_type
        logical :: initialised
        integer :: nel, nlevels, nmolec
        real(dp) :: xmy, suma, sumh, grph, summy, tmolim, molh
        type(element_type),    allocatable :: elements(:)
        type(level_type),      allocatable :: levels(:)
        type(molec_type),      allocatable :: molecules(:)
        type(transition_type), allocatable :: trans(:)
        real(dp), dimension(:,:,:,:), allocatable :: nhyd
        real(dp), dimension(:,:,:), allocatable :: temp, ne, rho
        real(dp), dimension(:,:,:), allocatable :: EK_over_temp, pe, dxi, fe, toth, saha
        real(dp), dimension(:,:,:), allocatable :: teta, eh, pg, E, xfakh, ro
        real(dp), dimension(:,:,:,:), allocatable :: abslte, sumsca
    contains
        procedure :: init => injon
        procedure :: init_molecules
        procedure :: load_state_variables
        procedure :: molecular_equilibrium
        procedure :: molecule_by_id
        procedure :: element_by_id
        procedure :: trans_by_id
        procedure :: ion_by_id
        procedure :: absorption
        procedure :: scattering
        procedure :: reassign_pointer
        procedure :: get_totn
        procedure :: get_nhyd
        procedure :: init_mol_dissociation
    end type

contains

    subroutine injon(self, absdat_file, abund_file)
        class(composition_type), target :: self
        character(len=80), intent(in) :: absdat_file, abund_file
        integer :: absdat, abund_unit
        integer :: i, j, k, m, nl, nk
        integer :: lvl_count, skip_back
        real(dp) :: abund, xion, g2, dump, tecor, decor
        character(len=3) :: id
        character(len=65) :: text
        integer :: nqfix, nqtemp, ifish, ncor
        real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
        integer :: nkomp,nextl,nutzl,nextt,nutzt,null
        integer :: ilogl, kvadl, minex, maxex, nnu
        integer :: ilogt, kvadt, minet, maxet, ntt, iteta
        real(dp), dimension(:), allocatable :: tdiff

        character(len=3), parameter :: & 
            element_ids(28) = (/'H  ', 'HE ', 'LI ', 'BE ', 'B  ', 'C  ', 'N  ', 'O  ', &
                                'F  ', 'NE ', 'NA ', 'MG ', 'AL ', 'SI ', 'P  ', 'S  ', &
                                'CL ', 'AR ', 'K  ', 'CA ', 'SC ', 'TI ', 'V  ', 'CR ', &
                                'MN ', 'FE ', 'CO ', 'NI '/)

        type(element_type),    pointer :: elements(:)
        type(element_type),    pointer :: element, H
        type(molec_type),      pointer :: molecule
        type(ion_type),        pointer :: ion, Hm, H0
        type(level_type),      pointer :: lvl
        type(transition_type), pointer :: trans

        call trace_begin ('composition_type%init')

        106 FORMAT(4(F10.3,F10.4))
        107 FORMAT(I10,5I5)

        open(file=absdat_file, newunit=absdat)

        read(absdat, *) self%nel
        allocate(self%elements(self%nel))

        elements => self%elements

        read(absdat, *) elements%id
        read(absdat, *) elements%awgt
        read(absdat, *) elements%nion

        do i = 1, self%nel
            element => elements(i)
            element%np = findloc(element_ids, element%id(:3), dim=1)

            if (i>1) then
                allocate(element%ions(element%nion))
            else
                allocate(element%ions(0:element%nion))
            endif
        end do

        open(file=abund_file, newunit=abund_unit)

        do i=1, self%nel
            read(abund_unit, *) id, abund
            elements(findloc(elements%id, id))%abund = abund
        end do

        close(abund_unit)

        if (elements(1)%abund == 1) then
            elements%abund = log10(elements%abund) + 12.0
        endif

        elements%abund = 10**(elements%abund - elements(1)%abund)

        ! XMY = GRAMS OF STELLAR MATTER / GRAMS OF HYDROGEN
        ! SUMH = NUMBER OF OTHER NUCLEI / NUMBER OF HYDROGEN NUCLEI.
        ! SUMM = NUMBER OF NUCLEI OTHER THAN H, C, N, O / NUMBER OF HYDROGEN NUCLEI
        self%xmy  = sum(elements%abund*elements%awgt) / elements(1)%awgt
        self%suma = sum(elements%abund)
        self%sumh = self%suma - 1
        !self%summ = sum(elements%abund) - sum(elements([1, 3, 4, 5])%abund)

        self%grph  = XMH*self%xmy
        self%summy = self%xmy / self%suma * 1.008


        ! First, we walk through the ions stored in absdat in order to get
        ! total number of energy levels (including fake levels for PFs)

        self%nlevels = 0
        do i=1, self%nel
            element => elements(i)

            do j=1, elements(i)%nion
                ion => element%ions(j)
                ion%j = j
                !IQFIX(I,J) = 2 MEANS THAT A 'FULL' PARTITION FUNCTION SHOULD BE COMPUTED
                ion%iqfix = 2

                !g0 = THE STATISTICAL WEIGHT OF THE GROUND LEVEL
                !nk = THE NUMBER OF ELECTRON CONFIGURATIONS CONSIDERED.

                read(absdat, *) ion%g0, nk
                ion%nk = nk

                skip_back = 0
                do k=1, nk
                    read(absdat, *) dump, dump, dump, nl
                    skip_back = skip_back + 1 + ceiling(nl / 4.0)

                    if (nl > 0) then
                        ion%nk = ion%nk + nl
                        read(absdat, 106) (xion, g2, m=1,nl)
                    endif

                end do

                do k=1, skip_back
                    backspace absdat
                end do

                allocate(ion%levels(ion%nk))
                ion%levels%xl=0

                lvl_count = 0
                do k=1, nk
                    lvl_count = lvl_count + 1
                    lvl => ion%levels(lvl_count)

                    !xion = THE IONIZATION ENERGY IN ELECTRON VOLTS
                    !g2   = THE STATISTICAL WEIGHT (2L+1)*(2J+1)
                    !xl   = THE LOWEST QUANTUM NUMBER OF THE ASYMPTOTIC
                    !       (HYDROGENIC) PARTOF THE PARTITION FUNCTION,
                    !nl   = THE NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSION
                    !       FOR THE 'MIDDLE PART' OF THE PARTITION FUNCTION ('QPRIME').

                    read(absdat, *) lvl%xion, lvl%g2, lvl%xl, nl

                    if (nl > 0) then
                        read(absdat, 106) (ion%levels(m)%xion, ion%levels(m)%g2, m=lvl_count+1, lvl_count+nl)
                        lvl_count = lvl_count + nl
                    endif
                end do

            end do
        end do

        ! Now we rewind the file and this time we actually store the
        ! information on the energy levels

        do i=1, self%nel
            element => elements(i)
            read(absdat, 107) element%active, (element%ions(j)%active, j=1, element%nion)
        end do

        read(absdat, *) nqfix ! number of ions with constant partition function

        do m=1, nqfix
            read(absdat, *) i, j, elements(i)%ions(j)%parco

            !IQFIX(I,J) = 0 MEANS THAT THE PARTITION FUNCTION IS CONSIDERED TO BE CONSTANT.
            elements(i)%ions(j)%iqfix = 0
        end do

        read(absdat, *) nqtemp ! number of ions with PF(temp) instead of PF(temp, pe)

        do m=1, nqtemp
            print *, 'tabulated partition functions to be implemented...'
            stop
        end do

        ! IFISH IS A PARAMETER FOR THE CHOICE OF THE ASYMPTOTIC PARTITION
        ! FUNCTION. IFISH=0 MEANS THAT THE ASYMPTOTIC PART WILL BE EVALU-
        ! ATED FOLLOWING BASCHEK ET AL., ABH. HAMB. VIII,26 (1966). IFISH
        ! =1 MEANS THAT IT WILL BE EVALUATED FOLLOWING FISCHEL AND SPARKS
        ! ASTROPHYS. J. 164, 356 (1971).

        read(absdat, *) ifish

        ! WHEN MOLH IS GREATER THAN ZERO THE MOLECULAR FORMATION WILL BE COMPUTED
        ! IN SUBR. MOLEQ (ONLY H2 AND H2+), ELSE MORE COMPLETE MOLECULAR
        ! FORMATION WILL BE EVALUATED IN SUBR. MOL.
        read(absdat, *) self%tmolim, self%molh

        read(absdat, *) ! skip comment

        ! READS CORONAL ION BALANCE PARAMETERS
        read(absdat, *) ncor, tecor, decor

        do m=1, ncor
            read(absdat, *) i, j, acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
            element => elements(findloc(elements%np, i, dim=1))

            if (j <= element%nion) then
                ion => element%ions(j)

                ion%acolc = acolc
                ion%tcolc = tcolc
                ion%aradc = aradc
                ion%adic  = adic
                ion%bdic  = bdic
                ion%t0c   = t0c
                ion%t1c   = t1c
            endif
        end do

        do i=1, self%nel
            element => elements(i)
            do j=lbound(element%ions, 1), ubound(element%ions, 1)
                ion => element%ions(j)

                if (j == 0) then
                    ion%id = element%id(:3) // '-'
                elseif (j < 4) then
                    ion%id = element%id(:3) // repeat('I', j)
                elseif (j == 4) then
                    ion%id = element%id(:3) // 'IV'
                elseif (j == 5) then
                    ion%id = element%id(:3) // 'V'
                endif

                if(j < ubound(element%ions, 1)) ion%above => element%ions(j+1)
                if(j > 1) ion%below => element%ions(j-1)
                ion%awgt = element%awgt
                ion%abund => element%abund
            end do
        end do

        call self%element_by_id('H', H)
        Hm => H%ions(0) ! H- is ion number 0
        H0 => H%ions(1)

        Hm%iqfix = 0
        Hm%parco = 1
        Hm%nk = 1
        Hm%j = -2
        allocate(Hm%levels(1))
        Hm%levels(1)%xion = 0.754

        call self%init_molecules()
        call self%init_mol_dissociation()

        101 FORMAT(7X,I3,5(9X,I1))
        999 FORMAT(4(9X,I1),8X,I2,9X,I1)

        997 FORMAT(6(F10.0)) ! lam format
        998 FORMAT(6(F10.5)) ! xkap format
        !102 FORMAT(6F10.4)
        !105 FORMAT(F6.3,F4.0,F5.1,I5)

        read(absdat,101) nkomp,nextl,nutzl,nextt,nutzt,null

        allocate(self%trans(nkomp))

        do i=1, nkomp
            trans => self%trans(i)
            trans%bound_free = .false.

            read(absdat,'(A65)') text
            trans%label = text(:20)
            
            if (text(1:1) == '/') then
                trans%key = text(2:6)
                id = text(2:3)
                call self%element_by_id(id, element)
                read(absdat, *) trans%j, xion, trans%level%g2
                trans%level%xion = xion * hc_e
                trans%ion => element%ions(trans%j)
                trans%bound_free = .true.
                trans%regular = .true.
            else
                trans%regular = .false.
                trans%key = text(1:5)
                
                select case(trim(trans%key))
                case('H- BF')
                    trans%bound_free = .true.
                    trans%ion => Hm
                case('H BF')
                    trans%bound_free = .true.
                    trans%ion => elements(1)%ions(1)
                    ! hydrogen energy level
                    read (text(8:10),*) trans%k
                    trans%level%xion = Ryd * (1. - 1./trans%k**2)
                    trans%level%g2 = 2 * trans%k**2
                case('H- FF', 'H FF', '(H+H)', 'H2+')
                    trans%ion => H0
                case('H2-')
                    call self%molecule_by_id('H2', molecule)
                    trans%ion => molecule
                case('HE-')
                    call self%element_by_id('HE', element)
                    trans%ion => element%ions(1)
                end select
            endif

            read(absdat,999) ilogl,kvadl,minex,maxex,nnu

            trans%nnu = nnu
            trans%ilogl = ilogl > 0
            trans%kvadl = kvadl > 0
            trans%minex = minex > 0
            trans%maxex = maxex > 0

            allocate(trans%lam(nnu))
            allocate(trans%ldisk(nnu))

            trans%ldisk = .false.

            read(absdat, *) trans%lam

            trans%minl = trans%lam(1)
            trans%maxl = trans%lam(nnu)

            do j=1, nnu-1
                if (trans%lam(j+1) - trans%lam(j) < 2) then
                    trans%ldisk(j) = .true.
                endif
            end do

            read(absdat,999) ilogt,kvadt,minet,maxet,ntt,iteta

            trans%ntt = ntt
            trans%ilogt = ilogt > 0
            trans%kvadt = kvadt > 0
            trans%minet = minet > 0
            trans%maxet = maxet > 0
            trans%iteta = iteta > 0

            if (ntt > 1) then
                allocate(trans%temp(ntt))
                allocate(trans%tt  (ntt))
                read(absdat, *) trans%temp

                if (iteta > 0) then
                    allocate(trans%teta(ntt))
                    trans%teta = trans%temp
                    trans%temp = temp_to_theta / trans%teta(ntt:1:-1)
                    trans%tt = trans%teta
                else
                    trans%tt = trans%temp
                endif

                trans%mint = trans%temp(1)
                trans%maxt = trans%temp(ntt)
            endif

            allocate(trans%xkap(nnu, ntt))

            do j=1, ntt
                read(absdat, '(6E10.3)') (trans%xkap(k,j),k=1,nnu)
            end do
            
            if (trans%ilogl .or. trans%ilogt) then
                trans%log_xkap = trans%xkap
                where (trans%log_xkap <= 0) trans%log_xkap = 1d-37
                trans%log_xkap = log(trans%log_xkap)
            endif
        end do

        close(absdat)

        self%initialised = .true.
        call trace_end ()
    end subroutine injon

    subroutine reassign_pointer(self)
        class(composition_type), target :: self
        type(element_type),     pointer :: element
        type(molec_type),       pointer :: molecule
        type(ion_type),         pointer :: ion
        type(level_type),       pointer :: lvl
        type(transition_type),  pointer :: trans
        integer :: i, j, m

        call trace_begin ('composition_type%reassign_pointer')
        
        do i=1, self%nel
            element => self%elements(i)

            do j=lbound(element%ions,1), ubound(element%ions,1)
                ion => element%ions(j)
                ion%abund => element%abund

                if (j > lbound(element%ions,1)) then
                    ion%below => element%ions(j-1)
                elseif(j < ubound(element%ions,1)) then
                    ion%above => element%ions(j+1)
                endif

                do m=1, ion%nmolec
                    call self%molecule_by_id(ion%molecules(m)%p%id, molecule)
                    ion%molecules(m)%p => molecule
                end do
            end do
        end do

        do i=1, size(self%molecules)
            molecule => self%molecules(i)

            do j=1, molecule%nions
                ion => molecule%ions(j)%p
                call self%ion_by_id(ion%id, ion)
                molecule%ions(j)%p => ion
            end do
        end do

        do i=1, size(self%trans)
            trans => self%trans(i)
            select type(ptr => trans%ion)
            class is (ion_type)
                ion => ptr
                call self%ion_by_id(ptr%id, ion)
                trans%ion => ion
            class is (molec_type)
                molecule => ptr
                call self%molecule_by_id(ptr%id, molecule)
                trans%ion => molecule
            end select
        end do
        call trace_end ()
    end subroutine reassign_pointer

    subroutine init_molecules(self)
        class(composition_type), target :: self
        type(molec_type),   pointer :: molecule
        type(element_type), pointer :: element
        type(ion_type),     pointer :: ion
        character(len=3) :: mol_ids(12)
        character(len=6) :: ion_id
        integer :: i, j

        call trace_begin ('composition_type%init_molecules')

        mol_ids = ['H2 ', 'H2p', 'H2O', 'OH ', 'CH ', 'CO ', 'CN ', 'C2 ', 'N2 ', 'O2 ', 'NO ', 'NH ']

        self%nmolec = size(mol_ids)

        allocate(self%molecules(self%nmolec))

        do i=1, size(mol_ids)
            molecule => self%molecules(i)
            molecule%nions = len(trim(mol_ids(i)))
            molecule%id = trim(mol_ids(i))

            allocate(molecule%ions(molecule%nions))

            do j=1, molecule%nions
                ion_id = molecule%id(j:j)

                select case (ion_id(:3))
                case ('H', 'C', 'N', 'O')
                    ion_id(4:4) = 'I'
                    call self%ion_by_id(ion_id, ion)
                case ('2')
                    molecule%ions(j)%p => ion ! use previous ion a second time
                    ion%nmolec = ion%nmolec + 1 
                    cycle
                case ('p')
                    molecule%nions = molecule%nions - 1
                    molecule%ions(j-1)%p => ion%above
                    cycle
                case default 
                    cycle
                end select

                molecule%ions(j)%p => ion
                ion%nmolec = ion%nmolec + 1 
            end do
        end do

        do i=1, self%nel
            element => self%elements(i)
            do j=1, element%nion
                ion => element%ions(j)
                allocate(ion%molecules(ion%nmolec))
                ion%mcount = 1
            end do
        end do

        do i=1, size(mol_ids)
            molecule => self%molecules(i)

            do j=1, molecule%nions
                ion_id = molecule%id(j:j)

                select case (ion_id)
                case ('H', 'C', 'N', 'O')
                    ion_id(4:4) = 'I'
                    call self%ion_by_id(ion_id, ion)
                case ('2')
                    ion => ion ! use previous ion a second time
                case ('p')
                    ion => ion%above
                case default 
                    cycle
                end select
                
                ion%molecules(ion%mcount)%p => molecule
                ion%mcount = ion%mcount + 1 
            end do
        end do
        call trace_end ()
    end subroutine

    subroutine init_mol_dissociation (self)
        class(composition_type), target :: self
        type(molec_type),   pointer :: molecule
        integer :: i

        do i=1, size(self%molecules)
            molecule => self%molecules(i)

            select case (trim(molecule%id))
            case ('H2', 'H2p')
                allocate(molecule%de(5))
            case default
                allocate(molecule%de(1))
            end select

            select case (trim(molecule%id))
            case ('H2')
                molecule%dc = [ 6.3021E-4, -1.4149E-2,  1.2572E-1, -5.1172E0,  12.739E0 ]
                molecule%de = [  0.00975,  -0.12427,   0.60602,   -1.4772,    2.6757 ]
            case ('H2p')
                molecule%dc = [ 0.0000e+00,  0.024790, -7.9196803E-2, -2.7942767E0,  11.206998E0 ]
                molecule%de = [  0.15253,  -0.82685,    1.7231,   -2.0036,    2.9216 ]
            case ('H2O')
                molecule%dc = [ 8.1730E-04, -1.8368E-02,  1.6939E-1, -10.522E0,  25.420E0 ]
                molecule%de = [       9.50 ]
            case ('OH')
                molecule%dc = [ 7.7224E-04, -1.6547E-02,  1.3822E-01, -5.0578E0,  12.371E0 ]
                molecule%de = [       4.38 ]
            case ('CH')
                molecule%dc = [ 7.2661E-04, -1.5473E-02,  1.2768E-01, -4.0760E0,  12.135E0 ]
                molecule%de = [       3.47 ]
            case ('CO')
                molecule%dc = [ 1.1349E-03, -2.2888E-02,  1.7217E-01, -11.795E0,  13.820E0 ]
                molecule%de = [      11.11 ]
            case ('CN')
                molecule%dc = [ 3.4666E-04, -7.3627E-03,  6.4162E-02, -8.2793E0,  12.805E0 ]
                molecule%de = [       7.90 ]
            case ('C2')
                molecule%dc = [ 6.2603E-04, -1.2739E-02,  9.7719E-02, -6.5178E0,  12.804E0 ]
                molecule%de = [       6.12 ]
            case ('N2')
                molecule%dc = [ 1.4993E-03, -2.9997E-02,  2.2067E-01, -10.585E0,  13.590E0 ]
                molecule%de = [       9.76 ]
            case ('O2')
                molecule%dc = [ 3.7970E-04, -8.1511E-03,  6.9935E-02, -5.5181E0,  13.228E0 ]
                molecule%de = [       5.12 ]
            case ('NO')
                molecule%dc = [ 1.1380E-03, -2.3065E-02,  1.7349E-01, -7.1964E0,  12.831E0 ]
                molecule%de = [       6.51 ]
            case ('NH')
                molecule%dc = [ 7.8691E-04, -1.6643E-02,  1.3629E-01, -3.8435E0,  12.033E0 ]
                molecule%de = [       3.21 ]
            end select
        end do
    end subroutine init_mol_dissociation

    subroutine load_state_variables(self, temp, ne, rho)
        class(composition_type), target :: self
        type(element_type), pointer :: element
        type(ion_type), pointer :: ion
        type(level_type), pointer :: ground_below
        type(element_type), pointer :: H
        type(ion_type), pointer :: Hm, H0
        real(dp), dimension(:,:,:), intent(in) :: temp, ne, rho
        real(dp), dimension(:,:,:), allocatable :: xnen, ejon, xnenh, xnecno
        integer :: i,j

        call trace_begin ('composition_type%load_state_variables')
        allocate(xnenh, xnecno, ejon, mold=temp)
        xnenh = 0
        xnecno = 0
        ejon = 0

        if (.not. self%initialised) then
            print *, 'composition has not been initialised yet'
            stop
        endif

        self%temp = temp
        self%ne   = ne
        self%rho  = rho
        self%EK_over_temp = EK/self%temp
        self%pe   = self%temp*BK*self%ne
        self%toth = self%rho / self%grph
        self%saha = 2 * op_saha_const * self%temp*sqrt(self%temp) / self%ne
        self%teta = temp_to_theta / self%temp

        ! -4.98e-4 = log10(XMH * op_saha_const)*(temp_to_theta*sqrt(BK))
        self%dxi = 4.98e-4*self%teta*sqrt(self%pe)

        do i=1, self%nel
            element => self%elements(i)

            do j=1, element%nion
                ion => element%ions(j)
                allocate(ion%f, ion%saha, ion%part, ion%anjon, mold=temp)
            end do

            call element%compute_ion_equilibrium(self%EK_over_temp, self%ne, self%saha, self%dxi)
            
            do j=1, element%nion
                ion => element%ions(j)
                ion%f = ion%anjon * element%abund

                ! number of electrons from ion
                xnen = ion%f  * (j-1)
                
                ! number of electrons from elements other than H
                if (i > 1) xnenh = xnenh + xnen
                
                ! number of electrons from elements C, N, O
                if (any(['C', 'N', 'O'] == trim(element%id))) then
                    xnecno = xnecno + xnen
                endif

                if (j > 1) then
                    ground_below => element%ions(j-1)%levels(1)
                    ejon = ejon + ion%f * (ground_below%xion - self%dxi*(j-1))
                endif
            end do
        end do

        ! number of electrons from elements other than H,C,N,O
        xnecno = xnenh - xnecno

        call self%element_by_id('H', H)
        Hm => H%ions(0)
        H0 => H%ions(1)

        Hm%saha = Hm%saha_eq(self%EK_over_temp, self%ne, self%saha, self%dxi, partition_ratio=2.d0) ! ratio H0/H-

        call self%molecular_equilibrium(self%temp, self%teta, self%pe, self%dxi, xnecno, xnenh, self%eh, self%pg, self%fe)

        do i=1, self%nmolec
            !self%molecules(i)%n = self%molecules(i)%f * self%ne / self%fe
            self%molecules(i)%rho_per_h = self%molecules(i)%f / self%grph
        end do

        do i=1, self%nel
            element => self%elements(i)

            element%ions(1)%rho_per_h = element%ions(1)%f / self%grph

            do j=2, element%nion
                ion => element%ions(j)
                ion%rho_per_h = ion%below%rho_per_h * ion%below%saha
            end do
        end do

        allocate (Hm%rho_per_h, self%xfakh, mold=temp)
        Hm%rho_per_h = H0%rho_per_h / Hm%saha

        where (temp < self%tmolim)
            ! keep this around for consistency with legacy code, its just a bunch of wrong constants.
            self%xfakh = self%ne*(BK/uop_BK) * H0%f / self%fe / self%rho * (1.008*UU / XMH)
        elsewhere
            self%xfakh = H0%rho_per_h
        end where

        !self%E = 1.5*self%pg/self%rho+(self%eh+ejon)*EE/self%grph

        ! legacy rho, should probably be replaced by actual rho
        self%ro = self%ne*(BK/uop_BK) * self%grph / self%fe
        call trace_end ()
    end subroutine load_state_variables

    subroutine absorption(self, nu, chi)
        !use phys_mod, only:planck
        class(composition_type), target :: self
        type(transition_type), pointer :: trans
        type(element_type), pointer :: H
        type(molec_type), pointer :: H2
        type(ion_type), pointer :: H0
        real(dp), dimension(:),     intent(in) :: nu
        real(4), dimension(:,:,:,:), intent(inout) :: chi
        real(dp), dimension(:,:,:), allocatable :: fakt, abso, hrest, abfakray, ab
        real(dp), dimension(:,:,:), allocatable :: umc, h1
        real(dp), dimension(:,:,:), allocatable :: fakray, h2ray, els, rayh, rayh2
        integer, dimension(:,:,:), allocatable :: nniv
        real(dp), dimension(:), allocatable :: lam
        real(dp) :: xray2
        integer :: i, f, nf, nx, ny, nz

        call trace_begin ('composition_type%absorption')

        nf = size(nu, 1)
        nx = size(self%temp, 1)
        ny = size(self%temp, 2)
        nz = size(self%temp, 3)

        lam = cc / nu * 1e8

        allocate(self%abslte(nx,ny,nz,nf))
        allocate(self%sumsca(nx,ny,nz,nf))
        allocate(hrest(nx,ny,nz), fakt(nx,ny,nz))

        self%abslte = 0
        self%sumsca = 0 ! line list not implemented yet

        do i=1, size(self%trans)
            trans => self%trans(i)
            select case (trim(trans%key))
            case ('H BF')
                call trans%add_absorption(self%abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%ro, self%dxi, self%xfakh, fakt)
            case ('H FF')
                call trans%add_absorption(self%abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%ro, self%dxi, self%xfakh, fakt)

                hrest = 0
                h1 = sqrt(13.595/self%dxi)
                nniv = min(int(h1), 15)
                umc = 2.3026 * self%dxi * self%teta
                
                where (h1 > 15.5)
                    hrest = fakt * (exp(31.30364*self%teta/((float(nniv)+0.5)**2)-umc)-1.)
                end where

                do f=1, nf
                    self%abslte(:,:,:,f) = self%abslte(:,:,:,f) + hrest
                end do
            end select
        end do

        do f=1, nf
            ! all continua above are multiplied by lam**3
            !------------------------------------------------------
            self%abslte(:,:,:,f) = self%abslte(:,:,:,f) * lam(f)**3
        end do

        do i=1, size(self%trans)
            trans => self%trans(i)

            select case (trim(trans%key))
            case ('H FF', 'H BF', 'H- FF')
                ! H BF and H FF have already been added
                ! H- FF still to be added after multiplication of stim
            case default
                call trans%add_absorption(self%abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%ro, self%dxi, self%xfakh, fakt)
            end select
        end do

        do f=1, nf
            !stim(:,:,:, f) = 1 - exp(- h_e * nu(f) * self%EK_over_temp)
            self%abslte(:,:,:, f) = self%abslte(:,:,:, f) * (1 - exp(-28556.*self%teta/(CC/nu(f)*1e8)))
        end do

        ! all continua above are multiplied by stim
        !--------------------------------------------------------------------------------

        ! H- FF absorption is added without stim factor
        call self%trans_by_id('H- FF', trans)
        call trans%add_absorption(self%abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                  self%pe, self%ro, self%dxi, self%xfakh, fakt)

        ! electron scattering
        els = 4.8206e-9 * self%pe / (self%temp * self%ro)

        ! rayleigh scattering
        call self%molecule_by_id('H2', H2)
        h2ray  = H2%rho_per_h

        call self%element_by_id('H', H)
        H0 => H%ions(1)
        fakray = 2 / H0%part * H0%rho_per_h

        do f=1, nf
            self%sumsca(:,:,:,f) = self%sumsca(:,:,:,f) + els

            if (lam(f) >= 911) then
                if (lam(f) < 1026) then
                    xray2 = 1 / 1026**2
                else
                    xray2 = 1 / lam(f)**2
                endif

                rayh  = xray2*xray2*(5.799e-13 + xray2*(1.422e-6 + xray2*2.784)) * fakray
                rayh2 = xray2*xray2*( 8.14e-13 + xray2*( 1.28e-6 + xray2*1.61 )) * h2ray
    
                ! add self%sumsca to self%sumsca from lines
                self%sumsca(:,:,:,f) = self%sumsca(:,:,:,f) + rayh + rayh2
            endif
        end do

        ! BACKGROUND OPACITY
        
        do f=1, nf
            chi(:,:,:,f) = real((self%abslte(:,:,:,f) + self%sumsca(:,:,:,f)) * self%rho)
        end do
        call trace_end ()
    end subroutine absorption

    subroutine scattering (self, nu, Jnu, eta)
        class(composition_type), target :: self
        real(dp), dimension(:),     intent(in) :: nu
        real(4), dimension(:,:,:,:), intent(in) :: Jnu
        real(4), dimension(:,:,:,:), intent(inout) :: eta
        integer :: f, nf

        nf = size(nu, 1)

        if(size(Jnu, 4) /= nf) then
            print *, 'nu and Jnu do not have the same number of frequencies. Stopping...'
            stop
        endif

        do f=1, nf
            eta(:,:,:,f) = real((self%abslte(:,:,:,f)*planck(nu(f), self%temp) + self%sumsca(:,:,:,f)*Jnu(:,:,:,f)) * self%rho)
        end do
    end subroutine scattering

    subroutine get_totn(self, id, abund, totn)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        real(dp), intent(in) :: abund
        real(dp), dimension(:,:,:), allocatable :: totn
        type(element_type), pointer :: element

        call trace_begin ('composition_type%get_totn')
        call self%element_by_id(id, element)

        if (associated(element)) then
            totn = element%ions(1)%f * (1 + element%ions(1)%saha) * self%toth
        else
            totn = 10.0**(abund-12.0) * self%rho / (self%xmy*1.008*uu)
        endif
        call trace_end ()
    end subroutine get_totn

    subroutine get_nhyd(self, nhyd)
        class(composition_type), target :: self
        real(dp), dimension(:,:,:,:), intent(out) :: nhyd
        real(dp), dimension(:,:,:), pointer :: nh, blzfac
        class(ion_type), pointer :: H0
        type(level_type) :: level
        integer :: k, nk

        call trace_begin ('composition_type%get_nhyd')
        nk = size(nhyd, 4)
        H0 => self%elements(1)%ions(1)

        allocate(nh, blzfac, mold=nhyd(:,:,:,1))

        !PRINT *, omp_get_thread_num(), allocated(self%fe), allocated(self%ne)
        !PRINT *, omp_get_thread_num(), shape(self%fe), 'fe'
        !PRINT *, omp_get_thread_num(), shape(self%ne), 'ne'
        !if (.not.allocated(self%ne)) then
        !  PRINT *, 'self%fe not allocated on thread', omp_get_thread_num()
        !end if
        nh = H0%f / self%fe * self%ne * (BK/uop_BK)

        do k=1, nk
            level%xion = Ryd * (1. - 1./k**2)
            blzfac = level%boltzmann(self%EK_over_temp, natural_log=.true.)
            nhyd(:,:,:,k) = k**2 * blzfac * nh
        end do

        deallocate (nh, blzfac)
        call trace_end ()
    end subroutine

    subroutine molecular_equilibrium(self, temp, teta, pe, dxi, xnecno, xnenh, eh, pg, fe)
        class(composition_type), target :: self
        real(dp), dimension(:,:,:), intent(in) :: temp, pe, teta, dxi, xnecno, xnenh
        real(dp), dimension(:,:,:), allocatable, intent(out) :: eh, pg, fe
        real(dp), dimension(:,:,:), allocatable :: r, s, b, pefe, p, q
        type(molec_type), pointer :: molecules(:)
        type(molec_type), pointer :: H2, CN, H2p
        type(ion_pointer_t) :: ions(4)
        type(ion_type), pointer :: H, C, N, O, ion, ioni, ionj
        type(ion_type), pointer :: Hm
        integer :: i, j, nx, ny, nz, ix, iy, iz, np1, iter, nions
        real(dp) :: new_diff
        real(dp), dimension(:,:,:), allocatable :: diff
        real(dp), dimension(:,:,:,:)  , allocatable, target :: Fvec
        real(dp), dimension(:,:,:,:,:), allocatable, target :: Amat
        real(dp), dimension(:,:,:), pointer :: Fi
        real(dp), dimension(:), allocatable :: Fold, Fnew
        logical, dimension(:,:,:), allocatable :: tmolim

        call trace_begin ('composition_type%molecular_equilibrium')

        nx = size(temp, 1)
        ny = size(temp, 2)
        nz = size(temp, 3)

        call self%ion_by_id('H  I', H)
        call self%ion_by_id('C  I', C)
        call self%ion_by_id('N  I', N)
        call self%ion_by_id('O  I', O)
        Hm => H%below

        molecules => self%molecules

        call self%molecule_by_id('H2' , H2)
        call self%molecule_by_id('CN' , CN)
        call self%molecule_by_id('H2p', H2p)
        
        ions(1)%p => H
        ions(2)%p => C
        ions(3)%p => N
        ions(4)%p => O
        nions = size(ions)
        np1 = nions + 1

        allocate(r, s, eh, pg, fe, diff, q, p, b, mold=temp)
        allocate(Fvec(np1, nx, ny, nz))
        allocate(Amat(np1, np1, nx, ny, nz))
        allocate(Fold(np1), Fnew(np1))
        allocate(H%fil, O%fil, C%fil, N%fil, Hm%pk, mold=temp)

        do i=1, self%nmolec
            allocate (molecules(i)%f, molecules(i)%e, molecules(i)%dis, mold=temp)
            molecules(i)%dis = molecules(i)%pk(teta)
        end do

        tmolim = temp < self%tmolim

        where (tmolim)
            Hm%pk = 1/Hm%saha
            H%fil = 1 + H%saha + Hm%pk
            fe = xnecno
            q = - fe/pe / H2%dis
            p = - H%fil * q
            H%f = solve_quadratic(p, q)

            ! number of free electrons donated by ions over toth
            ! fe = ne / toth
            fe = H%f * (H%saha - Hm%pk) + xnecno

            where (fe <= 0)
                r = H2%dis * Hm%saha * Hm%saha
                s = r - H%fil*Hm%saha
                p = (H%fil*Hm%saha * xnecno - 1 - 2*xnecno*r) / s
                q = r*xnecno*xnecno/s
                fe = solve_quadratic(p, q, avoid_negative=.true.)
                H%f = (xnecno - fe)*Hm%saha
            end where
        elsewhere
            fe = xnenh + self%elements(1)%ions(2)%anjon
        end where

        pefe = pe/fe

        do i=1, self%nmolec
            !partial electron pressures due to each molecule
            where (tmolim) molecules(i)%f = molecules(i)%calc_f(teta, pefe)
        end do

        where (tmolim)
            b = C%f_of_molecules_with(O) / C%f
            O%fil = 1 + O%saha + O%f_of_molecules_with(H) / O%f
            C%fil = 1 + C%saha + C%f_of_molecules_with(H) / C%f

            p = (C%abund-O%abund)/O%fil + C%fil/b
            q = - O%abund/O%fil * C%fil/b
            O%f = solve_quadratic(p, q)
            C%f = C%abund / (C%fil + b*O%f)

            CN%f = CN%calc_f(teta, pefe) ! has to be recalculated since C%f changed
            
            N%fil = 1 + N%saha + (N%f_of_molecules_with(H) + N%f_of_molecules_with(C)) / N%f
            b = N%f / N%f_of_molecules_with(N, count_double=.true.) 
            p = 2*N%fil*b
            q = - 2*N%abund*b
            N%f = solve_quadratic(p, q)
            where (N%f < 0) N%f = 0
        end where

        diff = 0

        do iter=1, 100
            Amat = 0
            do i=1, self%nmolec
                where (tmolim) molecules(i)%f = molecules(i)%calc_f(teta, pefe)
            end do

            do i=1, nions
                ioni => ions(i)%p
                Fi  => Fvec(i,:,:,:)

                where (tmolim) 
                    Fi = ioni%f_of_own_molecules()
                    Amat(i, np1,:,:,:) = - Fi / fe
                end where

                if (i==1) then
                    where (tmolim) 
                        Amat(i,  i,:,:,:) = 1 + H%saha + Hm%pk
                        Amat(np1,i,:,:,:) = H%saha + 2*H2p%f /(H%f * pefe) - Hm%pk
                        Fvec(np1,:,:,:) = H%f * (H%saha + H2p%f - Hm%pk) - fe + xnecno
                    end where
                else
                    where (tmolim) 
                        Amat(i,  i,:,:,:) = 1 + ioni%saha
                        Amat(np1,i,:,:,:) = ioni%saha
                        Fvec(np1,:,:,:) = Fvec(np1,:,:,:) + ioni%saha * ioni%f
                    end where
                endif

                where (tmolim) Fi = Amat(i,i,:,:,:)*ioni%f + Fi - ioni%abund

                do j=1, nions
                    ionj => ions(j)%p
                    where (tmolim) Amat(i,j,:,:,:) = Amat(i,j,:,:,:) + ionj%f_of_molecules_with(ioni, count_double=.true.) / ionj%f
                end do
            end do

            where (tmolim) Amat(np1, np1,:,:,:) = - H2p%f/pe - 1

            do iz=1, nz
            do iy=1, ny
            do ix=1, nx
                if (.not. tmolim(ix,iy,iz)) cycle

                if (iter > 2 .and. diff(ix,iy,iz) < 0.001) then
                    cycle
                endif

                Fold(1:nions) = [(ions(i)%p%f(ix, iy, iz), i=1, nions)]
                Fold(np1) = fe(ix,iy,iz)

                call AINV2(Amat(:,:,ix,iy,iz),np1)
                Fnew = Fold - matmul(Amat(:,:,ix,iy,iz),Fvec(:,ix,iy,iz))

                do i=1, nions
                    ions(i)%p%f(ix, iy, iz) = Fnew(i)
                end do

                fe(ix,iy,iz) = Fnew(np1)
                pefe(ix,iy,iz) = pe(ix,iy,iz) / fe(ix,iy,iz)

                diff(ix,iy,iz) = maxval(abs(1-Fold/Fnew))
            end do
            end do
            end do

            if (iter > 2 .and. all(diff < 0.001)) exit
        end do

        !print *, iter, 'iterations in molecular equilibrium calculations'
        where (tmolim)
            !eh = - H%f*((H%levels(1)%xion - dxi)*(1+Hm%pk) + (Hm%levels(1)%xion - 2*dxi)*Hm%pk)
            pg = fe + self%suma !sum(self%elements%abund)
        elsewhere
            !eh = (H%levels(1)%xion - dxi) * H%above%anjon
            pg = pe * (1 + self%suma/fe)
        end where

        do i=1, self%nmolec
            where (tmolim)
                molecules(i)%f = molecules(i)%calc_f(teta, pefe) ! has to be called first
                pg = pg + molecules(i)%f
                !molecules(i)%e = molecules(i)%calc_e(teta, temp, dxi)
                !eh = eh + molecules(i)%e
            end where
        end do

        do i=1, nions
            ion => ions(i)%p
            where (tmolim) pg = pg + ion%f*(1+ion%saha) - ion%abund
        end do

        where (tmolim) pg = pg * pefe

        call trace_end ()
    end subroutine molecular_equilibrium

    subroutine molecule_by_id(self, id, molecule)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        type(molec_type), pointer, intent(out) :: molecule

        call trace_begin ('composition_type%init')
        molecule => self%molecules(findloc(self%molecules%id, id, 1))
        call trace_end ()
    end subroutine molecule_by_id

    subroutine element_by_id(self, id, element)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        type(element_type), pointer, intent(out) :: element
        integer :: idx

        call trace_begin ('composition_type%init')
        idx = findloc(self%elements%id, uppercase(id, len(id)), 1)

        if (idx == 0) then
            !print *, 'Element of id= ', id, ' is not part of the opacity package.'
            element => null()
        else
            element => self%elements(idx)
        endif
        call trace_end ()
    end subroutine element_by_id

    subroutine ion_by_id(self, id, ion)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        type(ion_type), pointer, intent(out) :: ion
        type(element_type), pointer :: element
        integer :: idx

        call trace_begin ('composition_type%init')
        call self%element_by_id(id(1:3), element)
        idx = findloc(element%ions%id, id, 1)
        if (lbound(element%ions,1) /= 1) then
            idx = idx + lbound(element%ions,1) - 1
        endif
        ion => element%ions(idx)
        call trace_end ()
    end subroutine ion_by_id

    subroutine trans_by_id(self, id, trans)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        type(transition_type), pointer, intent(out) :: trans
        integer :: i, idx

        call trace_begin ('composition_type%init')
        do i=1, size(self%trans)
            if (trim(self%trans(i)%key) == id) then
                trans => self%trans(i)
                call trace_end ()
                return
            endif
        end do

        print *, 'transition key = ', id, ' not found. Stopping...'
        stop
    end subroutine trans_by_id

    !---------------------------------------------------------------------
    ! PE: 22-09-21: Convert to elemental function
    real(8) elemental function planck(nu,t) result(bnu)
        !---------------------------------------------------------------------
        !
        !  calculates planck function bnu at frequency nu, temp t
        !
        !---------------------------------------------------------------------
        use mconst
        implicit none
        real(8), intent(in)  :: nu
        real(8), intent(in)  :: t
        real(8) :: x, c1,c2

        c1=2.0*hh*nu**3/cc**2
        c2=hh*nu/bk

        x=c2/t

        if(x.lt.80.) then
            bnu=c1/(exp(x)-1.0)
        else
            bnu=c1*exp(-x)
        endif
    end function planck
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    pure function uppercase(string,ss)
        ! converts the first ss characters of string to uppercase and returns them as a new string.
        character(len=*),intent(in) :: string
        integer,intent(in) :: ss
        character(len=ss) :: uppercase
        integer :: i,k
        character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
        character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'

        uppercase(:)=' '
        uppercase(1:ss)=string(1:ss)
        do i=1,ss
        k=index(lcase,string(i:i))
        if (k > 0) then
            uppercase(i:i)=ucase(k:k)
        else
            uppercase(i:i)=string(i:i)
        end if
        end do

    end function uppercase
    !---------------------------------------------------------------------

end module composition_mod
