MODULE composition_mod
    use trace_mod
    use prec_mod
    use op_math_mod
    use element_mod
    use ion_mod
    use level_mod
    use bbtrans_mod, only:bbtrans_type
    use bftrans_mod, only:bftrans_type
    use molecule_mod, only:molec_type
    use osmet_mod
    use omp_lib
    use os_mod
    implicit none
    private

    type, public :: composition_type
        logical :: initialised, ll_initialised=.false., loaded=.false.
        integer :: nel, nmolec, verbose
        real(dp) :: xmy, suma, grph, tmolim, molh
        type(element_type), allocatable :: elements(:)
        type(molec_type),   allocatable :: molecules(:)
        type(bftrans_type), allocatable :: trans(:)
        type(osmet_type) :: osmet
        real(dp), dimension(:,:,:,:), allocatable :: nhyd
        real(dp), dimension(:,:,:), allocatable :: temp, ne, rho
        real(dp), dimension(:,:,:), allocatable :: EK_over_temp, pe, dxi, fe, toth, saha
        real(dp), dimension(:,:,:), allocatable :: teta, eh, pg, E, xfakh, ro, p6
    contains
        procedure :: init => injon
        procedure :: init_molecules
        procedure :: load_state_variables
        procedure :: load_state_variables_without_nel
        procedure :: molecular_equilibrium
        procedure :: molecule_by_id
        procedure :: element_by_id
        procedure :: trans_by_id
        procedure :: ion_by_id
        procedure :: absorption
        procedure :: add_hrest
        procedure :: electron_scattering
        procedure :: rayleigh_scattering
        procedure :: scattering
        procedure :: reassign_pointer
        procedure :: get_totn
        procedure :: get_nhyd
        procedure :: compute_grams_per_hydrogen
        procedure :: compute_p6
        procedure :: overwrite_abundance
        procedure :: init_linelist
        procedure :: print_ion_table
    end type

contains

    !===============================================================================
    SUBROUTINE injon(self, absdat_file, abund_file, absmet_file, absmet_big_end)
        class(composition_type), target :: self
        character(len=80), intent(in) :: absdat_file, abund_file
        character(len=80), intent(in), optional :: absmet_file
        logical, intent(in), optional :: absmet_big_end
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
        character(len=3), parameter :: & 
            element_ids(28) = (/'H  ', 'HE ', 'LI ', 'BE ', 'B  ', 'C  ', 'N  ', 'O  ', &
                                'F  ', 'NE ', 'NA ', 'MG ', 'AL ', 'SI ', 'P  ', 'S  ', &
                                'CL ', 'AR ', 'K  ', 'CA ', 'SC ', 'TI ', 'V  ', 'CR ', &
                                'MN ', 'FE ', 'CO ', 'NI '/)

        class(element_type),    pointer :: elements(:)
        class(element_type),    pointer :: element, H
        class(molec_type),      pointer :: molecule
        class(ion_type),        pointer :: ion, Hm, HI
        class(level_type),      pointer :: lvl
        class(bftrans_type), pointer :: trans
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')

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
        call self%compute_grams_per_hydrogen(self%suma, self%xmy, self%grph)

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

        do i=1, self%nel
            element => elements(i)
            read(absdat, 107) element%active, (element%ions(j)%active, j=1, element%nion)
        end do

        read(absdat, *) nqfix ! number of ions with constant partition FUNCTION

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
                    ion%id = trim(element%id(:3)) // ' -'
                elseif (j < 4) then
                    ion%id = trim(element%id(:3)) // ' ' // repeat('I', j)
                elseif (j == 4) then
                    ion%id = trim(element%id(:3)) // ' IV'
                elseif (j == 5) then
                    ion%id = trim(element%id(:3)) // ' V'
                endif

                if(j < ubound(element%ions, 1)) ion%above => element%ions(j+1)
                if(j > 1) ion%below => element%ions(j-1)
                ion%awgt = element%awgt
                ion%abund => element%abund
            end do
        end do

        call self%element_by_id('H', H)
        Hm => H%ions(0) ! H- is ion number 0
        HI => H%ions(1)
        HI%below => Hm

        Hm%iqfix = 0
        Hm%parco = 1
        Hm%nk = 1
        Hm%j = -2
        allocate(Hm%levels(1))
        Hm%levels(1)%xion = 0.754
        ! --- Different Hm ionisation in Oslo legacy code ---!
        !Hm%levels(1)%xion = 0.747
        ! --------------------------------------------------!
        RYD = HI%levels(1)%xion ! overwriting RYD constant
        ! --------------------------------------------------!

        call self%init_molecules()

        101 FORMAT(7X,I3,5(9X,I1))
        999 FORMAT(4(9X,I1),8X,I2,9X,I1)

        997 FORMAT(6(F10.0)) ! lam format
        998 FORMAT(6(F10.5)) ! xkap format
        !102 FORMAT(6F10.4)
        !105 FORMAT(F6.3,F4.0,F5.1,I5)

        read(absdat,101) nkomp,nextl,nutzl,nextt,nutzt,null

        allocate(self%trans(nkomp))
        self%trans%verbose = self%verbose

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
                    trans%level%xion = RYD * (1. - 1./trans%k**2)
                    trans%level%g2 = 2 * trans%k**2
                case('H- FF', 'H FF', '(H+H)', 'H2+')
                    trans%ion => HI
                    trans%level = HI%levels(1)
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
                allocate(trans%log_xkap, mold=trans%xkap)
                trans%log_xkap = trans%xkap
                where (trans%log_xkap <= 0) trans%log_xkap = 1d-37
                trans%log_xkap = log(trans%log_xkap)
            endif
        end do

        close(absdat)

        if (present(absmet_file)) then
            if (absmet_file /= '') then
                call self%osmet%init(absmet_file, absmet_big_end=absmet_big_end)
                self%osmet%eosmet = 0.5
            endif
        endif

        self%initialised = .true.
        call trace%end ()
    END SUBROUTINE injon

    !===============================================================================
    SUBROUTINE init_linelist (self, linelist_file)
        class(composition_type) :: self
        class(ion_type), pointer :: ion
        class(bbtrans_type), dimension(:), allocatable :: lines
        class(bbtrans_type), pointer :: line
        character(len=*) :: linelist_file
        character(20) :: mass_id
        character(20) :: ion_id
        real(dp) :: lam0, gv, gq, loggf, ent8, ent9, mass, cnt_ev
        character :: lower_orb, upper_orb
        character(len=150) :: label
        integer :: lunit, iostat, j, nline, i, prev, ntotal=0
        logical :: isotope
        !-----------------------------------------------------------------------------
        if (linelist_file == '') return

        open(newunit=lunit, file=linelist_file, iostat=iostat)
        do while (iostat==0)
            read(lunit,*,iostat=iostat) mass_id, j, nline
            if (iostat/=0) exit
            read(lunit,*,iostat=iostat) ion_id
            read(mass_id,*) mass
            ion_id = adjustl(trim(ion_id))
            call self%ion_by_id(ion_id(1:6), ion)
            isotope = mod(mass, 1.0) > 0

            !ignore isotopes for now
            if (.not. associated(ion) .or. isotope) then
                do i=1, nline
                    read(lunit, *, iostat=iostat) ! skipping this ion
                    if (iostat.ne.0) exit
                end do
                cycle
            endif

            if (ion%nline == 0) then
                allocate(ion%line(nline))
                ion%nline = nline
                prev = 0
            else
                allocate(lines(nline))
                lines = ion%line
                deallocate(ion%line)
                allocate(ion%line(ion%nline+nline))
                ion%line(1:ion%nline) = lines
                deallocate(lines)
                prev = ion%nline
                ion%nline = ion%nline + nline
            endif

            do i=1, nline
                line => ion%line(prev+i)
                ! Try to read with stark broadening parameter
                read(lunit,*,iostat=iostat) lam0, line%lower_ev, loggf, gv, line%upper_g, &
                    line%ga, gq, lower_orb, upper_orb, ent8, ent9, label

                if (iostat.ne.0) then
                    ! Try to read without stark broadening parameter
                    backspace lunit
                    read(lunit,*,iostat=iostat) lam0, line%lower_ev, loggf, gv, line%upper_g, &
                        line%ga, lower_orb, upper_orb, ent8, ent9, label

                    ! Stop if nothing has worked
                    if (iostat.ne.0) then
                        print *, "Failed to read a line in the linelist", lam0
                        stop
                    endif
                    gq = 0.0
                endif

                lam0 = air2vac(lam0) * 1e-8
                line%nu0 = CC / lam0
                line%Aji = f_to_Aji * 10**loggf / (line%upper_g * lam0**2)
                line%Bji = c2_2h / line%nu0**3 * line%Aji
                line%Bij = line%Bji * line%upper_g
                cnt_ev = ion%above%levels(1)%xion

                call line%prepare_Stark(gq)
                call line%prepare_VanderWaals(gv, j, cnt_ev)

                ntotal = ntotal + 1
            end do
        end do
        close(lunit)
        print *, ntotal, ' lines read from linelist ', linelist_file
        self%ll_initialised = .true.
    END SUBROUTINE init_linelist

    !===============================================================================
    SUBROUTINE reassign_pointer(self)
        class(composition_type), target :: self
        class(element_type),     pointer :: element
        class(molec_type),       pointer :: molecule
        class(ion_type),         pointer :: ion
        class(level_type),       pointer :: lvl
        class(bftrans_type),  pointer :: trans
        integer :: i, j, m
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%reassign_pointer')
        
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
        call trace%end ()
    END SUBROUTINE reassign_pointer

    !===============================================================================
    SUBROUTINE init_molecules(self)
        class(composition_type), target :: self
        class(molec_type),   pointer :: molecule
        class(element_type), pointer :: element
        class(ion_type),     pointer :: ion
        character(len=3) :: mol_ids(12)
        character(len=6) :: ion_id
        integer :: i, j
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init_molecules')

        mol_ids = ['H2 ', 'H2p', 'H2O', 'OH ', 'CH ', 'CO ', 'CN ', 'C2 ', 'N2 ', 'O2 ', 'NO ', 'NH ']

        self%nmolec = size(mol_ids)

        allocate(self%molecules(self%nmolec))

        do i=1, size(mol_ids)
            molecule => self%molecules(i)
            molecule%nions = len(trim(mol_ids(i)))
            molecule%id = trim(mol_ids(i))
            call molecule%init_mol_dissociation()

            allocate(molecule%ions(molecule%nions))

            do j=1, molecule%nions
                ion_id = molecule%id(j:j)

                select case (ion_id(:3))
                case ('H', 'C', 'N', 'O')
                    ion_id(3:3) = 'I'
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
                    ion_id(3:3) = 'I'
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
        call trace%end ()
    END SUBROUTINE init_molecules

    !===============================================================================
    SUBROUTINE load_state_variables(self, temp, ne, rho)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        class(molec_type), pointer :: molecule
        class(ion_type), pointer :: ion
        class(element_type), pointer :: H
        class(ion_type), pointer :: Hm, HI
        real(4), dimension(:,:,:), intent(in) :: temp, ne, rho
        real(dp), dimension(:,:,:), allocatable :: xnen, ejon, xnecno, nh0, xerg, ehm
        integer :: i,j, nx, ny, nz, lunit
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%load_state_variables')

        if (.not. self%initialised) then
            print *, 'composition has not been initialised yet'
            stop
        endif

        nx = size(temp,1)
        ny = size(temp,2)
        nz = size(temp,3)

        if (self%loaded) then
            deallocate(self%temp)
            deallocate(self%xfakh, self%ne,   self%rho,  self%pe, self%pg)
            deallocate(self%toth,  self%saha, self%teta, self%dxi, self%eh)
            deallocate(self%EK_over_temp, self%ro, self%p6, self%fe)
        endif
        
        allocate(self%temp(nx,ny,nz))
        allocate(self%xfakh, self%ne,   self%rho,  self%pe, self%pg, mold=self%temp)
        allocate(self%toth,  self%saha, self%teta, self%dxi, self%eh, mold=self%temp)
        allocate(self%EK_over_temp, self%ro, self%p6, self%fe, mold=self%temp)
        allocate(xnen, xnecno, ejon, nh0, xerg, ehm, mold=self%temp)

        self%temp = temp
        self%ne   = ne
        self%rho  = rho
        self%EK_over_temp = EK / self%temp
        self%pe   = self%temp*BK*self%ne
        self%toth = self%rho / self%grph
        self%saha = 2 * op_saha_const * self%temp*sqrt(self%temp) / self%ne
        self%teta = temp_to_theta / self%temp

        self%dxi = 4.98e-4*self%teta*sqrt(self%pe)
        ! --- this only shows up in the Oslo legacy code ---!
        self%dxi=min(self%dxi,0.747*0.5)
        ! --------------------------------------------------!

        self%fe = 0
        xnecno = 0
        ejon = 0

        do i=1, self%nel
            element => self%elements(i)
            element%nx = nx
            element%ny = ny
            element%nz = nz

            call element%compute_ion_equilibrium(self%EK_over_temp, self%ne, self%saha, self%dxi)
            
            do j=1, element%nion
                ion => element%ions(j)
                ion%f = ion%anjon * element%abund

                ! number of electrons from ion
                xnen = ion%f * (j-1)
                
                ! number of electrons from elements other than HI
                if (ion%id /= 'H I') self%fe = self%fe + xnen
                
                ! number of electrons from elements other than H, C, N, O
                if (all(trim(element%id) /= ['H', 'C', 'N', 'O'])) then
                    xnecno = xnecno + xnen
                endif

                if (j > 1) then
                    xerg = ion%below%levels(1)%xion - self%dxi*(j-1)

                    if (trim(element%id) == 'H') then
                        self%eh = ion%below%f * xerg
                    else
                        ! energy of ionization (ejon). hydrogen is not included.
                        ejon = ejon + ion%f * xerg
                    endif
                endif
            end do
        end do

        call self%element_by_id('H', H)
        Hm => H%ions(0)
        HI => H%ions(1)
        call Hm%allocate(nx,ny,nz)

        Hm%saha = Hm%saha_eq(self%EK_over_temp, self%ne, self%saha, self%dxi, partition_ratio=2.d0) ! ratio HI/H-
        ehm = (Hm%levels(1)%xion + HI%levels(1)%xion) * HI%f / Hm%saha

        self%pg = self%pe * (1 + self%suma/self%fe)

        do i=1, self%nmolec
            molecule => self%molecules(i)
            call molecule%allocate(nx,ny,nz)
            molecule%dis = molecule%pk(self%teta)
        end do

        call self%molecular_equilibrium(self%temp, self%teta, self%pe, self%dxi, xnecno, self%eh, self%pg, self%fe)

        do i=1, self%nel
            element => self%elements(i)
            element%ions(1)%rho_per_h = element%ions(1)%f / self%grph

            do j=2, element%nion
                ion => element%ions(j)
                ion%rho_per_h = ion%below%rho_per_h * ion%below%saha
            end do
        end do

        nh0 = HI%f * self%ne / self%fe * (BK/uop_BK)

        Hm%rho_per_h = HI%rho_per_h / Hm%saha

        where (temp < self%tmolim)
            ! keep this around for consistency with legacy code, its just a bunch of wrong constants.
            !self%xfakh = self%ne*(BK/uop_BK) * HI%f / (self%fe * self%rho) * (1.008*UU / XMH)
            self%xfakh = nh0 / self%rho * (1.008*UU / XMH)
        elsewhere
            self%xfakh = HI%rho_per_h
        end where

        ! legacy rho, should probably be replaced by actual rho
        self%ro = self%ne / self%fe * (BK/uop_BK) * self%grph

        self%E = 1.5*self%pg/self%ro + (ejon - self%eh)*EE/self%grph

        call self%compute_p6(self%p6)
        self%loaded = .true.
        call trace%end ()
    END SUBROUTINE load_state_variables

    !===============================================================================
    SUBROUTINE compute_p6(self, p6)
        class(composition_type), target :: self
        class(ion_type), pointer :: HeI, HI
        class(molec_type), pointer :: H2, H2p
        real(dp), dimension(:,:,:), intent(out) :: p6
        !---------------------------------------------------------------------------
        call self%molecule_by_id('H2',  H2)
        call self%molecule_by_id('H2p', H2p)
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('HE I', HeI)
        
        p6 = LOG10((HI%f + 0.42*HeI%f + 0.85*(H2%f + H2p%f)) * self%toth * BK * self%temp)
    END SUBROUTINE compute_p6

    !===============================================================================
    SUBROUTINE absorption(self, nu, chi, abslte, sumsca, include_absmet, absmet_big_end)
        !use phys_mod, only:planck
        class(composition_type), target :: self
        class(bftrans_type), pointer :: trans
        real(dp), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:,:), intent(inout) :: chi, abslte, sumsca
        real(dp), dimension(:,:,:), allocatable :: stim
        real(dp), dimension(:), allocatable :: lam
        logical, optional :: include_absmet, absmet_big_end
        real(dp) :: xray2, epsln
        real(4) :: linop
        integer :: i, j, k, f, nf
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%absorption')
        allocate(lam, mold=nu)

        nf = size(nu, 1)
        lam = cc / nu * 1e8
        abslte = 0

        do i=1, size(self%trans)
            trans => self%trans(i)
            select case (trim(trans%key))
            case ('H BF', 'H FF')
                call trans%add_absorption(abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%ro, self%dxi, self%xfakh)
            end select
        end do

        call self%add_hrest(abslte)

        do f=1, nf
            ! Hrest, H BF and H FF are multiplied by lam**3
            !------------------------------------------------------
            abslte(:,:,:,f) = abslte(:,:,:,f) * lam(f)**3
        end do

        !print *, 'after lam3', abslte(pp(1),pp(2),pp(3),pp(4)), lam(pp(4))**3

        do i=1, size(self%trans)
            trans => self%trans(i)
            select case (trim(trans%key))
            case ('H FF', 'H BF', 'H- FF')
                ! H BF and H FF have already been added
                ! H- FF still to be added after multiplication of stim
            case default
                call trans%add_absorption(abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%ro, self%dxi, self%xfakh)
            end select
        end do

        allocate(stim, mold=self%temp)
        do f=1, nf
            ! all continua above are multiplied by stim
            !--------------------------------------------------------------------------------
            !stim = 1. - exp(- h_e * nu(f) * self%EK_over_temp)
            stim = 1. - exp(-28556.*self%teta/lam(f))
            abslte(:,:,:, f) = abslte(:,:,:, f) * stim

            !if (f==pp(4)) then
            !    print *, 'stim', stim(pp(1),pp(2),pp(3)), lam(f)
            !endif
        end do

        ! H- FF absorption is added without stim factor
        call self%trans_by_id('H- FF', trans)
        call trans%add_absorption(abslte, lam, self%temp, self%EK_over_temp, self%teta, &
                                  self%pe, self%ro, self%dxi, self%xfakh)

        sumsca = 0
        call self%electron_scattering(sumsca)
        call self%rayleigh_scattering(lam, sumsca)

        if (present(include_absmet)) then
            if (self%osmet%initialised .and. include_absmet) then
                call self%osmet%add_line_opacities(nu, self%p6, self%temp, self%ne, abslte, sumsca)
            endif
        endif

        do f=1, nf
            ! BACKGROUND OPACITY
            chi(:,:,:,f) = real((abslte(:,:,:,f) + sumsca(:,:,:,f)) * self%rho)
        end do

        !print *, 'chi/rho', abslte(pp(1),pp(2),pp(3),pp(4)) + sumsca(pp(1),pp(2),pp(3),pp(4))
        !print *, 'abslte', abslte(pp(1),pp(2),pp(3),pp(4))
        !print *, 'sumsca', sumsca(pp(1),pp(2),pp(3),pp(4))

        call trace%end ()
    END SUBROUTINE absorption

    !===============================================================================
    SUBROUTINE scattering (self, nu, abslte, eta, Jnu, sumsca)
        class(composition_type), target :: self
        real(dp), dimension(:),     intent(in) :: nu
        real(4), dimension(:,:,:,:), intent(in) :: abslte
        real(4), dimension(:,:,:,:), intent(inout) :: eta
        real(4), dimension(:,:,:,:), intent(in), optional :: Jnu, sumsca
        integer :: f, nf
        !---------------------------------------------------------------------------
        nf = size(nu, 1)
        if (present(Jnu) .and. present(sumsca)) then
            do f=1, nf
                eta(:,:,:,f) = real((abslte(:,:,:,f)*planck(nu(f), self%temp) + sumsca(:,:,:,f)*Jnu(:,:,:,f)) * self%rho)
            end do
        else
            do f=1, nf
                eta(:,:,:,f) = real((abslte(:,:,:,f)*planck(nu(f), self%temp)) * self%rho)
            end do
        endif
    END SUBROUTINE scattering

    !===============================================================================
    SUBROUTINE get_totn(self, id, abund, totn)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        real(dp), intent(in) :: abund
        real(4), dimension(:,:,:) :: totn
        class(element_type), pointer :: element
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%get_totn')
        call self%element_by_id(trim(id), element)

        if (associated(element)) then
            totn = element%ions(1)%f * (1 + element%ions(1)%saha) * self%toth
        else
            totn = 10.0**(abund-12.0) * self%rho / (self%xmy*1.008*uu)
        endif
        call trace%end ()
    END SUBROUTINE get_totn

    !===============================================================================
    SUBROUTINE get_nhyd(self, nhyd)
        class(composition_type), target :: self
        real(4), dimension(:,:,:,:), intent(out) :: nhyd
        real(dp), dimension(:,:,:), pointer :: nh, blzfac
        class(ion_type), pointer :: HI
        type(level_type) :: level
        integer :: k, nk
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%get_nhyd')
        nk = size(nhyd, 4)
        HI => self%elements(1)%ions(1)

        allocate(nh, blzfac, mold=self%temp)
        nh = HI%f * self%ne / self%fe * (BK/uop_BK)

        do k=1, nk
            level%xion = RYD * (1. - 1./k**2)
            blzfac = level%boltzmann(self%EK_over_temp, natural_log=.true.)
            nhyd(:,:,:,k) = k**2 * blzfac * nh
        end do

        deallocate (nh, blzfac)
        call trace%end ()
    END SUBROUTINE get_nhyd

    !===============================================================================
    SUBROUTINE molecular_equilibrium(self, temp, teta, pe, dxi, xnecno, eh, pg, fe)
        class(composition_type), target :: self
        real(dp), dimension(:,:,:), intent(in) :: temp, pe, teta, dxi, xnecno
        real(dp), dimension(:,:,:), allocatable, intent(inout) :: eh, pg, fe
        real(dp), dimension(:,:,:), allocatable :: r, s, b, pefe, p, q
        class(molec_type), pointer :: molecules(:)
        class(molec_type), pointer :: H2, CN, H2p
        type(ion_pointer_t) :: ions(4)
        class(ion_type), pointer :: H, C, N, O, ion, ioni, ionj
        class(ion_type), pointer :: Hm
        integer :: i, j, nx, ny, nz, ix, iy, iz, np1, iter, nions
        real(dp) :: new_diff
        real(dp), dimension(:,:,:), allocatable :: diff
        real(dp), dimension(:,:,:,:)  , allocatable, target :: Fvec
        real(dp), dimension(:,:,:,:,:), allocatable, target :: Amat
        real(dp), dimension(:,:,:), pointer :: Fi
        real(dp), dimension(:), allocatable :: Fold, Fnew
        logical, dimension(:,:,:), allocatable :: tmolim
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%molecular_equilibrium')

        nx = size(temp, 1)
        ny = size(temp, 2)
        nz = size(temp, 3)

        call self%ion_by_id('H I', H)
        call self%ion_by_id('C I', C)
        call self%ion_by_id('N I', N)
        call self%ion_by_id('O I', O)
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

        allocate(r, s, diff, q, p, b, pefe, mold=temp)
        allocate(Fvec(np1, nx, ny, nz))
        allocate(Amat(np1, np1, nx, ny, nz))
        allocate(Fold(np1), Fnew(np1))
        allocate(tmolim(nx,ny,nz))

        if (allocated(Hm%pk)) then
            deallocate(Hm%pk, H%fil, O%fil, C%fil, N%fil)
        endif

        allocate(Hm%pk, H%fil, O%fil, C%fil, N%fil, mold=temp)

        tmolim = temp < self%tmolim

        where (tmolim)
            Hm%pk = 1/Hm%saha
            H%fil = 1 + H%saha + Hm%pk
            fe = xnecno
            q = - fe/pe / H2%dis
            p = - H%fil * q
            H%f = solve_quadratic(p, q)

            fe = H%f * (H%saha - Hm%pk) + xnecno

            where (fe <= 0)
                r = H2%dis * Hm%saha * Hm%saha
                s = r - H%fil*Hm%saha
                p = (H%fil*Hm%saha * xnecno - 1 - 2*xnecno*r) / s
                q = r*xnecno*xnecno/s
                fe = solve_quadratic(p, q, avoid_negative=.true.)
                H%f = (xnecno - fe)*Hm%saha
            end where
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
            eh = - H%f*((H%levels(1)%xion - dxi)*(1+Hm%pk) + (Hm%levels(1)%xion - 2*dxi)*Hm%pk)
            pg = fe + self%suma !sum(self%elements%abund)
        end where

        do i=1, self%nmolec
            where (tmolim)
                molecules(i)%f = molecules(i)%calc_f(teta, pefe) ! has to be called first
                molecules(i)%e = molecules(i)%calc_e(teta, temp, dxi)
                pg = pg + molecules(i)%f
                eh = eh + molecules(i)%e
            end where
        end do

        do i=1, nions
            ion => ions(i)%p
            where (tmolim) pg = pg + ion%f*(1+ion%saha) - ion%abund
        end do

        where (tmolim) pg = pg * pefe

        do i=1, self%nmolec
            !self%molecules(i)%n = self%molecules(i)%f * self%ne / self%fe
            self%molecules(i)%rho_per_h = self%molecules(i)%f / self%grph
        end do

        call trace%end ()
    END SUBROUTINE molecular_equilibrium

    !===============================================================================
    SUBROUTINE molecule_by_id(self, id, molecule)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(molec_type), pointer, intent(out) :: molecule
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        molecule => self%molecules(findloc(self%molecules%id, id, 1))
        call trace%end ()
    END SUBROUTINE molecule_by_id

    !===============================================================================
    SUBROUTINE element_by_id(self, id, element)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(element_type), pointer, intent(out) :: element
        integer :: idx
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        idx = findloc(self%elements%id, uppercase(id, len(id)), 1)

        if (idx == 0) then
            element => null()
            print *, 'element with id=' // id // ' was not found.'
        else
            element => self%elements(idx)
        endif
        call trace%end ()
    END SUBROUTINE element_by_id

    !===============================================================================
    SUBROUTINE ion_by_id(self, id, ion)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        character(3) :: eid
        class(ion_type), pointer, intent(out) :: ion
        class(element_type), pointer :: element
        integer :: idx, space
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        space = 0
        do while (id(space+1:space+1) /= ' ')
            space = space+1
        end do
        eid = '   '
        eid(1:space) = id(1:space)
        call self%element_by_id(eid, element)
        if (associated(element)) then
            idx = findloc(element%ions%id, uppercase(id, len(id)), 1)
            if (lbound(element%ions,1) /= 1) then
                idx = idx + lbound(element%ions,1) - 1
            endif
            if (idx < lbound(element%ions,1) .or. idx>ubound(element%ions,1)) then
                ion => null()
                print *, 'ion with id=' // id // ' was not found.'
            else
                ion => element%ions(idx)
            endif
        else
            ion => null()
        endif
        call trace%end ()
    END SUBROUTINE ion_by_id

    !===============================================================================
    SUBROUTINE trans_by_id(self, id, trans)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(bftrans_type), pointer, intent(out) :: trans
        integer :: i, idx
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        do i=1, size(self%trans)
            if (trim(self%trans(i)%key) == id) then
                trans => self%trans(i)
                call trace%end ()
                return
            endif
        end do

        print *, 'transition key = ', id, ' not found. Stopping...'
        stop
    END SUBROUTINE trans_by_id

    !===============================================================================
    ! PE: 22-09-21: Convert to elemental FUNCTION
    real(8) elemental FUNCTION planck(nu,t) result(bnu)
        !---------------------------------------------------------------------
        !
        !  calculates planck FUNCTION bnu at frequency nu, temp t
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
    end FUNCTION planck

    !================================================================================
    pure FUNCTION uppercase(string,ss)
        ! converts the first ss characters of string to uppercase and returns them as a new string.
        character(len=*),intent(in) :: string
        integer,intent(in) :: ss
        character(len=ss) :: uppercase
        integer :: i,k
        character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
        character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
        !---------------------------------------------------------------------------
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

    end FUNCTION uppercase

    !=========================================================================
    SUBROUTINE compute_grams_per_hydrogen(self, suma, xmy, grph)
        class(composition_type), target :: self
        real(dp), intent(out) :: suma, xmy, grph
        !---------------------------------------------------------------------
        ! Abundances should not be in log units when they enter here!
        ! XMY = GRAMS OF STELLAR MATTER / GRAMS OF HYDROGEN
        !---------------------------------------------------------------------
        suma = sum(self%elements%abund)
        xmy  = sum(self%elements%abund*self%elements%awgt) / self%elements(1)%awgt
        grph = XMH*xmy
    END SUBROUTINE compute_grams_per_hydrogen

    !=========================================================================
    SUBROUTINE overwrite_abundance(self, id, abnd, awgt)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        character(len=*), intent(in) :: id
        real(dp), intent(in) :: abnd, awgt
        real(dp) :: abund
        logical, save :: first_time=.true.
        !---------------------------------------------------------------------
        call self%element_by_id(trim(id), element)

        if (associated(element)) then
            abund = log10(element%abund) + 12.0
            if (abund /= abnd) then
                if (first_time) then
                    first_time = .false.
                    print '(A, A, A, F5.2, A, F5.2)', 'Changing abundance of ', trim(id), ' from ', abund, ' to ', abnd
                endif

                element%abund =  10**(abnd - 12.0)
                element%awgt  =  awgt / uu
                call self%compute_grams_per_hydrogen(self%suma, self%xmy, self%grph)
            endif
        else
            print *, 'overwrite_abundance: Element of id= ', id, ' is not part of the opacity package.'
        endif
    END SUBROUTINE overwrite_abundance

    !=========================================================================
    SUBROUTINE load_state_variables_without_nel(self,temp,rho,tol)
        ! find a state variables consistent with given temp & rho
        class(composition_type), target :: self
        real(4), dimension(:,:,:), intent(in) :: temp, rho
        real(4), intent(in) :: tol ! stop cryterium
        logical, dimension(:,:,:), allocatable :: ne_converged
        real(4), dimension(:,:,:), allocatable :: ne, a, b, x, xold, fx, fa, fb 
        integer :: i, nconverded, arr_size, arr_shape(3)
        integer, parameter :: nmax = 20 ! max number of iteration

        !self%temp = temp 
        !self%rho = rho
        
        arr_shape = shape(temp)
        arr_size = size(temp)

        allocate(ne, mold=temp)
        allocate(ne_converged(arr_shape(1),arr_shape(2),arr_shape(3))) !can use mold= on different kind=
        allocate(x,xold,a,b,fx,fa,fb, mold=temp)

        ne_converged(:,:,:) = .false.

        ! set initial value for electron density

        where (temp.gt.1.0d4) ! fully ionized gas
            ne =  rho/self%GRPH
        else where ! neutrals only
            ne = 1e-4*rho/self%GRPH
        end where

        a = log(ne)
        call self%load_state_variables(temp, ne, rho)

        fa = log(self%ro / rho)

        where (abs(fa) < tol)
            ne_converged = .true.
        endwhere

        b = a - 0.69*fa
        ne = exp(b)

        call self%load_state_variables(temp, ne, rho)
        fb = log(self%ro / rho)

        where (abs(fb) < tol)
            ne_converged = .true.
        endwhere

        nconverded = count(ne_converged)
        i = 0
        write(*,*)
        write(*,fmt="(5x,'it [',I2.2,']',1x,'n-converded:',1x,I3,1X,'/ total size: ',I3)") i, nconverded, arr_size

        if (nconverded < arr_size) then
            x=b
            do i=1,nmax
                xold = x
                where (.not.ne_converged)
                    x=a-(b-a)/(fb-fa)*fa
                end where
        
                where (abs(x-xold) < tol)
                    ne_converged = .true.
                end where
                nconverded = count(ne_converged)
                if (size(ne_converged) == count(ne_converged)) exit ! all good !

                where (.not.ne_converged)
                    ne = exp(x)
                endwhere

                call self%load_state_variables(temp, ne, rho)
                fx = log(self%ro / rho)
                where (abs(a-x).lt.abs(b-x)) 
                    b = x 
                    fb = fx
                else where 
                    a = x 
                    fa = fx 
                end where
            enddo
        endif

        write(*,fmt="(5x,'it [',I2.2,']',1x,'n-converded:',1x,I3,1X,'/ total size: ',I3)") i, nconverded, arr_size
        write(*,*)

        deallocate(ne)
        deallocate(ne_converged)
        deallocate(x,xold,a,b,fx,fa,fb)

    END SUBROUTINE load_state_variables_without_nel

    !=========================================================================
    SUBROUTINE add_hrest(self, abslte)
        class(composition_type), target :: self
        class(ion_type), pointer :: HI
        real(4), dimension(:,:,:,:), intent(inout) :: abslte
        real(dp), dimension(:,:,:), allocatable :: hrest, umc, h1, nniv, fakt
        integer :: f
        !---------------------------------------------------------------------
        allocate(hrest, fakt, umc, h1, nniv, mold=self%temp)
        call self%ion_by_id('H I', HI)

        fakt = HI%levels(1)%boltzmann(self%EK_over_temp, self%dxi, natural_log=.true.)
        fakt = 2.0898E-26 / (2*RYD*self%EK_over_temp) * fakt * HI%part * HI%rho_per_h
        h1 = sqrt(RYD/self%dxi)
        nniv = (min(int(h1), 15) + 0.5)**2
        
        hrest = 0
        where (h1 > 15.5)
            !umc = 2.3026 * self%dxi * self%teta
            !hrest = fakt * (exp(31.30364*self%teta/nniv-umc)-1.)
            hrest = exp(self%EK_over_temp*(HI%levels(1)%xion/nniv - self%dxi))
            hrest = fakt * (hrest - 1)
        end where

        do f=1, size(abslte, 4)
            abslte(:,:,:,f) = abslte(:,:,:,f) + hrest
        end do

        !print *, 'hrest', hrest(pp(1),pp(2),pp(3)), abslte(pp(1),pp(2),pp(3),pp(4))
    END SUBROUTINE add_hrest

    !=========================================================================
    SUBROUTINE electron_scattering(self, sumsca)
        class(composition_type), target :: self
        real(4), dimension(:,:,:,:), intent(inout) :: sumsca
        real(dp), dimension(:,:,:), allocatable :: els
        integer :: f
        !---------------------------------------------------------------------
        allocate(els, mold=self%temp)
        els = 4.8206e-9 * self%pe / (self%temp * self%ro)

        do f=1, size(sumsca, 4)
            sumsca(:,:,:,f) = sumsca(:,:,:,f) + els
        end do

        !print *, 'els', els(pp(1),pp(2),pp(3))
    END SUBROUTINE electron_scattering

    !=========================================================================
    SUBROUTINE rayleigh_scattering(self, lam, sumsca)
        class(composition_type), target :: self
        class(molec_type), pointer :: H2
        class(ion_type), pointer :: HI
        real(dp), dimension(:), intent(in) :: lam
        real(4), dimension(:,:,:,:), intent(inout) :: sumsca
        real(dp), dimension(:,:,:), allocatable :: h2ray, rayh, rayh2
        real(dp) :: xray2, ll
        integer :: f
        !---------------------------------------------------------------------
        allocate(h2ray, rayh, rayh2, mold=self%temp)

        ! rayleigh scattering
        call self%ion_by_id('H I', HI)
        call self%molecule_by_id('H2', H2)

        do f=1, size(lam)
            if (lam(f) >= 911) then
                ll = max(1026., lam(f))
                xray2 = 1 / ll**2

                rayh  = (5.799e-13 + xray2*(1.422e-6 + xray2*2.784)) * HI%rho_per_h
                rayh2 = ( 8.14e-13 + xray2*( 1.28e-6 + xray2*1.61 )) * H2%rho_per_h
    
                ! add sumsca to sumsca from lines
                sumsca(:,:,:,f) = sumsca(:,:,:,f) + (rayh + rayh2)*xray2*xray2

                !if (f==pp(4)) then
                !    print *, 'rayh', rayh(pp(1),pp(2),pp(3))*xray2**2, HI%rho_per_h(pp(1),pp(2),pp(3)), xray2
                !    print *, 'rayh2', rayh2(pp(1),pp(2),pp(3))*xray2**2
                !endif
            endif
        end do
    END SUBROUTINE rayleigh_scattering

    !=========================================================================
    SUBROUTINE print_ion_table(self)
        class(composition_type), target :: self
        integer :: i, lunit
        !---------------------------------------------------------------------
        open(newunit=lunit, file='./ion_info.txt')
        write(lunit, '(A8, es12.4)') 'temp=', self%temp(pp(1),pp(2),pp(3))
        write(lunit, '(A8, es12.4)') 'pe=', self%pe(pp(1),pp(2),pp(3))
        write(lunit, '(A8, es12.4)') 'ne=', self%ne(pp(1),pp(2),pp(3))
        write(lunit, '(A8, es12.4)') 'rho=', self%rho(pp(1),pp(2),pp(3))
        write(lunit, '(A8, es12.4)') 'ro=', self%ro(pp(1),pp(2),pp(3))
        write(lunit, '(A8, es12.4)') 'dxi=', self%dxi(pp(1),pp(2),pp(3))

        do i=1, self%nel
            call self%elements(i)%print_info(lunit, pp)
        end do
        close(lunit)
    END SUBROUTINE print_ion_table

END MODULE composition_mod
