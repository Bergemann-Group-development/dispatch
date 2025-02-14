MODULE composition_mod
    USE trace_mod
    USE prec_mod
    USE op_math_mod
    USE handyfuncs_mod
    USE element_mod
    USE ion_mod, only: ion_type, debye_shield
    USE level_mod, only: level_type, wcalc
    USE bbtrans_mod, only:bbtrans_type
    USE bftrans_mod, only:bftrans_type, H_occupation_probability
    USE molecule_mod, only:molec_type
    USE hprofl_mod, only:prepare_hydrogen_damping, add_lte_h_opacities
    USE osmet_mod, only:osmet
    USE omp_lib
    USE os_mod
    USE io_mod
    use, intrinsic :: ieee_arithmetic
    implicit none
    private

    type, public :: composition_type
        logical :: initialised, ll_initialised=.false., loaded=.false.
        integer :: nel, nmolec, nline, verbose
        real(dp) :: xmy, suma, grph, tmolim, molh, FeH
        type(element_type), allocatable :: elements(:)
        type(molec_type),   allocatable :: molecules(:)
        type(bftrans_type), allocatable :: trans(:)
        type(level_type), dimension(:), allocatable :: input_H_levels
        real(dp), dimension(:,:,:), allocatable :: temp, ne, rho, zero
        real(dp), dimension(:,:,:), allocatable :: EK_over_temp, pe, dxi, toth, saha
        real(dp), dimension(:,:,:), allocatable :: teta, pg, E, p6, shield, fo, vbarh, t3nhe
        real(4), dimension(:,:,:), allocatable :: chi500
    contains
        procedure :: init => injon
        procedure :: init_molecules
        procedure :: load_state_variables
        procedure :: iterate_thermal_eq
        procedure :: thermal_equilibrium
        procedure :: compute_energy_and_gas_pressure
        procedure :: set_aux_data
        procedure :: molecular_equilibrium
        procedure :: init_H_levels
        procedure :: Hoccupation_partition
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
        procedure :: compute_grams_per_hydrogen
        procedure :: compute_p6
        procedure :: overwrite_abundance
        procedure :: init_linelist
        procedure :: mask_linelist
        procedure :: add_linelist_opacities
        procedure :: print_ion_table
        procedure :: read_molec_data
        procedure :: print_opacity_contributions
    end type

    character(128), public, save :: exclude_elements=''
    logical, public :: scattering_absorption=.true., mhd_eos=.false.
    logical, public :: irwin_dis=.true., helium_scattering=.true.
    real(dp), public :: debug_opac=0
    integer, public :: switch=1
    real(dp), public :: ref_lam=5000.d0
    character(len=1024), public :: abundances
    logical, public :: relative_abundances=.false.

    character(len=3), save :: & 
        element_ids(92) = [ &
        'H ','He','Li','Be','B ','C ','N ','O ','F ', &
        'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', &
        'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co', &
        'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
        'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh', &
        'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
        'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu', &
        'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
        'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl', &
        'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
        'Pa','U ']

contains

    !===============================================================================
    SUBROUTINE injon(self, absdat_file, abund_file, FeH, &
        isotope_file, molec_file, tmolim, external_abund, alpha_enhancement, &
        r_process_enhancement)
        class(composition_type), target :: self
        character(len=*), intent(in) :: absdat_file, abund_file
        character(len=*), intent(in), optional :: isotope_file, molec_file
        real(4), intent(in), optional :: tmolim
        real(8), intent(in) :: FeH
        real(4), intent(in), optional :: external_abund(:), alpha_enhancement, r_process_enhancement
        integer :: absdat, lunit, iostat
        integer :: i, j, k, l, m, nl, nk, niso
        integer :: lvl_count, skip_back
        real(dp) :: abund, xion, g2, dump, tecor, decor
        character(len=3) :: id
        character(len=6) :: key
        character(len=65) :: text
        character(len=10) :: s10
        integer :: g0, nqfix, nqtemp, ifish, ncor
        real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c, scale
        integer :: nkomp,nextl,nutzl,nextt,nutzt,null
        integer :: ilogl, kvadl, minex, maxex, nnu, iscal
        integer :: ilogt, kvadt, minet, maxet, ntt, iteta, istim

        class(element_type),    pointer :: elements(:)
        class(element_type),    pointer :: element, H
        class(molec_type),      pointer :: molecule
        class(ion_type),        pointer :: ion, Hm, HI
        class(isotope_type),    pointer :: iso
        class(level_type),      pointer :: lvl
        class(bftrans_type), pointer :: trans
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')

        106 FORMAT(4(F10.3,F10.4))
        call os%check_file(absdat_file)
        open(file=absdat_file, newunit=absdat)

        read(absdat, *) self%nel
        allocate(self%elements(self%nel))

        elements => self%elements

        read(absdat, *) elements%id
        read(absdat, *) elements%awgt
        read(absdat, *) elements%nion

        do i=1, self%nel
            element => elements(i)
            element%np = findloc(element_ids, element%id(:3), dim=1)

            if (i>1) then
                allocate(element%ions(element%nion))
            else
                allocate(element%ions(0:element%nion))
            endif
        end do

        call os%check_file(abund_file)
        open(file=abund_file, newunit=lunit, iostat=iostat)
        call nextentry(lunit)
        do while (iostat==0)
            read(lunit, *, iostat=iostat) id, abund
            call self%element_by_id(id, element, verbose=.true.)
            if(associated(element)) element%abund = abund
        end do
        close(lunit)

        if (elements(1)%abund == 1) then
            elements%abund = log10(elements%abund) + 12.0
        endif

        self%FeH = merge(FeH, 0.0, FeH > -99)
        do i=1, self%nel
            element => elements(i)
            select case (element%id)
            case ('H', 'He', 'Li')
                ! do not scale with FeH
            case default
                element%abund = element%abund + self%FeH
            end select
            if (present(alpha_enhancement)) then
                select case (element%id)
                case ('O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti')
                    element%abund = element%abund + alpha_enhancement
                end select
            endif
            if (present(r_process_enhancement)) then
                select case (element%id)
                case ('Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', &
                      'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  &
                      'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', &
                      'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', &
                      'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'T', 'U')
                    element%abund = element%abund + r_process_enhancement
                end select
            endif
        end do

        if (present(external_abund)) then
            ! do not overwrite external abundances with FeH
            do i=1, size(external_abund)
                call self%element_by_id(element_ids(i), element, verbose=.true.)
                if(associated(element)) element%abund = external_abund(i)
            end do
        endif

        if (len_trim(abundances) > 0) then
            do i=1, self%nel
                element => elements(i)
                j = index(abundances, trim(element%id)//'=')
                if (j > 0) then
                    l = len(abundances)
                    l = j + index(abundances(j:l), ' ') - 1
                    j = j + len_trim(element%id) + 1
                    read(abundances(j:l), *) abund
                    if (relative_abundances) then
                        element%abund = element%abund + abund
                    else
                        element%abund = abund
                    endif
                    print *, 'Changing abundance of ', trim(element%id), ' to ', element%abund
                endif
            end do
        endif

        elements%abund = 10**(elements%abund - elements(1)%abund)
        call self%compute_grams_per_hydrogen()

        do i=1, self%nel
            element => elements(i)

            do j=1, elements(i)%nion
                ion => element%ions(j)
                ion%j = j
                ion%dxi_fac = j

                !g0 = THE STATISTICAL WEIGHT OF THE GROUND LEVEL
                !nk = THE NUMBER OF ELECTRON CONFIGURATIONS CONSIDERED.

                read(absdat, *, iostat=iostat) g0
                backspace(absdat)

                if (g0 > 0) then
                    read(absdat, *, iostat=iostat) g0, nk
                    ion%nk = nk+1
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
                    ion%levels(1)%g = g0
                    ion%levels(1)%xion = 0.d0

                    lvl_count = 1
                    do k=1, nk
                        lvl_count = lvl_count + 1
                        lvl => ion%levels(lvl_count)

                        !xion = THE IONIZATION ENERGY IN ELECTRON VOLTS
                        !g2   = THE STATISTICAL WEIGHT (2L+1)*(2J+1)
                        !xl   = THE LOWEST QUANTUM NUMBER OF THE ASYMPTOTIC
                        !       (HYDROGENIC) PARTOF THE PARTITION FUNCTION,
                        !nl   = THE NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSION
                        !       FOR THE 'MIDDLE PART' OF THE PARTITION FUNCTION ('QPRIME').

                        read(absdat, *) lvl%xion, lvl%g, lvl%xl, nl

                        if (nl > 0) then
                            read(absdat, 106) (ion%levels(m)%xion, ion%levels(m)%g, m=lvl_count+1, lvl_count+nl)
                            lvl_count = lvl_count + nl
                        endif
                    end do
                    ion%xion = ion%levels(2)%xion

                else ! IRWIN FORMAT
                    allocate(ion%irwin(6))
                    read(absdat, *) g0, ion%xion, ion%irwin
                    ion%irwin(:) = ion%irwin(6:1:-1)
                endif
            end do
        end do

        do i=1, self%nel
            element => elements(i)
            read(absdat, *) element%active, (element%ions(j)%active, j=1, element%nion)
        end do

        read(absdat, *) nqfix ! number of ions with constant partition FUNCTION

        do m=1, nqfix
            read(absdat, *) i, j, elements(i)%ions(j)%parco

            !IQFIX(I,J) = 0 MEANS THAT THE PARTITION FUNCTION IS CONSIDERED TO BE CONSTANT.
            elements(i)%ions(j)%fixed_part = .true.
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

        if (present(tmolim)) then
            if (tmolim > -98) then
                self%tmolim = tmolim
            endif
        endif

        read(absdat, '(A9)') text
        if (text(1:9)=='* CORONAL') then
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
        else
            backspace absdat
        endif

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
        HI%nk = 15

        Hm%fixed_part = .true.
        Hm%parco = 1
        Hm%xion = 0.754
        Hm%dxi_fac = 2
        ! --- Different Hm ionisation in Oslo legacy code ---!
        !Hm%xion = 0.747
        ! --------------------------------------------------!
        RYD = HI%xion ! overwriting RYD constant
        ! --------------------------------------------------!

        call self%init_molecules(molec_file)

        !101 FORMAT(7X,I3,5(9X,I1))
        !999 FORMAT(4(9X,I1),8X,I2,9X,I1)

        !997 FORMAT(6(F10.0)) ! lam format
        !998 FORMAT(6(F10.5)) ! xkap format
        !102 FORMAT(6F10.4)
        !105 FORMAT(F6.3,F4.0,F5.1,I5)

        !read(absdat,101) nkomp,nextl,nutzl,nextt,nutzt,null

        call read_fixed(absdat, '(A10)', nkomp)
        call read_fixed(absdat, '(A10)', nextl)
        call read_fixed(absdat, '(A10)', nutzl)
        call read_fixed(absdat, '(A10)', nextt)
        call read_fixed(absdat, '(A10)', nutzt)
        call read_fixed(absdat, '(A10)', null)
        read(absdat, *) ! advance

        allocate(self%trans(nkomp))
        self%trans%verbose = self%verbose

        self%elements(1)%ions%in_opac = .true.
        self%elements(2)%ions%in_opac = .true.
        call self%molecule_by_id('H2', molecule)
        molecule%in_opac = .true.

        do i=1, nkomp
            trans => self%trans(i)

            read(absdat,'(A65)') text
            text = uppercase(text)
            trans%label = text(:20)
            trans%key = text(1:5)
            trans%peach = index(trans%label, 'PEACH') > 0
            trans%lam3 = trans%peach
            
            select case(trim(trans%key))
            case('H- BF')
                trans%ion => Hm
            case('H BF')
                trans%ion => HI
                ! hydrogen energy level
                trans%k = extract_int(text(1:10))
                trans%lam3 = .true.
            case('H FF')
                trans%ion => HI
                trans%lam3 = .true.
            case('H- FF', '(H+H)', 'H2+')
                trans%ion => HI
            case('H2-')
                call self%molecule_by_id('H2', molecule)
                trans%ion => molecule
            case('HE-')
                call self%element_by_id('HE', element)
                trans%ion => element%ions(1)
            case ('H2O-')
                call self%molecule_by_id('H2O', molecule)
                trans%ion => molecule
                trans%molecular = .true.
            case default
                if (text(1:1) == '/') then
                    trans%key = text(2:6)
                    id = text(2:3)
                    call self%element_by_id(id, element)
                    read(absdat, *) j, xion, trans%level%g
                    trans%level%xion = xion * hc_e
                    trans%ion => element%ions(j)
                    trans%regular = .true.
                else
                    id = text(1:2)
                    j = scan(id,'_')
                    if (j > 1) id = id(1:j-1)
                    trans%key = id // text(3:4)
                    j = scan(trans%key, 'IV', back=.true.)
                    if (j > 0) trans%key = trans%key(1:j)
                    j = scan(trans%key, '-')
                    trans%regular = j == 0
                    if (j > 0) then
                        key = trans%key(1:j-1)
                    else
                        key = trans%key
                    endif
                    call self%molecule_by_id(key, molecule)
                    trans%molecular = associated(molecule)
                    trans%ion => molecule
                    if (trans%regular) then
                        if (.not. trans%molecular) then
                            call self%ion_by_id(trans%key, ion)
                            if (associated(ion)) then
                                trans%ion => ion
                                if (trans%peach) trans%ion => ion%above
                            endif
                        endif
                        trans%level%xion = 0
                        trans%level%g = 1
                    elseif (.not. trans%molecular) then
                        call self%ion_by_id(key//' I', ion)
                        trans%ion => ion
                    endif
                endif
            end select

            if (.not. associated(trans%ion)) then
                trans%active = .false.
                if (trans%key /= 'ME') then
                    print *, trans%key,  ' deactivated!'
                endif
            else
                trans%ion%in_opac = .true.
            endif

            !read(absdat,999) ilogl,kvadl,minex,maxex,nnu
            call read_fixed(absdat, '(A10)', ilogl)
            call read_fixed(absdat, '(A10)', kvadl)
            call read_fixed(absdat, '(A10)', minex)
            call read_fixed(absdat, '(A10)', maxex)
            call read_fixed(absdat, '(A10)', nnu)
            read(absdat, *) ! advance

            trans%nnu = nnu
            trans%ilogl = ilogl > 0
            trans%kvadl = kvadl > 0
            trans%minex = minex > 0
            trans%maxex = maxex > 0

            if (trans%key == 'H BF') then
                trans%maxex  = H_occupation_probability
                trans%loglam = H_occupation_probability
            endif

            allocate(trans%lam(nnu))
            allocate(trans%ldisk(nnu))

            trans%ldisk = .false.

            read(absdat, *) trans%lam

            if (trans%loglam) then
                allocate(trans%lnlam, source=log(trans%lam))
            endif

            trans%minl = trans%lam(1)
            trans%maxl = trans%lam(nnu)


            do j=1, nnu-1
                if (trans%lam(j+1) - trans%lam(j) <= 10.) then
                    trans%ldisk(j) = .true.
                endif
            end do

            !if (count(trans%ldisk)>0) then
            !    write(text, '(A,I0,A)') '(A10, ":[",', count(trans%ldisk), '(F14.1,","),F14.1,"],")'
            !    print text, "'"//trim(trans%key)//"'", pack(trans%lam, trans%ldisk), trans%maxl
            !else
            !    write(text, '(A,I0,A)') '(A10, ":[",F14.1,"],")'
            !    print text, "'"//trim(trans%key)//"'", trans%maxl
            !endif

            !read(absdat,999) ilogt,kvadt,minet,maxet,ntt,iteta
            call read_fixed(absdat, '(A10)', ilogt)
            call read_fixed(absdat, '(A10)', kvadt)
            call read_fixed(absdat, '(A10)', minet)
            call read_fixed(absdat, '(A10)', maxet)
            call read_fixed(absdat, '(A10)', ntt)
            call read_fixed(absdat, '(A10)', iteta)
            call read_fixed(absdat, '(A10)', iscal)
            call read_fixed(absdat, '(A10)', istim)
            read(absdat, *) ! advance

            if (iscal>0) read(absdat, *) scale

            trans%ntt = ntt
            trans%ilogt = ilogt > 0
            trans%kvadt = kvadt > 0
            trans%minet = minet > 0
            trans%maxet = maxet > 0
            trans%iteta = iteta > 0
            trans%pedep = iscal > 0
            trans%istim = istim > 0

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

            if (iscal > 0) then
                trans%xkap = trans%xkap * scale
            endif
            
            if (trans%ilogl .or. trans%ilogt) then
                allocate(trans%log_xkap, mold=trans%xkap)
                trans%log_xkap = trans%xkap
                where (trans%log_xkap <= 0) trans%log_xkap = 1d-37
                trans%log_xkap = log(trans%log_xkap)
            endif
        end do

        close(absdat)

        if (present(isotope_file)) then
          if (isotope_file /= '') then
            call os%check_file(isotope_file)
            open(file=isotope_file, newunit=lunit)
            read(lunit, *, iostat=iostat) text
            do while (iostat==0)
                do while (text(1:1) == '!')
                    read(lunit, *, iostat=iostat) text
                end do
                backspace lunit
                read(lunit, *, iostat=iostat) id, niso
                call self%element_by_id(id, element)
                if (associated(element)) then
                    element%niso = niso
                    allocate(element%isotope(element%niso))
                    do i=1, element%niso
                        iso => element%isotope(i)
                        read(lunit, *) iso%iso, iso%awgt, iso%frac
                    end do
                    do j=1, element%nion
                        element%ions(j)%isotope => element%isotope
                        element%ions(j)%niso = element%niso
                    end do
                endif
                read(lunit, *, iostat=iostat) text
            end do
            close(lunit)
            do i=1, self%nmolec
              molecule => self%molecules(i)
              call molecule%generate_isotopic_combinations()
            end do
          endif
        endif
    
        if (allocated(HI%levels)) then
            allocate(self%input_H_levels, source=HI%levels)
        endif

        self%initialised = .true.
        call trace%end ()

    contains
        !===========================================
        subroutine read_fixed(unit, fmt, val)
            integer :: i, unit, val
            character(len=*) :: fmt
            character(len=32) :: text
            !---------------------------------------
            read(unit, fmt, advance='no') text
            i = scan(text, '0123456789')
            read(text(i:len(text)), *) val
        end subroutine read_fixed
    END SUBROUTINE injon

    !===============================================================================
    SUBROUTINE init_linelist (self, linelist_file)
        class(composition_type) :: self
        class(ion_type), pointer :: ion
        class(molec_type), pointer :: mol
        class(particle_type), pointer :: p
        class(element_type), pointer :: element
        class(bbtrans_type), dimension(:), allocatable :: lines
        class(bbtrans_type), pointer :: line
        character(len=*) :: linelist_file
        character(20) :: mass_id
        character(20) :: ion_id
        character(len=512) :: line_text
        real(dp) :: lam0, ent8, ent9, cnt_ev, lam_vac, vbar2e4
        real(dp) :: upper_ev, upper_g, lower_g, gij, f
        character :: lower_orb, upper_orb
        character(len=150) :: label
        character(len=32) :: gq
        character :: char
        integer :: np, iso
        integer :: lunit, iostat, iostat2, j, nline, i, prev, iso_idx, nkh
        integer, save :: ntotal, nactive
        logical :: isotope, ok
        !-----------------------------------------------------------------------------
        if (linelist_file == '' .and. .not. self%ll_initialised) return
        self%ll_initialised = .true.
        ntotal = 0
        nactive = 0
        nkh = 0

        call os%check_file(linelist_file)
        open(newunit=lunit, file=linelist_file, iostat=iostat)
        do while (iostat==0)
            read(lunit,*,iostat=iostat) mass_id, j, nline
            if (iostat < 0) exit ! End-of-File
            call print_failed_line(lunit, iostat)
            read(lunit,'(A)',iostat=iostat) ion_id
            call print_failed_line(lunit, iostat)
            call replace_char(mass_id, '.', ' ')
            read(mass_id, *) np, iso
            isotope = iso > 0
            ion_id = trim_between_words(ion_id)
            nullify(mol)
            nullify(p)
            call self%ion_by_id(ion_id(1:6), ion)
            ok = associated(ion)
            
            if (.not. ok) then
                call self%molecule_by_id(ion_id(1:6), mol)
                ok = associated(mol)
            endif

            if (associated(ion)) then
                p => ion
            elseif (associated(mol)) then
                p => mol
            endif
            
            if (associated(p) .and. isotope) then
                if (associated(p%isotope)) then
                    iso_idx = findloc(p%isotope%iso, iso, 1)
                    ok = iso_idx > 0
                else
                    ok = iso == nint(p%awgt)
                    isotope = .false.
                endif
            endif

            if (.not. ok) then
                if (associated(p) .and. isotope) then
                    print *, p%id, ' isotope ', iso, ' not found!'
                endif
                do i=1, nline
                    read(lunit, *, iostat=iostat) ! skipping this ion
                    call print_failed_line(lunit, iostat)
                end do
                cycle
            endif
            
            p%in_linelist=.true.

            if (p%nline == 0) then
                allocate(p%line(nline))
                p%nline = nline
                prev = 0
            else
                allocate(lines(nline))
                lines = p%line
                deallocate(p%line)
                prev = p%nline
                p%nline = p%nline + nline
                allocate(p%line(p%nline))
                p%line(1:prev) = lines
                deallocate(lines)
            endif

            do i=1, nline
                line => p%line(prev+i)
                read(lunit,'(A)', iostat=iostat) line_text

                if (p%id == 'H I') then
                    read(line_text, *, iostat=iostat2) lam0, line%i, line%j, line%lower_ev, upper_ev, line%loggf, label, line%ga

                    ! Stop if nothing has worked
                    call print_failed_line(lunit, iostat2)

                    nkh = max(line%j, nkh)
                    lower_g = 2 * line%i**2
                    line%upper_g = 2 * line%j**2
                    lam_vac = air2vac(lam0) * 1e-8
                    line%nu0 = CC / lam_vac
                    gij = lower_g / line%upper_g
                    f = 10**line%loggf / lower_g
                    call line%init_einstein_coefficients(f, gij)
                    call prepare_hydrogen_damping(line) ! ga is overwritten
                else
                    ! It will try to read as many of these parameters as possible,
                    ! but only lam0 lower_ev and loggf are required.
                    ! Broadening parameters are zero by default if not overwritten here.
                    read(line_text, *, iostat=iostat2) lam0, line%lower_ev, line%loggf, line%gv, upper_g, line%ga, gq
                    call print_failed_line(lunit, iostat2)
                    ! read gq if it is there
                    read(gq, *, iostat=iostat2) line%gq
    
                    lam_vac = air2vac(lam0) * 1e-8
                    line%nu0 = CC / lam_vac

                    f = 10**line%loggf
                    vbar2e4 = p%compute_vbar2(10000.)
                    if (isotope) then
                        line%iso_idx = iso_idx
                        vbar2e4 = p%isotope(iso_idx)%compute_vbar2(10000.)
                        f = f * p%isotope(iso_idx)%frac
                    endif

                    call line%init_einstein_coefficients(f, 1.d0)
                    call line%prepare_Stark(line%gq)
                    call line%prepare_VanderWaals(line%gv, j, p%xion, vbar2e4)
                endif
                allocate(line%used, source=.false.)
                allocate(line%strength, source=0.d0)
                if (line%active) nactive = nactive + 1
                ntotal = ntotal + 1
            end do
        end do
        close(lunit)
        print '(I12, a, I12, a, a)', nactive, ' of ', ntotal, ' lines read from linelist ', trim(linelist_file)
        self%nline = nactive
        self%elements(1)%ions(1)%nk = max(self%elements(1)%ions(1)%nk, nkh)
        call self%mask_linelist()
    END SUBROUTINE init_linelist

    !===============================================================================
    SUBROUTINE mask_linelist (self)
      class(composition_type), target :: self
      class(element_type), pointer :: element
      class(molec_type), pointer :: mol
      integer :: n_excluded, i, ii, p
      character(len=3) :: id
      logical :: new_word
      ! ----------------------------------------------------------------------------
      call count_words(exclude_elements, n_excluded)
      p = 1
      do i=1, n_excluded
        id = ''
        do ii=1, 3
          new_word = .false.
          do while (exclude_elements(p:p) == ' ' .and. p < len(exclude_elements))
            new_word = .true.
            p = p + 1
          end do
          if (new_word) exit
          id(ii:ii) = exclude_elements(p:p)
          p = p + 1
        end do
        call self%element_by_id(id, element)
        if (associated(element)) then
          element%ions%in_linelist = .false.
          print *, 'removing ', trim(element%id), ' from linelist'
        else
          call self%molecule_by_id(id, mol)
          if (associated(mol)) then
            mol%in_linelist = .false.
            print *, 'removing ', trim(mol%id), ' from linelist'
          endif
        endif
      end do
    END SUBROUTINE mask_linelist

    !===============================================================================
    SUBROUTINE add_linelist_opacities(self, nu, all_chi, dlam, vproj, vmic, all_eta)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        class(ion_type), pointer :: ion, HeI
        class(molec_type), pointer :: mol
        class(isotope_type), pointer :: iso
        class(bbtrans_type), pointer :: line
        real(8), intent(in), dimension(:) :: nu
        real(4), intent(inout), dimension(:,:,:,:) :: all_chi
        real(4), intent(inout), dimension(:,:,:,:), optional :: all_eta
        real(4), intent(in), dimension(:,:,:), contiguous, target, optional :: vproj, vmic
        real(8), intent(in) :: dlam
        real(4), dimension(:,:,:), allocatable :: tothi, adamp, stark_pe, ne, temp, EK_over_temp, dxij, chi_fac, unsold
        real(4), dimension(:,:,:), allocatable, target :: inv_part, vbar2, inv_df, iso_vbar2, iso_inv_df, iso_inv_part
        real(4), dimension(:,:,:), pointer, contiguous :: vb2, idopf, ipart, vpr, fo, t3nhe, ni, nj
        real(4), dimension(:,:,:,:), pointer, contiguous :: eta, chi
        real(4), dimension(:,:,:,:), allocatable, target :: chi_l
        real(8), dimension(:), allocatable :: qq
        integer :: i, j, k, l, f, ired, iblue, nactive, prev_iso, m(3)
        real(8) :: cnt_ev, gf, awgt, iso_awgt, dl
        real(4) :: scalar, bnu, h_kT, hv_kT
        logical, dimension(:), allocatable :: fmask
        integer, save:: itimer=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%add_linelist_opacities', itimer=itimer)
        m = shape(self%ne)
        allocate(tothi(m(1),m(2),m(3)))
        allocate(adamp,        mold=tothi)
        allocate(stark_pe,     mold=tothi)
        allocate(ne,           mold=tothi)
        allocate(temp,         mold=tothi)
        allocate(EK_over_temp, mold=tothi)
        allocate(inv_df,       mold=tothi)
        allocate(iso_inv_df,   mold=tothi)
        allocate(inv_part,     mold=tothi)
        allocate(iso_inv_part, mold=tothi)
        allocate(vbar2,        mold=tothi)
        allocate(iso_vbar2,    mold=tothi)
        allocate(dxij,         mold=tothi)
        allocate(chi_fac,      mold=tothi)
        allocate(unsold,       mold=tothi)
        allocate(fo,           mold=tothi)
        allocate(t3nhe,        mold=tothi)
        allocate(ni,           mold=tothi)
        allocate(nj,           mold=tothi)
        allocate(fmask(size(nu)), source=.false.)
        allocate(qq(size(nu)))

        allocate(chi_l, mold=all_chi)
        chi_l = 0.

        if (present(vproj)) then
            vpr => vproj
        else
            allocate(vpr, mold=tothi)
            vpr = 0
        endif

        ne   = real(self%ne)
        temp = real(self%temp)
        EK_over_temp = real(self%EK_over_temp)
        tothi = real(self%elements(1)%ions(1)%n) ! total neutral hydrogen density
        stark_pe = ne * temp**one_sixth

        element => self%elements(1)
        ion => element%ions(1)
        if (ion%nline > 0) then
            fo = self%fo
            vbar2 = self%vbarh
            t3nhe = self%t3nhe
            inv_df = element%inverse_doppler_factor(temp, vmic)

            do l=1, ion%nline
                line => ion%line(l)
                if (.not. line%active) cycle
                call line%get_range(nu, dlam, ired, iblue)
                if (iblue == 0) cycle

                ni = ion%levels(line%i)%n
                nj = ion%levels(line%j)%n

                fmask(ired:iblue) = .true.
                call add_lte_h_opacities(line, ired, iblue, nu, vpr, inv_df, temp, &
                ni, nj, ne, tothi, vbar2, fo, t3nhe, chi_l, all_chi)
            end do
        endif

        do i=2, self%nel
            element => self%elements(i)
            if (.not. any(element%ions%in_linelist)) cycle
            do j=lbound(element%ions,1), ubound(element%ions,1)
                ion => element%ions(j)
                call particle_add_lte_opacities(ion)
            end do
        end do

        do i=1, self%nmolec
            mol => self%molecules(i)
            call particle_add_lte_opacities(mol)
        end do

        do k=1, m(3)
        do j=1, m(2)
        do i=1, m(1)
            do f=1, size(nu)
                if (fmask(f)) then
                    all_chi(f,i,j,k) = all_chi(f,i,j,k) + chi_l(f,i,j,k)
                endif
            end do
        end do
        end do
        end do

        if (present(all_eta)) then
            do k=1, m(3)
            do j=1, m(2)
            do i=1, m(1)
                h_kT = h_k / temp(i,j,k)
                do f=1, size(nu)
                    if (fmask(f)) then
                        ! computing boltzmann factor
                        hv_kT = nu(f) * h_kT
                        bnu = exp(-hv_kT)
                        if (bnu > 1e-10) then
                            bnu = bnu / (1.0 - bnu)
                        endif
                        all_eta(f,i,j,k) = all_eta(f,i,j,k) + chi_l(f,i,j,k) * bnu
                    endif
                end do
            end do
            end do
            end do
        endif
        if (.not. present(vproj)) deallocate(vpr)
        call trace%end(itimer)

    contains
        ! --------------------------------------------------------------------------
        subroutine particle_add_lte_opacities(ptc)
            class(particle_type), target :: ptc
            !---------------------------------------------------------------------------
            if (.not. ptc%in_linelist) return
            if (.not. allocated(ptc%n)) stop 'particle density not allocated'

            inv_df = ptc%inverse_doppler_factor(temp, vmic)
            vbar2 = ptc%compute_vbar2(temp)
            unsold = vbar2**0.3 * tothi
            prev_iso = 0

            ! ptc%part is actually the inverse partition
            inv_part = ptc%n * ptc%part
            dxij = real(self%dxi) * (ptc%j-1)

            do l=1, ptc%nline
                line => ptc%line(l)
                if (.not. line%active) cycle
                call line%get_range(nu, dlam, ired, iblue)
                if (iblue == 0) cycle

                fmask(ired:iblue) = .true.

                if (line%iso_idx == 0) then
                    vb2 => vbar2
                    idopf => inv_df
                elseif (line%iso_idx /= prev_iso) then
                    iso => ptc%isotope(line%iso_idx)
                    iso_inv_df = iso%inverse_doppler_factor(temp, vmic)
                    iso_vbar2 = iso%compute_vbar2(temp)
                    vb2 => iso_vbar2
                    idopf => iso_inv_df
                endif
                prev_iso = line%iso_idx

                call line%add_lte_opacities(ired, iblue, nu, qq, &
                    vpr, idopf, EK_over_temp, dxij, inv_part, ne, tothi, vb2, stark_pe, unsold, chi_l, all_chi)
            end do
        end subroutine particle_add_lte_opacities
    END SUBROUTINE add_linelist_opacities

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
                endif
                if(j < ubound(element%ions,1)) then
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
            if (associated(trans%ion)) then
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
            endif
        end do
        call trace%end ()
    END SUBROUTINE reassign_pointer

    !===============================================================================
    SUBROUTINE init_molecules(self, molec_file)
        class(composition_type), target :: self
        class(molec_type),   pointer :: molecule
        class(element_type), pointer :: element
        class(ion_type),     pointer :: ion
        character(len=*) :: molec_file
        character(len=3) :: mol_ids(12)
        character(len=6) :: ion_id
        integer :: i, j
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init_molecules')

        mol_ids = ['H2 ', 'H2+', 'H2O', 'OH ', 'CH ', 'CO ', 'CN ', 'C2 ', 'N2 ', 'O2 ', 'NO ', 'NH ']

        if (self%tmolim > 0) then
            self%nmolec = size(mol_ids)
        else
            self%nmolec = 0
        endif

        allocate(self%molecules(self%nmolec))

        do i=1, self%nmolec
            molecule => self%molecules(i)
            molecule%nions = len(trim(mol_ids(i)))
            molecule%id = trim(mol_ids(i))
            molecule%awgt = 0
            !call molecule%init_mol_dissociation()

            allocate(molecule%ions(molecule%nions))

            do j=1, molecule%nions
                ion_id = molecule%id(j:j)

                select case (ion_id(:3))
                case ('H', 'C', 'N', 'O')
                    ion_id(3:3) = 'I'
                    call self%ion_by_id(ion_id, ion)
                case ('2')
                    ion => ion ! use previous ion a second time
                case ('+')
                    molecule%nions = molecule%nions - 1
                    molecule%ions(j-1)%p => ion%above
                    cycle
                case default 
                    cycle
                end select

                molecule%ions(j)%p => ion
                ion%nmolec = ion%nmolec + 1 
                molecule%awgt = molecule%awgt + ion%awgt
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

        do i=1, self%nmolec
            molecule => self%molecules(i)

            do j=1, molecule%nions
                ion_id = molecule%id(j:j)

                select case (ion_id)
                case ('H', 'C', 'N', 'O')
                    ion_id(3:3) = 'I'
                    call self%ion_by_id(ion_id, ion)
                case ('2')
                    ion => ion ! use previous ion a second time
                case ('+')
                    ion => ion%above
                case default 
                    cycle
                end select
                
                ion%molecules(ion%mcount)%p => molecule
                ion%mcount = ion%mcount + 1 
            end do
        end do
        call self%read_molec_data(molec_file)
        call trace%end ()
    END SUBROUTINE init_molecules

    !===============================================================================
    SUBROUTINE load_state_variables(self, temp, ne, rho, use_rho, use_ne)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        class(ion_type), pointer :: ion
        class(molec_type), pointer :: mol
        real(dp), dimension(:,:,:), intent(in) :: temp, ne, rho
        logical, intent(in), optional :: use_rho, use_ne
        real(dp), dimension(:,:,:), allocatable :: ro, eh, charge
        real(4), dimension(:,:,:,:), pointer :: chi500
        real(dp) :: nu500(1), nudebug(1)
        integer :: i,j, nx, ny, nz, lunit, m(3)
        integer :: itimer=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%load_state_variables', itimer=itimer)

        if (.not. self%initialised) then
            print *, 'composition has not been initialised yet'
            stop
        endif

        if (self%loaded) then
            deallocate(self%temp)
            deallocate(self%ne, self%rho, self%pe, self%pg, self%E)
            deallocate(self%toth, self%saha, self%teta, self%dxi)
            deallocate(self%EK_over_temp, self%p6, self%zero)
        endif

        allocate(self%temp, self%ne, self%rho, self%pe, self%pg, self%E, mold=temp)
        allocate(self%toth, self%saha, self%teta, self%dxi, eh, mold=temp)
        allocate(self%EK_over_temp, self%p6, ro, self%zero, mold=temp)

        self%zero = 0.d0
        self%E = 0
        
        self%ne = ne
        if (option(use_ne, .true.)) then
            call self%thermal_equilibrium(temp, ne, self%rho)
        else
            call self%iterate_thermal_eq(temp, rho, self%ne, tol=1d-6)
            call self%thermal_equilibrium(temp, self%ne, self%rho)
        endif

        if (option(use_rho, .true.)) then
            self%rho = rho
        endif

        self%toth = self%rho / self%grph
        call self%set_aux_data(temp, self%ne)

        do i=1, self%nel
            element => self%elements(i)
            do j=lbound(element%ions,1), ubound(element%ions,1)
                ion => element%ions(j)
                ion%n = ion%f * self%toth
            end do
        end do

        do i=1, self%nmolec
            mol => self%molecules(i)
            mol%n = mol%f * self%toth
        end do

        call self%compute_energy_and_gas_pressure(eh, self%E, self%pg)
        call self%compute_p6(self%p6)
        call self%init_H_levels()

        do i=1, self%nel
            element => self%elements(i)
            do j=lbound(element%ions,1), ubound(element%ions,1)
                ion => element%ions(j)
                call ion%deallocate()
            end do
        end do

        do i=1, self%nmolec
            mol => self%molecules(i)
            call mol%deallocate()
        end do

        m = shape(self%temp)
        allocate(self%chi500(m(1),m(2),m(3)))
        chi500(1:1,1:m(1),1:m(2),1:m(3)) => self%chi500
        nu500 = 1d8 * CC / air2vac(ref_lam)
        call self%absorption(nu500, chi500, sorted=.true.)

        if (debug_opac > 0) then
            nudebug = 1d8 * CC / air2vac(debug_opac)
            call self%print_opacity_contributions(nudebug)
        endif

        self%loaded = .true.
        call trace%end(itimer)
    END SUBROUTINE load_state_variables

    !=========================================================================
    SUBROUTINE thermal_equilibrium(self, temp, ne, rho)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        class(ion_type), pointer :: ion, HI, HeI
        real(dp), intent(in),  dimension(:,:,:) :: temp, ne
        real(dp), intent(out), dimension(:,:,:) :: rho
        real(dp), dimension(:,:,:), allocatable :: fe, xnen, xnecno, toth, dxi, hpart
        integer :: i, j, itimer=0, itimer2=0, itimer3=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%thermal_equilibrium', itimer=itimer)
        allocate(fe, xnen, xnecno, toth, mold=temp)
        call self%set_aux_data(temp, ne)
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('He I', HeI)

        fe = 0
        xnecno = 0

        do i=1, self%nel
            element => self%elements(i)
            call trace%begin ('composition_type%partition_functions', itimer=itimer2)
            call element%compute_partition_functions(self%temp, self%EK_over_temp, self%ne)
            call trace%end(itimer2)

            call trace%begin ('composition_type%ion_equilibrium', itimer=itimer3)
            call element%compute_ion_equilibrium(self%EK_over_temp, self%saha, self%dxi)
            call trace%end(itimer3)
            
            do j=lbound(element%ions, 1), ubound(element%ions, 1)
                ion => element%ions(j)

                ! number of electrons from ion
                xnen = ion%f * (j-1)
                !if (element%id /= 'H') fe = fe + xnen

                select case (element%id)
                case ('H','He')
                  ! do not add
                case default
                  fe = fe + xnen
                end select

                select case (element%id)
                case ('H','C','N','O')
                    ! do not add
                case default
                    ! number of electrons from elements other than H, C, N, O
                    xnecno = xnecno + xnen
                end select
            end do
        end do

        if (mhd_eos) then
          allocate(dxi, hpart, mold=self%temp)
          do i=1, 5
            toth = ne / (fe + HI%above%f + HeI%above%f + 2*HeI%above%above%f)
            HI%n = HI%f * toth
            HeI%n = HeI%f * toth
            call self%Hoccupation_partition(HI%n, HeI%n, hpart)
            hpart = max(1e-40, hpart)
            dxi = - log(hpart*HI%part) / self%EK_over_temp
            call self%elements(1)%compute_ion_equilibrium(self%EK_over_temp, self%saha, dxi=dxi)
            call self%elements(2)%compute_ion_equilibrium(self%EK_over_temp, self%saha, dxi=dxi)
          end do
        endif

        fe = fe + HI%above%f + HeI%above%f + 2*HeI%above%above%f

        if (self%nmolec > 0) then
            call self%molecular_equilibrium(self%temp, self%teta, self%pe, xnecno, fe)
            where (fe < 0)
                fe = ieee_value(fe, ieee_quiet_nan)
            endwhere
        endif

        ! total hydrogen number density, computed from scratch
        toth = ne / fe
        ! mass density computed from scratch
        rho = toth * self%grph
        call trace%end(itimer)
    END SUBROUTINE thermal_equilibrium

    !=========================================================================
    SUBROUTINE set_aux_data(self, temp, ne)
        class(composition_type), target :: self
        real(dp), intent(in),  dimension(:,:,:) :: temp, ne
        !---------------------------------------------------------------------------
        self%ne = ne
        self%temp = temp
        self%teta = temp_to_theta / temp
        self%EK_over_temp = EK / temp
        self%pe = temp * BK * ne
        self%saha = 2 * op_saha_const * temp*sqrt(temp) / ne
        self%fo = 2.6031*EC * ne**two_third ! Holtsmark normal field strength
        
        if (debye_shield) then
          self%dxi = 4.98e-4*self%teta*sqrt(self%pe)
          ! --- this only shows up in the Oslo legacy code ---!
          !self%dxi=min(self%dxi,0.747*0.5)
          ! --------------------------------------------------!
        else
          self%dxi = 0
        endif
    END SUBROUTINE set_aux_data

    !=========================================================================
    ! find a state variables consistent with given temp & rho    
    !=========================================================================
    SUBROUTINE iterate_thermal_eq(self, temp, rho, ne, tol)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        real(dp), dimension(:,:,:), intent(in) :: temp, rho
        real(dp), dimension(:,:,:), intent(out) :: ne
        real(dp), intent(in) :: tol ! stop criterion
        logical, dimension(:,:,:), allocatable :: ne_converged
        real(dp), dimension(:,:,:), allocatable :: a, b, c, x, xold, fx, fa, fb, ro, pg
        real(dp), dimension(:,:,:), allocatable :: lrho
        integer :: i, nconverged, arr_size, m(3)
        integer, parameter :: nmax = 20 ! max number of iteration
        !---------------------------------------------------------------------------
        m = shape(temp)
        arr_size = size(temp)

        allocate(ne_converged(m(1),m(2),m(3)), source=.false.) ! cant use mold= on different kinds
        allocate(x, xold, a, b, c, fx, fa, fb, ro, pg, lrho, mold=temp)

        if (ne(1,1,1) == 0) then
            ! set initial value for electron density
            where (temp > 2.0d4) ! fully ionized gas
                ! threshold used to be 1d4, but led to convergence issues in some cases
                ne = rho / self%grph * 1e-1
            elsewhere ! neutrals only
                ne = rho / self%grph * 1e-4
            end where
        endif

        lrho = log(rho)

        call self%thermal_equilibrium(temp, ne, ro)
        fa = log(ro) - lrho
        ne_converged = abs(fa) < tol

        a = log(ne)
        b = a - 0.69*fa
        ne = exp(b)

        call self%thermal_equilibrium(temp, ne, ro)
        fb = log(ro) - lrho
        ne_converged = abs(fb) < tol

        nconverged = count(ne_converged)

        i = 0
        x = b
        !write(*,fmt="(5x,'it [',I2.2,']',1x,'n-converged:',1x,I8,1X,'/ total size: ',I8)") i, nconverged, arr_size

        do while (nconverged < arr_size .and. i <= nmax)
            xold = x
            where (.not. ne_converged)
                c = fb-fa ! log(ro2) - log(ro1)
                c = (b-a)/c
                x = a - c*fa
                x = max(x, min(a,b)-2.d0)
                x = min(x, max(a,b)+2.d0)
            end where
    
            ne_converged = abs(x-xold) < tol
            nconverged = count(ne_converged)
            if (arr_size == nconverged) exit ! all good !

            where (.not. ne_converged)
                ne = exp(x)
            endwhere

            call self%thermal_equilibrium(temp, ne, ro)
            fx = log(ro) - lrho
            where (abs(a-x) < abs(b-x)) 
                b = x 
                fb = fx
            else where 
                a = x 
                fa = fx 
            end where
            i = i+1
        end do

        if (nconverged /= arr_size) then
            write(*,fmt="(5x,'it [',I3.3,']',1x,'n-converged:',1x,I8,1X,'/ total size: ',I8)") i, nconverged, arr_size
        endif

        where (.not. ne_converged)
            self%ne = ieee_value(self%ne, ieee_quiet_nan)
            self%pg = ieee_value(self%pg, ieee_quiet_nan)
            self%E  = ieee_value(self%E,  ieee_quiet_nan)
        end where

        deallocate(ne_converged)
        deallocate(x,xold,a,b,fx,fa,fb)
    END SUBROUTINE iterate_thermal_eq

    !===============================================================================
    SUBROUTINE Hoccupation_partition(self, nh, nhe, part)
        class(composition_type), target :: self
        real(dp), intent(in), dimension(:,:,:) :: nh, nhe
        real(dp), intent(out), dimension(:,:,:) :: part
        class(ion_type), pointer :: HI
        class(level_type), pointer :: lvl
        integer :: i, itimer=0, itimer2=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%compute_occupation_weights', itimer=itimer)
        call self%ion_by_id('H I', HI)
        if (.not. allocated(HI%levels)) allocate(HI%levels(HI%nk))
        if (.not. allocated(self%shield)) then
            allocate(self%shield, mold=self%temp)
            self%shield = 0.08989 * self%ne**0.16667d0 / sqrt(self%temp) ! the shielding parameter
            self%shield = (1.d0 + self%shield)**3.15d0 / self%ne
        endif

        part = 0.
        do i=1, HI%nk
          lvl => HI%levels(i)
          if (.not. allocated(lvl%f)) then
            allocate(lvl%f, mold=self%temp)
            call os%mem_register(storage_size(lvl%f), product(HI%m), 'lvl%f')
            lvl%i = i
            lvl%xion = RYD * (1. - 1./i**2)
            lvl%g = 2 * i**2
            lvl%f = lvl%g * lvl%boltzmann(self%EK_over_temp, natural_log=.true.)
          endif
          if (.not. allocated(lvl%w)) then
            allocate(lvl%w, mold=self%temp)
            call os%mem_register(storage_size(lvl%w), product(HI%m), 'lvl%w')
          endif
          call trace%begin ('composition_type%compute_occupation_weights%wcalc', itimer=itimer2)
          call wcalc(nh, nhe, self%ne, self%shield, dble(i), lvl%w)
          call trace%end (itimer2)
          part = part + lvl%w * lvl%f
        end do
        call trace%end (itimer)
    END SUBROUTINE Hoccupation_partition

    !===============================================================================
    SUBROUTINE init_H_levels(self)
        class(composition_type), target :: self
        class(ion_type), pointer :: HI, HeI
        class(level_type), pointer :: lvl
        real(dp), dimension(:,:,:), allocatable :: h1, part
        integer, dimension(:,:,:), allocatable :: nniv
        integer :: i, itimer=0, itimer2=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init_H_levels', itimer=itimer)
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('HE I', HeI)

        if (allocated(self%shield)) deallocate(self%shield)
        allocate(self%shield, mold=self%temp)
        self%shield = 0.08989 * self%ne**0.16667d0 / sqrt(self%temp) ! the shielding parameter
        self%shield = (1.d0 + self%shield)**3.15d0 / self%ne

        !HI%nk = 15
        if (allocated(HI%levels)) deallocate(HI%levels)
        allocate(HI%levels(HI%nk))
        allocate(part, mold=self%temp)

        do i=1, HI%nk
          lvl => HI%levels(i)
          allocate(lvl%n, mold=self%temp)
          call os%mem_register(storage_size(lvl%n), product(HI%m), 'lvl%n')
          lvl%i = i
          lvl%xion = RYD * (1. - 1./i**2)
          lvl%g = 2 * i**2
          lvl%n = lvl%g * lvl%boltzmann(self%EK_over_temp, natural_log=.true.)
        end do

        part = HI%part
        if (H_occupation_probability) then
          part = 0.
          do i=1, HI%nk
            lvl => HI%levels(i)
            allocate(lvl%w, mold=self%temp)
            call os%mem_register(storage_size(lvl%w), product(HI%m), 'lvl%w')
            call trace%begin ('composition_type%init_H_levels%wcalc', itimer=itimer2)
            call wcalc(HI%n, HeI%n, self%ne, self%shield, dble(i), lvl%w)
            call trace%end (itimer2)
            part = part + lvl%w * lvl%n
          end do
          part = 1.d0 / part
        elseif (self%dxi(1,1,1) > 0) then
          allocate(h1,   source = sqrt(RYD/self%dxi) )
          allocate(nniv, source = int(h1)            )
          do i=1, HI%nk
              lvl => HI%levels(i)
              where (nniv < i)
                  lvl%n = 0
              elsewhere (nniv == i)
                  lvl%n = lvl%n * (h1 - nniv)
              end where
          end do
        endif
        do i=1, HI%nk
          lvl => HI%levels(i)
          lvl%n = HI%n * lvl%n * part
        end do
        call trace%end (itimer)
    END SUBROUTINE init_H_levels

    !===============================================================================
    SUBROUTINE compute_p6(self, p6)
        class(composition_type), target :: self
        class(ion_type), pointer :: HeI, HI
        class(molec_type), pointer :: H2, H2p
        real(dp), dimension(:,:,:), intent(out) :: p6
        integer :: x, y, z
        !---------------------------------------------------------------------------
        call self%molecule_by_id('H2',  H2)
        call self%molecule_by_id('H2+', H2p)
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('HE I', HeI)

        self%vbarh = HI%compute_vbar2(real(self%temp))
        self%t3nhe = (self%temp/10000.)**0.3 * HeI%n
        
        if (associated(H2) .and. associated(H2p)) then
            p6 = LOG10((HI%n + 0.41336*HeI%n + 0.85*(H2%n + H2p%n)) * BK * self%temp)
        else
            p6 = LOG10((HI%n + 0.41336*HeI%n) * BK * self%temp)
        endif
    END SUBROUTINE compute_p6

    !===============================================================================
    SUBROUTINE absorption(self, nu, chi, blz, include_absmet, sorted)
        !use phys_mod, only:planck
        class(composition_type), target :: self
        class(bftrans_type), pointer :: trans
        class(ion_type), pointer :: HeI
        real(dp), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
        real(4), dimension(:,:,:,:), contiguous, intent(in), optional :: blz
        real(dp), dimension(:), allocatable :: lam, hnu
        logical, optional :: include_absmet, sorted
        real(dp) :: xray2, epsln, stim
        integer :: i, x, y, z, f, nf, nx, ny, nz
        integer :: p(4)=pp, itimer=0, itimer2=0
        logical :: ok
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%absorption', itimer=itimer)
        allocate(lam, hnu, mold=nu)
        call self%ion_by_id('HE I', HeI)

        nx = size(self%temp, 1)
        ny = size(self%temp, 2)
        nz = size(self%temp, 3)
        nf = size(nu, 1)
        lam = cc / nu * 1e8
        hnu = h_e * nu
        chi = 0
        call self%add_hrest(chi, lam)

        ! add absorption from all tables, which do not include stimulated emission
        do i=1, size(self%trans)
            trans => self%trans(i)
            if (.not. trans%istim) then
                call trans%add_absorption(chi, lam, self%temp, self%EK_over_temp, self%teta, &
                                          self%pe, self%dxi, ne=self%ne, nhe=HeI%n, shield=self%shield, sorted=sorted)
            endif
        end do

        if (present(blz)) then
            do z=1, nz
            do y=1, ny
            do x=1, nx
            do f=1, nf
                ! all continua above are multiplied by stim
                !--------------------------------------------------------------------------------
                !stim = 1. - exp(-28556.*self%teta/lam(f))
                stim = 1. - blz(f,x,y,z)
                chi(f,x,y,z) = chi(f,x,y,z) * stim
            end do
            end do
            end do
            end do
        else
            do z=1, nz
            do y=1, ny
            do x=1, nx
            do f=1, nf
                stim = 1. - exp(-hnu(f) * self%EK_over_temp(x,y,z))
                chi(f,x,y,z) = chi(f,x,y,z) * stim
            end do
            end do
            end do
            end do
        endif

        do i=1, size(self%trans)
            trans => self%trans(i)
            if (trans%istim) then
                call trans%add_absorption(chi, lam, self%temp, self%EK_over_temp, self%teta, &
                self%pe, self%dxi, sorted=sorted)
            endif
        end do

        if (scattering_absorption) then
            call self%electron_scattering(chi)
            call self%rayleigh_scattering(lam, chi)
        endif

        if (present(include_absmet)) then
            if (osmet%initialised .and. include_absmet) then
                call trace%begin('osmet%add_line_opacities', itimer=itimer2)
                call osmet%add_line_opacities(nu, self%p6, self%temp, self%rho, chi)
                call trace%end(itimer2)
            endif
        endif

        call trace%end(itimer)
    END SUBROUTINE absorption

    !===============================================================================
    SUBROUTINE scattering (self, nu, eta, Jnu)
        class(composition_type), target :: self
        real(dp), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:,:), intent(inout) :: eta
        real(4), dimension(:,:,:,:), intent(in) :: Jnu
        real(4), dimension(:,:,:,:), allocatable :: sumsca
        real(4) :: bnu
        integer :: f, i, j, k, n(4)
        !---------------------------------------------------------------------------
        n = shape(Jnu)
        allocate(sumsca, mold=Jnu)
        sumsca = 0
        call self%electron_scattering(sumsca)
        call self%rayleigh_scattering(cc / nu * 1e8, sumsca)
        do k=1, n(4)
        do j=1, n(3)
        do i=1, n(2)
        do f=1, n(1)
            eta(f,i,j,k) = eta(f,i,j,k) + sumsca(f,i,j,k) * Jnu(f,i,j,k)
        end do
        end do
        end do
        end do
    END SUBROUTINE scattering

    !===============================================================================
    SUBROUTINE get_totn(self, id, abund, totn, el, mol)
        class(composition_type), target :: self
        class(element_type), intent(out), optional, pointer :: el
        class(molec_type), intent(out), optional, pointer :: mol
        class(element_type), pointer :: element
        class(molec_type), pointer :: molec
        character(len=*), intent(in) :: id
        real(dp), intent(in) :: abund
        real(dp), dimension(:,:,:) :: totn
        integer :: j
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%get_totn')
        element => null()
        molec   => null()

        call self%element_by_id(id, element)
        if (associated(element)) then
            totn = 0
            do j=1, element%nion
                totn = totn + element%ions(j)%n
            end do
        else
            call self%molecule_by_id(id, molec)
            if (associated(molec)) then
                totn = molec%n
            else
                totn = 10.0**(abund-12.0) * self%toth
            endif
        endif

        if (present(el) .and. present(mol)) then
            el => element
            mol => molec
        endif
        call trace%end ()
    END SUBROUTINE get_totn

    !===============================================================================
    SUBROUTINE molecular_equilibrium(self, temp, teta, pe, xnecno, fe)
        class(composition_type), target :: self
        real(dp), dimension(:,:,:), intent(in) :: temp, pe, teta, xnecno
        real(dp), dimension(:,:,:), intent(inout) :: fe
        real(dp), dimension(:,:,:), allocatable :: r, s, b, pefe, p, q
        class(molec_type), pointer :: molecules(:), mol
        class(molec_type), pointer :: H2, CN, H2p, CO, N2
        type(ion_pointer_t) :: ions(4)
        class(ion_type), pointer :: HI, CI, NI, OI, ion, ioni, ionj
        class(ion_type), pointer :: Hm
        integer :: i, j, nx, ny, nz, ix, iy, iz, np1, iter, nions
        real(dp) :: new_diff
        real(dp), dimension(:,:,:), allocatable :: diff, dis
        real(dp), dimension(:,:,:,:)  , allocatable, target :: Fvec
        real(dp), dimension(:,:,:,:,:), allocatable, target :: Amat
        real(dp), dimension(:), allocatable :: Fold, Fnew
        logical, dimension(:,:,:), allocatable :: tmolim
        integer :: itimer=0
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%molecular_equilibrium', itimer=itimer)

        nx = size(temp, 1)
        ny = size(temp, 2)
        nz = size(temp, 3)

        allocate(dis, mold=temp)
        do i=1, self%nmolec
            mol => self%molecules(i)
            call mol%allocate(shape(temp))
            mol%part = 1.0 / mol%irwin_partf(self%temp)
            dis = mol%get_dis(self%teta)
            if (irwin_dis .and. size(mol%de)==1 .and. mol%nions==2) then
                call mol%get_dis2(self%temp, mol%dis)
            else
               mol%dis = mol%get_dis(self%teta)
            endif
        end do

        call self%ion_by_id('H I', HI)
        call self%ion_by_id('C I', CI)
        call self%ion_by_id('N I', NI)
        call self%ion_by_id('O I', OI)
        Hm => HI%below

        molecules => self%molecules

        call self%molecule_by_id('H2' , H2)
        call self%molecule_by_id('N2' , N2)
        call self%molecule_by_id('CN' , CN)
        call self%molecule_by_id('CO' , CO)
        call self%molecule_by_id('H2+', H2p)
        
        ions(1)%p => HI
        ions(2)%p => CI
        ions(3)%p => NI
        ions(4)%p => OI
        nions = size(ions)
        np1 = nions + 1

        allocate(r, s, diff, q, p, b, pefe, mold=temp)
        allocate(Fvec(np1, nx, ny, nz))
        allocate(Amat(np1, np1, nx, ny, nz))
        allocate(Fold(np1), Fnew(np1))
        allocate(tmolim(nx,ny,nz))

        if (allocated(Hm%pk)) then
            deallocate(Hm%pk, HI%fil, OI%fil, CI%fil, NI%fil)
        endif

        allocate(Hm%pk, HI%fil, OI%fil, CI%fil, NI%fil, mold=temp)

        tmolim = temp < self%tmolim
        Hm%pk = 1/Hm%saha

        where (tmolim)
            fe = xnecno
            HI%fil = 1 + HI%saha + Hm%pk
            b = H2%dis * pe/fe * 2
            q = - HI%abund/b
            p = HI%fil/b
            HI%f = solve_quadratic(p, q)

            fe = HI%f * (HI%saha - Hm%pk) + xnecno

            where (fe <= 0)
                r = Hm%saha * pe * H2%dis
                s = HI%fil - r
                p = (HI%fil*xnecno - Hm%pk - 2*xnecno*r) / s
                q = r*xnecno*xnecno/s
                fe = solve_quadratic(p, q)
                HI%f = (xnecno - fe)*Hm%saha
            end where
        end where

        pefe = pe/fe

        do i=1, self%nmolec
            !partial electron pressures due to each molecule
            where (tmolim) molecules(i)%f = molecules(i)%calc_f(pefe)
        end do

        where (tmolim)
            b = CO%dis * pefe
            OI%fil = 1 + OI%saha + OI%f_of_molecules_with(HI) / OI%f
            CI%fil = 1 + CI%saha + CI%f_of_molecules_with(HI) / CI%f

            p = (CI%abund-OI%abund)/OI%fil + CI%fil/b
            q = - OI%abund/OI%fil * CI%fil/b
            OI%f = solve_quadratic(p, q)
            CI%f = CI%abund / (CI%fil + b*OI%f)

            CN%f = CN%calc_f(pefe) ! has to be recalculated since CI%f changed
            
            NI%fil = 1 + NI%saha + (NI%f_of_molecules_with(HI) + NI%f_of_molecules_with(CI)) / Ni%f
            b = 2*N2%dis * pefe
            q = -NI%abund/b
            p = NI%fil/b
            NI%f = solve_quadratic(p, q)
            where (NI%f < 0) NI%f = 0
        end where

        diff = 0
        do iter=1, 100
            Amat = 0
            do i=1, self%nmolec
                where (tmolim) molecules(i)%f = molecules(i)%calc_f(pefe)
            end do

            do i=1, nions
                ioni => ions(i)%p

                where (tmolim) 
                    Fvec(i,:,:,:) = ioni%f_of_own_molecules()
                    Amat(i, np1,:,:,:) = - Fvec(i,:,:,:) / fe
                end where

                if (i==1) then
                    where (tmolim) 
                        Amat(i,  i,:,:,:) = 1 + HI%saha + Hm%pk
                        Amat(np1,i,:,:,:) = HI%saha + 2*H2p%f /(HI%f * pefe) - Hm%pk
                        Fvec(np1,:,:,:) = HI%f * (HI%saha + H2p%f - Hm%pk) - fe + xnecno
                    end where
                else
                    where (tmolim) 
                        Amat(i,  i,:,:,:) = 1 + ioni%saha
                        Amat(np1,i,:,:,:) = ioni%saha
                        Fvec(np1,:,:,:) = Fvec(np1,:,:,:) + ioni%saha * ioni%f
                    end where
                endif

                where (tmolim) Fvec(i,:,:,:) = Fvec(i,:,:,:) + Amat(i,i,:,:,:)*ioni%f - ioni%abund

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

                call lineq_no_pivot(Amat(:,:,ix,iy,iz), Fvec(:,ix,iy,iz), np1)
                Fnew = Fold - Fvec(:,ix,iy,iz)
                diff(ix,iy,iz) = maxval(abs(1-Fold/Fnew))

                do i=1, nions
                    ions(i)%p%f(ix, iy, iz) = Fnew(i)
                end do

                fe(ix,iy,iz) = Fnew(np1)
                pefe(ix,iy,iz) = pe(ix,iy,iz) / fe(ix,iy,iz)
            end do
            end do
            end do

            if (iter > 2 .and. all(diff < 0.001)) exit
        end do

        Hm%f = HI%f * Hm%pk
        do i=1, nions
            ion => ions(i)%p
            do while (associated(ion%above))
                ion%above%f = ion%f * ion%saha
                ion => ion%above
            end do
        end do
        do i=1, self%nmolec
            where (tmolim)
                molecules(i)%f = molecules(i)%calc_f(pefe)
            end where
        end do
        call trace%end(itimer)
    END SUBROUTINE molecular_equilibrium

    !===============================================================================
    SUBROUTINE compute_energy_and_gas_pressure(self, eh, e, pg)
        class(composition_type), target :: self
        class(element_type), pointer :: element
        class(molec_type), pointer :: mol
        class(ion_type), pointer :: ion, HI, Hm
        real(dp), intent(out), dimension(:,:,:) :: eh, pg
        real(dp), intent(inout), dimension(:,:,:) :: e
        real(dp), dimension(:,:,:), allocatable :: xerg, ejon, emol
        integer :: i, j
        !---------------------------------------------------------------------------
        allocate(xerg, ejon, emol, mold=pg)
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('H -', Hm)

        ejon = 0
        pg = self%ne
        eh = - (HI%xion - self%dxi)*(HI%n + Hm%n) - (Hm%xion - 2*self%dxi)*Hm%n

        do i=1, self%nel
            element => self%elements(i)
            do j=1, element%nion
                ion => element%ions(j)
                pg = pg + ion%n
                if (i>1 .and. associated(ion%above)) then
                    xerg = ion%xion - self%dxi*j
                    ! energy of ionization (ejon). hydrogen is not included.
                    ejon = ejon + ion%above%n * xerg
                endif
            end do
        end do

        emol = 0
        do i=1, self%nmolec
            mol => self%molecules(i)
            pg = pg + mol%n
            emol = emol + mol%calc_e(self%teta, self%temp, self%dxi)
        end do

        ! The gas pressure is Pg = NkT, where N is (Ne + Ni + Nmol)
        pg = pg * BK * self%temp

        !ehm = (Hm%xion + HI%xion) * Hm%f
        where (.not.ieee_is_nan(E))
            E = EE*(ejon + eh + emol)
            E = (1.5*pg + E) / self%rho
        end where
    END SUBROUTINE compute_energy_and_gas_pressure

    !===============================================================================
    SUBROUTINE molecule_by_id(self, id, molecule, verbose)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(molec_type), pointer, intent(out) :: molecule
        integer :: idx
        logical, optional :: verbose
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        idx = findloc(self%molecules%id, uppercase(id), 1)
        if (idx == 0) then
            molecule => null()
            if (option(verbose)) then
                print *, 'molecule with id=' // id // ' was not found.'
            endif
        else
            molecule => self%molecules(idx)
        endif
        call trace%end ()
    END SUBROUTINE molecule_by_id

    !===============================================================================
    SUBROUTINE element_by_id(self, id, element, verbose)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(element_type), pointer, intent(out) :: element
        integer :: idx
        logical, optional :: verbose
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        idx = findloc(self%elements%id, capitalize(id), 1)

        element => null()
        if (idx == 0) then
            if (option(verbose)) then
                print *, 'Opacity package: element with id = ' // id // ' was not found.'
            endif
        else
            element => self%elements(idx)
        endif
        call trace%end ()
    END SUBROUTINE element_by_id

    !===============================================================================
    SUBROUTINE ion_by_id(self, id, ion, verbose, element_out)
        class(composition_type), target :: self
        character(len=*), intent(in) :: id
        class(ion_type), pointer, intent(out) :: ion
        class(element_type), intent(out), pointer, optional :: element_out
        class(element_type), pointer :: element
        logical, optional :: verbose
        character(3) :: eid, iid
        character(6) :: trim_id
        integer :: idx, iostat
        !---------------------------------------------------------------------------
        call trace%begin ('composition_type%init')
        read(id, *, iostat=iostat) eid, iid
        eid = capitalize(eid)
        iid = uppercase(iid)
        if (iostat /= 0) then
            if (option(verbose)) then
                print *, 'ion with id = ' // id // ' was not found.'
            endif
            ion => null()
        else
            trim_id = trim(eid)//' '//trim(iid)
            call self%element_by_id(eid, element)
            if (associated(element)) then
                idx = findloc(element%ions%id, trim_id, 1)
                if (lbound(element%ions,1) /= 1) then
                    idx = idx + lbound(element%ions,1) - 1
                endif
                if (idx < lbound(element%ions,1) .or. idx>ubound(element%ions,1)) then
                    ion => null()
                    print *, 'ion with id = ' // trim_id // ' was not found.'
                else
                    ion => element%ions(idx)
                endif
            else
                ion => null()
            endif
            if (present(element_out)) then
                element_out => element
            endif
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

    !=========================================================================
    SUBROUTINE compute_grams_per_hydrogen(self)
        class(composition_type), target :: self
        !---------------------------------------------------------------------
        ! Abundances should not be in log units when they enter here!
        ! XMY = GRAMS OF STELLAR MATTER / GRAMS OF HYDROGEN
        !---------------------------------------------------------------------
        self%suma = sum(self%elements%abund)
        self%xmy  = sum(self%elements%abund*self%elements%awgt) / self%elements(1)%awgt
        self%grph = XMH*self%xmy
    END SUBROUTINE compute_grams_per_hydrogen

    !=========================================================================
    SUBROUTINE overwrite_abundance(self, species, abnd, awgt, relative_abundance)
        class(composition_type), target :: self
        class(massive_type), pointer :: species
        class(element_type), pointer :: element
        class(molec_type), pointer :: mol
        class(ion_type), pointer :: ion
        real(dp), intent(in) :: abnd, awgt
        logical, intent(in) :: relative_abundance
        real(dp) :: abund, a
        integer :: j
        !---------------------------------------------------------------------
        select type (species)
        class is (element_type)
            element => species
        class is (molec_type)
            mol => species
            do j=1, mol%nions
                ion => mol%ions(j)%p
                if (ion%id(1:2) /= 'H ') then
                    call self%element_by_id(ion%id(1:2), element, verbose=.false.)
                    exit
                endif
            end do
        class default
            print *, 'Particle is not recognised by opacity package.'
            return
        end select

        abund = log10(element%abund) + 12.0
        if (relative_abundance) then
            a = abund + abnd
        else
            a = abnd
        endif
        if (abund /= a) then
            print '(A, A, A, F5.2, A, F5.2)', 'Changing abundance of ', trim(element%id), ' from ', abund, ' to ', a
            element%abund = 10.0**(a - 12.0)
            !if (.not. associated(mol)) element%awgt = awgt
            call self%compute_grams_per_hydrogen()
        endif
    END SUBROUTINE overwrite_abundance

    !=========================================================================
    SUBROUTINE add_hrest(self, chi, lam)
        class(composition_type), target :: self
        class(ion_type), pointer :: HI
        real(4), dimension(:,:,:,:), intent(inout) :: chi
        real(dp), dimension(:), intent(in) :: lam
        real(dp), dimension(:), allocatable :: lam3
        real(dp) :: hrest, umc, fakt, EK_over_temp, dxi
        real(dp), parameter :: nniv2 = 15.5**2
        integer :: f, x,y,z
        integer :: p(4)=pp
        !---------------------------------------------------------------------
        call self%ion_by_id('H I', HI)
        allocate(lam3, mold=lam)
        lam3 = lam**3
        do z=1, size(self%temp,3)
        do y=1, size(self%temp,2)
        do x=1, size(self%temp,1)
            EK_over_temp = self%EK_over_temp(x,y,z)
            dxi = self%dxi(x,y,z)
            fakt = HI%boltzmann(EK_over_temp, dxi, natural_log=.true.)
            fakt = 2.0898E-26 / (2*RYD*EK_over_temp) * fakt * HI%part(x,y,z) * HI%n(x,y,z)
            if (RYD >= nniv2 * dxi) then
                !NNIV=15
                !XNIV=15.
                !IF(H(1).LT.XNIV)NNIV=INT(H(1))
                !ADDF=EXP(TETA31/((FLOAT(NNIV)+0.5)**2)-UMC)-1.
                !IF(H(1).LT.XNIV+0.5)ADDF=0.
                !HREST(NTP)=EXPJ*ADDF
                hrest = exp(EK_over_temp*(RYD/nniv2 - dxi))
                hrest = fakt * (hrest - 1)
                chi(:,x,y,z) = chi(:,x,y,z) + hrest*lam3
            endif
        end do
        end do
        end do
    END SUBROUTINE add_hrest

    !=========================================================================
    SUBROUTINE electron_scattering(self, chi)
        class(composition_type), target :: self
        real(4), dimension(:,:,:,:), intent(inout) :: chi
        real(dp), parameter :: thomson = 8./3.*PI*(EC**2/(EM*CC**2))**2
        real(dp) :: els
        integer :: f, x, y, z
        !---------------------------------------------------------------------
        do z=1, size(chi, 4)
        do y=1, size(chi, 3)
        do x=1, size(chi, 2)
            !els = 4.8206e-9 * self%pe / self%temp
            els = thomson * self%ne(x,y,z)
            chi(:,x,y,z) = chi(:,x,y,z) + els
        end do
        end do
        end do
    END SUBROUTINE electron_scattering

    !=========================================================================
    SUBROUTINE rayleigh_scattering(self, lam, chi)
        class(composition_type), target :: self
        class(molec_type), pointer :: H2
        class(ion_type), pointer :: HI, HE
        real(dp), dimension(:), intent(in) :: lam
        real(4), dimension(:,:,:,:), intent(inout) :: chi
        real(dp), dimension(:), allocatable :: xray2, HeIscat
        real(dp) :: xr, ll, scalar, h2ray, rayh, rayh2, nHI, nH2, nHE
        logical :: do_H2
        integer :: f, x, y, z
        !---------------------------------------------------------------------
        allocate(xray2, HeIscat,  mold=lam)

        ! rayleigh scattering
        call self%ion_by_id('H I', HI)
        call self%ion_by_id('HE I', HE)
        call self%molecule_by_id('H2', H2)
        do_H2 = associated(H2)
        
        xray2 = 0.
        where (lam > 911.27) xray2 = 1 / max(1026., lam)**2

        if (helium_scattering) then
        ! Scattering by HeI in cm**2/atom, from Bues & Wehrse, 1976, A&A 51, 461
            HeIscat = 0.66520e-24*4.*(500./lam)**4
        else
            HeIscat = 0
        endif

        do z=1, size(chi, 4)
        do y=1, size(chi, 3)
        do x=1, size(chi, 2)
          nHI = HI%n(x,y,z)
          nHE = HE%n(x,y,z)
          if (do_H2) then
            nH2 = H2%n(x,y,z)
          else
            nH2 = 0
          endif

          do f=1, size(lam)
            xr = xray2(f)
            if (xr > 0.) then
              rayh  = (5.799e-13 + xr*(1.422e-6 + xr*2.784)) * nHI
              rayh2 = ( 8.14e-13 + xr*( 1.28e-6 + xr*1.61 )) * nH2
  
              ! add sumsca to sumsca from lines
              chi(f,x,y,z) = chi(f,x,y,z) + (rayh + rayh2)*xr*xr
            endif

            chi(f,x,y,z) = chi(f,x,y,z) + HeIscat(f) * nHE
          end do
        end do
        end do
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
        write(lunit, '(A8, es12.4)') 'dxi=', self%dxi(pp(1),pp(2),pp(3))

        do i=1, self%nel
            call self%elements(i)%print_info(lunit, pp)
        end do
        close(lunit)
    END SUBROUTINE print_ion_table

    !===============================================================================
    SUBROUTINE read_molec_data (self, molec_file)
        class(composition_type) :: self
        class(molec_type), pointer :: mol
        character(len=*) :: molec_file
        character(80) :: text
        character(3) :: id
        integer :: i, lunit, iostat, defix, n_irwin, n_de
        !---------------------------------------------------------------------------
        call os%check_file(molec_file)
        open(file=molec_file, newunit=lunit)
        read(lunit, *, iostat=iostat) text
        do while (iostat==0)
            do while (text(1:1) == '!')
                read(lunit, *, iostat=iostat) text
            end do
            backspace lunit
            read(lunit, *, iostat=iostat) id, defix, n_irwin
            call self%molecule_by_id(id, mol)
            if (associated(mol)) then
                n_de = merge(5, 1, defix==0)
                allocate(mol%de(n_de))
                allocate(mol%irwin(n_irwin))
                read(lunit, *) mol%de
                read(lunit, *) mol%dc ! always 5 numbers
                read(lunit, *) mol%irwin
                mol%irwin(:) = mol%irwin(n_irwin:1:-1)
                read(lunit, *) mol%xion
            endif
            read(lunit, *, iostat=iostat) text
        end do
        close(lunit)
    END SUBROUTINE read_molec_data

    !===============================================================================
    SUBROUTINE check_for_comment(lunit, comment, n_skipped)
        integer :: lunit, iostat, n_skipped
        character(len=1) :: comment, first_char
        !---------------------------------------------------------------------------
        iostat = 0
        n_skipped = 0
        do while (iostat==0)
            read(lunit, *, iostat=iostat) first_char
            if (first_char /= comment) then
                backspace(lunit)
                return
            endif
            n_skipped = n_skipped + 1
        end do
    END SUBROUTINE check_for_comment

    !===============================================================================
    SUBROUTINE print_opacity_contributions(self, nu)
        class(composition_type), target :: self
        class(bftrans_type), pointer :: trans
        class(ion_type), pointer :: HeI
        real(dp), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:,:), allocatable :: stim, tchi
        real(4), dimension(:,:,:,:,:), allocatable :: chi
        real(dp), dimension(:), allocatable :: lam, hnu
        real(4) :: ref
        character(len=1024) :: text
        integer, dimension(:), allocatable :: bin, nel, order
        character(6), dimension(:), allocatable :: key
        integer :: f, x, y, z, i, j, ib, nb, nh, nx, ny, nz, nf, nt, lunit
        logical :: first_time = .true.
        integer, save :: itimer=0
        !---------------------------------------------------------------------------
        call trace%begin('composition_type%print_opacity_contributions', itimer=itimer)
        nx = size(self%temp, 1)
        ny = size(self%temp, 2)
        nz = size(self%temp, 3)
        nf = size(nu, 1)
        nt = size(self%trans)
        allocate(stim(nf,nx,ny,nz))
        allocate(tchi(nf,nx,ny,nz))
        allocate(bin(nt))
        allocate(key(nt))
        allocate(nel(nt), order(nt))

        do i=1, nt
            trans=> self%trans(i)
            if (associated(trans%ion)) then
                nel(i) = findloc(element_ids, trans%ion%id(1:2), 1)
            else
                nel(i) = 0
            endif
            if (nel(i) == 0) nel(i) = nt+1
        end do
        call argsort_int(nel, order)

        nb = 3
        key(1:3) = ['Hrest', 'Escat', 'Rscat']
        do i=1, nt
            trans => self%trans(order(i))
            ib = findloc(key(1:nb), trim(trans%key), 1)
            if (ib == 0) then
                nb = nb + 1
                ib = nb
                key(ib) = trim(trans%key)
            endif
            bin(i) = ib
        end do
        
        allocate(chi(nf, nx, ny, nz, nb))
        call self%ion_by_id('HE I', HeI)
        lam = CC / nu * 1e8
        hnu = h_e * nu
        chi = 0

        do z=1, nz
        do y=1, ny
        do x=1, nx
        do f=1, nf
            stim(f,x,y,z) = 1. - exp(-hnu(f) * self%EK_over_temp(x,y,z))
        end do
        end do
        end do
        end do

        call self%add_hrest(chi(:,:,:,:,1), lam)
        chi(:,:,:,:,1) = chi(:,:,:,:,1) * stim

        if (scattering_absorption) then
            call self%electron_scattering(chi(:,:,:,:,2))
            call self%rayleigh_scattering(lam, chi(:,:,:,:,3))
        end if

        ! add absorption from all tables, which do not include stimulated emission
        do i=1, nt
            trans => self%trans(order(i))
            ib = bin(i)
            tchi = 0
            call trans%add_absorption(tchi, lam, self%temp, self%EK_over_temp, self%teta, &
                                        self%pe, self%dxi, ne=self%ne, nhe=HeI%n, shield=self%shield, sorted=.true.)
            if (.not. trans%istim) tchi = tchi * stim
            chi(:,:,:,:,ib) = chi(:,:,:,:,ib) + tchi
        end do

        call io%lock%set('opacity_contributions')
        open(newunit=lunit, file=trim(io%run_dir)//'opacity_contributions.txt', access='append')
        if (first_time) then
            first_time = .false.
            write(lunit, '(A8, A12, 300(A7))') 'Temp', 'Chi500/Rho', (trim(key(i)), i=1, nb)
        endif

        do z=1, nz
            associate(chi=>chi(1,1,1,z,:), rho=>self%rho(1,1,z), t=>self%temp(1,1,z))
                ref = sum(chi)
                write(text, '(F8.2, ES12.3, 300(F7.2))') t, ref/rho, (chi(i)/ref*100, i=1, nb)
            end associate
            call replace_substring(text, '-0.00', '     ')
            call replace_substring(text,  '0.00',  '    ')
            write(lunit, '(A)') text(1:18+7*nb)
        end do
        close(lunit)
        call io%lock%unset('opacity_contributions')
        call trace%end(itimer)
    END SUBROUTINE print_opacity_contributions
END MODULE composition_mod
