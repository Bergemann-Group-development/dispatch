module linelist_mod
    use io_mod, only: alloc_lu, messg
    use mconst
    use parallel_mod, only: xs,xe,ys,ye,zs,ze,zr,fs,fe,zs1,ze1,update_1_z
    use atom_var
    use mem_mod
    implicit none

    ! energy level type, can be filled with populations
    type level_lltype
        real(8) :: energy,ion_stage,g
    end type level_lltype

    ! Background line type for linelist.
    ! Extends the atom line type for turbospectrum line list
    ! The important parameters in the line type are filled.
    type,extends(line_type) :: line_lltype
        integer            :: ion_stage,lines_per_ion
        character(len=20)  :: mass_id
        character(len=20)  :: ion
        real(8)            :: wavelength,ion_pot_lo,loggf, &
                              fdamp,g_hi,raddamp,ent8,ent9, &
                              ch_lo,ch_hi,blue,red
        real(8)            :: stark_broad
        character          :: level_lo,level_hi
        character(len=150) :: level_info
        integer            :: active_lo,active_hi
        integer            :: il0,il1,ion_id
        contains 
            procedure :: get_boltzmann_factor
            procedure :: damping
    end type line_lltype

    ! ion type
    type ion_lltype
        ! ion name
        character(len=20)  :: name
        integer            :: ion_stage,n_levels
        real(8)            :: tol,ionization_energy,atomic_weight
        ! absdat id
        integer :: ele_id
        ! number of lines for this ions
        integer :: lines_per_ion
        ! Number of active lines
        integer,allocatable :: n_active(:)
        ! pointer to the list of lines
        type(line_lltype),pointer      :: lines(:)
        type(level_lltype),allocatable :: levels(:)
        integer,allocatable            :: i_active_lines(:,:)
        ! partition functions
        real(8),allocatable :: partition_function(:,:,:)
        contains
            procedure :: collect_levels
            procedure :: init_levels
            procedure :: init_part_func
            procedure :: init_active_lines
            procedure :: init_n_active
            procedure :: collect_partition_function
            procedure :: active_lines_of_ion
    end type ion_lltype

    ! element type
    type element_lltype
        ! Element type for element in linelist
        character(len=6) :: name
        ! ions
        integer :: n_ions
        ! pointer to the list of ions
        type(ion_lltype), pointer :: ions(:)
        contains
            procedure:: loc => loc_el
    end type element_lltype

    ! linelist_lltype type containing all information about 
    ! all ions in the Turbospectrum line list.
    type linelist_lltype
        integer :: n_lines,n_ions,n_elements     ! Number of lines, ions and elements
        integer :: n_io,n_io_tot                 ! Dummy integer (can be replaced in the future)
        real(8),allocatable :: frequency_grid(:) ! Frequency grid special for the background lines [Hz]
        integer :: grid_size                     ! Size of the Frequency grid
        logical :: active=.false.                ! Weather or not the linelist has been activated on a final grid
        character(len=2) :: grid_method = 'ex'   ! Method of frequency grid creation
        integer :: points_per_line = 9           ! Number of points per line
        real(8) :: line_influence  = 2.5         ! Size of the line influence region
        real(8) :: window(2)  = -99.0            ! Wavelength window (Angstrom) of increased precision
        integer :: cut_window = 1                ! Flag for the way of window cutting
        real(8) :: wvl_step   = 0.005            ! Minimal region that should only contain one point (Angstrom)
        type(line_lltype), allocatable    :: lines(:)
        type(element_lltype), allocatable :: elements(:)
        type(ion_lltype), allocatable     :: ions(:)
        contains 
            procedure :: load_from_file
            procedure :: read_lines_into_linelist
            procedure :: init_array
            procedure :: fill_ele_io_lists
            procedure :: init_eleio
            procedure :: activate
            procedure :: loc => loc_ll
            procedure :: populate_frequency_grid
            procedure :: init_frequency_grid
            procedure :: set_window
            procedure :: shrink_to_window
            procedure :: shrink_freq_grid
            procedure :: set_quadrature
    end type linelist_lltype

    type(linelist_lltype),target :: back_linelist             ! Global linelist object
    real(8),allocatable :: h_density(:,:,:),h2_density(:,:,:) ! Number density for H (needed for damping)

    ! Module procedures (bound and unbound)
    contains 
        ! Uppercase function
        function to_upper(strIn) result(strOut)
            character(len=*), intent(in) :: strIn
            character(len=len(strIn))    :: strOut
            integer :: i,j
            do i = 1, len(strIn)
                j = iachar(strIn(i:i))
                if (j>= iachar("a") .and. j<=iachar("z") ) then
                    strOut(i:i) = achar(iachar(strIn(i:i))-32)
                else
                    strOut(i:i) = strIn(i:i)
                end if
            end do    
        end function to_upper

        ! available elements subroutine
        function available_elements() result(av_ele_names)
            character(len=6),allocatable :: av_ele_names(:)
            character(len=6) :: element(120)
            real(8)          :: abundance_sol(120)
            integer          :: ios
            integer          :: lun,i
            ! Read the abund file
            call alloc_lu(lun)
            open(lun,file='abund',status='old',action='read',iostat=ios)
            do i = 1,120
                ! Read a line
                read(lun,*,iostat=ios) element(i),abundance_sol(i)
                if (ios.ne.0) exit
            enddo
            ! make the array for output
            allocate(av_ele_names(i-1))
            av_ele_names = element(1:i-1)
            ! close the file
            close(lun)
        end function available_elements

        logical function exclude_element(test_ele,list) result(exclude)
            character(len=*),intent(in) :: test_ele
            character(len=*),intent(in) :: list(:)
            ! Check if it is empty
            if (adjustl(trim(test_ele)).eq.'') then
                exclude = .true.
                return
            endif

            ! Check if is the foreground element
            if (adjustl(trim(test_ele(1:index(test_ele,' ')))).eq.to_upper(adjustl(trim(atom%id)))) then 
                exclude = .true.
                return
            endif

            ! Check if it is available
            if (any(list.eq. &
                        adjustl(trim(test_ele(1:index(test_ele,' ')))))) then
                exclude = .false.
            else
                exclude = .true.
            endif    
        end function exclude_element
        
        ! check if string is present in part of the list
        logical function string_in_list(string,list,subrange) result(found)
            character(len=*) :: string
            character(len=*) :: list(:)
            integer :: subrange,j
            ! Add the element, if not there
            found = .false.
            do j = 1, subrange + 1
                ! Check if it is in the list part that has 
                ! been added already
                if (index(list(j),string).ne.0) then
                    ! it has been found
                    found = .true.
                    exit
                endif
            enddo
        end function string_in_list
        
        subroutine strip_ion_name(name)
            character(len=*),intent(inout) :: name
            integer :: i
    
            name    = adjustl(trim(name))
            ! index of the first 'I'
            do i = index(name,' '),len(name)
                ! Check if it is either a space or an 'I'
                if ((.not.(name(i:i).eq.'I')).and.(.not.(name(i:i).eq.' '))) exit
            enddo
            name = adjustl(trim(name(:i-1)))

        end subroutine strip_ion_name

        ! Read the rest of the file (bound)
        subroutine read_lines_into_linelist(this,path)
            class(linelist_lltype),intent(inout) :: this
            character(len=*),intent(in)        :: path
            character(len=200)                 :: read_line
            integer                            :: ios,i,j,n
            character(20)                      :: mass_id
            integer                            :: ion_stage,lines_per_ion
            character(20)                      :: ion_name,ion_name_up
            integer                            :: lun
            character(len=6),allocatable       :: av_elements(:)
       
            ! check available elements
            av_elements = available_elements()

            ! Open the file
            call alloc_lu(lun)
            open(lun,file=path,status='old',action='read',iostat=ios)
            if (ios.ne.0) stop "Could not read linelist, stoping."

            n = 0 
            do i=1,this%n_io_tot
                ! Read the first line for this ion
                read(lun,*,iostat=ios) mass_id,ion_stage,lines_per_ion
                if (ios.ne.0) exit
                mass_id   = adjustl(trim(mass_id))
                
                ! Read the second line
                read(lun,*,iostat=ios) ion_name
                call strip_ion_name(ion_name)
                
                ! Upper case for checking element
                ion_name_up = to_upper(ion_name)

                ! check if the element should be excluded
                if (exclude_element(ion_name_up,av_elements)) then
                    do j=1,lines_per_ion
                        read(lun,'(A)',iostat=ios) read_line
                    enddo
                    cycle
                endif
                !write(*,*) ion_name, lines_per_ion, n
                ! Read the actual content
                do j=1,lines_per_ion
                    ! Fill one space in the linelist
                    n = n + 1
                    ! Try to read with stark broadening parameter
                    read(lun,*,iostat=ios) this%lines(n)%wavelength,this%lines(n)%ion_pot_lo,&
                                            this%lines(n)%loggf,this%lines(n)%fdamp,this%lines(n)%g_hi,&
                                            this%lines(n)%raddamp,this%lines(n)%stark_broad,&
                                            this%lines(n)%level_lo,this%lines(n)%level_hi,&
                                            this%lines(n)%ent8,this%lines(n)%ent9,this%lines(n)%level_info
                    if (ios.ne.0) then
                        !write(*,*) "Chosen linelist does not contain stark broadening."
                        ! Try the old version without this parameter
                        backspace lun
                        read(lun,*,iostat=ios) this%lines(n)%wavelength,this%lines(n)%ion_pot_lo,&
                                                this%lines(n)%loggf,this%lines(n)%fdamp,this%lines(n)%g_hi,&
                                                this%lines(n)%raddamp,&
                                                this%lines(n)%level_lo,this%lines(n)%level_hi,&
                                                this%lines(n)%ent8,this%lines(n)%ent9,this%lines(n)%level_info
                        ! Stop if nothing has worked
                        if (ios.ne.0) stop "(read_lines_into_linelist): Failed to read a line in the linelist"
                        this%lines(n)%stark_broad = 0.0
                    endif
                    ! Include the common header to each line
                    this%lines(n)%ion_stage         = ion_stage
                    this%lines(n)%lines_per_ion     = lines_per_ion
                    this%lines(n)%mass_id           = mass_id
                    this%lines(n)%ion               = ion_name

                    ! basic properties
                    this%lines(n)%ch_lo = this%lines(n)%ion_pot_lo 
                    this%lines(n)%ch_hi = ((this%lines(n)%ion_pot_lo * ee) & 
                                            + hc/(this%lines(n)%wavelength * 1.0d-8)) / ee
                    this%lines(n)%blue  = this%lines(n)%wavelength - this%line_influence
                    this%lines(n)%red   = this%lines(n)%wavelength + this%line_influence

                    ! Inherited properties from line_type
                    this%lines(n)%nu0 = cc/(this%lines(n)%wavelength * 1.0d-8) ! Freq. [Hz]
                    this%lines(n)%profile_type = 'voigt'
                    this%lines(n)%io           = 1
                enddo
            enddo
            ! Close the file
            close(lun)
        end subroutine read_lines_into_linelist

        ! Read linelist into linelist object (bound)
        subroutine load_from_file(this,path)
            class(linelist_lltype),intent(inout) :: this
            character(len=*),intent(in)          :: path
            integer :: n_lines
            integer :: n_io,n_io_tot
            ! Get the number of lines
            call number_of_lines(path,n_lines,n_io,n_io_tot)

            ! Number of lines in the file
            this%n_lines  = n_lines
            this%n_io     = n_io
            this%n_io_tot = n_io_tot

            ! Allocate the arrays to the right length
            call this%init_array    

            ! Read the content into the linelist_lltype object
            call this%read_lines_into_linelist(path)

        end subroutine load_from_file

        ! Initialize the linelist array (bound)
        subroutine init_array(this)
            class(linelist_lltype),intent(inout) :: this
            if (allocated(this%lines)) deallocate(this%lines)
            allocate(this%lines(this%n_lines))
        end subroutine init_array

        ! Get the number of lines from the linelist
        subroutine number_of_lines(path,n_lines,n_io,n_io_tot)
            character(len=*),intent(in)  :: path 
            integer,intent(out)          :: n_lines,n_io,n_io_tot
            character(len=200)           :: read_line
            integer                      :: ios
            character(20)                :: mass_id
            integer                      :: ion_stage,lines_per_ion,i,lun
            character(len=6),allocatable :: av_elements(:)
            character(len=20)            :: ion_name

            n_lines  = 0
            n_io     = 0
            n_io_tot = 0
            
            ! Available elements
            av_elements = available_elements()

            ! Open the file
            call alloc_lu(lun)
            open(lun,file=path,status='old',action='read',iostat=ios)
            if (ios.ne.0) stop "Could not read Linelist, stoping."
    
            do
                ! read first line
                read(lun,*,iostat=ios) mass_id,ion_stage,lines_per_ion
                if (ios.ne.0) exit

                ! Check if it is in available lines
                read(lun,*,iostat=ios) ion_name
                ion_name = to_upper(ion_name)
                if (.not.exclude_element(ion_name,av_elements)) then       
                    n_lines = n_lines + lines_per_ion
                    n_io    = n_io  + 1
                endif
                n_io_tot = n_io_tot + 1
                ! Read the rest
                do i=1,lines_per_ion 
                    read(lun,'(A)',iostat=ios) read_line
                enddo
            end do
    
            ! Close the file
            close(lun)
        end subroutine number_of_lines

        ! allocate space for elements and ion arrays
        subroutine init_eleio(this)
            class(linelist_lltype) :: this
            if (allocated(this%elements)) deallocate(this%elements)
            allocate(this%elements(this%n_elements))
            if (allocated(this%ions)) deallocate(this%ions)
            allocate(this%ions(this%n_ions))
        end subroutine init_eleio

        subroutine get_iotope_index(element,ions,strict,idx_start,idx_end)
            ! In/Out
            character(len=*),intent(in)  :: element
            character(len=*),intent(in)  :: ions(:)
            logical,intent(in)  :: strict
            integer,intent(out) :: idx_start,idx_end
            ! local variables
            integer :: i
            logical :: found
            character(len=len(element)) :: io_ele
            ! Go through the ion list and see if the element 
            ! exactly equal to the ion name
            idx_start = 1
            idx_end   = 1
            found     = .false.
            do i = 1, size(ions)
                ! If it is equal, there you go
                if (strict) then 
                    if ((element.eq.ions(i)).and.(.not.found)) then
                        idx_start = i
                        found     = .true.
                    else if ((element.ne.ions(i))) then
                        if (found) then
                            idx_end = i -1
                            exit
                        endif
                    else if (found) then
                        idx_end = i
                    endif
                else
                    ! Just the element name before the first I
                    io_ele = adjustl(trim(ions(i)(:index(ions(i)," "))))
                    if ((element.eq.io_ele).and.(.not.found)) then
                        idx_start = i
                        found     = .true.
                    else if ((element.ne.io_ele)) then
                        ! It is not the right ion, did we find it?
                        if (found) then
                            ! Yes we did, so this must be the end
                            idx_end = i -1
                            exit
                        endif
                    else if (found) then
                        idx_end = i
                    endif
                endif
            enddo
            if (found.and.(idx_end.eq.1)) idx_end = idx_start 
            
            ! If the end is still weird, there has been a problem
            ! Because all elements should have been read and ordered, we stop
            if (.not.found) then
                stop "There has been problem when reading the linelist." &
                   //"Please make sure that ions are ordered."
            endif

        end subroutine get_iotope_index

        ! Fill the elements and ion arrays
        subroutine fill_ele_io_lists(this)
            class(linelist_lltype), target :: this
            character(len=6)  :: ele_names(120)
            character(len=20) :: ion_names(this%n_lines)
            integer :: nele,nios,i,idx_start,idx_end,ele_len
            logical :: found
            
            ! initialize
            nele  = 0
            nios = 0
            ele_names = ""
            ion_names = ""
            found = .false.

            ! go through the lines and add elements/ions if they are new
            ! The dummy lists are too long, but will be sliced in the end
            do i = 1, this%n_lines
                ! Add the elements
                ele_len = index(this%lines(i)%ion," ") -1
                if (.not.string_in_list(this%lines(i)%ion(:ele_len),&
                                        ele_names,nele)) then
                    nele = nele + 1
                    ele_names(nele) = adjustl(trim(this%lines(i)%ion(:ele_len)))
                endif
                ! Same for ions
                if (.not.string_in_list(this%lines(i)%ion,&
                                        ion_names,nios)) then
                    nios = nios + 1
                    ion_names(nios) = this%lines(i)%ion
                endif
            enddo     
            
            ! allocate space for the element arrays
            this%n_elements = nele
            this%n_ions     = nios
            call this%init_eleio

            ! Prepare the pointers
            do i = 1,nios
                call get_iotope_index(ion_names(i),this%lines(:)%ion,.true., &
                                       idx_start,idx_end)
                !write(*,*) ion_names(i),idx_start,idx_end
                this%ions(i)%lines => this%lines(idx_start:idx_end)
                this%ions(i)%lines_per_ion = idx_end - idx_start + 1
                this%ions(i)%name      = ion_names(i)
                this%ions(i)%ion_stage = this%lines(idx_start)%ion_stage
            enddo
            do i = 1,nele
                ! Get beginning and end of this element in line list
                call get_iotope_index(ele_names(i),ion_names(:nios),.false., &
                                        idx_start,idx_end)
                ! Each element has different ions
                !write(*,*) this%ions(idx_start:idx_end)%name, ele_names(i)
                this%elements(i)%ions  => this%ions(idx_start:idx_end)
                this%elements(i)%name   = ele_names(i)
                this%elements(i)%n_ions = idx_end - idx_start + 1
            enddo
        end subroutine fill_ele_io_lists

        integer function loc_ll(this,element) result(idx)
            class(linelist_lltype) :: this
            character(len=*) :: element
            integer :: i

            ! Search the elements first
            idx = 0
            do i = 1, this%n_elements
                if (element.eq.trim(this%elements(i)%name)) then
                    idx = i
                    exit
                endif
            enddo
            if (idx.ne.0) return 

            ! Try search the ions
            do i = 1, this%n_ions
                if (element.eq.trim(this%ions(i)%name)) then
                    idx = i
                    exit
                endif
            enddo
            if (idx.eq.0) stop "Requested element could not be found."
        end function loc_ll
        
        integer function loc_el(this,element) result(idx)
            class(element_lltype) :: this
            character(len=*) :: element
            integer :: i
            ! Search the elements first
            idx = 0
            do i = 1, this%n_ions
                if (this%name//" "//element.eq.trim(this%ions(i)%name)) then
                    idx = i
                    exit
                endif
            enddo
            if (idx.eq.0) stop "Requested element could not be found."
        end function loc_el

        ! Turobespectrum version of the linelist reader function
        ! from MULTI3D (unbound, but could be bound, just the call would need to change)
        subroutine inlin_turbospec(linelist_object)
            use mcalin
            use uppsala_orig_var
            implicit none
           
            ! Linelist object to be used
            type(linelist_lltype) :: linelist_object
            integer :: i_ll,i_ion
            
            ! See if it has been loaded yet
            if (.not.allocated(linelist_object%lines)) then
                call init_linelist(linelist_object)
            endif

            ! Use the read data for filling the arrays
            nll = linelist_object%n_lines
            if (nll.gt.mll) stop "There are more lines in the line list than mll."
            do i_ll = 1,nll
                ! source description
                srcell(i_ll) = " "

                ! lower level info, log(gf) is taken as fll and gill is set to 1
                abnill(i_ll) = to_upper(linelist_object%lines(i_ll)%ion)
                
                ! This part is only needed for the old line list treatment
                !ionill(i_ll) = linelist_object%lines(i_ll)%ion_stage
                !chill(i_ll)  = linelist_object%lines(i_ll)%ion_pot_lo * ee/hc
                !gill(i_ll)   = 1.0

                !! upper level info
                !abnjll(i_ll) = to_upper(linelist_object%lines(i_ll)%ion)
                !ionjll(i_ll) = linelist_object%lines(i_ll)%ion_stage
                !chjll(i_ll)  = ((linelist_object%lines(i_ll)%ion_pot_lo * ee) & 
                !                + hc/(linelist_object%lines(i_ll)%wavelength * 1e-8)) / hc
                !gjll(i_ll)   = linelist_object%lines(i_ll)%g_hi

                !! General line info
                !alamll(i_ll) = linelist_object%lines(i_ll)%wavelength
                !bluell(i_ll) = linelist_object%lines(i_ll)%blue
                !redll(i_ll)  = linelist_object%lines(i_ll)%red
                !fll(i_ll)    = 10**(linelist_object%lines(i_ll)%loggf)
                !gall(i_ll)   = linelist_object%lines(i_ll)%raddamp
                !gqll(i_ll)   = linelist_object%lines(i_ll)%stark_broad
                !gwll(i_ll)   = linelist_object%lines(i_ll)%fdamp

                ! Fill the elmcll array with the element ids from ABSDAT
                call elmcll(i_ll)
                linelist_object%lines(i_ll)%ion_id = iellll(i_ll)
            
            do i_ion = 1,linelist_object%n_ions
                linelist_object%ions(i_ion)%ele_id = linelist_object%ions(i_ion)%lines(1)%ion_id
            enddo

            enddo 
        end subroutine inlin_turbospec

        subroutine init_linelist(linelist_object)
            class(linelist_lltype) :: linelist_object

            ! See if it has been loaded yet
            if (allocated(linelist_object%lines)) deallocate(linelist_object%lines)
            
            ! Read the line list into the the object
            call linelist_object%load_from_file('abslin')

            ! Cut the linelist to the spectral window if desired
            if ((linelist_object%window(1).gt.0).and.(linelist_object%cut_window.ne.0)) then
                call linelist_object%shrink_to_window
            endif

            ! Fill the isotopes and elements arrays
            call linelist_object%fill_ele_io_lists
        end subroutine init_linelist
    
        ! subroutine that says which lines are active at which frequency 
        ! and activates the line list for opacity considerations
        subroutine activate(this,nu)
            class(linelist_lltype),intent(inout) :: this
            real(8),intent(in),dimension(:)      :: nu
            integer :: i_loc
            ! Save a mask for every ion to save time later
            do i_loc=1,this%n_ions
                call this%ions(i_loc)%active_lines_of_ion(nu)
            enddo 
            ! Activate the linelist
            this%active = .true.
        end subroutine activate

        subroutine init_active_lines(this,max_active)
            class(ion_lltype) :: this 
            integer,intent(in) :: max_active
            if (allocated(this%i_active_lines)) deallocate(this%i_active_lines)
            allocate(this%i_active_lines(fs:fe,max_active))
        end subroutine init_active_lines

        subroutine init_n_active(this)
            class(ion_lltype) :: this 
            if (allocated(this%n_active)) deallocate(this%n_active)
            allocate(this%n_active(fs:fe))
        end subroutine init_n_active
         
        ! Get an array with active lines for every isotope
        subroutine active_lines_of_ion(this,nu)
            class(ion_lltype),intent(inout)             :: this      ! Ion
            logical,dimension(fs:fe,this%lines_per_ion) :: mask      ! Mask 
            integer                                     :: freq_loc  ! Current frequency idx
            real(8),intent(in),dimension(:)             :: nu        ! Frequency grid [Hz]
            integer :: i,j_loc,max_active

            call this%init_n_active

            do j_loc=fs,fe
                ! Get the active lines
                mask(j_loc,:) = (this%lines(:)%red.ge.(cc/nu(j_loc)*1.0d8)) .and. &
                                (this%lines(:)%blue.le.(cc/nu(j_loc)*1.0d8))
                ! Save the number of active lines
                this%n_active(j_loc) = count(mask(j_loc,:))
            enddo
            
            ! Maximum size
            max_active = maxval(this%n_active)
            
            ! Allocate the array
            call this%init_active_lines(max_active)

            ! Fill the array
            do j_loc=fs,fe
                this%i_active_lines(j_loc,1:this%n_active(j_loc)) = &
                            pack((/ (i, i=1,this%lines_per_ion)  /),mask(j_loc,:))
            enddo
        end subroutine active_lines_of_ion
        
        subroutine init_levels(this)
            class(ion_lltype),intent(inout) :: this
            if (allocated(this%levels)) deallocate(this%levels)
            allocate(this%levels(this%n_levels))
        end subroutine init_levels
        
        subroutine collect_levels(this,tol)
            class(ion_lltype),intent(inout) :: this
            type(level_lltype),dimension(2*this%lines_per_ion) :: dummy_level_list
            real(8),optional :: tol
            real(8) :: e_lo,e_hi,tol_loc
            integer :: i,n,ie
            tol_loc = 1.0d-3
            if (present(tol)) tol_loc = tol

            stop "Levels from background line list is not implemented yet."

            n = 1
            dummy_level_list(:)%energy = - 99.0
            ! Loop through the lines of this ion an collect the levels
            do i = 1,this%lines_per_ion
                e_lo = this%lines(i)%ch_lo
                e_hi = this%lines(i)%ch_hi
                
                if (.not.any(abs(dummy_level_list(:n)%energy - e_lo).lt.tol_loc)) then
                    ! This is a new level
                    dummy_level_list(n)%energy    = e_lo 
                    dummy_level_list(n)%ion_stage = this%lines(i)%ion_stage
                    dummy_level_list(n)%g         = 1.0
                    n = n + 1
                endif

                if (.not.any(abs(dummy_level_list(:n)%energy - e_hi).lt.tol_loc)) then
                    dummy_level_list(n)%energy    = e_hi 
                    dummy_level_list(n)%ion_stage = this%lines(i)%ion_stage
                    dummy_level_list(n)%g         = this%lines(i)%g_hi
                    n = n + 1
                else 
                    ! Maybe we got lucky and get the statistical weight
                    ie = minloc(abs(dummy_level_list(:n)%energy - e_hi),1)
                    if (abs(dummy_level_list(ie)%g - 1).lt.tol_loc) then
                        dummy_level_list(ie)%g = this%lines(i)%g_hi
                    endif
                endif
            enddo
            
            ! Allocate space for the number of levels found
            this%n_levels = n - 1
            this%tol      = tol_loc
            call this%init_levels

            ! Fill the level list
            this%levels           = dummy_level_list(:this%n_levels)
            this%levels(:)%energy = this%levels(:)%energy *ee/hc
            call quicksort_levels(this%levels,1,this%n_levels)

            ! see which levels correspond to which line
            ! TO BE IMPLEMENTED IF NEEDED
            
        end subroutine collect_levels

        ! Quicksort
        recursive subroutine quicksort_levels(a, first, last)
            type(level_lltype)  a(*),x, t
            integer first, last
            integer i, j
            x = a( (first+last) / 2 )
            i = first
            j = last
            do
                do while (a(i)%energy < x%energy)
                    i=i+1
                end do
                do while (x%energy < a(j)%energy)
                    j=j-1
                end do
                if (i >= j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                i=i+1
                j=j-1
            end do
            if (first < i-1) call quicksort_levels(a, first, i-1)
            if (j+1 < last)  call quicksort_levels(a, j+1, last)
        end subroutine quicksort_levels

        subroutine init_part_func(this)
            class(ion_lltype),intent(inout) :: this
            if (allocated(this%partition_function)) then
                call deallocate(this%partition_function,"linelist_mod")
            endif
            call allocate(this%partition_function,xs,xe,ys,ye,zs1,ze1,"linelist_mod")
        end subroutine init_part_func

        subroutine init_h_density
            if (allocated(h_density)) call deallocate(h_density,"linelist_mod")
            call allocate(h_density,xs,xe,ys,ye,zs1,ze1,"linelist_mod")
            if (allocated(h2_density)) call deallocate(h2_density,"linelist_mod")
            call allocate(h2_density,xs,xe,ys,ye,zs1,ze1,"linelist_mod")
        end subroutine init_h_density
        
        ! Read the partition function and save it 
        subroutine collect_partition_function(this)
            use atmos_var, only: temp,rho,ne
            use uppsala_orig_var

            class(ion_lltype),intent(inout) :: this
            integer :: i,j,iz,ie,it
            real(8) :: t_temp(zr),rho_temp(zr),pe_temp(zr),ggpgh
            real(8) :: roin(mt)
            real(8) :: pg,ro,e,hn
            logical :: save_h_density

            ! Variables for the common blocks
            real(8) :: xiong(mel,mj),ai(mel),abundm(mel),xmh,xmy,anjon(mel,mj),part(mel,mj)
            real(8) :: toth(mt),toth2(mt)
            real(8) :: abund(mel),h(5),dxi,f1,f2,f3,f4,f5,xkhm
            real(8) :: fl2(5),parco(mjmax),parq(4*mjmax),shxij(5),tparf(4),&
                        eev,enamn,sumh,xkbol,summ,nj(mel),nel
            real(8) :: tothi(mt),&
                        totch(mt),totco(mt),totcn(mt),&
                        totc(mt),totnit(mt),toto(mt),&
                        bhyd(mt,5),grph_bla
            
            ! This is a common block for jon, there the partition function is returned
            common/ci5/abund,anjon,h,part,dxi,&
                     f1,f2,f3,f4,f5,xkhm,xmh,xmy
            common/cabndm/ abundm
            common/ci1/fl2,parco,parq,shxij,tparf,&
                              xiong,eev,enamn,sumh,xkbol,summ,nj,nel
            common/ci9/ai
            common/cout/ toth,tothi,toth2,&
                     totch,totco,totcn,&
                     totc,totnit,toto,&
                     bhyd,grph_bla

            ! Allocate the array
            call this%init_part_func
            
            ! What element is this
            ie = this%ele_id

            ! save the ionization energy
            this%ionization_energy = xiong(ie,this%ion_stage)
            this%atomic_weight     = ai(ie)*uu

            ! check if H pop has been computed already
            if (allocated(h_density)) then
                save_h_density = .false.
            else
                save_h_density = .true.
                call init_h_density
            endif

            ! Go through the atmosphere column-wise like in absko and save the part. function
            yloop: do j=ys,ye
                xloop: do i=xs,xe
                    ! Temperature, density and electron pressure column
                    t_temp   = temp(i,j,zs:ze)
                    rho_temp = rho(i,j,zs:ze)
                    pe_temp  = ne(i,j,zs:ze)*bk*temp(i,j,zs:ze)
                    
                    ! Partition function & other stuff for this column
                    call tabs(zr,t_temp)

                    do it=1,zr
                        call jon(it,t_temp(it),pe_temp(it),1,pg,ro,e,0,rho_temp,ggpgh)
                        ! from detabs
                        hn=1./(xmh*xmy)
                        ! store the partition function for this ion
                        this%partition_function(i,j,zs+it-1) = abundm(ie) * hn * &
                                                anjon(ie,this%ion_stage) / part(ie,this%ion_stage)
                    enddo
                    ! Save the H density if it has not been saved
                    if (save_h_density) then
                        h_density(i,j,zs:ze)  = tothi(:zr)
                        h2_density(i,j,zs:ze) = toth2(:zr)
                    endif
                enddo xloop
            enddo yloop

            call update_1_z(this%partition_function)
            if (save_h_density) then
                call update_1_z(h_density)
                call update_1_z(h2_density)
            endif
        end subroutine collect_partition_function

        elemental function get_boltzmann_factor(this,t_loc,ch_loc_in) result(blz)
            class(line_lltype),intent(in) :: this
            real(8),intent(in) :: t_loc
            real(8),intent(in),optional ::ch_loc_in
            real(8) :: ch_loc
            real(8) :: blz
            ch_loc = this%ch_lo
            if (present(ch_loc_in)) ch_loc = ch_loc_in
            blz = exp(-e_k * ch_loc/t_loc)    
        end function get_boltzmann_factor
    
        ! Calculate the damping constant for this line at a given point in atmos
        ! Re-calculate them instead of pre-computing to save RAM
        ! If this turns out to be too slow, one can pre-tabulate them instead
        elemental function damping(this,weight,eion,t_loc,ne_loc,tothi_loc,toth2_loc) result(gamma_damp)
            class(line_lltype),intent(in) ::  this
            real(8),intent(in)    :: t_loc,ne_loc,toth2_loc,tothi_loc
            real(8),intent(in)    :: weight,eion
            real(8) :: sigma,alpha,gv,gx,gammaf,scr,c1
            real(8) :: gamma_damp
            real(8) :: doppler_factor_loc

            doppler_factor_loc = sqrt( 2.0*bk/weight * t_loc) / cc
            gamma_damp         = 0.0
            gv                 = 0.0

            ! Van der Waals part
            if (this%fdamp.ge.20.0) then
                gv    = this%fdamp
                sigma = int(gv)*2.80028d-21
                alpha = gv-int(gv)
                gx    = 1.0-alpha*0.5 
                gammaf=1.0+(-.5748646+(.9512363+(-.6998588+ &
                           (.4245549-.1010678*gx)*gx)*gx)*gx)*gx
                gv = (4.0/pi)**(alpha*0.5) * gammaf * 1.e4 * sigma
                c1 = 21173.*(1./1.008+1./(weight/uu))
                scr = sqrt(c1*t_loc)
                scr = gv*((scr/1.e4)**(1.-alpha))*tothi_loc*1.e6*2.
                gamma_damp = gamma_damp + scr
            elseif (this%fdamp.gt.0.0) then
                c1 = 1.283984d-12*this%ion_stage**0.8* &
                     (abs(1./(eion-this%ch_hi)**2-1./(eion-this%ch_lo)**2))**.4
                gamma_damp = gamma_damp + &
                            this%fdamp*8.411*(8.*bk*t_loc/pi*(1./(1.008*uu)+1./weight))**0.3 * &
                            tothi_loc*c1
            elseif (this%fdamp.lt.0.0) then
                ! this is copied from turobospectrum2019 Bpz 05/04-2013
                ! PE: Remark: Replace toth and toth2 + tothe with pressure from RH
                ! + phe(i)/ph(i)*0.41336 +
                gv = 10.0**this%fdamp * (t_loc/10000.)**0.3 * tothi_loc
                gv = gv * (1.0 + toth2_loc/tothi_loc*0.85) ! Missing He part, set it to 0 for now
                gamma_damp = gamma_damp + gv
            endif

            ! stark part
            if(this%stark_broad.ge.0.0) then
                gamma_damp = gamma_damp+this%stark_broad * ne_loc
            else
                ! formula from Gray
                gv = 19.4 + 2./3. * this%stark_broad + log10(ne_loc*bk) + log10(t_loc) / 6.
                gamma_damp = gamma_damp + 10**gv
            endif

            ! radiative part
            gamma_damp = gamma_damp + this%raddamp 
            gamma_damp = gamma_damp / (4.0*pi*(cc/(this%wavelength*1.0d-8))) / doppler_factor_loc
        end function damping

        ! build the line list frequency grid
        subroutine populate_frequency_grid(this)
            class(linelist_lltype),intent(inout) :: this
            integer :: i,j,if_loc,n_loc
            real(8) :: nucen_loc,spacing_loc,nublue_loc,nured_loc
            logical :: limit_window

            if (this%window(1).lt.0.0) then
                limit_window = .false.
            else 
                limit_window = .true.
            endif

            ! If it is not symmetric, add an additional point to make it symmetric
            if (modulo(this%points_per_line-1,2).ne.0) this%points_per_line = this%points_per_line + 1

            ! calculate the number of points that is desired
            this%grid_size = this%n_lines * this%points_per_line
            if (limit_window) then
                do i = 1,this%n_lines
                    if (.not.((this%lines(i)%wavelength.le.this%window(2))&
                    .and.(this%lines(i)%wavelength.ge.this%window(1)))) then
                         this%grid_size = this%grid_size - this%points_per_line
                    endif
                enddo
            endif
            call this%init_frequency_grid

            n_loc = (this%points_per_line-1)/2

            ! Loop through all lines and add frequency points at appropriate locations
            if_loc = 0
            do i = 1,this%n_lines
                if (limit_window .and. &
                    (.not.((this%lines(i)%wavelength.le.this%window(2))&
                    .and.(this%lines(i)%wavelength.ge.this%window(1))))) then
                    cycle
                endif
                ! populate the line from the center to the wings 
                select case (this%grid_method)
                    case ('eq')
                        ! populate evenly from red to blue
                        nublue_loc  = cc/(this%lines(i)%blue*1.0d-8)
                        nured_loc   = cc/(this%lines(i)%red*1.0d-8) 
                        spacing_loc = (nublue_loc - nured_loc)/(this%points_per_line-1)

                        do j = 1,this%points_per_line
                            if_loc = if_loc + 1
                            this%frequency_grid(if_loc) = nured_loc + spacing_loc*(j-1)
                        enddo

                    case ('ex')
                        ! exponentially spaced grid from blue to red
                        nucen_loc  = cc/(this%lines(i)%wavelength*1.0d-8)
                        nublue_loc = cc/(this%lines(i)%blue*1.0d-8)
                        nured_loc  = cc/(this%lines(i)%red*1.0d-8) 

                        ! For the blue part
                        spacing_loc = (log(nublue_loc) - log(nucen_loc))/n_loc
                        this%frequency_grid(if_loc+n_loc+2:if_loc+n_loc+1+n_loc) = &
                                            exp( (/ (log(nucen_loc)+spacing_loc*j,j=1,n_loc) /) )
                        
                        ! Mirror the red part
                        this%frequency_grid(if_loc+1:if_loc+n_loc) = &
                            2*nucen_loc - this%frequency_grid(if_loc+n_loc+2:if_loc+n_loc+2+n_loc)
                            
                        ! the central part
                        this%frequency_grid(if_loc+n_loc+1) = nucen_loc

                        if_loc = if_loc + this%points_per_line

                    case default
                        stop "LINELIST_MOD: The chosen frequency grid method is not available." 
                end select
            enddo

            ! Remove points in crowded regions 
            call this%shrink_freq_grid

        end subroutine populate_frequency_grid

        subroutine init_frequency_grid(this)
            class(linelist_lltype) :: this
            if (allocated(this%frequency_grid)) deallocate(this%frequency_grid)
            allocate(this%frequency_grid(this%grid_size))
        end subroutine init_frequency_grid

        subroutine read_injon
            use uppsala_orig_var
            use io_mod, only: close,alloc_lu,open
            use mclu
            integer iread, iwrit, islask, ireat,ireset, labsda
            common/utput/ iread, iwrit
            common/cfil/ireset,islask,ireat
            labsda = 0

            if(labsda.ne.0) call close(labsda)
            call alloc_lu(labsda)
            call open(labsda,'absdat',1,'old',.false.)
            iwrit=lout
            iread=labsda
            ireat=labsda
            call injon(1)
            call close(labsda)
        end subroutine
    
        subroutine set_window(this,w,cw)
            class(linelist_lltype),intent(inout) :: this
            real(8),dimension(2),intent(in) :: w
            integer,intent(in) :: cw
            this%window     = w
            this%cut_window = cw
        end subroutine set_window
    
        subroutine shrink_to_window(this)
            class(linelist_lltype),intent(inout) :: this
            integer :: i,n_lines_short
            logical,dimension(this%n_lines)  :: mask
            type(line_lltype),dimension(:),allocatable :: short_linelist

            ! go through the read line list and remove every line outside the window
            mask = (this%lines(:)%wavelength.ge.this%window(1)) .and. &
                    (this%lines(:)%wavelength.le.this%window(2))

            ! shrink down the line list
            n_lines_short = count(mask) 
            allocate(short_linelist(n_lines_short))
            short_linelist = pack(this%lines,mask)

            ! Reset the linelist
            this%n_lines = n_lines_short
            call this%init_array
            this%lines   = short_linelist

            deallocate(short_linelist)
        end subroutine shrink_to_window

        subroutine shrink_freq_grid(this)
            class(linelist_lltype),intent(inout) :: this
            real(8),dimension(:),allocatable  :: freq_new
            real(8),dimension(this%grid_size) :: wvl_arr_loc
            real(8) :: curr_point
            logical,dimension(this%grid_size) :: mask
            integer :: i

            print *,this%wvl_step
            wvl_arr_loc = cc/this%frequency_grid * 1e8

            do i = 1,this%grid_size
                if (wvl_arr_loc(i).le.0.0) cycle
                ! check the proximity of this point
                curr_point = wvl_arr_loc(i)
                where (abs(wvl_arr_loc - wvl_arr_loc(i)) .le. this%wvl_step)
                    wvl_arr_loc = -99999999.0
                end where
                wvl_arr_loc(i) = curr_point
            enddo

            ! Save only those that are positive
            mask = wvl_arr_loc.gt.0.0
            this%grid_size = count(mask)

            if (allocated(freq_new)) deallocate(freq_new)
            allocate(freq_new(this%grid_size))
            freq_new = pack(this%frequency_grid,mask)

            ! reallocate grid
            deallocate(this%frequency_grid)
            call this%init_frequency_grid
            this%frequency_grid = freq_new

            ! clean up
            deallocate(freq_new)

        end subroutine shrink_freq_grid

        subroutine set_quadrature(this,nq_background,fmethod_background,wavel_step)
            class(linelist_lltype),intent(inout) :: this
            integer,intent(in) :: nq_background
            real(8),intent(in) :: wavel_step
            character(len=2),intent(in) :: fmethod_background

            back_linelist%points_per_line = nq_background
            back_linelist%grid_method     = fmethod_background
            back_linelist%wvl_step        = wavel_step
        end subroutine set_quadrature
   
    end module linelist_mod