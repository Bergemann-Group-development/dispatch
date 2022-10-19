program main
    USE composition_mod
    USE atom_mod
    USE op_math_mod

    implicit none

    real(dp), dimension(1,1,1) :: temp, ne, rho
    real(dp), dimension(:,:,:), allocatable :: pe
    real(dp), dimension(:,:,:,:), allocatable :: Jnu, chi, eta
    integer :: i, j, dk, k
    type(ion_type), pointer :: ion
    type(element_type), pointer :: element
    type(atom_type) :: atom
    !type(composition_type) :: composition
    character(len=80) :: absdat_file, abund_file, atom_file
    real(dp), dimension(:), allocatable :: nu, xx, yy

    !test = [(i*i,i=1,5)]


    absdat_file = '../mul23/input/absdat'
    abund_file = '../mul23/input/abund'
    atom_file = '../mul23/input/atoms/atom.ba06'

    call composition%init(absdat_file, abund_file)

    temp = 9872.90039062500! 4072.5! 3.8751E+03
    ne = 6.153502517126205E+015! 39372336480.2160! 9.6924E+09
    rho = 3.195600032806396E-007! 1.045600045472383E-009! 2.2536249305856074e-10

    !allocate(pe, mold=temp)

    !pe = ne*BK*temp
    !pe = ne*1.38053E-16*temp

    allocate(Jnu(1,1,1,1))
    allocate(nu(1))

    !nu = cc / (5000*1e-8)
    nu = 461265747282872.
    Jnu = 0

    call composition%compute_abundances(temp, ne, rho)
    call composition%absorption(nu, temp, rho, Jnu, chi, eta)

    print *, 'chi = ', chi, 'eta', eta

    !do i=1, composition%nel
    !    element => composition%elements(i)
!
    !    select case (trim(element%id))
    !    case('H', 'C', 'N', 'O')
    !        print *, element%id , element%f * (1 + element%ions(1)%saha) * composition%toth
    !    end select
    !end do
    !call composition%element_by_id('BA', element)
    !print *, element%abund * composition%toth

    !ion =>composition%elements(1)%ions(1)

    !print *, pe, temp, pe / (1.38053E-16*temp)
    !print *, ion%levels%xl
    !print *, ion%levels(1)%fischel_sparks(1, temp, ne)

    !Scall atom%read(atom_file)

    !element => composition%elements(2)
    !ion => element%ions(2)


    !test = (/ ( 1*(i-0.5d0), i=1, 10 ) /)

    !print *, test
    !call element%ion_abundances(temp,ne,n)
    
    !print *, ion%id, ion%n

    !print *, element%ions(2)%levels%g2
    !print *, element%ions(2)%g0
    !print *, element%ions(2)%levels%fischel_sparks(temp, ne)
    !print *, element%ions(2)%levels%boltzmann(temp)

end program main