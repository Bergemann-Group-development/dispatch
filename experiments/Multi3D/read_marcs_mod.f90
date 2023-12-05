MODULE read_marcs_mod
  USE io_mod
  USE prec_mod, only: dp
  USE os_mod
  USE op_math_mod, only: interp1d, pchip_mono
  USE mconst
  implicit none

  integer, public :: zfac=4

contains

  !===============================================================================
  SUBROUTINE read_marcs (file, newn, ds, rho, temp, pe, FeH, vmic)
    character(len=*), intent(in) :: file
    integer, intent(in) :: newn 
    real(8), intent(out) :: ds(3), FeH
    real(8), intent(inout) :: vmic
    real(4), intent(out), dimension(:,:,:), pointer, contiguous :: rho, temp, pe
    integer :: lunit, iostat, ndepth, kk, k
    real(8) :: FeH_marcs
    real(4) :: teff, flux, vmic_Marcs, g, rmass, rmetal, alpha, radius, rluminosity
    real(4) :: xx, yy, zz, palfa, pny, py, pbeta, abund(92)
    real(4) :: lgTauR, lgTau5, pg, prad, pturb, xkapr, emu, vconv, fconv, rhox
    real(4), dimension(:), allocatable :: depth, in_temp, in_pe, in_depth, in_rho
    character(512) :: fileID, command
    !-----------------------------------------------------------------------------
    call os%check_file(file)
    open(newunit=lunit, file=file)
    read(lunit,*) fileID
    read(lunit,*) teff                  ! Teff [K].
    read(lunit,*) flux                  ! Flux [erg/cm2/s]
    read(lunit,*) g                     ! Surface gravity [cm/s2]
    read(lunit,*) vmic_Marcs            ! Microturbulence parameter [km/s]
    read(lunit,*) rmass                 ! No mass for plane-parallel models
    read(lunit,*) FeH_marcs, alpha      ! Metallicity [Fe/H] and [alpha/Fe]
    read(lunit,*) radius                ! 1 cm radius for plane-parallel models
    read(lunit,*) rluminosity           ! Luminosity [Lsun] FOR A RADIUS OF 1 cm!
    read(lunit,*) palfa, pny, py, pbeta ! Convection parameters: alpha, nu, y and beta
    read(lunit,*) xx, yy, zz            ! X, Y and Z, 12C/13C=89 (=solar)
    read(lunit,*) ! skip comment
    read(lunit,*) abund
    read(lunit,*) ndepth
    read(lunit,*) ! skip comment
    read(lunit,*) ! skip comment

    if (FeH /= FeH_marcs) then
      FeH = FeH_marcs
      write(stdout, '(A, F6.2)') 'FeH was overwritten by Marcs model: [Fe/H] =', FeH
    endif
    if (vmic < 0) then
      vmic = vmic_Marcs
      write(stdout, '(A, F6.2)') 'Vmic was overwritten by Marcs model: Vmic =', vmic
    endif
    
    allocate(in_depth(ndepth))
    allocate(in_temp (ndepth))
    allocate(in_rho  (ndepth))
    allocate(in_pe   (ndepth))
    
    do kk=1,ndepth
      read(lunit,*) k, lgTauR, lgTau5, in_depth(k), in_temp(k), in_pe(k), pg, prad, pturb
      if(kk /= k) stop 'model and loop indices differ'
    end do
    
    read(lunit,*) ! skip comment
    
    do kk=1,ndepth
      read(lunit,*) k, lgTauR, xkapr, in_rho(k), emu, vconv, fconv, rhox
      if(kk /= k) stop 'model and loop indices differ'
    end do
    close(lunit)

    allocate(depth(newn))
    allocate(temp (1,1,newn))
    allocate(rho  (1,1,newn))
    allocate(pe   (1,1,newn))

    ds(1:2) = 1d10
    ds(3) = (in_depth(ndepth) - in_depth(1)) / (newn - 1)
    depth = [( k*ds(3) + in_depth(1), k=0, newn-1 )]

    !temp(1,1,:) = interp1d(in_depth, in_temp, depth, order=1, sorted=.true.)
    !rho (1,1,:) = interp1d(in_depth, in_rho , depth, order=1, sorted=.true.)
    !pe  (1,1,:) = interp1d(in_depth, in_pe  , depth, order=1, sorted=.true.)

    call pchip_mono(in_depth, in_temp, depth, temp(1,1,:), sorted=.true.)
    call pchip_mono(in_depth, in_rho , depth, rho (1,1,:), sorted=.true.)
    call pchip_mono(in_depth, in_pe  , depth, pe  (1,1,:), sorted=.true.)

    write(command,*) 'rm -f ', trim(io%run_dir), 'atmos.marcs'
    call execute_command_line (command, exitstat=iostat)
    open(newunit=lunit, file=trim(io%run_dir) // 'atmos.marcs')
    do k=1, newn
      write(lunit, '(E17.8, E17.8, E17.8, F8.1, F8.1)') -depth(k)*1e-5, temp(1,1,k), pe(1,1,k)/(BK*temp(1,1,k)), 0., 0.
    end do
    close(lunit)
  END SUBROUTINE read_marcs
END MODULE read_marcs_mod