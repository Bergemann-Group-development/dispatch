!-------------------------------------------------------------------------------
MODULE eos_mod

  USE io_mod
  USE mpi_mod
  
  IMPLICIT none

  integer, parameter   :: mfile = 127
  character(len=mfile) :: eostabfile
  character(len=mfile) :: radtabfile
  logical              :: do_eos, do_table, do_cubic_eos
  logical              :: flag_scat

  real(kind=4), allocatable, dimension(:,:)   :: pbot0, dlnpdE_r, dlnpdlnr_E
  integer,      allocatable, dimension(:,:,:) :: irarr_e, iearr_e
  real,         allocatable, dimension(:,:,:) :: pxarr_e, pyarr_e
  integer,      allocatable, dimension(:,:,:) :: irarr_r, itarr_r
  real,         allocatable, dimension(:,:,:) :: pxarr_r, pyarr_r
  
  ! EOS table: parameters
  character(len=8)       :: dateeos, timeeos
  integer                :: NetabE, NrtabE, MaE, NaE, Ntbvar, NelemE
  integer, dimension(20) :: NxoptE, NxrevE
  real(kind=4)           :: lnrMinE, lnrMaxE, lneMinE, lneMaxE
  real(kind=4)           :: ulE, urE, utE
  real(kind=4),     allocatable, dimension(:)       :: abundE
  real(kind=4),     allocatable, dimension(:,:)     :: auxE
  real(kind=4),     allocatable, dimension(:,:,:,:) :: eostab
  character(len=4), allocatable, dimension(:)       :: ielE

  ! Radiation table: parameters
  character(len=8)       :: daterad,timerad
  integer                :: Nlam
  integer                :: NTtabR, NrtabR, MaR, NaR, NelemR
  integer, dimension(20) :: NxoptR, NxrevR
  real(kind=4)           :: lnrMinR, lnrMaxR, lntMinR, lntMaxR
  real(kind=4)           :: ulR, urR, utR
  real(kind=4),     allocatable, dimension(:)       :: abundR
  real(kind=4),     allocatable, dimension(:,:)     :: auxR
  real(kind=4),     allocatable, dimension(:,:,:,:) :: radtab
  character(len=4), allocatable, dimension(:)       :: ielR

  ! Maximum size for EOS and radiation tables
  integer, parameter :: maxtabe = 10000000
  integer, parameter :: maxtabr = 10000000
  integer            :: mtablee
  integer            :: mtabler

  ! Simulation units (cgs): length, time, density
  ! Simulation units (cgs): velocity, energy/mass, pressure, opacity,
  ! integrated Planck function (B*\Delta\lambda)
  real(kind=4) :: uur, uul, uut
  real(kind=4) :: uuv, uuee, uup, uuk, uubp

  ! EOS table units (cgs): velocity, energy/mass, pressure, opacity
  real(kind=4) :: uvE, ueeE, upE, ukE

  ! Radiation table units (cgs): velocity, energy/mass, pressure,
  ! opacity, integrated source function
  real(kind=4) :: uvR, ueeR, upR, ukR, ubpR

  ! EOS table: ln density and ln energy/mass increments
  ! Radiation table: ln density and ln temperature increments
  real(kind=4) :: dlnrE, dlneE
  real(kind=4) :: dlnrR, dlntR

  ! Switch: verbose output for debugging
  ! Switch: restrain density, temperature and energy/mass to within table limits
  logical :: verbose_tab
  logical :: do_limit_tab

  ! EOS and radiation tables: limits to density, energy/mass, and temperature
  real(kind=4) :: eelimMinE, eelimMaxE, rholimMinE, rholimMaxE
  real(kind=4) :: ttlimMinR, ttlimMaxR, rholimMinR, rholimMaxR

  ! Limits to density and internal energy (DEPRECATED)
  real(kind=4) :: eelimE, rholimE 

  ! Constants: Boltzmann constant (cgs)
  real(kind=4), parameter :: kBoltz = 1.306504e-16
  real(kind=4), parameter :: lnkBoltz = alog(kBoltz)

  ! NOTE: do_ionization is part of "params" module in Stagger
  logical :: do_ionization

  
  integer, parameter :: n_lambda=12

  TYPE, PUBLIC :: eos_t
     integer :: n_lambda = n_lambda
     real    :: w_lambda(n_lambda) = 1.0
     real(8) :: gamma=1.2d0       ! required for syntax in some solvers
     logical :: on = .false.
     logical :: do_rt = .false.
   CONTAINS
     procedure :: init => init_eos
     procedure :: lookup
  END TYPE eos_t

  type(eos_t), public :: eos_t

CONTAINS

!-------------------------------------------------------------------------------
  SUBROUTINE init
    
    ! This routine reads the EOS and radiation tables.
    !
    ! The EOS table contains the following variables stored in the
    ! following order:
    !
    ! - Pressure (gas plus radiation)
    ! - Rosseland exctinction coefficient (per unit length)
    ! - Temperature
    ! - Electron number density (OPTIONAL)
    ! - Planck-averaged extinction coefficient
    ! - Standard extinction coefficient at 500 nm
    !
    ! The variables are stored in ln scale and tabulated as a function of
    ! ln energy/mass and ln density
    ! 
    ! The EOS table also includes the derivatives of the ln variables with
    ! respect to ln energy/mass and ln density.
    !
    !
    ! The radiation table contains the following variables:
    !
    ! - 1st wavelength (bin):  ratio of monochromatic (bin) and
    !                          Rosseland opacities
    ! - 1st wavelength (bin):  integrated thermal emission
    !                          (epsilon*B*w_lambda)
    ! - 1st wavelength (bin):  destruction probability (epsilon)
    !
    ! - 2nd wavelength (bin):  same as above...
    ! ...
    ! ...
    ! last 3 entries: Rosseland extinction coefficient (per unit length)
    !                 Planck-averaged extinction coefficient,
    !                 extinction coefficient at 500 nm
    !
    ! All variables are stored in ln scale and tabulated as a function of
    ! ln temperature and ln density
    
    IMPLICIT none

    ! ID String
    character(len=32) :: IdStr = 'init'
    character(len=72) :: hline = repeat('-',72)
    
    ! EOS and radiation table file units
    integer :: lunE = 12
    integer :: lunR = 13

    ! String arrays for printout (logging)
    character(len=32), dimension (10) :: sopctype = (/ ('', i=1,10) /)
    character(len=32), dimension (4)  :: ssrctype = (/ ('', i=1, 4) /)
    character(len=32), dimension (4)  :: stautype = (/ ('', i=1, 4) /)
    character(len=32), dimension (4)  :: sxtrtype = (/ ('', i=1, 4) /)

    character(len=16) :: jtype
    character(len=8)  :: linedat

    integer :: i,j, rank

    logical :: flag_vanreg
    logical :: flag_bins
    logical :: flag_odf


    
    ! TO DO: change eqofst --> eos_params ?
    namelist /eqofst/ &
       do_eos, do_ionization, do_table, eostabfile, radtabfile, &
       do_cubic_eos, do_limit_tab, verbose_tab, uur, uul, uut, &
       eelimMinE, eelimMaxE, rholimMinE, rholimMaxE, &
       ttlimMinR, ttlimMaxR, rholimMinR, rholimMaxR, &
       eelimE, rholimE, do_rt

    ! Default values
    do_eos        = .true.
    do_rt         = .false.
    do_ionization = .true.
    do_table      = .false.
    eostabfile    = 'EOSrhoe.tab'
    radtabfile    = 'kaprhoT.tab'
    do_cubic_eos  = .true.
    do_limit_tab  = .false.
    verbose_tab   = .false.
    uur           = 1.0e-7
    uul           = 1.0e8
    uut           = 1.0e2
    eelimMinE     = 1.0e-4
    eelimMaxE     = 1.0e4
    rholimMinE    = 1.0e-9
    rholimMaxE    = 1.0e9
    ttlimMinR     = 1.0e1
    ttlimMaxR     = 1.0e7
    rholimMinR    = 1.0e-9
    rholimMaxR    = 1.0e9
    ! NOTE: eelimE, rholimE kept for backward compatibility
    eelimE        = 0.0       ! DEPRECATED
    rholimE       = 0.0       ! DEPRECATED

    ! Read namelist
    rewind (io%input); read (io%input, eqofst)

    ! Simulation units, derived (cgs)
    uuv  = uul / uut         ! velocity 
    uuee = uuv * uuv         ! energy/mass
    uup  = uur * uuee        ! pressure
    uuk  = 1.  / uul / uur   ! opacity
    uubp = uur * uuee * uuv  ! integrated Planck function ! [density] * [velocity]^3


    ! EOS table
    
    ! Read EOS table
    open (lunE, file=eostabfile, form='UNFORMATTED', status='OLD')
    read (lunE) NetabE, NrtabE, MaE, NaE, Ntbvar, NelemE, NxoptE, NxrevE
    allocate ( auxE(MaE, NaE) )
    allocate ( ielE(NelemE) )
    allocate ( abundE(NelemE) )
    allocate ( eostab(NetabE, NrtabE, 3, Ntbvar) )
    read (lunE) dateeos, timeeos, &
         lnrMinE, lnrMaxE, &
         lneMinE, lnEmaxE, &
         ulE, urE, utE, &
         iele(1:NelemE), &
         abundE(1:NelemE), &
         auxE(1:MaE,1:NaE)
    read (lunE) eostab(1:Netabe,1:Nrtabe,1:3,1:Ntbvar)
    close(lune) 

    ! EOS table units, derived:
    uvE  = ulE / utE         ! velocity
    ueeE = uvE * uvE         ! energy/mass
    upE  = urE * ueeE        ! pressure
    ukE  = 1.  / ulE / urE   ! opacity
    
    ! EOS table boundaries: convert to simulation units
    lneMaxE = lneMaxE + log(ueeE/uuee)
    lneMinE = lneMinE + log(ueeE/uuee)
    lnrMaxE = lnrMaxE + log(urE/uur)
    lnrMinE = lnrMinE + log(urE/uur)

    ! EOS table: density and internal energy increments
    dlnrE = (lnrMaxE-lnrMinE)/(real(NrtabE)-1.0)
    dlneE = (lneMaxE-lneMinE)/(real(NetabE)-1.0)

    ! EOS table: lower and upper limits to ln density and ln energy/mass
    eelimMinE  = max(eelimMinE, exp(lneMinE))
    eelimMaxE  = min(eelimMaxE, exp(lneMaxE))
    rholimMinE = max(rholimMinE,exp(lnrMinE))
    rholimMaxE = min(rholimMaxE,exp(lnrMaxE))

    ! DEPRECATED: eelimE and rholimE: lower limits to energy/mass and density
    if ( eelimE.ne.0.0 .or. rholimE.ne.0.0)  then
       print *, trim(IdStr)//': WARNING: eelimE and rholimE DEPRECATED. '// &
            'Please use eelimMinE and rholimMinE instead.'
       eelimMinE  = max(eelimMinE,  eelimE )
       rholimMinE = max(rholimMinE, rholimE)
    end if

    ! EOS table: print effective limits:
    if (do_limit_tab) then
       print '(A)', hline
       print '(A)', trim(IdStr)//&
            ': EOS table: effective lower and upper limits for ee and rho: '
       print '(3X,A,G10.2)', 'EELIMMINE  = ', eelimMinE
       print '(3X,A,G10.2)', 'EELIMMAXE  = ', eelimMaxE
       print '(3X,A,G10.2)', 'RHOLIMMINE = ', rholimMinE
       print '(3X,A,G10.2)', 'RHOLIMMAXE = ', rholimMaxE
       print '(A)', hline
    end if

    
    ! Radiation table

    ! Read radiation table
    open (lunR, file=radtabfile, form='UNFORMATTED', status='OLD')
    read (lunR) NTtabr, Nrtabr, Mar, Nar, Nlam, Nelemr, Nxoptr, Nxrevr
    allocate ( auxR(MaR, NaR) )
    allocate ( ielR(NelemR) )
    allocate ( abundR(NelemR) )
    allocate ( radtab(NTtabR,NrtabR,3,Nlam+1) )
    read (lunR) daterad, timerad, &
         lnrMinR, lnrMaxR, &
         lntMinR, lntMaxR, &
         ulR, urR, utR, &
         ielR(1:NelemR),&
         abundR(1:NelemR), auxR(1:MaR,1:NaR)
    do j=1,Nlam+1
       read (lunR) radtab( 1:NTtabR, 1:NrtabR, 1:3, j )
    enddo
    close(lunR)

    ! Radiation table, derived units:
    uvR  = ulR / utR           ! velocity
    ueeR = uvR * uvR           ! energy/mass
    upR  = urR * ueeR          ! pressure
    ukR  = 1.  / ulR  / urR    ! opacity
    ubpR = urR * ueeR * uvR    ! integrated Planck (B*\Delta\lambda) [density] * [velocity^3]

    ! Radiation table: ln density boundaries, convert to simulation units
    lnrMaxR = lnrMaxR + log(urR/uur)
    lnrMinR = lnrMinR + log(urR/uur)

    ! Radiation table: ln density and ln temperature steps
    dlnrR = (lnrMaxR-lnrMinR)/(real(NrtabR)-1.0)
    dlntR = (lntMaxR-lntMinR)/(real(NttabR)-1.0)

    
    ! EOS table: convert variables from table to simulation units
    ! 1: Pressure
    ! 2: Rosseland extinction coefficient:
    eostab(1:NetabE,1:NrtabE,1,1) = eostab(1:NetabE,1:NrtabE,1,1) + log(upE/uup)
    eostab(1:NetabE,1:NrtabE,1,2) = eostab(1:NetabE,1:NrtabE,1,2) - log(ulE/uul)
    
    select case(Ntbvar)
    case(5)
       ! 4: Planck mean extinction coefficient
       ! 5: Standard extinction coefficient at 500 nm
       eostab(1:NetabE,1:NrtabE,1,4) = eostab(1:NetabE,1:NrtabE,1,4) - log(ulE/uul)
       eostab(1:NetabE,1:NrtabE,1,5) = eostab(1:NetabE,1:NrtabE,1,5) - log(ulE/uul)
    case(6)
       ! 4: Electron number density
       ! 5: Planck mean extinction coefficient
       ! 6: Standard extinction coefficient at 500 nm
       eostab(1:NetabE,1:NrtabE,1,4) = eostab(1:NetabE,1:NrtabE,1,4) - log(ulE/uul)*3.0
       eostab(1:NetabE,1:NrtabE,1,5) = eostab(1:NetabE,1:NrtabE,1,5) - log(ulE/uul)
       eostab(1:NetabE,1:NrtabE,1,6) = eostab(1:NetabE,1:NrtabE,1,6) - log(ulE/uul)       
       ! TO DO: ensure that the ln electron number density does not exceed MAXEXPONENT in simulation units
    case default
       print '(A)', trim(IdStr)//': Unsupported table format!'
       stop
    end select

    ! Radiation table, convert from table to simulation units
   
    ! Integrated thermal emissions (\epsilon*B*\Delta\lambda)
    do j=1,Nlam
       radtab(1:NTtabR,1:NrtabR,2,j) = radtab(1:NTtabR,1:NrtabR,2,j) + log(ubpR/uubp)
    end do

    ! Rosseland, Planck-averaged, and 500 nm extinction coefficients
    ! NOTE: from revision 73 onwards, the table contains the extinction coefficients,
    ! not the opacities
    radtab(1:NTtabR,1:NrtabR,1:3,Nlam+1) = radtab(1:NTtabR,1:NrtabR,1:3,Nlam+1) - log(ulR/uul) 

    ! EOS and radiation tables: Check table size
    mtableE = NetabE * NrtabE * 3 * Ntbvar
    mtableR = NTtabR * NrtabR * 3 * (Nlam+1)
    if ( mtableE .gt. maxTabE ) then
       ! TO DO: should it be printed by master only?       
       print '(A)', trim(IdStr)//': EOS table exceeds size limits.'
       stop
    end if
    if ( mtableR .gt. maxTabR ) then
       ! TO DO: should it be printed by master only?       
       print '(A)',trim(IdStr)//': Radiation table exceeds size limits.'
       stop
    end if

    
  ! EOS and radiation tables: basic information

  ! flag_scat: consistent treatment of scattering (t) ?
  ! flag_bins: bins (t) or wavelengths (f) ?
  ! flag_odf:  ODFs (t) or opacity sampling (f) ?
  ! flag_vanreg: van Regemerter approximation for line scattering (t) ?
  select case(Nxoptr(2))
  case(1)
     flag_scat = .false.
  case(2)
     flag_scat = .true.
  case default
     ! TO DO: should it be printed by master only?       
     print '(A)', trim(IdStr)//': Unsupported source function type. STOP.'
     stop
  end select
  flag_odf    = Nxrevr(4).lt.70
  flag_bins   = Nxoptr(1).lt.6
  flag_vanreg = flag_scat.and.(Nxoptr(3).ne.0)

  if (flag_odf) then 
     linedat = 'odf.dat'
  else
     linedat = 'OS.dat'
  end if

  if (flag_bins) then 
     jtype='bins'
  else
     jtype='wavelengths'
  end if

  sopctype(1) = 'Opacity binning, scaled'
  sopctype(2) = 'Opacity binning'
  sopctype(6) = 'Opacity sampling'

  ssrctype(1) = 'Pure Planck function'
  ssrctype(2) = 'Scattering term included'

  stautype(1) = 'monochromatic opacity'
  stautype(2) = 'Rosseland optical depth'

  sxtrtype(1) = 'xcorr'
  sxtrtype(2) = 'Jlambda'
  sxtrtype(3) = 'T-scaled Jlambda'


  ! Print out information about EOS and radiation tables
  if (master) then
     print '(A)', hline
     print '(2A)', 'init_eos: basic information about EOS and RAD tables'
     print '(2A)',       ' (eos) table file           = ', eostabfile
     print '(A,2A10)',   ' (eos) date and time        = ', dateeos,timeeos
     print '(A,I4)',     ' (eos) SVN Rev              = ', NxrevE(2)
     print '(A,I2)',     ' (eos) Number of entries    = ', Ntbvar
     print '(A)', hl
     print '(2A)',       ' (rad) table file           = ', radtabfile
     print '(A,2A10)',   ' (rad) date and time        = ', daterad, timerad
     print '(A,I4)',     ' (rad) SVN Rev              = ', Nxrevr(2)
     print '(2A)',       ' (rad) opacity table type   = ', sopctype(Nxoptr(1))
     print '(2A)',       ' (rad) source function type = ', ssrctype(Nxoptr(2))
     print '(A,A10,A,I6)',' (rad) number of ', jtype,' = ', Nlam

     if (flag_bins) then     
        print '(A,A24,A)', &
             ' (rad) bin membership according to value of ', &
             stautype(Nxoptr(4)), ' where tau_lambda=1'
        print '(2A)',&
             ' (rad) type of extrapolation from 1D rad transf solution = ', &
             sxtrtype(Nxoptr(5))
        print '(A,F6.2,A)', &
             ' (rad) Bridging function    = exp(',&
             (-real(Nxoptr(6))*0.01),'*tau).'
     end if

     if (flag_scat) then
        if (flag_vanreg) then
           print '(A,F6.2)', &
                ' (rad) van Regemorter approx of line scattering; factor=', &
                (real(Nxoptr(3))*0.01)
        else
           print '(A)',' (rad) no van Regemorter approx of line scattering'
        end if
     end if

     print '(A,I4)',       ' (rad) EOS.tab version      = ',Nxrevr(3)
     print '(A,A10,A,I4)', ' (rad)',linedat,' version    = ',Nxrevr(4)
     print '(A,I4)',       ' (rad) subs.dat version     = ',Nxrevr(5)

     print '(A)', hline
  end if

 END SUBROUTINE init
 
!-------------------------------------------------------------------------------
SUBROUTINE lookup_eos(nx,ny,nz,rho,ee,ivar,var,do_indx,do_exp,do_cubic)

  USE params, only: ub,mid,mpi_rank

  IMPLICIT none

  character(len=32) :: IdStr = 'lookup_eos'
  character(len=72) :: hline = repeat('-',72)

  integer,                   intent(in)  :: nx, ny, nz
  integer,                   intent(in)  :: ivar
  logical,                   intent(in)  :: do_indx, do_exp, do_cubic
  real, dimension(nx,ny,nz)              :: rho, ee
  real, dimension(nx,ny,nz), intent(out) :: var

  integer  :: nx_chk, ny_chk, nz_chk
  integer  :: ie, ir
  real     :: px, py, qx, qy, pxqx, pyqy, dfx1, dfx2, dfy1, dfy2
  real     :: rmin, rmax, emin, emax, ri, ei

  integer  :: ix, iy, iz
  logical  :: is_outside, is_alloc, do_alloc


  ! Check whether arrays need to be (re-)allocated
  is_alloc = allocated(iearr_e)
  if (is_alloc) then
     nx_chk = size(iearr_e,dim=1)
     ny_chk = size(iearr_e,dim=2)
     nz_chk = size(iearr_e,dim=3)
     do_alloc = nx_chk.ne.nx .and. ny_chk.ne.ny .and. nz_chk.ne.nz
  else
     do_alloc = .true.
  end if
  if (do_alloc) then
     if (is_alloc) then
        deallocate( iearr_e, irarr_e, pxarr_e, pyarr_e )
        deallocate( pbot0, dlnpdE_r, dlnpdlnr_E )
     end if
     allocate( iearr_e(nx,ny,nz), irarr_e(nx,ny,nz) ) 
     allocate( pxarr_e(nx,ny,nz), pyarr_e(nx,ny,nz) )
     allocate( pbot0(nx,nz), dlnpdE_r(nx,nz), dlnpdlnr_E(nx,nz) )
  end if

  ! Calculate indexes for interpolation
  if (do_indx) then

     rmin = minval(rho(:,:,:))
     rmax = maxval(rho(:,:,:))
     emin = minval( ee(:,:,:))
     emax = maxval( ee(:,:,:))

     ! Check if outside table
     is_outside = &
          rmin.lt.rholimMinE.or.rmax.gt.rholimMaxE.or. &
          emin.lt. eelimMinE.or.emax.gt. eelimMaxE
     if (is_outside) then 
        if (do_limit_tab) then ! force rho and e to witihin table limits
           do iz=1,nz
              do iy=1,ny
                 do ix=1,nx
                    rho(ix,iy,iz) = min(max(rho(ix,iy,iz),rholimMinE),rholimMaxE)
                    ee(ix,iy,iz)  = min( max(ee(ix,iy,iz), eelimMinE), eelimMaxE)
                 enddo
              enddo
           enddo
           if (verbose_tab) then
              print '(A)', hline
              print '(A)', trim(IdStr)//': Density and energy forced to be within limits'
              print '(A,X,I6)', trim(IdStr)//': MPI rank =', mpi%rank
              print '(A)', hline
           end if
        else
           print '(A)', hline
           print '(A)',trim(IdStr)//': Outside table'
           print '(A)',trim(IdStr)//': MPI rank = ', mpi%rank
           print '(5X,A,X,4(1PG14.6))', 'Table: lnrMin, lnrMax, lneMin, lneMax =', &
                lnrMinE, lnrMaxE, lneMinE, lneMaxE
           print '(5X,A,X,4(1PG14.6))', 'Sim:   lnrMin, lnrMax, lneMin, lneMax =', &
                log(rMin), log(rMax), log(eMin), log(eMax)
           print '(A)', hline
           call mpi%abort (trim(IdStr)': Outside table')
           stop
        end if
     end if

     ! Calculate indexes
     do iz=1,nz
        do iy=1,ny
           do ix=1,nx

              ri = (log(rho(ix,iy,iz))-lnrMinE)/dlnrE + 1.0
              ir = max0(1,min0(Nrtabe-1,int(ri)))
              pxarr_e_alt(ix,iy,iz) = ri-ir
              irarr_e_alt(ix,iy,iz) = ir

              ei = (log(ee(ix,iy,iz))-lneMinE)/dlneE + 1.0
              ie = max0(1,min0(Netabe-1,int(ei)))
              pyarr_e_alt(ix,iy,iz) = ei-ie
              iearr_e_alt(ix,iy,iz) = ie

           end do
        end do
     end do

  endif

  ! Cubic interpolation ?
  if (do_cubic) then  

     do iz=1,nz
        do iy=1,ny
           do ix=1,nx

              ir = irarr_e(ix,iy,iz)
              ie = iearr_e(ix,iy,iz)
              px = pxarr_e(ix,iy,iz)
              py = pyarr_e(ix,iy,iz)
              qx   = 1. - px
              qy   = 1. - py
              pxqx = px*qx
              pyqy = py*qy

              dfx1 = eostab(ie  ,ir+1,1,ivar)-eostab(ie  ,ir  ,1,ivar)
              dfx2 = eostab(ie+1,ir+1,1,ivar)-eostab(ie+1,ir  ,1,ivar)
              dfy1 = eostab(ie+1,ir  ,1,ivar)-eostab(ie  ,ir  ,1,ivar)
              dfy2 = eostab(ie+1,ir+1,1,ivar)-eostab(ie  ,ir+1,1,ivar)

              var(ix,iy,iz) = qy*(                                           &
                   &             qx*(      eostab(ie  ,ir  ,1,ivar)          &
                   &               + pxqx*(eostab(ie  ,ir  ,3,ivar)-dfx1)    &
                   &               + pyqy*(eostab(ie  ,ir  ,2,ivar)-dfy1))   &
                   &           + px*(      eostab(ie  ,ir+1,1,ivar)          &
                   &               - pxqx*(eostab(ie  ,ir+1,3,ivar)-dfx1)    &
                   &               + pyqy*(eostab(ie  ,ir+1,2,ivar)-dfy2))   &
                   &                                                      )  &
                   &        + py*(                                           &
                   &             px*(      eostab(ie+1,ir+1,1,ivar)          &
                   &               - pxqx*(eostab(ie+1,ir+1,3,ivar)-dfx2)    &
                   &               - pyqy*(eostab(ie+1,ir+1,2,ivar)-dfy2))   &
                   &           + qx*(      eostab(ie+1,ir  ,1,ivar)          &
                   &               + pxqx*(eostab(ie+1,ir  ,3,ivar)-dfx2)    &
                   &               - pyqy*(eostab(ie+1,ir  ,2,ivar)-dfy1))   &
                   &                                                      )

           end do
        end do
     end do

  else
     
     ! Linear interpolation
     do iz=1,nz
        do iy=1,ny
           do ix=1,nx

              ir = irarr_e(ix,iy,iz)
              ie = iearr_e(ix,iy,iz)
              px = pxarr_e(ix,iy,iz)
              py = pyarr_e(ix,iy,iz)
              qx   = 1. - px
              qy   = 1. - py
              var(ix,iy,iz) = &
                   qy*(qx*eostab(ie  ,ir,1,ivar)+px*eostab(ie  ,ir+1,1,ivar)) + &
                   py*(qx*eostab(ie+1,ir,1,ivar)+px*eostab(ie+1,ir+1,1,ivar))

           enddo
        end do
     end do
  end if


  ! Pressure: calculate auxiliary arrays needed by boundary routines
  if (ivar.eq.1) then
     do iz=1,nz
        do ix=1,nx

           pbot0(ix,iz) = var(ix,ub,iz)
           ir = irarr_e(ix,ub,iz)
           ie = iearr_e(ix,ub,iz)
           px = pxarr_e(ix,ub,iz)
           py = pyarr_e(ix,ub,iz)
           qx   = 1. - px
           qy   = 1. - py

           ! Look up dlnpdlne
           ! NOTE:derivatives stored in the tables are multiplied by dlne
           dlnpdE_r(ix,iz) =  &
                qy*(qx*eostab(ie  ,ir,2,ivar)+px*eostab(ie  ,ir+1,2,ivar)) + &
                py*(qx*eostab(ie+1,ir,2,ivar)+px*eostab(ie+1,ir+1,2,ivar))

           ! dlnpde = dlnpdlne / ee   (in simulation units!)
           dlnpdE_r(ix,iz) = dlnpdE_r(ix,iz) / ee(ix,ub,iz) / dlneE


           ! Look up dlnpdlnr
           ! NOTE:derivatives stored in the tables are multiplied by dlnr
           dlnpdlnr_E_alt(ix,iz) =  &
                qy*(qx*eostab(ie  ,ir,3,ivar)+px*eostab(ie  ,ir+1,3,ivar)) + &
                py*(qx*eostab(ie+1,ir,3,ivar)+px*eostab(ie+1,ir+1,3,ivar))

           dlnpdlnr_E(ix,iz) = dlnpdlnr_E(ix,iz) / dlnrE
        enddo
     end do

     ! Exponentiate
     if (do_exp) then
        pbot0(:,:) = exp(pbot0(:,:))
     end if

  end if

  ! Exponentiate
  if (do_exp) then
     do iz=1,nz
        var(:,:,iz) = exp(var(:,:,iz))
     end do
  end if

END SUBROUTINE lookup_eos


!-------------------------------------------------------------------------------
SUBROUTINE lookup_rad(nx,ny,nz,rho,tt,ivar,ilam,var,do_indx,do_exp)

  IMPLICIT none

  character(len=32) :: IdStr = 'lookup_rad'
  character(len=72) :: hline = repeat('-',72)

  integer,                   intent(in)  :: nx,ny,nz
  integer,                   intent(in)  :: ivar, ilam
  logical,                   intent(in)  :: do_indx, do_exp
  real, dimension(nx,ny,nz)              :: rho,tt
  real, dimension(nx,ny,nz), intent(out) :: var

  integer  :: nx_chk, ny_chk, nz_chk
  integer  :: ir, itt
  real     :: px, py, qx, qy
  real     :: rmin,rmax,tmin,tmax,ri,ti
  integer  :: ix,iy,iz
  logical  :: is_outside, is_alloc, do_alloc

  ! (Re-)allocate arrays ?
  is_alloc = allocated(irarr_r)
  if (is_alloc) then
     nx_chk = size(irarr_r,dim=1)
     ny_chk = size(irarr_r,dim=2)
     nz_chk = size(irarr_r,dim=3)
     do_alloc = nx_chk.ne.nx .and. ny_chk.ne.ny .and. nz_chk.ne.nz
  else
     do_alloc = .true.
  end if
  if (do_alloc) then
     if (is_alloc) then
        deallocate( irarr_r, itarr_r, pxarr_r, pyarr_r )
     end if
     allocate( irarr_r(nx,ny,nz), itarr_r(nx,ny,nz) ) 
     allocate( pxarr_r(nx,ny,nz), pyarr_r(nx,ny,nz) )
  end if

  ! Calculate indexes for interpolation?
  if (do_indx) then

     rmin = minval(rho(:,:,:))
     rmax = maxval(rho(:,:,:))
     tmin = minval( tt(:,:,:))
     tmax = maxval( tt(:,:,:))

     is_outside = &
          rmin.lt.rholimMinR.or.rmax.gt.rholimMaxR.or.&
          tmin.lt. ttlimMinR.or.tmax.gt. ttlimMaxR
     if (is_outside) then
        if (do_limit_tab) then ! force rho and tt to witihin table limits
           do iz=1,nz
              do iy=1,ny
                 do ix=1,nx
                    rho(ix,iy,iz) = min(max(rho(ix,iy,iz),rholimMinR),rholimMaxR)
                    tt(ix,iy,iz)  = min( max(tt(ix,iy,iz), ttlimMinR), ttlimMaxR)
                 end do
              end do
           end do
           if (verbose_tab) then
              print '(A)', hline
              print '(A)', trim(IdStr)//': Density and temperature forced to be within limits'
              print '(A,X,I6)', trim(IdStr)//': MPI rank =', mpi%rank
              print '(A)', hline
           end if
        else   ! print error message and exit
           print '(A)', hline
           print '(A)',trim(IdStr)//': Outside table'
           print '(A)',trim(IdStr)//': MPI rank = ', mpi%rank
           print '(5X,A,X,4(1PG14.6))', 'Table: lnrMin, lnrMax, lnTMin, lnTMax =', &
                lnrMinR, lnrMaxR, lnTMinR, lnTMaxR
           print '(5X,A,X,4(1PG14.6))', 'Sim:   lnrMin, lnrMax, lnTMin, lnTMax =', &
                log(rMin), log(rMax), log(tMin), log(tMax)
           print '(A)', hline
           call mpi%abort (trim(IdStr)': Outside table')
           stop
        end if
     end if

     ! Calculate indexes
     do iz=1,nz
        do iy=1,ny
           do ix=1,nx

              ri = (log(rho(ix,iy,iz))-lnrMinR)/dlnrR + 1.0
              ir = max0(1,min0(Nrtabr-1,int(ri)))
              pxarr_r(ix,iy,iz) = ri-ir
              irarr_r(ix,iy,iz) = ir

              ti = (log(tt(ix,iy,iz))-lntMinR)/dlntR + 1.0
              itt = max0(1,min0(Nttabr-1,int(ti)))
              pyarr_r(ix,iy,iz) = ti-itt
              itarr_r(ix,iy,iz) = itt

           end do
        end do
     end do

  end if

  ! Linear interpolation
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           ir = irarr_r(ix,iy,iz)
           itt= itarr_r(ix,iy,iz)
           px = pxarr_r(ix,iy,iz)
           py = pyarr_r(ix,iy,iz)
           qx   = 1. - px
           qy   = 1. - py
           var(ix,iy,iz) = &
                qy*(qx*radtab(itt  ,ir,ivar,ilam)+px*radtab(itt  ,ir+1,ivar,ilam)) + &
                py*(qx*radtab(itt+1,ir,ivar,ilam)+px*radtab(itt+1,ir+1,ivar,ilam))
        end do
     end do
  end do

  ! Exponentiate
  if (do_exp) then
     do iz=1,nz
        var(:,:,iz) = exp(var(:,:,iz))
     end do
  end if

END SUBROUTINE lookup_rad

!-------------------------------------------------------------------------------
SUBROUTINE pressure(nx, ny, nz, rho, ee, pp)

  IMPLICIT none

  character(len=32) :: IdStr = 'pressure'

  integer,                   intent(in)    :: nx,ny,nz
  real, dimension(nx,ny,nz)                :: rho,ee
  real, dimension(nx,ny,nz), intent(out)   :: pp

  logical :: do_indx, do_exp

  do_indx = .true.
  do_exp  = .true.

  call lookup_eos(nx, ny, nz, rho, ee, 1, pp, do_indx ,do_exp, do_cubic)

END SUBROUTINE pressure

!-------------------------------------------------------------------------------
SUBROUTINE temperature(nx, ny, nz, rho, ee, tt)

  IMPLICIT none

  character(len=32) :: IdStr = 'temperature'

  integer,                   intent(in)    :: nx,ny,nz
  real, dimension(nx,ny,nz)                :: rho,ee
  real, dimension(nx,ny,nz), intent(out)   :: tt

  logical :: do_indx, do_exp

  do_indx = .true.
  do_exp  = .true.

  call lookup_eos(nx, ny, nz, rho, ee, 3, tt, do_indx, do_exp, do_cubic)

END SUBROUTINE temperature

!-------------------------------------------------------------------------------
SUBROUTINE rosseland(nx, ny, nz, rho, ee, lnrkRoss)

! Note: subroutine returns ln Rosseland extinction coefficient

  IMPLICIT none

  character(len=32) :: IdStr = 'rosseland'

  integer,                   intent(in)    :: nx, ny, nz
  real, dimension(nx,ny,nz)                :: rho, ee
  real, dimension(nx,ny,nz), intent(out)   :: lnrkRoss

  logical :: do_indx, do_exp

  do_indx = .true.
  do_exp  = .false.

  call lookup_eos(nx, ny, nz, rho, ee, 2, lnrkRoss, do_indx ,do_exp, do_cubic)
  !call lookup_rad(nx, ny, nz, rho, tt, 1, nlam+1, lnrkRoss, do_indx, do_exp)

END SUBROUTINE rosseland

!-------------------------------------------------------------------------------
SUBROUTINE elnumdens(nx, ny, nz, rho, tt, lnnel)

! Note: subroutine returns ln electron number density (in simulation units!)

  IMPLICIT none

  character(len=32) :: IdStr = 'elnumdens'

  integer,                   intent(in)    :: nx, ny, nz
  real, dimension(nx,ny,nz)                :: rho, tt
  real, dimension(nx,ny,nz), intent(out)   :: lnnel

  logical :: do_indx, do_exp
  logical :: has_elnumdens

  do_indx = .true.
  do_exp  = .false.

  has_elnumdens = Ntbvar.eq.6

  if (has_elnumdens) then
     call lookup_eos(nx, ny, nz, rho, ee, 4, lnnel, do_indx ,do_exp, do_cubic)
  else
     call mpi%abort (trim(IdStr)': EOS table does not contain electron number density')
     stop
  end if
  
END SUBROUTINE elnumdens

END MODULE eos_mod
!-------------------------------------------------------------------------------


