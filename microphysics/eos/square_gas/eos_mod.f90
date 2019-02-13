!> =============================================================================
!> Equation of state module. Imported from bifrost square_gas_mpi.f90
!> ---------------------------------------------------------------------
!> gamma: ratio of specific heats as in cs=sqrt(gamma*p/r)
!> p :    gas pressure [internal units]
!> tg:    gas temperature [K]
!> cs:    sound speed [internal units]
!> ee:    internal energy per mass unit [internal units]
!> ne:    electron density [internal units]
!> eostable: (:,:,:,4):    table with in this order
!>                         log Pg, tg, log ne, log KappaRos.
!>                         Units: [dyne/cm], [K], [cm-3],c[m2/g]
!> epstable: (:,:,:): table with log scattering probability per
!>                         radiation bin [dimensionless]
!> opatable: (:,:,:): table with log opacity per radiation bin [cm-1]
!> temtable: (:,:,:): table with log thermal emission per radiation bin
!>                         [erg / (s ster cm2)]
!> Lookup_tablse subroutine variables (all in code units):
!> dim :  dimensions of the array
!> e   :  internal energy
!> ee  :  e/d
!> d   :  density
!> lnd :  log(density)
!> lne :  log(e)
!> lnee:  log(ee)
!> pg  :  gas pressure
!> tt  :  gas temperature
!> ne  :  electron density
!> rk  :  opacity
!> src :  Thermal emission
!> scat:  scattering probability
!> ross:  Rosseland mean opacity
!> gamma: adiabatic index
!> ---------------------- Changelog ------------------------------------
!> 26/10/2018 Reader works. No NLTE for now.  -- AP
!> ---------------------------------------------------------------------
!> NOTES: for now, as tabparam.in is a separate file and not a header in
!> the tables, remember to add namelist flags (& and /) to the beginning
!> and the end of the file.
!> =====================================================================
MODULE eos_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE timer_mod
  implicit none
  private
  logical, save:: do_eos
  logical, save:: first_time=.true.
  type, public:: eos_t
    integer:: nd, ne
    integer:: nv      = 2
    integer:: i_pg    = 1
    integer:: i_tt    = 2
    integer:: i_kappa = 3
    integer:: i_src   = 4
    integer:: n_lambda = 4
    integer:: w_lambda = 1
    character(len=70) :: EOSTableFile = 'eostable.dat'
    character(len=70) :: RhoEiRadTableFile = 'rhoei_radtab.dat'
    character(len=70) ::NeTgRadTableFile = 'netg_radtab.dat'
    real:: lnd_min, lnd_max, lnd_stp
    real:: lne_min, lne_max, lne_stp
    real:: mu_l = 0.62    
    real:: gamma = 5.0e0/3.0e0
    real, dimension(:,:,:), pointer :: eostable,epstable
    real, dimension(:,:,:), pointer :: opatable,temtable
    real, dimension(:,:,:), pointer :: abl
    real, dimension(:,:,:), pointer :: bplp,epstab,opatab
    real, dimension(:,:,:,:), pointer :: ne_epstab,ne_opatab,ne_temtab
    logical :: do_rt     = .true.
    ! do_helium: swithces on non-equilibrium helium ionization if running helium exp.
    logical :: do_helium = .false.
    ! do_hion: if t then switches HION on, this requires the presence 
    !          of hion in the EXPERIMENT extras. It defaults to zero.
    logical :: do_hion   = .false.
    ! noneq: if T then use Ne-Tg radiation table in chromosphere 
    !        if F then use LTE based rho-ei radiation table in chromosphere
    logical :: noneq     = .false.                                      
    ! helium_exp: When running helium-experiment this variable is 
    ! derived from do_hion and do_helium. Controls incrad_euv that
    ! includes thin radiative losses pluss heating and cooling in the
    ! Lyman continuum and the two helium continuas.
    logical :: helium_exp=.false.
    ! rmfp: Rosseland mean free path [cm]. Used in determination of thin
    !       and thick weights of the ne-tg radiation tables
    real :: rmfp=1.5d7
    !allows to expand the EOS table using ideal gas
    logical :: expieos=.false.
    real :: maxlnelim=32.3
    real :: minlnrlim=-34.0
    ! arrays that contain the elements included in the OES table and
    ! their abundance
    integer :: nelmax=26
    integer :: nel
    character(len=2),dimension(:), pointer :: el
    real            ,dimension(:), pointer :: abnd,awght
    ! number of H atoms per gram
    real:: grph = 0.0
  contains
    procedure:: init
    procedure:: pressure
    procedure:: lookup
    procedure:: lookup_table
  end type
  type (eos_t), public:: eos
CONTAINS

!> =============================================================================
!> Initialize the EOS tables
!> =============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  integer, save:: nRhoBin = 300
  integer, save:: nEiBin  = 1000
  integer, save:: nRadBin = 4
  integer, save:: nRadBins = 4
  integer, save:: nel = 0
  real, save:: RhoMin = 1e-18
  real, save:: RhoMax = 7e-5
  real, save:: EiMin  = 2.00466172E+12
  real, save:: EiMax  = 2.42065997E+15
  real, save:: mu_l = 0.62
  character(len=2), dimension(30), save:: cel = ''
  real, dimension(3), save:: abund   = 0.0
  real, dimension(3), save:: aweight = 0.0
  character(len=70), save :: tabparam = 'tabparam.in'
  character(len=70), save :: EOSTableFile = 'eostable.dat'
  character(len=70), save :: RhoEiRadTableFile = 'rhoei_radtab.dat'
  character(len=70), save :: NeTgRadTableFile = 'netg_radtab.dat'
  character(len=70), save :: table_loc = '../../data/eos/square_gas'
  character(len=70)       :: file
  !.....................................................................
  integer:: i, recl
  real:: mudens, muntot
  !.....................................................................
  namelist /eos_params/ table_loc, tabparam
  namelist /tabpar/     EOSTableFile, RhoEiRadTableFile, NeTgRadTableFile, &
                        nRhoBin, RhoMin, RhoMax, nEiBin, EiMin, EiMax, nRadBins
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  
  !$omp critical (eos_cr)
  if (first_time) then
    first_time = .false.
    do_eos = .true.
    rewind (io%input); read (io%input,eos_params)
    if (io%master) write (*,eos_params)
    ! ------------------------------------------------------------------
    ! Set up table indexes. These are given in tabparam.in
    ! ------------------------------------------------------------------
    file = trim(table_loc)//'/'//trim(tabparam)
    open (io_unit%nml1, file=file, form='formatted', status='old')
    read (io_unit%nml1, tabpar)
    if (io%master) write (*,tabpar)
    close(io_unit%nml1)
    !
    self%nd = nrhobin
    self%ne = neibin
    self%n_lambda = nradbins
    self%nel = nel
    self%EOSTableFile      = EOSTableFile
    self%RhoEiRadTableFile = RhoEiRadTableFile
    self%NeTgRadTableFile  = NeTgRadTableFile
    ! set log mass density axis 
    self%lnd_min = log(rhomin/scaling%d)
    self%lnd_max = log(rhomax/scaling%d)
    self%lnd_stp = (self%lnd_max-self%lnd_min)/real(self%nd-1)
    ! set log internal energy axis
    self%lne_min = log(eimin/scaling%ee)
    self%lne_max = log(eimax/scaling%ee)
    self%lne_stp = (self%lne_max-self%lne_min)/real(self%ne-1)
    
    ! ------------------------------------------------------------------
    ! Allocate arrays
    ! ------------------------------------------------------------------
    allocate(self%eostable(self%ne,self%nd,4))                          ! eos table
    allocate(self%epstable(self%ne,self%nd,self%n_lambda))              ! scattering table
    allocate(self%opatable(self%ne,self%nd,self%n_lambda))              ! opacity table
    allocate(self%temtable(self%ne,self%nd,self%n_lambda))              ! thermal emission table
    !-------------------------------------------------------------------
    ! Read the table files
    !-------------------------------------------------------------------
    print *, 'reading tables..'
    recl = kind(1.0)*nrhobin*neibin
         
    ! read the EOS tables
    file = trim(table_loc)//'/'//trim(self%EOSTableFile)
    open  (io_unit%data, file=file, access='direct', form='unformatted', status='old', &
                                    recl=4*recl)
    !do i=1,4            
    !print *, i                    
      read  (io_unit%data,rec=1) self%eostable
    !end do
    close (io_unit%data)
    print *, 'EOS table read',trim(self%EOSTableFile)
    ! read radiation tables
    file = trim(table_loc)//'/'//trim(self%RhoEiRadTableFile)
    open  (io_unit%data, file=file, access='direct', form='unformatted', status='old', &
                                    recl=self%n_lambda*recl)
    read  (io_unit%data, rec=1) self%epstable
    read  (io_unit%data, rec=2) self%temtable
    read  (io_unit%data, rec=3) self%opatable
    close (io_unit%data)
    print *, 'radiation table',trim(self%RhoEiRadTableFile), 'read'
    ! ------------------------------------------------------------------
    ! Convert to code units
    ! ------------------------------------------------------------------
    self%eostable(:,:,1)=self%eostable(:,:,1)-log(scaling%p)            ! gas pressure
    self%eostable(:,:,2)=self%eostable(:,:,2)/scaling%temp              ! temperature
    self%eostable(:,:,3)=self%eostable(:,:,3)-3.0*log(scaling%l)        ! electron density
    self%eostable(:,:,4)=self%eostable(:,:,4)-log(scaling%kr)           ! rosseland opacity
    self%temtable(:,:,:)=self%temtable(:,:,:)-log(scaling%e/scaling%t*scaling%l) ! box thermal emission erg/(s ster cm2)]
    self%opatable(:,:,:)=self%opatable(:,:,:)+log(scaling%l)
    ! ------------------------------------------------------------------
    ! non-equilibrium stuff, H-ion, Heliu,, etc. should go here
    ! ------------------------------------------------------------------
    if (self%noneq.and.self%do_hion) &
      call io%abort("eos_mod::init: noneq and hion not set up in this eos yet")
    !
  end if
  !$omp end critical (eos_cr)
  !$omp barrier

  call trace%end()
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION pressure (self, lnd, lne) RESULT (pg)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: lnd, lne
  real, dimension(size(lnd,1),size(lnd,2),size(lnd,3)):: pg
  !-----------------------------------------------------------------------------
  call self%lookup ([size(lnd,1),size(lnd,2),size(lnd,3)], lnd, lne, pg=pg)
END FUNCTION pressure

!> =============================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =============================================================================
SUBROUTINE lookup (self, dim, lnd, lne, d, ee, pg, tt, ne, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, lne
  real, dimension(:,:,:),   optional :: d, ee, pg, tt, ne
  real, dimension(:,:,:,:), optional :: src, rk
  real, optional:: gamma
  !.............................................................................
  real, dimension(:,:,:), pointer:: lnd_loc, lne_loc
  real, dimension(dim(1)):: pr, pe, qr, qe
  real, dimension(dim(1),4):: g
  integer, dimension(dim(1)):: ir, ie
  integer:: ix, iy, iz, iv
  integer:: md, me, mv
  real:: clnd, clne
  logical:: skip
  call io%abort ("eos_mod::lookup is not set up for square_gas eos")
  if (present(d))  deallocate (lnd_loc)
  if (present(ee)) deallocate (lne_loc)
END SUBROUTINE lookup

!> =====================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =====================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne,lnee, pg, tt, ne, rk, src, scat, ross,gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:)  , pointer, optional :: lnd,d,lne,lnee,e,ee,tt,ne,ross
  real, dimension(:,:,:,:)         , optional :: src,rk,scat
  real, dimension(:,:,:)           , optional :: pg
  real, optional:: gamma
  ! --------------------------------------------------------------------
  real, dimension(:,:,:), pointer    :: lnd_loc, ee_loc
  integer,         dimension(dim(1)) :: il, jl
  real(kind=4),    dimension(dim(1)) :: f1, f2, f3, f4
  real    :: dmin,dmax,dd, emin,emax,de,fr, m1t, m1u, t, u
  integer :: ix,iy,iz, j
  integer :: mx,my,mz
  integer, save:: itimer=0
  ! ....................................................................
  call trace%begin ('eos_t%lookup_table', itimer=itimer)
  if (io%verbose>2) then
    print *, 'present(lnd) =', present(lnd)
    print *, 'present(ee)  =', present(ee)
    print *, 'present(pg)  =', present(pg)
    print *, 'present(src) =', present(src)
    print *, 'present(rk)  =', present(rk)
    print *, 'present(scat)  =', present(scat)
    print *, 'present(tt)  =', present(tt)
    print *, 'present(ne)  =', present(ne)
    print *, 'present(ross)  =', present(ross)
    print *, 'present(gamma)  =', present(gamma)
    print *, 'dim:', dim
  end if
  mx=dim(1)
  my=dim(2)
  mz=dim(3)
  ! --------------------------------------------------------------------
  ! Allocate temporary arrays, if needed
  ! --------------------------------------------------------------------
  ! For log(density)
  if (present(lnd)) then
    lnd_loc => lnd
  else if (present(d)) then
    allocate (lnd_loc(dim(1),dim(2),dim(3)))
    lnd_loc = log(d)
  else
    call io%abort ('eos_t: missing density')
  end if
  ! For log(internal energy)
  if (present(lnee)) then
    ee_loc => lnee
  else if (present(e)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = log(e/d)
  else if (present(lne)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = lne-lnd_loc
  else if (present(ee)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = log(ee)
  else
    call io%abort ('eos_t: missing internal energy')
  end if
  !
  dmin=self%lnd_min
  dmax=self%lnd_max
  dd  =self%lnd_stp
  !
  emin=self%lne_min
  emax=self%lne_max
  de  =self%lne_stp
  ! --------------------------------------------------------------------
  ! Loop over points
  ! --------------------------------------------------------------------
  do iz=1,mz
  do iy=1,my
    ! ------------------------------------------------------------------
    ! fractional indices
    ! ------------------------------------------------------------------
    do ix=1,mx
      ! log(rho)
      fr = (lnd_loc(ix,iy,iz)-dmin)/dd+1d0
      jl(ix) = aint(fr)
      u = mod(fr,1d0)
      m1u = 1d0-u
      ! log(ei)
      fr = (ee_loc(ix,iy,iz)-emin)/de+1d0
      il(ix) = aint(fr)
      t = mod(fr,1d0)
      m1t = 1d0-t
      !
      f1(ix)= m1t * m1u
      f2(ix)=   t * m1u
      f3(ix)=   t *   u
      f4(ix)= m1t *   u
    end do
    ! ------------------------------------------------------------------
    ! Check ranges
    ! ------------------------------------------------------------------
    if (any(jl>self%nd-1).or.any(jl<1)) &
      call io%abort("EOS outside of density range")
    if (any(il>self%ne-1).or.any(il<1)) &
      call io%abort("EOS outside of internal energy range")
    ! ------------------------------------------------------------------
    !  Loop over table entries, nvar = 4.  Would this be more efficent
    !  if split into one gather loop and one calculation loop?
    ! ------------------------------------------------------------------
    !
    ! gas pressure
    if (present(pg)) then
      do ix=1,mx
        pg(ix,iy,iz) = exp( &
          f1(ix) * self%eostable(il(ix)  ,jl(ix)  ,1) &
        + f2(ix) * self%eostable(il(ix)+1,jl(ix)  ,1) &
        + f3(ix) * self%eostable(il(ix)+1,jl(ix)+1,1) &
        + f4(ix) * self%eostable(il(ix)  ,jl(ix)+1,1))
      end do
    end if
    !-------------------------------------------------------------------
    ! temperature
    if (present(tt)) then
      do ix=1,mx
        tt(ix,iy,iz) = ( &
          f1(ix) * self%eostable(il(ix)  ,jl(ix)  ,2) &
        + f2(ix) * self%eostable(il(ix)+1,jl(ix)  ,2) &
        + f3(ix) * self%eostable(il(ix)+1,jl(ix)+1,2) &
        + f4(ix) * self%eostable(il(ix)  ,jl(ix)+1,2))
      end do
    end if
    !-------------------------------------------------------------------
    ! electron density
    if (present(ne)) then
      do ix=1,mx
        ne(ix,iy,iz) = exp( &
          f1(ix) * self%eostable(il(ix)  ,jl(ix)  ,3) &
        + f2(ix) * self%eostable(il(ix)+1,jl(ix)  ,3) &
        + f3(ix) * self%eostable(il(ix)+1,jl(ix)+1,3) &
        + f4(ix) * self%eostable(il(ix)  ,jl(ix)+1,3))
      end do
    end if
    !-------------------------------------------------------------------
    ! Rosseland mean opacity
    if (present(ross)) then
      do ix=1,mx
        ross(ix,iy,iz) = exp( &
          f1(ix) * self%eostable(il(ix)  ,jl(ix)  ,4) &
        + f2(ix) * self%eostable(il(ix)+1,jl(ix)  ,4) &
        + f3(ix) * self%eostable(il(ix)+1,jl(ix)+1,4) &
        + f4(ix) * self%eostable(il(ix)  ,jl(ix)+1,4))
      end do
    end if
    !-------------------------------------------------------------------
    ! Scattering probability
    if (present(scat)) then
      do j=1, self%n_lambda
        do ix=1,mx
          scat(ix,iy,iz,j) = ( &
            f1(ix) * self%epstable(il(ix)  ,jl(ix)  ,j) &
          + f2(ix) * self%epstable(il(ix)+1,jl(ix)  ,j) &
          + f3(ix) * self%epstable(il(ix)+1,jl(ix)+1,j) &
          + f4(ix) * self%epstable(il(ix)  ,jl(ix)+1,j))
        end do
      end do
    end if
    !-------------------------------------------------------------------
    ! Thermal emission
    if (present(src)) then
      do j=1, self%n_lambda
        do ix=1,mx
          src(ix,iy,iz,j) = exp( &
            f1(ix) * self%temtable(il(ix)  ,jl(ix)  ,j) &
          + f2(ix) * self%temtable(il(ix)+1,jl(ix)  ,j) &
          + f3(ix) * self%temtable(il(ix)+1,jl(ix)+1,j) &
          + f4(ix) * self%temtable(il(ix)  ,jl(ix)+1,j))
        end do
      end do
    end if
    !-------------------------------------------------------------------
    ! Opacity
    if (present(rk)) then
      do j=1, self%n_lambda
        do ix=1,mx
          rk(ix,iy,iz,j) = exp( &
            f1(ix) * self%opatable(il(ix)  ,jl(ix)  ,j) &
          + f2(ix) * self%opatable(il(ix)+1,jl(ix)  ,j) &
          + f3(ix) * self%opatable(il(ix)+1,jl(ix)+1,j) &
          + f4(ix) * self%opatable(il(ix)  ,jl(ix)+1,j))
        end do
      end do
    end if
    !-------------------------------------------------------------------
  end do
  end do
  !
  if (present(d)) then
    deallocate (lnd_loc)
    nullify (lnd_loc)
  end if
  if (present(e).or.present(lne).or.present(ee)) then
    deallocate (ee_loc)
    nullify (ee_loc)
  end if
  call trace%end (itimer)
END SUBROUTINE lookup_table

END MODULE eos_mod
