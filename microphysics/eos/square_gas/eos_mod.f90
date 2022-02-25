!===============================================================================
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
!=======================================================================
MODULE eos_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE timer_mod
  USE kinds_mod
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
    real:: lne_min, lne_max, lne_stp, ee_min
    real:: mu_l = 0.62    
    real(8):: gamma = 5.0e0/3.0e0
    real, dimension(:,:,:), pointer :: eostable,epstable
    real, dimension(:,:,:), pointer :: opatable,temtable
    real, dimension(:,:,:), pointer :: abl
    real, dimension(:,:,:), pointer :: bplp,epstab,opatab
    real, dimension(:,:,:,:), pointer :: ne_epstab,ne_opatab,ne_temtab
    logical :: initialized = .false.
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
    procedure:: ddeos
    procedure:: lookup
    procedure:: inverse_lookup
    procedure, private:: lookup_table
  end type
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> Initialize the EOS tables
!===============================================================================
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
  real,save:: gamma=5.0e0/3.0e0
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
  real:: mudens, muntot, lneir
  real, dimension(:),pointer:: lnei
  !.....................................................................
  namelist /eos_params/ table_loc, tabparam, gamma
  namelist /tabpar/     EOSTableFile, RhoEiRadTableFile, NeTgRadTableFile, &
                        nRhoBin, RhoMin, RhoMax, nEiBin, EiMin, EiMax, nRadBins
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call io%lock%set ('eos_t%init')
  if (first_time) then
    first_time = .false.
    do_eos = .true.
    call scaling%init
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
    self%gamma = gamma
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
    self%ee_min  = eimin/scaling%ee
    self%lne_max = log(eimax/scaling%ee)
    self%lne_stp = (log(eimax-eimin)-log(scaling%ee))/real(self%ne-1)
    allocate(lnEi(nEiBin))
    lnEir=log(EiMax)-log(EiMin)
    do i=1,nEiBin
      lnEi(i)=log(EiMin)+real(i-1)/real(nEiBin-1)*lnEir
    end do
    lnei=lnei-log(scaling%ee)
    self%lne_stp=lnEi(2)-lnEi(1)
    self%lne_min=lnei(1)
    self%lne_max=lnei(neibin)
    deallocate(lnei)
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
    if (io%master) print *, 'reading tables..'
    recl = kind(1.0)*nrhobin*neibin
    ! read the EOS tables
    file = trim(table_loc)//'/'//trim(self%EOSTableFile)
    open  (io_unit%data, file=file, access='direct', form='unformatted', status='old', &
                                    recl=4*recl)
    read  (io_unit%data,rec=1) self%eostable
    close (io_unit%data)
    if (io%master) print *, 'EOS table read',trim(self%EOSTableFile)
    ! read radiation tables
    file = trim(table_loc)//'/'//trim(self%RhoEiRadTableFile)
    open  (io_unit%data, file=file, access='direct', form='unformatted', status='old', &
                                    recl=self%n_lambda*recl)
    read  (io_unit%data, rec=1) self%epstable
    read  (io_unit%data, rec=2) self%temtable
    read  (io_unit%data, rec=3) self%opatable
    close (io_unit%data)
    if (io%master) print *, 'radiation table',trim(self%RhoEiRadTableFile), 'read'
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
  self%initialized=.true.
  call io%lock%unset ('eos_t%_init')
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, lnd, d, lne, lnee, e, ee, ss, s, pg, tt, ne, rk, ross, &
                   scat, src, dlnpdE_r, dlnpdlnr_E, gamma, id, e_below, label)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, d, lne, lnee, e, ee
  real, dimension(:,:,:), intent(in), pointer, optional :: ss, s
  real, dimension(:,:,:),                      optional :: pg, tt, ne, ross
  real, dimension(:,:,:),                      optional :: dlnpdE_r, dlnpdlnr_E
  real, dimension(:,:,:,:),                    optional :: src, rk, scat
  real, optional:: gamma
  integer, optional:: id
  logical, optional:: e_below
  character(len=*), optional::label
  !.............................................................................
  real, dimension(:,:,:), pointer :: lnd_L, ee_L
  real(8):: g1
  integer:: dim(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  !-----------------------------------------------------------------------------
  if (.not. self%initialized) call self%init
    if (io%verbose>2) then
    print *, 'present(lnd) =', present(lnd)
    print *, 'present(d) =', present(d)
    print *, 'present(ee)  =', present(ee)
    print *, 'present(e)  =', present(e)
    print *, 'present(lne)  =', present(lne)
    print *, 'present(lnee)  =', present(lnee)
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
  if (present(gamma)) then
    self%gamma = gamma
  end if
  g1 = merge(1d0, self%gamma-1d0, self%gamma==1d0)
  !-----------------------------------------------------------------------------
  ! Establish d_L as the density
  !-----------------------------------------------------------------------------
  nullify (lnd_L, ee_L)
  if (present(d)) then
    dim = shape(d)
    allocate (lnd_L(dim(1),dim(2),dim(3)))
    lnd_L = log(d)
  else if (present(lnd)) then
    dim = shape(lnd)
    lnd_L => lnd
  else
    call io%abort ('eos_t: missing 1st argument')
  end if
  !-----------------------------------------------------------------------------
  ! Establish ee_L as the internal energy per unit volume
  !-----------------------------------------------------------------------------
  if (present(ee)) then
    !--- energy per unit mass given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    ee_L = log(ee)
  else if (present(lnee)) then
    !--- log energy per unit mass given ---
    ee_L => lnee
  else if (present(lne)) then
    !--- log energy per unit volume given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    ee_L = lne - lnd_L
  else if (present(e)) then
    !--- energy per unit volume given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    if (.not. present(d)) then
      ee_L = log(e)-lnd_L
    else
      ee_L = log(e/d)
    end if
  else if (present(ss).or.present(s)) then
    call io%abort('square_gas EOS table not set up for entropy')
  else
    call io%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  ! Do the table lookup and deallocate locals
  !-----------------------------------------------------------------------------  
  call lookup_table (self, dim, lnd_L, ee_L, pg=pg, tt=tt, ne=ne, rk=rk, scat=scat, &
                     src=src, ross=ross, dlnpdE_r=dlnpdE_r, dlnpdlnr_E=dlnpdlnr_E, &
                     gamma=gamma, id=id, e_below=e_below, label=label)
  if (.not.present(lnee)) deallocate( ee_L)
  if (present(d)) deallocate(lnd_L)
  !-----------------------------------------------------------------------------
  call trace%end(itimer)  


END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, lnd, lnee, pg, tt, ne, rk, src, ross, scat,&
                         dlnpdE_r, dlnpdlnr_E, gamma, id, e_below, label)
  class (eos_t):: self
  integer:: dim(3)
  integer, optional:: id
  real, dimension(:,:,:)  , optional :: lnd, lnee, ne, ross
  real, dimension(:,:,:,:), optional :: src,rk,scat
  real, dimension(:,:,:)  , optional :: pg,tt
  real, dimension(:,:,:)  , optional:: dlnpdE_r, dlnpdlnr_E

  real, optional:: gamma
  logical, optional:: e_below
  character(len=*),optional::label
  ! --------------------------------------------------------------------
  integer,         dimension(dim(1)) :: il, jl
  real(kind=4),    dimension(dim(1)) :: f1, f2, f3, f4
  real        ,    dimension(dim(1)) :: fr, fe
  real, dimension(:), pointer :: pg_loc
  real    :: dmin,dmax,dd, emin,emax,de, m1t, m1u, t, u
  integer :: ix,iy,iz, j
  integer :: mx,my,mz
  integer, save:: itimer=0
  ! ....................................................................
  call trace%begin ('eos_t%lookup_table', itimer=itimer)
  mx=dim(1)
  my=dim(2)
  mz=dim(3)
  dmin=self%lnd_min
  dmax=self%lnd_max
  dd  =self%lnd_stp
  !
  emin=self%lne_min
  emax=self%lne_max
  de  =self%lne_stp
  ! ----------------------------------------------------------------------------
  ! Allocate scratches, if needed
  ! ----------------------------------------------------------------------------
  if (present(dlnpdE_r).or.present(dlnpdlnr_E)) then
    allocate(pg_loc(mx))
  end if
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
      fr(ix) = (lnd(ix,iy,iz)-dmin)/dd+1d0
      jl(ix) = aint(fr(ix))
      u = mod(fr(ix),1d0)
      m1u = 1d0-u
      ! log(ei)
      fe(ix) = (lnee(ix,iy,iz)-emin)/de+1d0
      il(ix) = aint(fe(ix))
      t = mod(fe(ix),1d0)
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
    if (any(jl>self%nd-1).or.any(jl<1)) then
      print *, 'minval, maxval lnd', minval(lnd), maxval(lnd)
      print *, 'table ranges', dmin, dmax
      if (present(label)) print *, 'table was called from ', trim(label)
      call io%abort("EOS outside of density range")
    end if
    if (any(il>self%ne-1).or.any(il<1)) then
      print *, il, 'ne total', self%ne
      print *, 'minval, maxval lnee', minval(lnee), exp(minval(lnee)), maxval(lnee), exp(maxval(lnee))
      print *, 'table ranges', emin, exp(emin),emax, exp(emax)
      print *, 'minloc ee', minloc(lnee)
      print *, 'maxloc ee', maxloc(lnee)
      print *, 'temperature would be approx., loc, table', exp(minval(lnee))*0.6666667, exp(emin)*0.6666667, exp(emax)*0.6666667
      if (present(label)) print *, 'table was called from ', trim(label)
      call io%abort("EOS outside of internal energy range")
    end if
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
    ! Derivatives
    if (present(dlnpdE_r).or.present(dlnpdlnr_E)) then
      ! gas pressure
      if (present(pg)) then
        pg_loc = pg(:,iy,iz)
      else
        do ix=1,mx
          pg_loc(ix) = exp( &
              f1(ix) * self%eostable(il(ix)  ,jl(ix)  ,1) &
            + f2(ix) * self%eostable(il(ix)+1,jl(ix)  ,1) &
            + f3(ix) * self%eostable(il(ix)+1,jl(ix)+1,1) &
            + f4(ix) * self%eostable(il(ix)  ,jl(ix)+1,1))
        end do
      end if
      do ix =1,mx
        u = fr(ix)-jl(ix)
        m1u = 1d0-u
        t = fr(ix)-il(ix)
        m1t = 1d0-t
        dlnpdE_r(ix,iy,iz)   = (m1u * self%eostable(il(ix)+1,jl(ix)  ,1) &
                             -  m1u * self%eostable(il(ix)  ,jl(ix)  ,1) &
                             +    u * self%eostable(il(ix)+1,jl(ix)+1,1) &
                             -    u * self%eostable(il(ix)  ,jl(ix)+1,1) ) &
                             /   de * pg_loc(ix)/exp(lnee(ix,iy,iz))

        dlnpdlnr_e(ix,iy,iz) = ( m1t * self%eostable(il(ix)  ,jl(ix)+1,1) &
                             -   m1t * self%eostable(il(ix)  ,jl(ix)  ,1) &
                             +     t * self%eostable(il(ix)+1,jl(ix)+1,1) &
                             -     t * self%eostable(il(ix)+1,jl(ix)  ,1) ) &
                             / dd
      end do
    end if
  end do
  end do
  !
  ! ----------------------------------------------------------------------------
  ! Cleanup
  ! ----------------------------------------------------------------------------
  if (associated(pg_loc)) deallocate(pg_loc)

  !
  call trace%end (itimer)
END SUBROUTINE lookup_table

!===============================================================================
!> Interface to call with ideal equation of state fall-back
!> FIXME: Check if this is actually used, and replace with generic call
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, v, name, ss, ee, label)
  class (eos_t):: self
  real, dimension(:,:,:), pointer :: d, lnd, v
  optional                        :: d, lnd, v
  real, dimension(:,:,:), pointer :: ss, ee
  optional                        :: ss, ee
  character(*), optional:: label, name
  !.............................................................................
  real, dimension(:,:,:), pointer :: x_L, y_L, lnd_l, ee_l
  integer:: dim(3)
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup', itimer=itimer)
  call io%abort ('eos_t%inverse_lookup not ready yet')
  call trace%end(itimer)
END SUBROUTINE inverse_lookup


!===============================================================================
!> Thermodynamic derivatives of the gas pressure
! ee - internal energy
! d  - density
! pg - gas pressure
! dpdecr - dpde at constant density
! dpdrce - dpdr at constant energy
! dpdrcs - dpdr at constant entropy
!===============================================================================
SUBROUTINE ddeos(self,ee,d,dlnpde_r,dlnpdlnr_E,dpdrcs)
  class(eos_t):: self
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee, d
  real(kind=KindScalarVar), dimension(:,:), pointer:: dlnpde_r, dlnpdlnr_E
  real(kind=KindScalarVar), dimension(:,:), pointer, optional ::dpdrcs
  ! ----------------------------------------------------------------------------
  real    :: dmin,dmax,dd, emin,emax,de,fr, m1t, m1u, t, u, pg
  real    :: f1, f2, f3, f4
  integer :: i,j, il, jl
  integer:: n(2)
  ! ............................................................................
  n = shape(d)

  do j=1,n(2)
      do i=1,n(1)
        ! ----------------------------------------------------------------------
        ! fractional indices
        ! ----------------------------------------------------------------------
        ! log(rho)
        fr = (log(d(i,j))-self%lnd_min)/self%lnd_stp+1d0
        jl = int(fr)
        u = fr-jl
        m1u = 1d0-u
        ! log(ei)
        fr = (log(ee(i,j))-self%lne_min)/self%lne_stp+1d0
        il = int(fr)
        t = fr-il
        m1t = 1d0-t
        !
        f1= m1t * m1u
        f2=   t * m1u
        f3=   t *   u
        f4= m1t *   u
        ! ------------------------------------------------------------------
        ! Check ranges
        ! ------------------------------------------------------------------
        if ((jl+1>self%nd).or.(jl<1)) then
          print *, jl+1, self%nd
          print *, log(d(i,j)), self%lnd_min, self%lnd_max
          call io%abort("ddEOS outside of density range")
        end if
        if ((il+1>self%ne).or.(il<1)) then
          print *, il+1, self%ne
          print *, ee(i,j), self%lnd_min, self%lnd_max
          call io%abort("ddEOS outside of internal energy range")
        end if

        if ((il+1.le.self%ne).and.(jl+1.le.self%nd)) then
          ! linear interpolation
          ! thermodynamic derivatives of the gas pressure
          pg = exp( &
               f1 * self%eostable(il  ,jl  ,1) &
             + f2 * self%eostable(il+1,jl  ,1) &
             + f3 * self%eostable(il+1,jl+1,1) &
             + f4 * self%eostable(il  ,jl+1,1))
          dlnpdE_r(i,j) = ( m1u * self%eostable(il+1,jl  ,1) &
                        -   m1u * self%eostable(il  ,jl  ,1) &
                        +     u * self%eostable(il+1,jl+1,1) &
                        -     u * self%eostable(il  ,jl+1,1) ) &
                        / self%lne_stp * pg/ee(i,j)

          dlnpdlnr_e(i,j) = ( m1t * self%eostable(il  ,jl+1,1) &
                          -   m1t * self%eostable(il  ,jl  ,1) &
                          +     t * self%eostable(il+1,jl+1,1) &
                          -     t * self%eostable(il+1,jl  ,1) ) &
                          / self%lnd_stp
        endif
      end do ! ~ end over i
    end do ! ~ end over j


END SUBROUTINE ddeos

END MODULE eos_mod
