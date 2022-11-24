MODULE coll_mod
    USE io_mod
    USE os_mod
    USE m3io_mod
    USE mmath,only:ei1,hermite_interpolation,hermite_derivative
    USE mconst
    USE prec_mod
    implicit none


   ! *******************************************************************

   type coll_type
      character(len=16) :: key ! transition type 
      integer :: nc ! size of coefficient grid
      integer :: ll,ul ! lower/upper level for transition
      integer :: itemp
      !real(4),dimension(:), pointer, contiguous :: c ! coefficients
      !real(4),dimension(:), pointer, contiguous :: dcdt ! derivatives
   end type coll_type

   ! *******************************************************************
 
    integer,parameter :: ncdi=5,nbadnell=2
    
    !********************************************************************
    real(8), parameter :: c0=((Eryd/sqrt(em)) * pi*Rbohr**2) * sqrt(8.0/(pi*bk))
 
    real(8), parameter :: sumscl=1.0d0
 
    ! the following paramters are for Burgess rates.
    ! actelc=1 to be replaced later by routine that computes number
    ! of active electrons from term information
    real(8), parameter :: cbar=2.3d0, x100=5.0d2
    real(8), parameter :: actelc=1.d0
 
    character(len=8),dimension(13),parameter :: theprocesses = &
          (/'OHMEGA  ','CI      ','CH      ','CE      ', &
             'FIXRAT  ','SHULL82 ','AR85-CEA',&
             'BURGESS ','AR85-CH ','AR85-CHE','AR85-CHP',&
             'CH0     ','CHI     '/)
    
    !integer :: n_coll_processes! number of keywords 
    type(coll_type),dimension(:),allocatable :: coll
     
 contains

   !====================================================================================
   SUBROUTINE gencol(coll,id,nlevel,g,ev,ion,cij,temp,ne,nstar,nh,tgrid,cgrid,dcdt)
      class(coll_type), dimension(:) :: coll
      character(len=*), intent(in) :: id
      integer, intent(in) :: nlevel
      real(8), intent(in), dimension(:) :: g, ev
      integer, intent(in), dimension(:) :: ion
      real(8), dimension(nlevel,nlevel),intent(out) :: cij
      real(8),intent(in) :: temp, ne
      real(8),dimension(:), intent(in) :: nh
      real(4),dimension(:), intent(in) :: nstar
  
      real(8) :: ct,cup,cdn,dekt,de_k,de
      real(8) :: nne,tg,stg,nsr
      real(8), dimension(:), allocatable :: nnh
      real(8), dimension(:,:) :: tgrid, cgrid, dcdt
  
      real(8) :: t4,summrs
      real(8) :: betab,ddekt,wlog,wb,dexin1,xj,fac,fxj
      real(8) :: dekti,qb,he1,he2
      real(8) :: saha, dxi, dev, tns, ek_over_temp
      integer :: ul,ll,i,j,k,z, nhlevels, ii, jj, nc, itemp
      !---------------------------------------------------------------------
      ! get local copies for ease of notation
      tg=temp
      stg=sqrt(tg)
      nne=ne
      nnh=nh(:)
      nhlevels = size(nh, 1)
      saha = 1.5*log(temp) - log(saha_const*ne)
      dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
      ek_over_temp = e_k / temp
  
      ! loop over all the transitions
      do i=1, size(coll, 1)
  
         cup=0.0d0
         cdn=0.0d0
         ll=coll(i)%ll
         ul=coll(i)%ul

         dev = ev(ul) - ev(ll)
         if (ion(ul)>ion(ll)) then
            ii = ion(ll) * (ion(ll)-1)
            jj = ion(ul) * (ion(ul)-1)
            dev = dev - dxi*0.5*(jj-ii)
         endif

         if (nstar(ll) > 0) then
            nsr = nstar(ul) / nstar(ll)
         else
            tns = - dev*ek_over_temp + (ion(ul) - ion(ll))*saha
            nsr = g(ul)/g(ll) * exp(tns)
         endif

         select case (trim(coll(i)%key))
         case ('OHMEGA','CE','CH','CI','FIXRAT','CH0','CHI')
            nc = coll(i)%nc
            itemp = coll(i)%itemp
            call hermite_interpolation(ct,tg,tgrid(1:nc,itemp),cgrid(1:nc,i),dcdt(1:nc,i))
         end select
  
         !********************************************************************
  
         select case (trim(coll(i)%key))
         case('OHMEGA')
            !  excitation of  electron collisons with positive ions
            cdn = c0 / g(ul) * ct  * nne / stg
            cup = cdn * nsr

         case('CE')
            ! excitation of neutrals by electron collisions 
            cdn = g(ll)/g(ul) * ct * nne * stg
            cup = cdn * nsr

         case('CH')
            ! collisions with neutral hydrogen)
            ! these are downward rates (cdn); cup are upward rates
            cdn = ct * nnh(1)
            cup = cdn * nsr

         case('CHI')
            ! b-f collisions with neutral hydrogen
            ! these are downward rates (cdn); cup are upward rates
            cup = ct * nnh(1)
            cdn = cup / nsr

         case('CI')
            ! ionization by electron collisions
            cup = ct * nne * stg * exp(-dev*ek_over_temp) 
            cdn = cup / nsr

         case('CH0')
            ! MB 5.04.2018 - charge transfer with H-, A+H <--> A+ + H-
            ! these are upwards rates (CT ionisation or ion-pair formation)
            ! also assumed for CI; but the difference is that CI downwards rates 
            ! are computed from detailed balance (they are bi-directional);
            ! CH0 are one-directional and the levels must be in the correct order
            !  in the atom (, i.e. ll  --> lu for ioninisation)
            ! in the QM papers, the inverse exothermic rates (neutralisation) 
            ! are derived from the detailed balance
            ! N_kj(T,E_k) = N_jk(t,E_k) * exp(-E_jk/ (k_B * T))
            cup = ct * nnh(1)
            cdn = cup / nsr
  
            !********************************************************************

         case('SHULL82')
            ! Coefficients for collisional ionization, radiative recombination, and dielectronic recombination
            ! following
            ! Shull & van Steenberg, Astrophysical Journal Supplement Series, vol. 48, Jan. 1982, p. 95-107
            ! ''The ionization equilibrium of astrophysically abundant elements''
  
            !  the idea is to reduce dielectronic recombination rate
            !  owing to collisional ionization of rydberg states
            !  following summers 1974, appleton lab rept.
            summrs=1.0d0
            ! to remove density dependence in dielectronic rate comment out the next line
            summrs=summers(id,ion,ul,nne)
  
            cdn = cgrid(3,i) * (tg/1.0d4)**(-cgrid(4,i)) &
                 +  summrs * cgrid(5,i) /tg/stg * exp(-cgrid(7,i)/tg) * &  
                 (1.d0 + cgrid(6,i) * (exp(- cgrid(8,i) /tg)))
            cdn = cdn * nne
  
            cup = cgrid(1,i) * stg * exp(-cgrid(2,i)/tg) / (1.d0 + 0.1d0 * tg / cgrid(2,i)) 
            cup = cup * nne
  
            !  3-body recombination (high density limit)
            cdn = cdn + cup / nsr
  
  
            !********************************************************************

         case('AR85-CDI')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            ! AR85-CDI is di-electronic recombination
  
            cup = 0.0d0
  
            do j=1,coll(i)%nc
  
               k = (j-1) * ncdi
  
               xj=cgrid(k+1,i) * ee / bk / tg
  
               fac=exp(-xj)*sqrt(xj)
  
               fxj=fac*( cgrid(k+2,i) + cgrid(k+3,i) * (1.d0+xj) &
                    +(cgrid(k+4,i) - xj*(cgrid(k+2,i) + cgrid(k+3,i) &
                    *(2.d0+xj)))*fone(xj)+ &
                    cgrid(k+5,i)*xj*ftwo(xj) )
  
               fac=6.69d-07/cgrid(k+1,i)/sqrt(cgrid(k+1,i))
  
               cup = cup + fac*fxj
  
            end do
  
            if (cup<0) cup=0
  
            cup=cup*nne
            cdn=cup / nsr
  
  
            !********************************************************************
  
         case('AR85-CEA')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"          
  
            ! AR85-CEA is collision auto-ionization
  
            call ar85cea(id, ion, tg,ll,ul,ct)  
            cup = cgrid(1,i) * ct * nne
  
            !********************************************************************
  
         case('AR85-CH')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            !AR85-CH is charge transfer with neutral hydrogen atom^+ + H -> atom + H^+
            if ( tg>=cgrid(1,i) .and. tg<=cgrid(2,i) )  then
               t4=tg/1.0d4
               cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * nnh(1)
            endif
  
            !********************************************************************  
  
         case('AR85-CH+')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            ! ar85-ch+ is  charge transfer with ionized hydrogen atom + H^+ -> atom^+ + H
            if ( tg>=cgrid(1,i) .and. tg<=cgrid(2,i) )  then
               t4=tg/1.0d4
               cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * nnh(nhlevels)
            endif
  
            !********************************************************************
  
         case('AR85-CHE')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            !AR85-CHE is charge transfer with neutral helium atom^+ + He -> atom + He^+
            if ( tg>=cgrid(1,i) .and. tg<=cgrid(2,i) ) then
               t4=tg/1.0d4
               call hepop(tg,sum(nnh),he1,he2)
               cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * he1 
            endif
  
            !********************************************************************
  
         case('AR85-CHP')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            ! ar85-chp or ar85-che+ is  charge transfer with ionized helium atom^+ + He^+ -> atom + He^2+
            if ( tg>=cgrid(1,i) .and. tg<=cgrid(2,i) ) then
               t4=tg/1.0d4
               call hepop(tg,sum(nnh),he1,he2)
               cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * he2 
            endif
  
            !********************************************************************
  
         case('BURGESS')
            ! Electron impact ionization of complex ions
            ! Burgess & Chidichimo, 1983, MNRAS, 203, 1269
            de=ev(ul)-ev(ll)
            z=ion(ll)-1
            betab = 0.25d0*(sqrt((1.d2*z +91.d0)/(4.d0*z+3.d0))-5.d0)
            dekt  = de*ee/bk/tg
            ddekt = -dekt/tg
            dekt  = min(x100,dekt)
            dekti = 1.d0/dekt
            wlog = log(1.d0 + dekti)
            wb = wlog**(betab/(1.d0+dekti))
            dexin1 = -exp(-dekt)/dekt*ddekt
            qb = 2.1715d-08 * cbar * (13.6d0/de)**1.5d0 * sqrt(dekt) * ei1(dekt) * wb
  
            cup = qb * actelc * nne * cgrid(1,i)
            cdn = cup/nsr
  
            !********************************************************************
  
         case('BADNELL')
            ! fits to dielectronic recombination from Badnell and co-workers.
            ! see for example Altun et al., 2006, A&A, 44 7, 1165
            cdn=0.0d0
            do j=1,coll(i)%nc
               cdn=cdn + cgrid(j+coll(i)%nc,i) * exp(-cgrid(j,i)/tg)
            end do
  
            summrs=sumscl*summers(id,ion,ul,nne) + (1.0d0-sumscl)
            cdn = cdn * nne * tg**(-1.5d0)* summrs
  
         end select
  
         !********************************************************************
  
         cij(ul,ll) = cij(ul,ll)  + cup ! upward rate Cll->ul 
         cij(ll,ul) = cij(ll,ul)  + cdn ! downward rate Cul->ll
  
      end do ! transition loop 
  
   END SUBROUTINE gencol


   !====================================================================================
   SUBROUTINE classic_col(g,ev,ion,temp,ne,nh,cdata,cij,ecol,hcol)
      real(8), intent(in), dimension(:) :: g, ev
      integer, intent(in), dimension(:) :: ion
      real(8), intent(in) :: temp, ne
      real(8), dimension(:), intent(in) :: nh
      real(4), dimension(:), intent(in) :: cdata
      real(8), dimension(:,:),intent(inout) :: cij
      logical, dimension(:),intent(inout) :: ecol, hcol
  
      real(8) :: c,cup,cdn, EK_over_temp
      real(8) :: nne,tg,stg,nsr
      real(8), dimension(:), allocatable :: nnh
  
      real(8) :: saha, dxi, dev, tns, e0, ex, u0, f, tmp, isqt, fac
      real(8) :: seaton, allen, regemorter
      integer :: nlevel, ul,ll,i,j,k,z, nhlevels, ii, jj, nc, icont, tri, n
      logical :: cont
      logical, save :: first_time=.true.
      !---------------------------------------------------------------------
      nlevel = size(ev,1)
      n = 0
      saha = 1.5*log(temp) - log(saha_const*ne)
      dxi  = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
      isqt = 1 / sqrt(temp)
      allen = 5.465e-11 * Rydk * isqt
      regemorter = 14.5 * allen * RydV
      seaton = 1.55e13 * isqt

      EK_over_temp = e_k / temp

      do i=1, nlevel
         do j=i+1, nlevel
            tri = (j-1)*(j-2)/2 + i
            fac = cdata(tri)
            cont = ion(j) > ion(i)
            if (.not. ecol(tri)) then
               e0 = ev(j) - ev(i)
               if (cont) then
                  ii = ion(j) * (ion(j)-1)
                  jj = ion(i) * (ion(i)-1)
                  e0 = e0 - dxi*0.5*(jj-ii)
                  u0 = e0 * EK_over_temp
                  ex = exp(-u0)
                  if (.not. ecol(tri)) then
                     ! ionization by electron collisions
                     tmp = min(0.1*ion(i), 0.3)
                     c = ne * fac * tmp / u0
                     cup = c * ex
                     cdn = c * g(i) / g(j)
                  endif
               else
                  u0 = e0 * EK_over_temp
                  ex = exp(-u0)
                  if (.not. ecol(tri)) then
                     ! excitation of neutrals by electron collisions 
                     if (fac > 0) then
                        ! Van Regemorter
                        f = max(fac, 1e-5)
                        tmp = max(0.2, 0.276*expin(u0)/ex)
                        c = ne * regemorter * f * tmp / e0
                     else
                        ! Allen
                        c = ne * allen / g(i)
                     endif
                     cup = c * ex
                     cdn = c * g(i) / g(j)
                  endif
               endif

               cij(j,i) = cij(j,i) + cup ! upward rate
               cij(i,j) = cij(i,j) + cdn ! downward rate
               n = n + 1
            endif
            if (cont) exit
         end do
      end do

      if (first_time) then
         first_time = .false.
         write(stdout, *) n, ' collisional transitions filled with classical recipes.'
      endif
   END SUBROUTINE classic_col
 
   !====================================================================================
   logical function is_valid_process(key)
    character(len=8),intent(in) :: key
    integer :: i
    !---------------------------------------------------------------------
    is_valid_process=.false.
    do i=1,size(theprocesses)
       if ( key(1:len_trim(key)) == trim(theprocesses(i)))  then
          is_valid_process=.true.
          exit
       end if
    end do
 
   end function is_valid_process
 
   !====================================================================================
   pure function summers(id, ion, j,edens)
      !*********************************************************************
      !
      !  DENSITY SENSITIVE DIELECTRONIC RECOMBINATION
      !
      !  the term adi may be multiplied by a density-sensitive factor
      !  if needed- this is crucial for Li and B-like ions colliding with
      !  impacting electrons.
      !  This simple formulation was derived from a study of the dependence of
      !  the dielectronic "bump" in the figures of Summers 1974 
      !  (Appleton Laboratory internal memo), and fitting according to the
      !  parameter Ne / z^7
      !  This should be accurate to typically +/- 0.1 in log in regions
      !  where it matters.  Worse case is e.g. C like Neon where it underestimates
      !  density factor by maybe 0.25 in log.
      !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  June 24, 2006 changes begin P.G.Judge 
      !  original (pre MAY 2006) code
      !         !  define rho = nne/ z^7 where z is charge on recombining ion
      !         rho=10.^(alog10(nne) - 7.* alog10(charge))
      !         rho0=2.e3
      !         if(isos eq 'LI' or isos eq 'NA' or isos eq 'K') then rho0 = 3.e1      
      !         ne_factor = 1./(1. + rho/rho0)^0.14
      !         print,'ne_factor',ne_factor
      !
      ! June 24, 2006, more accurate version. 
      ! get the row and column of the recombined ion's isoelectronic 
      ! sequence  in the periodic table
      ! the following parameters mimic the tables 1-19 of
      ! H. Summers' Appleton Lab Report 367, 1974
      ! 
      ! Mar  9 2012 jorritl: stolen from Phil Judge's DIPER IDL package
      !
      !*********************************************************************

      real(8) :: summers  
      character(len=*), intent(in) :: id
      integer, intent(in), dimension(:) :: ion
      integer, intent(in) :: j
      real(8), intent(in) :: edens

      real(8) :: z,rho0,rhoq,x,beta

      real(8),parameter :: ex=exp(1.0d0)

      integer :: anr,isoseq,row,col

      summers=1.0d0
      anr=atomic_number(id)   

      z=real(ion(j)-1) ! charge of recombining ion
      isoseq=anr-z+1 ! iso-electronic sequence of recombined ion
      call rowcol(isoseq,row,col) ! row and column in periodic table
      rhoq=edens/z**7
      x=one_third*(0.5*z + (col-1))*row
      beta= -0.2/log(x+ex)
      rho0= 30.0+50.0*x
      summers= (1.0d0 + rhoq/rho0)**beta

      return
   end function summers

   !*********************************************************************

   pure function fone(x)
      ! calculates f1(x) needed for collisional rates of arnaud and rothenflug 
      ! for large values of x, use asymptotic limit
      real(8) :: fone
      real(8),intent(in) :: x
      real(8),parameter :: split=50.0
      if (x<=split) then 
      fone = exp(x)*ei1(x)
      else
      fone= 1./x
      endif
      return
   end function fone

   !**********************************************************************

   function ftwo(x)
      USE mconst,only:pi,euler_const
      implicit none
      ! new routine 24-jan-1994  p. judge
      ! calculates f2(x) needed for collisional rates of arnaud and rothenflug 
      ! if argument x is < break, then use the eqns 5 and 6 of hummer (1983, jqsrt, 30 281)
      ! break given by hummer as 4.0
      !
      ! note that the suggested polynomial expansion of arnaud and rothenflug
      ! fails for arguments < break- see hummer's original paper.
      ! it is unclear if the tabulated ionization fractions of arnaud and rothenflug
      ! contain this error.   
      real(8) :: ftwo
      real(8),intent(in) :: x 

      real(8),parameter :: break=4.0
      real(8),dimension(15),parameter :: p= &
         (/1.0,2.1658e+02,2.0336e+04,1.0911e+06,3.7114e+07, &
         8.3963e+08,1.2889e+10,1.3449e+11,9.4002e+11,4.2571e+12,&
         1.1743e+13,1.7549e+13,1.0806e+13,4.9776e+11,0.0000/)

      real(8),dimension(15),parameter :: q= &
         (/1.0,2.1958e+02,2.0984e+04,1.1517e+06,4.0349e+07,&
         9.4900e+08,1.5345e+10,1.7182e+11,1.3249e+12,6.9071e+12,&
         2.3531e+13,4.9432e+13,5.7760e+13,3.0225e+13, 3.3641e+12/)

      real(8) :: px,xfact,qx,f0x,count,fact,term

      integer :: i

      if (x .ge. break) then 
      px=p(1)
      xfact=1.0
      do i=2,15
         xfact = xfact/x
         px=px+ p(i)*xfact
      end do
      qx=q(1)
      xfact=1.0
      do  i=2,15
         xfact = xfact/x
         qx=qx+ q(i)*xfact
      end do
      ftwo=px/qx/x/x
      else

      !  hummer's eqns 5 and 6
      f0x = pi*pi/12.
      term=1.0
      count=0.0
      fact=1.0
      xfact=1.0
      do while (abs(term/f0x) .gt. 1.e-8)
         count=count+1.0
         fact = fact*count
         xfact=xfact*(-x)
         term = xfact/count/count/fact
         f0x=f0x+term
         if(count .gt. 100.) then
            stop 'ftwo: too many iterations'
         endif
      enddo
      ftwo=exp(x)*((log(x)+euler_const)*(log(x)+euler_const)*0.5 + f0x)
      endif
      return
   end function ftwo


   !*********************************************************************

   pure integer function atomic_number(atomic_symbol)
      USE mconst,only:element_name
      implicit none
      CHARACTER(len=*), intent(in) :: atomic_symbol
      integer :: i

      atomic_number=-1
      do i=1,size(element_name)
         if (uppercase(trim(atomic_symbol),2)==element_name(i)) then
               atomic_number=i
               return
         end if
      end do
   end function atomic_number

   !*********************************************************************

   PURE SUBROUTINE rowcol(an,row,col)
      !**********************************************************************
      !
      ! returns the row and column of the element in
      ! the periodic table- the row is the "period", but 
      ! the column is just the order in which elements occur
      ! in the row.  it is not the "group"
      !
      ! Document name: rowcol.pro
      ! Created by:    judge, June 24, 2006
      !
      ! Last Modified: Sat Jun 24 13:08:00 2006 by judge (judge) on edlen.local
      !
      ! Mar  9 2012 jorritl: stolen from Phil Judge's DIPER IDL package.
      !
      !**********************************************************************

      integer, intent(in)  :: an
      integer, intent(out) :: row,col
      integer :: i
      integer,dimension(7),parameter :: istart=(/0,2,10,18,36,54,86/)

      do i=1,6
         if ( an>=istart(i) .and. an<=istart(i+1) ) then
               row=i
               exit
         endif
      end  do
      col=an-istart(row)
   END SUBROUTINE rowcol

   !**********************************************************************

   SUBROUTINE ar85cea(id, ion, tg,ilo,ihi,cup)
      !  new routine for computing collisional autoionization 
      !  rates using formalism and formulae from arnaud and rothenflug 1985
      ! 
      !  ar85cea 94-02-22  new routine: (philip judge)
      ! 
      !          96-03-07  modifications: (philip judge)
      !                    bug fixed: cup initialized to zero
      !
      character(len=*), intent(in) :: id
      integer, dimension(:), intent(in) :: ion
      integer,intent(in) :: ilo,ihi
      real(8),intent(in) :: tg
      real(8), intent(out) :: cup
   
      character(len=2) ::  isoseq_element, element
   
      integer :: z,atomic_charge
      real(8) :: zz,bkt,a,b,zeff,iea,y,f1y
   
      bkt=tg*bk/ee
   
      !  initialize output to zero 
      y=0.0
      f1y=0.0
      cup=0.0
   
      !  find element
      element = id(1:2)
      z=atomic_number(element)
      zz=real(z)
   
      if (z<1 .or. z>92) then 
         write(stdout, *) 'ar85-cea:  atomic number = ', z, ' for ', element
         write(stdout, *) 'ar85-cea:  no autoionization included'
         return
      endif
   
      !  find iso-electronic sequence
      atomic_charge=ion(ihi)-1
      isoseq_element=lowercase(element_name(z-atomic_charge+1),2)
   
      select case(isoseq_element)
      case('li') 
         ! lithium sequence
   
         b=1.0 / (1.0 + 2.0d-4*zz**3)
         zeff = zz-1.3 
         iea=13.6*((zz-0.835)**2 - 0.25*(zz-1.62)**2)
         y=iea/bkt
         f1y=fone(y)
         cup=1.60e-07*1.2*b/zeff**2/sqrt(bkt)*exp(-y)* &
               (2.22*f1y+0.67*(1.-y*f1y)+0.49*y*f1y+1.2*y*(1.-y*f1y))
   
         ! special cases:
         if(element == 'c ') cup = cup*0.6    ! C IV
         if(element == 'n ') cup = cup*0.8 ! N V
         if(element == 'o ') cup = cup*1.25   ! O VI
   
      case('na')
         ! sodium sequence
         if (z <= 16) then 
            a=2.8e-17*(zz-11.)**(-0.7)
            iea=26.*(zz-10.)
            y=iea/bkt
            f1y=fone(y)
            cup= 6.69e+7 * a *iea/sqrt(bkt) * exp(-y)*(1. - y*f1y)
         else if ( (z>=18) .and. (z<=28) ) then 
            a=1.3e-14*(zz-10.)**(-3.73)
            iea=11.*(zz-10.)*sqrt(zz-10.)
            y=iea/bkt
            f1y=fone(y)
            cup= 6.69e+7 * a *iea/sqrt(bkt) * exp(-y)  &
                  *(1. - 0.5*(y -y*y + y*y*y*f1y))
         else 
            return
         endif
   
      case('mg','al','si','p ','s ')
   
         if (isoseq_element=='mg') iea=10.3*(zz-10.)**1.52
         if (isoseq_element=='al') iea=18.0*(zz-11.)**1.33
         if (isoseq_element=='si') iea=18.4*(zz-12.)**1.36
         if (isoseq_element=='p ') iea=23.7*(zz-13.)**1.29
         if (isoseq_element=='s ') iea=40.1*(zz-14.)**1.1
         a=4.0e-13/zz**2/iea
         y=iea/bkt
         f1y=fone(y)
         cup= 6.69e+7 * a *iea/sqrt(bkt) * exp(-y)   &
               *(1. - 0.5*(y -y*y + y*y*y*f1y))
      case default
         ! do nothing
      end select
   
      !  special cases
   
      if(element=='ca' .and. atomic_charge==0) then 
         a = 9.8e-17
         iea= 25.
         b=1.12         
         cup=6.69e+7*a *iea/sqrt(bkt) * exp(-y)*(1. + b*f1y) 
      else if ( element=='ca' .and. atomic_charge==1 ) then 
         a = 6.0e-17
         iea= 25.
         b=1.12         
         cup=6.69e+7*a *iea/sqrt(bkt) * exp(-y)*(1. + b*f1y)
      else if (element=='fe' .and. atomic_charge==3) then 
         a = 1.8e-17
         iea= 60.
         b=1.0         
         cup=6.69e+7*a *iea/sqrt(bkt) * exp(-y)*(1. + b*f1y)
      else if ( element=='fe' .and. atomic_charge==4 ) then 
         a = 5.0e-17
         iea= 73.
         b=1.0         
         cup=6.69e+7*a *iea/sqrt(bkt) * exp(-y)*(1. + b*f1y)
      endif
   END SUBROUTINE ar85cea

   !*********************************************************************

   SUBROUTINE hepop(tg,hpop,he1,he2)
      USE mmath,only:linear1d
      implicit none
      !Ratio of Helium ionization fraction from table IV of Arnaud and Rothenflug
      !ABHE is abundance of Helium- assumed =0.1
      real(8),intent(in) :: tg,hpop
      real(8),intent(out) :: he1,he2

      real(8) :: tlog,dum

      real(8), parameter :: abhe=0.1
      integer, parameter :: nt=19

      real(8) ,dimension(nt),parameter :: &
            tt=(/3.50, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, &
            4.70, 4.80, 4.90, &
            5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70/)

      real(8),dimension(nt),parameter :: &
            ie1=(/0.0, 0.0 ,0.0 ,0.0 ,0.0 ,-0.07,-0.51,-1.33,-2.07, &
            -2.63,-3.20,-3.94,-4.67,-5.32,-5.90,-6.42, &
            -6.90,-7.33,-7.73/)

      real(8),dimension(nt),parameter :: &
            ie2=(/-20.0,-9.05,-6.10,-3.75,-2.12,-0.84,-0.16, &
            -0.02,-0.01,-0.05,-0.34,-0.96,-1.60,-2.16,-2.63, & 
            -3.03,-3.38,-3.68,-3.94/)

      tlog=log10(tg)
      call linear1d(nt,tt,ie1,tlog,dum,.true.)
      he1=hpop * abhe * 10.**dum!10.**splin(tlog,tt,ie1,nt,nt)
      call linear1d(nt,tt,ie1,tlog,dum,.true.)
      he2=hpop * abhe * 10.**dum    !10.**splin(tlog,tt,ie2,nt,nt)

   END SUBROUTINE hepop

   !*********************************************************************


   !====================================================================================
   ELEMENTAL FUNCTION expin (x) result(ex1)
      real(8), intent(in) :: x
      real(8) :: ex, ex1
      real(8), parameter, dimension(6) :: &
      a = [-44178.5471728217,57721.7247139444,9938.31388962037,    &
            1842.11088668000,101.093806161906,5.03416184097568]
      real(8), parameter, dimension(5) :: &
      b = [76537.3323337614,32597.1881290275,6106.10794245759,     &
            635.419418378382,37.2298352833327]
      real(8), parameter, dimension(7) :: &
      c = [ 4.65627107975096e-7,.999979577051595,9.04161556946329, &
           24.3784088791317,23.0192559391333,6.90522522784444,     &
             .430967839469389]
      real(8), parameter, dimension(6) :: &
      d = [10.0411643829054,32.4264210695138, 41.2807841891424,    &
           20.4494785013794,3.31909213593302, .103400130404874]
      real(8), parameter, dimension(7) :: &
      e = [-.999999999998447,-26.6271060431811,-241.055827097015,  &
           -895.927957772937,-1298.85688746484,-545.374158883133,  &
            -5.66575206533869]
      real(8), parameter, dimension(6) :: &
      f = [28.6271060422192,292.310039388533,1332.78537748257,     &
           2777.61949509163,2404.01713225909, 631.657483280800]
      !---------------------------------------------------------------------
      ex = exp(-x)
      if (x > 4) then
         ex1 = (ex+ex*(e(1)+(e(2)+(e(3)+(e(4)+(e(5)+(e(6)+e(7)/x)/x)/x)/x)/x)/x)/  &
               (x+f(1)+(f(2)+(f(3)+(f(4)+(f(5)+f(6)/x)/x)/x)/x)/x))/x
      elseif (x > 1) then
         ex1 = ex*(c(7)+(c(6)+(c(5)+(c(4)+(c(3)+(c(2)+c(1)*x)*x)*x)*x)*x)*x)/      &
               (d(6)+(d(5)+(d(4)+(d(3)+(d(2)+(d(1)+x)*x)*x)*x)*x)*x)
      elseif (x > 0) then
         ex1 = (a(1)+(a(2)+(a(3)+(a(4)+(a(5)+a(6)*x)*x)*x)*x)*x)/                  &
               (b(1)+(b(2)+(b(3)+(b(4)+(b(5)+x)*x)*x)*x)*x) - log(x)
      else
         ex1 = 0.d0
      endif
   END FUNCTION expin

 END MODULE coll_mod