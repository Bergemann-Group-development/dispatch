MODULE coll_mod
    USE io_mod
    USE os_mod
    USE mmath_mod,only:ei1,hermite_interpolation,hermite_derivative
    USE const_mod
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
   SUBROUTINE gencol(coll,id,nlevel,g,ev,ion,cij,temp,ne,nstar,nh1,tgrid,cgrid,dcdt)
      class(coll_type), dimension(:) :: coll
      character(len=*), intent(in) :: id
      integer, intent(in) :: nlevel
      real(8), intent(in), dimension(:) :: g, ev
      integer, intent(in), dimension(:) :: ion
      real(8), intent(in), dimension(:) :: nstar
      real(8), intent(out), dimension(nlevel,nlevel) :: cij
      real(8), intent(in) :: temp, ne, nh1
  
      real(8) :: ct,cup,cdn,dekt,de_k,de
      real(8) :: stg,nsr
      real(8), dimension(:,:) :: tgrid, cgrid, dcdt
  
      real(8) :: t4,summrs
      real(8) :: betab,ddekt,wlog,wb,dexin1,xj,fac,fxj
      real(8) :: dekti,qb,he1,he2
      real(8) :: saha, dxi, dev, tns, ek_over_temp
      integer :: ul,ll,i,j,k,z, nhlevels, ii, jj, nc, itemp
      !---------------------------------------------------------------------
      ! get local copies for ease of notation
      stg=sqrt(temp)
      saha = 1.5*log(temp) - log(saha_const*ne)
      dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
      ek_over_temp = e_k / temp
  
      ! loop over all the transitions
      do i=1, size(coll, 1)
         if (coll(i)%nc == 0) cycle
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
            call hermite_interpolation(ct,temp,tgrid(1:nc,itemp),cgrid(1:nc,i),dcdt(1:nc,i))
            ct = exp(ct)
         end select
  
         !********************************************************************
  
         select case (trim(coll(i)%key))
         case('OHMEGA')
            !  excitation of  electron collisons with positive ions
            cdn = c0 / g(ul) * ct * ne / stg
            cup = cdn * nsr

         case('CE')
            ! excitation of neutrals by electron collisions 
            cdn = g(ll)/g(ul) * ct * ne * stg
            cup = cdn * nsr

         case('CH')
            ! collisions with neutral hydrogen)
            ! these are downward rates (cdn); cup are upward rates
            cdn = ct * nh1
            cup = cdn * nsr

         case('CHI')
            ! b-f collisions with neutral hydrogen
            ! these are downward rates (cdn); cup are upward rates
            cup = ct * nh1
            cdn = cup / nsr

         case('CI')
            ! ionization by electron collisions
            cup = ct * ne * stg * exp(-dev*ek_over_temp) 
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
            cup = ct * nh1
            cdn = cup / nsr
  
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
  
               xj=cgrid(k+1,i) * ee / bk / temp
  
               fac=exp(-xj)*sqrt(xj)
  
               fxj=fac*( cgrid(k+2,i) + cgrid(k+3,i) * (1.d0+xj) &
                    +(cgrid(k+4,i) - xj*(cgrid(k+2,i) + cgrid(k+3,i) &
                    *(2.d0+xj)))*fone(xj)+ &
                    cgrid(k+5,i)*xj*ftwo(xj) )
  
               fac=6.69d-07/cgrid(k+1,i)/sqrt(cgrid(k+1,i))
  
               cup = cup + fac*fxj
  
            end do
  
            if (cup<0) cup=0
  
            cup=cup*ne
            cdn=cup / nsr
  
            !********************************************************************
  
         case('AR85-CH')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            !AR85-CH is charge transfer with neutral hydrogen atom^+ + H -> atom + H^+
            if ( temp>=cgrid(1,i) .and. temp<=cgrid(2,i) )  then
               t4=temp/1.0d4
               cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * nh1
            endif
  
            !********************************************************************  
  
         !case('AR85-CH+')
         
         !   ! The AR85**** rates are based on 
         !   ! Arnaud & Rothenflug
         !   ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
         !   ! "An updated evaluation of recombination and ionization rates"
         
         !   ! ar85-ch+ is  charge transfer with ionized hydrogen atom + H^+ -> atom^+ + H
         !   if ( temp>=cgrid(1,i) .and. temp<=cgrid(2,i) )  then
         !      t4=temp/1.0d4
         !      cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * nh(nhlevels)
         !   endif
  
            !********************************************************************
  
         !case('AR85-CHE')
  
         !   ! The AR85**** rates are based on 
         !   ! Arnaud & Rothenflug
         !   ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
         !   ! "An updated evaluation of recombination and ionization rates"
  
         !   !AR85-CHE is charge transfer with neutral helium atom^+ + He -> atom + He^+
         !   if ( temp>=cgrid(1,i) .and. temp<=cgrid(2,i) ) then
         !      t4=temp/1.0d4
         !      call hepop(temp,sum(nh),he1,he2)
         !      cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * he1 
         !   endif
  
            !********************************************************************
  
         !case('AR85-CHP')
         
         !   ! The AR85**** rates are based on 
         !   ! Arnaud & Rothenflug
         !   ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
         !   ! "An updated evaluation of recombination and ionization rates"
         
         !   ! ar85-chp or ar85-che+ is  charge transfer with ionized helium atom^+ + He^+ -> atom + He^2+
         !   if ( temp>=cgrid(1,i) .and. temp<=cgrid(2,i) ) then
         !      t4=temp/1.0d4
         !      call hepop(temp,sum(nh),he1,he2)
         !      cdn = cgrid(3,i) * 1d-9 * t4**cgrid(4,i) * (1.0d0+cgrid(5,i)*exp(cgrid(6,i)*t4)) * he2 
         !   endif
  
            !********************************************************************
  
         case('BURGESS')
            ! Electron impact ionization of complex ions
            ! Burgess & Chidichimo, 1983, MNRAS, 203, 1269
            de=ev(ul)-ev(ll)
            z=ion(ll)-1
            betab = 0.25d0*(sqrt((1.d2*z +91.d0)/(4.d0*z+3.d0))-5.d0)
            dekt  = de*ee/bk/temp
            ddekt = -dekt/temp
            dekt  = min(x100,dekt)
            dekti = 1.d0/dekt
            wlog = log(1.d0 + dekti)
            wb = wlog**(betab/(1.d0+dekti))
            dexin1 = -exp(-dekt)/dekt*ddekt
            qb = 2.1715d-08 * cbar * (13.6d0/de)**1.5d0 * sqrt(dekt) * ei1(dekt) * wb
  
            cup = qb * actelc * ne * cgrid(1,i)
            cdn = cup/nsr
  
         end select
  
         !********************************************************************
  
         cij(ul,ll) = cij(ul,ll)  + cup ! upward rate Cll->ul 
         cij(ll,ul) = cij(ll,ul)  + cdn ! downward rate Cul->ll
  
      end do ! transition loop 
  
   END SUBROUTINE gencol
 
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
      USE const_mod,only:pi,euler_const
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

   !*********************************************************************

   SUBROUTINE hepop(tg,hpop,he1,he2)
      USE mmath_mod,only:linear1d
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

   !====================================================================================
   ELEMENTAL FUNCTION gaunt (nu, n)
      real(8), intent(in) :: nu
      integer, intent(in) :: n
      real(8) :: gaunt
      real(8), parameter, dimension(6) :: &
      a1 = [ 1.078e00, 1.0926e0, 1.0983e0, 1.0954e0, 1.0912e0, 1.0925e0], &
      b1 = [-8.754e14,-2.019e14,-9.45e13, -5.188e13,-3.2e13,  -2.331e13], &
      c1 = [-1.791e29, 1.836e28, 9.177e27, 3.552e27, 1.576e27, 9.325e26], &
      a2 = [ 0.798e00, 0.768e00, 0.793e00, 0.831e00, 0.758e00, 0.79e00],  &
      b2 = [ 5.358e15, 6.242e15, 5.48e15,  4.094e15, 6.633e15, 5.808e15], &
      c2 = [-3.484e31,-3.208e31,-2.318e31,-1.43e31, -3.32e31, -2.844e31]
      !---------------------------------------------------------------------
      if (n > 6) then
         gaunt = 1
      else
         if (nu <= 1e16) then
            gaunt = a1(n) + (b1(n) + c1(n)/nu) / nu
         else
            gaunt = a2(n) + (b2(n) + c2(n)/nu) / nu
         endif
      endif
   END FUNCTION gaunt

   !====================================================================================
   ELEMENTAL FUNCTION Kramer (nu, dev) result(hyd)
      real(8), intent(in) :: nu, dev
      real(8) :: nl, hyd
      !---------------------------------------------------------------------
      if (dev > 0) then
         nl = sqrt(RydV/dev)
         hyd = 2.815e29 * gaunt(nu,ceiling(nl)) / (nl**5 * nu**3)
      else
         hyd = 0.
      endif
   END FUNCTION Kramer
 END MODULE coll_mod