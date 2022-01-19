module coll_mod
    use m3io_mod
    use mmath,only:ei1,hermite_interpolation,hermite_derivative
    use mconst
    implicit none


   ! *******************************************************************

   type coll_type
      character(len=16) :: key ! transition type 
      integer :: nc ! size of coefficient grid
      integer :: ll,ul ! lower/upper level for transition
      real(8),dimension(:), allocatable :: c ! coefficients
      real(8),dimension(:), allocatable :: dcdt ! derivatives
      real(8),dimension(:), allocatable :: t ! temperatures
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
 
    character(len=8),dimension(14),parameter :: theprocesses = &
          (/'OHMEGA  ','CI      ','CH      ','CE      ',&
             'FIXRAT  ','SHULL82 ','AR85-CDI','AR85-CEA',&
             'BURGESS ','AR85-CH ','AR85-CHE','AR85-CHP',&
             'BADNELL ','CH0     '/)
    
    !integer :: n_coll_processes! number of keywords 
    type(coll_type),dimension(:),allocatable :: coll
     
 contains
 
 
   subroutine read_collisions(file_unit, coll, max_n_temp, nlevel)
      class(coll_type), target :: coll(:)
      integer :: file_unit,ntemp, max_n_temp, nlevel
      character(len=8) :: key
      class(coll_type), pointer :: col
   
      real(8), dimension(:),allocatable :: tgrid, cgrid
      integer ::  j,i,u,l,ios
   
      allocate(tgrid(max_n_temp))
      allocate(cgrid(max_n_temp))
   
      i=0
      do while (i<=size(coll)) !le to catch the end keyword
   
         !  read the keyword 
         call nextentry(file_unit)
         read(file_unit,fmt='(a8)',iostat=ios) key
         if (ios<0) call stop('read_collisions: keyword "END" not found in collision file')
   
         key=uppercase(adjustl(key),len_trim(key))
         
         ! this is historic because multi only allowed 8 character process strings
         if (key(1:len_trim(key))=='AR85-CHE+') key(1:8)='AR85-CHP'
         
         if (is_valid_process(key)) then
            i=i+1
            col => coll(i)
         end if
   
         select case(key(1:len_trim(key)))             
         case('END')     
            write(*,*) "read_collisions: read all collisions."
            i=i+1
            exit
   
         case('TEMP')
            call nextentry(file_unit)
            read(file_unit,*) ntemp, (tgrid(j),j=1,ntemp)
            if (ntemp<2) call stop("read_collisions: ntemp<2")
   
         ! MB 5.04.2018 - added CH0 keyword
         case('OHMEGA','CE','CH','CI','FIXRAT','CH0')
            col%key=key
            call nextentry(file_unit)
            allocate(col%t(ntemp))
            allocate(col%c(ntemp))
            allocate(col%dcdt(ntemp))
            col%nc=ntemp    
            col%t = tgrid(1:ntemp)
            read(file_unit,*) col%ll,col%ul,(col%c(j),j=1,ntemp)  
   
            ! set derivatives for 3rd order hermite interpolation
            call hermite_derivative(col%dcdt, col%t, col%c)
   
         case('AR85-RR') 
            col%key=key
            call nextentry(file_unit)
            col%nc=2
            allocate(col%c(col%nc))
            read(file_unit,*) col%ll,col%ul,(col%c(j),j=1,2)
   
         case('AR85-CEA','BURGESS') 
            col%key=key
            call nextentry(file_unit)
            col%nc=1
            allocate(col%c(col%nc))
            read(file_unit,*) col%ll, col%ul, col%c(1)
   
         case('AR85-CH','AR85-CH+','AR85-CHE','AR85-CHP') 
            col%key=key
            call nextentry(file_unit)
            col%nc=6
            allocate(col%c(col%nc))
            read(file_unit,*) col%ll,col%ul,(col%c(j),j=1,col%nc)
   
         case('SHULL82') 
            col%key=key
            call nextentry(file_unit)
            col%nc=8
            allocate(col%c(col%nc))
            read(file_unit,*) col%ll,col%ul,(col%c(j),j=1,col%nc)
   
         case('BADNELL')
            col%key=key
            call nextentry(file_unit)
            read(file_unit,*) col%ll,col%ul, col%nc
            allocate(col%c(nbadnell*col%nc))
            call nextentry(file_unit)
            read(file_unit,*) (col%c(j),j=1,col%nc)
            call nextentry(file_unit)
            read(file_unit,*) (col%c(j+col%nc),j=1,col%nc)
   
         case('AR85-CDI') 
            col%key=key
            call nextentry(file_unit)
            read(file_unit,*)  col%ll,col%ul,col%nc
            allocate(col%c(ncdi*col%nc))
            do l=0,col%nc-1
               call nextentry(file_unit)
               read(file_unit,*) (col%c(j+ncdi*l),j=1,ncdi)
            end do
   
         case default
            if (.not.is_valid_process(key)) call stop('read_collisions: unrecognized keyword '//key)
         end select
   
         ! sanity check on levels
         if (is_valid_process(key)) then
            l = col%ll
            u = col%ul
            
            if (l==u) call stop("ll==ul for process " // col%key // " " // string3(l) // " " // string3(u))
   
            if (l>u) then
               col%ul=l
               col%ll=u
            endif
   
            ! within level range
            if(l<0.or.u>nlevel) then
               write(*,*) 'key=',col%key
               write(*,*) 'll=',l
               write(*,*) 'ul=',u
               write(*,*) 'nk=', nlevel
               call stop('read_collisions:  level index outside range (1,nk)')
            endif
         end if
   
      end do
   
   end subroutine read_collisions
   
   
   subroutine get_ncoll(file_unit, crout, n_coll_processes, max_n_temp)
      integer :: file_unit, ntemp
      character(len=8) :: key
      integer ::  ios, strl
      integer, intent(out) :: n_coll_processes, max_n_temp
      character(len=*) :: crout
      character(len=20) :: txt
      character(len=26),parameter :: upcase_alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
   
      n_coll_processes = 0
      max_n_temp = 0
   
      ! first run through to determine
      ! n_coll_processes=number of keywords
      ! mtgrd=size of the largest temperature grid
   
      do while (.true.)
         call nextentry(file_unit)
         read(file_unit,fmt='(a8)',iostat=ios) key
         if (ios<0) call stop('read_collisions: keyword "END" not found in collision file')
         key=uppercase(adjustl(key),len_trim(key))
   
         ! read next non-comment line if key doesnt start with a letter
         if (scan(key(1:1),upcase_alphabet)==0) cycle 
   
         ! we found a keyword, so interpret it
         select case(key(1:len_trim(key)))
         case('END')       
            exit
   
         case('TEMP')
            ! keep track of max size of temperature grid
            call nextentry(file_unit)
            read(file_unit,*) ntemp
            if(ntemp>max_n_temp) max_n_temp=ntemp
         case default
            if (.not.is_valid_process(key)) then
               call stop('read_collisions: unrecognized keyword '//key)
            end if
            n_coll_processes=n_coll_processes+1
         end select
      end do
   
      rewind(file_unit)
   
      do while (.true.)
         call nextentry(file_unit)
         read(file_unit,*) txt 
         txt=adjustl(txt)
         strl=len_trim(txt)
         txt=lowercase(txt,strl)
         if (txt==crout(1:strl)) exit
      enddo
   
   end subroutine get_ncoll


   !********************************************************************
   ! MOD PE: Add temp,ne,nh,nstar as input parameters
   subroutine gencol(coll,id,nlevel,g,ev,ion,cij,ix,iy,iz,temp,ne,nstar,nh)
    !use atmos_var
    !use atom_var
    use mconst
    use mclu
    use mem_mod
    implicit none

    class(coll_type), dimension(:) :: coll
    character(len=*), intent(in) :: id
    integer, intent(in) :: nlevel
    real(8), intent(in), dimension(:) :: g, ev
    integer, intent(in), dimension(:) :: ion
    real(8), dimension(nlevel,nlevel),intent(out) :: cij
    integer,intent(in) :: ix,iy,iz
    real(8),intent(in), dimension(:,:,:) :: temp, ne
    real(8),dimension(:,:,:,:), intent(in) :: nstar, nh

    real(8) :: ct,cup,cdn,dekt,de_k,de
    real(8) :: nne,tg,stg,nsr
    real(8), dimension(:), allocatable :: nnh

    real(8) :: t4,summrs
    real(8) :: betab,ddekt,wlog,wb,dexin1,xj,fac,fxj
    real(8) :: dekti,qb,he1,he2
    integer :: ul,ll,i,j,k,z, nhlevels

    character(len=6),parameter :: rname='gencol'

    cij=0.0d0

    ! get local copies for ease of notation
    tg=temp(ix,iy,iz)
    stg=sqrt(tg)
    nne=ne(ix,iy,iz)
    nnh=nh(ix,iy,iz,:)
    nhlevels = size(nh,4)

    ! loop over all the transitions
    do i=1, size(coll, 1)

       cup=0.0d0
       cdn=0.0d0
       ll=coll(i)%ll
       ul=coll(i)%ul
       nsr=nstar(ix,iy,iz,ul)/nstar(ix,iy,iz,ll) 

       !********************************************************************

       select case (trim(coll(i)%key))

          !********************************************************************

       case('OHMEGA')

          !  excitation of  electron collisons with positive ions
          call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
          cdn = c0 / g(ul) * ct  * nne / stg
          cup = cdn * nsr

          !********************************************************************

       case('CE')

          ! excitation of neutrals by electron collisions 
          call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
          cdn = g(ll)/g(ul) * ct * nne * stg
          cup = cdn * nsr

          !********************************************************************

       case('CH')

          ! collisions with neutral hydrogen)
          ! these are downward rates (cdn); cup are upward rates

          call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
          cdn = ct * nnh(1)
          cup = cdn * nsr  

          !********************************************************************

       case('CI')

          ! ionization by electron collisions
          call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
          de_k = (ev(ul)-ev(ll))*ee/bk
          cup = ct * nne * stg * exp(-de_k/tg) 
          cdn = cup / nsr

          !********************************************************************

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

          call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
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

          cdn = coll(i)%c(3) * (tg/1.0d4)**(-coll(i)%c(4)) &
               +  summrs * coll(i)%c(5) /tg/stg * exp(-coll(i)%c(7)/tg) * &  
               (1.d0 + coll(i)%c(6) * (exp(- coll(i)%c(8) /tg)))
          cdn = cdn * nne

          cup = coll(i)%c(1) * stg * exp(-coll(i)%c(2)/tg) / (1.d0 + 0.1d0 * tg / coll(i)%c(2)) 
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

             xj=coll(i)%c(k+1) * ee / bk / tg

             fac=exp(-xj)*sqrt(xj)

             fxj=fac*( coll(i)%c(k+2) + coll(i)%c(k+3) * (1.d0+xj) &
                  +(coll(i)%c(k+4) - xj*(coll(i)%c(k+2) + coll(i)%c(k+3) &
                  *(2.d0+xj)))*fone(xj)+ &
                  coll(i)%c(k+5)*xj*ftwo(xj) )

             fac=6.69d-07/coll(i)%c(k+1)/sqrt(coll(i)%c(k+1))

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
          cup = coll(i)%c(1) * ct * nne

          !********************************************************************

       case('AR85-CH')

          ! The AR85**** rates are based on 
          ! Arnaud & Rothenflug
          ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
          ! "An updated evaluation of recombination and ionization rates"

          !AR85-CH is charge transfer with neutral hydrogen atom^+ + H -> atom + H^+
          if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) )  then
             t4=tg/1.0d4
             cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * nnh(1)
          endif

          !********************************************************************  

       case('AR85-CH+')

          ! The AR85**** rates are based on 
          ! Arnaud & Rothenflug
          ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
          ! "An updated evaluation of recombination and ionization rates"

          ! ar85-ch+ is  charge transfer with ionized hydrogen atom + H^+ -> atom^+ + H
          if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) )  then
             t4=tg/1.0d4
             cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * nnh(nhlevels)
          endif

          !********************************************************************

       case('AR85-CHE')

          ! The AR85**** rates are based on 
          ! Arnaud & Rothenflug
          ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
          ! "An updated evaluation of recombination and ionization rates"

          !AR85-CHE is charge transfer with neutral helium atom^+ + He -> atom + He^+
          if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) ) then
             t4=tg/1.0d4
             call hepop(tg,sum(nnh),he1,he2)
             cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * he1 
          endif

          !********************************************************************

       case('AR85-CHP')

          ! The AR85**** rates are based on 
          ! Arnaud & Rothenflug
          ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
          ! "An updated evaluation of recombination and ionization rates"

          ! ar85-chp or ar85-che+ is  charge transfer with ionized helium atom^+ + He^+ -> atom + He^2+
          if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) ) then
             t4=tg/1.0d4
             call hepop(tg,sum(nnh),he1,he2)
             cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * he2 
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

          cup = qb * actelc * nne * coll(i)%c(1)
          cdn = cup/nsr

          !********************************************************************

       case('BADNELL')
          ! fits to dielectronic recombination from Badnell and co-workers.
          ! see for example Altun et al., 2006, A&A, 44 7, 1165
          cdn=0.0d0
          do j=1,coll(i)%nc
             cdn=cdn + coll(i)%c(j+coll(i)%nc) * exp(-coll(i)%c(j)/tg)
          end do

          summrs=sumscl*summers(id,ion,ul,nne) + (1.0d0-sumscl)
          cdn = cdn * nne * tg**(-1.5d0)* summrs

       end select

       !********************************************************************

       cij(ul,ll) = cij(ul,ll)  + cup ! upward rate Cll->ul 
       cij(ll,ul) = cij(ll,ul)  + cdn ! downward rate Cul->ll

    end do ! transition loop 

   end subroutine gencol


   subroutine gencol_new(coll,id,nlevel,g,ev,ion,cij,temp,ne,nstar,nh)
      use mconst
      use mclu
      use mem_mod
      implicit none
  
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
  
      real(8) :: t4,summrs
      real(8) :: betab,ddekt,wlog,wb,dexin1,xj,fac,fxj
      real(8) :: dekti,qb,he1,he2
      integer :: ul,ll,i,j,k,z, nhlevels
  
      character(len=6),parameter :: rname='gencol'
  
      cij=0.0d0
  
      ! get local copies for ease of notation
      tg=temp
      stg=sqrt(tg)
      nne=ne
      nnh=nh(:)
      nhlevels = size(nh, 1)
  
      ! loop over all the transitions
      do i=1, size(coll, 1)
  
         cup=0.0d0
         cdn=0.0d0
         ll=coll(i)%ll
         ul=coll(i)%ul
         nsr=nstar(ul)/nstar(ll) 
  
         !********************************************************************
  
         select case (trim(coll(i)%key))
  
            !********************************************************************
  
         case('OHMEGA')
  
            !  excitation of  electron collisons with positive ions
            call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
            cdn = c0 / g(ul) * ct  * nne / stg
            cup = cdn * nsr
  
            !********************************************************************
  
         case('CE')
  
            ! excitation of neutrals by electron collisions 
            call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
            cdn = g(ll)/g(ul) * ct * nne * stg
            cup = cdn * nsr
  
            !********************************************************************
  
         case('CH')
  
            ! collisions with neutral hydrogen)
            ! these are downward rates (cdn); cup are upward rates
  
            call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
            cdn = ct * nnh(1)
            cup = cdn * nsr  
  
            !********************************************************************
  
         case('CI')
  
            ! ionization by electron collisions
            call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
            de_k = (ev(ul)-ev(ll))*ee/bk
            cup = ct * nne * stg * exp(-de_k/tg) 
            cdn = cup / nsr
  
            !********************************************************************
  
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
  
            call hermite_interpolation(ct,tg,coll(i)%t,coll(i)%c,coll(i)%dcdt)
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
  
            cdn = coll(i)%c(3) * (tg/1.0d4)**(-coll(i)%c(4)) &
                 +  summrs * coll(i)%c(5) /tg/stg * exp(-coll(i)%c(7)/tg) * &  
                 (1.d0 + coll(i)%c(6) * (exp(- coll(i)%c(8) /tg)))
            cdn = cdn * nne
  
            cup = coll(i)%c(1) * stg * exp(-coll(i)%c(2)/tg) / (1.d0 + 0.1d0 * tg / coll(i)%c(2)) 
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
  
               xj=coll(i)%c(k+1) * ee / bk / tg
  
               fac=exp(-xj)*sqrt(xj)
  
               fxj=fac*( coll(i)%c(k+2) + coll(i)%c(k+3) * (1.d0+xj) &
                    +(coll(i)%c(k+4) - xj*(coll(i)%c(k+2) + coll(i)%c(k+3) &
                    *(2.d0+xj)))*fone(xj)+ &
                    coll(i)%c(k+5)*xj*ftwo(xj) )
  
               fac=6.69d-07/coll(i)%c(k+1)/sqrt(coll(i)%c(k+1))
  
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
            cup = coll(i)%c(1) * ct * nne
  
            !********************************************************************
  
         case('AR85-CH')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            !AR85-CH is charge transfer with neutral hydrogen atom^+ + H -> atom + H^+
            if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) )  then
               t4=tg/1.0d4
               cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * nnh(1)
            endif
  
            !********************************************************************  
  
         case('AR85-CH+')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            ! ar85-ch+ is  charge transfer with ionized hydrogen atom + H^+ -> atom^+ + H
            if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) )  then
               t4=tg/1.0d4
               cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * nnh(nhlevels)
            endif
  
            !********************************************************************
  
         case('AR85-CHE')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            !AR85-CHE is charge transfer with neutral helium atom^+ + He -> atom + He^+
            if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) ) then
               t4=tg/1.0d4
               call hepop(tg,sum(nnh),he1,he2)
               cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * he1 
            endif
  
            !********************************************************************
  
         case('AR85-CHP')
  
            ! The AR85**** rates are based on 
            ! Arnaud & Rothenflug
            ! Astronomy and Astrophysics Supplement Series (ISSN 0365-0138), vol. 60, June 1985, p. 425-457
            ! "An updated evaluation of recombination and ionization rates"
  
            ! ar85-chp or ar85-che+ is  charge transfer with ionized helium atom^+ + He^+ -> atom + He^2+
            if ( tg>=coll(i)%c(1) .and. tg<=coll(i)%c(2) ) then
               t4=tg/1.0d4
               call hepop(tg,sum(nnh),he1,he2)
               cdn = coll(i)%c(3) * 1d-9 * t4**coll(i)%c(4) * (1.0d0+coll(i)%c(5)*exp(coll(i)%c(6)*t4)) * he2 
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
  
            cup = qb * actelc * nne * coll(i)%c(1)
            cdn = cup/nsr
  
            !********************************************************************
  
         case('BADNELL')
            ! fits to dielectronic recombination from Badnell and co-workers.
            ! see for example Altun et al., 2006, A&A, 44 7, 1165
            cdn=0.0d0
            do j=1,coll(i)%nc
               cdn=cdn + coll(i)%c(j+coll(i)%nc) * exp(-coll(i)%c(j)/tg)
            end do
  
            summrs=sumscl*summers(id,ion,ul,nne) + (1.0d0-sumscl)
            cdn = cdn * nne * tg**(-1.5d0)* summrs
  
         end select
  
         !********************************************************************
  
         cij(ul,ll) = cij(ul,ll)  + cup ! upward rate Cll->ul 
         cij(ll,ul) = cij(ll,ul)  + cdn ! downward rate Cul->ll
  
      end do ! transition loop 
  
   end subroutine gencol_new
 
   !*********************************************************************
 
   logical function is_valid_process(key)
    character(len=8),intent(in) :: key
    integer :: i
 
    is_valid_process=.false.
    do i=1,size(theprocesses)
       if ( key(1:len_trim(key)) == trim(theprocesses(i)))  then
          is_valid_process=.true.
          exit
       end if
    end do
 
   end function is_valid_process
 
   ! *******************************************************************

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
      use mconst,only:pi,euler_const
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
         if(count .gt. 100.) call stop('ftwo: too many iterations')
      enddo
      ftwo=exp(x)*((log(x)+euler_const)*(log(x)+euler_const)*0.5 + f0x)
      endif
      return
   end function ftwo


   !*********************************************************************

   pure integer function atomic_number(atomic_symbol)
      use mconst,only:element_name
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

   pure subroutine rowcol(an,row,col)
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
   end subroutine rowcol

   !**********************************************************************

   subroutine ar85cea(id, ion, tg,ilo,ihi,cup)
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
         call messg('ar85-cea:  atomic number = ' //string3(z) // ' for ' // element)
         call messg('ar85-cea:  no autoionization included')
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
   end subroutine ar85cea

   !*********************************************************************

   subroutine hepop(tg,hpop,he1,he2)
      use mmath,only:linear1d
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

   end subroutine hepop

   !*********************************************************************

     
 end module coll_mod