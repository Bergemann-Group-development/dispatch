!===============================================================================
!> non-ideal MHD module
!===============================================================================
MODULE non_ideal_mod
  USE const
  USE io_mod
  USE io_unit_mod
  USE link_mod
  use hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE index_mod
  USE mesh_mod
  implicit none
  private
  type, public :: non_ideal_t
    real:: eta_Ohm, etaMD, gamma_AD, rho_ions
    real:: rho0 
    real:: crho 
    real:: T0   
    real:: cT   
    real:: B0   
    real:: cB   
    real:: E0
    real:: cE
    real:: mu
    real:: mp
    real:: kB
    real:: NA  
    real,dimension(200) :: rhoarr
    real,dimension(60)  :: Tarr
    real,dimension(150) :: Barr
    real,dimension(60) :: Earr
    real,dimension(200,60,150) :: eta_Ohm_tbl, eta_AD_tbl
    real,pointer,dimension(:,:,:,:,:) :: flux
    logical:: mhd_AD, mhd_Ohm, ntestDADM
    logical:: is_used
  contains
    procedure:: init
    procedure:: non_ideal_comp
    procedure:: flux_upd
    procedure:: non_ideal_emf_up
    procedure:: computejb2
    procedure:: computambip
    procedure:: computdifmag            
    procedure:: eta_table
    procedure:: eta_table_init
!    procedure:: comp_weight
!    procedure:: read_table
  end type
  type(non_ideal_t), public:: non_ideal
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE init (self)
  class (non_ideal_t):: self
  integer:: iostat
  real:: gamma_AD=75., rho_ions=1.0, eta_Ohm=0.0, etaMD=0.0
  logical:: mhd_AD=.false., mhd_Ohm=.false., ntestDADM=.true.
  logical, save:: first_time=.true.
  namelist /non_ideal_params/ gamma_AD, eta_Ohm, etaMD, rho_ions, mhd_AD, mhd_Ohm, ntestDADM
  !---------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, non_ideal_params,iostat=iostat)
    if (io%master) write (*, non_ideal_params)
  end if
  !$omp end critical (init_cr)
  self%mhd_Ohm  = mhd_Ohm
  self%eta_Ohm  = eta_Ohm
  self%mhd_AD   = mhd_AD
  self%gamma_AD = gamma_AD
  self%rho_ions = rho_ions
  self%etaMD    = etaMD
  self%is_used  = mhd_Ohm .or. mhd_AD
  self%ntestDADM=ntestDADM
  if ( (self%is_used .eqv. .true.) .and. (self%ntestDADM .eqv. .false.) ) then
    self%eta_Ohm_tbl = 0.d0
    self%eta_AD_tbl  = 0.d0
    self%kB  = 1.3807e-16                ! Boltzmann constant
    self%mu  = 2.37                      ! mean molecular weight
    self%mp  = 1.6737236e-24             ! Proton mass  
    self%NA  = 6.022140857e23            ! Avogadro constant
    call eta_table_init(self,self%eta_Ohm_tbl,self%eta_AD_tbl)
  endif

END SUBROUTINE init

subroutine non_ideal_comp(self,mesh,idx,bf,uin,qin,fluxambdiff,emfambdiff,u_max,ds,dt,gamma)
  class (non_ideal_t):: self
  class (mesh_t):: mesh(:)
  class (index_t):: idx
  integer:: n(4), l(3), u(3)
  real :: u_max
  real(8) :: gamma
  real(8) :: ds(3), dt
  real, dimension(:,:,:,:):: bf, uin, qin, fluxambdiff, emfambdiff
  ! Output electromotive force
  real(8), allocatable, dimension(:,:,:)::emfx
  real(8), allocatable, dimension(:,:,:)::emfy
  real(8), allocatable, dimension(:,:,:)::emfz

  ! Output courant vector in the cell
  real(8), allocatable, dimension(:,:,:,:,:)::bmagij
  real(8), allocatable, dimension(:,:,:)::jcentersquare,jxbsquare
  real(8), allocatable, dimension(:,:,:,:)::bemfx,bemfy,bemfz
  real(8), allocatable, dimension(:,:,:,:)::jemfx,jemfy,jemfz
  real(8), allocatable, dimension(:,:,:,:)::florentzx,florentzy,florentzz
  real(8), allocatable, dimension(:,:,:,:)::fluxmd,fluxh,fluxad
  real(8), allocatable, dimension(:,:,:,:)::emfohmdiss,fluxohm 

  n = shape(uin)
  allocate (emfx         (n(1),n(2),n(3)))
  allocate (emfy         (n(1),n(2),n(3)))
  allocate (emfz         (n(1),n(2),n(3)))
  allocate (bmagij       (n(1),n(2),n(3),3,3))
  allocate (jcentersquare(n(1),n(2),n(3)))
  allocate (jxbsquare    (n(1),n(2),n(3)))
  allocate (bemfx        (n(1),n(2),n(3),3))
  allocate (bemfy        (n(1),n(2),n(3),3))
  allocate (bemfz        (n(1),n(2),n(3),3))
  allocate (jemfx        (n(1),n(2),n(3),3))
  allocate (jemfy        (n(1),n(2),n(3),3))
  allocate (jemfz        (n(1),n(2),n(3),3))
  allocate (florentzx    (n(1),n(2),n(3),3))
  allocate (florentzy    (n(1),n(2),n(3),3))
  allocate (florentzz    (n(1),n(2),n(3),3))
  allocate (fluxmd       (n(1),n(2),n(3),3))
  allocate (fluxh        (n(1),n(2),n(3),3))
  allocate (fluxad       (n(1),n(2),n(3),3))
  allocate (emfohmdiss   (n(1),n(2),n(3),3))
  allocate (fluxohm      (n(1),n(2),n(3),3))

  emfx=0.d0
  emfy=0.d0
  emfz=0.d0
  bmagij=0.d0
  jcentersquare=0.d0
  jxbsquare=0d0
  bemfx=0.d0
  bemfy=0.d0
  bemfz=0.d0
  jemfx=0.d0
  jemfy=0.d0
  jemfz=0.d0
  florentzx=0.d0
  florentzy=0.d0
  florentzz=0.d0
  fluxmd=0.d0
  fluxh=0.d0
  fluxad=0.d0
  emfohmdiss=0.d0
  fluxohm=0.d0


  if((self%mhd_AD).or.(self%mhd_Ohm)) then
  !  compute Lorentz Force with current
     call computejb2(self,bf,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,ds)
  endif
  
  
  ! AMBIPOLAR DIFFUSION
  
  if(self%mhd_AD) then
     call computambip(self,uin,qin,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,ds,dt,gamma)
     l = mesh%lo
     u = mesh%uo
     u_max = max(u_max, maxval( &
       uin(l(1):u(1),l(2):u(2),l(3):u(3),idx%bx)**2 + &
       uin(l(1):u(1),l(2):u(2),l(3):u(3),idx%by)**2 + &
       uin(l(1):u(1),l(2):u(2),l(3):u(3),idx%bz)**2 / &
       (self%gamma_AD * self%rho_ions * uin(:,:,:,1) * minval(ds)) )  )
  endif
  
  
  ! OHMIC DISSIPATION
  
  if(self%mhd_Ohm) then
     call computdifmag(self,uin,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,ds,dt)
     u_max = max(u_max,self%eta_Ohm/ minval(ds))
  endif

  
  deallocate(emfx,emfy,emfz,bmagij,jcentersquare,jxbsquare,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,emfohmdiss,fluxohm)

end subroutine

subroutine flux_upd(self,flux,fluxambdiff, dx, dt, idim)
   class (non_ideal_t):: self

   integer:: n(5)
   integer :: i, j, k, idim
   real, dimension(:,:,:,:,:):: flux 
   real, dimension(:,:,:,:)::fluxambdiff

   real, allocatable, dimension(:,:,:,:)::fluxohm
   real(8) :: dx, dt


   n = shape(flux)
   allocate (fluxohm(n(1),n(2),n(3),3))

   fluxohm=0d0
   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then
      flux(:,:,:,5,idim) = flux(:,:,:,5,idim) + (fluxambdiff(:,:,:,idim)+fluxohm(:,:,:,idim) ) * dt / dx
   end if

   deallocate(fluxohm)

end subroutine


subroutine non_ideal_emf_up(self,emfup,emfambdiff, dx, dt, idim)
   class (non_ideal_t):: self
   integer:: n(4)
   integer :: idim
   real, dimension(:,:,:)::emfup
   real, dimension(:,:,:,:)::emfambdiff
   real, allocatable, dimension(:,:,:,:)::emfohmdiss
   real(8) :: dx, dt

   n = shape(emfambdiff)
   allocate (emfohmdiss(n(1),n(2),n(3),n(4)))
   emfohmdiss = 0d0

   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then

#if NDIM==1

   emfup(:,:,:) = ( emfambdiff(:,:,:,idim) + emfohmdiss(:,:,:,idim) ) * dt / dx
#else
   emfup(2:n(1)-1,2:n(2)-1,2:n(3)-1) = emfup(2:n(1)-1,2:n(2)-1,2:n(3)-1) + &
        ( emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,idim) + emfohmdiss(2:n(1)-1,2:n(2)-1,2:n(3)-1,idim) ) * dt / dx
#endif     
   end if

   deallocate(emfohmdiss)

end subroutine


subroutine computejb2(self,bf,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,ds)
  class (non_ideal_t):: self

  integer :: i, j, k,  ibegin, iend, jbegin, jend, kbegin, kend
  real, dimension(:,:,:,:)::bf 
  real, dimension(:,:,:,:)::q 


  real(8), dimension(:,:,:,:)::bemfx,bemfy,bemfz
  real(8), dimension(:,:,:,:)::jemfx,jemfy,jemfz
  real(8), dimension(:,:,:,:,:)::bmagij
  real(8), dimension(:,:,:,:)::florentzx,florentzy,florentzz
  real(8), dimension(:,:,:,:)::fluxmd,fluxh,fluxad
!  real,dimension(1:self%gn(1),1:self%gn(2),1:self%gn(3),1:3)::jcell


  ! declare local variables
  integer :: m

  integer:: n(4)
  real(8), allocatable, dimension(:,:,:,:)::bmagijbis
  real(8), allocatable, dimension(:,:,:,:,:)::jface
  real(8), allocatable, dimension(:,:,:,:)::bcenter
  real(8), allocatable, dimension(:,:,:,:,:)::fluxbis,fluxter,fluxquat
  real(8)::ds(3),dx,dy,dz

  n = shape(bemfx)
  allocate (bmagijbis(n(1),n(2),n(3),n(4)))
  allocate (jface    (n(1),n(2),n(3),n(4),3))
  allocate (bcenter  (n(1),n(2),n(3),n(4)))
  allocate (fluxbis  (n(1),n(2),n(3),n(4),3))
  allocate (fluxter  (n(1),n(2),n(3),n(4),3))
  allocate (fluxquat (n(1),n(2),n(3),n(4),3))

  bmagijbis = 0.d0
  jface     = 0.d0
  bcenter   = 0.d0
  fluxbis   = 0.d0
  fluxter   = 0.d0
  fluxquat  = 0.d0

  dx = ds(1) ; dy = ds(2) ; dz = ds(3)
   ! Compute before the loop because bmagijbis has to be computed with index i+1/j+1/k+1 in loop over i/j/k
   bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)= 0.25d0* (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) + bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) + bf(2:n(1)-1,1:n(2)-2,1:n(3)-2,1) )
   bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= 0.25d0* (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) + bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,2) + bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,2) + bf(1:n(1)-2,2:n(2)-1,1:n(3)-2,2) ) 
   bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= 0.25d0* (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) + bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,3) + bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,3) + bf(1:n(1)-2,1:n(2)-2,2:n(3)-1,3) ) 
   ! Define the 1 and 2 components of bcenter, since they are needed further below, when computing the current on the cell interfaces (jface) 
   ! magnetic field at center of cells
   bcenter(1:n(1),1:n(2),1:n(3),1:3)=q(1:n(1),1:n(2),1:n(3),6:8)

   !!!!!!!!!!!!!!!!!!
   ! EMF x
   !!!!!!!!!!!!!!!!!!
   bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=0.25d0* ( q (2:n(1)-1,2:n(2)-1,2:n(3)-1,6) + q (2:n(1)-1,1:n(2)-2,2:n(3)-1,6) + q(2:n(1)-1,2:n(2)-1,1:n(3)-2,6) + q(2:n(1)-1,1:n(2)-2,1:n(3)-2,6) )
   bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=0.5d0 * ( bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) + bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,2) )
   bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=0.5d0 * ( bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) + bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,3) )
   !!!!!!!!!!!!!!!!!!
   ! EMF y
   !!!!!!!!!!!!!!!!!!
   bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=0.5d0 * ( bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) )
   bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=0.25d0* ( q (2:n(1)-1,2:n(2)-1,2:n(3)-1,7) + q (1:n(1)-2,2:n(2)-1,2:n(3)-1,7) + q(2:n(1)-1,2:n(2)-1,1:n(3)-2,7) + q(1:n(1)-2,2:n(2)-1,1:n(3)-2,7) )
   bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=0.5d0 * ( bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,3) + bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) )
   !!!!!!!!!!!!!!!!!!
   ! EMF z
   !!!!!!!!!!!!!!!!!!
   bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=0.5d0 * ( bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) )
   bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=0.5d0 * ( bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) + bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,2) )
   bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=0.25d0* ( q (2:n(1)-1,2:n(2)-1,2:n(3)-1,8) + q (1:n(1)-2,2:n(2)-1,2:n(3)-1,8) + q(2:n(1)-1,1:n(2)-2,2:n(3)-1,8) + q(1:n(1)-2,1:n(2)-2,2:n(3)-1,8) )
   do m=1,3
   ! bmagij is the value of the magnetic field Bi where Bj 
   ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
   ! and we can write it Bx,y  
   !! m+5 mandatory cf Bx=uin(i,j,k,6)
      bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,m,m)=bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,m)
   end do
   ! case Bx,y
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,6) + q(2:n(1)-1,1:n(2)-2,2:n(3)-1,6) )
   ! case Bx,z
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,6) + q(2:n(1)-1,2:n(2)-1,1:n(3)-2,6) )
   ! case By,y
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,7) + q(1:n(1)-2,2:n(2)-1,2:n(3)-1,7) )
   ! case By,z
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,3) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,7) + q(2:n(1)-1,2:n(2)-1,1:n(3)-2,7) )
   ! case Bz,x
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,8) + q(1:n(1)-2,2:n(2)-1,2:n(3)-1,8) )
   ! case Bz,y
   bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) = 0.5d0 * (q(2:n(1)-1,2:n(2)-1,2:n(3)-1,8) + q(2:n(1)-1,1:n(2)-2,2:n(3)-1,8) )
   jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) = (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,3) ) / dy - (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,2) ) / dz 
   jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) = (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2) - bmagij(2:n(1)-1,2:n(2)-1,1:n(3)-2,1,2) ) / dz - (bmagijbis(1:n(1)-2,2:n(2)-1,2:n(3)-1,3) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dx
   jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) = (bmagijbis(1:n(1)-2,2:n(2)-1,2:n(3)-1,2) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) ) / dx - (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3)- bmagij(2:n(1)-1,1:n(2)-2,2:n(3)-1,1,3) ) / dy
   jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) = (bmagijbis(2:n(1)-1,3:n(2)  ,2:n(3)-1,3) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dy - (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1)- bmagij(2:n(1)-1,2:n(2)-1,1:n(3)-2,2,1) ) / dz
   jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) = (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - bf(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) ) / dz - (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,3) ) / dx
   jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) = (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,3) - bmagij(1:n(1)-2,2:n(2)-1,2:n(3)-1,2,3) ) / dx - (bmagijbis(2:n(1)-1,3:n(2)  ,2:n(3)-1,1) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dy
   jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) = (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1) - bmagij(2:n(1)-1,1:n(2)-2,2:n(3)-1,3,1) ) / dy - (bmagijbis(2:n(1)-1,2:n(2)-1,3:n(3)  ,2) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) ) / dz
   jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) = (bmagijbis(2:n(1)-1,2:n(2)-1,3:n(3)  ,1) - bmagijbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dz - (bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) - bmagij(1:n(1)-2,2:n(2)-1,2:n(3)-1,3,2) ) / dx
   jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) = (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - bf(1:n(1)-2,2:n(2)-1,2:n(3)-1,2) ) / dx - (bf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - bf(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) ) / dy


if(self%mhd_AD) then

! EMF x
   florentzx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
   florentzx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
   florentzx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
   florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
   florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
   florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
   florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
   florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
   florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
  

endif


! computation of current on faces

if((self%mhd_AD).or.(self%mhd_Ohm)) then

! face at i-1/2,j,k

   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,1)= (bemfz  (2:n(1)-1,3:n(2)  ,2:n(3)-1,3) - bemfz  (2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dy - (bemfy  (2:n(1)-1,2:n(2)-1,3:n(3)  ,2) - bemfy  (2:n(1)-1,2:n(2)-1,2:n(3)-1,2) ) / dz
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1)= (bemfy  (2:n(1)-1,2:n(2)-1,3:n(3)  ,1) - bemfy  (2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dz - (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - bcenter(1:n(1)-2,2:n(2)-1,2:n(3)-1,3) ) / dx
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1)= (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - bcenter(1:n(1)-2,2:n(2)-1,2:n(3)-1,2) ) / dx - (bemfz  (2:n(1)-1,3:n(2)  ,2:n(3)-1,1) - bemfz  (2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dy
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2)= (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - bcenter(2:n(1)-1,2:n(2)-1,1:n(3)-2,3) ) / dy - (bemfx  (2:n(1)-1,2:n(2)-1,3:n(3)  ,2) - bemfx  (2:n(1)-1,2:n(2)-1,2:n(3)-1,2) ) / dz
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,2)= (bemfx  (2:n(1)-1,2:n(2)-1,3:n(3)  ,1) - bemfx  (2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dz - (bemfz  (3:n(1)  ,2:n(2)-1,2:n(3)-1,3) - bemfz  (2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dx
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2)= (bemfz  (3:n(1)  ,2:n(2)-1,2:n(3)-1,2) - bemfz  (2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dx - (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - bcenter(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) ) / dy
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3)= (bemfx  (2:n(1)-1,3:n(2)  ,2:n(3)-1,3) - bemfx  (2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dy - (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - bcenter(2:n(1)-1,2:n(2)-1,1:n(3)-2,2) ) / dz             
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,3)= (bcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - bcenter(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) ) / dz - (bemfy  (3:n(1)  ,2:n(2)-1,2:n(3)-1,3) - bemfy  (2:n(1)-1,2:n(2)-1,2:n(3)-1,3) ) / dx             
   jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,3)= (bemfy  (3:n(1)  ,2:n(2)-1,2:n(3)-1,2) - bemfy  (2:n(1)-1,2:n(2)-1,2:n(3)-1,2) ) / dx - (bemfx  (2:n(1)-1,3:n(2)  ,2:n(3)-1,1) - bemfx  (2:n(1)-1,2:n(2)-1,2:n(3)-1,1) ) / dy            

   fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:) = jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:) - jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)
   fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:) = jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:) - jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)
   fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:) = jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:) - jface(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)

   fluxmd(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,1)
   fluxmd(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,2)
   fluxmd(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,3)


endif

if(self%mhd_AD) then

   fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:) - fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)
   fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:) - fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)
   fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)=fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:) - fluxbis(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)

   fluxh(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,1)
   fluxh(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,2)
   fluxh(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,3)

   fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:) - fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)
   fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:) - fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)
   fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,:)=fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:) - fluxter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,:)*bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,:)

   fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,1)
   fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,2)
   fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=fluxquat(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,3)


endif

  deallocate (bmagijbis, jface, bcenter, fluxbis, fluxter, fluxquat)



end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(self,u,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,ds,dt)
  class (non_ideal_t):: self

  integer :: i, j, k

  integer:: n(4)
  real, dimension(:,:,:,:)::u 
  real, dimension(:,:,:,:)::q 

  real(8), dimension(:,:,:,:)::bemfx,bemfy,bemfz
  real(8), dimension(:,:,:,:)::jemfx,jemfy,jemfz
  real(8), dimension(:,:,:,:,:)::bmagij
  real(8), dimension(:,:,:,:)::fluxmd

  real(8), dimension(:,:,:,:):: emfohmdiss,fluxohm 
  real(8), dimension(:,:,:)::jcentersquare

  ! declare local variables
    integer :: m, h


  ! WARNING following quantities defined with three components even
  ! if ndim<3 !
  real(8), allocatable, dimension(:,:,:,:)::jcenter
  real(8), allocatable, dimension(:,:,:,:,:)::jface
  real(8), allocatable, dimension(:,:,:,:)::jemf

  real(8)::rhocell,bcell,Cv

  real(8), allocatable, dimension(:,:,:)::bsquarex,bsquarey,bsquarez,dtlim
  real(8), allocatable, dimension(:,:,:)::rhox,rhoy,rhoz,epsx,epsy,epsz
  real(8), allocatable, dimension(:,:,:)::rhof,pf,bsqf,epsf

  integer, dimension(3) :: index_i,index_j,index_k

  real(8):: ds(3)
  real(8):: dx, dy, dz, dt

  n = shape(u)
  allocate (jcenter (n(1),n(2),n(3),3))
  allocate (jface   (n(1),n(2),n(3),3,3))
  allocate (jemf    (n(1),n(2),n(3),3))
  allocate (bsquarex(n(1)-2,n(2)-2,n(3)-2) )
  allocate (bsquarey(n(1)-2,n(2)-2,n(3)-2) )
  allocate (bsquarez(n(1)-2,n(2)-2,n(3)-2) )
  allocate (dtlim   (n(1)-2,n(2)-2,n(3)-2) )
  allocate (rhox    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (rhoy    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (rhoz    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (epsx    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (epsy    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (epsz    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (rhof    (n(1)-2,n(2)-2,n(3)-2) )
  allocate (pf      (n(1)-2,n(2)-2,n(3)-2) )
  allocate (bsqf    (n(1)-1,n(2)-1,n(3)-1) )
  allocate (epsf    (n(1)-1,n(2)-1,n(3)-1) )

  jcenter  = 0.d0
  jface    = 0.d0
  jemf     = 0.d0
  bsquarex = 0.d0
  bsquarey = 0.d0
  bsquarez = 0.d0
  dtlim    = 0.d0
  rhox     = 0.d0
  rhoy     = 0.d0
  rhoz     = 0.d0
  epsx     = 0.d0
  epsy     = 0.d0
  epsz     = 0.d0
  rhof     = 0.d0
  pf       = 0.d0
  bsqf     = 0.d0
  epsf     = 0.d0

  index_i = (/1,0,0/)
  index_j = (/0,1,0/)
  index_k = (/0,0,1/)

  dx = ds(1); dy = ds(2) ; dz = ds(3)
  dtlim = dt !neil

  jemf(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=jemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
  jemf(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=jemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
  jemf(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=jemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
  rhox=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) + u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) + u(2:n(1)-1,1:n(2)-2,1:n(3)-2,1))
  rhoy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) + u(1:n(1)-2,2:n(2)-1,1:n(3)-2,1))
  rhoz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) + u(1:n(1)-2,1:n(2)-2,2:n(3)-1,1))
  !epsx=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,j-1,2:ng-1,nvar) + u(2:ng-1,2:ng-1,1:ng-2,nvar) + u(2:ng-1,1:ng-2,1:ng-2,nvar))
  !epsy=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(i-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,2:ng-1,1:ng-2,nvar) + u(1:ng-2,2:ng-1,1:ng-2,nvar))
  !epsz=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(i-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,1:ng-2,2:ng-1,nvar) + u(1:ng-2,1:ng-2,2:ng-1,nvar))
   if(self%mhd_Ohm)then
     bsquarex=bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
     bsquarey=bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
     bsquarez=bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
   end if
              
  ! WARNING dB/dt=-curl(eta*J)
  emfohmdiss(2:n(1)-1,2:n(2)-1,2:n(3)-1,:)=-self%etaMD*jemf(2:n(1)-1,2:n(2)-1,2:n(3)-1,:)
  !!!!!!!!!!!!!!!!!!!!!!!
  !
  ! compute j at center of cells
  !
  ! mandatory for non isotherm case
  
  ! bmagij is the value of the magnetic field Bi where Bj 
  ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
  ! and we can write it Bx,y

  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)= (bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,3,2) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) ) / dz - (bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,2,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,3) ) / dy
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= (bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,1,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3) ) / dx - (bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,3,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1) ) / dz
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= (bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,2,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1) ) / dy - (bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,2,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1) ) / dx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  do h = 1,3
     
     rhof=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(2-index_i(h):n(1)-index_i(h),2-index_j(h):n(2)-index_j(h),2-index_k(h):n(3)-index_k(h),1))
     !epsf=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar)+u(2-index_i(h):ng-index_i(h),2-index_j(h):ng-index_j(h),2-index_k(h):ng-index_k(h),nvar))
     bsqf=bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,h)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,h)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,h)**2
          
     fluxohm(2:n(1)-1,2:n(2)-1,2:n(3)-1,h)=self%etaMD*fluxmd(2:n(1)-1,2:n(2)-1,2:n(3)-1,h)
     
  enddo

  deallocate(jcenter, jface, jemf, bsquarex, bsquarey, bsquarez, dtlim, rhox, rhoy, rhoz, epsx, epsy, epsz, rhof, pf, bsqf, epsf )
  
end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(self,u,q,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,ds,dt,gamma)
  class (non_ideal_t):: self
  
  real, dimension(:,:,:,:)::u 
  real, dimension(:,:,:,:)::q 
  
  real(8), dimension(:,:,:,:)::bemfx,bemfy,bemfz
  real(8), dimension(:,:,:,:)::florentzx,florentzy,florentzz
  real(8), dimension(:,:,:,:)::fluxad
  real(8), dimension(:,:,:,:,:)::bmagij
  real, dimension(:,:,:,:)::emfambdiff
  real, dimension(:,:,:,:)::fluxambdiff
  real(8), dimension(:,:,:)::jxbsquare

! declare local variables
  integer :: m, ivar
  integer:: n(4)

  real(8), allocatable, dimension(:,:,:)::rhocellmin
  real(8), allocatable, dimension(:,:,:)::rhofx,rhofy,rhofz
  real(8), allocatable, dimension(:,:,:)::bsquarex,bsquarey,bsquarez,bsquare
  real(8), allocatable, dimension(:,:,:)::bsquarexx,bsquareyy,bsquarezz
  real(8), allocatable, dimension(:,:,:)::betaad2,betaad,eta_AD,eta_Ohm
  real(8), allocatable, dimension(:,:,:)::rhox,rhoy,rhoz,rhocell,bcell,bcellold
  real(8)::dtlim,Cv

  real(8), allocatable, dimension(:,:,:,:)::florentz
  real(8), allocatable, dimension(:,:,:)  ::bsquaremax

  real(8), allocatable, dimension(:,:,:,:)::jcenter
  real(8), allocatable, dimension(:,:,:,:)::jxb

  real(8), allocatable, dimension(:,:,:)::uxfx,uxfy,uxfz,uyfx,uyfy,uyfz,uzfx,uzfy,uzfz,uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,Brms2,Eth,Brms,rho,B,P,T
  real(8), allocatable, dimension(:,:,:)::etfx,etfy,etfz,etx,ety,etz

  real(8) :: ds(3)
  real(8) :: dx, dy, dz, dt, c2_4pi, mp, kB, Bcode_to_Gauss, pi, scale_v, scale_d, scale_l, scale_t, scale_T2
  real(8) :: gamma 

  n = shape(u)

  allocate (rhocellmin(n(1),n(2),n(3)))
  allocate (rhofx     (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhofy     (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhofz     (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquarex  (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquarey  (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquarez  (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquarexx (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquareyy (n(1)-2,n(2)-2,n(3)-2))
  allocate (bsquarezz (n(1)-2,n(2)-2,n(3)-2))
  allocate (betaad    (n(1)-2,n(2)-2,n(3)-2))
  allocate (betaad2   (n(1)-2,n(2)-2,n(3)-2))
  allocate (eta_AD    (n(1)-2,n(2)-2,n(3)-2))
  allocate (eta_Ohm   (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhox      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhoy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhoz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhocell   (n(1)-2,n(2)-2,n(3)-2))
  allocate (bcell     (n(1)-2,n(2)-2,n(3)-2))
  allocate (bcellold  (n(1)-2,n(2)-2,n(3)-2))
  allocate (florentz  (n(1),n(2),n(3),3))
  allocate (bsquaremax(n(1),n(2),n(3)))
  allocate (jcenter   (n(1),n(2),n(3),3))
  allocate (jxb       (n(1),n(2),n(3),3))
  allocate (uxfx      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uxfy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uxfz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyfx      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyfy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyfz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzfx      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzfy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzfz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (uxx       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uxy       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uxz       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyx       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyy       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uyz       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzx       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzy       (n(1)-2,n(2)-2,n(3)-2))
  allocate (uzz       (n(1)-2,n(2)-2,n(3)-2))
  allocate (etfx      (n(1)-2,n(2)-2,n(3)-2))
  allocate (etfy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (etfz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (etx       (n(1)-2,n(2)-2,n(3)-2))
  allocate (ety       (n(1)-2,n(2)-2,n(3)-2))
  allocate (etz       (n(1)-2,n(2)-2,n(3)-2))
  allocate (Brms2     (n(1)-2,n(2)-2,n(3)-2))
  allocate (Eth       (n(1)-2,n(2)-2,n(3)-2))
  allocate (Brms      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rho       (n(1)-2,n(2)-2,n(3)-2))
  allocate (B         (n(1)-2,n(2)-2,n(3)-2))
  allocate (P         (n(1)-2,n(2)-2,n(3)-2))
  allocate (T         (n(1)-2,n(2)-2,n(3)-2))

  rhocellmin = 0.d0
  rhofx      = 0.d0
  rhofy      = 0.d0
  rhofz      = 0.d0
  bsquarex   = 0.d0
  bsquarey   = 0.d0
  bsquarez   = 0.d0
  bsquarexx  = 0.d0
  bsquareyy  = 0.d0
  bsquarezz  = 0.d0
  betaad     = 0.d0
  betaad2    = 0.d0
  eta_AD     = 0.d0
  eta_Ohm    = 0.d0
  rhox       = 0.d0
  rhoy       = 0.d0
  rhoz       = 0.d0
  rhocell    = 0.d0
  bcell      = 0.d0
  bcellold   = 0.d0
  florentz   = 0.d0
  bsquaremax = 0.d0
  jcenter    = 0.d0
  jxb        = 0.d0

  dx = ds(1) ; dy = ds(2) ; dz = ds(3)

  if ( self%ntestDADM .eqv. .false. ) then
    scale_v = 1e6
    scale_d = 1e-17
    scale_l = 5e15
    scale_t = scale_l/scale_v
    scale_T2 = (scale_v)**2 * self%mp / self%kB
    pi = 4.*atan(1.d0)
    c2_4pi = 7.152426325648638e+19
    Bcode_to_Gauss = (4.d0*pi*scale_d)**0.5*scale_v
  endif

  !dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)
  dtlim=dt!*coefalfven

  jxb=0.0d0
  jxbsquare=0.0d0
  jcenter=0.0d0

  rhox=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) + u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) + u(2:n(1)-1,1:n(2)-2,1:n(3)-2,1))
  rhoy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1) + u(1:n(1)-2,2:n(2)-1,1:n(3)-2,1))
  rhoz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) + u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1) + u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1) + u(1:n(1)-2,1:n(2)-2,2:n(3)-1,1))

  rhofx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1))
  rhofy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1))
  rhofz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1))
              
  rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)=min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)
  rhocell = u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
 
  bsquarex=bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
  bsquarey=bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
  bsquarez=bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2+bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2+bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2

  bsquarexx=bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,1)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1)**2
  bsquareyy=bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,2)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2)**2
  bsquarezz=bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,3)**2+bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,3)**2

  bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)=max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)
                 
  ! EMF x
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)= florentzx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) - florentzx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
  rhox=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1)+u(2:n(1)-1,1:n(2)-2,1:n(3)-2,1))
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  bcellold=bcell
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)

  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhox) 
  else
    ! Call resistivity table dependent on density, thermal energy and magnetic field strength
    ! Eth = Etot - kinetic energy - magnetic energy
    uxx=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,2)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,2)+u(2:n(1)-1,1:n(2)-2,1:n(3)-2,2))
    uyx=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,3)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,3)+u(2:n(1)-1,1:n(2)-2,1:n(3)-2,3))
    uzx=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,4)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,4)+u(2:n(1)-1,1:n(2)-2,1:n(3)-2,4))
    etx=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,5)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,5)+u(2:n(1)-1,1:n(2)-2,1:n(3)-2,5))
    Brms2 = bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2 + bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2 + bemfx(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2 
    Eth = etx - 0.5 / rhox * ( uxx**2 + uyx**2 + uzx**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5      
    rho  = rhox * scale_d   
    P    = Eth * (gamma-1.) 
    T    = P / rhox * scale_T2 * self%mu
    B    = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time

    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2  
  endif

  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*betaad2 

  ! EMF y
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
  rhoy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,1))
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  bcellold=bcell
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  ! comparison with hydro+idealMHD 

  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoy) 
  else
    ! Call resistivity table dependent on density, thermal energy and magnetic field strength
    ! Eth = Etot - kinetic energy - magnetic energy
    uxy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,2)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,2)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,2))
    uyy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,3)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,3)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,3))
    uzy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,4)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,4)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,4))
    ety=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,5)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,5)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,5))
    Brms2 = bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2 + bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2 + bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2
    Eth = ety - 0.5 / rhoy * ( uxy**2 + uyy**2 + uzy**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rho  = rhoy*scale_d 
    P    = Eth * (gamma-1.) 
    T    = P / rhoy * scale_T2 * self%mu
    B    = Brms * Bcode_to_Gauss   
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2
  endif
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*betaad2            
                    
  ! EMF z
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
  rhoz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,1))
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  bcellold=bcell
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)     

  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoz) 
  else
    ! Call resistivity table dependent on density, thermal energy and magnetic field strength
    ! Eth = Etot - kinetic energy - magnetic energy
    uxz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,2)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,2)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,2))
    uyz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,3)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,3)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,3))
    uzz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,4)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,4)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,4))
    etz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,5)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,5)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,5))
    Brms2 = bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)**2 + bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2 + bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2 
    Eth = etz - 0.5 / rhoz * ( uxz**2 + uyz**2 + uzz**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rho  = rhoz*scale_d 
    P    = Eth * (gamma-1.) 
    T    = P / rhoz * scale_T2 * self%mu
    B    = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2
  endif

  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*betaad2

  ! energy flux on faces
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1) 
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  rhofx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1))
  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofx) 
  else
    uxfx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,2))
    uyfx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,3))
    uzfx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,4))
    etfx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,5))
    Brms2 = bsquarex
    Eth = etfx - 0.5 / rhofx * ( uxfx**2 + uyfx**2 + uzfx**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rho  = rhofx*scale_d 
    P    = Eth * (gamma-1.) 
    T    = P / rhofx * scale_T2 *self%mu
    B    = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2  
  endif
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
 
  rhofy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1))
  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofy) 
  else
    uxfy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,2))
    uyfy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,3))
    uzfy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,4))
    etfy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,5))
    Brms2 = bsquarey
    Eth = etfy - 0.5 / rhofy * ( uxfy**2 + uyfy**2 + uzfy**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rho  = rhofy*scale_d 
    P    = Eth * (gamma-1.) 
    T    = P / rhofy * scale_T2 * self%mu
    B    = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2 
  endif
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
 
  rhofz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1))
  if (self%ntestDADM) then
    betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofz) 
  else
    uxfz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,2))
    uyfz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,3))
    uzfz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,4))
    etfz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,5))
    Brms2 = bsquarez
    Eth = etfz - 0.5 / rhofz * ( uxfz**2 + uyfz**2 + uzfz**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rho  = rhofz*scale_d 
    P    = Eth * (gamma-1.) 
    T    = P / rhofz * scale_T2 * self%mu
    B    = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad2 = eta_AD / Brms2 
  endif
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
 
  bcellold=bcell

  
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)= ( bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,3,2) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) ) / dz - ( bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,2,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) ) / dy
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= ( bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,1,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3) ) / dx - ( bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,3,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1) ) / dz
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= ( bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,2,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1) ) / dy - ( bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,1,2) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2) ) / dx
                
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,8) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,7)
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,6) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,8)
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,7) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,6)
 
  if (self%ntestDADM) then
    betaad=1.d0/(self%gamma_AD*self%rho_ions*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)) 
  else
    Brms2 = u(2:n(1)-1,2:n(2)-1,2:n(3)-1,6)**2 + u(2:n(1)-1,2:n(2)-1,2:n(3)-1,7)**2 + u(2:n(1)-1,2:n(2)-1,2:n(3)-1,8)**2 
    Eth = u(2:n(1)-1,2:n(2)-1,2:n(3)-1,5) - 0.5 / u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) * ( u(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)**2 + u(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)**2 + u(2:n(1)-1,2:n(2)-1,2:n(3)-1,4)**2 ) - &
          0.5 * Brms2
    Brms = Brms2**0.5
    rhocell = u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
    ! Convert variables to units used in the table
    ! rho in cm^(-3) ; T in K ; B in Gauss
    rho = rhocell * scale_d
    P   = Eth * (gamma-1.) 
    T   = P / rhocell * scale_T2 * self%mu
    B   = Brms * Bcode_to_Gauss
    call eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    ! Convert back to code units. Reminder: eta_Masson = c^2/(4*pi) * eta_Marchand
    ! from seconds to code time
    eta_AD = eta_AD * c2_4pi            ! account for necessary prefactor 
    eta_AD = eta_AD / scale_v**2 / scale_t ! transfer back to code units
    betaad = eta_AD / Brms2
  endif
  jxbsquare(2:n(1)-1,2:n(2)-1,2:n(3)-1)=(jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) +&
                                 & jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) +&
                                 & jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3))*&
                                 & betaad*dtlim

  deallocate (rhocellmin, rhofx, rhofy, rhofz, bsquarex, bsquarey, bsquarez, bsquarexx, bsquareyy, bsquarezz,  &
              betaad, betaad2, eta_AD, eta_Ohm, rhox, rhoy, rhoz, rhocell, bcell, bcellold, jcenter, jxb,  &
              uxfx, uxfy, uxfz, uyfx, uyfy, uyfz, uzfx, uzfy, uzfz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, &
              etfx, etfy, etfz, etx, ety, etz, Brms2, Eth, Brms, rho, B, P, T   )


end SUBROUTINE computambip
! fin modif nimhd


subroutine eta_table(self,eta_AD,eta_Ohm,rho,T,B)
    class (non_ideal_t):: self
    real(8), dimension(:,:,:) :: rho, T, B
    real(8), dimension(:,:,:) :: eta_AD, eta_Ohm
    real(8) :: rho0, T0, B0
    real(8) :: crho, cT, cB
    integer(8) :: rhoidx, Tidx, Bidx
    real :: w1(2), w2(2),w3(2)
    integer :: n(3),i,j,k,iw,jw,kw,irho,iT,iB

    n =shape(rho)

    crho = self%crho ; cT = self%cT ; cB = self%cB    

    eta_Ohm = 0.d0 ; eta_AD = 0.d0
    do i=1,n(1)
      do j=1,n(2)
        do k=1,n(3)
          call comp_weight_element(rho(i,j,k),self%rhoarr               ,crho               ,rhoidx,w1)
          call comp_weight_element(T(i,j,k)  ,self%Tarr                 ,cT                 ,Tidx  ,w2)
          call comp_weight_element(B(i,j,k)  ,self%Barr                 ,cB                 ,Bidx  ,w3)
          do iw=1,2
            do jw=1,2
              do kw=1,2
                if (self%mhd_Ohm) then
                  eta_Ohm(i,j,k) = eta_Ohm(i,j,k) + w1(iw)*w2(jw)*w3(kw) * self%eta_Ohm_tbl(rhoidx+iw, Tidx+jw, Bidx+kw)
                end if
                if (self%mhd_AD) then
                  eta_AD(i,j,k)  = eta_AD(i,j,k)  + w1(iw)*w2(jw)*w3(kw) * self%eta_AD_tbl (rhoidx+iw, Tidx+jw, Bidx+kw)
                end if
              end do
            end do
          end do
        end do
      end do
    end do
end subroutine eta_table

subroutine comp_weight_element(x,xarr,cx,idx,w)
    implicit none
    real(8)   :: x
    real, dimension(:) :: w
    real(8) :: x0,cx
    real, dimension(:) :: xarr
    real :: arr, xd, xu, idx_r
    real(8) :: mult
    integer(8) :: idx
    integer :: nxarr

    nxarr = size(xarr)
    idx = 0
    mult = 0.d0
    xd = 0.d0
    xu = 0.d0
    x0 = xarr(1)
    mult = x / x0 
    ! log rule: log_b(x) = log_d(x)/log_d(b)
    idx  = int( log10(mult) / log10(cx) )
    idx_r = real(idx)
 
    xd = x0 * cx**idx_r
    xu = x0 * cx**(idx_r+1.d0)
    ! weight factors
    w(1)  = ( x - xd ) / ( xu - xd ) 
    w(2)  = 1. - w(1)

end subroutine comp_weight_element

subroutine comp_weight(x,xarr,cx,idx,w)
    implicit none
    real(8), dimension(:,:,:)   :: x
    real, dimension(:,:,:,:) :: w
    real(8) :: x0,cx
    real, dimension(:) :: xarr
    real, allocatable, dimension(:,:,:) :: arr, xd, xu, idx_r
    real(8), allocatable, dimension(:,:,:) :: mult
    integer(8), dimension(:,:,:) :: idx
    integer :: n(3) 
    integer :: nxarr

    n = shape(x) 
    nxarr = size(xarr)
    allocate ( mult    (n(1), n(2), n(3) ) )
    allocate ( arr     (n(1), n(2), n(3) ) )
    allocate ( xd      (n(1), n(2), n(3) ) )
    allocate ( xu      (n(1), n(2), n(3) ) )
    allocate ( idx_r   (n(1), n(2), n(3) ) )
    arr = 1.d0
    idx = 0
    mult = 0.d0
    xd = 0.d0
    xu = 0.d0
    x0 = xarr(1)
    mult = x / (x0 * arr)
    ! log rule: log_b(x) = log_d(x)/log_d(b)
    idx  = int( log10(mult) / log10(cx*arr) )
    idx_r = real(idx)
 
    xd = x0*arr * (cx*arr)**idx_r
    xu = x0*arr * (cx*arr)**(idx_r+1.d0)
    ! weight factors
    w(:,:,:,1)  = ( x - xd ) / ( xu - xd ) 
    w(:,:,:,2)  = arr - w(:,:,:,1)

     
    deallocate(mult,arr,xd,xu,idx_r)

end subroutine comp_weight

subroutine eta_table_init(self,eta_ohm_tbl,eta_ad_tbl)
    class (non_ideal_t):: self
    real, dimension(:,:,:)    :: eta_Ohm_tbl, eta_AD_tbl
    real    :: rho0,rhoend,crho,T0,Tend,cT,E0,Eend,cE,B0,Bend,cB
    real, dimension(200) :: idxarr,rhoarr
    real, dimension(60)  :: Tarr
    real, dimension(150) :: Barr
    real, dimension(60)  :: Earr
    integer :: i,j
    real    :: nrho, nT, nB, nE, kB, mu, mp, NA

    kB = self%kB             ! Boltzmann constant
    mu = self%mu             ! mean molecular weight
    mp = self%mp             ! Proton mass  
    NA = self%NA             ! Avogadro constant

    nrho = 200. ; nT = 60. ; nE = nT ; nB = 150.
    rho0 = 300. * mu * mp ; rhoend = 1.7688710202180409e25 * mu * mp
    crho= (rhoend/rho0)**(1. / (nrho - 1.) ) 

    T0 = 10.    ; Tend   = 1.2449399509181750e5
    cT = (Tend/T0)**(1. / (nT - 1.) ) 

    E0 = 3./2 * kB * NA * T0    ; Eend   = 3./2 * kB * NA * Tend
    cE = (Eend/E0)**(1. / (nE - 1.) ) 

    B0 = 1e-10  ; Bend   = 1e10
    cB = (Bend/B0)**(1. / (nB - 1.) ) 

    do i=1,max(int(nrho),int(nT),int(nB))
      idxarr(i) = real(i-1)
    end do
    self%rhoarr = rho0 * crho**idxarr(1:int(nrho))
    self%Tarr   = T0   * cT**idxarr(1:int(nT))
    self%Earr   = E0   * cE**idxarr(1:int(nE))
    self%Barr   = B0   * cB**idxarr(1:int(nB))
    
    self%rho0 = rho0 
    self%crho = crho 
    self%T0   = T0   
    self%cT   = cT   
    self%E0   = E0   
    self%cE   = cE   
    self%B0   = B0   
    self%cB   = cB   

    print*,'call read table'    
    call read_table(eta_Ohm_tbl,eta_AD_tbl)
    print*,'read table was called!'    

end subroutine eta_table_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to read in resistivites from Marchand et al. 2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_table(eta_ohm_tbl, eta_ad_tbl) 
implicit none
integer   :: alpha           ! Number of bins
integer   :: Nvarchimie      ! Number of species
integer   :: nvar            ! Equal to Nvarchimie
integer   :: nion=9          ! Number of non-grain species

integer,parameter::dp=kind(1.0d0) ! default
!!!! Grains info
real(dp), allocatable, dimension(:) :: x_g        ! Abundance of grains per bins
real(dp), allocatable, dimension(:) :: r_g        ! Radius of grains per bins
real(dp), allocatable, dimension(:) :: m_g        ! Mass of grains per bins

!!!! Constants
real(dp), parameter :: pi=3.1415927410125732422_dp  
real(dp), parameter :: mp=1.6726d-24     ! Proton mass in g
real(dp), parameter :: me=9.1094d-28     ! Electron mass in g
real(dp), parameter :: e=4.803204d-10    ! Electron charge in cgs
real(dp), parameter :: Kb = 1.3807d-16   ! Boltzmann constant erg/K
real(dp), allocatable, dimension(:) :: q ! Electric charge 
real(dp), allocatable, dimension(:) :: m ! Mass

!!!! bins of grains
real(dp), parameter :: rho_s=2.3_dp        ! Bulk density g.cc
real(dp), parameter :: a_0=0.0375d-4       ! Reference radius cm
real(dp), parameter :: a_min=0.0181d-4     ! Minimum radius cm
real(dp), parameter :: a_max=0.9049d-4     ! Maximum radius cm
real(dp), parameter :: zeta=a_min/a_max    ! a_min/a_max
real(dp), parameter :: lambda_pow=-3.5d0   ! Coeff power law
real(dp)            :: rho_gtot            ! Grain density
real(dp) :: Lp1,Lp3,Lp4

! resistivites (cf Kunz & Mouschovias 2009, Marchand et al. 2016)
real(dp), allocatable, dimension(:)  :: sigma             
real(dp), allocatable, dimension(:)  :: zetas              
real(dp), allocatable, dimension(:)  :: phi
real(dp), allocatable, dimension(:)  :: tau_sn,tau_gp,tau_gm
real(dp), allocatable, dimension(:,:)  :: tau_inel
real(dp), allocatable, dimension(:)  :: gamma_zeta,gamma_omega
real(dp), allocatable, dimension(:)  :: omega,omega_bar
real(dp), parameter        :: c_l=299792458.d2     ! c in cm/s
real(dp)            :: sigv, muuu        
integer  :: iB,iH,iT,iX
real(dp) :: B,bmaxchimie,nH,T,xi,sigH,sigO,sigP

integer :: cmp_X

! For loops
integer :: tchimie   ! tchmie steps in temperature
real(dp) :: dtchimie   ! dtchmie step in temperature
real(dp) :: tminchimie   ! min tchmie 
integer :: nchimie   ! nchmie steps in density 
real(dp) :: dnchimie   ! dnchmie step in density 
real(dp) :: nminchimie   ! min nchmie  
integer :: bchimie   ! bchmie steps in B field
real(dp) :: dbchimie   ! dbchmie step in B field
real(dp) :: bminchimie   ! min bchmie 
integer :: xichimie ! Steps in ionisation rate
real(dp) :: dxichimie ! dxichimie step in ion rate
real(dp) :: ximinchimie ! min xichimie
integer :: nislin,tislin,xiislin ! Linear scale or not
real(dp),allocatable,dimension(:,:,:,:,:)::conductivities ! resistivities
real(dp),allocatable,dimension(:,:,:,:)::abundances ! abundances
integer :: i,j,k,ij


real(dp) :: sig_p, sig_perp, sig_h, eta_ohm, eta_hall, eta_ad

real, dimension(:,:,:) :: eta_Ohm_tbl, eta_AD_tbl

   open(42,file='../../tables/non-ideal_tbl/Table_abundances.dat', status='old')
   read(42,*) nchimie, tchimie, xichimie, nvar       ! Number of  steps in density, temperature and ionisation rate
   read(42,*) nislin, tislin, xiislin                      ! 1: Constant step in logscale, 0: Not constant step in logscale
   read(42,*)
   read(42,*)
   read(42,*)
   allocate(abundances(-2:nvar,nchimie,tchimie,xichimie)) 
   ! abundances(-2,:,:,:)=density
   ! abundances(-1,:,:,:)=temperature
   ! abundances(0,:,:,:)=ionisation rate
   ! abundances(n>°,:,:,:)=relative abundance of species n


   do ij=1,xichimie
     do i=1,tchimie
        do j=1,nchimie  
           read(42,*)abundances(-2:nvar,j,i,ij)
        end do
        read(42,*)
        read(42,*)
     end do
   end do
   close(42)
   nminchimie=(abundances(-2,1,1,1))   !min density
   tminchimie=(abundances(-1,1,1,1))   !min temperature
   ximinchimie=(abundances(0,1,1,1))   !min ionisation rate
   ! Density, temperature and IR steps in log (if nislin, tislin, xiislin = 1)
   dnchimie=(log10(abundances(-2,nchimie,1,1))-log10(abundances(-2,1,1,1)))/(nchimie-1)
   dtchimie=(log10(abundances(-1,1,tchimie,1))-log10(abundances(-1,1,1,1)))/(tchimie-1)
   dxichimie=(log10(abundances(0,1,1,xichimie))-log10(abundances(0,1,1,1)))/(xichimie-1)


   !Allocating arrays
   alpha=(nvar-nion)/3

   allocate(x_g(alpha))
   allocate(r_g(alpha))
   allocate(m_g(alpha))
   allocate(q(nvar))
   allocate(m(nvar))
   allocate(sigma(nvar))
   allocate(zetas(nvar))
   allocate(phi(nvar))
   allocate(tau_sn(nvar))
   allocate(tau_gp(nvar))
   allocate(tau_gm(nvar))
   allocate(gamma_zeta(nvar))
   allocate(gamma_omega(nvar))
   allocate(omega(nvar))
   allocate(omega_bar(nvar))
   allocate(tau_inel(nvar,nvar))
   

   ! Computation of grains radii and species mass
   Lp1=dble(lambda_pow+1)
   Lp3=dble(lambda_pow+3)
   Lp4=dble(lambda_pow+4)

    ! Grain radius
   if  (alpha==1) then
     r_g(1)=a_0
   else
     do  i=1,alpha    ! cf Marchand et al. 2016
       r_g(i)=a_min*zeta**(-dble(i)/dble(alpha)) * &
            & dsqrt( Lp1/Lp3* (1d0-zeta**(Lp3/dble(alpha)))/(1d0-zeta**(Lp1/dble(alpha))))
     end do
   end if


   q(:)=1.d0*e    ! cations
   q(1)=-1.d0*e   ! electrons
   do  i=nion+1,nvar
      if (mod(i-nion,3)==1) q(i)=1.d0*e   ! g+
      if (mod(i-nion,3)==2) q(i)=-1.d0*e  ! g-
      if (mod(i-nion,3)==0) q(i)=0.d0     ! g0
   end do
   m(:)=0.d0
   m(1) = me           ! e-
   m(2) = 23.5d0*mp    ! ions metalliques
   m(3) = 29.d0*mp     ! ions moleculaires
   m(4) = 3*mp         ! H3+
   m(5) = mp           ! H+
   m(6) = 12.d0*mp     ! C+
   m(7) = 4.d0*mp      ! He+
   m(8) = 39.098*mp    ! K+
   m(9) = 22.99d0*mp   ! Na+
   do  i=1,alpha       ! masse des grains
      m_g(i)=4.d0/3.d0*pi*r_g(i)**3*rho_s
      m(nion+3*(i-1)+1:nion+3*i)=m_g(i)
   end do


   ! values for Btable
   bminchimie=1.d-10
   bmaxchimie=1.d10               ! ok for first core in nimhd. maybe not enough for second core.
   bchimie=150
   dbchimie=(log10(bmaxchimie)-log10(bminchimie))/real((bchimie-1),dp)

   allocate(conductivities(0:3,1:nchimie,1:tchimie,1:xichimie,1:bchimie))
   !conductivities(1,:,:,:)= sigma_//
   !conductivities(2,:,:,:)= sigma_perp
   !conductivities(3,:,:,:)= sigma_hall
   !conductivities(0,:,:,:)= sigma_hall sign


  ! Computation of resistivities

   tau_sn      = 0.0_dp
   omega       = 0.0_dp
   sigma       = 0.0_dp
   phi         = 0.0_dp
   zetas       = 0.0_dp
   gamma_zeta  = 0.0_dp
   gamma_omega = 0.0_dp
   omega_bar   = 0.0_dp


do  iB=1,bchimie
do  iX=7,7 ! 7 : cosmic-ray ionisation rate of 1e-17 s^-1
!do  iX=1,xichimie
do  iT=1,tchimie
do  iH=1,nchimie

    nh=abundances(-2,iH,iT,iX)                ! density (.cc) of current point
    B =10.0d0**(log10(bminchimie)+dble(iB-1)*dbchimie)  ! Magnetic field
    T =abundances(-1,iH,iT,iX)                ! Temperature
    xi =abundances(0,iH,iT,iX)                ! Ionisation rate


    do i=1,nion
      if  (i==1) then  ! electron
        sigv=3.16d-11 * (dsqrt(8d0*kB*1d-7*T/(pi*me*1d-3))*1d-3)**1.3d0
        tau_sn(i) = 1.d0/1.16d0*(m(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*sigv)
      else if (i>=2 .and. i<=nion) then ! ions
        muuu=m(i)*2d0*mp/(m(i)+2d0*mp)
        if (i==2 .or. i==3) then
          sigv=2.4d-9 *(dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**0.6d0
        else if (i==4) then
          sigv=2d-9 * (dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**0.15d0
        else if (i==5) then
          sigv=3.89d-9 * (dsqrt(8d0*kB*1d-7*T/(pi*muuu*1d-3))*1d-3)**(-0.02d0)
        else
          sigv=1.69d-9
        end if
        tau_sn(i) = 1.d0/1.14d0*(m(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*sigv)
      end if
      omega(i) = q(i)*B/(m(i)*c_l)
      sigma(i) = abundances(i,iH,iT,iX)*nH*(q(i))**2*tau_sn(i)/m(i)
      phi(i) = 0.d0
      zetas(i) = 0.d0
      gamma_zeta(i) = 1.d0
      gamma_omega(i) = 1.d0
      omega_bar(i) = 0.d0
    end do
    
    do  i=1,alpha   ! grains
      tau_sn(nion+1+3*(i-1))=1.d0/1.28d0*(m_g(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*(pi*r_g(i)**2*(8.d0*Kb*T/(pi*2.d0*mp))**0.5))
      omega(nion+1+3*(i-1)) = q(nion+1+3*(i-1))*B/(m_g(i)*c_l)
      sigma(nion+1+3*(i-1)) = abundances(nion+1+3*(i-1),iH,iT,iX)*nH*(q(nion+1+3*(i-1)))**2*tau_sn(nion+1+3*(i-1))/m_g(i)
    
      tau_sn(nion+2+3*(i-1))=tau_sn(nion+1+3*(i-1))
      omega(nion+2+3*(i-1)) = q(nion+2+3*(i-1))*B/(m_g(i)*c_l)
      sigma(nion+2+3*(i-1)) = abundances(nion+2+3*(i-1),iH,iT,iX)*nH*(q(nion+2+3*(i-1)))**2*tau_sn(nion+2+3*(i-1))/m_g(i)
    
    end do
  
    sigP=0.d0
    sigO=0.d0
    sigH=0.d0

    do i=1,nvar
       sigP=sigP+sigma(i)
       sigO=sigO+sigma(i)/(1.d0+(omega(i)*tau_sn(i))**2)
       sigH=sigH-sigma(i)*omega(i)*tau_sn(i)/(1.d0+(omega(i)*tau_sn(i))**2)
    end do

    conductivities(1,iH,iT,iX,iB)=log10(sigP)
    conductivities(2,iH,iT,iX,iB)=log10(sigO)
    !conductivities(3,iH,iT,iX,iB)=log10(abs(sigH))
    !conductivities(0,iH,iT,iX,iB)=sign(1.0_dp,sigH)

    !conductivities(:,:,:,:) contains the log10(conductivities)


    eta_ohm  = 1d0/sigP
    eta_hall = sigH/(sigO**2 + sigH**2)
    eta_ad   = sigO/(sigO**2 + sigH**2)-1d0/sigP

    eta_Ohm_tbl(iH,iT,iB) = eta_ohm 
    eta_AD_tbl (iH,iT,iB) = eta_ad 
    if (eta_Ohm .lt. 0) eta_Ohm_tbl (iH,iT,iB) = 0d0 
    if (eta_ad .lt. 0)  eta_AD_tbl (iH,iT,iB) = 0d0 

end do
end do
end do
end do
!open(unit=3,file='eta_AD_tbl.dat',ACCESS='APPEND')
!write(3,*) eta_ad_tbl
!close(3)

end subroutine

end module
! modif nimhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
