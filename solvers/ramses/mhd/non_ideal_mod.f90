!> =============================================================================
!> non-ideal MHD module
!> =============================================================================
MODULE non_ideal_mod
  USE const
  USE io_mod
  USE io_unit_mod
  USE gpatch_mod
  USE link_mod
  use hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE force_mod
  USE index_mod
  implicit none
  private
  type, public :: non_ideal_t
    real:: eta_Ohm, etaMD, gamma_AD, rho_ions
    !real,pointer,dimension(:,:,:,:) :: uin, qin, fluxambdiff, emfambdiff
    real,pointer,dimension(:,:,:,:,:) :: flux
    logical:: mhd_AD, mhd_Ohm, ntestDADM
  contains
    procedure:: init
    procedure:: non_ideal_comp
    procedure:: flux_upd
    procedure:: non_ideal_emf_up
    procedure:: computejb2
    procedure:: computambip
    procedure:: computdifmag            
  end type
  type(non_ideal_t), public:: non_ideal
CONTAINS

!> =============================================================================
!> Read parameters
!> =============================================================================
SUBROUTINE init (self)
  class (non_ideal_t):: self
  integer:: iostat
  real:: gamma_AD=75., rho_ions=1.0, eta_Ohm=0.0, etaMD=0.0
  logical:: mhd_AD=.false., mhd_Ohm=.false., ntestDADM=.true.
  namelist /non_ideal_params/ gamma_AD, eta_Ohm, etaMD, rho_ions, mhd_AD, mhd_Ohm, ntestDADM
  !---------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, non_ideal_params,iostat=iostat)
  if (io%master) write (*, non_ideal_params)
  self%mhd_Ohm  = mhd_Ohm
  self%eta_Ohm  = eta_Ohm
  self%mhd_AD   = mhd_AD
  self%gamma_AD = gamma_AD
  self%rho_ions = rho_ions
  self%etaMD    = etaMD
END SUBROUTINE init

subroutine non_ideal_comp(self,bf,uin,qin,fluxambdiff,emfambdiff,u_max,dx,dt,ng)
  class (non_ideal_t):: self
  integer :: ng
  real :: u_max
  real(8) :: dx, dt
  real,dimension(1:ng,1:ng,1:ng,1:8)::uin, qin, fluxambdiff, emfambdiff
  real,dimension(1:ng,1:ng,1:ng,1:3)::bf
  ! Output electromotive force
  real(8),dimension(1:ng,1:ng,1:ng)::emfx
  real(8),dimension(1:ng,1:ng,1:ng)::emfy
  real(8),dimension(1:ng,1:ng,1:ng)::emfz

  ! Output courant vector in the cell
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::bmagij
  real(8),dimension(1:ng,1:ng,1:ng)::jcentersquare,jxbsquare
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::fluxmd,fluxh,fluxad
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::emfohmdiss,fluxohm 

  bmagij=0.d0
  emfohmdiss=0.d0
  fluxohm=0.d0
  jcentersquare=0.d0
  fluxmd=0.d0
  fluxh=0.d0
  fluxad=0.d0
  
  bemfx=0.d0
  bemfy=0.d0
  bemfz=0.d0
  jemfx=0.d0
  jemfy=0.d0
  jemfz=0.d0
  florentzx=0.d0
  florentzy=0.d0
  florentzz=0.d0

  if((self%mhd_AD).or.(self%mhd_Ohm)) then
  !  compute Lorentz Force with current
     call computejb2(self,bf,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,dx,ng)
  endif
  
  ! AMBIPOLAR DIFFUSION
  
  if(self%mhd_AD) then
     call computambip(self,uin,qin,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,dx,dt,ng)
     u_max = max(u_max, maxval(uin(:,:,:,6)**2+uin(:,:,:,7)**2+uin(:,:,:,8)**2 / (self%gamma_AD * self%rho_ions * uin(:,:,:,1) * dx)) )  
  endif
  
  
  ! OHMIC DISSIPATION
  
  if(self%mhd_Ohm) then
     call computdifmag(self,uin,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,dx,dt,ng)
     u_max = max(u_max,self%eta_Ohm/ dx)
  endif
  

end subroutine

subroutine flux_upd(self,flux,fluxambdiff, dx, dt, ndim, ng)
   class (non_ideal_t):: self
   
   integer :: i, j, k, ndim, ng
   real,dimension(1:ng,1:ng,1:ng,1:3)::fluxohm
   real,dimension(1:ng,1:ng,1:ng,1:3)::fluxambdiff
   real(8) :: dx, dt

   real, dimension(1:ng,1:ng,1:ng,1:8,1:3):: flux 

   fluxohm=0d0
   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then
      flux(:,:,:,5,ndim) = flux(:,:,:,5,ndim) + (fluxambdiff(:,:,:,ndim)+fluxohm(:,:,:,ndim) ) * dt / dx
   end if
end subroutine


subroutine non_ideal_emf_up(self,emfup,emfambdiff, dx, dt, ndim, ng)
   class (non_ideal_t):: self
   integer :: ndim, ng
   real,dimension(1:ng,1:ng,1:ng)::emfup
   real,dimension(1:ng,1:ng,1:ng,1:3)::emfohmdiss,emfambdiff
   real(8) :: dx, dt

   emfohmdiss = 0d0

   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then

#if NDIM==1

   emfup(:,:,:) = ( emfambdiff(:,:,:,ndim) + emfohmdiss(:,:,:,ndim) ) * dt / dx
#else
   emfup(2:ng-1,2:ng-1,2:ng-1) = emfup(2:ng-1,2:ng-1,2:ng-1) + &
        ( emfambdiff(2:ng-1,2:ng-1,2:ng-1,ndim) + emfohmdiss(2:ng-1,2:ng-1,2:ng-1,ndim) ) * dt / dx
#endif     
   end if
end subroutine


subroutine computejb2(self,bf,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,dx,ng)
  class (non_ideal_t):: self

  integer :: i, j, k,  ibegin, iend, jbegin, jend, kbegin, kend, ndim, ng
  real,dimension(:,:,:,:)::bf 
  real,dimension(:,:,:,:)::q 


  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::fluxmd,fluxh,fluxad
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::bmagij
!  real,dimension(1:self%gn(1),1:self%gn(2),1:self%gn(3),1:3)::jcell


  ! declare local variables
  integer :: m, n 

  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bmagijbis
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::jface
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bcenter
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::fluxbis,fluxter,fluxquat
  real(8)::dx,dy,dz

  dy = dx ; dz = dx
   ! Compute before the loop because bmagijbis has to be computed with index i+1/j+1/k+1 in loop over i/j/k
   bmagijbis(2:ng-1,2:ng-1,2:ng-1,1)= 0.25d0* (bf(2:ng-1,2:ng-1,2:ng-1,1) + bf(2:ng  ,1:ng-1,2:ng  ,1) + bf(2:ng,2:ng  ,1:ng-1,1) + bf(2:ng  ,1:ng-1,1:ng-1,1) )
   bmagijbis(2:ng-1,2:ng-1,2:ng-1,2)= 0.25d0* (bf(2:ng-1,2:ng-1,2:ng-1,2) + bf(1:ng-1,2:ng  ,2:ng  ,2) + bf(2:ng,2:ng  ,1:ng-1,2) + bf(1:ng-1,2:ng  ,1:ng-1,2) ) 
   bmagijbis(2:ng-1,2:ng-1,2:ng-1,3)= 0.25d0* (bf(2:ng-1,2:ng-1,2:ng-1,3) + bf(1:ng-1,2:ng  ,2:ng  ,3) + bf(2:ng,1:ng-1,2:ng  ,3) + bf(1:ng-1,1:ng-1,2:ng  ,3) ) 
   ! Define the 1 and 2 components of bcenter, since they are needed further below, when computing the current on the cell interfaces (jface) 
   ! magnetic field at center of cells
   bcenter(1:ng,1:ng,1:ng,1:3)=q(1:ng,1:ng,1:ng,6:8)

   !!!!!!!!!!!!!!!!!!
   ! EMF x
   !!!!!!!!!!!!!!!!!!
   bemfx(2:ng-1,2:ng-1,2:ng-1,1)=0.25d0* ( q (2:ng-1,2:ng-1,2:ng-1,6) + q (2:ng-1,1:ng-2,2:ng-1,6) + q(2:ng-1,2:ng-1,1:ng-2,6) + q(2:ng-1,1:ng-2,1:ng-2,6) )
   bemfx(2:ng-1,2:ng-1,2:ng-1,2)=0.5d0 * ( bf(2:ng-1,2:ng-1,2:ng-1,2) + bf(2:ng-1,2:ng-1,1:ng-2,2) )
   bemfx(2:ng-1,2:ng-1,2:ng-1,3)=0.5d0 * ( bf(2:ng-1,2:ng-1,2:ng-1,3) + bf(2:ng-1,1:ng-2,2:ng-1,3) )
   !!!!!!!!!!!!!!!!!!
   ! EMF y
   !!!!!!!!!!!!!!!!!!
   bemfy(2:ng-1,2:ng-1,2:ng-1,1)=0.5d0 * ( bf(2:ng-1,2:ng-1,2:ng-1,1) + bf(2:ng-1,2:ng-1,1:ng-2,1) )
   bemfy(2:ng-1,2:ng-1,2:ng-1,2)=0.25d0* ( q (2:ng-1,2:ng-1,2:ng-1,7) + q (1:ng-2,2:ng-1,2:ng-1,7) + q(2:ng-1,2:ng-1,1:ng-2,7) + q(1:ng-2,2:ng-1,1:ng-2,7) )
   bemfy(2:ng-1,2:ng-1,2:ng-1,3)=0.5d0 * ( bf(1:ng-2,2:ng-1,2:ng-1,3) + bf(2:ng-1,2:ng-1,2:ng-1,3) )
   !!!!!!!!!!!!!!!!!!
   ! EMF z
   !!!!!!!!!!!!!!!!!!
   bemfz(2:ng-1,2:ng-1,2:ng-1,1)=0.5d0 * ( bf(2:ng-1,2:ng-1,2:ng-1,1) + bf(2:ng-1,1:ng-2,2:ng-1,1) )
   bemfz(2:ng-1,2:ng-1,2:ng-1,2)=0.5d0 * ( bf(2:ng-1,2:ng-1,2:ng-1,2) + bf(1:ng-2,2:ng-1,2:ng-1,2) )
   bemfz(2:ng-1,2:ng-1,2:ng-1,3)=0.25d0* ( q (2:ng-1,2:ng-1,2:ng-1,8) + q (1:ng-2,2:ng-1,2:ng-1,8) + q(2:ng-1,1:ng-2,2:ng-1,8) + q(1:ng-2,1:ng-2,2:ng-1,8) )
   do m=1,3
   ! bmagij is the value of the magnetic field Bi where Bj 
   ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
   ! and we can write it Bx,y  
   !! m+5 mandatory cf Bx=uin(i,j,k,6)
      bmagij(2:ng-1,2:ng-1,2:ng-1,m,m)=bf(2:ng-1,2:ng-1,2:ng-1,m)
   end do
   ! case Bx,y
   bmagij(2:ng-1,2:ng-1,2:ng-1,1,2) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,6) + q(2:ng-1,1:ng-2,2:ng-1,6) )
   ! case Bx,z
   bmagij(2:ng-1,2:ng-1,2:ng-1,1,3) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,6) + q(2:ng-1,2:ng-1,1:ng-2,6) )
   ! case By,y
   bmagij(2:ng-1,2:ng-1,2:ng-1,2,1) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,7) + q(1:ng-2,2:ng-1,2:ng-1,7) )
   ! case By,z
   bmagij(2:ng-1,2:ng-1,2:ng-1,2,3) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,7) + q(2:ng-1,2:ng-1,1:ng-2,7) )
   ! case Bz,x
   bmagij(2:ng-1,2:ng-1,2:ng-1,3,1) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,8) + q(1:ng-2,2:ng-1,2:ng-1,8) )
   ! case Bz,y
   bmagij(2:ng-1,2:ng-1,2:ng-1,3,2) = 0.5d0 * (q(2:ng-1,2:ng-1,2:ng-1,8) + q(2:ng-1,1:ng-2,2:ng-1,8) )
   jemfx(2:ng-1,2:ng-1,2:ng-1,1) = (bf(2:ng-1,2:ng-1,2:ng-1,3) - bf(2:ng-1,1:ng-2,2:ng-1,3) ) / dy - (bf(2:ng-1,2:ng-1,2:ng-1,2) - bf(2:ng-1,2:ng-1,1:ng-2,2) ) / dz 
   jemfx(2:ng-1,2:ng-1,2:ng-1,2) = (bmagij(2:ng-1,2:ng-1,2:ng-1,1,2) - bmagij(2:ng-1,2:ng-1,1:ng-2,1,2) ) / dz - (bmagijbis(1:ng-2,2:ng-1,2:ng-1,3) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,3) ) / dx
   jemfx(2:ng-1,2:ng-1,2:ng-1,3) = (bmagijbis(1:ng-2,2:ng-1,2:ng-1,2) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,2) ) / dx - (bmagij(2:ng-1,2:ng-1,2:ng-1,1,3)- bmagij(2:ng-1,1:ng-2,2:ng-1,1,3) ) / dy
   jemfy(2:ng-1,2:ng-1,2:ng-1,1) = (bmagijbis(2:ng-1,3:ng  ,2:ng-1,3) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,3) ) / dy - (bmagij(2:ng-1,2:ng-1,2:ng-1,2,1) - bmagij(2:ng-1,2:ng-1,1:ng-2,2,1) ) / dz
   jemfy(2:ng-1,2:ng-1,2:ng-1,2) = (bf(2:ng-1,2:ng-1,2:ng-1,1) - bf(2:ng-1,2:ng-1,1:ng-2,1) ) / dz - (bf(2:ng-1,2:ng-1,2:ng-1,3) - bf(1:ng-2,2:ng-1,2:ng-1,3) ) / dx
   jemfy(2:ng-1,2:ng-1,2:ng-1,3) = (bmagij(2:ng-1,2:ng-1,2:ng-1,2,3) - bmagij(1:ng-2,2:ng-1,2:ng-1,2,3) ) / dx - (bmagijbis(2:ng-1,3:ng  ,2:ng-1,1) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,1) ) / dy
   jemfz(2:ng-1,2:ng-1,2:ng-1,1) = (bmagij(2:ng-1,2:ng-1,2:ng-1,3,1) - bmagij(2:ng-1,1:ng-2,2:ng-1,3,1) ) / dy - (bmagijbis(2:ng-1,2:ng-1,3:ng  ,2) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,2) ) / dz
   jemfz(2:ng-1,2:ng-1,2:ng-1,2) = (bmagijbis(2:ng-1,2:ng-1,3:ng  ,1) - bmagijbis(2:ng-1,2:ng-1,2:ng-1,1) ) / dz - (bmagij(2:ng-1,2:ng-1,2:ng-1,3,2) - bmagij(1:ng-2,2:ng-1,2:ng-1,3,2) ) / dx
   jemfz(2:ng-1,2:ng-1,2:ng-1,3) = (bf(2:ng-1,2:ng-1,2:ng-1,2) - bf(1:ng-2,2:ng-1,2:ng-1,2) ) / dx - (bf(2:ng-1,2:ng-1,2:ng-1,1) - bf(2:ng-1,1:ng-2,2:ng-1,1) ) / dy


if(self%mhd_AD) then

! EMF x
   florentzx(2:ng-1,2:ng-1,2:ng-1,1)=jemfx(2:ng-1,2:ng-1,2:ng-1,2)*bemfx(2:ng-1,2:ng-1,2:ng-1,3) - jemfx(2:ng-1,2:ng-1,2:ng-1,3)*bemfx(2:ng-1,2:ng-1,2:ng-1,2)
   florentzx(2:ng-1,2:ng-1,2:ng-1,2)=jemfx(2:ng-1,2:ng-1,2:ng-1,3)*bemfx(2:ng-1,2:ng-1,2:ng-1,1) - jemfx(2:ng-1,2:ng-1,2:ng-1,1)*bemfx(2:ng-1,2:ng-1,2:ng-1,3)
   florentzx(2:ng-1,2:ng-1,2:ng-1,3)=jemfx(2:ng-1,2:ng-1,2:ng-1,1)*bemfx(2:ng-1,2:ng-1,2:ng-1,2) - jemfx(2:ng-1,2:ng-1,2:ng-1,1)*bemfx(2:ng-1,2:ng-1,2:ng-1,2)
   florentzy(2:ng-1,2:ng-1,2:ng-1,1)=jemfy(2:ng-1,2:ng-1,2:ng-1,2)*bemfy(2:ng-1,2:ng-1,2:ng-1,3) - jemfy(2:ng-1,2:ng-1,2:ng-1,3)*bemfy(2:ng-1,2:ng-1,2:ng-1,2)
   florentzy(2:ng-1,2:ng-1,2:ng-1,2)=jemfy(2:ng-1,2:ng-1,2:ng-1,3)*bemfy(2:ng-1,2:ng-1,2:ng-1,1) - jemfy(2:ng-1,2:ng-1,2:ng-1,1)*bemfy(2:ng-1,2:ng-1,2:ng-1,3)
   florentzy(2:ng-1,2:ng-1,2:ng-1,3)=jemfy(2:ng-1,2:ng-1,2:ng-1,1)*bemfy(2:ng-1,2:ng-1,2:ng-1,2) - jemfy(2:ng-1,2:ng-1,2:ng-1,1)*bemfy(2:ng-1,2:ng-1,2:ng-1,2)
   florentzz(2:ng-1,2:ng-1,2:ng-1,1)=jemfz(2:ng-1,2:ng-1,2:ng-1,2)*bemfz(2:ng-1,2:ng-1,2:ng-1,3) - jemfz(2:ng-1,2:ng-1,2:ng-1,3)*bemfz(2:ng-1,2:ng-1,2:ng-1,2)
   florentzz(2:ng-1,2:ng-1,2:ng-1,2)=jemfz(2:ng-1,2:ng-1,2:ng-1,3)*bemfz(2:ng-1,2:ng-1,2:ng-1,1) - jemfz(2:ng-1,2:ng-1,2:ng-1,1)*bemfz(2:ng-1,2:ng-1,2:ng-1,3)
   florentzz(2:ng-1,2:ng-1,2:ng-1,3)=jemfz(2:ng-1,2:ng-1,2:ng-1,1)*bemfz(2:ng-1,2:ng-1,2:ng-1,2) - jemfz(2:ng-1,2:ng-1,2:ng-1,1)*bemfz(2:ng-1,2:ng-1,2:ng-1,2)
  

endif


! computation of current on faces

if((self%mhd_AD).or.(self%mhd_Ohm)) then

! face at i-1/2,j,k

   jface(2:ng-1,2:ng-1,2:ng-1,1,1)= (bemfz  (2:ng-1,3:ng  ,2:ng-1,3) - bemfz  (2:ng-1,2:ng-1,2:ng-1,3) ) / dy - (bemfy  (2:ng-1,2:ng-1,3:ng  ,2) - bemfy  (2:ng-1,2:ng-1,2:ng-1,2) ) / dz
   jface(2:ng-1,2:ng-1,2:ng-1,2,1)= (bemfy  (2:ng-1,2:ng-1,3:ng  ,1) - bemfy  (2:ng-1,2:ng-1,2:ng-1,1) ) / dz - (bcenter(2:ng-1,2:ng-1,2:ng-1,3) - bcenter(1:ng-3,2:ng-1,2:ng-1,3) ) / dx
   jface(2:ng-1,2:ng-1,2:ng-1,3,1)= (bcenter(2:ng-1,2:ng-1,2:ng-1,2) - bcenter(1:ng-3,2:ng-1,2:ng-1,2) ) / dx - (bemfz  (2:ng-1,3:ng  ,2:ng-1,1) - bemfz  (2:ng-1,2:ng-1,2:ng-1,1) ) / dy
   jface(2:ng-1,2:ng-1,2:ng-1,1,2)= (bcenter(2:ng-1,2:ng-1,2:ng-1,3) - bcenter(2:ng-1,2:ng-1,1:ng-3,3) ) / dy - (bemfx  (2:ng-1,2:ng-1,3:ng  ,2) - bemfx  (2:ng-1,2:ng-1,2:ng-1,2) ) / dz
   jface(2:ng-1,2:ng-1,2:ng-1,2,2)= (bemfx  (2:ng-1,2:ng-1,3:ng  ,1) - bemfx  (2:ng-1,2:ng-1,2:ng-1,1) ) / dz - (bemfz  (3:ng  ,2:ng-1,2:ng-1,3) - bemfz  (2:ng-1,2:ng-1,2:ng-1,3) ) / dx
   jface(2:ng-1,2:ng-1,2:ng-1,3,2)= (bemfz  (3:ng  ,2:ng-1,2:ng-1,2) - bemfz  (2:ng-1,2:ng-1,2:ng-1,3) ) / dx - (bcenter(2:ng-1,2:ng-1,2:ng-1,1) - bcenter(2:ng-1,1:ng-3,2:ng-1,1) ) / dy
   jface(2:ng-1,2:ng-1,2:ng-1,1,3)= (bemfx  (2:ng-1,3:ng  ,2:ng-1,3) - bemfx  (2:ng-1,2:ng-1,2:ng-1,3) ) / dy - (bcenter(2:ng-1,2:ng-1,2:ng-1,2) - bcenter(2:ng-1,2:ng-1,1:ng-3,2) ) / dz             
   jface(2:ng-1,2:ng-1,2:ng-1,2,3)= (bcenter(2:ng-1,2:ng-1,2:ng-1,1) - bcenter(2:ng-1,2:ng-1,1:ng-3,1) ) / dz - (bemfy  (3:ng  ,2:ng-1,2:ng-1,3) - bemfy  (2:ng-1,2:ng-1,2:ng-1,3) ) / dx             
   jface(2:ng-1,2:ng-1,2:ng-1,3,3)= (bemfy  (3:ng  ,2:ng-1,2:ng-1,2) - bemfy  (2:ng-1,2:ng-1,2:ng-1,2) ) / dx - (bemfx  (2:ng-1,3:ng  ,2:ng-1,1) - bemfx  (2:ng-1,2:ng-1,2:ng-1,1) ) / dy            
   do n=1,3
     fluxbis(2:ng-1,2:ng-1,2:ng-1,1,n) = jface(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n) - jface(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n)
     fluxbis(2:ng-1,2:ng-1,2:ng-1,2,n) = jface(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n) - jface(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n)
     fluxbis(2:ng-1,2:ng-1,2:ng-1,3,n) = jface(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n) - jface(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n)
   end do
   fluxmd(2:ng-1,2:ng-1,2:ng-1,1)=fluxbis(2:ng-1,2:ng-1,2:ng-1,1,1)
   fluxmd(2:ng-1,2:ng-1,2:ng-1,2)=fluxbis(2:ng-1,2:ng-1,2:ng-1,2,2)
   fluxmd(2:ng-1,2:ng-1,2:ng-1,3)=fluxbis(2:ng-1,2:ng-1,2:ng-1,3,3)


endif

if(self%mhd_AD) then

   do n=1,3
      fluxter(2:ng-1,2:ng-1,2:ng-1,1,n)=fluxbis(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n) - fluxbis(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n)
      fluxter(2:ng-1,2:ng-1,2:ng-1,2,n)=fluxbis(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n) - fluxbis(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n)
      fluxter(2:ng-1,2:ng-1,2:ng-1,3,n)=fluxbis(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n) - fluxbis(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n)
   end do
   fluxh(2:ng-1,2:ng-1,2:ng-1,1)=fluxter(2:ng-1,2:ng-1,2:ng-1,1,1)
   fluxh(2:ng-1,2:ng-1,2:ng-1,2)=fluxter(2:ng-1,2:ng-1,2:ng-1,2,2)
   fluxh(2:ng-1,2:ng-1,2:ng-1,3)=fluxter(2:ng-1,2:ng-1,2:ng-1,3,3)
   do n=1,3                 
      fluxquat(2:ng-1,2:ng-1,2:ng-1,1,n)=fluxter(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n) - fluxter(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n)
      fluxquat(2:ng-1,2:ng-1,2:ng-1,2,n)=fluxter(2:ng-1,2:ng-1,2:ng-1,3,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n) - fluxter(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,3,n)
      fluxquat(2:ng-1,2:ng-1,2:ng-1,3,n)=fluxter(2:ng-1,2:ng-1,2:ng-1,1,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,2,n) - fluxter(2:ng-1,2:ng-1,2:ng-1,2,n)*bmagij(2:ng-1,2:ng-1,2:ng-1,1,n)
   end do
   fluxad(2:ng-1,2:ng-1,2:ng-1,1)=fluxquat(2:ng-1,2:ng-1,2:ng-1,1,1)
   fluxad(2:ng-1,2:ng-1,2:ng-1,2)=fluxquat(2:ng-1,2:ng-1,2:ng-1,2,2)
   fluxad(2:ng-1,2:ng-1,2:ng-1,3)=fluxquat(2:ng-1,2:ng-1,2:ng-1,3,3)


endif


end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(self,u,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,dx,dt,ng)
  class (non_ideal_t):: self

  integer :: i, j, k, ng
  real,dimension(:,:,:,:)::u 
  real,dimension(:,:,:,:)::q 

  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::bmagij
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::fluxmd

  real(8),dimension(1:ng,1:ng,1:ng,1:3):: emfohmdiss,fluxohm 
  real(8),dimension(1:ng,1:ng,1:ng)::jcentersquare

  ! declare local variables
    integer :: m, n, h


  ! WARNING following quantities defined with three components even
  ! if ndim<3 !
  real(8),dimension(1:ng,1:ng,1:ng,3)::jcenter
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::jface
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::jemf

  real(8)::tcell,rhocell,bcell,ionisrate,Cv

  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::tcellx,tcelly,tcellz,bsquarex,bsquarey,bsquarez,dtlim
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::rhox,rhoy,rhoz,epsx,epsy,epsz
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::rhof,pf,bsqf,epsf,tcellf

  integer , dimension(1:3) :: index_i,index_j,index_k

  real(8):: dx, dy, dz, dt

  index_i = (/1,0,0/)
  index_j = (/0,1,0/)
  index_k = (/0,0,1/)

  dy = dx ; dz = dx
  dtlim = dt !neil

  jemf(2:ng-1,2:ng-1,2:ng-1,1)=jemfx(2:ng-1,2:ng-1,2:ng-1,1)
  jemf(2:ng-1,2:ng-1,2:ng-1,2)=jemfy(2:ng-1,2:ng-1,2:ng-1,2)
  jemf(2:ng-1,2:ng-1,2:ng-1,3)=jemfz(2:ng-1,2:ng-1,2:ng-1,3)
  rhox=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(2:ng-1,1:ng-2,2:ng-1,1) + u(2:ng-1,2:ng-1,1:ng-2,1) + u(2:ng-1,1:ng-2,1:ng-2,1))
  rhoy=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(1:ng-2,2:ng-1,2:ng-1,1) + u(2:ng-1,2:ng-1,1:ng-2,1) + u(1:ng-2,2:ng-1,1:ng-2,1))
  rhoz=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(1:ng-2,2:ng-1,2:ng-1,1) + u(2:ng-1,1:ng-2,2:ng-1,1) + u(1:ng-2,1:ng-2,2:ng-1,1))
  !epsx=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,j-1,2:ng-1,nvar) + u(2:ng-1,2:ng-1,1:ng-2,nvar) + u(2:ng-1,1:ng-2,1:ng-2,nvar))
  !epsy=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(i-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,2:ng-1,1:ng-2,nvar) + u(1:ng-2,2:ng-1,1:ng-2,nvar))
  !epsz=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar) + u(i-1,2:ng-1,2:ng-1,nvar) + u(2:ng-1,1:ng-2,2:ng-1,nvar) + u(1:ng-2,1:ng-2,2:ng-1,nvar))
   if(self%mhd_Ohm)then
     bsquarex=bemfx(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfx(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfx(2:ng-1,2:ng-1,2:ng-1,3)**2
     bsquarey=bemfy(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfy(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfy(2:ng-1,2:ng-1,2:ng-1,3)**2
     bsquarez=bemfz(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfz(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfz(2:ng-1,2:ng-1,2:ng-1,3)**2
   end if
  tcellx=1.0d0
  tcelly=1.0d0
  tcellz=1.0d0
  ionisrate=1 !default_ionisrate 
              
  ! WARNING dB/dt=-curl(eta*J)
  do n=1,3
    emfohmdiss(2:ng-1,2:ng-1,2:ng-1,n)=-self%etaMD*jemf(2:ng-1,2:ng-1,2:ng-1,n)
  end do
  !!!!!!!!!!!!!!!!!!!!!!!
  !
  ! compute j at center of cells
  !
  ! mandatory for non isotherm case
  
  ! bmagij is the value of the magnetic field Bi where Bj 
  ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
  ! and we can write it Bx,y

  jcenter(2:ng-1,2:ng-1,2:ng-1,1)= (bmagij(2:ng-1,3:ng  ,2:ng-1,3,2) - bmagij(2:ng-1,2:ng-1,2:ng-1,3,2) ) / dz - (bmagij(2:ng-1,2:ng-1,3:ng  ,2,3) - bmagij(2:ng-1,2:ng-1,2:ng-1,2,3) ) / dy
  jcenter(2:ng-1,2:ng-1,2:ng-1,2)= (bmagij(2:ng-1,2:ng-1,3:ng  ,1,3) - bmagij(2:ng-1,2:ng-1,2:ng-1,1,3) ) / dx - (bmagij(3:ng  ,2:ng-1,2:ng-1,3,1) - bmagij(2:ng-1,2:ng-1,2:ng-1,3,1) ) / dz
  jcenter(2:ng-1,2:ng-1,2:ng-1,3)= (bmagij(3:ng  ,2:ng-1,2:ng-1,2,1) - bmagij(2:ng-1,2:ng-1,2:ng-1,2,1) ) / dy - (bmagij(2:ng-1,3:ng  ,2:ng-1,2,3) - bmagij(2:ng-1,2:ng-1,2:ng-1,2,3) ) / dx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  do h = 1,3
     
     rhof=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(2-index_i(h):ng-index_i(h),2-index_j(h):ng-index_j(h),2-index_k(h):ng-index_k(h),1))
     !epsf=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,nvar)+u(2-index_i(h):ng-index_i(h),2-index_j(h):ng-index_j(h),2-index_k(h):ng-index_k(h),nvar))
     bsqf=bmagij(2:ng-1,2:ng-1,2:ng-1,1,h)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,2,h)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,3,h)**2
     
     tcellf=1.0d0
     
     fluxohm(2:ng-1,2:ng-1,2:ng-1,h)=self%etaMD*fluxmd(2:ng-1,2:ng-1,2:ng-1,h)
     
  enddo
  
end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(self,u,q,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,dx,dt,ng)
  class (non_ideal_t):: self
  
  real,dimension(:,:,:,:)::u 
  real,dimension(:,:,:,:)::q 
  
  integer :: ng
  
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ng,1:ng,1:ng)::rhocellmin
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::fluxad
  real(8),dimension(1:ng,1:ng,1:ng,1:3,1:3)::bmagij
  real,dimension(1:ng,1:ng,1:ng,1:3)::emfambdiff
  real,dimension(1:ng,1:ng,1:ng,1:3)::fluxambdiff

! declare local variables
  INTEGER :: m, n, ivar

  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::rhofx,rhofy,rhofz
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::bsquarex,bsquarey,bsquarez,bsquare
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::bsquarexx,bsquareyy,bsquarezz
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::betaad2,betaad
  real(8),dimension(1:ng-2,1:ng-2,1:ng-2)::rhox,rhoy,rhoz,rhocell,bcell,bcellold,tcell,ionisrate
  real(8)::dtlim,Cv

  real(8),dimension(1:ng,1:ng,1:ng,1:3)::florentz
  real(8),dimension(1:ng,1:ng,1:ng)  ::bsquaremax

  REAL(8),dimension(1:ng,1:ng,1:ng)::jxbsquare
  real(8),dimension(1:ng,1:ng,1:ng,3)::jcenter
  real(8),dimension(1:ng,1:ng,1:ng,1:3)::jxb

  real(8) :: dx, dy, dz, dt

  dy = dx ; dz = dx

  dtlim=dt!*coefalfven
  !dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

  jxb=0.0d0
  jxbsquare=0.0d0
  jcenter=0.0d0

  rhox=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(2:ng-1,1:ng-2,2:ng-1,1) + u(2:ng-1,2:ng-1,1:ng-2,1) + u(2:ng-1,1:ng-2,1:ng-2,1))
  rhoy=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(1:ng-2,2:ng-1,2:ng-1,1) + u(2:ng-1,2:ng-1,1:ng-2,1) + u(1:ng-2,2:ng-1,1:ng-2,1))
  rhoz=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1) + u(1:ng-2,2:ng-1,2:ng-1,1) + u(2:ng-1,1:ng-2,2:ng-1,1) + u(1:ng-2,1:ng-2,2:ng-1,1))

  rhofx=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(1:ng-2,2:ng-1,2:ng-1,1))
  rhofy=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(2:ng-1,1:ng-2,2:ng-1,1))
  rhofz=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(2:ng-1,2:ng-1,1:ng-2,1))
              
  rhocellmin(2:ng-1,2:ng-1,2:ng-1)=min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)
  rhocell = u(2:ng-1,2:ng-1,2:ng-1,1)

  ! Compute gas temperature in cgs
  tcell=1.0d0
 
  bsquarex=bemfx(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfx(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfx(2:ng-1,2:ng-1,2:ng-1,3)**2
  bsquarey=bemfy(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfy(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfy(2:ng-1,2:ng-1,2:ng-1,3)**2
  bsquarez=bemfz(2:ng-1,2:ng-1,2:ng-1,1)**2+bemfz(2:ng-1,2:ng-1,2:ng-1,2)**2+bemfz(2:ng-1,2:ng-1,2:ng-1,3)**2

  bsquarexx=bmagij(2:ng-1,2:ng-1,2:ng-1,1,1)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,2,1)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,3,1)**2
  bsquareyy=bmagij(2:ng-1,2:ng-1,2:ng-1,1,2)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,2,2)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,3,2)**2
  bsquarezz=bmagij(2:ng-1,2:ng-1,2:ng-1,1,3)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,2,3)**2+bmagij(2:ng-1,2:ng-1,2:ng-1,3,3)**2

  bsquaremax(2:ng-1,2:ng-1,2:ng-1)=max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)
                 
  ! EMF x
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,1)= florentzx(2:ng-1,2:ng-1,2:ng-1,2)*bemfx(2:ng-1,2:ng-1,2:ng-1,3) - florentzx(2:ng-1,2:ng-1,2:ng-1,3)*bemfx(2:ng-1,2:ng-1,2:ng-1,2)
  rhox=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(2:ng-1,1:ng-2,2:ng-1,1)+u(2:ng-1,2:ng-1,1:ng-2,1)+u(2:ng-1,1:ng-2,1:ng-2,1))
  bcell = bsquaremax(2:ng-1,2:ng-1,2:ng-1)
  bcellold=bcell
  rhocell = rhocellmin(2:ng-1,2:ng-1,2:ng-1)
  ionisrate=1!default_ionisrate
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhox) 
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,1)=emfambdiff(2:ng-1,2:ng-1,2:ng-1,1)*betaad2 

  ! EMF y
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,2)= florentzy(2:ng-1,2:ng-1,2:ng-1,3)*bemfy(2:ng-1,2:ng-1,2:ng-1,1) - florentzy(2:ng-1,2:ng-1,2:ng-1,1)*bemfy(2:ng-1,2:ng-1,2:ng-1,3)
  rhoy=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(1:ng-2,2:ng-1,2:ng-1,1)+u(2:ng-1,2:ng-1,1:ng-2,1)+u(1:ng-2,2:ng-1,1:ng-2,1))
  bcell = bsquaremax(2:ng-1,2:ng-1,2:ng-1)
  bcellold=bcell
  rhocell = rhocellmin(2:ng-1,2:ng-1,2:ng-1)
  ! comparison with hydro+idealMHD 
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoy) 
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,2)=emfambdiff(2:ng-1,2:ng-1,2:ng-1,2)*betaad2            
                    
  ! EMF z
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,3)= florentzz(2:ng-1,2:ng-1,2:ng-1,1)*bemfz(2:ng-1,2:ng-1,2:ng-1,2) - florentzz(2:ng-1,2:ng-1,2:ng-1,2)*bemfz(2:ng-1,2:ng-1,2:ng-1,1)
  rhoz=0.25d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(1:ng-2,2:ng-1,2:ng-1,1)+u(2:ng-1,1:ng-2,2:ng-1,1)+u(1:ng-2,1:ng-2,2:ng-1,1))
  bcell = bsquaremax(2:ng-1,2:ng-1,2:ng-1)
  bcellold=bcell
  rhocell = rhocellmin(2:ng-1,2:ng-1,2:ng-1)     
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoz) 
  emfambdiff(2:ng-1,2:ng-1,2:ng-1,3)=emfambdiff(2:ng-1,2:ng-1,2:ng-1,3)*betaad2

  ! energy flux on faces
  bcell = bsquaremax(2:ng-1,2:ng-1,2:ng-1) 
  rhocell = rhocellmin(2:ng-1,2:ng-1,2:ng-1)
  rhofx=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(1:ng-2,2:ng-1,2:ng-1,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofx) 
  fluxambdiff(2:ng-1,2:ng-1,2:ng-1,1)=-betaad2*fluxad(2:ng-1,2:ng-1,2:ng-1,1)
 
  rhofy=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(2:ng-1,1:ng-2,2:ng-1,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofy) 
  fluxambdiff(2:ng-1,2:ng-1,2:ng-1,2)=-betaad2*fluxad(2:ng-1,2:ng-1,2:ng-1,2)
 
  rhofz=0.5d0*(u(2:ng-1,2:ng-1,2:ng-1,1)+u(2:ng-1,2:ng-1,1:ng-2,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofz) 
  fluxambdiff(2:ng-1,2:ng-1,2:ng-1,3)=-betaad2*fluxad(2:ng-1,2:ng-1,2:ng-1,3)
 
  bcellold=bcell

  
  jcenter(2:ng-1,2:ng-1,2:ng-1,1)= ( bmagij(2:ng-1,3:ng  ,2:ng-1,3,2) - bmagij(2:ng-1,2:ng-1,2:ng-1,3,2) ) / dz - ( bmagij(2:ng-1,2:ng-1,3:ng  ,2,3) - bmagij(2:ng-1,2:ng-1,2:ng-1,3,2) ) / dy
  jcenter(2:ng-1,2:ng-1,2:ng-1,2)= ( bmagij(2:ng-1,2:ng-1,3:ng  ,1,3) - bmagij(2:ng-1,2:ng-1,2:ng-1,1,3) ) / dx - ( bmagij(3:ng  ,2:ng-1,2:ng-1,3,1) - bmagij(2:ng-1,2:ng-1,2:ng-1,3,1) ) / dz
  jcenter(2:ng-1,2:ng-1,2:ng-1,3)= ( bmagij(3:ng  ,2:ng-1,2:ng-1,2,1) - bmagij(2:ng-1,2:ng-1,2:ng-1,2,1) ) / dy - ( bmagij(2:ng-1,3:ng  ,2:ng-1,1,2) - bmagij(2:ng-1,2:ng-1,2:ng-1,1,2) ) / dx
                
  jxb(2:ng-1,2:ng-1,2:ng-1,1) = jcenter(2:ng-1,2:ng-1,2:ng-1,2)*u(2:ng-1,2:ng-1,2:ng-1,8) - jcenter(2:ng-1,2:ng-1,2:ng-1,3)*u(2:ng-1,2:ng-1,2:ng-1,7)
  jxb(2:ng-1,2:ng-1,2:ng-1,2) = jcenter(2:ng-1,2:ng-1,2:ng-1,3)*u(2:ng-1,2:ng-1,2:ng-1,6) - jcenter(2:ng-1,2:ng-1,2:ng-1,1)*u(2:ng-1,2:ng-1,2:ng-1,8)
  jxb(2:ng-1,2:ng-1,2:ng-1,3) = jcenter(2:ng-1,2:ng-1,2:ng-1,1)*u(2:ng-1,2:ng-1,2:ng-1,7) - jcenter(2:ng-1,2:ng-1,2:ng-1,2)*u(2:ng-1,2:ng-1,2:ng-1,6)
 
  betaad=1.d0/(self%gamma_AD*self%rho_ions*u(2:ng-1,2:ng-1,2:ng-1,1)) 
  jxbsquare(2:ng-1,2:ng-1,2:ng-1)=(jxb(2:ng-1,2:ng-1,2:ng-1,1)*jxb(2:ng-1,2:ng-1,2:ng-1,1) +&
                                 & jxb(2:ng-1,2:ng-1,2:ng-1,2)*jxb(2:ng-1,2:ng-1,2:ng-1,2) +&
                                 & jxb(2:ng-1,2:ng-1,2:ng-1,3)*jxb(2:ng-1,2:ng-1,2:ng-1,3))*&
                                 & betaad*dtlim

end SUBROUTINE computambip
! fin modif nimhd

end module
! modif nimhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
