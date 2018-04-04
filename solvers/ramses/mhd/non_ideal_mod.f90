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

subroutine non_ideal_comp(self,bf,uin,qin,fluxambdiff,emfambdiff,u_max,dx,dt,ngrid)
  class (non_ideal_t):: self
  integer :: ngrid
  real :: u_max
  real(8) :: dx, dt
  real,dimension(1:ngrid,1:ngrid,1:ngrid,1:8)::uin, qin, fluxambdiff, emfambdiff
  real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bf
  ! Output electromotive force
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::emfx
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::emfy
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::emfz

  ! Output courant vector in the cell
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::bmagij
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::jcentersquare,jxbsquare
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxmd,fluxh,fluxad
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::emfohmdiss,fluxohm 

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
     call computejb2(self,bf,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,dx,ngrid)
  endif
  
  ! AMBIPOLAR DIFFUSION
  
  if(self%mhd_AD) then
     call computambip(self,uin,qin,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,dx,dt,ngrid)
     u_max = max(u_max, maxval(uin(:,:,:,6)**2+uin(:,:,:,7)**2+uin(:,:,:,8)**2 / (self%gamma_AD * self%rho_ions * uin(:,:,:,1) * dx)) )  
  endif
  
  
  ! OHMIC DISSIPATION
  
  if(self%mhd_Ohm) then
     call computdifmag(self,uin,qin,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,dx,dt,ngrid)
     u_max = max(u_max,self%eta_Ohm/ dx)
  endif
  

end subroutine

subroutine flux_upd(self,flux,fluxambdiff, dx, dt, ndim, ngrid)
   class (non_ideal_t):: self
   
   integer :: i, j, k, ndim, ngrid
   real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxohm
   real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxambdiff
   real(8) :: dx, dt

   real, dimension(1:ngrid,1:ngrid,1:ngrid,1:8,1:3):: flux 

   fluxohm=0d0
   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then
     do k=1,ngrid
        do j=1,ngrid
           do i=1,ngrid
              flux(i,j,k,5,ndim) = flux(i,j,k,5,ndim) + (fluxambdiff(i,j,k,ndim)+fluxohm(i,j,k,ndim) ) * dt / dx
          end do
        end do
     end do
   end if
end subroutine


subroutine non_ideal_emf_up(self,emfup,emfambdiff, dx, dt, ndim, ngrid)
   class (non_ideal_t):: self
   integer :: i, j, k,  ibegin, iend, jbegin, jend, kbegin, kend, ndim, ngrid
   real,dimension(1:ngrid,1:ngrid,1:ngrid)::emfup
   real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::emfohmdiss,emfambdiff
   real(8) :: dx, dt

   emfohmdiss = 0d0

   if (self%mhd_AD .eqv. .true. .or. self%mhd_Ohm .eqv. .true.) then

#if NDIM==1

   do k=1,ngrid
      do j=1,ngrid
         do i=1,ngrid
             emfup(i,j,k) = ( emfambdiff(i,j,k,ndim) + emfohmdiss(i,j,k,ndim) ) * dt / dx
         end do 
      end do
   end do
#else
   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1
                 emfup(i,j,k) = emfup(i,j,k) + ( emfambdiff(i,j,k,ndim) + emfohmdiss(i,j,k,ndim) ) * dt / dx
           end do 
        end do
     end do
#endif     
   end if
end subroutine


subroutine computejb2(self,bf,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,dx,ngrid)
  class (non_ideal_t):: self

  integer :: i, j, k,  ibegin, iend, jbegin, jend, kbegin, kend, ndim, ngrid
  real,dimension(:,:,:,:)::bf 
  real,dimension(:,:,:,:)::q 


  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxmd,fluxh,fluxad
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::bmagij
!  real,dimension(1:self%gn(1),1:self%gn(2),1:self%gn(3),1:3)::jcell


  ! declare local variables
  integer :: m, n 

  real(8)::computdx,computdy,computdz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bmagijbis
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::jface
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bcenter
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::fluxbis,fluxter,fluxquat
  real(8)::computdxbis,computdybis,computdzbis
  real(8)::dx,dy,dz

  dy = dx ; dz = dx
   ! Compute before the loop because bmagijbis has to be computed with index i+1/j+1/k+1 in loop over i/j/k
   do k=2,ngrid
      do j=2,ngrid
         do i=2,ngrid
              bmagijbis(i,j,k,1)= 0.25d0* (bf(i  ,j  ,k  ,1) + bf(i  ,j-1,k  ,1) + bf(i  ,j  ,k-1,1) + bf(i  ,j-1,k-1,1) )
              bmagijbis(i,j,k,2)= 0.25d0* (bf(i  ,j  ,k  ,2) + bf(i-1,j  ,k  ,2) + bf(i  ,j  ,k-1,2) + bf(i-1,j  ,k-1,2) ) 
              bmagijbis(i,j,k,3)= 0.25d0* (bf(i  ,j  ,k  ,3) + bf(i-1,j  ,k  ,3) + bf(i  ,j-1,k  ,3) + bf(i-1,j-1,k  ,3) ) 
         end do
      end do
   end do

   ! Define the 1 and 2 components of bcenter, since they are needed further below, when computing the current on the cell interfaces (jface) 
   ! magnetic field at center of cells
   bcenter(1:ngrid,1:ngrid,1:ngrid,1:3)=q(1:ngrid,1:ngrid,1:ngrid,6:8)
   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1
              !!!!!!!!!!!!!!!!!!
              ! EMF x
              !!!!!!!!!!!!!!!!!!
              bemfx(i,j,k,1)=0.25d0* ( q(i  ,j  ,k  ,6) + q(i  ,j-1,k  ,6) + q(i  ,j  ,k-1,6) + q(i  ,j-1,k-1,6) )
              bemfx(i,j,k,2)=0.5d0 * ( bf(i  ,j  ,k  ,2) + bf(i  ,j  ,k-1,2) )
              bemfx(i,j,k,3)=0.5d0 * ( bf(i  ,j  ,k  ,3) + bf(i  ,j-1,k  ,3) )
              !!!!!!!!!!!!!!!!!!
              ! EMF y
              !!!!!!!!!!!!!!!!!!
              bemfy(i,j,k,1)=0.5d0 * ( bf(i  ,j  ,k  ,1) + bf(i  ,j  ,k-1,1) )
              bemfy(i,j,k,2)=0.25d0* ( q (i  ,j  ,k  ,7) + q (i-1,j  ,k  ,7) + q(i  ,j  ,k-1,7) + q(i-1,j  ,k-1,7) )
              bemfy(i,j,k,3)=0.5d0 * ( bf(i-1,j  ,k  ,3) + bf(i  ,j  ,k  ,3) )
              !!!!!!!!!!!!!!!!!!
              ! EMF z
              !!!!!!!!!!!!!!!!!!
              bemfz(i,j,k,1)=0.5d0 * ( bf(i  ,j  ,k  ,1) + bf(i  ,j-1,k  ,1) )
              bemfz(i,j,k,2)=0.5d0 * ( bf(i  ,j  ,k  ,2) + bf(i-1,j  ,k  ,2) )
              bemfz(i,j,k,3)=0.25d0* ( q (i  ,j  ,k  ,8) + q (i-1,j  ,k  ,8) + q(i  ,j-1,k  ,8) + q(i-1,j-1,k  ,8) )
              do m=1,3
              ! bmagij is the value of the magnetic field Bi where Bj 
              ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
              ! and we can write it Bx,y  
              !! m+5 mandatory cf Bx=uin(i,j,k,6)
                 bmagij(i,j,k,m,m)=bf(i,j,k,m)
              end do
              ! case Bx,y
              bmagij(i,j,k,1,2) = 0.5d0 * (q(i  ,j  ,k  ,6) + q(i  ,j-1,k  ,6) )
              ! case Bx,z
              bmagij(i,j,k,1,3) = 0.5d0 * (q(i  ,j  ,k  ,6) + q(i  ,j  ,k-1,6) )
              ! case By,y
              bmagij(i,j,k,2,1) = 0.5d0 * (q(i  ,j  ,k  ,7) + q(i-1,j  ,k  ,7) )
              ! case By,z
              bmagij(i,j,k,2,3) = 0.5d0 * (q(i  ,j  ,k  ,7) + q(i  ,j  ,k-1,7) )
              ! case Bz,x
              bmagij(i,j,k,3,1) = 0.5d0 * (q(i  ,j  ,k  ,8) + q(i-1,j  ,k  ,8) )
              ! case Bz,y
              bmagij(i,j,k,3,2) = 0.5d0 * (q(i  ,j  ,k  ,8) + q(i  ,j-1,k  ,8) )
              jemfx(i,j,k,1) = (bf(i  ,j  ,k  ,3) - bf(i  ,j-1,k  ,3) ) / dy - (bf(i  ,j  ,k  ,2) - bf(i  ,j  ,k-1,2) ) / dz 
              jemfx(i,j,k,2) = (bmagij(i  ,j  ,k  ,1,2) - bmagij(i  ,j  ,k-1,1,2) ) / dz - (bmagijbis(i+1,j  ,k  ,3) - bmagijbis(i  ,j  ,k  ,3) ) / dx
              jemfx(i,j,k,3) = (bmagijbis(i+1,j  ,k  ,2) - bmagijbis(i  ,j  ,k  ,2) ) / dx - (bmagij(i  ,j  ,k  ,1,3)- bmagij(i  ,j-1,k  ,1,3) ) / dy

              jemfy(i,j,k,1) = (bmagijbis(i  ,j+1,k  ,3) - bmagijbis(i  ,j  ,k  ,3) ) / dy - (bmagij(i  ,j  ,k  ,2,1) - bmagij(i  ,j  ,k-1,2,1) ) / dz
              jemfy(i,j,k,2) = (bf(i  ,j  ,k  ,1) - bf(i  ,j  ,k-1,1) ) / dz - (bf(i  ,j  ,k  ,3) - bf(i-1,j  ,k  ,3) ) / dx
              jemfy(i,j,k,3) = (bmagij(i  ,j  ,k  ,2,3) - bmagij(i-1,j  ,k  ,2,3) ) / dx - (bmagijbis(i  ,j+1,k  ,1) - bmagijbis(i  ,j  ,k  ,1) ) / dy

              jemfz(i,j,k,1) = (bmagij(i  ,j  ,k  ,3,1) - bmagij(i  ,j-1,k  ,3,1) ) / dy - (bmagijbis(i  ,j  ,k+1,2) - bmagijbis(i  ,j  ,k  ,2) ) / dz
              jemfz(i,j,k,2) = (bmagijbis(i  ,j  ,k+1,1) - bmagijbis(i  ,j  ,k  ,1) ) / dz - (bmagij(i  ,j  ,k  ,3,2) - bmagij(i-1,j  ,k  ,3,2) ) / dx
              jemfz(i,j,k,3) = (bf(i  ,j  ,k  ,2) - bf(i-1,j  ,k  ,2) ) / dx - (bf(i  ,j  ,k  ,1) - bf(i  ,j-1,k  ,1) ) / dy
         end do
      end do
   end do




if(self%mhd_AD) then

! EMF x
  
  do k=2,ngrid-1
     do j=2,ngrid-1
        do i=2,ngrid-1
              florentzx(i,j,k,1)=jemfx(i,j,k,2)*bemfx(i,j,k,3) - jemfx(i,j,k,3)*bemfx(i,j,k,2)
              florentzx(i,j,k,2)=jemfx(i,j,k,3)*bemfx(i,j,k,1) - jemfx(i,j,k,1)*bemfx(i,j,k,3)
              florentzx(i,j,k,3)=jemfx(i,j,k,1)*bemfx(i,j,k,2) - jemfx(i,j,k,1)*bemfx(i,j,k,2)

              florentzy(i,j,k,1)=jemfy(i,j,k,2)*bemfy(i,j,k,3) - jemfy(i,j,k,3)*bemfy(i,j,k,2)
              florentzy(i,j,k,2)=jemfy(i,j,k,3)*bemfy(i,j,k,1) - jemfy(i,j,k,1)*bemfy(i,j,k,3)
              florentzy(i,j,k,3)=jemfy(i,j,k,1)*bemfy(i,j,k,2) - jemfy(i,j,k,1)*bemfy(i,j,k,2)

              florentzz(i,j,k,1)=jemfz(i,j,k,2)*bemfz(i,j,k,3) - jemfz(i,j,k,3)*bemfz(i,j,k,2)
              florentzz(i,j,k,2)=jemfz(i,j,k,3)*bemfz(i,j,k,1) - jemfz(i,j,k,1)*bemfz(i,j,k,3)
              florentzz(i,j,k,3)=jemfz(i,j,k,1)*bemfz(i,j,k,2) - jemfz(i,j,k,1)*bemfz(i,j,k,2)
        end do
     end do
  end do
  

endif


! computation of current on faces

if((self%mhd_AD).or.(self%mhd_Ohm)) then

! face at i-1/2,j,k


   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1
              jface(i,j,k,1,1)=computdybis(bemfz,3,i,j,k,dy,ngrid)-computdzbis(bemfy,2,i,j,k,dz,ngrid)
              jface(i,j,k,2,1)=computdzbis(bemfy,1,i,j,k,dz,ngrid)-computdxbis(bcenter,3,i-1,j,k,dx,ngrid)
              jface(i,j,k,3,1)=computdxbis(bcenter,2,i-1,j,k,dx,ngrid)-computdybis(bemfz,1,i,j,k,dy,ngrid)
              jface(i,j,k,1,2)=computdybis(bcenter,3,i,j-1,k,dy,ngrid)-computdzbis(bemfx,2,i,j,k,dz,ngrid)
              jface(i,j,k,2,2)=computdzbis(bemfx,1,i,j,k,dz,ngrid)-computdxbis(bemfz,3,i,j,k,dx,ngrid)
              jface(i,j,k,3,2)=computdxbis(bemfz,2,i,j,k,dx,ngrid)-computdybis(bcenter,1,i,j-1,k,dy,ngrid)
              jface(i,j,k,1,3)=computdybis(bemfx,3,i,j,k,dy,ngrid)-computdzbis(bcenter,2,i,j,k-1,dz,ngrid)             
              jface(i,j,k,2,3)=computdzbis(bcenter,1,i,j,k-1,dz,ngrid)-computdxbis(bemfy,3,i,j,k,dx,ngrid)             
              jface(i,j,k,3,3)=computdxbis(bemfy,2,i,j,k,dx,ngrid)-computdybis(bemfx,1,i,j,k,dx,ngrid)            

              do n=1,3
                fluxbis(i,j,k,1,n) = jface(i,j,k,2,n)*bmagij(i,j,k,3,n) - jface(i,j,k,3,n)*bmagij(i,j,k,2,n)
                fluxbis(i,j,k,2,n) = jface(i,j,k,3,n)*bmagij(i,j,k,1,n) - jface(i,j,k,1,n)*bmagij(i,j,k,3,n)
                fluxbis(i,j,k,3,n) = jface(i,j,k,1,n)*bmagij(i,j,k,2,n) - jface(i,j,k,2,n)*bmagij(i,j,k,1,n)
              end do

              fluxmd(i,j,k,1)=fluxbis(i,j,k,1,1)
              fluxmd(i,j,k,2)=fluxbis(i,j,k,2,2)
              fluxmd(i,j,k,3)=fluxbis(i,j,k,3,3)
        end do
     end do
  end do


endif

if(self%mhd_AD) then

   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1

              do n=1,3
                 fluxter(i,j,k,1,n)=fluxbis(i,j,k,2,n)*bmagij(i,j,k,3,n) - fluxbis(i,j,k,3,n)*bmagij(i,j,k,2,n)
                 fluxter(i,j,k,2,n)=fluxbis(i,j,k,3,n)*bmagij(i,j,k,1,n) - fluxbis(i,j,k,1,n)*bmagij(i,j,k,3,n)
                 fluxter(i,j,k,3,n)=fluxbis(i,j,k,1,n)*bmagij(i,j,k,2,n) - fluxbis(i,j,k,2,n)*bmagij(i,j,k,1,n)
              end do

              fluxh(i,j,k,1)=fluxter(i,j,k,1,1)
              fluxh(i,j,k,2)=fluxter(i,j,k,2,2)
              fluxh(i,j,k,3)=fluxter(i,j,k,3,3)

              do n=1,3                 
                 fluxquat(i,j,k,1,n)=fluxter(i,j,k,2,n)*bmagij(i,j,k,3,n) - fluxter(i,j,k,3,n)*bmagij(i,j,k,2,n)
                 fluxquat(i,j,k,2,n)=fluxter(i,j,k,3,n)*bmagij(i,j,k,1,n) - fluxter(i,j,k,1,n)*bmagij(i,j,k,3,n)
                 fluxquat(i,j,k,3,n)=fluxter(i,j,k,1,n)*bmagij(i,j,k,2,n) - fluxter(i,j,k,2,n)*bmagij(i,j,k,1,n)
              end do
              fluxad(i,j,k,1)=fluxquat(i,j,k,1,1)
              fluxad(i,j,k,2)=fluxquat(i,j,k,2,2)
              fluxad(i,j,k,3)=fluxquat(i,j,k,3,3)

        end do
     end do
  end do

endif


end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(self,u,q,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare,dx,dt,ngrid)
  class (non_ideal_t):: self

  integer :: i, j, k, ngrid
  real,dimension(:,:,:,:)::u 
  real,dimension(:,:,:,:)::q 

  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::jemfx,jemfy,jemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::bmagij
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxmd

real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3):: emfohmdiss,fluxohm 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::jcentersquare

! declare local variables
  integer :: m, n, h


! WARNING following quantities defined with three components even
! if ndim<3 !
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,3)::jcenter
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::jface
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::jemf


real(8)::computdx,computdy,computdz,tcell,rhocell,bcell,ionisrate,Cv

real(8)::etaod2,tcellx,tcelly,tcellz,bsquarex,bsquarey,bsquarez,dtlim
real(8)::pressurex,pressurey,pressurez,rhox,rhoy,rhoz,epsx,epsy,epsz
real(8)::etaod2x,etaod2y,etaod2z,rhof,pf,bsqf,epsf,tcellf

integer , dimension(1:3) :: index_i,index_j,index_k

real(8):: dx, dt

index_i = (/1,0,0/)
index_j = (/0,1,0/)
index_k = (/0,0,1/)

dtlim = dt !neil

   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1
              jemf(i,j,k,1)=jemfx(i,j,k,1)
              jemf(i,j,k,2)=jemfy(i,j,k,2)
              jemf(i,j,k,3)=jemfz(i,j,k,3)


              rhox=0.25d0*(u(i,j,k,   1)+u(i  ,j-1,k,   1)+u(i,j  ,k-1,   1)+u(i  ,j-1,k-1,   1))
              rhoy=0.25d0*(u(i,j,k,   1)+u(i-1,j  ,k,   1)+u(i,j  ,k-1,   1)+u(i-1,j  ,k-1,   1))
              rhoz=0.25d0*(u(i,j,k,   1)+u(i-1,j  ,k,   1)+u(i,j-1,k  ,   1)+u(i-1,j-1,k  ,   1))
              epsx=0.25d0*(u(i,j,k,nvar)+u(i  ,j-1,k,nvar)+u(i,j  ,k-1,nvar)+u(i  ,j-1,k-1,nvar))
              epsy=0.25d0*(u(i,j,k,nvar)+u(i-1,j  ,k,nvar)+u(i,j  ,k-1,nvar)+u(i-1,j  ,k-1,nvar))
              epsz=0.25d0*(u(i,j,k,nvar)+u(i-1,j  ,k,nvar)+u(i,j-1,k  ,nvar)+u(i-1,j-1,k  ,nvar))
               if(self%mhd_Ohm)then
                 bsquarex=bemfx(i,j,k,1)**2+bemfx(i,j,k,2)**2+bemfx(i,j,k,3)**2
                 bsquarey=bemfy(i,j,k,1)**2+bemfy(i,j,k,2)**2+bemfy(i,j,k,3)**2
                 bsquarez=bemfz(i,j,k,1)**2+bemfz(i,j,k,2)**2+bemfz(i,j,k,3)**2
               end if

              tcellx=1.0d0
              tcelly=1.0d0
              tcellz=1.0d0

              ionisrate=1 !default_ionisrate 
              etaod2x=self%etaMD
              etaod2y=self%etaMD
              etaod2z=self%etaMD
              
! WARNING dB/dt=-curl(eta*J)
              emfohmdiss(i,j,k,1)=-etaod2x*jemf(i,j,k,1)
              emfohmdiss(i,j,k,2)=-etaod2y*jemf(i,j,k,2)
              emfohmdiss(i,j,k,3)=-etaod2z*jemf(i,j,k,3)

! !!!!!!!!!!!!!!!!!!!!!!!
! !
! ! compute j at center of cells
! !
! ! mandatory for non isotherm case
! 
! ! bmagij is the value of the magnetic field Bi where Bj 
! ! is naturally defined; Ex bmagij(i,j,k,1,2) is Bx at i,j-1/2,k
! ! and we can write it Bx,y

              jcenter(i,j,k,1)=computdy(bmagij,3,2,i,j,k,dx)-computdz(bmagij,2,3,i,j,k,dx)
              jcenter(i,j,k,2)=computdz(bmagij,1,3,i,j,k,dx)-computdx(bmagij,3,1,i,j,k,dx)
              jcenter(i,j,k,3)=computdx(bmagij,2,1,i,j,k,dx)-computdy(bmagij,1,2,i,j,k,dx)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
              do h = 1,3
                 
                 rhof=0.5d0*(u(i,j,k,   1)+u(i-index_i(h),j-index_j(h),k-index_k(h),   1))
                 epsf=0.5d0*(u(i,j,k,nvar)+u(i-index_i(h),j-index_j(h),k-index_k(h),nvar))
                 bsqf=bmagij(i,j,k,1,h)**2+bmagij(i,j,k,2,h)**2+bmagij(i,j,k,3,h)**2
                 
                 tcellf=1.0d0
                 
                 etaod2=self%etaMD
                 fluxohm(i,j,k,h)=etaod2*fluxmd(i,j,k,h)
                 
              enddo

        end do
     end do
  end do
  

end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(self,u,q,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,dx,dt,ngrid)
  class (non_ideal_t):: self
  
  real,dimension(:,:,:,:)::u 
  real,dimension(:,:,:,:)::q 
  
  integer :: i, j, k, ngrid
  
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::bemfx,bemfy,bemfz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)::rhocellmin
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::florentzx,florentzy,florentzz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxad
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::bmagij
  real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::emfambdiff
  real,dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::fluxambdiff

! declare local variables
  INTEGER :: m, n, ivar
  real(8)::computdx,computdy,computdz

  real(8)::rhofx,rhofy,rhofz
  real(8)::bsquarex,bsquarey,bsquarez,bsquare
  real(8)::bsquarexx,bsquareyy,bsquarezz
  real(8)::betaad2,betaad
  real(8)::rhox,rhoy,rhoz,rhocell,bcell,bcellold,tcell,ionisrate
  real(8)::dtlim,Cv,eps

  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::florentz
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid)  ::bsquaremax

  REAL(8),dimension(1:ngrid,1:ngrid,1:ngrid)::jxbsquare
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,3)::jcenter
  real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::jxb

  real(8) :: dx, dy, dz, dt

  dy = dx ; dz = dx

  dtlim=dt!*coefalfven
!dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

jxb=0.0d0

jxbsquare=0.0d0
jcenter=0.0d0

   do k=2,ngrid-1
      do j=2,ngrid-1
         do i=2,ngrid-1
              rhox=0.25d0*(u(i,j,k,1)+u(i,j-1,k,1)+u(i,j,k-1,1)+u(i,j-1,k-1,1))
              rhoy=0.25d0*(u(i,j,k,1)+u(i-1,j,k,1)+u(i,j,k-1,1)+u(i-1,j,k-1,1))
              rhoz=0.25d0*(u(i,j,k,1)+u(i-1,j,k,1)+u(i,j-1,k,1)+u(i-1,j-1,k,1))

              rhofx=0.5d0*(u(i,j,k,1)+u(i-1,j,k,1))
              rhofy=0.5d0*(u(i,j,k,1)+u(i,j-1,k,1))
              rhofz=0.5d0*(u(i,j,k,1)+u(i,j,k-1,1))
              
              rhocellmin(i,j,k)=min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)

              rhocell = u(i,j,k,1)

              ! Compute gas temperature in cgs
              tcell=1.0d0
             

              bsquarex=bemfx(i,j,k,1)**2+bemfx(i,j,k,2)**2+bemfx(i,j,k,3)**2
              bsquarey=bemfy(i,j,k,1)**2+bemfy(i,j,k,2)**2+bemfy(i,j,k,3)**2
              bsquarez=bemfz(i,j,k,1)**2+bemfz(i,j,k,2)**2+bemfz(i,j,k,3)**2


              bsquarexx=bmagij(i,j,k,1,1)**2+bmagij(i,j,k,2,1)**2+bmagij(i,j,k,3,1)**2
              bsquareyy=bmagij(i,j,k,1,2)**2+bmagij(i,j,k,2,2)**2+bmagij(i,j,k,3,2)**2
              bsquarezz=bmagij(i,j,k,1,3)**2+bmagij(i,j,k,2,3)**2+bmagij(i,j,k,3,3)**2

              bsquaremax(i,j,k)=max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)
                 
! EMF x
              emfambdiff(i,j,k,1)= florentzx(i,j,k,2)*bemfx(i,j,k,3) - florentzx(i,j,k,3)*bemfx(i,j,k,2)
              rhox=0.25d0*(u(i,j,k,1)+u(i,j-1,k,1)+u(i,j,k-1,1)+u(i,j-1,k-1,1))
              bcell = bsquaremax(i,j,k)
              bcellold=bcell

              rhocell = rhocellmin(i,j,k)

              ionisrate=1!default_ionisrate
              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhox) 
              emfambdiff(i,j,k,1)=emfambdiff(i,j,k,1)*betaad2 

! EMF y
              emfambdiff(i,j,k,2)= florentzy(i,j,k,3)*bemfy(i,j,k,1) - florentzy(i,j,k,1)*bemfy(i,j,k,3)
              rhoy=0.25d0*(u(i,j,k,1)+u(i-1,j,k,1)+u(i,j,k-1,1)+u(i-1,j,k-1,1))
              bcell = bsquaremax(i,j,k)
              bcellold=bcell

              rhocell = rhocellmin(i,j,k)
! comparison with hydro+idealMHD 

              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoy) 
              emfambdiff(i,j,k,2)=emfambdiff(i,j,k,2)*betaad2            
                    
! EMF z
              emfambdiff(i,j,k,3)= florentzz(i,j,k,1)*bemfz(i,j,k,2) - florentzz(i,j,k,2)*bemfz(i,j,k,1)
              rhoz=0.25d0*(u(i,j,k,1)+u(i-1,j,k,1)+u(i,j-1,k,1)+u(i-1,j-1,k,1))
              bcell = bsquaremax(i,j,k)
              bcellold=bcell
              rhocell = rhocellmin(i,j,k)
             
              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoz) 

              emfambdiff(i,j,k,3)=emfambdiff(i,j,k,3)*betaad2

! energy flux on faces

              bcell = bsquaremax(i,j,k)


              rhocell = rhocellmin(i,j,k)

              rhofx=0.5d0*(u(i,j,k,1)+u(i-1,j,k,1))
              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofx) 
              fluxambdiff(i,j,k,1)=-betaad2*fluxad(i,j,k,1)

              rhofy=0.5d0*(u(i,j,k,1)+u(i,j-1,k,1))
              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofy) 
              fluxambdiff(i,j,k,2)=-betaad2*fluxad(i,j,k,2)

              rhofz=0.5d0*(u(i,j,k,1)+u(i,j,k-1,1))
              betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofz) 
              fluxambdiff(i,j,k,3)=-betaad2*fluxad(i,j,k,3)

              bcellold=bcell

              jcenter(i,j,k,1)=computdy(bmagij,3,2,i,j,k,dx,ngrid)-computdz(bmagij,2,3,i,j,k,dx,ngrid)
              jcenter(i,j,k,2)=computdz(bmagij,1,3,i,j,k,dy,ngrid)-computdx(bmagij,3,1,i,j,k,dy,ngrid)
              jcenter(i,j,k,3)=computdx(bmagij,2,1,i,j,k,dz,ngrid)-computdy(bmagij,1,2,i,j,k,dz,ngrid)
              

              jxb(i,j,k,1) = jcenter(i,j,k,2)*u(i,j,k,8) - jcenter(i,j,k,3)*u(i,j,k,7)
              jxb(i,j,k,2) = jcenter(i,j,k,3)*u(i,j,k,6) - jcenter(i,j,k,1)*u(i,j,k,8)
              jxb(i,j,k,3) = jcenter(i,j,k,1)*u(i,j,k,2) - jcenter(i,j,k,1)*u(i,j,k,2)

              betaad=1.d0/(self%gamma_AD*self%rho_ions*u(i,j,k,1))

              jxbsquare(i,j,k)=(jxb(i,j,k,1)*jxb(i,j,k,1)+jxb(i,j,k,2)*jxb(i,j,k,2)+jxb(i,j,k,3)*jxb(i,j,k,3))*&
              & betaad*dtlim
        end do
     end do
  end do

end SUBROUTINE computambip
! fin modif nimhd

end module
! modif nimhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdx(vec,n2,n3,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::vec
real(8)::dx
integer::n2,n3,i,j,k, ngrid

computdx=(vec(i+1,j,k,n2,n3)-vec(i,j,k,n2,n3))/dx

end function computdx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdy(vec,n2,n3,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::vec
real(8)::dx
integer::n2,n3,i,j,k, ngrid

computdy=(vec(i,j+1,k,n2,n3)-vec(i,j,k,n2,n3))/dx

end function computdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdz(vec,n2,n3,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3,1:3)::vec
real(8)::dx
integer::n2,n3,i,j,k, ngrid

computdz=(vec(i,j,k+1,n2,n3)-vec(i,j,k,n2,n3))/dx

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdxbis(vec,n2,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
real(8)::dx
integer::n2,i,j,k, ngrid

computdxbis=(vec(i+1,j,k,n2)-vec(i,j,k,n2))/dx

end function computdxbis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdybis(vec,n2,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
real(8)::dx
integer::n2,i,j,k, ngrid

computdybis=(vec(i,j+1,k,n2)-vec(i,j,k,n2))/dx

end function computdybis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdzbis(vec,n2,i,j,k,dx,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none 
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
real(8)::dx
integer::n2,i,j,k, ngrid

computdzbis=(vec(i,j,k+1,n2)-vec(i,j,k,n2))/dx

end function computdzbis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdxvx(vec,i,j,k,dx,dy,dz,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
integer i,j,k, ngrid
real(8)::dx,dy,dz

computdxvx=(vec(i+1,j,k,1)-vec(i,j,k,1))/dx 


end function computdxvx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdyvy(vec,i,j,k,dx,dy,dz,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
integer i,j,k, ngrid
real(8)::dx,dy,dz

computdyvy=(vec(i,j+1,k,2)-vec(i,j,k,2))/dy 

end function computdyvy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdzvz(vec,i,j,k,dx,dy,dz,ngrid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none
real(8),dimension(1:ngrid,1:ngrid,1:ngrid,1:3)::vec
integer i,j,k, ngrid
real(8)::dx,dy,dz

computdzvz=(vec(i,j,k+1,3)-vec(i,j,k,3))/dz

end function computdzvz

! fin modif nimhd
