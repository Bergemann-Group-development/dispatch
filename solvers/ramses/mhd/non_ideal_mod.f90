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
    logical:: is_used
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
END SUBROUTINE init

subroutine non_ideal_comp(self,bf,uin,qin,fluxambdiff,emfambdiff,u_max,ds,dt)
  class (non_ideal_t):: self
  integer:: n(4)
  real :: u_max
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
     call computambip(self,uin,qin,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,ds,dt)
     u_max = max(u_max, maxval(uin(:,:,:,6)**2+uin(:,:,:,7)**2+uin(:,:,:,8)**2 / (self%gamma_AD * self%rho_ions * uin(:,:,:,1) * minval(ds)) )  )  
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

  real(8)::tcell,rhocell,bcell,ionisrate,Cv

  real(8), allocatable, dimension(:,:,:)::tcellx,tcelly,tcellz,bsquarex,bsquarey,bsquarez,dtlim
  real(8), allocatable, dimension(:,:,:)::rhox,rhoy,rhoz,epsx,epsy,epsz
  real(8), allocatable, dimension(:,:,:)::rhof,pf,bsqf,epsf,tcellf

  integer, dimension(3) :: index_i,index_j,index_k

  real(8):: ds(3)
  real(8):: dx, dy, dz, dt

  n = shape(u)
  allocate (jcenter (n(1),n(2),n(3),3))
  allocate (jface   (n(1),n(2),n(3),3,3))
  allocate (jemf    (n(1),n(2),n(3),3))
  allocate (tcellx  (n(1)-2,n(2)-2,n(3)-2) )
  allocate (tcelly  (n(1)-2,n(2)-2,n(3)-2) )
  allocate (tcellz  (n(1)-2,n(2)-2,n(3)-2) )
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
  allocate (tcellf  (n(1)-1,n(2)-1,n(3)-1) )

  jcenter  = 0.d0
  jface    = 0.d0
  jemf     = 0.d0
  tcellx   = 0.d0
  tcelly   = 0.d0
  tcellz   = 0.d0
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
  tcellf   = 0.d0

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
  tcellx=1.0d0
  tcelly=1.0d0
  tcellz=1.0d0
  ionisrate=1 !default_ionisrate 
              
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
     
     tcellf=1.0d0
     
     fluxohm(2:n(1)-1,2:n(2)-1,2:n(3)-1,h)=self%etaMD*fluxmd(2:n(1)-1,2:n(2)-1,2:n(3)-1,h)
     
  enddo

  deallocate(jcenter, jface, jemf, tcellx, tcelly, tcellz, bsquarex, bsquarey, bsquarez, dtlim, rhox, rhoy, rhoz, epsx, epsy, epsz, rhof, pf, bsqf, epsf, tcellf)
  
end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(self,u,q,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare,ds,dt)
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
  real(8), allocatable, dimension(:,:,:)::betaad2,betaad
  real(8), allocatable, dimension(:,:,:)::rhox,rhoy,rhoz,rhocell,bcell,bcellold,tcell,ionisrate
  real(8)::dtlim,Cv

  real(8), allocatable, dimension(:,:,:,:)::florentz
  real(8), allocatable, dimension(:,:,:)  ::bsquaremax

  real(8), allocatable, dimension(:,:,:,:)::jcenter
  real(8), allocatable, dimension(:,:,:,:)::jxb

  real(8) :: ds(3)
  real(8) :: dx, dy, dz, dt

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
  allocate (rhox      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhoy      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhoz      (n(1)-2,n(2)-2,n(3)-2))
  allocate (rhocell   (n(1)-2,n(2)-2,n(3)-2))
  allocate (bcell     (n(1)-2,n(2)-2,n(3)-2))
  allocate (bcellold  (n(1)-2,n(2)-2,n(3)-2))
  allocate (tcell     (n(1)-2,n(2)-2,n(3)-2))
  allocate (ionisrate (n(1)-2,n(2)-2,n(3)-2))
  allocate (florentz  (n(1),n(2),n(3),3))
  allocate (bsquaremax(n(1),n(2),n(3)))
  allocate (jcenter   (n(1),n(2),n(3),3))
  allocate (jxb       (n(1),n(2),n(3),3))

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
  rhox       = 0.d0
  rhoy       = 0.d0
  rhoz       = 0.d0
  rhocell    = 0.d0
  bcell      = 0.d0
  bcellold   = 0.d0
  tcell      = 0.d0
  ionisrate  = 0.d0
  florentz   = 0.d0
  bsquaremax = 0.d0
  jcenter    = 0.d0
  jxb        = 0.d0

  dx = ds(1) ; dy = ds(2) ; dz = ds(3)

  dtlim=dt!*coefalfven
  !dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

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

  ! Compute gas temperature in cgs
  tcell=1.0d0
 
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
  ionisrate=1!default_ionisrate
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhox) 
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*betaad2 

  ! EMF y
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) - florentzy(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfy(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
  rhoy=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1)+u(1:n(1)-2,2:n(2)-1,1:n(3)-2,1))
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  bcellold=bcell
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  ! comparison with hydro+idealMHD 
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoy) 
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*betaad2            
                    
  ! EMF z
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) - florentzz(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*bemfz(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
  rhoz=0.25d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1)+u(1:n(1)-2,1:n(2)-2,2:n(3)-1,1))
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  bcellold=bcell
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)     
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhoz) 
  emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=emfambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*betaad2

  ! energy flux on faces
  bcell = bsquaremax(2:n(1)-1,2:n(2)-1,2:n(3)-1) 
  rhocell = rhocellmin(2:n(1)-1,2:n(2)-1,2:n(3)-1)
  rhofx=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(1:n(1)-2,2:n(2)-1,2:n(3)-1,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofx) 
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)
 
  rhofy=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,1:n(2)-2,2:n(3)-1,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofy) 
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)
 
  rhofz=0.5d0*(u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)+u(2:n(1)-1,2:n(2)-1,1:n(3)-2,1))
  betaad2=1.d0/(self%gamma_AD*self%rho_ions*rhofz) 
  fluxambdiff(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)=-betaad2*fluxad(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)
 
  bcellold=bcell

  
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)= ( bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,3,2) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) ) / dz - ( bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,2,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,2) ) / dy
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)= ( bmagij(2:n(1)-1,2:n(2)-1,3:n(3)  ,1,3) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,3) ) / dx - ( bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,3,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,3,1) ) / dz
  jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)= ( bmagij(3:n(1)  ,2:n(2)-1,2:n(3)-1,2,1) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,2,1) ) / dy - ( bmagij(2:n(1)-1,3:n(2)  ,2:n(3)-1,1,2) - bmagij(2:n(1)-1,2:n(2)-1,2:n(3)-1,1,2) ) / dx
                
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,8) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,7)
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,6) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,8)
  jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3) = jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,7) - jcenter(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,6)
 
  betaad=1.d0/(self%gamma_AD*self%rho_ions*u(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)) 
  jxbsquare(2:n(1)-1,2:n(2)-1,2:n(3)-1)=(jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,1) +&
                                 & jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,2) +&
                                 & jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3)*jxb(2:n(1)-1,2:n(2)-1,2:n(3)-1,3))*&
                                 & betaad*dtlim

  deallocate (rhocellmin, rhofx, rhofy, rhofz, bsquarex, bsquarey, bsquarez, bsquarexx, bsquareyy, bsquarezz,  &
              betaad, betaad2, rhox, rhoy, rhoz, rhocell, bcell, bcellold, tcell, ionisrate, jcenter, jxb )


end SUBROUTINE computambip
! fin modif nimhd

end module
! modif nimhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
