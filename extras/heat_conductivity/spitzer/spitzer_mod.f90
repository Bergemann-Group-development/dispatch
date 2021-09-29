!=======================================================================
!> 
!=======================================================================
MODULE spitzer_mod
  USE aux_mod
  USE bits_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE io_unit_mod
  USE scalar_mod
  USE stagger_mod
  USE mesh_mod
  USE index_mod
  USE eos_mod
  USE vector_mod
  USE vector_ops
  USE scaling_mod
  USE gpatch_mod
  USE link_mod
  USE scratch_mod
  implicit none
  private
  logical, save  :: on = .false.
  logical, save  :: first_time = .true.
  character(len=8), save :: method='explicit'
  real(kind=KindScalarVar), save            :: sigma_e= 1.1e-6                  ! in CGS /ul/ul/(ue/ut)
  real(kind=KindScalarVar), parameter       :: csat   = 6.90256e12              ! in CGS ur/(ue/ut)
  real(kind=KindScalarVar), save            :: amp = 1.0                        ! spitzer amplification
  real(kind=KindScalarVar), save            :: theta_mg = 0.5                   ! should be ~0.5  for wave, how fast the wave propagates wrt alfvenic spd
  real(kind=KindScalarVar), save            :: bmin = 1e-10
  type, public:: spitzer_t
    real(8) :: sigma0, csat
    real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: bb     
    real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: bm     
    real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: fc                 ! conductive flux
    real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: dfc                ! conductive flux change
    real(kind=KindScalarVar), dimension(:,:,:),   pointer :: sigmac   
    real                                                  :: u_max = 0.0
    class(gpatch_t), pointer :: patch 
    logical:: on = .false.   
    logical:: first_time = .true.
    logical:: first_iter = .true.
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: compute_q
    procedure:: courant_condition
    procedure:: pre_update
    procedure:: post_update
    procedure:: pre_update_wave
    procedure:: post_update_explicit
    procedure:: compute_b_field_topology
    procedure:: conductive_flux
  end type
  type (spitzer_t), public:: spitzer_conductivity
CONTAINS

!=======================================================================
!> Read parameters
!=======================================================================
SUBROUTINE pre_init (self, gpatch)
  class (spitzer_t):: self
  class (gpatch_t), target :: gpatch
  ! --------------------------------------------------------------------
  integer:: iostat
  namelist /spitzer_conductivity_params/ on, amp, theta_mg, bmin, sigma_e, method
  !.....................................................................
  call io%lock%set ('spitzer_t%init')
  if (first_time) then
    rewind (io_unit%input)
    read (io_unit%input, spitzer_conductivity_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io_unit%input, spitzer_conductivity_params)
    end if
    write (io%namelist, spitzer_conductivity_params)
    first_time = .false.
  end if
  call io%lock%unset('spitzer_t%init')
  if (on) then
    self%patch => gpatch
    self%sigma0 = sigma_e / scaling%l/scaling%l/(scaling%e/scaling%t) * amp
    self%csat = csat * scaling%d/(scaling%e/scaling%t)
    self%on = on
    !---------------------------------------------------------------------------
    ! Request heating per unit volume for MHD
    !---------------------------------------------------------------------------
    
  end if
END SUBROUTINE pre_init

!=======================================================================
!> Allocate arrays
!=======================================================================
SUBROUTINE init (self)
  class (spitzer_t):: self
  ! --------------------------------------------------------------------
  integer:: iostat, gn(3)
  if (self%on) then
    gn = self%patch%mesh%gn
    call allocate_vectors(gn,self%fc, self%dfc)
    self%fc  = 0.0
    self%dfc = 0.0    
  end if
END SUBROUTINE init

!=======================================================================
!> Compute the initial conductive heating rate
!=======================================================================
SUBROUTINE compute_q(self, ff)
  class (spitzer_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), target :: ff
  ! --------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: d, e
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: tg, lntg
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: dtx, dty, dtz
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: bbx, bby, bbz
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: bmx, bmy, bmz
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: fcx, fcy, fcz
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: dtg
  real(kind=KindScalarVar) :: kappac, fsat, tmp
  real(8) :: ds(3)
  integer :: n(3), lb(3), ub(3), i, j, k
  ! ....................................................................
  if (.not.self%on) return
  do i=1,3
    lb(i) = self%patch%mesh(i)%lb
    ub(i) = self%patch%mesh(i)%ub
    ds(i) = self%patch%mesh(i)%d
  end do  
  d => ff(:,:,:,self%patch%idx%d)
  e => ff(:,:,:,self%patch%idx%e)
  n =  shape(d)
  ! allocate scratches
  call allocate_scalars (n, lntg, tg)
  call allocate_vectors (n, dtg, self%bb, self%bm)
  dtx => dtg(:,:,:,1)
  dty => dtg(:,:,:,2)
  dtz => dtg(:,:,:,3)
  fcx => self%fc(:,:,:,1)
  fcy => self%fc(:,:,:,2)
  fcz => self%fc(:,:,:,3)
  bbx => self%bb(:,:,:,1); bmx => self%bm(:,:,:,1)                      ! convenience
  bby => self%bb(:,:,:,2); bmy => self%bm(:,:,:,2)                      ! convenience
  bbz => self%bb(:,:,:,3); bmz => self%bm(:,:,:,3)                      ! convenience
  ! magnetic field topology
  call self%compute_b_field_topology(ff)
  ! gas temperature
  call eos%lookup_table (n, d=d, e=e, tt=tg)
  lntg = alog(tg)
  ! centered temperature gradient
  dtg=up(grad(self%patch%mesh,lntg))
  ! Spitzer wave propagation speed
    do k=lb(3),ub(3)
      do j=lb(2),ub(2)
        do i=lb(1),ub(1)
          kappac  = self%sigma0         *exp(3.5*lntg(i,j,k))
          fsat    = self%csat  *d(i,j,k)*exp(1.5*lntg(i,j,k))
          tmp = bbx(i,j,k)*dtx(i,j,k)+ &
                bby(i,j,k)*dty(i,j,k)+ &
                bbz(i,j,k)*dtz(i,j,k)
          fcx(i,j,k)=-kappac*(bbx(i,j,k)*tmp+bmx(i,j,k)*dtx(i,j,k))
          fcy(i,j,k)=-kappac*(bby(i,j,k)*tmp+bmy(i,j,k)*dty(i,j,k))
          fcz(i,j,k)=-kappac*(bbz(i,j,k)*tmp+bmz(i,j,k)*dtz(i,j,k))
        
          fcx(i,j,k)=fcx(i,j,k)/(1.+abs(fcx(i,j,k))/fsat)
          fcy(i,j,k)=fcy(i,j,k)/(1.+abs(fcy(i,j,k))/fsat)
          fcz(i,j,k)=fcz(i,j,k)/(1.+abs(fcz(i,j,k))/fsat) 
        end do
      end do
    end do
    call deallocate_scalars (lntg, tg)
    call deallocate_vectors (dtg, self%bb, self%bm)
END SUBROUTINE compute_q
  

!=======================================================================
!> Update Spitzer conductivity
!=======================================================================
SUBROUTINE pre_update(self)
  class (spitzer_t):: self
  ! --------------------------------------------------------------------
  if (.not.self%on) return
  !
  select case(trim(method))
    case('explicit')
    case ('wave')
      call self%pre_update_wave
    case default
      call io%abort("spitzer conductivity: method not chosen")
    end select
END SUBROUTINE pre_update


!=======================================================================
!> Update Spitzer conductivity.
!=======================================================================
SUBROUTINE post_update(self)
  class (spitzer_t):: self
  ! --------------------------------------------------------------------
  if (.not.self%on) return
  !
  select case(trim(method))
    case('explicit')
      call self%post_update_explicit
    case ('wave')
    case default
      call io%abort("spitzer conductivity: method not chosen")
    end select
END SUBROUTINE post_update

!=======================================================================
!> Update Spitzer conductivity -- wave formulation
!=======================================================================
SUBROUTINE pre_update_wave(self)
  class (spitzer_t):: self
  ! --------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: f
  real(8) :: ds(3), dtime
  ! ....................................................................
  if (.not.self%on) return
  !
  ds = self%patch%ds  
  !
  if (self%patch%dtime <= 1d-20) then 
    dtime = 1d-10
  else
    dtime = self%patch%dtime
  end if
  if (self%first_iter) then
    self%first_iter = .false.
    call self%compute_q(self%patch%mem(:,:,:,:,self%patch%it,1))
  else
    ! save conductive flux values from previous timestep
    self%fc(:,:,:,1) = self%fc(:,:,:,1)+self%dfc(:,:,:,1)*dtime
    self%fc(:,:,:,2) = self%fc(:,:,:,2)+self%dfc(:,:,:,2)*dtime
    self%fc(:,:,:,3) = self%fc(:,:,:,3)+self%dfc(:,:,:,3)*dtime
    ! Calculate conductive flux for this timestep
    f => self%patch%mem(:,:,:,:,self%patch%it,1)
    call self%conductive_flux(f,dtime)
  end if
  call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  scratch%heating_per_unit_volume = & 
  scratch%heating_per_unit_volume - div(ds,down(self%fc))

  !print *, 'max dedt before spitz', maxval(abs(ff(:,:,:,self%idx%e,it,2)))
  !print *, 'max spitz', maxval(abs(div(ds,down(self%fc))))
  !print *, 'max dspitz', maxval(abs(self%dfc))
  ! field aligned flux divergence
  
  !print *, 'max dedt after spitz', maxval(abs(ff(:,:,:,self%idx%e,it,2)))
  
END SUBROUTINE pre_update_wave

!===============================================================================
!> Update Spitzer conductivity -- explicit formulation
!===============================================================================
SUBROUTINE post_update_explicit(self)
  class (spitzer_t):: self
  ! ----------------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: dtg, B, Bc
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d, e, tg, lnd, lntg, dedt
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dtgx, dtgy, dtgz, dddt
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bcx, bcy, bcz, qspitz
  real(kind=KindScalarVar) ::  u_cond,lgtmin_loc,kappac,kappam,tmp
  real(8) :: ds(3), dsmax
  integer :: n(3), lb(3), ub(3), li(3), ui(3)
  
  real,save :: lgtmin = 0.0
  integer :: count,i,j,k
  ! ............................................................................
  do i=1,3
    lb(i) = self%patch%mesh(i)%lb
    ub(i) = self%patch%mesh(i)%ub
    li(i) = self%patch%mesh(i)%li
    ui(i) = self%patch%mesh(i)%ui
    ds(i) = self%patch%mesh(i)%d
     n(i) = self%patch%gn(i)
  end do
  dsmax = maxval(ds,n>1)
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  dedt => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,2)
  dddt => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,2)
  B => self%patch%mem(:,:,:,self%patch%idx%bx:self%patch%idx%bz,self%patch%it,1)
  n =  shape(d)
  ! allocate scratches
  call allocate_scalars (n, lntg, tg, qspitz)
  call allocate_vectors (n, dtg, Bc)
  ! gas temperature
  call eos%lookup_table (n, d=d, e=e, tt=tg)
  lntg = alog(tg)

  u_cond=8e-7*scaling%temp**3.5/scaling%l**2*(scaling%t/scaling%e)

  if (lgtmin == 0.0 ) then
    lgtmin=alog(8.0e3)
  endif

! ** face centered temperature
! scr7 - log temperature
! scr1 - dtx
! scr2 - dty
! scr3 - dtz
! scr4 - bxc
! scr5 - byc
! scr6 - bzc
  ! ----------------------------------------------------------------------------
  ! centered temperature gradient and magnetic field components
  ! ----------------------------------------------------------------------------
  dtg = up(grad(self%patch%mesh,lntg))
  Bc  = up(B)
  ! ----------------------------------------------------------------------------
  ! point for convenience
  ! ----------------------------------------------------------------------------
  dtgx => dtg(:,:,:,1)
  dtgy => dtg(:,:,:,2)
  dtgz => dtg(:,:,:,3)
  bcx  => bc (:,:,:,1)
  bcy  => bc (:,:,:,2)
  bcz  => bc (:,:,:,3)

  count=0
  do k=lb(3),ub(3)
    do j=lb(2),ub(2)
      do i=lb(1),ub(1)
        kappac  = u_cond*exp(3.5*lntg(i,j,k))
        kappam  = self%patch%u_max*e(i,j,k)
        tmp = min(kappac,kappam)*(dtgx(i,j,k)*bcx(i,j,k) + &
                                  dtgy(i,j,k)*bcy(i,j,k) + &
                                  dtgz(i,j,k)*bcz(i,j,k))/ &
                                 (bcx(i,j,k)*bcx(i,j,k) + &
                                  bcy(i,j,k)*bcy(i,j,k) + &
                                  bcz(i,j,k)*bcz(i,j,k) + 1e-20)
        bcx(i,j,k) = bcx(i,j,k)*tmp
        bcy(i,j,k) = bcy(i,j,k)*tmp
        bcz(i,j,k) = bcz(i,j,k)*tmp
        if ((kappac >= kappam).and. &
            (i >= li(1)).and.(i <= ui(1)).and. &
            (j >= li(2)).and.(j <= ui(2)).and. &
            (k >= li(3)).and.(k <= ui(3))) count=count+1
      end do
    end do
  end do

!call extrap_fullz_upper(scr1)
!call extrap_fullz_lower(scr2)

  ! ----------------------------------------------------------------------------
  ! divergence of a field aligned flux
  ! ----------------------------------------------------------------------------
  qspitz = div(ds,down(bc))


  do k=lb(3),ub(3)
    do j=lb(2),ub(2)
      do i=lb(1),ub(1)
        if (lntg(i,j,k) >= lgtmin) dedt(i,j,k)= dedt(i,j,k) + qspitz(i,j,k)
      end do
    end do
  end do

 ! print *, 'spitzer', self%patch%id, minval(qspitz), maxval(qspitz)
  self%u_max = 0.5*dsmax*self%patch%fmaxval(abs((dedt-e/d*dddt)/d))                   ! ops: 1 1 1
  call self%patch%aux%register ('qspitz', qspitz)
  call deallocate_scalars (lntg, tg, qspitz)
  call deallocate_vectors (dtg, Bc)

!print*,'number of points with lowered conduction: ',self%patch%id,(1.0*count)/(1.0*product(n))
  
END SUBROUTINE post_update_explicit

!=======================================================================
!>
!=======================================================================
SUBROUTINE conductive_flux(self, ff, dt)
  class (spitzer_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: ff
  real(8) :: dt
  ! --------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: dtg
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d, e, tg, lnd, lntg
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dtx, dty, dtz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dfcx, dfcy, dfcz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bbx, bby, bbz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bmx, bmy, bmz
  real(kind=KindScalarVar), dimension(:,:,:),   pointer ::   ux,   uy,   uz
  real(kind=KindScalarVar) :: fsat, tmp
  real(8) :: tau_spitz, u_cmax, dtheta, sigmac
  real(8) :: ds(3)
  integer :: n(3), lb(3), ub(3)
  integer :: i, j, k
  ! ....................................................................
  if (.not. self%on) return
    lb = self%patch%mesh%lb
    ub = self%patch%mesh%ub
    ds = self%patch%mesh%d
  d => ff(:,:,:,self%patch%idx%d)
  e => ff(:,:,:,self%patch%idx%e)
  n =  shape(d)
  call allocate_vectors (n, self%bb, self%bm, dtg)
  call allocate_scalars (n,ux, uy, uz)
  call allocate_scalars (n, lnd, self%sigmac, lntg, tg)
  ! gas temperature
  call eos%lookup_table (n, d=d, e=e, tt=tg)
  lntg = log(tg)
  ! gas velocities
  lnd = log(d)
  ux = ff(:,:,:,self%patch%idx%px)/exp(xdn(lnd))
  uy = ff(:,:,:,self%patch%idx%py)/exp(ydn(lnd))
  uz = ff(:,:,:,self%patch%idx%pz)/exp(zdn(lnd))
  dfcx => self%dfc(:,:,:,1)
  dfcy => self%dfc(:,:,:,2)
  dfcz => self%dfc(:,:,:,3)
  dtx => dtg(:,:,:,1)
  dty => dtg(:,:,:,2)
  dtz => dtg(:,:,:,3)
  bbx => self%bb(:,:,:,1); bmx => self%bm(:,:,:,1)                      ! convenience
  bby => self%bb(:,:,:,2); bmy => self%bm(:,:,:,2)                      ! convenience
  bbz => self%bb(:,:,:,3); bmz => self%bm(:,:,:,3)                      ! convenience
  ! gas velocities
  lnd = log(d)
  ux = ff(:,:,:,self%patch%idx%px)/exp(xdn(lnd))
  uy = ff(:,:,:,self%patch%idx%py)/exp(ydn(lnd))
  uz = ff(:,:,:,self%patch%idx%pz)/exp(zdn(lnd))

  ! magnetic field topology
  call self%compute_b_field_topology(ff)
  ! centered temperature gradient
  dtg=up(grad(self%patch%mesh,tg))
  ! Spitzer wave propagation speed
    dtheta = theta_mg*minval(ds)/dt
  !
  do k=lb(3),ub(3)
    do j=lb(2),ub(2)
      do i=lb(1),ub(1)
        sigmac  = self%sigma0*exp(3.5*tg(i,j,k))
        u_cmax = dtheta - max(abs(ux(i,j,k)),abs(uy(i,j,k)),abs(uz(i,j,k)))
        tau_spitz = max(dt, sigmac/e(i,j,k)/u_cmax**2) 
        fsat = self%csat*d(i,j,k)*exp(1.5*lntg(i,j,k))
        tmp  = bbx(i,j,k)*dtx(i,j,k)+ &
               bby(i,j,k)*dty(i,j,k)+ &
               bbz(i,j,k)*dtz(i,j,k)
        dfcx(i,j,k) = sigmac*(bbx(i,j,k)*tmp+bmx(i,j,k)*dtx(i,j,k))
        dfcy(i,j,k) = sigmac*(bby(i,j,k)*tmp+bmy(i,j,k)*dty(i,j,k))
        dfcz(i,j,k) = sigmac*(bbz(i,j,k)*tmp+bmz(i,j,k)*dtz(i,j,k))
        
        dfcx(i,j,k) = -(dfcx(i,j,k)/(1.+abs(dfcx(i,j,k))/fsat) + &
                        self%fc(i,j,k,1))/tau_spitz
        dfcy(i,j,k) = -(dfcy(i,j,k)/(1.+abs(dfcy(i,j,k))/fsat) + &
                        self%fc(i,j,k,2))/tau_spitz
        dfcz(i,j,k) = -(dfcz(i,j,k)/(1.+abs(dfcz(i,j,k))/fsat) + &
                        self%fc(i,j,k,3))/tau_spitz
      end do
    end do
  end do
  call deallocate_vectors (self%bb, self%bm,dtg)
  call deallocate_scalars (ux, uy, uz)
  call deallocate_scalars (lnd, self%sigmac, lntg, tg)
END SUBROUTINE conductive_flux

!=======================================================================
!> Compute \hat{b}, which is  \hat{b}= B/abs(B) the unit magnetic field
!> vector in the direction of the field 
!=======================================================================
SUBROUTINE compute_b_field_topology(self,ff)
  class (spitzer_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), target :: ff
  ! --------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bx, by, bz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bbx, bby, bbz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: bmx, bmy, bmz
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: scr1, scr2
  real(kind=KindScalarVar) :: a
  integer :: i,j,k, n(3), lb(3), ub(3)
  ! --------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin("spitzer_t::compute_b_field_topology")
  bx => ff(:,:,:,self%patch%idx%bx)
  by => ff(:,:,:,self%patch%idx%by)
  bz => ff(:,:,:,self%patch%idx%bz)
  bbx => self%bb(:,:,:,1); bmx => self%bm(:,:,:,1)
  bby => self%bb(:,:,:,2); bmy => self%bm(:,:,:,2)
  bbz => self%bb(:,:,:,3); bmz => self%bm(:,:,:,3)
  n=shape(bx)
  call allocate_scalars(n, scr1, scr2)
  do i=1,3
    lb(i) = self%patch%mesh(i)%lb
    ub(i) = self%patch%mesh(i)%ub
  end do
  !
  scr1=xdn(yup(by))
  scr2=xdn(zup(bz))
  if (n(1)>1) then
    do k=lb(3),ub(3)
      do j=lb(2),ub(2)
        do i=lb(1),ub(1)
          a = bx(i,j,k)**2+ scr1(i,j,k)**2+scr2(i,j,k)**2+bmin**2
          bbx(i,j,k) = bx(i,j,k)/sqrt(a)
          bmx(i,j,k) = bmin**2/a
        end do
      end do
    end do
    bbx = xup(bbx)
    bmx = xup(bmx)
  end if
  ! Fix boundary cells with zero gradient
  if (n(1)>1) bbx(ub(1),:,:) = bbx(ub(1)-1,:,:)
  if (n(2)>1) bbx(:,lb(2),:) = bbx(:,lb(2)+1,:)
  if (n(3)>1) bbx(:,:,lb(3)) = bbx(:,:,lb(3)+1)
  if (n(1)>1) bmx(ub(1),:,:) = bmx(ub(1)-1,:,:)
  if (n(2)>1) bmx(:,lb(2),:) = bmx(:,lb(2)+1,:)
  if (n(3)>1) bmx(:,:,lb(3)) = bmx(:,:,lb(3)+1)
  !
  if (n(2)>1) then 
    scr1=ydn(xup(bx))
    scr2=ydn(zup(bz))
    do k=lb(3),ub(3)
      do j=lb(2),ub(2)
        do i=lb(1),ub(1)
          a = by(i,j,k)**2+ scr1(i,j,k)**2+scr2(i,j,k)**2+bmin**2
          bby(i,j,k) = by(i,j,k)/sqrt(a)
          bmy(i,j,k) = bmin**2/a
        end do
      end do
    end do
    bby = yup(bby)
    bmy = yup(bmy)
  end if
  ! Fix boundary cells with zero gradient
  if (n(2)>1)bby(:,ub(2),:) = bby(:,ub(2)-1,:)
  if (n(1)>1)bby(lb(1),:,:) = bby(lb(1)+1,:,:)
  if (n(3)>1)  bby(:,:,lb(3)) = bby(:,:,lb(3)+1)
  if (n(2)>1)bmy(:,ub(2),:) = bmy(:,ub(2)-1,:)
  if (n(1)>1)bmy(lb(1),:,:) = bmy(lb(1)+1,:,:)
  if (n(3)>1)  bmy(:,:,lb(3)) = bmy(:,:,lb(3)+1)
  !
  if (n(3)>1) then
    scr1=zdn(xup(bx))
    scr2=zdn(yup(by))
    do k=lb(3),ub(3)
      do j=lb(2),ub(2)
        do i=lb(1),ub(1)
          a = bz(i,j,k)**2+ scr1(i,j,k)**2+scr2(i,j,k)**2+bmin**2
          bbz(i,j,k) = bz(i,j,k)/sqrt(a)
          bmz(i,j,k) = bmin**2/a
        end do
      end do
    end do
    bbz = zup(bbz)
    bmz = zup(bmz)
  end if
    ! Fix boundary cells with zero gradient
    if (n(3)>1)bbz(:,:,ub(3)) = bbz(:,:,ub(3)-1)
    if (n(2)>1)bbz(:,lb(2),:) = bbz(:,lb(2)+1,:)
    if (n(1)>1)bbz(lb(1),:,:) = bbz(lb(1)+1,:,:)
    if (n(3)>1)bmz(:,:,ub(3)) = bmz(:,:,ub(3)-1)
    if (n(2)>1)bmz(:,lb(2),:) = bmz(:,lb(2)+1,:)
    if (n(1)>1)bmz(lb(1),:,:) = bmz(lb(1)+1,:,:)
  !
  call deallocate_scalars(scr1,scr2)
  call trace%end()
END SUBROUTINE compute_b_field_topology

!===============================================================================
!> Intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(spitzer_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('spitzer_t%courant_condition')
  !print *, 'u_max(ee after spitzer)', self%patch%id, self%u_max, self%patch%u_max
  self%patch%u_max = max (self%patch%u_max, self%u_max)
  call trace%end()
END SUBROUTINE courant_condition


END MODULE spitzer_mod
