!> =============================================================================
!> RAMSES Godunov solvers
!>
!> godunov_fine
!>   unsplit
!>     ctoprim
!>     uslope
!>     trace3d
!>     cmpflxm
!>       riemann%hllc
!> =============================================================================
MODULE mhd_mod
  !.............................................................................
  USE io_unit_mod
  USE patch_mod
  USE riemann_mod
  use hydro_parameters
  use trace_mod
  implicit none
  PRIVATE
  type, public, extends(patch_t):: mhd_t
    logical:: mhd=.true.
    integer:: idim=3
    character(len=10):: riemann='hllc'
    real, pointer, dimension(:,:,:,:):: uin       ! conserved
    real, pointer, dimension(:,:,:,:):: unew      ! new conserved
    real, pointer, dimension(:,:,:,:):: q         ! primitive
    real, pointer, dimension(:,:,:,:):: fx        ! 1D fluxes
    real, pointer, dimension(:,:,:,:,:):: qp, qm  ! face centered
    real, pointer, dimension(:,:,:,:,:):: dq      ! slopes
    real, pointer, dimension(:,:,:,:,:):: fl      ! 3D fluxes
    real, pointer, dimension(:,:,:):: c, s        ! scratch
  contains
    procedure:: init
    procedure:: update
    procedure:: log_density
    procedure:: log_pressure
    procedure:: grav_potential
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  logical:: first_time=.true.
  character(len=16):: eqn='mhd'
  real(kind=8):: courant_number=0.2
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  character(len=10):: riemann='hllc'
  namelist /ramses_params/ gamma, riemann, slope_type, courant_factor, smallr, &
    smallc, do_isothermal, csound_isothermal
!...............................................................................
  call trace_begin('mhd_mod::init')
  call self%patch_t%init
  !$omp critical (pde)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, ramses_params)
    if (io_unit%master) write (*, ramses_params)
  end if
  !$omp end critical (pde)
  self%riemann = riemann
  self%courant = courant_factor
  self%gamma = gamma
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer:: n(3)
  !-----------------------------------------------------------------------------
  call trace_begin ('mhd_mod::update')
  self%uin => self%mem(:,:,:,:,self%it,1)
  self%unew => self%mem(:,:,:,:,self%new,1)
  n = self%mesh%gn
  allocate (self%q (n(1),n(2),n(3),self%nv))
  allocate (self%fx(n(1),n(2),n(3),self%nv))
  allocate (self%qp(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qm(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%dq(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%fl(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%c(n(1),n(2),n(3)))
  allocate (self%s(n(1),n(2),n(3)))
  call godunov (self)
  deallocate (self%q,self%fx,self%qp,self%qm,self%dq,self%fl,self%c,self%s)
  call trace_end
END SUBROUTINE update

!> =============================================================================
subroutine godunov (self)
  class(mhd_t)::self
  !-----------------------------------------------------------------------------
  ! This routine accesses the conserved variables in self. It then calls
  ! the Godunov solver that computes fluxes. These fluxes are zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array unew(:), both at the current level and at the
  ! coarser level if necessary.
  !-----------------------------------------------------------------------------
  integer::i0,j0,k0,i,j,k,idim,ivar,lo(3),uo(3)
  integer,save:: itimer=0
  associate (unew=>self%unew, uin=>self%uin, fl=>self%fl)
  !-----------------------------------------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------------------------------------
  call check_small (self, uin)
  call unsplit (self)
  !-----------------------------------------------------------------------------
  ! Conservative update, looping over internal cells only
  !-----------------------------------------------------------------------------
  call trace_begin ('mhd_mod::godunov', itimer=itimer)
  lo = self%mesh%lo
  uo = self%mesh%uo
  do idim=1,self%ndim
    i0 = merge(1,0,idim==1)
    j0 = merge(1,0,idim==2)
    k0 = merge(1,0,idim==3)
    do ivar=1,self%nv
      do k=lo(3),uo(3)
      do j=lo(2),uo(2)
      do i=lo(1),uo(1)
        ! Update conservative variables
        unew(i,j,k,ivar)=uin(i,j,k,ivar)+ &
          & (fl(i   ,j   ,k   ,ivar,idim) &
          & -fl(i+i0,j+j0,k+k0,ivar,idim))
      end do
      end do
      end do
    end do
  end do
  call trace_end (itimer)
  end associate
end subroutine godunov

!> =============================================================================
!> =============================================================================
subroutine check_small (self, uu)
  class(mhd_t):: self
  real(dp), dimension(:,:,:,:):: uu
  integer:: i,j,k,lo(3),uo(3),n_smalld,n_smallp
  real(dp):: c,d,u,v,w,p,ekin
  integer, save:: nprint=10
  !-----------------------------------------------------------------------------
  n_smalld = 0
  n_smallp = 0
  lo = self%mesh%li
  uo = self%mesh%ui
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  do i=lo(1),uo(1)
    if (uu(i,j,k,1) < smallr) then
      if (nprint>0) then
        nprint = nprint-1
        print '(i7,3i5,1p,2e12.3)',self%id,i,j,k,uu(i,j,k,1),smallr
      end if
      n_smalld = n_smalld+1
      uu(i,j,k,1) = smallr
    end if
    d =  uu(i,j,k,1)
    c = 1./d
    u =  uu(i,j,k,2)*c
    v =  uu(i,j,k,3)*c
    w =  uu(i,j,k,4)*c
    ekin = 0.5*d*(u**2+v**2+w**2)
    if (do_isothermal) then
      p = d*csound_isothermal**2
      uu(i,j,k,5) = p/(gamma-one)+ekin
    else
      p = (uu(i,j,k,5)-ekin)*(gamma-one)
      if (p < smallr*smallc**2) then
        p = smallr*smallc**2
        uu(i,j,k,5) = p/(gamma-one)+ekin
        n_smallp = n_smallp+1
      end if
    end if
  end do
  end do
  end do
  if (n_smalld+n_smallp > 0) then
    print*,'small: id,nr,np', self%id,n_smalld,n_smallp
  end if
end subroutine

!> =============================================================================
!>  UNSPLIT     Unsplit second order Godunov integrator for
!>              polytropic gas dynamics using MUSCL-HANCOCK
!> =============================================================================
subroutine unsplit(self)
  class(mhd_t):: self
  real(dp):: dtds(3)
  integer,save:: itimer=0
  associate (fl=>self%fl, fx=>self%fx)
  !-----------------------------------------------------------------------------
  call ctoprim(self)                  ! Primitive variables, sound speeds
  call self%courant_condition         ! compute dt
  dtds = self%dtime/self%mesh%d       ! needed here
  call uslope(self)                   ! TVD slopes
  call trace3d(self)                  ! predictor step
  call trace_begin ('mhd_mod::cmpflxm', itimer=itimer)
  call cmpflxm(self,2,3,4)            ! Solve for 1D flux in X direction
  fl(:,:,:,:,1)=fx(:,:,:,:)*dtds(1)   ! 1D to 3D flux
  call cmpflxm(self,3,2,4)            ! Solve for 1D flux in Y direction
  fl(:,:,:,:,2)=fx(:,:,:,:)*dtds(3)   ! 1D to 3D flux
  call cmpflxm(self,4,2,3)            ! Solve for 1D flux in Z direction
  fl(:,:,:,:,3)=fx(:,:,:,:)*dtds(3)   ! 1D to 3D flux
  call trace_end (itimer)
  end associate
end subroutine unsplit

!> =============================================================================
!> Conservative to primitive variables
!> =============================================================================
subroutine ctoprim (self)
  class(mhd_t):: self
  real(dp):: smallp, dtxhalf
  integer,save:: itimer=0
  associate (uin=>self%uin, q=>self%q, c=>self%c, s=>self%s)
  !-----------------------------------------------------------------------------
  call trace_begin ('mhd_mod::ctoprim', itimer=itimer)
  smallp = smallr*smallc**2/gamma
  q(:,:,:,1) = max(uin(:,:,:,1),smallr)
  s = 1.d0/q(:,:,:,1)
  q(:,:,:,2) = uin(:,:,:,2)*s
  q(:,:,:,3) = uin(:,:,:,3)*s
  q(:,:,:,4) = uin(:,:,:,4)*s
  s = q(:,:,:,2)**2+q(:,:,:,3)**2+q(:,:,:,4)**2
  q(:,:,:,5)=max((gamma-one)*(uin(:,:,:,5)-half*q(:,:,:,1)*s(:,:,:)),smallp)
  c(:,:,:)=sqrt(gamma*q(:,:,:,5)/q(:,:,:,1))
  self%u_max = self%fmaxval(c + abs(q(:,:,:,2)) &
                          + c + abs(q(:,:,:,3)) &
                          + c + abs(q(:,:,:,4)))
  call trace_end (itimer)
  end associate
end subroutine ctoprim

!> =============================================================================
!> Slope limiters
!> =============================================================================
subroutine uslope(self)
  class(mhd_t):: self
  integer::ngrid
  integer::i, j, k, n, l(3), u(3)
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, c1, c2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer,save:: itimer=0
  associate (q=>self%q, dq=>self%dq)
  !-----------------------------------------------------------------------------
  if (slope_type < 0) return
  call trace_begin ('mhd_mod::uslope', itimer=itimer)
  if(slope_type==0)then
     dq=0.0_dp
     return
  end if
  !
  l = self%mesh%lb+1
  u = self%mesh%ub-1
  if (slope_type==1) then  ! minmod
    do n = 1, self%nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = q(i  ,j,k,n) - q(i-1,j,k,n)
            drgt = q(i+1,j,k,n) - q(i  ,j,k,n)
            dq(i,j,k,n,1) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
            ! slopes in second coordinate direction
            dlft = q(i,j  ,k,n) - q(i,j-1,k,n)
            drgt = q(i,j+1,k,n) - q(i,j  ,k,n)
            dq(i,j,k,n,2) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
            ! slopes in third coordinate direction
            dlft = q(i,j,k  ,n) - q(i,j,k-1,n)
            drgt = q(i,j,k+1,n) - q(i,j,k  ,n)
            dq(i,j,k,n,3) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
          end do
        end do
      end do
    end do
  else if (slope_type==2) then ! moncen
    c1 = slope_type
    c2 = 0.5_dp/c1
    do n = 1, self%nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = c1*(q(i  ,j,k,n) - q(i-1,j,k,n))
            drgt = c1*(q(i+1,j,k,n) - q(i  ,j,k,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
          end do
          do i = l(1),u(1)
            ! slopes in second coordinate direction
            dlft = c1*(q(i,j  ,k,n) - q(i,j-1,k,n))
            drgt = c1*(q(i,j+1,k,n) - q(i,j  ,k,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one,dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
          end do
          do i = l(1),u(1)
            ! slopes in third coordinate direction
            dlft = c1*(q(i,j,k  ,n) - q(i,j,k-1,n))
            drgt = c1*(q(i,j,k+1,n) - q(i,j,k  ,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one,dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
          end do
        end do
      end do
    end do
  else if (slope_type==3) then ! positivity preserving 3d unsplit slope
    do n = 1, self%nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            dflll = q(i-1,j-1,k-1,n)-q(i,j,k,n)
            dflml = q(i-1,j  ,k-1,n)-q(i,j,k,n)
            dflrl = q(i-1,j+1,k-1,n)-q(i,j,k,n)
            dfmll = q(i  ,j-1,k-1,n)-q(i,j,k,n)
            dfmml = q(i  ,j  ,k-1,n)-q(i,j,k,n)
            dfmrl = q(i  ,j+1,k-1,n)-q(i,j,k,n)
            dfrll = q(i+1,j-1,k-1,n)-q(i,j,k,n)
            dfrml = q(i+1,j  ,k-1,n)-q(i,j,k,n)
            dfrrl = q(i+1,j+1,k-1,n)-q(i,j,k,n)
            !
            dfllm = q(i-1,j-1,k  ,n)-q(i,j,k,n)
            dflmm = q(i-1,j  ,k  ,n)-q(i,j,k,n)
            dflrm = q(i-1,j+1,k  ,n)-q(i,j,k,n)
            dfmlm = q(i  ,j-1,k  ,n)-q(i,j,k,n)
            dfmmm = q(i  ,j  ,k  ,n)-q(i,j,k,n)
            dfmrm = q(i  ,j+1,k  ,n)-q(i,j,k,n)
            dfrlm = q(i+1,j-1,k  ,n)-q(i,j,k,n)
            dfrmm = q(i+1,j  ,k  ,n)-q(i,j,k,n)
            dfrrm = q(i+1,j+1,k  ,n)-q(i,j,k,n)
            !
            dfllr = q(i-1,j-1,k+1,n)-q(i,j,k,n)
            dflmr = q(i-1,j  ,k+1,n)-q(i,j,k,n)
            dflrr = q(i-1,j+1,k+1,n)-q(i,j,k,n)
            dfmlr = q(i  ,j-1,k+1,n)-q(i,j,k,n)
            dfmmr = q(i  ,j  ,k+1,n)-q(i,j,k,n)
            dfmrr = q(i  ,j+1,k+1,n)-q(i,j,k,n)
            dfrlr = q(i+1,j-1,k+1,n)-q(i,j,k,n)
            dfrmr = q(i+1,j  ,k+1,n)-q(i,j,k,n)
            dfrrr = q(i+1,j+1,k+1,n)-q(i,j,k,n)
            !
            vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                 &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                 &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
            vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                 &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                 &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
            !
            dfx  = half*(q(i+1,j,k,n)-q(i-1,j,k,n))
            dfy  = half*(q(i,j+1,k,n)-q(i,j-1,k,n))
            dfz  = half*(q(i,j,k+1,n)-q(i,j,k-1,n))
            dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
            !
            slop = merge(min(one,min(abs(vmin),abs(vmax))/dff),one,dff>zero)
            dlim = slop
            !
            dq(i,j,k,n,1) = dlim*dfx
            dq(i,j,k,n,2) = dlim*dfy
            dq(i,j,k,n,3) = dlim*dfz
          end do
        end do
      end do
    end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
  call trace_end (itimer)
  end associate
end subroutine uslope

!> =============================================================================
!> Predictor step, using primitive variables
!> =============================================================================
subroutine trace3d(self)
  class(mhd_t):: self
  integer ::i, j, k, lb(3), ub(3)
  integer ::ir=1, iu=2, iv=3, iw=4, ip=5
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0
  integer,save:: itimer=0
  associate (q=>self%q, dq=>self%dq, qp=>self%qp, qm=>self%qm)
  !-----------------------------------------------------------------------------
  call trace_begin ('mhd_mod::trace3d', itimer=itimer)
  dtdx = self%dtime/self%mesh(1)%d
  dtdy = self%dtime/self%mesh(2)%d
  dtdz = self%dtime/self%mesh(3)%d
  !-----------------------------------------------------------------------------
  lb = self%mesh%lb
  ub = self%mesh%ub-1
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        ! Cell centered values
        r   =  q(i,j,k,ir)
        u   =  q(i,j,k,iu)
        v   =  q(i,j,k,iv)
        w   =  q(i,j,k,iw)
        p   =  q(i,j,k,ip)

        ! TVD slopes in all 3 directions
        drx = dq(i,j,k,ir,1)
        dpx = dq(i,j,k,ip,1)
        dux = dq(i,j,k,iu,1)
        dvx = dq(i,j,k,iv,1)
        dwx = dq(i,j,k,iw,1)

        dry = dq(i,j,k,ir,2)
        dpy = dq(i,j,k,ip,2)
        duy = dq(i,j,k,iu,2)
        dvy = dq(i,j,k,iv,2)
        dwy = dq(i,j,k,iw,2)

        drz = dq(i,j,k,ir,3)
        dpz = dq(i,j,k,ip,3)
        duz = dq(i,j,k,iu,3)
        dvz = dq(i,j,k,iv,3)
        dwz = dq(i,j,k,iw,3)

        ! Source terms (including transverse derivatives)
        sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
        sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*gamma*p
        su0 = -u*dux-v*duy-w*duz - (dpx        )/r
        sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
        sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r

        ! Right state at left interface
        qp(i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
        qp(i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
        qp(i,j,k,iu,1) = u - half*dux + su0*dtdx*half
        qp(i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
        qp(i,j,k,iw,1) = w - half*dwx + sw0*dtdx*half
        qp(i,j,k,ir,1) = max(smallr, qp(i,j,k,ir,1))

        ! Left state at right interface
        qm(i+1,j,k,ir,1) = r + half*drx + sr0*dtdx*half
        qm(i+1,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
        qm(i+1,j,k,iu,1) = u + half*dux + su0*dtdx*half
        qm(i+1,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
        qm(i+1,j,k,iw,1) = w + half*dwx + sw0*dtdx*half
        qm(i+1,j,k,ir,1) = max(smallr, qm(i+1,j,k,ir,1))

        ! Top state at bottom interface
        qp(i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
        qp(i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
        qp(i,j,k,iu,2) = u - half*duy + su0*dtdy*half
        qp(i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
        qp(i,j,k,iw,2) = w - half*dwy + sw0*dtdy*half
        qp(i,j,k,ir,2) = max(smallr, qp(i,j,k,ir,2))

        ! Bottom state at top interface
        qm(i,j+1,k,ir,2) = r + half*dry + sr0*dtdy*half
        qm(i,j+1,k,ip,2) = p + half*dpy + sp0*dtdy*half
        qm(i,j+1,k,iu,2) = u + half*duy + su0*dtdy*half
        qm(i,j+1,k,iv,2) = v + half*dvy + sv0*dtdy*half
        qm(i,j+1,k,iw,2) = w + half*dwy + sw0*dtdy*half
        qm(i,j+1,k,ir,2) = max(smallr, qm(i,j+1,k,ir,2))

        ! Back state at front interface
        qp(i,j,k,ir,3) = r - half*drz + sr0*dtdz*half
        qp(i,j,k,ip,3) = p - half*dpz + sp0*dtdz*half
        qp(i,j,k,iu,3) = u - half*duz + su0*dtdz*half
        qp(i,j,k,iv,3) = v - half*dvz + sv0*dtdz*half
        qp(i,j,k,iw,3) = w - half*dwz + sw0*dtdz*half
        qp(i,j,k,ir,3) = max(smallr, qp(i,j,k,ir,3))

        ! Front state at back interface
        qm(i,j,k+1,ir,3) = r + half*drz + sr0*dtdz*half
        qm(i,j,k+1,ip,3) = p + half*dpz + sp0*dtdz*half
        qm(i,j,k+1,iu,3) = u + half*duz + su0*dtdz*half
        qm(i,j,k+1,iv,3) = v + half*dvz + sv0*dtdz*half
        qm(i,j,k+1,iw,3) = w + half*dwz + sw0*dtdz*half
        qm(i,j,k+1,ir,3) = max(smallr, qm(i,j,k+1,ir,3))
      end do
    end do
  end do
  call trace_end (itimer)
  end associate
end subroutine trace3d

!> =============================================================================
!> Compute fluxes
!> =============================================================================
subroutine cmpflxm(self,ln,lt1,lt2)
  class(mhd_t):: self
  integer ::ln,lt1,lt2
  ! local variables
  integer ::i, j, k, n, l, idim, xdim, lb(3), ub(3), ngrid
  real(dp),dimension(self%mesh(1)%gn,1:self%nv)::qleft,qright,qgdnv,fgdnv
  integer,save:: itimer=0
  !-----------------------------------------------------------------------------
  associate (qp=>self%qp, qm=>self%qm, fx=>self%fx)
  xdim=ln-1
  ngrid = self%mesh(1)%gn
  !
  lb = self%mesh%lb
  ub = self%mesh%ub
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      ! Mass density
      qleft (:,1) = qm(:,j,k,1,xdim)
      qright(:,1) = qp(:,j,k,1,xdim)
      !
      ! Normal velocity
      qleft (:,2) = qm(:,j,k,ln,xdim)
      qright(:,2) = qp(:,j,k,ln,xdim)
      !
      ! Pressure
      qleft (:,3) = qm(:,j,k,5,xdim)
      qright(:,3) = qp(:,j,k,5,xdim)
      !
      ! Tangential velocity 1
      qleft (:,4) = qm(:,j,k,lt1,xdim)
      qright(:,4) = qp(:,j,k,lt1,xdim)
      !
      ! Tangential velocity 2
      qleft (:,5) = qm(:,j,k,lt2,xdim)
      qright(:,5) = qp(:,j,k,lt2,xdim)
      !
      ! Solve Riemann problem
      if(self%riemann.eq.'acoustic')then
        call riemann%acoustic(qleft,qright,qgdnv,fgdnv,ngrid,self%nv)
      else if (self%riemann.eq.'exact')then
        call riemann%approx  (qleft,qright,qgdnv,fgdnv,ngrid,self%nv)
      else if (self%riemann.eq.'llf')then
        call riemann%llf     (qleft,qright,qgdnv,fgdnv,ngrid,self%nv)
      else if (self%riemann.eq.'hllc')then
        call riemann%hllc    (qleft,qright,qgdnv,fgdnv,ngrid,self%nv)
      else if (self%riemann.eq.'hll')then
        call riemann%hll     (qleft,qright,qgdnv,fgdnv,ngrid,self%nv)
      else
        write(*,*)'unknown Riemann solver'
        stop
      end if
      !
      ! Compute fluxes
      fx(:,j,k,1)      = fgdnv(:,1)
      fx(:,j,k,ln)     = fgdnv(:,2)
      fx(:,j,k,lt1)    = fgdnv(:,4)
      fx(:,j,k,lt2)    = fgdnv(:,5)
      fx(:,j,k,ndim+2) = fgdnv(:,3)
    end do
  end do
  end associate
end subroutine cmpflxm

!*******************************************************************************
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE log_pressure (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE grav_potential (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE velocity_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!*******************************************************************************
SUBROUTINE magnetic_field_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

END MODULE
