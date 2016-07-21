!> =============================================================================
!> =============================================================================
!> RAMSES Godunov solvers
!>
!> godunov_fine
!>   unsplit              1   2   3      ui  uo  ub
!>     ctoprim            +---+---+......+---+---+
!>     uslope:dq              +---+......+---+
!>     trace3d:qp             +---+......+---+
!>     trace3d:qm                 +......+---+---+
!>     cmpflxm                    +......+---+
!>       riemann%hllc             +......+
!> =============================================================================
MODULE mhd_mod
  !.............................................................................
  USE io_unit_mod
  USE patch_mod
  USE link_mod
  USE riemann_mod
  use hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  implicit none
  PRIVATE
  type, public, extends(patch_t):: mhd_t
    logical:: mhd=.true.
    integer:: idim=3
    character(len=10):: riemann='hllc'
    real, pointer, dimension(:,:,:,:):: uin       ! conserved
    real, pointer, dimension(:,:,:,:):: unew      ! new conserved
    real, pointer, dimension(:,:,:,:):: q         ! primitive
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
    smallc, do_isothermal, csound
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
  allocate (self%qp(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qm(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%dq(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%fl(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%c(n(1),n(2),n(3)))
  allocate (self%s(n(1),n(2),n(3)))
  call godunov (self)
  deallocate (self%q,self%qp,self%qm,self%dq,self%fl,self%c,self%s)
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
  integer::i0,j0,k0,i,j,k,idim,ivar,li(3),ui(3)
  integer,save:: itimer=0
  associate (unew=>self%unew, uin=>self%uin, fl=>self%fl)
  !-----------------------------------------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------------------------------------
  li = self%mesh%li
  ui = self%mesh%ui
  call check_small (self, uin)
  call unsplit (self)
  !-----------------------------------------------------------------------------
  ! Conservative update, looping over internal cells only
  !-----------------------------------------------------------------------------
  call trace_begin ('mhd_mod::godunov', itimer=itimer)
  unew = uin
  do ivar=1,self%nv
    do idim=1,self%ndim
      i0 = merge(1,0,idim==1)
      j0 = merge(1,0,idim==2)
      k0 = merge(1,0,idim==3)
      do k=li(3),ui(3)
      do j=li(2),ui(2)
      do i=li(1),ui(1)
        ! Update conservative variables
        unew(i,j,k,ivar) = unew(i,j,k,ivar)+ &
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
  lo = self%mesh%lb
  uo = self%mesh%ub
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  do i=lo(1),uo(1)
    if (uu(i,j,k,1) < smallr) then
      if (nprint>0) then
        nprint = nprint-1
        write (io_unit%log,'(i7,3i5,1p,2e12.3)') self%id,i,j,k,uu(i,j,k,1),smallr
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
      p = d*csound**2
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
    write (io_unit%log,*) 'small: id,nr,np', self%id,n_smalld,n_smallp
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
  !-----------------------------------------------------------------------------
  call ctoprim(self)                  ! Primitive variables, sound speeds
  call self%courant_condition         ! compute dt
  call uslope(self)                   ! TVD slopes
  call trace3d(self)                  ! predictor step
  call trace_begin ('mhd_mod::cmpflxm', itimer=itimer)
  call cmpflxm(self)                  ! Solve for 1D flux in X direction
  call trace_end (itimer)
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
  call trace_begin ('mhd_mod::uslope', itimer=itimer)
  !
  l = self%mesh%lo
  u = self%mesh%uo
  if(slope_type==-1)then
    do n = 1, self%nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            dq(i,j,k,n,1) = half*(q(i+1,j,k,n) - q(i-1,j,k,n))
            dq(i,j,k,n,2) = half*(q(i,j+1,k,n) - q(i,j-1,k,n))
            dq(i,j,k,n,3) = half*(q(i,j,k+1,n) - q(i,j,k-1,n))
          end do
        end do
      end do
    end do
  else if(slope_type==0)then
     dq=0.0_dp
     return
  else if (slope_type==1) then  ! minmod
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
!          end do
!          do i = l(1),u(1)
            ! slopes in second coordinate direction
            dlft = c1*(q(i,j  ,k,n) - q(i,j-1,k,n))
            drgt = c1*(q(i,j+1,k,n) - q(i,j  ,k,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one,dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
!          end do
!          do i = l(1),u(1)
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
  lb = self%mesh%lo
  ub = self%mesh%uo
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
subroutine cmpflxm(self)
  class(mhd_t):: self
  ! local variables
  integer:: i, j, k, lb(3), ub(3)
  real(dp), dimension(self%mesh(1)%gn,self%nv):: qleft, qright, fgdnv
  real(dp):: dtds(3)
  logical:: error
  associate (qp=>self%qp, qm=>self%qm, fl=>self%fl)
  !-----------------------------------------------------------------------------
  error = .false.
  dtds = self%dtime/self%mesh%d
  lb = self%mesh%li
  ub = self%mesh%uo
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      qleft (lb(1):ub(1),1) = qm(lb(1):ub(1),j,k,1,1); qright(lb(1):ub(1),1) = qp(lb(1):ub(1),j,k,1,1)
      qleft (lb(1):ub(1),2) = qm(lb(1):ub(1),j,k,2,1); qright(lb(1):ub(1),2) = qp(lb(1):ub(1),j,k,2,1)
      qleft (lb(1):ub(1),3) = qm(lb(1):ub(1),j,k,5,1); qright(lb(1):ub(1),3) = qp(lb(1):ub(1),j,k,5,1)
      qleft (lb(1):ub(1),4) = qm(lb(1):ub(1),j,k,3,1); qright(lb(1):ub(1),4) = qp(lb(1):ub(1),j,k,3,1)
      qleft (lb(1):ub(1),5) = qm(lb(1):ub(1),j,k,4,1); qright(lb(1):ub(1),5) = qp(lb(1):ub(1),j,k,4,1)
      call solver (self,self%mesh(1)%gn,qleft,qright,fgdnv,error)
      if (error) call error_handler('x', j, k, lb(1), ub(1))
      fl(lb(1):ub(1),j,k,1,1) = fgdnv(lb(1):ub(1),1)*dtds(1)
      fl(lb(1):ub(1),j,k,2,1) = fgdnv(lb(1):ub(1),2)*dtds(1)
      fl(lb(1):ub(1),j,k,5,1) = fgdnv(lb(1):ub(1),3)*dtds(1)
      fl(lb(1):ub(1),j,k,3,1) = fgdnv(lb(1):ub(1),4)*dtds(1)
      fl(lb(1):ub(1),j,k,4,1) = fgdnv(lb(1):ub(1),5)*dtds(1)
    end do
  end do
  if (any(self%mesh%gn/=self%mesh(1)%gn)) then
    call cmpflxm_y
    call cmpflxm_z
  else
    do k = lb(3),ub(3)
      do i = lb(1),ub(1)
        qleft (lb(2):ub(2),1) = qm(i,lb(2):ub(2),k,1,2); qright(lb(2):ub(2),1) = qp(i,lb(2):ub(2),k,1,2)
        qleft (lb(2):ub(2),2) = qm(i,lb(2):ub(2),k,3,2); qright(lb(2):ub(2),2) = qp(i,lb(2):ub(2),k,3,2)
        qleft (lb(2):ub(2),3) = qm(i,lb(2):ub(2),k,5,2); qright(lb(2):ub(2),3) = qp(i,lb(2):ub(2),k,5,2)
        qleft (lb(2):ub(2),4) = qm(i,lb(2):ub(2),k,2,2); qright(lb(2):ub(2),4) = qp(i,lb(2):ub(2),k,2,2)
        qleft (lb(2):ub(2),5) = qm(i,lb(2):ub(2),k,4,2); qright(lb(2):ub(2),5) = qp(i,lb(2):ub(2),k,4,2)
        call solver (self,self%mesh(2)%gn,qleft,qright,fgdnv,error)
        if (error) call error_handler('y', i, k, lb(2), ub(2))
        fl(i,lb(2):ub(2),k,1,2) = fgdnv(lb(2):ub(2),1)*dtds(2)
        fl(i,lb(2):ub(2),k,3,2) = fgdnv(lb(2):ub(2),2)*dtds(2)
        fl(i,lb(2):ub(2),k,5,2) = fgdnv(lb(2):ub(2),3)*dtds(2)
        fl(i,lb(2):ub(2),k,2,2) = fgdnv(lb(2):ub(2),4)*dtds(2)
        fl(i,lb(2):ub(2),k,4,2) = fgdnv(lb(2):ub(2),5)*dtds(2)
      end do
    end do
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        qleft (lb(3):ub(3),1) = qm(i,j,lb(3):ub(3),1,3); qright(lb(3):ub(3),1) = qp(i,j,lb(3):ub(3),1,3)
        qleft (lb(3):ub(3),2) = qm(i,j,lb(3):ub(3),4,3); qright(lb(3):ub(3),2) = qp(i,j,lb(3):ub(3),4,3)
        qleft (lb(3):ub(3),3) = qm(i,j,lb(3):ub(3),5,3); qright(lb(3):ub(3),3) = qp(i,j,lb(3):ub(3),5,3)
        qleft (lb(3):ub(3),4) = qm(i,j,lb(3):ub(3),2,3); qright(lb(3):ub(3),4) = qp(i,j,lb(3):ub(3),2,3)
        qleft (lb(3):ub(3),5) = qm(i,j,lb(3):ub(3),3,3); qright(lb(3):ub(3),5) = qp(i,j,lb(3):ub(3),3,3)
        call solver (self,self%mesh(3)%gn,qleft,qright,fgdnv,error)
        if (error) call error_handler('z', i, j, lb(3), ub(3))
        fl(i,j,lb(3):ub(3),1,3) = fgdnv(lb(3):ub(3),1)*dtds(3)
        fl(i,j,lb(3):ub(3),4,3) = fgdnv(lb(3):ub(3),2)*dtds(3)
        fl(i,j,lb(3):ub(3),5,3) = fgdnv(lb(3):ub(3),3)*dtds(3)
        fl(i,j,lb(3):ub(3),2,3) = fgdnv(lb(3):ub(3),4)*dtds(3)
        fl(i,j,lb(3):ub(3),3,3) = fgdnv(lb(3):ub(3),5)*dtds(3)
      end do
    end do
  end if
  end associate
contains
  subroutine cmpflxm_y
  ! local variables
  real(dp), dimension(self%mesh(2)%gn,self%nv):: qleft, qright, fgdnv
  associate (qp=>self%qp, qm=>self%qm, fl=>self%fl)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  dtds = self%dtime/self%mesh%d       ! needed here
  lb = self%mesh%li
  ub = self%mesh%uo
  do k = lb(3),ub(3)
    do i = lb(1),ub(1)
      qleft (:,1) = qm(i,:,k,1,2); qright(:,1) = qp(i,:,k,1,2)
      qleft (:,2) = qm(i,:,k,3,2); qright(:,2) = qp(i,:,k,3,2)
      qleft (:,3) = qm(i,:,k,5,2); qright(:,3) = qp(i,:,k,5,2)
      qleft (:,4) = qm(i,:,k,2,2); qright(:,4) = qp(i,:,k,2,2)
      qleft (:,5) = qm(i,:,k,4,2); qright(:,5) = qp(i,:,k,4,2)
      call solver (self,self%mesh(2)%gn,qleft,qright,fgdnv,error)
      if (error) call error_handler('y', i, k, lb(2), ub(2))
      fl(i,:,k,1,2) = fgdnv(:,1)*dtds(2)
      fl(i,:,k,3,2) = fgdnv(:,2)*dtds(2)
      fl(i,:,k,5,2) = fgdnv(:,3)*dtds(2)
      fl(i,:,k,2,2) = fgdnv(:,4)*dtds(2)
      fl(i,:,k,4,2) = fgdnv(:,5)*dtds(2)
    end do
  end do
  end associate
  end subroutine cmpflxm_y
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine cmpflxm_z
  ! local variables
  integer:: i, j, k, lb(3), ub(3)
  real(dp), dimension(self%mesh(3)%gn,self%nv):: qleft, qright, fgdnv
  associate (qp=>self%qp, qm=>self%qm, fl=>self%fl)
  !-----------------------------------------------------------------------------
  dtds = self%dtime/self%mesh%d       ! needed here
  lb = self%mesh%li
  ub = self%mesh%uo
  do j = lb(2),ub(2)
    do i = lb(1),ub(1)
      qleft (:,1) = qm(i,j,:,1,3); qright(:,1) = qp(i,j,:,1,3)
      qleft (:,2) = qm(i,j,:,4,3); qright(:,2) = qp(i,j,:,4,3)
      qleft (:,3) = qm(i,j,:,5,3); qright(:,3) = qp(i,j,:,5,3)
      qleft (:,4) = qm(i,j,:,2,3); qright(:,4) = qp(i,j,:,2,3)
      qleft (:,5) = qm(i,j,:,3,3); qright(:,5) = qp(i,j,:,3,3)
      call solver (self,self%mesh(3)%gn,qleft,qright,fgdnv,error)
      if (error) call error_handler('z', i, j, lb(3), ub(3))
      fl(i,j,:,1,3) = fgdnv(:,1)*dtds(3)
      fl(i,j,:,4,3) = fgdnv(:,2)*dtds(3)
      fl(i,j,:,5,3) = fgdnv(:,3)*dtds(3)
      fl(i,j,:,2,3) = fgdnv(:,4)*dtds(3)
      fl(i,j,:,3,3) = fgdnv(:,5)*dtds(3)
    end do
  end do
  end associate
  end subroutine cmpflxm_z
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine error_handler (label, i1, i2, lb, ub)
  class(*), pointer:: link
  character(len=*):: label
  integer:: i1, i2, lb, ub
  class(link_t), pointer:: nbor
  !.............................................................................
  write (io_unit%log,*) 'solver error: id, i[12] =', self%id, i1, i2
  !-----------------------------------------------------------------------------
  ! Print the left and right states that give errors
  !-----------------------------------------------------------------------------
  do i=lb,ub
    write (io_unit%log,'(a,1p,2(5e12.3,2x))') &
      trim(label)//'solver error:', qleft(i,1:5), qright(i,1:5)
  end do
  !-----------------------------------------------------------------------------
  ! Print a list of neighbors and their time differences
  !-----------------------------------------------------------------------------
  link => self%link
  select type (link)
  class is (link_t)
    nbor => link%nbor
  class default
    nullify(nbor)
  end select
  do while (associated(nbor))
    write (io_unit%log,'(a,i9,f10.3,3x,3f6.2)') &
      trim(label)//'solver error: nbor, time diff', &
      nbor%task%id, (nbor%task%time-self%time)/self%dtime, &
      (nbor%task%position-self%position)/self%size
    nbor => nbor%next
  end do
  call mpi%abort ('HLLC')
  end subroutine
end subroutine cmpflxm

!> =============================================================================
!> Selection of 1-D solvers
!> =============================================================================
subroutine solver (self,ngrid,qleft,qright,fgdnv,error)
  class(mhd_t):: self
  integer:: ngrid
  real(dp),dimension(ngrid,self%nv):: qleft,qright,fgdnv
  logical:: error
  !.............................................................................
  if(self%riemann.eq.'acoustic')then
    call riemann%acoustic(qleft,qright,fgdnv,ngrid,self%nv)
  else if (self%riemann.eq.'exact')then
    call riemann%approx  (qleft,qright,fgdnv,ngrid,self%nv)
  else if (self%riemann.eq.'llf')then
    call riemann%llf     (qleft,qright,fgdnv,ngrid,self%nv)
  else if (self%riemann.eq.'hllc')then
    call riemann%hllc    (qleft,qright,fgdnv,ngrid,self%nv,error)
  else if (self%riemann.eq.'hllc_a')then
    call riemann%hllc_a  (qleft,qright,fgdnv,ngrid,self%nv)
  else if (self%riemann.eq.'hllc_v')then
    call riemann%hllc_v  (qleft,qright,fgdnv,ngrid,self%nv)
  else if (self%riemann.eq.'hll')then
    call riemann%hll     (qleft,qright,fgdnv,ngrid,self%nv)
  else
    write(*,*)'unknown Riemann solver'
    stop
  end if
end subroutine solver

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!> =============================================================================
SUBROUTINE log_pressure (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!> =============================================================================
SUBROUTINE grav_potential (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:):: v
  !-----------------------------------------------------------------------------
END SUBROUTINE

END MODULE
