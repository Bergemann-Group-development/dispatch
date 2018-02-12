!> =============================================================================
!> =============================================================================
!> RAMSES Godunov solvers
!>
!> godunov
!>   unsplit              1   2   3      ui  uo  ub
!>     ctoprim            +---+---+......+---+---+
!>     uslope:dq              +---+......+---+
!>     trace3d:qp             +---+......+---+
!>     trace3d:qm                 +......+---+---+
!>     cmpflxm                    +......+---+
!>       cmp_mag_flx
!>         cmp_mag_flx_hlld
!> =============================================================================
MODULE mhd_mod
  !.............................................................................
  USE const
  USE io_mod
  USE io_unit_mod
  USE gpatch_mod
  USE link_mod
  USE riemann_mod
  use hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE force_mod
  USE index_mod
  implicit none
  PRIVATE
  type, public, extends(gpatch_t):: mhd_t
    logical:: mhd=.true.
    integer:: idim=3
    character(len=10):: riemann='hlld'
    real, pointer, dimension(:,:,:,:):: q         ! primitive
    real, pointer, dimension(:,:,:,:):: bf, em    ! face centered B, emf
    real, pointer, dimension(:,:,:,:,:):: qp, qm  ! face centered
    real, pointer, dimension(:,:,:,:,:):: dq, dbf ! slopes
    real, pointer, dimension(:,:,:,:,:):: qLB,qRB,qLT,QRT
    real, pointer, dimension(:,:,:,:):: uin       ! conserved
    real, pointer, dimension(:,:,:,:):: unew      ! new conserved
    real, pointer, dimension(:,:,:):: Ex,Ey,Ez
    real, pointer, dimension(:,:,:):: c
    real, pointer, dimension(:,:,:,:):: ff        ! external force
    type(force_t):: force
  contains
    procedure:: init
    procedure:: update
    procedure:: gas_pressure
  end type
  logical:: first_time=.true.
  character(len=16):: eqn='mhd'
  real(kind=8):: courant_number=0.2
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  character(len=10):: riemann='hlld'
  type(index_t):: idx
  real, save:: csound=1.
  integer:: i
  namelist /ramses_params/ gamma, riemann, slope_type, courant_factor, smallr, &
    smallc, do_isothermal, csound
!...............................................................................
  call trace%begin('mhd_mod%init')

  self%nw = 1
  self%kind = 'ramses_mhd_patch'
  if (self%nv == 0) then
    self%nv = 5
    if (self%mhd) self%nv = self%nv + 3
  end if
  call self%gpatch_t%init
  !-----------------------------------------------------------------------------
  ! Read solver-relevant parameters and store in the task instance
  !-----------------------------------------------------------------------------
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
  self%staggered = .false.
  do i=1,3
    self%mesh(i)%h = 0d0
  end do
  call self%force%init(self%kind)

  call trace%end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  real, pointer, dimension(:,:,:,:):: em, p
  real, pointer, dimension(:,:,:,:,:):: fl
  real, pointer, dimension(:,:,:):: d, e, Ux=>null(), Uy=>null(), Uz=>null()
  real(dp):: dtdx, dtdy, dtdz
  integer:: n(3), nv3, i
  integer, parameter:: ndim=3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%update', itimer=itimer)
  nv3 = self%nv-3
  dtdx = self%dtime/self%ds(1)
  dtdy = self%dtime/self%ds(2)
  dtdz = self%dtime/self%ds(3)
  
  call self%gravity
  !----------------------------------------------------------------------------
  ! Checks whether to do output
  !----------------------------------------------------------------------------
  call self%output
  !-----------------------------------------------------------------------------
  ! Memory mapping and allocations.  All alloc and dealloc should be done here,
  ! even if it would be possible to do some of them locally, in the procedures
  ! being called.
  !-----------------------------------------------------------------------------
  self%uin  => self%mem(:,:,:,:,self%it ,1)
  self%unew => self%mem(:,:,:,:,self%new,1)
  n = self%mesh%gn
  allocate (self%c  (n(1),n(2),n(3)))
  allocate (self%Ex (n(1),n(2),n(3)))
  allocate (self%Ey (n(1),n(2),n(3)))
  allocate (self%Ez (n(1),n(2),n(3)))
  allocate (self%q  (n(1),n(2),n(3),self%nv))
  allocate (self%qp (n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qm (n(1),n(2),n(3),self%nv,ndim))
  allocate (self%dq (n(1),n(2),n(3),self%nv,ndim))
  allocate (self%dbf(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qLB(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qLT(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qRB(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qRT(n(1),n(2),n(3),self%nv,ndim))
  d => self%mem(:,:,:,          self%idx%d ,self%it,1)
  e => self%mem(:,:,:,          self%idx%e ,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  ! External force; note that Ux, Uy, Uz are empty pointers here, so they
  ! should not and cannot be used in the solver
  !-----------------------------------------------------------------------------
  allocate (self%ff(size(d,1),size(d,2),size(d,3),3))
  if (associated(self%force%selected)) then
    self%ff = self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  else
    self%ff = 0.0
  end if
  !-----------------------------------------------------------------------------
  ! Compute update using second-order Godunov method
  !-----------------------------------------------------------------------------
  call check_small (self)
  call ctoprim(self)                                    ! Primitive variables
  call self%courant_condition                           ! compute dt
  call uslope(self,self%q,self%dq,self%nv)              ! TVD slopes
  call uslope(self,self%bf,self%dbf,3)                  ! TVD slopes
  call trace3d(self)                                    ! predictor step

  allocate (fl(self%gn(1),self%gn(2),self%gn(3),self%nv,3))
  allocate (em(self%gn(1),self%gn(2),self%gn(3),3))

  call cmpflxm(self,self%gn(1),2,3,4,6,7,8,fl(:,:,:,:,1))
  call cmpflxm(self,self%gn(2),3,2,4,7,6,8,fl(:,:,:,:,2))
  call cmpflxm(self,self%gn(3),4,2,3,8,6,7,fl(:,:,:,:,3))

  call cmp_mag_flx(self,2,3,4,6,7,8,dtdz,em(:,:,:,3),self%gn(3),self%nv)
  call cmp_mag_flx(self,4,2,3,8,6,7,dtdy,em(:,:,:,2),self%gn(2),self%nv)
  call cmp_mag_flx(self,3,4,2,7,8,6,dtdx,em(:,:,:,1),self%gn(1),self%nv)
  !-----------------------------------------------------------------------------
  ! Update conservative variables
  !-----------------------------------------------------------------------------
  associate (i0=>self%mesh(1)%li, j0=>self%mesh(2)%li, k0=>self%mesh(3)%li, &
             i1=>self%mesh(1)%ui, j1=>self%mesh(2)%ui, k1=>self%mesh(3)%ui, &
             uin=>self%uin,unew=>self%unew)
  unew(1:i0-1,:,:,:) = uin(1:i0-1,:,:,:); unew(i1+1:,:,:,:) = uin(i1+1:,:,:,:)
  unew(:,1:j0-1,:,:) = uin(:,1:j0-1,:,:); unew(:,j1+1:,:,:) = uin(:,j1+1:,:,:)
  unew(:,:,1:k0-1,:) = uin(:,:,1:k0-1,:); unew(:,:,k1+1:,:) = uin(:,:,k1+1:,:)
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) = &
   uin(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,1) - &
    fl(i0+1:i1+1,j0  :j1  ,k0  :k1  ,1:nv3,1) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,2) - &
    fl(i0  :i1  ,j0+1:j1+1,k0  :k1  ,1:nv3,2) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,3) - &
    fl(i0  :i1  ,j0  :j1  ,k0+1:k1+1,1:nv3,3)
  unew = uin
  end associate
  !-----------------------------------------------------------------------------
  ! Add external force as source term
  !-----------------------------------------------------------------------------
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  do i=1,3
    p(:,:,:,i) = p(:,:,:,i) + self%dtime*self%ff(:,:,:,i)*d
  end do
  !-----------------------------------------------------------------------------
  ! Add gradient of grav. potential as source term (simple, 1st-order approach)
  !-----------------------------------------------------------------------------
  if (self%selfgravity) then
    call gravitational_force(self)
    call self%gpatch_t%update
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate
  !-----------------------------------------------------------------------------
  deallocate (fl,em,self%q,self%qp,self%qm,self%dq,self%dbf, &
              self%qLT,self%qLB,self%qRT,self%qRB,self%c, &
              self%Ex,self%Ey,self%Ez)
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
subroutine check_small (self)
  class(mhd_t):: self
  integer:: i,j,k,lo(3),uo(3),n_smalld,n_smallp
  real(dp):: c,d,u,v,w,p,ekin
  integer, save:: nprint=10
  real, save:: csound=1.
  integer, save:: itimer=0
  associate(uu=>self%uin)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_mod::check_small',itimer=itimer)
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
  end associate
  call trace%end (itimer)
end subroutine check_small

!> =============================================================================
!> Conservative to primitive variables
!> =============================================================================
subroutine ctoprim(self)
  class(mhd_t):: self
  integer::n(4)
  real(dp):: smallp, c
  integer, save:: itimer=0
  associate(uin=>self%uin,q=>self%q,bf=>self%bf,dt=>self%dtime,c=>self%c)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_mod::ctoprim',itimer=itimer)
  smallp = smallr*smallc**2/self%gamma
  n = shape(uin)-1

  self%bf => uin(:,:,:,6:8)                                          ! Face centered B
  q(:,:,:,1) = max(uin(:,:,:,1),smallr)                         ! Density
  q(:,:,:,2) = uin(:,:,:,2)/q(:,:,:,1)                          ! vx
  q(:,:,:,3) = uin(:,:,:,3)/q(:,:,:,1)                          ! vy
  q(:,:,:,4) = uin(:,:,:,4)/q(:,:,:,1)                          ! vz
  q(1:n(1),:,:,6) = (uin(1:n(1),:,:,6)+uin(2:,:,:,6))*half      ! Bx
  q(:,1:n(2),:,7) = (uin(:,1:n(2),:,7)+uin(:,2:,:,7))*half      ! By
  q(:,:,1:n(3),8) = (uin(:,:,1:n(3),8)+uin(:,:,2:,8))*half      ! Bz
  q(n(1)+1,:,:,6) = 0.0
  q(:,n(2)+1,:,7) = 0.0
  q(:,:,n(3)+1,8) = 0.0
  q(:,:,:,5) = MAX((self%gamma-one)*(uin(:,:,:,5) -&                 ! Gas pressure
    half*(q(:,:,:,6)**2+q(:,:,:,7)**2+q(:,:,:,8)**2)) - &       ! Magnetic pressure
    half*q(:,:,:,1)*(q(:,:,:,2)**2+q(:,:,:,3)**2+q(:,:,:,4)**2),smallp)
  ! Gravity predictor step
  !q(:,:,:,2:4) = q(:,:,:,idim+1) + Grho*gravin(:,:,:,1:3)*dt*half
  c = sqrt((q(:,:,:,5)*self%gamma+q(:,:,:,6)**2+q(:,:,:,7)**2+q(:,:,:,8)**2)/q(:,:,:,1))
  self%u_max = self%fmaxval(c + abs(q(:,:,:,2)) &
                          + c + abs(q(:,:,:,3)) &
                          + c + abs(q(:,:,:,4)))
  end associate
  call trace%end (itimer)
end subroutine ctoprim

!> =============================================================================
!> Slope limiters
!> =============================================================================
subroutine uslope(self,q,dq,nv)
  class(mhd_t):: self
  real, dimension(:,:,:,:):: q
  real, dimension(size(q,1),size(q,2),size(q,3),size(q,4),3):: dq
  integer::nv
  integer::i, j, k, n, l(3), u(3)
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, c1, c2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer,save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_mod::uslope', itimer=itimer)
  !
  l = self%mesh%lo
  u = self%mesh%uo
  if(slope_type==-1)then
    do n = 1, nv
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
    do n = 1, nv
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
    do n = 1, nv
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
  call trace%end (itimer)
end subroutine uslope

!> =============================================================================
!> Predictor step
!> =============================================================================
SUBROUTINE trace3d(self)
  IMPLICIT NONE
  class(mhd_t):: self
  real(dp):: dx, dy, dz, dt
  integer :: i, j, k, n
  integer, parameter ::ir=1, iu=2, iv=3, iw=4, ip=5, ia=6, ib=7, ic=8
  real(dp)::dtdx, dtdy, dtdz, smallp
  real(dp) :: invr
!  real(dp) :: r, u, v, w, p, A, B, C
!  real(dp) :: ELL, ELR, ERL, ERR
!  real(dp) :: FLL, FLR, FRL, FRR
!  real(dp) :: GLL, GLR, GRL, GRR
  real(dp) :: drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  real(dp) :: dry, duy, dvy, dwy, dpy, dAy, dbY, dCy
  real(dp) :: drz, duz, dvz, dwz, dpz, dAz, dBz, dCz
  real(dp) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
!  real(dp) :: AL, AR, BL, BR, CL, CR
  real(dp) :: dALy, dARy, dALz, dARz
  real(dp) :: dBLx, dBRx, dBLz, dBRz
  real(dp) :: dCLx, dCRx, dCLy, dCRy
  real(dp) :: sAL0, sAR0, sBL0, sBR0, sCL0, sCR0
  real(dp) :: eff_gamma
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  associate(q=>self%q,bf=>self%bf,dq=>self%dq,qm=>self%qm,qp=>self%qp, &
    dbf=>self%dbf,qRT=>self%qRT,qRB=>self%qRB,qLT=>self%qLT,qLB=>self%qLB, &
    Bx=>self%bf(:,:,:,1),By=>self%bf(:,:,:,2),Bz=>self%bf(:,:,:,3), &
    Ux=>self%q(:,:,:,2),Uy=>self%q(:,:,:,3),Uz=>self%q(:,:,:,4), &
    Ex=>self%Ex,Ey=>self%Ey,Ez=>self%Ez)
  !
  call trace%begin ('mhd_mod::trace3d', itimer=itimer)
  su0=0.0_dp; sv0=0.0_dp; sw0=0.0_dp
  !
  dtdx = self%dtime/self%ds(1)
  dtdy = self%dtime/self%ds(2)
  dtdz = self%dtime/self%ds(3)
  smallp = smallr*smallc**2/self%gamma
  !
  ! E = - UxB = BxU
  Ex = zdn(xdn(Uz))*zdn(By) - ydn(zdn(Uy))*ydn(Bz)
  Ey = xdn(ydn(Ux))*xdn(Bz) - zdn(xdn(Uz))*zdn(Bx)
  Ez = ydn(zdn(Uy))*ydn(Bx) - xdn(ydn(Ux))*xdn(By)
  !
  eff_gamma = self%gamma
  do k=1,self%gn(3)-1
  do j=1,self%gn(2)-1
  do i=1,self%gn(1)-1
    associate ( &
    r=>q(i,j,k,ir),&
    u=>q(i,j,k,iu),&
    v=>q(i,j,k,iv),&
    w=>q(i,j,k,iw),&
    p=>q(i,j,k,ip),&
    A=>q(i,j,k,iA),&
    B=>q(i,j,k,iB),&
    C=>q(i,j,k,iC),&
    ! Face centered variables
    AL => bf(i,j,k,1),&
    BL => bf(i,j,k,2),&
    CL => bf(i,j,k,3),&
    AR => bf(i+1,j  ,k  ,1),&
    BR => bf(i  ,j+1,k  ,2),&
    CR => bf(i  ,j  ,k+1,3),&
    ! Edge centered electric field in X, Y and Z directions
    ELL => Ex(i,j  ,k  ),&
    ELR => Ex(i,j  ,k+1),&
    ERL => Ex(i,j+1,k  ),&
    ERR => Ex(i,j+1,k+1),&
    FLL => Ey(i  ,j,k  ),&
    FLR => Ey(i  ,j,k+1),&
    FRL => Ey(i+1,j,k  ),&
    FRR => Ey(i+1,j,k+1),&
    GLL => Ez(i  ,j  ,k),&
    GLR => Ez(i  ,j+1,k),&
    GRL => Ez(i+1,j  ,k),&
    GRR => Ez(i+1,j+1,k))
    !
    ! Cell centered TVD slopes in X, Y and Z directions
    drx = half * dq(i,j,k,ir,1)
    dux = half * dq(i,j,k,iu,1)
    dvx = half * dq(i,j,k,iv,1)
    dwx = half * dq(i,j,k,iw,1)
    dpx = half * dq(i,j,k,ip,1)
    dBx = half * dq(i,j,k,iB,1)
    dCx = half * dq(i,j,k,iC,1)
    !
    dry = half * dq(i,j,k,ir,2)
    duy = half * dq(i,j,k,iu,2)
    dvy = half * dq(i,j,k,iv,2)
    dwy = half * dq(i,j,k,iw,2)
    dpy = half * dq(i,j,k,ip,2)
    duy = half * dq(i,j,k,iu,2)
    dAy = half * dq(i,j,k,iA,2)
    dCy = half * dq(i,j,k,iC,2)
    !
    drz = half * dq(i,j,k,ir,3)
    duz = half * dq(i,j,k,iu,3)
    dvz = half * dq(i,j,k,iv,3)
    dwz = half * dq(i,j,k,iw,3)
    dpz = half * dq(i,j,k,ip,3)
    dAz = half * dq(i,j,k,iA,3)
    dBz = half * dq(i,j,k,iB,3)
    !
    ! Face centered TVD slopes in transverse directions
    dALy = half * dbf(i  ,j  ,k  ,1,1)
    dARy = half * dbf(i+1,j  ,k  ,1,1)
    dALz = half * dbf(i  ,j  ,k  ,1,2)
    dARz = half * dbf(i+1,j  ,k  ,1,2)

    dBLx = half * dbf(i  ,j  ,k  ,2,1)
    dBRx = half * dbf(i  ,j+1,k  ,2,1)
    dBLz = half * dbf(i  ,j  ,k  ,2,2)
    dBRz = half * dbf(i  ,j+1,k  ,2,2)

    dCLx = half * dbf(i  ,j  ,k  ,3,1)
    dCRx = half * dbf(i  ,j  ,k+1,3,1)
    dCLy = half * dbf(i  ,j  ,k  ,3,2)
    dCRy = half * dbf(i  ,j  ,k+1,3,2)
    !
    ! Face centered TVD slopes in transverse directions
    dALy = half * dbf(i  ,j  ,k  ,1,1)
    dARy = half * dbf(i+1,j  ,k  ,1,1)
    dALz = half * dbf(i  ,j  ,k  ,1,2)
    dARz = half * dbf(i+1,j  ,k  ,1,2)
    !
    dBLx = half * dbf(i  ,j  ,k  ,2,1)
    dBRx = half * dbf(i  ,j+1,k  ,2,1)
    dBLz = half * dbf(i  ,j  ,k  ,2,2)
    dBRz = half * dbf(i  ,j+1,k  ,2,2)
    !
    dCLx = half * dbf(i  ,j  ,k  ,3,1)
    dCRx = half * dbf(i  ,j  ,k+1,3,1)
    dCLy = half * dbf(i  ,j  ,k  ,3,2)
    dCRy = half * dbf(i  ,j  ,k+1,3,2)
    !
    ! Face-centered predicted states
    sAL0 = +(GLR-GLL)*dtdy*half -(FLR-FLL)*dtdz*half
    sAR0 = +(GRR-GRL)*dtdy*half -(FRR-FRL)*dtdz*half
    sBL0 = -(GRL-GLL)*dtdx*half +(ELR-ELL)*dtdz*half
    sBR0 = -(GRR-GLR)*dtdx*half +(ERR-ERL)*dtdz*half
    sCL0 = +(FRL-FLL)*dtdx*half -(ERL-ELL)*dtdy*half
    sCR0 = +(FRR-FLR)*dtdx*half -(ERR-ELR)*dtdy*half
    !
    AL = AL + sAL0
    AR = AR + sAR0
    BL = BL + sBL0
    BR = BR + sBR0
    CL = CL + sCL0
    CR = CR + sCR0
    !
    ! Source terms (including transverse derivatives)
    invr = one / r
    sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy + (-w*drz-dwz*r)*dtdz
    su0 = (-u*dux-(dpx+B*dBx+C*dCx)*invr)*dtdx + (-v*duy+B*dAy*invr)*dtdy + (-w*duz+C*dAz*invr)*dtdz
    sv0 = (-u*dvx+A*dBx*invr)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)*invr)*dtdy + (-w*dvz+C*dBz*invr)*dtdz
    sw0 = (-u*dwx+A*dCx*invr)*dtdx + (-v*dwy+B*dCy*invr)*dtdy + (-w*dwz-(dpz+A*dAz+B*dBz)*invr)*dtdz
    !
    ! Evolve entropy instead of pressure
    sp0 = (-u*dpx-dux*eff_gamma*p)*dtdx + (-v*dpy-dvy*eff_gamma*p)*dtdy + (-w*dpz-dwz*eff_gamma*p)*dtdz
    !sS0 = (-u*dSx)*dtdx + (-v*dSy)*dtdy + (-w*dSz)*dtdz
    !
!    ! Save predicted cell-centered pseudo-entropy; we are borrowing space in q(:,:,:,:,ip),
!    ! which is not used below this point.  In the calling routine this is again copied over
!    ! into uin, and is thus made available to godfine1, through that array.
!    if (entropy_fix) then
!       if ((i==1.or.i==2).and.(j==1.or.j==2).and.(k==1.or.k==2)) then
!          q(:,i,j,k,ip) = S + 2.0*sS0
!          if (local_debug) then
!             do ig=1,ngrid
!                if (in_debug_region(ind_grid(ig))) then
!                   print'(3hMKE,i8,i3,3f10.6,1p,12e11.3)', &
!                     ind_grid(ig),(i-1)+(j-1)*2+(k-1)*4, &
!                     xg(ind_grid(ig),1)+(i-1.5)/32.0, &
!                     xg(ind_grid(ig),2)+(j-1.5)/32.0, &
!                     xg(ind_grid(ig),3)+(k-1.5)/32.0, &
!                     log(p(ig)-dpx(ig)) - gamma*log(r(ig)-drx(ig)), &
!                     log(p(ig))-gamma*log(r(ig)), &
!                     log(p(ig)+dpx(ig)) - gamma*log(r(ig)+drx(ig)), &
!                     q(ig,i,j,k,ip), &
!                     dSx(ig),dSy(ig),dSz(ig), &
!                     u(ig),v(ig),w(ig),dtdx
!                end if
!             end do
!          end if
!       end if
!    end if
    !
    ! Cell-centered predicted states
    r = r + sr0
    r = max(r,smallr)
    u = u + su0
    v = v + sv0
    w = w + sw0
    p = p + sp0
    p = max(p,smallp)
    A = half*(AL+AR)
    B = half*(BL+BR)
    C = half*(CL+CR)
    !
    ! Face averaged right state at left interface
    qp(i,j,k,ir,1) = r - drx
    qp(i,j,k,iu,1) = u - dux
    qp(i,j,k,iv,1) = v - dvx
    qp(i,j,k,iw,1) = w - dwx
    qp(i,j,k,ip,1) = p - dpx
    qp(i,j,k,iA,1) = AL
    qp(i,j,k,iB,1) = B - dBx
    qp(i,j,k,iC,1) = C - dCx
    !
    ! Face averaged left state at right interface
    qm(i,j,k,ir,1) = r + drx
    qm(i,j,k,iu,1) = u + dux
    qm(i,j,k,iv,1) = v + dvx
    qm(i,j,k,iw,1) = w + dwx
    qm(i,j,k,ip,1) = p + dpx
    qm(i,j,k,iA,1) = AR
    qm(i,j,k,iB,1) = B + dBx
    qm(i,j,k,iC,1) = C + dCx
    !
    ! Face averaged top state at bottom interface
    qp(i,j,k,ir,2) = r - dry
    qp(i,j,k,iu,2) = u - duy
    qp(i,j,k,iv,2) = v - dvy
    qp(i,j,k,iw,2) = w - dwy
    qp(i,j,k,ip,2) = p - dpy
    qp(i,j,k,iA,2) = A - dAy
    qp(i,j,k,iB,2) = BL
    qp(i,j,k,iC,2) = C - dCy
    !
    ! Face averaged bottom state at top interface
    qm(i,j,k,ir,2) = r + dry
    qm(i,j,k,iu,2) = u + duy
    qm(i,j,k,iv,2) = v + dvy
    qm(i,j,k,iw,2) = w + dwy
    qm(i,j,k,ip,2) = p + dpy
    qm(i,j,k,iA,2) = A + dAy
    qm(i,j,k,iB,2) = BR
    qm(i,j,k,iC,2) = C + dCy
    !
    ! Face averaged front state at back interface
    qp(i,j,k,ir,3) = r - drz
    qp(i,j,k,iu,3) = u - duz
    qp(i,j,k,iv,3) = v - dvz
    qp(i,j,k,iw,3) = w - dwz
    qp(i,j,k,ip,3) = p - dpz
    qp(i,j,k,iA,3) = A - dAz
    qp(i,j,k,iB,3) = B - dBz
    qp(i,j,k,iC,3) = CL
    !
    ! Face averaged back state at front interface
    qm(i,j,k,ir,3) = r + drz
    qm(i,j,k,iu,3) = u + duz
    qm(i,j,k,iv,3) = v + dvz
    qm(i,j,k,iw,3) = w + dwz
    qm(i,j,k,ip,3) = p + dpz
    qm(i,j,k,iA,3) = A + dAz
    qm(i,j,k,iB,3) = B + dBz
    qm(i,j,k,iC,3) = CR
    !
    ! X-edge averaged right-top corner state (RT->LL)
    qRT(i,j,k,ir,1) = r + (+dry+drz)
    qRT(i,j,k,iu,1) = u + (+duy+duz)
    qRT(i,j,k,iv,1) = v + (+dvy+dvz)
    qRT(i,j,k,iw,1) = w + (+dwy+dwz)
    qRT(i,j,k,ip,1) = p + (+dpy+dpz)
    qRT(i,j,k,iA,1) = A + (+dAy+dAz)
    qRT(i,j,k,iB,1) = BR+ (   +dBRz)
    qRT(i,j,k,iC,1) = CR+ (+dCRy   )
    !
    ! X-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,1) = r + (+dry-drz)
    qRB(i,j,k,iu,1) = u + (+duy-duz)
    qRB(i,j,k,iv,1) = v + (+dvy-dvz)
    qRB(i,j,k,iw,1) = w + (+dwy-dwz)
    qRB(i,j,k,ip,1) = p + (+dpy-dpz)
    qRB(i,j,k,iA,1) = A + (+dAy-dAz)
    qRB(i,j,k,iB,1) = BR+ (   -dBRz)
    qRB(i,j,k,iC,1) = CL+ (+dCLy   )
    !
    ! X-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,1) = r + (-dry+drz)
    qLT(i,j,k,iu,1) = u + (-duy+duz)
    qLT(i,j,k,iv,1) = v + (-dvy+dvz)
    qLT(i,j,k,iw,1) = w + (-dwy+dwz)
    qLT(i,j,k,ip,1) = p + (-dpy+dpz)
    qLT(i,j,k,iA,1) = A + (-dAy+dAz)
    qLT(i,j,k,iB,1) = BL+ (   +dBLz)
    qLT(i,j,k,iC,1) = CR+ (-dCRy   )
    !
    ! X-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,1) = r + (-dry-drz)
    qLB(i,j,k,iu,1) = u + (-duy-duz)
    qLB(i,j,k,iv,1) = v + (-dvy-dvz)
    qLB(i,j,k,iw,1) = w + (-dwy-dwz)
    qLB(i,j,k,ip,1) = p + (-dpy-dpz)
    qLB(i,j,k,iA,1) = A + (-dAy-dAz)
    qLB(i,j,k,iB,1) = BL+ (   -dBLz)
    qLB(i,j,k,iC,1) = CL+ (-dCLy   )
    !
    ! Y-edge averaged right-top corner state (RT->LL)
    qRT(i,j,k,ir,2) = r + (+drx+drz)
    qRT(i,j,k,iu,2) = u + (+dux+duz)
    qRT(i,j,k,iv,2) = v + (+dvx+dvz)
    qRT(i,j,k,iw,2) = w + (+dwx+dwz)
    qRT(i,j,k,ip,2) = p + (+dpx+dpz)
    qRT(i,j,k,iA,2) = AR+ (   +dARz)
    qRT(i,j,k,iB,2) = B + (+dBx+dBz)
    qRT(i,j,k,iC,2) = CR+ (+dCRx   )
    !
    ! Y-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,2) = r + (+drx-drz)
    qRB(i,j,k,iu,2) = u + (+dux-duz)
    qRB(i,j,k,iv,2) = v + (+dvx-dvz)
    qRB(i,j,k,iw,2) = w + (+dwx-dwz)
    qRB(i,j,k,ip,2) = p + (+dpx-dpz)
    qRB(i,j,k,iA,2) = AR+ (   -dARz)
    qRB(i,j,k,iB,2) = B + (+dBx-dBz)
    qRB(i,j,k,iC,2) = CL+ (+dCLx   )
    !
    ! Y-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,2) = r + (-drx+drz)
    qLT(i,j,k,iu,2) = u + (-dux+duz)
    qLT(i,j,k,iv,2) = v + (-dvx+dvz)
    qLT(i,j,k,iw,2) = w + (-dwx+dwz)
    qLT(i,j,k,ip,2) = p + (-dpx+dpz)
    qLT(i,j,k,iA,2) = AL+ (   +dALz)
    qLT(i,j,k,iB,2) = B + (-dBx+dBz)
    qLT(i,j,k,iC,2) = CR+ (-dCRx   )
    !
    ! Y-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,2) = r + (-drx-drz)
    qLB(i,j,k,iu,2) = u + (-dux-duz)
    qLB(i,j,k,iv,2) = v + (-dvx-dvz)
    qLB(i,j,k,iw,2) = w + (-dwx-dwz)
    qLB(i,j,k,ip,2) = p + (-dpx-dpz)
    qLB(i,j,k,iA,2) = AL+ (   -dALz)
    qLB(i,j,k,iB,2) = B + (-dBx-dBz)
    qLB(i,j,k,iC,2) = CL+ (-dCLx   )
    !
    ! Z-edge averaged right-top corner state (RT->LL)
    qRT(i,j,k,ir,3) = r + (+drx+dry)
    qRT(i,j,k,iu,3) = u + (+dux+duy)
    qRT(i,j,k,iv,3) = v + (+dvx+dvy)
    qRT(i,j,k,iw,3) = w + (+dwx+dwy)
    qRT(i,j,k,ip,3) = p + (+dpx+dpy)
    qRT(i,j,k,iA,3) = AR+ (   +dARy)
    qRT(i,j,k,iB,3) = BR+ (+dBRx   )
    qRT(i,j,k,iC,3) = C + (+dCx+dCy)
    !
    ! Z-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,3) = r + (+drx-dry)
    qRB(i,j,k,iu,3) = u + (+dux-duy)
    qRB(i,j,k,iv,3) = v + (+dvx-dvy)
    qRB(i,j,k,iw,3) = w + (+dwx-dwy)
    qRB(i,j,k,ip,3) = p + (+dpx-dpy)
    qRB(i,j,k,iA,3) = AR+ (   -dARy)
    qRB(i,j,k,iB,3) = BL+ (+dBLx   )
    qRB(i,j,k,iC,3) = C + (+dCx-dCy)
    !
    ! Z-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,3) = r + (-drx+dry)
    qLT(i,j,k,iu,3) = u + (-dux+duy)
    qLT(i,j,k,iv,3) = v + (-dvx+dvy)
    qLT(i,j,k,iw,3) = w + (-dwx+dwy)
    qLT(i,j,k,ip,3) = p + (-dpx+dpy)
    qLT(i,j,k,iA,3) = AL+ (   +dALy)
    qLT(i,j,k,iB,3) = BR+ (-dBRx   )
    qLT(i,j,k,iC,3) = C + (-dCx+dCy)
    !
    ! Z-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,3) = r + (-drx-dry)
    qLB(i,j,k,iu,3) = u + (-dux-duy)
    qLB(i,j,k,iv,3) = v + (-dvx-dvy)
    qLB(i,j,k,iw,3) = w + (-dwx-dwy)
    qLB(i,j,k,ip,3) = p + (-dpx-dpy)
    qLB(i,j,k,iA,3) = AL+ (   -dALy)
    qLB(i,j,k,iB,3) = BL+ (-dBLx   )
    qLB(i,j,k,iC,3) = C + (-dCx-dCy)
    !
    qRT(i,j,k,ir,1) = r + (+dry+drz)
    qRB(i,j,k,ir,1) = r + (+dry-drz)
    qLT(i,j,k,ir,1) = r + (-dry+drz)
    qLB(i,j,k,ir,1) = r + (-dry-drz)
    qRT(i,j,k,ir,2) = r + (+drx+drz)
    qRB(i,j,k,ir,2) = r + (+drx-drz)
    qLT(i,j,k,ir,2) = r + (-drx+drz)
    qLB(i,j,k,ir,2) = r + (-drx-drz)
    qRT(i,j,k,ir,3) = r + (+drx+dry)
    qRB(i,j,k,ir,3) = r + (+drx-dry)
    qLT(i,j,k,ir,3) = r + (-drx+dry)
    qLB(i,j,k,ir,3) = r + (-drx-dry)
    end associate
  end do
  end do
  end do
  n=size(q,1)
  qp(n,:,:,:,:)=qp(n-1,:,:,:,:)
  qm(n,:,:,:,:)=qm(n-1,:,:,:,:)
  n=size(q,2)
  qp(:,n,:,:,:)=qp(:,n-1,:,:,:)
  qm(:,n,:,:,:)=qm(:,n-1,:,:,:)
  n=size(q,3)
  qp(:,:,n,:,:)=qp(:,:,n-1,:,:)
  qm(:,:,n,:,:)=qm(:,:,n-1,:,:)
  !
  call trace%end (itimer)
  end associate
contains
  FUNCTION xdn(a) RESULT (b)
    real, dimension(:,:,:):: a
    real, dimension(size(a,1),size(a,2),size(a,3)):: b
    do k=1,size(a,3)
    do j=1,size(a,2)
    do i=2,size(a,1)
      b(i,j,k)=half*(a(i-1,j,k)+a(i,j,k))
    end do
    end do
    end do
    b(1,:,:) = a(1,:,:)
  END FUNCTION
  FUNCTION ydn(a) RESULT (b)
    real, dimension(:,:,:):: a
    real, dimension(size(a,1),size(a,2),size(a,3)):: b
    do k=1,size(a,3)
    do j=2,size(a,2)
    do i=1,size(a,1)
      b(i,j,k)=half*(a(i,j-1,k)+a(i,j,k))
    end do
    end do
    end do
    b(:,1,:) = a(:,1,:)
  END FUNCTION
  FUNCTION zdn(a) RESULT (b)
    real, dimension(:,:,:):: a
    real, dimension(size(a,1),size(a,2),size(a,3)):: b
    do k=2,size(a,3)
    do j=1,size(a,2)
    do i=1,size(a,1)
      b(i,j,k)=half*(a(i,j,k-1)+a(i,j,k))
    end do
    end do
    end do
    b(:,:,1) = a(:,:,1)
  END FUNCTION
END SUBROUTINE trace3d

!> =============================================================================
!> Compute fluxes
!> =============================================================================
subroutine cmpflxm(self, nx, ln, lt1, lt2, bn, bt1, bt2, fl)
  class(mhd_t):: self
  real(dp), dimension(:,:,:,:):: fl
  integer:: nx, ln, lt1, lt2, bn, bt1, bt2
  ! local variables
  integer:: i, j, k, lb(3), ub(3), xdim
  real(dp), dimension(nx,self%nv):: qleft, qright, fgdnv
  real(dp), dimension(nx):: bn_mean
  real(dp):: dtds(3)
  integer, save:: itimer=0
  associate (qp=>self%qp, qm=>self%qm)
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_mod::cmpflxm',itimer=itimer)
  dtds = self%dtime/self%mesh%d
  lb = self%mesh%li
  ub = self%mesh%uo
  xdim = ln-1
  !
  do k = lb(3),ub(3)
  do j = lb(2),ub(2)
    ! Enforce continuity for normal magnetic field
    bn_mean(1:nx) = half*(qm(1:nx,j,k,bn,xdim)+qp(1:nx,j,k,bn,xdim))
    ! Left state
    qleft (1:nx,1) = qm(1:nx,j,k,1  ,xdim) ! Mass density
    qleft (1:nx,2) = qm(1:nx,j,k,5  ,xdim) ! Pressure
    qleft (1:nx,3) = qm(1:nx,j,k,ln ,xdim) ! Normal velocity
    qleft (1:nx,4) = bn_mean(1:nx)         ! Normal magnetic field
    qleft (1:nx,5) = qm(1:nx,j,k,lt1,xdim) ! Tangential velocity 1
    qleft (1:nx,6) = qm(1:nx,j,k,bt1,xdim) ! Tangential magnetic field 1
    qleft (1:nx,7) = qm(1:nx,j,k,lt2,xdim) ! Tangential velocity 2
    qleft (1:nx,8) = qm(1:nx,j,k,bt2,xdim) ! Tangential magnetic field 2
    ! Right state
    qright(1:nx,1) = qp(1:nx,j,k,1  ,xdim) ! Mass density
    qright(1:nx,2) = qp(1:nx,j,k,5  ,xdim) ! Pressure
    qright(1:nx,3) = qp(1:nx,j,k,ln ,xdim) ! Normal velocity
    qright(1:nx,4) = bn_mean(1:nx)         ! Normal magnetic field
    qright(1:nx,5) = qp(1:nx,j,k,lt1,xdim) ! Tangential velocity 1
    qright(1:nx,6) = qp(1:nx,j,k,bt1,xdim) ! Tangential magnetic field 1
    qright(1:nx,7) = qp(1:nx,j,k,lt2,xdim) ! Tangential velocity 2
    qright(1:nx,8) = qp(1:nx,j,k,bt2,xdim) ! Tangential magnetic field 2
    call solver (self,qleft,qright,fgdnv)
    ! Output fluxes
    fl(1:nx,j,k,1  ) = fgdnv(1:nx,1)  ! Mass density
    fl(1:nx,j,k,5  ) = fgdnv(1:nx,2)  ! Total energy
    fl(1:nx,j,k,ln ) = fgdnv(1:nx,3)  ! Normal momentum
    fl(1:nx,j,k,bn ) = fgdnv(1:nx,4)  ! Normal magnetic field
    fl(1:nx,j,k,lt1) = fgdnv(1:nx,5)  ! Transverse momentum 1
    fl(1:nx,j,k,bt1) = fgdnv(1:nx,6)  ! Transverse magnetic field 1
    fl(1:nx,j,k,lt2) = fgdnv(1:nx,7)  ! Transverse momentum 2
    fl(1:nx,j,k,bt2) = fgdnv(1:nx,8)  ! Transverse magnetic field 2
  end do
  end do
  end associate
  call trace%end (itimer)
end subroutine cmpflxm

!> =============================================================================
!> 2D solver for EMF
!> =============================================================================
subroutine cmp_mag_flx (self, lp1 ,lp2 ,lor ,bp1 , bp2, bor, dtds, emf, ngrid, nvar)
  class(mhd_t):: self
  ! 2D Riemann solver to compute EMF at cell edges
  ! indices of the 2 planar velocity lp1 and lp2 and the orthogonal one,
  ! lor and idem for the magnetic field
  integer:: lp1, lp2, lor, bp1, bp2, bor, ngrid, nvar
  real(dp):: dtds
  real(dp),DIMENSION(:,:,:):: emf
  ! local variables
  integer:: j, k, n, xdim, l
  real(dp),DIMENSION(ngrid,nvar)::qLL,qRL,qLR,qRR
  real(dp):: smallp
  integer, save:: itimer=0
  call trace%begin ('mhd_mod::cmp_mag_flx',itimer=itimer)
  associate (qRT=>self%qRT,qLT=>self%qLT,qRB=>self%qRB,qLB=>self%qLB)
  !-----------------------------------------------------------------------------

  xdim = lor - 1
  smallp = smallr*smallc**2
  !
  DO k = 1,self%mesh(3)%ub
     DO j = 1,self%mesh(2)%ub
      qLL (1:ngrid,1) = max(qRT(1:ngrid,j,k,1,xdim),smallr)
      qRL (1:ngrid,1) = max(qLT(1:ngrid,j,k,1,xdim),smallr)
      qLR (1:ngrid,1) = max(qRB(1:ngrid,j,k,1,xdim),smallr)
      qRR (1:ngrid,1) = max(qLB(1:ngrid,j,k,1,xdim),smallr)
      qLL (1:ngrid,2) = max(qRT(1:ngrid,j,k,5,xdim),smallp)
      qRL (1:ngrid,2) = max(qLT(1:ngrid,j,k,5,xdim),smallp)
      qLR (1:ngrid,2) = max(qRB(1:ngrid,j,k,5,xdim),smallp)
      qRR (1:ngrid,2) = max(qLB(1:ngrid,j,k,5,xdim),smallp)
      qLL (1:ngrid,3) = qRT(1:ngrid,j,k,lp1,xdim)
      qRL (1:ngrid,3) = qLT(1:ngrid,j,k,lp1,xdim)
      qLR (1:ngrid,3) = qRB(1:ngrid,j,k,lp1,xdim)
      qRR (1:ngrid,3) = qLB(1:ngrid,j,k,lp1,xdim)
      qLL (1:ngrid,4) = qRT(1:ngrid,j,k,lp2,xdim)
      qRL (1:ngrid,4) = qLT(1:ngrid,j,k,lp2,xdim)
      qLR (1:ngrid,4) = qRB(1:ngrid,j,k,lp2,xdim)
      qRR (1:ngrid,4) = qLB(1:ngrid,j,k,lp2,xdim)
      qLL (1:ngrid,6) = half*(qRT(1:ngrid,j,k,bp1,xdim)+qLT(1:ngrid,j,k,bp1,xdim))
      qRL (1:ngrid,6) = half*(qRT(1:ngrid,j,k,bp1,xdim)+qLT(1:ngrid,j,k,bp1,xdim))
      qLR (1:ngrid,6) = half*(qRB(1:ngrid,j,k,bp1,xdim)+qLB(1:ngrid,j,k,bp1,xdim))
      qRR (1:ngrid,6) = half*(qRB(1:ngrid,j,k,bp1,xdim)+qLB(1:ngrid,j,k,bp1,xdim))
      qLL (1:ngrid,7) = half*(qRT(1:ngrid,j,k,bp2,xdim)+qRB(1:ngrid,j,k,bp2,xdim))
      qRL (1:ngrid,7) = half*(qLT(1:ngrid,j,k,bp2,xdim)+qLB(1:ngrid,j,k,bp2,xdim))
      qLR (1:ngrid,7) = half*(qRT(1:ngrid,j,k,bp2,xdim)+qRB(1:ngrid,j,k,bp2,xdim))
      qRR (1:ngrid,7) = half*(qLT(1:ngrid,j,k,bp2,xdim)+qLB(1:ngrid,j,k,bp2,xdim))
      qLL (1:ngrid,5) = qRT(1:ngrid,j,k,lor,xdim)
      qRL (1:ngrid,5) = qLT(1:ngrid,j,k,lor,xdim)
      qLR (1:ngrid,5) = qRB(1:ngrid,j,k,lor,xdim)
      qRR (1:ngrid,5) = qLB(1:ngrid,j,k,lor,xdim)
      qLL (1:ngrid,8) = qRT(1:ngrid,j,k,bor,xdim)
      qRL (1:ngrid,8) = qLT(1:ngrid,j,k,bor,xdim)
      qLR (1:ngrid,8) = qRB(1:ngrid,j,k,bor,xdim)
      qRR (1:ngrid,8) = qLB(1:ngrid,j,k,bor,xdim)
      call rieman%hlld_2d (qLL,qRL,qLR,qRR,emf,j,k,ngrid,nvar)
    END DO
  END DO
  call trace%end (itimer)
  end associate
END SUBROUTINE cmp_mag_flx

!> =============================================================================
!> Selection of 1-D solvers
!> =============================================================================
subroutine solver (self,qleft,qright,fgdnv)
  class(mhd_t):: self
  real(dp),dimension(:,:):: qleft,qright,fgdnv
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mhd_mod::solver', itimer=itimer)
  if (self%riemann.eq.'hlld')then
    call rieman%hlld_vec(qleft,qright,fgdnv)
  else
    write(*,*)'unknown Riemann solver'
    stop
  end if
  call trace%end (itimer)
end subroutine solver

!> =============================================================================
!> =============================================================================
FUNCTION gas_pressure (self)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: gas_pressure
  gas_pressure = self%mem(:,:,:,self%idx%e,self%it,1)*(self%gamma-1.0)
END FUNCTION gas_pressure

!> =============================================================================
!> =============================================================================
SUBROUTINE gravitational_force (self)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: phi, d, et, gradphi1, gradphi2, gradphi3
  real, dimension(:,:,:,:), pointer:: p
  integer:: i,j,k
  real, dimension(3):: wtp1, wtm1
  real:: etemp
  !.............................................................................
  associate (l => self%li, u => self%ui)
  wtp1 = 0.50 * 4.0 / 3.0 / self%ds
  wtm1 = 0.25 * 1.0 / 3.0 / self%ds
  phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  d   => self%mem(:,:,:,self%idx%d,self%it,1)
  et  => self%mem(:,:,:,self%idx%e,self%new,1)
  p   => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  allocate (gradphi1(self%gn(1),self%gn(2),self%gn(3)), &
            gradphi2(self%gn(1),self%gn(2),self%gn(3)), &
            gradphi3(self%gn(1),self%gn(2),self%gn(3)))

  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        ! calculate energy minus the kinetic energy
        etemp = et(i,j,k) - 0.5 * p(i,j,k,1)**2 / d(i,j,k) &
                          - 0.5 * p(i,j,k,2)**2 / d(i,j,k) &
                          - 0.5 * p(i,j,k,3)**2 / d(i,j,k)
        ! calculate gradient of phi
        gradphi1(i,j,k) = wtp1(1) * (phi(i+1,j,k) - phi(i-1,j,k)) &
                 - wtm1(1) * (phi(i+2,j,k) - phi(i-2,j,k))
        gradphi2(i,j,k) = wtp1(2) * (phi(i,j+1,k) - phi(i,j-1,k)) &
                 - wtm1(2) * (phi(i,j+2,k) - phi(i,j-2,k))
        gradphi3(i,j,k) = wtp1(3) * (phi(i,j,k+1) - phi(i,j,k-1)) &
                 - wtm1(3) * (phi(i,j,k+2) - phi(i,j,k-2))
        ! update momentum
        p(i,j,k,1) = p(i,j,k,1) - self%dtime * d(i,j,k) * gradphi1(i,j,k)
        p(i,j,k,2) = p(i,j,k,2) - self%dtime * d(i,j,k) * gradphi2(i,j,k)
        p(i,j,k,3) = p(i,j,k,3) - self%dtime * d(i,j,k) * gradphi3(i,j,k)
        ! update total energy with new momentum
        et(i,j,k) = etemp + 0.5 * p(i,j,k,1)**2 / d(i,j,k) &
                          + 0.5 * p(i,j,k,2)**2 / d(i,j,k) &
                          + 0.5 * p(i,j,k,3)**2 / d(i,j,k)
      end do
    end do
  end do

  deallocate (gradphi1, gradphi2, gradphi3)
  end associate
END SUBROUTINE gravitational_force

END MODULE
