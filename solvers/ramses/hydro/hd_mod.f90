!===============================================================================
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
!>
!> Note that it is NECESSARY to have a call to courant_condition here, after
!> ctoprim (which calculates an up-to-date self%u_max), and before trace3d,
!> which needs the new self%dtime
!===============================================================================
MODULE hd_mod
  !.............................................................................
  USE const
  USE io_mod
  USE io_unit_mod
  USE extras_mod
  USE link_mod
  USE riemann_mod
  USE hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE eos_mod
  USE index_mod
  USE patch_mod
  USE validate_mod
  implicit none
  PRIVATE
  type, public, extends(extras_t):: hd_t
    logical:: mhd=.false.
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
    procedure:: dealloc
    procedure:: update
    procedure:: gas_pressure
  end type
  integer, save:: courant_type=0
  logical, save:: debug=.false.
  logical, save:: first_time=.true.
  logical, save:: detailed_timer=.false.
  logical, save:: unsigned=.true.
  character(len=16):: eqn='mhd'
CONTAINS

!===============================================================================
!> Initialize the HD solver
!===============================================================================
SUBROUTINE init (self)
  class(hd_t):: self
  character(len=10):: riemann='hllc'
  real, save:: csound=1.
  integer:: i, iv
  namelist /ramses_params/ gamma, riemann, slope_type, courant_factor, smallr, &
    smallc, isothermal, csound, detailed_timer, unsigned, courant_type, debug
  type(index_t):: idx
  character(len=120):: id = &
   '$Id$ solvers/ramses/hydro/hd_mod.f90'
!...............................................................................
  call trace%begin('hd_mod%init')
  call trace%print_id (id)
  self%mhd = .false.
  self%nw = 1
  self%ng = 2
  self%kind = 'ramses_hd_patch'
  if (self%nv == 0) self%nv = 5
  call self%idx%init (5, self%mhd)
  call self%gpatch_t%init
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, ramses_params)
    if (mpi%rank==0) write (*, ramses_params)
  end if
  if (gamma == 1d0) &
    isothermal = .true.
  if (isothermal) &
    gamma = 1d0
  !$omp end critical (input_cr)
  self%riemann = riemann
  self%courant = courant_factor
  self%gamma = gamma
  self%csound = csound
  self%staggered = .false.
  do i=1,3
    self%mesh(i)%h = 0d0
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Deallocate permanently allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(hd_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%dealloc')
  call self%extras_t%dealloc
  call trace%end
END SUBROUTINE dealloc

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(hd_t):: self
  integer:: n(3), m, i, l(3), u(3)
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: p
  real, dimension(:,:,:), pointer:: d, e, Ux=>null(), Uy=>null(), Uz=>null()
  real, dimension(:,:,:), pointer:: tmp
  real:: error
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%update', itimer=itimer)
  if (io%smallr  /= 0.0) smallr  = io%smallr
  if (io%courant /= 0.0) self%courant = io%courant
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  !-----------------------------------------------------------------------------
  ! Memory mapping and allocations
  !-----------------------------------------------------------------------------
  self%uin => self%mem(:,:,:,:,self%it,1)
  self%unew => self%mem(:,:,:,:,self%new,1)
  n = self%mesh%gn
  m = n(1)
  allocate (self%q (n(1),n(2),n(3),self%nv))
  allocate (self%qp(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%qm(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%dq(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%fl(n(1),n(2),n(3),self%nv,ndim))
  allocate (self%c(n(1),n(2),n(3)))
  allocate (self%s(n(1),n(2),n(3)))
  !self%q=0.0; self%qp=0.0; self%qm=0.0; self%dq=0.0; self%fl=0.; self%c=0.0; self%s=0.0
  if (detailed_timer) call trace_end (itimer)
  if (debug) then
    if (self%id==229.or.self%id==231) then
      print *,self%id, 'mk1:px', self%mem(19-2:19+2,19,19,self%idx%px,self%it,1), self%idx%px
      print *,self%id, 'mk1: d', self%mem(19-2:19+2,19,19,self%idx%d ,self%it,1)
      print *,self%id, 'mk1: e', self%mem(19-2:19+2,19,19,self%idx%e ,self%it,1), self%idx%e
    end if
    if (self%id==230) then
      print *,self%id, 'mk1:px', self%mem( m-7:m-3 ,19,19,self%idx%px,self%it,1), self%idx%px
      print *,self%id, 'mk1: d', self%mem( m-7:m-3 ,19,19,self%idx%d ,self%it,1)
      print *,self%id, 'mk1: e', self%mem( m-7:m-3 ,19,19,self%idx%e ,self%it,1), self%idx%e
    end if
  end if
  d => self%mem(:,:,:,          self%idx%d ,self%it,1)
  e => self%mem(:,:,:,          self%idx%e ,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  ! Godunov solver
  !-----------------------------------------------------------------------------
  call godunov (self)
  if (isothermal) &
    self%mem(:,:,:,self%idx%e,self%new,1) = &
    self%mem(:,:,:,self%idx%d,self%new,1)*self%csound**2
  if (detailed_timer) call trace_begin ('hd_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Add the external force per unit mass as a source term
  !-----------------------------------------------------------------------------
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  if (allocated(self%force_per_unit_mass)) then
    do i=1,3
      p(:,:,:,i) = p(:,:,:,i) + self%dtime*self%force_per_unit_mass(:,:,:,i)*d
    end do
  end if
  if (allocated(self%force_per_unit_volume)) then
    p = p + self%dtime*self%force_per_unit_volume
  end if
  call validate%check (self%link, p, 'after force')
  if (omp_lock%tasks) then
    call self%lock%unset ('update')
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate
  !-----------------------------------------------------------------------------
  call self%counter_update
  deallocate (self%q,self%qp,self%qm,self%dq,self%fl,self%c,self%s)
  call trace%end (itimer)
END SUBROUTINE update

!=======================================================================
!> Get velocities from momentum
!=======================================================================
SUBROUTINE p2U(self, U, it)
  class(hd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d,           it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,size(p,4)
  U(:,:,:,i) = p(:,:,:,i)/d
  end do
  end associate
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momentum
!=======================================================================
SUBROUTINE U2p(self, U, it)
  class(hd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  integer:: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  do i=1,size(p,4)
  p(:,:,:,i)= U(:,:,:,i)*d
  end do
  end associate
END SUBROUTINE U2p

!=======================================================================
!> Get thermal energy from total/specific energy
!=======================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(hd_t) :: self
  real, dimension(:,:,:), pointer:: E_th
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: k
  integer :: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            e => self%mem(:,:,:,self%idx%e           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  E_th = (e-k)/d
  deallocate(k)
  end associate
END SUBROUTINE e2E_th

!=======================================================================
!> Get etropy per unit mass from total energy
!=======================================================================
SUBROUTINE e2S(self, s, it)
  class(hd_t) :: self
  real, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: k
  integer :: i
  !---------------------------------------------------------------------
  if (isothermal) then
    s = 0.0
    return
  end if
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            e => self%mem(:,:,:,self%idx%e           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  s = log((e-k)*(gamma-1d0)/d**gamma)/(gamma-1d0)
  deallocate(k)
  end associate
END SUBROUTINE e2S

!=======================================================================
!> Put etropy per unit mass to total energy
!=======================================================================
SUBROUTINE S2e (self, s, it)
  class(hd_t) :: self
  real, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: k
  integer :: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            e => self%mem(:,:,:,self%idx%e           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  if (isothermal) then
    e = d*self%csound**2
    return
  end if
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  e = d**gamma*exp(s*(gamma-1d0)) + k
  deallocate(k)
  end associate
END SUBROUTINE S2e

!=======================================================================
!> Put thermal energy to total/specific energy
!=======================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(hd_t) :: self
  integer:: it
  real, dimension(:,:,:), pointer:: E_th
  !.....................................................................
  real, dimension(:,:,:), pointer:: k
  integer :: i
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%idx%d           ,it,1), &
            e => self%mem(:,:,:,self%idx%e           ,it,1), &
            p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1))
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  e = E_th*d+k
  deallocate(k)
  end associate
END SUBROUTINE E_th2e

!===============================================================================
subroutine godunov (self)
  class(hd_t)::self
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
  if (.not.detailed_timer) call trace%begin ('hd_t%godunov')
  li = self%mesh%li
  ui = self%mesh%ui
  call check_small (self, uin)
  call unsplit (self)
  !-----------------------------------------------------------------------------
  ! Conservative update, looping over internal cells only
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('hd_t%godunov', itimer=itimer)
  if (omp_lock%tasks) then
    call self%lock%set ('update')
  end if
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
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end
  end if
  end associate
end subroutine godunov

!===============================================================================
!===============================================================================
subroutine check_small (self, uu)
  class(hd_t):: self
  real, dimension(:,:,:,:):: uu
  integer:: i,j,k,lo(3),uo(3),n_smalld,n_smallp
  real(dp):: c,d,u,v,w,p,ekin
  integer, save:: nprint=10
  integer, save:: itimer=0
  integer:: it, jt, iw
  logical:: panic
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('hd_t%check_small', itimer=itimer)
  n_smalld = 0
  n_smallp = 0
  panic = .false.
  lo = self%mesh%lb
  uo = self%mesh%ub
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  do i=lo(1),uo(1)
    if (uu(i,j,k,self%idx%d) < smallr) then
      panic = .true.
      if (nprint>0) then
        nprint = nprint-1
        write (io_unit%log,'(i7,3i5,1p,2e12.3)') self%id,i,j,k,uu(i,j,k,self%idx%d),smallr
      end if
      n_smalld = n_smalld+1
      uu(i,j,k,self%idx%d) = smallr
    end if
    d =  uu(i,j,k,self%idx%d)
    c = 1./d
    u =  uu(i,j,k,self%idx%px)*c
    v =  uu(i,j,k,self%idx%py)*c
    w =  uu(i,j,k,self%idx%pz)*c
    ekin = 0.5*d*(u**2+v**2+w**2)
    if (isothermal) then
      uu(i,j,k,self%idx%e) = max(uu(i,j,k,self%idx%e), smallr*smallc**2)
    else
      p = (uu(i,j,k,self%idx%e)-ekin)*(gamma-one)
      if (p < smallr*smallc**2) then
        p = smallr*smallc**2
        uu(i,j,k,self%idx%e) = p/(gamma-one)+ekin
        n_smallp = n_smallp+1
      end if
    end if
  end do
  end do
  end do
  if (n_smalld+n_smallp > 0) then
    write (io_unit%log,*) 'small: id,nr,np', self%id,n_smalld,n_smallp
  end if
  if (panic .and. self%n_dump > 0) then
    !$omp critical (panic_cr)
    self%n_dump = self%n_dump-1
    write (io_unit%dump) self%id, self%gn, self%nv, self%nt, self%nw, self%istep
    write (io_unit%dump) self%time, self%dtime
    do iw=1,self%nw
    do it=1,self%nt
      !$omp atomic read
      jt = self%iit(it)
      write (io_unit%dump) self%mem(:,:,:,:,jt,iw)
    end do
    end do
    write (io_unit%mpi,*) 'density < SMALLR, dumped', self%id, self%time
    !$omp end critical (panic_cr)
  end if
  if (detailed_timer) call trace%end (itimer)
end subroutine check_small

!===============================================================================
!>  UNSPLIT     Unsplit second order Godunov integrator for
!>              polytropic gas dynamics using MUSCL-HANCOCK
!===============================================================================
subroutine unsplit(self)
  class(hd_t):: self
  real(dp):: dtds(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%unsplit')
  call ctoprim(self)                  ! Primitive variables, sound speeds
  call self%courant_condition (detailed_timer)
  call uslope(self)                   ! TVD slopes
  call trace3d(self)                  ! predictor step
  call cmpflxm(self)                  ! Solve for 1D flux in X direction
  call trace%end
end subroutine unsplit

!===============================================================================
!> Conservative (uin) to primitive (q) variables, stored in the order
!> 1:d, 2:u, 3:v, 4:w, 5:p
!===============================================================================
subroutine ctoprim (self)
  class(hd_t):: self
  real(dp):: smallp, dtxhalf
  integer,save:: itimer=0
  integer:: l(3), u(3)
  real, dimension(:,:,:), pointer:: d, e, pg
  associate (uin=>self%uin, q=>self%q, c=>self%c, s=>self%s)
  !-----------------------------------------------------------------------------
  if (detailed_timer) then
    call trace%begin ('hd_mod::ctoprim', itimer=itimer)
  else
    call trace%begin ('hd_mod::ctoprim')
  end if
  smallp = smallr*smallc**2/gamma
  q(:,:,:,1) = max(uin(:,:,:,self%idx%d),smallr)
  s = 1.d0/q(:,:,:,1)
  q(:,:,:,2) = uin(:,:,:,self%idx%px)*s
  q(:,:,:,3) = uin(:,:,:,self%idx%py)*s
  q(:,:,:,4) = uin(:,:,:,self%idx%pz)*s
  s = q(:,:,:,2)**2+q(:,:,:,3)**2+q(:,:,:,4)**2
  !-----------------------------------------------------------------------------
  ! Equation-of-state: s is the internal energy, and pg is the gas pressure
  !-----------------------------------------------------------------------------
  s = max((uin(:,:,:,self%idx%e)-half*q(:,:,:,1)*s(:,:,:)),smallp)
  d => q(:,:,:,1)
  e => s
  pg => q(:,:,:,5)
  if (isothermal) then
    pg = d*self%csound**2
  else
    pg = e*(gamma-1d0)
  end if
  !call eos%lookup (shape(e), ee=e, d=d, pg=pg, gamma=gamma)
  !
  if (self%do_pebbles) then
    self%vgas(:,:,:,1:3,self%it) = q(:,:,:,2:4)
    self%pressure(:,:,:,  self%it) = pg
  end if
  c(:,:,:)=sqrt(gamma*q(:,:,:,5)/q(:,:,:,1))
  if (courant_type==1) then
    c = c + max(abs(q(:,:,:,2)), abs(q(:,:,:,3)), abs(q(:,:,:,4)))
  else
    c = c + (abs(q(:,:,:,2)) + abs(q(:,:,:,3)) + abs(q(:,:,:,4)))/3.0
  end if
  l = self%mesh%lo
  u = self%mesh%uo
  self%u_max = maxval(c(l(1):u(1),l(2):u(2),l(3):u(3)))
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end
  end if
  end associate
end subroutine ctoprim

!===============================================================================
!> Slope limiters
!===============================================================================
subroutine uslope(self)
  class(hd_t):: self
  integer::ngrid
  integer::i, j, k, n, l(3), u(3), islope_type
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, c1, c2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,idff,xslope_type
  integer,save:: itimer=0
  associate (q=>self%q, dq=>self%dq)
  !-----------------------------------------------------------------------------
  if (detailed_timer) then
    call trace%begin ('hd_mod::uslope', itimer=itimer)
  else
    call trace%begin ('hd_mod::uslope')
  end if
  !
  islope_type = int(slope_type)
  l = self%mesh%lo
  u = self%mesh%uo
  if(islope_type==-1)then
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
  else if(islope_type==0)then
    dq=0.0_dp
    if (detailed_timer) then
      call trace%end (itimer)
    else
      call trace%end
    end if
    return
  else if (islope_type==1) then  ! minmod
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
  else if (islope_type==2) then ! moncen
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
  else if (islope_type==3) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 2_dp,kind=dp)
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
            idff  = xslope_type / sqrt(dfx**2+dfy**2+dfz**2+tiny(1.0_dp)) ! xslope_type = 1 or 2 dep. on slope_type
            !
            ! vmin and vmax are the lower and upper limit of one-cell differences
            ! so a constant slope with step dq will give vmin=-dq and vmax=+dq,
            ! while dfx (dff for any direction) will also be = dq, unless we keep
            ! the 'half' in front of the sqrt.  With that kept, the dlim factor
            ! will remain unity until the slopes differ with a factor of 2. Then,
            ! if any of the 26 slopes is steeper than twice the gradient slope, 
            ! the gradient will be diminished.
            !
            ! idff > 0, so 1/dff = idff < 0.001 * huge is equivalent to dff > 0
            slop = merge(min(one,min(abs(vmin),abs(vmax))*idff), one, idff<0.001_dp * huge(1.0_dp))
            dlim = slop
            dq(i,j,k,n,1) = dlim*dfx
            dq(i,j,k,n,2) = dlim*dfy
            dq(i,j,k,n,3) = dlim*dfz
          end do
        end do
      end do
    end do
  else if (islope_type==4) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 3_dp,kind=dp)
    do k = l(3),u(3)
      do j = l(2),u(2)
        do i = l(1),u(1)
          block
          real(dp):: adfx(self%nv), adfy(self%nv), adfz(self%nv)
          dlim = 1.0
          do n = 1, self%nv
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
            idff  = xslope_type / sqrt(dfx**2+dfy**2+dfz**2+tiny(1.0_dp)) ! xslope_type = 1 or 2 dep. on slope_type
            !-------------------------------------------------------------------
            ! vmin and vmax are the lower and upper limit of one-cell differences
            ! so a constant slope with step dq will give vmin=-dq and vmax=+dq,
            ! while dfx (dff for any direction) will also be = dq, unless we keep
            ! the 'half' in front of the sqrt.  With that kept, the dlim factor
            ! will remain unity until the slopes differ with a factor of 2. Then,
            ! if any of the 26 slopes is steeper than twice the gradient slope, 
            ! the gradient will be diminished.
            !
            ! idff > 0, so 1/dff = idff < 0.001 * huge is equivalent to dff > 0
            !-------------------------------------------------------------------
            slop = merge(min(one,min(abs(vmin),abs(vmax))*idff), one, idff<0.001_dp * huge(1.0_dp))
            !-------------------------------------------------------------------
            ! Combine the dlim for different variables into one (experimental)
            !-------------------------------------------------------------------
            !dlim = slop
            !if (n>=2.and.n<=5) dlim = min(dlim,slop)
            dlim = min(dlim,slop)
            ! save fluxes for each variable
            adfx(n) = dfx
            adfy(n) = dfy
            adfz(n) = dfz
          end do
          do n = 1, self%nv
            ! use saved fluxes, with common dlim
            dq(i,j,k,n,1) = dlim*adfx(n)
            dq(i,j,k,n,2) = dlim*adfy(n)
            dq(i,j,k,n,3) = dlim*adfz(n)
          end do
          end block
        end do
      end do
    end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end
  end if
  end associate
end subroutine uslope

!===============================================================================
!> Predictor step, using primitive variables
!===============================================================================
subroutine trace3d(self)
  class(hd_t):: self
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
  if (detailed_timer) then
    call trace%begin ('hd_mod::trace3d', itimer=itimer)
  else
    call trace%begin ('hd_mod::trace3d')
  end if
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
        associate (ff => self%force_per_unit_mass)
        su0 = -u*dux-v*duy-w*duz - dpx/r + ff(i,j,k,1)*self%mesh(1)%d
        sv0 = -u*dvx-v*dvy-w*dvz - dpy/r + ff(i,j,k,2)*self%mesh(2)%d
        sw0 = -u*dwx-v*dwy-w*dwz - dpz/r + ff(i,j,k,3)*self%mesh(3)%d
        end associate

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
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end
  end if
  end associate
end subroutine trace3d

!===============================================================================
!> Compute fluxes
!===============================================================================
subroutine cmpflxm(self)
  class(hd_t):: self
  ! local variables
  integer:: i, j, k, lb(3), ub(3)
  real(dp), dimension(self%mesh(1)%gn,self%nv):: qleft, qright, fgdnv
  real(dp):: dtds(3)
  logical:: error
  integer,save:: itimer=0
  associate (qp=>self%qp, qm=>self%qm, fl=>self%fl)
  !-----------------------------------------------------------------------------
  if (detailed_timer) then
    call trace%begin ('hd_mod::cmpflxm', itimer=itimer)
  else
    call trace%begin ('hd_mod::cmpflxm')
  end if
  error = .false.
  dtds = self%dtime/self%mesh%d
  lb = self%mesh%li
  ub = self%mesh%uo
  qleft=0.0; qright=0.0
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        qleft (i,1) = qm(i,j,k,1,1); qright(i,1) = qp(i,j,k,1,1)
        qleft (i,2) = qm(i,j,k,2,1); qright(i,2) = qp(i,j,k,2,1)
        qleft (i,3) = qm(i,j,k,5,1); qright(i,3) = qp(i,j,k,5,1)
        qleft (i,4) = qm(i,j,k,3,1); qright(i,4) = qp(i,j,k,3,1)
        qleft (i,5) = qm(i,j,k,4,1); qright(i,5) = qp(i,j,k,4,1)
      end do
      if (detailed_timer) call trace_end (itimer)
      call solver (self,self%mesh(1)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (detailed_timer) call trace_begin ('hd_mod::cmpflxm', itimer=itimer)
      if (error) call error_handler('x', j, k, lb(1), ub(1))
      do i = lb(1),ub(1)
        fl(i,j,k,self%idx%d ,1) = fgdnv(i,1)*dtds(1)
        fl(i,j,k,self%idx%px,1) = fgdnv(i,2)*dtds(1)
        fl(i,j,k,self%idx%e ,1) = fgdnv(i,3)*dtds(1)
        fl(i,j,k,self%idx%py,1) = fgdnv(i,4)*dtds(1)
        fl(i,j,k,self%idx%pz,1) = fgdnv(i,5)*dtds(1)
      end do
    end do
  end do
  if (any(self%mesh%gn/=self%mesh(1)%gn)) then
    call cmpflxm_y
    call cmpflxm_z
  else
    do k = lb(3),ub(3)
      do i = lb(1),ub(1)
        do j = lb(2),ub(2)
          qleft (j,1) = qm(i,j,k,1,2); qright(j,1) = qp(i,j,k,1,2)
          qleft (j,2) = qm(i,j,k,3,2); qright(j,2) = qp(i,j,k,3,2)
          qleft (j,3) = qm(i,j,k,5,2); qright(j,3) = qp(i,j,k,5,2)
          qleft (j,4) = qm(i,j,k,2,2); qright(j,4) = qp(i,j,k,2,2)
          qleft (j,5) = qm(i,j,k,4,2); qright(j,5) = qp(i,j,k,4,2)
        end do
        if (detailed_timer) call trace_end (itimer)
        call solver (self,self%mesh(2)%gn,qleft,qright,fgdnv,error,detailed_timer)
        if (detailed_timer) call trace_begin ('hd_mod::cmpflxm', itimer=itimer)
        if (error) call error_handler('y', i, k, lb(2), ub(2))
        do j = lb(2),ub(2)
          fl(i,j,k,self%idx%d ,2) = fgdnv(j,1)*dtds(2)
          fl(i,j,k,self%idx%py,2) = fgdnv(j,2)*dtds(2)
          fl(i,j,k,self%idx%e ,2) = fgdnv(j,3)*dtds(2)
          fl(i,j,k,self%idx%px,2) = fgdnv(j,4)*dtds(2)
          fl(i,j,k,self%idx%pz,2) = fgdnv(j,5)*dtds(2)
        end do
      end do
    end do
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        do k = lb(3),ub(3)
          qleft (k,1) = qm(i,j,k,1,3); qright(k,1) = qp(i,j,k,1,3)
          qleft (k,2) = qm(i,j,k,4,3); qright(k,2) = qp(i,j,k,4,3)
          qleft (k,3) = qm(i,j,k,5,3); qright(k,3) = qp(i,j,k,5,3)
          qleft (k,4) = qm(i,j,k,2,3); qright(k,4) = qp(i,j,k,2,3)
          qleft (k,5) = qm(i,j,k,3,3); qright(k,5) = qp(i,j,k,3,3)
        end do
        if (detailed_timer) call trace_end (itimer)
        call solver (self,self%mesh(3)%gn,qleft,qright,fgdnv,error,detailed_timer)
        if (detailed_timer) call trace_begin ('hd_mod::cmpflxm', itimer=itimer)
        if (error) call error_handler('z', i, j, lb(3), ub(3))
        do k = lb(3),ub(3)
          fl(i,j,k,self%idx%d ,3) = fgdnv(k,1)*dtds(3)
          fl(i,j,k,self%idx%pz,3) = fgdnv(k,2)*dtds(3)
          fl(i,j,k,self%idx%e ,3) = fgdnv(k,3)*dtds(3)
          fl(i,j,k,self%idx%px,3) = fgdnv(k,4)*dtds(3)
          fl(i,j,k,self%idx%py,3) = fgdnv(k,5)*dtds(3)
        end do
      end do
    end do
  end if
  end associate
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end
  end if
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
  qleft=0.0; qright=0.0
  do k = lb(3),ub(3)
    do i = lb(1),ub(1)
      qleft (:,1) = qm(i,:,k,1,2); qright(:,1) = qp(i,:,k,1,2)
      qleft (:,2) = qm(i,:,k,3,2); qright(:,2) = qp(i,:,k,3,2)
      qleft (:,3) = qm(i,:,k,5,2); qright(:,3) = qp(i,:,k,5,2)
      qleft (:,4) = qm(i,:,k,2,2); qright(:,4) = qp(i,:,k,2,2)
      qleft (:,5) = qm(i,:,k,4,2); qright(:,5) = qp(i,:,k,4,2)
      if (detailed_timer) call trace_end (itimer)
      call solver (self,self%mesh(2)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (error) call error_handler('y', i, k, lb(2), ub(2))
      if (detailed_timer) call trace_begin ('hd_mod::cmpflxm', itimer=itimer)
      fl(i,:,k,self%idx%d ,2) = fgdnv(:,1)*dtds(2)
      fl(i,:,k,self%idx%py,2) = fgdnv(:,2)*dtds(2)
      fl(i,:,k,self%idx%e ,2) = fgdnv(:,3)*dtds(2)
      fl(i,:,k,self%idx%px,2) = fgdnv(:,4)*dtds(2)
      fl(i,:,k,self%idx%pz,2) = fgdnv(:,5)*dtds(2)
    end do
  end do
  end associate
  end subroutine cmpflxm_y
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine cmpflxm_z
  ! local variables
  real(dp), dimension(self%mesh(3)%gn,self%nv):: qleft, qright, fgdnv
  associate (qp=>self%qp, qm=>self%qm, fl=>self%fl)
  !-----------------------------------------------------------------------------
  dtds = self%dtime/self%mesh%d       ! needed here
  lb = self%mesh%li
  ub = self%mesh%uo
  qleft=0.0; qright=0.0
  do j = lb(2),ub(2)
    do i = lb(1),ub(1)
      qleft (:,1) = qm(i,j,:,1,3); qright(:,1) = qp(i,j,:,1,3)
      qleft (:,2) = qm(i,j,:,4,3); qright(:,2) = qp(i,j,:,4,3)
      qleft (:,3) = qm(i,j,:,5,3); qright(:,3) = qp(i,j,:,5,3)
      qleft (:,4) = qm(i,j,:,2,3); qright(:,4) = qp(i,j,:,2,3)
      qleft (:,5) = qm(i,j,:,3,3); qright(:,5) = qp(i,j,:,3,3)
      if (detailed_timer) call trace_end (itimer)
      call solver (self,self%mesh(3)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (error) call error_handler('z', i, j, lb(3), ub(3))
      if (detailed_timer) call trace_begin ('hd_mod::cmpflxm', itimer=itimer)
      fl(i,j,:,self%idx%d ,3) = fgdnv(:,1)*dtds(3)
      fl(i,j,:,self%idx%pz,3) = fgdnv(:,2)*dtds(3)
      fl(i,j,:,self%idx%e ,3) = fgdnv(:,3)*dtds(3)
      fl(i,j,:,self%idx%px,3) = fgdnv(:,4)*dtds(3)
      fl(i,j,:,self%idx%py,3) = fgdnv(:,5)*dtds(3)
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

!===============================================================================
!> Selection of 1-D solvers
!===============================================================================
subroutine solver (self,ngrid,qleft,qright,fgdnv,error,detailed_timer)
  class(hd_t):: self
  integer:: ngrid
  real(dp),dimension(:,:):: qleft,qright,fgdnv
  logical:: error, detailed_timer
  !.............................................................................
  if(self%riemann.eq.'acoustic')then
    call rieman%acoustic(qleft,qright,fgdnv,ngrid,self%idx%e)
  else if (self%riemann.eq.'exact')then
    call rieman%approx  (qleft,qright,fgdnv,ngrid,self%idx%e)
  else if (self%riemann.eq.'llf')then
    call rieman%llf     (qleft,qright,fgdnv,ngrid,self%idx%e)
  else if (self%riemann.eq.'hllc')then
    call rieman%hllc    (qleft,qright,fgdnv,ngrid,self%idx%e,isothermal,error,detailed_timer)
  else if (self%riemann.eq.'hllc_a')then
    call rieman%hllc_a  (qleft,qright,fgdnv,ngrid,self%idx%e)
  else if (self%riemann.eq.'hllc_v')then
    call rieman%hllc_v  (qleft,qright,fgdnv,ngrid,self%idx%e)
  else if (self%riemann.eq.'hll')then
    call rieman%hll     (qleft,qright,fgdnv,ngrid,self%idx%e)
  else
    write(*,*)'unknown Riemann solver'
    stop
  end if
end subroutine solver

!===============================================================================
!===============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(hd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  associate (d  => self%mem(:,:,:,self%idx%d ,self%it,1), &
             px => self%mem(:,:,:,self%idx%px,self%it,1), &
             py => self%mem(:,:,:,self%idx%py,self%it,1), &
             pz => self%mem(:,:,:,self%idx%pz,self%it,1), &
             e  => self%mem(:,:,:,self%idx%e ,self%it,1))
  if (isothermal) then
    pg = d*self%csound**2
  else
    pg = (gamma-1d0)*(e-0.5*((px**2+py**2+pz**2)/d))
  end if
  end associate
END FUNCTION gas_pressure

!===============================================================================
!===============================================================================
SUBROUTINE gravitational_force (self)
  class(hd_t):: self
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
