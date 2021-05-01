!===============================================================================
!> RAMSES Godunov solvers
!>
!> godunov
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
!>
!> Memory access and locks:  The sequence of access is
!>
!> update                            uin    unew    q
!>   unew = uin                      |      |
!>   godunov
!>     unsplit
!>       ctoprim                            |
!>       uslope                                     |
!>       trace3d                                    |
!>       cmpflxm                                    |
!>     unew .. flux update                  |
!>   unew ... forcing and heating           |
!>
!> By copying uin to unew early, that slot can be released, and only unew neeeds
!> to remain locked until the end of the update step
!===============================================================================
MODULE hd_mod
  !.............................................................................
  USE io_mod
  USE io_unit_mod
  USE task_mod
  USE extras_mod
  USE link_mod
  USE parameters_mod
  USE riemann_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE eos_mod
  USE index_mod
  USE patch_mod
  USE validate_mod
  USE scratch_mod
  implicit none
  PRIVATE
  type, public, extends(extras_t):: hd_t
    logical:: mhd=.false.
    character(len=10):: riemann='hllc'
    real, pointer, contiguous, dimension(:,:,:,:):: uin       ! conserved
    real, pointer, contiguous, dimension(:,:,:,:):: unew      ! new conserved
    real, pointer, contiguous, dimension(:,:,:,:):: q         ! primitive
    real, pointer, contiguous, dimension(:,:,:,:,:):: qp, qm  ! face centered
    real, pointer, contiguous, dimension(:,:,:,:,:):: dq      ! slopes
    real, pointer, contiguous, dimension(:,:,:,:,:):: fl      ! 3D fluxes
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: dealloc
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure:: nbor_relations
  end type
  integer, parameter:: ndim=3
  real, save:: slope_type=2.0
  real, save:: courant_factor=0.25
  integer, save:: verbose=0
  integer, save:: courant_type=1
  logical, save:: debug=.false.
  logical, save:: first_time=.true.
  logical, save:: detailed_timer=.false.
  logical, save:: unsigned=.true.
  character(len=16):: eqn='mhd'
CONTAINS

!===============================================================================
!> Pre-initialize the hydro solver, by indicating how many slots to use.
!> Other extensions, such as selfgravity, are then free to add to that
!===============================================================================
SUBROUTINE pre_init (self)
  class(hd_t), target:: self
  character(len=10):: riemann='hllc'
  real, save:: csound=1.
  integer:: i, iv, iostat
  namelist /ramses_params/ riemann, slope_type, courant_factor, smallr, &
    smallc, isothermal, csound, detailed_timer, unsigned, courant_type, debug, &
    verbose
  type(index_t):: idx
  !-----------------------------------------------------------------------------
  call trace%begin('hd_t%pre_init')
  !-----------------------------------------------------------------------------
  ! Make sure the extras modules have the correct starting point
  !-----------------------------------------------------------------------------
  self%nw = 1
  self%nv = 5
  self%ng = 2
  self%mhd = .false.
  self%kind = 'ramses_hd_patch'
  call self%idx%init (energy=5, mhd=self%mhd)
  call self%extras_t%pre_init
  !-----------------------------------------------------------------------------
  ! Make sure the extras modules do not mess this up
  !-----------------------------------------------------------------------------
  self%mhd = .false.
  self%kind = 'ramses_hd_patch'
  !-----------------------------------------------------------------------------
  ! Ramses specific parameters
  !-----------------------------------------------------------------------------
  call io%lock%set ('hd_t%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, ramses_params, iostat=iostat)
    call io%namelist_warning ('ramses_params', iostat)
    if (eos%gamma == 1d0) isothermal = .true.
    if (isothermal) eos%gamma = 1d0
    write (io%namelist, ramses_params)
    first_time = .false.
  end if
  call io%lock%unset ('hd_t%init')
  self%riemann = riemann
  self%courant = courant_factor
  self%csound = csound
  self%staggered = .false.
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the HD solver
!===============================================================================
SUBROUTINE init (self)
  class(hd_t):: self
  integer:: i
  character(len=120):: id = &
   '$Id$'
  !.............................................................................
  call trace%begin ('hd_t%init')
  call trace%print_id (id)
  call self%extras_t%init
  !-----------------------------------------------------------------------------
  ! Set the staggering, unsigned and pervolume flags before calling the rest
  ! of the cascade if %init(), where %gpatch_t%init may need them, in IC calls
  !-----------------------------------------------------------------------------
  do i=1,3
    self%mesh(i)%h = 0d0
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  !-----------------------------------------------------------------------------
  ! Initialize everything else
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Deallocate permanently allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(hd_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%dealloc')
  call self%extras_t%dealloc
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
! Preparation steps, of which calling extras_t%pre_update is mandatory.
! First, primitive variables, which include u_max, are computed and made
! available to the pre_update step of extras_t procedures, which may need to,
! or choose to, make use of the information to call their courant_conditions.
! The pre_update step of extras_t procedures may also need to know the final
! timestep, so it isn't called until after the final courant_condition call,
! via self%gpatch_t, has actually set self%dtime.
!
! Note that extras_t components that need to call their own courant_conditions
! only as part of their pre_update are free to do do, by not placing their
! calls inside the extras_t%courant_condition.
!===============================================================================
SUBROUTINE pre_update (self)
  class(hd_t), target:: self
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%pre_update')
  call self%extras_t%pre_update ()
  !-----------------------------------------------------------------------------
  ! Primitive variables may be needed in extras courant_conditions
  !-----------------------------------------------------------------------------
  if (io%smallr  /= 0.0) smallr = io%smallr
  if (io%courant /= 0.0) self%courant = io%courant
  call ctoprim(self)                  ! Primitive variables, sound speeds
  call self%extras_t%courant_condition (detailed_timer)
  !-----------------------------------------------------------------------------
  ! Pre-update extras components
  !-----------------------------------------------------------------------------
  call self%extras_t%pre_update
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(hd_t):: self
  real, dimension(:,:,:), pointer, contiguous:: d, e, px, py, pz
  real, dimension(:,:,:), allocatable:: ek
  real, dimension(:,:,:,:), pointer:: p
  real:: u0, u1
  type(lock_t), pointer:: lock
  integer:: n(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call trace_begin ('hd_t%update', itimer=itimer)
  call self%output
  !-----------------------------------------------------------------------------
  ! Memory allocations
  !-----------------------------------------------------------------------------
  n = self%mesh%gn
  allocate (self%qp(n(1),self%nv,ndim,n(2),n(3)))
  allocate (self%qm(n(1),self%nv,ndim,n(2),n(3)))
  allocate (self%dq(n(1),self%nv,ndim,n(2),n(3)))
  allocate (self%fl(n(1),n(2),n(3),self%nv,ndim))
  !-----------------------------------------------------------------------------
  ! Godunov solver
  !-----------------------------------------------------------------------------
  call godunov (self)
  !-----------------------------------------------------------------------------
  ! Save the current value of the kinetic energy
  !-----------------------------------------------------------------------------
  allocate (ek(n(1),n(2),n(3)))
  d  => self%mem(:,:,:,self%idx%d ,self%new,1)
  e  => self%mem(:,:,:,self%idx%e ,self%new,1)
  px => self%mem(:,:,:,self%idx%px,self%new,1)
  py => self%mem(:,:,:,self%idx%py,self%new,1)
  pz => self%mem(:,:,:,self%idx%pz,self%new,1)
  ek = 0.5*(px**2+py**2+pz**2)/d
  !-----------------------------------------------------------------------------
  ! Add the external force per unit mass as a source term
  !-----------------------------------------------------------------------------
                                                 lock => self%mem_lock(self%new)
                                                 call lock%set ('hd_t%update')
  if (associated(scratch%force_per_unit_mass)) then
    px = px + self%dtime*scratch%force_per_unit_mass(:,:,:,1)*d
    py = py + self%dtime*scratch%force_per_unit_mass(:,:,:,2)*d
    pz = pz + self%dtime*scratch%force_per_unit_mass(:,:,:,3)*d
  end if
  if (associated(scratch%force_per_unit_volume)) then
    px = px + self%dtime*scratch%force_per_unit_volume(:,:,:,1)
    py = py + self%dtime*scratch%force_per_unit_volume(:,:,:,2)
    pz = pz + self%dtime*scratch%force_per_unit_volume(:,:,:,3)
  end if
  !-----------------------------------------------------------------------------
  ! Adjust the total energy
  !-----------------------------------------------------------------------------
  e = e + (0.5*(px**2+py**2+pz**2)/d - ek)
  deallocate (ek)
  !call validate%check (self%link, p, 'after force')
  !-----------------------------------------------------------------------------
  ! Add the external heating per unit mass as a source term
  !-----------------------------------------------------------------------------
  if (associated(scratch%heating_per_unit_mass)) then
    e = e + self%dtime*scratch%heating_per_unit_mass*d
  end if
  if (associated(scratch%heating_per_unit_volume)) then
    e = e + self%dtime*scratch%heating_per_unit_volume
  end if
  !-----------------------------------------------------------------------------
  ! Reset the energy variable in isothermal cases
  !-----------------------------------------------------------------------------
  if (isothermal) &
    self%mem(:,:,:,self%idx%e,self%new,1) = &
    self%mem(:,:,:,self%idx%d,self%new,1)
                                                 call lock%unset ('hd_t%update')
  !.............................................................................
  if (verbose>0) then
    block
    integer:: i
    do i=1,self%nv
      d => self%mem(:,:,:,i,self%new,1)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') 'unew', self%id, i, &
      self%fminval(d), self%fminloc(d)+1, self%fmaxval(d), self%fmaxloc(d)+1
    end do
    do i=1,3
      d => scratch%force_per_unit_mass(:,:,:,i)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') 'fpm', self%id, i, &
      self%fminval(d), self%fminloc(d)+1, self%fmaxval(d), self%fmaxloc(d)+1
    end do
    do i=1,3
      d => self%fl(:,:,:,5,i)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') 'flux', self%id, i, &
        self%fminval(d), self%fminloc(d), self%fmaxval(d), self%fmaxloc(d)
    end do
    end block
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate
  !-----------------------------------------------------------------------------
  call self%counter_update
  deallocate (self%qp, self%qm, self%dq, self%fl, self%q)
  call trace%end (itimer)
END SUBROUTINE update

!=======================================================================
!> Measure time in extras_t%post_update
!=======================================================================
SUBROUTINE post_update (self)
  class(hd_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%post_update', itimer=itimer)
  call self%extras_t%post_update ()
  call trace%end (itimer)
END SUBROUTINE post_update

!=======================================================================
!> Get velocities from momentum
!=======================================================================
SUBROUTINE p2U(self, U, it)
  class(hd_t) :: self
  real, contiguous, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  integer:: i
  real, contiguous, dimension(:,:,:)  , pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,             it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
  do i=1,size(p,4)
  U(:,:,:,i) = p(:,:,:,i)/d
  end do
END SUBROUTINE p2U

!=======================================================================
!> Put velocities to momentum
!=======================================================================
SUBROUTINE U2p(self, U, it)
  class(hd_t) :: self
  real, contiguous, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  integer:: i
  real, contiguous, dimension(:,:,:)  , pointer:: d
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d             ,it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
  do i=1,size(p,4)
  p(:,:,:,i)= U(:,:,:,i)*d
  end do
END SUBROUTINE U2p

!=======================================================================
!> Get thermal energy from total/specific energy
!=======================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(hd_t) :: self
  real, contiguous, dimension(:,:,:), pointer:: E_th
  integer:: it
  !.....................................................................
  integer :: i
  real, contiguous, dimension(:,:,:), pointer:: k, d, e
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d             ,it,1)
  e => self%mem(:,:,:,self%idx%e             ,it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  E_th = (e-k)/d
  deallocate(k)
END SUBROUTINE e2E_th

!=======================================================================
!> Get etropy per unit mass from total energy
!=======================================================================
SUBROUTINE e2S(self, s, it)
  class(hd_t) :: self
  real, contiguous, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
  integer :: i
  real, contiguous, dimension(:,:,:), pointer:: k, d, e
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  if (isothermal) then
    s = 0.0
    return
  end if
  d => self%mem(:,:,:,self%idx%d             ,it,1)
  e => self%mem(:,:,:,self%idx%e             ,it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  s = log((e-k)*(eos%gamma-1d0)/d**eos%gamma)/(eos%gamma-1d0)
  deallocate(k)
END SUBROUTINE e2S

!=======================================================================
!> Put etropy per unit mass to total energy
!=======================================================================
SUBROUTINE S2e (self, s, it)
  class(hd_t) :: self
  real, contiguous, dimension(:,:,:), pointer:: s
  integer:: it
  !.....................................................................
  integer :: i
  real, contiguous, dimension(:,:,:), pointer:: k, d, e
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d             ,it,1)
  e => self%mem(:,:,:,self%idx%e             ,it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
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
  e = d**eos%gamma*exp(s*(eos%gamma-1d0)) + k
  deallocate(k)
END SUBROUTINE S2e

!=======================================================================
!> Put thermal energy to total/specific energy
!=======================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(hd_t) :: self
  integer:: it
  real, contiguous, dimension(:,:,:), pointer:: E_th
  !.....................................................................
  integer :: i
  real, contiguous, dimension(:,:,:), pointer:: k, d, e
  real, dimension(:,:,:,:), pointer:: p
  !---------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d             ,it,1)
  e => self%mem(:,:,:,self%idx%e             ,it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,it,1)
  allocate(k(size(d,1),size(d,2),size(d,3)))
  k = 0.0
  do i=1,size(p,4)
    k = k+p(:,:,:,i)**2
  end do
  k = 0.5*k/d
  e = E_th*d+k
  deallocate(k)
END SUBROUTINE E_th2e

!===============================================================================
SUBROUTINE godunov (self)
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
  real, contiguous, dimension(:,:,:,:)  , pointer:: unew
  real, contiguous, dimension(:,:,:,:,:), pointer:: fl
  !-----------------------------------------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------------------------------------
  unew => self%unew
  fl   => self%fl
  if (.not.detailed_timer) call trace%begin ('hd_t%godunov')
  li = self%mesh%li
  ui = self%mesh%ui
  call check_small (self, unew)
  call unsplit (self)
  !-----------------------------------------------------------------------------
  ! Conservative update, looping over internal cells only
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('hd_t%godunov', itimer=itimer)
  do ivar=1,self%nv
    if (ivar /= self%idx%phi .and. ivar /= self%idx%dphi) then
      do idim=1,ndim
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
        !.......................................................................
        if (verbose>0) then
          block
          real, pointer:: d(:,:,:)
          d => self%fl(:,:,:,ivar,idim)
          write (stdlog,'(a,i6,2i4,2(1p,e14.4,3i4))') &
          'fl', self%id, ivar, idim, &
          self%fminval(d), self%fminloc(d), self%fmaxval(d), self%fmaxloc(d)
          end block
        end if
        !.......................................................................
      end do
    end if
  end do
  call check_small (self, unew)
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end()
  end if
END SUBROUTINE godunov

!===============================================================================
!===============================================================================
SUBROUTINE check_small (self, uu)
  class(hd_t):: self
  real, dimension(:,:,:,:):: uu
  integer:: i,j,k,lo(3),uo(3),n_smalld,n_smallp
  real(dp):: c, d, u, v, w, p, ekin, eth
  integer, save:: nprint=10
  integer, save:: itimer=0
  integer:: it, jt, iw
  logical:: panic
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('hd_t%check_small', itimer=itimer)
  n_smalld = 0
  n_smallp = 0
  if (isothermal .or. eos%gamma==1d0) then
    self%ee_min = smallc**2
  else
    self%ee_min = smallc**2/(eos%gamma-1d0)
  end if
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
      uu(i,j,k,self%idx%e) = max(uu(i,j,k,self%idx%e), self%ee_min)
    else
      eth = uu(i,j,k,self%idx%e)-ekin
      if (eth < d*self%ee_min) then
        uu(i,j,k,self%idx%e) = d*self%ee_min + ekin
        n_smallp = n_smallp+1
      end if
    end if
  end do
  end do
  end do
  if (n_smalld+n_smallp > 0) then
    write (stdlog,'(a,i6,1p,e14.5,2i6)') &
      'small: id,nr,np', self%id, self%time, n_smalld, n_smallp
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
END SUBROUTINE check_small

!===============================================================================
!>  UNSPLIT     Unsplit second order Godunov integrator for
!>              polytropic gas dynamics using MUSCL-HANCOCK
!===============================================================================
SUBROUTINE unsplit(self)
  class(hd_t):: self
  real(dp):: dtds(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('hd_t%unsplit')
  call self%dump (1)                  ! dump, now that dtime is known
  call uslope(self)                   ! TVD slopes
  call trace3d(self)                  ! predictor step
  call cmpflxm(self)                  ! Solve for 1D flux in X direction
  call trace%end()
END SUBROUTINE unsplit

!===============================================================================
!> Conservative (uin) to primitive (q) variables, stored in the order
!> 1:d, 2:u, 3:v, 4:w, 5:p
!===============================================================================
SUBROUTINE ctoprim (self)
  class(hd_t):: self
  !.............................................................................
  real(dp):: smallp, dtds
  integer,save:: itimer=0
  integer:: l(3), u(3), n(3), i1, i2, i3
  real:: s
  real, contiguous, dimension(:,:,:)  , pointer:: d, e, pg, c
  real, contiguous, dimension(:,:,:,:), pointer:: unew, q
  !-----------------------------------------------------------------------------
  ! NOTE: Do NOT move this allocate to some other place; it needs to be here
  ! to allow the courant_condition to be evaluated in pre_update
  !-----------------------------------------------------------------------------
  n = self%mesh%gn
  allocate (self%q (n(1),n(2),n(3),self%nv))
  q    => self%q
  c    => self%u_phase
  unew => self%unew
  pg => q(:,:,:,5)
  if (detailed_timer) then
    call trace%begin ('hd_t%ctoprim', itimer=itimer)
  else
    call trace%begin ('hd_t%ctoprim')
  end if
  smallp = smallr*smallc**2/eos%gamma
  do i3=1,self%gn(3)
  do i2=1,self%gn(2)
  do i1=1,self%gn(1)
    q(i1,i2,i3,1) = max(unew(i1,i2,i3,self%idx%d),smallr)
    s = 1.0/q(i1,i2,i3,1)
    q(i1,i2,i3,2) = unew(i1,i2,i3,self%idx%px)*s
    q(i1,i2,i3,3) = unew(i1,i2,i3,self%idx%py)*s
    q(i1,i2,i3,4) = unew(i1,i2,i3,self%idx%pz)*s
    s = s*(q(i1,i2,i3,2)**2+q(i1,i2,i3,3)**2+q(i1,i2,i3,4)**2)
    pg(i1,i2,i3) = unew(i1,i2,i3,self%idx%e) - half*s
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Equation-of-state: s is the internal energy, and pg is the gas pressure
  !-----------------------------------------------------------------------------
  d => q(:,:,:,1)
  e => pg
  if (isothermal) then
    pg = d*self%csound**2
  else
    pg = e*(eos%gamma-1d0)
  end if
  !.............................................................................
  l = self%mesh%li
  u = self%mesh%ui
  if (any(pg(l(1):u(1),l(2):u(2),l(3):u(3)) < 0.0) .or. verbose > 0) then
    write (stdlog,'(a,i6,1p,2e11.2,3i4)') 'ctoprim: max(pg),min(pg),loc=', &
      self%id, self%fmaxval(pg), self%fminval(pg), self%fminloc(pg)
    pg = max(pg,smallp)
  end if
  !.............................................................................
  c(:,:,:)=sqrt(eos%gamma*q(:,:,:,5)/q(:,:,:,1))
  if (courant_type==1) then
    c = c + max(abs(q(:,:,:,2)), abs(q(:,:,:,3)), abs(q(:,:,:,4)))
  else
    c = c + (abs(q(:,:,:,2)) + abs(q(:,:,:,3)) + abs(q(:,:,:,4)))/3.0
  end if
  !-----------------------------------------------------------------------------
  ! With this solver we can take into account the phase speed also in the
  ! entire guard zone, since it is based on values valid in all points
  !-----------------------------------------------------------------------------
  self%u_max = maxval(c)
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end()
  end if
END SUBROUTINE ctoprim

!===============================================================================
!> Slope limiters
!===============================================================================
SUBROUTINE uslope(self)
  class(hd_t):: self
  integer::ngrid
  integer::i, j, k, n, l(3), u(3), islope_type, nv
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, c1, c2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,idff,xslope_type
  real, dimension(:,:,:,:)  , pointer, contiguous:: q
  real, dimension(:,:,:,:,:), pointer, contiguous:: dq
  integer,save:: itimer=0
  !-----------------------------------------------------------------------------
  q  => self%q
  dq => self%dq
  if (detailed_timer) then
    call trace%begin ('hd_t%uslope', itimer=itimer)
  else
    call trace%begin ('hd_t%uslope')
  end if
  !
  islope_type = int(slope_type)
  l = self%mesh%lo
  u = self%mesh%uo
  if(islope_type==-1)then
    do n = 1, self%nv
      if (n==self%idx%phi .or. n==self%idx%dphi) cycle
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd
          do i = l(1),u(1)
            dq(i,n,1,j,k) = half*(q(i+1,j,k,n) - q(i-1,j,k,n))
            dq(i,n,2,j,k) = half*(q(i,j+1,k,n) - q(i,j-1,k,n))
            dq(i,n,3,j,k) = half*(q(i,j,k+1,n) - q(i,j,k-1,n))
          end do
        end do
      end do
    end do
  else if(islope_type==0)then
    dq=0.0_dp
    if (detailed_timer) then
      call trace%end (itimer)
    else
      call trace%end()
    end if
    return
  else if (islope_type==1) then  ! minmod
    do n = 1, self%nv
      if (n==self%idx%phi .or. n==self%idx%dphi) cycle
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = q(i  ,j,k,n) - q(i-1,j,k,n)
            drgt = q(i+1,j,k,n) - q(i  ,j,k,n)
            dq(i,n,1,j,k) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
            ! slopes in second coordinate direction
            dlft = q(i,j  ,k,n) - q(i,j-1,k,n)
            drgt = q(i,j+1,k,n) - q(i,j  ,k,n)
            dq(i,n,2,j,k) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
            ! slopes in third coordinate direction
            dlft = q(i,j,k  ,n) - q(i,j,k-1,n)
            drgt = q(i,j,k+1,n) - q(i,j,k  ,n)
            dq(i,n,3,j,k) = merge(zero,merge(min(dlft,drgt),max(dlft,drgt),dlft>0),(dlft*drgt)<=zero)
          end do
        end do
      end do
    end do
  else if (islope_type==2) then ! moncen
    c1 = slope_type
    c2 = 0.5_dp/c1
    do n = 1, self%nv
      if (n==self%idx%phi .or. n==self%idx%dphi) cycle
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd private(dlft,drgt,dsgn,slop,dlim)
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = c1*(q(i  ,j,k,n) - q(i-1,j,k,n))
            drgt = c1*(q(i+1,j,k,n) - q(i  ,j,k,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,n,1,j,k) = dsgn*min(dlim,abs(dcen))
            ! slopes in second coordinate direction
            dlft = c1*(q(i,j  ,k,n) - q(i,j-1,k,n))
            drgt = c1*(q(i,j+1,k,n) - q(i,j  ,k,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one,dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,n,2,j,k) = dsgn*min(dlim,abs(dcen))
            ! slopes in third coordinate direction
            dlft = c1*(q(i,j,k  ,n) - q(i,j,k-1,n))
            drgt = c1*(q(i,j,k+1,n) - q(i,j,k  ,n))
            dcen = c2*(dlft+drgt)
            dsgn = sign(one,dcen)
            slop = min(abs(dlft),abs(drgt))
            dlim = merge(zero,slop,(dlft*drgt)<=zero)
            dq(i,n,3,j,k) = dsgn*min(dlim,abs(dcen))
          end do
        end do
      end do
    end do
  else if (islope_type==3) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 2_dp,kind=dp)
    do n = 1, self%nv
      if (n==self%idx%phi .or. n==self%idx%dphi) cycle
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ vector aligned
          !dir$ ivdep
          !$omp simd private(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml, &
          !$omp    & dfrrl,dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm, &
          !$omp    & dfrrm,dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr, &
          !$omp    & dfrrr,vmin,vmax,dfx,dfy,dfz,idff,slop,dlim)
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
            dlim = merge(min(one,min(abs(vmin),abs(vmax))*idff), one, idff<0.001_dp * huge(1.0_dp))
            !
            dq(i,n,1,j,k) = dlim*dfx
            dq(i,n,2,j,k) = dlim*dfy
            dq(i,n,3,j,k) = dlim*dfz
          end do
        end do
      end do
    end do
  else if (islope_type==4) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 3_dp,kind=dp)
    nv = self%nv
    if (self%idx%phi > 0) then
      nv = self%idx%phi-1
    else
      nv = self%nv
    end if
    do k = l(3),u(3)
      do j = l(2),u(2)
        do i = l(1),u(1)
          block
          real(dp):: adfx(self%nv), adfy(self%nv), adfz(self%nv)
          dlim = 1.0
          !dir$ ivdep
          !$omp simd private(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml, &
          !$omp    & dfrrl,dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm, &
          !$omp    & dfrrm,dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr, &
          !$omp    & dfrrr,vmin,vmax,dfx,dfy,dfz,idff,slop,dlim)
          do n = 1, nv
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
          !dir$ ivdep
          !$omp simd
          do n = 1, nv
            ! use saved fluxes, with common dlim
            dq(i,n,1,j,k) = dlim*adfx(n)
            dq(i,n,2,j,k) = dlim*adfy(n)
            dq(i,n,3,j,k) = dlim*adfz(n)
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
    call trace%end()
  end if
END SUBROUTINE uslope

!===============================================================================
!> Predictor step, using primitive variables
!===============================================================================
SUBROUTINE trace3d(self)
  class(hd_t):: self
  integer ::i, j, k, lb(3), ub(3)
  integer ::ir=1, iu=2, iv=3, iw=4, ip=5
  real(dp)::dtdx, dtdy, dtdz, smallp
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0
  real(dp)::ff(self%gn(1),3), qq(self%gn(1))
  real(dp), dimension(:,:,:)    , pointer, contiguous:: d
  real(dp), dimension(:,:,:,:)  , pointer, contiguous:: q
  real(dp), dimension(:,:,:,:,:), pointer, contiguous:: qp, qm, dq
  real(8):: dt2
  integer,save:: itimer=0
  !-----------------------------------------------------------------------------
  q  => self%q
  dq => self%dq
  qp => self%qp
  qm => self%qm
  if (detailed_timer) then
    call trace%begin ('hd_t%trace3d', itimer=itimer)
  else
    call trace%begin ('hd_t%trace3d')
  end if
  dtdx = self%dtime/self%mesh(1)%d
  dtdy = self%dtime/self%mesh(2)%d
  dtdz = self%dtime/self%mesh(3)%d
  smallp = smallr*smallc**2/eos%gamma
  !-----------------------------------------------------------------------------
  dt2 = 0.5d0*self%dtime
  lb = self%mesh%lo
  ub = self%mesh%uo
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      d => self%mem(:,:,:,self%idx%d,self%new,1)
      ff(:,:) = 0.0_dp
      if (associated(scratch%force_per_unit_mass)) then
        do i=lb(1),ub(1)
          ff(i,1) = ff(i,1) + scratch%force_per_unit_mass(i,j,k,1)*dt2
          ff(i,2) = ff(i,2) + scratch%force_per_unit_mass(i,j,k,2)*dt2
          ff(i,3) = ff(i,3) + scratch%force_per_unit_mass(i,j,k,3)*dt2
        end do
      end if
      if (associated(scratch%force_per_unit_volume)) then
        do i=lb(1),ub(1)
          ff(i,1) = ff(i,1) + scratch%force_per_unit_volume(i,j,k,1)*dt2/d(i,j,k)
          ff(i,2) = ff(i,2) + scratch%force_per_unit_volume(i,j,k,2)*dt2/d(i,j,k)
          ff(i,3) = ff(i,3) + scratch%force_per_unit_volume(i,j,k,3)*dt2/d(i,j,k)
        end do
      end if
      qq(:) = 0.0
      if (associated(scratch%heating_per_unit_volume)) then
        do i=lb(1),ub(1)
          qq(i) = qq(i) + scratch%heating_per_unit_volume(i,j,k)*dt2
        end do
      end if
      if (associated(scratch%heating_per_unit_mass)) then
        do i=lb(1),ub(1)
          qq(i) = qq(i) + scratch%heating_per_unit_mass(i,j,k)*d(i,j,k)*dt2
        end do
      end if
      !dir$ ivdep
      !$omp simd
      do i = lb(1),ub(1)
        ! Cell centered values
        r   =  q(i,j,k,ir)
        u   =  q(i,j,k,iu)
        v   =  q(i,j,k,iv)
        w   =  q(i,j,k,iw)
        p   =  q(i,j,k,ip)

        ! TVD slopes in all 3 directions
        drx = dq(i,ir,1,j,k)*half*dtdx
        dpx = dq(i,ip,1,j,k)*half*dtdx
        dux = dq(i,iu,1,j,k)*half*dtdx
        dvx = dq(i,iv,1,j,k)*half*dtdx
        dwx = dq(i,iw,1,j,k)*half*dtdx

        dry = dq(i,ir,2,j,k)*half*dtdy
        dpy = dq(i,ip,2,j,k)*half*dtdy
        duy = dq(i,iu,2,j,k)*half*dtdy
        dvy = dq(i,iv,2,j,k)*half*dtdy
        dwy = dq(i,iw,2,j,k)*half*dtdy

        drz = dq(i,ir,3,j,k)*half*dtdz
        dpz = dq(i,ip,3,j,k)*half*dtdz
        duz = dq(i,iu,3,j,k)*half*dtdz
        dvz = dq(i,iv,3,j,k)*half*dtdz
        dwz = dq(i,iw,3,j,k)*half*dtdz

        ! Source terms (including transverse derivatives)
        sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
        sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*eos%gamma*p + (eos%gamma-1d0)*qq(i)

        ! Uncomment to use unconstrained pressure gradients
        !dpx = half*(q(i+1,j,k,ip)-q(i-1,j,k,ip))*half*dtdx
        !dpy = half*(q(i,j+1,k,ip)-q(i,j-1,k,ip))*half*dtdy
        !dpz = half*(q(i,j,k+1,ip)-q(i,j,k-1,ip))*half*dtdz
        su0 = -u*dux-v*duy-w*duz - dpx/r + ff(i,1)
        sv0 = -u*dvx-v*dvy-w*dvz - dpy/r + ff(i,2)
        sw0 = -u*dwx-v*dwy-w*dwz - dpz/r + ff(i,3)

        ! Add predictor step in time
        r = r + sr0
        p = p + sp0
        u = u + su0
        v = v + sv0
        w = w + sw0

        ! Right state at left interface
        drx = dq(i,ir,1,j,k)*half
        dpx = dq(i,ip,1,j,k)*half
        dux = dq(i,iu,1,j,k)*half
        dvx = dq(i,iv,1,j,k)*half
        dwx = dq(i,iw,1,j,k)*half
        qp(i,ir,1,j,k) = r - drx
        qp(i,ip,1,j,k) = p - dpx
        qp(i,iu,1,j,k) = u - dux
        qp(i,iv,1,j,k) = v - dvx
        qp(i,iw,1,j,k) = w - dwx
        qp(i,ir,1,j,k) = max(smallr, qp(i,ir,1,j,k))
        qp(i,ip,1,j,k) = max(smallp, qp(i,ip,1,j,k))

        ! Left state at right interface
        qm(i+1,ir,1,j,k) = r + drx
        qm(i+1,ip,1,j,k) = p + dpx
        qm(i+1,iu,1,j,k) = u + dux
        qm(i+1,iv,1,j,k) = v + dvx
        qm(i+1,iw,1,j,k) = w + dwx
        qm(i+1,ir,1,j,k) = max(smallr, qm(i+1,ir,1,j,k))
        qm(i+1,ip,1,j,k) = max(smallp, qm(i+1,ip,1,j,k))

        ! Top state at bottom interface
        dry = dq(i,ir,2,j,k)*half
        dpy = dq(i,ip,2,j,k)*half
        duy = dq(i,iu,2,j,k)*half
        dvy = dq(i,iv,2,j,k)*half
        dwy = dq(i,iw,2,j,k)*half
        qp(i,ir,2,j,k) = r - dry
        qp(i,ip,2,j,k) = p - dpy
        qp(i,iu,2,j,k) = u - duy
        qp(i,iv,2,j,k) = v - dvy
        qp(i,iw,2,j,k) = w - dwy
        qp(i,ir,2,j,k) = max(smallr, qp(i,ir,2,j,k))
        qp(i,ip,2,j,k) = max(smallp, qp(i,ip,2,j,k))

        ! Bottom state at top interface
        qm(i,ir,2,j+1,k) = r + dry
        qm(i,ip,2,j+1,k) = p + dpy
        qm(i,iu,2,j+1,k) = u + duy
        qm(i,iv,2,j+1,k) = v + dvy
        qm(i,iw,2,j+1,k) = w + dwy
        qm(i,ir,2,j+1,k) = max(smallr, qm(i,ir,2,j+1,k))
        qm(i,ip,2,j+1,k) = max(smallp, qm(i,ip,2,j+1,k))

        ! Back state at front interface
        drz = dq(i,ir,3,j,k)*half
        dpz = dq(i,ip,3,j,k)*half
        duz = dq(i,iu,3,j,k)*half
        dvz = dq(i,iv,3,j,k)*half
        dwz = dq(i,iw,3,j,k)*half
        qp(i,ir,3,j,k) = r - drz
        qp(i,ip,3,j,k) = p - dpz
        qp(i,iu,3,j,k) = u - duz
        qp(i,iv,3,j,k) = v - dvz
        qp(i,iw,3,j,k) = w - dwz
        qp(i,ir,3,j,k) = max(smallr, qp(i,ir,3,j,k))
        qp(i,ip,3,j,k) = max(smallp, qp(i,ip,3,j,k))

        ! Front state at back interface
        qm(i,ir,3,j,k+1) = r + drz
        qm(i,ip,3,j,k+1) = p + dpz
        qm(i,iu,3,j,k+1) = u + duz
        qm(i,iv,3,j,k+1) = v + dvz
        qm(i,iw,3,j,k+1) = w + dwz
        qm(i,ir,3,j,k+1) = max(smallr, qm(i,ir,3,j,k+1))
        qm(i,ip,3,j,k+1) = max(smallp, qm(i,ip,3,j,k+1))
      end do
    end do
  end do
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end()
  end if
END SUBROUTINE trace3d

!===============================================================================
!> Compute fluxes
!===============================================================================
SUBROUTINE cmpflxm(self)
  class(hd_t):: self
  ! local variables
  integer:: i, j, k, lb(3), ub(3)
  real(dp), dimension(self%mesh(1)%gn,self%nv):: qleft, qright, fgdnv
  real(dp):: dtds(3)
  logical:: error
  integer,save:: itimer=0
  real(dp), dimension(:,:,:,:,:), pointer, contiguous:: qp, qm, fl
  !-----------------------------------------------------------------------------
  qp => self%qp
  qm => self%qm
  fl => self%fl
  if (detailed_timer) then
    call trace%begin ('hd_t%cmpflxm', itimer=itimer)
  else
    call trace%begin ('hd_t%cmpflxm')
  end if
  error = .false.
  dtds = self%dtime/self%mesh%d
  lb = self%mesh%li
  ub = self%mesh%uo
  qleft=0.0; qright=0.0
  do k = lb(3),ub(3)
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        qleft (i,1) = qm(i,1,1,j,k); qright(i,1) = qp(i,1,1,j,k)
        qleft (i,2) = qm(i,2,1,j,k); qright(i,2) = qp(i,2,1,j,k)
        qleft (i,3) = qm(i,5,1,j,k); qright(i,3) = qp(i,5,1,j,k)
        qleft (i,4) = qm(i,3,1,j,k); qright(i,4) = qp(i,3,1,j,k)
        qleft (i,5) = qm(i,4,1,j,k); qright(i,5) = qp(i,4,1,j,k)
      end do
      call solver (self,self%mesh(1)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (error) call error_handler('x', j, k, lb(1), ub(1))
      do i = lb(1),ub(1)
        fl(i,j,k,self%idx%d ,1) = fgdnv(i,1)*dtds(1)
        fl(i,j,k,self%idx%px,1) = fgdnv(i,2)*dtds(1)
        fl(i,j,k,self%idx%e ,1) = fgdnv(i,3)*dtds(1)
        fl(i,j,k,self%idx%py,1) = fgdnv(i,4)*dtds(1)
        fl(i,j,k,self%idx%pz,1) = fgdnv(i,5)*dtds(1)
      end do
      do i = lb(1),ub(1)
        qleft (i,1) = qm(i,1,2,j,k); qright(i,1) = qp(i,1,2,j,k)
        qleft (i,2) = qm(i,3,2,j,k); qright(i,2) = qp(i,3,2,j,k)
        qleft (i,3) = qm(i,5,2,j,k); qright(i,3) = qp(i,5,2,j,k)
        qleft (i,4) = qm(i,2,2,j,k); qright(i,4) = qp(i,2,2,j,k)
        qleft (i,5) = qm(i,4,2,j,k); qright(i,5) = qp(i,4,2,j,k)
      end do
      call solver (self,self%mesh(1)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (error) call error_handler('y', j, k, lb(1), ub(1))
      do i = lb(1),ub(1)
        fl(i,j,k,self%idx%d ,2) = fgdnv(i,1)*dtds(2)
        fl(i,j,k,self%idx%py,2) = fgdnv(i,2)*dtds(2)
        fl(i,j,k,self%idx%e ,2) = fgdnv(i,3)*dtds(2)
        fl(i,j,k,self%idx%px,2) = fgdnv(i,4)*dtds(2)
        fl(i,j,k,self%idx%pz,2) = fgdnv(i,5)*dtds(2)
      end do
      do i = lb(1),ub(1)
        qleft (i,1) = qm(i,1,3,j,k); qright(i,1) = qp(i,1,3,j,k)
        qleft (i,2) = qm(i,4,3,j,k); qright(i,2) = qp(i,4,3,j,k)
        qleft (i,3) = qm(i,5,3,j,k); qright(i,3) = qp(i,5,3,j,k)
        qleft (i,4) = qm(i,2,3,j,k); qright(i,4) = qp(i,2,3,j,k)
        qleft (i,5) = qm(i,3,3,j,k); qright(i,5) = qp(i,3,3,j,k)
      end do
      call solver (self,self%mesh(1)%gn,qleft,qright,fgdnv,error,detailed_timer)
      if (error) call error_handler('z', j, k, lb(1), ub(1))
      do i = lb(1),ub(1)
        fl(i,j,k,self%idx%d ,3) = fgdnv(i,1)*dtds(3)
        fl(i,j,k,self%idx%pz,3) = fgdnv(i,2)*dtds(3)
        fl(i,j,k,self%idx%e ,3) = fgdnv(i,3)*dtds(3)
        fl(i,j,k,self%idx%px,3) = fgdnv(i,4)*dtds(3)
        fl(i,j,k,self%idx%py,3) = fgdnv(i,5)*dtds(3)
      end do
    end do
  end do
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end()
  end if
contains
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine error_handler (label, i1, i2, lb, ub)
  class(*), pointer:: link
  character(len=*):: label
  integer:: i1, i2, lb, ub
  class(nbor_t), pointer:: nbor
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
  END SUBROUTINE
END SUBROUTINE cmpflxm

!===============================================================================
!> Selection of 1-D solvers
!===============================================================================
SUBROUTINE solver (self,ngrid,qleft,qright,fgdnv,error,detailed_timer)
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
END SUBROUTINE solver

!===============================================================================
!===============================================================================
SUBROUTINE gravitational_force (self)
  class(hd_t):: self
  real, contiguous, dimension(:,:,:), pointer:: phi, d, et, gradphi1, gradphi2, gradphi3
  real, dimension(:,:,:,:), pointer:: p
  integer:: i,j,k
  real, dimension(3):: wtp1, wtm1
  real:: etemp
  !.............................................................................
  associate (l => self%li, u => self%ui)
  wtp1 = 0.50 * 4.0 / 3.0 / self%ds
  wtm1 = 0.25 * 1.0 / 3.0 / self%ds
  phi => self%mem(:,:,:,self%idx%phi,self%new,1)
  d   => self%mem(:,:,:,self%idx%d  ,self%new,1)
  et  => self%mem(:,:,:,self%idx%e  ,self%new,1)
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

!===============================================================================
!> Set up nbor relations between MHD tasks
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(hd_t):: self
  class(task_t), pointer:: another, task
  logical:: needed, needs_me, download, overlaps
  ! ----------------------------------------------------------------------------
  integer, save:: itimer=0
  class(hd_t),  pointer:: nb_hd
  ! ----------------------------------------------------------------------------
  task => another
  call self%task_t%nbor_relations(task,needed,needs_me,download)
  select type (another)
  class is (hd_t)
    ! do nothing, all is done in the base task_t layer
  class default
    needed=.false.; needs_me=.false.; download=.false.
  end select

END SUBROUTINE nbor_relations

END MODULE
