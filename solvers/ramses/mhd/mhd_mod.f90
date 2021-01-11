!===============================================================================
!> RAMSES Godunov solvers, use of guard zones; specifically in HLLD
!>
!> godunov
!>   unsplit              1   2   3   4     ui  uo      ub
!>     ctoprim:uin        +---+---+---+......+---+---+---+    ! in, with guards
!>     ctoprim:bx       +---+---+---+......+---+---+---+      ! down-staggered
!>     ctoprim:q          +---+---+---+......+---+---+        ! up-staggered in b
!>     uslope:dq              +---+---+......+---+            ! both ends loose
!>     trace3d:qm           +---+---+......+---+              ! down-staggered
!>     trace3d:qp               +---+......+---+---+          ! up-shifted
!>     cmpflxm                  +---+......+---+              ! need qm & qp
!>       cmp_mag_flx            +---+......+---+              ! is down-staggered
!>     unew                       +---+......+                ! up-staggered
!===============================================================================
MODULE mhd_mod
  !.............................................................................
  !USE const
  USE io_mod
  USE io_unit_mod
  USE omp_lock_mod
  USE extras_mod
  USE link_mod
  USE riemann_mod
  use parameters_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE index_mod
  USE non_ideal_mod
  USE flux_mod
  USE omp_timer_mod
  USE timer_mod
  USE scratch_mod
  implicit none
  PRIVATE
  type, public, extends(extras_t):: mhd_t
    logical:: mhd=.true.
    integer:: idim=3
    character(len=10):: riemann='hlld'
    real, contiguous, pointer, dimension(:,:,:,:):: q         ! primitive
    real, contiguous, pointer, dimension(:,:,:,:):: bf, em    ! face centered B, emf
    real, contiguous, pointer, dimension(:,:,:,:,:):: qp, qm  ! face centered
    real, contiguous, pointer, dimension(:,:,:,:,:):: dq, dq2, dbf ! slopes
    real, contiguous, pointer, dimension(:,:,:,:,:):: qLB,qRB,qLT,qRT
    real, contiguous, pointer, dimension(:,:,:,:):: uin       ! conserved
    real, contiguous, pointer, dimension(:,:,:,:):: unew      ! new conserved
    real, contiguous, pointer, dimension(:,:,:):: Ex,Ey,Ez
    real, contiguous, pointer, dimension(:,:,:):: c
    !dir$ attributes align:64 :: q,bf,em,qp,qm,dq,dq2,dbf,qLB,qRB,qLT,qRT,Ex,Ey,Ez,c
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: update
  end type
  real, save:: slope_type=2.0
  real, save:: courant_factor=0.25
  logical, save:: first_time=.true.
  logical, save:: detailed_timer=.false.
  logical, save:: unsigned=.true.
  logical, save:: do_temperature=.false.
  integer, save:: verbose=0
  integer, save:: divb_cleaner=-1
  integer, save:: divb_iter=4
  integer, save:: b_slope_type=0
  integer, save:: u_slope_type=0
  integer, save:: courant_type=0
  character(len=16):: eqn='mhd'
CONTAINS

!===============================================================================
!> Pre-init; read parameters, etc
!===============================================================================
SUBROUTINE pre_init (self)
  class(mhd_t), target:: self
  character(len=10):: riemann='hlld'
  real, save:: max_dlogd=20.
  integer, save:: n_dump=0
  real:: csound=1.
  integer:: iostat
  namelist /ramses_params/ gamma, riemann, slope_type, b_slope_type, &
    u_slope_type, courant_factor, smallr, smallc, isothermal, csound, &
    detailed_timer, unsigned, courant_type, max_dlogd, n_dump, &
    divb_cleaner, divb_iter, do_temperature, verbose
  character(len=120):: id = &
   '$Id$'
!...............................................................................
  call trace%begin ('mhd_mod%pre_init')
  call trace%print_id (id)
  self%nw = 1
  self%ng = 3
  self%nv = 8
  call self%extras_t%pre_init
  self%kind = 'ramses_mhd_patch'
  call self%idx%init (ie=5, mhd=self%mhd)
  if (do_temperature) then
    call self%idx%request_index (self%idx%tt, self%nv)
  end if
  !-----------------------------------------------------------------------------
  ! Read solver-relevant parameters and store in the task instance
  !-----------------------------------------------------------------------------
  call io%lock%set ('mhd_mod%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, ramses_params, iostat=iostat)
    call io%namelist_warning ('ramses_params', iostat)
    if (gamma == 1d0) isothermal = .true.
    if (isothermal)   gamma = 1d0
    write (io%namelist, ramses_params)
    first_time = .false.
  end if
  call io%lock%unset ('mhd_mod%pre_init')
  rieman%max_dlogd = max_dlogd
  self%n_dump  = n_dump
  self%riemann = riemann
  self%courant = courant_factor
  self%gamma = gamma
  self%csound = csound
  self%staggered = .false.
  !-----------------------------------------------------------------------------
  call non_ideal%init  
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  integer:: i
!...............................................................................
  call trace%begin('mhd_mod%init')
  call self%extras_t%init
  if (minval(self%mesh%n) < 3) then
    call io%abort ('SOLVER=ramses/mhd needs at least 3 mesh points in any direction')
  end if
  do i=1,3
    self%mesh(i)%h = 0d0
    self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  call trace%end()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE pre_update (self)
  class(mhd_t), target:: self
  integer:: n(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%pre_update')
  !-----------------------------------------------------------------------------
  ! Primitive variables may be needed in extras courant_conditions
  !-----------------------------------------------------------------------------
  if (io%smallr  /= 0.0) smallr = io%smallr
  if (io%courant /= 0.0) self%courant = io%courant
  call ctoprim(self)                  ! Primitive variables, sound speeds
  call self%extras_t%courant_condition (detailed_timer)
  call self%gpatch_t%courant_condition (detailed_timer)
  !-----------------------------------------------------------------------------
  ! Pre-update extras components
  !-----------------------------------------------------------------------------
  call self%extras_t%pre_update
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  real, pointer, dimension(:,:,:,:):: p
  real, allocatable, dimension(:,:,:,:):: em
  real, pointer, contiguous, dimension(:,:,:,:,:):: fl
  real, allocatable, dimension(:,:,:,:)::fluxambdiff,emfambdiff,fluxohmdiss,emfohmdiss
  real, pointer, dimension(:,:,:):: d, e, Ux=>null(), Uy=>null(), Uz=>null()
  real(dp):: dtdx, dtdy, dtdz
  real:: u0, u1
  integer:: n(3), nv3, i, j, k, idir
  type(lock_t), pointer:: lock
  integer, parameter:: ndim=3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%update', itimer=itimer)
  n = self%mesh%ui-self%mesh%li+1
  if (io%smallr > 0) smallr = io%smallr
  if (io%courant > 0) self%courant = io%courant
  !----------------------------------------------------------------------------
  ! Checks whether to do output
  !----------------------------------------------------------------------------
  call self%output
  !-----------------------------------------------------------------------------
  ! Memory mapping and allocations.  All alloc and dealloc should be done here,
  ! even if it would be possible to do some of them locally, in the procedures
  ! being called.  Set all values to zero initially, in case some loop extends
  ! beyond defined values.
  !-----------------------------------------------------------------------------

  n = self%mesh%gn
  allocate (self%Ex (n(1),n(2),n(3)))                           ! self%Ex  = 0_dp
  allocate (self%Ey (n(1),n(2),n(3)))                           ! self%Ey  = 0_dp
  allocate (self%Ez (n(1),n(2),n(3)))                           ! self%Ez  = 0_dp
  allocate (self%qp (n(1),self%nv,ndim,n(2),n(3)))              ! self%qp  = 0_dp
  allocate (self%qm (n(1),self%nv,ndim,n(2),n(3)))              ! self%qm  = 0_dp
  allocate (self%dq (n(1),self%nv,ndim,n(2),n(3)))              ! self%dq  = 0_dp
  allocate (self%dq2(n(1),self%nv,ndim,n(2),n(3)))              ! self%dq2 = 0_dp
  allocate (self%dbf(n(1),3,2,n(2),n(3)))                       ! self%dbf = 0_dp
  allocate (self%qLB(n(1),self%nv,ndim,n(2),n(3)))              ! self%qLB = 0_dp
  allocate (self%qLT(n(1),self%nv,ndim,n(2),n(3)))              ! self%qLT = 0_dp
  allocate (self%qRB(n(1),self%nv,ndim,n(2),n(3)))              ! self%qRB = 0_dp
  allocate (self%qRT(n(1),self%nv,ndim,n(2),n(3)))              ! self%qRT = 0_dp
  d => self%mem(:,:,:,          self%idx%d ,self%new,1)
  e => self%mem(:,:,:,          self%idx%e ,self%new,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  !-----------------------------------------------------------------------------
  ! Compute update using second-order Godunov method
  !-----------------------------------------------------------------------------
  !self%dq = 0d0 ; self%dq2 = 0d0 ; self%dbf = 0d0
  call uslope(self,self%bf,self%q,self%dq,self%dbf,slope_type)   ! TVD slopes
  if (verbose > 2) write (stdlog,*) &
   'uslope', self%id, minval(self%dbf), maxval(self%dbf)
  !-----------------------------------------------------------------------------
  ! For slope types in the interval 3-4, use the conservative 3.0 for corners
  ! (cf. RAMSES/CPH)
  !-----------------------------------------------------------------------------
  if (slope_type <= 3.0) then
    self%dq2=self%dq
  else
    call uslope(self,self%bf,self%q,self%dq2,self%dbf,3.0)       ! TVD slopes
  endif

  if (verbose > 2) then
    write (stdlog,*) 'uslope', self%id, minval(self%dq), maxval(self%dq)
    write (stdlog,*) 'uslope', self%id, minval(self%dq2), maxval(self%dq2)
  end if

  if (non_ideal%is_used) then
    allocate (fluxambdiff(self%gn(1),self%gn(2),self%gn(3),3))
    allocate (emfambdiff(self%gn(1),self%gn(2),self%gn(3),3))
    allocate (fluxohmdiss(self%gn(1),self%gn(2),self%gn(3),3))
    allocate (emfohmdiss(self%gn(1),self%gn(2),self%gn(3),3))
    fluxambdiff = 0d0; emfambdiff = 0d0; fluxohmdiss = 0d0; emfohmdiss = 0d0;   
    call non_ideal%non_ideal_comp(self%mesh,self%idx,self%bf,self%uin,self%q,&
      fluxambdiff,emfambdiff,fluxohmdiss,emfohmdiss,self%u_max,self%ds,self%dtime,self%gamma)
  end if
  call self%courant_condition                           ! compute dt

  call trace3d(self)                                    ! predictor step

  allocate (fl(self%gn(1),self%gn(2),self%gn(3),self%nv,3))
  allocate (em(self%gn(1),self%gn(2),self%gn(3),3))
  fl = 0.0d0
  em = 0.0d0
  call cmpflxm(self,self%gn(1),2,3,4,6,7,8,fl(:,:,:,:,1))
  call cmpflxm(self,self%gn(2),3,4,2,7,8,6,fl(:,:,:,:,2))
  call cmpflxm(self,self%gn(3),4,2,3,8,6,7,fl(:,:,:,:,3))

  call cmp_mag_flx(self, self%qRT, self%qRB, self%qLT, self%qLB,2,3,4,6,7,8,em(:,:,:,3),self%gn(1),self%nv)
  call cmp_mag_flx(self, self%qRT, self%qLT, self%qRB, self%qLB,4,2,3,8,6,7,em(:,:,:,2),self%gn(1),self%nv)
  call cmp_mag_flx(self, self%qRT, self%qRB, self%qLT, self%qLB,3,4,2,7,8,6,em(:,:,:,1),self%gn(1),self%nv)

  dtdx = self%dtime/self%ds(1) 
  dtdy = self%dtime/self%ds(2) 
  dtdz = self%dtime/self%ds(3) 
  em(:,:,:,1) = em(:,:,:,1)*dtdy
  em(:,:,:,2) = em(:,:,:,2)*dtdz
  em(:,:,:,3) = em(:,:,:,3)*dtdx
  if (self%id==io%id_debug) then
    write (io%output,'("X6",f10.6,i4,1p,10e12.3)') self%time, 6, em(4:13,8,8,1)
    write (io%output,'("X7",f10.6,i4,1p,10e12.3)') self%time, 7, em(4:13,8,8,2)
    write (io%output,'("X8",f10.6,i4,1p,10e12.3)') self%time, 8, em(4:13,8,8,3)
  end if

  if (non_ideal%is_used) then
    call non_ideal%flux_upd(fl,fluxambdiff,fluxohmdiss,self%ds(1),self%dtime,1)
    call non_ideal%flux_upd(fl,fluxambdiff,fluxohmdiss,self%ds(2),self%dtime,2)
    call non_ideal%flux_upd(fl,fluxambdiff,fluxohmdiss,self%ds(3),self%dtime,3)
    call non_ideal%non_ideal_emf_up(em(:,:,:,3),emfambdiff,emfohmdiss,self%ds(3),self%dtime,3)
    call non_ideal%non_ideal_emf_up(em(:,:,:,2),emfambdiff,emfohmdiss,self%ds(2),self%dtime,2)
    call non_ideal%non_ideal_emf_up(em(:,:,:,1),emfambdiff,emfohmdiss,self%ds(1),self%dtime,1)
  end if
  ! Compute fluxes with appropriate time step
  fl(:,:,:,:,1) = fl(:,:,:,:,1) * dtdx
  fl(:,:,:,:,2) = fl(:,:,:,:,2) * dtdy
  fl(:,:,:,:,3) = fl(:,:,:,:,3) * dtdz
  if (self%id==io%id_debug) then
    do i=1,5
      write (io%output,'("X3",f10.6,i4,1p,10e12.3)') self%time, i, fl(4:13,8,8,i,1)
      write (io%output,'("X4",f10.6,i4,1p,10e12.3)') self%time, i, fl(4:13,8,8,i,2)
      write (io%output,'("X5",f10.6,i4,1p,10e12.3)') self%time, i, fl(4:13,8,8,i,3)
    end do
  end if

  !-----------------------------------------------------------------------------
  ! Update conservative variables
  !-----------------------------------------------------------------------------
  associate (i0=>self%mesh(1)%li, j0=>self%mesh(2)%li, k0=>self%mesh(3)%li, &
             i1=>self%mesh(1)%ui, j1=>self%mesh(2)%ui, k1=>self%mesh(3)%ui, &
             uin=>self%uin,unew=>self%unew)
  lock => self%mem_lock(self%new)
                                                        call lock%set ('update')
  nv3 = 5
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,1) - &
    fl(i0+1:i1+1,j0  :j1  ,k0  :k1  ,1:nv3,1) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,2) - &
    fl(i0  :i1  ,j0+1:j1+1,k0  :k1  ,1:nv3,2) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,3) - &
    fl(i0  :i1  ,j0  :j1  ,k0+1:k1+1,1:nv3,3)   

  if (self%id == io%id_debug) then
    write (io%output,'("M1",f10.6,i4,1p,10e12.3)') self%time, i, em(4:13,8,8,1)
    write (io%output,'("M2",f10.6,i4,1p,10e12.3)') self%time, i, em(4:13,8,8,2)
    write (io%output,'("M3",f10.6,i4,1p,10e12.3)') self%time, i, em(4:13,8,8,3)
  end if
  !-----------------------------------------------------------------------------
  ! Update Bx 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,6) =  &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,6) + &
   (( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,2)   - &
      em(i0  :i1  ,j0  :j1  ,k0+1:k1+1,2) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,3)   - &
      em(i0  :i1  ,j0+1:j1+1,k0  :k1  ,3) ))
  !-----------------------------------------------------------------------------
  ! Update By 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,7) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,7) + &
   (( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,3)   - &
      em(i0+1:i1+1,j0  :j1  ,k0  :k1  ,3) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,1)   - &
      em(i0  :i1  ,j0  :j1  ,k0+1:k1+1,1) ))
  !-----------------------------------------------------------------------------
  ! Update Bz 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,8) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,8) + &
   (( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,1)   - &
      em(i0  :i1  ,j0+1:j1+1,k0  :k1  ,1) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,2)   - &
      em(i0+1:i1+1,j0  :j1  ,k0  :k1  ,2) ))
  end associate
  !-----------------------------------------------------------------------------
  ! Add the external force per unit mass as a source term, and add the
  ! corresponding total energy
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%new,1)
  e => self%mem(:,:,:,self%idx%e,self%new,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  do idir=1,3
  do k=self%li(3),self%ui(3)
  do j=self%li(2),self%ui(2)
    if (associated(scratch%force_per_unit_mass)) then
      do i=self%li(1),self%ui(1)
        u0 = p(i,j,k,idir)/d(i,j,k)
        u1 = u0 + self%dtime*scratch%force_per_unit_mass(i,j,k,idir)
        p(i,j,k,idir) = p(i,j,k,idir) + (u1-u0)*d(i,j,k)
        e(i,j,k) = e(i,j,k) + 0.5*(u1**2-u0**2)*d(i,j,k)
      end do
    end if
    if (associated(scratch%force_per_unit_volume)) then
      do i=self%li(1),self%ui(1)
        u0 = p(i,j,k,idir)/d(i,j,k)
        u1 = u0 + self%dtime*scratch%force_per_unit_mass(i,j,k,idir)/d(i,j,k)
        p(i,j,k,idir) = p(i,j,k,idir) + (u1-u0)*d(i,j,k)
        e(i,j,k) = e(i,j,k) + 0.5*(u1**2-u0**2)*d(i,j,k)
      end do
    end if
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Reset the energy variable in isothermal cases
  !-----------------------------------------------------------------------------
  if (isothermal) &
    self%mem(:,:,:,self%idx%e,self%new,1) = &
    self%mem(:,:,:,self%idx%d,self%new,1)
  !.............................................................................
  if (verbose > 1) then
    block
    integer:: i
    do i=1,5 !self%nv
      d => self%mem(:,:,:,i,self%new,1)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') 'unew', self%id, i, &
      self%fminval(d), self%fminloc(d), self%fmaxval(d), self%fmaxloc(d)
    end do
    do i=1,3
      d => scratch%force_per_unit_mass(:,:,:,i)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') ' fpm', self%id, i, &
      self%fminval(d), self%fminloc(d), self%fmaxval(d), self%fmaxloc(d)
    end do
    do i=1,3
      d => fl(:,:,:,5,i)
      write (stdlog,'(a,i6,i4,2(1p,e14.4,3i4))') 'flux', self%id, i, &
        self%fminval(d), self%fminloc(d), self%fmaxval(d), self%fmaxloc(d)
    end do
    end block
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate
  !-----------------------------------------------------------------------------
  call self%counter_update
  if (non_ideal%is_used) &
    deallocate (fluxambdiff,emfambdiff,fluxohmdiss,emfohmdiss)
  deallocate (fl,em,self%q,self%qp,self%qm,self%dq,self%dq2,self%bf,self%dbf, &
              self%qLT,self%qLB,self%qRT,self%qRB,self%c, &
              self%Ex,self%Ey,self%Ez)
  !-----------------------------------------------------------------------------
  ! Make flux corrections, based on flux information from nbors.
  !-----------------------------------------------------------------------------
  call flux%corrections (self)
  !-----------------------------------------------------------------------------
  ! Use one of several methods to clean div(B) from guard zones
  !-----------------------------------------------------------------------------
  select case (divb_cleaner)
  case (0)
  call divb_clean0 (self)
  case (1)
  call divb_clean1 (self)
  case (2)
  call divb_clean2 (self)
  end select
                                                      call lock%unset ('update')
  !.............................................................................
  call trace%end (itimer)
  if (rieman%ndiff > 0) then
    write(io_unit%mpi,*) 'id,ndiff', self%id, rieman%ndiff, rieman%nsolu
    !$omp atomic write
    rieman%ndiff = 0
    !$omp atomic write
    rieman%nsolu = 0
  end if
END SUBROUTINE update

!===============================================================================
!> Determine the guard zone values of the magnetic field components from the
!> div(B)=0 condition.  As appreciated from the expressions below, the only
!> external knowledge that is required are the perpendicular components of the
!> magnetic field along cube edges, and values interpolated from nbors are used
!> there.
!===============================================================================
SUBROUTINE divb_clean0 (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: b1, b2, b3
  integer:: iter, i1, i2, i3, li(3), ui(3), lb(3), ub(3)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean0', itimer=itimer)
  li = self%mesh%li
  ui = self%mesh%ui
  lb = self%mesh%lb
  ub = self%mesh%ub
  b1 => self%mem(:,:,:,self%idx%bx,self%new,1)
  b2 => self%mem(:,:,:,self%idx%by,self%new,1)
  b3 => self%mem(:,:,:,self%idx%bz,self%new,1)
  do iter=1,divb_iter
    call iter_x
    call iter_y
    call iter_z
    call iter_y
  end do
  call trace%end (itimer)
contains
  subroutine iter_x
    !---------------------------------------------------------------------------
    ! Left x-face.  Initially un-constrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=li(1),lb(1),-1
    do i3=lb(3),ub(3)-1
    do i2=lb(2),ub(2)-1
      b1(i1,i2,i3) = b1(i1+1,i2,i3)                 &
                   + b2(i1,i2+1,i3) - b2(i1,i2,i3)  &
                   + b3(i1,i2,i3+1) - b3(i1,i2,i3)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right x-face.  Initially unconstrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=ui(1),ub(1)-1
    do i3=lb(3),ub(3)-1
    do i2=lb(2),ub(2)-1
      b1(i1+1,i2,i3) =                + b1(i1,i2,i3) &
                     - b2(i1,i2+1,i3) + b2(i1,i2,i3) &
                     - b3(i1,i2,i3+1) + b3(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_x
  subroutine iter_y
    !---------------------------------------------------------------------------
    ! Left y-face.  .  Initially unconstrained: b3(y) along edges
    !---------------------------------------------------------------------------
    do i2=li(2),lb(2),-1
    do i3=lb(3),ub(3)-1
    do i1=lb(1),ub(1)-1
      b2(i1,i2,i3) = b2(i1,i2+1,i3)                &
                   + b3(i1,i2,i3+1) - b3(i1,i2,i3) &
                   + b1(i1+1,i2,i3) - b1(i1,i2,i3)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right y-face.
    !---------------------------------------------------------------------------
    do i2=ui(2),ub(2)-1
    do i3=lb(3),ub(3)-1
    do i1=lb(1),ub(1)-1
      b2(i1,i2+1,i3) =                + b2(i1,i2,i3) &
                     - b3(i1,i2,i3+1) + b3(i1,i2,i3) &
                     - b1(i1+1,i2,i3) + b1(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_y
  subroutine iter_z
    !---------------------------------------------------------------------------
    ! Left z-face.  Not known: b3(:,:,l3)
    !---------------------------------------------------------------------------
    do i3=li(3),lb(3),-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3) = b3(i1,i2,i3+1)                &
                   + b1(i1+1,i2,i3) - b1(i1,i2,i3) &
                   + b2(i1,i2+1,i3) - b2(i1,i2,i3)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right z-face: Not known, b3(i1,i2,i3+1)
    !---------------------------------------------------------------------------
    do i3=ui(3),ub(3)-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3+1) =                + b3(i1,i2,i3) &
                     - b1(i1+1,i2,i3) + b1(i1,i2,i3) &
                     - b2(i1,i2+1,i3) + b2(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_z
END SUBROUTINE divb_clean0

!===============================================================================
!> Perform a single iteration of a div(B) clean.  The limit of stability is
!> eps=1./6.  This is intended to prevent build-up of layers of large div(B)
!> between patches. An exact solution would be to solve an internal Laplace
!> equation, with boundary conditiosn = the boundary value of div(B), and then
!> subtract the gradient of the resulting potential.
!>
!> The problem can also be greatly reduced by using a higher order time
!> interpoation scheme, thus making the electric field computed from the
!> guard zone values of B and U a better match to the neighbors E-fields.
!===============================================================================
SUBROUTINE divb_clean1 (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: divb
  real, parameter:: eps=0.11
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean1', itimer=itimer)
  bx => self%mem(:,:,:,self%idx%bx,self%new,1)
  by => self%mem(:,:,:,self%idx%by,self%new,1)
  bz => self%mem(:,:,:,self%idx%bz,self%new,1)
  allocate (divb(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  divb = ddup(bx,1) + ddup(by,2) + ddup(bz,3)
  bx = bx + eps*dddn(divb,1)
  by = by + eps*dddn(divb,2)
  bz = bz + eps*dddn(divb,3)
  deallocate (divb)
  call trace%end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  function ddup(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = cshift(f,1,i) - f
  end function
  !-----------------------------------------------------------------------------
  function dddn(f,i) result (g)
  real:: f(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = f - cshift(f,-1,i)
  end function
END SUBROUTINE divb_clean1

!===============================================================================
!> Enforce div(B)=0 in the guard zones, by succesively imposing the condition,
!> layer-by-layer, on the component perpendicular to each of the patch faces.
!===============================================================================
SUBROUTINE divb_clean2 (self)
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: divb
  real, parameter:: eps=0.11
  integer:: l(3), u(3), ix, iy, iz
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%divb_clean2', itimer=itimer)
  bx => self%mem(:,:,:,self%idx%bx,self%new,1)
  by => self%mem(:,:,:,self%idx%by,self%new,1)
  bz => self%mem(:,:,:,self%idx%bz,self%new,1)
  l = self%mesh%lb
  u = self%mesh%ub-1
  !--- x-direction -------------------------------------------------------------
  do ix=self%mesh(1)%lo,self%mesh(1)%lb,-1
    bx(ix  ,l(2):u(2),l(3):u(3)) = &
    bx(ix+1,l(2):u(2),l(3):u(3)) + ( &
      by(ix,l(2)+1:u(2)+1,l(3)  :u(3)  ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)  :u(2)  ,l(3)+1:u(3)+1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  do ix=self%mesh(1)%ui,self%mesh(1)%ub-1
    bx(ix+1,l(2):u(2),l(3):u(3)) = &
    bx(ix  ,l(2):u(2),l(3):u(3)) - ( &
      by(ix,l(2)+1:u(2)+1,l(3)  :u(3)  ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)  :u(2)  ,l(3)+1:u(3)+1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  !--- y-direction -------------------------------------------------------------
  do iy=self%mesh(2)%lo,self%mesh(2)%lb,-1
    by(l(1):u(1),iy  ,l(3):u(3)) = &
    by(l(1):u(1),iy+1,l(3):u(3)) + ( &
      bx(l(1)+1:u(1)+1,iy,l(3)  :u(3)  ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)  :u(1)  ,iy,l(3)+1:u(3)+1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  do iy=self%mesh(2)%ui,self%mesh(2)%ub-1
    by(l(1):u(1),iy+1,l(3):u(3)) = &
    by(l(1):u(1),iy  ,l(3):u(3)) - ( &
      bx(l(1)+1:u(1)+1,iy,l(3)  :u(3)  ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)  :u(1)  ,iy,l(3)+1:u(3)+1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  !--- z-direction -------------------------------------------------------------
  do iz=self%mesh(3)%lo,self%mesh(3)%lb,-1
    bz(l(1):u(1),l(2):u(2),iz  ) = &
    bz(l(1):u(1),l(2):u(2),iz+1) + ( &
      bx(l(1)+1:u(1)+1,l(2)  :u(2)  ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)  :u(1)  ,l(2)+1:u(2)+1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  do iz=self%mesh(3)%ui,self%mesh(3)%ub-1
    bz(l(1):u(1),l(2):u(2),iz+1) = &
    bz(l(1):u(1),l(2):u(2),iz  ) - ( &
      bx(l(1)+1:u(1)+1,l(2)  :u(2)  ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)  :u(1)  ,l(2)+1:u(2)+1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  call trace%end (itimer)
END SUBROUTINE divb_clean2

!===============================================================================
!===============================================================================
subroutine check_small (self)
  class(mhd_t):: self
  integer:: i,j,k,lo(3),uo(3),n_smalld,n_smallp
  real(dp):: c,d,u,v,w,p,ekin,B1,B2,B3,emag
  integer, save:: nprint=10
  integer, save:: itimer=0
  integer:: it, jt, iw
  logical:: panic
  associate(uu=>self%uin)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%check_small', itimer=itimer, detailed_timer=detailed_timer)
  n_smalld = 0
  n_smallp = 0
  panic = .false.
  lo = self%mesh%lb
  uo = self%mesh%ub-1
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  do i=lo(1),uo(1)
    if (uu(i,j,k,1) < smallr) then
      panic = .true.
      if (nprint>0) then
        nprint = nprint-1
        write (stdlog,'(i7,3i5,1p,2e12.3)') self%id,i,j,k,uu(i,j,k,1),smallr
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
    B1 = 0.5*(uu(i,j,k,6)+uu(i+1,j,k,6))
    B2 = 0.5*(uu(i,j,k,7)+uu(i,j+1,k,7))
    B3 = 0.5*(uu(i,j,k,8)+uu(i,j,k+1,8))
    emag = 0.5*(B1**2+B2**2+B3**2)
    if (isothermal) then
      p = d*self%csound**2
      uu(i,j,k,5) = p+ekin+emag
    else
      p = (uu(i,j,k,5)-ekin)*(gamma-one)
      if (p < smallr*smallc**2) then
        p = smallr*smallc**2
        uu(i,j,k,5) = p/(gamma-one)+ekin+emag
        n_smallp = n_smallp+1
      end if
    end if
  end do
  end do
  end do
  end associate
  !.............................................................................
  if (n_smalld+n_smallp > 0) then
    write (stdlog,*) 'small: id,nr,np', &
      wallclock(),self%id,self%time,n_smalld,n_smallp
  end if
  !.............................................................................
  if (panic .and. self%n_dump > 0) then
    !$omp critical (panic_cr)
    self%n_dump = self%n_dump-1
    write (io_unit%dump) self%id, self%gn, self%nv, self%nt, self%nw, self%istep
    write (io_unit%dump) self%time, self%dtime
    do iw=1,self%nw
    do it=1,self%nt
      jt = self%iit(it)
      write (io_unit%dump) self%mem(:,:,:,:,jt,iw)
    end do
    end do
    write (io_unit%mpi,*) 'density < SMALLR, dumped', self%id, self%time
    !$omp end critical (panic_cr)
  end if
  !.............................................................................
  call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine check_small

!===============================================================================
!> Conservative to primitive variables.  Using scalar scratch variables reduces
!> teh need for a large (non-default) OMP_STACKSIZE
!===============================================================================
subroutine ctoprim(self)
  class(mhd_t):: self
  integer::n(4)
  real(dp):: smallp, c, entho
  real(dp):: eken, emag, etot, eint, B1, B2, B3
  integer, save:: itimer=0
  integer:: i, j, k, loc(3), l(3), u(3), m(3)
  real, pointer:: pg(:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%ctoprim', itimer=itimer, detailed_timer=detailed_timer)
  !-----------------------------------------------------------------------------
  self%uin  => self%mem(:,:,:,:,self%new,1)
  self%unew => self%mem(:,:,:,:,self%new,1)
  m = self%mesh%gn
  allocate (self%c  (m(1),m(2),m(3)))
  allocate (self%q  (m(1),m(2),m(3),self%nv))
  allocate (self%bf (m(1),m(2),m(3),3))
  call check_small (self)
  associate(uin=>self%uin,q=>self%q,bf=>self%bf,dt=>self%dtime,c=>self%c)
  !-----------------------------------------------------------------------------
  if (isothermal) then
    entho = one
  else
    entho = one/max(gamma-one,1e-6)
  end if
  smallp = smallr*smallc**2/self%gamma
  n = shape(uin)-1
  !
  !self%bf => uin(:,:,:,6:8)                                     ! Left  face centered B
  uin(:,:,:,1) = max(uin(:,:,:,1),smallr)                       ! Density
  q(:,:,:,1) = uin(:,:,:,1)                                     ! Density
  q(:,:,:,2) = uin(:,:,:,2)/q(:,:,:,1)                          ! vx
  q(:,:,:,3) = uin(:,:,:,3)/q(:,:,:,1)                          ! vy
  q(:,:,:,4) = uin(:,:,:,4)/q(:,:,:,1)                          ! vz
  q(:,:,:,5) = 0.0                                              ! p_t; initialised for safety
  q(:,:,:,6) = uin(:,:,:,6)                                     ! Bx
  q(:,:,:,7) = uin(:,:,:,7)                                     ! By
  q(:,:,:,8) = uin(:,:,:,8)                                     ! Bz
  bf(:,:,:,1) = uin(:,:,:,6)
  bf(:,:,:,2) = uin(:,:,:,7)
  bf(:,:,:,3) = uin(:,:,:,8)
  if (verbose > 1) then
    write (stdlog,*) 'ctoprim:rh', self%id, minval(q(:,:,:,1)), maxval(q(:,:,:,1))
    write (stdlog,*) 'ctoprim:bf', self%id, minval(bf), maxval(bf)
  end if
  !
  do k=1,self%gn(3)-1
  do j=1,self%gn(2)-1
  do i=1,self%gn(1)-1
    eken = half*(q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4)**2)*q(i,j,k,1)
    B1 = half*(q(i,j,k,6)+q(i+1,j,k,6))
    B2 = half*(q(i,j,k,7)+q(i,j+1,k,7))
    B3 = half*(q(i,j,k,8)+q(i,j,k+1,8))
    emag = half*(B1**2 + B2**2 + B3**2)
    q(i,j,k,5) = (uin(i,j,k,5) - emag - eken)*(gamma-1d0)
  end do
  end do
  end do
  !.............................................................................
  pg => q(:,:,:,5)
  l = self%mesh%li
  u = self%mesh%ui
  if (any(pg(l(1):u(1),l(2):u(2),l(3):u(3)) < 0.0) .or. verbose > 0) then
    write (stdlog,'(a,i6,1p,2e11.2,3i4)') 'ctoprim: max(pg),min(pg),loc=', &
      self%id, self%fmaxval(pg), self%fminval(pg), self%fminloc(pg)
  end if
  !.............................................................................
  do k=1,self%gn(3)-1
  do j=1,self%gn(2)-1
  do i=1,self%gn(1)-1
    q(i,j,k,5) = merge(q(i,j,k,1)*self%csound**2, max(q(i,j,k,5),smallp), isothermal)
    c(i,j,k) = sqrt((q(i,j,k,5)*gamma + 2.0*emag)/q(i,j,k,1))
  end do
  end do
  end do
  if (self%id == io%id_debug) then
    do i=1,8
      write (io%output,'("X1",f10.6,i4,1p,10e12.3)') self%time, i, uin(4:13,8,8,i)
      write (io%output,'("X2",f10.6,i4,1p,10e12.3)') self%time, i,   q(4:13,8,8,i)
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Alternative Courant condition formulations
  !-----------------------------------------------------------------------------
  if (courant_type==1) then
    c(:,:,:) = c(:,:,:) + max(abs(q(:,:,:,2)),  abs(q(:,:,:,3)),  abs(q(:,:,:,4)))
  else
    c(:,:,:) = c(:,:,:) +    (abs(q(:,:,:,2)) + abs(q(:,:,:,3)) + abs(q(:,:,:,4)))/3.0
  end if
  l = self%mesh%lo
  u = self%mesh%uo
  if (any(l < 0 .or. u > self%gn)) then
    write (stderr,*) 'ERROR: lo, uo =', self%id, l, u
    call io%abort ('ctoprim: lo or uo inconsistent')
  end if
  self%u_max = maxval(c(l(1):u(1),l(2):u(2),l(3):u(3)))
  if (self%get_umax_location) then
    loc = maxloc(c(l(1):u(1),l(2):u(2),l(3):u(3))) + l-1
    write (stdlog,'(a,i6,f12.6,2(1x,3i4),3f9.4,2(1x,a,1p,4g10.2))') 'u_max:', &
      self%id, self%time, self%ipos, loc, &
      (self%mesh(1)%p + self%mesh(1)%r(loc(1))-self%mesh(2)%llc_cart)/self%mesh(1)%b, &
      (self%mesh(2)%p + self%mesh(2)%r(loc(2))-self%mesh(2)%llc_cart)/self%mesh(2)%b, &
      (self%mesh(3)%p + self%mesh(3)%r(loc(3))-self%mesh(2)%llc_cart)/self%mesh(3)%b, &
      'P,B:', q(loc(1),loc(2),loc(3),5:8), &
      'D,U:', q(loc(1),loc(2),loc(3),1:4)
    flush (stdlog)
  end if
  end associate
  call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine ctoprim

!===============================================================================
!> Slope limiters
!===============================================================================
subroutine uslope(self,bf,q,dq,dbf,rslope_type)
  class(mhd_t):: self
  real, contiguous, dimension(:,:,:,:):: q
  real, contiguous, dimension(:,:,:,:):: bf
  real, contiguous, dimension(:,:,:,:,:):: dq
  real, contiguous, dimension(:,:,:,:,:):: dbf
  real:: rslope_type
  integer::nv, islope_type
  integer::i, j, k, n, l(3), u(3)
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, xslope_type
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,idff
  integer,save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%uslope', itimer=itimer, detailed_timer=detailed_timer)
  !
  islope_type = int(rslope_type)
  l = self%mesh%lb+1
  u = self%mesh%ub-1
  nv=self%nv
  if(slope_type==-1)then
    do n = 1, nv
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
     return
  else if (islope_type==1) then  ! minmod
    do n = 1, nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd private(dlft,drgt)
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
    xslope_type= 2 !MIN(slope_type,2)
    do n = 1, nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd private(dlft,drgt,dcen,dsgn,slop)
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = (q(i  ,j,k,n) - q(i-1,j,k,n))          ! = dq
            drgt = (q(i+1,j,k,n) - q(i  ,j,k,n))          ! = dq
            dcen = half*(dlft+drgt)                       ! = dq
            dsgn = sign(one, dcen)                        ! sign(slope)
            slop = xslope_type*min(abs(dlft),abs(drgt))   ! slope_type=1: slope=dq, slope_type=2: slope=2*dq 
            slop = merge(slop,zero,(dlft*drgt)>zero)      ! if both are non-zero, use the smaller one, or times 2
            dq(i,n,1,j,k) = dsgn*min(slop,abs(dcen))      ! the smaller of the twice the smaller slope, and the central one            ! slopes in second coordinate direction
            dlft = (q(i,j  ,k,n) - q(i,j-1,k,n))
            drgt = (q(i,j+1,k,n) - q(i,j  ,k,n))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dq(i,n,2,j,k) = dsgn*min(slop,abs(dcen))
            ! slopes in third coordinate direction
            dlft = (q(i,j,k  ,n) - q(i,j,k-1,n))
            drgt = (q(i,j,k+1,n) - q(i,j,k  ,n))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dq(i,n,3,j,k) = dsgn*min(slop,abs(dcen))
          end do
        end do
      end do
    end do

    do k = l(3),u(3)
      do j = l(2),u(2)
        !dir$ vector aligned
        !dir$ ivdep
        !$omp simd private(dlft,drgt,dcen,dsgn,slop)
        do i = l(1),u(1)
          ! Bx
          dlft = (bf(i,j  ,k,1) - bf(i,j-1,k,1))
          drgt = (bf(i,j+1,k,1) - bf(i,j  ,k,1))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,1,1,j,k) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j,k  ,1) - bf(i,j,k-1,1))
          drgt = (bf(i,j,k+1,1) - bf(i,j,k  ,1))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,1,2,j,k) = dsgn*min(slop,abs(dcen))
          ! By
          ! slopes in first coordinate direction
          dlft = (bf(i  ,j,k,2) - bf(i-1,j,k,2))
          drgt = (bf(i+1,j,k,2) - bf(i  ,j,k,2))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,2,1,j,k) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j,k  ,2) - bf(i,j,k-1,2))
          drgt = (bf(i,j,k+1,2) - bf(i,j,k  ,2))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,2,2,j,k) = dsgn*min(slop,abs(dcen))
          ! Bz
          ! slopes in first coordinate direction
          dlft = (bf(i  ,j,k,3) - bf(i-1,j,k,3))
          drgt = (bf(i+1,j,k,3) - bf(i  ,j,k,3))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,3,1,j,k) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j  ,k,3) - bf(i,j-1,k,3))
          drgt = (bf(i,j+1,k,3) - bf(i,j  ,k,3))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,3,2,j,k) = dsgn*min(slop,abs(dcen))
        end do
      end do
    end do

  else if (islope_type==3) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 2_dp,kind=dp)
    do n = 1, nv
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
      do k = l(3),u(3)
        do j = l(2),u(2)
          !dir$ ivdep
          !$omp simd private(dlft,drgt,dcen,dsgn,slop)
          do i = l(1),u(1)
            ! Bx
            dlft = (bf(i,j  ,k,1) - bf(i,j-1,k,1))
            drgt = (bf(i,j+1,k,1) - bf(i,j  ,k,1))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,1,1,j,k) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j,k  ,1) - bf(i,j,k-1,1))
            drgt = (bf(i,j,k+1,1) - bf(i,j,k  ,1))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,1,2,j,k) = dsgn*min(slop,abs(dcen))
            !By
            ! slopes in first coordinate direction
            dlft = (bf(i  ,j,k,2) - bf(i-1,j,k,2))
            drgt = (bf(i+1,j,k,2) - bf(i  ,j,k,2))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,2,1,j,k) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j,k  ,2) - bf(i,j,k-1,2))
            drgt = (bf(i,j,k+1,2) - bf(i,j,k  ,2))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,2,2,j,k) = dsgn*min(slop,abs(dcen))
            ! Bz
            ! slopes in first coordinate direction
            dlft = (bf(i  ,j,k,3) - bf(i-1,j,k,3))
            drgt = (bf(i+1,j,k,3) - bf(i  ,j,k,3))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,3,1,j,k) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j  ,k,3) - bf(i,j-1,k,3))
            drgt = (bf(i,j+1,k,3) - bf(i,j  ,k,3))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,3,2,j,k) = dsgn*min(slop,abs(dcen))
          end do
        end do
     end do

  else
     write(*,*)'Unknown slope type'
     stop
  endif
  call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine uslope

!===============================================================================
!> Predictor step
!===============================================================================
SUBROUTINE trace3d(self)
  IMPLICIT NONE
  class(mhd_t):: self
  real(dp):: dx, dy, dz, dt
  integer :: i, j, k, n
  integer, parameter ::ir=1, iu=2, iv=3, iw=4, ip=5, ia=6, ib=7, ic=8
  real(dp)::dtdx, dtdy, dtdz, smallp
  real(dp) :: invr
  real(dp) :: r, u, v, w, p, A, B, C
  real(dp) :: ELL, ELR, ERL, ERR
  real(dp) :: FLL, FLR, FRL, FRR
  real(dp) :: GLL, GLR, GRL, GRR
  real(dp) :: drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  real(dp) :: dry, duy, dvy, dwy, dpy, dAy, dbY, dCy
  real(dp) :: drz, duz, dvz, dwz, dpz, dAz, dBz, dCz
  real(dp) :: drx2,dpx2,dry2,dpy2,drz2,dpz2
  real(dp) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
  real(dp) :: AL, AR, BL, BR, CL, CR
  real(dp) :: dALy, dARy, dALz, dARz
  real(dp) :: dBLx, dBRx, dBLz, dBRz
  real(dp) :: dCLx, dCRx, dCLy, dCRy
  real(dp) :: sAL0, sAR0, sBL0, sBR0, sCL0, sCR0
  real(dp), dimension(:,:,:), allocatable:: rho_t, p_t
  real(dp) :: eff_gamma, off
  real(dp) :: ff(self%gn(1),3), qq(self%gn(1))
  real(8):: dt2
  integer:: lo(3), uo(3)
  integer, save:: itimer(4)=0
  !-----------------------------------------------------------------------------
  real(dp), dimension(:,:,:)    , pointer, contiguous:: Bx, By, Bz, Ux, Uy, Uz
  real(dp), dimension(:,:,:)    , pointer, contiguous:: Ex, Ey, Ez, d
  real(dp), dimension(:,:,:,:)  , pointer, contiguous:: q, bf, em, uin, unew
  real(dp), dimension(:,:,:,:,:), pointer, contiguous:: qp, qm, dq, dq2, dbf
  real(dp), dimension(:,:,:,:,:), pointer, contiguous:: qLB, qRB, qLT, qRT
  !-----------------------------------------------------------------------------
  !associate(q=>self%q,bf=>self%bf,dq=>self%dq,dq2=>self%dq2,qm=>self%qm,qp=>self%qp, &
  !  dbf=>self%dbf,qRT=>self%qRT,qRB=>self%qRB,qLT=>self%qLT,qLB=>self%qLB, &
  !  Bx=>self%bf(:,:,:,1),By=>self%bf(:,:,:,2),Bz=>self%bf(:,:,:,3), &
  !  Ux=>self%q(:,:,:,2),Uy=>self%q(:,:,:,3),Uz=>self%q(:,:,:,4), &
  !  Ex=>self%Ex,Ey=>self%Ey,Ez=>self%Ez)
  !-----------------------------------------------------------------------------
  uin => self%uin
  unew=> self%unew
  q   => self%q
  bf  => self%bf
  dq  => self%dq
  dq2 => self%dq2
  qm  => self%qm
  qp  => self%qp
  Bx  => self%bf(:,:,:,1)
  By  => self%bf(:,:,:,2)
  Bz  => self%bf(:,:,:,3)
  Ux  => self%q(:,:,:,2)
  Uy  => self%q(:,:,:,3)
  Uz  => self%q(:,:,:,4)
  Ex  => self%Ex
  Ey  => self%Ey
  Ez  => self%Ez
  dbf => self%dbf
  qRT => self%qRT
  qRB => self%qRB
  qLT => self%qLT
  qLB => self%qLB
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%trace3d', itimer=itimer(1), detailed_timer=detailed_timer)
  su0=0.0_dp; sv0=0.0_dp; sw0=0.0_dp
  !
  dtdx = self%dtime/self%ds(1)
  dtdy = self%dtime/self%ds(2)
  dtdz = self%dtime/self%ds(3)
  smallp = smallr*smallc**2/gamma
  !
  ! E = - UxB = BxU
  Ex = ydn(zdn(Uy))*ydn(Bz) - zdn(ydn(Uz))*zdn(By) 
  Ey = zdn(xdn(Uz))*zdn(Bx) - xdn(zdn(Ux))*xdn(Bz) 
  Ez = xdn(ydn(Ux))*xdn(By) - ydn(xdn(Uy))*ydn(Bx) 
  !
  !call trace%end (itimer(1))
  !call trace%begin ('mhd_t%trace3d(2)', itimer=itimer(2))
  eff_gamma = gamma
  !
  allocate (rho_t(1:self%gn(1),1:self%gn(2),1:self%gn(3)))
  allocate (  p_t(1:self%gn(1),1:self%gn(2),1:self%gn(3)))
  !
  ! Fill rho-t & p_t with non-zero values at guard zone edges,
  ! to ! avoid arithmetic errors (internal values are set in loop)
  rho_t(1,:,:) = q(1,:,:,ir)
  rho_t(:,1,:) = q(:,1,:,ir)
  rho_t(:,:,1) = q(:,:,1,ir)
  p_t  (1,:,:)   = q(1,:,:,ip)
  p_t  (:,1,:)   = q(:,1,:,ip)
  p_t  (:,:,1)   = q(:,:,1,ip)
  n = self%gn(1)
  rho_t(n,:,:) = q(n,:,:,ir)
  p_t  (n,:,:)   = q(n,:,:,ip)
  n = self%gn(2)
  rho_t(:,n,:) = q(:,n,:,ir)
  p_t  (:,n,:)   = q(:,n,:,ip)
  n = self%gn(3)
  rho_t(:,:,n) = q(:,:,n,ir)
  p_t  (:,:,n)   = q(:,:,n,ip)
  !-----------------------------------------------------------------------------
  ! This loop needs to produce data in the interval (lo-1,uo+1), since the loop
  ! below accesses {i,j,k}+-1, and needs to be good in (lo,uo).  With ng=3, the
  ! loop limits are (lo-1,uo+1) = (2,gn-1)
  !-----------------------------------------------------------------------------
  dt2 = 0.5d0*self%dtime
  lo = self%mesh%lo-1
  uo = self%mesh%uo+1
  d => self%mem(:,:,:,self%idx%d,self%new,1)
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  ff(:,:) = 0.0_dp
  if (associated(scratch%force_per_unit_mass)) then
    !dir$ simd
    do i=2,self%gn(2)-1
      ff(i,1) = ff(i,1) + scratch%force_per_unit_mass(i,j,k,1)*dt2
      ff(i,2) = ff(i,2) + scratch%force_per_unit_mass(i,j,k,2)*dt2
      ff(i,3) = ff(i,3) + scratch%force_per_unit_mass(i,j,k,3)*dt2
    end do
  end if
  if (associated(scratch%force_per_unit_volume)) then
    !dir$ simd
    do i=lo(1),uo(1)
      ff(i,1) = ff(i,1) + scratch%force_per_unit_mass(i,j,k,1)*dt2/d(i,j,k)
      ff(i,2) = ff(i,2) + scratch%force_per_unit_mass(i,j,k,2)*dt2/d(i,j,k)
      ff(i,3) = ff(i,3) + scratch%force_per_unit_mass(i,j,k,3)*dt2/d(i,j,k)
    end do
  end if
  qq(:) = 0.0
  if (associated(scratch%heating_per_unit_volume)) then
    !dir$ simd
    do i=lo(1),uo(1)
      qq(i) = qq(i) + scratch%heating_per_unit_volume(i,j,k)
    end do
  end if
  if (associated(scratch%heating_per_unit_mass)) then
    !dir$ simd
    do i=lo(1),uo(1)
      qq(i) = qq(i) + scratch%heating_per_unit_mass(i,j,k)*d(i,j,k)
    end do
  end if
  !dir$ ivdep
  !$omp simd private( &
  !$omp & A,AL,AR,B,BL,BR,C,CL,CR,ELL,ELR,ERL,ERR,FLL,FLR,FRL,FRR,GLL,GLR,GRL,GRR, &
  !$omp & dALy,dALz,dARy,dARz,dAy,dAz,dBLx,dBLz,dBRx,dBRz,dBx,dBz,dCLx,dCLy,dCRx, &
  !$omp & dCRy,dCx,dCy,dpx,dpx2,dpy,dpy2,dpz,dpz2,drx,drx2,dry,dry2,drz,drz2,dux, &
  !$omp & duy,duz,dvx,dvy,dvz,dwx,dwy,dwz,invr,p,r,sAL0,sAR0,sBL0,sBR0,sCL0,sCR0, &
  !$omp & sp0,sr0,su0,sv0,sw0,u,v,w)
  do i=lo(1),uo(1)
    r=q(i,j,k,ir)
    u=q(i,j,k,iu)
    v=q(i,j,k,iv)
    w=q(i,j,k,iw)
    p=q(i,j,k,ip)
    A=q(i,j,k,iA)
    B=q(i,j,k,iB)
    C=q(i,j,k,iC)
    ! Face centered variables
    AL = bf(i,j,k,1)
    BL = bf(i,j,k,2)
    CL = bf(i,j,k,3)
    AR = bf(i+1,j  ,k  ,1)
    BR = bf(i  ,j+1,k  ,2)
    CR = bf(i  ,j  ,k+1,3)
    ! Edge centered electric field in X, Y and Z directions
    ELL = Ex(i,j  ,k  )
    ELR = Ex(i,j  ,k+1)
    ERL = Ex(i,j+1,k  )
    ERR = Ex(i,j+1,k+1)
    FLL = Ey(i  ,j,k  )
    FLR = Ey(i  ,j,k+1)
    FRL = Ey(i+1,j,k  )
    FRR = Ey(i+1,j,k+1)
    GLL = Ez(i  ,j  ,k)
    GLR = Ez(i  ,j+1,k)
    GRL = Ez(i+1,j  ,k)
    GRR = Ez(i+1,j+1,k)
    !
    ! Cell centered TVD slopes in X, Y and Z directions
    drx = half * dq(i,ir,1,j,k)
    drx2= half *dq2(i,ir,1,j,k)
    dux = half * dq(i,iu,1,j,k)
    dvx = half * dq(i,iv,1,j,k)
    dwx = half * dq(i,iw,1,j,k)
    dpx = half * dq(i,ip,1,j,k)
    dpx2= half *dq2(i,ip,1,j,k)
    dBx = half * dq(i,iB,1,j,k)
    dCx = half * dq(i,iC,1,j,k)
    !
    dry = half * dq(i,ir,2,j,k)
    dry2= half *dq2(i,ir,2,j,k)
    duy = half * dq(i,iu,2,j,k)
    dvy = half * dq(i,iv,2,j,k)
    dwy = half * dq(i,iw,2,j,k)
    dpy = half * dq(i,ip,2,j,k)
    dpy2= half *dq2(i,ip,2,j,k)
    duy = half * dq(i,iu,2,j,k)
    dAy = half * dq(i,iA,2,j,k)
    dCy = half * dq(i,iC,2,j,k)
    !
    drz = half * dq(i,ir,3,j,k)
    drz2= half *dq2(i,ir,3,j,k)
    duz = half * dq(i,iu,3,j,k)
    dvz = half * dq(i,iv,3,j,k)
    dwz = half * dq(i,iw,3,j,k)
    dpz = half * dq(i,ip,3,j,k)
    dpz2= half *dq2(i,ip,3,j,k)
    dAz = half * dq(i,iA,3,j,k)
    dBz = half * dq(i,iB,3,j,k)
    !
    ! Face centered TVD slopes in transverse directions
    dALy = half * dbf(i  ,1,1,j  ,k  )
    dARy = half * dbf(i+1,1,1,j  ,k  )
    dALz = half * dbf(i  ,1,2,j  ,k  )
    dARz = half * dbf(i+1,1,2,j  ,k  )
    !
    dBLx = half * dbf(i  ,2,1,j  ,k  )
    dBRx = half * dbf(i  ,2,1,j+1,k  )
    dBLz = half * dbf(i  ,2,2,j  ,k  )
    dBRz = half * dbf(i  ,2,2,j+1,k  )
    !
    dCLx = half * dbf(i  ,3,1,j  ,k  )
    dCRx = half * dbf(i  ,3,1,j  ,k+1)
    dCLy = half * dbf(i  ,3,2,j  ,k  )
    dCRy = half * dbf(i  ,3,2,j  ,k+1)
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
    sp0 = (-u*dpx-dux*eff_gamma*p)*dtdx &
        + (-v*dpy-dvy*eff_gamma*p)*dtdy &
        + (-w*dpz-dwz*eff_gamma*p)*dtdz
    !
    ! Add external force
    su0 = su0 + ff(i,1)
    sv0 = sv0 + ff(i,2)
    sw0 = sw0 + ff(i,3)
    !
    ! Add external heating
    sp0 = sp0  + (eff_gamma-1d0)*qq(i)*half*self%dtime
    !
    ! Evolve entropy instead of pressure
    !sS0 = (-u*dSx)*dtdx + (-v*dSy)*dtdy + (-w*dSz)*dtdz
    !
!    ! Save predicted cell-centered pseudo-entropy; we are borrowing space in q(:,:,:,:,ip),
!    ! which is not used below this point.  In the calling routine this is again copied over
!    ! into uin, and is thus made available to godfine1, through that array.
!    if (entropy_fix) then
!       if ((i==1.or.i==2).and.(j==1.or.j==2).and.(k==1.or.k==2)) then
!          q(:,i,j,k,ip) = S + 2.0*sS0
!          if (verbose > 0) then
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
    rho_t(i,j,k) = r
    u = u + su0
    v = v + sv0
    w = w + sw0
    p = p + sp0
    p = max(p,smallp)
    p_t(i,j,k) = p
    A = half*(AL+AR)
    B = half*(BL+BR)
    C = half*(CL+CR)
    !
    ! Face averaged right state at left interface
    qp(i,ir,1,j,k) = r - drx
    qp(i,iu,1,j,k) = u - dux
    qp(i,iv,1,j,k) = v - dvx
    qp(i,iw,1,j,k) = w - dwx
    qp(i,ip,1,j,k) = p - dpx
    qp(i,iA,1,j,k) = AL
    qp(i,iB,1,j,k) = B - dBx
    qp(i,iC,1,j,k) = C - dCx
    !
    ! Face averaged left state at right interface
    qm(i,ir,1,j,k) = r + drx
    qm(i,iu,1,j,k) = u + dux
    qm(i,iv,1,j,k) = v + dvx
    qm(i,iw,1,j,k) = w + dwx
    qm(i,ip,1,j,k) = p + dpx
    qm(i,iA,1,j,k) = AR
    qm(i,iB,1,j,k) = B + dBx
    qm(i,iC,1,j,k) = C + dCx
    !
    ! Face averaged top state at bottom interface
    qp(i,ir,2,j,k) = r - dry
    qp(i,iu,2,j,k) = u - duy
    qp(i,iv,2,j,k) = v - dvy
    qp(i,iw,2,j,k) = w - dwy
    qp(i,ip,2,j,k) = p - dpy
    qp(i,iA,2,j,k) = A - dAy
    qp(i,iB,2,j,k) = BL
    qp(i,iC,2,j,k) = C - dCy
    !
    ! Face averaged bottom state at top interface
    qm(i,ir,2,j,k) = r + dry
    qm(i,iu,2,j,k) = u + duy
    qm(i,iv,2,j,k) = v + dvy
    qm(i,iw,2,j,k) = w + dwy
    qm(i,ip,2,j,k) = p + dpy
    qm(i,iA,2,j,k) = A + dAy
    qm(i,iB,2,j,k) = BR
    qm(i,iC,2,j,k) = C + dCy
    !
    ! Face averaged front state at back interface
    qp(i,ir,3,j,k) = r - drz
    qp(i,iu,3,j,k) = u - duz
    qp(i,iv,3,j,k) = v - dvz
    qp(i,iw,3,j,k) = w - dwz
    qp(i,ip,3,j,k) = p - dpz
    qp(i,iA,3,j,k) = A - dAz
    qp(i,iB,3,j,k) = B - dBz
    qp(i,iC,3,j,k) = CL
    !
    ! Face averaged back state at front interface
    qm(i,ir,3,j,k) = r + drz
    qm(i,iu,3,j,k) = u + duz
    qm(i,iv,3,j,k) = v + dvz
    qm(i,iw,3,j,k) = w + dwz
    qm(i,ip,3,j,k) = p + dpz
    qm(i,iA,3,j,k) = A + dAz
    qm(i,iB,3,j,k) = B + dBz
    qm(i,iC,3,j,k) = CR
    !
    ! X-edge averaged right-top corner state (RT->LL)
    qRT(i,ir,1,j,k) = r + (+dry2+drz2)
    qRT(i,iu,1,j,k) = u + (+duy+duz)
    qRT(i,iv,1,j,k) = v + (+dvy+dvz)
    qRT(i,iw,1,j,k) = w + (+dwy+dwz)
    qRT(i,ip,1,j,k) = p + (+dpy2+dpz2)
    qRT(i,iA,1,j,k) = A + (+dAy+dAz)
    qRT(i,iB,1,j,k) = BR+ (   +dBRz)
    qRT(i,iC,1,j,k) = CR+ (+dCRy   )
    !
    ! X-edge averaged right-bottom corner state (RB->LR)
    qRB(i,ir,1,j,k) = r + (+dry2-drz2)
    qRB(i,iu,1,j,k) = u + (+duy-duz)
    qRB(i,iv,1,j,k) = v + (+dvy-dvz)
    qRB(i,iw,1,j,k) = w + (+dwy-dwz)
    qRB(i,ip,1,j,k) = p + (+dpy2-dpz2)
    qRB(i,iA,1,j,k) = A + (+dAy-dAz)
    qRB(i,iB,1,j,k) = BR+ (   -dBRz)
    qRB(i,iC,1,j,k) = CL+ (+dCLy   )
    !
    ! X-edge averaged left-top corner state (LT->RL)
    qLT(i,ir,1,j,k) = r + (-dry2+drz2)
    qLT(i,iu,1,j,k) = u + (-duy+duz)
    qLT(i,iv,1,j,k) = v + (-dvy+dvz)
    qLT(i,iw,1,j,k) = w + (-dwy+dwz)
    qLT(i,ip,1,j,k) = p + (-dpy2+dpz2)
    qLT(i,iA,1,j,k) = A + (-dAy+dAz)
    qLT(i,iB,1,j,k) = BL+ (   +dBLz)
    qLT(i,iC,1,j,k) = CR+ (-dCRy   )
    !
    ! X-edge averaged left-bottom corner state (LB->RR)
    qLB(i,ir,1,j,k) = r + (-dry2-drz2)
    qLB(i,iu,1,j,k) = u + (-duy-duz)
    qLB(i,iv,1,j,k) = v + (-dvy-dvz)
    qLB(i,iw,1,j,k) = w + (-dwy-dwz)
    qLB(i,ip,1,j,k) = p + (-dpy2-dpz2)
    qLB(i,iA,1,j,k) = A + (-dAy-dAz)
    qLB(i,iB,1,j,k) = BL+ (   -dBLz)
    qLB(i,iC,1,j,k) = CL+ (-dCLy   )
    !
    ! Y-edge averaged right-top corner state (RT->LL)
    qRT(i,ir,2,j,k) = r + (+drx2+drz2)
    qRT(i,iu,2,j,k) = u + (+dux+duz)
    qRT(i,iv,2,j,k) = v + (+dvx+dvz)
    qRT(i,iw,2,j,k) = w + (+dwx+dwz)
    qRT(i,ip,2,j,k) = p + (+dpx2+dpz2)
    qRT(i,iA,2,j,k) = AR+ (   +dARz)
    qRT(i,iB,2,j,k) = B + (+dBx+dBz)
    qRT(i,iC,2,j,k) = CR+ (+dCRx   )
    !
    ! Y-edge averaged right-bottom corner state (RB->LR)
    qRB(i,ir,2,j,k) = r + (+drx2-drz2)
    qRB(i,iu,2,j,k) = u + (+dux-duz)
    qRB(i,iv,2,j,k) = v + (+dvx-dvz)
    qRB(i,iw,2,j,k) = w + (+dwx-dwz)
    qRB(i,ip,2,j,k) = p + (+dpx2-dpz2)
    qRB(i,iA,2,j,k) = AR+ (   -dARz)
    qRB(i,iB,2,j,k) = B + (+dBx-dBz)
    qRB(i,iC,2,j,k) = CL+ (+dCLx   )
    !
    ! Y-edge averaged left-top corner state (LT->RL)
    qLT(i,ir,2,j,k) = r + (-drx2+drz2)
    qLT(i,iu,2,j,k) = u + (-dux+duz)
    qLT(i,iv,2,j,k) = v + (-dvx+dvz)
    qLT(i,iw,2,j,k) = w + (-dwx+dwz)
    qLT(i,ip,2,j,k) = p + (-dpx2+dpz2)
    qLT(i,iA,2,j,k) = AL+ (   +dALz)
    qLT(i,iB,2,j,k) = B + (-dBx+dBz)
    qLT(i,iC,2,j,k) = CR+ (-dCRx   )
    !
    ! Y-edge averaged left-bottom corner state (LB->RR)
    qLB(i,ir,2,j,k) = r + (-drx2-drz2)
    qLB(i,iu,2,j,k) = u + (-dux-duz)
    qLB(i,iv,2,j,k) = v + (-dvx-dvz)
    qLB(i,iw,2,j,k) = w + (-dwx-dwz)
    qLB(i,ip,2,j,k) = p + (-dpx2-dpz2)
    qLB(i,iA,2,j,k) = AL+ (   -dALz)
    qLB(i,iB,2,j,k) = B + (-dBx-dBz)
    qLB(i,iC,2,j,k) = CL+ (-dCLx   )
    !
    ! Z-edge averaged right-top corner state (RT->LL)
    qRT(i,ir,3,j,k) = r + (+drx2+dry2)
    qRT(i,iu,3,j,k) = u + (+dux+duy)
    qRT(i,iv,3,j,k) = v + (+dvx+dvy)
    qRT(i,iw,3,j,k) = w + (+dwx+dwy)
    qRT(i,ip,3,j,k) = p + (+dpx2+dpy2)
    qRT(i,iA,3,j,k) = AR+ (   +dARy)
    qRT(i,iB,3,j,k) = BR+ (+dBRx   )
    qRT(i,iC,3,j,k) = C + (+dCx+dCy)
    !
    ! Z-edge averaged right-bottom corner state (RB->LR)
    qRB(i,ir,3,j,k) = r + (+drx2-dry2)
    qRB(i,iu,3,j,k) = u + (+dux-duy)
    qRB(i,iv,3,j,k) = v + (+dvx-dvy)
    qRB(i,iw,3,j,k) = w + (+dwx-dwy)
    qRB(i,ip,3,j,k) = p + (+dpx2-dpy2)
    qRB(i,iA,3,j,k) = AR+ (   -dARy)
    qRB(i,iB,3,j,k) = BL+ (+dBLx   )
    qRB(i,iC,3,j,k) = C + (+dCx-dCy)
    !
    ! Z-edge averaged left-top corner state (LT->RL)
    qLT(i,ir,3,j,k) = r + (-drx2+dry2)
    qLT(i,iu,3,j,k) = u + (-dux+duy)
    qLT(i,iv,3,j,k) = v + (-dvx+dvy)
    qLT(i,iw,3,j,k) = w + (-dwx+dwy)
    qLT(i,ip,3,j,k) = p + (-dpx2+dpy2)
    qLT(i,iA,3,j,k) = AL+ (   +dALy)
    qLT(i,iB,3,j,k) = BR+ (-dBRx   )
    qLT(i,iC,3,j,k) = C + (-dCx+dCy)
    !
    ! Z-edge averaged left-bottom corner state (LB->RR)
    qLB(i,ir,3,j,k) = r + (-drx2-dry2)
    qLB(i,iu,3,j,k) = u + (-dux-duy)
    qLB(i,iv,3,j,k) = v + (-dvx-dvy)
    qLB(i,iw,3,j,k) = w + (-dwx-dwy)
    qLB(i,ip,3,j,k) = p + (-dpx2-dpy2)
    qLB(i,iA,3,j,k) = AL+ (   -dALy)
    qLB(i,iB,3,j,k) = BL+ (-dBLx   )
    qLB(i,iC,3,j,k) = C + (-dCx-dCy)
  end do
  end do
  end do
  !
  !call trace%end (itimer(2))
  !call trace%begin ('mhd_t%trace3d(3)', itimer=itimer(3))
  !-----------------------------------------------------------------------------
  ! This loop needs to produce data in the interval (lo,uo), since the loop
  ! in cmp_mag_flux uses {i,j,k}-1
  !-----------------------------------------------------------------------------
  lo = self%mesh%lo-1
  uo = self%mesh%uo+1
  do k=lo(3),uo(3)
  do j=lo(2),uo(2)
  !$omp simd private(drx,dry,drz,dpx,dpy,dpz)
  !dir$ ivdep
  do i=lo(1),uo(1)
    r = rho_t(i,j,k)
    p = p_t(i,j,k)
    !
    ! X-edge averaged right-top corner state (RT->LL)
    !qRT(i,ir,1,j,k) = r + (+dry+drz)
    drx = sqrt(sqrt(r*rho_t(i,j+1,k)*rho_t(i,j,k+1)*rho_t(i,j+1,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j+1,k)*  p_t(i,j,k+1)*  p_t(i,j+1,k+1))) - p
    dry = qRT(i,ir,1,j,k) - r
    dpy = qRT(i,ip,1,j,k) - p
    qRT(i,ir,1,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,ip,1,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,ir,1,j,k) = r + (+dry-drz)
    drx = sqrt(sqrt(r*rho_t(i,j+1,k)*rho_t(i,j,k-1)*rho_t(i,j+1,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j+1,k)*  p_t(i,j,k-1)*  p_t(i,j+1,k-1))) - p
    dry = qRB(i,ir,1,j,k) - r
    dpy = qRB(i,ip,1,j,k) - p
    qRB(i,ir,1,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,ip,1,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged left-top corner state (LT->RL)
    !qLT(i,ir,1,j,k) = r + (-dry+drz)
    drx = sqrt(sqrt(r*rho_t(i,j-1,k)*rho_t(i,j,k+1)*rho_t(i,j-1,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j-1,k)*  p_t(i,j,k+1)*  p_t(i,j-1,k+1))) - p
    dry = qLT(i,ir,1,j,k) - r
    dpy = qLT(i,ip,1,j,k) - p
    qLT(i,ir,1,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,ip,1,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,ir,1,j,k) = r + (-dry-drz)
    drx = sqrt(sqrt(r*rho_t(i,j-1,k)*rho_t(i,j,k-1)*rho_t(i,j-1,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j-1,k)*  p_t(i,j,k-1)*  p_t(i,j-1,k-1))) - p
    dry = qLB(i,ir,1,j,k) - r
    dpy = qLB(i,ip,1,j,k) - p
    qLB(i,ir,1,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,ip,1,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged right-top corner state (RT->LL)
    !qRT(i,ir,2,j,k) = r + (+drx+drz)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j,k+1)*rho_t(i+1,j,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j,k+1)*  p_t(i+1,j,k+1))) - p
    dry = qRT(i,ir,2,j,k) - r
    dpy = qRT(i,ip,2,j,k) - p
    qRT(i,ir,2,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,ip,2,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,ir,2,j,k) = r + (+drx-drz)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j,k-1)*rho_t(i+1,j,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j,k-1)*  p_t(i+1,j,k-1))) - p
    dry = qRB(i,ir,2,j,k) - r
    dpy = qRB(i,ip,2,j,k) - p
    qRB(i,ir,2,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,ip,2,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged left-top corner state (LT->RL)
    !qLT(i,ir,2,j,k) = r + (-drx+drz)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j,k+1)*rho_t(i-1,j,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j,k+1)*  p_t(i-1,j,k+1))) - p
    dry = qLT(i,ir,2,j,k) - r
    dpy = qLT(i,ip,2,j,k) - p
    qLT(i,ir,2,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,ip,2,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,ir,2,j,k) = r + (-drx-drz)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j,k-1)*rho_t(i-1,j,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j,k-1)*  p_t(i-1,j,k-1))) - p
    dry = qLB(i,ir,2,j,k) - r
    dpy = qLB(i,ip,2,j,k) - p
    qLB(i,ir,2,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,ip,2,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged right-top corner state (RT->LL)
    !qRT(i,ir,3,j,k) = r + (+drx+dry)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j+1,k)*rho_t(i+1,j+1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j+1,k)*  p_t(i+1,j+1,k))) - p
    dry = qRT(i,ir,3,j,k) - r
    dpy = qRT(i,ip,3,j,k) - p
    qRT(i,ir,3,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,ip,3,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,ir,3,j,k) = r + (+drx-dry)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j-1,k)*rho_t(i+1,j-1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j-1,k)*  p_t(i+1,j-1,k))) - p
    dry = qRB(i,ir,3,j,k) - r
    dpy = qRB(i,ip,3,j,k) - p
    qRB(i,ir,3,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,ip,3,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged left-top corner state (LT->RL)
    !qLT(i,ir,3,j,k) = r + (-drx+dry)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j+1,k)*rho_t(i-1,j+1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j+1,k)*  p_t(i-1,j+1,k))) - p
    dry = qLT(i,ir,3,j,k) - r
    dpy = qLT(i,ip,3,j,k) - p
    qLT(i,ir,3,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,ip,3,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,ir,3,j,k) = r + (-drx-dry)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j-1,k)*rho_t(i-1,j-1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j-1,k)*  p_t(i-1,j-1,k))) - p
    dry = qLB(i,ir,3,j,k) - r
    dpy = qLB(i,ip,3,j,k) - p
    qLB(i,ir,3,j,k) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,ip,3,j,k) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
  end do
  end do
  end do
  !call trace%end (itimer(3))
  deallocate (rho_t, p_t)
  call trace%end (itimer(1), detailed_timer=detailed_timer)
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

!===============================================================================
!> Compute fluxes
!===============================================================================
SUBROUTINE cmpflxm (self, nx, ln, lt1, lt2, bn, bt1, bt2, fl)
  class(mhd_t):: self
  real(dp), dimension(:,:,:,:):: fl
  integer:: nx, ln, lt1, lt2, bn, bt1, bt2
  ! local variables
  integer:: i, j, k, lb(3), ub(3), xdim, l(3), idim
  real(dp), dimension(:,:), allocatable:: qleft, qright, fgdnv
  real(dp):: dtds(3)
  integer, save:: itimer=0
  integer:: n, nvar, i0, i1 
  associate (qp=>self%qp, qm=>self%qm)
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%cmpflxm', itimer=itimer, detailed_timer=detailed_timer)
  lb = self%mesh%li
  ub = self%mesh%uo
  nvar = self%nv
  xdim = ln-1
  if (verbose > 2) write (stdlog,'(a,2(2x,3i4))') &
   'mhd_t%cmpflxm: lb, ub=', lb, ub
  !
  ! Avoid using spurious values from allocated, but non-defined entries
  !
  if (xdim == 1) then
    n = ub(1)-lb(1)+1
    allocate (qleft (n,nvar))
    allocate (qright(n,nvar))
    allocate (fgdnv (n,nvar))
    do k = lb(3),ub(3)
      do j = lb(2),ub(2)
        ! Enforce continuity for normal magnetic field
        i0 = lb(1)-2
        i1 = lb(1)-1
        !dir$ simd
        !dir$ ivdep
        do i=1,n
          ! Left state
          qleft (i,1) = qm(i0+i,1  ,xdim,j,k) ! Mass density
          qleft (i,2) = qm(i0+i,5  ,xdim,j,k) ! Pressure
          qleft (i,3) = qm(i0+i,ln ,xdim,j,k) ! Normal velocity
          qleft (i,4) = qm(i0+i,bn ,xdim,j,k) ! Normal magnetic field
          qleft (i,5) = qm(i0+i,lt1,xdim,j,k) ! Tangential velocity 1
          qleft (i,6) = qm(i0+i,bt1,xdim,j,k) ! Tangential magnetic field 1
          qleft (i,7) = qm(i0+i,lt2,xdim,j,k) ! Tangential velocity 2
          qleft (i,8) = qm(i0+i,bt2,xdim,j,k) ! Tangential magnetic field 2
          ! Right state
          qright(i,1) = qp(i1+i,1  ,xdim,j,k) ! Mass density
          qright(i,2) = qp(i1+i,5  ,xdim,j,k) ! Pressure
          qright(i,3) = qp(i1+i,ln ,xdim,j,k) ! Normal velocity
          qright(i,4) = qp(i1+i,bn ,xdim,j,k) ! Normal magnetic field
          qright(i,5) = qp(i1+i,lt1,xdim,j,k) ! Tangential velocity 1
          qright(i,6) = qp(i1+i,bt1,xdim,j,k) ! Tangential magnetic field 1
          qright(i,7) = qp(i1+i,lt2,xdim,j,k) ! Tangential velocity 2
          qright(i,8) = qp(i1+i,bt2,xdim,j,k) ! Tangential magnetic field 2
        end do
        call solver (self,qleft,qright,fgdnv,n,nvar)
        ! Output fluxes
        !dir$ ivdep
        !dir$ simd
        do i=1,n
          fl(i1+i,j,k,1  ) = fgdnv(i,1)  ! Mass density
          fl(i1+i,j,k,5  ) = fgdnv(i,2)  ! Total energy
          fl(i1+i,j,k,ln ) = fgdnv(i,3)  ! Normal momentum
          fl(i1+i,j,k,bn ) = fgdnv(i,4)  ! Normal magnetic field
          fl(i1+i,j,k,lt1) = fgdnv(i,5)  ! Transverse momentum 1
          fl(i1+i,j,k,bt1) = fgdnv(i,6)  ! Transverse magnetic field 1
          fl(i1+i,j,k,lt2) = fgdnv(i,7)  ! Transverse momentum 2
          fl(i1+i,j,k,bt2) = fgdnv(i,8)  ! Transverse magnetic field 2
        end do
      end do
    end do
    deallocate (qleft, qright, fgdnv)
  else if (xdim == 2) then
    n = ub(1)-lb(1)+1
    allocate (qleft (n,nvar))
    allocate (qright(n,nvar))
    allocate (fgdnv (n,nvar))
    do k = lb(3),ub(3)
      do j = lb(2),ub(2)
        ! Enforce continuity for normal magnetic field
        i1 = lb(1)-1
        !dir$ ivdep
        !dir$ simd
        do i=1,n
          ! Left state
          qleft (i,1) = qm(i+i1,1  ,xdim,j-1,k) ! Mass density
          qleft (i,2) = qm(i+i1,5  ,xdim,j-1,k) ! Pressure
          qleft (i,3) = qm(i+i1,ln ,xdim,j-1,k) ! Normal velocity
          qleft (i,4) = qm(i+i1,bn ,xdim,j-1,k) ! Normal magnetic field
          qleft (i,5) = qm(i+i1,lt1,xdim,j-1,k) ! Tangential velocity 1
          qleft (i,6) = qm(i+i1,bt1,xdim,j-1,k) ! Tangential magnetic field 1
          qleft (i,7) = qm(i+i1,lt2,xdim,j-1,k) ! Tangential velocity 2
          qleft (i,8) = qm(i+i1,bt2,xdim,j-1,k) ! Tangential magnetic field 2
          ! Right state
          qright(i,1) = qp(i+i1,1  ,xdim,j  ,k) ! Mass density
          qright(i,2) = qp(i+i1,5  ,xdim,j  ,k) ! Pressure
          qright(i,3) = qp(i+i1,ln ,xdim,j  ,k) ! Normal velocity
          qright(i,4) = qp(i+i1,bn ,xdim,j  ,k) ! Normal magnetic field
          qright(i,5) = qp(i+i1,lt1,xdim,j  ,k) ! Tangential velocity 1
          qright(i,6) = qp(i+i1,bt1,xdim,j  ,k) ! Tangential magnetic field 1
          qright(i,7) = qp(i+i1,lt2,xdim,j  ,k) ! Tangential velocity 2
          qright(i,8) = qp(i+i1,bt2,xdim,j  ,k) ! Tangential magnetic field 2
        end do
        call solver (self,qleft,qright,fgdnv,n,nvar)
        !dir$ ivdep
        !dir$ simd
        do i=1,n
          ! Output fluxes
          fl(i+i1,j,k,1  ) = fgdnv(i,1)  ! Mass density
          fl(i+i1,j,k,5  ) = fgdnv(i,2)  ! Total energy
          fl(i+i1,j,k,ln ) = fgdnv(i,3)  ! Normal momentum
          fl(i+i1,j,k,bn ) = fgdnv(i,4)  ! Normal magnetic field
          fl(i+i1,j,k,lt1) = fgdnv(i,5)  ! Transverse momentum 1
          fl(i+i1,j,k,bt1) = fgdnv(i,6)  ! Transverse magnetic field 1
          fl(i+i1,j,k,lt2) = fgdnv(i,7)  ! Transverse momentum 2
          fl(i+i1,j,k,bt2) = fgdnv(i,8)  ! Transverse magnetic field 2
        end do
      end do
    end do
    deallocate (qleft, qright, fgdnv)
  else if (xdim == 3) then
    n = ub(1)-lb(1)+1
    allocate (qleft (n,nvar))
    allocate (qright(n,nvar))
    allocate (fgdnv (n,nvar))
    do k = lb(3),ub(3)
      do j = lb(2),ub(2)
        ! Enforce continuity for normal magnetic field
        i1 = lb(1)-1
        !dir$ ivdep
        !dir$ simd
        do i=1,n
          ! Left state
          qleft (i,1) = qm(i+i1,1  ,xdim,j,k-1) ! Mass density
          qleft (i,2) = qm(i+i1,5  ,xdim,j,k-1) ! Pressure
          qleft (i,3) = qm(i+i1,ln ,xdim,j,k-1) ! Normal velocity
          qleft (i,4) = qm(i+i1,bn ,xdim,j,k-1) ! Normal magnetic field
          qleft (i,5) = qm(i+i1,lt1,xdim,j,k-1) ! Tangential velocity 1
          qleft (i,6) = qm(i+i1,bt1,xdim,j,k-1) ! Tangential magnetic field 1
          qleft (i,7) = qm(i+i1,lt2,xdim,j,k-1) ! Tangential velocity 2
          qleft (i,8) = qm(i+i1,bt2,xdim,j,k-1) ! Tangential magnetic field 2
          ! Right state
          qright(i,1) = qp(i+i1,1  ,xdim,j,k  ) ! Mass density
          qright(i,2) = qp(i+i1,5  ,xdim,j,k  ) ! Pressure
          qright(i,3) = qp(i+i1,ln ,xdim,j,k  ) ! Normal velocity
          qright(i,4) = qp(i+i1,bn ,xdim,j,k  ) ! Normal magnetic field
          qright(i,5) = qp(i+i1,lt1,xdim,j,k  ) ! Tangential velocity 1
          qright(i,6) = qp(i+i1,bt1,xdim,j,k  ) ! Tangential magnetic field 1
          qright(i,7) = qp(i+i1,lt2,xdim,j,k  ) ! Tangential velocity 2
          qright(i,8) = qp(i+i1,bt2,xdim,j,k  ) ! Tangential magnetic field 2
        end do
        call solver (self,qleft,qright,fgdnv,n,nvar)
        !dir$ ivdep
        !dir$ simd
        do i=1,n
          ! Output fluxes
          fl(i+i1,j,k,1  ) = fgdnv(i,1)  ! Mass density
          fl(i+i1,j,k,5  ) = fgdnv(i,2)  ! Total energy
          fl(i+i1,j,k,ln ) = fgdnv(i,3)  ! Normal momentum
          fl(i+i1,j,k,bn ) = fgdnv(i,4)  ! Normal magnetic field
          fl(i+i1,j,k,lt1) = fgdnv(i,5)  ! Transverse momentum 1
          fl(i+i1,j,k,bt1) = fgdnv(i,6)  ! Transverse magnetic field 1
          fl(i+i1,j,k,lt2) = fgdnv(i,7)  ! Transverse momentum 2
          fl(i+i1,j,k,bt2) = fgdnv(i,8)  ! Transverse magnetic field 2
        end do
      end do 
    end do
    deallocate (qleft, qright, fgdnv)
  endif
  !
  end associate
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE cmpflxm

!===============================================================================
!> Interface to 2D solver for EMF
!===============================================================================
SUBROUTINE cmp_mag_flx (self, qRT, qRB, qLT, qLB, lp1 ,lp2 ,lor ,bp1 , bp2, bor, emf, ngrid, nvar)
  class(mhd_t):: self
  ! 2D Riemann solver to compute EMF at cell edges
  ! indices of the 2 planar velocity lp1 and lp2 and the orthogonal one,
  ! lor and idem for the magnetic field
  integer:: lp1, lp2, lor, bp1, bp2, bor, ngrid, nvar
  integer:: iRT, jRT, kRT, iRB, jRB, kRB, iLT, jLT, kLT, iLB, jLB, kLB
  real(dp),DIMENSION(:,:,:):: emf
  real(dp),DIMENSION(:,:,:,:,:):: qRT, qLT, qRB, qLB
  ! local variables
  integer:: i, j, k, n, xdim, l, u, lo(3), uo(3)
  real(dp),DIMENSION(1:ngrid,nvar)::qLL,qRL,qLR,qRR
  real(dp):: smallp
  integer, save:: itimer=0
  !associate (qRT=>self%qRT,qLT=>self%qLT,qRB=>self%qRB,qLB=>self%qLB)
  call trace%begin ('mhd_t%cmp_mag_flx', itimer=itimer, detailed_timer=detailed_timer)
  !-----------------------------------------------------------------------------
  xdim = lor - 1
!  if (xdim == 3 .or. xdim == 1) then
!    associate (qRT=>self%qRT,qLT=>self%qLT,qRB=>self%qRB,qLB=>self%qLB)
!  else if (xdim == 2) then
!    associate (qRT=>self%qRT,qLT=>self%qRB,qRB=>self%qLT,qLB=>self%qLB) 
!  endif 
  !
  smallp = smallr*smallc**2
  !
  qLL = 0d0
  qRL = 0d0
  qLR = 0d0
  qRR = 0d0
  !
  if (xdim == 3) then
    iRT = 1; jRT = 1; kRT = 0
    iRB = 1; jRB = 0; kRB = 0
    iLT = 0; jLT = 1; kLT = 0
    iLB = 0; jLB = 0; kLB = 0
  else if (xdim == 2) then
    iRT = 1; jRT = 0; kRT = 1
    iRB = 0; jRB = 0; kRB = 1
    iLT = 1; jLT = 0; kLT = 0
    iLB = 0; jLB = 0; kLB = 0
  else if (xdim == 1) then
    iRT = 0; jRT = 1; kRT = 1
    iRB = 0; jRB = 1; kRB = 0
    iLT = 0; jLT = 0; kLT = 1
    iLB = 0; jLB = 0; kLB = 0
  endif 
  !-----------------------------------------------------------------------------
  ! This loop produces the EMFs, which are accessed at {i,j,k}+1, so it needs to
  ! use the interval (li,uo).  It accesses input arrays at {i,j,k}-1, so the
  ! input data is needed in the interval (lo,uo).
  !-----------------------------------------------------------------------------
  lo = self%mesh%li
  uo = self%mesh%uo
  l = lo(1)
  u = uo(1)
  do k = lo(3),uo(3)
    do j = lo(2),uo(2)
     !dir$ ivdep
     !dir$ simd
     do i = l,u
      qLL (i,1) = max(qRT(i-iRT,  1, xdim, j-jRT, k-kRT), smallr)
      qRL (i,1) = max(qLT(i-iLT,  1, xdim, j-jLT, k-kLT), smallr)
      qLR (i,1) = max(qRB(i-iRB,  1, xdim, j-jRB, k-kRB), smallr)
      qRR (i,1) = max(qLB(i-iLB,  1, xdim, j-jLB, k-kLB), smallr)
      qLL (i,2) = max(qRT(i-iRT,  5, xdim, j-jRT, k-kRT), smallp)
      qRL (i,2) = max(qLT(i-iLT,  5, xdim, j-jLT, k-kLT), smallp)
      qLR (i,2) = max(qRB(i-iRB,  5, xdim, j-jRB, k-kRB), smallp)
      qRR (i,2) = max(qLB(i-iLB,  5, xdim, j-jLB, k-kLB), smallp)
      qLL (i,3) = qRT(i-iRT,  lp1, xdim, j-jRT, k-kRT)
      qRL (i,3) = qLT(i-iLT,  lp1, xdim, j-jLT, k-kLT)
      qLR (i,3) = qRB(i-iRB,  lp1, xdim, j-jRB, k-kRB)
      qRR (i,3) = qLB(i-iLB,  lp1, xdim, j-jLB, k-kLB)
      qLL (i,4) = qRT(i-iRT,  lp2, xdim, j-jRT, k-kRT)
      qRL (i,4) = qLT(i-iLT,  lp2, xdim, j-jLT, k-kLT)
      qLR (i,4) = qRB(i-iRB,  lp2, xdim, j-jRB, k-kRB)
      qRR (i,4) = qLB(i-iLB,  lp2, xdim, j-jLB, k-kLB)
      qLL (i,6) = half*(qRT(i-iRT,  bp1, xdim, j-jRT, k-kRT) + qLT(i-iLT, bp1, xdim, j-jLT, k-kLT))
      qRL (i,6) = half*(qRT(i-iRT,  bp1, xdim, j-jRT, k-kRT) + qLT(i-iLT, bp1, xdim, j-jLT, k-kLT))
      qLR (i,6) = half*(qRB(i-iRB,  bp1, xdim, j-jRB, k-kRB) + qLB(i-iLB, bp1, xdim, j-jLB, k-kLB))
      qRR (i,6) = half*(qRB(i-iRB,  bp1, xdim, j-jRB, k-kRB) + qLB(i-iLB, bp1, xdim, j-jLB, k-kLB))
      qLL (i,7) = half*(qRT(i-iRT,  bp2, xdim, j-jRT, k-kRT) + qRB(i-iRB, bp2, xdim, j-jRB, k-kRB))
      qRL (i,7) = half*(qLT(i-iLT,  bp2, xdim, j-jLT, k-kLT) + qLB(i-iLB, bp2, xdim, j-jLB, k-kLB))
      qLR (i,7) = half*(qRT(i-iRT,  bp2, xdim, j-jRT, k-kRT) + qRB(i-iRB, bp2, xdim, j-jRB, k-kRB))
      qRR (i,7) = half*(qLT(i-iLT,  bp2, xdim, j-jLT, k-kLT) + qLB(i-iLB, bp2, xdim, j-jLB, k-kLB))
      qLL (i,5) = qRT(i-iRT,  lor, xdim, j-jRT, k-kRT)
      qRL (i,5) = qLT(i-iLT,  lor, xdim, j-jLT, k-kLT)
      qLR (i,5) = qRB(i-iRB,  lor, xdim, j-jRB, k-kRB)
      qRR (i,5) = qLB(i-iLB,  lor, xdim, j-jLB, k-kLB)
      qLL (i,8) = qRT(i-iRT,  bor, xdim, j-jRT, k-kRT)
      qRL (i,8) = qLT(i-iLT,  bor, xdim, j-jLT, k-kLT)
      qLR (i,8) = qRB(i-iRB,  bor, xdim, j-jRB, k-kRB)
      qRR (i,8) = qLB(i-iLB,  bor, xdim, j-jLB, k-kLB)
     end do
     call solver_2d (self,qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer=detailed_timer)
    end do
  end do
  call trace%end (itimer, detailed_timer=detailed_timer)
!  end associate
END SUBROUTINE cmp_mag_flx
!===============================================================================
!> Selection of 1-D solvers
!===============================================================================
SUBROUTINE solver (self, qleft, qright, fgdnv, n, nvar)
  class(mhd_t):: self
  real(dp),dimension(:,:):: qleft,qright,fgdnv
  integer :: n, nvar, n1
  !.............................................................................
  n1 = size(qleft,1)
  if (self%riemann == 'hlld') then
    call rieman%hlld (qleft, qright, fgdnv, n1, detailed_timer)
  else if (self%riemann == 'hlld_vec') then
    call rieman%hlld_vec (qleft, qright, fgdnv, n1, nvar, detailed_timer)
  else if (self%riemann == 'llf')then
    call rieman%llf (qleft, qright, fgdnv, n1, detailed_timer)
  else
    write(*,*)'unknown Riemann solver'
    stop
  end if
END SUBROUTINE solver

!===============================================================================
!> Selection of 2-D solvers
!===============================================================================
SUBROUTINE solver_2d (self,qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer)
  class(mhd_t):: self
  real(DP),dimension(:,:),intent(in):: qLL,qRL,qLR,qRR
  real(DP),dimension(:,:,:):: emf
  integer, intent(in):: j,k,l,u,nvar,xdim
  logical, optional:: detailed_timer
  !.............................................................................
  if (self%riemann == 'hlld') then
    call rieman%hlld_2d (qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer=detailed_timer)
  else if (self%riemann == 'hlld_vec') then
    call rieman%hlld_2d (qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer=detailed_timer)
  else if (self%riemann == 'llf')then
    call rieman%llf_2d (qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer=detailed_timer)
  else
    write(*,*)'unknown Riemann 2D solver'
    stop
  end if
END SUBROUTINE solver_2d

!===============================================================================
FUNCTION up(f,i) RESULT (g)
  real, dimension(:,:,:), intent(in):: f
  real, dimension(size(f,1),size(f,2),size(f,3)):: g
  integer:: i
  !.............................................................................
  g = 0.5*(cshift(f,1,i)+f)
END FUNCTION

!===============================================================================
!> Add gradient of phi gravitational force.  FIXME: This should be in extras
!===============================================================================
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
  phi => self%mem(:,:,:,self%idx%phi,self%new,1)
  d   => self%mem(:,:,:,self%idx%d,self%new,1)
  et  => self%mem(:,:,:,self%idx%e,self%new,1)
  p   => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
  allocate (gradphi1(self%gn(1),self%gn(2),self%gn(3)), &
            gradphi2(self%gn(1),self%gn(2),self%gn(3)), &
            gradphi3(self%gn(1),self%gn(2),self%gn(3)))
  !
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
