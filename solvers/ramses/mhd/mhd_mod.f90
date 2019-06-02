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
  USE const
  USE io_mod
  USE io_unit_mod
  USE extras_mod
  USE link_mod
  USE riemann_mod
  use hydro_parameters
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE index_mod
  USE non_ideal_mod
  USE omp_timer_mod
  implicit none
  PRIVATE
  type, public, extends(extras_t):: mhd_t
    logical:: mhd=.true.
    integer:: idim=3
    character(len=10):: riemann='hlld'
    real, pointer, dimension(:,:,:,:):: q         ! primitive
    real, pointer, dimension(:,:,:,:):: bf, em    ! face centered B, emf
    real, pointer, dimension(:,:,:,:,:):: qp, qm  ! face centered
    real, pointer, dimension(:,:,:,:,:):: dq, dq2, dbf ! slopes
    real, pointer, dimension(:,:,:,:,:):: qLB,qRB,qLT,qRT
    real, pointer, dimension(:,:,:,:):: uin       ! conserved
    real, pointer, dimension(:,:,:,:):: unew      ! new conserved
    real, pointer, dimension(:,:,:):: Ex,Ey,Ez
    real, pointer, dimension(:,:,:):: c
  contains
    procedure:: init
    procedure:: update
    procedure:: gas_pressure
  end type
  logical, save:: first_time=.true.
  logical, save:: detailed_timer=.false.
  logical, save:: unsigned=.true.
  integer, save:: verbose=0
  integer, save:: divb_cleaner=2
  integer, save:: b_slope_type=0
  integer, save:: u_slope_type=0
  character(len=16):: eqn='mhd'
CONTAINS

!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  character(len=10):: riemann='hlld'
  real, save:: max_dlogd=20.
  integer, save:: n_dump=0
  real:: csound=1.
  integer:: i, iv
  namelist /ramses_params/ gamma, riemann, slope_type, b_slope_type, &
    u_slope_type, courant_factor, smallr, smallc, isothermal, csound, &
    verbose, detailed_timer, unsigned, courant_type, max_dlogd, n_dump, &
    divb_cleaner
  character(len=120):: id = &
   '$Id$ solvers/ramses/mhd/mhd_mod.f90'
!...............................................................................
  call trace%begin('mhd_mod%init')
  call trace%print_id (id)
  self%nw = 1
  self%ng = 3
  self%nv = 8
  self%kind = 'ramses_mhd_patch'
  call self%idx%init (5, self%mhd)
  call self%gpatch_t%init
  !-----------------------------------------------------------------------------
  ! Read solver-relevant parameters and store in the task instance
  !-----------------------------------------------------------------------------
  !$omp critical (pde)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, ramses_params)
    isothermal = isothermal .or. gamma==1d0
    !---------------------------------------------------------------------------
    ! The HLLD solver throws a divide fault if gamma==1, hence this workaround:
    !---------------------------------------------------------------------------
    if (isothermal) gamma=1.0d0
    if (mpi%master) write (*, ramses_params)
  end if
  !$omp end critical (pde)
  rieman%max_dlogd = max_dlogd
  self%n_dump = n_dump
  self%riemann = riemann
  self%courant = courant_factor
  self%gamma = gamma
  self%csound = csound
  self%staggered = .false.
  do i=1,3
    self%mesh(i)%h = 0d0
    self%mesh(i)%h(self%idx%bx+i-1) = -0.5
  end do
  call non_ideal%init  
  self%unsigned(self%idx%d) = unsigned
  self%unsigned(self%idx%e) = unsigned
  self%pervolume(self%idx%px) = unsigned
  self%pervolume(self%idx%py) = unsigned
  self%pervolume(self%idx%pz) = unsigned
  call trace%end
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  real, pointer, dimension(:,:,:,:):: p
  real, allocatable, dimension(:,:,:,:):: em
  real, allocatable, dimension(:,:,:,:,:):: fl
  real, allocatable, dimension(:,:,:,:)::fluxambdiff,emfambdiff
  real, pointer, dimension(:,:,:):: d, e, Ux=>null(), Uy=>null(), Uz=>null()
  real(dp):: dtdx, dtdy, dtdz
  integer:: n(3), nv3, i
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
  self%uin  => self%mem(:,:,:,:,self%it ,1)
  self%unew => self%mem(:,:,:,:,self%new,1)

  n = self%mesh%gn
  allocate (self%c  (n(1),n(2),n(3)))                           ! self%c   = 0_dp
  allocate (self%Ex (n(1),n(2),n(3)))                           ! self%Ex  = 0_dp
  allocate (self%Ey (n(1),n(2),n(3)))                           ! self%Ey  = 0_dp
  allocate (self%Ez (n(1),n(2),n(3)))                           ! self%Ez  = 0_dp
  allocate (self%q  (n(1),n(2),n(3),self%nv))                   ! self%q   = 0_dp
  allocate (self%qp (n(1),n(2),n(3),self%nv,ndim))              ! self%qp  = 0_dp
  allocate (self%qm (n(1),n(2),n(3),self%nv,ndim))              ! self%qm  = 0_dp
  allocate (self%dq (n(1),n(2),n(3),self%nv,ndim))              ! self%dq  = 0_dp
  allocate (self%dq2(n(1),n(2),n(3),self%nv,ndim))              ! self%dq2 = 0_dp
  allocate (self%bf (n(1),n(2),n(3),3))                         ! self%bf  = 0_dp
  allocate (self%dbf(n(1),n(2),n(3),3,2))                       ! self%dbf = 0_dp
  allocate (self%qLB(n(1),n(2),n(3),self%nv,ndim))              ! self%qLB = 0_dp
  allocate (self%qLT(n(1),n(2),n(3),self%nv,ndim))              ! self%qLT = 0_dp
  allocate (self%qRB(n(1),n(2),n(3),self%nv,ndim))              ! self%qRB = 0_dp
  allocate (self%qRT(n(1),n(2),n(3),self%nv,ndim))              ! self%qRT = 0_dp
  d => self%mem(:,:,:,          self%idx%d ,self%it,1)
  e => self%mem(:,:,:,          self%idx%e ,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  ! Compute update using second-order Godunov method
  !-----------------------------------------------------------------------------
  call check_small (self)                               ! Check smallr, smallp
  call ctoprim(self)                                    ! Primitive variables
  !self%dq = 0d0 ; self%dq2 = 0d0 ; self%dbf = 0d0
  call uslope(self,self%bf,self%q,self%dq,self%dbf,slope_type)   ! TVD slopes
  if (verbose>0) print *, 'uslope', self%id, minval(self%dbf), maxval(self%dbf)
  !-----------------------------------------------------------------------------
  ! For slope types in the interval 3-4, use the conservative 3.0 for corners
  ! (cf. RAMSES/CPH)
  !-----------------------------------------------------------------------------
  if (slope_type <= 3.0) then
    self%dq2=self%dq
  else
    call uslope(self,self%bf,self%q,self%dq2,self%dbf,3.0)       ! TVD slopes
  endif

  if (verbose>0) print *, 'uslope', self%id, minval(self%dq), maxval(self%dq)
  if (verbose>0) print *, 'uslope', self%id, minval(self%dq2), maxval(self%dq2)

  if (non_ideal%is_used) then
    allocate (fluxambdiff(self%gn(1),self%gn(2),self%gn(3),3))
    allocate (emfambdiff(self%gn(1),self%gn(2),self%gn(3),3))
    call non_ideal%non_ideal_comp(self%mesh,self%idx,self%bf,self%uin,self%q,&
      fluxambdiff,emfambdiff,self%u_max,self%ds,self%dtime,self%gamma)
  end if
  call self%courant_condition                           ! compute dt

  call trace3d(self)                                    ! predictor step

  allocate (fl(self%gn(1),self%gn(2),self%gn(3),self%nv,3))
  allocate (em(self%gn(1),self%gn(2),self%gn(3),3))
  fl = 0.0d0
  em = 0.0d0
  call cmpflxm(self,self%gn(1),2,3,4,6,7,8,fl(:,:,:,:,1))
  call cmpflxm(self,self%gn(2),3,2,4,7,6,8,fl(:,:,:,:,2))
  call cmpflxm(self,self%gn(3),4,2,3,8,6,7,fl(:,:,:,:,3))

  call cmp_mag_flx(self, self%qRT, self%qRB, self%qLT, self%qLB,2,3,4,6,7,8,em(:,:,:,3),self%gn(3)-1,self%nv)
  call cmp_mag_flx(self, self%qRT, self%qLT, self%qRB, self%qLB,4,2,3,8,6,7,em(:,:,:,2),self%gn(2)-1,self%nv)
  call cmp_mag_flx(self, self%qRT, self%qRB, self%qLT, self%qLB,3,4,2,7,8,6,em(:,:,:,1),self%gn(1)-1,self%nv)

  dtdx = self%dtime/self%ds(1) 
  dtdy = self%dtime/self%ds(2) 
  dtdz = self%dtime/self%ds(3) 
  em(:,:,:,1) = em(:,:,:,1)*dtdy
  em(:,:,:,2) = em(:,:,:,2)*dtdz
  em(:,:,:,3) = em(:,:,:,3)*dtdx

  if (non_ideal%is_used) then
    call non_ideal%flux_upd(fl,fluxambdiff,self%ds(1),self%dtime,1)
    call non_ideal%flux_upd(fl,fluxambdiff,self%ds(2),self%dtime,2)
    call non_ideal%flux_upd(fl,fluxambdiff,self%ds(3),self%dtime,3)
    call non_ideal%non_ideal_emf_up(em(:,:,:,3),emfambdiff,self%ds(3),self%dtime,3)
    call non_ideal%non_ideal_emf_up(em(:,:,:,2),emfambdiff,self%ds(2),self%dtime,2)
    call non_ideal%non_ideal_emf_up(em(:,:,:,1),emfambdiff,self%ds(1),self%dtime,1)
  end if
  ! Compute fluxes with appropriate time step
  fl(:,:,:,:,1) = fl(:,:,:,:,1) * dtdx
  fl(:,:,:,:,2) = fl(:,:,:,:,2) * dtdy
  fl(:,:,:,:,3) = fl(:,:,:,:,3) * dtdz

  !-----------------------------------------------------------------------------
  ! Update conservative variables
  !-----------------------------------------------------------------------------
  associate (i0=>self%mesh(1)%li, j0=>self%mesh(2)%li, k0=>self%mesh(3)%li, &
             i1=>self%mesh(1)%ui, j1=>self%mesh(2)%ui, k1=>self%mesh(3)%ui, &
             uin=>self%uin,unew=>self%unew)

  call self%lock%set ('update')
  unew = uin
  nv3 = 5
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3  ) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,1) - &
    fl(i0+1:i1+1,j0  :j1  ,k0  :k1  ,1:nv3,1) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,2) - &
    fl(i0  :i1  ,j0+1:j1+1,k0  :k1  ,1:nv3,2) + &
    fl(i0  :i1  ,j0  :j1  ,k0  :k1  ,1:nv3,3) - &
    fl(i0  :i1  ,j0  :j1  ,k0+1:k1+1,1:nv3,3)   

  !!!! TO DO !!!! 
  !!! Instead of applying the update in time explicitly in mhd_mod
  !!! only compute the fluxes and provide them to call timestep%update
  !!! analogous to what is done in stagger2/mhd_mod.f90  

  ! Compute emf fluxes with appropriate time steps


  !-----------------------------------------------------------------------------
  ! Update Bx 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,6) =  &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,6) + &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,2)   - &
      em(i0  :i1  ,j0  :j1  ,k0+1:k1+1,2) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,3)   - &
      em(i0  :i1  ,j0+1:j1+1,k0  :k1  ,3) )

  !-----------------------------------------------------------------------------
  ! Update By 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,7) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,7) + &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,3)   - &
      em(i0+1:i1+1,j0  :j1  ,k0  :k1  ,3) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,1)   - &
      em(i0  :i1  ,j0  :j1  ,k0+1:k1+1,1) )

  !-----------------------------------------------------------------------------
  ! Update Bz 
  !-----------------------------------------------------------------------------
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,8) = &
  unew(i0  :i1  ,j0  :j1  ,k0  :k1  ,8) + &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,1)   - &
      em(i0  :i1  ,j0+1:j1+1,k0  :k1  ,1) ) - &
    ( em(i0  :i1  ,j0  :j1  ,k0  :k1  ,2)   - &
      em(i0+1:i1+1,j0  :j1  ,k0  :k1  ,2) )

  if (verbose>0) print *, 'unew', self%id, minval(unew), maxval(unew)
  end associate
  !-----------------------------------------------------------------------------
  ! Reset the energy variable in isothermal cases
  !-----------------------------------------------------------------------------
  if (isothermal) &
    self%mem(:,:,:,self%idx%e,self%new,1) = &
    self%mem(:,:,:,self%idx%d,self%new,1)
  !-----------------------------------------------------------------------------
  ! Add external force as source term
  !-----------------------------------------------------------------------------
  if (allocated(self%force_per_unit_mass)) then
    p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%new,1)
    do i=1,3
      p(:,:,:,i) = p(:,:,:,i) + self%dtime*self%force_per_unit_mass(:,:,:,i)*d
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate
  !-----------------------------------------------------------------------------
  call self%counter_update
  if (non_ideal%is_used) &
    deallocate (fluxambdiff,emfambdiff)
  deallocate (fl,em,self%q,self%qp,self%qm,self%dq,self%dq2,self%bf,self%dbf, &
              self%qLT,self%qLB,self%qRT,self%qRB,self%c, &
              self%Ex,self%Ey,self%Ez)
  !----------------------------------------------------------------------------
  ! Use one of several methods to clean div(B) from guard zones
  !----------------------------------------------------------------------------
  select case (divb_cleaner)
  case (1)
  call divb_clean1 (self)
  case (2)
  call divb_clean2 (self)
  end select
  !.............................................................................
  call self%lock%unset ('update')
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
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
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
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
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
  if (n_smalld+n_smallp > 0) then
    write (io_unit%log,*) 'small: id,nr,np', &
      wallclock(),self%id,self%time,n_smalld,n_smallp
  end if
  end associate
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
  integer:: i, j, k, loc(3), l(3), u(3)
  associate(uin=>self%uin,q=>self%q,bf=>self%bf,dt=>self%dtime,c=>self%c)
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%ctoprim', itimer=itimer, detailed_timer=detailed_timer)
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
  if (verbose>0) print *,  'ctoprim:rh', self%id, minval(q(:,:,:,1)), maxval(q(:,:,:,1))
  if (verbose>0) print *,  'ctoprim:bf', self%id, minval(bf), maxval(bf)

  do k=1,self%gn(3)-1
  do j=1,self%gn(2)-1
  do i=1,self%gn(1)-1
    eken = half*(q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4)**2)*q(i,j,k,1)
    B1 = half*(q(i,j,k,6)+q(i+1,j,k,6))
    B2 = half*(q(i,j,k,7)+q(i,j+1,k,7))
    B3 = half*(q(i,j,k,8)+q(i,j,k+1,8))
    emag = half*(B1**2 + B2**2 + B3**2)
    q(i,j,k,5) = (uin(i,j,k,5) - emag - eken)*(gamma-1d0)
    q(i,j,k,5) = merge(q(i,j,k,1)*self%csound**2, max(q(i,j,k,5),smallp), isothermal)
    c(i,j,k) = sqrt((q(i,j,k,5)*gamma + 2.0*emag)/q(i,j,k,1))
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Alternative Courant condition formulations
  !-----------------------------------------------------------------------------
  if (courant_type==1) then
    c(:,:,:) = c(:,:,:) + max(abs(q(:,:,:,2)),  abs(q(:,:,:,3)),  abs(q(:,:,:,4)))
  else
    c(:,:,:) = c(:,:,:) +    (abs(q(:,:,:,2)) + abs(q(:,:,:,3)) + abs(q(:,:,:,4)))/3.0
  end if
  ! Gravity predictor step
  !q(:,:,:,2:4) = q(:,:,:,idim+1) + Grho*gravin(:,:,:,1:3)*dt*half
  l = self%mesh%lo
  u = self%mesh%uo
  self%u_max = maxval(c(l(1):u(1),l(2):u(2),l(3):u(3)))
  if (self%get_umax_location) then
    loc = maxloc(c(l(1):u(1),l(2):u(2),l(3):u(3))) + l-1
    write (io_unit%log,'(a,i6,f12.6,2(1x,3i4),3f9.4,2(1x,a,1p,4g10.2))') 'u_max:', &
      self%id, self%time, self%ipos, loc, &
      (self%mesh(1)%p + self%mesh(1)%r(loc(1))-self%mesh(2)%llc_cart)/self%mesh(1)%b, &
      (self%mesh(2)%p + self%mesh(2)%r(loc(2))-self%mesh(2)%llc_cart)/self%mesh(2)%b, &
      (self%mesh(3)%p + self%mesh(3)%r(loc(3))-self%mesh(2)%llc_cart)/self%mesh(3)%b, &
      'P,B:', q(loc(1),loc(2),loc(3),5:8), &
      'D,U:', q(loc(1),loc(2),loc(3),1:4)
    flush (io_unit%log)
  end if
  end associate
  call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine ctoprim

!===============================================================================
!> Slope limiters
!===============================================================================
subroutine uslope(self,bf,q,dq,dbf,rslope_type)
  class(mhd_t):: self
  real, dimension(:,:,:,:):: q
  real, dimension(:,:,:,:):: bf
  real, dimension(:,:,:,:,:):: dq
  real, dimension(:,:,:,:,:):: dbf
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
     return
  else if (islope_type==1) then  ! minmod
    do n = 1, nv
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
    xslope_type= 2 !MIN(slope_type,2)
    do n = 1, nv
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            ! slopes in first coordinate direction
            dlft = (q(i  ,j,k,n) - q(i-1,j,k,n))          ! = dq
            drgt = (q(i+1,j,k,n) - q(i  ,j,k,n))          ! = dq
            dcen = half*(dlft+drgt)                       ! = dq
            dsgn = sign(one, dcen)                        ! sign(slope)
            slop = xslope_type*min(abs(dlft),abs(drgt))   ! slope_type=1: slope=dq, slope_type=2: slope=2*dq 
            slop = merge(slop,zero,(dlft*drgt)>zero)      ! if both are non-zero, use the smaller one, or times 2
            dq(i,j,k,n,1) = dsgn*min(slop,abs(dcen))      ! the smaller of the twice the smaller slope, and the central one            ! slopes in second coordinate direction
            dlft = (q(i,j  ,k,n) - q(i,j-1,k,n))
            drgt = (q(i,j+1,k,n) - q(i,j  ,k,n))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dq(i,j,k,n,2) = dsgn*min(slop,abs(dcen))
            ! slopes in third coordinate direction
            dlft = (q(i,j,k  ,n) - q(i,j,k-1,n))
            drgt = (q(i,j,k+1,n) - q(i,j,k  ,n))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dq(i,j,k,n,3) = dsgn*min(slop,abs(dcen))
          end do
        end do
      end do
    end do

    do k = l(3),u(3)
      do j = l(2),u(2)
        do i = l(1),u(1)
          ! Bx
          dlft = (bf(i,j  ,k,1) - bf(i,j-1,k,1))
          drgt = (bf(i,j+1,k,1) - bf(i,j  ,k,1))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,1,1) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j,k  ,1) - bf(i,j,k-1,1))
          drgt = (bf(i,j,k+1,1) - bf(i,j,k  ,1))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,1,2) = dsgn*min(slop,abs(dcen))
          ! By
          ! slopes in first coordinate direction
          dlft = (bf(i  ,j,k,2) - bf(i-1,j,k,2))
          drgt = (bf(i+1,j,k,2) - bf(i  ,j,k,2))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,2,1) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j,k  ,2) - bf(i,j,k-1,2))
          drgt = (bf(i,j,k+1,2) - bf(i,j,k  ,2))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,2,2) = dsgn*min(slop,abs(dcen))
          ! Bz
          ! slopes in first coordinate direction
          dlft = (bf(i  ,j,k,3) - bf(i-1,j,k,3))
          drgt = (bf(i+1,j,k,3) - bf(i  ,j,k,3))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,3,1) = dsgn*min(slop,abs(dcen))
          ! slopes in second coordinate direction
          dlft = (bf(i,j  ,k,3) - bf(i,j-1,k,3))
          drgt = (bf(i,j+1,k,3) - bf(i,j  ,k,3))
          dcen = half*(dlft+drgt)
          dsgn = sign(one, dcen)
          slop = xslope_type*min(abs(dlft),abs(drgt))
          slop = merge(slop,zero,(dlft*drgt)>zero)
          dbf(i,j,k,3,2) = dsgn*min(slop,abs(dcen))
        end do
      end do
    end do

  else if (islope_type==3) then ! positivity preserving 3d unsplit slope
    xslope_type = real(slope_type - 2_dp,kind=dp)
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
            dq(i,j,k,n,1) = dlim*dfx
            dq(i,j,k,n,2) = dlim*dfy
            dq(i,j,k,n,3) = dlim*dfz
          end do
        end do
      end do
    end do
      do k = l(3),u(3)
        do j = l(2),u(2)
          do i = l(1),u(1)
            ! Bx
            dlft = (bf(i,j  ,k,1) - bf(i,j-1,k,1))
            drgt = (bf(i,j+1,k,1) - bf(i,j  ,k,1))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,1,1) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j,k  ,1) - bf(i,j,k-1,1))
            drgt = (bf(i,j,k+1,1) - bf(i,j,k  ,1))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,1,2) = dsgn*min(slop,abs(dcen))
            !By
            ! slopes in first coordinate direction
            dlft = (bf(i  ,j,k,2) - bf(i-1,j,k,2))
            drgt = (bf(i+1,j,k,2) - bf(i  ,j,k,2))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,2,1) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j,k  ,2) - bf(i,j,k-1,2))
            drgt = (bf(i,j,k+1,2) - bf(i,j,k  ,2))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,2,2) = dsgn*min(slop,abs(dcen))
            ! Bz
            ! slopes in first coordinate direction
            dlft = (bf(i  ,j,k,3) - bf(i-1,j,k,3))
            drgt = (bf(i+1,j,k,3) - bf(i  ,j,k,3))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,3,1) = dsgn*min(slop,abs(dcen))
            ! slopes in second coordinate direction
            dlft = (bf(i,j  ,k,3) - bf(i,j-1,k,3))
            drgt = (bf(i,j+1,k,3) - bf(i,j  ,k,3))
            dcen = half*(dlft+drgt)
            dsgn = sign(one, dcen)
            slop = xslope_type*min(abs(dlft),abs(drgt))
            slop = merge(slop,zero,(dlft*drgt)>zero)
            dbf(i,j,k,3,2) = dsgn*min(slop,abs(dcen))
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
  real(dp), dimension(1:self%gn(1),1:self%gn(2),1:self%gn(3)):: rho_t, p_t
  real(dp) :: eff_gamma, off
  integer, save:: itimer(4)=0
  real, allocatable:: ff(:,:,:,:)
  !-----------------------------------------------------------------------------
  associate(q=>self%q,bf=>self%bf,dq=>self%dq,dq2=>self%dq2,qm=>self%qm,qp=>self%qp, &
    dbf=>self%dbf,qRT=>self%qRT,qRB=>self%qRB,qLT=>self%qLT,qLB=>self%qLB, &
    Bx=>self%bf(:,:,:,1),By=>self%bf(:,:,:,2),Bz=>self%bf(:,:,:,3), &
    Ux=>self%q(:,:,:,2),Uy=>self%q(:,:,:,3),Uz=>self%q(:,:,:,4), &
    Ex=>self%Ex,Ey=>self%Ey,Ez=>self%Ez)
  !
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

  ! fill positive definite variables with reasonable values in "ghost-zones"
  rho_t(:,:,:) = q(:,:,:,ir)
  p_t(:,:,:)   = q(:,:,:,ip)

  !
  !call trace%end (itimer(1))
  !call trace%begin ('mhd_t%trace3d(2)', itimer=itimer(2))
  eff_gamma = gamma
  do k=2,self%gn(3)-1
  do j=2,self%gn(2)-1
  !!dir! simd
  do i=2,self%gn(1)-1
!    associate ( &
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
    drx = half * dq(i,j,k,ir,1)
    drx2= half *dq2(i,j,k,ir,1)
    dux = half * dq(i,j,k,iu,1)
    dvx = half * dq(i,j,k,iv,1)
    dwx = half * dq(i,j,k,iw,1)
    dpx = half * dq(i,j,k,ip,1)
    dpx2= half *dq2(i,j,k,ip,1)
    dBx = half * dq(i,j,k,iB,1)
    dCx = half * dq(i,j,k,iC,1)
    !
    dry = half * dq(i,j,k,ir,2)
    dry2= half *dq2(i,j,k,ir,2)
    duy = half * dq(i,j,k,iu,2)
    dvy = half * dq(i,j,k,iv,2)
    dwy = half * dq(i,j,k,iw,2)
    dpy = half * dq(i,j,k,ip,2)
    dpy2= half *dq2(i,j,k,ip,2)
    duy = half * dq(i,j,k,iu,2)
    dAy = half * dq(i,j,k,iA,2)
    dCy = half * dq(i,j,k,iC,2)
    !
    drz = half * dq(i,j,k,ir,3)
    drz2= half *dq2(i,j,k,ir,3)
    duz = half * dq(i,j,k,iu,3)
    dvz = half * dq(i,j,k,iv,3)
    dwz = half * dq(i,j,k,iw,3)
    dpz = half * dq(i,j,k,ip,3)
    dpz2= half *dq2(i,j,k,ip,3)
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
    sp0 = (-u*dpx-dux*eff_gamma*p)*dtdx + (-v*dpy-dvy*eff_gamma*p)*dtdy + (-w*dpz-dwz*eff_gamma*p)*dtdz
    !
    ! Add external force
    associate (ff => self%force_per_unit_mass)
    su0 = su0 + ff(i,j,k,1)*half*self%dtime
    sv0 = sv0 + ff(i,j,k,2)*half*self%dtime
    sw0 = sw0 + ff(i,j,k,3)*half*self%dtime
    end associate
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
!          if (verbose>0) then
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
    qRT(i,j,k,ir,1) = r + (+dry2+drz2)
    qRT(i,j,k,iu,1) = u + (+duy+duz)
    qRT(i,j,k,iv,1) = v + (+dvy+dvz)
    qRT(i,j,k,iw,1) = w + (+dwy+dwz)
    qRT(i,j,k,ip,1) = p + (+dpy2+dpz2)
    qRT(i,j,k,iA,1) = A + (+dAy+dAz)
    qRT(i,j,k,iB,1) = BR+ (   +dBRz)
    qRT(i,j,k,iC,1) = CR+ (+dCRy   )
    !
    ! X-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,1) = r + (+dry2-drz2)
    qRB(i,j,k,iu,1) = u + (+duy-duz)
    qRB(i,j,k,iv,1) = v + (+dvy-dvz)
    qRB(i,j,k,iw,1) = w + (+dwy-dwz)
    qRB(i,j,k,ip,1) = p + (+dpy2-dpz2)
    qRB(i,j,k,iA,1) = A + (+dAy-dAz)
    qRB(i,j,k,iB,1) = BR+ (   -dBRz)
    qRB(i,j,k,iC,1) = CL+ (+dCLy   )
    !
    ! X-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,1) = r + (-dry2+drz2)
    qLT(i,j,k,iu,1) = u + (-duy+duz)
    qLT(i,j,k,iv,1) = v + (-dvy+dvz)
    qLT(i,j,k,iw,1) = w + (-dwy+dwz)
    qLT(i,j,k,ip,1) = p + (-dpy2+dpz2)
    qLT(i,j,k,iA,1) = A + (-dAy+dAz)
    qLT(i,j,k,iB,1) = BL+ (   +dBLz)
    qLT(i,j,k,iC,1) = CR+ (-dCRy   )
    !
    ! X-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,1) = r + (-dry2-drz2)
    qLB(i,j,k,iu,1) = u + (-duy-duz)
    qLB(i,j,k,iv,1) = v + (-dvy-dvz)
    qLB(i,j,k,iw,1) = w + (-dwy-dwz)
    qLB(i,j,k,ip,1) = p + (-dpy2-dpz2)
    qLB(i,j,k,iA,1) = A + (-dAy-dAz)
    qLB(i,j,k,iB,1) = BL+ (   -dBLz)
    qLB(i,j,k,iC,1) = CL+ (-dCLy   )
    !
    ! Y-edge averaged right-top corner state (RT->LL)
    qRT(i,j,k,ir,2) = r + (+drx2+drz2)
    qRT(i,j,k,iu,2) = u + (+dux+duz)
    qRT(i,j,k,iv,2) = v + (+dvx+dvz)
    qRT(i,j,k,iw,2) = w + (+dwx+dwz)
    qRT(i,j,k,ip,2) = p + (+dpx2+dpz2)
    qRT(i,j,k,iA,2) = AR+ (   +dARz)
    qRT(i,j,k,iB,2) = B + (+dBx+dBz)
    qRT(i,j,k,iC,2) = CR+ (+dCRx   )
    !
    ! Y-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,2) = r + (+drx2-drz2)
    qRB(i,j,k,iu,2) = u + (+dux-duz)
    qRB(i,j,k,iv,2) = v + (+dvx-dvz)
    qRB(i,j,k,iw,2) = w + (+dwx-dwz)
    qRB(i,j,k,ip,2) = p + (+dpx2-dpz2)
    qRB(i,j,k,iA,2) = AR+ (   -dARz)
    qRB(i,j,k,iB,2) = B + (+dBx-dBz)
    qRB(i,j,k,iC,2) = CL+ (+dCLx   )
    !
    ! Y-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,2) = r + (-drx2+drz2)
    qLT(i,j,k,iu,2) = u + (-dux+duz)
    qLT(i,j,k,iv,2) = v + (-dvx+dvz)
    qLT(i,j,k,iw,2) = w + (-dwx+dwz)
    qLT(i,j,k,ip,2) = p + (-dpx2+dpz2)
    qLT(i,j,k,iA,2) = AL+ (   +dALz)
    qLT(i,j,k,iB,2) = B + (-dBx+dBz)
    qLT(i,j,k,iC,2) = CR+ (-dCRx   )
    !
    ! Y-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,2) = r + (-drx2-drz2)
    qLB(i,j,k,iu,2) = u + (-dux-duz)
    qLB(i,j,k,iv,2) = v + (-dvx-dvz)
    qLB(i,j,k,iw,2) = w + (-dwx-dwz)
    qLB(i,j,k,ip,2) = p + (-dpx2-dpz2)
    qLB(i,j,k,iA,2) = AL+ (   -dALz)
    qLB(i,j,k,iB,2) = B + (-dBx-dBz)
    qLB(i,j,k,iC,2) = CL+ (-dCLx   )
    !
    ! Z-edge averaged right-top corner state (RT->LL)
    qRT(i,j,k,ir,3) = r + (+drx2+dry2)
    qRT(i,j,k,iu,3) = u + (+dux+duy)
    qRT(i,j,k,iv,3) = v + (+dvx+dvy)
    qRT(i,j,k,iw,3) = w + (+dwx+dwy)
    qRT(i,j,k,ip,3) = p + (+dpx2+dpy2)
    qRT(i,j,k,iA,3) = AR+ (   +dARy)
    qRT(i,j,k,iB,3) = BR+ (+dBRx   )
    qRT(i,j,k,iC,3) = C + (+dCx+dCy)
    !
    ! Z-edge averaged right-bottom corner state (RB->LR)
    qRB(i,j,k,ir,3) = r + (+drx2-dry2)
    qRB(i,j,k,iu,3) = u + (+dux-duy)
    qRB(i,j,k,iv,3) = v + (+dvx-dvy)
    qRB(i,j,k,iw,3) = w + (+dwx-dwy)
    qRB(i,j,k,ip,3) = p + (+dpx2-dpy2)
    qRB(i,j,k,iA,3) = AR+ (   -dARy)
    qRB(i,j,k,iB,3) = BL+ (+dBLx   )
    qRB(i,j,k,iC,3) = C + (+dCx-dCy)
    !
    ! Z-edge averaged left-top corner state (LT->RL)
    qLT(i,j,k,ir,3) = r + (-drx2+dry2)
    qLT(i,j,k,iu,3) = u + (-dux+duy)
    qLT(i,j,k,iv,3) = v + (-dvx+dvy)
    qLT(i,j,k,iw,3) = w + (-dwx+dwy)
    qLT(i,j,k,ip,3) = p + (-dpx2+dpy2)
    qLT(i,j,k,iA,3) = AL+ (   +dALy)
    qLT(i,j,k,iB,3) = BR+ (-dBRx   )
    qLT(i,j,k,iC,3) = C + (-dCx+dCy)
    !
    ! Z-edge averaged left-bottom corner state (LB->RR)
    qLB(i,j,k,ir,3) = r + (-drx2-dry2)
    qLB(i,j,k,iu,3) = u + (-dux-duy)
    qLB(i,j,k,iv,3) = v + (-dvx-dvy)
    qLB(i,j,k,iw,3) = w + (-dwx-dwy)
    qLB(i,j,k,ip,3) = p + (-dpx2-dpy2)
    qLB(i,j,k,iA,3) = AL+ (   -dALy)
    qLB(i,j,k,iB,3) = BL+ (-dBLx   )
    qLB(i,j,k,iC,3) = C + (-dCx-dCy)
  end do
  end do
  end do
    !
  !call trace%end (itimer(2))
  !call trace%begin ('mhd_t%trace3d(3)', itimer=itimer(3))
  do k=2,self%gn(3)-1
  do j=2,self%gn(2)-1
  do i=2,self%gn(1)-1
    r = rho_t(i,j,k)
    p = p_t(i,j,k)

    ! X-edge averaged right-top corner state (RT->LL)
    !qRT(i,j,k,ir,1) = r + (+dry+drz)
    drx = sqrt(sqrt(r*rho_t(i,j+1,k)*rho_t(i,j,k+1)*rho_t(i,j+1,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j+1,k)*  p_t(i,j,k+1)*  p_t(i,j+1,k+1))) - p
    dry = qRT(i,j,k,ir,1) - r
    dpy = qRT(i,j,k,ip,1) - p
    qRT(i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,j,k,ir,1) = r + (+dry-drz)
    drx = sqrt(sqrt(r*rho_t(i,j+1,k)*rho_t(i,j,k-1)*rho_t(i,j+1,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j+1,k)*  p_t(i,j,k-1)*  p_t(i,j+1,k-1))) - p
    dry = qRB(i,j,k,ir,1) - r
    dpy = qRB(i,j,k,ip,1) - p
    qRB(i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged left-top corner state (LT->RL)
    !qLT(i,j,k,ir,1) = r + (-dry+drz)
    drx = sqrt(sqrt(r*rho_t(i,j-1,k)*rho_t(i,j,k+1)*rho_t(i,j-1,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j-1,k)*  p_t(i,j,k+1)*  p_t(i,j-1,k+1))) - p
    dry = qLT(i,j,k,ir,1) - r
    dpy = qLT(i,j,k,ip,1) - p
    qLT(i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! X-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,j,k,ir,1) = r + (-dry-drz)
    drx = sqrt(sqrt(r*rho_t(i,j-1,k)*rho_t(i,j,k-1)*rho_t(i,j-1,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i,j-1,k)*  p_t(i,j,k-1)*  p_t(i,j-1,k-1))) - p
    dry = qLB(i,j,k,ir,1) - r
    dpy = qLB(i,j,k,ip,1) - p
    qLB(i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged right-top corner state (RT->LL)
    !qRT(i,j,k,ir,2) = r + (+drx+drz)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j,k+1)*rho_t(i+1,j,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j,k+1)*  p_t(i+1,j,k+1))) - p
    dry = qRT(i,j,k,ir,2) - r
    dpy = qRT(i,j,k,ip,2) - p
    qRT(i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,j,k,ir,2) = r + (+drx-drz)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j,k-1)*rho_t(i+1,j,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j,k-1)*  p_t(i+1,j,k-1))) - p
    dry = qRB(i,j,k,ir,2) - r
    dpy = qRB(i,j,k,ip,2) - p
    qRB(i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

    ! Y-edge averaged left-top corner state (LT->RL)
    !qLT(i,j,k,ir,2) = r + (-drx+drz)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j,k+1)*rho_t(i-1,j,k+1))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j,k+1)*  p_t(i-1,j,k+1))) - p
    dry = qLT(i,j,k,ir,2) - r
    dpy = qLT(i,j,k,ip,2) - p
    qLT(i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Y-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,j,k,ir,2) = r + (-drx-drz)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j,k-1)*rho_t(i-1,j,k-1))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j,k-1)*  p_t(i-1,j,k-1))) - p
    dry = qLB(i,j,k,ir,2) - r
    dpy = qLB(i,j,k,ip,2) - p
    qLB(i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged right-top corner state (RT->LL)
    !qRT(i,j,k,ir,3) = r + (+drx+dry)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j+1,k)*rho_t(i+1,j+1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j+1,k)*  p_t(i+1,j+1,k))) - p
    dry = qRT(i,j,k,ir,3) - r
    dpy = qRT(i,j,k,ip,3) - p
    qRT(i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRT(i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged right-bottom corner state (RB->LR)
    !qRB(i,j,k,ir,3) = r + (+drx-dry)
    drx = sqrt(sqrt(r*rho_t(i+1,j,k)*rho_t(i,j-1,k)*rho_t(i+1,j-1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i+1,j,k)*  p_t(i,j-1,k)*  p_t(i+1,j-1,k))) - p
    dry = qRB(i,j,k,ir,3) - r
    dpy = qRB(i,j,k,ip,3) - p
    qRB(i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qRB(i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged left-top corner state (LT->RL)
    !qLT(i,j,k,ir,3) = r + (-drx+dry)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j+1,k)*rho_t(i-1,j+1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j+1,k)*  p_t(i-1,j+1,k))) - p
    dry = qLT(i,j,k,ir,3) - r
    dpy = qLT(i,j,k,ip,3) - p
    qLT(i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLT(i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p
    ! Z-edge averaged left-bottom corner state (LB->RR)
    !qLB(i,j,k,ir,3) = r + (-drx-dry)
    drx = sqrt(sqrt(r*rho_t(i-1,j,k)*rho_t(i,j-1,k)*rho_t(i-1,j-1,k))) - r
    dpx = sqrt(sqrt(p*  p_t(i-1,j,k)*  p_t(i,j-1,k)*  p_t(i-1,j-1,k))) - p
    dry = qLB(i,j,k,ir,3) - r
    dpy = qLB(i,j,k,ip,3) - p
    qLB(i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
    qLB(i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

!    end associate
  end do
  end do
  end do
  !call trace%end (itimer(3))
  call trace%end (itimer(1), detailed_timer=detailed_timer)
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

!===============================================================================
!> Compute fluxes
!===============================================================================
subroutine cmpflxm(self, nx, ln, lt1, lt2, bn, bt1, bt2, fl)
  class(mhd_t):: self
  real(dp), dimension(:,:,:,:):: fl
  integer:: nx, ln, lt1, lt2, bn, bt1, bt2
  ! local variables
  integer:: i, j, k, lb(3), ub(3), xdim, l(3), idim
  real(dp), dimension(nx,self%nv):: qleft, qright, fgdnv
  real(dp), dimension(nx):: bn_mean
  real(dp):: dtds(3)
  integer, save:: itimer=0
  integer:: n
  associate (qp=>self%qp, qm=>self%qm)
  !-----------------------------------------------------------------------------
  call trace%begin ('mhd_t%cmpflxm', itimer=itimer, detailed_timer=detailed_timer)
  lb = self%mesh%lb+2
  ub = self%mesh%ub-1
  xdim = ln-1
  if (verbose>0) print '(a,2(2x,3i4))', 'mhd_t%cmpflxm: lb, ub=', lb, ub
  !
  ! Avoid using spurious values from allocated, but non-defined entries
  bn_mean=0d0

  if (xdim == 1) then
    do k = lb(3),ub(3)
      do j = lb(2),ub(2)
        ! Enforce continuity for normal magnetic field
        n = ub(1)-lb(1)+1
        bn_mean(1:n) = half*(qm(lb(1)-1:ub(1)-1,j,k,bn,xdim)+qp(lb(1):ub(1),j,k,bn,xdim))
        ! Left state
        qleft (1:n,1) = qm(lb(1)-1:ub(1)-1,j,k,1  ,xdim) ! Mass density
        qleft (1:n,2) = qm(lb(1)-1:ub(1)-1,j,k,5  ,xdim) ! Pressure
        qleft (1:n,3) = qm(lb(1)-1:ub(1)-1,j,k,ln ,xdim) ! Normal velocity
        qleft (1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qleft (1:n,5) = qm(lb(1)-1:ub(1)-1,j,k,lt1,xdim) ! Tangential velocity 1
        qleft (1:n,6) = qm(lb(1)-1:ub(1)-1,j,k,bt1,xdim) ! Tangential magnetic field 1
        qleft (1:n,7) = qm(lb(1)-1:ub(1)-1,j,k,lt2,xdim) ! Tangential velocity 2
        qleft (1:n,8) = qm(lb(1)-1:ub(1)-1,j,k,bt2,xdim) ! Tangential magnetic field 2
        ! Right state
        qright(1:n,1) = qp(lb(1):ub(1),j,k,1  ,xdim) ! Mass density
        qright(1:n,2) = qp(lb(1):ub(1),j,k,5  ,xdim) ! Pressure
        qright(1:n,3) = qp(lb(1):ub(1),j,k,ln ,xdim) ! Normal velocity
        qright(1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qright(1:n,5) = qp(lb(1):ub(1),j,k,lt1,xdim) ! Tangential velocity 1
        qright(1:n,6) = qp(lb(1):ub(1),j,k,bt1,xdim) ! Tangential magnetic field 1
        qright(1:n,7) = qp(lb(1):ub(1),j,k,lt2,xdim) ! Tangential velocity 2
        qright(1:n,8) = qp(lb(1):ub(1),j,k,bt2,xdim) ! Tangential magnetic field 2
        call solver (self,qleft,qright,fgdnv,n)
        ! Output fluxes
        fl(lb(1):ub(1),j,k,1  ) = fgdnv(1:n,1)  ! Mass density
        fl(lb(1):ub(1),j,k,5  ) = fgdnv(1:n,2)  ! Total energy
        fl(lb(1):ub(1),j,k,ln ) = fgdnv(1:n,3)  ! Normal momentum
        fl(lb(1):ub(1),j,k,bn ) = fgdnv(1:n,4)  ! Normal magnetic field
        fl(lb(1):ub(1),j,k,lt1) = fgdnv(1:n,5)  ! Transverse momentum 1
        fl(lb(1):ub(1),j,k,bt1) = fgdnv(1:n,6)  ! Transverse magnetic field 1
        fl(lb(1):ub(1),j,k,lt2) = fgdnv(1:n,7)  ! Transverse momentum 2
        fl(lb(1):ub(1),j,k,bt2) = fgdnv(1:n,8)  ! Transverse magnetic field 2
      end do
    end do
  else if (xdim == 2) then
    do i = lb(1),ub(1)
      do k = lb(3),ub(3)
        ! Enforce continuity for normal magnetic field
        n = ub(2)-lb(2)+1
        bn_mean(1:n) = half*(qm(i,lb(2)-1:ub(2)-1,k,bn,xdim)+qp(i,lb(2):ub(2),k,bn,xdim))
        ! Left state
        qleft (1:n,1) = qm(i,lb(2)-1:ub(2)-1,k,1  ,xdim) ! Mass density
        qleft (1:n,2) = qm(i,lb(2)-1:ub(2)-1,k,5  ,xdim) ! Pressure
        qleft (1:n,3) = qm(i,lb(2)-1:ub(2)-1,k,ln ,xdim) ! Normal velocity
        qleft (1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qleft (1:n,5) = qm(i,lb(2)-1:ub(2)-1,k,lt1,xdim) ! Tangential velocity 1
        qleft (1:n,6) = qm(i,lb(2)-1:ub(2)-1,k,bt1,xdim) ! Tangential magnetic field 1
        qleft (1:n,7) = qm(i,lb(2)-1:ub(2)-1,k,lt2,xdim) ! Tangential velocity 2
        qleft (1:n,8) = qm(i,lb(2)-1:ub(2)-1,k,bt2,xdim) ! Tangential magnetic field 2
        ! Right state
        qright(1:n,1) = qp(i,lb(2):ub(2),k,1  ,xdim) ! Mass density
        qright(1:n,2) = qp(i,lb(2):ub(2),k,5  ,xdim) ! Pressure
        qright(1:n,3) = qp(i,lb(2):ub(2),k,ln ,xdim) ! Normal velocity
        qright(1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qright(1:n,5) = qp(i,lb(2):ub(2),k,lt1,xdim) ! Tangential velocity 1
        qright(1:n,6) = qp(i,lb(2):ub(2),k,bt1,xdim) ! Tangential magnetic field 1
        qright(1:n,7) = qp(i,lb(2):ub(2),k,lt2,xdim) ! Tangential velocity 2
        qright(1:n,8) = qp(i,lb(2):ub(2),k,bt2,xdim) ! Tangential magnetic field 2
        call solver (self,qleft,qright,fgdnv,n)
        ! Output fluxes
        fl(i,lb(2):ub(2),k,1  ) = fgdnv(1:n,1)  ! Mass density
        fl(i,lb(2):ub(2),k,5  ) = fgdnv(1:n,2)  ! Total energy
        fl(i,lb(2):ub(2),k,ln ) = fgdnv(1:n,3)  ! Normal momentum
        fl(i,lb(2):ub(2),k,bn ) = fgdnv(1:n,4)  ! Normal magnetic field
        fl(i,lb(2):ub(2),k,lt1) = fgdnv(1:n,5)  ! Transverse momentum 1
        fl(i,lb(2):ub(2),k,bt1) = fgdnv(1:n,6)  ! Transverse magnetic field 1
        fl(i,lb(2):ub(2),k,lt2) = fgdnv(1:n,7)  ! Transverse momentum 2
        fl(i,lb(2):ub(2),k,bt2) = fgdnv(1:n,8)  ! Transverse magnetic field 2
      end do
    end do
  else if (xdim == 3) then
    do j = lb(2),ub(2)
      do i = lb(1),ub(1)
        ! Enforce continuity for normal magnetic field
        n = ub(3)-lb(3)+1
        bn_mean(1:n) = half*(qm(i,j,lb(3)-1:ub(3)-1,bn,xdim)+qp(i,j,lb(3):ub(3),bn,xdim))
        ! Left state
        qleft (1:n,1) = qm(i,j,lb(3)-1:ub(3)-1,1  ,xdim) ! Mass density
        qleft (1:n,2) = qm(i,j,lb(3)-1:ub(3)-1,5  ,xdim) ! Pressure
        qleft (1:n,3) = qm(i,j,lb(3)-1:ub(3)-1,ln ,xdim) ! Normal velocity
        qleft (1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qleft (1:n,5) = qm(i,j,lb(3)-1:ub(3)-1,lt1,xdim) ! Tangential velocity 1
        qleft (1:n,6) = qm(i,j,lb(3)-1:ub(3)-1,bt1,xdim) ! Tangential magnetic field 1
        qleft (1:n,7) = qm(i,j,lb(3)-1:ub(3)-1,lt2,xdim) ! Tangential velocity 2
        qleft (1:n,8) = qm(i,j,lb(3)-1:ub(3)-1,bt2,xdim) ! Tangential magnetic field 2
        ! Right state
        qright(1:n,1) = qp(i,j,lb(3):ub(3),1  ,xdim) ! Mass density
        qright(1:n,2) = qp(i,j,lb(3):ub(3),5  ,xdim) ! Pressure
        qright(1:n,3) = qp(i,j,lb(3):ub(3),ln ,xdim) ! Normal velocity
        qright(1:n,4) = bn_mean(1:n)         ! Normal magnetic field
        qright(1:n,5) = qp(i,j,lb(3):ub(3),lt1,xdim) ! Tangential velocity 1
        qright(1:n,6) = qp(i,j,lb(3):ub(3),bt1,xdim) ! Tangential magnetic field 1
        qright(1:n,7) = qp(i,j,lb(3):ub(3),lt2,xdim) ! Tangential velocity 2
        qright(1:n,8) = qp(i,j,lb(3):ub(3),bt2,xdim) ! Tangential magnetic field 2
        call solver (self,qleft,qright,fgdnv,n)
        ! Output fluxes
        fl(i,j,lb(3):ub(3),1  ) = fgdnv(1:n,1)  ! Mass density
        fl(i,j,lb(3):ub(3),5  ) = fgdnv(1:n,2)  ! Total energy
        fl(i,j,lb(3):ub(3),ln ) = fgdnv(1:n,3)  ! Normal momentum
        fl(i,j,lb(3):ub(3),bn ) = fgdnv(1:n,4)  ! Normal magnetic field
        fl(i,j,lb(3):ub(3),lt1) = fgdnv(1:n,5)  ! Transverse momentum 1
        fl(i,j,lb(3):ub(3),bt1) = fgdnv(1:n,6)  ! Transverse magnetic field 1
        fl(i,j,lb(3):ub(3),lt2) = fgdnv(1:n,7)  ! Transverse momentum 2
        fl(i,j,lb(3):ub(3),bt2) = fgdnv(1:n,8)  ! Transverse magnetic field 2
      end do 
    end do
      
  endif

  end associate
  call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine cmpflxm

!===============================================================================
!> 2D solver for EMF
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
  integer:: i, j, k, n, xdim, l, u
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
  !
  if (xdim == 1) then
    do k = self%mesh(3)%li,self%mesh(3)%uo
      do j = self%mesh(2)%li,self%mesh(2)%uo
       l = self%mesh(1)%li
       u = self%mesh(1)%ui
       do i = l,u
        qLL (i,1) = max(qRT(i-iRT, j-jRT, k-kRT, 1, xdim), smallr)
        qRL (i,1) = max(qLT(i-iLT, j-jLT, k-kLT, 1, xdim), smallr)
        qLR (i,1) = max(qRB(i-iRB, j-jRB, k-kRB, 1, xdim), smallr)
        qRR (i,1) = max(qLB(i-iLB, j-jLB, k-kLB, 1, xdim), smallr)
        qLL (i,2) = max(qRT(i-iRT, j-jRT, k-kRT, 5, xdim), smallp)
        qRL (i,2) = max(qLT(i-iLT, j-jLT, k-kLT, 5, xdim), smallp)
        qLR (i,2) = max(qRB(i-iRB, j-jRB, k-kRB, 5, xdim), smallp)
        qRR (i,2) = max(qLB(i-iLB, j-jLB, k-kLB, 5, xdim), smallp)
        qLL (i,3) = qRT(i-iRT, j-jRT, k-kRT, lp1, xdim)
        qRL (i,3) = qLT(i-iLT, j-jLT, k-kLT, lp1, xdim)
        qLR (i,3) = qRB(i-iRB, j-jRB, k-kRB, lp1, xdim)
        qRR (i,3) = qLB(i-iLB, j-jLB, k-kLB, lp1, xdim)
        qLL (i,4) = qRT(i-iRT, j-jRT, k-kRT, lp2, xdim)
        qRL (i,4) = qLT(i-iLT, j-jLT, k-kLT, lp2, xdim)
        qLR (i,4) = qRB(i-iRB, j-jRB, k-kRB, lp2, xdim)
        qRR (i,4) = qLB(i-iLB, j-jLB, k-kLB, lp2, xdim)
        qLL (i,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qRL (i,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qLR (i,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qRR (i,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qLL (i,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRL (i,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLR (i,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRR (i,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLL (i,5) = qRT(i-iRT, j-jRT, k-kRT, lor, xdim)
        qRL (i,5) = qLT(i-iLT, j-jLT, k-kLT, lor, xdim)
        qLR (i,5) = qRB(i-iRB, j-jRB, k-kRB, lor, xdim)
        qRR (i,5) = qLB(i-iLB, j-jLB, k-kLB, lor, xdim)
        qLL (i,8) = qRT(i-iRT, j-jRT, k-kRT, bor, xdim)
        qRL (i,8) = qLT(i-iLT, j-jLT, k-kLT, bor, xdim)
        qLR (i,8) = qRB(i-iRB, j-jRB, k-kRB, bor, xdim)
        qRR (i,8) = qLB(i-iLB, j-jLB, k-kLB, bor, xdim)
       end do
       call solver_2d (self,qLL,qRL,qLR,qRR,emf,j,k,l,u,nvar,xdim,detailed_timer=detailed_timer)
      end do
    end do
  else if (xdim == 2) then
    do i = self%mesh(1)%li,self%mesh(1)%uo
      do k = self%mesh(3)%li,self%mesh(3)%uo
       l = self%mesh(2)%li
       u = self%mesh(2)%ui
       do j = l,u
        qLL (j,1) = max(qRT(i-iRT, j-jRT, k-kRT, 1, xdim), smallr)
        qRL (j,1) = max(qLT(i-iLT, j-jLT, k-kLT, 1, xdim), smallr)
        qLR (j,1) = max(qRB(i-iRB, j-jRB, k-kRB, 1, xdim), smallr)
        qRR (j,1) = max(qLB(i-iLB, j-jLB, k-kLB, 1, xdim), smallr)
        qLL (j,2) = max(qRT(i-iRT, j-jRT, k-kRT, 5, xdim), smallp)
        qRL (j,2) = max(qLT(i-iLT, j-jLT, k-kLT, 5, xdim), smallp)
        qLR (j,2) = max(qRB(i-iRB, j-jRB, k-kRB, 5, xdim), smallp)
        qRR (j,2) = max(qLB(i-iLB, j-jLB, k-kLB, 5, xdim), smallp)
        qLL (j,3) = qRT(i-iRT, j-jRT, k-kRT, lp1, xdim)
        qRL (j,3) = qLT(i-iLT, j-jLT, k-kLT, lp1, xdim)
        qLR (j,3) = qRB(i-iRB, j-jRB, k-kRB, lp1, xdim)
        qRR (j,3) = qLB(i-iLB, j-jLB, k-kLB, lp1, xdim)
        qLL (j,4) = qRT(i-iRT, j-jRT, k-kRT, lp2, xdim)
        qRL (j,4) = qLT(i-iLT, j-jLT, k-kLT, lp2, xdim)
        qLR (j,4) = qRB(i-iRB, j-jRB, k-kRB, lp2, xdim)
        qRR (j,4) = qLB(i-iLB, j-jLB, k-kLB, lp2, xdim)
        qLL (j,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qRL (j,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qLR (j,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qRR (j,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qLL (j,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRL (j,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLR (j,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRR (j,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLL (j,5) = qRT(i-iRT, j-jRT, k-kRT, lor, xdim)
        qRL (j,5) = qLT(i-iLT, j-jLT, k-kLT, lor, xdim)
        qLR (j,5) = qRB(i-iRB, j-jRB, k-kRB, lor, xdim)
        qRR (j,5) = qLB(i-iLB, j-jLB, k-kLB, lor, xdim)
        qLL (j,8) = qRT(i-iRT, j-jRT, k-kRT, bor, xdim)
        qRL (j,8) = qLT(i-iLT, j-jLT, k-kLT, bor, xdim)
        qLR (j,8) = qRB(i-iRB, j-jRB, k-kRB, bor, xdim)
        qRR (j,8) = qLB(i-iLB, j-jLB, k-kLB, bor, xdim)
       end do
       call solver_2d (self,qLL,qRL,qLR,qRR,emf,k,i,l,u,nvar,xdim,detailed_timer=detailed_timer)
      end do
    end do
  else if (xdim == 3) then    
    do j = self%mesh(2)%li,self%mesh(2)%uo
      do i = self%mesh(1)%li,self%mesh(1)%uo
       l = self%mesh(3)%li
       u = self%mesh(3)%ui
       do k = l,u
        qLL (k,1) = max(qRT(i-iRT, j-jRT, k-kRT, 1, xdim), smallr)
        qRL (k,1) = max(qLT(i-iLT, j-jLT, k-kLT, 1, xdim), smallr)
        qLR (k,1) = max(qRB(i-iRB, j-jRB, k-kRB, 1, xdim), smallr)
        qRR (k,1) = max(qLB(i-iLB, j-jLB, k-kLB, 1, xdim), smallr)
        qLL (k,2) = max(qRT(i-iRT, j-jRT, k-kRT, 5, xdim), smallp)
        qRL (k,2) = max(qLT(i-iLT, j-jLT, k-kLT, 5, xdim), smallp)
        qLR (k,2) = max(qRB(i-iRB, j-jRB, k-kRB, 5, xdim), smallp)
        qRR (k,2) = max(qLB(i-iLB, j-jLB, k-kLB, 5, xdim), smallp)
        qLL (k,3) = qRT(i-iRT, j-jRT, k-kRT, lp1, xdim)
        qRL (k,3) = qLT(i-iLT, j-jLT, k-kLT, lp1, xdim)
        qLR (k,3) = qRB(i-iRB, j-jRB, k-kRB, lp1, xdim)
        qRR (k,3) = qLB(i-iLB, j-jLB, k-kLB, lp1, xdim)
        qLL (k,4) = qRT(i-iRT, j-jRT, k-kRT, lp2, xdim)
        qRL (k,4) = qLT(i-iLT, j-jLT, k-kLT, lp2, xdim)
        qLR (k,4) = qRB(i-iRB, j-jRB, k-kRB, lp2, xdim)
        qRR (k,4) = qLB(i-iLB, j-jLB, k-kLB, lp2, xdim)
        qLL (k,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qRL (k,6) = half*(qRT(i-iRT, j-jRT, k-kRT, bp1, xdim) + qLT(i-iLT, j-jLT, k-kLT, bp1, xdim))
        qLR (k,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qRR (k,6) = half*(qRB(i-iRB, j-jRB, k-kRB, bp1, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp1, xdim))
        qLL (k,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRL (k,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLR (k,7) = half*(qRT(i-iRT, j-jRT, k-kRT, bp2, xdim) + qRB(i-iRB, j-jRB, k-kRB, bp2, xdim))
        qRR (k,7) = half*(qLT(i-iLT, j-jLT, k-kLT, bp2, xdim) + qLB(i-iLB, j-jLB, k-kLB, bp2, xdim))
        qLL (k,5) = qRT(i-iRT, j-jRT, k-kRT, lor, xdim)
        qRL (k,5) = qLT(i-iLT, j-jLT, k-kLT, lor, xdim)
        qLR (k,5) = qRB(i-iRB, j-jRB, k-kRB, lor, xdim)
        qRR (k,5) = qLB(i-iLB, j-jLB, k-kLB, lor, xdim)
        qLL (k,8) = qRT(i-iRT, j-jRT, k-kRT, bor, xdim)
        qRL (k,8) = qLT(i-iLT, j-jLT, k-kLT, bor, xdim)
        qLR (k,8) = qRB(i-iRB, j-jRB, k-kRB, bor, xdim)
        qRR (k,8) = qLB(i-iLB, j-jLB, k-kLB, bor, xdim)
       end do
       call solver_2d (self,qLL,qRL,qLR,qRR,emf,i,j,l,u,nvar,xdim,detailed_timer=detailed_timer)
      end do
    end do
  endif
  call trace%end (itimer, detailed_timer=detailed_timer)
!  end associate
END SUBROUTINE cmp_mag_flx

!===============================================================================
!> Selection of 1-D solvers
!===============================================================================
SUBROUTINE solver (self, qleft, qright, fgdnv, ub)
  class(mhd_t):: self
  real(dp),dimension(:,:):: qleft,qright,fgdnv
  integer :: lb, ub
  !.............................................................................
  if (self%riemann == 'hlld') then
    call rieman%hlld (qleft, qright, fgdnv, ub)
  else if (self%riemann == 'llf')then
    call rieman%llf (qleft, qright, fgdnv, ub)
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
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  integer :: l(3), u(3)
  !-----------------------------------------------------------------------------
  associate (d  => self%mem(:,:,:,self%idx%d ,self%it,1), &
             px => self%mem(:,:,:,self%idx%px,self%it,1), &
             py => self%mem(:,:,:,self%idx%py,self%it,1), &
             pz => self%mem(:,:,:,self%idx%pz,self%it,1), &
             bx => self%mem(:,:,:,self%idx%bx,self%it,1), &
             by => self%mem(:,:,:,self%idx%by,self%it,1), &
             bz => self%mem(:,:,:,self%idx%bz,self%it,1), &
             e  => self%mem(:,:,:,self%idx%e ,self%it,1))
  if (isothermal) then
    pg = d*self%csound**2
  else
    pg = (gamma-1d0)* &
      (e-0.5*(up(bx,1)**2+up(by,2)**2+up(bz,3)**2+(px**2+py**2+pz**2)/d))
  end if
  end associate
END FUNCTION gas_pressure

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
