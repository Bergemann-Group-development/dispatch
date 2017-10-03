!> =============================================================================
!> Staggered mesh solver using entropy per unit volume as the energy variable.
!>
!> Flop counting is enable via annotation at the end of lines. To count flops,
!> use a Linux pipeline such as
!>
!> % grep ops: ../../../solvers/stagger2/mhd_mod.f90 | sed 's/.*ops://' | sumcol 2
!>
!> where 'sumcol' is any script that sums up values in column 2
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE mpi_mod
  USE trace_mod
  USE kinds_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE vector_ops
  USE stagger_mod
  USE gpatch_mod
  USE force_mod
  USE index_mod
  USE timestep_mod
  USE eos_mod
  USE non_ideal_mod
  USE initial_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, s, dsdt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: csound
    real:: nu(6)
    real:: cdtd
    logical:: mhd=.false.
    logical:: first_time=.true.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: pde
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
  end type
  logical, save:: flop_count=.false.
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  type(index_t):: idx
  integer:: i
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  logical, save:: first_time=.true.
  character(len=16), save:: eos
  namelist /mhd_params/ nu, gamma, csound, cdtd, hardwire, eos, flop_count
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    eos = self%eos
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  !$omp end critical (input_cr)
  self%eos = eos
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%cdtd = cdtd
  self%initial%mhd = self%mhd
  call self%initial%init (self%kind, self%gamma)
  self%mhd = self%initial%mhd
  if (self%mhd) then
    self%kind = 'stagger2_mhd_patch'
  else
    self%kind = 'stagger2_hd_patch'
  end if
  if (self%mhd) then
    self%nv = 8
  else
    self%nv = 5
  end if
  !-----------------------------------------------------------------------------
  ! Allocate memory and mesh, etc
  !-----------------------------------------------------------------------------
  call self%gpatch_t%init
  self%i_d  = idx%d
  self%i_e  = idx%e
  self%i_s  = idx%s
  self%i_px = idx%px
  self%i_py = idx%py
  self%i_pz = idx%pz
  self%i_bx = idx%bx
  self%i_by = idx%by
  self%i_bz = idx%bz
  self%type = 'mhd_t'
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    if (self%mhd) self%mesh(i)%h(5+i) = -0.5
  end do
  self%unsigned(idx%d)  = .true.
  self%pervolume(idx%s) = .true.
  call self%force%init (self%kind)
  !$omp critical (input_cr)
  if (first_time) then
    call stagger_test
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call trace%end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, emf, J, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb, Q
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz, phi
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   pointer:: Txx, Txy
  real, dimension(:,:,:),   pointer::      Tyy, Tyz
  real, dimension(:,:,:),   pointer:: Tzx,      Tzz
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true.
  !.............................................................................
  real(8):: ds(3)
  integer, save:: itimer=0
  type(index_t):: i
  real:: u_max
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%pde',itimer=itimer)
  self%d => self%mem(:,:,:,i%d      ,self%it,1); self%dddt => self%mem(:,:,:,i%d      ,self%it,2)
  self%s => self%mem(:,:,:,i%s      ,self%it,1); self%dsdt => self%mem(:,:,:,i%s      ,self%it,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%it,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  d=>self%d; dddt=>self%dddt; s=>self%s; dsdt=>self%dsdt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%it,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  call allocate_vectors (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pa, du, ss, cs, pp, pb, Q)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux=>U(:,:,:,1); Uy=>U(:,:,:,2); Uz=>U(:,:,:,3)
  Ex=>emf(:,:,:,1); Ey=>emf(:,:,:,2); Ez=>emf(:,:,:,3)
  Jx=>J(:,:,:,1); Jy=>J(:,:,:,2); Jz=>J(:,:,:,3)
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2)
                    fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1);                   fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)                                                                  ! ops: 0 0 0 0 3
  ld = down(lnd)
  dd = exp(ld)                                                                  ! ops: 0 0 0 0 3
  U = p/dd                                                                      ! ops: 0 0 3
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
  !-----------------------------------------------------------------------------
  ss = s/d                                                                      ! ops: 0 0 1
  pg = self%gas_pressure()
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  ds = self%ds
  du = -maxval(ds)*div(ds,U)                            ! convergence rate      ! ops: 0 1 0
  du = max(du,0.0)                                      ! positive part
  self%u_max = 0.0
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  ! The diagonal fxx is valid in (2:gn-1)
  !-----------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txx, Tyy, Tzz)
  Txx = ddxup(ds,Ux)
  Tyy = ddyup(ds,Uy)
  Tzz = ddzup(ds,Uz)
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = xup(Bx)**2                                                             ! mhd: 0 1 0
    Ey = yup(By)**2                                                             ! mhd: 0 1 0
    Ez = zup(Bz)**2                                                             ! mhd: 0 1 0
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)                                                         ! mhd: 2 1 0
    ! Fast mode speed + velocity
    pa = self%nu(2)*d*du**2
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)                                    ! mhd: 1 1 0
    if (io%verbose>2) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      if (io%verbose>2) print 1, 'u_max(ca)  ', u_max, self%u_max
    1 format(1x,a,1p,2g12.3)
      u_max = self%fmaxval(cs + sqrt(pp))
      self%u_max = max(self%u_max,u_max)
      if (io%verbose>2) print 1, 'u_max(ca+u)', u_max, self%u_max
    end if
    u_max = self%fmaxval(cs + sqrt(pp) + du)
    self%u_max = max(self%u_max,u_max)
    if (io%verbose>2) print 1, 'u_max(ca+u)', u_max, self%u_max
    ! Isotropic pressure part
    pp = pg + pa + pb                                                           ! mhd: 1 0 0
    ! Total diagonal stress
    du = (maxval(ds)*self%nu(1))*cs
    fxx = pp + d*(fxx - du*Txx) - Ex                                            ! mhd: 1 0 0
    fyy = pp + d*(fyy - du*Tyy) - Ey                                            ! mhd: 1 0 0
    fzz = pp + d*(fzz - du*Tzz) - Ez                                            ! mhd: 1 0 0
  else
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2                                                            ! ops: 0 1 0
    fyy = yup(Uy)**2                                                            ! ops: 0 1 0
    fzz = zup(Uz)**2                                                            ! ops: 0 1 0
    pp = (fxx+fyy+fzz)                                                          ! ops: 2 0 0
    ! Isotropic pressure part
    pa = self%nu(2)*d*du**2                          ! dimension: pressure      ! ops: 0 3 0
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)                                             ! ops: 1 1 1 1 0
    if (io%verbose>2) then
      u_max = self%fmaxval(cs)
      self%u_max = max(self%u_max,u_max)
      if (io%verbose>2) print 1, 'u_max(cs)  ', u_max, self%u_max
      u_max = self%fmaxval(cs + sqrt(pp))
      self%u_max = max(self%u_max,u_max)
      if (io%verbose>2) print 1, 'u_max(cs+u)', u_max, self%u_max
    end if
    u_max = self%fmaxval(cs + sqrt(pp) + du)
    self%u_max = max(self%u_max,u_max)
    if (io%verbose>2) print 1, 'u_max(cs+u)', u_max, self%u_max
    pp = pg + pa                                                                ! ops: 1 0 0
    ! Total diagonal stress, with density factor
    du = (maxval(ds)*self%nu(1))*cs                                             ! ops: 0 1 0
    !---------------------------------------------------------------------------
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = pp + d*(fxx - du*Txx)                                                 ! ops: 2 2 0
    fyy = pp + d*(fyy - du*Tyy)                                                 ! ops: 2 2 0
    fzz = pp + d*(fzz - du*Tzz)                                                 ! ops: 2 2 0
  end if
  !-----------------------------------------------------------------------------
  ! du has units area per unit time, which gives velocity squared times mass
  ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)                                           ! ops: 2 5 0
  call deallocate_scalars (Txx, Tyy, Tzz)
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux   ! ops: 0 9 0
  fm = p+fd                                             ! mass flux             ! ops: 3 0 0
  dddt = -div(ds,fm)                                    ! density derivative
  u_max = self%cdtd*maxval(ds)*self%fmaxval(abs(dddt/d))                        ! ops: 0 0 1
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>2) print 1, 'u_max(dddt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !-----------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)                                                         ! ops: 0 1 0
  fyz = zdn(Uy)*ydn(Uz)                                                         ! ops: 0 1 0
  fzx = xdn(Uz)*zdn(Ux)                                                         ! ops: 0 1 0
  !-----------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !-----------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txy, Tyz, Tzx)
  Txy = ddxdn(ds,Uy)+ddydn(ds,Ux)                                               ! ops: 0 1 0
  Tyz = ddydn(ds,Uz)+ddzdn(ds,Uy)                                               ! ops: 0 1 0
  Tzx = ddzdn(ds,Ux)+ddxdn(ds,Uz)                                               ! ops: 0 1 0
  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy                                 ! ops: 1 3 0
  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz                                 ! ops: 1 3 0
  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx                                 ! ops: 1 3 0
  Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &                             ! ops: 2 5 0
                               yup(zup(Tyz))**2 + &                             ! ops: 1 1 0
                               zup(xup(Tzx))**2)                                ! ops: 1 0 0
  call deallocate_scalars (Txy, Tyz, Tzx)
  !-----------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !-----------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy                                                       ! ops: 0 1 0 0 1
  fyz = exp(zdn(ldy))*fyz                                                       ! ops: 0 1 0 0 1
  fzx = exp(xdn(ldz))*fzx                                                       ! ops: 0 1 0 0 1
  !-----------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - ydn(Bx)*xdn(By)                                                 ! mhd: 1 1 0
    fyz = fyz - zdn(By)*ydn(Bz)                                                 ! mhd: 1 1 0
    fzx = fzx - xdn(Bz)*zdn(Bx)                                                 ! mhd: 1 1 0
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)                                                        ! ops: 0 2 0
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
    fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)                                        ! ops: 1 2 0
    fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)                                        ! ops: 1 2 0
    fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)                                        ! ops: 1 2 0
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))                           ! ops: 2 2 0
    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))                           ! ops: 2 2 0
    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))                           ! ops: 2 2 0
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)                    ! ops: 2 0 0
  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)                    ! ops: 2 0 0
  dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)                    ! ops: 2 0 0
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
   dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)! ops: 1 0 0
  if (self%selfgravity) then
    phi => self%mem(:,:,:,self%i_phi,self%it,1)
    dpdtx = dpdtx - dd(:,:,:,1)*ddxdn(ds,phi)
    dpdty = dpdty - dd(:,:,:,2)*ddydn(ds,phi)
    dpdtz = dpdtz - dd(:,:,:,3)*ddzdn(ds,phi)
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
    fm = fm*down(ss) - (self%nu(5)*dd*down(du))*ddown(ds,ss)                    ! ops: 1 4 0
    !---------------------------------------------------------------------------
    ! Q is heating per unit volume (energy per unit volume and time), which when
    ! divided by presssure (energy per unit volume) becomes entropy (unitless)
    ! per unit mass.  Multiplied by mass density it becomes mass-entropy per
    ! unit volume, which is our units for s.
    !---------------------------------------------------------------------------
    dsdt = -div(ds,fm) + Q*d/pg                                                 ! ops: 1 1 1
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    emf = self%nu(6)*edge(du)*J                                                 ! mhd: 0 6 0
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))                     ! mhd: 2 3 0
    if (self%gamma/=1.0) &
      dsdt = dsdt + Q*d/pg                                                      ! mhd: 1 1 1
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)                                 ! mhd: 2 2 0
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)                                 ! mhd: 2 2 0
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)                                 ! mhd: 2 2 0
    call non_ideal%update (self%gn, minval(ds), d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  ! Time step limitation for rate of change dS/dt = (ds/dt-(s/d)*dddt)/d
  !-----------------------------------------------------------------------------
  u_max = self%cdtd*maxval(ds)*self%fmaxval(abs((dsdt-s*dddt/d)/d))             ! ops: 1 1 2
  self%u_max = max(self%u_max,u_max)
  if (io%verbose>2) print 1, 'u_max(dsdt)', u_max, self%u_max
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ss, cs, pp, pb, Q)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace%end (itimer)
  if (flop_count) then
    print *, stagger%count, ' calls', stagger%flops,' flops'
    stagger%count=0
    stagger%flops=0
  end if
END SUBROUTINE pde

!> =============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod, where boundary conditions may be inserted
!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Prepare the external forces for the pde call, in particular selfgravity. At
  ! this point the nbors are declared up-to-date with respect to everything,
  ! including their potentials, which we need for guard zones.  They may in fact
  ! be providing only their preliminary potentials in their (nt-1) time slots,
  ! while their older time slots should by now contain the converged potential.
  !-----------------------------------------------------------------------------
  call self%gravity
  !-----------------------------------------------------------------------------
  ! Output the state after guard zone downloads, selfgravity, and boundary
  ! conditions. If another call is placed at an earlier point in the task handling
  ! the call here will have no effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  call self%pde
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  !---------------------------------------------------------------------------
  ! Having now advanced the MHD state to the new time, and before allowing
  ! task_list_t%update to rotate the time slots, we make a preliminary
  ! call to the Poisson solver, based on extrapolation of the nbor potentials.
  ! We then provide this preliminary potential to nbors, until we have a
  ! better offer.
  !---------------------------------------------------------------------------
  if (self%selfgravity) then
    call self%dnload (only=self%i_phi)
    call self%gravity (only=self%i_phi)
  end if
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(mhd_t):: self
  character(len=64):: experiment_name, filename
  optional:: experiment_name
  logical:: exists
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%time < self%out_next .or. .not.io%do_output .or. .not.io%do_legacy) return
  call trace%begin('mhd_t%output',itimer=itimer)
  call self%gpatch_t%output (experiment_name)
  !$omp critical (out_cr)
  write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
    self%id, self%iout
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (out_cr)
  call trace%end (itimer)
END SUBROUTINE output

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: d, s
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (trim(self%eos)=='ideal') then
    d => self%mem(:,:,:,self%i_d,self%it,1)
    s => self%mem(:,:,:,self%i_s,self%it,1)
    if (self%gamma==1.0) then
      pg = d*self%csound**2
    else
      s => self%mem(:,:,:,self%i_s,self%it,1)
      pg = d**self%gamma*exp(s/d*(self%gamma-1))
    end if
  else
    call eos%lookup (shape(d), x=d, y=s, pg=pg)
  end if
END FUNCTION gas_pressure

END MODULE mhd_mod
