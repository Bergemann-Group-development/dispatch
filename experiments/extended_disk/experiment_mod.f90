!===============================================================================
!> Experiment initialization and update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE bits_mod
  USE boundaries_mod
  USE mesh_mod
  USE scaling_mod
  USE scratch_mod
  USE eos_mod
  USE units_mod
  USE vector_mod, only: cross
  USE disk_model_mod
  USE disk_particles_mod
  USE task_mod
  USE validate_mod
  USE omp_lock_mod
  USE mpi_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    real:: mass=1.0
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: initial_condition
    procedure:: update
    procedure:: post_update
    procedure:: force
    procedure:: bndries
    procedure:: boundary
    procedure:: update_position
    procedure:: get_position
    procedure:: velocity_shift
    procedure:: write_nml
  end type                                       ! profile
  integer      :: verbose     = 0
  real         :: f_outflow   = 0.01
  logical      :: with_planet = .false.
CONTAINS

!===============================================================================
!> Pre-allocate extras
!===============================================================================
SUBROUTINE pre_init (self)
  class(experiment_t), target:: self
  !-----------------------------------------------------------------------------
  allocate (disk_particles_t:: self%particle_solver)
  call self%solver_t%pre_init ()
  self%gn = self%n + 2*self%ng
END SUBROUTINE pre_init

!===============================================================================
!> Read parameters and set up initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  logical, save:: first_time=.true.
  namelist /experiment_params/ verbose, f_outflow, with_planet
  !.............................................................................
  integer:: l(3), u(3), iostat
  real(8):: r, z, v_k, v_phi
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  !----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, experiment_params)
    end if
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    io%nml_version = 5
  end if
  call self%lock%init
  call io%lock%unset
  !-----------------------------------------------------------------------------
  r = self%radius                                        ! radius in code units
  r = max(r,0.1)  ! FIXME
  z = self%position(3)                                   ! height
  v_k = disk_model%v_kepler (r)                          ! Kepler speed
  self%omega = [0d0, 0d0, v_k/r]                         ! angular velocity
  v_phi = v_k+disk_model%v_head(r,z)                     ! azimuthal speed
  self%velocity = [0d0, v_phi, 0d0]                      ! initial velocity
  self%csound = disk_model%c_sound(r)/sqrt(eos%gamma)    ! isothermal sound
                                                         ! speed
  eos%csound = self%csound                               ! important!
  !----------------------------------------------------------------------------
  if (verbose > 1) then
    l = self%mesh%li
    u = self%mesh%ui
    d  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d ,self%it,1)
    e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,self%it,1)
    px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,self%it,1)
    py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,self%it,1)
    pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,self%it,1)
    print '(i6,1p,4g12.3)', self%id, minval(d), maxval(d), maxval(abs(px/d)), &
      maxval(abs(py/d))
  end if
  call trace%end ()
END SUBROUTINE init

SUBROUTINE post_update (self)
  class(experiment_t):: self
  call trace%begin ('experiment_t%post_update')
  call self%solver_t%post_update ()
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Note, the force should be moved to pre_update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer:: id
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  id = self%id
  !-----------------------------------------------------------------------------
  ! Allocate (first time) and reset the force per unit mass -- extras may add
  ! through the call to solver_t%pre_update()
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%force_per_unit_mass, [self%gn,3])
  scratch%force_per_unit_mass = 0.0
  self%particle_solver%v_patch => self%v_patch
  if (self%id==1) &
    call validate%check (self, self%mem(:,:,:,:,self%it,1), 'upd1')
  !-----------------------------------------------------------------------------
  ! When computing forces already here, we need to do the force reset here also
  !-----------------------------------------------------------------------------
  call self%gpatch_t%pre_update ()
  call self%force ()
  if (self%id==1) &
    call validate%check (self, scratch%force_per_unit_mass(:,:,:,:), 'ff')
  call self%solver_t%update ()
  if (self%id==1) &
    call validate%check (self, self%mem(:,:,:,:,self%it,1), 'upd2')
  !-----------------------------------------------------------------------------
  ! Time slots have been rotated in the solver_t%update()
  !-----------------------------------------------------------------------------
  call self%boundary (self%it)
  if (self%id==1) &
    call validate%check (self, self%mem(:,:,:,:,self%it,1), 'bdr2')
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Disk with power-law mid-plane density and temperature profiles
!===============================================================================
SUBROUTINE initial_condition (self, l, u, it)
  class(experiment_t):: self
  integer            :: l(3), u(3)
  integer, optional  :: it
  !.............................................................................
  integer:: i1, i2, i3
  real(8):: r, z, csound, h, dmid, v_k, v_phi, v1, v2, p(3), g1
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  logical, save :: first_time = .true.
  !-----------------------------------------------------------------------------
  if (present(it)) then
    d  => self%mem(:,:,:,self%idx%d ,it,1)
    e  => self%mem(:,:,:,self%idx%e ,it,1)
    px => self%mem(:,:,:,self%idx%px,it,1)
    py => self%mem(:,:,:,self%idx%py,it,1)
    pz => self%mem(:,:,:,self%idx%pz,it,1)
  else
    d  => self%mem(:,:,:,self%idx%d ,self%it,1)
    e  => self%mem(:,:,:,self%idx%e ,self%it,1)
    px => self%mem(:,:,:,self%idx%px,self%it,1)
    py => self%mem(:,:,:,self%idx%py,self%it,1)
    pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  end if
  v_k = disk_model%v_kepler (self%radius)
  g1 = merge (1d0, eos%gamma-1d0, eos%gamma==1d0)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  !!omp parallel do default(none) private(i1,p,r,z,csound,h,dmid,e,d,v_phi) &
  !!omp private(v1,v2) shared(i2,i3,disk_model,eos,self,px,py,pz,v_k,l,u)
  do i1=l(1),u(1)
    p(1) = self%mesh(1)%r(i1) + self%radius
    p(2) = self%mesh(2)%r(i2)
    p(3) = self%mesh(3)%r(i3) + self%position(3)
    !---------------------------------------------------------------------------
    ! Gas density and thermal energy
    !---------------------------------------------------------------------------
    r = sqrt(p(1)**2 + p(2)**2)                       ! local radius
    z = p(3)                                          ! height
    csound = disk_model%c_sound(r)                    ! sound speed
    h = disk_model%h_scale (r)                        ! scale height
    dmid = disk_model%d_midplane (r)                  ! mid-plane density
    d (i1,i2,i3) = dmid*exp(-(z/h)**2/2d0)            ! Gaussian profile
    e (i1,i2,i3) = d(i1,i2,i3)*csound**2 / &          ! adiabatic
      (eos%gamma*g1)
    !---------------------------------------------------------------------------
    ! Initial head wind, relative to the patch
    !---------------------------------------------------------------------------
    v_phi = disk_model%v_head(r,z)                        ! azimuthal speed
    v1 = v_phi*self%phi_hat(i1,i2,i3,1)                   ! local x1
    v2 = v_phi*self%phi_hat(i1,i2,i3,2)                   ! local x2
    v1 = v1 + self%v_patch(i1,i2,i3,1)
    v2 = v2 + self%v_patch(i1,i2,i3,2) - v_k
    px(i1,i2,i3) = v1*d(i1,i2,i3)                         ! radial momentum
    py(i1,i2,i3) = v2*d(i1,i2,i3)                         ! azimuthal momentum
    pz(i1,i2,i3) = 0.0                                    ! vertical momentum
  end do
  end do
  end do
  if (self%kind(1:7)=='ramses_' .and. eos%gamma > 1.0) then
    e = e + 0.5*(px**2+py**2+pz**2)/d
  end if
  call self%particle_solver%initial_condition
END SUBROUTINE initial_condition

!===============================================================================
!> Boundary condition call from halo_t%update
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  !.............................................................................
  call boundary (self, self%it)
END SUBROUTINE bndries

!===============================================================================
!> Boundary conditions
!===============================================================================
SUBROUTINE boundary (self, it)
  class(experiment_t):: self
  integer :: it
  real, dimension(:,:,:), contiguous, pointer :: d, e, px, py, pz, ui
  integer:: l(3), u(3), i1, i2, i3
  real(8):: p(3), r, z, v_phi
  type(lock_t), pointer:: lock
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%boundary', itimer=itimer)
  lock => self%mem_lock(self%it)
  call lock%set ('experiment_t%boundary')
  !-----------------------------------------------------------------------------
  ! Inner radial boundary
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,it,1)
  px => self%mem(:,:,:,self%idx%px,it,1)
  py => self%mem(:,:,:,self%idx%py,it,1)
  pz => self%mem(:,:,:,self%idx%pz,it,1)
  allocate (ui(size(d,1),size(d,2),size(d,3)))
  allocate (e (size(d,1),size(d,2),size(d,3)))
  if (self%boundaries%is_set (bits%xl)) then
    call self%e2E_th (e,it)
    !---------------------------------------------------------------------------
    ! Apply zero gradient conditions to all variables
    !---------------------------------------------------------------------------
    call self%boundaries%zero_gradient (d, self%idx%d)
    call self%boundaries%zero_gradient (e, self%idx%e)
    call self%boundaries%zero_gradient (px, self%idx%px)
    call self%boundaries%zero_gradient (py, self%idx%py)
    call self%boundaries%zero_gradient (pz, self%idx%pz)
    !---------------------------------------------------------------------------
    ! Subtract headwind velocity at innermost non-guard cell from radial and
    ! azimuthal velocity because of lacking pressure gradient in zero gradient
    ! conditions
    !---------------------------------------------------------------------------
    l = [self%mesh(1)%lb, self%mesh(2)%lb, self%mesh(3)%lb]
    u = [self%mesh(1)%lo, self%mesh(2)%ub, self%mesh(3)%ub]
    call cache_v_head (self)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          p(1) = self%mesh(1)%r(self%mesh(1)%li) + self%radius
          p(2) = self%mesh(2)%r(i2)
          p(3) = self%mesh(3)%r(i3) + self%position(3)
          r = sqrt(p(1)**2 + p(2)**2)
          z = p(3)
          v_phi = self%v_head(i1,i2,i3)
          px(i1,i2,i3) = px(i1,i2,i3) - &
            v_phi*self%phi_hat(i1,i2,i3,1)*d(i1,i2,i3)
          py(i1,i2,i3) = py(i1,i2,i3) - &
            v_phi*self%phi_hat(i1,i2,i3,2)*d(i1,i2,i3)
        end do
        if (associated(self%filled)) then
          do i1=l(1),u(1)
            self%filled(i1,i2,i3) = self%level
          end do
        end if
      end do
    end do
    call self%E_th2e (e,it)
  end if
  !-----------------------------------------------------------------------------
  ! Outer radial boundary
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%xu)) then
    !---------------------------------------------------------------------------
    ! Apply initial conditions to all variables
    !---------------------------------------------------------------------------
    l = [self%mesh(1)%uo, self%mesh(2)%lb, self%mesh(3)%lb]
    u = [self%mesh(1)%ub, self%mesh(2)%ub, self%mesh(3)%ub]
    call initial_condition (self, l, u, it)
    if (associated(self%filled)) then
      self%filled(l(1):u(1),l(2):u(2),l(3):u(3)) = self%level
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Lower vertical boundary
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zl)) then
    call vertical_boundary (self, self%mesh(3)%lo, it)
  end if
  !-----------------------------------------------------------------------------
  ! Upper vertical boundary
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zu)) then
    call vertical_boundary (self, self%mesh(3)%uo, it)
  end if
  !-----------------------------------------------------------------------------
  ! Restore total energy
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zl+bits%zu)) then
    call self%e2E_th (e,it)
    call self%boundaries%zero_gradient (e, self%idx%e)
    call self%E_th2e (e,it)
  end if
  deallocate (ui,e)
  call lock%unset ('experiment_t%boundary')
  call trace%end (itimer)
END SUBROUTINE boundary

!===============================================================================
!===============================================================================
SUBROUTINE cache_v_head (self)
  class(experiment_t):: self
  integer:: i1, i2, i3
  real(8):: r, z, p(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Cache expensive calculation
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%cache_v_head', itimer=itimer)
  if (.not.associated(self%v_head)) then
    call self%lock%set ('cache_v_head')
    if (.not.associated(self%v_head)) then
      allocate (self%v_head(self%gn(1),self%gn(2),self%gn(3)))
      do i3=self%mesh(3)%lb,self%mesh(3)%ub
        do i2=self%mesh(2)%lb,self%mesh(2)%ub
          do i1=self%mesh(1)%lb,self%mesh(1)%ub
            p(1) = self%mesh(1)%r(i1) + self%radius
            p(2) = self%mesh(2)%r(i2)
            p(3) = self%mesh(3)%r(i3) + self%position(3)
            r = sqrt(p(1)**2 + p(2)**2)
            z = p(3)
            self%v_head(i1,i2,i3) = disk_model%v_head(r,z)
          end do
        end do
      end do
    end if
    call self%lock%unset ('cache_v_head')
  end if
  call trace%end (itimer)
END SUBROUTINE cache_v_head

!===============================================================================
!> Vertical boundary conditions
!===============================================================================
SUBROUTINE vertical_boundary (self, o, it)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer :: d, e, px, py, pz
  integer:: o, it, b, i, step, i1, i2, i3, l(3), u(3)
  real(8):: sign, zi, p(3), r, z, h, d_ratio, dv_phi, v, v_h, v_k
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%vertical_boundary', itimer=itimer)
  d  => self%mem(:,:,:,self%idx%d ,it,1)
  e  => self%mem(:,:,:,self%idx%e ,it,1)
  px => self%mem(:,:,:,self%idx%px,it,1)
  py => self%mem(:,:,:,self%idx%py,it,1)
  pz => self%mem(:,:,:,self%idx%pz,it,1)
  l = self%mesh%lb
  u = self%mesh%ub
  if (o == self%mesh(3)%lo) then
    sign = -1.0
    l(3) = self%mesh(3)%lb
    u(3) = self%mesh(3)%lo
    i = self%mesh(3)%li
  else
    sign = +1.0
    l(3) = self%mesh(3)%uo
    u(3) = self%mesh(3)%ub
    i = self%mesh(3)%ui
  end if
  v_k = disk_model%v_kepler (self%radius)
  zi = self%mesh(3)%r(i) + self%position(3)
  !-----------------------------------------------------------------------------
  call cache_v_head (self)
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        p(1) = self%mesh(1)%r(i1) + self%radius
        p(2) = self%mesh(2)%r(i2)
        p(3) = self%mesh(3)%r(i3) + self%position(3)
        r = sqrt(p(1)**2 + p(2)**2)
        z = p(3)
        !-----------------------------------------------------------------------
        ! Correct for difference in Gaussian density profile
        !-----------------------------------------------------------------------
        h = disk_model%h_scale (r)
        d_ratio = exp(-(z**2-zi**2)/(2.0*h**2))
        d(i1,i2,i3) = d(i1,i2,i) * d_ratio
        e(i1,i2,i3) = e(i1,i2,i) * d_ratio
        !-----------------------------------------------------------------------
        ! Outflow conditions that pull angular momentum (and mass) out
        !-----------------------------------------------------------------------
        px(i1,i2,i3) = px(i1,i2,i) * d_ratio
        v = f_outflow*self%v_patch(i1,i2,i3,2)
        v_h = self%v_head(i1,i2,i3)
        v_h = v_h + self%v_patch(i1,i2,i3,2) - v_k
        py(i1,i2,i3) = d(i1,i2,i3)*(- v + v_h)
        pz(i1,i2,i3) = d(i1,i2,i3)*(+ v*sign)
      end do
      if (associated(self%filled)) then
        do i1=l(1),u(1)
          self%filled(i1,i2,i3) = self%level
        end do
      end if
    end do
  end do
  if (verbose > 0) then
    write (stdout,'(a,i6,":",i2.2,3i4,2(2x,3i4))') &
      'boundary: target,l,u=', self%id, self%key, l, u
  end if
call trace%end (itimer)
END SUBROUTINE vertical_boundary

!===============================================================================
!> Classical Coriolis force and similar force terms
!===============================================================================
SUBROUTINE force (self)
  class(experiment_t):: self
  !.............................................................................
  integer:: i1, i2, i3, iz, j, m, n(3), l(3), u(3)
  real(8):: v_phi
  real:: u1, u2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%force', itimer=itimer)
  m = self%mesh(1)%gn
  l = self%mesh(:)%li
  u = self%mesh(:)%ui
  iz = self%key(4)
  if (self%kind(1:7) /= 'ramses_') then
    call io%abort('experiment_t%force')
  end if
  call scratch%allocate (scratch%force_per_unit_mass, [m,m,m,3])
  !-----------------------------------------------------------------------------
  ! Compute the centripetal force and add to the force of gravity.  If the
  ! local velocity field is zero v_cart becomes equal to the self%v_patch
  ! velocity field, which is parallel to phi_hat, with velocity v_phi = r*omega
  !-----------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !$omp simd
    do i1=l(1),u(1)
      u1 =               self%mem(i1,i2,i3,self%idx%px,self%it,1) / &
                         self%mem(i1,i2,i3,self%idx%d ,self%it,1)
      u2 =               self%mem(i1,i2,i3,self%idx%py,self%it,1) / &
                         self%mem(i1,i2,i3,self%idx%d ,self%it,1)
      v_phi =        self%phi_hat(i1,i2,i3,1)*(self%v_patch(i1,i2,i3,1) + u1) &
                   + self%phi_hat(i1,i2,i3,2)*(self%v_patch(i1,i2,i3,2) + u2)
      scratch%force_per_unit_mass(i1,i2,i3,1)    = &
                      self%f_grav(i1,i2,i3,1,iz) + &
              v_phi**2*self%r_inv(i1,i2,i3,1,iz) + 2.*self%omega(3)*u2
      scratch%force_per_unit_mass(i1,i2,i3,2)    = &
                      self%f_grav(i1,i2,i3,2,iz) + &
              v_phi**2*self%r_inv(i1,i2,i3,2,iz) - 2.*self%omega(3)*u1
      scratch%force_per_unit_mass(i1,i2,i3,3)    = &
                      self%f_grav(i1,i2,i3,3,iz)
    end do
  end do
  end do
  call trace%end (itimer)
END SUBROUTINE force

!===============================================================================
!> Update the position, velocity, and unit vectors.  Each update moves the
!> patch a courant_factor fraction of a cell forward.  The motion relative
!> to the next or previous cylinder is a fraction 'dlnr' of that.  Hence a
!> cylinder needs to move 1/dlnr full cells (= 1 radian), before a difference
!> of one cell has developed.
!===============================================================================
SUBROUTINE update_position (self)
  class(experiment_t):: self
  real(8):: sinphi, cosphi
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update_position')
  self%centre_nat(2) = &
  self%centre_nat(2) + self%dtime*self%omega(3)
  cosphi = cos(self%centre_nat(2))
  sinphi = sin(self%centre_nat(2))
  self%position = [+self%radius*cosphi, &
                   +self%radius*sinphi, 0d0]
  self%velocity = [-self%radius*sinphi*self%omega(3), &
                   +self%radius*cosphi*self%omega(3), 0d0]
  self%erot(:,1) = [ cosphi, sinphi, 0d0]
  self%erot(:,2) = [-sinphi, cosphi, 0d0]
  self%erot(:,3) = [    0d0,    0d0, 1d0]
  if (verbose > 0) then
    write (stdlog, '(a,i6,f10.2)') &
      'experiment_t%update_position: id,cells=', &
      self%id, self%omega(3)*self%time/self%ds(2)
  end if
  call trace%end()
END SUBROUTINE update_position

!===============================================================================
!> Return the position and velocity of the source in the target coordinates at
!> the target time, and the relative orientation of the source at that time.
!===============================================================================
SUBROUTINE get_position (self, source, p, v, e)
  class(experiment_t):: self
  class(task_t):: source
  real(8):: time, p(3), v(3), e(3,3)
  real(8):: r_src, r_tgt, phi_src, phi_tgt, v_src, v_tgt, z_src, z_tgt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%get_position', itimer=itimer)
  select type (source)
  class is (experiment_t)
    r_src   = source%centre_nat(1)
    r_tgt   =   self%centre_nat(1)
    phi_src = source%centre_nat(2) + (self%time-source%time)*source%omega(3)
    phi_tgt =   self%centre_nat(2)
    z_src   = source%centre_nat(3)
    z_tgt   =   self%centre_nat(3)
    v_src   = r_src*source%omega(3)
    v_tgt   = r_tgt*  self%omega(3)
    !---------------------------------------------------------------------------
    ! Relative position and velocity in absolute coordinates
    !---------------------------------------------------------------------------
    p = [ r_src*cos(phi_src), r_src*sin(phi_src), z_src] &
      - [ r_tgt*cos(phi_tgt), r_tgt*sin(phi_tgt), z_tgt]
    v = [-v_src*sin(phi_src), v_src*cos(phi_src), 0d0] &
      - [-v_tgt*sin(phi_tgt), v_tgt*cos(phi_tgt), 0d0]
    !---------------------------------------------------------------------------
    ! Transform from abolute coordinates to local target coordinates
    !---------------------------------------------------------------------------
    p = matmul (p , self%erot)
    v = matmul (v, self%erot)
    !---------------------------------------------------------------------------
    ! Return the relative rotation matrix, at the given target time
    !---------------------------------------------------------------------------
    e(:,1) = [ cos(phi_src-phi_tgt), sin(phi_src-phi_tgt), 0d0]
    e(:,2) = [-sin(phi_src-phi_tgt), cos(phi_src-phi_tgt), 0d0]
    e(:,3) = [                  0d0,                  0d0, 1d0]
  end select
  call trace%end (itimer)
END SUBROUTINE get_position

!===============================================================================
!> Return the velocity shift from source to target, expressed in the local
!> coordinate system.
!===============================================================================
FUNCTION velocity_shift (self, source) RESULT (v)
  class(experiment_t):: self
  class(task_t):: source
  real(8), dimension(3):: v, v_src, v_tgt
  !-----------------------------------------------------------------------------
  ! Velocities in the local coordinate syhstems
  !-----------------------------------------------------------------------------
  select type (source)
  class is (experiment_t)
  v_src = [0d0, disk_model%v_kepler(source%radius), 0d0]
  v_tgt = [0d0, disk_model%v_kepler(  self%radius), 0d0]
  !-----------------------------------------------------------------------------
  ! Transform to the target system
  !-----------------------------------------------------------------------------
  v = matmul (matmul (source%erot, v_src), self%erot) - v_tgt
  end select
END FUNCTION velocity_shift

!===============================================================================
!> Write metadata namelist
!===============================================================================
SUBROUTINE write_nml (self, unit)
  class(experiment_t):: self
  integer:: unit
  namelist /experiment_nml/ omega, centre_nat
  real(8):: omega(3), centre_nat(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%write_nml')
  call self%solver_t%write_nml (unit)
  omega      = self%omega
  centre_nat = self%centre_nat
  write (unit, experiment_nml)
  call trace%end ()
END SUBROUTINE write_nml

END MODULE experiment_mod
