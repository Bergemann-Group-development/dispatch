MODULE disk_particles_mod
  USE io_mod
  USE scaling_mod
  USE task_mod
  USE particle_solver_mod
  USE particles_mod
  USE gpatch_mod
  USE disk_model_mod
  USE vector_ops
  USE link_mod
  USE math_mod
  implicit none
  private
  type, public, extends(particle_solver_t):: disk_particles_t
    real :: m_star
    real :: omega
    real :: gm
    logical :: ICs_done = .false.
    real(8), pointer, contiguous:: v_patch(:,:,:,:) => null()
  contains
    procedure:: init
    procedure:: initial_condition
    procedure:: set_attributes
    procedure:: particle_force
    procedure:: import
    procedure:: refine_Flash
  end type
CONTAINS

!===============================================================================
!> Initialize -- a chance to add parameters, or just pass the call on
!===============================================================================
SUBROUTINE init (self)
  class(disk_particles_t):: self
  !.............................................................................
  if (.not.self%on) return
  call self%set_attributes
  call self%particle_solver_t%init
END SUBROUTINE init

!===============================================================================
!> Set attributes by call rather then accessing the attributes directly
!===============================================================================
SUBROUTINE set_attributes (self)
  class(disk_particles_t):: self
  !.............................................................................
  self%m_star = disk_model%m_star
  !FIXME: self%omega = self%patch%omega(3)
  self%gm = scaling%grav*self%m_star
END SUBROUTINE set_attributes

!===============================================================================
!> Particle initial conditions: power-law mid-plane density profile, Gaussian
!> vertical density profile, Keplerian velocities
!===============================================================================
SUBROUTINE initial_condition (self)
  class(disk_particles_t)     :: self
  class(particles_t), pointer :: particles
  class(particle_t), pointer  :: p
  type(particle_t), pointer   :: pn
  class(gpatch_t), pointer    :: patch
  integer                     :: m(3), i1, i2, i3, i4, j, n, nv, nb, id(5)
  real(8)                     :: dg, w, w1, vg(3), chance, h, dmid, v_k
  real(8)                     :: x(3), r, z
  logical                     :: draw
  integer, save               :: itimer=0, extra_bits
  integer, save               :: lid = 0
  !$omp threadprivate (lid)
  !-----------------------------------------------------------------------------
  ! Set flags indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  if (.not. self%on .or. self%ICs_done) return
  self%ICs_done = .true.
  call trace%begin ('disk_particles_t%initial_condition', itimer=itimer)
  patch => self%patch
  extra_bits = 200
  if (self%patch%id == 1) then
    allocate (p)
    write (stdout,'(a,4i6)') 'particle storage_size:', &
      storage_size(p)
    write (stdout,'(a,6i6)') 'particle storage_size:', &
      storage_size(p%t), storage_size(p%r), &
      storage_size(p%q), storage_size(p%v), storage_size(p%w), &
      storage_size(p%t) + 3*storage_size(p%r) + 4*storage_size(p%q) + &
      storage_size(p%w) + 3*storage_size(p%v) + 2*64
    call p%init
    nb = (storage_size(p)+extra_bits)/8
    print '(a,i3,a,f18.3,a)', ' estimated memory for particles (', &
      nb,' bytes):',nb*self%n_bin*self%particles_per_bin/1024.**3 &
      *product(self%patch%mesh%b/self%patch%mesh%d), ' GB'
    call pool%append (p)
  end if
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
#ifdef __GFORTRAN__
  !$omp critical (alloc_cr)
#endif
  if (self%output_every > 0) then
    chance = 1.0/self%output_every
  else
    chance = 0.0
  end if
  self%np = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
    n = 0
    do i1=self%lb(1),self%ub(1)
      associate (r1=>patch%mesh(1)%r, r2=>patch%mesh(2)%r, r3=>patch%mesh(3)%r)
      x(1) = r1(i1) + patch%radius
      x(2) = r2(i2)
      x(3) = r3(i3) + self%position(3)
      end associate
      r = sqrt(x(1)**2 + x(2)**2)                       ! radius
      z = x(3)                                          ! height
      h = disk_model%h_scale_dust (r)                   ! scale height
      dmid = disk_model%d_midplane_dust (r)             ! mid-plane density
      w = dmid*exp(-(z/h)**2/2d0) / &                   ! Gaussian profile
        (self%n_bin*self%particles_per_bin)
      dg = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)   ! gas density
      vg = patch%mem(i1,i2,i3,patch%idx%px:patch%idx%pz,patch%it,1) / &
          (patch%ds*dg)                                 ! gas velocity
        do i4=1,self%n_bin
        particles => self%particles(i1,i2,i3,i4)
        id = [self%patch%id,i1,i2,i3,i4]
        call particles%init (id, 0)
        w1 = w*(self%s_cgs(i4)/self%s_cgs(1))**(self%alpha_exponent+1.0)
        n = n+2
        nv = 0
        do j=1,self%particles_per_bin
          n = n+1
          if (self%output_every > 0) then
            draw = self%random%ran1() < chance
          else
            draw = .false.
          end if
          if (draw) then  ! every output_every is for I/O
            allocate(pn)
            p => pn
            call p%init_vparticle
            p%cell%id = 1 + omp%thread_global + (lid+nv)*omp%nthreads_global
            nv = nv+1
            p%cell%dg = dg
            p%cell%vg = vg
          else
            allocate(p)
            call p%init
          end if
          p%q = [i1,i2,i3,i4]
          p%r = self%random%uniform (3)
          p%v = vg
          p%w = w1
          p%t = 0
          call particles%append (p)
        end do
        self%np = self%np + particles%n
      end do
      lid = lid+nv
    end do
    call io%bits_mem (storage_size(p)+extra_bits, n, 'p')
  end do
  end do
#ifdef __GFORTRAN__
  !$omp end critical (alloc_cr)
#endif
  !.............................................................................
  if (self%verbose > 0) then
    if (self%verbose > 1) then
      write (io%output,1) &
        'init: id, box, pos, size, lower, upper, n =', &
        patch%id, patch%box, patch%position, patch%size, &
        self%lower, self%upper, &
        real(self%np)
      1 format(a,i4,3(1x,3f7.3),1x,2(1x,3l1),1p,e11.2)
    end if
    write (stdlog,1) &
      'init: id, box, pos, size, lower, upper, n =', &
      patch%id, patch%box, patch%position, patch%size, &
      self%lower, self%upper, &
      real(self%np)
  else
    write(stdlog,'(a,i6,2(3f8.3,2x),1x,2(1x,3l1),1p,e11.2)') &
     'disk_particles_t%init: id, pos, box, lower, upper, n =', &
      patch%id, patch%position, patch%box, self%lower, self%upper, real(self%np)
  end if
  call trace%end (itimer)
END SUBROUTINE initial_condition

!===============================================================================
!> External force on a particle
!===============================================================================
FUNCTION particle_force (self, p, vp, grad_phi, v) RESULT (a)
  class(disk_particles_t)                  :: self
  class(gpatch_t), pointer                 :: patch
  real, dimension(3)                       :: p, v, a, grad_phi
  real, dimension(2)                       :: vp
  real                                     :: ri, v_phi, ds, cs
  real, dimension(:,:,:), pointer          :: phi
  !-----------------------------------------------------------------------------
  patch => self%patch
  ds = patch%ds(1)
  cs = 1.0/ds
  !-----------------------------------------------------------------------------
  ! Correct particle position
  !-----------------------------------------------------------------------------
  p(1) = p(1) + patch%radius
  p(3) = p(3) + patch%position(3)
  ri = 1d0/sqrt(p(1)**2+p(2)**2)
  !-----------------------------------------------------------------------------
  ! Stellar gravity, centripetal and Coriolis forces
  !-----------------------------------------------------------------------------
  v_phi = (-(v(1)+vp(1))*p(2)*ri + (v(2) + vp(2))*p(1)*ri)*ds
  a(1) = (- self%gm*ri + v_phi**2)*ri*ri*p(1)*cs + 2d0*self%omega*v(2)
  a(2) = (- self%gm*ri + v_phi**2)*ri*ri*p(2)*cs - 2d0*self%omega*v(1)
  a(3) =  - self%gm*ri            *ri*ri*p(3)*cs
  !-----------------------------------------------------------------------------
  ! Self-gravity
  !-----------------------------------------------------------------------------
  a(1) = a(1) + grad_phi(1)*cs
  a(2) = a(2) + grad_phi(2)*cs
  a(3) = a(3) + grad_phi(3)*cs
END FUNCTION particle_force

!===============================================================================
!> Import particles from nbors export lists
!> Verbatim copy from particle_solver_mod, with coordinate and velocity
!> transformations added
!===============================================================================
SUBROUTINE import (self, count)
  class(disk_particles_t)    :: self
  logical, optional          :: count
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particles_t), pointer:: particles
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3), dist(3), pt(3)
  integer                    :: ni, j, bin, n(4), n_before, n_after
  integer(kind=2)            :: q(3)
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('disk_particles_t%import',itimer=itimer)
  self%ni = 0
  n_before = 0
  n_after = 0
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      call nbtask%lock%set ('import')
      dist = self%patch%distance_to(nbtask%patch) / self%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      ni = 0
      n = 0
      do bin=1,self%n_bin
        p => nbtask%exports(bin)%head%next
        n(1) = nbtask%exports(bin)%n
        n_before = n_before + nbtask%exports(bin)%n
        n(3) = self%np
        do j=1,nbtask%exports(bin)%n
          if (abs(maxval(p%q)) > 100) write (stdlog,*) 'import error', j, p%q
          next => p%next
          where (p%r==1d0) p%r=0.9999d0
          ! Transform particle coordinates from nbor to self
          pt = matmul (matmul (nbtask%patch%erot, &
            real(p%r,kind=8)+real(p%q(1:3),kind=8)), self%patch%erot)
          r = dist + pt
          q = floor(r)
          if (all(q >= self%lb .and. q <= self%ub)) then
            p%r = r-q
            p%q(1:3) = q
            if (self%verbose > 2) then
              write(stdlog,*) 'import: q=', self%patch%id, p%q
            end if
            call nbtask%exports(bin)%remove (p)
            particles => self%particles(p%q(1),p%q(2),p%q(3),p%q(4))
            ! Transform particle velocity from nbor to self
            p%v = matmul (matmul (nbtask%patch%erot, p%v), self%patch%erot) - &
               self%patch%velocity_shift (nbtask)
            call particles%append (p)
            ni = ni + 1
          else if (self%verbose > 2) then
            write (stdlog,'(2i6,3f10.3,2x,a,4(2x,3i6))') &
              self%patch%id,nbtask%id,dist,'rejected',p%q,q,self%lb,self%ub
          end if
          p => next
        end do
        n_after = n_after + nbtask%exports(bin)%n
        self%np = self%np + ni
        n(2) = nbtask%exports(bin)%n
        n(4) = self%np
        if (self%verbose > 1) then
          if (ni > 0) then
            write(stdlog,'(2(i6,2x,a),i6,4i8)') &
              self%patch%id, 'picked up', ni, 'from', nbtask%id, n
          end if
        end if
        if ((ni > 0 .and. self%verbose > 1) .or. self%verbose > 3) then
          write (stdlog,'(a,2i6,2i10,2x,3i4)') 'imported: id, nbor, ni, offset =', &
            self%patch%id, nbtask%patch%id, nbtask%exports(bin)%n, ni, &
            nint(dist)
        end if
      end do
      self%ni = self%ni + ni
      call nbtask%lock%unset ('import')
    end select
    nbor => nbor%next
  end do
  self%np = self%np + self%ni
  if ((self%ni > 0 .and. self%verbose > 1) .or. self%verbose > 3  ) then
    write (stdlog,'(a,i6,i10)') 'imported: id, self%ni =', &
      self%patch%id, self%ni
    flush (stdlog)
  end if
  call trace%end (itimer)
END SUBROUTINE import

!===============================================================================
!> Adoption of the refinement criterion of the FLASH simulations presented in
!> Schäfer et al. (2020). In their simulations, the resolution is increased if
!> the number of dust particles in a cell exceeds 10 and decreased if it is less
!> than 1.
!> We consider their simulations spanning 10 au to 100 au in the radial
!> dimension, but assume that the vertical gas density profile is given by
!> rho_g(z) = rho_g(z=0)*exp(-z^2/2H_g^2) rather than their Eq. 6.
!> Under these conditions, the mass of each dust particle (their Eq. 15) is
!> equal to m_d = 1.359e26*dust_to_gas g, where dust_to_gas is the dust-to-gas
!> surface density ratio.
!> Their simulations are two-dimensional and do not include the azimuthal
!> dimension. Each particle therefore represents the entire mass in a ring with
!> the circumference 2*pi*r, where r is the radius. To adopt their refinement
!> criterion for three-dimensional DISPATCH simulations, we thus multiply the
!> above mass threshold by the ratio of azimuthal cell size to 2*pi*r.
!===============================================================================
FUNCTION refine_Flash (self) RESULT(refine)
  integer :: refine
  !.............................................................................
  class(disk_particles_t)               :: self
  class(gpatch_t), pointer              :: patch
  class(mesh_t), pointer, dimension(:)  :: mesh
  integer                               :: l(3), u(3), ix, iy, iz
  real(8)                               :: r, mass_thres, mass
  !-----------------------------------------------------------------------------
  call trace%begin ('disk_particles_t%refine_Flash')
  patch => self%patch
  mesh => self%patch%mesh
  call self%density (self%dd)
  l = mesh%li
  u = mesh%ui
  do iz = l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        r = sqrt((patch%radius+mesh(1)%r(ix))**2 + mesh(2)%r(iy)**2)
        mass_thres = 1.359e26*disk_model%sigma_ratio*mesh(2)%d / &
          (2d0*math%pi*r*scaling%m)
        mass = sum (self%dd(ix,iy,iz,:))*mesh(1)%d*mesh(2)%d*mesh(3)%d
        if (mass >= 10d0 * mass_thres) then
          patch%frefine(ix,iy,iz) = 1.
        else if (mass < mass_thres) then
          patch%frefine(ix,iy,iz) = -1.
        end if
      end do
    end do
  end do
  refine = 0
  if (maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3))) == 1.) then
    refine = 1
  else if (any(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)) == -1.)) then
    refine = -1
  end if
  call trace%end ()
END FUNCTION refine_Flash

END MODULE disk_particles_mod
