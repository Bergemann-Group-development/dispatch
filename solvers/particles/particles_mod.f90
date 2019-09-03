!===============================================================================
!> Particle list, extending a doubly-linked list.  Each particle maintains 
!> arrays with previous positions, velocities, and times, which are brought
!> along if/when it changes owner patch or rank.  The memory footprint of a
!> particle is 12 words for position, 12 words for velocity, and 8 words for
!> time, all together 32 words = 128 bytes (plus a few words for id and weight).
!> This could be reduced to half, by keeping only two previous positions and
!> velocities in the data type.
!===============================================================================
MODULE particles_mod
  USE io_unit_mod
  USE dll_mod
  USE omp_timer_mod
  USE particle_mod
  USE patch_mod
  USE patch_force_mod
  USE patch_forces_mod
  USE hash_table_mod
  USE trace_mod
  USE scaling_mod
  USE lagrange_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: particles_t
    real:: rsoft=1.0
    real(8):: courant=0.2
    real(8):: rsoft2=tiny(1d0)**(2./3.)
    integer:: verbose=0
    integer:: nt=4
    type(hash_table_t):: ht
  contains
    procedure:: init
    procedure:: force
    procedure:: force_field
    procedure:: predict
    procedure:: courant_time
    procedure:: test
    procedure:: print => print_list
  end type
  type(particles_t), public:: particles
  !-----------------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE init (self, name)
  class (particles_t):: self
  character(len=*), optional:: name
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%init')
  call self%dll_t%init (name)
  !call self%ht%init                            ! either this ..
  !call hash_table%init                         ! .. or this
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Append a new particle_t data type to the particle list
!===============================================================================
SUBROUTINE append (self, new)
  class (particles_t):: self
  class (particle_t), pointer:: new
  class (dll_node_t), pointer:: new_dll
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%append')
  new_dll => new
  call self%dll_t%append (new_dll)
  call trace%end ()
END SUBROUTINE append

!===============================================================================
!> Compute force from other particles by direct summation and add forces from
!> patches
!===============================================================================
FUNCTION force (self, p, it)
  class(particles_t):: self
  class(particle_t), pointer:: p
  integer:: it
  real(8):: force(3)
  !.............................................................................
  class(dll_node_t), pointer:: o
  real(8):: time, ro(3), rp(3), r(3), r1, r2, r3, r2s
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%force')
  r2s = self%rsoft2
  force = 0.0_8
  !-----------------------------------------------------------------------------
  ! Copy relevant particle attributes atomically
  !-----------------------------------------------------------------------------
  !$omp atomic read
  time = p%t(it)
  do i=1,3
    !$omp atomic read
    rp(i) = p%r(i,it)
  end do
  o => self%head
  do while (associated(o))
    if (.not.associated (o, p)) then
      select type (o)
      class is (particle_t)
      !-------------------------------------------------------------------------
      ! Interpolate the other (o) particle position to time = p%t(it)
      !-------------------------------------------------------------------------
      ro = o%interpolate (time, o%r)
      r = ro-rp
      r2 = sum(r**2)
      !-------------------------------------------------------------------------
      ! r is the vector from p to o, so the force is in the same direction
      !-------------------------------------------------------------------------
      r2 = max(r2,r2s)
      r1 = 1.0/sqrt(r2)
      r3 = r1**3
      force = force + o%mass*r*r3
      if (particles%verbose > 1) &
        write (stdout,'(a,2i2,5(1x,3f8.5))') &
          'force: id,oid,ro,rp,rm,r,f =', &
          p%id, o%id, ro, rp, (ro+rp)/2d0, r, force
      end select
    end if
    o => o%next
  end do
  !-----------------------------------------------------------------------------
  ! Scale by the gravitational constant in code units
  !-----------------------------------------------------------------------------
  force = force*scaling%grav
  if (particles%verbose > 0) &
    print *,'sink_patch_t%force: id, f(1) =', p%id, force
  !-----------------------------------------------------------------------------
  ! Add the force from patches, which is summed up by p%patch_forces; this
  ! is already scaled by the gravitational constant when first evaluated
  !-----------------------------------------------------------------------------
  !force = force + p%patch_forces%force (time)
  if (particles%verbose > 0) &
    print *,'sink_patch_t%force: id, f(2) =', p%id, force
  call trace%end ()
END FUNCTION force

!===============================================================================
!> Compute the force field (action) in the cells of a patch, and while doing so,
!> also add up the opposite force (reaction) onto the particles from the cells.
!===============================================================================
SUBROUTINE force_field (self, patch)
  class(particles_t):: self
  class(patch_t):: patch
  !.............................................................................
  class(dll_node_t), pointer:: p
  real(8), allocatable:: rp(:,:), m(:), f(:,:)
  real(8):: force(3), pt, qt, rc(3), r(3), r1, r3, dm, r2, rsoft2, f2, f2max
  integer:: i, it, i0, i1, ip, ix, iy, iz
  class(*), pointer:: ptr
  class(dll_node_t), pointer:: dll_node
  type(patch_force_t), pointer:: patch_force
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Cache all particle position at a given time into an array, in order to be
  ! able to use a vectorized loop over particles
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%force_field', itimer=itimer)
  allocate (rp(self%n,3), m(self%n), f(self%n,3))
  ip = 0
  p => self%head
  do while (associated(p))
    ip = ip+1
    select type (p)
    class is (particle_t)
    !$omp atomic read
    it = p%it
    do i=2,self%nt
      i1 = mod(it+i-1,p%nt)+1
      if (p%t(i1) > patch%time) exit
    end do
    i0 = modulo(i1-2,p%nt)+1
    pt = (patch%time - p%t(i0))/max(p%t(i1)-p%t(i0),tiny(1d0))
    qt = 1d0-pt
    rp(ip,:) = qt*p%r(:,i0) + pt*p%r(:,i1)
    m(ip) = p%mass
    f(ip,:) = 0.0
    end select
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! For each cell, loop over particles and accumulate force in double precision
  !-----------------------------------------------------------------------------
  if (.not.allocated(patch%force_per_unit_mass)) then
    associate (m=>patch%gn)
    allocate (patch%force_per_unit_mass(m(1),m(2),m(3),3))
    end associate
  end if
  associate (m1=>patch%mesh(1), m2=>patch%mesh(2), m3=>patch%mesh(3))
  rsoft2 = (self%rsoft*minval(patch%ds))**2
  f2max = 0.0
  do iz=m3%li,m3%ui
    rc(3) = m3%p + m3%r(iz)
    do iy=m2%li,m2%ui
      rc(2) = m2%p + m2%r(iy)
      do ix=m1%li,m1%ui
        rc(1) = m1%p + m1%r(ix)
        dm = product(patch%ds)*patch%mem(ix,iy,iz,patch%idx%d,patch%it,1)
        force = 0d0
        do ip=1,self%n
          r = rp(ip,:) - rc
          r2 = max(r(1)**2 + r(2)**2 + r(3)**2, rsoft2)
          r1 = 1d0/sqrt(r2)
          r3 = r1*r1*r1
          r = r*r3
          !---------------------------------------------------------------------
          ! Force per unit mass on the cell from the sink particle
          !---------------------------------------------------------------------
          force = force - m(ip)*r
          !---------------------------------------------------------------------
          ! Force per unit mass on the sink particle from the cell
          !---------------------------------------------------------------------
          f(ip,:) = f(ip,:) + dm*r
!print '(a,i6,2x,3i4,i5,1p,2(2x,3e11.2))', &
!'particles_t: ix,iy,iz,ip,f(ip,:),force =', &
!ix,iy,iz,ip,f(ip,:),force
        end do
        force = force*scaling%grav
        if (particles%verbose > 3) then
          f2 = sum(force**2)
          if (f2 > f2max) then
            f2max = f2
            write (stdout,'(a,i6,2x,3i4,2x,1p,3e11.2)') &
              'particles_t: id,ix,iy,iz,force =', &
              patch%id,ix,iy,iz,force
          end if
        end if
        patch%force_per_unit_mass(ix,iy,iz,:) = &
        patch%force_per_unit_mass(ix,iy,iz,:) + force
      end do
    end do
  end do
  end associate
  deallocate (rp)
  !-----------------------------------------------------------------------------
  ! Store the patch forces in the sink particle patch force slots, using a
  ! pointer to the force obtained from a hash table lookup
  !-----------------------------------------------------------------------------
  ip = 0
  p => self%head
  do while (associated(p))
    ip = ip+1
    select type (p)
    class is (particle_t)
      call hash_table%get ([p%id, patch%id], ptr)
      select type (ptr)
      class is (patch_force_t)
        f(ip,:) = f(ip,:)*scaling%grav
        ptr%f(ptr%it,:) = f(ip,:)
        if (particles%verbose > 1) then
          write(stdout,'(a,2(2i4,2x),1p,e12.4,2x,3e10.2)') &
            'force_field: key,ip,it,time,f =', &
            p%id, patch%id, ip, ptr%it, p%time, f(ip,:)
        end if
      end select
    end select
    p => p%next
  end do
  call trace%end (itimer)
END SUBROUTINE force_field

!===============================================================================
!> Estimate the position and velocity at a future time -- by default the
!> approximately reflexive half time step.
!> 
!> FIXME:  This could/should be used to put the other particles at preliminary
!> positions during the first time step, when extrapolation of positions is not
!> possible, since there are no particle histories.
!===============================================================================
SUBROUTINE predict (self, p, dt)
  class(particles_t):: self
  class(particle_t), pointer:: p
  real(8), optional:: dt
  !.............................................................................
  integer:: prv, now, new, nstep, i
  real(8):: x0, x(p%nt), a(3), dt_l
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_t%predict')
  !$omp atomic read
  nstep = p%nstep
  !$omp atomic read
  now = p%it
  new = mod(now,p%nt) + 1
  !-----------------------------------------------------------------------------
  ! The very first time there is no history to extrapolate from, so we need
  ! a direct estimate, first of the timestep, then of the half-way position
  !-----------------------------------------------------------------------------
  if (nstep==1) then
    if (present(dt)) then
      dt_l = dt
    else
      dt_l = 0.5d0*self%courant_time (p)
    end if
    a = self%force (p, now)
    p%r(:,new) = p%r(:,now) + dt_l*p%v(:,now) + 0.5d0*dt_l**2*a
    p%v(:,new) = p%v(:,now) + dt_l*a
    p%t(  new) = p%t(  now) + dt_l
  else
    prv = modulo(now-2,p%nt) + 1
    p%r(:,new) = p%r(:,now)*1.5_8 - p%r(:,prv)*0.5_8
    p%v(:,new) = p%v(:,now)*1.5_8 - p%v(:,prv)*0.5_8
    p%t(  new) = p%t(  now)*1.5_8 - p%t(  prv)*0.5_8
  end if
  call trace%end()
END SUBROUTINE predict

!===============================================================================
!> Courant conditions:  Estimate, as inexpensively as possible, an approximately
!> reflexive time step. If an estimate of the acceleration results in it not
!> being important for the timestep size, then choose a fraction of the cell
!> travel time (whatever that means in this context)
!===============================================================================
FUNCTION courant_time (self, p) RESULT (dt)
  class(particles_t):: self
  class(particle_t), pointer:: p
  real(8):: dt
  !.............................................................................
  real(8):: a(3), ac, vc, dt1, dt2
  integer:: prv, now, new
  !-----------------------------------------------------------------------------
  call trace%begin('particle_solver_t%courant_time')
  if (self%courant < 0.0) then
    dt = -self%courant
  else
    now = p%it
    new = mod   (now  ,p%nt) + 1
    prv = modulo(now-2,p%nt) + 1
    call p%predict()
    call courant_expressions
    if (p%nstep==1) then
      dt = dt*0.5d0
      p%r(:,new) = p%r(:,now) + dt*p%v(:,now) + 0.5d0*dt**2*a
      call courant_expressions
    end if
    !-----------------------------------------------------------------------------
    if (self%verbose > 0) &
      print '(a,i4,1p,4e12.3)', &
        'courant_time: id,dt1,dt2,ac,vc =', p%id, dt1, dt2, ac, vc
  end if
  call trace%end()
contains
  subroutine courant_expressions
    associate (v=>p%v(:,new))
    a = self%force(p,new)
    ac = 1.0/sqrt(sum(a**2)+tiny(1.0_8))
    vc = 1.0/sqrt(sum(v**2)+tiny(1.0_8))
    dt1 = p%ds*vc*self%courant                ! one  cell  per step at courant=1      
    dt2 = (ac/vc)*self%courant*0.2            ! 0.2 radian per step at courant=1
    dt = min(dt1,dt2)
    end associate
  end subroutine courant_expressions
END FUNCTION courant_time

!===============================================================================
!> Test particle list functionality
!===============================================================================
SUBROUTINE test (self)
  class(particles_t):: self
  integer:: id, n
  class(particle_t), pointer:: p
  class(dll_node_t), pointer :: node, next, new
  real(8)::used
  logical:: flip
  !-----------------------------------------------------------------------------
  call self%init
  !-----------------------------------------------------------------------------
  ! Make a list with 100000 particles
  !-----------------------------------------------------------------------------
  n = 100000
  used = wallclock()
  do id=1,n
    allocate (p)
    call p%init
    p%r  = 0.5_8
    p%ds = 1.0_8
    p%v  = 1.0_8
    node => p
    call self%append (node)
  end do
  used = wallclock()-used
if (particles%verbose > 0) &
print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Remove every 2nd particle
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    select type (node)
    type is (particle_t)
    p => node
    end select
    if (flip) then
      call self%remove (node)
      deallocate (node)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
if (particles%verbose > 0) &
print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Add one new particle for every 2nd existing one
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    allocate (p)
    new => p
    if (flip) then
      call self%insert_before (node, new)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
if (particles%verbose > 0) &
print *, 'particle list n =', self%n, used/self%n
END SUBROUTINE test

!===============================================================================
!> Print a particle list
!===============================================================================
SUBROUTINE print_list (self)
  class(particles_t):: self
  class(dll_node_t), pointer:: p
  p => self%head
  do while (associated(p))
    select type (p)
    class is (particle_t)
    print *, p%id, p%r(:,p%it)
    end select
    p => p%next
  end do
END SUBROUTINE print_list

END MODULE particles_mod
