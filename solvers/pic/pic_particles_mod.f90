!===============================================================================
!> Particle solver for DISPATCH/PIC
!>
!> To increase vectorization, compute the E and B field from the list, while
!> while saving particle data in arrays, conditioned on that they still need to
!> be moved.  Then perform one particle move based on the list, put data back in
!> list, while checking if they still are inside, and repeat until all particles
!> have been moved.
!===============================================================================
MODULE pic_particles_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE dll_mod
  implicit none
  private
  integer, parameter:: nt=5
  !-----------------------------------------------------------------------------
  ! PIC particle data type
  !-----------------------------------------------------------------------------
  type, public:: pic_particle_t
    type(pic_particle_t), pointer:: prev, next
    real(8)         :: t(nt)
    real(8)         :: dt
    real            :: r(3,nt)
    real            :: p(3)
    real            :: e(3)
    integer(kind=2) :: q(3,nt)
    integer(kind=2) :: it=1
    integer         :: n=0
  end type
  !-----------------------------------------------------------------------------
  ! List of PIC particles
  !-----------------------------------------------------------------------------
  type:: pic_list_t
    class(pic_particle_t), pointer:: head => null()
    class(pic_particle_t), pointer:: tail => null()
    integer:: n=0
  contains
    procedure:: init => init_list
    procedure:: append
    procedure:: remove
    procedure:: delete
  end type
  !-----------------------------------------------------------------------------
  ! The particles_t data type contains everything needed by a particle solver.
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: pic_particles_t
    type(pic_list_t)        :: particle_list
    type(pic_list_t)        :: export_list
    class(patch_t), pointer :: patch => null()
    real                    :: mass, charge
    integer                 :: np
  contains
    procedure:: init
    procedure:: update
    procedure:: test
  end type
  real:: dt=2.0
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Initialize particle solver
!===============================================================================
SUBROUTINE init (self)
  class(pic_particles_t):: self
  integer:: iostat
  real:: time=1.0
  namelist /pic_particles_params/ verbose, time, dt
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%init')
  call io%lock%set
  if (first_time) then
    rewind (io_unit%input)
    read (io_unit%input, pic_particles_params, iostat=iostat)
    !call io%namelist_warning ('pic_particles_params', iostat > 0)
    write (io_unit%output, pic_particles_params)
    first_time = .false.
  end if
  self%time = time
  call io%lock%unset
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(pic_particles_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%update')
  call import_particles (self)
  call   move_particles (self)
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Move the particles inside the solver past self%time
!===============================================================================
SUBROUTINE move_particles_nv (self)
  USE units, ONLY: c, elm
  class(pic_particles_t)            :: self
  real, dimension(3)                :: r, dr, p
  integer, dimension(3)             :: q, ir
  real   , dimension(3)             :: ods ,EE, BB, t, p1
  real                              :: hE, hB, o, t2, s, ut
  real                              :: og, p2, pp
  real                              :: mass, com
  integer                           :: i, ip, isp, it, it0, it1, np, nmv
  integer(8)                        :: iter
  real, parameter                   :: sqrt2 = sqrt(2.)
  real, allocatable                 :: p_EE(:,:), p_BB(:,:)
  integer(kind=2), allocatable      :: p_q(:,:)
  class(pic_particle_t), pointer    :: pa, next
  class(*), pointer                 :: car
  real(8)                           :: time(2), tmin, tmax
  logical                           :: outside
!-------------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%move_particles')
  ods = 1.0/self%patch%ds
  np = self%particle_list%n
  allocate (p_EE(np,3),p_BB(np,3))
  time(1) = wallclock()
  !-----------------------------------------------------------------------------
  !> Non-vectorized Vay pusher
  !-----------------------------------------------------------------------------
  ip = 0
  iter = 0
  com = self%charge/self%mass
  pa => self%particle_list%head
  do while (associated(pa))
    next => pa%next
      !-------------------------------------------------------------------------
      ! Loop until the particle has moved enough (could be zero times)
      !-------------------------------------------------------------------------
      nmv = 0
      !do while (pa%t(pa%it) < self%time)
        call scatter_fields2 (pa, EE, BB)                     ! E&B at particle
        it0 = pa%it                                           ! time slot
        it1 = 1+mod(it0,nt)
        hE = pa%dt*com                                        ! time centered
        hB = hE*elm%kf*0.5                                    ! B prefactor
        r(1) = pa%r(1,it); r(2) = pa%r(2,it); r(3) = pa%r(3,it)
        q(1) = pa%q(1,it); q(2) = pa%q(2,it); q(3) = pa%q(3,it)
        p(1) = pa%p(1)   ; p(2) = pa%p(2)   ; p(3) = pa%p(3)
        o = 1./sqrt(1.+sum(p**2)*c%oc2)
        t(1) = hB*BB(1); t(2) = hB*BB(2); t(3) = hB*BB(3)
        !t2 = t(1)**2 + t(2)**2 + t(3)**2
        t2 = sum(t**2)
        p1(1) = p(1) + 2.*hE*EE(1) + o*(p(2)*t(3)-p(3)*t(2))
        p1(2) = p(2) + 2.*hE*EE(2) + o*(p(3)*t(1)-p(1)*t(3))
        p1(3) = p(3) + 2.*hE*EE(3) + o*(p(1)*t(2)-p(2)*t(1))
        !s = 1. + (p1(1)**2+p1(2)**2+p1(3)**2)*c%oc2 - t2
        s = 1. + sum(p1**2)*c%oc2 - t2
        pp = sum(p1*t)*c%oc
        o = sqrt2 * 1./sqrt(s + sqrt(s**2 + 4.*(t2 + pp**2)))
        t(1) = t(1)*o; t(2) = t(2)*o; t(3) = t(3)*o
        s = 1./(1 + t2*o**2)
        ut = sum(p1*t)
        p(1) = s * (p1(1) + ut*t(1) + p1(2)*t(3) - p1(3)*t(2))
        p(2) = s * (p1(2) + ut*t(2) + p1(3)*t(1) - p1(1)*t(3))
        p(3) = s * (p1(3) + ut*t(3) + p1(1)*t(2) - p1(2)*t(1))
        p2 = sum(p**2)
        og = 1./sqrt(1.+p2*c%oc2)                             ! 1./gamma
        pa%p(1) = p(1)                                        ! v*gamma
        pa%p(2) = p(2)                                        ! v*gamma
        pa%p(3) = p(3)                                        ! v*gamma
        pa%e = self%mass*p2*og/(1.+og)                        ! kinetic energy
        og=og*pa%dt
        dr(1) = p(1)*og*ods(1)
        dr(2) = p(2)*og*ods(2)
        dr(3) = p(3)*og*ods(3)
        r(1) = pa%r(1,it0) + dr(1)
        r(2) = pa%r(2,it0) + dr(2)
        r(3) = pa%r(3,it0) + dr(3)
        q(1) = floor(r(1))
        q(2) = floor(r(2))
        q(3) = floor(r(3))
        pa%r(1,it1) = r(1) - q(1)
        pa%r(2,it1) = r(2) - q(2)
        pa%r(3,it1) = r(3) - q(3)
        pa%q(1,it1) = pa%q(1,it0) + q(1)
        pa%q(2,it1) = pa%q(2,it0) + q(2)
        pa%q(3,it1) = pa%q(3,it0) + q(3)
        pa%t(it1) = pa%t(it0) + pa%dt
        pa%it = it1
        iter = iter+1
        nmv = nmv+1
        if (verbose > 2) then
          write (stdout,'(a,i6,2(2x,f12.6,1x,3f12.6))') &
            'move_particles_nv: it, after, before', pa%it, &
            pa%t(it1)      , pa%q(:,it1)+pa%r(:,it1), &
            sqrt(sum(p**2)), pa%q(:,it0)+pa%r(:,it0)
        end if
      !end do
      outside = is_outside (pa, self%patch)
    if (verbose > 1) then
      write (stdout,'("move_particles_nv: ip, nmv, nexp, n, outside =",4i6,l4)') &
        ip, nmv, self%export_list%n, self%particle_list%n, outside
    end if
    !---------------------------------------------------------------------------
    ! Move the particle to the export list if it falls outside
    !---------------------------------------------------------------------------
    if (outside) then
      call self%export_list%append (pa)
      call self%particle_list%remove (pa)
    end if
    pa => next
    ip = ip+1
  end do
  if (verbose > 0) then
    time(2) = wallclock()-time(1)
    write (stdout,'(i6,a,9p,f8.1,a)') &
      ip, ' updates', time(2)/iter, ' ns/update non-vectorized'
  end if
  deallocate (p_EE, p_BB)
  call trace%end ()
END SUBROUTINE move_particles_nv

!===============================================================================
!> Move the particles inside the solver past %time
!===============================================================================
SUBROUTINE move_particles_v (self)
  USE units, ONLY: c, elm
  type particle_ptr
    class(pic_particle_t), pointer :: p
  end type
  class(pic_particles_t)            :: self
  real, dimension(3)                :: r, dr, p
  integer, dimension(3)             :: q, ir
  real   , dimension(3)             :: ods ,EE, BB, t, p1
  real                              :: hE, hB, o, t2, s, ut
  real                              :: og, p2, pp
  real                              :: mass
  integer                           :: i, ip, isp, it, it0, it1, np, iter
  real, parameter                   :: sqrt2 = sqrt(2.)
  real, allocatable                 :: p_r(:,:), p_p(:,:), p_t(:), p_dt(:)
  real, allocatable                 :: p_EE(:,:), p_BB(:,:)
  integer(kind=2), allocatable      :: p_q(:,:)
  type(particle_ptr), allocatable   :: p_ptr(:)
  class(pic_particle_t), pointer        :: pa, next, pa_dll
  class(*), pointer                 :: car
  real(8)                           :: time(9), tmin, tmax
  logical                           :: outside
!-------------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%move_particles')
  ods = 1.0/self%patch%ds
  np = self%particle_list%n
  allocate (p_r(np,3),p_q(np,3),p_p(np,3),p_t(np),p_dt(np),p_EE(np,3),p_BB(np,3))
  allocate (p_ptr(np))
  !-----------------------------------------------------------------------------
  ! Loop as long as there still remains particles to be moved
  !-----------------------------------------------------------------------------
  iter = 0
  time(5) = wallclock()
  !do while (np > 0)
    time(1) = wallclock()
    !---------------------------------------------------------------------------
    ! Make arrays of values for particles that need update
    !---------------------------------------------------------------------------
    ip = 0
    pa => self%particle_list%head                           ! particle list
    do while (associated(pa))
      it = pa%it                                        ! current slot
      if (pa%t(it) < self%time) then
        !call scatter_fields2 (pa, EE, BB)               ! E&B at particle
        ip = ip+1                                       ! count active partd
        p_t (ip)   = pa%t(it)                           ! time
        p_ptr(ip)%p => pa                               ! reverse ptr
        p_dt(ip)   = pa%dt                              ! timestep
        !do i=1,3
          p_EE(ip,:) = 0.0                              ! E-field
          p_BB(ip,:) = 0.0                              ! B-field
         !p_EE(ip,:) = EE(:)                            ! E-field
         !p_BB(ip,:) = BB(:)                            ! B-field
          p_p (ip,:) = pa%p(:)                          ! momentum
          p_r (ip,:) = pa%r(:,it)                       ! float position
          p_q (ip,:) = pa%q(:,it)                       ! int position
        !end do
      end if
      pa => pa%next
    end do
    np = ip
    time(2) = wallclock()
    !-----------------------------------------------------------------------------
    !> Vectorized Vay pusher
    !-----------------------------------------------------------------------------
    do ip=1,np
      hE = 0.5*p_dt(ip)*self%charge/self%mass                 ! time centered
      hB = hE*elm%kf                                          ! B prefactor
      EE(1) = p_EE(ip,1); EE(2) = p_EE(ip,2); EE(3) = p_EE(ip,3)
      BB(1) = p_BB(ip,1); BB(2) = p_BB(ip,2); BB(3) = p_BB(ip,3)
      r(1) = p_r(ip,1); r(2) = p_r(ip,2); r(3) = p_r(ip,3)
      q(1) = p_q(ip,1); q(2) = p_q(ip,2); q(3) = p_q(ip,3)
      p(1) = p_p(ip,1); p(2) = p_p(ip,2); p(3) = p_p(ip,3)
      o = 1./sqrt(1.+(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))*c%oc2)
      t(1) = hB * BB(1); t(2) = hB * BB(2); t(3) = hB * BB(3)
      t2 = t(1)*t(1)+t(2)*t(2)+t(3)*t(3)
      p1(1) = p(1) + 2.*he*EE(1) + o*(p(2)*t(3)-p(3)*t(2))
      p1(2) = p(2) + 2.*he*EE(2) + o*(p(3)*t(1)-p(1)*t(3))
      p1(3) = p(3) + 2.*he*EE(3) + o*(p(1)*t(2)-p(2)*t(1))
      s = 1. + (p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))*c%oc2 - t2
      pp = (p1(1)*t(1) + p1(2)*t(2) + p1(3)*t(3))*c%oc
      o = sqrt2 * 1./sqrt(s + sqrt(s*s + 4.*(t2 + pp*pp)))
      t(1) = t(1)*o; t(2) = t(2)*o; t(3) = t(3)*o
      s = 1./(1 + t2*o*o)
      ut = p1(1)*t(1) + p1(2)*t(2) + p1(3)*t(3)
      p(1) = s * (p1(1) + ut*t(1) + p1(2)*t(3) - p1(3)*t(2))
      p(2) = s * (p1(2) + ut*t(2) + p1(3)*t(1) - p1(1)*t(3))
      p(3) = s * (p1(3) + ut*t(3) + p1(1)*t(2) - p1(2)*t(1))
      !pa%p = p
      p2= p(1)*p(1) + &
          p(2)*p(2) + &
          p(3)*p(3)
      og = 1./sqrt(1.+p2*c%oc2)                             ! 1./gamma
      !pa%e = self%mass*p2*og/(1.+og)                        ! kinetic energy
      og=og*p_dt(ip)*1e-10
      dr(1) = p(1)*og*ods(1)
      dr(2) = p(2)*og*ods(2)
      dr(3) = p(3)*og*ods(3)
      !drmax = max(drmax,dr(1),dr(2),dr(3))
      r(1) = p_r(ip,1) + dr(1)
      r(2) = p_r(ip,2) + dr(2)
      r(3) = p_r(ip,3) + dr(3)
      q(1) = floor(r(1))
      q(2) = floor(r(2))
      q(3) = floor(r(3))
      p_r(ip,1) = r(1) - q(1)
      p_r(ip,2) = r(2) - q(2)
      p_r(ip,3) = r(3) - q(3)
      p_q(ip,1) = p_q(ip,1) + q(1)
      p_q(ip,2) = p_q(ip,2) + q(2)
      p_q(ip,3) = p_q(ip,3) + q(3)
      !p_r(ip,1) = 0.5
      !p_r(ip,2) = 0.5
      !p_r(ip,3) = 0.5
      !p_q(ip,1) = 0
      !p_q(ip,2) = 0
      !p_q(ip,3) = 0
      p_t(ip) = p_t(ip) + p_dt(ip)
    end do
    time(3) = wallclock()
    !---------------------------------------------------------------------------
    ! Put results back in list, and prune.  The DLL list corresponds one-to-one
    ! to the array, but particles that have reach self%time, or gone out-of-
    ! bounds will be removed
    !---------------------------------------------------------------------------
    do ip=1,np
      pa => p_ptr(ip)%p
      it = 1 + mod(pa%it,nt)
      pa%r(:,it) = p_r(ip,:)
      pa%q(:,it) = p_q(ip,:)
      pa%p       = p_p(ip,:)
      pa%t(it)   = p_t(ip)
      pa%it = it
      outside = is_outside (pa, self%patch)
      if (outside) then
        call self%export_list%append (pa)
        call self%particle_list%remove (pa)
      end if
    end do
    time(4) = wallclock()
  !end do ! while (np > 0)
  time(6) = wallclock() - time(5)
  time(1:3) = time(2:4) - time(1:3)
  iter = np
  if (verbose > 1) then
    write (stdout,'(i6,a,9p,4f8.1,1x,a)') &
      iter, ' updates', time(1:3)/iter, time(6)/iter, 'ns/upd vectorized'
  else if (verbose > 0) then
    write (stdout,'(i6,a,9p,f8.1,a)') &
      iter, ' updates', time(6)/iter, ' ns/update vectorized'
  end if
  deallocate (p_r, p_p, p_q, p_t, p_dt, p_EE, p_BB, p_ptr)
  call trace%end ()
END SUBROUTINE move_particles_v

!===============================================================================
!> Import particles from nbors export lists
!===============================================================================
SUBROUTINE import_particles (self)
  class(pic_particles_t)     :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: task, nbtask
  class(patch_t)   , pointer :: nbpatch
  class(pic_particle_t), pointer :: pa, next
  real(8)                    :: p(3)
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%import_particles')
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (pic_particles_t)
      nbtask => nbor%task
      select type (nbtask)
      class is (pic_particles_t)
        nbpatch => nbtask%patch
        !-----------------------------------------------------------------------
        ! Loop over particles in export list and grab those that belong here
        !-----------------------------------------------------------------------
        pa => nbtask%export_list%head
        do while (associated(pa))
          next => pa%next
          select type (pa)
          class is (pic_particle_t)
            p = nbpatch%llc_cart + (pa%r(:,pa%it) + pa%q(:,pa%it))*nbpatch%ds
          end select
          if (self%patch%contains(p)) then
            call self%particle_list%append (pa)
            call task%export_list%remove (pa)
          end if
          pa => next
        end do
      end select
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE import_particles

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_outside (self, patch)
  class(pic_particle_t):: self
  class(patch_t):: patch
  !.............................................................................
  is_outside = any(abs((self%r + self%q) - 0.5) > 0.5)
END FUNCTION is_outside

!===============================================================================
!> Test subroutine
!===============================================================================
SUBROUTINE scatter_fields2 (pa, EE, BB)
  class(pic_particle_t):: pa
  !real, dimension(:,:,:):: Ex, Ey, Ez, Bx, By, Bz
  real, dimension(3):: EE, BB
  call trace%begin ('pic_particles_t%scatter_fields2')
  EE = 0.0
  BB = 0.0
  call trace%end ()
END SUBROUTINE scatter_fields2

!===============================================================================
!> Initialize a doubly-linked-list (DLL), with anonymous carry
!===============================================================================
SUBROUTINE init_list (self)
  class(pic_list_t):: self
  !.............................................................................
  nullify (self%head, self%tail)
  !$omp atomic write
  self%n = 0
END SUBROUTINE init_list

!===============================================================================
!> Add an item at the tail
!===============================================================================
SUBROUTINE append (self, new)
  class(pic_list_t)              :: self
  class(pic_particle_t), pointer :: new
  !.............................................................................
  new%prev => self%tail
  nullify(new%next)
  if (.not.associated(self%head)) self%head => new
  if (associated(self%tail)) self%tail%next => new
  self%tail => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE append

!===============================================================================
!> Remove an item from the DLL; note that the item is not deallocated!
!===============================================================================
SUBROUTINE remove (self, old)
  class(pic_list_t):: self
  class(pic_particle_t), pointer:: old
  !.............................................................................
  if (associated(old%prev)) old%prev%next => old%next
  if (associated(old%next)) old%next%prev => old%prev
  if (associated(old, self%head)) then
    self%head => old%next
  end if
  if (associated(old, self%tail)) then
    self%tail => old%prev
  end if
  !$omp atomic
  self%n = self%n-1
END SUBROUTINE remove

!===============================================================================
!> Delete a double-linked list, including carried content
!===============================================================================
SUBROUTINE delete (self)
  class(pic_list_t):: self
  !.............................................................................
  class (pic_particle_t), pointer:: item, next
  !-----------------------------------------------------------------------------
  item => self%head
  do while (associated(item))
    next => item%next
    deallocate (item)
    item => next
  end do
  nullify (self%head)
  nullify (self%tail)
  self%n = 0
  call self%init
END SUBROUTINE delete

!===============================================================================
!> Test data
!===============================================================================
SUBROUTINE init_test (self)
  class(pic_particles_t)     :: self
  integer                    :: np, ip
  real                       :: r(3), dtmin
  class(pic_particle_t), pointer :: pa
  !-----------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%test')
  np = 16**3
  self%mass = 1.0
  self%charge = 0.0
  if (.not.associated(self%patch)) then
    allocate (self%patch)
    allocate (self%patch%link)
  end if
  call self%particle_list%init
  call self%export_list%init
  do ip=1,np
    allocate (pa)
    call self%particle_list%append (pa)
  end do
  self%patch%size = 1.0
  self%patch%n = 16
  self%patch%ds = self%patch%size/self%patch%n
  self%patch%position = 0.5
  self%patch%llc_cart = 0.0
  dtmin = 1.0
  ip = 0
  pa => self%particle_list%head
  do while (associated(pa))
    call random_number (pa%r(:,1))
    call random_number (pa%p(:))
    pa%q = 0
    pa%t = 0
    pa%it = 1
    call random_number (pa%dt)
    pa%dt = dt*10.**(-1.0+pa%dt)
    dtmin = min(dtmin,pa%dt)
    ip = ip+1
    pa => pa%next
  end do
  if (verbose > 1) &
    print *, 'np:', np, self%particle_list%n, dtmin
  call trace%end ()
END SUBROUTINE init_test

!===============================================================================
!> Test
!===============================================================================
SUBROUTINE test (self)
  class(pic_particles_t)         :: self
  integer                        :: np, ip
  real                           :: r(3), dtmin
  class(pic_particle_t), pointer :: pa
  !-----------------------------------------------------------------------------
  call trace%begin ('pic_particles_t%test')
  call self%particle_list%init
  !-----------------------------------------------------------------------------
  call init_test (self)
  call move_particles_nv (self)
  !-----------------------------------------------------------------------------
  call init_test (self)
  call move_particles_v (self)
  !-----------------------------------------------------------------------------
  call self%export_list%delete
  call self%particle_list%delete
  if (verbose > 0) &
    stop
  call trace%end ()
END SUBROUTINE test

END MODULE pic_particles_mod
