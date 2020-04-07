!===============================================================================
!> Data type for dust-gas interaction:
!>
!> The dust evolution is simulated using a fairly large number of particles oer
!> cell, e.g. of order 10-30 per size bin, or of the order of 100 in total
!>
!> The particles are moved using a particle solver, which operates on lists of
!> particles, which have attributes suitable to describe dust properties.  The
!> most important of these are
!>
!>    particle_t%bin  : a bin index specifying size and mass
!>    particle_t%w    : a weight that expresses number density per particle
!>    particle_t%q    : an index triple that giving the integer cell coordinates
!>    particle_t%r    : a fractional index giving the position in the cell
!>
!> The main procedures that contribute to evolving the dust are
!>
!>    dust_t%gather   : gathers distributions of particles inside each cell
!>    dust_t%evolve   : decides how to evolve the dust distribution
!>    dust_t%scatter  : updates the particle properties accordingly
!>
!> The dust evolution does not necessarily have to be done after each update of
!> particle positions; it could be done at regular intervals, semi-independently.
!>
!> The dust module could be called as a sub-process, or as a separate task.
!> Initially it is most likely the simplest and most transparent organization
!> to have the experiment_t level be the task, from the point of view of the
!> framework, and let it call the various sub-solvers at various points in time
!>
!> TODO:
!>   * Arrange to refer to patch velocity rather than momentum
!>   * Possibly correct drag for supersonic relative motion
!>   * Handle empty bins gracefully
!===============================================================================
!define DBG
MODULE dust_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE scaling_mod
  USE random_mod
  USE mesh_mod
  USE task_mod
  USE gpatch_mod
  USE list_mod
  USE particles_mod
  USE aux_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Dust data type
  !-----------------------------------------------------------------------------
  type, public:: dust_t
    real(8):: time = 0.0                      ! dust age
    real(8):: dtime = 1d0                     ! dust evolution cadence
    real:: solid_density                      ! mass per unit volume, code units
    real(8), allocatable:: m_cgs(:)           ! bin mass
    real(8), allocatable:: s_cgs(:)           ! particle size
    real(8), allocatable:: dv(:)              ! average speed difference
    real      , pointer, contiguous:: n (:,:,:,:) => null()  ! number density per bin
    integer(2), pointer, contiguous:: na(:,:,:,:) => null()  ! added per bin
    integer(2), pointer, contiguous:: np(:,:,:,:) => null()  ! number per bin
    real      , pointer, contiguous:: k (:,:)     => null()  ! kernel
    class(gpatch_t), pointer:: patch          ! pointer to the MHD patch
    type(particles_t), pointer:: particles    ! pointer to particles
    real, pointer, contiguous:: d(:,:,:)      ! density
    real, pointer, contiguous:: p(:,:,:,:)    ! momentum
    type(random_t):: random                   ! local random number generator
    logical, dimension(3):: lower             ! has lower boundary
    logical, dimension(3):: upper             ! has upper boundary
    type(aux_t):: aux                         ! auxiliary output
    real(8):: next_time                       ! next evolve time
    procedure(test_kernel), pointer:: kernel  ! kernel procedure
    class(list_t), pointer:: task_list        ! task list
    integer:: verbose
    integer:: n_bin
    logical:: detailed_timer = .false.
  contains
    procedure:: init
    procedure:: initial_condition
    procedure:: update
    procedure:: drag
    procedure:: output
    procedure:: check_particles
    procedure, private:: gather
    procedure, private:: evolve
    procedure, private:: scatter
    procedure, private:: test_kernel
    procedure, private:: renormalize
  end type
  integer:: verbose = 0                       ! information level
  integer:: cadence = 10                      ! particle steps per evolve step
  integer:: seed = 11                         ! random seed
  integer:: particles_per_bin = 10            ! particles per cell-bin
  real(8):: dust_time = 0.01                  ! dust evolution timestep
  real:: min_fraction = 0.7                   ! particles per bin (max decrease)
  real:: max_fraction = 1.4                   ! particles per bin (max increase)
  real:: d_min_cgs = 1e-15                    ! lowest density for renormalizing
  real:: dust_to_gas = 0.01                   ! dust-to-gas ratio  
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: solid_density_cgs = 2.0              ! g/cm**3
  real:: density_yr_cgs = 1e-12               ! volume density giving 1 yr stop
  real:: drag_coeff = 0.0                     ! constant in front of drag force
  real:: growth_rate = 1e-2                   ! per unit code time
  integer:: n_bin = 15                        ! 5 bins per power of ten
  logical:: particle_io = .false.             ! output particle positions?
  type(mpi_counter_t):: counter               ! counter returning uniq id
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self, particles)
  class(dust_t):: self
  type(particles_t), pointer:: particles
  !.............................................................................
  integer:: iostat, bin
  real:: dln_m
  logical:: first_time = .true.
  namelist /dust_params/ n_bin, cadence, seed, particles_per_bin, particle_io, &
    min_size_cgs, max_size_cgs, solid_density_cgs, density_yr_cgs, dust_time, &
    growth_rate, min_fraction, max_fraction, d_min_cgs, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%init')
  self%particles => particles
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_params, iostat=iostat)
    call io%namelist_warning ('dust_params', iostat)
    write (io%output, dust_params)
    if (io%master) write (io%namelist, dust_params)
    call counter%init (0_8)
    n_bin = max(1,n_bin)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  self%kernel => test_kernel
  self%n_bin = n_bin
  self%verbose = verbose
  !-----------------------------------------------------------------------------
  ! Scaling to code units (keep dust size and mass in cgs)
  !-----------------------------------------------------------------------------
  self%solid_density = solid_density_cgs/(scaling%m/scaling%l**3)
  allocate (self%s_cgs(n_bin), self%m_cgs(n_bin))
  dln_m = log(max_size_cgs/min_size_cgs)/n_bin
  do bin=1,n_bin
    self%s_cgs(bin) = exp(log(min_size_cgs) + (bin-0.5)*dln_m)
    self%m_cgs(bin) = self%solid_density*self%s_cgs(bin)**3
    write(io_unit%log,'(a,i4,1p,2e11.2)') &
      'dust_t%init: bin, s_cgs, m_cgs =', bin, self%s_cgs(bin), self%m_cgs(bin)
  end do
  !-----------------------------------------------------------------------------
  ! Pre-factor for drag force.  The input parameter density_yr_cgs is the
  ! density (in CGS) for which the stopping time in 1 yr for a 1 cm particle.
  ! Since we keep size in cm, we just need to multiply with the density scaling
  ! constant, and divide with the value of one year in code units.
  !-----------------------------------------------------------------------------
  drag_coeff = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize an initial particle distribution, with weights equal to the local
!> mass density in code units (this must be taken into account in the kernel)
!===============================================================================
SUBROUTINE initial_condition (self, l, u)
  class(dust_t):: self
  integer:: l(3), u(3)
  !.............................................................................
  class(gpatch_t), pointer:: patch
  class(particle_t), pointer:: p
  real, pointer, contiguous:: d(:,:,:)
  integer:: i1, i2, i3, bin, j
  real:: c
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%initial_condition', itimer=itimer)
  call self%random%init
  patch => self%patch
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    do bin=1,n_bin
      do j=1,particles_per_bin
        allocate (p)
        p%bin = bin
        p%id = counter%update (1_8)
        p%q = [i1,i2,i3]
        p%r = self%random%uniform (3)
        p%v = 0.0
        p%c = drag_coeff/self%s_cgs(bin)
        call self%particles%append (p)
        p%w = d(i1,i2,i3)*dust_to_gas
        if (verbose > 2) &
          write (io_unit%log,'(a,i6,2x,3i3,1x,3f5.2,1x,1p,2e11.2)') &
          'particle_t%IC: id, q, r, w, c =', p%id, p%q, p%r, p%w, p%c
      end do
    end do
  end do
  end do
  end do
  !.............................................................................
  if (verbose > 1) &
    write (io_unit%log,'(a,i6,1p,2e11.2)') &
    'particle_t%IC: id, n =', patch%id, real(counter%update(0)), &
    real(product(u-l+1)*n_bin*particles_per_bin)
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE initial_condition

!===============================================================================
!> Update the dust state
!===============================================================================
SUBROUTINE update (self, particles)
  class(dust_t):: self
  type(particles_t), pointer:: particles
  logical:: do_dust
  integer:: n(3), l(3), u(3)
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('dust_t%update', itimer=itimer)
  self%particles => particles
  !-----------------------------------------------------------------------------
  ! Check if any cell is running low on macro-particles
  !-----------------------------------------------------------------------------
  do_dust = cadence > 0 .and. n_bin > 1
  if (do_dust) then
    do_dust = mod(self%patch%istep,cadence)==0
  end if
  !...........................................................................
  if (verbose > 0) then
    write (io_unit%log,'(a,i6,1p,3e11.2,l3)') &
      'dust_t%update: id, patch%time, time, dtime =', &
      self%patch%id, self%patch%time, self%time, self%dtime, do_dust
  end if
  !...........................................................................
  n = self%patch%mesh%gn                        ! cells in patch
  l = self%patch%mesh%li                        ! lower index limit
  u = self%patch%mesh%ui                        ! upper index limit
  allocate (self%n (n(1),n(2),n(3),n_bin))      ! work array for histograms
  allocate (self%na(n(1),n(2),n(3),n_bin))      ! number of added particles
  allocate (self%dv(n_bin))                     ! average speed difference
  call self%gather
  call self%renormalize
  call self%output
  if (do_dust) then
    allocate (self%k(n_bin,n_bin))              ! kernel expressing exchange
    call self%evolve
    call self%scatter
    deallocate (self%k)
  end if
  deallocate (self%n, self%na, self%dv)
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Gather information about the current state of the dust and gas
!===============================================================================
SUBROUTINE gather (self)
  class(dust_t):: self
  !.............................................................................
  integer:: bin, it, i(3), n(3), np, l(3), u(3)
  type(particle_t), pointer:: p
  real(8):: dv(n_bin)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%gather', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Compute histograms of particle number densities in each cell
  !-----------------------------------------------------------------------------
  self%n = 0.0                                  ! prepare to accumulate
  self%na = 0                                   ! initialize
  np = 0                                        ! number of updates
  dv = 0d0                                      ! average speed difference
  if (self%lower(3)) l(3)=self%patch%mesh(3)%lb ! adjust at boundaries
  if (self%upper(3)) u(3)=self%patch%mesh(3)%ub ! adjust at boundaries
  p => self%particles%head                      ! first particle
  !.............................................................................
  if (verbose > 0) &
    write (io_unit%log,*) 'GATHER'
  !.............................................................................
  do while (associated(p))                      ! loop over particle
    bin = p%bin                                 ! particle bin
    i = p%q                                     ! particle space index
    self%n(i(1),i(2),i(3),bin) = &              ! accumulate
    self%n(i(1),i(2),i(3),bin) + p%w            ! into self%n
    np = np+1
    !...........................................................................
    if (verbose > 2) then
      write(io_unit%log,'(a,i6,2i8,i3,1x,3i3,1p,2e11.2)') &
        'gather: patch%id, pa%id, np, bin, i, w, n =', &
        self%patch%id, p%id, np, bin, i, p%w, self%n(i(1),i(2),i(3),bin)
    end if
    !...........................................................................
    dv(bin) = dv(bin) + sqrt(sum(p%dv**2))      ! sum up speed differences
    p => p%next                                 ! next particle
  end do
  self%dv = dv/np                               ! average speed difference
  !.............................................................................
  if (verbose > 1) then
    write (io_unit%log,'(a,i6,1p,10e11.2)') &
    'gather: id, dv =', self%patch%id, dv
  end if
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE gather

!===============================================================================
!> Based on the current state of gas and dust, decide how to evolve it, by
!> applying a kernel to the distribution over bins and adding up the results.
!> The coefficient to apply is the ratio of the old to the new number density.
!===============================================================================
SUBROUTINE evolve (self)
  class(dust_t):: self
  !.............................................................................
  real:: n(n_bin), delta_n
  integer:: i1, i2, i3, i, j
  class(mesh_t), pointer:: mesh(:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%evolve', itimer=itimer)
  if (self%particles%n > 0) then
    mesh  => self%patch%mesh
    do i3=mesh(3)%li,mesh(3)%ui
    do i2=mesh(2)%li,mesh(2)%ui
    do i1=mesh(1)%li,mesh(1)%ui
      call self%kernel
      n(:) = self%n(i1,i2,i3,:)
      !-------------------------------------------------------------------------
      ! In order to conserve mass the number density added to bin j must be
      ! normalized by the ratio of bin masses.
      !-------------------------------------------------------------------------
      do i=1,n_bin
        do j=i+1,n_bin
          delta_n = self%n(i1,i2,i3,i)*self%k(i,j)*abs(self%dv(j)-self%dv(i))
          n(i) = n(i) - delta_n
          n(j) = n(j) + delta_n*self%m_cgs(i)/self%m_cgs(j)
        end do
      end do
      !-------------------------------------------------------------------------
      ! Guard against negative densities.  Since sum(k(i,:)) should be less than
      ! one, it should not be possible to come out with a negative n, so the
      ! test is a sanity test
      !-------------------------------------------------------------------------
      if (minval(n) <= 0.0) then
        !.......................................................................
        if (verbose > 0) then
          write (io_unit%log,'(a,3i5,1p,100e10.2)') &
            'evolve WARNING: i1, i2, i3, n =', i1, i2, i3, n
          write (stderr,'(a,3i5,1p,100e10.2)') &
            'evolve WARNING: i1, i2, i3, n =', i1, i2, i3, n
          flush (stderr)
        end if
        !.......................................................................
        n = max(0.0,n)
      end if
      if (verbose > 2) then
        write(io_unit%log,'(a,3i5,1p,100e10.2)') &
        'evolve: i1,i2,i3,n =', i1, i2, i3, n
      end if
      !-------------------------------------------------------------------------
      ! The final coefficient of change is the new n divided by the old n.  If
      ! the bin would appear to become empty, maintain current weights.
      !-------------------------------------------------------------------------
      do i=1,n_bin
        if (n(i) == 0.0 .or. self%n(i1,i2,i3,i) == 0.0) then
          self%n(i1,i2,i3,i) = 1.0
        else
          self%n(i1,i2,i3,i) = n(i)/self%n(i1,i2,i3,i)
        end if
      end do
    end do
    end do
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE evolve

!===============================================================================
!> Construct a kernel k(i,j), where for j>i the k(i,j) element is the fraction
!> of bin i given or taken to/from slot j; therefore, sum(k(i,:)) must be < 1
!===============================================================================
SUBROUTINE test_kernel (self)
  class(dust_t):: self
  !.............................................................................
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (verbose > 0) call trace%begin ('dust_t%test_kernel', itimer=itimer)
  self%k = 0.0
  do i=1,n_bin-1
    !---------------------------------------------------------------------------
    ! Simple kernel for test, assuming bin populations move to the next bin
    ! at a fixed rate per unit code time
    !---------------------------------------------------------------------------
    self%k(i,i+1) = growth_rate*self%dtime
  end do
  if (verbose > 0) call trace%end (itimer)
END SUBROUTINE test_kernel

!===============================================================================
!> Update the dust state in the task, based on the desired evolution
!===============================================================================
SUBROUTINE scatter (self)
  class(dust_t):: self
  !.............................................................................
  integer:: bin, it, i(3), n(3)
  type(particle_t), pointer:: pa
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%scatter', itimer=itimer)
  pa => self%particles%head                     ! first particle
  do while (associated(pa))                     ! loop over particle
    bin = pa%bin                                ! particle bin
    i = pa%q                                    ! particle space index
    pa%w = pa%w*self%n(i(1),i(2),i(2),bin)      ! renormalize
    pa => pa%next                               ! next particle
  end do
  call trace%end (itimer)
END SUBROUTINE scatter

!===============================================================================
!> Compute drag force on a particle, and add opposite force to the gas.  Since
!> velocity is measured in cells per unit time, the value is divided by the cell
!> size.
!===============================================================================
SUBROUTINE drag (self, p)
  class(dust_t):: self
  type(particle_t):: p
  integer:: i(3), bin
  real:: d, s, v(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%detailed_timer) &
    call trace%begin ('particle_solver_t%drag',itimer=itimer)
  i = p%q
  bin = p%bin
  s = self%s_cgs(bin)
  i = min(max(i,1),self%patch%gn)
  d = self%d(i(1),i(2),i(3))
  v = self%p(i(1),i(2),i(3),:)/d 
  p%dv = v - p%v
  p%c = drag_coeff*d/s            !*sqrt((speed2+sound2)/sound2)
  p%d = p%c*p%dv
!  self%patch%force_per_unit_volume(i(1),i(2),i(3),:) = &
!  self%patch%force_per_unit_volume(i(1),i(2),i(3),:) - p%d*p%w
  !.............................................................................
  if (verbose > 4) then
    write (io_unit%log,'(a,i8,1p,3e9.2,0p,3i3,1x,3f6.2,4(1x,1p,3e9.2))') &
      'drag: id, s, d, c, q, r, v, dv, d =', &
      p%id, s, d, p%c, p%q, p%r, p%v, p%dv, p%d
  end if
  !.............................................................................
  if (self%detailed_timer) call trace%end (itimer)
END SUBROUTINE drag

!===============================================================================
!> Add or remove macro-particles, to maintain good statistics
!===============================================================================
SUBROUTINE renormalize (self)
  class(dust_t):: self
  !.............................................................................
  integer:: i(3), l(3), u(3), bin, i1, i2, i3, count, unit
  class(particle_t), pointer:: p, next
  !integer(kind=2), allocatable:: n(:,:,:,:)
  logical:: add, rem, alloc
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%renormalize', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check if they are within limits, and adjust if not
  !-----------------------------------------------------------------------------
  rem = any(self%np > particles_per_bin*max_fraction)
  add = any(self%np < particles_per_bin*min_fraction)
  if (add .or. rem) then
    l = self%patch%mesh%li                        ! lower index limit
    u = self%patch%mesh%ui                        ! upper index limit
    !---------------------------------------------------------------------------
    if (verbose > 0) then
      unit = io_unit%log
      call print_before
      if (add) call add_particles
      if (rem) call remove_particles
      call print_after
    else
      if (add) call add_particles
      if (rem) call remove_particles
    end if
  end if
  !.............................................................................
  call trace%end (itimer)
contains
  !=============================================================================
  ! The particle add procedure it triggered when a cell (any cell and bin)
  ! has less than particles_per_bin*min_fraction particles, but while running
  ! brings the number of particles per cell up to at most particles_per_bin in
  ! all cells, provided the dust density is at least d_min_cgs in CGS units. 
  ! Afterwards the particle weight in the cells are renormalized, so mass is
  ! conserved.
  !=============================================================================
  subroutine add_particles
    integer:: np, n, na
    real:: w, c, d_cgs
    integer, save:: itimer1=0
    !---------------------------------------------------------------------------
    call trace%begin ('dust_t%add_particles', itimer=itimer1)
    na = 0
    p => self%particles%head
    do while (associated(p))
      next => p%next
      bin = p%bin                                 ! particle bin
      i = p%q                                     ! particle space index
      if (all(i >= l .and. i <= u)) then
        !-----------------------------------------------------------------------
        d_cgs = self%n(i(1),i(2),i(3),bin)*scaling%d
        if (d_cgs > d_min_cgs) then
          np = self%np(i(1),i(2),i(3),bin) + self%na(i(1),i(2),i(3),bin)
          !.....................................................................
          if (verbose > 1) then
            if (i(1)==3 .and. i(2)==3 .and. bin==3) &
            write (io_unit%log,'(a,4i4,2x,3i6,1p,2e11.2,2(2x,3i4))') &
            'add_part',i,bin,self%np(i(1),i(2),i(3),bin), &
            self%na(i(1),i(2),i(3),bin), &
            particles_per_bin, d_cgs, d_min_cgs, l, u
          endif
          !.....................................................................
          if (np < particles_per_bin) then
            c = drag_coeff/self%s_cgs(bin)
            do while (np < particles_per_bin)
              allocate (p)
              p%id = counter%update (1_8)
              p%q = i
              p%bin = bin
              p%r = self%random%uniform (3)
              p%v = 0.0
              p%c = c
              p%w = 0.0
              call self%particles%append (p)
              self%na(i(1),i(2),i(3),bin) = self%na(i(1),i(2),i(3),bin) + 1
              np = self%np(i(1),i(2),i(3),bin) + self%na(i(1),i(2),i(3),bin)
              na = na + 1
            end do
            !...................................................................
            if (verbose > 1) then
              write (io_unit%log,*) 'ADDED', n, np, particles_per_bin, &
                self%na(i(1),i(2),i(3),bin), self%np(i(1),i(2),i(3),bin)
            end if
            !...................................................................
          end if
        end if
      end if
      p => next
    end do
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,*) 'add_particles: added', na
    end if
    !...........................................................................
    call normalize_weights
    call trace%end (itimer1)
  end subroutine add_particles
  !=============================================================================
  ! The particle removal procedure it triggered when a cell (any cell and bin)
  ! has more than particles_per_bin*max_fraction particles, but while running
  ! brings the number of particles per cell down to at most particles_per_bin
  ! per cell and bin.  Afterwards the particle weight in the cells are
  ! renormalized, so mass is conserved.
  !=============================================================================
  subroutine remove_particles
    integer:: np, nr
    integer, save:: itimer2=0
    !---------------------------------------------------------------------------
    call trace%begin ('dust_t%remove_particles', itimer=itimer2)
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,*) 'REMOVE'
    endif
    !...........................................................................
    nr = 0
    p => self%particles%head
    do while (associated(p))
      next => p%next
      bin = p%bin                                 ! particle bin
      i = p%q                                     ! particle space index
      if (all(i >= l .and. i <= u)) then
        np = self%np(i(1),i(2),i(3),bin) + self%na(i(1),i(2),i(3),bin) 
        if (np > particles_per_bin) then
          !.....................................................................
          if (verbose > 0) then
            write (io_unit%log,*) 'remove', p%id, p%q, p%bin
          endif
          !.....................................................................
          call self%particles%remove (p)
          deallocate (p)
          self%na(i(1),i(2),i(3),bin) = self%na(i(1),i(2),i(3),bin) - 1
          nr = nr + 1
        end if
      end if
      p => next
    end do
    if (verbose > 0) then
      write (io_unit%log,*) 'remove_particles: removed', nr
    end if
    call trace%end (itimer2)
    call normalize_weights
  end subroutine remove_particles
  !=============================================================================
  !> Renormalize weights, so the new set of particles in a cell carry the same
  !> mass as the previous set did.  The total weight is saved in self%n(), 
  !> and the number of macro-particles is np = self%np()+self%na(), so the
  !> new particle weight should be set to self%n()/np.
  !=============================================================================
  subroutine normalize_weights
    integer:: np
    integer, save:: itimer3=0
    !---------------------------------------------------------------------------
    call trace%begin ('dust_t%normalize_weights', itimer=itimer3)
    !---------------------------------------------------------------------------
    ! Renormalize weights to get the same total as before
    !---------------------------------------------------------------------------
    p => self%particles%head
    do while (associated(p))
      i= p%q
      bin = p%bin
      np = self%np(i(1),i(2),i(3),bin) + self%na(i(1),i(2),i(3),bin)
      if (np > 0) then
        p%w = self%n(i(1),i(2),i(3),bin)/np
      else
        p%w = 0.0
      end if
      p => p%next
    end do
    call trace%end (itimer3)
  end subroutine normalize_weights
  !=============================================================================
  ! Compute the sum of particle weights and the number of particles per cell-bin
  !=============================================================================
  subroutine get_densities
    integer:: i(3), bin
    integer, save:: itimer4=0
    !---------------------------------------------------------------------------
    call trace%begin ('dust_t%get_densities', itimer=itimer4)
    !---------------------------------------------------------------------------
    ! Compute total weights of new particle set
    !---------------------------------------------------------------------------
    self%n  = 0.0
    self%na = 0.0
    self%np = 0
    p => self%particles%head
    do while (associated(p))
      i = p%q                                     ! particle space index
      bin = p%bin                                 ! particle bin
      self%n(i(1),i(2),i(3),bin) = &
      self%n(i(1),i(2),i(3),bin) + p%w
      self%np(i(1),i(2),i(3),bin) = &
      self%np(i(1),i(2),i(3),bin) + 1
      p => p%next
    end do
    call trace%end (itimer4)
  end subroutine get_densities
  !=============================================================================
  subroutine print_before
    do bin=1,n_bin
      write(unit,'(a,i3,2x,20i4)') 'particles before:', bin, self%np(3,3,:,bin)
    end do
    do bin=1,n_bin
      write(unit,'(a,i3,2x,1p,10e11.2)') 'bins before:', bin, self%n(3,3,1:10,bin)
    end do
    if (verbose > 1) then
      call get_densities
      do bin=1,n_bin
        write(unit,'(a,i3,2x,1p,10e11.2)') 'bins redone:', bin, self%n(3,3,1:10,bin)
      end do
    end if
    call normalize_weights
  end subroutine
  !=============================================================================
  subroutine print_after
    do bin=1,n_bin
      write(unit,'(a,i3,2x,20i4)') 'added particles:', bin, self%na(3,3,:,bin)
    end do
    if (verbose > 1) then
      call get_densities
      do bin=1,n_bin
        write(unit,'(a,i3,2x,1p,10e11.2)') 'bins  after:', bin, self%n(3,3,1:10,bin)
      end do
    end if
    do bin=1,n_bin
      write(unit,'(a,i3,2x,20i4)') 'particles  after:', bin, &
        self%np(3,3,:,bin) + self%na(3,3,:,bin)
    end do
  end subroutine
END SUBROUTINE renormalize

!===============================================================================
!> Output dust densities and (optionally) dust particles
!===============================================================================
SUBROUTINE output (self)
  class(dust_t):: self
  !.............................................................................
  class(particle_t), pointer:: p
  class(mesh_t), pointer:: mesh(:)
  real, pointer, contiguous:: tmp(:,:)
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%output')
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_d', self%n)
    !...........................................................................
    if (verbose > 1) then
      do i=1,n_bin
        write (io_unit%log, '(a,i3,2x,1p,10e11.2)') &
          'bins aux IO:', i, self%n(3,3,1:10,i)
      end do
    end if
    !...........................................................................
    if (particle_io) then
      allocate (tmp(3,self%particles%n))
      i = 1
      mesh => self%patch%mesh
      p => self%particles%head
      do while (associated(p))
        tmp(:,i) = mesh%p + (p%q + p%r - mesh%o)*mesh%d
        i = i+1
        p => p%next
      end do
      call self%aux%register ('dust_p', tmp)
      call self%aux%output_sequence (self%patch%id, self%time)
      deallocate (tmp)
    else
      call self%aux%output_sequence (self%patch%id, self%time)
    end if
  end if
  call trace%end ()
END SUBROUTINE output

SUBROUTINE check_particles (self, label)
  class(dust_t):: self
  character(len=*):: label
#ifdef DBG
  class(particle_t), pointer:: p
  write (io_unit%log,*) 'check_particles: '//label
  p => self%particles%head
  do while (associated(p))
    if (any(p%q <= 0)) write (io_unit%log,*) 'ERROR: '//label, p%id
    p => p%next
  end do
#endif
END SUBROUTINE

END MODULE dust_mod
