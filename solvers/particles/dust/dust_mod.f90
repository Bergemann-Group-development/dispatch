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
  USE units_mod
  USE random_mod
  USE mesh_mod
  USE task_mod
  USE gpatch_mod
  USE list_mod
  USE particles_mod
  USE aux_mod
  USE omp_timer_mod
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
    real, pointer, contiguous:: s_cgs(:)      ! particle size
    real(8), allocatable:: dv(:)              ! average speed difference
    real      , pointer, contiguous:: dd(:,:,:,:) => null()  ! number density per bin
    integer(2), pointer, contiguous:: na(:,:,:,:) => null()  ! added per bin
    integer(2), pointer, contiguous:: np(:,:,:,:) => null()  ! number per bin
    real      , pointer, contiguous:: k (:,:)     => null()  ! kernel
    real      , pointer, contiguous:: vg(:,:,:,:) => null()  ! gas velocity
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
    integer, dimension(3):: lb=0, ub=0        ! lower and upper index bounds
    integer, dimension(5):: i_dbg=0           ! debug selection
    integer:: verbose
    integer:: n_bin, pp_bin
    integer:: id_dbg=0, ip_dbg=0
    logical:: do_sort = .false.
    logical:: detailed_timer = .false.
    real:: drag_coeff=0.0                     ! constant in front of drag force
    real:: n_min, n_max
  contains
    procedure:: init
    procedure:: init_task_list
    procedure:: initial_condition
    procedure:: update
    procedure:: output
    procedure, private:: evolve
    procedure, private:: scatter
    procedure, private:: test_kernel
    procedure:: consistency
    procedure:: consistency2
  end type
  integer:: verbose = 0                       ! information level
  integer:: cadence = 10                      ! particle steps per evolve step
  integer:: seed = 11                         ! random seed
  integer:: particles_per_bin = 10            ! particles per cell-bin
  integer:: output_every=1000                 ! output only that fraction
  real(8):: dust_time = 0.01                  ! dust evolution timestep
  real:: min_fraction = 0.7                   ! particles per bin (max decrease)
  real:: max_fraction = 1.4                   ! particles per bin (max increase)
  real:: dust_to_gas = 0.01                   ! dust-to-gas ratio  
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: d_min_cgs = 1e-24                    ! lowest density for renormalizing
  real:: solid_density_cgs = 2.0              ! g/cm**3
  real:: density_yr_cgs = 1e-12               ! volume density giving 1 yr stop
  real:: drag_coeff = 0.0                     ! constant in front of drag force
  real:: growth_rate = 1e-2                   ! per unit code time
  integer:: n_bin = 15                        ! 5 bins per power of ten
  logical:: particle_io = .false.             ! output particle positions?
  logical:: validate = .false.                ! check consistency
  public:: n_bin
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self, particles)
  class(dust_t):: self
  type(particles_t), pointer:: particles
  !.............................................................................
  integer:: iostat, bin, m(3)
  real:: dln_m
  logical:: first_time = .true.
  namelist /dust_params/ n_bin, cadence, seed, particles_per_bin, particle_io, &
    min_size_cgs, max_size_cgs, solid_density_cgs, density_yr_cgs, dust_time, &
    growth_rate, min_fraction, max_fraction, d_min_cgs, output_every, validate, &
    verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%init')
  self%particles => particles
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_params, iostat=iostat)
    call io%namelist_warning ('dust_params', iostat)
    write (io%namelist, dust_params)
    n_bin = max(1,n_bin)
    first_time = .false.
  end if
  call io%lock%unset
  self%kernel => test_kernel
  self%n_bin = n_bin
  self%pp_bin = particles_per_bin
  self%n_min = particles_per_bin*min_fraction
  self%n_max = particles_per_bin*max_fraction
  self%verbose = verbose
  call self%aux%init
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
  m = self%patch%gn
  allocate (self%dd(m(1),m(2),m(3),n_bin))
  !-----------------------------------------------------------------------------
  ! Pre-factor for drag force.  The input parameter density_yr_cgs is the
  ! density (in CGS) for which the stopping time in 1 yr for a 1 cm particle.
  ! Since we keep size in cm, we just need to multiply with the density scaling
  ! constant, and divide with the value of one year in code units.
  !-----------------------------------------------------------------------------
  drag_coeff = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  self%drag_coeff = drag_coeff
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after nbor lists are created, making it possible
!> to create a copy of the patch nbor list, with tasks pointing to self
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(dust_t):: self
  class(list_t), pointer:: task_list
  type(particle_t):: particle
  integer:: bytes
  real:: gb
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%init_task_list')
  if (self%patch%id == 1) then
    bytes = storage_size(particle)/8
    write (stdout,'(" bytes per particle:",i5)') bytes
    gb = io%ntotal*product(self%patch%n) &
      *real(particles_per_bin*n_bin)*bytes/1024.**3
    write (stdout,'(" estimated memory for particles:",f5.1," --", f5.1," GB")') &
      gb, max_fraction*gb
    gb = io%ntotal*product(shape(self%patch%mem))*4./1024**3
    write (stdout,'(" estimated memory for cells:",f6.2," GB")') gb
  end if
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Initialize an initial particle distribution, with weights equal to the local
!> mass density in code units (this must be taken into account in the kernel)
!===============================================================================
SUBROUTINE initial_condition (self)
  class(dust_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  class(particle_t), pointer:: p
  real, pointer, contiguous:: d(:,:,:)
  integer:: i1, i2, i3, bin, j
  integer:: l(3), u(3), i
  real:: c, w, vg(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%initial_condition', itimer=itimer)
  call self%random%init
  patch => self%patch
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  l = self%lb
  u = self%ub
  i = 0
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    w = d(i1,i2,i3)*dust_to_gas/(n_bin*particles_per_bin)
    vg = patch%mem(i1,i2,i3,patch%idx%px:patch%idx%pz,patch%it,1) &
       / (patch%ds*patch%mem(i1,i2,i3,patch%idx%d,patch%it,1))
    do bin=1,n_bin
      self%dd(i1,i2,i3,bin) = 0.0
      do j=1,particles_per_bin
        i = i+1
        allocate (p)
        call p%init
        if (mod(i,output_every)==0) then
          p%flag = 1   ! select for tracking
        else
          p%flag = 0
        end if
        p%bin = bin
        p%q = [i1,i2,i3]
        p%r = self%random%uniform (3)
        p%v = vg
        p%w = w
        self%dd(i1,i2,i3,bin) = self%dd(i1,i2,i3,bin) + w
        call self%particles%append (p)
        !.......................................................................
        if (verbose > 2 .or. p%id == self%ip_dbg) &
          write (io_unit%log,'(a,i8,2x,3i3,1x,3f5.2,1x,1p,2e11.2)') &
          'particle_t%IC: id, q, r, w, c =', p%id, p%q, p%r, p%w
      end do
    end do
  end do
  end do
  end do
  !.............................................................................
  if (verbose > 1) &
    write (io_unit%log,'(a,i8,1p,2e11.2)') &
    'particle_t%IC: id, n =', patch%id, &
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
  integer:: n(3), i
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
  !allocate (self%na(n(1),n(2),n(3),n_bin))      ! number of added particles
  !allocate (self%dv(n_bin))                     ! average speed difference
  !allocate (self%np(n(1),n(2),n(3),n_bin))      ! work array for pp-cell-bin
  !call self%gather
  !call self%renormalize
  call self%output
  if (do_dust) then
    allocate (self%k(n_bin,n_bin))              ! kernel expressing exchange
    call self%evolve
    call self%scatter
    deallocate (self%k)
  end if
  !deallocate (self%na, self%dv, self%np)
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Based on the current state of gas and dust, decide how to evolve it, by
!> applying a kernel to the distribution over bins and adding up the results.
!> The coefficient to apply is the ratio of the old to the new number density.
!===============================================================================
SUBROUTINE evolve (self)
  class(dust_t):: self
  !.............................................................................
  real:: n(n_bin), delta_n
  integer:: i1, i2, i3, i, j, l(3), u(3)
  class(mesh_t), pointer:: mesh(:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%evolve', itimer=itimer)
  if (self%particles%n > 0) then
    l = self%lb
    u = self%ub
    mesh  => self%patch%mesh
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      call self%kernel
      n(:) = self%dd(i1,i2,i3,:)
      !-------------------------------------------------------------------------
      ! In order to conserve mass the number density added to bin j must be
      ! normalized by the ratio of bin masses.
      !-------------------------------------------------------------------------
      do i=1,n_bin
        do j=i+1,n_bin
          delta_n = self%dd(i1,i2,i3,i)*self%k(i,j)*abs(self%dv(j)-self%dv(i))
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
        if (n(i) == 0.0 .or. self%dd(i1,i2,i3,i) == 0.0) then
          self%dd(i1,i2,i3,i) = 1.0
        else
          self%dd(i1,i2,i3,i) = n(i)/self%dd(i1,i2,i3,i)
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
    pa%w = pa%w*self%dd(i(1),i(2),i(2),bin)      ! renormalize
    pa => pa%next                               ! next particle
  end do
  call trace%end (itimer)
END SUBROUTINE scatter

!===============================================================================
!> Output dust densities and (optionally) dust particles
!===============================================================================
SUBROUTINE output (self)
  class(dust_t):: self
  !.............................................................................
  class(particle_t), pointer:: p
  class(mesh_t), pointer:: mesh(:)
  real, pointer, contiguous:: tmp_r(:,:)
  integer, pointer, contiguous:: tmp_i(:)
  integer:: i, j, n
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%output')
  if (verbose > 1) then
    write (io_unit%log, '(a,i6,2f8.4)') &
      'dust_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_s', self%s_cgs)
    call self%aux%register ('dust_d', self%dd)
    !...........................................................................
    if (verbose > 1) then
      do i=1,n_bin
        write (io_unit%log, '(a,i6, i3,2x,1p,10e11.2)') &
          'dust_t%output: id, bin, n =', self%patch%id, i, self%dd(3,3,1:10,i)
      end do
    end if
    !...........................................................................
    if (particle_io) then
      mesh => self%patch%mesh
      n = 0
      p => self%particles%head
      do while (associated(p))
        if (p%flag > 0) then
          n = n+1
        end if
        p => p%next
      end do
      allocate (tmp_r(5,n))
      allocate (tmp_i(  n))
      i = 0
      p => self%particles%head
      do while (associated(p))
        if (p%flag > 0) then
          i = i+1
          tmp_r(1:3,i) = mesh%p + (p%q + p%r - mesh%o)*mesh%d
          tmp_r(  4,i) = p%w
          tmp_r(  5,i) = self%s_cgs(p%bin)
          tmp_i(    i) = p%id
        end if
        p => p%next
      end do
      call self%aux%register ('dust_p' , tmp_r)
      call self%aux%register ('dust_id', tmp_i)
      call self%aux%output_sequence (self%patch%id, self%time)
      if (verbose > 0) then
        write (io_unit%log, '(a,i6,f8.4,i3,2i8)') &
         'dust_t%output: id, particle array =', &
          self%patch%id, self%time, shape(tmp_r), i
      end if
      deallocate (tmp_r, tmp_i)
    else
      call self%aux%output_sequence (self%patch%id, self%time)
    end if
    if (verbose > 1) then
      write (io_unit%log, '(a,i6,1p,10e11.2)') &
       'dust_t%output: self%aux%aux_next =', &
        self%patch%id, self%aux%aux_next
    end if
  end if
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!===============================================================================
SUBROUTINE consistency (self, label)
  class(dust_t):: self
  character(len=*):: label
  class(sorted_t), pointer:: s
  class(particle_t), pointer:: p
  integer:: n, l(3), u(3), j
  integer(kind=2):: i1, i2, i3, bin
  !-----------------------------------------------------------------------------
  if (.not.validate) return
  print *,label
  l = self%patch%mesh%lb
  u = self%patch%mesh%ub
  p => self%particles%head
  do while (associated(p))
    if (any(p%q < l .or. p%q > u)) then
      print *,'outside:', p%id,p%q,p%bin
    end if
    p => p%next
  end do
  n = 0
  do bin=1,n_bin
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        call self%consistency2 (label, i1, i2, i3, bin)
        s => self%particles%sorted(i1,i2,i3,bin)
        n = n + s%n
      end do
    end do
  end do
  end do
  if (n /= self%particles%n) then
    write (stdout,2) 'mismatch in total n:', n, self%particles%n,i1,i2,i3,label
    2 format(a,2i9,2x,3i4,2x,a)
  end if
END SUBROUTINE consistency

!===============================================================================
!> Check that the head to tail range is consistent with the cell count
!===============================================================================
SUBROUTINE consistency2 (self, label, i1, i2, i3, bin)
  class(dust_t):: self
  character(len=*):: label
  integer(kind=2):: i1, i2, i3, bin
  class(sorted_t), pointer:: s
  class(particle_t), pointer:: p
  integer:: n, j
  !-----------------------------------------------------------------------------
  if (.not.validate) return
  s => self%particles%sorted(i1,i2,i3,bin)
  if (s%n > 0) then
    p => s%head
    do j=1,s%n
      if (j < s%n) then
        if (associated(p,s%tail)) then
          write (stdout,1) 'mismatch1 in local n:',j,s%n,i1,i2,i3,label
          write (io_unit%log,1) 'mismatch1 in local n:',j,s%n,i1,i2,i3,label
        1 format(a,2i4,2x,3i4,2x,a)
          exit
        end if
      else if (.not. associated(p,s%tail)) then
        write (stdout,1) 'mismatch2 in local n:',j,s%n,i1,i2,i3,label
        write (io_unit%log,1) 'mismatch2 in local n:',j,s%n,i1,i2,i3,label
        exit
      end if
      p => p%next
    end do
  end if
END SUBROUTINE consistency2

END MODULE dust_mod
