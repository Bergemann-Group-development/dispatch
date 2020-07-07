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
    real(8), allocatable:: dv(:)              ! average speed difference
    real      , pointer, contiguous:: dd(:,:,:,:) => null()  ! number density per bin
    integer(2), pointer, contiguous:: na(:,:,:,:) => null()  ! added per bin
    integer(2), pointer, contiguous:: np(:,:,:,:) => null()  ! number per bin
    real      , pointer, contiguous:: k (:,:)     => null()  ! kernel
    real      , pointer, contiguous:: vg(:,:,:,:) => null()  ! gas velocity
    class(gpatch_t), pointer:: patch          ! pointer to the MHD patch
    type(particles_t), pointer:: particles(:)
    type(particles_t), pointer:: sorted(:,:,:,:)
    real, pointer, contiguous:: m_cgs(:)      ! particle size
    real, pointer, contiguous:: d(:,:,:)      ! density
    real, pointer            :: p(:,:,:,:)    ! momentum
    logical, dimension(3):: lower             ! has lower boundary
    logical, dimension(3):: upper             ! has upper boundary
    type(aux_t):: aux                         ! auxiliary output
    real(8):: next_time                       ! next evolve time
    procedure(test_kernel), pointer:: kernel  ! kernel procedure
    class(list_t), pointer:: task_list        ! task list
    integer, dimension(3):: lb=0, ub=0        ! lower and upper index bounds
    integer, dimension(5):: i_dbg=0           ! debug selection
    integer:: verbose
    integer:: id_dbg=0, ip_dbg=0
    logical:: on = .false.
    logical:: detailed_timer = .false.
    real:: n_min, n_max
  contains
    procedure:: init
    procedure:: update
    procedure, private:: evolve
    procedure, private:: scatter
    procedure, private:: test_kernel
    procedure:: consistency
    procedure:: consistency2
  end type
  integer:: verbose = 0                       ! information level
  integer:: cadence = 10                      ! particle steps per evolve step
  real(8):: dust_time = 0.01                  ! dust evolution timestep
  real:: growth_rate = 1e-2                   ! per unit code time
  integer:: n_bin = 15                        ! 5 bins per power of ten
  logical:: validate = .false.                ! check consistency
  logical:: on = .false.
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self, nbin, m_cgs, dd)
  class(dust_t):: self
  integer:: nbin
  real, pointer, contiguous:: m_cgs(:), dd(:,:,:,:)
  !.............................................................................
  integer:: iostat, bin, m(3)
  real:: dln_m
  logical:: first_time = .true.
  namelist /dust_params/ on, cadence, dust_time, growth_rate, validate, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%init')
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_params, iostat=iostat)
    call io%namelist_warning ('dust_params', iostat)
    write (io%namelist, dust_params)
    n_bin = nbin
    first_time = .false.
  end if
  call io%lock%unset
  self%m_cgs => m_cgs
  self%dd => dd
  self%on = on
  self%kernel => test_kernel
  self%verbose = verbose
  m = self%patch%gn
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state
!===============================================================================
SUBROUTINE update (self, particles)
  class(dust_t):: self
  type(particles_t), pointer:: particles(:)
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
  allocate (self%dv(n_bin))                     ! average speed difference
  !allocate (self%na(n(1),n(2),n(3),n_bin))      ! number of added particles
  !allocate (self%np(n(1),n(2),n(3),n_bin))      ! work array for pp-cell-bin
  !call self%gather
  !call self%renormalize
  if (do_dust) then
    allocate (self%k(n_bin,n_bin))              ! kernel expressing exchange
    call io%bits_mem (storage_size(self%k), product(shape(self%k)), 'k')
    call self%evolve
    call self%scatter
    call io%bits_mem (storage_size(self%k), -product(shape(self%k)), '-k')
    deallocate (self%k)
  end if
  !deallocate (self%na, self%dv, self%np)
  deallocate (self%dv)
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
  !pa => self%particles%head                     ! first particle
  !do while (associated(pa))                     ! loop over particle
  !  bin = pa%bin                                ! particle bin
  !  i = pa%q                                    ! particle space index
  !  pa%w = pa%w*self%dd(i(1),i(2),i(2),bin)      ! renormalize
  !  pa => pa%next                               ! next particle
  !end do
  call trace%end (itimer)
END SUBROUTINE scatter

!===============================================================================
!===============================================================================
SUBROUTINE consistency (self, label)
  class(dust_t):: self
  character(len=*):: label
  class(particles_t), pointer:: s
  class(particle_t), pointer:: p
  integer:: n, l(3), u(3), j
  integer(kind=2):: i1, i2, i3, bin
  !-----------------------------------------------------------------------------
  if (.not.validate) return
  print *,label
  l = self%patch%mesh%lb
  u = self%patch%mesh%ub
  p => self%particles(1)%head%next
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
        s => self%sorted(i1,i2,i3,bin)
        n = n + s%n
      end do
    end do
  end do
  end do
END SUBROUTINE consistency

!===============================================================================
!> Check that the head to tail range is consistent with the cell count
!===============================================================================
SUBROUTINE consistency2 (self, label, i1, i2, i3, bin)
  class(dust_t):: self
  character(len=*):: label
  integer(kind=2):: i1, i2, i3, bin
  class(particles_t), pointer:: s
  class(particle_t), pointer:: p
  integer:: n, j
  !-----------------------------------------------------------------------------
  if (.not.validate) return
  s => self%sorted(i1,i2,i3,bin)
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
