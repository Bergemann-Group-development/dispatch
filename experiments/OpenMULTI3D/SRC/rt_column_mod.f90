!===============================================================================
!> One RT task per frequency window and ray direction.  An rt_column_t sweeps
!> the depth slabs (rt_patch_t) of its window in the direction of its rays,
!> one slab at a time, and hands the exit plane of each slab to the next.  The
!> boundary buffer therefore exists only while the column runs, and at most one
!> per busy thread.
!>
!> Within a slab the column runs the update sequence a slab task had
!> (pre_update, update, rotate), so each slab keeps its own iteration count,
!> time and files.  With premature_finish, slabs freeze from the bottom once
!> their atoms and those of all deeper slabs have converged.  In an up column
!> the highest frozen slab retains its last exit plane, in memory or with
!> rt_to_file on disk, and the column passes that plane on in place of a new
!> one.  The column finishes with its last slab.
!===============================================================================
MODULE rt_column_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE mpi_mod, only: mpi
  USE trace_mod
  USE bits_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t, link_t
  USE list_mod, only: list_t
  USE shared_mod
  USE console_mod
  USE m3d_patch_mod, only: m3d_patch_t, hash_string, rotate_atmos, mem_limit_gb
  USE rt_params_mod, only: rt_out
  USE rt_patch_mod, only: rt_patch_t, rt_patch_ptr_t
  USE composition_patch_mod, only: composition_patch_t
  USE nlte_patch_mod, only: nlte_patch_t
  implicit none
  private

  type, public, extends(m3d_patch_t) :: rt_column_t
    type(rt_patch_ptr_t), allocatable :: slab(:)   ! in sweep order
  contains
    procedure :: hash_key
    procedure :: execute
    procedure :: finish
    procedure :: nbor_relations
    procedure :: init_nbors_hash
  end type rt_column_t

  public :: init_rt_columns, report_memory_estimate, check_memory_limit

contains

!===============================================================================
!> Group the RT slab states of all composition patches into columns, one per
!> (i_nu, i_mu), and append the columns to the task list.  Called once all
!> slabs exist and before the nbor lists are built.
!===============================================================================
SUBROUTINE init_rt_columns (task_list)
  class(list_t) :: task_list
  class(link_t), pointer :: link
  type(rt_patch_ptr_t), allocatable :: slabs(:), grown(:)
  class(rt_column_t), pointer :: column
  logical, allocatable :: grouped(:)
  integer, allocatable :: member(:)
  integer :: i, j, k, tmp
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_column_mod%init_rt_columns')
  ! slabs and member are reallocated explicitly: the optimized builds use
  ! -assume norealloc_lhs
  allocate(slabs(0))
  link => task_list%head
  do while (associated(link))
    select type (task => link%task)
    class is (composition_patch_t)
      if (allocated(task%rt_slabs)) then
        call set_twins(task%rt_slabs)
        allocate(grown, source=[slabs, task%rt_slabs])
        call move_alloc(grown, slabs)
      endif
    end select
    link => link%next
  end do
  if (size(slabs) > 0 .and. mpi%size > 1) then
    call console_error ('RT columns need all depth slabs on one MPI rank')
    error stop 'MPI runs are not supported'
  endif

  allocate(grouped(size(slabs)), source=.false.)
  do i=1, size(slabs)
    call slabs(i)%p%init_unique (same=.true.)
    if (grouped(i)) cycle
    if (allocated(member)) deallocate(member)
    allocate(member, source=pack([(j, j=1, size(slabs))], [(slabs(j)%p%i_nu == slabs(i)%p%i_nu .and. &
                                 slabs(j)%p%i_mu == slabs(i)%p%i_mu, j=1, size(slabs))]))
    grouped(member) = .true.
    ! sweep order: along the rays, from the slab where they enter
    do j=2, size(member)
      k = j
      do while (k > 1)
        if (slabs(member(k-1))%p%ipos(3) <= slabs(member(k))%p%ipos(3)) exit
        tmp = member(k); member(k) = member(k-1); member(k-1) = tmp
        k = k - 1
      end do
    end do
    if (slabs(i)%p%dir < 0) member = member(size(member):1:-1)
    allocate(column)
    allocate(column%slab(size(member)))
    do j=1, size(member)
      column%slab(j)%p => slabs(member(j))%p
    end do
    call init_column (column)
  end do
  call trace%end()
END SUBROUTINE init_rt_columns

!===============================================================================
!> The i_mu twins of one depth slab exchange their hemispheres
!> (combine_twin_hemispheres), and the final pass reads the output pair.
!===============================================================================
SUBROUTINE set_twins (slabs)
  type(rt_patch_ptr_t), intent(in) :: slabs(:)
  class(rt_patch_t), pointer :: s, other
  integer :: i, j
  !-----------------------------------------------------------------------------
  do i=1, size(slabs)
    s => slabs(i)%p
    if (s%has_twin()) s%up_twin => s
    do j=1, size(slabs)
      other => slabs(j)%p
      if (other%i_nu /= s%i_nu .or. other%do_final) cycle
      if (s%do_final .and. rt_out%maxiter > 0) then
        if (other%i_mu == 1) s%up_twin => other
        if (other%i_mu == 2) s%dn_twin => other
      elseif (s%has_twin() .and. other%i_mu == 3 - s%i_mu) then
        s%dn_twin => other
      endif
    end do
  end do
END SUBROUTINE set_twins

!===============================================================================
!> A column takes its time model, flags and maxiter from the kind of its slabs,
!> its level from the lowest slab and its extent from all of them.  It has no
!> mem and nothing to download.
!===============================================================================
SUBROUTINE init_column (self)
  class(rt_column_t), target :: self
  class(rt_patch_t), pointer :: s
  real(8) :: lower, upper
  integer :: j
  !-----------------------------------------------------------------------------
  s => self%slab(1)%p
  self%kind = 'rt_column'
  self%i_nu = s%i_nu
  self%i_mu = s%i_mu
  self%dir = s%dir
  self%patch => s%patch
  do j=2, size(self%slab)
    if (self%slab(j)%p%ipos(3) < self%patch%ipos(3)) self%patch => self%slab(j)%p%patch
  end do
  call self%m3d_patch_t%pre_init()
  ! pre_init runs on the parent component, and so sets the parent's key
  self%key = self%hash_key()
  self%initialized = .true.
  ! the extent of the whole column, so that nbor lists (which require overlap)
  ! can hold the tasks of every slab
  lower = minval([(self%slab(j)%p%position(3) - 0.5d0*self%slab(j)%p%size(3), j=1, size(self%slab))])
  upper = maxval([(self%slab(j)%p%position(3) + 0.5d0*self%slab(j)%p%size(3), j=1, size(self%slab))])
  self%size(3) = upper - lower
  self%position(3) = 0.5d0*(lower + upper)
  self%n(3) = sum([(self%slab(j)%p%n(3), j=1, size(self%slab))])
  do j=1, size(self%slab)
    ! a composition patch checks the tasks on its key list when it updates
    call self%slab(j)%p%m3d_parent%add_to_key_list(self)
    self%slab(j)%p%column_end = j == size(self%slab) .and. &
      self%slab(j)%p%patch%boundaries%is_set(merge(bits%zu, bits%zl, self%dir > 0))
  end do
  call self%post_pre_init()
END SUBROUTINE init_column

!===============================================================================
FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
  class(rt_column_t):: self
  real(8), dimension(3), optional:: position, origin, size, box
  integer, optional:: level
  integer:: key(4)
  !-----------------------------------------------------------------------------
  key = [0, self%i_nu, self%i_mu, hash_string('rt_column')]
END FUNCTION hash_key

!===============================================================================
!> Sweep the slabs.  plane is the exit plane handed on by the slab before s.
!> Before the sweep, mark the slabs that may freeze (settled, see
!> rt_patch_pre_update): those whose atoms, and those of every deeper slab, have
!> converged.  The populations of this iteration are final, so both columns of
!> a window mark the same slabs.  A slab that freezes hands its exit plane on
!> as usual if the slab after it freezes too, and otherwise retains it for the
!> sweeps to come, in which the slab after it receives a copy.
!===============================================================================
SUBROUTINE execute (self)
  class(rt_column_t) :: self
  class(rt_patch_t), pointer :: s, prev
  real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: plane
  logical :: settled
  integer :: i, j, k, n
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_column_t%execute', itimer=itimer)
  n = size(self%slab)
  settled = .true.
  do k=1, n
    j = merge(k, n+1-k, self%dir > 0)
    s => self%slab(j)%p
    if (s%has_atom) settled = settled .and. all([(s%atoms(i)%data%converged, i=1, s%natom)])
    s%settled = settled .and. s%has_atom
  end do
  nullify (plane, prev)
  do j=1, n
    s => self%slab(j)%p
    if (s%is_finished()) then
      ! a finished slab no longer reads the plane before it
      call free (plane, prev)
      if (associated(prev)) then
        if (prev%static_mem) call prev%release_exit_plane()
      endif
    else
      if (associated(plane)) then
        s%mem => plane
        nullify (plane)
      elseif (associated(prev)) then
        if (.not. prev%static_mem) &
          call shared%abort ('rt_column_t%execute: no boundary for a slab inside the column')
        call prev%restore_exit_plane(s)
      endif
      s%upstream => prev
      call s%pre_update()
      call s%update()
      call s%rotate()
      if (s%static_mem .and. j < n) then
        if (self%slab(j+1)%p%settled) then
          plane => s%mem
          nullify (s%mem)
        else
          call s%retain_exit_plane()
        endif
      else
        plane => s%mem
        nullify (s%mem)
      endif
    endif
    prev => s
  end do
  call free (plane, prev)
  if (all([(self%slab(j)%p%is_finished(), j=1, n)])) call self%finish()
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  !> Free a plane that no slab reads any more, through the slab that made it.
  !-----------------------------------------------------------------------------
  subroutine free (plane, prev)
    real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: plane
    class(rt_patch_t), pointer :: prev
    !---------------------------------------------------------------------------
    if (.not. associated(plane)) return
    prev%mem => plane
    call prev%deallocate_mem()
    nullify (plane)
  end subroutine free
END SUBROUTINE execute

!===============================================================================
!> A column finishes when its iterations are done, when the global NLTE
!> iteration has converged (m3d_patch_t%pre_update) or when all of its slabs
!> have.  Slabs that have not finished themselves by then do so here, which
!> frees their twins' mean intensities.
!===============================================================================
SUBROUTINE finish (self)
  class(rt_column_t) :: self
  integer :: j
  !-----------------------------------------------------------------------------
  call self%m3d_patch_t%finish()
  do j=1, size(self%slab)
    if (.not. self%slab(j)%p%is_finished()) call self%slab(j)%p%finish()
    if (self%slab(j)%p%static_mem) call self%slab(j)%p%release_exit_plane()
  end do
END SUBROUTINE finish

!===============================================================================
!> Dependencies of a column:
!>   nlte_patch   needed and needs me: iteration k of an iteration column needs
!>                update k of every slab's populations, and update k+1 needs the
!>                rates of every column's iteration k.
!>   composition  needed: the slab's composition must exist.
!>   rt_column    the i_mu twin of a window, needed and needs me: iteration k+1
!>                reads the hemisphere the twin combined at the end of k; the
!>                final column needs the output scattering pair of its i_nu.
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(rt_column_t) :: self
  class(task_t), pointer :: another
  logical :: needed, needs_me, download
  !-----------------------------------------------------------------------------
  needed = .false.; needs_me = .false.; download = .false.
  select case (another%kind)
  case ('nlte_patch')
    needed = .true.; needs_me = .true.
  case ('composition')
    needed = .true.
  case ('rt_column')
    select type (another)
    class is (rt_column_t)
      if (another%i_nu == self%i_nu .and. another%i_mu /= self%i_mu) then
        if (self%do_final) then
          needed = .true.
        elseif (another%do_final) then
          needs_me = .true.
        else
          needed = .true.; needs_me = .true.
        endif
      endif
    end select
  end select
END SUBROUTINE nbor_relations

!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(rt_column_t):: self
  class(nbor_t), pointer:: nbors
  integer :: j, key(4)
  !-----------------------------------------------------------------------------
  do j=1, size(self%slab)
    key = [self%slab(j)%p%key(1), 0, 0, hash_string('nlte_patch')]
    call self%check_nbor(key, nbors)
    key(4) = hash_string('composition')
    call self%check_nbor(key, nbors)
  end do
  key = self%key
  if (self%do_final) then
    if (rt_out%maxiter > 0) then
      key(3) = 1
      call self%check_nbor(key, nbors)
      key(3) = 2
      call self%check_nbor(key, nbors)
    endif
  else
    key(3) = 3 - self%i_mu
    call self%check_nbor(key, nbors)
    if (self%do_out_sp) then
      key(3) = -1
      call self%check_nbor(key, nbors)
    endif
  endif
END SUBROUTINE init_nbors_hash

!===============================================================================
!> Per-rank estimate of the peak registered allocation.  The EOS setup, the NLTE
!> iteration and the output pass run one after the other (output tasks wait for
!> nlte_finished), so the estimate is the startup allocation and the retained
!> EOS species plus the largest phase.  A phase holds, per busy thread, the
!> workspace of one slab and the boundary buffer of its column, plus the
!> retained RT fields of all slabs.  Allocations outside Dispatch's registry,
!> opacity caches, solver temporaries and, without rt_to_file, the exit plane
!> an up column's frozen slabs retain are excluded, so this is an estimate, not an
!> upper bound.  Slabs
!> allocate u_phase when they first execute and keep it, so both RT phases
!> count it for every slab.
!===============================================================================
SUBROUTINE report_memory_estimate (list, estimate)
  class(list_t), intent(in) :: list
  real(8), intent(out) :: estimate
  class(link_t), pointer :: link
  real(8), parameter :: gib=1024d0**3
  real(8) :: work, retained, boundary, persistent, transient, eos_retained, eos_work, rt_persistent
  real(8) :: max_work(2), sum_retained(2), max_boundary(2), phase(3), atom(2), atom_nlte, atom_out
  integer :: n_column(2), n_comp, workers(2), eos_workers, j, k
  character(len=6), parameter :: label(2) = ['NLTE  ', 'output']
  !-----------------------------------------------------------------------------
  n_column = 0; n_comp = 0
  max_work = 0d0; sum_retained = 0d0; max_boundary = 0d0
  eos_work = 0d0; eos_retained = 0d0; rt_persistent = 0d0; atom = 0d0
  link => list%head
  do while (associated(link))
    if (link%task%is_clear(bits%remove+bits%external+bits%virtual)) then
      select type (task => link%task)
      class is (rt_column_t)
        k = merge(2, 1, task%do_out_sp)
        n_column(k) = n_column(k) + 1
        do j=1, size(task%slab)
          call task%slab(j)%p%estimate_memory(work, retained, boundary, persistent)
          rt_persistent = rt_persistent + persistent
          max_work(k) = max(max_work(k), work)
          sum_retained(k) = sum_retained(k) + retained
          max_boundary(k) = max(max_boundary(k), boundary)
        end do
      class is (composition_patch_t)
        call task%estimate_eos_memory(transient, retained)
        n_comp = n_comp + 1
        ! the patch being loaded holds its own retained arrays inside transient
        eos_work = max(eos_work, transient - retained/merge(1+task%n_out, 1, rotate_atmos))
        eos_retained = eos_retained + retained
      class is (nlte_patch_t)
        call task%estimate_atom_memory(atom_nlte, atom_out)
        atom = atom + [atom_nlte, atom_out]
      end select
    endif
    link => link%next
  end do
  workers = min(omp%nthreads, n_column)
  eos_workers = min(omp%nthreads, n_comp)
  phase(1) = eos_workers*eos_work
  phase(2:3) = workers*(max_work + max_boundary) + sum_retained + rt_persistent + atom
  estimate = os%gb + (eos_retained + maxval(phase))/gib
  if (.not. console_on(ch_m3d, 1, all_ranks=.true.)) return
  write(console_unit,'(a,i0,a)') 'Memory estimate rank ', mpi%rank, ' [GiB]:'
  write(console_unit,'(a,f10.4)') '  startup registered allocation ', os%gb
  write(console_unit,'(a,f10.4)') '  retained EOS species          ', eos_retained/gib
  write(console_unit,'(a,f10.4,a,i0,a)') '  EOS phase                     ', phase(1)/gib, &
    '  (', eos_workers, ' concurrent composition patches)'
  do k=1, 2
    if (n_column(k) == 0) cycle
    write(console_unit,'(2x,a6,a,f10.4,a,i0,a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a)') label(k), &
      ' phase                  ', phase(k+1)/gib, '  (', workers(k), ' of ', n_column(k), &
      ' columns x (', max_work(k)/gib, ' work + ', max_boundary(k)/gib, ' buffer), ', &
      sum_retained(k)/gib, ' retained, ', rt_persistent/gib, ' u_phase, ', atom(k)/gib, ' atom arrays)'
  end do
  write(console_unit,'(a,f10.4)') '  estimated peak allocation     ', estimate
  flush(console_unit)
END SUBROUTINE report_memory_estimate

!===============================================================================
!> Stop before the dispatcher starts when the estimated peak exceeds the memory
!> the run may use: mem_limit_gb if positive, else the physical memory of the
!> machine.  The estimate leaves out some allocations, so a run that passes may
!> still run out of memory, but one that fails would have.
!===============================================================================
SUBROUTINE check_memory_limit (estimate)
  real(8), intent(in) :: estimate
  real(8) :: limit
  character(len=160) :: msg
  !-----------------------------------------------------------------------------
  if (mem_limit_gb < 0d0) return
  limit = mem_limit_gb
  if (limit == 0d0) limit = physical_memory_gib()
  if (limit <= 0d0) then
    call console_error ('Physical memory unknown; set mem_limit_gb in &run_params to enable the memory check')
    return
  endif
  if (estimate <= limit) return
  write(msg,'(a,f0.1,a,f0.1,a)') 'Estimated peak memory ', estimate, ' GiB exceeds the limit of ', &
    limit, ' GiB (mem_limit_gb in &run_params; a negative value skips the check)'
  call console_error (trim(msg))
  error stop 'insufficient memory'
END SUBROUTINE check_memory_limit

!===============================================================================
!> Physical memory of this machine [GiB], 0 if it cannot be determined.  Linux
!> lists it in /proc/meminfo; elsewhere sysctl hw.memsize gives it in bytes.
!===============================================================================
FUNCTION physical_memory_gib () RESULT (gib)
  real(8) :: gib
  character(len=128) :: line
  integer(8) :: value
  integer :: unit, iostat
  !-----------------------------------------------------------------------------
  gib = 0d0
  open (newunit=unit, file='/proc/meminfo', status='old', action='read', iostat=iostat)
  if (iostat == 0) then
    do
      read (unit, '(a)', iostat=iostat) line
      if (iostat /= 0) exit
      if (line(1:9) == 'MemTotal:') then
        read (line(10:), *) value             ! kB
        gib = real(value,8)/1024d0**2
        exit
      endif
    end do
    close (unit)
    return
  endif
  call execute_command_line ('sysctl -n hw.memsize > .hw_memsize 2> /dev/null', exitstat=iostat)
  if (iostat /= 0) return
  open (newunit=unit, file='.hw_memsize', status='old', action='read', iostat=iostat)
  if (iostat /= 0) return
  read (unit, *, iostat=iostat) value
  if (iostat == 0) gib = real(value,8)/1024d0**3
  close (unit, status='delete')
END FUNCTION physical_memory_gib

END MODULE rt_column_mod
