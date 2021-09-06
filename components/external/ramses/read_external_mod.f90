!===============================================================================
!> Read a RAMSES snapshot and convert to DISPATCH data structures.
!>
!> An efficient method would be to read each RAMSES patch (oct), creat it as
!> a temporary, and check if it overlaps with the cut-out-region.  If it does,
!> then search for patches it overlaps with, and if the max level of these is
!> smaller than the patch level, then create a higher level DISPATCH patch there.
!> The lowest level DISPATCH does not have to be fully populated (or even at all
!> populated -- we could have an irregular collection.
!>
!> There are three overlap decisions that need to be consistent:
!>
!> 1) First the choice of RAMSES octs to add to the oct_list
!> 2) The criterion for adding RAMSES octs to the DISPATCH nbor lists
!> 3) The criterion for creating DISPATCH patches
!>
!> In the end, we should have
!> A) some patches with resolution up to level=maxlvl
!> B) Coverage up to that level within a cube of size cut_out_size
!> C) Successively larger support levels doen to a frozen level=minlvl
!>
!>                                                        cell level  RAMSES
!>                   ......|......              32k              24          21
!>                   ......|......                .               .
!>                   ......|......               64               .
!>                   ......|......                8              20
!>                   ......|......                1              19          16
!>            .............|...............                      18          15
!> ........................|...............................      17          14
!>
!> If, for example, we want full coverage of the innermost 120 AU, which
!> corresponds to a level 16 RAMSES cell, then that level must contain at least
!> one DISPATCH 16x16x16 patch, which means a DISPATCH level 20 patch is needed.
!> Level 28 would then contain up to 256 patches, which is fine for laptop work.
!>
!> To include structures up to 1000 AU, we should have a single DISPATCH patch
!> at levels 19 (240 AU), 18 (480 UA), and 17 (960 AU).  We could conveniently
!> go to level 16 DISPATCH, at 2000 AU, and freeze that one patch.  The size of
!> each level should then be constant (120 AU) from level 20 and up, with the
!> formula size=0.5**min(level-4,16).  One should sample RAMSES grids from level
!> 13 (where the cell size is large enough to over a level 16 DISPATCH patch. So,
!> at level
!> 13: sample RAMSES grids at level 13, out to distance 1 cell
!> 14:                              14                  1 cell
!> 17:                              17                  2 cells
!> 20:                              20                 16 cells = 1 patch
!> 24                               24                256 cells =
!>
!> minlvl=17 input should mean the frozen cells are at level 17, and there are 16
!> of them, corresponding to one level 13 cell (= 1000 AU).  To ensure support,
!> level 18, 19, and 20 should be a factor of 2 smaller per level, so the size
!> should contain the expression 16*min(level, minlvl+3), or 8*min() for half-size.
!> init
!>   read_snapshot                  ! read and amass RAMSES grids
!>   create_cartesian               ! create Cartesian patches
!>       cartesian_init             ! make one DISPATCH patch
!>   download_values                ! pick up values from RAMSES
!>   remove_ramses_octs             ! deallocate RAMSES oct list
!>   init_all_nbors
!>   reset_status
!>   do while (associated(link))
!>     patch%dnload
!>     patch%update
!>     patch%info
!>   end do
!===============================================================================
MODULE read_external_mod
  USE io_mod
  USE omp_mod
  USE trace_mod
  USE bits_mod
  USE timer_mod
  USE mpi_mod
  USE task_mod
  USE link_mod
  USE patch_mod
  USE gpatch_mod
  USE extras_mod
  USE mesh_mod
  USE mhd_mod
  USE scaling_mod
  USE units_mod
  USE star_mod
  USE list_mod
  USE task_list_mod
  USE solver_mod
  USE experiment_mod
  USE download_mod
  USE index_mod
  USE dll_mod
  USE hash_table_mod
  USE omp_timer_mod
  USE shared_mod
  USE stagger_mod
  USE eos_mod
  implicit none
  private
  integer, parameter :: dp=kind(1.0D0)                        ! double precision
  integer, parameter :: d1=kind(kind(1d0))                    ! kind of a kind
  integer, parameter :: nmetal=2
  type :: sink_t                        ! star data type
    SEQUENCE
    real(dp) x,y,z                      ! positions
    real(dp) px,py,pz                   ! momenta
    real(dp) ux,uy,uz                   ! velocities
    real(dp) dpx,dpy,dpz                ! momenta
    real(dp) mass                       ! mass
    real(dp) dmass                      ! mass
    real(dp) potential                  ! gravitational potential
    real(dp) density                    ! neigborhood density
    real(dp) t_create                   ! creation time
    real(dp) t_explode                  ! explosion time
    real(dp) metal(nmetal)              ! metal density
    real(dp) dmetal(nmetal)             ! metal density
    real(dp) beta_create                ! Beta at creation
    real(dp) dtold                      ! previous dt
    integer gid                         ! pointer to grid cell
    integer father                      ! pointer to father cell
    integer ilevel                      ! level of AMR hierachy where star lives
    integer old_ilevel                  ! level of hierachy star used to have, if changed
    integer owner                       ! owner thread
    integer has_moved
  end type
  !-----------------------------------------------------------------------------
  ! Light-weight RAMSES oct, using 4*(2+4+3+8*12) = 4*105 = 420 bytes, which
  ! means we can read in up to about 500 million octs on a frontend
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: oct_t
    real(kind=4):: u_max=0.0, cost=0.0
    integer:: key(4)=0, ii(3)=0
    real, allocatable:: var(:,:)
  end type
  !-----------------------------------------------------------------------------
  ! List of RAMSES octs
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: oct_list_t
  end type
  !-----------------------------------------------------------------------------
  ! RAMSES data type
  !-----------------------------------------------------------------------------
  type, public:: read_external_t
    logical:: big_endian, double_precision, quad_precision
    integer:: rank
    integer:: nout, minlvl, maxlvl, rootlvl
    integer:: npatch, nmerged_max
    real(kind=4):: scale_l, scale_d, scale_t
    real(kind=8):: rank_load
    real(kind=8):: xmin, xmax, ymin, ymax, zmin, zmax         ! corners of the cutout
    real(kind=8):: size, d_min, tt_min
    real(kind=4):: u_min=huge(1.0), u_max=0.0
    integer:: n_bins = 200
    integer, allocatable:: u_bins(:)
    logical :: centre_on_sink
    character(len=128):: path
    integer n_cell                      ! number of cells per direction
    type(oct_list_t), pointer:: oct_list
    class(task_list_t), pointer:: task_list
    type(task_list_t), pointer:: dispatch_task_list
    type(sink_t), dimension(:), pointer:: sink
    type(hash_table_t):: u_max_table
    type(hash_table_t):: patch_table
  contains
    procedure:: init
    procedure:: read_snapshot
    procedure:: phase_speed
    procedure:: create_patches
    procedure:: add_to
    procedure:: add_to_task_list
    procedure:: u_max_table_get
    procedure:: count_costs
    procedure:: energy
    procedure:: read_sink
    procedure:: evaluate
    procedure:: cartesian_init
    procedure:: initialize
    procedure:: init_all_nbors
    procedure:: init_nbors
    procedure:: remove_ramses_octs
    procedure:: check_support
    procedure:: prolong
    procedure:: dnload
    procedure:: download_link
    procedure:: check_filled
    procedure:: output
    procedure:: total_energy
  end type
  logical, save:: make_cartesian        ! replace RAMSES patches
  logical, save:: centre_on_sink        ! centre on the last sink
  logical, save:: do_evaluate           ! evaluate costs
  logical, save:: remove_b=.false.      ! remove B-field, for tests
  logical, save:: do_memtest=.false.    ! remove B-field, for tests
  logical, save:: direct_load=.false.
  logical, save:: detailed_timer=.false.
  logical, save:: fast_nbors=.true.
  integer, save:: oct_size=700          ! oct_t overhead (ifort)
  integer, save:: verbose=0             ! verbosity
  real(8), save:: cut_out_size=1d0      ! full size of cutout
  real(8), save:: start=0d0
  integer, allocatable:: level_cost(:)
  type(read_external_t), public:: read_external
CONTAINS

!===============================================================================
!> Facilitate call to %read_snapshot().  All heavy lifting is done there
!===============================================================================
SUBROUTINE init (self, task_list, name)
  class(read_external_t):: self
  class(task_list_t), pointer:: task_list
  integer:: iostat
  character(len=*), optional:: name
  logical, save:: first_time=.true.
  integer, save:: hash_table_size=5000000, n_cell=16
  integer, save:: nout, levelmin, levelmax, nmerged_max, patch_n
  real, save:: xmin, xmax, ymin, ymax, zmin, zmax, rank_load
  real, save:: d_min, tt_min, size, size_au
  character(len=64), save:: path='.'
  logical:: keep, load
  namelist /ramses_reader_params/ nout, xmin, xmax, ymin, ymax, zmin, zmax, &
    size, n_cell, nmerged_max, levelmin, levelmax, rank_load, path, centre_on_sink, &
    make_cartesian, do_evaluate, patch_n, remove_b, cut_out_size, d_min, &
    detailed_timer, direct_load, do_memtest, fast_nbors, verbose
  !.............................................................................
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%init ('ramses snapshot reader')
  call io%header ('ramses snapshot reader:')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and saves them
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time=.false.
    ! default values
    nout = -1
    xmin = 0.0; xmax = 1.0
    ymin = 0.0; ymax = 1.0
    zmin = 0.0; zmax = 1.0
    levelmin = shared%levelmin
    levelmax = shared%levelmax
    rank_load = 1.0
    d_min = 1d-18 ! cgs (1e6 code units)
    tt_min = 5.   ! K
    path = './'
    size = 1d0
    nmerged_max = 16
    make_cartesian = .true.
    centre_on_sink = .false.
    do_evaluate = .true.
    remove_b = .false.
    patch_n = 16
    !--- read namelist ---
    rewind (io%input)
    read (io%input, ramses_reader_params, iostat=iostat)
    call io%namelist_warning ('ramses_reader_params', iostat)
    write(io%namelist, ramses_reader_params)
    shared%levelmin = levelmin
    shared%levelmax = levelmax
    io%nml_version = 4
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! All threads save the parameters in their private parts, including the
  ! scaling information.  The selection of RAMSES data we read is based on
  ! the parameters defined here.
  !-----------------------------------------------------------------------------
  call scaling%init()
  self%n_cell    = n_cell
  self%d_min     = d_min/scaling%d
  self%tt_min    = tt_min/scaling%temp
  self%xmin      = xmin
  self%xmax      = xmax
  self%ymin      = ymin
  self%ymax      = ymax
  self%zmin      = zmin
  self%zmax      = zmax
  self%minlvl    = levelmin
  self%maxlvl    = levelmax
  self%path      = path
  self%npatch    = 1
  self%rank_load = rank_load
  self%nout      = nout
  self%size      = size
  self%centre_on_sink = centre_on_sink
  self%nmerged_max = nmerged_max
  !-----------------------------------------------------------------------------
  ! Each rank reads the RAMSES grids belonging to it, filling in the data,
  ! and initializing the patch structure of every grid.
  !-----------------------------------------------------------------------------
  allocate (task_list, self%oct_list)
  self%task_list => task_list
  call self%task_list%init ('ramses_mod')
  call self%oct_list%init ('ramses_snapshot')
  if (do_memtest) &
    call mem_test
  !-----------------------------------------------------------------------------
  ! If there is enough memory, read and keep the oct%var values and create
  ! patches based on that
  !-----------------------------------------------------------------------------
  if (direct_load) then
    call self%read_snapshot (keep=.true., load=.false.)
    if (do_evaluate) then
      call self%count_costs
    end if
    !call self%create_patches
  !-----------------------------------------------------------------------------
  ! If not, read once to decide load distribution, and once again to load data
  !-----------------------------------------------------------------------------
  else
    call self%read_snapshot (keep=.false., load=.false.)
    if (do_evaluate) then
      call self%count_costs
    end if
    call self%read_snapshot (keep=.false., load=.true.)
  end if
  !-----------------------------------------------------------------------------
  ! Delete the RAMSES octs and check patch support
  !-----------------------------------------------------------------------------
  call self%remove_ramses_octs
  call self%check_support
  !-----------------------------------------------------------------------------
  ! Initialize nbor lists, which as a byproduct adds tasks to the ready queue
  !-----------------------------------------------------------------------------
  if (fast_nbors) then
    call self%init_all_nbors
  else
    call self%task_list%init_all_nbors
  end if
  if (verbose > 0) &
    call log_nbor_lists
  !-----------------------------------------------------------------------------
  ! Prolong to unfilled interior points, dnload guard zones and check that
  ! everything is filled
  !-----------------------------------------------------------------------------
  call self%prolong
  call self%dnload
  call self%check_filled
  call self%total_energy
  call self%output
  !-----------------------------------------------------------------------------
  ! Add all tasks to the ready queue, preventing a slow default add_by_time
  !-----------------------------------------------------------------------------
  call self%task_list%queue_all
  call self%task_list%print_queue
  call trace%end()
contains

!===============================================================================
subroutine log_nbor_lists
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    write (stdlog,'(a,i6,":",i2.2)') 'log_nbor_lists: link=', &
      link%task%id, link%task%level
    nbor => link%nbor
    do while (associated(nbor))
      write (stdlog,'(2x,a,i6,":",i2.2)') 'log_nbor_lists: nbor=', &
        nbor%task%id, nbor%task%level
      nbor => nbor%next
    end do
    link => link%next
  end do
end subroutine log_nbor_lists
END SUBROUTINE init

!===============================================================================
!> Memory use test
!===============================================================================
SUBROUTINE mem_test
  integer:: n=1000000
  real:: ms=6e3
  type(oct_t)       , allocatable:: oct(:)
  type(task_t)      , allocatable:: task(:)
  type(patch_t)     , allocatable:: patch(:)
  type(gpatch_t)    , allocatable:: gpatch(:)
  type(extras_t)    , allocatable:: extras(:)
  type(mhd_t)       , allocatable:: mhd(:)
  type(solver_t)    , allocatable:: solver(:)
  type(experiment_t), allocatable:: experiment(:)
  !-----------------------------------------------------------------------------
  print*,'MEMORY TEST: Top GB gives kB per instance'  ; call omp%delay (ms)
  print*,'oct'    ; allocate (oct(n))    ; call omp%delay(ms); deallocate(oct)
  print*,'task'   ; allocate (task(n))   ; call omp%delay(ms); deallocate(task)
  print*,'patch'  ; allocate (patch(n))  ; call omp%delay(ms); deallocate(patch)
  print*,'extras' ; allocate (extras(n)) ; call omp%delay(ms); deallocate(extras)
  print*,'mhd'    ; allocate (mhd(n))    ; call omp%delay(ms); deallocate(mhd)
  print*,'solver' ; allocate (solver(n)) ; call omp%delay(ms); deallocate(solver)
END SUBROUTINE mem_test

!===============================================================================
!> Read a RAMSES snapshot (likely split across multiple files), and create a
!> patch for each "grid".
!===============================================================================
SUBROUTINE read_snapshot (self, keep, load)
  class(read_external_t):: self
  logical:: keep, load
  !
  !---- variables for reading the snapshot ----
  ! scalars
  integer(kind=4):: ilun=10, ilun2=11, nlevelmin, icpu, isentinel, ncpu, ndims, &
                    nx, ny, nz, nlevelmax, nboundary, ngrids, noutput, iout, ifout, &
                    overload, ndomain, ilevel, i, j, k, l, ivar, ioerr, &
                    bit_length, maxdom, imin, imax, jmin, jmax, kmin, kmax, ngrid_act, &
                    nx_full, ny_full, nz_full, ncpu_read, ix, iy, iz, ncache, ndom, &
                    ncells, max_dom, nvar_mhd, twotondim, ngrids_max, ngrids_k, it
  real(kind=8):: boxlen, dx, time, gamma_ramses, tmp_load, tot_load
  real(kind=16):: dmax, dkey, order_min, order_min_1d(1)
  real(kind=8):: r, rr(3)
  ! arrays
  integer(kind=4), dimension(1:8):: idom, jdom, kdom, cpu_min, cpu_max
  integer(kind=4), allocatable, dimension(:):: cpu_domain, isub_domain, index_grid, iig, &
                                               father, next, prev, cpu_list
  integer(kind=4), allocatable, dimension(:,:):: head_list, tail_list, numb_list, numb_tot, &
                                                 head_bndy, tail_bndy, numb_bndy, map_domain, &
                                                 son, neighbour, ngrid_bound, ngrid_level, &
                                                 ngrid_file
  logical, allocatable, dimension(:):: ref, cpu_read
  real(kind=4), allocatable, dimension(:):: fill
  real(kind=8), allocatable, dimension(:):: dtold, dtnew, bound_key_dp
  real(kind=16), dimension(1:8):: bounding_min, bounding_max
  real(kind=16), allocatable, dimension(:):: bound_key
  real(kind=8), allocatable, dimension(:,:):: xc, xg, x
  real(kind=8), allocatable, dimension(:,:):: bound_key_dpx2
  real(kind=8), allocatable, dimension(:,:,:):: var
  real(kind=8), dimension(3):: pmin=1d30, pmax=-1d30
  real(kind=8), pointer:: cost
  character(len=80):: ordering
  character(len=128):: filebase, filename, conv
  !
  ! derived types
  type:: level_t
    integer(kind=4):: ngrid, imin, imax, jmin, jmax, kmin, kmax
  end type level_t
  type(level_t), allocatable, dimension(:):: levels
  !
  ! diagnostic variables
  real:: d_min, d_min_level, u_max
  integer, allocatable, dimension(:,:):: ncellsum
  integer, save:: nprint_f=1000
  integer(8), save:: n_passive=0, n_active=0, n_onrank=0, n_offrank=0, &
    n_d_min=0, n_tt_min=0, n_frozen=0, n_grid
  !
  !---- dispatch variables ----
  class(task_t), pointer:: task, parent
  class(experiment_t), pointer:: patch
  class(dll_node_t), pointer:: dll_node
  type(oct_t), pointer:: oct
  real(8):: cut, ds, size, position(3), total_cost, weight, patch_cost
  integer, allocatable:: n_level(:)
  integer:: hash_table_size
  logical, save:: print_once=.true.
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin('read_snapshot', itimer=itimer)
  start = wallclock()
  !-----------------------------------------------------------------------------
  ! In a second sweep we alread know the number of octs
  !-----------------------------------------------------------------------------
  hash_table_size = 2000000
  if (direct_load) then
    call self%u_max_table%init (hash_table_size, ndim=4)
    call self%patch_table%init (hash_table_size, ndim=4)
    !hash_table_size=self%oct_list%n/50
  end if
  self%rank = 0
  self%rootlvl = floor(log(real(self%n_cell))/log(2.0))
  !
  write(filebase, '(a,"/output_",i5.5,"/amr_",i5.5,".out")') trim(self%path), self%nout, self%nout
  !
  ! We need some basic info about the snapshot; let's read the `info_` and `amr_` icpu = 1 files.
  !
  !$omp master
    !------ read snapshot `info_` file ------
    ! note that the `info_` and `amr_` files contain some redundant information.
    write(filename, '(a,"/output_",i5.5,"/info_",i5.5,".txt")') trim(self%path), self%nout, self%nout
    open(unit=ilun,file=filename,status='old',form='formatted')
    read(ilun,*) ! SKIP ncpu
    read(ilun,*) ! SKIP ndims
    read(ilun,'(13x,i11)') nlevelmin
    read(ilun,*) ! SKIP nlevelmax
    read(ilun,*) ! SKIP ngridmax
    read(ilun,*) ! SKIP nstep_coarse
    read(ilun,*) ! SKIP blank line
    read(ilun,*) ! SKIP boxlen
    read(ilun,*) ! SKIP time
    read(ilun,*) ! SKIP aexp
    read(ilun,*) ! SKIP H0
    read(ilun,*) ! SKIP omega_m
    read(ilun,*) ! SKIP omega_l
    read(ilun,*) ! SKIP omega_k
    read(ilun,*) ! SKIP omega_b
    read(ilun,'(13x,e23.15)') self%scale_l ! length scale factor
    read(ilun,'(13x,e23.15)') self%scale_d ! density scale factor
    read(ilun,'(13x,e23.15)') self%scale_t ! time scale factor
    close(ilun)
    !
    !------ read `amr_` file for icpu = 1 ------
    icpu = 1
    write(filename, '(a,i5.5)') trim(filebase), icpu
    !
#if __INTEL_COMPILER || defined(__GFORTRAN__)
    open(unit=ilun,file=filename,status='old',form='unformatted',access='stream')
    read(ilun) isentinel
    self%big_endian = (isentinel .ne. 4)  ! the data is in big endian format.
    close(ilun)
    !
    if (self%big_endian) then
      open(unit=ilun,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
    else
      open(unit=ilun,file=filename,status='old',form='unformatted')
    endif
#else
    open(unit=ilun,file=filename,status='old',form='unformatted')
#endif
    !---------------------------------------------------------------------------
    ! Establish if the file has quad precision Hilbert keys, assuming "new" fmt
    !---------------------------------------------------------------------------
    do i=1,25
      read(ilun)
    end do
    read(ilun) overload, ndomain
    backspace(ilun)
    allocate (bound_key(0:ndomain), bound_key_dp(0:ndomain), &
              isub_domain(1:ndomain), cpu_domain(1:ndomain), &
              map_domain(1:overload,1:ncpu))
    read(ilun,iostat=ioerr) &
      overload, ndomain, bound_key, isub_domain, cpu_domain, map_domain
    self%quad_precision = ioerr==0
    !---------------------------------------------------------------------------
    ! Rewind and read, without iostat
    !---------------------------------------------------------------------------
    rewind (ilun)
    read(ilun) ncpu
    read(ilun) ndims
    if (ndims < 3) then
      call mpi%abort('Only 3-D data is supported at this time.')
      close(ilun)
      stop
    endif
    read(ilun) nx, ny, nz
    read(ilun) nlevelmax
    read(ilun) ngrids_max
    read(ilun) nboundary
    read(ilun) ngrids
    !
    ! Try to read `r_boxlength` as a double precision variable; if an error occurs,
    ! then the data is single precision.
    read(ilun,iostat=ioerr) boxlen
    self%double_precision = (ioerr == 0)
    !
    ! At this time, only double precision RAMSES data is supported.
    if (.not. self%double_precision) then
      call mpi%abort('Only double-precision RAMSES data is supported at this time.')
      close(ilun)
      stop
    endif
    !
    read(ilun) noutput, iout, ifout
    if (verbose > 2) &
      write(stdlog,*)'noutput, iout, ifout:', noutput, iout, ifout
    read(ilun) ! SKIP tout
    read(ilun) ! SKIP aout
    read(ilun) time
    allocate(dtold(1:nlevelmax), dtnew(1:nlevelmax))
    read(ilun) dtold
    read(ilun) dtnew
    read(ilun) ! SKIP nstep, nstep_coarse
    read(ilun) ! SKIP const, mass_tot, rho_tot
    read(ilun) ! SKIP omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini
    read(ilun) ! SKIP aexp, hexp, aexp_old, epot_tot_ini, epot_tot_old
    read(ilun) ! SKIP mass_sph
    !
    ! levels variables
    allocate(head_list(1:ncpu,1:nlevelmax), tail_list(1:ncpu,1:nlevelmax), &
             numb_list(1:ncpu,1:nlevelmax), numb_tot(1:10,1:nlevelmax))
    read(ilun) head_list
    read(ilun) tail_list
    read(ilun) numb_list
    read(ilun) numb_tot
    !
    ! boundary linked lists
    if (nboundary > 0) then
      allocate(head_bndy(1:nboundary,1:nlevelmax), &
               tail_bndy(1:nboundary,1:nlevelmax), &
               numb_bndy(1:nboundary,1:nlevelmax))
      read(ilun) head_bndy
      read(ilun) tail_bndy
      read(ilun) numb_bndy
    endif
    !
    read(ilun) ! SKIP headf, tailf, numbf, used_mem, used_mem_tot
    read(ilun) ordering
    if (ordering == 'bisection') then
      call mpi%abort('Only Peano-Hilbert ordered RAMSES data can be read')
      close(ilun)
      stop
    endif
    !
    if (self%quad_precision) then
      write (stdout,*)  'Reading quad precision Hilbert keys'
      read(ilun) &
        overload, ndomain, bound_key, isub_domain, cpu_domain, map_domain
    else
      write (stdout,*)  'Reading double precision Hilbert keys'
      read(ilun) &
        overload, ndomain, bound_key_dp, isub_domain, cpu_domain, map_domain
      bound_key = bound_key_dp ! put the key into a quad precision variable.
    endif
    write (stdout,*) 'Using new format for key storage'
    deallocate (bound_key_dp, isub_domain, map_domain)
    ! ------ finished reading Peano-Hilbert keys ------
    !
    ! we have extracted enough information; close the file (for the time being).
    close(ilun)
    !
    ! output info about snapshot to terminal
    call io%header ('-RAMSES parameters from snapshot file')
    write(stdout,*) 'Data is '//merge('double','single',self%double_precision)//' precision with '// &
             trim(merge('quad  ','double',self%quad_precision))//' precision hilbert keys'
    write(stdout,*) '      ncpu =',ncpu
    write(stdout,*) ' nlevelmin =',nlevelmin
    write(stdout,*) ' nlevelmax =',nlevelmax
    allocate (n_level(1:nlevelmax))
    if (.not.allocated(level_cost)) &
      allocate (level_cost(nlevelmin:nlevelmax))
    !
    !------ allocate arrays ------
    ! now that we have more information, we can set and allocate a few required variables.
    twotondim = 2**ndims
    self%maxlvl = min(self%maxlvl, nlevelmax)
    if (self%maxlvl == 0) self%maxlvl = nlevelmax
    ncells = nx * ny * nz + twotondim * ngrids
    allocate(ngrid_file(1:ncpu+nboundary,1:nlevelmax), ngrid_level(1:ncpu,1:nlevelmax), &
             cpu_list(1:ncpu), cpu_read(1:ncpu), xc(1:twotondim,1:3), &
             father(1:ngrids_max), neighbour(1:ngrids_max,1:2*ndims), levels(1:nlevelmax), &
             ncellsum(1:nlevelmax,1:ncpu))
    if (nboundary > 0) allocate(ngrid_bound(1:nboundary,1:nlevelmax))
    !
    !------ given limits, make preparations to cut a sub-volume from the box ------
    ! uses the Hilbert-Peano curve to select which CPUs contain the relevant data
    dmax = max(self%xmax - self%xmin, self%ymax - self%ymin, self%zmax - self%zmin)
    ! for a given cut, what is the bare minimum level we need to consider.
    do ilevel=1,self%maxlvl
      dx = 0.5_8**ilevel
      if (dx < dmax) exit
    enddo
    bit_length = ilevel - 1
    if (verbose > 0) &
      write(stdout,*) 'need to go to at least level', bit_length
    !self%minlvl = bit_length
    max_dom = 2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    ndom = 1
    dkey = (2.0_16**(nlevelmax + 1) / real(max_dom,kind=16))**ndims
    if (bit_length > 0) then
      imin = int(self%xmin * dble(max_dom))
      imax = imin + 1
      jmin = int(self%ymin * dble(max_dom))
      jmax = jmin + 1
      kmin = int(self%zmin * dble(max_dom))
      kmax = kmin + 1
      ndom = 8
    end if
    !
    idom(1) = imin; idom(2) = imax
    idom(3) = imin; idom(4) = imax
    idom(5) = imin; idom(6) = imax
    idom(7) = imin; idom(8) = imax
    jdom(1) = jmin; jdom(2) = jmin
    jdom(3) = jmax; jdom(4) = jmax
    jdom(5) = jmin; jdom(6) = jmin
    jdom(7) = jmax; jdom(8) = jmax
    kdom(1) = kmin; kdom(2) = kmin
    kdom(3) = kmin; kdom(4) = kmin
    kdom(5) = kmax; kdom(6) = kmax
    kdom(7) = kmax; kdom(8) = kmax
    !
    if (bit_length > 0) then
      do i=1,ndom
        call hilbert3d(idom(i),jdom(i),kdom(i),order_min_1d,bit_length,1)
        order_min = order_min_1d(1)
        bounding_min(i) = (order_min) * dkey
        bounding_max(i) = (order_min + 1.0_16) * dkey
      enddo
    else
      bounding_min = 0.0_16
      bounding_max = dkey
    endif
    ! generate list of CPUs/output files to read.
    cpu_min = 0; cpu_max = 0
    cpu_read = .false.
    do j=1,ndomain
      do i=1,ndom
        if ( (bound_key(j-1) <= bounding_min(i)) .and. &
             (bound_key(j  ) >  bounding_min(i)) ) &
          cpu_min(i) = j
        if ( (bound_key(j-1) <  bounding_max(i)) .and. &
             (bound_key(j  ) >= bounding_max(i)) ) &
          cpu_max(i) = j
      enddo
    enddo
    deallocate (bound_key)
    ncpu_read = 0
    do i=1,ndom
      do j=cpu_min(i),cpu_max(i)
        if(.not. cpu_read(cpu_domain(j))) then
           ncpu_read = ncpu_read + 1
           cpu_list(ncpu_read) = cpu_domain(j)
           cpu_read(cpu_domain(j)) = .true.
        endif
      enddo
    enddo
    deallocate (cpu_domain)
    !
    ! calculate the full hierarchy
    do ilevel=1,nlevelmax
      nx_full = 2**ilevel
      ny_full = 2**ilevel
      nz_full = 2**ilevel
      levels(ilevel)%imin = int(self%xmin * dble(nx_full)) + 1
      levels(ilevel)%imax = int(self%xmax * dble(nx_full)) + 1
      levels(ilevel)%jmin = int(self%ymin * dble(ny_full)) + 1
      levels(ilevel)%jmax = int(self%ymax * dble(ny_full)) + 1
      levels(ilevel)%kmin = int(self%zmin * dble(nz_full)) + 1
      levels(ilevel)%kmax = int(self%zmax * dble(nz_full)) + 1
    enddo
    !
    ! ------ loop over CPU files ------
    call mpi%barrier ('ncpu begin')
    if (verbose > 0) then
      write(stdout,*) 'reading these files:'
      do k=1,ncpu_read
        icpu = cpu_list(k)
        write(filename, '(a,i5.5)') trim(filebase), icpu
        write(stdout,*) k, trim(filename)
      end do
    end if
    u_max = 0.0
    do k=1,ncpu_read
    !do k=1,32
      icpu = cpu_list(k)
      ! open `amr_` file
      write(filename, '(a,i5.5)') trim(filebase), icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (self%big_endian) then
        open(unit=ilun,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
      else
        open(unit=ilun,file=filename,status='old',form='unformatted')
      endif
#else
      open(unit=ilun,file=filename,status='old',form='unformatted')
#endif
      if (verbose > 2) &
        write(stdout,*) 'Reading file '//trim(filename)
      !
      ! skip header
      do i=1,6
        read(ilun) ! SKIP
      enddo
      read(ilun) ngrids_k
      do i=8,21
        read(ilun) ! SKIP
      enddo
      ! read grid numbers
      read(ilun) ngrid_level ! aka `numb_list`
      ngrid_file(1:ncpu,1:nlevelmax) = ngrid_level
      if (verbose > 3) &
        write (stdlog, '(a,i4.1,a,/,40i6)') &
        'ngrids(cpu=',k,')',ngrid_level(k,1:nlevelmax)
      !
      read(ilun) ! SKIP numb_tot
      if (nboundary > 0) then
        read(ilun) ! SKIP headb
        read(ilun) ! SKIP tailb
        read(ilun) ngrid_bound ! ala `numb_bndy`
        ngrid_file(ncpu+1:ncpu+nboundary,1:nlevelmax) = ngrid_bound
      endif
      read(ilun) ! SKIP headf, tailf, numbf, used_mem, used_mem_tot
      read(ilun) ! SKIP ordering
      if(TRIM(ordering) == 'bisection') then
        read(ilun) ! SKIP bisec_wall
        read(ilun) ! SKIP bisec_next
        read(ilun) ! SKIP bisec_indx
        read(ilun) ! SKIP bisec_cpubox_min
        read(ilun) ! SKIP bisec_cpubox_max
      else
        read(ilun) ! SKIP overload, ndomain, etc.
      endif
      read(ilun) ! SKIP son
      read(ilun) ! SKIP flag1
      read(ilun) ! SKIP cpu_map
      !
      ! now we open the `hydro_` file to get access to the hydro data
      write(filename, '(a,"/output_",i5.5,"/hydro_",i5.5,".out",i5.5)') trim(self%path), &
            self%nout, self%nout, icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (self%big_endian) then
        open(unit=ilun2,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
      else
        open(unit=ilun2,file=filename,status='old',form='unformatted')
      endif
#else
      open(unit=ilun2,file=filename,status='old',form='unformatted')
#endif
      read(ilun2) ! SKIP ncpu_dump
      read(ilun2) nvar_mhd
      read(ilun2) ! SKIP ndims
      read(ilun2) ! SKIP nlevelmax
      read(ilun2) ! SKIP nboundary
      read(ilun2) gamma_ramses
      eos%gamma = gamma_ramses
      !
      if (print_once) then
        print_once = .false.
        write(stdout,*) '      nvar =', nvar_mhd
        write(stdout,*) '    boxlen =', boxlen
        write(stdout,*) '      time =', time
        write(stdout,*) '     gamma =', eos%gamma
      end if
      !
      ! loop over levels for current file
      d_min = 1e30
      do ilevel=1,self%maxlvl
        n_grid = 0
        d_min_level = 1e30
        ! set cell geometry info
        dx = 0.5_8**ilevel
        nx_full = 2**ilevel
        ny_full = 2**ilevel
        nz_full = 2**ilevel
        do i=1,twotondim
          iz = (i - 1) / 4
          iy = (i - 1 - 4 * iz) / 2
          ix = (i - 1 - 2 * iy - 4 * iz)
          xc(i,1) = (dble(ix) - 0.5_8) * dx
          xc(i,2) = (dble(iy) - 0.5_8) * dx
          xc(i,3) = (dble(iz) - 0.5_8) * dx
        enddo
        !
        ! allocate worker/scratch arrays
        ngrid_act = ngrid_file(icpu,ilevel)
        levels(ilevel)%ngrid = ngrid_act
        if (ngrid_act > 0) allocate(xg(1:ngrid_act,1:ndims), son(1:ngrid_act,1:twotondim), &
                                    x(1:ndims,1:ngrid_act), var(1:ngrid_act,1:twotondim,1:nvar_mhd+1), ref(1:ngrid_act))
        !
        ! loop over domains
        do j=1,nboundary+ncpu
          ncache = ngrid_file(j,ilevel)
          !
          ! read relevant AMR tree data
          if (ncache > 0) then
            allocate(index_grid(1:ncache), iig(1:ncache))
            if (j == icpu) then
              read(ilun) index_grid ! read grid index
            else
              read(ilun) ! SKIP
            endif
            read(ilun) ! SKIP next index
            read(ilun) ! SKIP prev index
            ! read grid centre
            do i=1,ndims
              if (j == icpu) then
                read(ilun) xg(:,i)
              else
                read(ilun) ! SKIP
              endif
            enddo
            if (j == icpu) then
              read(ilun) iig ! read father index
              father(index_grid(1:ncache)) = iig(1:ncache)
            else
              read(ilun) ! SKIP
            endif
            do i=1,2*ndims
              if (j == icpu) then
                read(ilun) iig ! read neighbour index
                neighbour(index_grid(1:ncache),i) = iig(1:ncache)
              else
                read(ilun) ! SKIP
              endif
            enddo
            do i=1,twotondim
              if (j == icpu) then
                read(ilun) son(:,i) ! read son index
              else
                read(ilun) ! SKIP
              endif
            enddo
            do i=1,twotondim
              read(ilun) ! SKIP cpu map
            enddo
            do i=1,twotondim
              read(ilun) ! SKIP refinement map
            enddo
            deallocate(index_grid, iig)
          endif
          !
          ! read HYDRO data
          read(ilun2) ! SKIP ilevel
          read(ilun2) ! SKIP ncache
          if (ncache > 0) then
            do i=1,twotondim
              do ivar=1,nvar_mhd
                if (j == icpu) then
                  read(ilun2) var(:,i,ivar)
                else
                  read(ilun2) ! SKIP
                endif
              enddo
            enddo
          endif
        enddo ! end loop j=1,nboundary+ncpu
        ncellsum(ilevel,k) = 0
        position(:) = 0.0
        if (ngrid_act > 0) then
          if (ilevel >= self%rootlvl .and. ilevel <= self%maxlvl) then
            do l=1,ngrid_act
              self%npatch = self%npatch+1
              n_level(ilevel) = n_level(ilevel)+1
              !-----------------------------------------------------------------
              ! Save oct info in oct_list, optionally including oct data
              !-----------------------------------------------------------------
              allocate (oct)
              call io%bits_mem (8, oct_size, 'oct')
              position(:) = xg(l,:)
              ds = 0.5d0**ilevel
              size = self%n_cell*ds
              oct%key(1) = ilevel
              oct%key(2:4) = position/size
              oct%ii = floor((position-ds-oct%key(2:4)*size)/ds)
              if (verbose > 2) then
                write (stdlog,'(a,3f10.4,4i4)') &
                  'oct: position,key=', position, oct%key
              end if
              if ((keep .or. load) .and. ilevel <= shared%levelmax) then
                if (.not.allocated(oct%var)) then
                  allocate (oct%var(8,nvar_mhd+1))
                  call io%bits_mem (&
                    storage_size(oct%var), product(shape(oct%var)), 'oct')
                  oct%var = var(l,:,:)
                  call self%cartesian_init (oct%key, patch)
                  call self%add_to (patch, oct)
                  oct%u_max = self%phase_speed (var(l,:,:))
                  u_max = max(u_max,oct%u_max)
                  self%u_min = min(self%u_min,oct%u_max)
                  dll_node => oct
                  call self%oct_list%append (dll_node)
                  deallocate (oct%var)
                end if
              end if
              !-----------------------------------------------------------------
              ! Keep only octs at or below the LEVELMAX input parameter
              !-----------------------------------------------------------------
              !if (.not.keep .and. allocated(oct%var)) then
                !deallocate (oct%var)
              !end if
              !-----------------------------------------------------------------
              n_grid = n_grid + 1
            end do
            if (verbose > 1) &
              write(stdlog,'(1x,a,i10,i6,i4,1p,e12.3)') &
                'n_oct, rank, level, u_max:',self%npatch, k, ilevel, u_max
          endif
          deallocate(xg,son,x,var,ref)
        endif
      end do ! end loop ilevel=1,self%maxlvl
      self%u_max = u_max
      if (verbose > 0) &
        write(stdlog,'(1x,a,i10,i6,1p,e12.3)') &
          'n_oct, rank(in), u_max:',self%npatch, k, u_max
      close(ilun)  ! close `amr_` file
      close(ilun2) ! close `hydro_` file
      if (verbose > 2) &
        write (stdlog,'(4(a,i8,3x),a,1p,e9.2)') &
        ' n_active:',n_active,'n_passive:',n_passive,'n_onrank:',n_onrank, &
        ' n_offrank:',n_offrank
      !
      call mpi%barrier ('ncpu loop')
    enddo ! end loop k=1,ncpu_read
    write (stdout,'(1x,a,1p,e12.3)')   'u_max =', u_max
    write (stdout,'(1x,a,i12.1)')      'n_oct =', self%npatch
    if (verbose > 1) then
      write (stdout,'(a,2(i6.1,1x))') ' n_d_min, n_tt_min =',n_d_min,n_tt_min
      write (stdout,'(a,i11.1)')      ' n_frozen = ', n_frozen
    end if
  !$omp end master
  !$omp barrier
  !
  ! set scale factors (other scale factors will be calculated)
  scaling%l = self%scale_l
  scaling%t = self%scale_t
  scaling%d = self%scale_d
  !
  ! output some diagnostics to the terminal.
  allocate (fill(self%maxlvl))
  do ilevel=nlevelmin,self%maxlvl
    fill(ilevel)=sum(ncellsum(ilevel,:))/2.0**(3*(ilevel-1))
  enddo
  if (verbose > 0) then
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  end if
  !-----------------------------------------------------------------------------
  total_cost=0.0
  self%minlvl = nlevelmin
  self%maxlvl = nlevelmax
  do ilevel=self%minlvl,self%maxlvl
    ds = 0.5_8**ilevel
    weight = 8./ds
    patch_cost = n_level(ilevel)*weight
    total_cost = total_cost + patch_cost
  end do
  write(stdout,*) ' level  n_patch  cost(%)'
  do ilevel=self%minlvl,self%maxlvl
    ds = 0.5_8**ilevel
    weight = 8./ds
    patch_cost = n_level(ilevel)*weight
    write(stdout,'(i7,i9,f9.2)') &
      ilevel, n_level(ilevel), patch_cost*100./total_cost
  end do
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    if (cut_out_size > 0 .and. self%centre_on_sink) then
      write (*,'(a,1p,3(1e14.7,2x))') 'pmin (AU) =', &
        (pmin-star%position)*scaling%l/cgs%au
      write (*,'(a,1p,3(1e14.7,2x))')'pmax (AU) =', &
        (pmax-star%position)*scaling%l/cgs%au
    end if
    write (*,'(" fill (by level):",/,(1p,10g11.2))') fill(nlevelmin:self%maxlvl)
  end if
  deallocate (fill)
  deallocate(ngrid_file, ngrid_level, cpu_list, cpu_read, xc, father, &
             neighbour, levels, ncellsum)
  !-----------------------------------------------------------------------------
  ! If not doing direct_load, then setup the hash tables now
  !-----------------------------------------------------------------------------
  if (.not.direct_load) then
    !hash_table_size=self%oct_list%n/50
    call self%u_max_table%init (hash_table_size, ndim=4)
    call self%patch_table%init (hash_table_size, ndim=4)
  end if
  call memtest_delay
  call trace%end (itimer)
END SUBROUTINE read_snapshot

!===============================================================================
!===============================================================================
SUBROUTINE memtest_delay
  if (do_memtest) then
    write (stdout,*) 'create_patches: exit in 5 sec'
    call omp%delay (5e3)
  end if
END SUBROUTINE memtest_delay

!===============================================================================
!> Lifted directly from `amr2cell.f90`/`hilbert3d_qdp` (in the Copenhagen
!> RAMSES repository).
!===============================================================================
SUBROUTINE hilbert3d(x,y,z,order,bit_length,npoint)
  integer     , intent(in):: bit_length, npoint
  integer     , intent(in), dimension(1:npoint):: x, y, z
  integer:: i, ip, cstate, nstate, b0, b1, b2, sdigit, hdigit
  real(kind=16),intent(out), dimension(1:npoint):: order
  logical, dimension(0:3*bit_length-1):: i_bit_mask
  logical, dimension(0:1*bit_length-1):: x_bit_mask, y_bit_mask, z_bit_mask
  integer, dimension(0:7,0:1,0:11), parameter:: &
   state_diagram= reshape( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )
  !.............................................................................
  if(bit_length > bit_size(bit_length)) then
     write (stdout,*) 'Maximum bit length = ', bit_size(bit_length)
     write (stdout,*) 'ERROR in hilbert3d'
     stop
  endif
  do ip=1,npoint
     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i) = btest(x(ip), i)
        y_bit_mask(i) = btest(y(ip), i)
        z_bit_mask(i) = btest(z(ip), i)
     enddo
     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2) = x_bit_mask(i)
        i_bit_mask(3*i+1) = y_bit_mask(i)
        i_bit_mask(3*i  ) = z_bit_mask(i)
     enddo
     ! build Hilbert ordering using state diagram
     cstate = 0
     do i=bit_length-1,0,-1
        b2 = merge(1, 0, i_bit_mask(3*i+2))
        b1 = merge(1, 0, i_bit_mask(3*i+1))
        b0 = merge(1, 0, i_bit_mask(3*i  ))
        sdigit = b2 * 4 + b1 * 2 + b0
        nstate = state_diagram(sdigit, 0, cstate)
        hdigit = state_diagram(sdigit, 1, cstate)
        i_bit_mask(3*i+2) = btest(hdigit,2)
        i_bit_mask(3*i+1) = btest(hdigit,1)
        i_bit_mask(3*i  ) = btest(hdigit,0)
        cstate = nstate
     enddo
     ! save Hilbert key as double/quad precision real
     order(ip) = 0.0_16
     do i=0,3*bit_length-1
        order(ip) = order(ip) + merge(1.0_16, 0.0_16, i_bit_mask(i)) * 2.0_16**i
     enddo
  enddo
END SUBROUTINE hilbert3d

!===============================================================================
!> Compute max phase speed directly from RAMSES vars
!===============================================================================
FUNCTION phase_speed (self, var) RESULT (umax)
  class(read_external_t):: self
  real(kind=8):: var(:,:)
  real(kind=8):: d, ux, uy, uz, bx, by, bz, pg, cf, ua
  real:: umax
  integer:: i
  umax = 0.0
  do i=1,size(var,1)
    d  = var(i,1)
    ux = var(i,2)
    uy = var(i,3)
    uz = var(i,4)
    bx = 0.5*(var(i,5)+var(i,8))
    by = 0.5*(var(i,6)+var(i,9))
    bz = 0.5*(var(i,7)+var(i,10))
    pg = var(i,11)
    cf = sqrt((bx**2+by**2+bz**2+eos%gamma*pg)/d)
    ua = sqrt(ux**2+uy**2+uz**2)
    umax = max(umax,ua+cf)
  end do
END FUNCTION phase_speed

!===============================================================================
!> Run through all RAMSES oct and insert their values into patches.  For each
!> patch introduced at level L, there must be one also at the corresponding
!> location at level L-1.  Carrying this through down to the minlvl, one can
!> be sure to have support for prolongation, if only one proceeds from the
!> bottom up with prolongation.
!===============================================================================
SUBROUTINE create_patches (self)
  class(read_external_t):: self
  class(dll_node_t), pointer:: oct
  class(experiment_t), pointer:: patch
  class(*), pointer:: ptr
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%create_patches', itimer=itimer)
  call io%header ('-create_patches: begin', left=25)
  oct => self%oct_list%head
  do while (associated(oct))
    select type (oct)
    class is (oct_t)
    call self%cartesian_init (oct%key, patch)
    call self%add_to (patch, oct)
    end select
    oct => oct%next
  end do
  if (verbose > 0) then
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  end if
  call memtest_delay
  call io%header ('-create_patches: end', left=25)
  call trace%end (itimer)
END SUBROUTINE create_patches

!===============================================================================
!> Add oct values to a patch, located via a hash table.  If the patch does not
!> yet exist, create and initialize it, and add the address into the table.
!===============================================================================
SUBROUTINE add_to (self, patch, oct)
  class(read_external_t):: self
  class(experiment_t), pointer:: patch
  class(oct_t):: oct
  !.............................................................................
  real:: enth, d, ux, uy, uz, bx, by, bz, e
  integer:: i1, i2, i3, j1, j2, j3, j, k1, k2, k3, level
  logical:: new
  !-----------------------------------------------------------------------------
  ! Insert oct values into patch%mem.  FIXME: HLLD solver assumed
  !-----------------------------------------------------------------------------
  enth = 1d0/(eos%gamma-1d0)
  i1 = patch%li(1) + oct%ii(1)
  i2 = patch%li(2) + oct%ii(2)
  i3 = patch%li(3) + oct%ii(3)
  k1 =           1 + oct%ii(1)
  k2 =           1 + oct%ii(2)
  k3 =           1 + oct%ii(3)
  level = oct%key(1)
  j = 1
  do j3=0,1
  do j2=0,1
  do j1=0,1
    d  = oct%var(j,1)
    ux = oct%var(j,2)
    uy = oct%var(j,3)
    uz = oct%var(j,4)
    !bx = (oct%var(j,5) + oct%var(j,8 ))*0.5
    !by = (oct%var(j,6) + oct%var(j,9 ))*0.5
    !bz = (oct%var(j,7) + oct%var(j,10))*0.5
    e  = oct%var(j,11)*enth
    !e  = e + 0.5*(d*(ux**2+uy**2+uz**2) + bx**2 + by**2 + bz**2)
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%d ,1,1) = d
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%e ,1,1) = e
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%px,1,1) = d*ux
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%py,1,1) = d*uy
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%pz,1,1) = d*uz
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%bx,1,1) = oct%var(j,5)
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%by,1,1) = oct%var(j,6)
    patch%mem(i1+j1,i2+j2,i3+j3,patch%idx%bz,1,1) = oct%var(j,7)
    patch%filled(i1+j1,i2+j2,i3+j3) = level
    if (verbose > 3) &
      write(stdlog,'("add_to:",i7,i4,2x,3i4,1p,e11.2)') &
        patch%id, patch%level, i1+j1, i2+j2, i3+j3, d
    j = j+1
  end do
  end do
  end do
END SUBROUTINE add_to

!===============================================================================
!> Add a patch to the task list
!===============================================================================
SUBROUTINE add_to_task_list (self, patch)
  class(read_external_t):: self
  class(experiment_t), pointer:: patch
  class(link_t), pointer:: link
  integer, save:: i=0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('add_to_task_list', itimer=itimer)
  allocate (link)
  link%task => patch
  patch%link => link
  call self%task_list%prepend_link (link)
  if (verbose > 3) then
    write (stdlog,'(a,i7,2x,4i6)') &
      'add_to_task_list: id, key =', patch%id, patch%key
  else if (verbose > 0) then
    i = i+1
    if (mod(i,5000)==0) &
      write (stdlog,'(a,i7,2x,4i6)') &
        'add_to_task_list: id, key =', patch%id, patch%key
  end if
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE add_to_task_list

!===============================================================================
!> Count the costs of RAMSES and DISPATCH
!===============================================================================
SUBROUTINE count_costs (self)
  class(read_external_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: task
  class(experiment_t), pointer:: patch
  real(kind=8):: total_cost, max_cost, total_weight
  real(kind=4), pointer:: value
  real(kind=4):: ratio, weight, u_max, ds, oct_cost, patch_cost, c_bin
  integer, allocatable:: n_level(:)
  logical:: new
  integer:: n(2), level, bin, i
  !-----------------------------------------------------------------------------
  call io%header ('-count_costs: begin', left=25)
  allocate (self%u_bins(self%n_bins))
  self%u_bins = 0
  n = 0
  u_max = 0.0
  total_weight = 0d0
  !-----------------------------------------------------------------------------
  ! First loop over octs, to create hash table entries with patch costs
  !-----------------------------------------------------------------------------
  c_bin = 1.0/alog(self%u_max/self%u_min)*(self%n_bins-1)
  task => self%oct_list%head
  do while (associated(task))
    select type (task)
    class is (oct_t)
      !-------------------------------------------------------------------------
      ! Search for max phase speed over all octs, while summing up weights
      !-------------------------------------------------------------------------
      u_max = max(u_max,task%u_max)
      bin = floor(alog(task%u_max/self%u_min)*c_bin) + 1
      self%u_bins(bin) = self%u_bins(bin) + 1
      level = task%key(1)
      ds = 0.5_8**level
      weight = 2.0**3/ds
      total_weight = total_weight + weight
      oct_cost = task%u_max*weight
      call self%u_max_table_get (task%key, value, new)
      n(1) = n(1)+1
      if (mod(n(1),10000000)==0) then
        write(stdout,*) n(1)
      end if
      if (new) then
        n(2) = n(2)+1
      end if
      value = max(value,task%u_max)
    end select
    task => task%next
  end do
  if (verbose > 0) &
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  open (io_unit%tmp, file=trim(io%run_dir)//'u_bins.dat', form='unformatted', &
    status='unknown')
  write (io_unit%tmp) self%n_bins
  write (io_unit%tmp) alog10(self%u_min), alog10(self%u_max)
  write (io_unit%tmp) self%u_bins
  close (io_unit%tmp)
  write (stdout,*) &
    'read_external_t%count_costs: n, n(1)/n(2) =', n, n(1)/n(2)
  max_cost = u_max*total_weight
  !-----------------------------------------------------------------------------
  ! Second loop over octs, to sum up total patch costs
  !-----------------------------------------------------------------------------
  total_cost = 0d0
  allocate (n_level(self%rootlvl:self%maxlvl))
  n_level = 0
  task => self%oct_list%head
  do while (associated(task))
    select type (task)
    class is (oct_t)
      level = task%key(1)
      if (level >= self%rootlvl) then
        ds = 0.5_8**level
        weight = self%n_cell**3/ds
        call self%u_max_table_get (task%key, value, new)
        if (mod(n(1),10000000)==0) then
          write(stdout,*) n(1)
        end if
        if (value > 0.0) then
          total_cost = total_cost + value*weight
          n_level(level) = n_level(level)+1
          value = 0.0
        end if
      end if
    end select
    task => task%next
  end do
  ratio = max_cost/total_cost
  if (verbose > 0) &
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  !-----------------------------------------------------------------------------
  ! Print a table with level costs for octs, and for patches
  !-----------------------------------------------------------------------------
  total_cost=0.0
  do level=self%rootlvl,self%maxlvl
    ds = 0.5_8**level
    weight = self%n_cell**3/ds
    patch_cost = n_level(level)*weight
    total_cost = total_cost + patch_cost
  end do
  write(stdout,*) ' level  n_patch  cost(%)'
  do level=self%minlvl,self%maxlvl
    ds = 0.5_8**level
    weight = self%n_cell**3/ds
    patch_cost = n_level(level)*weight
    write(stdout,'(i7,i9,f9.2)') level, n_level(level), patch_cost*100./total_cost
  end do
  write (stdout,*) 'read_external_t%count_costs:        ratio =', ratio
  if (verbose > 0) & 
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  call io%header ('-count_costs: end', left=25)
  !-----------------------------------------------------------------------------
END SUBROUTINE count_costs

!===============================================================================
!> Get a hash table pointer to a double precision cost variable
!===============================================================================
SUBROUTINE u_max_table_get (self, key, value, new)
  class(read_external_t):: self
  integer:: key(4)
  real(kind=4), pointer:: value
  class(*), pointer:: ptr
  logical, optional:: new
  !-----------------------------------------------------------------------------
  call self%u_max_table%get_ptr (key, ptr)
  if (.not. associated(ptr)) then
    allocate (value)
    ptr => value
    call self%u_max_table%set_ptr (key, ptr)
    value = 0d0
    if (present(new)) then
      new = .true.
    end if
  else if (present(new)) then
    select type (ptr)
    type is (real)
    value => ptr
    end select
    new = .false.
  end if
END SUBROUTINE u_max_table_get

!===============================================================================
!> Convenience function to calculate the entropy (because the STAGGER solver in
!> DISPATCH uses the entropy rather than the internal energy/pressure or total
!> energy.)
!===============================================================================
FUNCTION energy (self, vars, patch) RESULT (e)
  class(read_external_t):: self
  real(kind=8):: e, bx, by, bz
  real(kind=8), dimension(:):: vars
  class(experiment_t), pointer:: patch
  !.............................................................................
  if (eos%gamma==1d0) then
    e = vars(11)
  else if (patch%kind(1:6)=='ramses') then
    bx = 0.5_8*(vars(5)+vars(8))
    by = 0.5_8*(vars(6)+vars(9))
    bz = 0.5_8*(vars(7)+vars(10))
    e = vars(11)/(eos%gamma-1d0) &
      + 0.5d0*vars(1)*(vars(2)**2+vars(3)**2+vars(4)**2) &
      + 0.5d0        *(bx**2+by**2+bz**2)
  else if (patch%kind(1:9)=='stagger2e') then
    e = vars(1)*(log(vars(11))-log(vars(1))*eos%gamma)/(eos%gamma-1.0)
  else
    e = vars(11)/(eos%gamma-1d0)
  end if
END FUNCTION energy

!===============================================================================
!> Check, from the bottom up that every patch has support, on all sides; i.e.,
!> if no level L nbor exists, then a level L-1 must exist, and will be created
!> if it doesn't
!===============================================================================
SUBROUTINE check_support (self)
  class(read_external_t):: self
  integer:: level, i1, i2, i3, key(4), n
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task, nbtask
  class(experiment_t), pointer:: patch
  class(*), pointer:: ptr
  real(8):: position(3), size(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%check_support', itimer=itimer)
  call io%header ('-check_support: begin', left=25)
  call timer%tic
  do level=self%rootlvl+1,self%maxlvl
    link => self%task_list%head
    do while (associated(link))
      if (link%task%level == level) then
        task => link%task
        select type (task)
        class is (patch_t)
          size = task%size
          do i3=-1,1
          do i2=-1,1
          do i1=-1,1
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            !-------------------------------------------------------------------
            ! Check if a level L nbor exists at this point
            !-------------------------------------------------------------------
            position = task%position + [i1,i2,i3]*size
            key = task_key (task%level, task%box, position, size)
            call self%patch_table%get_ptr (key, ptr)
            if (associated(ptr)) then
              if (verbose > 0) then
                select type (ptr)
                class is (experiment_t)
                  if (verbose > 0) &
                    write (stdlog,'(a,2(i6,":",i2.2),2x,4i3,2x,3i3)') &
                    'read_external_t%check_support: id,nbor=', &
                    task%id, task%level, ptr%id, ptr%level, key, [i1,i2,i3]
                end select
              end if
            else
              !-----------------------------------------------------------------
              ! Check if level L-1 nbor exist that covers the place, and create
              ! one if it does not exist
              !-----------------------------------------------------------------
              key = task_key (task%level-1, task%box, position, 2d0*size)
              call self%patch_table%get_ptr (key, ptr)
              if (associated(ptr)) then
                if (verbose > 0) then
                  select type (ptr)
                  class is (experiment_t)
                    if (verbose > 0) &
                      write (stdlog,'(a,2(i6,":",i2.2),2x,4i3,2x,3i3)') &
                      'read_external_t%check_support: id,nbor=', &
                      task%id, task%level, ptr%id, ptr%level, key, [i1,i2,i3]
                  end select
                end if
              else
                call self%cartesian_init (key, patch)
                if (verbose > 1) &
                  write (stdlog,'(a,2(i6,":",i2.2),2x,4i3,2x,3i3)') &
                    'read_external_t%check_support: id,new_nbor=', &
                    task%id, task%level, patch%id, patch%level, key, [i1,i2,i3]
              end if
            end if
          end do
          end do
          end do
        end select
      end if
      link => link%next
    end do
  end do
  if (verbose > 0) &
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  call timer%toc ('time to check_support', 1)
  call io%header ('-check_support: end', left=25)
  call trace%end (itimer)
END SUBROUTINE check_support

!===============================================================================
!> Compute the hash key, given the task leve, periodic box size, position and
!> size, with periodic wrapping
!===============================================================================
FUNCTION task_key (level, box, position, size) RESULT (key)
  integer:: level
  real(8), dimension(3):: box, position, size
  integer:: n(3), key(4)
  !-----------------------------------------------------------------------------
  key(1) = level
  key(2:4) = floor(position/size)
  n = nint(box/size)
  key(2:4) = modulo(key(2:4),n)
END FUNCTION task_key

!===============================================================================
!> Starting one level about levelmin, prolong into those patches that are not
!> properly filled by oct values
!===============================================================================
SUBROUTINE prolong (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  class(*), pointer:: patch
  real:: dmin, dmin_all
  real(8):: size(3)
  integer:: level, n, n0, i, key(4), l(3), u(3), n_task
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%prolong')
  call io%header ('-prolong: begin', left=25)
  call timer%tic
  n = 0; n0 = 0
  n_task = 0
  do level=self%minlvl+1,self%maxlvl
    dmin_all=huge(1.0)
    link => self%task_list%head
    do while (associated(link))
      if (link%task%level == level) then
        associate (task => link%task)
        select type (task)
        class is (patch_t)
          l = task%mesh%li
          u = task%mesh%ui
          if (any(task%filled(l(1):u(1),l(2):u(2),l(3):u(3)) == -1)) then
            n = n+1
            call self%download_link (link)
            task%filled = task%level-1
          end if
          dmin = task%fminval (task%idx%d)
          dmin_all = min(dmin,dmin_all)
          if (dmin <= 0.0) then
            n0 = n0+1
            write(stderr,*) 'prolong: ERROR, id, dmin, n=', task%id, dmin, n0
          end if
          if (verbose > 1 .or. (verbose > 0 .and. mod(n,5000)==0)) &
            write(stdlog,'(a,i7,i4,1p,2e11.2)') &
            'prolong: id, level, dmin, dmin_all =', task%id, task%level, dmin, dmin_all
        end select
        end associate
        n_task = n_task+1
      end if
      !if (download%n_tot > 1e7) then
      !  write (stdout,'(a,2i10,2p,f6.1,"%")') &
      !    'prolong:', download%n_tot, download%n_skip, real(n_task)/self%task_list%n
      !  download%n_tot = 0
      !  download%n_skip = 0
      !end if
      link => link%next
    end do
  end do
  if (verbose > 0) then
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  end if
  call timer%toc ('time to fill missing cells', 1)
  call io%header ('-prolong: end', left=25)
  call trace%end ()
END SUBROUTINE prolong

!===============================================================================
!> Fill guard zones for all tasks
!===============================================================================
SUBROUTINE dnload (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%dnload')
  call io%header ('-dnload: begin', left=25)
  call timer%tic ()
  link => self%task_list%head
  do while (associated(link))
    call link%task%dnload
    link => link%next
  end do
  call timer%toc ('time to dnload', 1)
  call io%header ('-dnload: end', left=25)
  call trace%end ()
END SUBROUTINE dnload

!===============================================================================
!> Find the patch that contains the link task and prolong from it
!===============================================================================
SUBROUTINE download_link (self, link)
  class(read_external_t):: self
  class(link_t), pointer:: link
  class(patch_t), pointer:: target, source
  class(*), pointer:: ptr
  integer:: jt(2), key(4)
  real:: pt(2)
  logical:: ok
  !-----------------------------------------------------------------------------
  associate (task=>link%task)
  select type(task)
  class is (experiment_t)
    target => task
    key = task_key (task%level-1, task%box, task%position, 2d0*task%size)
    call self%patch_table%get_ptr (key, ptr)
    ok = .false.
    if (associated(ptr)) then
      select type(ptr)
      class is (experiment_t)
        source => ptr
        jt = [1,1]
        pt = [1.0,0.0]
        download%check_filled = .true.
        call download%different (target, source, jt, pt, all_cells=.true.)
        ok = .true.
      end select
    else
      write (stderr,'(a,4i5)') 'read_external_t%download_link: key not found =', key
    end if
    if (verbose > 3) then
      write (stdlog,'(a,4i5,l3)') 'read_external_t%download_link: key =', key, ok
    end if
  end select
  end associate
END SUBROUTINE download_link

!===============================================================================
!===============================================================================
SUBROUTINE download_link2 (self, link)
  class(read_external_t):: self
  class(nbor_t), pointer:: nbor
  class(link_t), pointer:: link
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real:: pt(2)
  !-----------------------------------------------------------------------------
  nbor => link%nbor
  associate (task=>link%task)
  select type(task)
  class is (patch_t)
    target => task
    do while (associated(nbor))
      if (nbor%task%level == link%task%level-1) then
        associate(nbtask => nbor%task)
        select type (nbtask)
        class is (patch_t)
          source => nbtask
          jt = [1,1]
          pt = [1.0,0.0]
          if (nbtask%contains(task%position)) then
            call download%different (target, source, jt, pt, all_cells=.true.)
          end if
        end select
        end associate
      end if
      nbor => nbor%next
    end do
  end select
  end associate
END SUBROUTINE download_link2

!===============================================================================
!===============================================================================
SUBROUTINE check_filled (self)
  class(read_external_t):: self
  integer:: level, gap, nok, l(3), u(3)
  class(nbor_t), pointer:: nbor
  class(link_t), pointer:: link
  class(task_t), pointer:: task, nbtask
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%check_filled', itimer=itimer)
  call io%header ('-check_filled: begin', left=25)
  nok = 0
  do level=self%maxlvl,self%minlvl+1,-1
    link => self%task_list%head
    do while (associated(link))
      if (link%task%level == level) then
        task => link%task
        select type (task)
        class is (patch_t)
        l = task%mesh%li
        u = task%mesh%ui
        if (any(task%filled(l(1):u(1),l(2):u(2),l(3):u(3)) == -1)) then
          gap = 99
          nbor => link%nbor
          do while (associated(nbor))
            nbtask => nbor%task
            select type (nbtask)
            class is (patch_t)
            if (nbtask%level < task%level .and. &
                nbtask%contains (task%position)) then
              gap = min(gap,task%level-nbtask%level)
            end if
            end select
            nbor => nbor%next
          end do
          if (gap > 1) then
            nok = nok+1
            if (verbose > 0 .and. mod(nok,5000)==0) &
              write (stdout,'(a,i4,i7,i4,i7)') &
                'read_external_t%check_filled: not filled level, id, gap, nok =', &
                  task%level, task%id, gap, nok
            call create_support (task)
          end if
        end if
        end select
      end if
      link => link%next
    end do
  end do
  call self%task_list%count_status()
  if (verbose > 0) then
    write(stdout,*) 'check_filled: nok =', nok
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  end if
  call io%header ('-check_filled: end', left=25)
  call trace%end (itimer)
contains
  !=============================================================================
  !> Create a new patch as support to task
  !=============================================================================
  subroutine create_support (task)
    class(patch_t):: task
    class(experiment_t), pointer:: new
    integer:: key(4)
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    !...........................................................................
    call trace%begin ('create_support', itimer=itimer)
    !---------------------------------------------------------------------------
    ! Create a task at the overlapping location on level L-1
    !---------------------------------------------------------------------------
    key = task_key (task%level-1, task%box, task%position, 2d0*task%size)
    call self%cartesian_init (key, new)
    nbor => new%link%nbor
    if (verbose > 1) then
      do while (associated(nbor))
        write (stdout,'(a,i6,i4,i6,i4)') 'create_support: id,lelvel,nbor,level', &
          new%link%task%id, new%link%task%level, nbor%task%id, nbor%task%level
        nbor => nbor%next
      end do
    end if
    call download%download_link (new%link, all_cells=.true.)
    new%filled = new%level-1
    call trace%end (itimer)
  end subroutine create_support
END SUBROUTINE check_filled

!===============================================================================
!> Read sink particle information; the zoom-in sink is the last one, generally
!===============================================================================
 SUBROUTINE read_sink (self)
  class(read_external_t):: self
  integer, parameter:: sink_unit=9
  character(len=80):: fname
  logical:: do_translate
  integer:: center_star, version, nstar, szstar, sizeofstar, i, nstars, nmet, dp1
  real(4):: t,beta_create,dtold
  real(8):: t8,tflush,r_translate(3),v_translate(3),x_refine,y_refine,z_refine
  type(sink_t), pointer:: sink
  !.............................................................................
  call trace%begin('read_external_t%read_sink')
  write(stdout,*) '=================================================='
  write (fname,'(a,a,i5.5,a)') trim(self%path), '/output_', self%nout, '/stars_output.dat'
  write(stdout,*) fname
  open (sink_unit, file=fname,form='unformatted', status='unknown')
  read (sink_unit) version,nstars,t,szstar,sizeofstar
  write(stdout,*) 'versions =',version
  read (sink_unit) dp1,nmet,t8
  write(stdout,*) 'dp, nmetal =',dp1,nmet
  write(stdout,*) 't8 =',t8
  write(stdout,*) 'nstars =',nstars
  read (sink_unit) do_translate, center_star, r_translate, v_translate,  &
                   x_refine, y_refine, z_refine
  write(stdout,*) 'center_star',center_star
  allocate (self%sink(nstars))
  do i=1,nstars
    sink => self%sink(i)
    read (sink_unit) sink%x,sink%y,sink%z,sink%px,sink%py,&
                     sink%pz,sink%dpx,sink%dpy,sink%dpz,sink%mass,&
                     sink%dmass,sink%potential,sink%density,sink%t_create,&
                     sink%t_explode,sink%metal(1:nmetal),sink%dmetal(1:nmetal),&
                     sink%gid,sink%father,sink%ilevel,sink%old_ilevel, &
                     sink%owner,sink%has_moved,sink%beta_create,sink%dtold
  enddo
  close(sink_unit)
  sink => self%sink(center_star)
  star%position = [sink%x,sink%y,sink%z]
  star%mass = sink%mass
  write(stdout,*) 'type(sink_t) :: star => center_star'
  write(stdout,'(a,"["2(f14.9,","),f14.9,"]",2x,a,1i3.1)') &
    ' pos = ',star%position,', level = ',sink%ilevel
  write(stdout,*) '=================================================='
  call trace%end()
END SUBROUTINE read_sink

!===============================================================================
!> Evaluate update cost for all patches
!===============================================================================
SUBROUTINE evaluate (self, task_list)
  class (read_external_t):: self
  type (task_list_t), pointer:: task_list
  type (link_t), pointer:: link
  class (task_t), pointer:: task
  real, dimension(:,:,:,:), pointer:: mem
  integer:: l(3), u(3), i, ix, iy, iz, nskip, level, n, key(4), n_cost
  real:: ds, d, e, pg, ux, uy, uz, umax, bx, by, bz, p(3), cs, cf, ua, ub, ctmax, dsmin, ct
  real(8):: cost1, cost2, cost3, max_cost, sum_cost, cost1_l, cost2_l
  real(4), pointer:: cost
  !.............................................................................
  call trace%begin('read_external_t%evaluate')
  call io%header ('-evaluate: begin', left=25)
  cost1 = 0d0
  cost2 = 0d0
  cost3 = 0d0
  ctmax = 0d0
  level = 0
  dsmin = 0.5d0**self%maxlvl
  link => task_list%head
  i = 0
  do while (associated(link))
    task => link%task
    i = i+1
    if (mod(i,100000)==0) &
      write (stdout,*) i
    select type (task)
    class is (experiment_t)
      call self%u_max_table_get (task%key, cost)
      if (verbose > 1) &
        write (stdout,'(i6,o4,4(2x,3i4))') task%id, task%is_clear(bits%frozen), &
          task%mesh%li, task%mesh%ui, task%mesh%lo, task%mesh%uo
      if (task%is_clear(bits%frozen)) then
        if (verbose > 3) &
          write (stdout,*) 'PATCH:', task%id, &
            task%level, task%mesh%n, task%mem(3,3,3,1,task%it,1), task%it
        mem => task%mem(:,:,:,:,task%it,1)
        l = task%mesh%li
        u = task%mesh%ui
        ds = task%mesh(1)%d
        d = sum(mem(l(1):u(1),l(2):u(2),l(3):u(3),1))
        if (d <= 0.0) then
          write (stdout,*) 'SKIPPING', task%id
          cycle
        end if
        do i=1,3
          p(i) = sum(mem(l(1):u(1),l(2):u(2),l(3):u(3),task%idx%px+i-1))/d
        end do
        !-------------------------------------------------------------------------
        ! Evaluate the Courant condition for the oct and compute the cost
        !-------------------------------------------------------------------------
        nskip = 0
        cost1_l = 0.0
        cost2_l = 0.0
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              d =  mem(ix,iy,iz,task%idx%d)
              e =  mem(ix,iy,iz,task%idx%e)
              if (d <= 0.0) then
                nskip = nskip+1
                cycle
              end if
              ux = mem(ix,iy,iz,task%idx%px)/d
              uy = mem(ix,iy,iz,task%idx%py)/d
              uz = mem(ix,iy,iz,task%idx%pz)/d
              bx = mem(ix,iy,iz,task%idx%bx)
              by = mem(ix,iy,iz,task%idx%by)
              bz = mem(ix,iy,iz,task%idx%bz)
              if (task%kind(1:6)=='ramses') then
                !pg = (e-0.5*d*((ux**2+uy**2+uz**2) &
                !              +(bx**2+by**2+bz**2)))*(eos%gamma-1d0)
                pg = (e-0.5*d*((ux**2+uy**2+uz**2)))*(eos%gamma-1d0)
              else if (task%kind(1:9)=='stagger2e') then
                pg = d**eos%gamma*exp(e/d*(eos%gamma-1.0))
              else
                pg = e*(eos%gamma-1d0)
              end if
              cs = sqrt(eos%gamma*pg/d)
              cf = sqrt((bx**2+by**2+bz**2+eos%gamma*pg)/d)
              ua = sqrt(ux**2+uy**2+uz**2)
              ux = ux-p(1)
              uy = uy-p(2)
              uz = uz-p(3)
              ub = sqrt(ux**2+uy**2+uz**2)
              ct = (ua+cf)/dsmin
              cost1_l = max(cost1_l,(ua+cf)/ds)
              cost2_l = max(cost1_l,(ub+cf)/ds)
              cost3 = cost3 + 1d0/ds
              if (ct > ctmax) then
                ctmax = ct
                umax = ua+cf
                level = task%level
              end if
            end do
          end do
        end do
        if (nskip > 0) then
          write (stdout,*) 'SKIPPING PATCH', task%id, nskip, product(u-l+1)
        else
          cost  = cost  + cost1_l*product(u-l+1)
          cost1 = cost1 + cost1_l*product(u-l+1)
          cost2 = cost2 + cost2_l*product(u-l+1)
        end if
      end if
    end select
    link => link%next
  end do
  cost3 = cost3*umax
  write(stdout,'(a,i4,1p,e12.4,i8)') ' level, umax:', level, umax
  write(stdout,'(a,1p,3e12.4)') ' costs:', cost1, cost2, cost3
  write(stdout,'(a,1p,3e12.4)') ' costs ratios:     ', &
    cost3/max(tiny(1e0),cost1), cost3/max(tiny(1e0),cost2)
  !-----------------------------------------------------------------------------
  ! Sum up the cost over the larger, 16x16x16 DISPATCH patches by running threw
  ! the task list, making sure to count each cost only once
  !-----------------------------------------------------------------------------
  link => task_list%head
  level_cost = 0.0
  max_cost = 0.0
  sum_cost = 0.0
  n_cost = 0
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
      call self%u_max_table_get (task%key, cost)
      if (cost > 0.0) then
        !level_cost(task%level) = level_cost(task%level) + cost
        max_cost = max(cost,max_cost)
        sum_cost = sum_cost+cost
        cost = 0.0
        n_cost = n_cost+1
      end if
    end select
    link => link%next
  end do
  cost3 = n_cost*max_cost
  write(stdout,'(a,1p,3e12.4)') ' costs ratios:     ', cost3/max(tiny(1e0),sum_cost)
  call io%header ('-evaluate: end', left=25)
  call trace%end()
END SUBROUTINE evaluate

!===============================================================================
!> Initialize a minimum amount of patch information in a Cartesian arrangement
!> This assumes that positions lie btw 0.0 and self%size
!===============================================================================
SUBROUTINE cartesian_init (self, key, patch)
  class (read_external_t):: self
  integer:: key(4)
  class (experiment_t), pointer:: patch
  !.............................................................................
  integer:: level, ix, iy, iz
  class (task_t), pointer:: task
  class(*), pointer:: ptr
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%cartesian_init')
  !-----------------------------------------------------------------------------
  ! Check if a patch already exists at the key
  !-----------------------------------------------------------------------------
  call self%patch_table%get_ptr (key, ptr)
  if (associated(ptr)) then
    select type (ptr)
    class is (experiment_t)
    patch => ptr
    end select
    if (verbose > 3) then
      write (stdlog,1) 'cartesian_init: patch already exists: id =', &
        patch%id, key
      1 format(a,i7,2x,4i5,3f10.4)
    end if
  else
    allocate (patch)
    ptr => patch
    call self%patch_table%set_ptr (key, ptr)
    level = key(1)
    ix = key(2)
    iy = key(3)
    iz = key(4)
    patch%key = key
    patch%level = level
    patch%n = 16
    self%size = 1d0
    patch%ds = self%size*0.5d0**level
    patch%size = patch%ds*patch%n
    patch%position = [(ix+0.5d0)*patch%size(1), &
                      (iy+0.5d0)*patch%size(2), &
                      (iz+0.5d0)*patch%size(3)]
    patch%ng = 3
    patch%nv = 8
    patch%nw = 1
    call patch%init
    patch%mhd = .true.
    patch%periodic = .true.
    call self%initialize (patch)  
    !---------------------------------------------------------------------------
    ! Make sure interpolations are done optimally
    !---------------------------------------------------------------------------
    patch%unsigned(patch%idx%d) = .true.
    patch%unsigned(patch%idx%e) = .true.
    patch%pervolume(patch%idx%px) = .true.
    patch%pervolume(patch%idx%py) = .true.
    patch%pervolume(patch%idx%pz) = .true.
    !---------------------------------------------------------------------------
    ! Patches below levelmin are needed for selfgravity
    !---------------------------------------------------------------------------
    if (patch%level < self%rootlvl) then
      call patch%set (bits%frozen)
    else
      call patch%set (bits%has_output)
    end if
    !---------------------------------------------------------------------------
    allocate (patch%filled(patch%gn(1),patch%gn(2),patch%gn(3)))
    call io%bits_mem (storage_size(patch%filled), &
                     product(shape(patch%filled)), 'fill')
    patch%filled = -1
    !---------------------------------------------------------------------------
    ! Add the patch to the task list -- this creates the patch%link
    !---------------------------------------------------------------------------
    call self%add_to_task_list (patch)
    ptr => patch
    if (verbose > 1) then
      write (stdlog,1) &
        'cartesian_init: new patch created: id, key =', &
        patch%id, key, patch%position
    end if
  end if
  call trace%end ()
END SUBROUTINE cartesian_init

!===============================================================================
!> Initialize reasonable default values, in case download has gaps
!===============================================================================
SUBROUTINE initialize (self, patch)
  class(read_external_t):: self
  class(patch_t):: patch
  integer:: iv
  !.............................................................................
  do iv=1,patch%nv
    if      (iv==patch%idx%d) then
      patch%mem(:,:,:,iv,:,:) = self%d_min
    else if (iv==patch%idx%e) then
      if (eos%gamma == 1d0) then
        patch%mem(:,:,:,iv,:,:) = self%d_min
      else
        patch%mem(:,:,:,iv,:,:) = self%d_min*self%tt_min/(eos%gamma-1d0)
      end if
    else
      patch%mem(:,:,:,iv,:,:) = 0.0
    end if
  end do
END SUBROUTINE initialize

!===============================================================================
!> Add nbor lists to all tasks
!===============================================================================
SUBROUTINE init_all_nbors (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%init_all_nbors', itimer=itimer)
  call io%header ('-init_all_nbors: begin', left=25)
  call timer%tic
  call self%task_list%reset_status                 ! update task bits
    link => self%task_list%head
    do while (associated(link))
      call self%init_nbors (link)
      link => link%next
    end do
  call timer%toc ('time to initialize nbor lists', 1)
  call io%header ('-init_all_nbors: end', left=25)
  call trace%end (itimer)
END SUBROUTINE init_all_nbors

!===============================================================================
!> Add nbor list to one task
!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(read_external_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer:: i1, i2, i3, n
  real(8):: size(3), position(3)
  type(hash_table_t):: nbor_table
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%init_nbors', itimer=itimer)
  call nbor_table%init (200, ndim=4)
  n = 0
  task => link%task
  select type (task)
  class is (experiment_t)
    !---------------------------------------------------------------------------
    ! For a level above levelmin, check all L-1 positions
    !---------------------------------------------------------------------------
    if (task%level > self%rootlvl) then
      size = task%size*2d0
      do i3=-1,1
      do i2=-1,1
      do i1=-1,1
        !-------------------------------------------------------------------
        ! Position of level L nbors, and key for a L-1 patches there
        !-------------------------------------------------------------------
        position = task%position + [i1,i2,i3]*size
        call check_nbor (task%level-1, position, size)
      end do
      end do
      end do
    end if
    !---------------------------------------------------------------------------
    ! At the same level, check all level L positions except self
    !---------------------------------------------------------------------------
    size = task%size
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
      if (i1==0 .and. i2==0 .and. i3==0) cycle
      position = task%position + [i1,i2,i3]*size
      call check_nbor (task%level, position, size)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! For a level below levelmax, check all L+1 positions
    !---------------------------------------------------------------------------
    if (task%level < shared%levelmax) then
      size = task%size/2d0
      do i3=0,3
      do i2=0,3
      do i1=0,3
        position = task%position + ([i1,i2,i3]-1.5d0)*size
        call check_nbor (task%level+1, position, size)
      end do
      end do
      end do
    end if
  end select
  link%nbors_by_level => link%nbor
  if (verbose > 1) &
    write (stdout,*) link%task%id, ' added', n, 'nbors', associated(link%nbor)
  call nbor_table%dealloc
  call trace%end (itimer)
contains
!===============================================================================
!> Check if nbor task exists, and if it overlaps -- if so add to nbor list
!===============================================================================
subroutine check_nbor (level, position, size)
  integer:: level
  real(8):: position(3), size(3)
  integer:: key(4)
  class(*), pointer:: nbtask, ptr
  integer, save, target:: c=1
  !-----------------------------------------------------------------------------
  ! Level and positions key, with periodic wrapping
  !-----------------------------------------------------------------------------
  key = task_key (level, task%box, position, size)
  call nbor_table%get_ptr (key, ptr)
  if (.not.associated(ptr)) then
    call self%patch_table%get_ptr (key, nbtask)
    if (associated(nbtask)) then
      ptr => c
      call nbor_table%set_ptr (key, ptr)
      select type(nbtask)
      class is (experiment_t)
        if (nbtask%overlaps (task)) then
          nbtask%quality = nbtask%level
          if (verbose > 2) &
            write (stdout,'(a,2(i5,":",i2.2),2x,4i3)') &
              'read_external_t%check_nbor: task,nbor:', &
              link%task%id, link%task%level, nbtask%id, nbtask%level, key
          n = n+1
          call link%add_nbor (nbtask%link)
        end if
      end select
    end if
  end if
end subroutine check_nbor
END SUBROUTINE init_nbors

!===============================================================================
!> Add all patches with level <= L+1 and within a good margin of overlap
!===============================================================================
SUBROUTINE init_nbors2 (list, link)
  type(task_list_t), pointer:: list
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: scan, nbor
  real(8):: dist(3)
  integer:: n
  !-----------------------------------------------------------------------------
  nullify (link%nbor)
  link%task%n_nbors = 0
  n = 0
  scan => list%head
  do while (associated(scan))
    !---------------------------------------------------------------------------
    ! Include only patches at lower levels, for prolongation to initial values
    !---------------------------------------------------------------------------
    if (abs(scan%task%level-link%task%level) < 2) then
      n = n+1
      dist = abs(scan%task%position-link%task%position)
      if (all(dist <= 0.5*(scan%task%size + link%task%size))) then
        call link%add_nbor (scan)
      end if
    end if
    scan => scan%next
  end do
  if (verbose > 1) &
    write(stdout,*) 'id, level, n, n_nbors =', &
      link%task%id, link%task%level, n, link%task%n_nbors
END SUBROUTINE init_nbors2

!===============================================================================
!> Remove the initial RAMSES octs
!===============================================================================
SUBROUTINE remove_ramses_octs (self)
  class (read_external_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: oct, next
  integer:: id
  integer, save:: itimer=0, nprint=100
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%remove_ramses_octs', itimer=itimer)
  call io%header ('-remove_ramses_octs: begin', left=25)
  oct => self%oct_list%head
  do while (associated(oct))
    next => oct%next
    select type (oct)
    class is (oct_t)
    if (allocated(oct%var)) then
      call io%bits_mem (-storage_size(oct%var), product(shape(oct%var)), 'oct')
      deallocate (oct%var)
    end if
    end select
    call io%bits_mem (-8, oct_size, 'oct')
    deallocate (oct)
    oct => next
  end do
  if (verbose > 0) then
    write (stdout,'(" elapsed:",f8.1,"s")') wallclock()-start
  end if
  call memtest_delay
  call io%header ('-remove_ramses_octs: end', left=25)
  call trace%end (itimer)
END SUBROUTINE remove_ramses_octs

!===============================================================================
!> Initial output of all tasks, in case no task updates follow
!===============================================================================
SUBROUTINE output (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call self%task_list%reset_status
  link => self%task_list%head
  do while (associated(link))
    call link%task%output
    link%task%out_next = 0d0
    link => link%next
  end do
END SUBROUTINE output

!===============================================================================
!> Convert the energy variable to total energy
!===============================================================================
SUBROUTINE total_energy (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  real, dimension(:,:,:), allocatable:: cx, cy, cz
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz
  integer:: m(3)
  real:: e_min
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%total_energy', itimer=itimer)
  call io%header ('-total_energy: begin', left=25)
  call timer%tic
  e_min = huge(1.0)
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      m = task%mesh%gn
      allocate (cx(m(1),m(2),m(3)))
      allocate (cy(m(1),m(2),m(3)))
      allocate (cz(m(1),m(2),m(3)))
      d  => task%mem(:,:,:,task%idx%d ,task%it,1)
      e  => task%mem(:,:,:,task%idx%e ,task%it,1)
      px => task%mem(:,:,:,task%idx%px,task%it,1)
      py => task%mem(:,:,:,task%idx%py,task%it,1)
      pz => task%mem(:,:,:,task%idx%pz,task%it,1)
      bx => task%mem(:,:,:,task%idx%bx,task%it,1)
      by => task%mem(:,:,:,task%idx%by,task%it,1)
      bz => task%mem(:,:,:,task%idx%bz,task%it,1)
      cx = stagger%xup(bx)
      cy = stagger%yup(by)
      cz = stagger%zup(bz)
      e_min = min(e_min, minval(e))
      e = e + 0.5*((px**2+py**2+pz**2)/d + (cx**2+cy**2+cz**2))
      deallocate (cx, cy, cz)
    end select
    link => link%next
  end do
  call timer%toc ('time to compute total energy', 1)
  call io%header ('-total_energy: end', left=25)
  call trace%end (itimer)
END SUBROUTINE total_energy

END MODULE read_external_mod
