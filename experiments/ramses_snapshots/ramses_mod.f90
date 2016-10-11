!> =============================================================================
!> Read a RAMSES snapshot and convert to DISPATCH data structures.
!>
!> In this first version, each 2**3 block in turned into 1 DISPATCH
!> patch.
!>
!> =============================================================================
MODULE ramses_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE mpi_mod
  USE task_mod
  USE patch_mod
  USE mesh_mod
  USE mhd_mod
  USE scaling_mod
  USE units_mod
  USE star_mod
  USE list_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  private

  integer, parameter :: dp=kind(1.0D0)                        ! double precision
  integer, parameter :: d1=kind(kind(1d0))                    ! kind of a kind

  integer, parameter :: nmetal=2

  type, public:: ramses_t
    logical:: big_endian, double_precision, quad_precision
    integer:: rank
    integer:: nout, minlvl, maxlvl
    integer:: npatch, nmerged_max
    real(kind=4):: scale_l, scale_d, scale_t
    real(kind=8):: gamma
    real(kind=8):: rank_load
    real(kind=8):: xmin, xmax, ymin, ymax, zmin, zmax         ! corners of the cutout
    real(kind=8):: size, d_min, tt_min
    logical :: centre_on_sink, merge_octs
    character(len=128):: path
    type(task_list_t):: task_list
  contains
    procedure :: init
    procedure :: read_snapshot
    procedure :: entropy
    procedure :: read_sink
    procedure :: merge_ramses_octs
  end type

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
    integer gid                         ! pointer to grid cell
    integer father                      ! pointer to father cell
    integer ilevel                      ! level of AMR hierachy where star lives
    integer old_ilevel                  ! level of hierachy star used to have, if changed
    integer owner                       ! owner thread
    integer has_moved
    real(dp) dtold
  end type
  type(sink_t):: sink

CONTAINS
!> =============================================================================
!> Facilitate call to `ramses%read`.  All heavy lifting is done there
!> =============================================================================
SUBROUTINE init (self, name)
  class(ramses_t):: self
  character(len=*), optional:: name
  logical, save:: first_time=.true.
  integer, save:: nout, minlvl, maxlvl, nmerged_max
  real, save:: xmin, xmax, ymin, ymax, zmin, zmax, rank_load
  real, save:: d_min, tt_min, size, size_au
  logical, save :: centre_on_sink, merge_octs
  character(len=64), save:: path='.'
  namelist / init_params / nout, xmin, xmax, ymin, ymax, zmin, zmax, size, &
                           nmerged_max, minlvl, maxlvl, rank_load, path, &
                           centre_on_sink, merge_octs
  !.............................................................................
  call trace_begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and save them
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time=.false.
    ! default values
    nout = -1
    xmin = 0.0; xmax = 1.0
    ymin = 0.0; ymax = 1.0
    zmin = 0.0; zmax = 1.0
    minlvl = 1
    maxlvl = 1
    rank_load = 1.0
    d_min = 1d-30 ! cgs
    tt_min = 5.   ! K
    path = './'
    size = -1d0
    centre_on_sink = .false.
    merge_octs = .true.
    nmerged_max = 16
    ! read namelist
    rewind (io%input); read (io%input,init_params)
    if (io%master) write(*, init_params)
  end if
  !$omp end critical (init_cr)
  !-----------------------------------------------------------------------------
  ! All threads save the parameters in their private parts, including the 
  ! scaling information.  The selection of RAMSES data we read is based on
  ! the parameters defined here.
  !-----------------------------------------------------------------------------
  call scaling%init()
  self%d_min     = d_min/scaling%d
  self%tt_min    = tt_min/scaling%temp
  self%xmin      = xmin
  self%xmax      = xmax
  self%ymin      = ymin
  self%ymax      = ymax
  self%zmin      = zmin
  self%zmax      = zmax
  self%minlvl    = minlvl
  self%maxlvl    = maxlvl
  self%path      = path
  self%npatch    = 1
  self%rank_load = rank_load
  self%nout      = nout
  self%size      = size
  self%centre_on_sink = centre_on_sink
  self%nmerged_max = nmerged_max
  self%merge_octs = merge_octs

  !-----------------------------------------------------------------------------
  ! Each rank reads the RAMSES grids belonging to it, filling in the data,
  ! and initializing the patch structure of every grid.
  !-----------------------------------------------------------------------------
  call self%read_snapshot()

  ! attempt to merge RAMSES 2**3 octs into larger patches
  if (self%merge_octs) call self%merge_ramses_octs(nmerged_max)

  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Read a RAMSES snapshot (likely split across multiple files), and create a
!> patch for each "grid".
!> =============================================================================
SUBROUTINE read_snapshot (self)
  class(ramses_t):: self

  !---- variables for reading the snapshot ----
  ! scalars
  integer(kind=4):: ilun=10, ilun2=11, ilun3=12, nlevelmin, icpu, isentinel, ncpu, ndims, &
                    nx, ny, nz, nlevelmax, nboundary, ngrids, noutput, iout, ifout, &
                    overload, ndomain, ilevel, i, j, k, l, ivar, ioerr, &
                    bit_length, maxdom, imin, imax, jmin, jmax, kmin, kmax, ngrid_act, &
                    nx_full, ny_full, nz_full, ncpu_read, ix, iy, iz, ncache, ndom, &
                    ncells, max_dom, nvar_mhd, twotondim, ngrids_max, ngrids_k, it
  real(kind=8):: boxlen, dx, dxnext, time, gamma_ramses, tmp_load, tot_load
  real(kind=16):: dmax, dkey, order_min, order_min_1d(1)
  real(kind=8):: xmn, xmx, ymn, ymx, zmn, zmx, r

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
  real(kind=8), dimension(3):: pmin=1d30, pmax=-1d30, xg_plus

  character(len=80):: ordering
  character(len=128):: filebase, filename, conv

  ! derived types
  type level
    integer(kind=4):: ngrid, imin, imax, jmin, jmax, kmin, kmax
  end type level
  type(level), allocatable, dimension(:):: levels

  ! diagnostic variables
  integer, allocatable, dimension(:,:):: ncellsum
  integer, save:: nprint_f=1000
  integer(8), save:: n_passive=0, n_active=0, n_onrank=0, n_offrank=0, &
    n_d_min=0, n_tt_min=0, n_frozen=0

  !---- dispatch variables ----
  class(task_t), pointer:: task, parent
  class(experiment_t), pointer:: patch

  !.............................................................................
  call trace_begin('read_ramses_snapshot')

  if (self%centre_on_sink .and. self%size > 0) then
    call self%read_sink
    xmn = star%position(1)-self%size*0.5d0
    xmx = star%position(1)+self%size*0.5d0
    ymn = star%position(2)-self%size*0.5d0
    ymx = star%position(2)+self%size*0.5d0
    zmn = star%position(3)-self%size*0.5d0
    zmx = star%position(3)+self%size*0.5d0
  end if
  self%rank = 0
  tmp_load = 0d0
  tot_load = 0d0

  write(filebase, '(a,"/output_",i5.5,"/amr_",i5.5,".out")') trim(self%path), self%nout, self%nout

  ! We need some basic info about the snapshot; let's read the `info_` and `amr_` icpu = 1 files.

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

    !------ read `amr_` file for icpu = 1 ------
    icpu = 1
    write(filename, '(a,i5.5)') trim(filebase), icpu

#if __INTEL_COMPILER || defined(__GFORTRAN__)
    open(unit=ilun,file=filename,status='old',form='unformatted',access='stream')
    read(ilun) isentinel
    self%big_endian = (isentinel .ne. 4)  ! the data is in big endian format.
    close(ilun)

    if (self%big_endian) then
      open(unit=ilun,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
    else
      open(unit=ilun,file=filename,status='old',form='unformatted')
    endif
#else
    open(unit=ilun,file=filename,status='old',form='unformatted')
#endif

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

    ! Try to read `r_boxlength` as a double precision variable; if an error occurs,
    ! then the data is single precision.
    read(ilun,iostat=ioerr) boxlen
    self%double_precision = (ioerr == 0)

    ! At this time, only double precision RAMSES data is supported.
    if (.not. self%double_precision) then
      call mpi%abort('Only double-precision RAMSES data is supported at this time.')
      close(ilun)
      stop
    endif

    read(ilun) noutput, iout, ifout
    print*,'noutput, iout, ifout:', noutput, iout, ifout
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

    ! levels variables
    allocate(head_list(1:ncpu,1:nlevelmax), tail_list(1:ncpu,1:nlevelmax), &
             numb_list(1:ncpu,1:nlevelmax), numb_tot(1:10,1:nlevelmax))
    read(ilun) head_list
    read(ilun) tail_list
    read(ilun) numb_list
    read(ilun) numb_tot

    ! boundary linked lists
    if (nboundary > 0) then
      allocate(head_bndy(1:nboundary,1:nlevelmax), tail_bndy(1:nboundary,1:nlevelmax), &
               numb_bndy(1:nboundary,1:nlevelmax))
      read(ilun) head_bndy
      read(ilun) tail_bndy
      read(ilun) numb_bndy
    endif

    read(ilun) ! SKIP headf, tailf, numbf, used_mem, used_mem_tot
    read(ilun) ordering
    if (ordering == 'bisection') then
      call mpi%abort('Only Peano-Hilbert ordered RAMSES data can be read at this time.')
      close(ilun)
      stop
    endif
    read(ilun) overload, ndomain
    backspace(ilun) ! rewind one entry (why?)

    ! ------ read (double or quad precision) Peano-Hilbert keys ------
    ! includes precision and endianness checking.
    if (ndomain /= ncpu * overload .or. ndomain == 0 .or. overload == 0) then
      ! "old" format
      overload = 1
      ndomain = ncpu
      print*,'ndomain,ncpu:',ndomain,ncpu
      allocate(bound_key(0:ndomain), bound_key_dp(0:ndomain), isub_domain(1:ndomain), &
             cpu_domain(1:ndomain), map_domain(1:overload,1:ncpu))
      read(ilun, iostat=ioerr) bound_key
      if (ioerr /= 0) then
        ! the keys were stored at double precision
        backspace(ilun)
        read(ilun,iostat=ioerr) bound_key_dp
        bound_key = bound_key_dp ! put the key into a quad precision variable.
        self%quad_precision = .false.
      else
        backspace(ilun)
        self%quad_precision = .true.
        allocate(bound_key_dpx2(2,0:ndomain))
        read(ilun,iostat=ioerr) bound_key_dpx2
        if (any( (abs( real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16) - bound_key) &
               / (abs(bound_key) + 1e-45_16)) > 1e-15_16) ) then
#ifdef __xlc__
          ! We are on BlueGene/Q and using the XLF90 compiler. These two numbers should be identical.
          ! Therefore we are trying to read little endian data, but will never have any luck doing that,
          ! since the little-endian quad-precision format is not supported with XLF.
          call mpi%abort('We are on an IBM big-endian platform trying to read little endian quad-precision data. This is not possible. Giving up.')
          close(ilun)
          stop
#endif
#if __INTEL_COMPILER || defined(__GFORTRAN__)
          ! We should test if we are byteswapping or not. If we are reading big-endian, we should use
          ! the double-double format (because then the data comes from JUQUEEN), else we should use
          ! quad-precision data, because it comes from a x86 run.
          inquire(unit=ilun,convert=conv)
          if (conv=='BIG_ENDIAN') then
            bound_key = real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16)
            print *, 'Data is from a BIG-ENDIAN platform. Assuming from Blue-Gene.'
          endif
#endif
        endif
        deallocate(bound_key_dpx2)
      endif
      isub_domain = 1
      do i=1,ndomain
        cpu_domain(i) = i
        map_domain(1,i)= i
      enddo
      print *,'Using old format for key storage.'

    else
      ! "new" format
      print*,'ndomain,ncpu:',ndomain,ncpu
      allocate(bound_key(0:ndomain), bound_key_dp(0:ndomain), isub_domain(1:ndomain), &
               cpu_domain(1:ndomain), map_domain(1:overload,1:ncpu))
      read(ilun,iostat=ioerr) overload, ndomain, bound_key, isub_domain, cpu_domain, map_domain
      if (ioerr /= 0) then
        ! the keys were stored at double precision
        backspace(ilun)
        read(ilun,iostat=ioerr) overload, ndomain, bound_key_dp, isub_domain, cpu_domain, map_domain
        bound_key = bound_key_dp ! put the key into a quad precision variable.
        self%quad_precision = .false.
      else
        backspace(ilun)
        self%quad_precision = .true.
        allocate(bound_key_dpx2(2,0:ndomain))
        read(ilun,iostat=ioerr) overload, ndomain, bound_key_dpx2, isub_domain, cpu_domain, map_domain
        if (any( (abs( real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16) - bound_key) &
               / (abs(bound_key) + 1e-45_16)) > 1e-15_16) ) then
#ifdef __xlc__
          ! We are on BlueGene/Q and using the XLF90 compiler. These two numbers should be identical.
          ! Therefore we are trying to read little endian data, but will never have any luck doing that,
          ! since the little-endian quad-precision format is not supported with XLF.
          call mpi%abort('We are on an IBM big-endian platform trying to read little endian quad-precision data. This is not possible. Giving up.')
          close(ilun)
          stop
#endif
#if __INTEL_COMPILER || defined(__GFORTRAN__)
          ! We should test if we are byteswapping or not. If we are reading big-endian, we should use
          ! the double-double format (because then the data comes from JUQUEEN), else we should use
          ! quad-precision data, because it comes from a x86 run.
          inquire(unit=ilun,convert=conv)
          if (conv=='BIG_ENDIAN') then
            bound_key = real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16)
            print *, 'Data is from a BIG-ENDIAN platform. Assuming from Blue-Gene.'
          endif
#endif
        endif
        deallocate(bound_key_dpx2)
      endif

      print *,'Using new format for key storage.'
    endif
    deallocate(bound_key_dp)
    ! ------ finished reading Peano-Hilbert keys ------

    ! we have extracted enough information; close the file (for the time being).
    close(ilun)

    ! output info about snapshot to terminal
    print *, 'Data is '//merge('double','single',self%double_precision)//' precision with '// &
             trim(merge('quad  ','double',self%quad_precision))//' precision hilbert keys'
    print *, 'RAMSES snapshot time = ',time
    print *, 'ncpu = ',ncpu
    print *, 'nlevelmin = ',nlevelmin
    print *, 'nlevelmax = ',nlevelmax
    print *, 'boxlen = ',boxlen
    print *, 'time = ',time

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

    !------ given limits, make preparations to cut a sub-volume from the box ------
    ! uses the Hilbert-Peano curve to select which CPUs contain the relevant data
    dmax = max(self%xmax - self%xmin, self%ymax - self%ymin, self%zmax - self%zmin)
    ! for a given cut, what is the bare minimum level we need to consider.
    do ilevel=1,self%maxlvl
      dx = 0.5_8**ilevel
      if (dx < dmax) exit
    enddo
    bit_length = ilevel - 1
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

    ! ------ loop over CPU files ------
    call mpi%barrier ('ncpu begin')
    do k=1,ncpu_read
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
      if (io%master) print *, 'Reading file '//trim(filename)

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
      if (self%rank==mpi%rank) &
        print'(a,i4.1,a,/,40i6)','ngrids(cpu=',k,')',ngrid_level(k,1:nlevelmax)
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
      self%gamma = gamma_ramses

      ! open the `grav_` file to get access to the hydro data
      write(filename, '(a,"/output_",i5.5,"/grav_",i5.5,".out",i5.5)') trim(self%path), &
            self%nout, self%nout, icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (self%big_endian) then
        open(unit=ilun3,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
      else
        open(unit=ilun3,file=filename,status='old',form='unformatted')
      endif
#else
      open(unit=ilun3,file=filename,status='old',form='unformatted')
#endif
      read(ilun3) ! SKIP ncpu_dump
      read(ilun3) ! SKIP ndims
      read(ilun3) ! SKIP nlevelmax
      read(ilun3) ! SKIP nboundary

      ! loop over levels for current file
      do ilevel=1,self%maxlvl
        ! set cell geometry info
        dx = 0.5_8**ilevel
        dxnext = 0.5_8 * dx
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

        ! allocate worker/scratch arrays
        ngrid_act = ngrid_file(icpu,ilevel)
        levels(ilevel)%ngrid = ngrid_act
        if (ngrid_act > 0) allocate(xg(1:ngrid_act,1:ndims), son(1:ngrid_act,1:twotondim), &
                                    x(1:ndims,1:ngrid_act), var(1:ngrid_act,1:twotondim,1:nvar_mhd+1), ref(1:ngrid_act))

        ! loop over domains
        do j=1,nboundary+ncpu
          ncache = ngrid_file(j,ilevel)

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

          ! read GRAV data
          read(ilun3) ! SKIP ilevel
          read(ilun3) ! SKIP ncache
          if (ncache > 0) then
            do i=1,twotondim
              do ivar=1,4
                if (j == icpu .and. i==4) then
                  read(ilun3) var(:,i,nvar_mhd+1)
                else
                  read(ilun3) ! SKIP
                endif
              enddo
            enddo
          endif

        enddo ! end loop j=1,nboundary+ncpu

        ! compute cut-out
        ncellsum(ilevel,k) = 0
        if (ngrid_act > 0 .and. ilevel >= max(self%minlvl,nlevelmin)) then
          do l=1,ngrid_act
            ! the grid centre is given by the `xg` variable.
            ! check if cell is in the cut-out region.
            if (self%size > 0d0 .and. self%centre_on_sink) then
              r = sqrt(sum((xg(l,:)-star%position)**2))
              ref(l) = r < self%size
            else
              ref(l) = (xg(l,1) > self%xmin .and. xg(l,1) < self%xmax .and. &
                        xg(l,2) > self%ymin .and. xg(l,2) < self%ymax .and. &
                        xg(l,3) > self%zmin .and. xg(l,3) < self%zmax)
            end if
            pmin = min(pmin,xg(l,:))
            pmax = max(pmax,xg(l,:))
            if (ref(l)) then
              ! create a new DISPATCH patch/task for the current grid and fill with the
              ! cells with data if we own the patch
              self%npatch = self%npatch+1
              allocate (experiment_t :: task)
              select type (task)
              class is (experiment_t)
                patch => task
              end select
              patch%nv = 8
              patch%rank = self%rank
              patch%level = ilevel
              patch%position(:) = xg(l,:)
              patch%n = 2
              patch%size = 2.0 * dx
              call patch%init()
              ! if the patch is at the very boundary, freeze it!
              xg_plus = xg(l,:) + 2d0 * dx
              if (self%size > 0d0 .and. self%centre_on_sink) then
                if (r + 2d0 * dx > self%size) then
                  patch%status = ior(patch%status,bits%frozen)
                  n_frozen = n_frozen + 1
                  if (nprint_f > 0) then
                    nprint_f = nprint_f - 1
                    print '(a,i7.1,a,1e16.9,1x,i3.1)', 'patch', patch%id, ' frozen', r, ilevel
                  end if
                end if
              else if (xg_plus(1) < self%xmin .or. xg_plus(1) > self%xmax .or. &
                       xg_plus(2) < self%ymin .or. xg_plus(2) > self%ymax .or. &
                       xg_plus(3) < self%zmin .or. xg_plus(3) > self%zmax) then
                patch%status = ior(patch%status,bits%frozen)
                n_frozen = n_frozen + 1
                if (nprint_f > 0) then
                  nprint_f = nprint_f - 1
                  print '(a,i7.1,a,1p,3e16.9,i3.1)', 'patch', patch%id, ' frozen', xg(l,:), ilevel
                end if
              end if
              patch%time = time
              patch%out_next = 0.999 * time ! output to file immediately
              patch%velocity = 0.0
              call patch%set(bits%static)
              tmp_load = tmp_load + 1d-6*product(patch%n)*0.5**(nlevelmax-ilevel)
              tot_load = tot_load + 1d-6*product(patch%n)*0.5**(nlevelmax-ilevel)
              if (patch%rank==mpi%rank) then
                call patch%init()
                if (mod(self%npatch,1000)==0) &
                  print'(a,2i7,3x,a,1p,2e10.2)','npatch:',self%npatch,n_active,'load:',tot_load,tmp_load
                !!!patch%offset = patch%position ! why did I put this here? it causes issues later...
                ! loop over cells
                j = 1
                ! RAMSES variable mapping:
                ! 1 = density; 2:4 = zone-centred velocity components;
                ! 5:7 = "left" magnetic field components; 8:10 = "right" magnetic field components;
                ! 11 = thermal pressure; 12:13 = passive scalars.
                ! DISPATCH variable mapping:
                ! 1 = density; 2 = entropy; 3 = grav. potential; 4:6 = momentum components
                ! 7:9 = magnetic field components
                it = patch%it
                ! the RAMSES snapshot does not contain any information to fill boundary
                ! cells; this needs to be done by DISPATCH. To prevent the boundary cells
                ! from containing garbage, initialise them.
                ! Note: should `patch%init()` do this by default?
                patch%mem(:,:,:,:,:,:) = 0.0

                do iz=patch%mesh(3)%li,patch%mesh(3)%ui
                  do iy=patch%mesh(2)%li,patch%mesh(2)%ui
                    do ix=patch%mesh(1)%li,patch%mesh(1)%ui
                      ! note that ghost cells are left *empty* at this stage.
                      if (var(l,j,1) < self%d_min) then
                        var(l,j,1) = self%d_min
                        n_d_min = n_d_min+1
                      end if
                      if (var(l,j,11) < self%tt_min*var(l,j,1)) then
                        var(l,j,11) = self%tt_min*var(l,j,1)
                        n_tt_min = n_tt_min+1
                      end if
                      patch%mem(ix,iy,iz,1,patch%it,1) = var(l,j,1) ! density
                      patch%mem(ix,iy,iz,2,patch%it,1) = self%entropy(var(l,j,:)) ! entropy
                      patch%mem(ix,iy,iz,3,patch%it,1) = var(l,j,1) * var(l,j,2) ! convert to momentum
                      patch%mem(ix,iy,iz,4,patch%it,1) = var(l,j,1) * var(l,j,3)
                      patch%mem(ix,iy,iz,5,patch%it,1) = var(l,j,1) * var(l,j,4)
                      patch%mem(ix,iy,iz,6,patch%it,1) = var(l,j,5) ! only the "left" fields are kept
                      patch%mem(ix,iy,iz,7,patch%it,1) = var(l,j,6)
                      patch%mem(ix,iy,iz,8,patch%it,1) = var(l,j,7)
                      !!!patch%mem(ix,iy,iz,3,patch%it,1) = var(l,j,nvar_mhd+1) ! grav. pot; SKIP!
                      j = j + 1
                    enddo
                  enddo
                enddo
                ! note that neighbour lists still need to be determined (see below).
                call self%task_list%append (task)

                ! diagnostic info
                ncellsum(ilevel,k) = ncellsum(ilevel,k) + twotondim
                n_onrank = n_onrank + 1
              else
                n_offrank = n_offrank + 1
              endif
              n_active = n_active + 1
            else
              n_passive = n_passive + 1
            endif
          enddo
        endif ! end ngrid_act > 0 .and. ilevel >= nlevelmin
        if (ngrid_act > 0) deallocate(son, xg, var, ref, x)

      enddo ! end loop ilevel=1,self%maxlvl

      close(ilun) ! close `amr_` file
      close(ilun2) ! close `hydro_` file
      if (self%rank==mpi%rank) print'(4(a,i8,3x),a,1p,e9.2)', &
        'n_active:',n_active,'n_passive:',n_passive,'n_onrank:',n_onrank, &
        'n_offrank:',n_offrank,'tmp_load:',tmp_load
      !
      if (tmp_load > self%rank_load) then
        if (self%rank<mpi%size-1) then
          if (self%rank==mpi%rank) &
            print*,'bumping rank, with tmp_load =', tmp_load, n_active
          self%rank = self%rank + 1
          n_active=0
          n_passive=0
          n_onrank=0
          n_offrank=0
          tmp_load = 0d0
        end if
      end if
      call mpi%barrier ('ncpu loop')
    enddo ! end loop k=1,ncpu_read
    write (*,'(a)') '=================================================='
    write (*,'(a,2(i6.1,1x))') 'n_d_min, n_tt_min =',n_d_min,n_tt_min
    write (*,'(a,i11.1)') 'n_frozen = ', n_frozen
    write (*,'(a,i11.1)') 'n_patch = ',self%npatch

  !$omp end master
  !$omp barrier

  ! set scale factors (other scale factors will be calculated)
  scaling%l = self%scale_l
  scaling%t = self%scale_t
  scaling%d = self%scale_d

  ! output some diagnostics to the terminal.
  allocate (fill(self%maxlvl))
  do ilevel=nlevelmin,self%maxlvl
    fill(ilevel)=sum(ncellsum(ilevel,:))/2.0**(3*(ilevel-1))
  enddo
  write (*,'(a)') '=================================================='
  if (self%size > 0 .and. self%centre_on_sink) then
    write (*,'(a,1p,3(1e14.7,2x))') 'pmin (AU) =', &
      (pmin-star%position)*scaling%l/cgs%au
    write (*,'(a,1p,3(1e14.7,2x))')'pmax (AU) =', &
      (pmax-star%position)*scaling%l/cgs%au
  end if
  write (*,'("fill (by level):",30f7.3)') fill(nlevelmin:self%maxlvl)
  deallocate (fill)

  deallocate(ngrid_file, ngrid_level, cpu_list, cpu_read, xc, father, &
             neighbour, levels, ncellsum)

  ! finish task list initialisation
  print*,'Finishing task list initialisation (including neighbour determination).'
  print*,'This could take a while...'
  call self%task_list%init ('ramses_snapshot')
  call self%task_list%init_all_nbors

  call trace_end
END SUBROUTINE read_snapshot

!> =============================================================================
!> Lifted directly from `amr2cell.f90`/`hilbert3d_qdp` (in the Copenhagen
!> RAMSES repository).
!> =============================================================================
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
     print *,'Maximum bit length = ', bit_size(bit_length)
     print *,'ERROR in hilbert3d'
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

!> =============================================================================
!> Convenience function to calculate the entropy (because the STAGGER solver in
!> DISPATCH uses the entropy rather than the internal energy/pressure or total
!> energy.)
!> =============================================================================
FUNCTION entropy (self, vars) RESULT (s)
  class(ramses_t):: self
  real(kind=8):: s
  real(kind=8), dimension(:):: vars
  !.............................................................................
  if (self%gamma==1d0) then
    s = vars(11)/vars(1)                ! temperature
  else
    s = vars(1)*(log(vars(11))-log(vars(1))*self%gamma)/(self%gamma-1.0)
  end if
END FUNCTION entropy

!> =============================================================================
!> Read sink particle information; the zoom-in sink is the last one, generally
!> =============================================================================
 SUBROUTINE read_sink (self)
  class(ramses_t):: self
  integer, parameter:: sink_unit=9
  character(len=80):: fname
  logical:: do_translate
  integer:: center_star, version, nstar, szstar, sizeofstar, i, nstars, nmet, dp1
  real(4):: t,beta_create,dtold
  real(8):: t8,tflush,r_translate(3),v_translate(3),x_refine,y_refine,z_refine
  !.............................................................................
  call trace_begin('ramses_t%read_sink')

  print*,'=================================================='
  write (fname,'(a,a,i5.5,a)') trim(self%path), '/output_', self%nout, '/stars_output.dat'
  print*,fname
  open (sink_unit, file=fname,form='unformatted', status='unknown')
  read (sink_unit) version,nstars,t,szstar,sizeofstar
  print*,'versions =',version
  read (sink_unit) dp1,nmet,t8
  print*,'dp, nmetal =',dp1,nmet
  print*,'t8 =',t8
  print*,'nstars =',nstars
  read (sink_unit) do_translate, center_star, r_translate, v_translate,  &
                   x_refine, y_refine, z_refine
  do i=1,nstars
    read (sink_unit) sink%x,sink%y,sink%z,sink%px,sink%py,&
                     sink%pz,sink%dpx,sink%dpy,sink%dpz,sink%mass,&
                     sink%dmass,sink%potential,sink%density,sink%t_create,&
                     sink%t_explode,sink%metal(1:nmetal),sink%dmetal(1:nmetal),&
                     sink%gid,sink%father,sink%ilevel,sink%old_ilevel, &
                     sink%owner,sink%has_moved,sink%beta_create,sink%dtold
  enddo
  close(sink_unit)
  star%position = [sink%x,sink%y,sink%z]
  star%mass = sink%mass
  write (*,*) 'type(sink_t) :: star => last sink in file'
  write (*,'(a,3(1e11.5,1x),a,1i3.1)') ' pos = ',star%position,'; level = ',sink%ilevel
  print*,'=================================================='
  
  call trace_end
END SUBROUTINE read_sink

!> =============================================================================
!> Attempt to merge 2**3 RAMSES octs into large patches (up to nmax**3 cells).
!>
!> - This subroutine is specific to RAMSES data; it does not cover the possibility
!>   of overlapping patches.
!> - Only existing octs are used; no interpolation is performed.
!> - It is assumed that the neighbour list only contains octs which abut the
!>   current patch.
!> - Even with `merge_nbor_lists`, the neighbour list will still need to be
!>   reconstructed at the end of this procedure to ensure the new, merged patches
!>   are properly accounted for.
!> =============================================================================
SUBROUTINE merge_ramses_octs (self, n_max)
  class(ramses_t) :: self
  integer, intent(in) :: n_max
  class(task_t), pointer :: task, nbor_task, merged_task
  class(patch_t), pointer :: merged_patch
  class(link_t), pointer :: link, nbor, old_link, new_link
  integer :: n_now(3), n_nbor(3), n_updated(3), iterations, adir, n_merges=1, &
             merge_count, patch_count, diff_count
  real :: dpos(3)
  logical :: merged
  !.............................................................................
  call trace_begin('ramses_mod::merge_ramses_octs')

  patch_count = 0
  merge_count = 0
  iterations = 0
  
  do while (n_merges > 0)
    n_merges = 0
    link => self%task_list%head
    do while (associated(link))
      merged = .false.
      task => link%task

      ! all tasks should have their ready status reset to ensure the list and
      ! queue can be properly reconstructed at the end of the merger process.
      if (task%is_set(bits%ready)) then
        call task%clear(bits%ready)
        self%task_list%n_ready = self%task_list%n_ready - 1
      end if

      if (task%is_set(bits%remove)) then
        if (io%verbose > 1) print*,'Skipping link ID = ',task%id
        link => link%next
        cycle
      end if
    
      select type (task)
      class is (patch_t)
        n_now = task%n
      end select
      patch_count = patch_count + 1
    
      ! only proceed if no. of patch cells is less than `n_max`.
      if (any(n_now >= n_max)) then
        link => link%next
        cycle
      end if
    
      nbor => link%nbor
      do while (associated(nbor))
        nbor_task => nbor%task
        if (nbor_task%is_set(bits%remove)) then
          if (io%verbose > 1) print*,'Skipping nbor ID = ',nbor_task%id
          nbor => nbor%next
          cycle
        end if
    
        ! skip if neighbour is not at the same level
        if (task%level /= nbor%task%level) then
          nbor => nbor%next
          cycle
        end if
    
        ! only allow merges between octs that share a face.
        dpos = task%position - nbor%task%position
        if (count(abs(dpos) > 0.0) > 1) then
          nbor => nbor%next
          cycle
        end if
        adir = maxloc(abs(dpos),dim=1)
    
        select type (nbor_task)
        class is (patch_t)
          n_nbor = nbor_task%n
        end select

        diff_count = count(abs(n_nbor - n_now) > 0)
        if (n_nbor(adir) - n_now(adir) > 0) diff_count = diff_count - 1
        ! don't permit merges of patches with different sizes (except along
        ! the direction of the proposed merge, that's okay).
        if (diff_count > 0) then
          nbor => nbor%next
          cycle
        end if

        if (all(n_nbor < n_max)) then
          ! try to merge two patches; includes inheritence of existing neighbours
          merged_patch => merge_neighbours(link, nbor, dpos)
          ! add new merged patch to the *end* of the list
          merged_task => merged_patch
          call self%task_list%append(merged_task)
          ! merge the neighbour lists of the merged patches (and give it the
          ! result to the new task).
          associate (back_link=>merged_patch%link)
          select type (back_link)
          type is (link_t)
            new_link => back_link
          end select
          !!!call merge_nbor_lists(link, nbor, new_link)
          end associate
    
          ! prepare to remove neighbour and current patch from task list
          merged = .true.
          n_merges = n_merges + 1
          merge_count = merge_count + 1
          n_updated = merged_patch%n
          write (*,'(a,i10.1,a,i10.1,a,i10.1,a,3(i3.1,1x),a)') 'patches merged: ids = ', &
            link%task%id, ', ', nbor%task%id,' => ', merged_task%id, ' (', n_updated, ')'
    
          ! let's be conservative: only allow one merge at a given time; this will
          ! be slower, but safer.
          call nbor%task%set(bits%remove)  ! flag task for removal.
          exit ! exit do while over nbor
        end if
        nbor => nbor%next
      end do
      old_link => link
      link => link%next
      ! now it's safe to remove current link
      if (merged) then
        call old_link%task%set(bits%remove)  ! flag task for removal.
      end if
    end do
    
    iterations = iterations + 1
    write (*,'(a,i3.1)')  'iteration = ', iterations
    write (*,'(a,i10.1)') 'no. of merges = ', n_merges
    write (*,'(a,i10.1)') 'total no. of merges = ', merge_count
    write (*,'(a,i10.1)') 'no. of patches examined = ', patch_count
    
    ! Last but not least, tidy up. To ensure the merged patches have the
    ! correct neighbours, we need to search for neighbours again.
    call tidy_up(self%task_list)
    write (*,'(a,i10.1)') 'current no. of patches = ', self%task_list%n
    call self%task_list%init_all_nbors()
  end do

  ! Finished. Check the ready state of all patches.
  call self%task_list%check_all()

  call trace_end
END SUBROUTINE merge_ramses_octs

!> =============================================================================
!> Perform the actual merger of two patches; does not carry over the neighbour
!> lists of the two patches.
!> =============================================================================
FUNCTION merge_neighbours(link, nbor, dpos) RESULT(merger)
  class(link_t), pointer :: link, nbor
  class(patch_t), pointer :: patch, nbor_patch, merger
  class(mesh_t), pointer :: pm, nm, mm
  real, intent(in) :: dpos(3)
  integer :: it, ix, iy, iz, iv, idir, adir, dum(1), o(3)
  !.............................................................................
  patch => task2patch(link%task)
  nbor_patch => task2patch(nbor%task)

  adir = maxloc(abs(dpos),dim=1)
  idir = int(sign(float(adir), dpos(adir)))

  ! create new patch
  allocate (merger, source=patch)
  merger%mem_allocated = .false.
  nullify (merger%mem)
  ! the merged patch only grows in the `idir`-direction
  merger%size = patch%size
  merger%size(adir) = patch%size(adir) + nbor_patch%size(adir)
  merger%position = (patch%position + nbor_patch%position) / 2
  merger%n = patch%n
  merger%n(adir) = patch%n(adir) + nbor_patch%n(adir)
  merger%id = 0 ! will result in the assignment of an ID
  call merger%clear(bits%ready) ! clear ready bit so that checks are made by `task_list%check_all`.
  call merger%init()

  ! fill variable array
  it = merger%it
  pm => patch%mesh(adir)
  nm => nbor_patch%mesh(adir)
  mm => merger%mesh(adir)
  ! first copy over `patch`
  ! this offset formula reasonably assumes edges of RAMSES octs are co-spatial.
  o = 0
  o(adir) = int(((pm%p + pm%r(pm%li)) - (mm%p + mm%r(mm%li)))/mm%d)
  merger%mem(:,:,:,:,:,:) = 0.0
  do iz=patch%mesh(3)%li,patch%mesh(3)%ui
    do iy=patch%mesh(2)%li,patch%mesh(2)%ui
      do ix=patch%mesh(1)%li,patch%mesh(1)%ui
        do iv=1,patch%nv
          merger%mem(ix+o(1),iy+o(2),iz+o(3),iv,it,1) = patch%mem(ix,iy,iz,iv,it,1)
        end do
      end do
    end do
  end do

  ! then copy `nbor_patch`
  o = 0
  o(adir) = int(((nm%p + nm%r(nm%li)) - (mm%p + mm%r(mm%li)))/mm%d)
  do iz=nbor_patch%mesh(3)%li,nbor_patch%mesh(3)%ui
    do iy=nbor_patch%mesh(2)%li,nbor_patch%mesh(2)%ui
      do ix=nbor_patch%mesh(1)%li,nbor_patch%mesh(1)%ui
        do iv=1,patch%nv
          merger%mem(ix+o(1),iy+o(2),iz+o(3),iv,it,1) = nbor_patch%mem(ix,iy,iz,iv,it,1)
        end do
      end do
    end do
  end do
  ! done.
  nullify (pm, nm, mm)

END FUNCTION merge_neighbours

!> =============================================================================
!> Walk through the task list, deallocating tasks which are flagged for removal
!> (using `bits%remove`) while also removing the associated link in the task list.
!> =============================================================================
SUBROUTINE tidy_up (task_list)
  class(task_list_t) :: task_list
  class(task_t), pointer :: task
  class(link_t), pointer :: link, prev_link
  integer :: n_removed
  !.............................................................................
  call trace_begin('ramses_mod::tidy_up')

  ! deallocate tasks (and associated links) that are marked for removal.
  n_removed = 0
  link => task_list%head
  do while (associated(link))
    task => link%task
    prev_link => link
    link => link%next
    if (task%is_set(bits%remove)) then
      if (io%verbose > 1) print*,'removing task id ',task%id
      call task_list%remove_link(prev_link)
      deallocate(task)
      n_removed = n_removed + 1
    end if
  end do
  ! we want to remake the queue from scratch.
  nullify(task_list%queue)
  task_list%nq = 0

  write (*,'(a,i8.1)') 'no. of tasks removed = ', n_removed

  call trace_end
END SUBROUTINE tidy_up

END MODULE ramses_mod
