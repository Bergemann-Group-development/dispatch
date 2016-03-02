!> =============================================================================
!> Read a RAMSES snapshot and convert to DISPATCH data structures.
!>
!> In this first version, each 2**3 block in turned into 1 DISPATCH
!> patch.
!>
!> =============================================================================
MODULE read_ramses_mod
  USE mpi_base,    only: master
  USE plist_mod,   only: patch_list, wait_queue, x_list, add_patch_by_x, &
                         add_patch_by_quality, append_to_plist, check_quality_order
  USE patch_mod,   only: patch_t, init_patch_task, npatch, patch_table
  USE task_mod,    only: task_t, t_now
  USE state_mod,   only: gamma
  USE scalar_mod,  only: scalar_stats
  USE scaling_mod, only: code, cgs
  USE star_mod,    only: star
  USE mpi_base,    only: mpi_size
  USE io_unig,     only: unig
  USE io
  implicit none
  logical:: big_endian, double_precision, quad_precision
  integer:: rank
  real(kind=4), save:: scale_l, scale_d, scale_t
  real(kind=8), save:: rank_load
  integer, save:: nout, minlvl, maxlvl
  real(kind=8), save:: size_au, size, d_min, tt_min
  character(len=128):: path
PRIVATE
PUBLIC read_ramses
CONTAINS

SUBROUTINE read_ramses
  namelist / init_params / nout, size_au, &
                           minlvl, maxlvl, rank_load, path, d_min, tt_min
  logical:: error
  type(patch_t), pointer:: patch
  !.............................................................................
  ! default values
  nout = -1
  minlvl = 1
  maxlvl = 1
  d_min = 1d-30 ! cgs
  tt_min = 5.   ! K
  rank_load = 1.0
  path = './'
  size = 1d0
  size_au = size*code%l/cgs%au
  ! read namelist
  rewind (input); read (input,init_params); if (master) write(*, init_params)
  size = size_au*cgs%au/code%l
  npatch = 1
  d_min = d_min/code%d
  tt_min = tt_min/code%temp
  call read_ramses_snapshot (path, nout, size, &
                             minlvl, maxlvl, d_min, tt_min)
END SUBROUTINE read_ramses

!> =============================================================================
!> Read a RAMSES snapshot (likely split across multiple files), and create a
!> patch for each "grid".
!> =============================================================================
SUBROUTINE read_ramses_snapshot (datadir, nout, size, minlvl, maxlvl, &
                                 d_min, tt_min)

  !---- variables for reading the snapshot ----
  ! scalars
  integer(kind=4):: ilun=10, ilun2=11, ilun3=12, nlevelmin, icpu, isentinel, ncpu, ndims, &
                    nx, ny, nz, nlevelmax, nboundary, ngrids, noutput, iout, ifout, &
                    overload, ndomain, minlvl, maxlvl, ilevel, i, j, k, l, ivar, ioerr, &
                    bit_length, maxdom, imin, imax, jmin, jmax, kmin, kmax, ngrid_act, &
                    nx_full, ny_full, nz_full, ncpu_read, ix, iy, iz, nout, ncache, ndom, &
                    ncells, max_dom, nvar_mhd, twotondim, ngrids_max, ngrids_k, it
  real(kind=8):: boxlen, dx, dxnext, time, gamma_ramses, tmp_load, tot_load
  real(kind=16):: dmax, dkey, order_min, order_min_1d(1)
  real(kind=8):: xmin, xmax, ymin, ymax, zmin, zmax, d_min, tt_min
  real(kind=8):: size, xmn, xmx, ymn, ymx, zmn, zmx, r

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

  character(len=80):: ordering
  character(len=128):: datadir, filebase, filename, conv
  integer(8), save:: n_passive=0, n_active=0, n_onrank=0, n_offrank=0, &
    n_d_min=0, n_tt_min=0

  ! derived types
  type level
    integer(kind=4):: ngrid, imin, imax, jmin, jmax, kmin, kmax
  end type level
  type(level), allocatable, dimension(:):: levels

  ! diagnostic variables
  integer, allocatable, dimension(:,:):: ncellsum
  integer, save:: nprint=100, nprint_f=100

  !---- dispatch variables ----
  type(patch_t), pointer:: patch
  type(task_t), pointer:: task, parent

  !.............................................................................
  call trace_begin('read_ramses_snapshot')
  call read_sink (datadir, nout)
  xmin = 0.0; xmax = 1.0
  ymin = 0.0; ymax = 1.0
  zmin = 0.0; zmax = 1.0
  xmn = star%position(1)-size*0.5d0
  xmx = star%position(1)+size*0.5d0
  ymn = star%position(2)-size*0.5d0
  ymx = star%position(2)+size*0.5d0
  zmn = star%position(3)-size*0.5d0
  zmx = star%position(3)+size*0.5d0
  rank = 0
  tmp_load = 0d0
  tot_load = 0d0

  write(filebase, '(a,"/output_",i5.5,"/amr_",i5.5,".out")') trim(datadir), nout, nout

  ! We need some basic info about the snapshot; let's read the `info_` and `amr_` icpu = 1 files.

  !$omp master
  !if (master) then
    !---- read snapshot `info_` file ----
    ! note that the `info_` and `amr_` files contain some redundant information.
    write(filename, '(a,"/output_",i5.5,"/info_",i5.5,".txt")') trim(datadir), nout, nout
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
    read(ilun,'(13x,e23.15)') scale_l ! length scale factor
    read(ilun,'(13x,e23.15)') scale_d ! density scale factor
    read(ilun,'(13x,e23.15)') scale_t ! time scale factor
    close(ilun)

    !---- read `amr_` file for icpu = 1 ----
    icpu = 1
    write(filename, '(a,i5.5)') trim(filebase), icpu

#if __INTEL_COMPILER || defined(__GFORTRAN__)
    open(unit=ilun,file=filename,status='old',form='unformatted',access='stream')
    read(ilun) isentinel
    big_endian = (isentinel .ne. 4)  ! the data is in big endian format.
    close(ilun)

    if (big_endian) then
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
      call die('Only 3-D data is supported at this time.')
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
    double_precision = (ioerr == 0)

    ! At this time, only double precision RAMSES data is supported.
    if (.not. double_precision) then
      call die('Only double-precision RAMSES data is supported at this time.')
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
      call die('Only Peano-Hilbert ordered RAMSES data can be read at this time.')
      close(ilun)
      stop
    endif
    read(ilun) overload, ndomain
    backspace(ilun) ! rewind one entry (why?)

    ! read (double or quad precision) Peano-Hilbert keys;
    ! incl. precision and endianness checking.
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
        quad_precision = .false.
      else
        backspace(ilun)
        quad_precision = .true.
        allocate(bound_key_dpx2(2,0:ndomain))
        read(ilun,iostat=ioerr) bound_key_dpx2
        if (any( (abs( real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16) - bound_key) &
               / (abs(bound_key) + 1e-45_16)) > 1e-15_16) ) then
#ifdef __xlc__
          ! We are on BlueGene/Q and using the XLF90 compiler. These two numbers should be identical.
          ! Therefore we are trying to read little endian data, but will never have any luck doing that,
          ! since the little-endian quad-precision format is not supported with XLF.
          call die('We are on an IBM big-endian platform trying to read little endian quad-precision data. This is not possible. Giving up.')
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
        quad_precision = .false.
      else
        backspace(ilun)
        quad_precision = .true.
        allocate(bound_key_dpx2(2,0:ndomain))
        read(ilun,iostat=ioerr) overload, ndomain, bound_key_dpx2, isub_domain, cpu_domain, map_domain
        if (any( (abs( real(bound_key_dpx2(1,:),kind=16) + real(bound_key_dpx2(2,:),kind=16) - bound_key) &
               / (abs(bound_key) + 1e-45_16)) > 1e-15_16) ) then
#ifdef __xlc__
          ! We are on BlueGene/Q and using the XLF90 compiler. These two numbers should be identical.
          ! Therefore we are trying to read little endian data, but will never have any luck doing that,
          ! since the little-endian quad-precision format is not supported with XLF.
          call die('We are on an IBM big-endian platform trying to read little endian quad-precision data. This is not possible. Giving up.')
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

    ! we have extracted enough information; close the file (for the time being).
    close(ilun)

    ! output info about snapshot to terminal
    print *, 'Data is '//merge('double','single',double_precision)//' precision with '// &
             trim(merge('quad  ','double',quad_precision))//' precision hilbert keys'
    print *, 'RAMSES snapshot time = ',time
    print *, 'ncpu = ',ncpu
    print *, 'nlevelmin = ',nlevelmin
    print *, 'nlevelmax = ',nlevelmax
    print *, 'boxlen = ',boxlen
    t_now = time
    print *, 'time = ',time
    unig%t_next = time

    !---- allocate arrays ----
    ! now that we have more information, we can set and allocate a few required variables.
    twotondim = 2**ndims
    maxlvl = min(maxlvl, nlevelmax)
    if (maxlvl == 0) maxlvl = nlevelmax
    ncells = nx * ny * nz + twotondim * ngrids
    allocate(ngrid_file(1:ncpu+nboundary,1:nlevelmax), ngrid_level(1:ncpu,1:nlevelmax), &
             cpu_list(1:ncpu), cpu_read(1:ncpu), xc(1:twotondim,1:3), &
             father(1:ngrids_max), neighbour(1:ngrids_max,1:2*ndims), levels(1:nlevelmax), &
             ncellsum(1:nlevelmax,1:ncpu))
    if (nboundary > 0) allocate(ngrid_bound(1:nboundary,1:nlevelmax))

    !---- given limits, make preparations to cut a sub-volume from the box ----
    ! uses the Hilbert-Peano curve to select which CPUs contain the relevant data
    dmax = max(xmax - xmin, ymax - ymin, zmax - zmin)
    ! for a given cut, what is the bare minimum level we need to consider.
    do ilevel=1,maxlvl
      dx = 0.5_8**ilevel
      if (dx < dmax) exit
    enddo
    bit_length = ilevel - 1
    max_dom = 2**bit_length
    imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
    ndom = 1
    dkey = (2.0_16**(nlevelmax + 1) / real(max_dom,kind=16))**ndims
    if (bit_length > 0) then
      imin = int(xmin * dble(max_dom))
      imax = imin + 1
      jmin = int(ymin * dble(max_dom))
      jmax = jmin + 1
      kmin = int(zmin * dble(max_dom))
      kmax = kmin + 1
      ndom = 8
    endif

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
      levels(ilevel)%imin = int(xmin * dble(nx_full)) + 1
      levels(ilevel)%imax = int(xmax * dble(nx_full)) + 1
      levels(ilevel)%jmin = int(ymin * dble(ny_full)) + 1
      levels(ilevel)%jmax = int(ymax * dble(ny_full)) + 1
      levels(ilevel)%kmin = int(zmin * dble(nz_full)) + 1
      levels(ilevel)%kmax = int(zmax * dble(nz_full)) + 1
    enddo

    !---- create root grid ----
    ! at this point, we also have enough information to create the *root* grid, but
    ! *not* enough to fill the MHD variables with data; we want data to be downloaded
    ! from `nlevelmin` patches (`download` is called in `worker`).
    !call get_new_patch (patch)
    !patch%size = boxlen
    !patch%pos = 0.5 * boxlen
    !patch%n = 2**nlevelmin
    !patch%level = nlevelmin
    !patch%time = time
    !patch%velocity = 0.0
    !call init_patch_task (patch, task)
    !parent => task
    ! leave field variables empty for the time being; will download data later
    !call add_patch_by_quality (patch_list, patch)
    !call append_to_plist (wait_queue, patch)

    ! loop over CPU files
    call barrier_mpi ('ncpu begin')
    do k=1,ncpu_read
      icpu = cpu_list(k)
      ! open `amr_` file
      write(filename, '(a,i5.5)') trim(filebase), icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (big_endian) then
        open(unit=ilun,file=filename,status='old',form='unformatted',convert='BIG_ENDIAN')
      else
        open(unit=ilun,file=filename,status='old',form='unformatted')
      endif
#else
      open(unit=ilun,file=filename,status='old',form='unformatted')
#endif
      if (master) print *, 'Reading file '//trim(filename)

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
      if (rank==mpi_rank) &
        print'(a,i3,30i6)','ngr',k,ngrid_level(k,nlevelmin:nlevelmax)
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
      write(filename, '(a,"/output_",i5.5,"/hydro_",i5.5,".out",i5.5)') trim(datadir), &
            nout, nout, icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (big_endian) then
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
      gamma = gamma_ramses

      ! open the `grav_` file to get access to the hydro data
      write(filename, '(a,"/output_",i5.5,"/grav_",i5.5,".out",i5.5)') trim(datadir), &
            nout, nout, icpu
#if __INTEL_COMPILER || defined(__GFORTRAN__)
      if (big_endian) then
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
      do ilevel=1,maxlvl
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
                                    x(1:ndims,1:ngrid_act), var(1:ngrid_act,1:twotondim,1:nvar_mhd+1), &
                                    ref(1:ngrid_act))

        ! loop over domains
        do j=1,nboundary+ncpu
          ncache = ngrid_file(j,ilevel)

          ! read relevant AMR tree data
          if (ncache > 0) then
            allocate(index_grid(1:ncache), iig(1:ncache))
            read(ilun) index_grid ! read grid index
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
            read(ilun) iig ! read father index
            father(index_grid(1:ncache)) = iig(1:ncache)
            do i=1,2*ndims
              read(ilun) iig ! read neighbour index
              neighbour(index_grid(1:ncache),i) = iig(1:ncache)
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
          if (ngrid_file(j,ilevel) > 0) then
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
          if (ngrid_file(j,ilevel) > 0) then
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
        if (ngrid_act > 0 .and. ilevel >= max(minlvl,nlevelmin)) then
          do l=1,ngrid_act
            ! the grid centre is given by the `xg` variable.
            ! check if cell is in the cut-out region.
            r = sqrt(sum((xg(l,:)-star%position)**2))
            ref(l) = r < size
            pmin = min(pmin,xg(l,:))
            pmax = max(pmax,xg(l,:))
            if (ref(l)) then
              ! create a new DISPATCH patch/task for the current grid and fill with the
              ! cells with data if we own the patch
              npatch = npatch+1
              patch => patch_table(npatch)
              patch%rank = rank
              patch%level = ilevel
              patch%pos(:) = xg(l,:)
              patch%n = 2
              patch%size = 2.0 * dx
              !> if the patch is at the very boundary, freeze it!
              if (r+2d0*dx > size) then
                patch%status = ior(patch%status,frozen_bit)
                if (nprint_f > 0) then
                  nprint_f = nprint_f-1
                  print '(a,i7,a,1p,3e16.9)', &
                    'patch', patch%id, ' frozen', r, dx, size
                end if
              end if
              patch%time = time
              patch%velocity = 0.0
              tmp_load = tmp_load + 1d-6*product(patch%n)*0.5**(nlevelmax-ilevel)
              tot_load = tot_load + 1d-6*product(patch%n)*0.5**(nlevelmax-ilevel)
              if (patch%rank==mpi_rank) then
                call init_patch_task(patch, task)
                if (nprint>0) then
                  nprint = nprint-1
                  print'("task:",i6,i3,1p,3(3e13.5,2x))', &
                    patch%id, patch%level, patch%pos, task%m%p, task%m%s
                end if
                if (mod(npatch,1000)==0) &
                  print'(a,2i7,3x,a,1p,2e10.2)','npatch:',npatch,n_active,'load:',tot_load,tmp_load
                task%offset = patch%pos
                task%parent => parent
                ! loop over cells
                j = 1
                ! RAMSES variable mapping:
                ! 1 = density; 2:4 = zone-centred velocity components;
                ! 5:7 = "left" magnetic field components; 8:10 = "right" magnetic field components;
                ! 11 = thermal pressure; 12:13 = passive scalars.
                ! DISPATCH variable mapping:
                ! 1 = density; 2 = entropy; 3 = grav. potential; 4:6 = momentum components
                ! 7:9 = magnetic field components
                task%d_average = 0.0
                it = task%it
                do iz=task%m(3)%li,task%m(3)%ui
                  do iy=task%m(2)%li,task%m(2)%ui
                    do ix=task%m(1)%li,task%m(1)%ui
                      ! note that ghost cells are left *empty* at this stage.
                      if (var(l,j,1) < d_min) then
                        var(l,j,1) = d_min
                        n_d_min = n_d_min+1
                      end if
                      if (var(l,j,11) < tt_min*var(l,j,1)) then
                        var(l,j,11) = tt_min*var(l,j,1)
                        n_tt_min = n_tt_min+1
                      end if
                      task%mem(ix,iy,iz,1,task%it,1) = var(l,j,1)
                      task%mem(ix,iy,iz,2,task%it,1) = calculate_entropy(gamma_ramses, var(l,j,:))
                      task%mem(ix,iy,iz,3,task%it,1) = var(l,j,nvar_mhd+1)
                      !task%mem(ix,iy,iz,3,task%it,1) = 0.0
                      task%mem(ix,iy,iz,4,task%it,1) = var(l,j,1) * var(l,j,2) ! convert to momentum
                      task%mem(ix,iy,iz,5,task%it,1) = var(l,j,1) * var(l,j,3)
                      task%mem(ix,iy,iz,6,task%it,1) = var(l,j,1) * var(l,j,4)
                      task%mem(ix,iy,iz,7,task%it,1) = var(l,j,5) ! only the "left" fields are used
                      task%mem(ix,iy,iz,8,task%it,1) = var(l,j,6)
                      task%mem(ix,iy,iz,9,task%it,1) = var(l,j,7)
                      task%d_average = task%d_average + var(l,j,1)
                      task%v(:,it) = task%v(:,it) + task%mem(ix,iy,iz,4:6,task%it,1)
                      j = j + 1
                    enddo
                  enddo
                enddo
                task%v(:,it) = task%v(:,it)/task%d_average
                task%d_average = task%d_average/j
                !call scalar_stats (task%m, task%f%d, 'd')
                call add_patch_by_quality(patch_list, patch)
                !call add_patch_by_x(x_list, patch)
                call append_to_plist(wait_queue, patch)

                ! diagnostic info
                ncellsum(ilevel,k) = ncellsum(ilevel,k) + twotondim
                n_onrank=n_onrank+1
              else
                n_offrank=n_offrank+1
              endif
              n_active=n_active+1
            else
              n_passive=n_passive+1
            endif
          enddo
        endif ! end ngrid_act > 0 .and. ilevel >= nlevelmin
        if (ngrid_act > 0) deallocate(son, xg, var, ref, x)

      enddo ! end loop ilevel=1,maxlvl
      call check_quality_order (patch_list, 'loops')

      close(ilun) ! close `amr_` file
      close(ilun2) ! close `hydro_` file
      if (rank==mpi_rank) print'(4(a,i8,3x),a,1p,e9.2)', &
        'n_active:',n_active,'n_passive:',n_passive,'n_onrank:',n_onrank, &
        'n_offrank:',n_offrank,'tmp_load:',tmp_load
      !
      if (tmp_load > rank_load) then
        if (rank<mpi_size-1) then
          if (rank==mpi_rank) &
            print*,'bumping rank, with tmp_load =', tmp_load, n_active
          rank = rank+1
          n_active=0
          n_passive=0
          n_onrank=0
          n_offrank=0
          tmp_load = 0d0
        end if
      end if
      call barrier_mpi ('ncpu loop')
    enddo ! end loop k=1,ncpu_read
    print*,'n_d_min, n_tt_min =',n_d_min,n_tt_min

  !endif ! if (master)
  !$omp end master
  !$omp barrier

  ! set scale factors (other scale factors will be calculated)
  code%l = scale_l
  code%t = scale_t
  code%d = scale_d

  ! output some diagnostics to the terminal.
  allocate (fill(maxlvl))
  do ilevel=nlevelmin,maxlvl
    fill(ilevel)=sum(ncellsum(ilevel,:))/2.0**(3*(ilevel-1))
  enddo
  write (*,'(a,1p,2(3e10.2,2x))') 'pmin, pmax (AU) =', &
    (pmin-star%position)*code%l/cgs%au, (pmax-star%position)*code%l/cgs%au
  write (*,'("fill:",30f7.3)') fill(nlevelmin:maxlvl)
  deallocate (fill)

  deallocate(ngrid_file, ngrid_level, cpu_list, cpu_read, xc, father, &
             neighbour, levels, ncellsum)

  call trace_end
END SUBROUTINE read_ramses_snapshot

!> =============================================================================
!> Lifted directly from `amr2cell.f90`/`hilbert3d_qdp` (in the Copenhagen
!> RAMSES repository.
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
!> Convenience function to calculate the entropy (because DISPATCH tracks the
!> entropy rather than the internal energy/pressure or total energy.)
!> =============================================================================
FUNCTION calculate_entropy (gamma_ramses, vars) RESULT (s)
  implicit none
  real(kind=8):: gamma_ramses, s
  real(kind=8), dimension(:):: vars
  !.............................................................................
  ! vars(0) = rho
  ! vars(10) = pressure
  if (gamma==1d0) then
    s = vars(11)/vars(1)                ! temperature
  else
    s = vars(1)*(log(vars(11))-log(vars(1))*gamma_ramses)/(gamma_ramses-1.0)
  end if
END FUNCTION calculate_entropy

!> =============================================================================
!> Read sink particle information; the zoom-in sink is the last one, generally
!> =============================================================================
 SUBROUTINE read_sink (datadir,nout)
  implicit none
  integer, parameter:: nmetal=2, sink_unit=9
  character(len=*):: datadir
  integer:: nout
  !.............................................................................
  integer, parameter::p8=kind(1.0D0)                                            ! double precision
  integer, parameter::dp=kind(1.0D0)                                            ! double precision
  integer, parameter:: d1=kind(kind(1d0))                                       ! kind of a kind
  type :: sink_t                                                                ! star data type
    SEQUENCE
    real(p8) x,y,z                                                              ! positions
    real(dp) px,py,pz                                                           ! momenta
    real(dp) ux,uy,uz                                                           ! velocities
    real(dp) dpx,dpy,dpz                                                        ! momenta
    real(dp) mass                                                               ! mass
    real(dp) dmass                                                              ! mass
    real(dp) potential                                                          ! gravitational potential
    real(dp) density                                                            ! neigborhood density
    real(dp) t_create                                                           ! creation time
    real(dp) t_explode                                                          ! explosion time
    real(dp) metal(nmetal)                                                      ! metal density
    real(dp) dmetal(nmetal)                                                     ! metal density
    real(dp) beta_create                                                        ! Beta at creation
    integer gid                                                                 ! pointer to grid cell
    integer father                                                              ! pointer to father cell
    integer ilevel                                                              ! level of AMR hierachy where star lives
    integer old_ilevel                                                          ! level of hierachy star used to have, if changed
    integer owner                                                               ! owner thread
    integer has_moved
    real(dp) dtold
  end type
  type(sink_t):: sink
  character(len=80):: fname
  logical:: do_translate
  integer:: center_star, version, nstar, szstar, sizeofstar, i, nstars, nmet, dp1
  real(4):: t,beta_create,dtold
  real(8):: t8,tflush,r_translate(3),v_translate(3),x_refine,y_refine,z_refine
  !.............................................................................
  write (fname,'(a,a,i5.5,a)') trim(datadir), '/output_', nout, '/stars_output.dat'
  print*,fname
  open (sink_unit, file=fname,form='unformatted', status='unknown')
  read (sink_unit) version,nstars,t,szstar,sizeofstar
  print*,'versions =',version
  read (sink_unit) dp1,nmet,t8
  print*,'dp, nmetal =',dp1,nmet
  print*,'t8 =',t8
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
  allocate (star)
  star%position = [sink%x,sink%y,sink%z]
  star%mass = sink%mass
END SUBROUTINE read_sink

END MODULE read_ramses_mod
