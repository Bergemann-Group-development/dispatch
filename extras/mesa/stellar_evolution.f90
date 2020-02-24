!###############################################################################
! Module to provide glue between sink particles in DISPATCH and MESA
MODULE stellar_evolution_mod
  USE amr_commons,    only : myid, ncpu
  USE mesa_star_mod
  implicit none
  !
  integer, dimension(:), allocatable :: owner                                   ! rank that evolves the MESA model for star i
  character(len=80) :: mesadir, inlist
  character(len=255) :: output_dir
  real(kind=8) :: dt_min_for_evolve
  real(kind=8) :: u_t, u_m, u_e
CONTAINS
  !###############################################################################
  ! Init hooks for the stellar evolution code. Only done once.
  SUBROUTINE init_stellar_evolution
    USE amr_parameters, only : dp, datadir
    logical, save :: first_call = .true.
    real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

    if (.not. first_call) return
    first_call = .true.
    mesadir = trim(datadir)//'/mesa'
    !namelist /mesa/ inlist, dt_min_for_evolve
    !inlist = 'inlist'                                                           ! Name of MESA namelist file
    dt_min_for_evolve = 3600.*24.*365.25                                        ! Minimum timestep taken for evolving the stellar structure model
 !..............................................................................
    !rewind(1); read (1,mesa)                                                    ! namelist parameters
    !if (myid==1) write(*,mesa)

    ! ############################
    ! Setup units for interfacing between numerical and cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    u_t = scale_t
    u_m = real(scale_d,kind=8) * real(scale_l,kind=8)**3
    u_e = u_m * real(scale_v,kind=8)**2

  END SUBROUTINE init_stellar_evolution
  !
END MODULE
!###############################################################################
! Evolve all stellar models one time-step forward in time
SUBROUTINE evolve_stellar_models(time,mass,nstars,force_evolution)
  USE stellar_evolution_mod
  implicit none
  integer,      intent(in) :: nstars                                            ! number of sink particles in total
  real(kind=8), intent(in) :: time                                              ! time over which star is evolved
  real(kind=8), intent(in), dimension(nstars) :: mass                           ! current mass
  logical,      intent(in) :: force_evolution                                   ! Force the models to evolve no matter what -- useful for sync of outputs
  !
  real(kind=8) :: dt, dm, mdot, Lacc, Lacc_in, Lacc_out
  integer :: istar, isink

  ! Loop over all sinks, and evolve the sinks this rank owns
  istar = 0

  ! ############################
  do isink=1,nstars
    if (owner(isink)==myid) then
      ! ############################
      istar = istar + 1
      s => star_handles(istar) ! Get star pointer
      print *, 'MESA star_age before evolve', s% star_age*secyer
      dt = (time * u_t) - s% x_ctrl(26) - s% star_age*secyer ! In seconds (time in ramses - birth time in ramses - MESA evolved time ) in seconds
      dm = (mass(isink) * u_m) - s%mstar !  [g], check how much the star has accreted
      mdot = dm / dt ! Average mdot in time since last evolve call
      print *, ' ID: ' , myid
      !############################
      if (s% x_logical_ctrl(1) .eqv. .false.) print *, 'not evolving anymore...'
      print *, 'mdot, dt, dt(years):', dm, dt, dt / secyer
      !############################
      if (dm > 0 .and. s% x_logical_ctrl(1)) then ! if dm < 0 then second core mass not accreted yet
      ! x_logical_ctrl logical for skipping stellar evolution that can be set in mesa_star_mod.f90
      	if (dt > dt_min_for_evolve .or. force_evolution) then
      	  print *, 'calling step_MESA'
      	  call step_MESA(mdot, dt, istar) ! Evolve the model forward in time
      	  !! If evolve_stellar_models is called from write_models force_evolution = .true.
      	  !! Write photos to output_dir
      	endif
      endif
      ! ############################
      if (force_evolution) then ! If true write photo and store in DISPATCH output directory.
      	 print *, 'calling write_photos, istar, isink:', istar, isink
      	 call write_photos(output_dir, istar)
      endif
      ! ############################
    endif
  enddo
  ! ############################
END SUBROUTINE evolve_stellar_models
!###############################################################################
! Create a new stellar model on the least occupied MPI rank
! This model is coupled to the sink index "isink"
SUBROUTINE create_stellar_model(time, isink)
  USE stellar_evolution_mod
  implicit none
  integer, intent(in) :: isink
  integer, allocatable, dimension(:) :: owner_p
  real(kind=8), intent(inout) :: time

  call init_stellar_evolution

  if (allocated(owner)) then
    owner_p = owner
    deallocate(owner)
  endif
  allocate(owner(isink))
  if (allocated(owner_p)) then
    owner(1:isink-1) = owner_p
    deallocate(owner_p)
  endif
  ! ### Now decide which rank should own this sink...
  !
  owner(isink) = modulo(isink-1,ncpu)+1
  !
  if (owner(isink).ne.myid) return
  ! ### Create new stellar model
  print *, myid, ' calling create_MESA_star ...'
  call create_MESA_star(myid, time*u_t)

END SUBROUTINE create_stellar_model
!###############################################################################
! Write MESA restart snapshots
SUBROUTINE write_stellar_models(nstars,dir,time,mass)
  USE stellar_evolution_mod
  implicit none
  integer, intent(in) :: nstars
  character(len=*) :: dir
  character(len=255) :: make_MESA_dir
  real(kind=8), intent(in) :: time
  real(kind=8), intent(in), dimension(nstars) :: mass
  logical :: lexist
  ! ############################
  if (dir=='') then
    print *, 'Called from MESA module - nothing done in for this call at the moment'
  ! ############################
  else
    print *, 'Called from sinkparticles routine - will save latest MESA photo.'
    print *, 'Directory has to be string manipulated: '//trim(dir(:13))

    output_dir = dir(:13)
    make_MESA_dir = 'mkdir -p ./data/' // trim(output_dir) // 'MESA' ! command to make MESA dir in ouput directory

    call execute_command_line(make_MESA_dir) ! make MESA dir in current output directory
    call evolve_stellar_models(time,mass,nstars,.true.) ! Force stellar evolution and move latest photo to MESA directory in output directory.
  endif
  ! ############################

END SUBROUTINE write_stellar_models
!###############################################################################
! Read data for stellar models for restarting code
SUBROUTINE read_stellar_models(time, dir)
  USE stellar_evolution_mod
  implicit none
  integer :: i, reason, nfiles
  real :: r
  character(len=*) :: dir
  real(kind=8), intent(inout) :: time
  character(len=255), dimension(:), allocatable :: photo_list
  character(len=255) :: photo
  integer, allocatable, dimension(:) :: owner_p
  ! ############################
  call init_stellar_evolution

  ! ############################
  ! ### Find MESA photos, count and store in an array
  call execute_command_line('ls '// trim(dir) // 'photo* > list_restartphotos.txt')
  open(22, file = 'list_restartphotos.txt', action="read")
  i = 0
  do
    read(22, fmt = '(a)', iostat = reason) r
    if (reason/=0) exit
    i = i+1
  enddo

  nfiles = i
  if (myid == 1) print *, 'Number of stars in output dir: ', i
  allocate(photo_list(nfiles))
  rewind(22)
  do i = 1, nfiles
    read(22, '(a)') photo_list(i)
  enddo
  close(22)

  ! ############################
  ! ### Assign owners
  if (allocated(owner)) then
    owner_p = owner
    deallocate(owner)
  endif
  allocate(owner(nfiles))
  if (allocated(owner_p)) then
    owner = owner_p
    deallocate(owner_p)
  endif

  do i = 1, nfiles
    if (myid == modulo(i-1,ncpu)+1) then
      owner(i) = myid
      ! Restart stellar model
      print *, myid, ' calling restart_MESA_star ... '
      photo = photo_list(i)
      print *, 'photo: ', photo
      call restart_MESA_star(myid, time*u_t, photo)
    else
      continue
    endif
  enddo

END SUBROUTINE read_stellar_models
!###############################################################################
! Return something about star #i based on text string prop
! Examples could be
!   property = 'radius', 'temperature', 'luminosity'
!
REAL*8 FUNCTION get_stellar_property(property, isink, nstars)
  USE stellar_evolution_mod
  implicit none
  integer :: i, istar, isink, nstars
  character(len=*) :: property
  ! Loop over all sinks, and find the sink (isink)
  istar = 0
  do i=1,nstars
    if (owner(i)==myid) then
      istar = istar + 1
      if (istar == isink) then
      	s => star_handles(istar)
      	if (property == 'mass' .or. property == 'm' .or. property == 'Mass' .or. property == 'M') get_stellar_property = s%mstar
      	if (property == 'radius' .or. property == 'r' .or. property == 'Radius' .or. property == 'R') get_stellar_property = s%photosphere_r
      	if (property == 'Teff' .or. property == 'temperature' .or. property == 'teff') get_stellar_property = s%Teff
      	if (property == 'Lstar') get_stellar_property = s%L_phot
      	return
      endif
    endif
  enddo

END FUNCTION get_stellar_property
