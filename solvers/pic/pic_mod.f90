!===============================================================================
!> The data part of the pic_t derived type must contain all patch specific
!> dynamically  changing data for a patch.  This is on the one hand everything
!> that is permanently allocated via allocate commands in the PPcode, and on the
!> other hand scalars and arrays that change with time, or are location specific.

!> The PIC patch is 1 guard cell smaller, thant the MHD patch, i.e. ng = ng_mhd-1
!> This is necessary, as J and emf have to be computed for the guard zones.
!===============================================================================
MODULE pic_mod
  USE gpatch_mod
  USE io_mod
  USE link_mod 
  USE rand_mod,        ONLY: rand_seed
  USE trace_mod
  USE mesh_mod
  USE boundaries_mod
  USE bits_mod
  USE kinds_mod
  USE vector_mod
  USE vector_ops
  USE scalar_mod
  USE eos_mod
  USE scaling_mod
  ! PIC modules
  USE species_m,       ONLY: spec, sp, field, fields, vthb, nspecies, particle !, counters, n
  USE params,          ONLY: time, dt, tend, master, mpi, mcoord, &
                             nomp, romp, datadir, periodic
  USE pic_io,          ONLY: toutf, toutp
  USE grid_m,          ONLY: g, mesh_type, grid_type
  USE grid_m,          ONLY: Bx, By, Bz, Ex, Ey, Ez, fmhd, B_bnd, E_bnd, J_bnd
  USE maxwell_sources, ONLY: jx,jy,jz, calc_current, calc_charge
  USE random,          ONLY: wi,mt, mti
  implicit none
  private
  ! Linked-list for particle trade
  type, public:: part_list_t
    integer :: np = 0
    class(particle),    pointer :: particle => null()
    class(part_list_t), pointer :: next     => null()
    class(part_list_t), pointer :: prev     => null()
    class(part_list_t), pointer :: tail     => null()
  end type
  !
  type, public:: pic_t
    class(gpatch_t), pointer:: patch
    type(boundaries_t) :: boundaries
    type(pic_scaling_t):: scaling                                               ! quantity conversion between MHD and PIC
    real:: llc(3), extent(3)                                                    ! PIC box
    real(8) :: start_time                                                       ! time when pic should be activated
    type(field),        pointer, dimension(:,:,:,:)      :: fields              ! species on mesh
    real(kind=KindScalarVar), pointer:: mem(:,:,:,:,:)   ! store mhd-pic data 
    real,               pointer, dimension(:,:,:,:)      :: E                   ! store electric field
    real,               pointer, dimension(:,:,:,:)      :: E_old               ! point to it time slot of E for BCs
    real,               pointer, dimension(:,:,:,:)      :: B                   ! point to B
    real,               pointer, dimension(:,:,:,:)      :: B_old               ! point to it time slot of B for BCs
    real,               pointer, dimension(:,:,:,:)      :: J_bnd               ! BCs for J
    type(spec),         pointer, dimension(:)            :: sp                  ! species.f90
    type(part_list_t),  pointer, dimension(:)            :: exported            ! dim(nspecies) 
    real,               pointer, dimension(:,:,:,:,:)    :: vthb                ! species.f90
    real,               pointer, dimension(:,:,:)        :: fmhd                ! grid.f90
    type(grid_type),    pointer                          :: g=>null()           ! grid.f90
    integer(kind = wi), pointer, dimension(:)            :: mti
    integer(kind = wi), pointer, dimension(:,:)          :: mt
    integer,            pointer, dimension(:,:,:)        :: rand_seed
    integer,                     dimension(3)            :: p_lb, p_ub, p_gn    ! PIC domain WRT MHD domain
    real(kind=8),       pointer                          :: dt1
    logical::        on= .false.
    logical:: do_pairs = .true.
    logical:: do_mhd   = .true.
    logical:: first_pic= .true.
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    !
    procedure:: connect
    procedure:: connect_eb
    !procedure:: check_export
    procedure:: export_particle
    procedure:: export_particles
    procedure:: import_particles
    procedure:: needed
    procedure:: impose_particle_boundaries
    procedure:: mhd_to_pic
    procedure:: pic_to_mhd
    procedure:: convert_to_pic
    procedure:: convert_to_mhd
    procedure:: print_histogram
    procedure:: remove_particles
    procedure:: remove_particles2
    procedure:: match_units
  end type
  type(pic_t), public:: pic
  logical, save :: periodic_bcs(3)=.false.
CONTAINS

!===============================================================================
!> Request additional memory slots for the patch
!===============================================================================
SUBROUTINE pre_init (self, link)
  class(pic_t):: self
  class(link_t), pointer :: link
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%pre_init')
  !-----------------------------------------------------------------------------
  ! Make sure there is enough memory slots for the electro-magnetic field
  ! before calling the generic patch initialization
  !-----------------------------------------------------------------------------
  self%patch => gpatch%cast2gpatch(link%task)
  if (self%patch%idx%ex <1) then
    call self%patch%idx%request_index(self%patch%idx%ex, self%patch%nv)
    call self%patch%idx%request_index(self%patch%idx%ey, self%patch%nv)
    call self%patch%idx%request_index(self%patch%idx%ez, self%patch%nv)
    call self%patch%idx%request_index(self%patch%idx%jx, self%patch%nv)
    call self%patch%idx%request_index(self%patch%idx%jy, self%patch%nv)
    call self%patch%idx%request_index(self%patch%idx%jz, self%patch%nv)
  end if
  !-----------------------------------------------------------------------------
  call trace_end
END SUBROUTINE pre_init

!===============================================================================
!> Allocate and initialize the memory and paramters belonging to a PIC patch.
!> This procedure is called once for every patch.  It calls init_species to
!> initialize everything that belongs to the patch being initialized, and the
!> dynamically changing part of this info needs to be copied to the patch
!> data type, to be recovered when a patch is to be updated.
!===============================================================================
SUBROUTINE init (self)
  USE units,     ONLY: c

  class(pic_t):: self
  ! ............................................................................
  class(gpatch_t), pointer :: patch
  real, dimension(:,:,:), pointer:: d, e
  integer, save:: itimer=0
  integer:: isp, l(3), u(3), m(3), lb(3), ub(3), n(3)
  real, save::region_llc(3)=0.0, region_extent(3)=-1.0
  real(8), save:: start_time = 0.0
  logical, pointer:: m2p_mask(:,:,:)
  logical, save:: do_pairs=.true., do_tests=.true., on=.false., do_mhd=.true.
  logical, save:: first_time=.true., first_pic=.true.
  namelist /pic_params/ on, do_pairs, do_mhd, periodic_bcs, do_tests, region_llc, region_extent, start_time
  ! ----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Initialize new PIC data for a patch
  !-----------------------------------------------------------------------------
  call trace_begin('pic_t%pic_mod%init', itimer=itimer)
  self%first_pic = first_time
  call io%lock%set
  if (first_time) then
    stdin = io%input
    rewind (io%input)
    read (io%input, pic_params)                                                 ! rewind file descriptor and read input from file
    if (io%master) write (io%output, pic_params)
    call init_mpi                                                               ! MPI/OMP setup
    nomp = 1
    romp = 1
    call init_io                                                                ! input/output setup
    call init_params                                                            ! simulation parameters
    call init_units                                                             ! physical scaling units
    datadir = io%outputname
    toutf = io%out_time
    toutp = io%out_time
    tend  = io%end_time
    mpi%lb = .true.
    mpi%ub = .true.
    mpi%interior = .false.
    periodic = .false.
  end if
  call io%lock%unset
  self%do_pairs = do_pairs
  self%llc = region_llc
  self%extent = region_extent
  self%start_time = start_time
  self%on = on
  patch => self%patch
  ! ----------------------------------------------------------------------------
  ! Set emf index offsets in patch%mem()
  ! ----------------------------------------------------------------------------
  if (patch%idx%ex >0) then
    if (trim(patch%kind) == 'stagger2e_pic') then
      patch%mesh(1)%h(patch%idx%ex) = +0.5
      patch%mesh(2)%h(patch%idx%ey) = +0.5
      patch%mesh(3)%h(patch%idx%ez) = +0.5
    else
      patch%mesh(2)%h(patch%idx%ex) = -0.5
      patch%mesh(3)%h(patch%idx%ex) = -0.5
      patch%mesh(1)%h(patch%idx%ey) = -0.5
      patch%mesh(3)%h(patch%idx%ey) = -0.5
      patch%mesh(1)%h(patch%idx%ez) = -0.5
      patch%mesh(2)%h(patch%idx%ez) = -0.5
    end if
  end if
  ! ----------------------------------------------------------------------------
  ! Check if PIC is needed for this patch
  !-----------------------------------------------------------------------------
  if (self%on) self%on = self%needed()
  ! ----------------------------------------------------------------------------
  ! If PIC is needed, finalize the initialization
  ! ----------------------------------------------------------------------------
  if (self%on) then
    self%do_mhd = do_mhd                                                        ! MHD should be done in the background?
    !patch%do_mhd = self%do_mhd
    !---------------------------------------------------------------------------
    ! Allocate the PIC grid and match it to the MHD grid.
    !---------------------------------------------------------------------------
    allocate(g)
    self%g => g
    !---------------------------------------------------------------------
    ! Minimum boundaries needed for interpolation operator and field solver
    !---------------------------------------------------------------------
    call interpolation_boundaries
    call field_solver_boundaries
    print *, 'pic needs lb, ub', g%lb, g%ub
    !
    !self%p_lb = patch%ng-g%lb+2
    !self%p_ub = patch%gn-patch%ng+g%ub!patch%ui+g%ub+1
    g%lb = max(g%lb,patch%ng-1)
    g%ub = max(g%ub,patch%ng-2)
    g%gn = patch%mesh%n+3                                                       ! number of active zones
    g%n  = patch%mesh%n+3
    self%p_lb = patch%ng-g%lb
    self%p_ub = self%p_lb+g%lb+g%gn+g%ub-1
    self%p_gn = g%lb+g%gn+g%ub
    !self%p_ub = patch%gn-patch%ng+g%ub!patch%ui+g%ub+1
    !where(patch%n>2)
    !  self%p_lb = patch%mesh%lb+2
    !  self%p_ub = patch%mesh%ub-2
    !endwhere
    g%ds  = patch%ds                                                          ! PIC domain is of the same size as MHD patch
    
  
    ! --------------------------------------------------------------------------
    ! the pic can to go from m%lb+2 : m%ub-3, if m%ng>=4; Additionally, PIC
    ! domain loops go from g%lb to g%ub-1, thus g%ub needs one zone more to
    ! cover the exactly the same domain. The grid initialization calculates
    ! the domains as
    !  g%n(i) = g%n(i) + g%lb(i) + g%ub(i)                                      ! points in domain
    !  g%ub(i) = g%n(i) - g%ub(i)                                               ! upper bdry index
    ! thus initially g%lb should be m%ng-2 and g%ub=m%ng-2
    ! --------------------------------------------------------------------------
    g%grlb = patch%mesh%llc_cart+((self%p_lb-1.0)-0.5_8)*patch%mesh%d
    ! --------------------------------------------------------------------------
    ! Now we have enought parameters to set up the PIC grid.
    ! --------------------------------------------------------------------------
    call init_grid                                                              ! computational mesh
    !
    print *, 'grids.....'
    print *, 'p_lb, ub, ng', self%p_lb, self%p_ub, self%p_gn
    print *, 'rg',g%lb, g%ub
    print *, 'gn',g%n, g%gn
    print *, 'mesh'
    print *,patch%mesh%lb, patch%mesh%ub
    print *,patch%mesh%li, patch%mesh%ui
    print *, 'ds', g%ds, patch%mesh%d
    !
    mpi%offset = -g%lb
    !
    !if (patch%idx%ex<1) then
    n = self%p_gn
    allocate(self%mem(n(1),n(2),n(3), 9,2))
    call io%bits_mem(storage_size(self%mem),product(shape(self%mem)),'pic_mem')
    call allocate_vectors(n,self%J_bnd)
    call io%bits_mem(storage_size(self%J_bnd),product(shape(self%J_bnd)),'pic_j_bnd')
    !end if
    !call self%match_units
    call self%convert_to_pic(patch%it)
     
    call init_species                                                             ! define particle species
    !
    call init_distributions                                                    ! particle distibutions
    call init_random                                                              ! random numbers
    call init_switches                                                          ! execution control
    dt   = 0.0                                                                    ! trigger dt=dtnew
    ! --------------------------------------------------------------------------
    ! Initialize the first PIC state by spawning the whole domain with particles
    ! --------------------------------------------------------------------------
    l = 1; u = g%n
    allocate(m2p_mask(g%n(1), g%n(2),g%n(3)))
    m2p_mask = .true.
    call self%mhd_to_pic (m2p_mask)
    deallocate(m2p_mask)
    !-----------------------------------------------------------------------------
    ! Initialize the PIC state from MHD for patches of PIC type
    !-----------------------------------------------------------------------------
    call init_maxwell                                                           ! initialize maxwell eq
    call init_diagnostics                                                       ! energetics etc
    call write_log                                                              ! initialize log output
    self%fields    => fields
    self%vthb      => vthb                                                      ! species.f90
    self%sp        => sp                                                        ! species.f90
    self%g         => g                                                         ! grid.f90
    self%fmhd      => fmhd                                                      ! grid.f90
    self%mti       => mti
    self%mt        => mt
    self%rand_seed => rand_seed
    self%dt1       => dt
    !call write_io_pic (-1)
    ! --------------------------------------------------------------------------
    ! Allocate export arrays and export particles
    !---------------------------------------------------------------------------
    allocate(self%exported(nspecies))
    !call self%export_particles
    !
    call self%print_histogram ('after init')
       !
    call allocate_sources
    call sort
    call calc_current
    call calc_charge
    call self%convert_to_mhd(patch%it)
    call deallocate_sources
  end if    ! self%on
  if (first_time) first_time=.false.
  !
  call trace%end(itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Currently no pre-update procedures are necessary for pic.
!===============================================================================
SUBROUTINE pre_update (self)
  class(pic_t), target:: self
  ! ----------------------------------------------------------------------------
  
END SUBROUTINE pre_update

!===============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!===============================================================================
SUBROUTINE post_update (self)
  class(pic_t):: self
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  class(gpatch_t), pointer :: patch
  real, dimension(:,:,:,:), pointer :: mem
  real, dimension(:,:,:), pointer :: aux
  integer:: iv, ix, iy, iz, l(3), u(3)
  integer, save:: itimer(2)=0
  ! ............................................................................
  if (.not.self%on) return
  !
  patch => self%patch
  !mem   => patch%mem(:,:,:,:,patch%new,1)
  !-----------------------------------------------------------------------------
  ! Prepare for update, by pointing to the patch data and allocating scratch
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::post_update', itimer=itimer(1))
  !---------------------------------------------------------------------------
  ! Copy mem arrays to new time slot
  !---------------------------------------------------------------------------
  !mem(:,:,:,:) = patch%mem(:,:,:,:,patch%it,1)
  call self%convert_to_pic(patch%it)
  !
  call self%connect
  call self%import_particles
  print *, 'updating pic', self%patch%time
  !-----------------------------------------------------------------------------
  ! Output the state 
  !-----------------------------------------------------------------------------
  call self%output
  !-----------------------------------------------------------------------------
  ! Do the PIC update, using calls to the PPcode routines
  !-----------------------------------------------------------------------------
  if (self%patch%dt_fixed>0.0) then
    self%patch%dtime=self%patch%dt_fixed
    dt = self%patch%dt_fixed
  else
    call set_dt                                                                 ! set timestep dt
    self%patch%dtime = dt                                                       ! update dtime
  end if
  call self%print_histogram ('before move')
  call move_particles                                                           ! move particles
  !call send_particles                                                           ! across patch-bdries   -- it is just a mock subroutine now!!
  call self%print_histogram ('after move')
  call self%impose_particle_boundaries                                          ! boundary conditions
  !
  call self%print_histogram ('after BCs')
  call sort                                                                     ! within patches
  call trace_end (itimer(1))
  call trace_begin ('pic_t%solve_maxwell', itimer=itimer(2))
  call solve_maxwell                                                            ! evolve B and E-fields
  call self%print_histogram ('after after maxwell')
  time = time + dt                                                              ! time advance
  call trace_end (itimer(2)) 
  call self%convert_to_mhd(patch%new)                                           ! compute MHD values
  call self%export_particles
  !
  if (io%debug(2)) then
    associate (l => patch%mesh%li+1, u => patch%mesh%ui-1)
    block
    real, dimension(:,:,:), pointer:: d, e
    d => self%patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%new,1)
    e => self%patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%e,patch%new,1)
    print*,'average and minimum mass density   after pic_to_mhd:',&
            sum(d)/product(u-l+1), minval(d), patch%it
    print*,'average and minimum thermal energy after pic_to_mhd:',&
            sum(e)/product(u-l+1), minval(e), patch%it
    end block
    end associate
  end if
  print *, 'pic dtime set', dt
  !
  if (self%patch%dt_fixed>0) then
    self%patch%dtime=self%patch%dt_fixed
  else
    self%patch%dtime = dt                                                       ! update dtime
  end if
  !
  ! ----------------------------------------------------------------------------
  ! Auxilary variables
  ! ----------------------------------------------------------------------------
  !call allocate_scalars(self%patch%gn,aux)
  !l=self%p_lb
  !u=self%p_ub
  !aux=0.0
  print *, 'minmax jx', minval(jx), maxval(jx)
  print *, 'minmax jy', minval(jy), maxval(jy)
  print *, 'minmax jz', minval(jz), maxval(jz)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = ex
  !call self%patch%aux%register ('ex', aux)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = ey
  !call self%patch%aux%register ('ey', aux)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = ez
  !call self%patch%aux%register ('ez', aux)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = jx
  !call self%patch%aux%register ('jx', aux)
  !print *, 'jx', jx(6:10,10,10)
  !print *, 'jy', jy(6:10,10,10)
  !print *, 'jz', jz(6:10,10,10)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = jy
  !call self%patch%aux%register ('jy', aux)
  !aux(l(1):u(1),l(2):u(2),l(3):u(3)) = jz
  !call self%patch%aux%register ('jz', aux)
  !call deallocate_scalars(aux)
  ! ----------------------------------------------------------------------------
  call deallocate_sources
END SUBROUTINE post_update

!===============================================================================
!> Prepare the state of the PIC patch to take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!===============================================================================
SUBROUTINE connect (self)
  class(pic_t):: self
  integer     :: l(3), u(3)
  !-----------------------------------------------------------------------------
  ! Prepare for access, by pointing to the patch data
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::connect')
  fields => self%fields                                 ! species.f90
  vthb   => self%vthb                                   ! species.f90
  sp     => self%sp                                     ! species.f90
  g    => self%g                                        ! grid.f90
  !-----------------------------------------------------------------------------
  ! Copy mem arrays to new time slot and update pointers
  !-----------------------------------------------------------------------------
  fmhd => self%fmhd                                                             ! grid.f90
  mti  => self%mti
  mt   => self%mt
  rand_seed => self%rand_seed
  dt => self%dt1
  call allocate_sources
  !-----------------------------------------------------------------------------
  ! A few scalar data items need to be copied, to avoid having to make pointers
  !-----------------------------------------------------------------------------
  time   = self%patch%time
  call trace_end
END SUBROUTINE connect

!===============================================================================
!>
!===============================================================================
SUBROUTINE connect_eb(self, it, it2)
  class(pic_t):: self
  integer:: it, it2                            ! new mem slot, old mem slot for boundaries
  ! ............................................................................
  class(gpatch_t), pointer :: patch
  integer:: l(3), u(3)
  ! ----------------------------------------------------------------------------
  patch => self%patch
  l = self%p_lb
  u = self%p_ub
  self%B = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%bx:patch%idx%bz,it,1) &
         * self%scaling%convert%b
  Bx => self%B(:,:,:,1)
  By => self%B(:,:,:,2)
  Bz => self%B(:,:,:,3)
  !
  self%E = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%ex:patch%idx%ez,it,1) &
           * self%scaling%convert%emf
  Ex => self%E(:,:,:,1)
  Ey => self%E(:,:,:,2)
  Ez => self%E(:,:,:,3)
  ! Connect it time slot to use for boundaries
  self%B_old = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%bx:patch%idx%bz,it2,1) &
         * self%scaling%convert%b
  self%E_old = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%ex:patch%idx%ez,it2,1) &
         * self%scaling%convert%emf
  B_bnd => self%B_old
  E_bnd => self%E_old
  J_bnd => self%J_bnd
END SUBROUTINE connect_eb

!===============================================================================
!>
!===============================================================================
SUBROUTINE convert_to_pic(self, it)
  class(pic_t):: self
  integer:: it                                                                  ! new mem slot
  ! ............................................................................
  class(gpatch_t), pointer :: patch
  real, dimension(:,:,:,:), pointer:: p, b, emf, sub, du
  real, dimension(:,:,:), pointer:: d, e
  integer:: l(3), u(3)
  ! ----------------------------------------------------------------------------
  patch => self%patch
  l = self%p_lb
  u = self%p_ub
  !
  d   => patch%mem(:,:,:,patch%idx%d,it,1)
  e   => patch%mem(:,:,:,patch%idx%e,it,1)
  p   => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,it,1)
  b   => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,it,1)
  emf => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,it,1)
  !
  call allocate_vectors(self%p_gn,du)
  du = up1(d)
  select case(trim(patch%kind))
  case('stagger2e_pic')
    self%mem(:,:,:,1:3,1) =  p(l(1):u(1),l(2):u(2),l(3):u(3),:) / &
                            du(l(1):u(1),l(2):u(2),l(3):u(3),:)
    self%mem(:,:,:,4:6,1) = b(l(1):u(1),l(2):u(2),l(3):u(3),:)
    self%mem(:,:,:,7:9,1) = emf(l(1):u(1),l(2):u(2),l(3):u(3),:)
    !
    self%mem(:,:,:,:,2)   = self%mem(:,:,:,:,1)
  case default
    call allocate_vectors(patch%gn,sub)
    sub = super_up1(p)
    self%mem(:,:,:,1:3,1) = sub(l(1):u(1),l(2):u(2),l(3):u(3),:) / &
                             du(l(1):u(1),l(2):u(2),l(3):u(3),:)
    sub = super_up1(b)
    self%mem(:,:,:,4:6,1) = sub(l(1):u(1),l(2):u(2),l(3):u(3),:)
    sub = super_up1(emf)
    self%mem(:,:,:,7:9,1) = sub(l(1):u(1),l(2):u(2),l(3):u(3),:)
    !
    self%mem(:,:,:,:,2)   = self%mem(:,:,:,:,1)
    call deallocate_vectors(sub)
  end select
  Bx => self%mem(:,:,:,4,1)
  By => self%mem(:,:,:,5,1)
  Bz => self%mem(:,:,:,6,1)
  B_bnd => self%mem(:,:,:,4:6,2)
  Ex => self%mem(:,:,:,7,1)
  Ey => self%mem(:,:,:,8,1)
  Ez => self%mem(:,:,:,9,1)
  E_bnd => self%mem(:,:,:,7:9,2)
  J_bnd => self%J_bnd
  call deallocate_vectors(du)
END SUBROUTINE convert_to_pic

!===============================================================================
!>
!===============================================================================
SUBROUTINE convert_to_mhd(self, it)
  class(pic_t):: self
  integer:: it                                                                  ! new mem slot
  ! ............................................................................
  class(gpatch_t), pointer :: patch
  real, dimension(:,:,:,:), pointer:: p, b, emf, sub, du
  real, dimension(:,:,:), pointer:: d,e
  integer:: l(3), u(3)
  ! ----------------------------------------------------------------------------
  call self%pic_to_mhd(it)
  !
  patch => self%patch
  l = self%p_lb
  u = self%p_ub
  !
  d   => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,it,1)
  p   => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%px:patch%idx%pz,it,1)
  b   => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%bx:patch%idx%bz,it,1)
  emf => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%ex:patch%idx%ez,it,1)
  call allocate_vectors(self%p_gn,du)
  du = up1(d)
  select case(trim(patch%kind))
  case('stagger2e_pic')
    p   = self%mem(:,:,:,1:3,1)*d
    b   = self%mem(:,:,:,4:6,1)
    emf = self%mem(:,:,:,7:9,1)
  case default
    call allocate_vectors(self%p_gn,sub)
    sub = super_down1(self%mem(:,:,:,1:3,1)*d)
    p   = sub
    sub = super_down1(self%mem(:,:,:,4:6,1))
    b   = sub
    sub = super_down1(self%mem(:,:,:,7:9,1))
    emf = sub
    call deallocate_vectors(sub)
  end select
  call deallocate_vectors(du)
  !
END SUBROUTINE convert_to_mhd

!===============================================================================
!> Organize disk output.
!===============================================================================
SUBROUTINE output (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! Connect to the patch data and do I/O, using calls to the PPcode routine
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::output')
  if (associated(self%g)) then
  !  call write_io_pic (self%patch%iout)                                         ! conditionally, do io
  end if
  call trace_end
END SUBROUTINE output

!===============================================================================
!> Check if PIC is needed for this patch. Should be expanded with more conditions.
!> Also adds the boundary bits if they are needed
!===============================================================================
logical FUNCTION needed(self)
  class(pic_t):: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3), pic_urc(3), eps
  integer:: i, lbit(3), ubit(3)
  ! ............................................................................
  eps=1e-5
  llc = self%patch%position-0.5_8*self%patch%size+eps
  urc = self%patch%position+0.5_8*self%patch%size-eps
  pic_urc = self%llc+self%extent
  needed = .true.
  !
  lbit = [bits%xl, bits%yl, bits%zl]                                            ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]                                            ! upper boundary bits -- convenience
  !  
  do i = 1,3
    eps = self%patch%ds(i)
    if (self%extent(i) > 0.) then
      if ((urc(i) < self%llc(i))) then
        needed = .false.
      else if (llc(i) <= self%llc(i)+eps) then
        call self%boundaries%set (lbit(i))
      end if
      !
      if ((llc(i) > pic_urc(i))) then
        needed = .false.
      else if (urc(i) >= pic_urc(i)-eps) then
        call self%boundaries%set (ubit(i))
      end if
    end if
  end do

END FUNCTION needed

!===============================================================================
!> Print pic historgram
!===============================================================================
SUBROUTINE print_histogram (self, label)
  class(pic_t)                            :: self
  character(len=*)                        :: label
  integer                                 :: ip, isp, iz, l(2), u(2)
  type(particle), pointer                 :: pa
  integer, dimension(:,:), pointer        :: n
  real, dimension(:,:,:), pointer         :: v
  real                                    :: d(2)
  !-----------------------------------------------------------------------------
  !if (.not.io%debug(2)) return
  allocate (n(0:g%n(3)+1,nspecies))
  allocate (v(0:g%n(3)+1,3,nspecies))
  n = 0
  v = 0.0
  do isp=1,nspecies
    do ip=1,sp(isp)%np
      pa => sp(isp)%particle(ip)
      if (pa%w==0.0) cycle
      iz = pa%q(3) - mpi%offset(3)
      iz = max(0,min(iz,g%n(3)+1))
      n(iz,isp) = n(iz,isp) + 1
      v(iz,:,isp) = v(iz,:,isp) + pa%p
    end do
  end do
  print *, 'hist: '//trim(label)
  l = g%lb(1:2)
  u = g%ub(1:2)
  do iz=1,g%n(3)
    do isp=1,nspecies
      if (n(iz,isp)>0) v(iz,:,isp) = v(iz,:,isp)/n(iz,isp)
      if (iz>0.and.iz<=self%patch%gn(3)) then
        if (isp==1) then
          d(isp) = sum(bz(l(1):u(1),l(2):u(2),iz)**2)/product(u-l+1)
        else
          d(isp) = sum(ey(l(1):u(1),l(2):u(2),iz)**2)/product(u-l+1)
        end if
      else
        d(isp) = 0.0
      end if
    end do
    print '(i4,2i8,2x,1p,2e12.3,0p,2(2x,1p,3f9.3),3x,a)', &
      iz, n(iz,:), d, v(iz,:,:), ' hist: '//trim(label)
  end do
  deallocate (n,v)
END SUBROUTINE print_histogram


!=====================================================================!
! Non periodic boundaries
!=====================================================================!
SUBROUTINE impose_particle_boundaries(self)
  class(pic_t)                            :: self
  ! ----------------------------------------------------------------------------
  class(mesh_t), pointer                  :: m(:)
  integer                                 :: l(3),u(3)
  integer                                 :: lr(3),ur(3),la(3),ua(3)
  integer                                 :: lo(3), uo(3), li(3), ui(3)
  integer:: itimer=0
  logical, pointer:: m2p_mask(:,:,:)
  ! ............................................................................
  call trace_begin ('pic_t%impose_particle_boundaries',itimer=itimer)
  !
  m=>self%patch%mesh
  lr = 0;    ur = g%n+1
  la = 1;    ua = g%n
  lo=g%lb-1; uo = g%ub
  ! ----------------------------------------------------------------------------
  ! create a mask where boundary conditions should apply and remove particles
  ! in these areas
  ! ----------------------------------------------------------------------------
  allocate(m2p_mask(g%n(1), g%n(2),g%n(3)))
  m2p_mask = .false.
  if (self%boundaries%is_set(bits%xl)) then
    l = [lr(1), lr(2), lr(3)]
    u = [lo(1), ur(2), ur(3)]
    call self%remove_particles (l, u)
    l = [la(1), la(2), la(3)]
    u = [lo(1), ua(2), ua(3)]
    m2p_mask(la(1):lo(1),la(2):ua(2),la(3):ua(3)) = .true.
  end if
  !  
  if (self%boundaries%is_set(bits%xu)) then
    l = [uo(1), lr(2), lr(3)]
    u = [ur(1), ur(2), ur(3)]
    call self%remove_particles (l, u)
    l = [uo(1), la(2), la(3)]
    u = [ua(1), ua(2), ua(3)]
    m2p_mask(uo(1):ua(1),la(2):ua(2),la(3):ua(3)) = .true.
  end if
  !----------------------------------------------------------------------!
  if (self%boundaries%is_set(bits%yl)) then
    l = [lr(1), lr(2), lr(3)]
    u = [ur(1), lo(2), ur(3)]
    call self%remove_particles (l, u)
    l = [la(1), la(2), la(3)]
    u = [ua(1), lo(2), ua(3)]
    m2p_mask(la(1):ua(1),la(2):lo(2),la(3):ua(3)) = .true.
  end if
  !
  if (self%boundaries%is_set(bits%yu)) then    
    l = [lr(1), uo(2), lr(3)]
    u = [ur(1), ur(2), ur(3)]
    call self%remove_particles (l, u)
    l = [la(1), uo(2), la(3)]
    u = [ua(1), ua(2), ua(3)]
    m2p_mask(la(1):ua(1),uo(2):ua(2),la(3):ua(3)) = .true.
  end if
  !----------------------------------------------------------------------!
  if (self%boundaries%is_set(bits%zl)) then
    l = [lr(1), lr(2), lr(3)]
    u = [ur(1), ur(2), lo(3)]
    call self%remove_particles (l, u)
    l = [la(1), la(2), la(3)]
    u = [ua(1), ua(2), lo(3)]
    m2p_mask(la(1):ua(1),la(2):ua(2),la(3):lo(3)) = .true.
  end if
  !....................................................................!
  if (self%boundaries%is_set(bits%zu)) then
    l = [lr(1), lr(2), uo(3)]
    u = [ur(1), ur(2), ur(3)]
    call self%remove_particles (l, u)
    l = [la(1), la(2), uo(3)]
    u = [ua(1), ua(2), ua(3)]
    m2p_mask(la(1):ua(1),la(2):ua(2),uo(3):ua(3)) = .true.
  end if
  ! ----------------------------------------------------------------------------
  ! Impose the BCs where needed
  ! ----------------------------------------------------------------------------
  call self%mhd_to_pic(m2p_mask)
  deallocate(m2p_mask)
  call trace_end (itimer)
END SUBROUTINE impose_particle_boundaries

!===============================================================================
!> Convert pic variables to mhd and store them in the self%new time slot
!===============================================================================
SUBROUTINE pic_to_mhd (self, it)
  class(pic_t):: self
  integer :: it
  ! ----------------------------------------------------------------------------
  class(gpatch_t), pointer :: patch
  class(mesh_t),   pointer :: m(:)
  integer:: isp, i, ix, iy, iz, l(3), u(3), li(3), lb(3), ui(3), ub(3)
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:,:), pointer:: du
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace_begin ('pic_t%pic_to_mhd', itimer=itimer)
  !
  patch => self%patch
  m => self%patch%mesh
  !-----------------------------------------------------------------------------
  ! Copy over fields (B is omitted since is using pointers to self%mem)
  !-----------------------------------------------------------------------------
  call allocate_scalars(g%n, d, e, px, py, pz)
  call allocate_vectors(g%n, du)
  
  d = 0.0
  e = 0.0
  px=0.0
  py=0.0
  pz=0.0
  !
  do isp=1,nspecies
    do iz=1,g%n(3)
      do iy=1,g%n(2)
        do ix=1,g%n(1)
           d(ix,iy,iz) =  d(ix,iy,iz) + sp(isp)%mass*fields(ix,iy,iz,isp)%d
           e(ix,iy,iz) =  e(ix,iy,iz) + sp(isp)%mass*fields(ix,iy,iz,isp)%d &
                                      * fields(ix,iy,iz,isp)%vth**2 * 0.5
          px(ix,iy,iz) = px(ix,iy,iz) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(1)
          py(ix,iy,iz) = py(ix,iy,iz) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(2)
          pz(ix,iy,iz) = pz(ix,iy,iz) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(3)
        end do
      end do
    end do
  end do

  du = up1(d)
  !
  l = self%p_lb
  u = self%p_ub
  !
  patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d ,it,1) = d
  patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%e ,it,1) = e
  !
  self%mem(:,:,:,1,1) = px/du(:,:,:,1)
  self%mem(:,:,:,2,1) = py/du(:,:,:,2)
  self%mem(:,:,:,3,1) = pz/du(:,:,:,3)
  !
  patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%jx,it,1) = jx  
  patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%jy,it,1) = jy 
  patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%jz,it,1) = jz 
  !
  call deallocate_scalars(d, e, px, py, pz)
  call deallocate_vectors(du)
  !-----------------------------------------------------------------------------
  call trace_end (itimer)
END SUBROUTINE pic_to_mhd

!===============================================================================
!> Convert mhd variables to pic
!> For the conversion, we will need to compute J and emf, thus full size grid
!> of B as well as cell-centred u.
!===============================================================================
SUBROUTINE mhd_to_pic (self, m2p_mask)
  ! check these
  USE units,              ONLY: c
  USE particleoperations, ONLY: make_particles_with_temperature_3d, &
                                make_particle_pairs_with_temperature_3d, &
                                make_boundary_particles_with_temperature, &
                                make_boundary_particle_pairs_with_temperature, &
                                make_nonrelativistic_particle
  ! ----------------------------------------------------------------------------
  class(pic_t):: self
  logical:: m2p_mask(:,:,:)
  ! ............................................................................
  class(gpatch_t), pointer:: patch
  class(mesh_t),   pointer:: m(:)
  real :: msum, vex, vey, vez, weight, vg(3)
  integer :: isp, ix, iy, iz, i, l(3), u(3), n(3)
  real :: mean_mass
  real :: xupd, yupd, zupd
  real, dimension(:,:,:), pointer :: d, nd, e, ux, uy, uz
  real, dimension(:,:,:), pointer :: xdu, ydu, zdu, jxs, jys, jzs
  real, dimension(:,:,:), pointer :: ex_l, ey_l, ez_l
  real, dimension(:,:,:), pointer :: bx_l, by_l, bz_l
  real, dimension(:,:,:), pointer :: d_pic, temp
  real, dimension(:,:,:,:), pointer :: P, emf, sub
  real, dimension(:,:,:,:), pointer :: du, B, J, emf_pic, U_pic, B_pic
  real, dimension(:,:,:,:), pointer :: B_loc, emf_loc, U_loc, J_loc
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace_begin ('pic_t%mhd_to_pic', itimer=itimer)
  patch => self%patch
  m     => patch%mesh
  !
  l = self%p_lb; u = self%p_ub
  ! ----------------------------------------------------------------------------
  ! Allocate scratches and point 
  ! ----------------------------------------------------------------------------
  !call allocate_vectors(self%patch%gn, du, uu, J, B_pic, emf)                   ! patch-size
  !call allocate_scalars(g%n,d_pic, e_pic,temp,jxs,jys,jzs)                      ! PIC size
  d   => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  e   => patch%mem(:,:,:,patch%idx%e,patch%it,1)
  p   => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  b   => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%it,1)
  emf => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,patch%it,1)
  call allocate_vectors(patch%gn, du, U_loc, B_loc, emf_loc, J)
  call allocate_scalars(g%n, temp)
  du = up1(d)
  select case(trim(patch%kind))
  case('stagger2e_pic')
    U_loc   = P
    B_loc   = B
    emf_loc = emf
  case default
    call allocate_vectors(patch%gn,sub)
    sub     = super_up1(p)
    U_loc   = sub
    sub     = super_up1(b)
    b_loc   = sub
    sub     = super_up1(emf)
    emf_loc = sub
    call deallocate_vectors(sub)
  end select
  ! ----------------------------------------------------------------------------
  ! Convert momentum to velocity
  ! ----------------------------------------------------------------------------
  U_loc = U_loc / du
  ! ----------------------------------------------------------------------------
  ! sum of different particle masses and mean mass
  ! ----------------------------------------------------------------------------
  msum      = sum(sp(1:nspecies)%mass)                                          
  mean_mass = msum/2.
  ! ----------------------------------------------------------------------------
  ! density
  ! ----------------------------------------------------------------------------
  d_pic => d(l(1):u(1),l(2):u(2),l(3):u(3))
  do isp=1, nspecies
    do iz = l(3), u(3)
      do iy = l(2), u(2)
        do ix = l(1), u(1)
          if (m2p_mask(ix,iy,iz)) &
             fields(ix,iy,iz,isp)%d = d_pic(ix,iy,iz)/msum
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Calculation of electric current, which is needed to get the relative speeds
  ! of electrons and ions, and the particle temperature
  !-----------------------------------------------------------------------------
  J = curl_dn(self%patch%ds,B_loc)/c%mu0
  ! ----------------------------------------------------------------------------
  ! Velocity, up-staggered density and temperature
  ! ----------------------------------------------------------------------------
  ux   => U_loc(l(1):u(1),l(2):u(2),l(3):u(3),1)
  uy   => U_loc(l(1):u(1),l(2):u(2),l(3):u(3),2)
  uz   => U_loc(l(1):u(1),l(2):u(2),l(3):u(3),3)
  temp = mean_mass*e(l(1):u(1),l(2):u(2),l(3):u(3)) * &
         (self%patch%gamma-1.0)/(c%kb*d_pic)
  ! ----------------------------------------------------------------------------
  ! electric field
  ! ----------------------------------------------------------------------------
  !
  !emf(:,:,:,1) = xUP(zDN(Uz))*zDN(By_l) - xUP(yDN(Uy))*yDN(Bz_l)                    
  !emf(:,:,:,2) = yUP(xDN(Ux))*xDN(Bz_l) - yUP(zDN(Uz))*zDN(Bx_l)                    
  !emf(:,:,:,3) = zUP(yDN(Uy))*yDN(Bx_l) - zUP(xDN(Ux))*xDN(By_l)
  !Ex_l => emf(l(1):u(1),l(2):u(2),l(3):u(3),1)
  !Ey_l => emf(l(1):u(1),l(2):u(2),l(3):u(3),2)
  !Ez_l => emf(l(1):u(1),l(2):u(2),l(3):u(3),3)
  !Ex(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = Ex_l(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))                    
  !Ey(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = Ey_l(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))                    
  !Ez(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = Ez_l(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))
  !E_bnd(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),1) = Ex(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))                   
  !E_bnd(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),2) = Ey(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))                   
  !E_bnd(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),3) = Ez(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))   
  !
  do i = 1,3
    self%mem(:,:,:,7+i-1,1) = merge(emf_loc(l(1):u(1),l(2):u(2),l(3):u(3),i), &
    self%mem(:,:,:,7+i-1,1),m2p_mask)
  end do                
  ! ----------------------------------------------------------------------------
  ! point sub-arrays
  ! ----------------------------------------------------------------------------
  xdu   => du(l(1):u(1),l(2):u(2),l(3):u(3),1)
  ydu   => du(l(1):u(1),l(2):u(2),l(3):u(3),2)
  zdu   => du(l(1):u(1),l(2):u(2),l(3):u(3),3)
  J_loc =>  J(l(1):u(1),l(2):u(2),l(3):u(3),:)
  ! ----------------------------------------------------------------------------
  ! set J boundary
  ! ----------------------------------------------------------------------------
  self%J_bnd = J_loc
  !-----------------------------------------------------------------------------
  ! Prepare calls to make_particles_with_temperature_3d for each species, or
  ! the corresponding one for pairs.  Each particle is assigned a weight equal
  ! to fields(ix,iy,iz,isp)%d / sp(isp)%ntarget inside those subroutines.
  !-----------------------------------------------------------------------------
  do isp=1, nspecies
    !---------------------------------------------------------------------------
    ! To conserve the linear momentum when adding the electic current velocities
    ! they have to be weighted by the (normalized) inverse mass. If, for example,
    ! the electrons are much lighter than the ions, they will provide most of the
    ! electic current, while the ions provide the opposite momentum to cancel the
    ! one from the electrons.  This also follows from how they are accelerated by
    ! an electric field.
    !---------------------------------------------------------------------------
    weight = 1.0/sp(isp)%mass / sum(1.0/sp%mass)
    if (sp(isp)%charge/=0) then
      do iz=l(3), u(3)
        do iy=l(2), u(2)
          do ix=l(1), u(1)
            if (m2p_mask(ix,iy,iz)) then
              vex = weight*j_loc(ix,iy,iz,1)*(msum/(xdu(ix,iy,iz)*sp(isp)%charge))   ! electric current vel
              vey = weight*j_loc(ix,iy,iz,2)*(msum/(ydu(ix,iy,iz)*sp(isp)%charge))
              vez = weight*j_loc(ix,iy,iz,3)*(msum/(zdu(ix,iy,iz)*sp(isp)%charge))
              ! --- total velocity ---
              fields(ix,iy,iz,isp)%v(1) = (ux(ix,iy,iz) + vex)
              fields(ix,iy,iz,isp)%v(2) = (uy(ix,iy,iz) + vey)
              fields(ix,iy,iz,isp)%v(3) = (uz(ix,iy,iz) + vez)
            end if
          end do
        end do
      end do
    else
      do iz=l(3), u(3)
        do iy=l(2), u(2)
          do ix=l(1), u(1)
            if (m2p_mask(ix,iy,iz)) then
              ! --- velocity---
              fields(ix,iy,iz,isp)%v(1) = ux(ix,iy,iz)
              fields(ix,iy,iz,isp)%v(2) = uy(ix,iy,iz) 
              fields(ix,iy,iz,isp)%v(3) = uz(ix,iy,iz)
            end if
          end do
        end do
      end do
    end if
  end do
  !----------------------------------------------------------------------------!
  ! Since make_particle_pairs needs velocitites for both species, the calls 
  ! below must be in a separate loop over species
  !----------------------------------------------------------------------------!
  do isp=1, nspecies
    if (self%do_pairs) then
      if (mod(isp, 2) /= 0) then
        call make_particle_pairs_with_temperature_3d(isp, isp+1, temp, m2p_mask)
      end if
    else
      call make_particles_with_temperature_3d (isp, temp, m2p_mask)
    end if
  end do
  !
  call deallocate_vectors(du, U_loc, B_loc, emf_loc, J)
  call deallocate_scalars(temp)                  
  call trace_end (itimer)
END SUBROUTINE mhd_to_pic

!===============================================================================
!> Place a particle to an export list
!===============================================================================
SUBROUTINE export_particle(self, p, isp)
  class(pic_t) :: self
  class(particle):: p
  integer:: isp                                                                 ! species index
  logical :: export
  ! ----------------------------------------------------------------------------
  class(particle),    pointer:: ep
  class(part_list_t), pointer:: lp, this
  ! ............................................................................
  !  
  ! TODO: This takes only the main data. Particle history has to be manually added.
  allocate(ep, source=p)
  !
  ! append the particle
  if (associated(self%exported(isp)%particle)) then
    allocate(lp)
    lp%particle => ep
    if (associated(self%exported(isp)%tail)) then
      lp%prev => self%exported(isp)%tail
      lp%prev%next => lp
    else
      lp%prev => self%exported(isp)
      lp%prev%next => lp
    end if
    self%exported(isp)%tail => lp
    self%exported(isp)%np = self%exported(isp)%np + 1
  ! otherwise start the linked list
  else
    self%exported(isp)%particle => ep
    self%exported(isp)%np = self%exported(isp)%np + 1
  end if
  ! mark particle as removed
  p%w = 0

END SUBROUTINE export_particle

!===============================================================================
!> Remove particles in 3-D region
!===============================================================================
SUBROUTINE export_particles (self)
  class(pic_t):: self
  ! ----------------------------------------------------------------------------
  integer:: q(3), ip, isp, lr(3),ur(3),la(3),ua(3), l(3), u(3), lo(3), uo(3)
  integer:: minq(3), maxq(3), minq0(3), maxq0(3), exported
  class(particle), pointer, dimension(:)   :: pa
  class(mesh_t), pointer:: m(:)
  ! ............................................................................
  !
  m => self%patch%mesh
  lo = m%lo;  uo = m%uo
  lr = m%lb;  ur = m%ub;  la = m%lb+1;  ua = m%ub-1  
  ! ----------------------------------------------------------------------------
  ! Export particles wherever needed
  ! ----------------------------------------------------------------------------
  ! X boundaries
  if (self%boundaries%is_clear(bits%xl)) then
  print *, 'export xl'
    l = [lr(1), m(2)%li, m(3)%li]
    u = [lo(1), m(2)%ui, m(3)%ui]
    call export(l,u)
  end if
  if (self%boundaries%is_clear(bits%xu)) then
  print *, 'export xu'
    l = [uo(1), m(2)%li, m(3)%li]
    u = [ur(1), m(2)%ui, m(3)%ui]
    call export(l, u)
  end if
  ! Y boundaries
  if (self%boundaries%is_clear(bits%yl)) then
  print *, 'export yl'
    l = [m(1)%li, lr(2), m(3)%li]
    u = [m(1)%ui, lo(2), m(3)%ui]
    call export(l,u)
  end if
  if (self%boundaries%is_clear(bits%yu)) then
  print *, 'export yu'
    l = [m(1)%li, uo(2), m(3)%li]
    u = [m(1)%ui, ur(2), m(3)%ui]
    call export(l, u)
  end if
  ! Z boundaries
  if (self%boundaries%is_clear(bits%zl)) then
  print *, 'export zl'
    l = [m(1)%li, m(2)%li, lr(3)]
    u = [m(1)%ui, m(2)%ui, lo(3)]
    call export(l,u)
  end if
  if (self%boundaries%is_clear(bits%zu)) then
  print *, 'export zu'
    l = [m(1)%li, m(2)%li, uo(3)]
    u = [m(1)%ui, m(2)%ui, ua(3)]
    call export (l, u)
  end if
   
  contains
  subroutine export (ll,uu)
    integer:: ll(3), uu(3), i
    real(8):: lc(3), uc(3), ppos(3)
  do i=1,3
    lc(i) = self%patch%mesh(i)%r(ll(i))-0.5_8*self%patch%ds(i)
    uc(i) = self%patch%mesh(i)%r(uu(i))+0.5_8*self%patch%ds(i)
  end do
  exported = 0
  minq = 2**20
  maxq = -minq
  minq0 = 2**20
  maxq0 = -minq
  do isp = 1,nspecies
    pa => sp(isp)%particle
    do ip=1,sp(isp)%np
      !q = pa(ip)%q - mpi%offset
      do i=1,3
        ppos = self%patch%mesh(i)%r(pa(ip)%q(i) - mpi%offset(i))+0.5_8*self%patch%ds(i) + pa(ip)%r(i)* self%patch%ds(i)
      end do
      !if (all(q >= ll .and. q <= uu)) then
      if (all(ppos >= lc .and. ppos <= uc)) then
        call self%export_particle(pa(ip), isp)
        exported = exported + 1
      else if (pa(ip)%w > 0.0) then
        minq = min(minq,q)
        maxq = max(maxq,q)
      else
        minq0 = min(minq0,q)
        maxq0 = max(maxq0,q)
      end if
    end do
  end do
    print*, "Particles exported, in export queue, remaining ", exported, self%exported%np, sp%np
  end subroutine export
END SUBROUTINE export_particles

!===============================================================================
!> Import particles from nbors, if they have active PIC
!===============================================================================
SUBROUTINE import_particles(self)
  class(pic_t) :: self
  ! ----------------------------------------------------------------------------
  class(nbor_t), pointer:: nb_link
  class(gpatch_t), pointer:: nbor
  class(*), pointer:: nbp
  class(pic_t), pointer:: nb_pic
  class(particle),    pointer:: ep
  class(part_list_t), pointer:: plist, remove
  integer:: isp, i, imported
  real(8):: ppos(3), r(3)
  integer:: cpos(3)
  ! ............................................................................
  !
  print *, "Particles before import", sp%np
  nullify(remove)
  imported = 0
  !
  nb_link => self%patch%link%nbor
  do while (associated(nb_link))
    nbor => gpatch%cast2gpatch(nb_link%task)
    if (associated(nbor%connect%pic_patch)) then
      nbp => nbor%connect%pic_patch
      select type(nbp)
        class is (pic_t)
        nb_pic => nbp
      end select
      if (nb_pic%on) then
        do isp = 1, nspecies
          plist => nb_pic%exported(isp)
          do while(associated(plist))
            if (associated(plist%particle)) then
              ep => plist%particle
              if (ep%w /= 0) then
                do i=1,3
                  ppos(i) = nbor%position(i)+nbor%mesh(i)%r(ep%q(i)- mpi%offset(i))- &
                            nbor%ds(i)*0.5_8+ep%r(i)*nbor%ds(i)                 ! check the absolute particle position
                end do ! i=1,3
                if (self%patch%contains(ppos)) then                               ! this particle is inside the patch, thus has to be imported
                  !$omp atomic
                  sp(isp)%np = sp(isp)%np+1
                  sp(isp)%particle(sp(isp)%np) = ep
                  r = ppos - self%patch%position/self%patch%ds+self%patch%offset
                  cpos = nint(r)
                  r = r - real(cpos,kind=8) + 0.5_8*self%patch%ds
                  sp(isp)%particle(sp(isp)%np)%r = r
                  sp(isp)%particle(sp(isp)%np)%q = cpos+mpi%offset
                  imported = imported +1
                  nb_pic%exported(isp)%np = nb_pic%exported(isp)%np - 1
                  deallocate(ep)
                  ! remove the linked list entry
                  if (.not.associated(plist,nb_pic%exported(isp))) then
                    remove => plist
                    if(associated(remove%prev)) remove%prev%next => remove%next
                    if(associated(remove%next)) remove%next%prev => remove%prev
                    if (associated(remove,nb_pic%exported(isp)%tail)) &
                      nb_pic%exported(isp)%tail => remove%prev
                  end if
                else
                  !print *, 'not imported', self%patch%contains(ppos)
                  !print *, 'particle pos', ppos
                  !print *, 'llc, urc', self%patch%position - 0.5*self%patch%size, self%patch%position + 0.5*self%patch%size
                  !print *, 'contains', abs(ppos-self%patch%position) < 0.5_8*self%patch%size
                  !print *, 'r, q',ep%r, ep%q
                  !print *, 'particle weight', ep%w
                end if ! (self%patch%contains(ppos))
              end if
            end if
            plist => plist%next
            if(associated(remove)) deallocate(remove)
          end do ! while(associated(plist))
        end do ! isp = 1, nspecies
      end if ! (nb_pic%on)
    end if ! (associated(nbor%connect%pic_patch))
    nb_link => nb_link%next
  end do ! do while (associated(nb_link))
  print *, "Particles imported, total", imported, sp%np
END SUBROUTINE import_particles

!===============================================================================
!> Remove particles in 3-D region
!===============================================================================
SUBROUTINE remove_particles (self, l, u)
  class(pic_t):: self
  integer:: l(3), u(3)
  logical:: export
  ! ----------------------------------------------------------------------------
  logical:: periodic_dir
  integer:: q(3), ip, isp
  integer:: minq(3), maxq(3), minq0(3), maxq0(3), removed
  type(particle), pointer, dimension(:)   :: pa
  !--------------------------------------------------------------------------!
  print '(a,2(3i7,2x),3l3)', "Removal region:", l, u, periodic
  print*, "Total number of particles before remove", sp%np
  !
  removed = 0
  minq = 2**20
  maxq = -minq
  minq0 = 2**20
  maxq0 = -minq
  do isp = 1,nspecies
    pa => sp(isp)%particle
    do ip=1,sp(isp)%np
      q = pa(ip)%q - mpi%offset
      if (all(q >= l .and. q <= u)) then
        pa(ip)%w = 0
        removed = removed + 1
      else if (pa(ip)%w > 0.0) then
        minq = min(minq,q)
        maxq = max(maxq,q)
      else
        minq0 = min(minq0,q)
        maxq0 = max(maxq0,q)
      end if
    end do
  end do
  !if (io%debug(1)) then
    print*, "Particles removed, remaining ", removed, sp%np
  !end if
  !if (io%debug(2))
  print *,'minq,maxq', minq,maxq
END SUBROUTINE remove_particles


!===============================================================================
!> Remove particles in 3-D region
!===============================================================================
SUBROUTINE remove_particles2 (self, l, u)
  class(pic_t):: self
  integer:: l(3), u(3)
  logical:: export
  ! ----------------------------------------------------------------------------
  logical:: periodic_dir
  integer:: q(3), ip, isp
  integer:: minq(3), maxq(3), minq0(3), maxq0(3), removed
  real(8):: lp(3), up(3), pp(3)
  type(particle), pointer, dimension(:)   :: pa
  !--------------------------------------------------------------------------!
  print '(a,2(3i7,2x),3l3)', "Removal region:", l, u, periodic
  print*, "Total number of particles before remove", sp%np
  lp = l - 0.5_8*g%ds
  up = u + 0.5_8*g%ds
  removed = 0
  minq = 2**20
  maxq = -minq
  minq0 = 2**20
  maxq0 = -minq
  do isp = 1,nspecies
    pa => sp(isp)%particle
    do ip=1,sp(isp)%np
      pp = pa(ip)%q - mpi%offset + pa(ip)%r
      if (all(pp >= lp .and. pp <= up)) then
        pa(ip)%w = 0
        removed = removed + 1
      else if (pa(ip)%w > 0.0) then
        minq = min(minq,pa(ip)%q-mpi%offset)
        maxq = max(maxq,pa(ip)%q-mpi%offset)
      else
        minq0 = min(minq0,pa(ip)%q-mpi%offset)
        maxq0 = max(maxq0,pa(ip)%q-mpi%offset)
      end if
    end do
  end do
  !if (io%debug(1)) then
    print*, "Particles removed, remaining ", removed, sp%np
  !end if
  !if (io%debug(2))
  print *,'minq,maxq', minq,maxq
END SUBROUTINE remove_particles2

!===============================================================================
!> Remove particles in 3-D region
!===============================================================================
SUBROUTINE match_units(self)
  USE units_mod,       ONLY: cgs, si
  USE scaling_mod,     ONLY: scaling
  USE units,           ONLY: u,c,elm
  class(pic_t):: self
  real(8):: ds
  !--------------------------------------------------------------------------!
  ds = minval(self%patch%ds,self%patch%n > 1)
  ! SI
  !self%scaling%q_scale = 1d-8                                                ! electron charge scaling
  !self%scaling%e_mass_scale = 1836/18                                           ! electron mass  scaling
  !self%scaling%mu0_scale    = 316                                                            ! vacuum permeability
  !self%scaling%eps_scale   = 100                                                            ! vacuum permittivity
  ! Gaussian
  self%scaling%q_scale = 1d-16                                                ! electron charge scaling
  self%scaling%e_mass_scale = 1836/18                                           ! electron mass  scaling
  self%scaling%mu0_scale    = 1/cgs%mu_0!100000                                                            ! vacuum permeability
  self%scaling%eps_scale   = 1/cgs%eps_0!1000                                                            ! vacuum permittivity

  call self%scaling%init_pic(ds)
  call self%scaling%convert%init(self%scaling)
  c%c = 1!self%scaling%c
  c%c2= c%c**2
  c%oc        = 1./c%c
  c%oc2       = c%oc**2
  c%pi = cgs%pi
  !c%elm_q = self%scaling%q
  c%twopi = 2*c%pi
  c%sqrtpi = sqrt(c%pi)
  !c%fourpi = 4*c%pi
  !c%evK =
  !c%kb = self%scaling%k_b
  !c%mu = self%scaling%m_u
  !c%me = self%scaling%m_e
  c%gamma_gas = 1.666667
  c%mc = c%me*c%c
  c%mc2 = c%mc**2

  !elm%ke      = self%scaling%ke
  !elm%kb      = self%scaling%kb
  !elm%kf      = self%scaling%kf

    c%fourpi    = 1.
  !elm%ke      = 1./c%fourpi
  !elm%kb      = 1./(c%c*c%fourpi)
  !elm%kf      = 1./c%c
  elm%ke      = 1.
  elm%kb      = 1./c%c
  elm%kf      = 1./c%c
 
  c%c2        = c%c**2
  c%oc        = 1./c%c
  c%oc2       = c%oc**2

END SUBROUTINE match_units

END MODULE pic_mod
