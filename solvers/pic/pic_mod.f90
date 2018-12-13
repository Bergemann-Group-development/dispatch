!> =============================================================================
!> The data part of the pic_t derived type must contain all patch specific
!> dynamically  changing data for a patch.  This is on the one hand everything
!> that is permanently allocated via allocate commands in the PPcode, and on the
!> other hand scalars and arrays that change with time, or are location specific.
!> =============================================================================
MODULE pic_mod
  USE io_mod,          ONLY: io
  USE trace_mod,       ONLY: trace_begin, trace_end
  USE solver_mod,      ONLY: solver_t
  USE grid_m,          ONLY: g, mesh_type, grid_type
  USE grid_m,          ONLY: Bx, By, Bz, Ex, Ey, Ez, fmhd
  USE params,          ONLY: time, dt, it, tend, master, mpi, mcoord, stdin, &
                             nomp, romp, datadir, periodic
  USE species_m,       ONLY: spec, sp, field, fields, vthb, nspecies, particle !, counters, n
  USE maxwell_sources, ONLY: jx,jy,jz
  USE random,          ONLY: wi,mt, mti
  USE rand_mod,        ONLY: rand_seed
  USE pic_io,          ONLY: toutf, toutp
  USE pic_stagger,     ONLY: xup, yup, zup, xdn, ydn, zdn
  implicit none
  private

  type, public, extends(solver_t):: pic_t
    integer                                              :: i_ex=9, i_ey=10, i_ez=11
    type(field), pointer, dimension(:,:,:,:)             :: fields       ! species on mesh
    type(spec), pointer, dimension(:)                    :: sp           ! species.f90
    real, pointer, dimension(:,:,:,:,:)                  :: vthb         ! species.f90
    type(mesh_type), pointer, dimension(:,:,:)           :: em           ! grid.f90
    real, pointer, dimension(:,:,:)                      :: fmhd         ! grid.f90
    type(grid_type), pointer                             :: g=>null()    ! grid.f90
    integer(kind = wi), pointer, dimension(:)            :: mti
    integer(kind = wi), pointer, dimension(:,:)          :: mt
    integer, pointer, dimension(:,:,:)                   :: rand_seed
    real(kind=8), pointer                                :: dt1
    logical                                              :: do_pairs
  contains
    procedure :: print_patch_type
    procedure:: pre_init
    procedure:: post_init
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: connect
    procedure:: copy
    procedure:: info
    procedure:: mhd_to_pic
    procedure:: pic_to_mhd
    procedure:: impose_particle_boundaries_np
    procedure:: print_histogram
    procedure:: custom_refine
  end type
  logical, save :: firstcall = .true.
  logical, save :: periodic_bcs(3)=.false.
CONTAINS

!> =============================================================================
!> Initialize the PIC properties of a patch
!> =============================================================================
SUBROUTINE init (self)
 implicit none
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%init')
  !-----------------------------------------------------------------------------
  ! Make sure there is enough memory slots for the electro-magnetic field
  ! before calling the generic patch initialization
  !-----------------------------------------------------------------------------
  self%nv = 11
  call self%mhd_t%init
  self%idx%ex=9
  self%idx%ey=10
  self%idx%ez=11
  !-----------------------------------------------------------------------------
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Allocate and initialize the memory and paramters belonging to a PIC patch.
!> This procedure is called once for every patch.  It calls init_species to
!> initialize everything that belongs to the patch being initialized, and the
!> dynamically changing part of this info needs to be copied to the patch
!> data type, to be recovered when a patch is to be updated.
!> =============================================================================
SUBROUTINE pre_init (self)
 implicit none
  class(pic_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e
  integer:: isp, l(3), u(3), m(3), lb(3), ub(3)
  real:: msum, c
  logical, save:: do_pairs=.true., do_tests=.true.
  namelist /init/ do_pairs, periodic_bcs, do_tests
  !-----------------------------------------------------------------------------
  ! Initialize new PIC data for a patch
  !-----------------------------------------------------------------------------
  call trace_begin('pic_t%pre_init', itimer=itimer)
  !$omp critical (input_cr)
  if (firstcall) then
     stdin = io%input
     
     rewind(stdin); read (stdin,init)              ! rewind file descriptor and read input from file
     write (*,init)
     if (io%master) write (*, init)
     call init_mpi                                                      ! MPI/OMP setup
     nomp = 1
     romp = 1
     call init_io                                                       ! input/output setup
     call init_params                                                   ! simulation parameters
     call init_units                                                    ! physical scaling units
     datadir = io%outputname
     toutf = io%out_time
     toutp = io%out_time
     tend  = io%end_time
     mpi%lb = .true.
     mpi%ub = .true.
     mpi%interior = .false.
     periodic = periodic_bcs
  endif
  !$omp end critical (input_cr)
  self%do_pairs = do_pairs
  !-----------------------------------------------------------------------------
  ! Allocate PIC-related per-patch memory and copy static properties to PIC
  !-----------------------------------------------------------------------------
  Bx => self%mem(:,:,:,self%idx%bx,self%it,1)                             ! pic_mod.f90
  By => self%mem(:,:,:,self%idx%by,self%it,1)                             ! pic_mod.f90
  Bz => self%mem(:,:,:,self%idx%bz,self%it,1)                             ! pic_mod.f90
  Ex => self%mem(:,:,:,self%idx%ex,self%it,1)                             ! pic_mod.f90
  Ey => self%mem(:,:,:,self%idx%ey,self%it,1)                             ! pic_mod.f90
  Ez => self%mem(:,:,:,self%idx%ez,self%it,1)                             ! pic_mod.f90
  Ex=0.0; Ey=0.0; Ez=0.0
  !-----------------------------------------------------------------------------
  ! Additional initialization of the PIC code.  The grid data type 'g' is
  ! allocated in init_boundaries, and we pick up the location here.
  !-----------------------------------------------------------------------------
  call init_boundaries                                                  ! define boundaries
  self%g => g
  g%gs  = self%size
  g%ds = 0.0
  g%gn = self%mesh%n
  g%lb = self%mesh%ng                                                   ! first point inside
  g%ub = self%mesh%ng                                                   ! last point inside
  call init_grid                                                        ! computational mesh
  g%lb = self%mesh%li                                                   ! first point inside
  g%ub = self%mesh%uo                                                   ! last point inside
  mpi%offset = -g%lb
  call init_species                                                     ! define particle species
  call init_random                                                      ! random numbers
  if (firstcall) then
     call init_distributions                                            ! particle distibutions
     call init_switches                                                 ! execution control
  end if

  dt   = 0.0                                                            ! trigger dt=dtnew

  !-----------------------------------------------------------------------------
  ! Initialize the PIC state from MHD for patches of PIC type
  !-----------------------------------------------------------------------------
  lb = merge(self%mesh%li,self%mesh%li,periodic_bcs)
  ub = merge(self%mesh%ui,self%mesh%ui,periodic_bcs)
  call allocate_sources
  if (do_tests) then
     !call test_gather (self)
     dt = 0.5*minval(self%ds)
     !call test_scatter (self)
     !--------------------------------------------------------------------------
     ! Sanity check first; use one cell less that interior, since the last
     ! point receives only half a contribution
     !--------------------------------------------------------------------------
     l = self%mesh%li+1
     u = self%mesh%ui-1

     if (io%debug(2)) then
       d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1)
       e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%it,1)
       print*,'average and minimum mass density   before:', sum(d)/product(u-l+1), minval(d), self%it
       print*,'average and minimum thermal energy before:', sum(e)/product(u-l+1), minval(e), self%it
     end if

     call self%mhd_to_pic (lb, ub)
     call self%print_histogram ('after mhd_to_pic')
     !--------------------------------------------------------------------------
     ! Sort triggers a counting of particles per cell and index computations
     ! needed for mesh_sources to work, calling gather to compute vth
     !--------------------------------------------------------------------------
     call sort
     call self%print_histogram ('after sort')
     if (io%debug(2)) then
       print*,'number of particles expected', product(u-l+1)*sp(1)%ntarget
       print*,'number of particles after sort isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%n), &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%n)
       print*,'average density after sort isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%d)*sp(1)%mass/product(u-l+1)+ &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%d)*sp(2)%mass/product(u-l+1)
     end if

     call self%impose_particle_boundaries_np               ! boundary conditions
     call self%print_histogram ('after BCs')

     call sort
     call mesh_sources
     if (io%debug(2)) then
       do isp=1,nspecies
         print *,'average vth', isp, sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),isp)%vth)/product(u-l+1)
       end do
       print*,'average density after mesh_sources isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%d)*sp(1)%mass/product(u-l+1) + &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%d)*sp(2)%mass/product(u-l+1)
     end if

     call self%pic_to_mhd (self%it)
     if (io%debug(2)) then
       d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%new,1)
       e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%new,1)
       print*,'average and minimum mass density    after:', sum(d)/product(u-l+1), minval(d), self%it
       print*,'average and minimum thermal energy  after:', sum(e)/product(u-l+1), minval(e), self%it
     end if
     firstcall = .false.
  else
     call self%mhd_to_pic (lb, ub)
     call self%impose_particle_boundaries_np
  endif
  call deallocate_sources
  call trace_end (itimer=itimer)
END SUBROUTINE pre_init

!> =============================================================================
!> After allocating PIC-related arrays during initialization, save their
!> memory locations in the patch data type,  so they can be restored in the
!> connect procedure
!> =============================================================================
SUBROUTINE post_init (self)
  class(pic_t):: self
  call trace_begin ('pic_t%post_init')
  call init_maxwell                                     ! initialize maxwell eq
  if (io%debug(1)) &
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -1
  call init_diagnostics                                 ! energetics etc
  if (io%debug(1)) &
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), 0
  call write_log                                        ! initialize log output
  call write_io_pic (-1)                                ! write initial snapshot
  !self%out_time = tend
  !-----------------------------------------------------------------------------
  ! Save the location of dynamic PIC data for this patch
  !-----------------------------------------------------------------------------
  Bx => self%mem(:,:,:,self%idx%bx,self%it,1)           ! pic_mod.f90
  By => self%mem(:,:,:,self%idx%by,self%it,1)             ! pic_mod.f90
  Bz => self%mem(:,:,:,self%idx%bz,self%it,1)             ! pic_mod.f90
  Ex => self%mem(:,:,:,self%idx%ex,self%it,1)             ! pic_mod.f90
  Ey => self%mem(:,:,:,self%idx%ey,self%it,1)             ! pic_mod.f90
  Ez => self%mem(:,:,:,self%idx%ez,self%it,1)             ! pic_mod.f90
  self%fields    => fields
  self%vthb      => vthb                                ! species.f90
  self%sp        => sp                                  ! species.f90
  self%g         => g                                   ! grid.f90
  self%fmhd      => fmhd                                ! grid.f90
  self%mti       => mti
  self%mt        => mt
  self%rand_seed => rand_seed
  self%dt1       => dt
  call trace_end
END SUBROUTINE post_init

!> =============================================================================
!> Prepare the state of the PIC patch to take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!> =============================================================================
SUBROUTINE connect (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! Prepare for access, by pointing to the patch data
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::connect')
  fields => self%fields                                 ! species.f90
  vthb   => self%vthb                                   ! species.f90
  sp     => self%sp                                     ! species.f90
  !n      => self%n                                     ! species.f90
  g    => self%g                                        ! grid.f90
  Bx   => self%mem(:,:,:,self%idx%bx,self%it,1)           ! pic_mod.f90
  By   => self%mem(:,:,:,self%idx%by,self%it,1)           ! pic_mod.f90
  Bz   => self%mem(:,:,:,self%idx%bz,self%it,1)           ! pic_mod.f90
  Ex   => self%mem(:,:,:,self%idx%ex,self%it,1)           ! pic_mod.f90
  Ey   => self%mem(:,:,:,self%idx%ey,self%it,1)           ! pic_mod.f90
  Ez   => self%mem(:,:,:,self%idx%ez,self%it,1)           ! pic_mod.f90
  fmhd => self%fmhd                                     ! grid.f90
  mti  => self%mti
  mt   => self%mt
  rand_seed => self%rand_seed
  dt => self%dt1
  call allocate_sources
  !-----------------------------------------------------------------------------
  ! A few scalar data items need to be copied, to avoid having to make pointers
  !-----------------------------------------------------------------------------
  time   = self%time
  it     = self%it
  call trace_end
END SUBROUTINE connect

!> =============================================================================
!> Copy scalars back to the patch data type
!> =============================================================================
SUBROUTINE copy (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! A few scalar data items need to be copied, to avoid having to make pointers
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::copy')
  self%dtime  = dt
  call trace_end
END SUBROUTINE copy

!> =============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!> =============================================================================
SUBROUTINE update (self)
  class(pic_t):: self
  integer:: iv, ix, iy, iz
  integer, save:: itimer(2)=0
  !-----------------------------------------------------------------------------
  ! Prepare for update, by pointing to the patch data and allocating scratch
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::update', itimer=itimer(1))
  call self%connect
  !
  !-----------------------------------------------------------------------------
  ! Copy mem arrays to new time slot and update pointers
  !-----------------------------------------------------------------------------
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  if (periodic_bcs(1)) then
    do iv=1,self%nv
      call overlap_x(self%mem(:,:,:,iv,self%new,1))
    end do
  end if
  if (periodic_bcs(2)) then
    do iv=1,self%nv
      call overlap_y(self%mem(:,:,:,iv,self%new,1))
    end do
  end if
  if (periodic_bcs(3)) then
    do iv=1,self%nv
      call overlap_z(self%mem(:,:,:,iv,self%new,1))
    end do
  end if
  Bx => self%mem(:,:,:,self%idx%bx,self%new,1)            ! pic_mod.f90
  By => self%mem(:,:,:,self%idx%by,self%new,1)            ! pic_mod.f90
  Bz => self%mem(:,:,:,self%idx%bz,self%new,1)            ! pic_mod.f90
  Ex => self%mem(:,:,:,self%idx%ex,self%new,1)            ! pic_mod.f90
  Ey => self%mem(:,:,:,self%idx%ey,self%new,1)            ! pic_mod.f90
  Ez => self%mem(:,:,:,self%idx%ez,self%new,1)            ! pic_mod.f90
  !-----------------------------------------------------------------------------
  ! Do the PIC update, using calls to the PPcode routines
  !-----------------------------------------------------------------------------
  if (io%debug(1)) &
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),1
  call set_dt                                           ! set timestep dt
  call move_particles                                   ! move particles
  call send_particles                                   ! across patch-bdries
  call self%print_histogram ('after move')
  call self%impose_particle_boundaries_np               ! boundary conditions
  call self%print_histogram ('after BCs')
  call sort                                             ! within patches
  call trace_end (itimer(1))
  call solve_maxwell                                    ! evolve B and E-fields
  call trace_begin ('pic_t%solve_maxwell', itimer=itimer(2))
  time = time + dt                                      ! time advance
  call trace_end (itimer(2))
  call self%pic_to_mhd (self%new)                       ! compute MHD values
  if (io%debug(2)) then
    associate (l => self%mesh%li+1, u => self%mesh%ui-1)
    block
    real, dimension(:,:,:), pointer:: d, e
    d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%new,1)
    e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,self%new,1)
    print*,'average and minimum mass density   after pic_to_mhd:', sum(d)/product(u-l+1), minval(d), self%it
    print*,'average and minimum thermal energy after pic_to_mhd:', sum(e)/product(u-l+1), minval(e), self%it
    end block
    end associate
  end if
  call self%copy
  call deallocate_sources
END SUBROUTINE update

!> =============================================================================
!> Organize disk output.
!> =============================================================================
SUBROUTINE output (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! Connect to the patch data and do I/O, using calls to the PPcode routine
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::output')
  if (associated(self%g)) then
    call write_io_pic (self%iout)                             ! conditionally, do io
  end if
  call self%mhd_t%output
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> Print information to log output
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(pic_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::info')
  if (associated(self%g)) then
    call write_log                                        ! write log to stdout
  else
    call self%mhd_t%patch_t%info (nq, ntask)
  end if
  call trace_end
END SUBROUTINE info

!> =============================================================================
!> Convert mhd variables to pic
!> =============================================================================
SUBROUTINE mhd_to_pic (self, lb, ub)
  USE pic_stagger,        ONLY: ddxdn, ddydn, ddzdn, ddxdna, ddydna, ddzdna
  USE units,              ONLY: c
  USE particleoperations, ONLY: make_particles_with_temperature_3d, &
                                make_particle_pairs_with_temperature_3d
  implicit none
  class(pic_t):: self
  integer, dimension(3) :: lb, ub
  real :: msum, vex, vey, vez, weight, vg(3)
  integer :: isp, ix, iy, iz, i, l(3), u(3)
  real :: mean_mass, jmax
  real :: vemax, xupd, yupd, zupd
  real, dimension(:,:,:), pointer :: d, e, px, py, pz
  real, dimension(:,:,:), pointer :: temp
  integer, save:: itimer=0
  real:: ux, uy, uz
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%mhd_to_pic', itimer=itimer)
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  msum = sum(sp(1:nspecies)%mass)              ! sum of different particle masses
  do isp=1, nspecies
    fields(:,:,:,isp)%d = d/msum
    if (io%debug(2)) print*, 'msum, max(d)', msum, maxval(d)/msum
  end do
  !-----------------------------------------------------------------------------
  ! Calculation of electric current, which is needed to get the relative speeds
  ! of electrons and ions, and the particle temperature
  !-----------------------------------------------------------------------------
  allocate (temp(g%n(1),g%n(2),g%n(3)))
  mean_mass = msum/2.
  jmax = 0.0
  
  do iz=lb(3), ub(3)
    do iy=lb(2), ub(2)
      do ix=lb(1), ub(1)
        jx(ix,iy,iz) = (ddydn(bz,ix,max(iy,2),iz) - ddzdn(by,ix,iy,max(iz,2)))/c%mu0
        jy(ix,iy,iz) = (ddzdn(bx,ix,iy,max(iz,2)) - ddxdn(bz,max(ix,2),iy,iz))/c%mu0
        jz(ix,iy,iz) = (ddxdn(by,max(ix,2),iy,iz) - ddydn(bx,ix,max(iy,2),iz))/c%mu0
        jmax = max(jmax,jx(ix,iy,iz)**2+jy(ix,iy,iz)**2+jz(ix,iy,iz)**2)
        temp(ix,iy,iz) = mean_mass*e(ix,iy,iz)*(self%gamma-1.0)/(c%kb*d(ix,iy,iz))
      end do
    end do
  end do
  if (io%debug(1)) print *, 'jmax =', sqrt(jmax)
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
    vemax = 0.0
    do iz=lb(3), ub(3)
      do iy=lb(2), ub(2)
        do ix=lb(1), ub(1)
          xupd = xup(d,min(ix,self%mesh(1)%ub-1),iy,iz)
          yupd = yup(d,ix,min(iy,self%mesh(2)%ub-1),iz)
          zupd = zup(d,ix,iy,min(iz,self%mesh(3)%ub-1))
          vex = weight*jx(ix,iy,iz)*msum/(xupd*sp(isp)%charge) ! electric current vel
          vey = weight*jy(ix,iy,iz)*msum/(yupd*sp(isp)%charge)
          vez = weight*jz(ix,iy,iz)*msum/(zupd*sp(isp)%charge)
          vemax = max(vemax,vex**2+vey**2+vez**2)
          ux = px(ix,iy,iz)/xupd
          uy = py(ix,iy,iz)/yupd
          uz = pz(ix,iy,iz)/zupd
          ! --- total velocity ---
          fields(ix,iy,iz,isp)%v(1) = ux + vex
          fields(ix,iy,iz,isp)%v(2) = uy + vey
          fields(ix,iy,iz,isp)%v(3) = uz + vez
          ! --- electric field - U x B ---
          ex(ix,iy,iz) = - uy*ydn(bz,ix,max(iy,2),iz) + uz*zdn(by,ix,iy,max(iz,2))
          ey(ix,iy,iz) = - uz*zdn(bx,ix,iy,max(iz,2)) + ux*xdn(bz,max(ix,2),iy,iz)
          ez(ix,iy,iz) = - ux*xdn(by,max(ix,2),iy,iz) + uy*ydn(bx,ix,max(iy,2),iz)
        end do
      end do
    end do
    if (io%debug(2)) print *, 'max current speed', isp, sqrt(vemax)
  end do
  !----------------------------------------------------------------------------!
  ! Since make_particle_pairs needs velocitites for both species, the calls 
  ! below must be in a separate loop over species
  !----------------------------------------------------------------------------!
  do isp=1, nspecies
    if (self%do_pairs) then
      if (mod(isp, 2) /= 0) then
        call make_particle_pairs_with_temperature_3d (isp, isp+1, temp, lb, ub+1)
      end if
    else
      call make_particles_with_temperature_3d (isp, temp, lb, ub+1)
    end if
  end do
  deallocate (temp)
  call trace_end (itimer)
END SUBROUTINE mhd_to_pic

!> =============================================================================
!> Convert mhd variables to pic
!> =============================================================================
SUBROUTINE test_gather (self)
  USE particleoperations, ONLY: make_particles_with_temperature_3d
  implicit none
  class(pic_t):: self
  integer, dimension(3) :: lb, ub
  integer :: isp, ix, iy, iz, m(3)
  real, dimension(:,:,:), pointer :: temp
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%test_gather')
  !-----------------------------------------------------------------------------
  allocate (temp(g%n(1),g%n(2),g%n(3)))
  isp=1
  lb = self%mesh%lb
  ub = self%mesh%ub
  m = (lb+ub)/2
  do iz=lb(3), ub(3)
    do iy=lb(2), ub(2)
      do ix=lb(1), ub(1)
        if (ix==m(1).and.iy==m(2).and.iz==m(3)) then
          fields(ix,iy,iz,isp)%v = [1.,2.,3.]
          fields(ix,iy,iz,isp)%d = 1.0
          ex(ix,iy,iz) = 0.0
          ey(ix,iy,iz) = 0.0
          ez(ix,iy,iz) = 0.0
        else
          fields(ix,iy,iz,isp)%v = 0.0
          fields(ix,iy,iz,isp)%d = 0.0
          ex(ix,iy,iz) = 0.0
          ey(ix,iy,iz) = 0.0
          ez(ix,iy,iz) = 0.0
        end if
        temp(ix,iy,iz) = 1e-30
      end do
    end do
  end do
  call make_particles_with_temperature_3d (isp, temp, lb, ub+1)
  deallocate (temp)
  !--------------------------------------------------------------------------!
  ! Now check what we got back -- should be smooth charge and current
  !--------------------------------------------------------------------------!
  call sort
  call mesh_sources
  write (90) fields(:,:,:,isp)%d
  write (90) fields(:,:,:,isp)%v(1)
  write (90) fields(:,:,:,isp)%v(2)
  write (90) fields(:,:,:,isp)%v(3)
  flush (90)
  call remove_particles (self, lb, ub)
  !--------------------------------------------------------------------------!
  call trace_end
END SUBROUTINE test_gather

!> =============================================================================
!> Convert mhd variables to pic
!> =============================================================================
SUBROUTINE test_scatter (self)
  USE particleoperations, ONLY: make_particles_with_temperature_3d
  implicit none
  class(pic_t):: self
  integer, dimension(3) :: lb, ub
  integer :: isp, ix, iy, iz, m(3)
  real, dimension(:,:,:), pointer :: temp
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%test_scatter')
  !-----------------------------------------------------------------------------
  allocate (temp(g%n(1),g%n(2),g%n(3)))
  isp=1
  lb = self%mesh%lb
  ub = self%mesh%ub
  m = (lb+ub)/2
  do iz=lb(3), ub(3)
    do iy=lb(2), ub(2)
      do ix=lb(1), ub(1)
        fields(ix,iy,iz,isp)%d = 1.0
        fields(ix,iy,iz,isp)%v = 0.0
        temp(ix,iy,iz) = 1e-30
        ex(ix,iy,iz) = 0.0
        ey(ix,iy,iz) = 0.0
        ez(ix,iy,iz) = 0.0
        bx(ix,iy,iz) = 0.0
        by(ix,iy,iz) = 0.0
        bz(ix,iy,iz) = 0.0
        if (ix==m(1).and.iy==m(2).and.iz==m(3)) then
          ex(ix,iy,iz) = sp(isp)%mass/sp(isp)%charge
        else  if (ix==m(1)-10.and.iy==m(2).and.iz==m(3)) then
          ey(ix,iy,iz) = sp(isp)%mass/sp(isp)%charge
        else  if (ix==m(1)+10.and.iy==m(2).and.iz==m(3)) then
          ez(ix,iy,iz) = sp(isp)%mass/sp(isp)%charge
        end if
      end do
    end do
  end do
  call make_particles_with_temperature_3d (isp, temp, lb, ub+1)
  deallocate (temp)
  !--------------------------------------------------------------------------!
  ! Now check what we got back -- should be smooth charge and current
  !--------------------------------------------------------------------------!
  call sort
  call move_particles
  call sort
  call mesh_sources
  write (91) ex
  write (91) ey
  write (91) ez
  write (91) fields(:,:,:,isp)%v(1)
  write (91) fields(:,:,:,isp)%v(2)
  write (91) fields(:,:,:,isp)%v(3)
  flush (91)
  call remove_particles (self, lb, ub)
  !stop
  !--------------------------------------------------------------------------!
  call trace_end
END SUBROUTINE test_scatter

!----------------------------------------------------------------------------!
! Remove particles in 3-D region
SUBROUTINE remove_particles (self, l, u)
  implicit none
  class(pic_t):: self
  logical:: periodic_dir
  integer:: l(3), u(3), q(3), ip, isp
  integer:: minq(3), maxq(3), minq0(3), maxq0(3), removed
  type(particle), pointer, dimension(:)   :: pa
  !--------------------------------------------------------------------------!
  if (io%debug(1)) print '(a,2(3i7,2x),3l3)', "Removal region:", l, u, periodic
  if (io%debug(2)) print*, "Total number of particles before remove", sp%np
  removed = 0
  minq = 2**20
  maxq = -minq
  minq0 = 2**20
  maxq0 = -minq
  do isp = 1,nspecies
    pa => sp(isp)%particle
    do ip=1,sp(isp)%np
      q = pa(ip)%q - mpi%offset
      where (periodic) q = modulo(q-g%lb,g%gn) + g%lb
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
  if (io%debug(1)) then
    print*, "Particles removed, remaining ", removed, sp%np
  end if
  if (io%debug(2)) print '(a,4(3i4,2x))', 'minq,maxq', minq,maxq
END SUBROUTINE remove_particles

!> =============================================================================
!> Convert pic variables to mhd and store them in the self%new time slot
!> =============================================================================
SUBROUTINE pic_to_mhd (self, it)
  implicit none
  integer:: it
  class(pic_t):: self
  integer:: isp, i, ix, iy, iz
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%pic_to_mhd', itimer=itimer)
  if (io%debug(2)) &
  print *,'PIC_TO_MHD', self%mesh%li, self%mesh%ui
  !-----------------------------------------------------------------------------
  ! Copy over fields (B is omitted since is using pointers to self%mem)
  !-----------------------------------------------------------------------------
  associate (m => self%mesh)
  do iz=m(3)%lb,m(3)%ub
  do iy=m(2)%lb,m(2)%ub
  do ix=m(1)%lb,m(1)%ub
    self%mem(ix,iy,iz,self%idx%d ,it,1) = 0.0
    self%mem(ix,iy,iz,self%idx%e ,it,1) = 0.0
    self%mem(ix,iy,iz,self%idx%px,it,1) = 0.0
    self%mem(ix,iy,iz,self%idx%py,it,1) = 0.0
    self%mem(ix,iy,iz,self%idx%pz,it,1) = 0.0
  end do
  end do
  end do
  do isp=1,nspecies
    do iz=m(3)%li-1,m(3)%ui+2
    do iy=m(2)%li-1,m(2)%ui+2
    do ix=m(1)%li-1,m(1)%ui+2
      self%mem(ix,iy,iz,self%idx%d ,it,1) = &
      self%mem(ix,iy,iz,self%idx%d ,it,1) + sp(isp)%mass*fields(ix,iy,iz,isp)%d
      self%mem(ix,iy,iz,self%idx%e ,it,1) = &
      self%mem(ix,iy,iz,self%idx%e ,it,1) + sp(isp)%mass*fields(ix,iy,iz,isp)%d* &
                                                         0.5*fields(ix,iy,iz,isp)%vth**2
      self%mem(ix,iy,iz,self%idx%px,it,1) = &
      self%mem(ix,iy,iz,self%idx%px,it,1) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(1)
      self%mem(ix,iy,iz,self%idx%py,it,1) = &
      self%mem(ix,iy,iz,self%idx%py,it,1) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(2)
      self%mem(ix,iy,iz,self%idx%pz,it,1) = &
      self%mem(ix,iy,iz,self%idx%pz,it,1) + sp(isp)%mass*fields(ix,iy,iz,isp)%v(3)
  end do
  end do
  end do
  end do
  do i=1,5
    call overlap_add (self%mem(:,:,:,i,it,1))
  end do
  end associate
  !-----------------------------------------------------------------------------
  call trace_end (itimer)
END SUBROUTINE pic_to_mhd

!=====================================================================!
! Non periodic boundaries
!=====================================================================!
SUBROUTINE print_histogram (self, label)
  implicit none
  class(pic_t)                            :: self
  character(len=*)                        :: label
  integer                                 :: ip, isp, iz, l(2), u(2)
  type(particle), pointer                 :: pa
  integer, dimension(:,:), pointer        :: n
  real, dimension(:,:,:), pointer         :: v
  real                                    :: d(2)
  !-------------------------------------------------------------------!
  if (.not.io%debug(2)) return
  allocate (n(0:self%gn(3)+1,nspecies))
  allocate (v(0:self%gn(3)+1,3,nspecies))
  n = 0
  v = 0.0
  do isp=1,nspecies
    do ip=1,sp(isp)%np
      pa => sp(isp)%particle(ip)
      if (pa%w==0.0) cycle
      iz = pa%q(3) - mpi%offset(3)
      iz = max(0,min(iz,self%gn(3)+1))
      n(iz,isp) = n(iz,isp) + 1
      v(iz,:,isp) = v(iz,:,isp) + pa%p
    end do
  end do
  print *, 'hist: '//trim(label)
  l = self%mesh(1:2)%li
  u = self%mesh(1:2)%ui
  do iz=0,self%gn(3)+1
    do isp=1,nspecies
      if (n(iz,isp)>0) v(iz,:,isp) = v(iz,:,isp)/n(iz,isp)
      if (iz>0.and.iz<=self%gn(3)) then
        if (isp==1) then
          d(isp) = sum(bz(l(1):u(1),l(2):u(2),iz)**2)/product(u-l+1)
        else
          d(isp) = sum(ex(l(1):u(1),l(2):u(2),iz)**2)/product(u-l+1)
        end if
      else
        d(isp) = 0.0
      end if
    end do
    print '(i4,2i8,2x,1p,2e12.3,0p,2(2x,1p,3f9.3),3x,a)', &
      iz, n(iz,:), d, v(iz,:,:), ' hist: '//trim(label)
  end do
  deallocate (n)
END SUBROUTINE print_histogram

!=====================================================================!
! Non periodic boundaries
!=====================================================================!
SUBROUTINE impose_particle_boundaries_np(self)
  implicit none
  class(pic_t)                            :: self
  integer                                 :: ip,isp,idim,ii(3)
  integer                                 :: i,ix,iy,iz,l(3),u(3)
  integer                                 :: lr(3),ur(3),la(3),ua(3),lo(3),uo(3)
  real, dimension(mcoord)                 :: velb
  type(particle), pointer, dimension(:)   :: pa
  integer                                 :: below_lb, above_ub
  integer, save:: itimer=0
  !----------------------------------------------------------------------!
  ! If alĺ boundaries are periodic we don't need guard zone particles
  if (all(periodic_bcs)) return
  call trace_begin ('pic_t%impose_particle_boundaries_np',itimer=itimer)
  !----------------------------------------------------------------------!
  ! Only in non-periodic directions:
  !   First remove any particles in guard zones
  !   Then add new particles, avoiding the last layers
  !----------------------------------------------------------------------!
  associate (m=>self%mesh,lo=>self%mesh%lo,uo=>self%mesh%uo)
  lr = merge(m%lb  ,m%lb  ,periodic)
  ur = merge(m%ub  ,m%ub  ,periodic)
  la = merge(m%li  ,m%lb+1,periodic)
  ua = merge(m%ui  ,m%ub-1,periodic)
  !do i=1,3
  !  lo(i) = m(i)%lo
  !  uo(i) = m(i)%uo
  !  lr(i) = m(i)%lb
  !  ur(i) = m(i)%ub
  !  if (periodic(i)) then
  !    la(i) = m(i)%li
  !    ua(i) = m(i)%ui
  !  else
  !    la(i) = m(i)%li+1
  !    ua(i) = m(i)%ub-1
  !  end if
  !end do
  if (.not.periodic_bcs(1)) then
    l = [lr(1), lr(2), lr(3)]
    u = [lo(1), ur(2), ur(3)]
    call remove_particles (self, l, u)
    l = [la(1), la(2), la(3)]
    u = [lo(1), ua(2), ua(3)]
    if (io%debug(2)) print 1, 'Add particles lower x:', l , u
    1 format(a,2(2x,3i4))
    call self%mhd_to_pic (l, u)
    !....................................................................!
    l = [uo(1), lr(2), lr(3)]
    u = [ur(1), ur(2), ur(3)]
    call remove_particles (self, l, u)
    l = [uo(1), la(2), la(3)]
    u = [ua(1), ua(2), ua(3)]
    if (io%debug(2)) print 1, 'Add particles upper x:', l, u
    call self%mhd_to_pic (l, u)
    lr(1) = m(1)%li       ! shrink removal boundary
    ur(1) = m(1)%ui       ! shrink removal boundary
    la(1) = m(1)%li       ! shrink addition boundary
    ua(1) = m(1)%ui       ! shrink addition boundary
  end if
  !----------------------------------------------------------------------!
  if (.not.periodic_bcs(2)) then
    l = [lr(1), lr(2), lr(3)]
    u = [ur(1), lo(2), ur(3)]
    call remove_particles (self, l, u)
    l = [la(1), la(2), la(3)]
    u = [ua(1), lo(2), ua(3)]
    if (io%debug(2)) print 1, 'Add particles lower y:', l , u
    call self%mhd_to_pic (l, u)
    !....................................................................!
    l = [lr(1), uo(2), lr(3)]
    u = [ur(1), ur(2), ur(3)]
    call remove_particles (self, l, u)
    l = [la(1), uo(2), la(3)]
    u = [ua(1), ua(2), ua(3)]
    if (io%debug(2)) print 1, 'Add particles upper y:', l, u
    call self%mhd_to_pic (l, u)
    lr(2) = m(2)%li       ! shrink removal boundary
    ur(2) = m(2)%ui       ! shrink removal boundary
    la(2) = m(2)%li       ! shrink addition boundary
    ua(2) = m(2)%ui       ! shrink addition boundary
  end if
  !----------------------------------------------------------------------!
  if (.not.periodic_bcs(3)) then
    l = [lr(1), lr(2), lr(3)]
    u = [ur(1), ur(2), lo(3)]
    call remove_particles (self, l, u)
    l = [la(1), la(2), la(3)]
    u = [ua(1), ua(2), lo(3)]
    if (io%debug(2)) print 1, 'Add particles lower z:', l , u
    call self%mhd_to_pic (l, u)
    !....................................................................!
    l = [lr(1), lr(2), uo(3)]
    u = [ur(1), ur(2), ur(3)]
    call remove_particles (self, l, u)
    l = [la(1), la(2), uo(3)]
    u = [ua(1), ua(2), ua(3)]
    if (io%debug(2)) print 1, 'Add particles upper z:', l, u
    call self%mhd_to_pic (l, u)
  end if
  end associate
  call trace_end (itimer)
END SUBROUTINE impose_particle_boundaries_np

!> =============================================================================
SUBROUTINE print_patch_type (self, patch_type)
  class(pic_t) :: self
  logical :: patch_type
  print*, '==============================================================================='
  if (patch_type) then
    print*, '| Initializing PIC patch                                                      |'
  else
    print*, '| Initializing MHD patch                                                      |'
  end if
  print*, '==============================================================================='
END SUBROUTINE print_patch_type

!> =============================================================================
SUBROUTINE custom_refine (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  if (self%initial%first_pic) then
      self%initial%pic_needed = .true.
      call self%print_patch_type(self%initial%pic_needed)
      call self%patch_t%init
      call self%pre_init
      call self%post_init
      self%initial%first_pic = .false.
  end if
END SUBROUTINE custom_refine

END MODULE pic_mod
