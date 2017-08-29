!> =============================================================================
!> The data part of the pic_t derived type must contain all patch specific
!> dynamically  changing data for a patch.  This is on the one hand everything
!> that is permanently allocated via allocate commands in the PPcode, and on the
!> other hand scalars and arrays that change with time, or are location specific.
!> =============================================================================
MODULE pic_mod
  USE io_mod,          ONLY: io
  USE trace_mod,       ONLY: trace_begin, trace_end
  USE mhd_mod,         ONLY: mhd_t
  USE grid_m,          ONLY: g, mesh_type, grid_type
  USE grid_m,          ONLY: Bx, By, Bz, Ex, Ey, Ez, fmhd
  USE params,          ONLY: time, dt, it, tend, master, mpi, mcoord, stdin, datadir
  USE species_m,       ONLY: spec, sp, field, fields, vthb, nspecies, particle !, counters, n
  USE maxwell_sources, ONLY: jx,jy,jz,sx,sy,sz
  USE random,          ONLY: wi,mt, mti
  USE rand_mod,        ONLY: rand_seed
  USE pic_io,          ONLY: toutf, toutp
  USE pic_stagger,     ONLY: xup, yup, zup
  implicit none
  private
  type:: scratch_t
    real, dimension(:,:,:), pointer:: jx, jy, jz
    real, dimension(:,:,:), pointer:: sx, sy, sz
  end type
    !integer,  parameter :: intg = selected_int_kind( 9 )
    !integer,  parameter :: wi = intg

  type, public, extends(mhd_t):: pic_t
    integer                                              :: i_ex=9, i_ey=10, i_ez=11
    !type(counters), pointer                             :: n            ! species.f90
    type(spec), pointer, dimension(:)                    :: sp           ! species.f90
    real, pointer, dimension(:,:,:,:,:)                  :: vthb         ! species.f90
    type(mesh_type), pointer, dimension(:,:,:)           :: em           ! grid.f90
    !real, pointer, dimension(:,:,:)                      :: Bx, By, Bz   ! grid.f90
    !real, pointer, dimension(:,:,:)                      :: Ex, Ey, Ez   ! grid.f90
    real, pointer, dimension(:,:,:)                      :: fmhd         ! grid.f90
    type(grid_type), pointer                             :: g=>null()    ! grid.f90
    type(scratch_t)                                      :: scratch      ! maxwell_sources.f90
    integer(kind = wi), pointer, dimension(:)            :: mti
    integer(kind = wi), pointer, dimension(:,:)          :: mt
    integer, pointer, dimension(:,:,:)                   :: rand_seed
    real(kind=8), pointer                                :: dt1
    real, pointer, dimension(:,:,:)                      :: temp
    logical                                              :: do_pairs
    !type(patch_t), pointer, dimension(:,:,:,:,:,:)      :: mem1         ! patch_mod.f90
  contains
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
    procedure, nopass:: print_patch_type
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
  integer:: isp, l(3), u(3), m(3)
  real:: msum, c
  !-----------------------------------------------------------------------------
  ! Initialize new PIC data for a patch
  !-----------------------------------------------------------------------------
  call trace_begin('pic_t%pre_init', itimer=itimer)
  if (firstcall) then
     call init_mpi                                                      ! MPI/OMP setup
     call init_io                                                       ! input/output setup
     call init_params                                                   ! simulation parameters
     call init_units                                                    ! physical scaling units
     datadir = io%outputname
  endif
  !-----------------------------------------------------------------------------
  ! Allocate PIC-related per-patch memory and copy static properties to PIC
  !-----------------------------------------------------------------------------
  m = self%mesh%gn
  print *,'nspecies',m,nspecies,self%nt
  Bx => self%mem(:,:,:,self%i_bx,self%it,1)                             ! pic_mod.f90
  By => self%mem(:,:,:,self%i_by,self%it,1)                             ! pic_mod.f90
  Bz => self%mem(:,:,:,self%i_bz,self%it,1)                             ! pic_mod.f90
  Ex => self%mem(:,:,:,self%i_ex,self%it,1)                             ! pic_mod.f90
  Ey => self%mem(:,:,:,self%i_ey,self%it,1)                             ! pic_mod.f90
  Ez => self%mem(:,:,:,self%i_ez,self%it,1)                             ! pic_mod.f90
  !-----------------------------------------------------------------------------
  ! Additional initialization of the PIC code.  The grid data type 'g' is
  ! allocated in init_boundaries, and we pick up the location here.
  !-----------------------------------------------------------------------------
  call init_boundaries                                                  ! define boundaries
  self%g => g
  g%s  = self%size
  g%ds = self%ds
  g%gn = self%mesh%n
  call init_grid                                                        ! computational mesh
  Ex=0; Ey=0; Ez=0
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -5
  call init_species                                                     ! define particle species
  if (firstcall) then
     call init_distributions                                            ! particle distibutions
  end if
  call init_random                                                      ! random numbers

  dt   = 0.0                                                            ! trigger dt=dtnew
  toutf = io%out_time
  toutp = io%out_time

  if (firstcall) then
     call init_switches                                                 ! execution control
     !--------------------------------------------------------------------------
     ! Sanity check first; use one cell less that interior, since the last
     ! point receives only half a contribution (FIXME)
     !--------------------------------------------------------------------------
     l = self%mesh%li+1
     u = self%mesh%ui-1

     if (io%debug(2)) then
       d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%i_d,self%it,1)
       e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%i_e,self%it,1)
       print*,'average and minimum mass density   before:', sum(d)/product(u-l+1), minval(d), self%it
       print*,'average and minimum thermal energy before:', sum(e)/product(u-l+1), minval(e), self%it
     end if

     call self%mhd_to_pic (self%mesh%lb, self%mesh%ub-1)                ! MHD->PIC in full patch
     !--------------------------------------------------------------------------
     ! Sort triggers a counting of particles per cell and index computations
     ! needed for mesh_sources to work, calling gather to compute vth
     !--------------------------------------------------------------------------
     call sort
     if (io%debug(2)) then
       print*,'number of particles expected', product(u-l+1)*sp(1)%ntarget
       print*,'number of particles after sort isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%n), &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%n)
       print*,'average density after sort isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%d)*sp(1)%mass/product(u-l+1)+ &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%d)*sp(2)%mass/product(u-l+1)
     end if

     call mesh_sources
     if (io%debug(2)) then
       do isp=1,nspecies
         print *,'average vth', isp, sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),isp)%vth)/product(u-l+1)
       end do
       print*,'average density after mesh_sources isp1, isp2:', &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),1)%d)*sp(1)%mass/product(u-l+1) + &
         sum(fields(l(1):u(1),l(2):u(2),l(3):u(3),2)%d)*sp(2)%mass/product(u-l+1)
     end if

     call self%pic_to_mhd
     if (io%debug(2)) then
       d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%i_d,self%new,1)
       e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%i_e,self%new,1)
       print*,'average and minimum mass density    after:', sum(d)/product(u-l+1), minval(d), self%it
       print*,'average and minimum thermal energy  after:', sum(e)/product(u-l+1), minval(e), self%it
     end if
     !--------------------------------------------------------------------------
     ! Reset the IC again
     !--------------------------------------------------------------------------
     firstcall = .false.
  endif
  !-----------------------------------------------------------------------------
  ! Initialize the PIC state from MHD for patches of PIC type
  !-----------------------------------------------------------------------------
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -4
  call self%mhd_to_pic (self%mesh%lb, self%mesh%ub)                  ! MHD->PIC in full patch
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -3
  Ex=0.0; Ey=0.0; Ez=0.0
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -2

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
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), -1
  call init_diagnostics                                 ! energetics etc
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2), 0
  call write_log                                        ! initialize log output
  call write_io                                         ! write initial snapshot
  !io%out_time = tend
  !-----------------------------------------------------------------------------
  ! Save the location of dynamic PIC data for this patch
  !-----------------------------------------------------------------------------
  Bx => self%mem(:,:,:,self%i_bx,self%it,1)           ! pic_mod.f90
  By => self%mem(:,:,:,self%i_by,self%it,1)             ! pic_mod.f90
  Bz => self%mem(:,:,:,self%i_bz,self%it,1)             ! pic_mod.f90
  Ex => self%mem(:,:,:,self%i_ex,self%it,1)             ! pic_mod.f90
  Ey => self%mem(:,:,:,self%i_ey,self%it,1)             ! pic_mod.f90
  Ez => self%mem(:,:,:,self%i_ez,self%it,1)             ! pic_mod.f90
  self%vthb      => vthb                                ! species.f90
  self%sp        => sp                                  ! species.f90
  self%g         => g                                   ! grid.f90
  self%fmhd      => fmhd                                ! grid.f90
  self%mti       => mti
  self%mt        => mt
  self%rand_seed => rand_seed
  self%dt1       => dt

  !self%mem1 => mem
  !self%time => time
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
  vthb   => self%vthb                                   ! species.f90
  sp     => self%sp                                     ! species.f90
  !n      => self%n                                     ! species.f90
  g    => self%g                                        ! grid.f90
  Bx   => self%mem(:,:,:,self%i_bx,self%it,1)           ! pic_mod.f90
  By   => self%mem(:,:,:,self%i_by,self%it,1)           ! pic_mod.f90
  Bz   => self%mem(:,:,:,self%i_bz,self%it,1)           ! pic_mod.f90
  Ex   => self%mem(:,:,:,self%i_ex,self%it,1)           ! pic_mod.f90
  Ey   => self%mem(:,:,:,self%i_ey,self%it,1)           ! pic_mod.f90
  Ez   => self%mem(:,:,:,self%i_ez,self%it,1)           ! pic_mod.f90
  fmhd => self%fmhd                                     ! grid.f90
  mti  => self%mti
  mt   => self%mt
  rand_seed => self%rand_seed
  dt => self%dt1
  !mem => self%mem1
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
  self%time   = time
  self%dtime  = dt
  !self%it     = it
  !self%dt1    = dt
  call trace_end
END SUBROUTINE copy

!> =============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!> =============================================================================
SUBROUTINE update (self)
  class(pic_t):: self
  integer:: iv
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Prepare for update, by pointing to the patch data and allocating scratch
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::update', itimer=itimer)
  call self%connect
  !
  associate (scratch=>self%scratch)
  allocate(scratch%jx(g%n(1),g%n(2),g%n(3)), &
           scratch%jy(g%n(1),g%n(2),g%n(3)), &
           scratch%jz(g%n(1),g%n(2),g%n(3)))
  allocate(scratch%sx(g%n(1),g%n(2),g%n(3)), &
           scratch%sy(g%n(1),g%n(2),g%n(3)), &
           scratch%sz(g%n(1),g%n(2),g%n(3)))
  jx => scratch%jx
  jy => scratch%jy
  jz => scratch%jz
  sx => scratch%sx
  sy => scratch%sy
  sz => scratch%sz
  end associate
  !-----------------------------------------------------------------------------
  ! Copy mem arrays to new time slot and update pointers
  !-----------------------------------------------------------------------------
  do iv=self%i_bx,self%i_ez
    self%mem(:,:,:,iv,self%new,1) = self%mem(:,:,:,iv,self%it,1)
  end do
  Bx => self%mem(:,:,:,self%i_bx,self%new,1)            ! pic_mod.f90
  By => self%mem(:,:,:,self%i_by,self%new,1)            ! pic_mod.f90
  Bz => self%mem(:,:,:,self%i_bz,self%new,1)            ! pic_mod.f90
  Ex => self%mem(:,:,:,self%i_ex,self%new,1)            ! pic_mod.f90
  Ey => self%mem(:,:,:,self%i_ey,self%new,1)            ! pic_mod.f90
  Ez => self%mem(:,:,:,self%i_ez,self%new,1)            ! pic_mod.f90
  !-----------------------------------------------------------------------------
  ! Do the PIC update, using calls to the PPcode routines
  !-----------------------------------------------------------------------------
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),1
  call self%impose_particle_boundaries_np               ! boundary conditions
  call set_dt                                           ! set timestep dt
  call move_particles                                   ! move particles
  call send_particles                                   ! across patch-bdries
  call self%impose_particle_boundaries_np               ! boundary conditions
  call sort                                             ! within patches
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),2
  print *,'g:',g%ds,g%n,g%lb,g%ub
  call solve_maxwell                                    ! evolve B and E-fields
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),3
  !it   = it + 1                                        ! time step counter
  time = time + dt                                      ! time advance
  call self%pic_to_mhd                                  ! compute MHD values
  call self%copy
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),4
  call deallocate_sources
  print *,'Bz', sum(Ex**2), sum(Ey**2), sum(Bz**2),5
  call trace_end (itimer=itimer)
END SUBROUTINE update

!> =============================================================================
!> Organize disk output.
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(pic_t):: self
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  ! Connect to the patch data and do I/O, using calls to the PPcode routine
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::output')
  if (associated(self%g)) then
    call self%connect
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
    call self%connect
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
  integer :: isp, ix, iy, iz, i
  real :: mean_mass
  real, dimension(:,:,:), pointer :: d, e, px, py, pz, bx, by, bz, jx, jy, jz
  !real, dimension(:,:,:), pointer :: temp
  integer, save:: itimer=0
  logical, save:: do_pairs, first_time=.true.
  namelist /init/ do_pairs, periodic_bcs
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%mhd_to_pic', itimer=itimer)

  !$omp critical (input_cr)
  if (first_time) then
     rewind(stdin); read (stdin,init)              ! rewind file descriptor and read input from file
     first_time = .false.
     write (*,init)
  end if
  !$omp end critical (input_cr)


  self%do_pairs = do_pairs
  d  => self%mem(:,:,:,self%i_d ,self%it,1)
  e  => self%mem(:,:,:,self%i_e ,self%it,1)
  px => self%mem(:,:,:,self%i_px,self%it,1)
  py => self%mem(:,:,:,self%i_py,self%it,1)
  pz => self%mem(:,:,:,self%i_pz,self%it,1)
  bx => self%mem(:,:,:,self%i_bx,self%it,1)
  by => self%mem(:,:,:,self%i_by,self%it,1)
  bz => self%mem(:,:,:,self%i_bz,self%it,1)

  !print*, "sp(1)%particle%w =", sp(1)%particle(1)%w
  msum = sum(sp(1:nspecies)%mass)              ! sum of different particle masses
  do isp=1, nspecies
    fields(:,:,:,isp)%d = d/msum            ! number density of each species (patch+gz)
    if (io%debug(2)) print*, 'msum, max(d)', msum, maxval(d)/msum
  end do


  !if (io%debug(2)) print*,'mhd_to_pic number of particles isp1, isp2:', sum(fields(:,:,:,1)%d), sum(fields(:,:,:,2)%d)
  !-----------------------------------------------------------------------------
  ! Calculation of electric current, which is needed to get the relative speeds
  ! of electrons and ions
  !-----------------------------------------------------------------------------
  allocate(jx(ub(1),ub(1),ub(1)), jy(ub(2),ub(2),ub(2)), jz(ub(3),ub(3),ub(3)))
  jx(lb(1),:,:)=0.0; jx(:,lb(2),:)=0.0; jx(:,:,lb(3))=0.0
  jy(lb(1),:,:)=0.0; jy(:,lb(2),:)=0.0; jy(:,:,lb(3))=0.0
  jz(lb(1),:,:)=0.0; jz(:,lb(2),:)=0.0; jz(:,:,lb(3))=0.0
  do iz=lb(3)+1, ub(3)
    do iy=lb(2)+1, ub(2)
      do ix=lb(1)+1, ub(1)
        jx(ix,iy,iz) = (ddydn(bz,ix,iy,iz) - ddzdn(by,ix,iy,iz))/c%mu0
        jy(ix,iy,iz) = (ddzdn(bx,ix,iy,iz) - ddxdn(bz,ix,iy,iz))/c%mu0
        jz(ix,iy,iz) = (ddxdn(by,ix,iy,iz) - ddydn(bx,ix,iy,iz))/c%mu0
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Prepare calls to make_particles_with_temperature_3d for each species, or
  ! the corresponding one for pairs.  Each particle is assigned a weight equal
  ! to fields(ix,iy,iz,isp)%d / sp(isp)%ntarget inside those subroutines.
  !-----------------------------------------------------------------------------
  allocate (self%temp(g%n(1),g%n(2),g%n(3)))
  mean_mass = msum/2.
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
    do iz=lb(3), ub(3)
      do iy=lb(2), ub(2)
        do ix=lb(1), ub(1)
          vex = weight*jx(ix,iy,iz)/sp(isp)%charge ! electric current vel
          vey = weight*jy(ix,iy,iz)/sp(isp)%charge
          vez = weight*jz(ix,iy,iz)/sp(isp)%charge
          fields(ix,iy,iz,isp)%v(1) = px(ix,iy,iz)/xup(d,ix,iy,iz) + vex  ! total vel
          fields(ix,iy,iz,isp)%v(2) = py(ix,iy,iz)/yup(d,ix,iy,iz) + vey
          fields(ix,iy,iz,isp)%v(3) = pz(ix,iy,iz)/zup(d,ix,iy,iz) + vez
          self%temp(ix,iy,iz) = mean_mass*e(ix,iy,iz)*(self%gamma-1.0)/(c%kb*d(ix,iy,iz))
        end do
      end do
    end do
    if (io%debug(2)) print*,'average temp:', sum(self%temp)/product(ub-lb+1)
    if (self%do_pairs) then
      if (mod(isp, 2) /= 0) then
        call make_particle_pairs_with_temperature_3d (isp, isp+1, self%temp)
      end if
    else
      call make_particles_with_temperature_3d (isp, self%temp, self%li, self%ui+1)
    end if
  end do
  deallocate (jx, jy, jz)
  call trace_end (itimer)
END SUBROUTINE mhd_to_pic

!> =============================================================================
!> Convert pic variables to mhd and store them in the self%new time slot
!> =============================================================================
SUBROUTINE pic_to_mhd (self)
  implicit none
  class(pic_t):: self
  integer:: isp, i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_t%pic_to_mhd', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Copy over fields (B is omitted since is using pointers to self%mem)
  !-----------------------------------------------------------------------------
  self%mem(:,:,:,self%i_d ,self%new,1) = 0.0
  self%mem(:,:,:,self%i_e ,self%new,1) = 0.0
  self%mem(:,:,:,self%i_px,self%new,1) = 0.0
  self%mem(:,:,:,self%i_py,self%new,1) = 0.0
  self%mem(:,:,:,self%i_pz,self%new,1) = 0.0

  do isp=1,nspecies
    self%mem(:,:,:,self%i_d ,self%new,1) = &
    self%mem(:,:,:,self%i_d ,self%new,1) + sp(isp)%mass*fields(:,:,:,isp)%d
    self%mem(:,:,:,self%i_e ,self%new,1) = &
    self%mem(:,:,:,self%i_e ,self%new,1) + sp(isp)%mass*fields(:,:,:,isp)%d* &
                                           0.5*fields(:,:,:,isp)%vth**2
    self%mem(:,:,:,self%i_px,self%new,1) = &
    self%mem(:,:,:,self%i_px,self%new,1) + sp(isp)%mass*fields(:,:,:,isp)%v(1)
    self%mem(:,:,:,self%i_py,self%new,1) = &
    self%mem(:,:,:,self%i_py,self%new,1) + sp(isp)%mass*fields(:,:,:,isp)%v(2)
    self%mem(:,:,:,self%i_pz,self%new,1) = &
    self%mem(:,:,:,self%i_pz,self%new,1) + sp(isp)%mass*fields(:,:,:,isp)%v(3)
  end do
  !-----------------------------------------------------------------------------
  call trace_end (itimer)
END SUBROUTINE pic_to_mhd
!=====================================================================!
! Non periodic boundaries
!=====================================================================!
SUBROUTINE impose_particle_boundaries_np(self)
  implicit none
  class(pic_t)                            :: self
  integer                                 :: ip,isp,idim,ii(3)
  integer                                 :: i,ix,iy,iz,lb(3),ub(3)
  real, dimension(mcoord)                 :: velb
  type(particle), pointer, dimension(:)   :: pa
  integer                                 :: below_lb, above_ub
  !----------------------------------------------------------------------!
  ! If alĺ boundaries are periodic we don't need guard zone particles
  if (all(periodic_bcs)) return
  call trace_begin ('pic_t%impose_particle_boundaries_np')
  !----------------------------------------------------------------------!
  ! First remove any particles outside the boundaries
  ! Replace particles in boundaries
  associate (m => self%mesh)
  call remove_boundary_particles
  !----------------------------------------------------------------------!
  ! Set limits depending on if a direction is periodic or not
  do i=1,3
    if (periodic_bcs(i)) then
      lb(i) = m(i)%li
      ub(i) = m(i)%ui
    else
      lb(i) = m(i)%lb
      ub(i) = m(i)%ub
    end if
  end do
  !----------------------------------------------------------------------!
  ! Add particles only in non-periodic directions
  if (.not.periodic_bcs(3)) then
    if (io%debug(2)) print 1, 'Add particles lower z:', [lb(1), lb(2), m(3)%lb], [ub(1), ub(2), m(3)%lo]
    1 format(1x,a,2(2x,3i4))
    call add_boundary_particles ([lb(1), lb(2), m(3)%lb], [ub(1), ub(2), m(3)%lo])
    if (io%debug(2)) print 1, 'Add particles upper z:', [lb(1), lb(2), m(3)%uo], [ub(1), ub(2), m(3)%ub]
    call add_boundary_particles ([lb(1), lb(2), m(3)%uo], [ub(1), ub(2), m(3)%ub])
    lb(3) = m(3)%li
    ub(3) = m(3)%ui
  end if
  if (.not.periodic_bcs(2)) then
    if (io%debug(2)) print 1, 'Add particles lower y:', [lb(1), m(2)%lb, lb(3)], [ub(1), m(2)%lo, ub(3)]
    call add_boundary_particles ([lb(1), m(2)%lb, lb(3)], [ub(1), m(2)%lo, ub(3)])
    if (io%debug(2)) print 1, 'Add particles upper y:', [lb(1), m(2)%uo, lb(3)], [ub(1), m(2)%ub, ub(3)]
    call add_boundary_particles ([lb(1), m(2)%uo, lb(3)], [ub(1), m(2)%ub, ub(3)])
    lb(2) = m(2)%li
    ub(2) = m(2)%ui
  end if
  if (.not.periodic_bcs(1)) then
    if (io%debug(2)) print 1, 'Add particles lower x:', [m(1)%lb, lb(2), lb(3)], [m(1)%lo, ub(2), ub(3)]
    call add_boundary_particles ([m(1)%lb, lb(2), lb(3)], [m(1)%lo, ub(2), ub(3)])
    if (io%debug(2)) print 1, 'Add particles upper x:', [m(1)%uo, lb(2), lb(3)], [m(1)%ub, ub(2), ub(3)]
    call add_boundary_particles ([m(1)%uo, lb(2), lb(3)], [m(1)%ub, ub(2), ub(3)])
  end if
  end associate
  call trace_end

  contains
  !----------------------------------------------------------------------------!
  ! Remove particles in guard zones
  SUBROUTINE remove_boundary_particles
    implicit none
    real:: vl(3), vu(3), vi(3)
    integer:: inside, offset
    !--------------------------------------------------------------------------!
    below_lb = 0; above_ub = 0     ! Removed particle counters
    vl=0.0; vu=0.0; vi=0.0; inside=0
    if (io%debug(1)) print*, "Total number of particles", sp%np
    do idim = 1,3
      if (periodic_bcs(idim)) cycle
      offset = mpi%offset(idim) + self%mesh(idim)%li
      if (io%debug(2)) print *,'offset',offset,mpi%offset(idim)
      do isp=1,nspecies
        pa => sp(isp)%particle
        do ip=1,sp(isp)%np
          if (pa(ip)%q(idim)-offset < 0) then
            if (pa(ip)%w > 0.0) then
              below_lb = below_lb + 1
              pa(ip)%w = 0.
              vl = vl+pa(ip)%p
            end if
          else if (pa(ip)%q(idim)-offset >= g%gn(idim)) then
            if (pa(ip)%w > 0.0) then
              above_ub = above_ub + 1
              pa(ip)%w = 0.
              vu = vu+pa(ip)%p
            end if
          else
            vi = vi + pa(ip)%p
            inside = inside + 1
          end if
        end do
      end do
    end do
    call sort
    if (io%debug(1)) then
      print*, "Particles removed below lower boundary", below_lb, vl/max(1,below_lb)
      print*, "Particles removed above upper boundary", above_ub, vu/max(1,above_ub)
      if (io%debug(2)) print*, "Particles inside both boundaries      ", inside  , vi/max(1,inside)
    end if
  END SUBROUTINE remove_boundary_particles
  !----------------------------------------------------------------------------!
  ! Add particles in region specified by (l, u)
  SUBROUTINE add_boundary_particles (l, u)
    USE particleoperations, ONLY: make_boundary_particle_pairs_with_temperature, &
                                  make_boundary_particles_with_temperature
    implicit none
    integer                   :: l(3), u(3)
    integer                   :: isp
    real                      :: c
    c = 1./sum(sp(1:nspecies)%mass)
    if (io%debug(2)) print *,'add_boundary_particles: minval(d) =',minval(self%mem(:,:,:,1,self%it,1))
    do isp=1,nspecies
      do iz=lb(3), ub(3)
        do iy=lb(2), ub(2)
          do ix=lb(1), ub(1)
            fields(ix,iy,iz,isp)%d = self%mem(ix,iy,iz,1,self%it,1)*c
          end do
        end do
      end do
      if (self%do_pairs) then
        if (mod(isp, 2) /= 0) then
          call make_boundary_particle_pairs_with_temperature (isp, isp+1, self%temp, l, u, io%debug(1))
        end if
      else
        call make_boundary_particles_with_temperature (isp, self%temp, l, u, io%debug(1))
      end if
    end do
  END SUBROUTINE add_boundary_particles
END SUBROUTINE impose_particle_boundaries_np

!> =============================================================================
SUBROUTINE print_patch_type (patch_type)
  implicit none
  logical:: patch_type
  print*, '==============================================================================='
  if (patch_type) then
    print*, '| Initializing PIC patch                                                      |'
  else
    print*, '| Initializing MHD patch                                                      |'
  end if
  print*, '==============================================================================='
END SUBROUTINE print_patch_type

END MODULE pic_mod
