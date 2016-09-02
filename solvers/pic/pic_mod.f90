!> =============================================================================
!> The data part of the pic_t derived type must contain all patch specific
!> dynamically  changing data for a patch.  This is on the one hand everything
!> that is permanently allocated via allocate commands in the PPcode, and on the
!> other hand scalars and arrays that change with time, or are location specific.
!> =============================================================================
MODULE pic_mod
  USE io_mod
  USE trace_mod
  USE mhd_mod
  USE grid_m, ONLY: g, mesh_type, grid_type
  USE grid_m, ONLY: Bx, By, Bz, Ex, Ey, Ez, fmhd
  USE params, ONLY: time, dt, it
  USE species_m, ONLY: spec, sp, field, fields, vthb !, counters, n
  USE maxwell_sources, ONLY: jx,jy,jz,sx,sy,sz
  implicit none
  private
  type:: scratch_t
    real, dimension(:,:,:), pointer:: jx, jy, jz
    real, dimension(:,:,:), pointer:: sx, sy, sz
  end type
  type, public, extends(mhd_t):: pic_t
    !type(counters), pointer:: n                         ! species.f90
    type(spec), pointer, dimension(:):: sp              ! species.f90
    type(field), pointer, dimension(:,:,:,:):: fields   ! species.f90
    real, pointer, dimension(:,:,:,:,:):: vthb          ! species.f90
    type(mesh_type), pointer, dimension(:,:,:):: em     ! grid.f90
    real, pointer, dimension(:,:,:):: Bx, By, Bz        ! grid.f90
    real, pointer, dimension(:,:,:):: Ex, Ey, Ez        ! grid.f90
    real, pointer, dimension(:,:,:):: fmhd              ! grid.f90
    type(grid_type), pointer:: g                        ! grid.f90
    type(scratch_t):: scratch                           ! maxwell_sources.f90
  contains
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: connect
    procedure:: copy
    procedure:: info
  end type
CONTAINS

!> =============================================================================
!> Allocate and initialize the memory and paramters belonging to a PIC patch.
!> This procedure is called once for every patch.  It calls init_species to
!> initialize everything that belongs to the patch being initialized, and the
!> dynamically changing part of this info needs to be copied to the patch
!> data type, to be recovered when a patch is to be updated.
!> =============================================================================
SUBROUTINE init (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! Initialize new PIC data for a patch
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::init')
  call init_grid                                                        ! computational mesh
  call init_species                                                     ! define particle species
  call init_distributions                                               ! particle distibutions
  call init_random                                                      ! random numbers
  call init_pic_experiment                                              ! initialize the experiment
  call init_maxwell                                                     ! initialize maxwell eq
  call init_diagnostics                                                 ! energetics etc
  !-----------------------------------------------------------------------------
  ! Save the location of dynamic PIC data for this patch
  !-----------------------------------------------------------------------------
  self%fields => fields                                 ! species.f90
  self%vthb   => vthb                                   ! species.f90
  self%sp     => sp                                     ! species.f90
  !self%n      => n                                      ! species.f90
  self%g    => g                                        ! grid.f90
  self%Bx   => Bx                                       ! grid.f90
  self%By   => By                                       ! grid.f90
  self%Bz   => Bz                                       ! grid.f90
  self%Ex   => Ex                                       ! grid.f90
  self%Ey   => Ey                                       ! grid.f90
  self%Ez   => Ez                                       ! grid.f90
  self%fmhd => fmhd                                     ! grid.f90
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
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
  !n      => self%n                                      ! species.f90
  g    => self%g                                        ! grid.f90
  Bx   => self%Bx                                       ! grid.f90
  By   => self%By                                       ! grid.f90
  Bz   => self%Bz                                       ! grid.f90
  Ex   => self%Ex                                       ! grid.f90
  Ey   => self%Ey                                       ! grid.f90
  Ez   => self%Ez                                       ! grid.f90
  fmhd => self%fmhd                                     ! grid.f90
  !-----------------------------------------------------------------------------
  ! A few scalar data items need to be copied, to avoid having to make pointers
  !-----------------------------------------------------------------------------
  time   = self%time
  dt     = self%dtime
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
  self%it     = it
  call trace_end
END SUBROUTINE copy

!> =============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!> =============================================================================
SUBROUTINE update (self)
  class(pic_t):: self
  associate (scratch=>self%scratch)
  !-----------------------------------------------------------------------------
  ! Prepare for update, by pointing to the patch data and allocating scratch
  !-----------------------------------------------------------------------------
  call trace_begin ('pic_mod::update')
  call self%connect
  allocate(scratch%jx(g%gn(1),g%gn(2),g%gn(3)), &
           scratch%jy(g%gn(1),g%gn(2),g%gn(3)), &
           scratch%jz(g%gn(1),g%gn(2),g%gn(3)))
  allocate(scratch%sx(g%gn(1),g%gn(2),g%gn(3)), &
           scratch%sy(g%gn(1),g%gn(2),g%gn(3)), &
           scratch%sz(g%gn(1),g%gn(2),g%gn(3)))
  jx => scratch%jx
  jy => scratch%jy
  jz => scratch%jz
  sx => scratch%sx
  sy => scratch%sy
  sz => scratch%sz
  !-----------------------------------------------------------------------------
  ! Do the PIC update, using calls to the PPcode routines
  !-----------------------------------------------------------------------------
  call set_dt                                           ! set timestep dt
  call move_particles                                   ! move particles
  call impose_particle_boundaries                       ! boundary conditions
  call send_particles                                   ! across MPI-bdries
  call sort                                             ! within MPI sub-domains
  call solve_maxwell                                    ! evolve B and E-fields
  it   = it + 1                                         ! time step counter
  time = time + dt                                      ! time advance
  call self%copy
  deallocate(scratch%jx,scratch%jy,scratch%jz)
  deallocate(scratch%sx,scratch%sy,scratch%sz)
  end associate
  call trace_end
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
  call self%connect
  call write_io                                         ! conditionally, do io
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
  call self%connect
  call write_log                                        ! write log to stdout
  call trace_end
END SUBROUTINE info

END MODULE pic_mod
