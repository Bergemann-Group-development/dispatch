!> =============================================================================
!> The data part of the pic_t derived type must contain all patch specific
!> dynamically  changing data for a patch.  This is on the one hand everything
!> that is permanently allocated via allocate commands in the PPcode, and on the
!> other hand scalars and arrays that change with time, or are location specific.
!> =============================================================================
MODULE pic_mod
  USE io_mod
  USE mhd_mod
  USE grid_m, only: g, mesh_type, grid_type
  USE grid_m, only: Bx, By, Bz, Ex, Ey, Ez, fmhd
  USE params, only: time, dt, it
  USE species, only: spec, sp, field, fields, vthb !, counters, n
  implicit none
  private
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
  call init_species
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
END SUBROUTINE connect

!> =============================================================================
!> Copy scalars back to the patch data type
!> =============================================================================
SUBROUTINE copy (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! A few scalar data items need to be copied, to avoid having to make pointers
  !-----------------------------------------------------------------------------
  self%time   = time
  self%dtime  = dt
  self%it     = it
END SUBROUTINE copy

!> =============================================================================
!> Update the state of the PIC patch -- take a timestep.  To be able to use the
!> PIC update procedures on the patch specific data, we need to first point the
!> dynamically changing part of this info to the PIC patch values.
!> =============================================================================
SUBROUTINE update (self)
  class(pic_t):: self
  !-----------------------------------------------------------------------------
  ! Prepare for update, by pointing to the patch data
  !-----------------------------------------------------------------------------
  call self%connect
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
END SUBROUTINE update

!> =============================================================================
!> Organize disk output.
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(pic_t):: self
  character(len=*), optional:: experiment_name
  !-----------------------------------------------------------------------------
  ! Connect to the patch data and do I/O, using calls to the PPcode routine
  !-----------------------------------------------------------------------------
  call self%connect
  call write_io                                         ! conditionally, do io
END SUBROUTINE output

!> =============================================================================
!> Print information to log output
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(pic_t):: self
  integer, optional:: nq, ntask
  character(len=*), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call self%connect
  call write_log                                        ! write log to stdout
END SUBROUTINE info

END MODULE pic_mod
