!===============================================================================
!> Data type for dust-gas interaction:
!>
!> The dust evolution is simulated using a fairly large number of particles oer
!> cell, e.g. of order 10-30 per size bin, or of the order of 100 in total
!>
!> The particles are moved using a particle solver, which operates on lists of
!> particles, which have attributes suitable to describe dust properties.  The
!> most important of these are
!>
!>    particle_t%w    : a weight that expresses number density per particle
!>    particle_t%d    : the particle diameter (size)
!>    particle_t%m    : the particle_mass
!>    particle_t%q    : an index triple that giving the integer cell coordinates
!>    particle_t%r    : a fractional index giving the position in the cell
!>
!> The main procedures that contribute to evolving the dust are
!>
!>    dust_t%gather   : gathers distributions of particles inside each cell
!>    dust_t%evolve   : decides how to evolve the dust distribution
!>    dust_t%scatter  : updates the particle properties accordingly
!>
!> The dust evolution does not necessarily have to be done after each update of
!> particle positions; it could be done at regular intervals, semi-independently.
!>
!> The dust module could be called as a sub-process, or as a separate task.
!> Initially it is most likely the simplest and most transparent organization
!> to have the experiment_t level be the task, from the point of view of the
!> framework, and let it call the various sub-solvers at various points in time
!===============================================================================
MODULE dust_mod
  USE io_mod
  USE trace_mod
  USE particle_solver_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Dust data type
  !-----------------------------------------------------------------------------
  type, public:: dust_t
    real(8):: time = 0.0
    type(particle_solver_t):: particle_solver
  contains
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure, private:: gather
    procedure, private:: evolve
    procedure, private:: scatter
  end type
  integer:: verbose=0
  real:: min_size_cgs = 3e-4     ! 3 micron
  real:: max_size_cgs = 3000e-4  ! 3000 micron = 3 mm
  integer:: n_bin = 6            ! two bins per power of ten
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self)
  class(dust_t):: self
  !.............................................................................
  integer:: iostat
  logical:: first_time = .true.
  namelist /dust_parameters/ verbose, min_size_cgs, max_size_cgs, n_bin
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%update')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_parameters, iostat=iostat)
    call io%namelist_warning ('dust_parameters', iostat)
    write (io%output, dust_parameters)
    if (io%master) write (io%namelist, dust_parameters)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state
!===============================================================================
SUBROUTINE pre_update (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%pre_update')
  call self%gather
  call self%evolve
  call self%scatter
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Do whatever should be done after the HD/MHD update
!===============================================================================
SUBROUTINE post_update (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%post_update')
  !...
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Gather information about the current state of the dust and gas
!===============================================================================
SUBROUTINE gather (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%gather')
  call trace%end ()
END SUBROUTINE gather

!===============================================================================
!> Based on the current state of gas and dust, decide how to evolve it
!===============================================================================
SUBROUTINE evolve (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%evolve')
  call trace%end ()
END SUBROUTINE evolve

!===============================================================================
!> Update the dust state in the task, based on the desired evolution
!===============================================================================
SUBROUTINE scatter (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%scatter')
  call trace%end ()
END SUBROUTINE scatter

END MODULE dust_mod
