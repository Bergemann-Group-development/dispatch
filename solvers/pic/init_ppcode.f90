! main.f90 $Id$
! vim: nowrap
!-----------------------------------------------------------------------
SUBROUTINE init_ppcode
  USE params, ONLY: it, time, dt
  implicit none
  logical, external :: check_stop
!.......................................................................
  !call init_mpi                                                         ! MPI/OMP setup
  !call init_io                                                          ! input/output setup
  !call init_params                                                      ! simulation parameters
  !call init_units                                                       ! physical scaling units
  !call init_boundaries                                                  ! define boundaries
  !call init_switches                                                    ! execution control
  !call init_distributions                                               ! particle distibutions
  !call init_random                                                      ! random numbers
  !call init_pic_experiment                                              ! initialize the experiment
  !call init_maxwell                                                     ! initialize maxwell eq
  !call init_diagnostics                                                 ! energetics etc

  !call write_log                                                        ! initialize log output
  !call write_io                                                         ! write initial snapshot

END SUBROUTINE init_ppcode

! These calls can be triggered by touching the file reread.flag
!-----------------------------------------------------------------------
SUBROUTINE reread
implicit none
  call read_run_info                                                    ! new info params
  call read_params                                                      ! new simulation params
  call read_debug                                                       ! debuggin parameters
  call read_units                                                       ! adjusted scaling
  call read_species                                                     ! new species properties
  call read_grid                                                        ! new grid properties
  call read_switches                                                    ! new switch values
END SUBROUTINE reread

! This entry is called to make an orderly stop, also if errors occur
!-----------------------------------------------------------------------
SUBROUTINE stop_program (success)
  implicit none
  logical success
  call finalize_timer                                                   ! print final timers
  call finalize_compare                                                 ! final comparisons
  call finalize_mpi (success)                                           ! MPI close down
END SUBROUTINE stop_program
