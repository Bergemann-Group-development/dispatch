!===============================================================================
!> Port and test dust evolution from Podolak_1.py
!===============================================================================
MODULE dust_evolve_mod
  USE omp_lib
  USE units_mod
  USE math_mod
  USE trace_mod
  USE io_mod
  USE scaling_mod
  USE particles_mod
  implicit none
  private
  type, public:: dust_evolve_t
    logical:: on=.false.
    integer:: n_bin
  contains
    procedure:: init
    procedure:: update
    procedure:: test
  end type
  integer:: verbose=0
  logical:: on=.false.
  type(dust_evolve_t), public:: dust_evolve
CONTAINS

!===============================================================================
!> Read parameters and initialize -- add more parameters, if they are known
!> in the particle solver
!===============================================================================
SUBROUTINE init (self, n_bin, m_cgs)
  class(dust_evolve_t) :: self
  integer              :: n_bin               ! number of bins
  real                 :: m_cgs(:)            ! particle masses
  !.............................................................................
  integer              :: iostat              ! missing(<0), syntax errors (>0)
  logical, save        :: first_time=.true.   ! shared memory
  !--- add parameters here -----------------------------------------------------
  namelist /dust_evolve_params/ on, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%init')
  !-----------------------------------------------------------------------------
  ! Make sure that only one thread reads the namelist
  !-----------------------------------------------------------------------------
  call io%lock%set ('dust_evolve_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_evolve_params, iostat=iostat)
    if (iostat < 0) then                       ! if syntax errors
      rewind (io%input)                        ! rewind to top of file
      read (io%input, dust_evolve_params)      ! read againg to get more info
    end if
    call io%namelist_warning ('dust_evolve_params', iostat)
    write (io%namelist, dust_evolve_params)    ! write to data/run/params.nml
    call self%test                             ! run a self%test once
    first_time = .false.
  end if
  call io%lock%set ('dust_evolve_t%init')
  self%on = on                                 ! task private copy
  self%n_bin = n_bin
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state for one cell, with _in and _out arguments
!> holding particle mass density dd(n_bin) and velocity vd(n_bin,3)
!===============================================================================
SUBROUTINE update (self, dd_in, vd_in, dd_out, vd_out)
  class(dust_evolve_t):: self
  real:: dd_in(:), dd_out(:), vd_in(:,:), vd_out(:,:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%update', itimer=itimer)
  ! ....= ... dd_in
  ! ....= vd_in(i) ...
  ! ... dd_out = ...
  ! ... vd_out(i) = ...
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(dust_evolve_t):: self
  call trace%begin ('dust_evolve_t%test')
  ! ...
  call trace%end ()
END SUBROUTINE test

END MODULE dust_evolve_mod
