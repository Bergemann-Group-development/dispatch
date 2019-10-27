!==============================================================================
!==============================================================================
MODULE cooling_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  implicit none
  integer, parameter:: dp=4
  real, parameter:: T_MC=10.0
  private
  type, public:: cooling_t
    procedure(cooling_Casiana), pointer:: function => null()
  contains
    procedure:: init
    procedure:: pre_update
  end type
  logical:: on=.false.
  integer:: type=0, verbose=0
  type(cooling_t), public:: cooling
CONTAINS

!==============================================================================
!==============================================================================
SUBROUTINE init (self)
 class (cooling_t):: self
 integer:: iostat
 logical, save:: first_time=.true.
 namelist /cooling_params/ on, verbose, type
 !------------------------------------------------------------------------------
 call trace%begin ('simple/cooling_t%init')
 !$omp critical (input_cr)
 if (first_time) then
   first_time = .false.
   rewind (io%input)
   read (io%input, cooling_params, iostat=iostat)
   write (io%output, cooling_params)
 end if
 !$omp end critical (input_cr)
 select case (type)
 case (1)
   self%function => cooling_Dalgarno
 case (2)
   self%function => cooling_Gnedin
 case (3)
   self%function => cooling_Rosen_etal
 case (4)
   self%function => cooling_Shure
 case default
   self%function => cooling_Casiana
 end select
 call trace%end
END SUBROUTINE init

!==============================================================================
!> Implement actual cooling
!==============================================================================
SUBROUTINE pre_update (self, patch)
 class (cooling_t):: self
 class (gpatch_t), pointer:: patch
 if (on) then
   call trace%begin ('simple/cooling_t%pre_update')
   call trace%end
 end if
END SUBROUTINE pre_update

!==============================================================================
!==============================================================================
include "cooling_Casiana.f90"
include "cooling_Dalgarno.f90"
include "cooling_Gnedin.f90"
include "cooling_Rosen_etal.f90"
include "cooling_Shure.f90"

END MODULE cooling_mod
