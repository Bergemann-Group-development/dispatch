!> =============================================================================
!> =============================================================================
MODULE sinkparticle_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE particle_history_mod
  USE download_mod
  USE data_hub_mod
  implicit none
  private
  type, public, extends (particle_t):: sinkparticle_t
    type (patch_t):: patch
  contains
    procedure:: init
    procedure:: dnload
    procedure:: update
    procedure:: accrete
    procedure:: create_sink
  end type
  real:: accretion_radius=8.
  real:: exclusion_radius=16.
  real:: jeans_resolution=2.
  type(data_hub_t):: data_hub
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class (sinkparticle_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /sinkparticle_params/ accretion_radius, exclusion_radius, jeans_resolution
  !-----------------------------------------------------------------------------
  !$omp critical (sinkparticle_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, sinkparticle_params, iostat=iostat)
    if (io%master) write (io%output, sinkparticle_params)
  end if
  !$omp end critical (sinkparticle_cr)
  self%patch%n = 2*ceiling(accretion_radius)
  self%patch%ng = 0
  self%patch%nv = 4
  self%patch%nw = 2
  call self%patch%init
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE dnload (self, only)
  class (sinkparticle_t):: self
  integer, optional:: only
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /sinkparticle_params/ accretion_radius, exclusion_radius, jeans_resolution
  !-----------------------------------------------------------------------------
  if (self%patch%use_data_hub) then
    call data_hub%dnload (only=self%patch%i_d)
  else
    call download%download_link (self%patch%link, only=self%patch%i_d)
  end if
END SUBROUTINE dnload

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class (sinkparticle_t):: self
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
SUBROUTINE accrete (self)
  class (sinkparticle_t):: self
END SUBROUTINE accrete

!> =============================================================================
!> =============================================================================
SUBROUTINE create_sink (self)
  class (sinkparticle_t):: self
END SUBROUTINE create_sink

END MODULE sinkparticle_mod
