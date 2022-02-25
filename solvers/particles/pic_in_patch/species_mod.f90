!===============================================================================
!===============================================================================
MODULE species_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE patch_mod
  USE particles_mod
  implicit none
  private
  type, public:: species_t
    integer :: n = 2
    integer :: per_cell = 10
    integer :: lb(3), ub(3)
    class(patch_t), pointer:: patch
    class(particles_t), pointer:: particles(:)
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
  end type
  integer:: n = 2
  integer:: per_cell = 10
  integer:: verbose = 1
  logical:: first_time = .true.
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(species_t):: self
  class(patch_t), pointer:: patch
  integer:: iostat, i
  namelist /species_params/ n, per_cell, verbose
  !-----------------------------------------------------------------------------
  self%patch => patch
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read   (io%input, species_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read   (io%input, species_params)
    end if
    call io%namelist_warning ('species_params', iostat)
    write (io%namelist, species_params)
  end if
  call io%lock%unset
  self%n = n
  !-----------------------------------------------------------------------------
  allocate (self%particles(n))
  do i=1,n
    self%particles(i)%lb = self%lb
    self%particles(i)%ub = self%ub
    call self%particles(i)%pre_init (patch, per_cell)
  end do
END SUBROUTINE pre_init

!===============================================================================
!===============================================================================
SUBROUTINE init (self, patch)
  class(species_t):: self
  class(patch_t), pointer:: patch
  integer:: i, n
  !-----------------------------------------------------------------------------
  do i=1,n
    call self%particles(i)%init (patch, self%per_cell)
  end do
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self, B, E, J)
  class(species_t):: self
  real(var_k), contiguous, dimension(:,:,:,:):: B, E, J
  integer:: i
  !-----------------------------------------------------------------------------
  do i=1,n
    call self%particles(i)%update (B, E, J)
  end do
END SUBROUTINE update

END MODULE species_mod
