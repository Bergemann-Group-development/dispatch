!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE standard_boundaries_mod
  USE patch_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
  ! --------------------------------------------------------------------------
    real, pointer, contiguous:: grav_bot(:,:,:)
    real, pointer, contiguous:: grav_top(:,:,:)
    real, pointer, contiguous:: p_bb(:,:), d_bb(:,:)
    real, pointer, contiguous:: p_ub(:,:), d_ub(:,:)
    integer:: l_z, u_z
    logical:: surface=.false.
    logical:: initialized=.false.
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t)     :: self
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  integer:: iostat
  logical, save:: set_bc_bits=.false.
  logical, save:: first_time=.true.
  namelist /boundary_params/ set_bc_bits
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  patch => self%patch
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, boundary_params)
    end if
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  self%set_bc_bits = set_bc_bits
  call self%standard_boundaries_t%init ()
  self%initialized=.true.
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
END SUBROUTINE derivative

END MODULE boundaries_mod
