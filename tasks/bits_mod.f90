!> =============================================================================
!> $Id$
!> Status bits for tasks
!> =============================================================================
MODULE bits_mod
  integer, parameter::       static_bit=1
  integer, parameter::         leaf_bit=2
  integer, parameter::        ready_bit=2**2
  integer, parameter::     internal_bit=2**3
  integer, parameter::      virtual_bit=2**4
  integer, parameter::     boundary_bit=2**5
  integer, parameter::     external_bit=2**6
  integer, parameter::       frozen_bit=2**7
  integer, parameter:: swap_request_bit=2**8
  integer, parameter::    root_grid_bit=2**9
  integer, parameter::         star_bit=2**10
  integer, parameter::       planet_bit=2**11
  integer, parameter::         disk_bit=2**12
  integer, parameter::     particle_bit=2**13
  integer, parameter::        retry_bit=2**14
  integer, parameter::       remove_bit=2**15
  type bits_t
    integer:: remove       = remove_bit
    integer:: retry        = retry_bit
    integer:: ready        = ready_bit
    integer:: internal     = internal_bit
    integer:: virtual      = virtual_bit
    integer:: boundary     = boundary_bit
    integer:: external     = external_bit
    integer:: root_grid    = root_grid_bit
    integer:: not_leaf     = leaf_bit
    integer:: frozen       = frozen_bit
    integer:: swap_request = swap_request_bit
    integer:: particle     = particle_bit
    integer:: star         = star_bit
    integer:: planet       = planet_bit
    integer:: static       = static_bit
  end type
  type(bits_t):: bits
END MODULE bits_mod
