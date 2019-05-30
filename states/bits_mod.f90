!===============================================================================
!> $Id$
!> Status bits for tasks
!===============================================================================
MODULE bits_mod
  integer, parameter::       static_bit=1
  integer, parameter::         leaf_bit=2
  integer, parameter::        ready_bit=2**2
  integer, parameter::     internal_bit=2**3
  integer, parameter::     boundary_bit=2**4
  integer, parameter::      virtual_bit=2**5
  integer, parameter::     external_bit=2**6
  integer, parameter:: swap_request_bit=2**7
  integer, parameter::      support_bit=2**8
  integer, parameter::    root_grid_bit=2**9
  integer, parameter::         star_bit=2**10
  integer, parameter::       planet_bit=2**11
  integer, parameter::         disk_bit=2**12
  integer, parameter::     particle_bit=2**13
  integer, parameter::        retry_bit=2**14
  integer, parameter::       remove_bit=2**15
  integer, parameter::   top_bottom_bit=2**16
! NOTE: the first (or last) 16 bits are masked out in some contexts; please CHECK
  integer, parameter::       has_BC_bit=2**30
  integer, parameter::        trace_bit=2**29
  integer, parameter::       frozen_bit=2**28
  integer, parameter::         busy_bit=2**27
  integer, parameter::  no_download_bit=2**26
  integer, parameter::   no_density_bit=2**25
  integer, parameter::        no_rt_bit=2**24
  integer, parameter::        nbors_bit=2**23   ! trigger init_nbors on vnbor
  integer, parameter::      counter_bit=2**20
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
    integer:: support      = support_bit
    integer:: particle     = particle_bit
    integer:: star         = star_bit
    integer:: planet       = planet_bit
    integer:: static       = static_bit
    integer:: top_bottom   = top_bottom_bit
    integer:: has_BC       = has_BC_bit
    integer:: trace        = trace_bit
    integer:: busy         = busy_bit
    integer:: no_download  = no_download_bit
    integer:: no_density   = no_density_bit
    integer:: no_rt        = no_rt_bit
    integer:: init_nbors   = nbors_bit
    integer:: has_counter  = counter_bit
    integer:: hands_off    = ready_bit+busy_bit+virtual_bit+frozen_bit+external_bit
    !-------- BOUNDARY CONDITION BITS ------------------------------------------
    integer:: xl           = 2**1
    integer:: xu           = 2**2
    integer:: yl           = 2**3
    integer:: yu           = 2**4
    integer:: zl           = 2**5
    integer:: zu           = 2**6
    integer:: spherical    = 2**7
    integer:: symmetric    = 2**8
    integer:: antisymmetric= 2**9
    integer:: extrapolate  = 2**10
  end type
  type(bits_t):: bits
END MODULE bits_mod
