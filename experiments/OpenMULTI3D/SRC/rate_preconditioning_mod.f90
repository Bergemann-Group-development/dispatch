!===============================================================================
!> Keep the largest term of a running sum apart from the rest, so that one
!> term can later be removed from the total without a cancelling subtraction.
!===============================================================================
MODULE rate_preconditioning_mod
  implicit none
  private
  public :: add_largest, other_contribution
CONTAINS
  !=============================================================================
  !> Keep the largest-magnitude term separate as contributions arrive.  The
  !> remainder never contains that term, even when a later term replaces it.
  !=============================================================================
  pure subroutine add_largest (value, index, largest, rest, largest_index)
    real(4), intent(in) :: value
    integer, intent(in) :: index
    real(4), intent(inout) :: largest, rest
    integer, intent(inout) :: largest_index
    if (abs(value) > abs(largest)) then
      rest = rest + largest
      largest = value
      largest_index = index
    else
      rest = rest + value
    endif
  end subroutine add_largest

  pure function other_contribution (value, index, largest, rest, largest_index) result(other)
    real(4), intent(in) :: value, largest, rest
    integer, intent(in) :: index, largest_index
    real(4) :: other
    if (index == largest_index) then
      other = rest
    else
      other = largest + (rest-value)
    endif
  end function other_contribution

END MODULE rate_preconditioning_mod
