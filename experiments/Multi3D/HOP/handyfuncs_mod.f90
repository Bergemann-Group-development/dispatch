MODULE handyfuncs_mod
  USE op_math_mod
  implicit none
  
contains

!===============================================================================
SUBROUTINE count_words (string, nwords)
  character(len=*), intent(in) :: string
  integer, intent(out) :: nwords
  integer :: i, n
  logical :: new_word=.true.
  !-----------------------------------------------------------------------------
  n = len(string)
  nwords = 0
  do i=1, n
  if (string(i:i) == ' ') then
      new_word = .true.
  elseif(new_word) then
      nwords = nwords + 1
      new_word = .false.
  endif
  end do
END SUBROUTINE count_words

!===============================================================================
FUNCTION trim_between_words(string) result(res)
    character(len=*) :: string
    character(len=len(string)) :: res
    character :: char
    integer :: i, j
    !---------------------------------------------------------------------------
    res = ''
    char = ' '
    j = 1
    do i=1, len(string)
        if (string(i:i) == ' ') then
            if (char == ' ') then
                cycle
            else
                res(j:j) = string(i:i)
                j = j + 1
            endif
        else
            res(j:j) = string(i:i)
            j = j + 1
        endif
        char = string(i:i)
    end do
END FUNCTION

!===============================================================================
FUNCTION string_in_text (text, str) result(in_text)
  character(len=*) :: str
  character(len=*) :: text
  character(len=128) :: dummy
  integer :: i, j 
  logical :: in_text
  ! ----------------------------------------------------------------------------
  in_text = .false.
  dummy = ''
  j = 0
  do i=1, len(text)
    if (text(i:i) == ' ') then
      if (dummy(1:j) == str) exit
      j = 0
    else
      j = j + 1
      dummy(j:j) = text(i:i)
    endif
  end do
  if (dummy(1:j) == str) in_text = .true.
END FUNCTION

!===============================================================================
FUNCTION extract_int (str) result(val)
  character(len=*) :: str
  character(len=128) :: dummy
  integer :: i, o, val
  ! ----------------------------------------------------------------------------
  i = scan(str, '0123456789')
  o = scan(str, '0123456789', back=.true.)
  read(str(i:o), *) val
END FUNCTION

!===============================================================================
FUNCTION hash_string(text) result(hashed) 
  character(*), intent(in) :: text 
  integer :: hashed ! 32-bit integers assumed 
  integer, parameter :: magic_numb = int(Z'5D7A9F43', kind=kind(magic_numb))  ! hexadecimal
  integer :: i, j 
  !-----------------------------------------------------------------------------
  hashed = 0
  do i = 1, len_trim(text) 
    j = mod(i-1, 4) * 8 
    hashed = IEOR( hashed, ISHFT( ICHAR( text(i:i) ), j ) ) 
  end do
  hashed = abs( IEOR( hashed, magic_numb ) )
END FUNCTION hash_string 

!===============================================================================
FUNCTION thread2place(thread, ncores, nthreads) result(place)
  integer, intent(in) :: thread, ncores, nthreads
  integer :: place
  ! ----------------------------------------------------------------------------
  place = nint(thread*ncores/real(nthreads))
END FUNCTION thread2place

!===============================================================================
SUBROUTINE crop1d8(arr, i, o)
  real(8), intent(inout), allocatable :: arr(:)
  real(8), allocatable :: tmp(:)
  integer, intent(in) :: i, o
  integer :: n
  ! ----------------------------------------------------------------------------
  n = o-i+1
  allocate(tmp(n))
  tmp = arr(i:o)
  deallocate(arr)
  allocate(arr(n))
  arr = tmp
  deallocate(tmp)
END SUBROUTINE crop1d8

!====================================================================================
PURE FUNCTION lowercase(string)
  character(len=*),intent(in) :: string
  character(len=len(string)) :: lowercase
  character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
  integer :: i,k
  !---------------------------------------------------------------------------------
  lowercase=' '
  lowercase=string
  do i=1, len(string)
     k = index(ucase,string(i:i))
     if (k > 0) then
        lowercase(i:i)=lcase(k:k)
     else
        lowercase(i:i)=string(i:i)
     end if
  end do
END FUNCTION lowercase

!====================================================================================
PURE FUNCTION uppercase(string)
  ! converts the first ss characters of string to uppercase and returns them as a new string.
  implicit none
  character(len=*),intent(in) :: string
  character(len=len(string)) :: uppercase
  character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
  integer :: i,k
  !---------------------------------------------------------------------------------
  uppercase=' '
  uppercase=string
  do i=1, len(string)
     k=index(lcase,string(i:i))
     if (k > 0) then
        uppercase(i:i)=ucase(k:k)
     else
        uppercase(i:i)=string(i:i)
     end if
  end do
END FUNCTION uppercase

!====================================================================================
PURE FUNCTION capitalize(string)
  character(len=*),intent(in) :: string
  character(len=len(string)) :: capitalize
  integer :: i,k
  !---------------------------------------------------------------------------------
  capitalize = lowercase(trim(string))
  capitalize(1:1) = uppercase(capitalize(1:1))
END FUNCTION capitalize

!====================================================================================
SUBROUTINE sort_array(array, indices)
  real(8), intent(inout) :: array(:)
  integer, intent(out) :: indices(size(array))
  real(8) :: tmp
  integer :: i, j, min_idx, idx
  !---------------------------------------------------------------------------------
  ! initialize indices to the identity permutation
  indices = [(i, i=1, size(array))]

  ! selection sort algorithm
  do i = 1, size(array)-1
    min_idx = i
    do j = i+1, size(array)
      if (array(j) < array(min_idx)) then
        min_idx = j
      end if
    end do
    if (min_idx /= i) then
      ! swap elements in array and indices
      tmp = array(i)
      array(i) = array(min_idx)
      array(min_idx) = tmp
      idx = indices(i)
      indices(i) = indices(min_idx)
      indices(min_idx) = idx
    end if
  end do
END SUBROUTINE sort_array

!====================================================================================
SUBROUTINE refine_array(array, n, fac)
  real(8), intent(inout), dimension(:), pointer, contiguous :: array
  integer, intent(inout) :: n
  integer, intent(in) :: fac
  real(8), dimension(:), pointer, contiguous :: tmp
  integer :: i, nnew
  !---------------------------------------------------------------------------------
  nnew = (n-1)*fac + 1
  allocate(tmp(nnew))
  tmp = interp1d([(dble(i), i=0, n-1)], array, [(dble(i)/fac, i=0, (n-1)*fac)], order=1, sorted=.true.)
  deallocate(array)
  array => tmp
  n = nnew
END SUBROUTINE refine_array

!===============================================================================
SUBROUTINE remove_close_entries(array, threshold)
  real(8), intent(inout), allocatable :: array(:)
  real(8), intent(in) :: threshold
  real(8), allocatable :: scratch(:)
  integer :: i, new_size
  !-----------------------------------------------------------------------------
  new_size = size(array)
  i = 1
  do while (i < new_size)
    if (abs(array(i+1) - array(i)) < threshold) then
      array(i+1:new_size-1) = array(i+2:new_size)
      new_size = new_size - 1
    else
      i = i + 1
    end if
  end do
  
  allocate(scratch, source=array(1:new_size))
  deallocate(array)
  allocate(array, source=scratch)
  deallocate(scratch)
END SUBROUTINE remove_close_entries

!===============================================================================
elemental FUNCTION air2vac (air) result(vac)
  real(8), intent(in) :: air
  real(8) :: vac, n, s, s2
  !-----------------------------------------------------------------------------
  s = 1d4 / air
  s2 = s*s
  n = 1d0 + 0.00008336624212083d0 + 0.02408926869968d0 / (130.1065924522d0 - s2) &
        + 0.0001599740894897d0 / (38.92568793293d0 - s2)
  vac = air*n
END FUNCTION air2vac

END MODULE handyfuncs_mod