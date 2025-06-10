MODULE handyfuncs_mod
  USE op_math_mod
  implicit none

  interface argsort
     module procedure argsort_int, argsort_r4, argsort_r8
  end interface argsort

  abstract interface
    function generic_func(input) result(f)
      real(8), intent(in) :: input
      real(8) :: f
    end function generic_func
  end interface
  
contains

!====================================================================================
SUBROUTINE print_failed_line(lunit, iostat)
  integer, intent(in) :: lunit, iostat
  character(len=512) :: line_text, filename
  ! ----------------------------------------------------------------------
  if (iostat > 0) then
      inquire(unit=lunit, name=filename)
      backspace lunit
      read(lunit,'(A)') line_text
      print '(I4, 1x, A)', iostat, "Failed to read a line in file:"
      print '(A)', trim(filename)
      print '(A)', line_text
      stop
  endif
END SUBROUTINE print_failed_line

!====================================================================================
SUBROUTINE NextEntry(fileunit, comment)
  ! ***********************************************************************
  ! --- Reads an input line, returns when it is not a 
  ! --- comment or an empty line.
  ! --- Feb 21 2005: 
  ! ---   Created by Jorrit Leenaarts (j.leenaarts@astro.uu.nl)
  ! --- May 28 2008 :
  !       jorritl: updated, will now also check for empty lines,
  !       and allow data in the first column of the file.
  ! ***********************************************************************
  integer,intent(in) :: fileunit
  character(len=132) :: instring
  character(len=1)   :: token1, cmt
  character(len=*),intent(in), optional :: comment
  !---------------------------------------------------------------------------------
  if (present(comment)) then
    cmt = comment
  else
    cmt = '*'
  endif
  do 
     ! read the line
     read(unit=fileunit,fmt='(A72)') instring
     instring=adjustl(instring)
     token1=instring(1:1)

     ! exit if line is not a comment and not empty
     if (token1 .ne. cmt .and. token1 .ne. ' ') then
        backspace(fileunit)
        exit
     endif
  end do
END SUBROUTINE NextEntry

!===============================================================================
SUBROUTINE unequal_intervals(n, i, ni, start, end)
  integer, intent(in) :: n, i, ni
  integer, intent(out) :: start, end
  integer :: nv, md, j
  !-----------------------------------------------------------------------------
  nv = n / ni
  md = mod(n, ni)
  j = mod(abs(i)-1, ni)+1
  start = (j-1)*nv + 1 + min(md, j-1)
  if (j <= md) then
    nv = nv + 1
  endif
  end = start + nv - 1
END SUBROUTINE unequal_intervals

!===============================================================================
SUBROUTINE replace_char (string, out, in)
  character(len=*), intent(inout) :: string
  character, intent(in) :: out, in
  integer :: i, n
  !-----------------------------------------------------------------------------
  n = len(string)
  do i=1, n
  if (string(i:i) == out) then
    string(i:i) = in
  endif
  end do
END SUBROUTINE replace_char

!===============================================================================
SUBROUTINE replace_substring (string, out, in)
  character(len=*), intent(inout) :: string
  character(len=*), intent(in) :: out, in
  character(:), allocatable :: tmp
  integer :: i, j, n, no, ni
  !-----------------------------------------------------------------------------
  n = len(string)
  no = len(out)
  ni = len(in)
  allocate(tmp, mold=string)
  i = 1
  j = 1
  do while (i<=n-no+1)
      if (string(i:i+no-1) == out) then
        tmp(j:j+ni-1) = in
        j = j + ni
        i = i + no
      else
        tmp(j:j) = string(i:i)
        i = i + 1
        j = j + 1
      endif
  end do
  string = tmp
END SUBROUTINE replace_substring

!===============================================================================
SUBROUTINE count_words (string, nwords, sep)
  character(len=*), intent(in) :: string
  integer, intent(out) :: nwords
  character, intent(in), optional :: sep
  character :: seperator
  integer :: i, n
  logical :: new_word
  !-----------------------------------------------------------------------------
  if (present(sep)) then
    seperator = sep
  else
    seperator = ' '
  endif

  n = len(string)
  nwords = 0
  new_word = .true.
  do i=1, n
  if (string(i:i) == seperator) then
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
      if (string(i:i) == "'") then
        string(i:i) = ' '
      endif
    end do
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
  character(len=256) :: dummy
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
  character(len=256) :: dummy
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
SUBROUTINE argsort_int(a, order)
  integer, intent(in) :: a(:)
  integer, intent(out) :: order(:)
  integer, allocatable :: scr(:)
  integer :: i, j, n
  !-----------------------------------------------------------------------------
  allocate(scr, source=a)
  n = size(a, 1)
  order = 1

  do i = 2, n
    j = i - 1
    do while (j >= 1)
      if (scr(j) <= a(i)) exit
      scr(j + 1) = scr(j)
      order(j+1) = order(j)
      j = j - 1
    end do
    scr(j + 1) = a(i)
    order(j+1) = i
  end do
END SUBROUTINE argsort_int

!===============================================================================
SUBROUTINE argsort_r4(a, order)
  real(4), intent(in) :: a(:)
  integer, intent(out) :: order(:)
  real(4), allocatable :: scr(:)
  integer :: i, j, n
  !-----------------------------------------------------------------------------
  allocate(scr, source=a)
  n = size(a, 1)
  order = 1

  do i = 2, n
    j = i - 1
    do while (j >= 1)
      if (scr(j) <= a(i)) exit
      scr(j + 1) = scr(j)
      order(j+1) = order(j)
      j = j - 1
    end do
    scr(j + 1) = a(i)
    order(j+1) = i
  end do
END SUBROUTINE argsort_r4

!===============================================================================
SUBROUTINE argsort_r8(a, order)
  real(8), intent(in) :: a(:)
  integer, intent(out) :: order(:)
  real(8), allocatable :: scr(:)
  integer :: i, j, n
  !-----------------------------------------------------------------------------
  allocate(scr, source=a)
  n = size(a, 1)
  order = 1

  do i = 2, n
    j = i - 1
    do while (j >= 1)
      if (scr(j) <= a(i)) exit
      scr(j + 1) = scr(j)
      order(j+1) = order(j)
      j = j - 1
    end do
    scr(j + 1) = a(i)
    order(j+1) = i
  end do
END SUBROUTINE argsort_r8

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

!===============================================================================
FUNCTION find_closest_from_center(arr, tar) result(closest_index)
  real(8), intent(in) :: arr(:)
  real(8), intent(in) :: tar
  integer :: closest_index
  integer :: n, i, center
  real(8) :: min_diff, current_diff
  !-----------------------------------------------------------------------------
  ! Start from the center
  n = size(arr)
  center = (n + 1) / 2
  closest_index = center
  min_diff = abs(arr(center) - tar)

  ! Move right from center
  do i = center + 1, n
      current_diff = abs(arr(i) - tar)
      if (current_diff < min_diff) then
          min_diff = current_diff
          closest_index = i
      else
          exit  ! Stop if values start diverging
      end if
  end do

  ! Move left from updated closest_index
  do i = closest_index - 1, 1, -1
      current_diff = abs(arr(i) - tar)
      if (current_diff < min_diff) then
          min_diff = current_diff
          closest_index = i
      else
          exit
      end if
  end do
END FUNCTION find_closest_from_center

!===============================================================================
FUNCTION bisect(f_ptr, a, b, tol) result(c)
  procedure(generic_func), pointer :: f_ptr
  real(8), intent(in) :: a, b, tol
  real(8) :: c
  real(8) :: fa, fb, fc
  real(8) :: left, right
  !-----------------------------------------------------------------------------
  left = a
  right = b
  fa = f_ptr(left)
  fb = f_ptr(right)

  if (fa * fb > 0.0) then
     print *, "Root not bracketed"
     c = -1.0
     return
  elseif (fa == 0.0) then
     c = left
     return
  elseif (fb == 0.0) then
     c = right
     return
  end if

  do
     c = 0.5 * (left + right)
     fc = f_ptr(c)
     if (abs(fc) < tol .or. abs(right - left) < tol) exit
     if (fa * fc < 0.0) then
        right = c
        fb = fc
     else
        left = c
        fa = fc
     end if
  end do
END FUNCTION bisect

END MODULE handyfuncs_mod