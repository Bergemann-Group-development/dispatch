MODULE handyfuncs_mod
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
FUNCTION hash_string(text) result(hashed) 
  character(*), intent(in) :: text 
  integer :: hashed ! 32-bit integers assumed 
  integer, parameter :: magic_numb = Z'5D7A9F43'  ! hexadecimal
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
pure function lowercase(string)
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
end function lowercase

!====================================================================================
pure function uppercase(string)
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
end function uppercase

!====================================================================================
pure function option(bool) result(out)
  logical, intent(in), optional :: bool
  logical :: out
  !---------------------------------------------------------------------------------
  if (present(bool)) then
    out = bool
  else
    out = .false.
  endif
end function option

END MODULE handyfuncs_mod