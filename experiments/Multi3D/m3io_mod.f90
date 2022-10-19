MODULE m3io_mod 
  implicit none

  integer, parameter,private :: maxlu=500
  logical, dimension(-maxlu:maxlu), private :: lu=.false.
  character(len=128) :: error_mesg

  interface NextEntry
     module procedure NextEntry1, NextEntry2
  end interface NextEntry

contains

  !====================================================================================
  SUBROUTINE NextEntry1(fileunit)
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
    character(len=1)   :: token1
    character(len=1),parameter :: comment='*',empty=' '
    !---------------------------------------------------------------------------------
    do 

       ! read the line
       read(unit=fileunit,fmt='(A72)') instring
       instring=adjustl(instring)
       token1=instring(1:1)

       ! exit if line is not a comment and not empty
       if (token1.ne.comment.and.token1.ne.empty) then
          backspace(fileunit)
          exit
       endif

    end do
  END SUBROUTINE NextEntry1

  !====================================================================================
  SUBROUTINE NextEntry2(fileunit,comment)
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
    character(len=1)   :: token1
    character(len=1),parameter :: empty=' '
    character(len=*),intent(in) :: comment
    !---------------------------------------------------------------------------------
    do 
       ! read the line
       read(unit=fileunit,fmt='(A72)') instring
       instring=adjustl(instring)
       token1=instring(1:1)

       ! exit if line is not a comment and not empty
       if (token1.ne.comment.and.token1.ne.empty) then
          backspace(fileunit)
          exit
       endif
    end do
  END SUBROUTINE NextEntry2

  !====================================================================================
  pure function lowercase(string,ss)
    character(len=*),intent(in) :: string
    integer,intent(in) :: ss

    integer :: i,k
    character(len=ss) ::lowercase
    character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
    !---------------------------------------------------------------------------------
    lowercase(:)=' '
    lowercase(1:ss)=string(1:ss)
    do i=1,ss
       k=index(ucase,string(i:i))
       if (k > 0) then
          lowercase(i:i)=lcase(k:k)
       else
          lowercase(i:i)=string(i:i)
       end if
    end do
  end function lowercase
  
  !====================================================================================
  pure function uppercase(string,ss)
    ! converts the first ss characters of string to uppercase and returns them as a new string.
    implicit none
    character(len=*),intent(in) :: string
    integer,intent(in) :: ss
    character(len=ss) :: uppercase
    integer :: i,k
    character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
    !---------------------------------------------------------------------------------
    uppercase(:)=' '
    uppercase(1:ss)=string(1:ss)
    do i=1,ss
       k=index(lcase,string(i:i))
       if (k > 0) then
          uppercase(i:i)=ucase(k:k)
       else
          uppercase(i:i)=string(i:i)
       end if
    end do
  end function uppercase
END MODULE m3io_mod
