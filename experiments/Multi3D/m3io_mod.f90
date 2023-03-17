MODULE m3io_mod 
  implicit none

  integer, parameter,private :: maxlu=500
  logical, dimension(-maxlu:maxlu), private :: lu=.false.
  character(len=128) :: error_mesg

contains

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
END MODULE m3io_mod
