PROGRAM mesa_code
  USE io_unit_mod
  implicit none
  integer:: iostat, i
  character(len=32):: line
  !-----------------------------------------------------------------------------
  open (io_unit%mesa1, file='mesa.pipe1', form='formatted', status='unknown', action='write')
  open (io_unit%mesa2, file='mesa.pipe2', form='formatted', status='unknown', action='read')
  print '(40x,a)', 'MESA code ready'
  !-----------------------------------------------------------------------------
  do
    read (11,*) line
    print '(40x,a)', 'MESA code  got this: '//line
    read (line,*) i
    i = i+1
    if (i > 10) then
      exit
    end if
    write (line,*) i
    write (10,*) line
    print '(40x,a)', 'MESA code sent this: '//line
  end do
  print *,'MESA exiting normally, after read, before write'
END
