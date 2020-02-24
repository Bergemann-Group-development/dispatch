PROGRAM mesa_code
  implicit none
  integer:: iostat, i
  character(len=32):: line
  !-----------------------------------------------------------------------------
  open (10, file='mesa.pipe1', form='formatted', status='unknown', action='write')
  open (11, file='mesa.pipe2', form='formatted', status='unknown', action='read')
  print '(40x,a)', 'MESA code ready'
  !-----------------------------------------------------------------------------
  do
    read (11,*) line
    print '(40x,a)', 'MESA code  got this: '//line
    read (line,*) i
    i = i+1
    write (line,*) i
    write (10,*) line
    print '(40x,a)', 'MESA code sent this: '//line
  end do
END
