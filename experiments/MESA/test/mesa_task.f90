PROGRAM mesa_task
  implicit none
  integer:: iostat, i
  character(len=32):: line
  !-----------------------------------------------------------------------------
  print *, 'MESA task starting'
  call sleep(2)
  !-----------------------------------------------------------------------------
  print *, 'MESA task making pipes'
  call system("rm -f mesa.pipe1 mesa.pipe2")
  call system("mknod mesa.pipe1 p")
  call system("mknod mesa.pipe2 p")
  call sleep(1)
  !-----------------------------------------------------------------------------
  print *, 'MESA task starting MESA code'
  call system("./mesa.sh")
  !-----------------------------------------------------------------------------
  open (10, file='mesa.pipe1', form='formatted', status='unknown', action='read')
  open (11, file='mesa.pipe2', form='formatted', status='unknown', action='write')
  !-----------------------------------------------------------------------------
  i = 0
  do
    call sleep(5)
    write (line,*) i
    write (11,*) line
    print *, 'MESA task sent this: ', trim(line)
    read (10,*, iostat=iostat) line
    if (iostat /= 0) exit
    print *, 'MESA task  got this: ', trim(line)
    i = i+5
  end do
  print *, 'exiting with iostat =', iostat
END
