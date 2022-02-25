
PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: unit, i, rec
  integer, parameter:: n=32, nrec=1000
  real:: buf(n,n,n)
  character(len=32) filename

  call init_omp

  call tic
  !$omp parallel private(filename, unit,i, rec)
  write (filename,'("data/",i4.4,".dat")') unit
  !$omp critical
  filename="data/crit.dat"
  unit=10+omp_mythread()
  open (unit=unit,file=filename,form='unformatted', &
    access='direct',recl=sizeof(buf),status='unknown')
  do i=1,nrec
    rec=i+nrec*omp_mythread()
    write (unit,rec=rec) buf
  end do
  close (unit)
  !$omp end critical
  !$omp end parallel
  call timer%toc ('critical')

  call tic
  !$omp parallel private(filename, unit,i, rec)
  filename="data/shared.dat"
  do i=1,nrec
    rec=i+nrec*omp_mythread()
    unit=10+omp_mythread()
    !$omp critical
    open (unit=unit,file=filename,form='unformatted', &
      access='direct',recl=sizeof(buf),status='unknown')
    write (unit,rec=rec) buf
    close (unit)
    !$omp end critical
  end do
  !$omp end parallel
  call timer%toc ('shared')

END PROGRAM
