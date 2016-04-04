
SUBROUTINE subr (ndata, data)
  integer:: ndata, bytes_per_word=4
  integer, dimension(ndata):: data
  open (10, file='test.dat', form='unformatted', access='direct', &
    recl=size(data)*bytes_per_word)
  write (10, rec=1) data
  close (10)
END SUBROUTINE
