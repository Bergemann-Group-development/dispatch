PROGRAM test
  USE io_mod
  USE output_mod
  real, dimension(:), pointer:: data
  integer:: ndata=10
  !-----------------------------------------------------------------------------
  call io%init ('test')
  allocate (data(ndata))
  print*, size(data), sizeof(data)
  call subr (ndata, data)
  call output%test
END PROGRAM test
