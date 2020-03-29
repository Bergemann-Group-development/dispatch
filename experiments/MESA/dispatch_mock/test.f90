program test
    use mesa_mod
    implicit none
    real, parameter :: secyer = 3.1558149984e7 ! seconds per year
    real, parameter :: msun   = 1.9892e33      ! solar mass (g)
  
    type(mesa_sink) :: sink
    sink = create_sink(0.)
    call submit_rate(sink, 1e-6*msun/secyer, 1.)
    call submit_rate(sink, 1e-6*msun/secyer, 10*secyer)
    write(*,*) "hi"
  
  
end program test