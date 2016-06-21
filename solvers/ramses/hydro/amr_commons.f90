module amr_commons
  use amr_parameters

  ! Define default MPI types for transferring kind=dp and kind=qdp data
#ifndef WITHOUTMPI
  use mpi, only : MPI_REAL4, MPI_REAL8, MPI_REAL16
  private :: MPI_REAL4, MPI_REAL8, MPI_REAL16
#ifndef NPRE
  integer, parameter :: DPMPI_REAL=MPI_REAL4 ! default
#else
#if NPRE==4
  integer, parameter :: DPMPI_REAL=MPI_REAL4 ! real*4
#else
  integer, parameter :: DPMPI_REAL=MPI_REAL8 ! real*8
#endif
#endif
#ifdef QUADHILBERT
  integer, parameter :: QDPMPI_REAL=MPI_REAL16 ! real*16
#else
  integer, parameter :: QDPMPI_REAL=MPI_REAL8 ! real*8
#endif
#endif
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level

end module
