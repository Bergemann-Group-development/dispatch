!===============================================================================
!> Program to test packing of particle data types.
!===============================================================================

PROGRAM test
  USE io_mod
  USE trace_mod
  USE mpi_buffer_mod
  USE particles_mod
  USE omp_lib
  implicit none
  type(particle_t), allocatable:: p(:)
  type(particle_t), pointer:: p1
  type(particles_t), allocatable:: pp(:)
  type(mpi_buffer_t):: mpi_buffer
  real(8):: used
  integer:: n, i, j, ppc, npp, id(5)=0
  integer(8):: m
  !-----------------------------------------------------------------------------
  ! Test packing individual particles
  !-----------------------------------------------------------------------------
  n=2**20
  allocate (p(n))
  m = n*(2+3+3+2+1)+1
  !      t r v q w
  call mpi_buffer%init (m)
  call mpi_buffer%append (n)
  used = omp_get_wtime()
  do i=1,n
    call p(i)%pack(mpi_buffer)
  end do
  used = omp_get_wtime()-used
  print *,used,used/n
  !-----------------------------------------------------------------------------
  ! Test packing particle bunches
  !-----------------------------------------------------------------------------
  ppc = 128
  npp = n/ppc
  allocate (pp(npp))
  do i=1,npp
    call pp(i)%init (id)
    do j=1,ppc
      allocate (p1)
      call pp(i)%append (p1)
    end do
  end do
  call mpi_buffer%reset ()
  m = n*(2+3+3+2+1)+1+npp*3+n*(3 +1 +2)   ! = 11*n + 3*np + 6*n
  !      t r v q w             vg d id
  call mpi_buffer%init (m)
  used = omp_get_wtime()
  print *,npp
  do i=1,npp
    if (i < 10)  print *,i,npp, mpi_buffer%position
    call mpi_buffer%append (ppc)
    call pp(i)%pack (mpi_buffer)
  end do
  used = omp_get_wtime()-used
  print *,used,used/n
END PROGRAM
