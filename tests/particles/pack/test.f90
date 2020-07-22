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
  type(particle_t), allocatable:: pn(:), pv(:)       ! particle arrays
  class(particle_t), pointer:: pn1, pv1      ! viz particle
  type(particles_t), allocatable:: cells(:)          ! cell info
  type(mpi_buffer_t):: mpi_buffer, viz_buffer        ! separate buffers
  real(8):: used
  integer:: nn, nv, nc, np, i, j, ppc, npp, id(5)=0, every=1000, sz(2)
  integer(8):: mp, mv
  !-----------------------------------------------------------------------------
  ! Find the sizes of particles in the buffers
  !-----------------------------------------------------------------------------
  nn = 2**21                      ! number of MPI-only particles
  print '(f6.1,a)', nn/1e6, ' million particles'
  mp = 1000
  call mpi_buffer%init (mp)
  call viz_buffer%init (mp)
  allocate (particle_t:: pn1)
  allocate (vparticle_t:: pv1)
  call pn1%pack (mpi_buffer)
  call pv1%pack (viz_buffer)
  sz(1) = mpi_buffer%position
  sz(2) = viz_buffer%position
  print *, 'sizes of single particles:', sz
  call mpi_buffer%reset
  call viz_buffer%reset
  !-----------------------------------------------------------------------------
  ! Test packing individual particles
  !-----------------------------------------------------------------------------
  nv = nn/every                   ! number of viz particles
  np = nn + nv                    ! total number of particles
  allocate (pn(nn), pv(nv))
  mp = np*sz(1)+2                 ! buffer size, all particles
  mv = nv*sz(2)+2                 ! buffer size, only viz particles
  !        t r v q w d vg id 
  call mpi_buffer%init (mp)       ! buffer for MPI-only
  call viz_buffer%init (mv)       ! buffer for viz
  call mpi_buffer%append (mp)     ! size of MPI buffer
  call mpi_buffer%append (mv)     ! size of viz buffer
  used = omp_get_wtime()
  do i=1,nn
    call pn(i)%pack(mpi_buffer)   ! pack particle info into MPI-buffer
  end do
  do i=1,nv
    !print*,i,nv
    call pv(i)%pack(mpi_buffer)   ! pack viz info into MPI-buffer
    call pv(i)%pack(viz_buffer)   ! pack viz info into viz-buffer
  end do
  used = omp_get_wtime()-used
  print 1, used, used/np
  1 format(f6.2,' sec',9p,f8.1,' ns/p')
  !-----------------------------------------------------------------------------
  ! Divide the particles among a number of cells, with 128 in each
  !-----------------------------------------------------------------------------
  ppc = 128                           ! particles per cell
  nc = np/ppc                         ! number of cells
  allocate (cells(nc))                ! allocate cells
  do i=1,nc
    call cells(i)%init (id)           ! initialize cell
    do j=1,ppc
      allocate (particle_t:: pn1)     ! make a normal particle
      call cells(i)%append (pn1)      ! append to cell
    end do
    if (mod(i,every)==1) then         ! every so often also add viz particles
      do j=1,ppc
        allocate (vparticle_t:: pv1)  ! allocate viz particle
        call cells(i)%append (pv1)    ! append to cell
      end do
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Find the sizes of particles in the buffers
  !-----------------------------------------------------------------------------
  mp = 10000
  call mpi_buffer%init (mp)
  call viz_buffer%init (mp)
  call cells(1)%pack (mpi_buffer, viz_buffer)
  sz(1) = mpi_buffer%position
  sz(2) = viz_buffer%position
  print *, 'sizes of cell buffers:', sz
  call mpi_buffer%reset
  call viz_buffer%reset
  !-----------------------------------------------------------------------------
  ! Test particles_t%pack()
  !-----------------------------------------------------------------------------
  call mpi_buffer%reset ()
  call viz_buffer%reset ()
  mp = nc*sz(1) + 2
  mv = nc*sz(2) + 2
  call mpi_buffer%init (mp)
  call viz_buffer%init (mv)
  used = omp_get_wtime()
  do i=1,nc
    call cells(i)%pack (mpi_buffer, viz_buffer)
  end do
  used = omp_get_wtime()-used
  print 1, used, used/np
  !-----------------------------------------------------------------------------
  ! Write out the buffer (into fort.1), for Python development
  !-----------------------------------------------------------------------------
  write (1) viz_buffer%buffer
END PROGRAM