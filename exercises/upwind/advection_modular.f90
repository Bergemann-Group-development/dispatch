!*******************************************************************************
PROGRAM advection_test
  USE mesh                                                                      ! mesh module
  USE operators                                                                 ! div() operator
  USE boundaries                                                                ! ghost zone handling
  USE mpi_coords                                                                ! MPI dimensions
  USE mpi_reduce                                                                ! max, min, average
  USE io                                                                        ! input/output module
  implicit none
  character(len=mch), save:: id= &
    'main.f90 $Id$'
  integer :: istep, iout
  real, allocatable, dimension(:,:,:):: rho                                     ! density field
  real, allocatable, dimension(:,:,:) :: rho_interior                           ! without ghost zones
  real:: dt                                                                     ! time step
  integer:: nstep=100                                                           ! number of time steps
  integer:: fout                                                                ! output frequency
  integer:: ix, iy, iz                                                          ! loop indices
  real:: u(3) = (/1.0,0.5,0.0/)                                                 ! velocity components
  real:: width=10.0                                                             ! Gaussian width
  real:: courant=0.3                                                            ! Courant number
  namelist /advection/ u, width, courant, nstep, fout                           ! input namelist
!-------------------------------------------------------------------------------
! Initialize mesh
!-------------------------------------------------------------------------------
  call init_mpi                                                                 ! start mpi
  call print_id(id)                                                             ! identify version
  call init_io                                                                  ! setup input/output
  call init_mesh                                                                ! setup mesh
  call init_boundary
!-------------------------------------------------------------------------------
! Setup initial conditions
!-------------------------------------------------------------------------------
  rewind (input); read (input,advection); if (master) write(*,advection)        ! read namelist input
  allocate (rho(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub))               ! allocate density
  do iz=0,m(3)%n-1; do iy=0,m(2)%n-1; do ix=0,m(1)%n-1                          ! loop over 3-D
    rho(ix,iy,iz) = exp(-(m(1)%r(ix)**2+m(2)%r(iy)**2+m(3)%r(iz)**2)/width**2)  ! Gaussian, centered
  end do; end do; end do
  call boundary (rho)                                                           ! fill ghost zones
!-------------------------------------------------------------------------------
! Setup output
!-------------------------------------------------------------------------------
  allocate(rho_interior(m(1)%n,m(2)%n,m(3)%n))                                  ! buffer for I/O
  iout = 0                                                                      ! zero based snapshot
  call file_openw_mpi ('snapshot.dat', m%gn, m%n, m%offset)                     ! open data file
  call write_field(rho)                                                         ! dump data
!-------------------------------------------------------------------------------
! Evolve the solution in time
!-------------------------------------------------------------------------------
  dt = courant*minval(m%d)/sqrt(sum(u**2))                                      ! dt = C dx/u
  do istep = 1, nstep                                                           ! time step loop

    ! Change here to use upwind divergence instead
    rho = rho - dt*div(rho*u(1),rho*u(2),rho*u(3))                              ! advection equation


    if (modulo(istep,fout)==0) call write_field(rho)                            ! dump data
  enddo
!-------------------------------------------------------------------------------
! Clean up
!-------------------------------------------------------------------------------
  call file_close_mpi                                                           ! close snapshot file
  deallocate (rho)                                                              ! dellocate density
  call end_mpi                                                                  ! end mpi
CONTAINS

!*******************************************************************************
SUBROUTINE write_field(f)
  implicit none
  real f(m(1)%lb:m(1)%ub, m(2)%lb:m(2)%ub, m(3)%lb:m(3)%ub)
  real:: mn, av, mx
!...............................................................................
  rho_interior = f(0:m(1)%n-1,0:m(2)%n-1,0:m(3)%n-1)                            ! skip ghost zones
  call file_write_mpi (rho_interior, m%n, iout)                                 ! collective write
  mn = minval(rho_interior); call min_real_mpi(mn)                              ! global min
  av = sum(rho_interior)/product(m%gn); call sum_real_mpi(av)                   ! global average
  mx = maxval(rho_interior); call max_real_mpi(mx)                              ! global max
  if (master) print '(a,i5,3e12.4)', 'snap, min, aver, max :', iout, mn, av, mx ! printout
  iout = iout + 1                                                               ! increment snapshot
END SUBROUTINE write_field
END PROGRAM
