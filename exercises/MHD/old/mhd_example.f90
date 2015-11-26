!*******************************************************************************
PROGRAM mhd_test
  USE mhd_mod                                                                   ! MHD equations
  USE vectors                                                                   ! vector operators
  USE mesh                                                                      ! mesh module
  USE boundaries                                                                ! ghost zone handling
  USE mpi_coords                                                                ! MPI dimensions
  USE mpi_reduce                                                                ! max, min, average
  USE io                                                                        ! input/output module
  USE velocity_fields                                                           ! initial velocities
  USE timestep                                                                  ! time step routine
  implicit none
  character(len=mch), save:: id= &
    'main.f90 $Id$'
  integer :: istep, iout                                                        ! step and output
  type(mhd_t):: v                                                               ! variables
  integer:: nstep=100                                                           ! number of time steps
  integer:: fout=10                                                             ! output frequency
  integer:: ix, iy, iz                                                          ! loop indices
  real:: width=10.0, density=2.0, entropy=-1.0                                  ! Gaussian, density, entropy
  namelist /example/ width, nstep, fout, density, entropy                       ! input namelist
  real(8):: wc                                                                  ! wall clock
!-------------------------------------------------------------------------------
! Initialize
!-------------------------------------------------------------------------------
  !do_trace = .true.                                                             ! trace execution
  !do_debug = .true.                                                             ! debug printout
  call init_mpi                                                                 ! start mpi
  call print_id(id)                                                             ! identify version
  call init_io                                                                  ! setup input/output
  call init_mesh                                                                ! setup mesh
  call init_boundary                                                            ! ghost zones
  call init_timestep                                                            ! time integration
  call init_pde                                                                 ! PDE variables
!-------------------------------------------------------------------------------
! Setup initial conditions
!-------------------------------------------------------------------------------
  rewind (input); read (input,example); if (master) write(*,example)            ! read namelist input
  call allocate_mhd (v)                                                         ! allocate variables
  call init_velocity  (v%p)
  do iz=0,m(3)%n-1; do iy=0,m(2)%n-1; do ix=0,m(1)%n-1                          ! loop over 3-D
      v%d(ix,iy,iz) = density                                                   ! density
      v%s(ix,iy,iz) = density*entropy                                           ! entropy per unit vol
    v%p%x(ix,iy,iz) = density*v%p%x(ix,iy,iz)                                   ! momentum
    v%p%y(ix,iy,iz) = density*v%p%y(ix,iy,iz)
    v%p%z(ix,iy,iz) = density*v%p%z(ix,iy,iz)
    v%B%z(ix,iy,iz) = exp(-(m(1)%r(ix)**2+m(2)%r(iy)**2)/width**2)              ! Gaussian Bz(x,y), centered
  end do; end do; end do
  call boundary       (v%d)
  call boundary       (v%s)
  call boundary_vector(v%p)
  call boundary_vector(v%B)
!-------------------------------------------------------------------------------
! Setup output
!-------------------------------------------------------------------------------
  iout = 0                                                                      ! zero based snapshot
  call file_openw_mpi ('snapshot.dat', m%gn, m%n, m%offset)                     ! open data file
  call write_mhd(v)                                                             ! dump data
!-------------------------------------------------------------------------------
! Evolve the solution f in time
!-------------------------------------------------------------------------------
  wc = wallclock()                                                              ! initial wall time
  do istep = 1, nstep                                                           ! time step loop
    call time_step (v)                                                          ! step forward
    if (modulo(istep,fout)==0) call write_mhd(v)                                ! dump data
  enddo
  wc = wallclock()-wc                                                           ! final wall time
  if (master) print'(a,f6.1,a,f6.1)', &                                         ! nano-seconds per step
    'time:',wc,' s,  ns/p:',wc*1e9/(nstep*product(real(m%n)))
!-------------------------------------------------------------------------------
! Clean up
!-------------------------------------------------------------------------------
  call file_close_mpi                                                           ! close snapshot file
  call deallocate_mhd (v)                                                       ! kill variables
  call end_mpi                                                                  ! end mpi
CONTAINS

!*******************************************************************************
SUBROUTINE write_mhd (v)
  implicit none
  type(mhd_t):: v
!...............................................................................
  call write_field (v%d, 'd')
  call write_field (v%s, 's')
  call write_field (v%p%x/v%d, 'Ux')
  call write_field (v%p%y/v%d, 'Uy')
  call write_field (v%p%z/v%d, 'Uz')
  call write_field (v%B%x, 'Bx')
  call write_field (v%B%y, 'By')
  call write_field (v%B%z, 'Bz')
END SUBROUTINE write_mhd

!*******************************************************************************
SUBROUTINE write_field (f, c)
  implicit none
  character(len=*):: c
  real, dimension(m(1)%lb:m(1)%ub, m(2)%lb:m(2)%ub, m(3)%lb:m(3)%ub):: f
  real, allocatable, dimension(:,:,:) :: interior                               ! without ghost zones
  real:: mn, av, mx
!...............................................................................
  allocate(interior(m(1)%n,m(2)%n,m(3)%n))                                      ! buffer for I/O
  interior = f(0:m(1)%n-1,0:m(2)%n-1,0:m(3)%n-1)                                ! skip ghost zones
  call file_write_mpi (interior, m%n, iout)                                     ! collective write
  mn = minval(interior);              call min_real_mpi(mn)                     ! global min
  av = sum(interior)/product(m%gn);   call sum_real_mpi(av)                     ! global average
  mx = maxval(interior);              call max_real_mpi(mx)                     ! global max
  if (master) print '(i5,3x,a,1p,e11.4,3x,a,3e11.3,2x,a)', &
    iout/8, 't =',t,' min, aver, max :', mn, av, mx, c                          ! printout
  iout = iout + 1                                                               ! increment snapshot
  deallocate (interior)
END SUBROUTINE write_field
END PROGRAM
