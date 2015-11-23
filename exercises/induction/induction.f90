!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the induction problem in 3-D, for now with the veloocity field constant
!> in space
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE induction_mod
  USE io_mod                                            ! input/output module
  USE mpi_mod                                           ! MPI dimensions
  USE mpi_coords_mod                                    ! MPI dimensions
  USE mpi_reduce_mod                                    ! max, min, average
  USE mpi_send_mod                                      ! send / recv
  USE mpi_par_io_mod                                    ! MPI parallel I/O
  implicit none
  private
  character(len=mch), save:: id= &
    'main.f90 $Id$'
  type, public:: induction_t
    real, dimension(3):: position, u, size, ds
    real:: time, dtime, end_time=1.0, out_time=0.1, out_next
    real:: courant, w
    logical:: do_upwind
    integer:: n(3), gn(3), iout=0, ns=0
    real, dimension(:,:,:), pointer:: Bx, Ex, dBxdt
    real, dimension(:,:,:), pointer:: By, Ey, dBydt
    real, dimension(:,:,:), pointer:: Bz, Ez, dBzdt
  contains
    procedure:: init
    procedure:: gaussian
    procedure:: run
    procedure:: output
    procedure:: print
    procedure:: finish
  end type
CONTAINS

!> =============================================================================
!> Initialize the problem
!> =============================================================================
SUBROUTINE init (self)
  class(induction_t):: self
  real, dimension(3):: u=1.0, center
  real:: width=10., size=1.0, courant=0.2, end_time=1.0, out_time=0.1
  integer:: n=100, dims(3), offset(3)=0, coords(3)
  character(len=64):: file='induction.dat'
  logical:: do_upwind=.true.
  namelist /induction_params/ n, u, size, width, courant, end_time, out_time, &
                              do_upwind
  !-----------------------------------------------------------------------------
  ! Initialize MPI and mesh
  !-----------------------------------------------------------------------------
  call mpi%init                                         ! start mpi
  call io%init                                          ! setup input/output
  call print_id(id)                                     ! identify version
  call mpi_send%init
  call mpi_par_io%init                                  ! MPI parallel I/O
  !-----------------------------------------------------------------------------
  ! Setup initial conditions
  !-----------------------------------------------------------------------------
  size = n                                              ! default ds=1
  call mpi%init                                         ! start MPI
  rewind (io%input); read (io%input,induction_params)   ! read namelist
  if (io%master) write(*,induction_params)              ! namelist output
  call mpi_coords%init (dims=[n,n,n])                   ! MPI geometry
  self%gn = n                                           ! global dimensions
  self%ds = size/n                                      ! cell size
  self%n  = n/mpi_coords%dims                           ! local dimensions
  coords  = mpi_coords%coords                           ! integer MPI coords
  center  = (mpi_coords%dims-1)*0.5                     ! center of ditto
  self%size = size/mpi_coords%dims                      ! size of sub-domains
  self%position = self%size*(coords-center)             ! real position
  self%u = u                                            ! velocity
  self%w = width                                        ! Gaussian width
  self%courant = courant                                ! Courant number
  self%end_time = end_time                              ! run end time
  self%out_time = out_time                              ! file output interval
  self%out_next = out_time                              ! next file output
  self%do_upwind = do_upwind                            ! choose derivative type
  call mpi_reduce%init                                  ! reduction operators
  call mpi_par_io%init                                  ! parallel I/O
  call mpi_par_io%openw (file, self%gn, self%n, offset) ! open output file
  allocate (self%Bx(self%n(1),self%n(2),self%n(3)), &   ! allocate magnetic field
            self%By(self%n(1),self%n(2),self%n(3)), &   !
            self%Bz(self%n(1),self%n(2),self%n(3)), &   !
            self%Ex(self%n(1),self%n(2),self%n(3)), &   ! allocate electric field
            self%Ey(self%n(1),self%n(2),self%n(3)), &   !
            self%Ez(self%n(1),self%n(2),self%n(3)), &   !
         self%dBxdt(self%n(1),self%n(2),self%n(3)), &   ! allocate B time derivative
         self%dBydt(self%n(1),self%n(2),self%n(3)), &   !
         self%dBzdt(self%n(1),self%n(2),self%n(3)))     !
  call self%gaussian (self%bz)
  call self%output (self%bz)
END SUBROUTINE

!> =============================================================================
!> Gaussian centered on self%position, for initial conditions and comparisons.
!> The width of the Gaussian is in grid units, so the coordinates are also.
!> =============================================================================
SUBROUTINE gaussian (self, f)
  class(induction_t):: self
  real, dimension(:,:,:), pointer:: f
  integer:: i, j, k
  real:: x, y, z
  !.............................................................................
  do k=1,self%n(3)
  do j=1,self%n(2)
  do i=1,self%n(1)
    x = (i-(self%n(1)+1)*0.5) + self%position(1)/self%ds(1)
    y = (j-(self%n(2)+1)*0.5) + self%position(2)/self%ds(2)
    z = (k-(self%n(3)+1)*0.5) + self%position(3)/self%ds(3)
    f(i,j,k) = exp(-(x**2+y**2+z**2)/self%w**2)
  end do
  end do
  end do
END SUBROUTINE

!> =============================================================================
!> Evolve the solution in time.  The associate / end associate is effectively a
!> region with aliasing, which makes expressions there more easily readable.
!> =============================================================================
SUBROUTINE run (self)
  USE derivatives_mod
  USE central_mod
  USE upwind_mod
  class(induction_t):: self
  integer:: istep
  type(central_t), target:: c
  type(upwind_t), target:: u
  class(derivatives_t), pointer:: d
  !.............................................................................
  if (self%do_upwind) then
    d => u
  else
    d => c
  end if
  associate (Bx=>self%bx, Ex=>self%ex, dbxdt=>self%dbxdt, &
             By=>self%by, Ey=>self%ey, dbydt=>self%dbydt, &
             Bz=>self%bz, Ez=>self%ez, dbzdt=>self%dbzdt, &
             u=> self%u, ds=>self%ds, n=>self%n, ns=>self%ns, &
             t=>self%time, dt=>self%dtime, courant=>self%courant)
  dt = courant*minval(ds)/sqrt(sum(u**2))               ! dt = C ds/u
  do while (t < self%end_time)
    Ex = -u(2)*bz + u(3)*by                             ! E = -(U x B)
    Ey = -u(3)*bx + u(1)*bz
    Ez = -u(1)*by + u(2)*bx
    d%ds = ds
    d%n  = n
    dbxdt = -d%dy(Ez) + d%dz(Ey)                        ! dB/dt = -curl(E)
    dbydt = -d%dz(Ex) + d%dx(Ez)
    dbzdt = -d%dx(Ey) + d%dy(Ex)
    bx = bx - dt*dbxdt                                  ! B = B + dt*dB/dt
    by = by - dt*dbydt
    bz = bz - dt*dbzdt
    t  = t  + dt                                        ! increment time
    ns = ns + 1
    call self%output (bz)                               ! write?
  enddo
  end associate
END SUBROUTINE run

!> =============================================================================
!> File I/O
!> =============================================================================
SUBROUTINE output (self, f)
  implicit none
  class(induction_t):: self
  real f(self%n(1), self%n(2), self%n(3))
  real:: mn, av, mx
!...............................................................................
  if (self%time < self%out_next) return                 ! next output?
  self%out_next = self%out_next + self%out_time         ! cadence
  call mpi_par_io%write (self%bz, self%iout)            ! collective write
  mn = minval(self%bz)
  mx = maxval(self%bz)
  av = sum(self%bz)/product(self%gn)
  call mpi_reduce%min_real (mn)                         ! global min
  call mpi_reduce%max_real (mx)                         ! global min
  call mpi_reduce%sum_real (av)                         ! global sum
!  print*,self%gn, size(self%bz), av
  if (io%master) print '(i5,3x,a,f10.6,3x,a,1p,3e12.4)', self%iout, 't:', &
    self%time, 'min, aver, max:', mn, av, mx            ! printout
  self%iout = self%iout + 1                             ! increment snapshot
END SUBROUTINE output

!> =============================================================================
!> Print summary, with error measure
!> =============================================================================
SUBROUTINE print (self)
  class(induction_t):: self
  integer:: istep, omp_get_thread_num
  real, dimension(:,:,:), pointer:: tmp
  real:: sqsum(2)
  !.............................................................................
  self%position = self%position + self%time*self%u
  allocate (tmp(self%n(1), self%n(2), self%n(3)))
  call self%gaussian (tmp)
  sqsum(1) = sum((self%bz-tmp)**2)
  sqsum(2) = sum(tmp**2)
  deallocate (tmp)
  call mpi_reduce%sum_real (sqsum(1))
  call mpi_reduce%sum_real (sqsum(2))
  !call mpi_reduce%sum_reals (2, sqsum)
  if (io%master) &
    print'(1x,2(a,i5,3x),a,f6.3,3x,a,1p,2e10.2)', &
    'thread:', omp_get_thread_num(), &
    'n_step:', self%ns, &
    'Courant:', self%courant, &
    'rms:', sqrt(sqsum(1)/sqsum(2))
END SUBROUTINE

!> =============================================================================
!> Clean up
!> =============================================================================
SUBROUTINE finish (self)
  class(induction_t):: self
  !.............................................................................
  call self%print
  deallocate (self%Bx, self%Ex, self%dBxdt)
  deallocate (self%By, self%Ey, self%dBydt)
  deallocate (self%Bz, self%Ez, self%dBzdt)
  call mpi%end
END SUBROUTINE finish
END MODULE induction_mod

!> =============================================================================
!> =============================================================================
PROGRAM induction_prog
  USE induction_mod
  implicit none
  type(induction_t):: induction
  call induction%init
  call induction%run
  call induction%finish
END PROGRAM induction_prog
