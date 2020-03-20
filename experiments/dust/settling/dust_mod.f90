!===============================================================================
!> Data type for dust-gas interaction:
!>
!> The dust evolution is simulated using a fairly large number of particles oer
!> cell, e.g. of order 10-30 per size bin, or of the order of 100 in total
!>
!> The particles are moved using a particle solver, which operates on lists of
!> particles, which have attributes suitable to describe dust properties.  The
!> most important of these are
!>
!>    particle_t%w    : a weight that expresses number density per particle
!>    particle_t%d    : the particle diameter (size)
!>    particle_t%m    : the particle_mass
!>    particle_t%q    : an index triple that giving the integer cell coordinates
!>    particle_t%r    : a fractional index giving the position in the cell
!>
!> The main procedures that contribute to evolving the dust are
!>
!>    dust_t%gather   : gathers distributions of particles inside each cell
!>    dust_t%evolve   : decides how to evolve the dust distribution
!>    dust_t%scatter  : updates the particle properties accordingly
!>
!> The dust evolution does not necessarily have to be done after each update of
!> particle positions; it could be done at regular intervals, semi-independently.
!>
!> The dust module could be called as a sub-process, or as a separate task.
!> Initially it is most likely the simplest and most transparent organization
!> to have the experiment_t level be the task, from the point of view of the
!> framework, and let it call the various sub-solvers at various points in time
!===============================================================================
MODULE dust_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE mesh_mod
  USE gpatch_mod
  USE particles_mod
  USE particle_solver_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Dust data type
  !-----------------------------------------------------------------------------
  type, public:: dust_t
    real(8):: time = 0.0                      ! dust age
    real(8):: dtime = 1d0                     ! cadence
    real:: min_size, max_size                 ! min and max size in code units
    real:: specific_mass                      ! mass per unit volume, code units
    real, allocatable:: m_bin(:)              ! bin mass
    real, allocatable:: s_bin(:)              ! particle size
    real, allocatable:: m(:,:,:,:)            ! mass per bin
    real, allocatable:: kernel(:,:)           ! k(i,j) = increase at j from i
    class(gpatch_t), pointer:: patch          ! pointer to the MHD patch
    type(particle_solver_t):: particle_solver ! has the particles
  contains
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure, private:: gather
    procedure, private:: evolve
    procedure, private:: scatter
    procedure, private:: make_kernel
  end type
  integer:: verbose = 0                       ! information level
  integer:: cadence                           ! particle steps per evolve step
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: specific_mass_cgs = 2.0              ! g/cm**3
  integer:: n_bin = 15                        ! 5 bins per power of ten
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self)
  class(dust_t):: self
  !.............................................................................
  integer:: iostat, bin
  logical:: first_time = .true.
  namelist /dust_parameters/ verbose, min_size_cgs, max_size_cgs, n_bin, &
    cadence
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%update')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_parameters, iostat=iostat)
    call io%namelist_warning ('dust_parameters', iostat)
    write (io%output, dust_parameters)
    if (io%master) write (io%namelist, dust_parameters)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call self%particle_solver%init                ! initialize particle solver
  self%particle_solver%patch => self%patch      ! make the patch known to it
  !-----------------------------------------------------------------------------
  ! Scaling to code units
  !-----------------------------------------------------------------------------
  self%min_size = min_size_cgs/scaling%l        ! scale to code units
  self%max_size = max_size_cgs/scaling%l
  self%specific_mass = specific_mass_cgs/(scaling%m/scaling%l**3)
  allocate (self%s_bin(n_bin), self%m_bin(n_bin))
  do bin=1,n_bin
    self%s_bin(bin) = exp(log(self%min_size) &
                + (bin-1.0)/(n_bin-1.0)*log(self%max_size/self%min_size))
    self%m_bin(bin) = self%specific_mass*self%s_bin(bin)**3
  end do
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state
!===============================================================================
SUBROUTINE pre_update (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%pre_update')
  call self%particle_solver%update
  self%dtime = cadence*self%particle_solver%dtime
  if (self%patch%time > self%time) then
    call self%gather
    call self%evolve
    call self%scatter
    self%time = self%time + self%dtime
  end if
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Do whatever should be done after the HD/MHD update
!===============================================================================
SUBROUTINE post_update (self)
  class(dust_t):: self
  !.............................................................................
  call trace%begin ('dust_t%post_update')
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Gather information about the current state of the dust and gas
!===============================================================================
SUBROUTINE gather (self)
  class(dust_t):: self
  !.............................................................................
  integer:: bin, it, i(3), n(3)
  type(particle_t), pointer:: pa
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%gather')
  !-----------------------------------------------------------------------------
  ! Compute histograms of particle number densities in each cell
  !-----------------------------------------------------------------------------
  n = self%patch%mesh%gn                        ! cells in patch
  allocate (self%m(n(1),n(2),n(3),n_bin))       ! work array for histograms
  allocate (self%kernel(n_bin,n_bin))           ! kernel expressing exchange
  self%m = 0.0                                  ! prepare to accumulate
  pa => self%particle_solver%particles%head     ! first particle
  do while (associated(pa))                     ! loop over particle
    bin = pa%bin                                ! particle bin
    it = pa%it                                  ! particle time slice
    i = pa%q(:,it)                              ! particle space index
    self%m(i(1),i(2),i(3),bin) = &              ! accumulate
    self%m(i(1),i(2),i(3),bin) + &              ! into self%m
    pa%w*self%m_bin(bin)                        ! particle mass
    pa => pa%next                               ! next particle
  end do
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE gather

!===============================================================================
!> Based on the current state of gas and dust, decide how to evolve it, by
!> applying a kernel to the distribution over bins and adding up the results.
!> The coefficient to apply is the ratio of the old to the new number density.
!===============================================================================
SUBROUTINE evolve (self)
  class(dust_t):: self
  !.............................................................................
  real:: m(n_bin)
  integer:: i1, i2, i3, i, j
  class(mesh_t), pointer:: mesh(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%evolve')
  if (self%particle_solver%particles%n > 0) then
    mesh  => self%patch%mesh
    do i3=mesh(3)%li,mesh(3)%ui
    do i2=mesh(2)%li,mesh(2)%ui
    do i1=mesh(1)%li,mesh(1)%ui
      call self%make_kernel 
      m(:) = self%m(i1,i2,i3,:)    
      do i=1,n_bin
        do j=i+1,n_bin
          m(i) = m(i) - self%kernel(i,j)*self%m(i1,i2,i3,j)
          m(j) = m(j) + self%kernel(i,j)*self%m(i1,i2,i3,j)     
        end do
      end do
      !-------------------------------------------------------------------------
      ! The final coefficient of change is the new mass divided by the old
      !-------------------------------------------------------------------------
      do i=1,n_bin
        self%m(i1,i2,i3,i) = m(i)/self%m(i1,i2,i3,i)
      end do
    end do
    end do
    end do
  end if
  call trace%end ()
END SUBROUTINE evolve

!===============================================================================
!> Construct a kernel k(i,j), where for j>i the k(i,j) element is how much bin i
!> gives or takes to/from slot j, while k(j,i) is how much bin j looses or gains
!> because of that.  Because the weights count number densities, the two coeffs
!> are related to each other by k(j,i)*m(j) = k(i,j)*m(i), where m(i) is the bin
!> mass.  This ensures mass conservation.
!===============================================================================
SUBROUTINE make_kernel (self)
  class(dust_t):: self
  !.............................................................................
  integer:: i, j
  !-----------------------------------------------------------------------------
  do i=1,n_bin
    do j=i+1,n_bin
      !self%kernel(i,j) = ... FIXME: add code
      self%kernel(i,j) = 0.0
    end do
  end do
END SUBROUTINE make_kernel

!===============================================================================
!> Update the dust state in the task, based on the desired evolution
!===============================================================================
SUBROUTINE scatter (self)
  class(dust_t):: self
  !.............................................................................
  integer:: bin, it, i(3), n(3)
  type(particle_t), pointer:: pa
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%scatter')
  pa => self%particle_solver%particles%head     ! first particle
  do while (associated(pa))                     ! loop over particle
    bin = pa%bin                                ! particle bin
    it = pa%it                                  ! particle time slice
    i = pa%q(:,it)                              ! particle space index
    pa%w = pa%w*self%m(i(1),i(2),i(2),bin)      ! renormalize
    pa => pa%next                               ! next particle
  end do
  deallocate (self%m, self%kernel)
  call trace%end ()
END SUBROUTINE scatter

END MODULE dust_mod
