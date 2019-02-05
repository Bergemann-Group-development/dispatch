!> =============================================================================
!> Trace particle module.  Two classes of trace particles will be supported
!> (this is not yet fully implemented):
!>
!> 1) The bulk of the particles are used to accurate trace the motion of mass,
!> with a fair number of particles per cell, constantly renormalized to accurately
!> represent the mass in the patch.
!>
!> 2) A small fraction of trace particles, intended for graphics, are flagged
!> with a bit that results in a higher cadence of saving to disk, with particles
!> moving only a few cells between snapshots.
!>
!> Disk files should be organized with one file per snapshot for the bulk
!> particles, while the graphics particles should be accessible in one file
!> per rank, since it needs to be written too frequently for collective MPI
!> calls.  When reading and plotting, one can open these rank files, and make
!> a selection, based on for example their first location.  Subsequently, as
!> particles transition to other ranks, this has to be handled by the reader.
!>
!> One example:  For superonic turbulence, one may choose to save only two
!> snapshots per dynamic time L/M, while that would correspond to at least N/2
!> cell transitions. Hence, the number of particles needs to be of order N
!> smaller in the graphics files, if disks files should have about the same
!> size. In principle, only positions and particle IDs need to be save in the
!> bulk snapshots, which means of order 16 bytes per particle, and one can then
!> allow a few particles per cell, w/o a major increase in data volume.
!> 
!> Conversely, one can allow of order N times fewer particles in the graphics
!> files, or about one particle per 1000 cells, or one particle per radius 10.
!> One could choose to make the graphics particles have a constant fraction of
!> the mass, especially if one ensures that their paths accurately reflect the
!> the path of the mass. For such a small fraction it does not make sense to
!> even consider mass conservation, and one should instead extra precision in
!> their motion.
!>
!> The particle ID and a few flags may be encoded in a 64 bit integer, with the
!> upper 4 bits reserved for flags, being masked out when considering IDs.
!>
!> In the case of the bulk particles, we are interested in two partially
!> conflicting requirements: 1) they should be renomalized frequently, to keep
!> the weight of each particle from becoming too extreme (small or large)
!> compared with the amount of mass in the cell (or in the patch). 2) particle
!> heritage should be accessible, so one can identify exactly where a particle
!> comes from.  This could be achieved by reserving the lowermost No bits for
!> their "original" ID.
!> =============================================================================
MODULE trace_particles_mod
  USE io_mod
  USE trace_mod
  USE dll_mod
  USE particle_list_mod
  USE omp_timer_mod
  USE random_mod
  USE mesh_mod
  implicit none
  private
  integer, parameter:: nt=4
  type, public, extends(particle_list_t):: trace_particles_t
    real:: per_cell
    real:: lf(3), uf(3)
    integer:: id, it, nt, time_order, n_trade, verbose
    integer:: offset(3)
    !real(8), allocatable:: r(:,:,:), v(:,:,:), w(:)
    real(8), allocatable:: wt(:)
    real(8):: position(3)
    logical:: on
    type(particle_list_t):: export, import
  contains
    procedure:: init_nt
    procedure:: add
    procedure:: update
    procedure:: interpolate
    procedure:: trade
  end type
  type(random_t):: random
CONTAINS

!> =============================================================================
!> Initialize trace particles
!> =============================================================================
SUBROUTINE init_nt (self, nt)
  class(trace_particles_t):: self
  integer:: nt, iostat, np
  real, save:: per_cell=0.0
  integer, save:: time_order=3, verbose=0
  logical, save:: on=.false.
  logical, save:: first_time=.true.
  namelist /trace_particle_params/ on, verbose, per_cell, time_order
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, trace_particle_params, iostat=iostat)
    if (io%master) &
      write (io%output, trace_particle_params)
    call random%init
  end if
  !$omp end critical (input_cr)
  self%on = on
  self%nt = nt
  self%per_cell = per_cell
  self%time_order = time_order
  self%verbose = verbose
  if (.not.on) return
  call self%particle_list_t%init
  !-----------------------------------------------------------------------------
  ! Initialize Adams-Bashforth integration; the concept of symplectic integration
  ! does not apply when particles are just folling the gas, so there is no point
  ! in doing kick-drift-kick for trace particles.  One could possibly improve
  ! accuracy by using the fact that one has both location and velocity history,
  ! but a higher order cell velocity interpolation may be more important.
  !-----------------------------------------------------------------------------
  time_order = min(time_order,nt)
  allocate (self%wt(time_order))
  select case (time_order)
  case (2)
    self%wt(1:2) = [1.5, -0.5]
  case (3)
    self%wt(1:3) = [23./12., -16./12., 5./12.]
  case default
  end select
END SUBROUTINE init_nt

!> =============================================================================
!> Add trace particles, in proportion to density, and with velocity histories
!> intialized to a single time.
!> =============================================================================
SUBROUTINE add (self, mesh, t, dt, d, v)
  class(trace_particles_t):: self
  class(mesh_t), pointer:: mesh(:)
  real:: d(:,:,:), v(:,:,:,:)
  real(8):: t, dt
  !.............................................................................
  integer:: ip, np, ix, iy, iz
  class(particle_t), pointer:: p
  class(dll_node_t), pointer:: node
  !-----------------------------------------------------------------------------
  ! Add particles in proportion to the density, relative to the average
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('trace_particles_t%add')
  call self%particle_list_t%init ('trace')
  self%id = mesh(1)%id
  if (self%per_cell > 0.0) then
    np = ifix(self%per_cell)
    do iz=mesh(3)%li,mesh(3)%ui
    do iy=mesh(2)%li,mesh(2)%ui
    do ix=mesh(3)%li,mesh(1)%ui
      do ip=1,ifix(self%per_cell)
       call add1
      end do
      if (random%ran1() < (self%per_cell-np)) then
       call add1
      end if
    end do
    end do
    end do
  end if
  call trace%end
contains
  !-----------------------------------------------------------------------------
  real function ran()
    ran = random%ran1() - 0.5
  end function
  !-----------------------------------------------------------------------------
  subroutine add1
    integer:: it
    allocate (p)
    call p%init
    do it=1,nt
      p%r(:,it) = [ix+ran(),iy+ran(),iz+ran()]
      p%dr(:,it) = v(ix,iy,iz,:)*dt
      p%t(it) = t
    end do
    node => p
    call self%append (node)
  end subroutine add1
  !-----------------------------------------------------------------------------
END SUBROUTINE add

!> =============================================================================
!> Update particle positions.  To optimize, particle positions are already
!> normalized to floating point array coordinates.  This allows using single
!> precision, also on disk, where the coordinates just need to be complemented
!> with the double precision patch position.
!> =============================================================================
SUBROUTINE update (self, mesh, time, dtime, d, v)
  class(trace_particles_t):: self
  class(mesh_t), pointer:: mesh(:)
  real(8):: time, dtime, dt
  real:: d(:,:,:), v(:,:,:,:)
  !.............................................................................
  integer:: it, o1, o2, o3, i ,j, nrem
  class(dll_node_t), pointer:: node, next
  class(particle_t), pointer:: p
  class(*), pointer:: car
  real, pointer:: r(:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('trace_particles_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Import any particles exported by neighbors
  !-----------------------------------------------------------------------------
  if (self%import%n > 0) then
    node => self%import%head
    do while (associated(node))
      next => node%next
      call self%import%remove (node)
      call self%append (node)
      node => next
    end do
    if (self%verbose > 0) &
      write (io%output,*) self%id, ' has', self%n, ' particles'
  end if
  self%it = o1
  !-----------------------------------------------------------------------------
  ! Direct, non-vectorized update
  !-----------------------------------------------------------------------------
  self%lf = mesh%li - 0.5
  self%uf = mesh%ui + 0.5
  call self%export%init ('export')
  call self%import%init ('import')
  call self%interpolate (v)
  node => self%head
  nrem = 0
  do while (associated(node))
    select type (node)
    class is (particle_t)
    p => node
    end select
    it = 1+mod(p%it,nt)
    o1 = 1+modulo(it-2,nt)
    o2 = 1+modulo(it-3,nt)
    dt = min(dtime,time-p%t(o1))
    p%dr(:,o1) = p%dr(:,o1)*dt
    if (self%time_order == 3) then
      o3 = 1+modulo(it-4,nt)
      p%r(:,it) = &
      p%r(:,o1) + &
      p%dr(:,o1)*self%wt(1) + &
      p%dr(:,o2)*self%wt(2) + &
      p%dr(:,o3)*self%wt(3)
    else
      p%r(:,it) = &
      p%r(:,o1) + &
      p%dr(:,o1)*self%wt(1) + &
      p%dr(:,o2)*self%wt(2)
    end if
    p%it = it
    p%t(it) = time
    !---------------------------------------------------------------------------
    ! Particles outside the patch interior are appended to an export list and
    ! removed from the patch particle list
    !---------------------------------------------------------------------------
    r => p%r(:,it)
    if (any(r < self%lf .or. r > self%uf)) then
      call self%remove (node)
      call self%export%append (node)
      nrem = nrem+1
    end if
    node => node%next
  end do
  if (self%verbose > 0 .and. nrem>0) &
    write (io%output,*) self%id, ' has', self%n, ' particles', nrem
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Tri-linear interpolate in a velocity array.  The particle positions within
!> a patch are maintained in units of mesh points, both to simplify interpolation
!> and updates, and to improve precision, while staying with singe precision
!> in coordinates and velocities.
!> =============================================================================
SUBROUTINE interpolate (self, v)
  class(trace_particles_t):: self
  integer:: it
  real:: v(:,:,:,:)
  !.............................................................................
  integer:: sz(4), ip(3), i, j
  real:: p(3), q(3)
  class(dll_node_t), pointer:: node
  class(particle_t), pointer:: part
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Straight loops.  This could be optimized further, by using arrays instead
  ! of lists, and accpting "holes" in the arrays up to a certain fraction.  But
  ! this is already fast enough, given the disk space contsraints on the number
  ! of particles.
  !-----------------------------------------------------------------------------
  call trace%begin ('trace_particles_t%interpolate', itimer=itimer)
  sz = size(v)
  node => self%head
  do while (associated(node))
    select type (node)
    type is (particle_t)
    part => node
    it = part%it
    end select
    p = 1.0 + part%r(:,it)
    ip = floor(p)
    ip = max(1,min(ip,sz(1:3)))
    p = p-ip
    do j=1,3
      part%dr(j,it) = &
        q(3)*(q(2)*(q(1)*v(ip(1)  ,ip(2)  ,ip(3)  ,j)   &
                  + p(1)*v(ip(1)+1,ip(2)  ,ip(3)  ,j))  &
            + p(2)*(q(1)*v(ip(1)  ,ip(2)+1,ip(3)  ,j)   &
                  + p(1)*v(ip(1)+1,ip(2)+1,ip(3)  ,j))) &
      + p(3)*(q(2)*(q(1)*v(ip(1)  ,ip(2)  ,ip(3)+1,j)   &
                  + p(1)*v(ip(1)+1,ip(2)  ,ip(3)+1,j))  &
             +p(2)*(q(1)*v(ip(1)  ,ip(2)+1,ip(3)+1,j)   &
                  + p(1)*v(ip(1)+1,ip(2)+1,ip(3)+1,j)))
    end do
    node => node%next
  end do
  call trace%end (itimer)
END SUBROUTINE interpolate

!> =============================================================================
!> Trade particles between patches.  
!> =============================================================================
SUBROUTINE trade (self, nbor)
  class(trace_particles_t):: self
  class(trace_particles_t), pointer:: nbor
  class(dll_node_t), pointer:: node, next
  class(particle_t), pointer:: part
  real:: r(3)
  !-----------------------------------------------------------------------------
  self%n_trade = 0
  if (self%export%n == 0) return
  node => self%export%head
  do while (associated(node))
    next => node%next
    select type (node)
    class is (particle_t)
    part => node
    end select
    r = part%r(:,part%it) + nbor%offset
    !---------------------------------------------------------------------------
    ! Particles in the interior of an nbor node are removed from the export list
    ! of this patch, and put on the import list of that patch, with adjusted r
    !---------------------------------------------------------------------------
    if (all(r >= self%lf .and. r <= self%uf)) then
      if (self%verbose > 1) &
        write (io%output,*) nbor%id, self%id, part%id, part%r(:,part%it), nbor%offset
      part%r(:,part%it) = r
      call self%export%remove (node)
      call nbor%import%append (node)
      self%n_trade = self%n_trade + 1
    end if
    node => node%next
  end do
END SUBROUTINE trade

END MODULE trace_particles_mod
