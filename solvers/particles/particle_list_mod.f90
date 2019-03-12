!===============================================================================
!> Particle list, extending a doubly-linked list.  Each particle maintains 
!> arrays with previous positions, velocities, and times, which are brought
!> along if/when it changes owner patch or rank.  The memory footprint of a
!> particle is 12 words for position, 12 words for velocity, and 8 words for
!> time, all together 32 words = 128 bytes (plus a few words for id and weight).
!> This could be reduced to half, by keeping only two previous positions and
!> velocities in the data type.
!===============================================================================
MODULE particle_list_mod
  USE io_unit_mod
  USE dll_mod
  USE omp_timer_mod
  USE particle_mod
  USE patch_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: particle_list_t
    integer:: verbose=0
    integer:: nt=4
  contains
    procedure:: init
    procedure:: force
    procedure:: force_field
    procedure:: test
    procedure:: print
  end type
  !-----------------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE init (self, name)
  class (particle_list_t):: self
  character(len=*), optional:: name
  !-----------------------------------------------------------------------------
  call self%dll_t%init (name)
  call particle%init
END SUBROUTINE init

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE append (self, new)
  class (particle_list_t):: self
  class (particle_t), pointer:: new
  class (dll_node_t), pointer:: new_dll
  !-----------------------------------------------------------------------------
  new_dll => new
  call self%dll_t%append (new_dll)
END SUBROUTINE append

!===============================================================================
!> Compute force from other particles by direct summation
!===============================================================================
FUNCTION force (self, p)
  class(particle_list_t):: self
  class(particle_t), pointer:: p
  real(8):: force(3)
  !.............................................................................
  class(dll_node_t), pointer:: o
  integer:: i, i0, i1, j
  real(8):: rp(3), pt, qt, r(3), r1, r2, r3
  !-----------------------------------------------------------------------------
  force = 0.0_8
  rp = p%r(:,p%it)
  o => self%head
  do while (associated(o))
    if (.not.associated (o, p)) then
      select type (o)
      class is (particle_t)
      do i=2,self%nt-1
        i1 = o%iit(i)
        if (o%t(i1) > p%time) exit
      end do
      i0 = o%iit(i-1)
      pt = (p%time-o%t(i0))/max(o%t(i1)-o%t(i0),tiny(1d0))
      qt = 1d0-pt
      r2 = tiny(1d0)
      do j=1,3
        r(j) = (qt*o%r(j,i0) + pt*o%r(j,i1)) - rp(j)
        r2 = r2 + r(j)**2
      end do
      r1 = 1.0/sqrt(r2)
      r3 = r1**3
      force = force - o%mass*r*r3
      end select
    end if
    o => o%next
  end do
END FUNCTION force

!===============================================================================
!===============================================================================
SUBROUTINE force_field (self, patch)
  class(particle_list_t):: self
  class(patch_t):: patch
  !.............................................................................
  class(dll_node_t), pointer:: p
  real(8), allocatable:: rp(:,:), m(:)
  real(8):: force(3), pt, qt, rc(3), r(3), r1, r3
  integer:: i, i0, i1, ip, ix, iy, iz
  !-----------------------------------------------------------------------------
  ! Cache all particle position into array, in order to be able to loop over
  ! particles
  !-----------------------------------------------------------------------------
  allocate (rp(self%n,3), m(self%n))
  ip = 0
  p => self%head
  do while (associated(p))
    ip = ip+1
    select type (p)
    class is (particle_t)
    do i=2,self%nt-1
      i0 = p%iit(i-1)
      i1 = p%iit(i)
      if (p%t(i1) > patch%time) exit
    end do
    pt = (patch%time - p%t(i0))/max(p%t(i1)-p%t(i0),tiny(1d0))
    qt = 1d0-pt
    rp(:,ip) = qt*p%r(:,i0) + pt*p%r(:,i1)
    m(ip) = p%mass
    end select
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! For each cell, loop over particles and accumulate force in double precision
  !-----------------------------------------------------------------------------
  associate (m1=>patch%mesh(1), m2=>patch%mesh(2), m3=>patch%mesh(3))
  do iz=m3%li,m3%ui
    rc(3) = m3%p + m3%r(iz)
    do iy=m2%li,m2%ui
      rc(2) = m2%p + m2%r(iy)
      do ix=m1%li,m1%ui
        rc(1) = m1%p + m1%r(ix)
        force = 0d0
        do ip=1,self%n
          r = rp(:,ip) - rc
          r1 = 1d0/sqrt(r(1)**2 + r(2)**2 + r(3)**2 + tiny(1d0))
          r3 = r1*r1*r1
          force = force - m(ip)*r*r3
        end do
        patch%force_per_unit_mass(ix,iy,iz,:) = &
        patch%force_per_unit_mass(ix,iy,iz,:) + force
      end do
    end do
  end do
  end associate
  deallocate (rp)
END SUBROUTINE force_field

!===============================================================================
!> Test particle list functionality
!===============================================================================
SUBROUTINE test (self)
  class(particle_list_t):: self
  integer:: id, n
  class(particle_t), pointer:: p
  class(dll_node_t), pointer :: node, next, new
  real(8)::used
  logical:: flip
  !-----------------------------------------------------------------------------
  call self%init
  !-----------------------------------------------------------------------------
  ! Make a list with 100000 particles
  !-----------------------------------------------------------------------------
  n = 100000
  used = wallclock()
  do id=1,n
    allocate (p)
    call p%init
    p%r  = 0.5_8
    p%ds = 1.0_8
    p%v  = 1.0_8
    node => p
    call self%append (node)
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Remove every 2nd particle
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    select type (node)
    type is (particle_t)
    p => node
    end select
    if (flip) then
      call self%remove (node)
      deallocate (node)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Add one new particle for every 2nd existing one
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    allocate (p)
    new => p
    if (flip) then
      call self%insert_before (node, new)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
END SUBROUTINE test

!===============================================================================
!> Print a particle list
!===============================================================================
SUBROUTINE print (self)
  class(particle_list_t):: self
  class(dll_node_t), pointer:: p
  p => self%head
  do while (associated(p))
    select type (p)
    class is (particle_t)
    print *, p%id, p%r(:,p%it)
    end select
    p => p%next
  end do
END SUBROUTINE print

END MODULE particle_list_mod
