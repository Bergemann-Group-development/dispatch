!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE patch_mod
  USE eos_mod
  !
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real::   theta=10.0
    real::   mpoly=1.0
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t)     :: self
  class(patch_t),target          :: patch
  !.............................................................................
  call self%standard_boundaries_t%init (patch)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=var_k), dimension(:,:,:,:), pointer :: mem
  real(kind=var_k), dimension(:,:,:), pointer:: eth
  real(kind=var_k), dimension(:,:,:), pointer:: ux, uy, uz
  real(kind=var_k), dimension(:,:,:), pointer:: d0, eth0
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  integer:: ix, iy, iz, n(3)
  real(kind=var_k):: pg, rhov2, d, e
  real(8):: z, zz
  ! ............................................................................
   call trace_begin('boundaries_t%primary')
  m => self%mesh
  ! ----------------------------------------------------------------------------
  ! Allocate scratches
  ! ----------------------------------------------------------------------------
  n = m%gn
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (d0(n(1),n(2),n(3)))
  allocate (eth(n(1),n(2),n(3)))
  allocate (eth0(n(1),n(2),n(3)))
  !
  ! ----------------------------------------------------------------------------
  ! Fill in scracthes with values
  ! ----------------------------------------------------------------------------
  ux = mem(:,:,:,idx%px) / mem(:,:,:,idx%d)
  uy = mem(:,:,:,idx%py) / mem(:,:,:,idx%d)
  uz = mem(:,:,:,idx%pz) / mem(:,:,:,idx%d)
  ! ----------------------------------------------------------------------------
  ! Get thermal energy
  ! ----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb, m(1)%ub
        rhov2 = (ux(ix,iy,iz)**2 + &
                 uy(ix,iy,iz)**2 + &
                 uz(ix,iy,iz)**2)*mem(ix,iy,iz,idx%d) 
        eth(ix,iy,iz) = mem(ix,iy,iz,idx%e) - 0.5*rhov2
      end do
    end do
  end do
  ! ----------------------------------------------------------------------------
  ! Estimate the fluctuations for density and thermal energy
  ! ----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    z = max(m(3)%p + m(3)%r(iz),0.0)
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb, m(1)%ub
        d0(ix,iy,iz)   = mem(ix,iy,iz,idx%d) - (1.0+self%theta*z)**self%mpoly
      end do
    end do
  end do
  !
  do iz=m(3)%lb,m(3)%ub
    z = max(m(3)%p + m(3)%r(iz),0.0)
    pg = (1.0+self%theta*z)**(self%mpoly+1.0)
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb, m(1)%ub
        eth0(ix,iy,iz) = eth(ix,iy,iz) - pg/(eos%gamma-1.0)
      end do
    end do
  end do
  ! ----------------------------------------------------------------------------
  ! Both top and bottom boundary are identical and are an "inpenetrable wall"
  ! Boundary conditions should be:
  ! 1) antisymmetric for vertical velocity
  ! 2) symmetric for horizontal velocities
  ! 3) symmetric for logarithmic fluctuations of density and energy
  ! ----------------------------------------------------------------------------
  call self%symmetric(ux,idx%px)
  call self%symmetric(uy,idx%py)
  call self%antisymmetric(uz,idx%pz)
  call self%symmetric(d0,idx%d)
  call self%symmetric(eth0,idx%e)
  ! ----------------------------------------------------------------------------
  ! Put the quantities back
  ! ----------------------------------------------------------------------------
  if (m(3)%lower_boundary) then
    do iz=m(3)%lo,m(3)%lb,-1
      z = max(m(3)%p + m(3)%r(iz),0.0)
      !z=1.0-zz
      pg = (1.0+self%theta*z)**(self%mpoly+1.0)
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb, m(1)%ub
          d = (1.0+self%theta*z)**self%mpoly + d0(ix,iy,iz)
          e = pg/(eos%gamma-1.0) + eth0(ix,iy,iz)
          mem(ix,iy,iz,idx%d) = merge(d,(1.0+self%theta*z)**self%mpoly, d>1e-10)
          mem(ix,iy,iz,idx%e) = merge(e,pg/(eos%gamma-1.0),e >1e-10)
          mem(ix,iy,iz,idx%px) = ux(ix,iy,iz) * mem(ix,iy,iz,idx%d)
          mem(ix,iy,iz,idx%py) = uy(ix,iy,iz) * mem(ix,iy,iz,idx%d)
          mem(ix,iy,iz,idx%pz) = uz(ix,iy,iz) * mem(ix,iy,iz,idx%d)

          rhov2 = (ux(ix,iy,iz)**2 + &
                   uy(ix,iy,iz)**2 + &
                   uz(ix,iy,iz)**2)*mem(ix,iy,iz,idx%d) 
          mem(ix,iy,iz,idx%e) = mem(ix,iy,iz,idx%e) + 0.5*rhov2
        end do
      end do
    end do
  end if
  !
  if (m(3)%upper_boundary) then
    do iz=m(3)%uo,m(3)%ub
      z = max(m(3)%p + m(3)%r(iz),0.0)
      !z=1.0-zz
      pg = (1.0+self%theta*z)**(self%mpoly+1.0)
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb, m(1)%ub
          d = (1.0+self%theta*z)**self%mpoly + d0(ix,iy,iz)
          e = pg/(eos%gamma-1.0) + eth0(ix,iy,iz)
          mem(ix,iy,iz,idx%d) = merge(d,(1.0+self%theta*z)**self%mpoly, d>1e-10)
          mem(ix,iy,iz,idx%e) = merge(e,pg/(eos%gamma-1.0),e >1e-10)
          mem(ix,iy,iz,idx%px) = ux(ix,iy,iz) * mem(ix,iy,iz,idx%d)
          mem(ix,iy,iz,idx%py) = uy(ix,iy,iz) * mem(ix,iy,iz,idx%d)
          mem(ix,iy,iz,idx%pz) = uz(ix,iy,iz) * mem(ix,iy,iz,idx%d)

          rhov2 = (ux(ix,iy,iz)**2 + &
                   uy(ix,iy,iz)**2 + &
                   uz(ix,iy,iz)**2)*mem(ix,iy,iz,idx%d) 
          mem(ix,iy,iz,idx%e) = mem(ix,iy,iz,idx%e) + 0.5*rhov2
        end do
      end do
    end do
  end if
  if (any(mem(:,:,:,idx%d)<=0)) print *, 'low density', minval(mem(:,:,:,idx%d)), minloc(mem(:,:,:,idx%d))
  if (any(mem(:,:,:,idx%e)<=0)) then
    n = minloc(mem(:,:,:,idx%e))
    z = m(3)%p + m(3)%r(1)
      !z=1.0-zz
      pg = (1.0+self%theta*z)**(self%mpoly+1.0)
    print *, 'low energy', minval(mem(:,:,:,idx%e)), minloc(mem(:,:,:,idx%e)), eth0(n(1),n(2),n(3))
    print *, pg/(eos%gamma-1.0), z
  end if
  deallocate(ux,uy,uz,d0, eth, eth0)
  call trace%end()
END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)

END SUBROUTINE derivative


END MODULE boundaries_mod
