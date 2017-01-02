!> =============================================================================
!> HD module for the Sedov blast wave exercise. This implements the hydrodynamic
!> equations in conservative form on a staggered mesh, where density and entropy
!> are centered at cell-centers, linear momenta are centered on cell faces (walls),
!> and momentum fluxes are centered on cell edges.
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5                     ! memory slots
  !-----------------------------------------------------------------------------
  type, public:: initial_t
    !procedure, pointer:: conditions => advection
    logical:: mhd
    real(8):: time
  contains
    procedure:: init
    procedure:: cleanup
    procedure:: conditions
    procedure:: sedov
    procedure:: advection
  end type
  !-----------------------------------------------------------------------------
  type(initial_t), public:: initial
  real, save:: velocity(2) = [0.2,0.1]
  real, save:: w_gauss=3.0, p_gauss=4.0, &
               amplitude_s=3.0, amplitude_d=3., gamma=5./3.
  logical, save:: first_time=.true.
  character(len=16), save:: type='sedov'
  real, dimension(:,:,:,:), pointer, save:: p
  real, dimension(:,:,:),   pointer, save:: d, s
CONTAINS

!> =============================================================================
!> Generic initial conditions
!> =============================================================================
SUBROUTINE init (self)
  class(initial_t):: self
  !.............................................................................
END SUBROUTINE init

!> =============================================================================
!> Generic initial conditions
!> =============================================================================
SUBROUTINE cleanup (self)
  class(initial_t):: self
  !.............................................................................
END SUBROUTINE cleanup

!> =============================================================================
!> Generic initial conditions
!> =============================================================================
SUBROUTINE conditions (self, mesh, mem)
  class(initial_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  integer:: ix, iostat
  namelist /initial_params/ type, velocity, w_gauss, p_gauss, amplitude_s, &
    amplitude_d, gamma
  !-----------------------------------------------------------------------------
  ! Read and print the namelist only once
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_t%conditions')
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, initial_params,iostat=iostat)
    if (io%master) write (*,initial_params)
  end if
  d => mem(:,:,:,i_d)
  s => mem(:,:,:,i_s)
  p => mem(:,:,:,i_px:i_pz)
  if (trim(type)=='sedov') then
    call self%sedov (mesh, mem)
  else if (trim(type)=='advection') then
    call self%advection (mesh, mem)
  end if
END SUBROUTINE conditions

!> =============================================================================
!> Sedov blast wave initial conditions
!> =============================================================================
SUBROUTINE sedov (self, mesh, mem)
  class(initial_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  integer:: ix, iy, iz
  real:: x, y, z, f
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy if at time=0
  !-----------------------------------------------------------------------------
  do iz=1,mesh(3)%ub
    z = mesh(3)%p + mesh(3)%r(iz)
    do iy=1,mesh(2)%ub
      y = mesh(2)%p + mesh(2)%r(iy)
      do ix=1,mesh(1)%ub
        !   patch position + mesh offset        - gaussian position
        x = mesh(1)%p + mesh(1)%r(ix)
        if (mesh(3)%ub > 1) then
          f = exp(-((x**2+y**2+z**2)/w_gauss**2)**(p_gauss/2.0))
        else
          f = exp(-((x**2+y**2     )/w_gauss**2)**(p_gauss/2.0))
          !f = cos(2.*pi*x*8.)
        end if
        ! ln = gaussian
        d(ix,iy,iz) = exp(amplitude_d*f)
        s(ix,iy,iz) = amplitude_s*f*d(ix,iy,iz)
        !if (ix==127.and.iy==127) print*,amplitude_s,w_gauss,p_gauss,f,s(ix,iy,iz)
        p(ix,iy,iz,:) = 0.0
      end do
    end do
  end do
  call trace%end
END SUBROUTINE sedov

!> =============================================================================
!> Map variables to the current memory slot and initialize the deriv operators.
!> pde_t%init is called from timestep_t%init, again called from task_list_t%init.
!> =============================================================================
SUBROUTINE advection (self, mesh, mem)
  class(initial_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  integer:: ix, iy, iz
  real:: x, y, z, f
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Constant advection speed
  !-----------------------------------------------------------------------------
  do iz=1,mesh(3)%ub
  do iy=1,mesh(2)%ub
  do ix=1,mesh(1)%ub
    !   patch position + mesh offset - gaussian offset
    x = mesh(1)%p + mesh(1)%r(ix)
    y = mesh(2)%p + mesh(2)%r(iy)
    z = mesh(3)%p + mesh(3)%r(iz)
    if (mesh(3)%ub > 1) then
      f = exp(-((x**2+y**2+z**2)/w_gauss**2)**(p_gauss/2.0))
    else
      f = exp(-((x**2+y**2     )/w_gauss**2)**(p_gauss/2.0))
    end if
    ! ln = gaussian
    d(ix,iy,iz) = exp(amplitude_d*f)
    s(ix,iy,iz) = log(1./d(ix,iy,iz)**gamma)/(gamma-1.0)
    p(ix,iy,iz,1) = velocity(1)
    p(ix,iy,iz,2) = velocity(2)
    p(ix,iy,iz,3) = 0.0
  end do
  end do
  end do
  call trace%end
END SUBROUTINE advection

END MODULE initial_mod
