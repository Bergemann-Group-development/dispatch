!===============================================================================
!> HD module for the Sedov blast wave exercise. This implements the hydrodynamic
!> equations in conservative form on a staggered mesh, where density and entropy
!> are centered at cell-centers, linear momenta are centered on cell faces (walls),
!> and momentum fluxes are centered on cell edges.
!===============================================================================
MODULE initial_mod
  USE mpi_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public:: initial_t
    procedure(advection), pointer:: condition=>null()
    logical:: mhd
    real(8):: time
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: cleanup
    procedure:: conditions
    procedure:: sedov
    procedure:: advection
  end type
  !-----------------------------------------------------------------------------
  type(initial_t), public:: initial
  real, dimension(:,:,:,:), pointer, save:: p
  real, dimension(:,:,:),   pointer, save:: d, s
CONTAINS

!===============================================================================
!> Generic initial conditions
!===============================================================================
SUBROUTINE init (self, solver)
  class(initial_t):: self
  character(len=64):: solver
  !.............................................................................
  self%solver = solver
END SUBROUTINE init

!===============================================================================
!> Generic initial conditions
!===============================================================================
SUBROUTINE cleanup (self)
  class(initial_t):: self
  !.............................................................................
END SUBROUTINE cleanup

!===============================================================================
!> Generic initial conditions
!===============================================================================
SUBROUTINE conditions (self, mesh, mem, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  integer:: ix, iostat
  logical, save:: first_time=.true.
  !.............................................................................
  character(len=16), save:: type='sedov'
  namelist /IC_params/ type
  !-----------------------------------------------------------------------------
  ! Read and print the namelist only once
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_t%conditions')
  call io%lock%set
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, initial_params,iostat=iostat)
    if (io%master) write (*,initial_params)
  end if
  call io%lock%unset
  d => mem(:,:,:,idx%d)
  s => mem(:,:,:,idx%s)
  p => mem(:,:,:,idx%px:idx%pz)
  select case (type)
  case ('sedov')
    call self%sedov (mesh, mem)
   case ('advection')
    call self%advection (mesh, mem)
  case default
    print *, 'UNKNOWN INITIAL CONDITION TYPE'
    call mpi%abort ('initial_mod')
  end select
  call trace%end()
END SUBROUTINE conditions

!===============================================================================
!> Sedov blast wave initial conditions
!===============================================================================
SUBROUTINE sedov (self, mesh, mem, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  real(8), dimension(:), pointer:: rx, ry, rz
  integer:: ix, iy, iz, iostat
  real:: x, y, z, f
  logical, save:: first_time=.true.
  !.............................................................................
  real, save:: w_gauss=3.0, p_gauss=4.0, amplitude_s=3.0, amplitude_d=3., gamma=5./3.
  namelist /sedov_params/ w_gauss, p_gauss, amplitude_s, amplitude_d, gamma
  !-----------------------------------------------------------------------------
  ! Read and print the namelist only once
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_t%conditions')
  call io%lock%set
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, sedov_params,iostat=iostat)
    if (io%master) write (*,sedov_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat in entropy if at time=0
  !-----------------------------------------------------------------------------
  rx => mesh(1)%r
  ry => mesh(2)%r
  rz => mesh(3)%r
  do iz=1,size(mem,3)
    z = mesh(3)%p + rz(iz)
    do iy=1,size(mem,2)
      y = mesh(2)%p + ry(iy)
      do ix=1,size(mem,1)
        !   patch position + mesh offset        - gaussian position
        x = mesh(1)%p + rx(ix)
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
  call trace%end()
END SUBROUTINE sedov

!===============================================================================
!> Map variables to the current memory slot and initialize the deriv operators.
!> mhd_t%init is called from timestep_t%init, again called from task_list_t%init.
!===============================================================================
SUBROUTINE advection (self, mesh, mem, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:,:), pointer:: mem
  !.............................................................................
  real(8), dimension(:), pointer:: rx, ry, rz
  integer:: ix, iy, iz, iostat
  real:: x, y, z, f
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  real, save:: velocity(2) = [0.2,0.1]
  real, save:: w_gauss=3.0, p_gauss=4.0, amplitude_s=1.0, amplitude_d=1., gamma=5./3.
  namelist /advection_params/velocity, w_gauss, p_gauss,  amplitude_s, &
                                       amplitude_d, gamma
  ! Read and print the namelist only once
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_t%conditions')
  call io%lock%set
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, advection_params,iostat=iostat)
    if (io%master) write (*,advection_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Constant advection speed
  !-----------------------------------------------------------------------------
  rx => mesh(1)%r
  ry => mesh(2)%r
  rz => mesh(3)%r
  do iz=1,size(mem,3)
    z = mesh(3)%p + rz(iz)
    do iy=1,size(mem,2)
      y = mesh(2)%p + ry(iy)
      do ix=1,size(mem,1)
        !   patch position + mesh offset - gaussian offset
        x = mesh(1)%p + rx(ix)
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
  call trace%end()
END SUBROUTINE advection

END MODULE initial_mod
