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
  USE timer_mod
  USE index_mod
  USE math_mod
  USE patch_mod
  USE standard_boundaries_mod
  USE eos_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real(kind=KindScalarVar)       :: gravity = -2.74
    real(kind=KindScalarVar)       :: d_bot, e_bot, ee_top
    real(kind=KindScalarVar), dimension(:,:,:,:), pointer::ic_copy => null()
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
    procedure:: lower_y
    procedure:: upper_y
    procedure:: lower_z
    procedure:: upper_z
    procedure:: initial_conditions_copy
    procedure:: static_values
    procedure:: add_flow
    procedure:: smooth_reduce
  end type
  integer, save:: verbose=1
  logical, save:: smudge_boundary = .false.
  logical, save:: scale = .true.
  logical, save:: deriv_bnd = .false.
  real, save:: d_bot=0.0, e_bot=0.0, ee_top=0.0
  real, save:: bottom_bc_rate=0.01
  real, save:: damp_rate=0.01
  real, save:: r_sponge = 0.5
  real, save:: pb_fact = 1.0
  integer, save:: divb_iter = 4
  ! --- lower_damp parameters
  real, save:: beta_lim   = 1e+2
  real, save:: spot_lower = 0.0
  real, save:: pb_spot    = 0.0
  real, save:: t_bdry     = 3e+0
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t)              :: self
  class(patch_t), target:: patch
  integer:: id
  !.............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /boundary_params/ d_bot, e_bot, ee_top, bottom_bc_rate, verbose, &
           divb_iter, scale, t_bdry, r_sponge, damp_rate, deriv_bnd, pb_fact
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  !call io%lock%set('boundaries_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (io%master) write (*, boundary_params)
    !d_bot = d_bot/scaling%d
    !e_bot = e_bot/scaling%e
    first_time = .false.
    !self%d_bot = d_bot
  end if
  !call io%lock%unset('boundaries_t%init')
  call self%standard_boundaries_t%init(patch)

  call self%initial_conditions_copy(self%patch%mesh, self%patch%idx)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%mesh(2)%lower_boundary) call self%lower_y(mem)
  if (self%mesh(2)%upper_boundary) call self%upper_y(mem)
  if (self%mesh(3)%lower_boundary) call self%lower_z(mem)
  if (self%mesh(3)%upper_boundary) call self%upper_z(mem)

END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)
  !if (self%mesh(3)%lower_boundary) call self%ddt_boundary_bot(mem)
END SUBROUTINE derivative

SUBROUTINE lower_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d ,idx%d , select=bits%yl)
  call self%symmetric(e ,idx%e , select=bits%yl)
  call self%symmetric(px,idx%px, select=bits%yl)
  call self%symmetric(py,idx%py, select=bits%yl)
  call self%symmetric(pz,idx%pz, select=bits%yl)
  call self%symmetric(bx,idx%bx, select=bits%yl)
  call self%symmetric(by,idx%by, select=bits%yl)
  call self%symmetric(bz,idx%bz, select=bits%yl)
END SUBROUTINE lower_y

SUBROUTINE upper_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d ,idx%d , select=bits%yu)
  call self%symmetric(e ,idx%e , select=bits%yu)
  call self%symmetric(px,idx%px, select=bits%yu)
  call self%symmetric(py,idx%py, select=bits%yu)
  call self%symmetric(pz,idx%pz, select=bits%yu)
  call self%symmetric(bx,idx%bx, select=bits%yu)
  call self%symmetric(by,idx%by, select=bits%yu)
  call self%symmetric(bz,idx%bz, select=bits%yu)
END SUBROUTINE upper_y

SUBROUTINE lower_z(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem!, ic_copy
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  integer                        :: lb, lo
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  
  if(self%patch%id==11) print *,'pre-bn', minval(d)
  call self%symmetric(d,idx%d, select=bits%zl)
  if(self%patch%id==11) print *,'post-bn', minval(d)
  call self%symmetric(e,idx%e, select=bits%zl)

  call self%symmetric(px,idx%px, select=bits%zl)
  call self%smooth_reduce(py,idx%py, select=bits%zl)
  call self%smooth_reduce(pz,idx%pz, select=bits%zl)

  call self%static_values(bx,idx%bx, select=bits%zl)
  call self%static_values(by,idx%by, select=bits%zl)
  call self%static_values(bz,idx%bz, select=bits%zl)

END SUBROUTINE lower_z

SUBROUTINE upper_z(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  call self%symmetric(d,idx%d, select=bits%zu)
  call self%symmetric(e,idx%e, select=bits%zu)

  call self%antisymmetric(px,idx%px, select=bits%zu)
  call self%symmetric(py,idx%py, select=bits%zu)
  call self%antisymmetric(pz,idx%pz, select=bits%zu)

  call self%extrapolate(bx,idx%bx, select=bits%zu)
  call self%extrapolate(by,idx%by, select=bits%zu)
  call self%extrapolate(bz,idx%bz, select=bits%zu)

  call self%add_flow(py, idx%py, select=bits%zu)

END SUBROUTINE upper_z

SUBROUTINE initial_conditions_copy(self, m, idx)
  class(boundaries_t):: self
  type(index_t)      :: idx
  class(mesh_t), pointer, dimension(:):: m
  ! ----------------------------------------------------------------------------
  integer:: ix, iy, iz, n
  real :: entropy, pg, pb, bb, bx_sum, d, bx, by, bz, bxp, byp, bzp
  real(8) :: xs, ys, zs, xp, yp, zp
  logical, save                  :: first_time_ic=.true.
  real                           :: gamma=1.6666666667
  character(len=20)               :: solver
  real                           :: d0, beta, alpha, idir
  namelist /initial_params/ d0, beta, alpha, idir
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !-----------------------------------------------------------------------------

  allocate(self%ic_copy(m(1)%gn,m(2)%gn,m(3)%gn,8))

  solver='stagger2e_bifrost'

  d0=1.0; beta=0.02; alpha=2.0; idir=1

!  pg = 0.01
!  bb = sqrt(2*pg/beta/4/math%pi) !mu_0 possibly
  bb = 3.!1./sqrt(4./math%pi) !mu_0 possibly
  d = 20.

  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz)
    zs = zp - 0.5*m(3)%d
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix)
        xs = xp - 0.5*m(1)%d
!        bx_sum = 0.
!        do n=1,25
!          bx_sum=bx_sum+(((1-(yp**2+zp**2)/(d**2*n**2))*((8*zp)/(d**2*n**2)- &
!                        (4*zp*(1-(yp**2+zp**2)/(d**2*n**2)))/(d**2*n**2)))/  &
!                        (n**2*((4*zp**2)/(d**2*n**2)+(1-(yp**2+zp**2)/       &
!                        (d**2*n**2))**2)**2))+(2*zp)/(d**2*n**4*((4*zp**2)/  &
!                        (d**2*n**2)+(1-(yp**2+zp**2)/(d**2*n**2))**2))
!        end do 
        bx = 0.0 !alpha*bb*yp*(2/d)**2*bx_sum
        by = 4*alpha*bb*ys*zp/(ys**2+zp**2)**2
        bz = -bb*(1+2*alpha*(yp**2-zs**2)/(yp**2+zs**2)**2)
        bxp = 0.0 !alpha*bb*yp*(2/d)**2*bx_sum
        byp = 4*alpha*bb*yp*zp/(yp**2+zp**2)**2
        bzp = -bb*(1+2*alpha*(yp**2-zp**2)/(yp**2+zp**2)**2)
        self%ic_copy(ix,iy,iz,idx%d)=d0
        self%ic_copy(ix,iy,iz,idx%px)=0.0
        self%ic_copy(ix,iy,iz,idx%py)=0.0
        self%ic_copy(ix,iy,iz,idx%pz)=0.0
        self%ic_copy(ix,iy,iz,idx%bx)=bx !0.0!alpha*bb*yp*(2/d)**2*bx_sum 
        self%ic_copy(ix,iy,iz,idx%by)=by !4*alpha*bb*ys*zp/(ys**2+zp**2)**2
        self%ic_copy(ix,iy,iz,idx%bz)=bz !-bb*(1+2*alpha*(yp**2-zs**2)/(yp**2+zs**2)**2)
        !pb = 0.5*(bx**2+by**2+bz**2)
        !pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = beta*bb
        if (gamma==1.0) then
          self%ic_copy(ix,iy,iz,idx%s) = 1.0
        else if (solver(1:14) == 'zeus_mhd_patch' .or. &
                 solver(1:10) == 'stagger2e_') then 
          self%ic_copy(ix,iy,iz,idx%e) = pg / (gamma - 1.0)
        else if (solver(1:9) == 'stagger2_') then
          self%ic_copy(ix,iy,iz,idx%s) = &
              (log(pg)-log(d0)*gamma)/(gamma-1.0)*d0
        else if (solver(1:6) == 'ramses')  then
          self%ic_copy(ix,iy,iz,idx%e) = pg / (gamma - 1.0)
        end if
      end do
    end do
  end do
  end associate

END SUBROUTINE initial_conditions_copy

SUBROUTINE static_values(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  integer :: lb, lo, iv
  integer, optional :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  if (self%is_set(bits%zl)) then
    lb = m(3)%lb
    lo = m(3)%lo
    mem(:,:,lb:lo) = self%ic_copy(:,:,lb:lo, iv)
  end if
END SUBROUTINE static_values

SUBROUTINE add_flow(self, mem, iv, select)
  class(boundaries_t)     :: self
  class(mesh_t), pointer  :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer                 :: ub, uo, iv
  real                    :: t0, v0, flow, t
  integer, optional       :: select
  !-----------------------------------------------------------------------------
  t0 = 100.0
  v0 = 0.45
  m => self%mesh
  t = timer%time

  if (self%is_set(bits%zu)) then
    uo = m(3)%uo
    ub = m(3)%ub
    if (t < t0) then
      flow = v0/2.*(1.-cos(math%pi*t/t0)) !-v0/2.*(1.-cos(math%pi*(t-(t-1.))/t0))
      mem(:,:,uo:ub) = flow 
    else
      mem(:,:,uo:ub) =  v0
    end if
  end if

END SUBROUTINE add_flow

SUBROUTINE smooth_reduce(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  integer :: lb, lo, iv, z
  integer, optional :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  if (self%is_set(bits%zl)) then
    lb = m(3)%lb
    lo = m(3)%lo
    do z=lo,lb,-1
      mem(:,:,z) = mem(:,:,z+1)*0.7
    end do
  end if
END SUBROUTINE smooth_reduce

END MODULE boundaries_mod
