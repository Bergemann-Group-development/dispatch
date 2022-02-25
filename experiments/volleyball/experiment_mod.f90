!===============================================================================
!> Minimal experiment initialization -- all other calls fall through to solver
!>
!> TODO:
!> o) Coriolis forces
!> o) Force of gravity
!> o) Newton cooling
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE omp_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE link_mod
  USE task_mod
  USE solver_mod
  USE model_1d_mod
  USE scaling_mod
  USE units_mod
  USE region_mod
  USE sphere_point_mod
  USE interpolate_mod
  USE lagrange_mod
  USE timer_mod
  USE eos_mod
  USE os_mod
  USE bits_mod
  USE random_mod
  USE newton_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(sphere_point_t), pointer:: point
    real(8):: r_min, r_max
    integer:: layer
  contains
    procedure:: init
    procedure:: init_level
    procedure:: update
    procedure:: initial_conditions
    procedure:: bndries
    procedure:: loc2glob
    procedure:: info_e
    procedure:: prepare_cooling
    procedure:: prepare_gravity
    procedure:: perturbations
  end type
  integer:: verbose=0, dbglog=0
  logical:: detailed_timer=.false.
  logical:: initial_bc=.false.
  logical:: generic_IC=.false.
  logical:: perturb=.false.
  real:: a_perturb=1e-1, r_perturb=660., w_perturb=1., f_newton=0.9
CONTAINS

!===============================================================================
!> Initialize parameters
!> Generate 3-D values
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  integer:: iostat, i
  namelist /experiment_params/ verbose, detailed_timer, initial_bc, &
    generic_IC, perturb, r_perturb, w_perturb, a_perturb, f_newton
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init', itimer=itimer)
  call io%lock%set ('experiment_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, experiment_params)
    end if
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    call model_1d%init ()
    dbglog = merge(stderr, stdlog, verbose > 1)
  end if
  call io%lock%unset ('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! Rescale the two parameters, position and size, which uniquely determine the
  ! mesh and related attributes, and then initialize
  !-----------------------------------------------------------------------------
  self%position = self%point%p *cgs%r_sun/scaling%l
  self%llc_nat  = self%llc_nat *cgs%r_sun/scaling%l
  self%size     = self%size    *cgs%r_sun/scaling%l
  self%box      = self%box     *cgs%r_sun/scaling%l
  self%point%p  = self%position
  self%r_max = self%r_max*cgs%r_sun/scaling%l
  self%r_min = self%r_min*cgs%r_sun/scaling%l
  call self%point%to_spherical ()
  ! redo unit vectors
  call self%point%unit_vectors ()
  self%erot = self%point%e
  !-----------------------------------------------------------------------------
  ! This generates correctly scaled mesh%r, but mesh%p is set = self%position
  !-----------------------------------------------------------------------------
  call self%solver_t%init ()
  !-----------------------------------------------------------------------------
  ! Make sure that centre_nat agrees, and set mesh%p to the layer radius
  !-----------------------------------------------------------------------------
  if (any(self%mesh%centre_nat /= self%position)) then
    call io%warning ('volleyball_t%inconsitent positions')
    self%mesh%centre_nat  = self%position
  end if
  self%mesh%p = sqrt(sum(self%position**2))
  self%periodic = .false.
  !-----------------------------------------------------------------------------
  ! Assign regions of authority and interest
  !-----------------------------------------------------------------------------
  self%authority%l = -self%size/2d0
  self%authority%u = +self%size/2d0
  self%interest%l = -self%size/2d0 - self%mesh%ng*self%mesh%d
  self%interest%u = +self%size/2d0 + self%mesh%ng*self%mesh%d
  if (verbose > 2) then
    do i=1,3
      call self%mesh(i)%print
    end do
    call self%point%info ('experiment_t%init')
  end if
  !-----------------------------------------------------------------------------
  ! Set up initial conditions
  !-----------------------------------------------------------------------------
  call self%initial_conditions()
  if (perturb) call self%perturbations()
  !-----------------------------------------------------------------------------
  ! Prepare all extra physics
  !-----------------------------------------------------------------------------
  call self%prepare_cooling()
  call self%prepare_gravity()
  !
  call trace_end (itimer)
END SUBROUTINE init

!===============================================================================
!> Set a fixed initial level for this experiment
!===============================================================================
SUBROUTINE init_level (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  self%level = 10
END SUBROUTINE init_level

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call self%solver_t%update ()
END SUBROUTINE update

!===============================================================================
!> Perturbations intended to start convection near the surface
!===============================================================================
SUBROUTINE perturbations(self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:), pointer:: d,e,ee
  real, parameter:: pi = 4.*atan(1.)
  integer:: i, j, k, k1, k2, n(3)
  real:: a, kx, ky, d1, f, a1, a2
  real(8), pointer:: x(:), y(:), z(:)
  !-----------------------------------------------------------------------------
  m => self%mesh
  d => self%mem(:,:,:,self%idx%d,1,1)
  e => self%mem(:,:,:,self%idx%e,1,1)
  n=shape(d)
  allocate (x (n(1)),y(n(2)),z(n(3)))
  allocate (ee(n(1),n(2),n(3)))
  ee=e/d
  x = m(1)%r/m(1)%s
  y = m(2)%r/m(2)%s
  z = m(3)%r+m(3)%p
  n = shape(d)
  do k=1,n(3)
    d1 = d(1,1,k)
    do j=1,n(2)
    do i=1,n(1)
      f = 1.0
      do k1=1,4
      do k2=1,4
        if (sqrt(real(k1)**2+real(k2)**2) <= 4.) then
          a1 = self%random%ran1()
          a2 = self%random%ran1()
          kx = 2.*pi*k1*x(i)
          ky = 2.*pi*k2*y(j)
          a = a_perturb/(1+a1*(k1-1)+a2*(k2-1))
          f = f + a*cos(kx)*cos(ky)*exp(-((z(k)-r_perturb)/w_perturb)**2)
        end if
      end do
      end do
      d(i,j,k) = d1*f
    end do
    end do
  end do
  e=ee*d
  deallocate(x, y, z,ee)
END SUBROUTINE perturbations

!===============================================================================
! Set up the initial HSE from JCD Model-S
!===============================================================================
SUBROUTINE initial_conditions(self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  real(var_k), dimension(:,:,:), pointer:: d, tg, px, py, pz, bx, by, bz, e
  real(var_k), pointer:: ff(:,:,:,:)
  integer :: ix, iy, iz, itg
  real(8) :: r_max, xp, yp, zp, pos(3), r, d1, tg1, e1
  ! ............................................................................ 
  m => self%mesh
  if (generic_IC) then
    call self%allocate_generic()
    ff => self%generic
    itg=max(self%idx%s, self%idx%e)
    tg => ff(:,:,:,itg )
  else
    ff => self%mem(:,:,:,:,self%it,1)
    if (trim(self%etype)=='thermal') e  => ff(:,:,:,self%idx%e)
  end if
  d  => ff(:,:,:,self%idx%d )
  px => ff(:,:,:,self%idx%px)
  py => ff(:,:,:,self%idx%py)
  pz => ff(:,:,:,self%idx%pz)
  if (self%idx%bx>0 ) then
    bx => ff(:,:,:,self%idx%bx)
    by => ff(:,:,:,self%idx%by)
    bz => ff(:,:,:,self%idx%bz)
  end if
  r_max=-1d2
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%r(iz)
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%r(iy)
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%r(ix)
        pos=self%loc2glob([xp,yp,zp])
        r = sqrt(pos(3)**2+pos(2)**2+pos(1)**2)
        r_max=max(r_max,r)
        !-----------------------------------------------------------------------
        ! Setting index l=0 triggers computing and reusing l and p
        !-----------------------------------------------------------------------
        !call model_1d%interpolate (r, d1, e1)
        d1 = interpolate%one_d (r, model_1d%r, model_1d%d )
        if (generic_IC) then
          tg1= interpolate%one_d (r, model_1d%r, model_1d%tg)
          tg(ix,iy,iz) = tg1
        else
          e1= interpolate%one_d (r, model_1d%r, model_1d%e)
          e(ix,iy,iz) = e1
        end if
        d (ix,iy,iz) = d1
        px(ix,iy,iz) = 0.0
        py(ix,iy,iz) = 0.0
        pz(ix,iy,iz) = 0.0
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! If the top boundary is at the Solar surface, we need to use a different BC
  !-----------------------------------------------------------------------------
  if (r_max>=cgs%r_sun/scaling%l) then
    self%boundaries%surface=.true.
  end if
  if (generic_IC) then
    call self%from_generic()
    call self%deallocate_generic()
  end if
END SUBROUTINE initial_conditions

!===============================================================================
!> Place the cooling point at the center of the first patch, which by
!> constuction is in the outermost layer
!===============================================================================
SUBROUTINE prepare_cooling(self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  real:: e, d
  real(8):: pos, r_pos, ee_top, ee_slope
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  if (.not. self%newton%on) return
  call trace%begin ('experiment_t%prepare_cooling')
  !-----------------------------------------------------------------------------
  ! Get the chosen Newton cooling position, and interpolate to find the
  ! energy per unit mass there and the slope
  !-----------------------------------------------------------------------------
  call self%newton%get (position=pos)
  d = interpolate%one_d (pos, model_1d%r, model_1d%d)
  e = interpolate%one_d (pos, model_1d%r, model_1d%e)
  ee_top = e/d
  d = interpolate%one_d (pos-1., model_1d%r, model_1d%d)
  e = interpolate%one_d (pos-1., model_1d%r, model_1d%e)
  ee_slope = (ee_top-e/d)
  !-----------------------------------------------------------------------------
  ! Reduce the energy to cool to with f_newton
  !-----------------------------------------------------------------------------
  ee_top = ee_top*f_newton
  !-----------------------------------------------------------------------------
  ! Set the property via a %set() function, not directly (transparency)
  !-----------------------------------------------------------------------------
  call self%newton%set (ee0=ee_top, ee1=ee_slope)
  if (self%id == 1) &
    write (stdout,*) &
    'experiment_t%prepare_cooling: pos,ee_top,ee_slope=', pos, ee_top, ee_slope
  call trace%end ()
END SUBROUTINE prepare_cooling

SUBROUTINE prepare_gravity(self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  real(8):: pos(3), rp, rhat(3), g, xp, yp, zp, vect(3)
  integer:: n(3), ix, iy, iz
  real(8), pointer:: e(:,:)
  real(8):: r
  !-----------------------------------------------------------------------------
  if (.not. self%gravity%on) return

  n = self%gn
  allocate (self%gravity%vector_field(n(1),n(2),n(3),3))
  self%gravity%vector_field = 0.0
  m => self%mesh
  e => self%point%e
  r = sqrt(sum(self%position**2))
  do iz=m(3)%lb, m(3)%ub
    zp = m(3)%r(iz) + r
    do iy=m(2)%lb, m(2)%ub
      yp = m(2)%r(iy)
      do ix=m(1)%lb, m(1)%ub
        xp = m(1)%r(ix)
            rp=sqrt(xp**2+yp**2+zp**2)
        !-----------------------------------------------------------------------
        ! interpolate table values and store them
        !-----------------------------------------------------------------------
        g = lagrange%interpolate1d (rp, model_1d%r, model_1d%g)
        vect=[xp,yp,zp]
        self%gravity%vector_field(ix,iy,iz,:) = -g*vect/sqrt(sum(vect**2))
      end do
    end do
  end do
END SUBROUTINE prepare_gravity

FUNCTION loc2glob(self,pos) RESULT (out)
  class(experiment_t):: self
  real(8):: pos(3), out(3)
  !-----------------------------------------------------------------------------
  out = matmul(self%point%e,pos)+self%position
END FUNCTION loc2glob

!===============================================================================
!> Return the previous index corresponding to a point in my coordinates, and
!> the fractional remaining float index
!===============================================================================
SUBROUTINE my_index (self, p_in, i4, ii, p, valid)
  class(experiment_t)    :: self
  real(8), dimension(3)  :: p_in, p
  integer                :: ii(3), i4
  class(mesh_t), pointer :: m1, m2, m3
  logical                :: valid
  real, parameter        :: margin=0.5d0
  !.............................................................................
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  p  = p_in/self%mesh%d + self%mesh%o - [m1%h(i4),m2%h(i4),m3%h(i4)]
  ii = floor(p)
  ii = min (max (ii, self%mesh%li), self%mesh%ui-1)
  p  = p-ii
  valid = all(p > -margin .and. p < 1d0+margin)
END SUBROUTINE my_index

!===============================================================================
!> Return the local position, given indices and component staggering
!===============================================================================
FUNCTION my_position (self, i1, i2, i3, i4) RESULT (p)
  class(experiment_t)    :: self
  integer                :: i1, i2, i3, i4
  real(8)                :: p(3)
  !-----------------------------------------------------------------------------
  p(1) = self%mesh(1)%r(i1) + self%mesh(1)%h(i4)*self%mesh(1)%d
  p(2) = self%mesh(2)%r(i2) + self%mesh(2)%h(i4)*self%mesh(2)%d
  p(3) = self%mesh(3)%r(i3) + self%mesh(3)%h(i4)*self%mesh(3)%d
END FUNCTION my_position

!===============================================================================
!> Write out info on the position, size, and orientation of self
!===============================================================================
SUBROUTINE info_e (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  if (verbose > 2) then
    write (dbglog,1) 'experiment_t%info: position=', self%position
    write (dbglog,1) 'experiment_t%info:     size=', self%size
    write (dbglog,1) 'experiment_t%info:       e1=', self%point%e(:,1)
    write (dbglog,1) 'experiment_t%info:       e2=', self%point%e(:,2)
    write (dbglog,1) 'experiment_t%info:       e3=', self%point%e(:,3)
    1 format(a,3f12.6)
  end if
END SUBROUTINE info_e

!===============================================================================
!> Utility functions for single 3-D vector operations
!===============================================================================
FUNCTION dot (a, b) RESULT (c)
  real(8):: a(3), b(3), c
  c = sum(a*b)
END FUNCTION dot
!===============================================================================
FUNCTION normalize (a) RESULT (b)
  real(8):: a(3), b(3)
  b = a/sqrt(sum(a**2))
END FUNCTION normalize
!===============================================================================
FUNCTION cross (a, b) RESULT (c)
  real(8):: a(3), b(3), c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION cross

!===============================================================================
!===============================================================================
SUBROUTINE bndries(self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zl)) then
    call lower_bc (self)
  end if
  if (self%boundaries%is_set (bits%zu)) then
    call upper_bc (self)
  end if
END SUBROUTINE bndries

!===============================================================================
!===============================================================================
SUBROUTINE upper_bc (self)
  class(experiment_t):: self
  !.............................................................................
  real(var_k),                   contiguous, pointer :: uu(:,:,:,:)
  real,                          contiguous, pointer :: g(:,:,:)
  real(var_k), dimension(:,:,:), contiguous, pointer :: d, e, s, p, ux, uy, uz
  real(var_k), dimension(:,:,:), contiguous, pointer :: px, py, pz, eth
  class(mesh_t), pointer :: m
  real(8)                :: dz
  integer                :: iz, jz, it, ix, iy, l(3), u(3), n(3)
  integer                :: zi
  integer, save          :: verbose_l=1
  !-----------------------------------------------------------------------------
  if (verbose_l > 1) then
    PRINT *, self%id, self%it, self%idx%d, self%idx%e, self%idx%s
  end if
  it = self%it
  m => self%mesh(3)
  l=self%mesh%lb
  u=self%mesh%ub
  l(3)=self%mesh(3)%ui
  n = self%gn
  n(3) = self%ng(3)+1
  d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,it,1)
  g => self%gravity%vector_field(l(1):u(1),l(2):u(2),l(3):u(3),3)
  call os%allocate (p, n)
  if (trim(self%etype)=='thermal') then
    eth => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e,it,1)
    call eos%lookup (d=d, e=eth, pg=p)
  else if (trim(self%etype)=='total') then
    call os%allocate (eth, n)
    e  => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e ,it,1)
    px => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px,it,1)
    py => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%py,it,1)
    pz => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%pz,it,1)
    do iz=1,n(3)
      do iy=1,n(2)
        do ix=1,n(1)
          eth(ix,iy,iz) =  e (ix,iy,iz) - &
                     0.5* (px(ix,iy,iz)**2 + &
                           py(ix,iy,iz)**2 + &
                           pz(ix,iy,iz)**2)/ &
                            d(ix,iy,iz)
        end do
      end do
    end do
    call eos%lookup (d=d, e=eth, pg=p)
  else if (trim(self%etype)=='entropy') then
    s => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%s,it,1)
    call eos%lookup (d=d, s=s, pg=p)
  end if
  !-----------------------------------------------------------------------------
  ! Density BC
  !-----------------------------------------------------------------------------
  jz=m%ui
  do iz=2,n(3)
    dz=m%r(jz+iz-1)-m%r(jz)
    d(:,:,iz) = d(:,:,1)*exp(dz*g(:,:,1)*d(:,:,1)/p(:,:,1))
    !---------------------------------------------------------------------------
    ! Energy or entropy BC: Constant energy or entropy per unit mass
    !---------------------------------------------------------------------------
    if (self%idx%e > 0) then
      eth(:,:,iz) = eth(:,:,1)/d(:,:,1)*d(:,:,iz)
    else if (self%idx%s > 0) then
      s(:,:,iz) = s(:,:,1)/d(:,:,1)*d(:,:,iz)
    end if
  end do
  call os%deallocate (p)
  !-----------------------------------------------------------------------------
  ! Velocity boundary conditions: Symmtric across boundary
  !-----------------------------------------------------------------------------
  call os%allocate (uu, [self%gn,3])
  call self%p2U (uu, self%it)
  ux => uu(:,:,:,1); uy => uu(:,:,:,2); uz => uu(:,:,:,3)
  call self%boundaries%symmetric (ux, self%idx%px)
  call self%boundaries%symmetric (uy, self%idx%py)
  call self%boundaries%symmetric (uz, self%idx%pz)
  call self%U2p (uu, self%it)
  call os%deallocate (uu)
  ! ----------------------------------------------------------------------------
  ! If total energy is used, convert back to it.
  ! ----------------------------------------------------------------------------
  if (trim(self%etype)=='total') then
    do iz=2,n(3)
      do iy=1,n(2)
        do ix=1,n(1)
          e(ix,iy,iz) = eth(ix,iy,iz) + &
                   0.5* (px(ix,iy,iz)**2 + &
                         py(ix,iy,iz)**2 + &
                         pz(ix,iy,iz)**2)/ &
                          d(ix,iy,iz)
        end do
      end do
    end do
    call os%deallocate(eth)
  end if
  !-----------------------------------------------------------------------------
  if (associated(self%filled)) then
    self%filled(:,:,m%uo:m%ub) = self%level
  end if
END SUBROUTINE upper_bc

!===============================================================================
!===============================================================================
SUBROUTINE lower_bc (self)
  class(experiment_t):: self
  !.............................................................................
  real(var_k),                   contiguous, pointer :: uu(:,:,:,:)
  real(var_k), dimension(:,:,:), contiguous, pointer :: d, e, eth, s, ux, uy, uz
  real(var_k), dimension(:,:,:), contiguous, pointer :: px,py, pz
  real(var_k), dimension(:,:,:), contiguous, pointer :: lnf
  integer :: it, ix, iy, iz, l(3), u(3)
  !-----------------------------------------------------------------------------
  it = self%it
  l=self%mesh%lb
  u=self%mesh%ub
  d => self%mem(:,:,:,self%idx%d,it,1)
  if (trim(self%etype)=='thermal') then
    eth => self%mem(:,:,:,self%idx%e,it,1)
  else if (trim(self%etype)=='total') then
    call os%allocate (eth, self%gn)
    e  => self%mem(:,:,:,self%idx%e ,it,1)
    px => self%mem(:,:,:,self%idx%px,it,1)
    py => self%mem(:,:,:,self%idx%py,it,1)
    pz => self%mem(:,:,:,self%idx%pz,it,1)
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          eth(ix,iy,iz) =  e (ix,iy,iz) - &
                     0.5* (px(ix,iy,iz)**2 + &
                           py(ix,iy,iz)**2 + &
                           pz(ix,iy,iz)**2)/ &
                            d(ix,iy,iz)
        end do
      end do
    end do
  else if (trim(self%etype)=='entropy') then
    s => self%mem(:,:,:,self%idx%s,it,1)
  end if
  !-----------------------------------------------------------------------------
  ! Density and energy BCs
  !-----------------------------------------------------------------------------
  call os%allocate (lnf, self%gn)
  lnf = log(d)
  call self%boundaries%extrapolate (lnf, self%idx%d)
  if (trim(self%etype)=='entropy') then
    lnf = log(s)
    call self%boundaries%extrapolate (lnf, self%idx%s)
  else
    lnf = log(eth)
    call self%boundaries%extrapolate (lnf, self%idx%e)
  end if
  call os%deallocate (lnf)
  !-----------------------------------------------------------------------------
  ! Velocity boundary conditions: No-shear ux & uy, vanishing uz
  !-----------------------------------------------------------------------------
  call os%allocate (uu, [self%gn,3])
  call self%p2U (uu, it)
  ux => uu(:,:,:,1); uy => uu(:,:,:,2); uz => uu(:,:,:,3)
  call self%boundaries%symmetric (ux, self%idx%px)
  call self%boundaries%symmetric (uy, self%idx%py)
  call self%boundaries%antisymmetric (uz, self%idx%pz)
  call self%U2p (uu, it)
  call os%deallocate (uu)
  ! ----------------------------------------------------------------------------
  ! If total energy is used, convert back to it.
  ! ----------------------------------------------------------------------------
  if (trim(self%etype)=='total') then
    call os%allocate (eth, self%gn)
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          e(ix,iy,iz) = eth(ix,iy,iz) + &
                   0.5* (px(ix,iy,iz)**2 + &
                         py(ix,iy,iz)**2 + &
                         pz(ix,iy,iz)**2)/ &
                          d(ix,iy,iz)
        end do
      end do
    end do
    call os%deallocate(eth)
  end if
  !-----------------------------------------------------------------------------
  if (associated(self%filled)) then
    self%filled(:,:,self%mesh(3)%lb:self%mesh(3)%li) = self%level
  end if
END SUBROUTINE lower_bc

END MODULE experiment_mod
