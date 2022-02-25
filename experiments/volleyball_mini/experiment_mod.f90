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
  USE timer_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(sphere_point_t), pointer:: point
    real(8):: r_min, r_max
  contains
    procedure:: init
    procedure:: initial_conditions
    procedure:: bndries
    procedure:: get_position
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
  real:: amplitude=1e-1
CONTAINS

!===============================================================================
!> Initialize parameters
!> Generate 3-D values
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  integer:: iostat, i
  namelist /experiment_params/ verbose, detailed_timer, initial_bc, generic_IC, perturb
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init', itimer=itimer)
  call io%lock%set ('experiment_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    call model_1d%init ()
    dbglog = merge(stderr, stdlog, verbose > 1)
  end if
  call io%lock%unset ('experiment_t%init')
  self%mhd=.false.
  !-----------------------------------------------------------------------------
  ! Rescale the two parameters, position and size, which uniquely determine the
  ! mesh and related attributes, and then initialize
  !-----------------------------------------------------------------------------
  !self%position = self%point%p *cgs%r_sun/scaling%l
  !self%llc_nat  = self%llc_nat *cgs%r_sun/scaling%l
  !self%size     = self%size    *cgs%r_sun/scaling%l
  !self%box      = self%box     *cgs%r_sun/scaling%l
  allocate(self%point)
  self%point%p  = self%position
  call self%point%to_spherical ()
  ! redo unit vectors
  !call self%point%unit_vectors ()
  self.point.e(:,1) = [1.0,0.0,0.0]
  self.point.e(:,2) = [0.0,1.0,0.0]
  self.point.e(:,3) = [0.0,0.0,1.0]
  !allocate(self%erot(3,3))
  !self%erot = self%point%e
  !
  call self%solver_t%init ()
  !-----------------------------------------------------------------------------
  ! Make darn sure that all position information agrees
  !-----------------------------------------------------------------------------
  if (any(self%mesh%p          /= self%position .or. &
          self%mesh%centre_nat /= self%position)) then
    call io%warning ('volleyball_t%inconsitent positions')
    self%mesh%p           = self%position
    self%mesh%centre_nat  = self%position
  end if
  !self%periodic = .false.
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

  ! ----------------------------------------------------------------------------
  ! Set up initial conditions
  ! ----------------------------------------------------------------------------
  call self%initial_conditions()
  if (perturb) call self%perturbations()
  ! ----------------------------------------------------------------------------
  ! Prepare all extra physics
  ! ----------------------------------------------------------------------------
  call self%prepare_cooling()
  call self%prepare_gravity()


  call trace_end (itimer)
END SUBROUTINE init

SUBROUTINE perturbations(self)
  class(experiment_t):: self
  ! ----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:), pointer:: d,e,ee
  real, parameter:: pi = 4.*atan(1.)
  integer:: i, j, k, k1, k2, n(3)
  real:: a, kx, ky, d1, f
  real(8), pointer:: x(:), y(:), z(:)
  !-----------------------------------------------------------------------------
  m => self%mesh
  d => self%mem(:,:,:,self%idx%d,1,1)
  e => self%mem(:,:,:,self%idx%e,1,1)
  n=shape(d)
  allocate (x (n(1)),y(n(2)),z(n(3)))
  allocate (ee(n(1),n(2),n(3)))
  ee=e/d
  x = m(1)%r+m(1)%p
  y = m(2)%r+m(2)%p
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
          kx = 2.*pi/4.*k1*x(i) - 11.3*k1**0.5
          ky = 2.*pi/4.*k2*y(j) - 50.1*k2**0.2
          a = amplitude/(1+(k1-1)+3*(k2-1))
          f = f + a*sin(kx)*sin(ky)*exp(-(z(k)-cgs%r_sun/scaling%l)**2/0.01)
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


! ------------------------------------------------------------------------------
! Set up the initial HSE from JCD Model-S
! ------------------------------------------------------------------------------
SUBROUTINE initial_conditions(self)
  class(experiment_t):: self
  ! ----------------------------------------------------------------------------
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
        r = sqrt(pos(3)**2)
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

  ! ----------------------------------------------------------------------------
  ! If the top boundary is at the Solar surface, we need to use a different BC
  ! ----------------------------------------------------------------------------
  if (r_max>=cgs%r_sun/scaling%l) then
    self%boundaries%surface=.true.
  end if
  if (generic_IC) then
    call self%from_generic()
    call self%deallocate_generic()
  end if

END SUBROUTINE initial_conditions

SUBROUTINE prepare_cooling(self)
  class(experiment_t):: self
  ! ----------------------------------------------------------------------------
  real:: e, d
  real(8):: pos(3)
  ! ----------------------------------------------------------------------------
  if (.not. self%newton%on) return

  ! ----------------------------------------------------------------------------
  ! Set the cooling to be at the surface of the Sun. We need to orientate
  ! the position to be in the patch reference frame. 
  ! To do that, we use matmul(pos,e)
  ! ----------------------------------------------------------------------------
  pos = [0d0,0d0,cgs%r_sun/scaling%l]
  pos = matmul(pos,self%point%e)
  self%newton%position= pos(3)
  ! ----------------------------------------------------------------------------
  ! Interpolate the specific energy at the surface from JCD model
  ! ----------------------------------------------------------------------------
  d = interpolate%one_d (pos(3), model_1d%r, model_1d%d)
  e = interpolate%one_d (pos(3), model_1d%r, model_1d%e)
  self%newton%ee_top = e/d

END SUBROUTINE prepare_cooling

SUBROUTINE prepare_gravity(self)
  class(experiment_t):: self
  ! ----------------------------------------------------------------------------
  class(mesh_t), pointer:: m(:)
  real(8):: pos(3), rp, rhat(3), g, xp, yp, zp, vect(3)
  integer:: n(3), ix, iy, iz
  real(8), pointer:: e(:,:)
  ! ----------------------------------------------------------------------------
  if (.not. self%gravity%on) return

  n = self%gn
  allocate (self%gravity%vector_field(n(1),n(2),n(3),3))
  self%gravity%vector_field = 0.0
  m => self%mesh
  e => self%point%e
  do iz=m(3)%lb, m(3)%ub
    zp = m(3)%r(iz)
    do iy=m(2)%lb, m(2)%ub
      yp = m(2)%r(iy)
      do ix=m(1)%lb, m(1)%ub
        xp = m(1)%r(ix)
        pos = self%loc2glob([xp,yp,zp])
        rp=sqrt(pos(3)**2)          ! needs to be all components in full VB
        rp=max(rp,self%gravity%minrad)
        rhat = pos/rp
        ! ----------------------------------------------------------------------
        ! interpolate table values and store them
        ! ----------------------------------------------------------------------
        g = interpolate%one_d (rp, model_1d%r, model_1d%g)
        vect=matmul(-g*rhat,self%point%e)
        self%gravity%vector_field(ix,iy,iz,3)=vect(3)  ! needs to be all components in full VB
      end do
    end do
  end do

END SUBROUTINE prepare_gravity

SUBROUTINE bndries(self)
  class(experiment_t):: self
  ! ----------------------------------------------------------------------------
  class(mesh_t), pointer:: m1, m2, m3
  real(var_k), dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz
  real(var_k), dimension(:,:,:,:), pointer:: mem
  integer:: i, i1, i2, i3, i0, iostat, id
  real(8):: x1, x2, x3, pos(3), r, d1, e1, max_r
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  if (self%idx%bx>0) then
    bx => self%mem(:,:,:,self%idx%bx,self%it,1)
    by => self%mem(:,:,:,self%idx%by,self%it,1)
    bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  end if

  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  !if (self%id==1) print *, 'minmax d,e before bnd', minval(d), maxval(d)
  !if (self%id==1) print *, 'minmax d,e before bnd',minval(e), maxval(e)
  !
  if (initial_bc) then
    if (self%boundaries%is_set (bits%zu)) then    
      do i3=m3%uo,m3%ub
        x3 = m3%r(i3)
        do i2=m2%lb,m2%ub
          x2 = m2%r(i2)
          do i1=m1%lb,m1%ub
            x1 = m1%r(i1)
            pos=self%loc2glob([x1,x2,x3])
            r = sqrt(pos(3)**2)
            d1 = interpolate%one_d (r, model_1d%r, model_1d%d)
            e1 = interpolate%one_d (r, model_1d%r, model_1d%e)
            d(i1,i2,i3) = d1
            e(i1,i2,i3) = e1
            px(i1,i2,i3) = 0.0
            py(i1,i2,i3) = 0.0
            pz(i1,i2,i3) = 0.0
            if (self%idx%bx>0) then
              bx(i1,i2,i3) =0.0
              by(i1,i2,i3) =0.0
              bz(i1,i2,i3) =0.0
            end if
            if (associated(self%filled)) self%filled(i1,i2,i3) = self%level
          end do
        end do
      end do  
    end if
    if (self%boundaries%is_set (bits%zl)) then
      do i3=m3%lb,m3%lo
        x3 = m3%r(i3)
        do i2=m2%lb,m2%ub
          x2 = m2%r(i2)
          do i1=m1%lb,m1%ub
            x1 = m1%r(i1)
            pos=self%loc2glob([x1,x2,x3])
            r = sqrt(pos(3)**2)
            d1 = interpolate%one_d (r, model_1d%r, model_1d%d)
            e1 = interpolate%one_d (r, model_1d%r, model_1d%e)
            d(i1,i2,i3) = d1
            e(i1,i2,i3) = e1
            px(i1,i2,i3) = 0.0
            py(i1,i2,i3) = 0.0
            pz(i1,i2,i3) = 0.0
            if (self%idx%bx>0) then
              bx(i1,i2,i3) =0.0
              by(i1,i2,i3) =0.0
              bz(i1,i2,i3) =0.0
            end if
            if (associated(self%filled)) self%filled(i1,i2,i3) = self%level
          end do
        end do
      end do
    end if
  else
    mem => self%mem(:,:,:,:,self%it,1)
    call self%boundaries%primary(mem)
  end if
END SUBROUTINE bndries

FUNCTION loc2glob(self,pos) RESULT (out)
  class(experiment_t):: self
  real(8):: pos(3), out(3)
  ! ----------------------------------------------------------------------------
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
!> Return the position and velocity of the source in the target coordinates at
!> the target time, and the relative orientation of the source at that time.
!===============================================================================
SUBROUTINE get_position (self, source, p, v, e)
  class(experiment_t):: self
  class(task_t):: source
  real(8):: time, p(3), v(3), e(3,3)
  optional:: p, v, e
  real(8):: crn(8,3)
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('experiment_t%get_position', itimer=itimer)
  select type (source)
  class is (experiment_t)    
    !---------------------------------------------------------------------------
    ! Relative position and velocity in target coordinates
    !---------------------------------------------------------------------------
    if (present(p)) then
      p = source%position-self%position
      p = matmul (p , self%point%e)
    end if
    !---------------------------------------------------------------------------
    ! Return the relative rotation matrix, at the given target time
    !---------------------------------------------------------------------------
    if (present(e)) then
      e=matmul(transpose(source%point%e),self%point%e)
    end if
  end select
  if (timer%detailed) &
    call trace%end (itimer)
END SUBROUTINE get_position

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

END MODULE experiment_mod
