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
  USE standard_boundaries_mod
  USE eos_mod
  USE patch_mod
  USE stagger_mod
  USE vector_ops
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real(kind=KindScalarVar) :: t0, v0, vanishing_factor, beta, alpha, d0, d, bb
    logical :: vecpot
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
    procedure:: lower_x
    procedure:: upper_x
    procedure:: lower_y
    procedure:: upper_y
    procedure:: lower_z
    procedure:: upper_z
    procedure:: static_values
    procedure:: add_flow
    procedure:: vanishing
  end type
  integer, save:: verbose=1
  ! --- add_flow parameters
  real, save:: t0 = 100.
  real, save:: v0 = 0.045
  ! --- vanishing parameter
  ! --- for constant boundary
  real, save:: d0 = 1.0
  real, save:: alpha = 2.0
  real, save:: beta = 0.01
  real, save:: d = 20.0
  real, save:: bb = 3.0
  logical, save:: vecpot = .false.
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t)   :: self
  !.............................................................................
  integer:: iostat
  real, save:: t0, v0
  real, save:: vanishing_factor = 0.7
  logical, save:: first_time=.true.
  namelist /boundary_params/ verbose, t0, v0, vanishing_factor, d0, beta, &
    alpha, bb, d, vecpot
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (io%master) write (*, boundary_params)
    first_time = .false.
  end if
  self%t0 = t0
  self%v0 = v0
  self%vanishing_factor = vanishing_factor
  self%d0 = d0
  self%alpha = alpha
  self%beta = beta
  self%d = d
  self%bb = bb
  self%vecpot = vecpot

  call io%lock%unset('boundaries_t%init')
  call self%standard_boundaries_t%init ()
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%mesh(2)%lower_boundary) call self%lower_y(mem)
  if (self%mesh(2)%upper_boundary) call self%upper_y(mem)
  if (self%mesh(1)%lower_boundary) call self%lower_x(mem)
  if (self%mesh(1)%upper_boundary) call self%upper_x(mem)
  if (self%mesh(3)%lower_boundary) call self%lower_z(mem)
  if (self%mesh(3)%upper_boundary) call self%upper_z(mem)
  call self%set_filled ()
END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
END SUBROUTINE derivative

SUBROUTINE lower_x(self, mem)
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

  call self%antisymmetric(px,idx%px, select=bits%xl)
  call self%antisymmetric(py,idx%py, select=bits%xl)
  call self%antisymmetric(pz,idx%pz, select=bits%xl)
  call self%static_values(d ,idx%d , select=bits%xl)
  call self%static_values(e ,idx%e , select=bits%xl)
  call self%static_values(bx,idx%bx, select=bits%xl)
  call self%static_values(by,idx%by, select=bits%xl)
  call self%static_values(bz,idx%bz, select=bits%xl)
END SUBROUTINE lower_x  

SUBROUTINE upper_x(self, mem)
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

  call self%symmetric(d,idx%d, select=bits%xu)
  call self%symmetric(e,idx%e, select=bits%xu)
  call self%antisymmetric(px,idx%px, select=bits%xu)
  call self%antisymmetric(py,idx%py, select=bits%xu)
  call self%antisymmetric(pz,idx%pz, select=bits%xu)
  call self%extrapolate(bx,idx%bx, select=bits%xu)
  call self%extrapolate(by,idx%by, select=bits%xu)
  call self%extrapolate(bz,idx%bz, select=bits%xu)
  
END SUBROUTINE upper_x

SUBROUTINE lower_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), pointer:: ux, uy, uz
  real, dimension(:,:,:), allocatable:: xdnd, ydnd, zdnd
  integer:: l(3),u(3),li,i,j,k
  !-----------------------------------------------------------------------------
  l=self%mesh%lb
  u=self%mesh%ub
  li=self%mesh(2)%li
  u(2)=self%mesh(2)%lo
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  allocate(xdnd,ydnd,zdnd, mold=d)
  allocate(ux,uy,uz, mold=d)
  xdnd=0; ydnd=0; zdnd=0
  ! Convert to velocity
  xdnd = xdn(d); ydnd = ydn(d); zdnd = zdn(d)
  ux = px/xdnd; uy = py/ydnd; uz = pz/zdnd
!  call self%antisymmetric(ux,idx%px, select=bits%yl)
  call self%vanishing(ux,idx%px, select=bits%yl)
  call self%vanishing(uy,idx%py, select=bits%yl)
  call self%vanishing(uz,idx%pz, select=bits%yl)
  ! Back to momentum
  call self%static_values(d ,idx%d , select=bits%yl)
  call self%static_values(e ,idx%e , select=bits%yl)
  call self%static_values(bx,idx%bx, select=bits%yl)
  call self%static_values(by,idx%by, select=bits%yl)
  call self%static_values(bz,idx%bz, select=bits%yl)
  ! down-stagger updated density
  xdnd = xdn(d); ydnd = ydn(d); zdnd = zdn(d)
  ! Back to momentum
  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        px(i,j,k) = ux(i,j,k)*xdnd(i,j,k)
        py(i,j,k) = uy(i,j,k)*ydnd(i,j,k)
        pz(i,j,k) = uz(i,j,k)*zdnd(i,j,k)
      end do
    end do
  end do

  deallocate(xdnd,ydnd,zdnd)
  deallocate(ux,uy,uz)

END SUBROUTINE lower_y

SUBROUTINE upper_y(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), pointer:: ux
  real, dimension(:,:,:), allocatable:: xdnd
  integer:: l(3),u(3),li,i,j,k
  !-----------------------------------------------------------------------------
  l=self%mesh%lb
  u=self%mesh%ub
  l(2)=self%mesh(2)%uo
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)

  allocate(ux, xdnd, mold=d)
  xdnd=0
  ! Convert to velocity
  xdnd = xdn(d)
  ux = px/xdnd
!  call self%add_flow(ux, select=bits%yu)

!  call self%symmetric(d,idx%d, select=bits%yu)
!  call self%symmetric(e,idx%e, select=bits%yu)
  call self%antisymmetric(px,idx%px, select=bits%yu)

!  ! Convert to velocity
!  xdnd = xdn(d)
!  ux = px/xdnd
  call self%add_flow(ux, select=bits%yu)

  call self%symmetric(d,idx%d, select=bits%yu)
  call self%symmetric(e,idx%e, select=bits%yu)

  call self%antisymmetric(py,idx%py, select=bits%yu)
  call self%antisymmetric(pz,idx%pz, select=bits%yu)
  call self%zero_gradient(bx,idx%bx, select=bits%yu)
!  call self%symmetric(bx,idx%bx, select=bits%yu)
  call self%symmetric(by,idx%by, select=bits%yu)
  call self%symmetric(bz,idx%bz, select=bits%yu)

  ! down-stagger updated denxsity
  xdnd = xdn(d)
  ! Back to momentum
  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        px(i,j,k) = ux(i,j,k)*xdnd(i,j,k)
      end do
    end do
  end do
  deallocate(ux,xdnd)
END SUBROUTINE upper_y

SUBROUTINE lower_z(self, mem)
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
  
  call self%antisymmetric(px,idx%px, select=bits%zl)
  call self%antisymmetric(py,idx%py, select=bits%zl)
  call self%antisymmetric(pz,idx%pz, select=bits%zl)
  call self%static_values(d ,idx%d , select=bits%zl)
  call self%static_values(e ,idx%e , select=bits%zl)
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
  call self%antisymmetric(py,idx%py, select=bits%zu)
  call self%antisymmetric(pz,idx%pz, select=bits%zu)
  call self%extrapolate(bx,idx%bx, select=bits%zu)
  call self%extrapolate(by,idx%by, select=bits%zu)
  call self%extrapolate(bz,idx%bz, select=bits%zu)

END SUBROUTINE upper_z
!===============================================================================
!> Static boundary condition
!===============================================================================
SUBROUTINE static_values(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  real(kind=KindScalarVar), dimension(:,:,:,:), allocatable:: ff
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: Ax2, Ax2_xdnyup, Ax2_xupydn
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: Ax0, Ax0_xdnyup, Ax0_xupydn
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: A_vec, B_vec
  integer :: lb, lo, uo, ub, iv, ix, iy, iz, n
  real :: bb, d, alpha, beta, d0, p_smoother, s_smoother, Ax0_sum
  real :: bx_sum, by_sum1, by_sum2 
  real :: xp, xs, yp, ys, zp, zs
  real :: bx, by, bz, bxp, byp, bzp, pg, pb
  real(8):: ds(3), pi = asin(1.0)*2
  integer, optional :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  ds = self%patch%ds
  allocate(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%lo,m(3)%lb:m(3)%ub,8))

  bb = self%bb/sqrt(4*pi); d = self%d
  beta = self%beta; alpha = self%alpha; d0 = self%d0

  if (self%vecpot) then
  allocate(Ax2, Ax2_xdnyup, Ax2_xupydn, mold=ff(:,:,:,1))
  allocate(Ax0, Ax0_xdnyup, Ax0_xupydn, mold=ff(:,:,:,1))
  allocate(A_vec, B_vec, mold=ff(:,:,:,1:3))
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + m(3)%r(iz)
    zs = zp - 0.5*m(3)%d
    do iy=m(2)%lb,m(2)%lo+8
      yp = m(2)%p + m(2)%r(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + m(1)%r(ix)
        xs = xp - 0.5*m(1)%d
        Ax2(ix,iy,iz) = bb*(xp-2*alpha*xp)/(xp**2+yp**2)
        Ax0_sum = 0.0
        do n=1,25
          Ax0_sum = Ax0_sum + (1/n**2)*(1-1/n**2*(xp**2+yp**2)/  &
                    d**2)/((1-1/n**2*(xp**2+yp**2)/d**2)**2+     &
                    (1/n*2*yp/d)**2)
        end do
        Ax0(ix,iy,iz) = alpha*bb*xp*(2/d)**2*Ax0_sum
        !bz = 0.0!self%alpha*bb*xp*(2/d)**2*bz_sum
        !s_smoother = abs(tanh(xs+d/2)+tanh(xs-d/2))
        !p_smoother = abs(tanh(xp+d/2)+tanh(xp-d/2))
        !bx =  4*self%alpha*bb*xs*yp/(xs**2+yp**2)**2
        !by = -bb*(1+2*self%alpha*(xp**2-ys**2)/(xp**2+ys**2)**2)
        !bzp = 0.0
        !bxp = 4*self%alpha*bb*xp*yp/(xp**2+yp**2)**2
        !byp = -bb*(1+2*self%alpha*(xp**2-yp**2)/(xp**2+yp**2)**2)
        ff(ix,iy,iz,idx%d)=self%d0
        ff(ix,iy,iz,idx%px)=0.0
        ff(ix,iy,iz,idx%py)=0.0
        ff(ix,iy,iz,idx%pz)=0.0
        !ff(ix,iy,iz,idx%bx)=bx*(1-s_smoother)
        !ff(ix,iy,iz,idx%by)=by
        !ff(ix,iy,iz,idx%bz)=bz
        !bxp = bxp*(1-p_smoother)
        !pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = 0.01
        !pg = self%beta*bb-pb
        !pg = 30 -pb
        if (idx%e>0) then
          ff(ix,iy,iz,idx%e) = pg / (eos%gamma - 1.0)
          if (self%patch%kind(1:6) == 'ramses')  then
            ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + pb
          end if
        else if (idx%s>0) then
          if (eos%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else
            ff(ix,iy,iz,idx%s) = &
                (log(pg)-log(self%d0)*eos%gamma)/(eos%gamma-1.0)*self%d0
          end if
        end if
      end do
    end do
  end do
  A_vec(:,:,:,1) = Ax0
  A_vec(:,:,:,2) = 0.0
  A_vec(:,:,:,3) = Ax2
  B_vec = curl_dn(ds, A_vec)
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%lo
      do ix=m(1)%lb,m(1)%ub 
        ff(ix,iy,iz,idx%bx) = B_vec(ix,iy,iz,1)
        ff(ix,iy,iz,idx%by) = B_vec(ix,iy,iz,2)
        ff(ix,iy,iz,idx%bz) = B_vec(ix,iy,iz,3)
      end do
    end do
  end do
  deallocate(Ax2, Ax2_xdnyup, Ax2_xupydn)
  deallocate(Ax0, Ax0_xdnyup, Ax0_xupydn)
  deallocate(A_vec, B_vec)
  else

  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + m(3)%r(iz)
    zs = zp - 0.5*m(3)%d
    do iy=m(2)%lb,m(2)%lo
      yp = m(2)%p + m(2)%r(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + m(1)%r(ix)
        xs = xp - 0.5*m(1)%d
        bx_sum = 0.
        by_sum1 = 0.
        by_sum2 = 0.
        do n=1,25
          bx_sum = bx_sum + &
                 ((1 - (xs**2 + yp**2)/(d**2*n**2))*(4*yp*(1-(xs**2+yp**2)/ &
                 (d**2*n**2))/(d**2*n**2) - 8*yp/(d**2*n**2))/(n**2*((1 -   &
                 (xs**2 + yp**2)/(d**2*n**2))**2+4*yp**2/(d**2*n**2))**2) - &
                 2*yp/(d**2*n**4*((1 - (xs**2+yp**2)/(d**2*n**2))**2 +      &
                 4*yp**2/(d**2*n**2))))
          by_sum1 = by_sum1 + &
                 (4*xp*(1 -(xp**2 + ys**2)/(d**2*n**2))**2/(d**2*n**4*((1 - &
                 (xp**2 + ys**2)/(d**2*n**2))**2+4*ys**2/(d**2*n**2))**2) - &
                 2*xp/(d**2*n**4*((1 - (xp**2 + ys**2)/(d**2*n**2))**2 +    &
                 4*ys**2/(d**2*n**2))))
          by_sum2 = by_sum2 + &
                 ((1 - (xp**2 + ys**2)/(d**2*n**2))/(n**2*((1 - (xp**2 +    &
                 ys**2)/(d**2*n**2))**2 + 4*ys**2/(d**2*n**2))))
        end do
!        bz_sum = 0.
!        do n=1,25
!          bz_sum=bz_sum+(((1-(xp**2+yp**2)/(d**2*n**2))*((8*yp)/(d**2*n**2)- &
!                        (4*yp*(1-(xp**2+yp**2)/(d**2*n**2)))/(d**2*n**2)))/  &
!                        (n**2*((4*yp**2)/(d**2*n**2)+(1-(xp**2+yp**2)/       &
!                        (d**2*n**2))**2)**2))+(2*yp)/(d**2*n**4*((4*yp**2)/  &
!                        (d**2*n**2)+(1-(xp**2+yp**2)/(d**2*n**2))**2))
!        end do
        bz = 0.0
!        s_smoother = abs(tanh(xs+d/2)+tanh(xs-d/2))
!        p_smoother = abs(tanh(xp+d/2)+tanh(xp-d/2))
        bx = 4*self%alpha*bb*xs*yp/(xs**2+yp**2)**2 + &
             4*bb*self%alpha*xs*bx_sum/d**2
        by = -bb*(1+2*self%alpha*(xp**2-ys**2)/(xp**2+ys**2)**2) - &
             4*bb*self%alpha*xp*by_sum1/d**2-4*bb*self%alpha*by_sum2/d**2
        bzp = 0.0
        bxp = 4*self%alpha*bb*xp*yp/(xp**2+yp**2)**2
        byp = -bb*(1+2*self%alpha*(xp**2-yp**2)/(xp**2+yp**2)**2)
        ff(ix,iy,iz,idx%d)=self%d0
        ff(ix,iy,iz,idx%px)=0.0
        ff(ix,iy,iz,idx%py)=0.0
        ff(ix,iy,iz,idx%pz)=0.0
        ff(ix,iy,iz,idx%bx)=bx!*(1-s_smoother)
        ff(ix,iy,iz,idx%by)=by
        ff(ix,iy,iz,idx%bz)=bz
        !bxp = bxp*(1-p_smoother)
        pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = 0.01
        !pg = self%beta*bb-pb
        !pg = 30 -pb
!        s_smoother = abs(tanh(xs+d/2)+tanh(xs-d/2))
!        p_smoother = abs(tanh(xp+d/2)+tanh(xp-d/2))
!        bz = 0.0
!        bx = 4*alpha*bb*xs*yp/(xs**2+yp**2)**2
!        by = -bb*(1+2*alpha*(xp**2-ys**2)/(xp**2+ys**2)**2)
!        bzp = 0.0
!        bxp = 4*alpha*bb*xp*yp/(xp**2+yp**2)**2
!        byp = -bb*(1+2*alpha*(xp**2-yp**2)/(xp**2+yp**2)**2)
!        ff(ix,iy,iz,1)=d0
!        ff(ix,iy,iz,2)=0.0
!        ff(ix,iy,iz,3)=0.0
!        ff(ix,iy,iz,4)=0.0
!        ff(ix,iy,iz,6)=bx*(1-s_smoother)
!        ff(ix,iy,iz,7)=by
!        ff(ix,iy,iz,8)=bz
!        bxp = bxp*(1-p_smoother)
!        pb = 0.5*(bxp**2+byp**2+bzp**2) 
!        pg = 0.01 !beta*bb
        if (idx%e>0) then 
          ff(ix,iy,iz,5) = pg / (eos%gamma - 1.0)
          if (self%patch%kind(1:6) == 'ramses')  then
            ff(ix,iy,iz,5) = ff(ix,iy,iz,5) + pb
          end if
        else if (idx%s>0) then
          if (eos%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else
            ff(ix,iy,iz,idx%s) = &
                (log(pg)-log(d0)*eos%gamma)/(eos%gamma-1.0)*d0
          end if
        end if
      end do
    end do
  end do
  end if
  mem(:,m(2)%lb:m(2)%lo,:) = ff(:,:,:,iv)
  deallocate(ff)
END SUBROUTINE static_values
!===============================================================================
!> Vanishing boundary condition, more gentle than antisymmetric
!===============================================================================
SUBROUTINE vanishing(self, mem, iv, select)
  class(boundaries_t)      :: self
  class(mesh_t), pointer   :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: mem
  integer :: lb, lo, uo, ub, iv, ix, iy, iz
  integer, optional :: select
  real :: factor
  !-----------------------------------------------------------------------------
  m => self%mesh
  factor = self%vanishing_factor
  if (self%is_set(bits%xl)) then
    lb = m(1)%lb
    lo = m(1)%lo
    do ix=lo,lb,-1
      mem(ix,:,:) = mem(ix+1,:,:)*factor
    end do
    self%was_filled(1,1) = .true.
  else if (self%is_set(bits%xu)) then
    ub = m(1)%ub
    uo = m(1)%uo
    do ix=uo,ub
      mem(ix,:,:) = mem(ix-1,:,:)*factor
    end do
    self%was_filled(1,2) = .true.
  else if (self%is_set(bits%yl)) then
    lb = m(2)%lb
    lo = m(2)%lo
    do iy=lo,lb,-1
      mem(:,iy,:) = mem(:,iy+1,:)*factor
    end do
    self%was_filled(2,1) = .true.
  else if (self%is_set(bits%yu)) then
    ub = m(2)%ub
    uo = m(2)%uo
    do iy=uo,ub
      mem(:,iy,:) = mem(:,iy-1,:)*factor
    end do
    self%was_filled(2,2) = .true.
  else if (self%is_set(bits%zl)) then
    lb = m(3)%lb
    lo = m(3)%lo
    do iz=lo,lb,-1
      mem(:,:,iz) = mem(:,:,iz+1)*factor
    end do
    self%was_filled(3,1) = .true.
  else if (self%is_set(bits%zu)) then
    ub = m(3)%ub
    uo = m(3)%uo
    do iz=uo,ub
      mem(:,:,iz) = mem(:,:,iz-1)*factor
    end do
    self%was_filled(3,2) = .true.
  end if
END SUBROUTINE vanishing
!===============================================================================
!> Add a flow to the boundary
!===============================================================================
SUBROUTINE add_flow(self, mem, select)
  class(boundaries_t)     :: self
  class(mesh_t), pointer  :: m(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer                 :: ub, uo, lb, lo
  real                    :: t_start, t_end, v0, flow, t
  integer, optional       :: select
  !-----------------------------------------------------------------------------
  m => self%mesh
  t = self%patch%time
  t_end = self%t0
  t_start = t_end-100
  v0 = self%v0
  if (t<t_start) then
    flow = 0.0
  else
    if (t<t_end) then
      flow = v0/2.*(1.-cos(math%pi*(t-t_start)/(t_end-t_start)))
    else
      flow = v0
    end if
  end if

  if (self%is_set(bits%yl)) then
    lo = m(2)%lo
    lb = m(2)%lb
    mem(:,lb:lo,:) = -flow
  else if (self%is_set(bits%yu)) then
    uo = m(2)%uo
    ub = m(2)%ub
    mem(:,uo:ub,:) = flow
  end if
END SUBROUTINE add_flow

END MODULE boundaries_mod
