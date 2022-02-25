!===============================================================================
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE eos_mod
  USE gpatch_mod
  USE vector_ops
  implicit none
  private
  type, public :: initial_conditions_t
    integer :: idir
    logical :: mhd=.true., vecpot=.false.
    real    :: d0=1.0, beta=0.1, alpha=0.2, d=20.0, bb=3.
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
  end type
CONTAINS

!===============================================================================
!> Pre-initialize the initial conditions.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer,save:: idir
  real, save :: d0, beta, alpha, d, bb
  logical, save :: first_time=.true., vecpot=.false.
  namelist /IC_params/ d0, beta, alpha, d, bb, idir, vecpot
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_conditions_t%init')
  self%patch => patch

  call io%lock%set('initial_conditions_t%pre_init')
  if (first_time) then
    ! default values for namelist parameters
!    d0    = 1.0
!   beta  = 0.1
!    alpha = 0.2
!    d     = 20.0
!    bb    = 3.0
!    idir  = 3
    rewind (io%input)
    read (io%input, IC_params)
    write (io%namelist, IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_conditions_t%pre_init')

  self%d0     = d0
  self%beta   = beta
  self%alpha  = alpha
  self%d      = d
  self%bb     = bb
  self%idir   = idir
  self%vecpot = vecpot
  !
  call trace%end()
END subroutine pre_init
!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  ! ----------------------------------------------------------------------------
  type(index_t), pointer:: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: ff, ff_copy, A_vec, B_vec
  !.............................................................................
  integer:: ix, iy, iz, n                                  ! loop indices
  real:: entropy, pg, pb, bb, alpha, bz_sum, d, bx, by, bz, bxp, byp, bzp
  real, dimension(:,:,:), pointer:: Ax2, Ax2_xdnyup, Ax2_xupydn
  real, dimension(:,:,:), pointer:: Ax0, Ax0_xdnyup, Ax0_xupydn
  real:: s_smoother, p_smoother, Ax0_sum, bx_sum, by_sum1, by_sum2
  real(8):: ds(3), xs, ys, zs, xp, yp, zp
  real(8):: pi = asin(1.0)*2
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')

  m => self%patch%mesh
  ds = self%patch%ds
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)

  alpha = self%alpha
  bb = self%bb/sqrt(4*pi)
  d = self%d

  if (self%vecpot) then
  allocate(Ax2, Ax2_xdnyup, Ax2_xupydn, mold=ff(:,:,:,1))
  allocate(Ax0, Ax0_xdnyup, Ax0_xupydn, mold=ff(:,:,:,1))
  allocate(A_vec, B_vec, mold=ff(:,:,:,1:3))
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz)
    zs = zp - 0.5*m(3)%d
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix)
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
  ff(:,:,:,idx%bx) = B_vec(:,:,:,1)
  ff(:,:,:,idx%by) = B_vec(:,:,:,2)
  ff(:,:,:,idx%bz) = B_vec(:,:,:,3)
  deallocate(Ax2, Ax2_xdnyup, Ax2_xupydn)
  deallocate(Ax0, Ax0_xdnyup, Ax0_xupydn)
  deallocate(A_vec, B_vec)
  else

  if (self%idir==1) then ! dimensions are YZ
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
        s_smoother = abs(tanh(ys+d/2)+tanh(ys-d/2))
        p_smoother = abs(tanh(yp+d/2)+tanh(yp-d/2))
        bx = 0.0
        by = 4*self%alpha*bb*ys*zp/(ys**2+zp**2)**2
        bz = -bb*(1+2*self%alpha*(yp**2-zs**2)/(yp**2+zs**2)**2)
        bxp = 0.0
        byp = 4*self%alpha*bb*yp*zp/(yp**2+zp**2)**2
        bzp = -bb*(1+2*self%alpha*(yp**2-zp**2)/(yp**2+zp**2)**2)
        ff(ix,iy,iz,idx%d)=self%d0
        ff(ix,iy,iz,idx%px)=0.0
        ff(ix,iy,iz,idx%py)=0.0
        ff(ix,iy,iz,idx%pz)=0.0
        ff(ix,iy,iz,idx%bx)=bx
        ff(ix,iy,iz,idx%by)=by*(1-s_smoother)
        ff(ix,iy,iz,idx%bz)=bz
        byp = byp*(1-p_smoother)
        pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = self%beta*bb
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

  else if (self%idir==2) then ! dimensions are XZ
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
        s_smoother = abs(tanh(zs+d/2)+tanh(zs-d/2))
        p_smoother = abs(tanh(zp+d/2)+tanh(zp-d/2))
        by = 0.0
        bz = 4*self%alpha*bb*zs*xp/(zs**2+xp**2)**2
        bx = -bb*(1+2*self%alpha*(zp**2-xs**2)/(zp**2+xs**2)**2)
        byp = 0.0
        bxp = 4*self%alpha*bb*zp*xp/(zp**2+xp**2)**2
        bzp = -bb*(1+2*self%alpha*(zp**2-xp**2)/(zp**2+xp**2)**2)
        ff(ix,iy,iz,idx%d)=self%d0
        ff(ix,iy,iz,idx%px)=0.0
        ff(ix,iy,iz,idx%py)=0.0
        ff(ix,iy,iz,idx%pz)=0.0
        ff(ix,iy,iz,idx%bx)=bx
        ff(ix,iy,iz,idx%by)=by
        ff(ix,iy,iz,idx%bz)=bz*(1-s_smoother)
        bxp = bxp*(1-p_smoother)
        pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = self%beta*bb
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

  else if (self%idir==3) then ! dimensions are XY
  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz)
    zs = zp - 0.5*m(3)%d
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix)
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
  end if
  end if
  call trace_end
  end associate
END SUBROUTINE init

END MODULE initial_conditions_mod
