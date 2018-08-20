!> =============================================================================
!> Given a direction Omega, which may align mainly with the x-. y-, or z-axis,
!> solve the integral form of the RT on short characteristics, extending across
!> only three layers perpendicular to the main direction.
!>
!> This is done by transposing the incoming problem, defined by rk, src, and q,
!> in connection with doing bilinear interpolation in the planes perpendicular
!> to the main RT axis. 
!> =============================================================================
MODULE rt_integral_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  type, public:: rt_integral_t
    real:: mu(3)
  contains
    procedure:: solve
    procedure, nopass:: diagnostics
  end type
  integer:: verbose=0
  real:: dtau_min=1e-3, dtau_dif=1e3, dtau_max=1e4
  integer(8), save:: n_cases(4)=0
  type(rt_integral_t), public:: rt_integral
CONTAINS

!> =============================================================================
!> Solve the radiative transfer for all bins, along the closest axis direction
!> =============================================================================
SUBROUTINE solve (self, k, s, q, mesh, mu3, phi3)
  class(rt_integral_t):: self
  real, dimension(:,:,:,:), pointer:: k, s, q
  class(mesh_t), pointer:: mesh(:)
  real:: mu3, phi3
  !.............................................................................
  real, dimension(:,:,:,:), allocatable:: k3, s3, q3
  real:: dr
  integer:: m(4), loc(1), dir, i, i1, i2, i3, ib, l(3), u(3), ii(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_integral_t%solve')
  m = shape(k)
  !-----------------------------------------------------------------------------
  ! Determine closest axis direction
  !-----------------------------------------------------------------------------
  self%mu(1) = sin(acos(mu3))*cos(phi3)
  self%mu(2) = sin(acos(mu3))*sin(phi3)
  self%mu(3) = mu3
  loc = maxloc(abs(self%mu))
  dir = loc(1)
  dr = mesh(dir)%d/abs(self%mu(dir))
  !-----------------------------------------------------------------------------
  ! Transpose dimensions to solver index coordinates and allocate
  !-----------------------------------------------------------------------------
  if (dir==1) then
    ii = [2,3,1]
  else if (dir==2) then
    ii = [1,3,2]
  else
    ii = [1,2,3]
  end if
  m(1:3) = m(ii)
  allocate (k3(m(1),m(2),m(4),3))
  allocate (s3(m(1),m(2),m(4),3))
  allocate (q3(m(1),m(2),m(4),3))
  !-----------------------------------------------------------------------------
  ! Determine loop limits and loop direction
  !-----------------------------------------------------------------------------
  l = mesh(ii)%li
  u = mesh(ii)%ui
  if (self%mu(dir) > 0.0) then
    i1=l(3); i2=u(3); i3=+1
  else
    i1=u(3); i2=l(3); i3=-1
  end if
  if (verbose > 0) &
    print *, 'rt_integral_t%solve: mu, phi, dir, dr =', self%mu, dir, i1, i2, i3
  !-----------------------------------------------------------------------------
  ! Prepare layers for the single layer solver and call it
  !-----------------------------------------------------------------------------
  do i=i1,i2,i3
    call prepare_layers (self, dir, mesh, k, k3, i, i3)
    call prepare_layers (self, dir, mesh, s, s3, i, i3)
    call prepare_layers (self, dir, mesh, q, q3, i, i3)
    if (dir==1) then
      call solve1 (k3, s3, q3, q(i,:,:,:), dr, l, u)
    else if (dir==2) then
      call solve1 (k3, s3, q3, q(:,i,:,:), dr, l, u)
    else
      call solve1 (k3, s3, q3, q(:,:,i,:), dr, l, u)
    end if
  end do
  deallocate (k3, s3, q3)
  if (verbose > 0) then
    do ib=1,m(4)
      print '(a,i4,1p,10e10.2)', ' rt_integral_t%solve: q =', ib, &
        q(l(1):l(1)+9,(l(2)+u(2))/2,(l(3)+u(3))/2,ib)
    end do
  end if
  !call diagnostics
  call trace%end
END SUBROUTINE solve

!> =============================================================================
!> Prepare three layers; one centered on the mesh points, one before and one
!> after, in the RT direction.  The indices in the layers are (i1,i2,ib,1:3), 
!> where i1-i3 are 3D coords, ib is bin index, and the last is layer index
!> =============================================================================
SUBROUTINE prepare_layers (self, dir, mesh, in, out, i, i3)
  class(rt_integral_t):: self
  integer:: dir, i, i3
  class(mesh_t), pointer:: mesh(:)
  real, dimension(:,:,:,:):: in, out
  integer:: o(2,2), l(2), u(2)
  real:: p(2,2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_integral_t%prepare_layers', itimer=itimer)
  !-----------------------------------------------------------------------------
  call calculate_offsets (self, dir, mesh, o, p)
  if (dir==1) then
    l = [mesh(2)%li, mesh(3)%li]
    u = [mesh(2)%ui, mesh(3)%ui]
    call bilinear (in(i-i3,:,:,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(i,:,:,:)
    call bilinear (in(i+i3,:,:,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  else if (dir==2) then
    l = [mesh(1)%li, mesh(3)%li]
    u = [mesh(1)%ui, mesh(3)%ui]
    call bilinear (in(:,i-i3,:,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(:,i,:,:)
    call bilinear (in(:,i+i3,:,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  else if (dir==3) then
    l = [mesh(1)%li, mesh(2)%li]
    u = [mesh(1)%ui, mesh(2)%ui]
    call bilinear (in(:,:,i-i3,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(:,:,i,:)
    call bilinear (in(:,:,i+i3,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  end if
  if (verbose > 1) &
    print *, 'rt_integral_t%prepare_layers: dir, l, u =', dir, l, u
  call trace%end (itimer)
END SUBROUTINE prepare_layers

!> =============================================================================
!> Calculate perpendicular index offsets from one layer to the next
!> =============================================================================
SUBROUTINE calculate_offsets (self, dir, mesh, o, p)
  class(rt_integral_t):: self
  class(mesh_t), pointer:: mesh(:)
  integer:: dir, o(2,2)
  real:: p(2,2)
  !.............................................................................
  real(8), pointer, dimension(:):: x, y, z
  real:: dpdi(2)
  integer:: j
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_integral_t%calculate_offsets')
  ! -- Fractional steps in perpendicular indices per parallel index
  if (dir==1) then
    dpdi(1) = self%mu(2)/self%mu(1)*mesh(1)%d/mesh(2)%d
    dpdi(2) = self%mu(3)/self%mu(1)*mesh(1)%d/mesh(3)%d
  else if (dir==2) then
    dpdi(1) = self%mu(1)/self%mu(2)*mesh(2)%d/mesh(1)%d
    dpdi(2) = self%mu(3)/self%mu(2)*mesh(2)%d/mesh(3)%d
  else
    dpdi(1) = self%mu(1)/self%mu(3)*mesh(3)%d/mesh(1)%d
    dpdi(2) = self%mu(2)/self%mu(3)*mesh(3)%d/mesh(2)%d
  end if
  ! -- backward and forwards one index
  do j=1,2
    p(:,j) = (2*j-3)*dpdi
  end do
  o = floor(p)
  p = p - o
  if (verbose > 1) &
    print '(a,2x,3f7.3,2x,2f7.3,2(2x,2i3),0p,2(2x,2f7.3))', &
      ' rt_integral_t%calculate_offsets: mu, dpdi, o, p =', &
      self%mu, dpdi, o, p
  call trace%end
END SUBROUTINE calculate_offsets

!> =============================================================================
!> Bilinear interpolation in a plane
!> =============================================================================
SUBROUTINE bilinear (in, out, o, p, l, u)
  real, dimension(:,:,:), intent(in)  :: in
  real, dimension(:,:,:), intent(out) :: out
  real, dimension(2):: p, q
  integer, dimension(2):: o, l, u
  !.............................................................................
  integer:: i1, i2, ib, m(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_integral_t%bilinear', 2, itimer=itimer)
  m = shape(in)
  q = 1.0-p
  !-----------------------------------------------------------------------------
  ! Purely vertical case
  !-----------------------------------------------------------------------------
  if (all(abs(o+p) < 1e-5)) then
    out = in
    if (verbose > 1) &
      print *, 'rt_integral_t%bilinear: no interpolation'
  !-----------------------------------------------------------------------------
  ! Interpolation only in y
  !-----------------------------------------------------------------------------
  else if (abs(o(1)+p(1)) < 1e-5) then
    if (verbose > 1) &
      print *, 'rt_integral_t%bilinear: 2nd index, o, p =', &
        o(2), p(2)
    do ib=1,m(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,ib) = &
        q(2)*in(i1,i2+o(2),ib) + p(2)*in(i1,i2+o(2)+1,ib)
    end do
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Interpolation only in x
  !-----------------------------------------------------------------------------
  else if (abs(o(2)+p(2)) < 1e-5) then
    if (verbose > 1) &
      print *, 'rt_integral_t%bilinear: 1st index, o, p =', &
        o(1), p(1)
    do ib=1,m(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,ib) = &
        q(1)*in(i1+o(1),i2,ib) + p(1)*in(i1+o(1)+1,i2,ib)
    end do
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Interpolation in x and y
  !-----------------------------------------------------------------------------
  else
    if (verbose > 1) &
      print *, 'rt_integral_t%bilinear: o, p =', &
        o, p
    do ib=1,m(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,ib) = &
          q(1)*q(2)*in(i1+o(1)  ,i2+o(2)  ,ib) + &
          p(1)*q(2)*in(i1+o(1)+1,i2+o(2)  ,ib) + &
          q(1)*p(2)*in(i1+o(1)+1,i2+o(2)+1,ib) + &
          p(1)*p(2)*in(i1+o(1)+1,i2+o(2)+1,ib)
    end do
    end do
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE bilinear

!> =============================================================================
!> Integral solver for one layer.  The exact integral solution Q=I-S for a 2nd
!> order polynomial has coefficient for the upstream Q, and for the 1st and 2nd
!> derivative of S.  However, if one is doing a solution both forwards and 
!> backward, the term in the 1st derivative cancels.
!>
!> The order of indices is in general permuted in the call, with 
!> =============================================================================
SUBROUTINE solve1 (k, s, qin, qout, dr, l, u)
  real, dimension(:,:,:,:) :: k, s, qin
  real, dimension(:,:,:) :: qout
  integer:: l(3), u(3)
  real:: dr
  !.............................................................................
  integer:: i1, i2, ib, m(4)
  real:: dsdtau1, dsdtau2, d2sdtau2, ex0, ex2, dtaumin, dtaumax
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  real, dimension(size(k,1)):: dtau1, dtau2
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_integral_t%solve1', 1, itimer=itimer)
  m = shape(k)
  do ib=1,m(4)
    do i2=l(2),u(2)
      dtau1(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,1))
      dtau2(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,3))
      dtaumin = minval(dtau1(l(1):u(1)))
      dtaumax = maxval(dtau1(l(1):u(1)))
      if (dtaumin > dtau_max) then
        do i1=l(1),u(1)
          qout(i1,i2,ib) = 0.0
        end do
        !$omp atomic
        n_cases(4) = n_cases(4)+1
      else if (dtaumin > dtau_dif) then
        do i1=l(1),u(1)
          dsdtau1  =     (s(i1,i2,ib,2) - s(i1,i2,ib,1))/dtau1(i1)
          dsdtau2  =     (s(i1,i2,ib,3) - s(i1,i2,ib,2))/dtau2(i1)
          d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
          qout(i1,i2,ib) = d2sdtau2
        end do
        !$omp atomic
        n_cases(3) = n_cases(3)+1
      else if (dtaumax < dtau_min) then
        do i1=l(1),u(1)
          ex0      = 1.0-dtau1(i1)
          qout(i1,i2,ib) = qin(i1,i2,ib,1)*ex0
        end do
        !$omp atomic
        n_cases(1) = n_cases(1)+1
      else 
        do i1=l(1),u(1)
          dsdtau1  =     (s(i1,i2,ib,2) - s(i1,i2,ib,1))/dtau1(i1)
          dsdtau2  =     (s(i1,i2,ib,3) - s(i1,i2,ib,2))/dtau2(i1)
          d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
          ex0      = exp(-dtau1(i1))
          ex2      = (1.0-ex0) - dtau1(i1)*ex0
          qout(i1,i2,ib) = qin(i1,i2,1,ib)*ex0 &
                         + d2sdtau2       *ex2
        end do
        !$omp atomic
        n_cases(2) = n_cases(2)+1
      end if
    end do
    if (verbose > 1) then
      print '(a,i4,1p,10e10.2)', ' rt_integral_t%solve1: ib, q =', ib, &
        qout(l(1):l(1)+9,(l(2)+u(2))/2,ib)
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE solve1

!> =============================================================================
!> =============================================================================
SUBROUTINE diagnostics (label)
  character(len=*):: label
  real:: f
  f = 1./max(sum(n_cases),1_8)
  print *, 'fractions of thin, full, diff, and zero cases: ', &
    trim(label), n_cases*f
END SUBROUTINE diagnostics 

END MODULE rt_integral_mod
