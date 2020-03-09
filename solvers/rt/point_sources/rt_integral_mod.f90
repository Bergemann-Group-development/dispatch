!===============================================================================
!> Given a direction Omega, which may align mainly with the x-. y-, or z-axis,
!> solve the integral form of the RT on short characteristics, extending across
!> only three layers perpendicular to the main direction.
!>
!> This is done by transposing the incoming problem, defined by rk, src, and q,
!> in connection with doing bilinear interpolation in the planes perpendicular
!> to the main RT axis. 
!===============================================================================
MODULE rt_integral_mod
  USE timer_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  implicit none
  private
  type, public:: rt_integral_t
    procedure(solve1), pointer:: solver => null()
  contains
    procedure:: init
    procedure:: solve
    procedure:: prepare_layers
    procedure:: bilinear
    procedure:: calculate_offsets
    procedure:: solve1
    procedure:: solve2
    procedure, nopass:: diagnostics
  end type
  integer, save:: verbose=0, order=3, solver=2
  logical:: detailed_timer=.false.
  real:: dtau_min=1e-3, dtau_dif=1e3, dtau_max=1e4
  integer(8), save:: n_cases(4)=0
  type(rt_integral_t), public:: rt_integral
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(rt_integral_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /rt_integral_params/ verbose, solver, order, dtau_min, dtau_dif, &
    dtau_max, detailed_timer
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time=.false.
      rewind (io_unit%input)
      read (io_unit%input, rt_integral_params, iostat=iostat)
      write (io_unit%output, rt_integral_params)
    end if
    !$omp end critical (input_cr)
    select case (solver)
    case (1)
      self%solver => solve1
    case (2)
      self%solver => solve2
    case default
      write (stderr,*) solver
      flush (stderr)
      call io%abort ('rt_integral%init: wrong solver')
    end select
  end if
END SUBROUTINE init

!===============================================================================
!> Solve the radiative transfer for all bins, along the closest axis direction
!===============================================================================
SUBROUTINE solve (self, k, s, q, mesh, dir, mu)
  class(rt_integral_t):: self
  real, dimension(:,:,:,:), pointer:: k, s, q
  class(mesh_t), pointer:: mesh(:)
  real:: mu3, phi3, mu(3)
  !.............................................................................
  real, dimension(:,:,:,:), allocatable:: k3, s3, q3
  real:: dr
  logical:: debug
  integer:: m(4), loc(1), dir, i, i1, i2, i3, ib, l(3), u(3), ii(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trc_begin ('rt_integral_t%solve', itimer=itimer)
  debug = .false.
  m = shape(k)
  !-----------------------------------------------------------------------------
  dr = 0.5*mesh(dir)%d/abs(mu(dir))
  !-----------------------------------------------------------------------------
  ! Transpose dimensions to solver index coordinates and allocate
  !-----------------------------------------------------------------------------
  if (dir==1) then
    ii = [2,3,1]
  else if (dir==2) then
    ii = [3,1,2]
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
  if (mu(dir) > 0.0) then
    i1=l(3); i2=u(3); i3=+1
  else
    i1=u(3); i2=l(3); i3=-1
  end if
  if (verbose > 1) &
    write(io_unit%log,*) 'rt_integral_t%solve: mu, phi, dir, dr =', mu, dir, i1, i2, i3
  !-----------------------------------------------------------------------------
  ! Prepare layers for the single layer solver and call it
  !-----------------------------------------------------------------------------
  do i=i1,i2,i3
    call self%prepare_layers (mu, dir, mesh, k, k3, i, i3)
    call self%prepare_layers (mu, dir, mesh, s, s3, i, i3)
    call self%prepare_layers (mu, dir, mesh, q, q3, i, i3)
    if (dir==1) then
      call self%solver (k3, s3, q3, q(i,:,:,:), dr, l, u, debug)
    else if (dir==2) then
      call self%solver (k3, s3, q3, q(:,i,:,:), dr, l, u, debug)
    else
      call self%solver (k3, s3, q3, q(:,:,i,:), dr, l, u, debug)
    end if
  end do
  deallocate (k3, s3, q3)
  if (verbose > 1) then
    do ib=1,m(4)
      write(io_unit%log,'(a,i4,1p,10e10.2)') ' rt_integral_t%solve: q =', &
        ib, q(l(1):l(1)+9,(l(2)+u(2))/2,(l(3)+u(3))/2,ib)
    end do
  end if
  if (verbose > 0) &
    call diagnostics
  !-----------------------------------------------------------------------------
  call trc_end (itimer)
END SUBROUTINE solve

!===============================================================================
!> Prepare three layers; one centered on the mesh points, one before and one
!> after, in the RT direction.  The indices in the layers are (i1,i2,ib,1:3), 
!> where i1-i3 are 3D coords, ib is bin index, and the last is layer index
!===============================================================================
SUBROUTINE prepare_layers (self, mu, dir, mesh, in, out, i, i3)
  class(rt_integral_t):: self
  real:: mu(3)
  integer:: dir, i, i3
  class(mesh_t), pointer:: mesh(:)
  real, dimension(:,:,:,:):: in, out
  integer:: o(2,2), l(2), u(2)
  real:: p(2,2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trc_begin ('rt_integral_t%prepare_layers', itimer=itimer)
  !-----------------------------------------------------------------------------
  call self%calculate_offsets (mu, dir, mesh, o, p)
  if (dir==1) then
    l = [mesh(2)%li, mesh(3)%li]
    u = [mesh(2)%ui, mesh(3)%ui]
    call self%bilinear (in(i-i3,:,:,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(i,:,:,:)
    call self%bilinear (in(i+i3,:,:,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  else if (dir==2) then
    l = [mesh(1)%li, mesh(3)%li]
    u = [mesh(1)%ui, mesh(3)%ui]
    call self%bilinear (in(:,i-i3,:,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(:,i,:,:)
    call self%bilinear (in(:,i+i3,:,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  else if (dir==3) then
    l = [mesh(1)%li, mesh(2)%li]
    u = [mesh(1)%ui, mesh(2)%ui]
    call self%bilinear (in(:,:,i-i3,:), out(:,:,:,1), o(:,1), p(:,1), l, u)
    out(:,:,:,2) = in(:,:,i,:)
    call self%bilinear (in(:,:,i+i3,:), out(:,:,:,3), o(:,2), p(:,2), l, u)
  end if
  if (verbose > 1) &
    write(io_unit%log,*) 'rt_integral_t%prepare_layers: dir, l, u =', dir, l, u
  call trc_end (itimer)
END SUBROUTINE prepare_layers

!===============================================================================
!> Calculate perpendicular index offsets from one layer to the next
!===============================================================================
SUBROUTINE calculate_offsets (self, mu, dir, mesh, o, p)
  class(rt_integral_t):: self
  class(mesh_t), pointer:: mesh(:)
  integer:: dir, o(2,2)
  real:: mu(3), p(2,2)
  !.............................................................................
  real(8), pointer, dimension(:):: x, y, z
  real:: dpdi(2)
  integer:: j
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! -- Fractional steps in perpendicular indices per parallel index
  if (dir==1) then
    ! -- when dir==1, the 1st index is y, the 2nd is z
    dpdi(1) = mu(2)/abs(mu(1))*mesh(1)%d/mesh(2)%d
    dpdi(2) = mu(3)/abs(mu(1))*mesh(1)%d/mesh(3)%d
  else if (dir==2) then
    ! -- when dir==2, the 1st index is z, the 2nd is x
    dpdi(1) = mu(3)/abs(mu(2))*mesh(2)%d/mesh(3)%d
    dpdi(2) = mu(1)/abs(mu(2))*mesh(2)%d/mesh(1)%d
  else
    ! -- when dir==3, the 1st index is x, the 2nd is y
    dpdi(1) = mu(1)/abs(mu(3))*mesh(3)%d/mesh(1)%d
    dpdi(2) = mu(2)/abs(mu(3))*mesh(3)%d/mesh(2)%d
  end if
  ! -- backward and forwards one index
  do j=1,2
    p(:,j) = (2*j-3)*dpdi
  end do
  ! -- prevent o=1 p=0.0 from happening
  o = max(-1,min(0,floor(p)))
  p = p - o
  if (verbose > 1) then
    write(io_unit%log,'(i2,2x,a,2x,3f7.3,2x,2f7.3,2(2x,2i3),0p,2(2x,2f7.3))') &
      omp%thread, 'rt_integral_t%calculate_offsets: mu, dpdi, o, p =', &
      mu, dpdi, o, p
  end if
END SUBROUTINE calculate_offsets

!===============================================================================
!> Bilinear interpolation in a plane
!===============================================================================
SUBROUTINE bilinear (self, in, out, o, p, l, u)
  class(rt_integral_t):: self
  real, dimension(:,:,:), intent(in)  :: in
  real, dimension(:,:,:), intent(out) :: out
  real, dimension(2):: p, q
  integer, dimension(2):: o, l, u
  !.............................................................................
  integer:: i1, i2, ib, m(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trc_begin ('rt_integral_t%bilinear', itimer=itimer)
  m = shape(in)
  q = 1.0-p
  !-----------------------------------------------------------------------------
  ! Purely vertical case
  !-----------------------------------------------------------------------------
  if (all(abs(o+p) < 1e-5)) then
    out = in
    if (verbose > 1) &
      write(io_unit%log,*) 'rt_integral_t%bilinear: no interpolation'
  !-----------------------------------------------------------------------------
  ! Interpolation only in y
  !-----------------------------------------------------------------------------
  else if (abs(o(1)+p(1)) < 1e-5) then
    if (verbose > 1) &
      write(stderr,*) 'rt_integral_t%bilinear: 2nd index, o, p =', &
        o(2), p(2)
    call check_index_2 ('(1)')
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
      write(stderr,*) 'rt_integral_t%bilinear: 1st index, o, p =', &
        o(1), p(1)
    call check_index_1 ('(1)')
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
      write(stderr,*) 'rt_integral_t%bilinear: o, p =', &
        o, p
    call check_index_1 ('(2)')
    call check_index_2 ('(2)')
    do ib=1,m(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,ib) = &
          q(1)*q(2)*in(i1+o(1)  ,i2+o(2)  ,ib) + &
          p(1)*q(2)*in(i1+o(1)+1,i2+o(2)  ,ib) + &
          q(1)*p(2)*in(i1+o(1)  ,i2+o(2)+1,ib) + &
          p(1)*p(2)*in(i1+o(1)+1,i2+o(2)+1,ib)
    end do
    end do
    end do
  end if
  call trc_end (itimer)
contains
!===============================================================================
subroutine check_index_1 (label)
  character(len=*):: label
  do i1=l(1),u(1),u(1)-l(1)
    if (i1+o(1) < 1 .or. i1+o(1)+1 > size(in,1)) then
      write(stderr,*) 'index 1 outside range '//label, i1, size(in,1), o(1), p(1)
      call io%abort ('rt_integral_t%bilinear: index outside range')
    end if
  end do
end subroutine check_index_1
!===============================================================================
subroutine check_index_2 (label)
  character(len=*):: label
  do i2=l(2),u(2),u(2)-l(2)
    if (i2+o(2) < 1 .or. i2+o(2)+1 > size(in,2)) then
      write(stderr,*) 'index 2 outside range '//label, i2, size(in,2), o(2), p(2)
      call io%abort ('rt_integral_t%bilinear: index outside range')
    end if
  end do
end subroutine check_index_2
END SUBROUTINE bilinear

!===============================================================================
!> Integral solver for one layer.  The exact integral solution Q=I-S for a 2nd
!> order polynomial has coefficient for the upstream Q, and for the 1st and 2nd
!> derivative of S.  However, if one is doing a solution both forwards and 
!> backward, the term in the 1st derivative cancels.
!>
!> The order of indices is in general permuted in the call, with 
!===============================================================================
SUBROUTINE solve1 (self, k, s, qin, qout, dr, l, u, debug)
  class(rt_integral_t):: self
  real, dimension(:,:,:,:) :: k, s, qin
  real, dimension(:,:,:) :: qout
  integer:: l(3), u(3)
  real:: dr
  logical:: debug
  !.............................................................................
  integer:: i1, i2, ib, m(4), icase
  real:: dsdtau1, dsdtau2, d2sdtau2, ex0, ex2, dtaumin, dtaumax
  real(8), pointer:: z(:)
  real, dimension(size(k,1)):: dtau1, dtau2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trc_begin ('rt_integral_t%solve1', itimer=itimer)
  m = shape(k)
  do ib=1,m(3)
    do i2=l(2),u(2)
      dtau1(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,1))
      dtau2(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,3))
      dtaumin = minval(dtau1(l(1):u(1)))
      dtaumax = maxval(dtau1(l(1):u(1)))
      !-------------------------------------------------------------------------
      ! For optical depth increments larger than dtau_max for all points: Q=0.0
      !-------------------------------------------------------------------------
      if (dtaumin > dtau_max) then
        do i1=l(1),u(1)
          qout(i1,i2,ib) = 0.0
        end do
        icase = 4
      !-------------------------------------------------------------------------
      ! For optical depth increments larger than dtau_dif: Q = diffusion approx
      !-------------------------------------------------------------------------
      else if (dtaumin > dtau_dif) then
        do i1=l(1),u(1)
          dsdtau1  =     (s(i1,i2,ib,2) - s(i1,i2,ib,1))/dtau1(i1)
          dsdtau2  =     (s(i1,i2,ib,3) - s(i1,i2,ib,2))/dtau2(i1)
          d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
          qout(i1,i2,ib) = d2sdtau2
        end do
        icase = 3
      !-------------------------------------------------------------------------
      ! For optical depth increments all maller than dtau_min: linear absorption
      !-------------------------------------------------------------------------
      else if (dtaumax < dtau_min) then
        do i1=l(1),u(1)
          ex0      = 1.0-dtau1(i1)
          qout(i1,i2,ib) = qin(i1,i2,ib,1)*ex0
        end do
        icase = 1
      !-------------------------------------------------------------------------
      ! For all other cases, use the full expression (assume sum of +/- directions)
      !-------------------------------------------------------------------------
      else 
        do i1=l(1),u(1)
          dsdtau1  =     (s(i1,i2,ib,2) - s(i1,i2,ib,1))/dtau1(i1)
          dsdtau2  =     (s(i1,i2,ib,3) - s(i1,i2,ib,2))/dtau2(i1)
          d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
          ex0      = exp(-dtau1(i1))
          ex2      = (1.0-ex0) - dtau1(i1)*ex0
          qout(i1,i2,ib) = qin(i1,i2,ib,1)*ex0 &
                         + d2sdtau2       *ex2
        end do
        icase = 2
      end if
      !$omp atomic
      n_cases(icase) = n_cases(icase)+1
    end do
    if (verbose > 1 .and. debug) then
      write(io_unit%log,1) ' rt_integral_t%solve1: ib, dt =', &
        ib, icase,  dtau1(l(1):l(1)+9)
      write(io_unit%log,1) ' rt_integral_t%solve1: ib, qi =', &
        ib, icase,  qin(l(1):l(1)+9,(l(2)+u(2))/2,ib,1)
      write(io_unit%log,1) ' rt_integral_t%solve1: ib, qo =', &
        ib, icase, qout(l(1):l(1)+9,(l(2)+u(2))/2,ib)
      1 format(a,i4,i2,1p,10e10.2)
    end if
  end do
  call trc_end (itimer)
END SUBROUTINE solve1

SUBROUTINE solve2 (self, k, s, qin, qout, dr, l, u, debug)
  class(rt_integral_t):: self
  real, dimension(:,:,:,:) :: k, s, qin
  real, dimension(:,:,:) :: qout
  integer:: l(3), u(3)
  real:: dr
  logical:: debug
  !.............................................................................
  integer:: i1, i2, ib, m(4), icase
  real:: dsdtau1, dsdtau2, d2sdtau1, d2sdtau2, ex0, ex1, ex2, dtaumin, dtaumax
  real(8), pointer:: z(:)
  real, dimension(size(k,1)):: dtau1, dtau2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trc_begin ('rt_integral_t%solve2', itimer=itimer)
  m = shape(k)
  do ib=1,m(3)
    do i2=l(2),u(2)
      dtau1(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,1))
      dtau2(l(1):u(1)) = dr*(k(l(1):u(1),i2,ib,2)+k(l(1):u(1),i2,ib,3))
        do i1=l(1),u(1)
          ex0 = exp(-dtau1(i1))
          ex1 = 1.-ex0
          ex2 = ex1 - dtau1(i1)*ex0
          dsdtau1  =     (s(i1,i2,ib,2) - s(i1,i2,ib,1))/dtau1(i1)
          dsdtau2  =     (s(i1,i2,ib,3) - s(i1,i2,ib,2))/dtau2(i1)
          d2sdtau1 = (dsdtau2 * dtau1(i1) + dsdtau1 * dtau2(i1))/(dtau1(i1)+dtau2(i1))
          d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1(i1)+dtau2(i1))
          
          qout(i1,i2,ib) = qin(i1,i2,ib,1)*ex0 &
                         - d2sdtau1       *ex1 &
                   + merge(d2sdtau2       *ex2, 0.0, dtau1(i1) > dtau_min) 
        end do
    end do
    if (verbose > 1 .and. debug) then
      write(io_unit%log,1) ' rt_integral_t%solve2: ib, dt =', &
        ib, icase,  dtau1(l(1):l(1)+9)
      write(io_unit%log,1) ' rt_integral_t%solve2: ib, qi =', &
        ib, icase,  qin(l(1):l(1)+9,(l(2)+u(2))/2,ib,1)
      write(io_unit%log,1) ' rt_integral_t%solve2: ib, qo =', &
        ib, icase, qout(l(1):l(1)+9,(l(2)+u(2))/2,ib)
      1 format(a,i4,i2,1p,10e10.2)
    end if
  end do
  call trc_end (itimer)
END SUBROUTINE solve2

!===============================================================================
SUBROUTINE diagnostics (label)
  character(len=*), optional:: label
  real:: f
  f = 1./max(sum(n_cases),1_8)
  if (present(label)) then
    write(io_unit%log,*) 'fractions of thin, full, diff, and zero cases: ', &
      n_cases*f, trim(label)
  else
    write(io_unit%log,*) 'fractions of thin, full, diff, and zero cases: ', &
    n_cases*f
  end if
END SUBROUTINE diagnostics 

!===============================================================================
SUBROUTINE trc_begin (label, itimer)
  character(len=*):: label
  integer, optional:: itimer
  if (detailed_timer) then
    call timer%begin (label, itimer)
  end if
END SUBROUTINE

!===============================================================================
SUBROUTINE trc_end (itimer)
  integer, optional:: itimer
  if (detailed_timer) then
    call timer%end (itimer)
  end if
END SUBROUTINE

END MODULE rt_integral_mod
