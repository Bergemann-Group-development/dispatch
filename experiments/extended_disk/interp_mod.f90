!===============================================================================
!> Handles intepolations at the face between two tasks.  This module
!> services guards_mod, and is replaceable via the makefile macro INTERP.
!===============================================================================
MODULE interp_mod
  USE trace_mod
  USE io_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE mesh_mod
  implicit none
  private
  type, public:: meshp_t
    class(mesh_t), pointer:: p => null()
  end type
  type, public:: interp_t
    class(task_t), pointer:: source => null()    ! source patch
    class(task_t), pointer:: target => null()    ! target patch
    real, allocatable:: r(:,:,:,:)
    type(lock_t):: lock
  contains
    procedure:: init
    procedure:: from_source
    procedure:: indices
    procedure, nopass:: in_t
    procedure, nopass:: in_x
    procedure, nopass:: in_y
    procedure, nopass:: in_z
  end type
  integer:: verbose = 0
  logical:: detailed_timer = .false.
  type(interp_t), public:: interp
CONTAINS

!===============================================================================
!> Initialize interpolation module
!===============================================================================
SUBROUTINE init (self, nbtask)
  class(interp_t):: self
  class(task_t), pointer, optional:: nbtask
  !.............................................................................
  integer:: i1, i2, i3, m(3), iostat
  logical, save:: first_time = .true.
  namelist /interp_params/ verbose, detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('interp_t%init')
  call io%lock%set ('interp_t%init')
  if (first_time) then
    !---------------------------------------------------------------------------
    ! Read namelist parameters
    !---------------------------------------------------------------------------
    rewind (io%input)
    read (io%input, interp_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, interp_params)
    end if
    call io%namelist_warning ('interp_params', iostat)
    write (io%namelist, interp_params)
    first_time = .false.
  end if
  call self%lock%init
  call io%lock%unset ('interp_t%init')
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Interpolate in 3-D index space
!===============================================================================
SUBROUTINE from_source (self, tgt, lt, ut, src)
  class(interp_t)           :: self
  class(patch_t)            :: tgt, src
  integer, dimension(3)     :: lt, ut
  !.............................................................................
  real, allocatable         :: p(:,:,:,:)
  real, pointer, contiguous :: tm(:,:,:,:), sm(:,:,:,:)
  real, allocatable         :: sp(:,:,:,:), tp(:,:,:,:)
  real                      :: pt(2)
  real(8)                   :: p1, p2, q1, q2, v_shift(3)
  logical                   :: outside(3)
  integer                   :: i1, i2, i3, j1 ,j2, iv, j, it(2)
  integer, dimension(3)     :: l, u, m(4)
  integer, save             :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%from_source', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Interpolate in time, in per unit mass variables
  !-----------------------------------------------------------------------------
  tm => tgt%mem(:,:,:,:,tgt%it,1)
  sm => src%mem(:,:,:,:,src%it,1)
  m = shape(sm)
  allocate (tp(m(1),m(2),m(3),m(4)))
  allocate (sp(m(1),m(2),m(3),m(4)))
  call src%time_interval (tgt%time, it, pt)
  call to_per_unit_mass (src, sp, it, pt)
  !-----------------------------------------------------------------------------
  ! Interpolation coordinates
  !-----------------------------------------------------------------------------
  allocate (p(src%gn(1),src%gn(3),src%gn(3),3))
  call interp%indices (tgt, lt, ut, src, p)
  !-----------------------------------------------------------------------------
  ! Bi-linear interpolation in per-unit-mass variables
  !-----------------------------------------------------------------------------
  v_shift = tgt%velocity_shift (src)
  l = src%mesh%li
  u = src%mesh%ui
  do iv=1,tgt%nv
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          !---------------------------------------------------------------------
          ! Skip interior points
          !---------------------------------------------------------------------
          j1 = p(i1,i2,i3,1)
          j2 = p(i1,i2,i3,2)
          j1 = min(max(j1,l(1)),u(1)-1)
          j2 = min(max(j2,l(2)),u(2)-1)
          p1 = p(i1,i2,i3,1) - j1
          p2 = p(i1,i2,i3,2) - j2
          q1 = 1.0 - p1
          q2 = 1.0 - p2
          !...................................................................
          tp(i1,i2,i3,iv) = &
            q2*(q1*sp(j1,j2  ,i3  ,iv) + p1*sp(j1+1,j2  ,i3  ,iv)) + &
            p2*(q1*sp(j1,j2+1,i3  ,iv) + p1*sp(j1+1,j2+1,i3  ,iv))
        end do
      end do
    end do
  end do
  l = src%mesh%li
  u = src%mesh%ui
  do i3=lt(3),ut(3)
    outside(3) = i3 < tgt%mesh(3)%li .or. i3 > tgt%mesh(3)%ui
    do i2=lt(2),ut(2)
      outside(2) = i2 < tgt%mesh(2)%li .or. i2 > tgt%mesh(2)%ui
      do i1=lt(1),ut(1)
        outside(1) = i1 < tgt%mesh(1)%li .or. i1 > tgt%mesh(1)%ui
        !-----------------------------------------------------------------------
        ! Skip interior points
        !-----------------------------------------------------------------------
        if (any(outside)) then
          tm(i1,i2,i3,tgt%idx%px) = (tp(i1,i2,i3,tgt%idx%d))
          !tm(i1,i2,i3,tgt%idx%px) = exp(tp(i1,i2,i3,tgt%idx%d))
          !---------------------------------------------------------------------
          ! Add velocity shift
          !---------------------------------------------------------------------
          tp(i1,i2,i3,tgt%idx%px) = tp(i1,i2,i3,tgt%idx%px) + v_shift(1)
          tp(i1,i2,i3,tgt%idx%py) = tp(i1,i2,i3,tgt%idx%py) + v_shift(2)
          tp(i1,i2,i3,tgt%idx%pz) = tp(i1,i2,i3,tgt%idx%pz) + v_shift(3)
          !---------------------------------------------------------------------
          ! Linear momenta
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%px) = &
          tp(i1,i2,i3,tgt%idx%px) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%py) = &
          tp(i1,i2,i3,tgt%idx%py) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%pz) = &
          tp(i1,i2,i3,tgt%idx%pz) * tm(i1,i2,i3,tgt%idx%d)
          !---------------------------------------------------------------------
          ! Total energy
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%e ) = &
          tp(i1,i2,i3,tgt%idx%e ) + ( &
          tp(i1,i2,i3,tgt%idx%px)**2 + &
          tp(i1,i2,i3,tgt%idx%py)**2 + &
          tp(i1,i2,i3,tgt%idx%pz)**2) * 0.5*tm(i1,i2,i3,tgt%idx%d)
        end if
      end do
    end do
  end do
  deallocate (tp, sp)
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE from_source

SUBROUTINE to_per_unit_mass (src, sp, it, pt)
  class(patch_t):: src
  real:: sp(:,:,:,:)
  integer:: it(2)
  real:: pt(2)
  !.............................................................................
  real, dimension(:,:,:), contiguous, pointer:: d1, e1, px1, py1, pz1
  real, dimension(:,:,:), contiguous, pointer:: d2, e2, px2, py2, pz2
  real:: di1, di2
  integer:: i1, i2, i3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%to_per_unit_mass', itimer=itimer)
  d1  => src%mem(:,:,:,src%idx%d ,it(1),1)
  e1  => src%mem(:,:,:,src%idx%e ,it(1),1)
  px1 => src%mem(:,:,:,src%idx%px,it(1),1)
  py1 => src%mem(:,:,:,src%idx%py,it(1),1)
  pz1 => src%mem(:,:,:,src%idx%pz,it(1),1)
  !.............................................................................
  d2  => src%mem(:,:,:,src%idx%d ,it(2),1)
  e2  => src%mem(:,:,:,src%idx%e ,it(2),1)
  px2 => src%mem(:,:,:,src%idx%px,it(2),1)
  py2 => src%mem(:,:,:,src%idx%py,it(2),1)
  pz2 => src%mem(:,:,:,src%idx%pz,it(2),1)
  do i3=1,src%gn(3)
  do i2=1,src%gn(2)
    !$omp simd
    do i1=1,src%gn(1)
      di1 = 1.0/d1(i1,i2,i3)
      di2 = 1.0/d2(i1,i2,i3)
      sp(i1,i2,i3,src%idx%d ) = pt(1)*d1(i1,i2,i3) &
                              + pt(2)*d2(i1,i2,i3)
      sp(i1,i2,i3,src%idx%px) = pt(1)*px1(i1,i2,i3)*di1 &
                              + pt(2)*px2(i1,i2,i3)*di2
      sp(i1,i2,i3,src%idx%py) = pt(1)*py1(i1,i2,i3)*di1 &
                              + pt(2)*py2(i1,i2,i3)*di2
      sp(i1,i2,i3,src%idx%pz) = pt(1)*pz1(i1,i2,i3)*di1 &
                              + pt(2)*pz2(i1,i2,i3)*di2
      sp(i1,i2,i3,src%idx%e ) = pt(1)*(e1(i1,i2,i3) - ( &
                                      px1(i1,i2,i3)**2 + &
                                      py1(i1,i2,i3)**2 + &
                                      pz1(i1,i2,i3)**2)*0.5*di1)*di1 &
                              + pt(2)*(e2(i1,i2,i3) - ( &
                                      px2(i1,i2,i3)**2 + &
                                      py2(i1,i2,i3)**2 + &
                                      pz2(i1,i2,i3)**2)*0.5*di2)*di2
    end do
  end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE to_per_unit_mass

!===============================================================================
!> Vectorized transformation from target index space to source index space.
!> The r(m,m,m,1:3) array holds the floating index distances from the patch
!> position, in units of cells.  The patch%e(j,k) array holds the absolute
!> Cartesian displacement in the j-direction for one cell along the local
!> unit vector in the k-direction.
!===============================================================================
SUBROUTINE indices (self, tgt, lt, ut, src, ps)
  class(interp_t):: self
  class(patch_t):: tgt, src
  integer, dimension(3):: lt, ut, ls, us
  real, contiguous:: ps(:,:,:,:)
  !.............................................................................
  real, allocatable:: pt(:,:)
  real(8):: t2s(3,3), tc(3,3), sc(3,3)
  integer:: i1, i2, i3, i, j, k
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%indices', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Cache the index space cell coordinates -- common to all patches
  !-----------------------------------------------------------------------------
  call self%lock%set ('indices')
  if (allocated(self%r)) then
    if (any(shape(self%r) /= [tgt%gn(1),tgt%gn(2),tgt%gn(3),3])) then
      deallocate (self%r)
    end if
  end if
  if (.not. allocated(self%r)) then
    allocate (self%r(tgt%gn(1),tgt%gn(2),tgt%gn(3),3))
    do i3=1,tgt%gn(3)
    do i2=1,tgt%gn(2)
    do i1=1,tgt%gn(1)
      self%r(i1,i2,i3,1) = (i1-tgt%mesh(1)%o)
      self%r(i1,i2,i3,2) = (i2-tgt%mesh(2)%o)
      self%r(i1,i2,i3,3) = (i3-tgt%mesh(3)%o)
    end do
    end do
    end do
  end if
  call self%lock%unset ('indices')
  !-----------------------------------------------------------------------------
  ! Cache target and source coefficient
  !-----------------------------------------------------------------------------
  do j=1,3
  do k=1,3
    tc(j,k) = tgt%e(j,k)*tgt%ds(k)
    sc(j,k) = src%e(j,k)/src%ds(k)
  end do
  end do
  allocate (pt(size(ps,1),3))
  ls = src%mesh%li
  us = src%mesh%ui
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      !-------------------------------------------------------------------------
      ! pt(:,j) is the displacement in the common Cartesian coordinate
      ! j-direction for each index triplet (i1,i2,i3) in the target index space,
      ! where the three indices give the location in the local unit vector
      ! coordinates. self%r(i1,i2,i3) is in units of cells, and tgt%e(j,k) is
      ! in code units, and hence p(:,j) is in code units.
      !-------------------------------------------------------------------------
      do j=1,3
        do i1=lt(1),ut(1)
          pt(i1,j) = tgt%position(j) &
                   - src%position(j) &
                   + self%r(i1,i2,i3,1)*tc(j,1) &
                   + self%r(i1,i2,i3,2)*tc(j,2) &
                   + self%r(i1,i2,i3,3)*tc(j,3)
        end do
      end do
      !-------------------------------------------------------------------------
      ! ps(:,:,:,j) are the final index space coordinates in the local source
      ! coordinate system, with a center at the src%position, and with the
      ! t2s(i,j) matrix transforming from the local target to the local source
      ! coordinates.
      !-------------------------------------------------------------------------
!define DEBUG
#ifdef DEBUG
      do i1=lt(1),ut(1)
        do j=1,3
          ps(i1,i2,i3,j) = src%mesh(j)%o &
                         + pt(i1,1)*sc(j,1) &
                         + pt(i1,2)*sc(j,2) &
                         + pt(i1,3)*sc(j,3)
        end do
        write(stdlog,'(a,2(i6,1x,4i3),1x,4i3,3f8.3)') &
          'interp_t%indices: tgt,src,ii,p=', &
          tgt%id, tgt%key, src%id, src%key, i1,i2,i3,j, ps(i1,i2,i3,:)
        block
          real(8),dimension(3):: pos, p, p2
          integer:: i
          pos = tgt%position - src%position + &
                (i1-tgt%mesh(1)%o)*tgt%ds(1)*tgt%e(:,1) + &
                (i2-tgt%mesh(2)%o)*tgt%ds(2)*tgt%e(:,2) + &
                (i3-tgt%mesh(3)%o)*tgt%ds(3)*tgt%e(:,3)
          p(1) = sum(pos*src%e(:,1))/src%ds(1) + src%mesh(1)%o
          p(2) = sum(pos*src%e(:,2))/src%ds(2) + src%mesh(2)%o
          p(3) = sum(pos*src%e(:,3))/src%ds(3) + src%mesh(3)%o
          if (any(abs(p-ps(i1,i2,i3,:)) > 0.001)) then
            write(stdlog,1) pt(i1,:), 'pnt src position in tgt coords'
            write(stdlog,1) src%position, 'src position'
            write(stdlog,1) tgt%position, 'tgt position'
            write(stdlog,1) pos, 'pnt distance from src%pos'
            write(stdlog,1) p, 'pnt src indices'
            write(stdlog,1) ps(i1,i2,i3,:), ' ps src indices'
          1 format(3f10.6,2x,a)
            flush (stdlog)
          end if
        end block
      end do
#else
      do j=1,3
        do i1=lt(1),ut(1)
          ps(i1,i2,i3,j) = src%mesh(j)%o &
                         + pt(i1,1)*sc(j,1) &
                         + pt(i1,2)*sc(j,2) &
                         + pt(i1,3)*sc(j,3)
        end do
      end do
#endif
    end do
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE indices

!===============================================================================
!> Interpolate in time over the relevant source region
!===============================================================================
SUBROUTINE in_t (target, source, ls, us, in, out, lt, ut, all, only, extrap)
  class(patch_t)                  :: source, target
  integer, dimension(3)           :: ls, us
  real                            :: in(:,:,:,:,:,:)
  real                            :: out(:,:,:,:)
!  intent(in)                      :: in
!  contiguous                      :: in, out
  integer, dimension(3), optional :: lt, ut
  logical, optional               :: all, extrap
  integer, optional               :: only
  !.............................................................................
  integer, dimension(3)           :: l, u, o
  integer                         :: i1, i2, i3, i4, jt(2), nv, n, li, ui
  real                            :: pt(2), v_min, d_min
  real, allocatable               :: p(:)
  real(8)                         :: wc, e(3,3), c12
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: lock1, lock2
  integer, save                   :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_t', itimer=itimer)
  if (verbose > 1) then
    wc = wallclock()
    if (ntimer == 0) then
      wtimer = 0d0
    end if
    if (verbose > 2) then
      call check_positive (source, &
        in(:,:,:,1,source%it,1), ls, us, 'src')
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If the the unit vector matrix exists the patches are tiltet slightly
  !-----------------------------------------------------------------------------
  if (.false.) then
    allocate(p(target%n(1)))
    ui = source%mesh(1)%ui
    write (stdout,'(a,3(/,3f12.6))') 'in_t: source%e =', source%e
    write (stdout,'(a,3(/,3f12.6))') 'in_t: target%e =', target%e
    e = matmul (target%e, source%einv)
    write (stdout,'(a,3(/,3f12.6))') 'in_t:        e =', e
    e = matmul (source%einv, target%e)
    write (stdout,'(a,3(/,3f12.6))') 'in_t:        e =', e
    !--- radial coordinates in
    p = (target%mesh(1)%r - target%mesh(1)%r(ui))
    write (stdout,'(a,/,16f8.4)')    'in_t:        r =', p
    p = p*e(1,2)/target%ds(2)
    write (stdout,'(a,/,16f8.4)')    'in_t:        p =', p
    !--- target distance along 1-vector
    p = (target%mesh(1)%r - target%mesh(1)%r(ui))
    !--- target 1-vector projection onto source 2-vector
    print *, e(1,2), sum(target%e(:,1)*source%e(2,:))
    !--- displcement along source 2-vector, in units of source cells
    c12 = sum(target%e(:,1)*source%e(2,:))/source%ds(2)
    p = p*c12
    write (stdout,'(a,/,16f8.4)')    'in_t:        p =', p
  end if
  !-----------------------------------------------------------------------------
  ! If the src is at the same time and level, just copy
  !-----------------------------------------------------------------------------
  if (source%time==target%time .and. source%level==target%level) then
    if (present(lt))   then
      l = lt
      u = ut
    else
      l = ls
      u = us
    end if
    o = 0
    !...........................................................................
    if (verbose > 1) then
      write (stdlog,'(a,2i6,4(2x,3i4))') &
      'interp_t%in_t: copy ls,us=', target%id, source%id, ls, us, l, u
      flush (stdlog)
    end if
    !...........................................................................
    out(l (1):u (1),l (2):u (2),l (3):u (3),:) = &
    in(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,1,1)
    n = product(us-ls+1)
    !...........................................................................
    if (verbose > 1) then
      wc = wallclock()-wc
      n = product(us-ls+1)
      write (stdlog,'(a,f6.3,1p,4(2x,3e11.2))') &
      'interp_t%in_t(1): mus/p,min,max,min,max=', wc/n, &
      minval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,1,1)), &
      maxval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,1,1)), &
      minval(out(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
      maxval(out(ls(1):us(1),ls(2):us(2),ls(3):us(3),1))
    end if
  !-----------------------------------------------------------------------------
  ! Interpolate in time
  !-----------------------------------------------------------------------------
  else
    call source%time_interval (target%time, jt, pt, id=target%id,extrap=extrap)
    !...........................................................................
    if (verbose > 1) then
      write (stdlog,'(a,6i3,3x,4i3)') &
      'interp_t%in_t: shapes=', shape(in), shape(out)
    end if
    !---------------------------------------------------------------------------
    ! Lock access to source while interpolating
    !---------------------------------------------------------------------------
    if (associated(source%lock)) then
      lock1 => source%mem_lock(jt(1))
      lock2 => source%mem_lock(jt(2))
      !call lock1%set ('interp_t%in_t')
      !call lock2%set ('interp_t%in_t')
    end if
    !---------------------------------------------------------------------------
    ! Shifted from source index space
    !---------------------------------------------------------------------------
    if (present(lt)) then
      if (verbose > 0) then
        write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,8l1)') &
          wallclock(), pad_trace(1:indent)// &
          'interp_t%in_t: tgt,src,jt,pt,ls,us,lt,ut,usg=', &
          target%id, source%id, jt, pt, ls, us, lt, ut, &
          source%unsigned
      end if
      !.........................................................................
      if (verbose > 2) then
        write (stdlog,'(a,1p,4(2x,2e11.2))') &
          'interp_t%in_t: min/max=', &
          minval( in(:,:,:,1,jt(1),1)), &
          maxval( in(:,:,:,1,jt(1),1)), &
          minval( in(:,:,:,1,jt(2),1)), &
          maxval( in(:,:,:,1,jt(2),1)), &
          minval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(1),1)), &
          maxval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(1),1)), &
          minval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(2),1)), &
          maxval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(2),1))
        flush (stdlog)
      end if
      !.........................................................................
      o = lt-ls
      do i4=1,target%nv
        if (skip(i4,only)) &
          cycle
        do i3=ls(3),us(3)
          !---------------------------------------------------------------------
          ! Use log interpolation only when needed
          !---------------------------------------------------------------------
          if (target%unsigned(i4) .and. product(pt) < 0.0) then
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                out(i1+o(1),i2+o(2),i3+o(3),i4) = &
                exp(log(in(i1,i2,i3,i4,jt(1),1))*pt(1) &
                  + log(in(i1,i2,i3,i4,jt(2),1))*pt(2))
              end do
            end do
          !---------------------------------------------------------------------
          ! Signed variables, or safe positive, shifted in index
          !---------------------------------------------------------------------
          else
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                out(i1+o(1),i2+o(2),i3+o(3),i4) = in(i1,i2,i3,i4,jt(1),1)*pt(1) &
                                                + in(i1,i2,i3,i4,jt(2),1)*pt(2)
              end do
            end do
          end if
        end do
      end do
      n = product(us-ls+1)
    !---------------------------------------------------------------------------
    ! Unshifted interpolation in time
    !---------------------------------------------------------------------------
    else
      o = 0
      !.........................................................................
      if (verbose > 0) then
        write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,3i3),1x,8l1)') &
          wallclock(), &
          pad_trace(1:indent)//'interp_t%in_t: tgt,src,jt,pt,ls,us,usg=', &
          target%id, source%id, jt, pt, ls, us, &
          source%unsigned
      end if
      !.........................................................................
      do i4=1,target%nv
        if (skip(i4,only)) &
          cycle
        v_min = huge(1.)
        do i3=ls(3),us(3)
          if (source%unsigned(i4)) then
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                out(i1,i2,i3,i4) = &
                  exp(log(in(i1,i2,i3,i4,jt(1),1))*pt(1) &
                    + log(in(i1,i2,i3,i4,jt(2),1))*pt(2))
                v_min = min(v_min,out(i1,i2,i3,i4))
                !...............................................................
                !write (stdlog,'(a,1(2x,3i3),2x,2i3,2f7.3,2x,1p,3e11.2)') &
                !  'in_t(2):', i1, i2, i3, jt, pt, &
                !  in(i1,i2,i3,i4,jt(1),1), in(i1,i2,i3,i4,jt(2),1), &
                !  out(i1,i2,i3,i4)
                !...............................................................
              end do
            end do
          else
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                out(i1,i2,i3,i4) = in(i1,i2,i3,i4,jt(1),1)*pt(1) &
                                 + in(i1,i2,i3,i4,jt(2),1)*pt(2)
                v_min = min(v_min,out(i1,i2,i3,i4))
              end do
            end do
          end if
        end do
        if (i4==1) d_min = v_min
      end do
      n = product(us-ls+1)
    end if
    if (associated(source%lock)) then
      !call lock2%unset ('interp_t%in_t')
      !call lock1%unset ('interp_t%in_t')
    end if
    !...........................................................................
    if (verbose > 1) then
      wc = wallclock()-wc
      n = product(us-ls+1)
      lt = ls+o
      write (stdlog,'(a,f6.3,1p,6(2x,3e11.2))') &
      'interp_t%in_t(2): mus/p,min,max,min,max=', wc/n, &
      minval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(1),1)), &
      maxval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(1),1)), &
      minval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(2),1)), &
      maxval( in(ls(1):us(1),ls(2):us(2),ls(3):us(3),1,jt(2),1)), &
      minval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),1)), &
      maxval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),1))
      if (verbose > 2) then
        call check_positive (target, out(:,:,:,1), lt, ut, 'in_t: out')
      end if
    end if
  end if
  !.............................................................................
  if (verbose > 1) then
    ntimer = ntimer + n
    wtimer = wtimer + wc
    if (ntimer > 1e6) then
      write (stdlog,'(a,6p,f7.3)') 'in_t: average mus/pt=', wtimer/ntimer
      ntimer = 0
    end if
  end if
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_t

!===============================================================================
!> Interpolate in x over the relevant target region.  The cell fraction
!> difference in x should be small -- of 2nd order in the tilt.  The exact
!> displacement in x may be found by taking the target position and ask for
!> the index position in the source patch.
!==============================================5================================
SUBROUTINE in_x (tgt, lt, ut, src, ls, us, om, in, out, all, only)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  class(meshp_t)           :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), nf(2), ub
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: x, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  integer, optional        :: only
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_x', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! If the source level is lower than the target level (prolongation), then
  ! ignore the fact that the target%filled may be at the target level already.
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Choose lower and upper index limits based on source and target ranges;
  ! where they already agree the values have already been interpolated
  !-----------------------------------------------------------------------------
  si = us - ls + 1                                      ! size per direction
  st = ut - lt + 1                                      ! tgt size
  where (si==st)                                        ! where size matches tgt
    l = lt                                              ! assume tgt
    u = ut                                              ! assume tgt
  else where                                            ! .. else ..
    l = ls                                              ! assume src
    u = us                                              ! assume src
  end where
  l(1) = lt(1)                                          ! new target range
  u(1) = ut(1)                                          ! new target range
  !-----------------------------------------------------------------------------
  ! The om() array holds the current interpolation state, and is updated here
  !-----------------------------------------------------------------------------
  om(1)%p => tgt%mesh(1)                                ! interpolation mesh
  tm => tgt%mesh(1)                                     ! target mesh
  sm => src%mesh(1)                                     ! source mesh
  d = tgt%distance_to (src)
  sm_p = d(1)                                           ! x-distance to src
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(a,2i6,3f6.2,3(2x,2(1x,3i3)),l2)') &
      'interp_t%in_x: src,tgt,dist,ls,us,lt,ut,l,u,all=', &
      src%id, tgt%id, d/tgt%size, ls, us, lt, ut, l, u, all
    if (verbose > 2) then
      write (stdlog,'(a,2(i6,2x,3f9.4))') &
        'interp_t%in_x: tm,sm n,p,d,o=', &
        tm%n, tm%p, tm%d, tm%o, &
        sm%n, sm_p, sm%d, sm%o
    end if
    flush (stdlog)
  end if
  !.............................................................................
  do i4=1,tgt%nv                                        ! loop over variables
    if (skip(i4,only)) &
      cycle
    !---------------------------------------------------------------------------
    ! Using the full valid range
    !---------------------------------------------------------------------------
    nf = 0
    ub = sm%ub
    if (all) then
      do i1=l(1),u(1)
        x = tm%r(i1)                                    ! x-position
        p = (x-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        j = min(j,ub-1)
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        do i3=l(3),u(3)
          !$omp simd
          do i2=l(2),u(2)
            out(i1,i2,i3,i4) = in(j  ,i2,i3,i4)*q &
                             + in(j+1,i2,i3,i4)*p
          end do
        end do
        !.......................................................................
        if (verbose > 1 .and. i4 == 1) then
          if (any(out(i1,l(2):u(2),l(3):u(3),i4) <= 0.0)) then
            write (stdlog,'(a,i6,2(i4,f6.4))') &
              'negative density in interp_t%in_x: id,i1,j,p=', &
              src%id, i1, x, j, p
            do i3=l(3),u(3)
              do i2=l(2),u(2)
                write (stdlog,'(3i3,1p,e11.2,i3)') &
                  i1,i2,i3,out(i1,i2,i3,i4), tgt%filled(i1,i2,i3)
              end do
            end do
            call io%abort ('non-positive densities(1) in interp_mod%in_x')
          end if
        end if
        !.......................................................................
      end do
      if (verbose > 1 .and. i4 == 1) then
        d_min = minval(out(l(1):u(1),l(2):u(2),l(3):u(3),1))
      end if
      n = product(u-l+1)
    !---------------------------------------------------------------------------
    ! Skipping the interior
    !---------------------------------------------------------------------------
    else
      n = 0
      do i1=l(1),u(1)
        outside(1) = (i1 < om(1)%p%li .or. i1 > om(1)%p%ui)
        x = tm%r(i1)                                    ! x-position
        p = (x-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        if (verbose > 2) &
          write (stdlog,*) i1, j, p
        !.......................................................................
        v_min = huge(1.0)
        do i3=l(3),u(3)
          outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
          do i2=l(2),u(2)
            outside(2) = (i2 < om(2)%p%li .or. i2 > om(2)%p%ui)
            if (any(outside)) then
              out(i1,i2,i3,i4) = in(j  ,i2,i3,i4)*q &
                               + in(j+1,i2,i3,i4)*p
#ifdef DEBUG
              !.................................................................
              write (stdlog,'(a,3(2x,3i3),2x,3l1,i4,f7.3,2x,1p,3e11.2)') &
                'in_x:', i1, i2, i3, li, ui, outside, j, p, &
                in(j,i2,i3,i4), in(j+1,i2,  i3,i4), out(i1,i2,i3,i4)
              !.................................................................
#endif
              n = n+1
              v_min = min(v_min,out(i1,i2,i3,i4))
            end if
          end do
        end do
      end do
      if (i4 == 1) d_min = v_min
    end if
  end do
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()-wc
    write (stdlog,'(a,3pf7.3,6p,f6.2,1p,2x,e11.2,2i5)') &
      'in_x: ms,mus/pt,min,max,nf=', wc, wc/n, d_min, nf
    !if (verbose > 2) then
    !  call check_positive (tgt,  in(:,:,:,1), ls, us, 'in_x in')
    !  call check_positive (tgt, out(:,:,:,1),  l,  u, 'in_x out')
    !end if
  end if
  !.............................................................................
  ls = l
  us = u
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_x

!===============================================================================
!> Interpolate in x over the relevant target region
!===============================================================================
SUBROUTINE in_y (tgt, lt, ut, src, ls, us, om, in, out, all, only)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  class(meshp_t)           :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4)
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: y, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  integer, optional        :: only
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_y', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Choose lower and upper index limits based on source and target ranges;
  ! where they already agree the values have already been interpolated
  !-----------------------------------------------------------------------------
  si = us - ls + 1                                      ! size per direction
  st = ut - lt + 1                                      ! tgt size
  where (si==st)                                        ! where size matches tgt
    l = lt                                              ! assume tgt
    u = ut                                              ! assume tgt
  else where                                            ! .. else ..
    l = ls                                              ! assume src
    u = us                                              ! assume src
  end where
  l(2) = lt(2)                                          ! new target range
  u(2) = ut(2)                                          ! new target range
  !-----------------------------------------------------------------------------
  ! The om() array holds the current interpolation state, and is updated here
  !-----------------------------------------------------------------------------
  om(2)%p => tgt%mesh(2)                                ! interpolation mesh
  tm => tgt%mesh(2)                                     ! target mesh
  sm => src%mesh(2)                                     ! source mesh
  d = tgt%distance_to (src)
  sm_p = d(2)                                           ! y-distance to src
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(a,2i6,3f6.2,3(2x,2(1x,3i3)),l2)') &
      'interp_t%in_y: src,tgt,dist,ls,us,lt,ut,l,u,all=', &
      src%id, tgt%id, d/tgt%size, ls, us, lt, ut, l, u, all
    if (verbose > 2) then
      write (stdlog,'(a,2(i6,2x,3f9.4))') &
        'interp_t%in_y: tm,sm n,p,d,o=', &
        tm%n, tm%p, tm%d, tm%o, &
        sm%n, sm_p, sm%d, sm%o
    end if
    flush (stdlog)
  end if
  !.............................................................................
  do i4=1,tgt%nv                                        ! loop over variables
    if (skip(i4,only)) &
      cycle
!---------------------------------------------------------------------------
    ! Using the full valid range
    !---------------------------------------------------------------------------
    if (all) then
      do i2=l(2),u(2)
        y = tm%r(i2)                                    ! y-position
        p = (y-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        do i3=l(3),u(3)
          !$omp simd
          do i1=l(1),u(1)
              out(i1,i2,i3,i4) = in(i1,j  ,i3,i4)*q &
                               + in(i1,j+1,i3,i4)*p
          end do
        end do
        !.......................................................................
        if (verbose > 2 .and. i4 == 1) then
          if (any(out(l(1):u(1),i2,l(3):u(3),i4) <= 0.0)) then
            write (stdlog,'(a,i6,2(i4,f7.4))') &
              'non-positive density(1) in interp_t%in_y: id,i3,j,p=', &
              src%id, i2, y, j, p
            call io%abort ('non-positive densities(1) in interp_mod%in_y')
          end if
        end if
        !.......................................................................
      end do
      n = product(u-l+1)
      if (verbose > 1 .and. i4 == 1) then
        d_min = minval(out(l(1):u(1),l(2):u(2),l(3):u(3),1))
        if (verbose > 2) then
          call check_positive (tgt, out(:,:,:,1), l, u, 'interp_mod%in_y')
        end if
      end if
    !---------------------------------------------------------------------------
    ! Skipping the interior
    !---------------------------------------------------------------------------
    else
      n = 0
      do i2=l(2),u(2)
        outside(2) = (i2 < om(2)%p%li .or. i2 > om(2)%p%ui)
        y = tm%r(i2)                                    ! y-position
        p = (y-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        v_min = huge(1.0)
        if (verbose > 0) then
          do i3=l(3),u(3)
            outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
            do i1=l(1),u(1)
              outside(1) = (i1 < om(1)%p%li .or. i1 > om(1)%p%ui)
              if (any(outside)) then
                out(i1,i2,i3,i4) = in(i1,j  ,i3,i4)*q &
                                 + in(i1,j+1,i3,i4)*p
                n = n+1
                v_min = min(v_min,out(i1,i2,i3,i4))
              end if
            end do
          end do
        else
          !---------------------------------------------------------------------
          ! Optimized into vectorizablae x-loops
          !---------------------------------------------------------------------
          do i3=l(3),u(3)
            outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
            if (any(outside(2:3))) then
              !$omp simd
              do i1=l(1),u(1)
                out(i1,i2,i3,i4) = in(i1,j  ,i3,i4)*q &
                                 + in(i1,j+1,i3,i4)*p
              end do
            else
              !$omp simd
              do i1=l(1),om(1)%p%li-1
                out(i1,i2,i3,i4) = in(i1,j  ,i3,i4)*q &
                                 + in(i1,j+1,i3,i4)*p
              end do
              !$omp simd
              do i1=om(1)%p%ui+1,u(1)
                out(i1,i2,i3,i4) = in(i1,j  ,i3,i4)*q &
                                 + in(i1,j+1,i3,i4)*p
              end do
            end if
          end do
        end if
      end do
      if (i4 == 1) d_min = v_min
    end if
  end do
  !............................................................................
  if (verbose > 1) then
    wc = wallclock()-wc
    write (stdlog,'(a,3pf7.3,6p,f6.2,1p,2x,2e11.2)') &
      'in_y: ms,mus/pt,min,max=', wc, wc/n, d_min
    write (stdlog,*) 'd_min=', d_min
  end if
  !.............................................................................
  ls = l
  us = u
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_y

!===============================================================================
!> Interpolate in x over the relevant target region
!===============================================================================
SUBROUTINE in_z (tgt, lt, ut, src, ls, us, om, in, out, all, only)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  class(meshp_t)           :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: z, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  integer, optional        :: only
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_z', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Choose lower and upper index limits based on source and target ranges;
  ! where they already agree the values have already been interpolated
  !-----------------------------------------------------------------------------
  si = us - ls + 1                                      ! size per direction
  st = ut - lt + 1                                      ! tgt size
  where (si==st)                                        ! where size matches tgt
    l = lt                                              ! assume tgt
    u = ut                                              ! assume tgt
  else where                                            ! .. else ..
    l = ls                                              ! assume src
    u = us                                              ! assume src
  end where
  l(3) = lt(3)                                          ! new target range
  u(3) = ut(3)                                          ! new target range
  !-----------------------------------------------------------------------------
  ! The om() array holds the current interpolation state, and is updated here
  !-----------------------------------------------------------------------------
  tm => tgt%mesh(3)                                     ! target mesh
  sm => src%mesh(3)                                     ! source mesh
  d = tgt%distance_to (src)                             ! relative to target
  sm_p = d(3)                                           ! z-displacement
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(a,2i6,3f6.2,3(2x,2(1x,3i3)),l2)') &
      'interp_t%in_z: src,tgt,dist,ls,us,lt,ut,l,u,all=', &
      src%id, tgt%id, d/tgt%size, ls, us, lt, ut, l, u, all
    if (verbose > 2) then
      write (stdlog,'(a,2(i6,2x,3f9.4))') &
        'interp_t%in_z: tm,sm n,p,d,o=', &
        tm%n, tm%p, tm%d, tm%o, &
        sm%n, sm_p, sm%d, sm%o
    end if
    flush (stdlog)
  end if
  !.............................................................................
  do i4=1,tgt%nv                                        ! loop over variables
    if (skip(i4,only)) &
      cycle
    !---------------------------------------------------------------------------
    ! Using the full valid range, not excluding the interior
    !---------------------------------------------------------------------------
    if (all) then
      do i3=l(3),u(3)
        z = tm%r(i3)                                    ! z-position
        p = (z-sm_p)/sm%d + sm%o                        ! float src index
        j = floor(p)                                    ! integer src index
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        do i2=l(2),u(2)
          !$omp simd
          do i1=l(1),u(1)
              out(i1,i2,i3,i4) = in(i1,i2,j  ,i4)*q &
                               + in(i1,i2,j+1,i4)*p
          end do
        end do
        !.......................................................................
      end do
      n = product(u-l+1)
      if (verbose > 1 .and. i4 == 1) then
        d_min = minval(out(l(1):u(1),l(2):u(2),l(3):u(3),1))
        if (verbose > 3) then
          call check_positive (tgt, out(:,:,:,1), l, u, 'interp_mod%in_z')
        end if
      end if
    !---------------------------------------------------------------------------
    ! Skipping the interior
    !---------------------------------------------------------------------------
    else
      n = 0
      do i3=l(3),u(3)
        outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
        z = tm%r(i3)                                    ! z-position
        p = (z-sm_p)/sm%d + sm%o                        ! float src index
        j = floor(p)                                    ! integer src index
        p = p - j                                       ! 2nd fraction
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        v_min = huge(1.0)
        do i2=l(2),u(2)
          outside(2) = (i2 < om(2)%p%li .or. i2 > om(2)%p%ui)
          do i1=l(1),u(1)
            outside(1) = (i1 < om(1)%p%li .or. i1 > om(1)%p%ui)
            if (any(outside)) then
              out(i1,i2,i3,i4) = in(i1,i2,j  ,i4)*q &
                               + in(i1,i2,j+1,i4)*p
#ifdef DEBUG
              !.................................................................
              write (stdlog,'(a,3(2x,3i3),2x,3l1,i4,f7.3,2x,1p,3e11.2)') &
                'in_z:', i1, i2, i3, tm%li, tm%ui, outside, j, p, &
                in(i1,i2,j,i4), in(i1,i2,j+1,i4), out(i1,i2,i3,i4)
              !.................................................................
#endif
              n = n+1
              v_min = min(v_min,out(i1,i2,i3,i4))
            end if
          end do
        end do
        !.......................................................................
      end do ! i3
      if (i4 == 1) d_min = v_min
    end if ! .not. all  => no interior
  end do ! i4
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()-wc
    write (stdlog,'(a,3pf7.3,6p,f6.2,1p,2x,2e11.2)') &
      'in_z: ms,mus/pt,min=', wc, wc/n, d_min
  end if
  !-----------------------------------------------------------------------------
  ! Update the source index range
  !-----------------------------------------------------------------------------
  ls = l
  us = u
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_z

!===============================================================================
!> Check positivity of array
!===============================================================================
SUBROUTINE check_positive (task, v, l, u, label)
  class(task_t)          :: task
  integer, dimension(3)  :: l, u
  real                   :: v(:,:,:)
  !.............................................................................
  character(len=*)       :: label
  integer                :: i1, i2, i3
  integer, save          :: itimer=0, max_error=30
  !-----------------------------------------------------------------------------
  call trace%begin ('interp_t%check_positive', itimer=itimer)
  write (stdlog,'(a,2(2x,3i3),2x,1p,2e12.3)') &
    'check_positive: l,u,v=', l, u, v(l(1),l(2),l(3)), v(u(1),u(2),u(3))
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    if (v(i1,i2,i3) <= 0.0) then
      write (stdlog,'(a,i6,2x,3i4,1p,e11.2,2x,a)') &
        'interp_t%check_positive:', &
        task%id, i1, i2, i3, v(i1,i2,i3), label
      !$omp atomic update
      max_error = max_error-1
      call io%assert (max_error >= 0, 'interpt%check_positive failed')
    end if
  end do
  end do
  end do
  call trace%end (itimer)
END SUBROUTINE check_positive

!===============================================================================
!> Skip an index, or not
!===============================================================================
LOGICAL FUNCTION skip (i, only)
  integer:: i, only
  optional:: only
  if (present(only)) then
    skip = i /= only
  else
    skip = .false.
  end if
END FUNCTION skip

END MODULE interp_mod
