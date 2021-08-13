!===============================================================================
!> Handles intepolations at the face between two tasks.  This module is
!> services face_mod, and is replaceable via the makefile macro INTERP
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
  contains
    procedure:: init
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
  call io%lock%unset ('interp_t%init')
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

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
  integer                         :: i1, i2, i3, i4, jt(2), nv, n, it
  integer                         :: unit, iv0, iv1
  real                            :: pt(2), v_min, d_min
  real(8)                         :: wc
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: lock1, lock2
  integer, save                   :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_t', itimer=itimer)
  unit = merge (stdout, stdlog, verbose > 2 .or. target%id==io%id_debug)
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
    it = source%it
    !...........................................................................
    if (verbose > 1) then
      write (stdlog,'(a,2i6,4(2x,3i4))') &
      'interp_t%in_t: copy ls,us=', target%id, source%id, ls, us, l, u
      flush (stdlog)
    end if
    !...........................................................................
    out(l (1):u (1),l (2):u (2),l (3):u (3),:) = &
     in(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,it,1)
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
      if (unit == stdout) then
        write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,8l1)') &
          'interp_t%in_t(1): tgt,src,jt,pt,ls,us,lt,ut,usg=', &
          target%id, target%level, target%time, &
          source%id, source%level, source%time, &
          jt, pt, ls, us, lt, ut, &
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
      if (present(only)) then
        iv0 = only
        iv1 = only
      else
        iv0 = 1
        iv1 = target%nv
      end if
      do i4=iv0,iv1                                         ! loop over variables
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
      if (unit == stdout) then
        write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,3i3),1x,8l1)') &
          'interp_t%in_t(2): tgt,src,jt,pt,ls,us,lt,ut,usg=', &
          target%id, target%level, target%time, &
          source%id, source%level, source%time, &
          jt, pt, ls, us, &
          source%unsigned
      end if
      !.........................................................................
      if (present(only)) then
        iv0 = only
        iv1 = only
      else
        iv0 = 1
        iv1 = target%nv
      end if
      do i4=iv0,iv1                                         ! loop over variables
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
                ! !...............................................................
                ! write (stdlog,'(a,1(2x,3i3),2x,2i3,2f7.3,2x,1p,3e11.2)') &
                !   'in_t(2):', i1, i2, i3, jt, pt, &
                !   in(i1,i2,i3,i4,jt(1),1), in(i1,i2,i3,i4,jt(2),1), &
                !   out(i1,i2,i3,i4)
                ! !...............................................................
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
  !...........................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_t

!===============================================================================
!> Interpolate in x over t  he relevant target region
!==============================================5=================================
SUBROUTINE in_x (tgt, lt, ut, src, ls, us, om, in, out, all, only)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  class(meshp_t)           :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), nf(2), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: x, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer)  call trace%begin ('interp_t%in_x', itimer=itimer)
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
  d = tgt%distance_to (src)                             ! src%pos - tgt%pos
  sm_p = tm%p + d(1)
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
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = tgt%nv
  end if
  do i4=iv0,iv1                                         ! loop over variables
    !---------------------------------------------------------------------------
    ! Using the full valid range
    !---------------------------------------------------------------------------
    nf = 0
    if (all) then
      do i1=l(1),u(1)
        x = tm%p + tm%r(i1)                             ! x-position
        p = (x-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(1)),us(1)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        do i3=l(3),u(3)
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
        x = tm%p + tm%r(i1)                             ! x-position
        p = (x-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(1)),us(1)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
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
        !.......................................................................
        if (verbose > 2 .and. i4 == 1) then
          if (v_min <= 0.0) then
            write (stdlog,'(a,i6,2(i4,f6.4))') &
              'non-positive density in interp_t%in_x: id,i3,j,p=', &
              src%id, i1, x, j, p
            if (verbose > 3) then
              do i3=l(3),u(3)
                 outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
                do i2=l(2),u(2)
                  outside(2) = (i2 < om(2)%p%li .or. i2 > om(2)%p%ui)
                  if (any(outside) .and. out(i1,i2,i3,i4) <= 0.0) then
                    write (stdlog,'(a,3i4,1p,3e11.2)') &
                      'in_x non-positive(2):', i1, i2, i3, &
                      in (j  ,i2,i3,i4), &
                      in (j+1,i2,i3,i4), &
                      out(i1  ,i2,i3,i4)
                  end if
                end do
              end do
            end if
            call io%abort ('non-positive densities(2) in interp_mod%in_x')
          end if
        end if
        !.......................................................................
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
  if (detailed_timer)  call trace%end (itimer)
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
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: y, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer)  call trace%begin ('interp_t%in_y', itimer=itimer)
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
  d = tgt%distance_to (src)                             ! src%pos - tgt%pos
  sm_p = tm%p + d(2)
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
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = tgt%nv
  end if
  do i4=iv0,iv1                                         ! loop over variables
!---------------------------------------------------------------------------
    ! Using the full valid range
    !---------------------------------------------------------------------------
    if (all) then
      do i2=l(2),u(2)
        y = tm%p + tm%r(i2)                             ! x-position
        p = (y-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(2)),us(2)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
        q = 1.0 - p                                     ! 1st fraction
        !.......................................................................
        do i3=l(3),u(3)
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
        y = tm%p + tm%r(i2)                             ! x-position
        p = (y-sm_p)/sm%d + sm%o                        ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(2)),us(2)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
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
        !.......................................................................
        if (verbose > 2 .and. i4 == 1) then
          if (v_min <= 0.0) then
            write (stdlog,'(a,i6,2(i4,f7.4))') &
            'non-positive density(2) in interp_t%in_y: id,i3,j,p=', &
            src%id, i2, y, j, p
            if (verbose > 3) then
              do i3=l(3),u(3)
              outside(3) = (i3 < om(3)%p%li .or. i3 > om(3)%p%ui)
                do i1=l(1),u(1)
                outside(1) = (i1 < om(1)%p%li .or. i1 > om(1)%p%ui)
                  if (any(outside) .and. out(i1,i2,i3,i4) <= 0.0) then
                    write (stdlog,'(a,3i4,1p,3e11.2)') &
                      'in_y non-positive(2):', i1, i2, i3, &
                      in (i1,j  ,i3,i4), &
                      in (i1,j+1,i3,i4), &
                      out(i1,i2 ,i3 ,i4)
                  end if
                end do
              end do
            end if
            call io%abort ('non-positive densities(2) in interp_mod%in_y')
          end if
        end if
        !.......................................................................
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
  if (detailed_timer)  call trace%end (itimer)
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
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: z, wc, wtimer, sm_p, d(3)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer)  call trace%begin ('interp_t%in_z', itimer=itimer)
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
  d = tgt%distance_to (src)                             ! src%p - tgt%p
  sm_p = tm%p + d(3)
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
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = tgt%nv
  end if
  do i4=iv0,iv1                                         ! loop over variables
    !---------------------------------------------------------------------------
    ! Using the full valid range, not excluding the interior
    !---------------------------------------------------------------------------
    if (all) then
      do i3=l(3),u(3)
        z = tm%p + tm%r(i3)                             ! x-position
        p = (z-sm_p)/sm%d + sm%o                        ! float src index
        j = floor(p)                                    ! integer src index
        j = min(max(j,ls(3)),us(3)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
        q = 1.0 - p                                     ! 1st fraction
        do i2=l(2),u(2)
          do i1=l(1),u(1)
              out(i1,i2,i3,i4) = in(i1,i2,j  ,i4)*q &
                               + in(i1,i2,j+1,i4)*p
          end do
        end do
        !.......................................................................
        if (verbose > 2 .and. i4 == 1) then
          if (any(out(l(1):u(1),l(2):u(2),i3,i4) <= 0.0)) then
            write (stdlog,'(a,i6,2(i4,f7.3))') &
            'negative density in interp_t%in_z: id,i3,z,j,p=', &
            src%id, i3, z, j, p
            if (verbose > 3) then
              do i2=l(2),u(2)
                do i1=l(1),u(1)
                  if (out(i1,i2,i3,i4) <= 0.0) then
                    write (stdlog,'(a,3i4,1p,3e11.2)') &
                      'non-positive(1):', i1, i2, i3, &
                      in (i2,i2,j  ,i4), &
                      in (i2,i2,j+1,i4), &
                      out(i2,i2,i3 ,i4)
                  end if
                end do
              end do
            end if
            call io%abort ('negative densities(1)')
          end if
        end if
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
        z = tm%p + tm%r(i3)                             ! x-position
        p = (z-sm_p)/sm%d + sm%o                        ! float src index
        j = floor(p)                                    ! integer src index
        j = min(max(j,ls(3)),us(3)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
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
        if (verbose > 2 .and. i4 == 1) then
          if (v_min <= 0.0) then
            write (stdlog,'(a,i6,2(i4,f7.3))') &
            'negative density in interp_t%in_z: id,i3,z,j,p=', &
            src%id, i3, z, j, p
            if (verbose > 3) then
              do i2=l(2),u(2)
                outside(2) = (i2 < om(2)%p%li .or. i2 > om(2)%p%ui)
                do i1=l(1),u(1)
                  outside(1) = (i1 < om(1)%p%li .or. i1 > om(1)%p%ui)
                  if (any(outside) .and. out(i1,i2,i3,i4) <= 0.0) then
                    write (stdlog,'(a,3i4,1p,3e11.2)') &
                      'non-positive(2):', i1, i2, i3, &
                      in (i2,i2,j  ,i4), &
                      in (i2,i2,j+1,i4), &
                      out(i2,i2,i3 ,i4)
                  end if
                end do ! i1
              end do ! i2
            end if ! verbose > 3
            call io%abort ('negative densities(2)')
          end if ! v_min <= 0.0
        end if ! versose > 1 .and. i4 == 1
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
  if (detailed_timer)  call trace%end (itimer)
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
