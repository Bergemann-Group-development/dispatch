!===============================================================================
!> Handles intepolations at the face between two tasks.  This module
!> services guards_mod, and is replaceable via the makefile macro INTERP.
!===============================================================================
MODULE interp_mod
  USE iso_fortran_env, only: int8
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
    real, pointer, contiguous:: r(:,:,:,:)
    type(lock_t):: lock
  contains
    procedure:: init
    procedure:: from_source
    procedure:: indices
    procedure, nopass:: in_t
    procedure, nopass:: no_t
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
  !-----------------------------------------------------------------------------
  ! Read namelist inside double-if, since this is called repeatedly
  !-----------------------------------------------------------------------------
  if (first_time) then
    call io%lock%set ('interp_t%init')
    if (first_time) then
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
  end if
  call self%lock%init
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Interpolate in 3-D index space
!===============================================================================
SUBROUTINE from_source (self, tgt, lt, ut, src, filled, all_src)
  class(interp_t)             :: self
  class(patch_t)              :: tgt, src
  integer, dimension(3)       :: lt, ut
  integer(kind=int8), pointer :: filled(:,:,:)
  optional                    :: filled
  logical, optional           :: all_src
  !.............................................................................
  real, pointer, contiguous   :: p(:,:,:,:)
  real(8), pointer, contiguous   :: tr1(:), tr2(:), tr3(:)
  real, pointer, contiguous   :: tm(:,:,:,:), sm(:,:,:,:)
  real, pointer, contiguous   :: sp(:,:,:,:), tp(:,:,:,:)
  real                        :: pt(2), r1, r2, r
  real(8)                     :: p1, p2, q1, q2, p3, q3, v(3), v_shift(3), d_omega
  real(8)                     :: d(3), e(3,3)
  logical                     :: outside(3)
  integer                     :: i1, i2, i3, j1 ,j2, j3, iv, j, it(2), idx_es
  integer, dimension(3)       :: l, u, m(4)
  logical, pointer            :: available(:,:,:)
  logical                     :: all_src_loc
  integer                     :: m_extrap
  integer, save               :: itimer=0, n_extrap=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%from_source', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Interpolate in time, in per unit mass variables, leaving source values in
  ! sp(:,:,:,:)
  !-----------------------------------------------------------------------------
  tm => tgt%mem(:,:,:,:,tgt%it,1)
  sm => src%mem(:,:,:,:,src%it,1)
  m = shape(sm)
  if (present(all_src)) then
    all_src_loc = all_src
  else
    all_src_loc = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Combine handling of energy and entropy per unit mass
  !-----------------------------------------------------------------------------
  idx_es = max(src%idx%e,src%idx%s)  
  !-----------------------------------------------------------------------------
  ! the available grid holds logicals, telling us, which cells are available to
  ! be filled from source.
  !-----------------------------------------------------------------------------
  allocate (available(m(1),m(2),m(3)))
  available(:,:,:) = .false.
  !
  allocate (tp(m(1),m(2),m(3),m(4)))
  allocate (sp(m(1),m(2),m(3),m(4)))
  call src%time_interval (tgt, tgt%time, it, pt)
  call to_per_unit_mass (src, sp, it, pt)
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from tgt to src
  !-----------------------------------------------------------------------------
  call src%get_position (tgt, p=d, e=e)
  !-----------------------------------------------------------------------------
  ! Interpolation coordinates and weights p(:,:,:,1:3)
  !-----------------------------------------------------------------------------
  allocate (p(src%gn(1),src%gn(3),src%gn(3),3))
  call interp%indices (tgt, lt, ut, src, p, available, d, e)
  !-----------------------------------------------------------------------------
  ! Bi-linear interpolation in per-unit-mass variables
  !-----------------------------------------------------------------------------
  l = src%mesh%li
  u = src%mesh%ui
  do iv=1,tgt%nv
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        if (all_src_loc) then
          do i1=lt(1),ut(1)
            if (.not.available(i1,i2,i3)) &
              cycle
            !-------------------------------------------------------------------
            ! Skip interior points
            !-------------------------------------------------------------------
            j1 = p(i1,i2,i3,1)
            j2 = p(i1,i2,i3,2)
            j3 = p(i1,i2,i3,3)
            j1 = min(max(j1,l(1)),u(1)-1)
            j2 = min(max(j2,l(2)),u(2)-1)
            j3 = min(max(j3,l(3)),u(3)-1)
            p1 = p(i1,i2,i3,1) - j1
            p2 = p(i1,i2,i3,2) - j2
            p3 = p(i1,i2,i3,3) - j3
            if (any([p1,p2,p2] > 1.5 .or. [p1,p2,p3] < -0.5)) then
              ! PRINT '(2i6,4x,3f8.3,2x,a)', &
              ! tgt%id, src%id, j1+p1,j2+p2,j3+p3,'extrapolation prevented'
              !$omp atomic capture
              n_extrap = n_extrap+1
              m_extrap = n_extrap
              !$omp end atomic
              if (m_extrap == 100000) then
                call io%warning ("100,000 extrapolations prevented")
                !$omp atomic write
                n_extrap = 0
              end if
            end if
            p1 = min(max(p1,0.0),1.0)
            p2 = min(max(p2,0.0),1.0)
            p3 = min(max(p3,0.0),1.0)
            q1 = 1.0 - p1
            q2 = 1.0 - p2
            q3 = 1.0 - p3
            !...................................................................
            if (filled(i1,i2,i3)==0) then
              tp(i1,i2,i3,iv) = &
                q3*(q2*(q1*sp(j1,j2  ,j3  ,iv) + p1*sp(j1+1,j2  ,j3  ,iv)) + &
                    p2*(q1*sp(j1,j2+1,j3  ,iv) + p1*sp(j1+1,j2+1,j3  ,iv))) + &
                p3*(q2*(q1*sp(j1,j2  ,j3+1,iv) + p1*sp(j1+1,j2  ,j3+1,iv)) + &
                    p2*(q1*sp(j1,j2+1,j3+1,iv) + p1*sp(j1+1,j2+1,j3+1,iv)))
            end if
          end do
        else
          do i1=lt(1),ut(1)
            if (.not.available(i1,i2,i3)) &
              cycle
            !-------------------------------------------------------------------
            ! Skip interior points
            !-------------------------------------------------------------------
            j1 = p(i1,i2,i3,1)
            j2 = p(i1,i2,i3,2)
            j3 = p(i1,i2,i3,3)
            j1 = min(max(j1,l(1)),u(1)-1)
            j2 = min(max(j2,l(2)),u(2)-1)
            j3 = min(max(j3,l(3)),u(3)-1)
            p1 = p(i1,i2,i3,1) - j1
            p2 = p(i1,i2,i3,2) - j2
            p3 = p(i1,i2,i3,3) - j3
            if (any([p1,p2,p2] > 1.5 .or. [p1,p2,p3] < -0.5)) then
              ! PRINT '(2i6,4x,3f8.3,2x,a)', &
              ! tgt%id, src%id, j1+p1,j2+p2,j3+p3,'extrapolation prevented'
              !$omp atomic capture
              n_extrap = n_extrap+1
              m_extrap = n_extrap
              !$omp end atomic
              if (m_extrap == 100000) then
                call io%warning ("100,000 extrapolations prevented")
                !$omp atomic write
                n_extrap = 0
              end if
            end if
            p1 = min(max(p1,0.0),1.0)
            p2 = min(max(p2,0.0),1.0)
            p3 = min(max(p3,0.0),1.0)
            q1 = 1.0 - p1
            q2 = 1.0 - p2
            q3 = 1.0 - p3
            !...................................................................
            tp(i1,i2,i3,iv) = &
                q3*(q2*(q1*sp(j1,j2  ,j3  ,iv) + p1*sp(j1+1,j2  ,j3  ,iv)) + &
                    p2*(q1*sp(j1,j2+1,j3  ,iv) + p1*sp(j1+1,j2+1,j3  ,iv))) + &
                p3*(q2*(q1*sp(j1,j2  ,j3+1,iv) + p1*sp(j1+1,j2  ,j3+1,iv)) + &
                    p2*(q1*sp(j1,j2+1,j3+1,iv) + p1*sp(j1+1,j2+1,j3+1,iv)))
          end do
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from src to tgt
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, p=d, e=e)
  !-----------------------------------------------------------------------------
  ! Deposit interpolated values in target mem, taking into account that the
  ! tp values are "per unit mass"
  !-----------------------------------------------------------------------------
  l = src%mesh%li
  u = src%mesh%ui
  tr1 => tgt%mesh(1)%r
  tr2 => tgt%mesh(2)%r
  do i3=lt(3),ut(3)
    outside(3) = i3 < tgt%mesh(3)%li .or. i3 > tgt%mesh(3)%ui
    do i2=lt(2),ut(2)
      outside(2) = i2 < tgt%mesh(2)%li .or. i2 > tgt%mesh(2)%ui
      do i1=lt(1),ut(1)
        outside(1) = i1 < tgt%mesh(1)%li .or. i1 > tgt%mesh(1)%ui
        !-----------------------------------------------------------------------
        ! Skip interior points
        !-----------------------------------------------------------------------
        if (any(outside) .and. available(i1,i2,i3)) then
          !---------------------------------------------------------------------
          ! Mass density
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%d) = tp(i1,i2,i3,tgt%idx%d)
          !---------------------------------------------------------------------
          ! The coordinate velocity of each cell in a patch is omega*r, where
          ! r is the local radius, and the velocity shift to add from source
          ! to target is the difference from target to source.  The direction
          ! of that velocity is along phi_hat
          !---------------------------------------------------------------------
          r1 = tr1(i1)*tgt%ds(1) + tgt%centre_nat(1)
          r2 = tr2(i2)*tgt%ds(2)
          r = sqrt(r1**2+r2**2)
          d_omega = (src%omega(3)-tgt%omega(3))
          v_shift = d_omega*[-r2,+r1,0.0]
          !v_shift = d_omega*[r2,+r1,0.0]
          !---------------------------------------------------------------------
          ! Get the src velocity interpolated to the tgt points
          !---------------------------------------------------------------------
          v(1) = tp(i1,i2,i3,tgt%idx%px)
          v(2) = tp(i1,i2,i3,tgt%idx%py)
          v(3) = tp(i1,i2,i3,tgt%idx%pz)
          !---------------------------------------------------------------------
          ! Transform the velocity components to the tgt system.  Since e is
          ! essentially the rotation matrix of src (minus that of tgt), the
          ! matmul() operatin to go from src to tgt should be matmul(e,v):
          ! The full operation, via Eulerian, would be:
          !   1) get the local src velocity
          !   2) transform components to the Eulerian
          !   3) add the velocity difference btw patches in the E-frame
          !   4) transform to tgt coordinates
          ! But the velocity difference is already expressed in the tgt system,
          ! so this should be fine:
          !---------------------------------------------------------------------
          v = matmul (v, e)
          !---------------------------------------------------------------------
          ! Add velocity shift, which is already in the tgt system
          !---------------------------------------------------------------------
          v = v + v_shift
          tm(i1,i2,i3,tgt%idx%px) = v(1)
          tm(i1,i2,i3,tgt%idx%py) = v(2)
          tm(i1,i2,i3,tgt%idx%pz) = v(3)
          !---------------------------------------------------------------------
          ! Total energy in the target coordinate system.  The tp(...idx%e) is
          ! the thermal energy per unit mass, interpolated to the target point,
          ! to which needs to be added the kinetic energy in the target frame,
          ! keeping in mind that the tm(...%px) values are still velocities
          !---------------------------------------------------------------------
         ! print*, 'd,e', tm(i1,i2,i3,tgt%idx%d), tp(i1,i2,i3,tgt%idx%e )
          if (tgt%etype == 'total') then
            tm(i1,i2,i3,tgt%idx%e ) = tm(i1,i2,i3,tgt%idx%d)*( &
            tp(i1,i2,i3,tgt%idx%e ) + 0.5*( &
            tm(i1,i2,i3,tgt%idx%px)**2 + &
            tm(i1,i2,i3,tgt%idx%py)**2 + &
            tm(i1,i2,i3,tgt%idx%pz)**2))
          else
            tm(i1,i2,i3,    idx_es) = tm(i1,i2,i3,tgt%idx%d)* &
                                      tp(i1,i2,i3,   idx_es)
          end if
         ! print*, 'tot e', tm(i1,i2,i3,tgt%idx%e ), i1,i2,i3, p(i1,i2,i3,:)
          !---------------------------------------------------------------------
          ! Convert from velocities to linear momenta
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%px) = &
          tm(i1,i2,i3,tgt%idx%px) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%py) = &
          tm(i1,i2,i3,tgt%idx%py) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%pz) = &
          tm(i1,i2,i3,tgt%idx%pz) * tm(i1,i2,i3,tgt%idx%d)
          if (tgt%id == io%id_debug) then
            PRINT 1,src%id,i1,i2,i3,tm(i1,i2,i3,tgt%idx%d)
            1 format("from:source",i6,3i3,1p,g12.3)
          end if
          ! --------------------------------------------------------------------
          ! Fill in the filled grid
          ! --------------------------------------------------------------------
          if (associated(filled)) filled(i1,i2,i3)=max(src%level,tgt%level)
        end if
      end do
    end do
  end do
  deallocate (p, tp, sp, available)
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE from_source

!===============================================================================
!> Interpolate in time and convert to velocity and internal energy per mass
!===============================================================================
SUBROUTINE to_per_unit_mass (src, sp, it, pt)
  class(patch_t):: src
  real:: sp(:,:,:,:)
  integer:: it(2), l(3), u(3)
  real:: pt(2)
  !.............................................................................
  real, dimension(:,:,:), contiguous, pointer:: d1, e1, px1, py1, pz1
  real, dimension(:,:,:), contiguous, pointer:: d2, e2, px2, py2, pz2
  real:: di1, di2
  integer:: i1, i2, i3, idx_es
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%to_per_unit_mass', itimer=itimer)
  idx_es = max(src%idx%e,src%idx%s)  
  d1  => src%mem(:,:,:,src%idx%d ,it(1),1)
  d2  => src%mem(:,:,:,src%idx%d ,it(2),1)
  px1 => src%mem(:,:,:,src%idx%px,it(1),1)
  px2 => src%mem(:,:,:,src%idx%px,it(2),1)
  py1 => src%mem(:,:,:,src%idx%py,it(1),1)
  py2 => src%mem(:,:,:,src%idx%py,it(2),1)
  pz1 => src%mem(:,:,:,src%idx%pz,it(1),1)
  pz2 => src%mem(:,:,:,src%idx%pz,it(2),1)
  e1  => src%mem(:,:,:,    idx_es,it(1),1)
  e2  => src%mem(:,:,:,    idx_es,it(2),1)
  !.............................................................................
  l = src%mesh%lo
  u = src%mesh%uo
  if (trim(src%etype)=='total') then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
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
        sp(i1,i2,i3,    idx_es) = pt(1)*(e1(i1,i2,i3) - ( &
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
  else
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
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
        sp(i1,i2,i3,    idx_es) = pt(1)* e1(i1,i2,i3)*di1 &
                                + pt(2)* e2(i1,i2,i3)*di2
      end do
    end do
    end do
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE to_per_unit_mass

!===============================================================================
!> Vectorized transformation from target index space to source index space.
!> The r(m,m,m,1:3) array holds the floating index distances from the patch
!> position, in units of cells.  The patch%erot(j,k) array holds the absolute
!> Cartesian displacement in the j-direction for one cell along the local
!> unit vector in the k-direction.
!===============================================================================
SUBROUTINE indices (self, tgt, lt, ut, src, ps, available, d, e)
  class(interp_t) :: self
  class(patch_t)  :: tgt, src
  integer         :: lt(3), ut(3)
  real(8)         :: d(3), e(3,3), c(3), o(3), v(3)
  real, contiguous:: ps(:,:,:,:)
  logical, pointer:: available(:,:,:)
  !.............................................................................
  real, contiguous, pointer:: r(:,:,:,:)
  integer:: i1, i2, i3, i, j, k, ls(3), us(3)
  integer, save:: itimer=0
  real(8), parameter:: eps=0.5001d0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%indices', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Cache the index space cell coordinates -- common to all patches
  !-----------------------------------------------------------------------------
  allocate (r(tgt%gn(1),tgt%gn(2),tgt%gn(3),3))
  do i3=1,tgt%gn(3)
    do i2=1,tgt%gn(2)
      do i1=1,tgt%gn(1)
        r(i1,i2,i3,1) = tgt%mesh(1)%r(i1)
        r(i1,i2,i3,2) = tgt%mesh(2)%r(i2)
        r(i1,i2,i3,3) = tgt%mesh(3)%r(i3)
      end do
    end do
  end do
  c = 1d0/src%mesh%d
  o = src%mesh%o
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      !-------------------------------------------------------------------------
      ! ps(:,j) are the integer index space coordinates in the local source
      ! coordinate system.  We get this by rotating the target points into the
      ! source coordinate system, and rescaling them to integer index space.
      !-------------------------------------------------------------------------
      do i1=lt(1),ut(1)
        ! v = matmul(r(i1,i2,i3,:), tgt%erot)                ! tgt to global
        ! ps(i1,i2,i3,:) = (matmul(src%erot,v) + d)*c + o    ! global to src
        !-----------------------------------------------------------------------
        ! The e-matrix may be used to go from source to target coordinates, but
        ! here we need to go from target to source, so we should multiply from
        ! the left!
        !-----------------------------------------------------------------------
        ps(i1,i2,i3,:) = (matmul(r(i1,i2,i3,:),e)+d)*c + o
        !ps(i1,i2,i3,:) = (matmul(e,r(i1,i2,i3,:))+d)*c + o
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! available(:,:,:) is an array, indicating if the ps() is in the valid range
  ! inside the source li and ui
  !-----------------------------------------------------------------------------
  ls = src%mesh%li
  us = src%mesh%ui
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        available(i1,i2,i3) = (all(int(ps(i1,i2,i3,:)-eps)<us).and. &
                               all(int(ps(i1,i2,i3,:)+eps)>=ls))
      end do
    end do
  end do
  deallocate (r)
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
  real, pointer, contiguous       :: in1(:,:,:,:), in2(:,:,:,:)
  integer, dimension(3)           :: l, u, o
  integer                         :: i1, i2, i3, i4, jt(2), nv, n, it, id, m(4)
  integer                         :: unit, iv0, iv1
  real                            :: pt(2), v_min, d_min
  real, pointer, contiguous       :: d(:,:,:)
  real(8)                         :: wc
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: lock1, lock2, lock3
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
  call source%time_interval (target, target%time, jt, pt, extrap=extrap)
  !-----------------------------------------------------------------------------
  ! Lock access to source while interpolating
  !-----------------------------------------------------------------------------
  m = shape(out)
  allocate (in1(m(1),m(2),m(3),m(4)))
  allocate (in2(m(1),m(2),m(3),m(4)))
  !
  lock1 => source%mem_lock(jt(1))
  lock2 => source%mem_lock(jt(2))
  lock3 => target%mem_lock(target%it)
                                       call lock1%set   ('interp_t%in_t src(1)')
  do i4=1,source%nv
  do i3=ls(3),us(3)
  do i2=ls(2),us(2)
  do i1=ls(1),us(1)
    in1(i1,i2,i3,i4) = in(i1,i2,i3,i4,jt(1),1)
  end do
  end do
  end do
  end do
                                       call lock1%unset ('interp_t%in_t src(1)')
                                       call lock2%set   ('interp_t%in_t src(2)')
  do i4=1,source%nv
  do i3=ls(3),us(3)
  do i2=ls(2),us(2)
  do i1=ls(1),us(1)
    in2(i1,i2,i3,i4) = in(i1,i2,i3,i4,jt(2),1)
  end do
  end do
  end do
  end do
                                       call lock2%unset ('interp_t%in_t src(2)')
  !...........................................................................
  if (target%id == io%id_debug) then
    print '("debug:",2(i6,":",i2.2,i3),f10.6,2(i4,f10.6,f7.3))', &
      target%id, target%level, target%istep, &
      source%id, source%istep, source%level, &
      target%time, &
      jt(1), source%t(jt(1)), pt(1), &
      jt(2), source%t(jt(2)), pt(2)
  end if
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
                                                   call lock3%set ('interp tgt')
  if (present(lt)) then
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,9l1)') &
        wallclock(), pad_trace(1:indent)// &
        'interp_t%in_t: tgt,src,jt,pt,ls,us,lt,ut,usg=', &
        target%id, source%id, jt, pt, ls, us, lt, ut, &
        source%unsigned
    end if
    if (unit == stdout) then
      write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,9l1)') &
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
        minval( in1(:,:,:,1)), &
        maxval( in1(:,:,:,1)), &
        minval( in2(:,:,:,1)), &
        maxval( in2(:,:,:,1)), &
        minval( in1(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
        maxval( in1(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
        minval( in2(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
        maxval( in2(ls(1):us(1),ls(2):us(2),ls(3):us(3),1))
      flush (stdlog)
    end if
    !.........................................................................
    id = target%idx%d
    if (any(target%pervolume)) then
      allocate (d(target%gn(1),target%gn(2),target%gn(3)))
    end if
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
        if (target%unsigned(i4)) then
          do i2=ls(2),us(2)
            !$omp simd
            do i1=ls(1),us(1)
              out(i1+o(1),i2+o(2),i3+o(3),i4) = &
              exp(log(in1(i1,i2,i3,i4))*pt(1) &
                + log(in2(i1,i2,i3,i4))*pt(2))
            end do
          end do
          !-------------------------------------------------------------------
          ! If per-volume variables are used, save the density
          !-------------------------------------------------------------------
          if (i4==id .and. any(target%pervolume)) then
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                  d(i1+o(1),i2+o(2),i3+o(3)) = &
                out(i1+o(1),i2+o(2),i3+o(3),i4)
              end do
            end do
          end if
        !---------------------------------------------------------------------
        ! Use per-volume interpolation
        !---------------------------------------------------------------------
        else if (target%pervolume(i4)) then
          do i2=ls(2),us(2)
            !$omp simd
            do i1=ls(1),us(1)
              out(i1+o(1),i2+o(2),i3+o(3),i4) = &
                d(i1+o(1),i2+o(2),i3+o(3))* &
                (in1(i1,i2,i3,i4)/in1(i1,i2,i3,id)*pt(1) &
               + in2(i1,i2,i3,i4)/in2(i1,i2,i3,id)*pt(2))
            end do
          end do
        !---------------------------------------------------------------------
        ! Signed variables, or safe positive, shifted in index
        !---------------------------------------------------------------------
        else
          do i2=ls(2),us(2)
            !$omp simd
            do i1=ls(1),us(1)
              out(i1+o(1),i2+o(2),i3+o(3),i4) = in1(i1,i2,i3,i4)*pt(1) &
                                              + in2(i1,i2,i3,i4)*pt(2)
            end do
          end do
        end if
      end do
    end do
    if (any(target%pervolume)) then
      deallocate (d)
    end if
    n = product(us-ls+1)
  !---------------------------------------------------------------------------
  ! Unshifted interpolation in time
  !---------------------------------------------------------------------------
  else
    o = 0
    !.........................................................................
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,3i3),1x,9l1)') &
        wallclock(), &
        pad_trace(1:indent)//'interp_t%in_t: tgt,src,jt,pt,ls,us,usg=', &
        target%id, source%id, jt, pt, ls, us, &
        source%unsigned
    end if
    if (unit == stdout) then
      write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,3i3),1x,9l1)') &
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
            if (verbose > 0) then
              do i1=ls(1),us(1)
                if (in1(i1,i2,i3,i4) <= 0.0 .or. &
                    in2(i1,i2,i3,i4) <= 0.0) then
                  write (stdlog,'(a,2(i6,":",i2.2,3i4))') &
                    'interp_t%in_t: tgt,src=', &
                    target%id, target%key, source%id, source%key
                end if
                out(i1,i2,i3,i4) = &
                  exp(log(in1(i1,i2,i3,i4))*pt(1) &
                    + log(in2(i1,i2,i3,i4))*pt(2))
                v_min = min(v_min,out(i1,i2,i3,i4))
              end do
            else
              !$omp simd
              do i1=ls(1),us(1)
                out(i1,i2,i3,i4) = &
                  exp(log(in1(i1,i2,i3,i4))*pt(1) &
                    + log(in2(i1,i2,i3,i4))*pt(2))
                v_min = min(v_min,out(i1,i2,i3,i4))
              end do
            end if
          end do
        else
          do i2=ls(2),us(2)
            !$omp simd
            do i1=ls(1),us(1)
              out(i1,i2,i3,i4) = in1(i1,i2,i3,i4)*pt(1) &
                               + in2(i1,i2,i3,i4)*pt(2)
              v_min = min(v_min,out(i1,i2,i3,i4))
            end do
          end do
        end if
      end do
      if (i4==1) d_min = v_min
    end do
    n = product(us-ls+1)
  end if
                                                 call lock3%unset ('interp tgt')
  !...........................................................................
  if (verbose > 1 .and. present(lt)) then
    wc = wallclock()-wc
    n = product(us-ls+1)
    !lt = ls+o
    write (stdlog,'(a,f6.3,1p,6(2x,3e11.2))') &
    'interp_t%in_t(2): mus/p,min,max,min,max=', wc/n, &
    minval(in1(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
    maxval(in1(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
    minval(in2(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
    maxval(in2(ls(1):us(1),ls(2):us(2),ls(3):us(3),1)), &
    minval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),1)), &
    maxval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),1))
    if (verbose > 2) then
      call check_positive (target, out(:,:,:,1), lt, ut, 'in_t: out')
    end if
  end if
  deallocate (in1, in2)
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
!> Interpolate in x over the relevant target region.
!==============================================5================================
SUBROUTINE in_x (tgt, lt, ut, src, ls, us, om, in, out, all, only)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), nf(2), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: x, wc, wtimer, sm_p, d(3)
  real(8), allocatable     :: dst(:), dss(:)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_x', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Staggering offsets
  !-----------------------------------------------------------------------------
  allocate (dst(tgt%mv), dss(src%mv))
  dst = 0
  dss = 0
  dst(1:tgt%nv) = tgt%mesh(1)%h*tgt%mesh(1)%d
  dss(1:tgt%nv) = src%mesh(1)%h*src%mesh(1)%d
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
        x = tm%p + tm%r(i1) + dst(i4)                   ! x-position
        p = (x-sm_p-dss(i4))/sm%d + sm%o                ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(1)),us(1)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
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
        x = tm%p + tm%r(i1) + dst(i4)                   ! x-position
        p = (x-sm_p-dss(i4))/sm%d + sm%o                ! float source index
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
  if (verbose > 1 .and. i4==1) then
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
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: y, wc, wtimer, sm_p, d(3)
  real(8), allocatable     :: dst(:), dss(:)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_y', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Staggering offsets
  !-----------------------------------------------------------------------------
  allocate (dst(tgt%mv), dss(src%mv))
  dst = 0
  dss = 0
  dst(1:tgt%nv) = tgt%mesh(2)%h*tgt%mesh(2)%d
  dss(1:tgt%nv) = src%mesh(2)%h*src%mesh(2)%d
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
        y = tm%p + tm%r(i2) + dst(i4)                   ! y-position
        p = (y-sm_p-dss(i4))/sm%d + sm%o                ! float source index
        j = floor(p)                                    ! integer source index
        j = min(max(j,ls(2)),us(2)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
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
        y = tm%p + tm%r(i2) + dst(i4)                   ! y-position
        p = (y-sm_p-dss(i4))/sm%d + sm%o                ! float source index
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
  if (verbose > 1 .and. i4==1) then
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
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: z, wc, wtimer, sm_p, d(3)
  real(8), allocatable     :: dst(:), dss(:)
  class(mesh_t), pointer   :: tm, sm
  logical, dimension(3)    :: outside
  integer, save            :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_z', itimer=itimer)
  !.............................................................................
  if (verbose > 1) then
    wc = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Staggering offsets
  !-----------------------------------------------------------------------------
  allocate (dst(tgt%mv), dss(src%mv))
  dst = 0
  dss = 0
  dst(1:tgt%nv) = tgt%mesh(3)%h*tgt%mesh(3)%d
  dss(1:tgt%nv) = src%mesh(3)%h*src%mesh(3)%d
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
        z = tm%p + tm%r(i3) + dst(i4)                   ! z-position
        p = (z-sm_p-dss(i4))/sm%d + sm%o                ! float source index
        j = floor(p)                                    ! integer src index
        j = min(max(j,ls(3)),us(3)-1)                   ! make sure inside src
        p = p - j                                       ! 2nd fraction
        p = min(max(p,0.0),1.0)                         ! no extrapolation
        q = 1.0 - p                                     ! 1st fraction
        do i2=l(2),u(2)
          !$omp simd
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
        z = tm%p + tm%r(i3) + dst(i4)                   ! z-position
        p = (z-sm_p-dss(i4))/sm%d + sm%o                ! float source index
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
        end if ! verbose > 1 .and. i4 == 1
        !.......................................................................
      end do ! i3
      if (i4 == 1) then 
        d_min = v_min 
      else 
        d_min = minval(d)
      end if 
    end if ! .not. all  => no interior
  end do ! i4
  !.............................................................................
  if (verbose > 1 .and. i4==1) then
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
!> Not interpolating in time, but simply copying over the relevant source region
!===============================================================================
SUBROUTINE no_t (target, source, ls, us, in, out, lt, ut, only)
  class(patch_t)                  :: source, target
  integer, dimension(3)           :: ls, us
  real                            :: in (:,:,:,:,:,:)
  real                            :: out(:,:,:,:,:,:)
!  intent(in)                      :: in
!  contiguous                      :: in, out
  integer, dimension(3), optional :: lt, ut
  integer, optional               :: only
  !.............................................................................
  integer, dimension(3)           :: l, u, o
  integer                         :: i1, i2, i3, i4, i5, n, it, id, m(5)
  integer                         :: iv0, iv1
  real(8)                         :: wc
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: lock1, lock2
  integer, save                   :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%no_t', itimer=itimer)
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
  ! Lock access to source while copying
  !-----------------------------------------------------------------------------
  lock1 => source%mem_lock(1)
  lock2 => target%mem_lock(1)
  !...........................................................................
  if (target%id == io%id_debug) then
    print '("debug:",2(i6,":",i2.2,i3),f10.6,2(i4,f10.6,f7.3))', &
      target%id, target%level, target%istep, &
      source%id, source%istep, source%level, &
      target%time, &
      1, source%time, 1
  end if
  !...........................................................................
  if (verbose > 1) then
    write (stdlog,'(a,6i3,3x,4i3)') &
    'interp_t%no_t: shapes=', shape(in), shape(out)
  end if
  !---------------------------------------------------------------------------
  ! Lock access to source while interpolating
  !---------------------------------------------------------------------------
  if (associated(source%lock)) then
    lock1 => source%mem_lock(1)
  end if
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = target%nv
  end if
  !---------------------------------------------------------------------------
  ! Shifted from source index space
  !---------------------------------------------------------------------------
                                                   call lock2%set ('interp tgt')
  if (present(lt)) then
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2(1x,2(1x,3i3)),1x,9l1)') &
        wallclock(), pad_trace(1:indent)// &
        'interp_t%no_t: tgt,src,ls,us,lt,ut,usg=', &
        target%id, source%id, ls, us, lt, ut, &
        source%unsigned
    end if
    !.........................................................................
    if (verbose > 2) then
      write (stdlog,'(a,1p,4(2x,2e11.2))') &
        'interp_t%no_t: min/max=', &
        minval(in(:,:,:,:,1,:)), &
        maxval(in(:,:,:,:,1,:)), &
        minval(in(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,1,:)), &
        maxval(in(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,1,:))
      flush (stdlog)
    end if
    !.........................................................................
    id = target%idx%d
    o = lt-ls
    do i5=1,source%nw
    do i4=iv0,iv1                                         ! loop over variables
    do i3=ls(3),us(3)
    do i2=ls(2),us(2)
    !$omp simd
    do i1=ls(1),us(1)
      out(i1+o(1),i2+o(2),i3+o(3),i4,1,i5) = in(i1,i2,i3,i4,1,i5)
    end do
    end do
    end do
    end do
    end do
    n = product(us-ls+1)
  !---------------------------------------------------------------------------
  ! Unshifted copy in time
  !---------------------------------------------------------------------------
  else
    o = 0
    !.........................................................................
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2(1x,3i3),1x,9l1)') &
        wallclock(), &
        pad_trace(1:indent)//'interp_t%no_t: tgt,src,ls,us,usg=', &
        target%id, source%id, ls, us, &
        source%unsigned
    end if
    !.........................................................................
    do i5=1,source%nw
    do i4=iv0,iv1                                         ! loop over variables
    do i3=ls(3),us(3)
    do i2=ls(2),us(2)
    !$omp simd
    do i1=ls(1),us(1)
      out(i1,i2,i3,i4,1,i5) = in(i1,i2,i3,i4,1,i5)
    end do
    end do
    end do
    end do
    end do
    n = product(us-ls+1)
  end if
                                                 call lock2%unset ('interp tgt')
  !...........................................................................
  if (verbose > 1 .and. present(lt)) then
    wc = wallclock()-wc
    !lt = ls+o
    write (stdlog,'(a,f6.3,1p,6(2x,3e11.2))') &
    'interp_t%in_t(2): mus/p,min,max,min,max=', wc/n, &
    minval(in (ls(1):us(1),ls(2):us(2),ls(3):us(3),:,1,:)), &
    maxval(in (ls(1):us(1),ls(2):us(2),ls(3):us(3),:,1,:)), &
    minval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),:,1,:)), &
    maxval(out(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3),:,1,:))
  end if
  !.............................................................................
  if (verbose > 1) then
    ntimer = ntimer + n
    wtimer = wtimer + wc
    if (ntimer > 1e6) then
      write (stdlog,'(a,6p,f7.3)') 'no_t: average mus/pt=', wtimer/ntimer
      ntimer = 0
    end if
  end if
  !.............................................................................
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE no_t

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

END MODULE interp_mod
