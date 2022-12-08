!===============================================================================
!> Handles intepolations at the face between two tasks.  This module
!> services guards_mod, and is replaceable via the makefile macro INTERP.
!define DEBUG
!> The DEBUG macro turns on code sections that may slow down normal operation
!===============================================================================
MODULE interp_mod
  USE iso_fortran_env, only: int8
  USE trace_mod
  USE io_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE task_mod
  USE patch_mod
  USE mesh_mod
  USE shared_mod
  USE stagger_mod
  USE os_mod
  USE vector_ops
  USE bits_mod
  implicit none
  private
  intrinsic :: log
  type, public:: meshp_t
    class(mesh_t), pointer:: p => null()
  end type
  type, public:: interp_t
    class(task_t), pointer:: source => null()    ! source patch
    class(task_t), pointer:: target => null()    ! target patch
    real, pointer, contiguous:: r(:,:,:,:)
    type(lock_t):: lock
    logical:: use_staggered=.false.
    procedure(indices1), nopass, pointer:: indices => null()
  contains
    procedure:: init
    procedure:: deposit_to_target
    procedure:: deposit_to_target_stag
    procedure:: from_source
    procedure:: from_source_stag
    procedure:: from_only
    procedure, nopass :: indices1
    procedure, nopass :: indices2
    procedure:: indices_stag
    procedure:: to_target
    procedure:: to_target_stag
    procedure, nopass:: in_t
    procedure, nopass:: no_t
    procedure, nopass:: in_x
    procedure, nopass:: in_y
    procedure, nopass:: in_z
  end type
  integer:: verbose = 0
  logical:: detailed_timer = .false.
  real(8), parameter:: eps=0.5001d0
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
  logical, save:: use_staggered=.false.
  integer, save:: index_pick = 2
  namelist /interp_params/ verbose, detailed_timer, use_staggered, index_pick
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
  !$omp atomic write
  self%use_staggered = use_staggered
  select case(index_pick)
  case(1)
    self%indices => indices1
  case(2)
    self%indices => indices2
  case default
    self%indices => indices2
  end select   
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
  integer(2), pointer         :: filled(:,:,:)
  logical, optional           :: all_src
  !.............................................................................
  real(8), pointer, contiguous   :: p(:,:,:,:)
  real(8), pointer, contiguous:: tr1(:), tr2(:), tr3(:)
  real, pointer, contiguous   :: tm(:,:,:,:), sm(:,:,:,:)
  real, pointer, contiguous   :: sp(:,:,:,:), tp(:,:,:,:)
  real, pointer, contiguous   :: v1(:,:,:), v2(:,:,:), v3(:,:,:)
  real                        :: pt(2), level
  real(8)                     :: p1, p2, q1, q2, p3, q3, v(3), v_shift, d_omega
  real(8)                     :: r1, r2, r, d(3), e(3,3), r_tgt(3), vr
  real(8)                     :: p_src(3), p_tgt(3)
  integer                     :: i1, i2, i3, j1 ,j2, j3, iv, j, it(2), idx_es
  integer, dimension(3)       :: l, u, ls, us, m(4)
  logical, pointer, contiguous:: available(:,:,:)
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
  all_src_loc = shared%option (all_src)
  !-----------------------------------------------------------------------------
  ! Combine handling of energy and entropy per unit mass
  !-----------------------------------------------------------------------------
  idx_es = max(src%idx%e,src%idx%s)  
  !-----------------------------------------------------------------------------
  ! the available grid holds logicals, telling us which cells are available to
  ! be filled from source
  !-----------------------------------------------------------------------------
  call os%allocate (available, m(1:3), 'available')
  available(:,:,:) = .false.
  !-----------------------------------------------------------------------------
  ! Allocate an array holding the source values in per-unit-mass, and an array
  ! for the resulting target values
  !-----------------------------------------------------------------------------
  call os%allocate (tp, m, 'tp')
  call os%allocate (sp, m, 'sp')
  call src%time_interval (tgt, tgt%time, it, pt)
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from tgt to src.  On return d is the
  ! distance from tgt to src, expressed in tgt coordinates, while the matrix
  ! e has the unit vectors of the source coordinate system, expressed in the
  ! target coordinate system, in the columns, and the target unit vectors
  ! expressed in the source coordinate system as rows.  To tranform a vector
  ! from the tgt to the src system, one should therefore multiply the matrix
  ! from the right onto the vectors expressed in target coordinates.
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, p=d, e=e)
  ! call                                src%lock%set_read ('interp_t%from_source')
  !-----------------------------------------------------------------------------
  ! Floating interpolation coordinates p(:,:,:,1:3) -- also returns the logical
  ! mask that exclude points that are too far outside
  !-----------------------------------------------------------------------------
  call os%allocate (p, [src%gn,3], 'p')
  call self%indices (tgt, lt, ut, src, p, available, d, e)
  call source_indices (lt, ut, src, p, available, all_src_loc, ls, us)
  ! call                              src%lock%unset_read ('interp_t%from_source')
  !-----------------------------------------------------------------------------
  ! Convert to per unit mass source values
  !-----------------------------------------------------------------------------
  call to_per_unit_mass (src, ls, us, sp, it, pt)
  !-----------------------------------------------------------------------------
  ! Transform velocities to {r,phi,z} components
  !-----------------------------------------------------------------------------
  v1 => sp(:,:,:,src%idx%px)
  v2 => sp(:,:,:,src%idx%py)
  v3 => sp(:,:,:,src%idx%pz)
  call src%to_radial (ls, us, v1, v2, v3)
  !-----------------------------------------------------------------------------
  ! Handle all_src cases
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if
  !-----------------------------------------------------------------------------
  ! Tri-linear interpolation in per-unit-mass variables
  !-----------------------------------------------------------------------------
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        !-----------------------------------------------------------------------
        ! The available mask allows only points at most 0.5 cell outside ROA
        !-----------------------------------------------------------------------
        if (.not.available(i1,i2,i3)) &
          cycle
        j1 = p(i1,i2,i3,1)
        j2 = p(i1,i2,i3,2)
        j3 = p(i1,i2,i3,3)
        j1 = min(max(j1,l(1)),u(1)-1)
        j2 = min(max(j2,l(2)),u(2)-1)
        j3 = min(max(j3,l(3)),u(3)-1)
        p1 = p(i1,i2,i3,1) - j1
        p2 = p(i1,i2,i3,2) - j2
        p3 = p(i1,i2,i3,3) - j3
        q1 = 1. - p1
        q2 = 1. - p2
        q3 = 1. - p3
        level = min(p1,p2,p3,q1,q2,q3,0.0_8)
        !-----------------------------------------------------------------------
        ! Among remaining candidate source points, keep the one (or the last
        ! one), that has either the smallest or no extrapolation
        !-----------------------------------------------------------------------
        level = 10*(level + src%level)
#ifdef DEBUG
        if (verbose>0 .and. i1>=26 .and. i2<=2 .and. i3==14) then
          write (stdout,'(a,3f7.2,2(2x,3i3),f8.2,i5,3f7.3)') &
            'debug:', p1,p2,p3,i1,i2,i3,j1,j2,j3,level,filled(i1,i2,i3),d
          !--- source point in source coordinates ---
          p_src = [src%mesh(1)%r(j1)+p1*src%ds(1) + src%centre_nat(1), &
                   src%mesh(2)%r(j2)+p2*src%ds(2), &
                   src%mesh(3)%r(j3)+p3*src%ds(2)]
          !--- scalar products with tgt unit vectors in src coordinates  ---
          p_tgt = matmul (e, p_src)
          !--- subtract radius and transform to index space ---
          p_tgt(1) = p_tgt(1) - tgt%centre_nat(1)
          !--- should be equal to [i1,i2,i3]
          p_tgt = p_tgt/tgt%mesh%d + tgt%mesh%o
        end if
#endif
        if (filled(i1,i2,i3) < int(level)) then
#ifdef DEBUG
          !.....................................................................
          if (verbose > 0 .and. filled(i1,i2,i3) > 0.0) &
            write (stdout,'(3i3,i5,f6.1,2x,a)') &
              i1,i2,i3,filled(i1,i2,i3),level,'replace'
          !.....................................................................
#endif
          filled(i1,i2,i3) = int(level)
          !$omp simd
          do iv=1,tgt%nv
            tp(i1,i2,i3,iv) = &
              q3*(q2*(q1*sp(j1,j2  ,j3  ,iv) + p1*sp(j1+1,j2  ,j3  ,iv)) + &
                  p2*(q1*sp(j1,j2+1,j3  ,iv) + p1*sp(j1+1,j2+1,j3  ,iv))) + &
              p3*(q2*(q1*sp(j1,j2  ,j3+1,iv) + p1*sp(j1+1,j2  ,j3+1,iv)) + &
                  p2*(q1*sp(j1,j2+1,j3+1,iv) + p1*sp(j1+1,j2+1,j3+1,iv)))
          end do
        else
          available(i1,i2,i3)=.false.
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Add the velocity difference to the phi-component.  This loop is faster
  ! vectorized than if using available conditions, since vectorization gives
  ! much more speedup than the fraction of internal cells that could be skipped.
  !-----------------------------------------------------------------------------
  v1 => tp(:,:,:,tgt%idx%px)
  v2 => tp(:,:,:,tgt%idx%py)
  v3 => tp(:,:,:,tgt%idx%pz)
  tr1 => tgt%mesh(1)%r
  tr2 => tgt%mesh(2)%r
  tr3 => tgt%mesh(3)%r
  if (src%omega(3) /= tgt%omega(3)) then
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          !---------------------------------------------------------------------
          ! The coordinate velocity of each cell in a patch is omega*r
          !---------------------------------------------------------------------
          r1 = tr1(i1) + tgt%centre_nat(1)
          r2 = tr2(i2)
          v2(i1,i2,i3) = &
          v2(i1,i2,i3) + (src%omega(3)-tgt%omega(3))*sqrt(r1**2+r2**2)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Convert back to cartesian velocity components and conservative variables.
  ! The code in the loop is conditional so does not vectorize.
  !-----------------------------------------------------------------------------
  call tgt%to_cartesian (lt, ut, v1, v2, v3)
  !.............................................................................
  if (tgt%id == io%id_debug) &
    write (stdlog,'(a,2(i6,":",i2.2,3i4,i5,f8.4,2(1x,3i3),2(1x,3f8.4)),f8.3)') &
      'from_source:', &
      tgt%id, tgt%key, tgt%istep, tgt%time, lt, ut, tgt%position, tgt%centre_nat, &
      src%id, src%key, src%istep, src%time, ls, us, src%position, src%centre_nat, &
      maxval(v1(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)))
  !.............................................................................
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        if (tgt%outside(i1,i2,i3) .and. available(i1,i2,i3)) then
          !---------------------------------------------------------------------
          ! Mass density
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%d) = tp(i1,i2,i3,tgt%idx%d)
          !---------------------------------------------------------------------
          ! Convert from velocities to linear momenta
          !---------------------------------------------------------------------
          tm(i1,i2,i3,tgt%idx%px) = &
          tp(i1,i2,i3,tgt%idx%px) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%py) = &
          tp(i1,i2,i3,tgt%idx%py) * tm(i1,i2,i3,tgt%idx%d)
          tm(i1,i2,i3,tgt%idx%pz) = &
          tp(i1,i2,i3,tgt%idx%pz) * tm(i1,i2,i3,tgt%idx%d)
          !---------------------------------------------------------------------
          ! Energy variable in the target coordinate system.
          !---------------------------------------------------------------------
          if (tgt%etype == 'Entropy') then
            tm(i1,i2,i3, tgt%idx%s) = tp(i1,i2,i3,tgt%idx%s)
          else
            tm(i1,i2,i3,    idx_es) = tm(i1,i2,i3,tgt%idx%d)* &
                                      tp(i1,i2,i3,   idx_es)
          end if
        end if
      end do
    end do
  end do
  call os%deallocate (p, 'p')
  call os%deallocate (tp, 'tp')
  call os%deallocate (sp, 'sp')
  call os%deallocate (available, 'available')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE from_source

!===============================================================================
!> Interpolate in 3-D index space
!===============================================================================
SUBROUTINE from_source_stag (self, tgt, lt, ut, src, filled, all_src)
  class(interp_t)             :: self
  class(patch_t)              :: tgt, src
  integer, dimension(3)       :: lt, ut
  integer(2), pointer         :: filled(:,:,:)
  optional                    :: filled
  logical, optional           :: all_src
  !.............................................................................
  real(8), pointer, contiguous   :: p(:,:,:,:), p_stag(:,:,:,:,:)
  real, pointer, contiguous   :: tm(:,:,:,:), sm(:,:,:,:)
  real, pointer, contiguous   :: sp(:,:,:,:), tp(:,:,:,:)
  real, pointer, contiguous   :: vel_stag(:,:,:,:,:), B_stag(:,:,:,:,:)
  real                        :: pt(2)
  real(8)                     :: d(3), e(3,3)
  logical                     :: outside(3)
  integer                     :: i1, i2, j1 ,j2, iv, j, it(2), idx_es
  integer, dimension(3)       :: l, u, idx, ls, us, m(4)
  integer(2), pointer, contiguous :: filled_stag(:,:,:)
  logical, pointer, contiguous:: available(:,:,:), available_stag(:,:,:,:)
  logical, pointer, contiguous:: stag(:)
  logical                     :: all_src_loc, stag_vel, stag_B
  integer, pointer, contiguous:: idx_cent(:), idx_stag(:)
  integer, save               :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%from_source_stag', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Interpolate in time, in per unit mass variables, leaving source values in
  ! sp(:,:,:,:)
  !-----------------------------------------------------------------------------
  tm => tgt%mem(:,:,:,:,tgt%it,1)
  sm => src%mem(:,:,:,:,src%it,1)
  m = shape(sm)
  call os%allocate (sp, m, 'sp') ! array, holding per unit mass quantities in source
  all_src_loc = shared%option (all_src)
  !-----------------------------------------------------------------------------
  ! Combine handling of energy and entropy per unit mass
  !-----------------------------------------------------------------------------
  idx_es = max(src%idx%e,src%idx%s)  
  ! ----------------------------------------------------------------------------
  ! Book-keeping necessities.
  ! Go through nv and identify centered and staggered quantities
  ! borrowing iX and jX here.
  ! j1 - number of centered quantities
  ! j2 - number of staggered quantities
  ! ----------------------------------------------------------------------------
  call os%allocate(stag, src%nv)
  stag     = .false.
  stag_vel = .false.
  stag_B   = .false.
  j1=0; j2=0
  do i1=1,src%nv
    do i2=1,3
      if (tgt%mesh(i2)%h(i1)/=0.0) stag(i1)=.true.
    end do
    if (stag(i1)) then
      j2=j2+1
    else
      j1=j1+1
    end if
  end do
  call os%allocate(idx_cent, j1, 'idx_cent')
  call os%allocate(idx_stag, j2, 'idx_stag')
  j1=1; j2=1
  do i1=1,src%nv
    if (stag(i1)) then
      idx_stag(j2)=i1
      j2 = j2+1
      if (i1==src%idx%px) stag_vel = .true.
      if (i1==src%idx%bx) stag_B   = .true.
    else
      idx_cent(j1)=i1
      j1 = j1+1
    end if
  end do
  !-----------------------------------------------------------------------------
  ! the available grid holds logicals, telling us, which cells are available to
  ! be filled from source, for cell-centred quantities.
  !-----------------------------------------------------------------------------
  call os%allocate (available, m(1:3), 'available')
  available(:,:,:) = .false.
  
  ! ----------------------------------------------------------------------------
  ! arrays, holding quantities in target's frame of reference. Vectors are not
  ! de-rotated in these arrays.
  ! ----------------------------------------------------------------------------
  call os%allocate (tp, [m(1:3),size(idx_cent,1)], 'tp')
  if (stag_vel) call os%allocate (vel_stag, [m(1:3),3,3], 'vel_stag')
  if (stag_B)   call os%allocate (  B_stag, [m(1:3),3,3], 'B_stag')
  !-----------------------------------------------------------------------------
  ! Duplicate filled array for consistency. The filled array will still be
  ! controled by the cell-centered information
  !-----------------------------------------------------------------------------
  if (stag_B.or.stag_vel) then
    if (present(filled)) then
      call os%allocate (filled_stag, shape(filled), 'filled_stag')
    end if
    filled_stag=filled
  end if
  ! ----------------------------------------------------------------------------
  ! Get interpolation in time slots and weights
  ! ----------------------------------------------------------------------------
  call src%time_interval (tgt, tgt%time, it, pt)
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from tgt to src
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, p=d, e=e)
  !-----------------------------------------------------------------------------
  ! Interpolation coordinates and weights p(:,:,:,1:3)
  !-----------------------------------------------------------------------------
  call os%allocate (p, [m(1:3),3], 'p')
  call self%indices (tgt, lt, ut, src, p, available, d, e)
  call source_indices (lt, ut, src, p, available, all_src_loc, ls, us)
  ! ----------------------------------------------------------------------------
  ! the available_stag and p_stag holds corresponding logicals 
  ! and interpolation coordinates/weights for staggered quantities.
  ! ----------------------------------------------------------------------------
  if (stag_vel .or. stag_B) then
    call os%allocate (available_stag, [m(1:3),3], 'avail_stag')
    call os%allocate (p_stag, [m(1:3),3,3], 'p_stag')
    available_stag=.false.
    call self%indices_stag (tgt, lt, ut, src, p_stag, available_stag, d, e)
    call source_indices_stag (lt, ut, src, p_stag, available_stag, &
                              all_src_loc, ls, us)
  end if
  ! ----------------------------------------------------------------------------
  ! convert source mem to per unit mass. 
  ! NB!: there are specific cases, where this conversion is not done
  ! ----------------------------------------------------------------------------
  call to_per_unit_mass (src, ls, us, sp, it, pt)
  ! ----------------------------------------------------------------------------
  ! Interpolate to target reference frame
  ! ----------------------------------------------------------------------------
  call self%to_target(src, tp, sp, lt, ut, p, available, &
                      idx_cent, all_src_loc, filled)
  if (stag_vel) then
    idx(1)=src%idx%px
    idx(2)=src%idx%py
    idx(3)=src%idx%pz
    call self%to_target_stag(src, vel_stag, sp, lt, ut, p_stag, available_stag, &
                             idx, all_src_loc, filled_stag)
  end if
  if (stag_B) then
    idx(1)=src%idx%bx
    idx(2)=src%idx%by
    idx(3)=src%idx%bz
    call self%to_target_stag(src, B_stag, sp, lt, ut, p_stag, available_stag, &
                             idx, all_src_loc, filled_stag)
  end if

  !-----------------------------------------------------------------------------
  ! Deposit interpolated values in target mem, taking into account that the
  ! tp values are "per unit mass".
  ! once again, staggered vectors are done separately.
  !  We pass the filled arrays only to centered subroutine.
  !-----------------------------------------------------------------------------
  call self%deposit_to_target(tgt, src, tp, lt, ut, idx_cent, &
                available, all_src_Loc)
  if (stag_vel) then
    idx(1)=tgt%idx%px
    idx(2)=tgt%idx%py
    idx(3)=tgt%idx%pz
    call self%deposit_to_target_stag(tgt, src, vel_stag, lt, ut, idx, &
                available_stag, all_src_loc)
  end if
  if (stag_B) then
    idx(1)=tgt%idx%bx
    idx(2)=tgt%idx%by
    idx(3)=tgt%idx%bz
    call self%deposit_to_target_stag(tgt, src, B_stag, lt, ut, idx, &
                available_stag, all_src_loc)
  end if
  !
  call os%deallocate (p, 'p')
  call os%deallocate (tp, 'tp')
  call os%deallocate (sp, 'sp')
  call os%deallocate (available, 'available')
  call os%deallocate (stag, 'stag')
  call os%deallocate(idx_stag, 'idx_stag')
  call os%deallocate(idx_cent, 'idx_cent')
  if (stag_vel .or. stag_B) then
    call os%deallocate (available_stag, 'avail_stag')
    call os%deallocate (p_stag, 'p_stag')
    if (stag_vel) call os%deallocate (vel_stag, 'vel_stag')
    if (stag_B)   call os%deallocate (B_stag, 'B_stag')
    if (present(filled)) then
      call os%deallocate (filled_stag, 'filled_stag')
    end if
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE from_source_stag

!===============================================================================
!> Interpolate to target reference frame, using tri-linear interpolation
!===============================================================================
SUBROUTINE to_target(self, src, tp, sp, lt, ut, p, available, &
                      idx_cent, all_src, filled)
  class(interp_t)           :: self
  class(patch_t)            :: src
  real, pointer, contiguous :: sp(:,:,:,:), tp(:,:,:,:)
  real(8), pointer, contiguous:: p(:,:,:,:)
  integer, dimension(3)     :: lt, ut
  logical, pointer          :: available(:,:,:)
  integer, pointer          :: idx_cent(:)
  logical                   :: all_src
  integer(2), pointer, optional :: filled(:,:,:)
  ! ----------------------------------------------------------------------------
  integer       :: l(3), u(3), iv, nv, idx
  integer       :: i1, i2, i3
  integer       :: j1, j2, j3
  real          :: level
  real(8)       :: p1, q1, p2, q2, p3, q3
  integer       :: m_extrap
  integer, save :: n_extrap=0
  ! ............................................................................

  !-----------------------------------------------------------------------------
  ! Handle all_src cases
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if
  nv=size(idx_cent,1)
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
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
        q1 = 1.0 - p1
        q2 = 1.0 - p2
        q3 = 1.0 - p3
        level = min(p1,p2,p3,q1,q2,q3,0.0_8)
        !---------------------------------------------------------------------
        ! Among remaining candidate source points, keep the one (or the last
        ! one), that has either the smallest or no extrapolation
        !---------------------------------------------------------------------
        level = 10*(level + src%level)
        if (filled(i1,i2,i3) < int(level)) then
          filled(i1,i2,i3) = int(level)
          do iv=1,nv
            idx=idx_cent(iv)
            tp(i1,i2,i3,iv) = &
              q3*(q2*(q1*sp(j1,j2  ,j3  ,idx) + p1*sp(j1+1,j2  ,j3  ,idx)) + &
                  p2*(q1*sp(j1,j2+1,j3  ,idx) + p1*sp(j1+1,j2+1,j3  ,idx))) + &
              p3*(q2*(q1*sp(j1,j2  ,j3+1,idx) + p1*sp(j1+1,j2  ,j3+1,idx)) + &
                  p2*(q1*sp(j1,j2+1,j3+1,idx) + p1*sp(j1+1,j2+1,j3+1,idx)))
          end do
        else
          available(i1,i2,i3)=.false.
        end if
      end do
    end do
  end do

END SUBROUTINE to_target

!===============================================================================
!> Interpolate to target reference frame, using tri-linear interpolation
!> for staggered quantities.
!> FIXME: perhaps don't prevent extrapolation within small range...
!===============================================================================
SUBROUTINE to_target_stag(self, src, tp, sp, lt, ut, p, available, &
                      idx, all_src, filled)
  class(interp_t)           :: self
  class(patch_t)            :: src
  real, pointer, contiguous :: sp(:,:,:,:), tp(:,:,:,:,:)
  real(8), pointer, contiguous:: p(:,:,:,:,:)
  integer, dimension(3)     :: lt, ut
  logical, pointer          :: available(:,:,:,:)
  integer                   :: idx(3)
  logical                   :: all_src
  integer(2), pointer, optional :: filled(:,:,:)
  ! ----------------------------------------------------------------------------
  real, pointer, contiguous :: vect(:,:,:,:), vect_cent(:,:,:,:), vect_restag(:,:,:,:,:)
  real, pointer, dimension(:,:,:), contiguous :: Vx_y, Vx_z, Vy_x, Vy_z, Vz_x, Vz_y
  real, pointer, dimension(:,:,:), contiguous :: x,y,z
  integer       :: l(3), u(3), iv, n(4), m(5)
  integer       :: i1, i2, i3
  integer       :: j1, j2, j3
  real          :: level
  real(8)       :: p1, q1, p2, q2, p3, q3
  integer       :: m_extrap
  integer, save :: n_extrap=0
  ! ............................................................................
  ! ----------------------------------------------------------------------------
  ! Point to the vector and allocate a dataset, where the staggered vector will
  ! be spread. We don't hardcode the choice of index in this subroutine because
  ! the subroutine doesn't really know what vector it is dealing with.
  ! ----------------------------------------------------------------------------
  vect => sp(:,:,:,idx(1):idx(3))
  n=shape(vect)
  m=shape(tp)
  call os%allocate(vect_cent, n, 'vect_cent')
  call os%allocate(vect_restag, m, 'vect_restag')
  ! ----------------------------------------------------------------------------
  ! Cell-center the vectors
  ! ----------------------------------------------------------------------------
  vect_cent=up(vect)
  x => vect_cent(:,:,:,1)
  y => vect_cent(:,:,:,2)
  z => vect_cent(:,:,:,3)
  ! ----------------------------------------------------------------------------
  ! Re-stagger staggered vectors to have them available for rotation. E.g.:
  ! bx_x = xup(bx))
  ! bx_y = ydn(bx_x)
  ! bx_z = zdn(bx_x)
  ! by_y = yup(by)
  ! by_x = xdn(by_y)
  ! by_z = zdn(by_y)
  ! bz_z = zup(bz)
  ! bz_x = xdn(bz_z)
  ! bz_y = ydn(bz_z)
  ! vect  has dimension (x|y|z|3)
  ! tp will have dimension (x|y|z|3|3), where last dimension is the X, Y Z vector 
  ! component at a given staggered original vector field
  ! ----------------------------------------------------------------------------
  vect_restag(:,:,:,1,1)=vect(:,:,:,1)
  vect_restag(:,:,:,2,2)=vect(:,:,:,2)
  vect_restag(:,:,:,3,3)=vect(:,:,:,3)
  Vx_y => vect_restag(:,:,:,1,2)
  Vx_z => vect_restag(:,:,:,1,3)
  Vy_x => vect_restag(:,:,:,2,1)
  Vy_z => vect_restag(:,:,:,2,3)
  Vz_x => vect_restag(:,:,:,3,1)
  Vz_y => vect_restag(:,:,:,3,2)
  !
  Vx_y=ydn(x)
  Vx_z=zdn(x)
  !
  Vy_x=xdn(y)
  Vy_z=zdn(y)
  !
  Vz_x=xdn(z)
  Vz_y=ydn(z)
  !-----------------------------------------------------------------------------
  ! Handle all_src cases
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if
  ! ----------------------------------------------------------------------------
  ! Reinterpolate to target mesh
  ! ----------------------------------------------------------------------------
  do iv=1,3
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          if (.not.available(i1,i2,i3,iv)) &
            cycle
          !-------------------------------------------------------------------
          ! Skip interior points
          !-------------------------------------------------------------------
          j1 = p(i1,i2,i3,iv,1)
          j2 = p(i1,i2,i3,iv,2)
          j3 = p(i1,i2,i3,iv,3)
          j1 = min(max(j1,l(1)),u(1)-1)
          j2 = min(max(j2,l(2)),u(2)-1)
          j3 = min(max(j3,l(3)),u(3)-1)
          p1 = p(i1,i2,i3,iv,1) - j1
          p2 = p(i1,i2,i3,iv,2) - j2
          p3 = p(i1,i2,i3,iv,3) - j3
          q1 = 1.0 - p1
          q2 = 1.0 - p2
          q3 = 1.0 - p3
          level = min(p1,p2,p3,q1,q2,q3,0.0_8)
          !---------------------------------------------------------------------
          ! Among remaining candidate source points, keep the one (or the last
          ! one), that has either the smallest or no extrapolation
          !---------------------------------------------------------------------
          level = 10*(level + src%level)
          !...................................................................
          if (filled(i1,i2,i3) < int(level)) then
            ! this is commented out as the real filled block is filled
            ! in "to_target" subroutine. Ideally we might want an additional
            ! dimension for `filled`, so different IVs would have their own
            ! values
            !filled(i1,i2,i3) = level
              tp(i1,i2,i3,iv,:) = &
                q3*(q2*(q1*vect_restag(j1,j2  ,j3  ,:,iv) + p1*vect_restag(j1+1,j2  ,j3  ,:,iv)) + &
                    p2*(q1*vect_restag(j1,j2+1,j3  ,:,iv) + p1*vect_restag(j1+1,j2+1,j3  ,:,iv))) + &
                p3*(q2*(q1*vect_restag(j1,j2  ,j3+1,:,iv) + p1*vect_restag(j1+1,j2  ,j3+1,:,iv)) + &
                    p2*(q1*vect_restag(j1,j2+1,j3+1,:,iv) + p1*vect_restag(j1+1,j2+1,j3+1,:,iv)))
            !...................................................................
            if (verbose > 0 .and. filled(i1,i2,i3) > 0.0) &
              write (stdlog,*) i1,i2,i3,filled(i1,i2,i3),level,'replace'
          else
            available(i1,i2,i3,iv)=.false.
          end if
        end do
      end do
    end do
  end do

  ! ----------------------------------------------------------------------------
  ! cleanup
  ! ----------------------------------------------------------------------------
  call os%deallocate (vect_cent, 'vect_cent')
  call os%deallocate (vect_restag, 'vect_restag')
END SUBROUTINE to_target_stag

!===============================================================================
!> Deposits the cell-centered quantities to target. Keeps in mind, that some 
!> quantities are in per unit mass.
!===============================================================================
SUBROUTINE deposit_to_target(self, tgt, src, tp, lt, ut, idx_cent, available, all_src)
  class(interp_t)           :: self
  class(patch_t)            :: tgt, src
  real, pointer, contiguous :: tp(:,:,:,:)
  integer, dimension(3)     :: lt, ut
  logical, pointer          :: available(:,:,:)
  integer, pointer          :: idx_cent(:)
  logical:: all_src
  ! ----------------------------------------------------------------------------
  logical:: outside(3)
  integer:: l(3), u(3), idv, ndv
  integer:: i1, i2, i3, idx_es
  real   :: r1, r2, r
  real(8):: e(3,3), dist(3)
  real(8):: v(3), v_shift(3), d_omega
  real(8), pointer, contiguous:: tr1(:), tr2(:), tr3(:)
  real, pointer, contiguous   :: tm(:,:,:,:)
  logical:: cent_vel
  ! ............................................................................
  tm => tgt%mem(:,:,:,:,tgt%it,1)
  !-----------------------------------------------------------------------------
  ! Combine handling of energy and entropy per unit mass
  !-----------------------------------------------------------------------------
  idx_es = max(src%idx%e,src%idx%s)  
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from src to tgt
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, p=dist, e=e)
  !
  ! ----------------------------------------------------------------------------
  ! Check if velocities are cell-centered
  ! ----------------------------------------------------------------------------
  ndv=size(idx_cent,1)
  if (any(idx_cent==tgt%idx%px)) then
    cent_vel=.true.
  else
    cent_vel=.false.
  end if
  !-----------------------------------------------------------------------------
  ! Handle all_src cases
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if
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
          do idv=1,ndv
            if (idx_cent(idv)==tgt%idx%d) exit
          end do
          tm(i1,i2,i3,tgt%idx%d) = tp(i1,i2,i3,idv)
          !---------------------------------------------------------------------
          ! Cell-centered velocity, if applicable
          !
          ! The coordinate velocity of each cell in a patch is omega*r, where
          ! r is the local radius, and the velocity shift to add from source
          ! to target is the difference from target to source.  The direction
          ! of that velocity is along phi_hat
          !---------------------------------------------------------------------
          do idv=1,ndv
            if (idx_cent(idv)==tgt%idx%px) exit
          end do
          if (cent_vel) then
            if (tgt%is_clear(bits%static)) then
            r1 = tr1(i1)*tgt%ds(1) + tgt%centre_nat(1)
            r2 = tr2(i2)*tgt%ds(2)
            r = sqrt(r1**2+r2**2)
            d_omega = (src%omega(3)-tgt%omega(3))
            v_shift = d_omega*[-r2,+r1,0.0]
            else
              v_shift= 0.0
            end if
            !---------------------------------------------------------------------
            ! Get the src velocity interpolated to the tgt points
            !---------------------------------------------------------------------
            v(1) = tp(i1,i2,i3,idv  )
            v(2) = tp(i1,i2,i3,idv+1)
            v(3) = tp(i1,i2,i3,idv+2)
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
            v = matmul (e, v)
            !---------------------------------------------------------------------
            ! Add velocity shift, which is already in the tgt system
            !---------------------------------------------------------------------
            v = v + v_shift
            tm(i1,i2,i3,tgt%idx%px) = v(1)
            tm(i1,i2,i3,tgt%idx%py) = v(2)
            tm(i1,i2,i3,tgt%idx%pz) = v(3)
          end if
          !---------------------------------------------------------------------
          ! Total energy in the target coordinate system.  The tp(...idx%e) is
          ! the thermal energy per unit mass, interpolated to the target point,
          ! to which needs to be added the kinetic energy in the target frame,
          ! keeping in mind that the tm(...%px) values are still velocities
          !---------------------------------------------------------------------
          do idv=1,ndv
              if (idx_cent(idv)==idx_es) exit
          end do
          if (tgt%etype == 'total') then
            tm(i1,i2,i3,tgt%idx%e ) = tm(i1,i2,i3,tgt%idx%d)*( &
            tp(i1,i2,i3,      idv ) + 0.5*( &
            tm(i1,i2,i3,tgt%idx%px)**2 + &
            tm(i1,i2,i3,tgt%idx%py)**2 + &
            tm(i1,i2,i3,tgt%idx%pz)**2))
          else if (tgt%etype == 'Entropy') then
            tm(i1,i2,i3, tgt%idx%s) = tp(i1,i2,i3,idv)
          else
            tm(i1,i2,i3,    idx_es) = tm(i1,i2,i3,tgt%idx%d)* &
                                      tp(i1,i2,i3,      idv)
          end if
          !---------------------------------------------------------------------
          ! Convert from velocities to linear momenta
          !---------------------------------------------------------------------
          if (cent_vel) then
            tm(i1,i2,i3,tgt%idx%px) = &
            tm(i1,i2,i3,tgt%idx%px) * tm(i1,i2,i3,tgt%idx%d)
            tm(i1,i2,i3,tgt%idx%py) = &
            tm(i1,i2,i3,tgt%idx%py) * tm(i1,i2,i3,tgt%idx%d)
            tm(i1,i2,i3,tgt%idx%pz) = &
            tm(i1,i2,i3,tgt%idx%pz) * tm(i1,i2,i3,tgt%idx%d)
          end if
          !
        end if
      end do
    end do
  end do

END SUBROUTINE deposit_to_target

!===============================================================================
!> Deposits the staggered quantities to target. 
!===============================================================================
SUBROUTINE deposit_to_target_stag(self, tgt, src, vector, lt, ut, idx, available, all_src)
  class(interp_t)           :: self
  class(patch_t)            :: tgt, src
  real, pointer, contiguous :: vector(:,:,:,:,:)
  integer, dimension(3)     :: lt, ut, idx
  logical, pointer          :: available(:,:,:,:)
  logical :: all_src
  ! ----------------------------------------------------------------------------
  logical:: outside(3)
  integer:: l(3), u(3), iv
  integer:: i1, i2, i3
  real   :: r1, r2, r
  real(8):: e(3,3), dist(3)
  real(8):: v(3)
  real, pointer, contiguous   :: tm(:,:,:,:)
  ! ............................................................................
  tm => tgt%mem(:,:,:,:,tgt%it,1)
  !-----------------------------------------------------------------------------
  ! Get the distance and rotation matrix from src to tgt
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, p=dist, e=e)
  !
  !-----------------------------------------------------------------------------
  ! Handle all_src cases
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if

  do iv=1,3
    do i3=lt(3),ut(3)
      outside(3) = i3 < tgt%mesh(3)%li .or. i3 > tgt%mesh(3)%ui
      do i2=lt(2),ut(2)
        outside(2) = i2 < tgt%mesh(2)%li .or. i2 > tgt%mesh(2)%ui
        do i1=lt(1),ut(1)
          outside(1) = i1 < tgt%mesh(1)%li .or. i1 > tgt%mesh(1)%ui
          !-----------------------------------------------------------------------
          ! Skip interior points
          !-----------------------------------------------------------------------
          if (any(outside) .and. available(i1,i2,i3,iv)) then
            !---------------------------------------------------------------------
            ! pick up the vector component
            !---------------------------------------------------------------------
            v(1) = vector(i1,i2,i3,iv,1)
            v(2) = vector(i1,i2,i3,iv,2)
            v(3) = vector(i1,i2,i3,iv,3)
            !---------------------------------------------------------------------
            ! Transform the vector components to the tgt system.  Since e is
            ! essentially the rotation matrix of src (minus that of tgt), the
            ! matmul() operation to go from src to tgt should be matmul(e,v):
            ! The full operation, via Eulerian, would be:
            !   1) get the local src velocity
            !   2) transform components to the Eulerian
            !   3) add the velocity difference btw patches in the E-frame
            !   4) transform to tgt coordinates
            ! But the velocity difference is already expressed in the tgt system,
            ! so this should be fine:
            !---------------------------------------------------------------------
            v = matmul (e, v)
            tm(i1,i2,i3,idx(iv)) = v(iv)
          end if
        end do
      end do
    end do
  end do
END SUBROUTINE deposit_to_target_stag

!===============================================================================
!> Interpolate in 3-D index space and time for a single variable
!===============================================================================
SUBROUTINE from_only (self, tgt, lt, ut, src, only, filled, all_src)
  class(interp_t)             :: self
  class(patch_t)              :: tgt, src
  integer, dimension(3)       :: lt, ut
  integer                     :: only
  integer(2), pointer         :: filled(:,:,:)
  optional                    :: filled
  logical, optional           :: all_src
  !.............................................................................
  real(8), pointer, contiguous:: p(:,:,:,:)
  real, pointer, contiguous   :: tm(:,:,:), s1(:,:,:)
  real                        :: pt(2), level
  real(8)                     :: p1, p2, q1, q2, p3, q3
  real(8)                     :: d(3), e(3,3)
  integer                     :: i1, i2, i3, j1 ,j2, j3, it(2)
  integer, dimension(3)       :: l, u, ls, us
  logical, pointer, contiguous:: available(:,:,:)
  logical                     :: all_src_loc
  integer                     :: m_extrap, j
  integer, save               :: itimer=0, n_extrap=0
  class(lock_t), pointer      :: lock
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%from_only', itimer=itimer)
  all_src_loc = shared%option (all_src)
  !-----------------------------------------------------------------------------
  ! Interpolate in time, in per unit mass variables, leaving source values in
  ! sp(:,:,:,:)
  !-----------------------------------------------------------------------------
  call src%time_interval (tgt, tgt%time, it, pt)
  tm => tgt%mem(:,:,:,only,tgt%it,1)
  !-----------------------------------------------------------------------------
  ! Interpolation coordinates and weights p(:,:,:,1:3)
  !-----------------------------------------------------------------------------
  call os%allocate (p, [src%gn,3], 'p')
  call os%allocate (available, src%gn, 'avail')
  call tgt%get_position (src, p=d, e=e)
  call self%indices (tgt, lt, ut, src, p, available, d, e)
  call source_indices (lt, ut, src, p, available, all_src_loc, ls, us)
  !-----------------------------------------------------------------------------
  ! Tri-linear interpolation + time
  !-----------------------------------------------------------------------------
  if (shared%option(all_src)) then
    l = src%mesh%lb
    u = src%mesh%ub
  else
    l = src%mesh%li
    u = src%mesh%ui
  end if
  tm = 0.0
  do j=1,2
    s1 => src%mem(:,:,:,only,it(j),1)
    lock => src%mem_lock(it(j))
    call                                             lock%set_read ('from_only')
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          !---------------------------------------------------------------------
          ! The available mask allows only points at most 0.5 cell outside ROA
          !---------------------------------------------------------------------
          if (.not.available(i1,i2,i3)) &
            cycle
          j1 = p(i1,i2,i3,1)
          j2 = p(i1,i2,i3,2)
          j3 = p(i1,i2,i3,3)
          j1 = min(max(j1,l(1)),u(1)-1)
          j2 = min(max(j2,l(2)),u(2)-1)
          j3 = min(max(j3,l(3)),u(3)-1)
          p1 = p(i1,i2,i3,1) - j1
          p2 = p(i1,i2,i3,2) - j2
          p3 = p(i1,i2,i3,3) - j3
          q1 = 1. - p1
          q2 = 1. - p2
          q3 = 1. - p3
          level = min(p1,p2,p3,q1,q2,q3,0.0_8)
          !---------------------------------------------------------------------
          ! Among remaining candidate source points, keep the one (or the last
          ! one), that has either the smallest or no extrapolation
          !---------------------------------------------------------------------
          level = 10*(level + src%level)
          if (filled(i1,i2,i3) < int(level)) then
            if (j == 2) &
              filled(i1,i2,i3) = int(level)
            tm(i1,i2,i3) = tm(i1,i2,i3) + pt(j)*( &
              q3*(q2*(q1*s1(j1,j2  ,j3  ) + p1*s1(j1+1,j2  ,j3  )) + &
                  p2*(q1*s1(j1,j2+1,j3  ) + p1*s1(j1+1,j2+1,j3  ))) + &
              p3*(q2*(q1*s1(j1,j2  ,j3+1) + p1*s1(j1+1,j2  ,j3+1)) + &
                  p2*(q1*s1(j1,j2+1,j3+1) + p1*s1(j1+1,j2+1,j3+1)))) 
            !.....................................................................
            if (verbose > 0 .and. filled(i1,i2,i3) > 0.0) &
              write (stdlog,*) i1,i2,i3,filled(i1,i2,i3),level,'replace'
          end if
        end do
      end do
    end do
    call                                           lock%unset_read ('from_only')
  end do
  call os%deallocate (p, 'p')
  call os%deallocate (available, 'avail')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE from_only

!===============================================================================
!> Interpolate in time and convert to velocity and internal energy per mass
!===============================================================================
SUBROUTINE to_per_unit_mass (src, ls, us, sp, it, pt)
  class(patch_t):: src
  integer:: ls(3), us(3)
  real, contiguous:: sp(:,:,:,:)
  integer:: it(2), l(3), u(3)
  real:: pt(2)
  !.............................................................................
  class(lock_t), pointer:: lock
  real, contiguous, pointer:: mem(:,:,:,:)
  integer:: i1, i2, i3, idx_es, j
  real:: di, ptj
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%to_per_unit_mass', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Point pointer, initialize
  !-----------------------------------------------------------------------------
  idx_es = max(src%idx%e,src%idx%s)  
  sp = 0.0
  !-----------------------------------------------------------------------------
  ! Index range
  !-----------------------------------------------------------------------------
  ! l = src%mesh%lo
  ! u = src%mesh%uo
  l = ls
  u = us
  do j=1,2
    lock => src%mem_lock(it(j))
    call                                      lock%set_read ('to_per_unit_mass')
    !---------------------------------------------------------------------------
    ! Density
    !---------------------------------------------------------------------------
    mem => src%mem(:,:,:,:,it(j),1)
    ptj = pt(j)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd
        do i1=l(1),u(1)
          sp(i1,i2,i3,src%idx%d) = &
          sp(i1,i2,i3,src%idx%d) + ptj*mem(i1,i2,i3,src%idx%d)
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    ! Next -  conversion from momenta to velocities. These can be either
    ! cell-centered or staggered. If the quantity is staggered, it is most likely
    ! stagger solver, and since kinetic energy is not needed, we, for simplicity,
    ! use momenta.
    !---------------------------------------------------------------------------
    if (src%mesh(1)%h(src%idx%px) /= 0.0) then
      do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd
        do i1=l(1),u(1)
          sp(i1,i2,i3,src%idx%px) = &
          sp(i1,i2,i3,src%idx%px) + ptj*mem(i1,i2,i3,src%idx%px)
          sp(i1,i2,i3,src%idx%py) = &
          sp(i1,i2,i3,src%idx%py) + ptj*mem(i1,i2,i3,src%idx%py)
          sp(i1,i2,i3,src%idx%pz) = &
          sp(i1,i2,i3,src%idx%pz) + ptj*mem(i1,i2,i3,src%idx%pz)
        end do
      end do
      end do
    else
      do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd private(di)
        do i1=l(1),u(1)
          di = 1.0/mem(i1,i2,i3,src%idx%d)
          sp(i1,i2,i3,src%idx%px) = &
          sp(i1,i2,i3,src%idx%px) + ptj*mem(i1,i2,i3,src%idx%px)*di
          sp(i1,i2,i3,src%idx%py) = &
          sp(i1,i2,i3,src%idx%py) + ptj*mem(i1,i2,i3,src%idx%py)*di
          sp(i1,i2,i3,src%idx%pz) = &
          sp(i1,i2,i3,src%idx%pz) + ptj*mem(i1,i2,i3,src%idx%pz)*di
        end do
      end do
      end do
    end if
    !---------------------------------------------------------------------------
    ! magnetic fields, if applicable
    !---------------------------------------------------------------------------
    if (src%idx%bx > 0) then
      do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd
        do i1=l(1),u(1)
          sp(i1,i2,i3,src%idx%bx) = &
          sp(i1,i2,i3,src%idx%bx) + ptj*mem(i1,i2,i3,src%idx%bx)
          sp(i1,i2,i3,src%idx%by) = &
          sp(i1,i2,i3,src%idx%by) + ptj*mem(i1,i2,i3,src%idx%by)
          sp(i1,i2,i3,src%idx%bz) = &
          sp(i1,i2,i3,src%idx%bz) + ptj*mem(i1,i2,i3,src%idx%bz)
        end do
      end do
      end do
    end if
    !---------------------------------------------------------------------------
    ! Energy of entropy. Quantities with capital letters are already per unit
    ! mass, so just interpolation in time is needed.
    !---------------------------------------------------------------------------
    if (src%etype == 'Entropy') then
      do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd
        do i1=l(1),u(1)
          sp(i1,i2,i3,idx_es) = &
          sp(i1,i2,i3,idx_es) + ptj*mem(i1,i2,i3,idx_es)
        end do
      end do
      end do
    else
      if (src%etype == 'total' .and. verbose > 0 .and. src%idx%bx>0) &
        print *, 'interp_mod: using magnetic fields with total energy in tilted setup might give errors!'
      do i3=l(3),u(3)
      do i2=l(2),u(2)
        !$omp simd
        do i1=l(1),u(1)
          sp(i1,i2,i3,idx_es) = &
          sp(i1,i2,i3,idx_es) + ptj*mem(i1,i2,i3,   idx_es) &
                                   /mem(i1,i2,i3,src%idx%d)
        end do
      end do
      end do
    end if
    call                                    lock%unset_read ('to_per_unit_mass')
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE to_per_unit_mass

!===============================================================================
!> Vectorized transformation from target index space to source index space.
!> For cell-centered quantities.
!> The r(m,m,m,1:3) array holds the floating index distances from the patch
!> position, in units of cells.  The patch%erot(j,k) array holds the absolute
!> Cartesian displacement in the j-direction for one cell along the local
!> unit vector in the k-direction.
!===============================================================================
SUBROUTINE indices1 (tgt, lt, ut, src, ps, available, d, e)
  class(patch_t)       :: tgt, src
  integer, intent(in)  :: lt(3), ut(3)
  real(8), contiguous  :: ps(:,:,:,:)
  logical, pointer     :: available(:,:,:)
  real(8)              :: d(3), e(3,3)
  !.............................................................................
  real(8)                      :: c(3), o(3), v(3)
  integer                      :: i1, i2, i3, i, j, k, ls(3), us(3)
  real, contiguous,    pointer :: r(:,:,:,:)
  integer,                save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%indices1', itimer=itimer)

  !-----------------------------------------------------------------------------
  ! Cache the index space cell coordinates -- common to all patches
  !-----------------------------------------------------------------------------
  call os%allocate (r, [tgt%gn,3], 'r')
  do i3=1,tgt%gn(3)
    do i2=1,tgt%gn(2)
      do i1=1,tgt%gn(1)
        r(i1,i2,i3,1) = tgt%mesh(1)%r(i1)
        r(i1,i2,i3,2) = tgt%mesh(2)%r(i2)
        r(i1,i2,i3,3) = tgt%mesh(3)%r(i3)
      end do
    end do
  end do  
  !-----------------------------------------------------------------------------
  ! To get the coordinates in the src system, we need to take the scalar product
  ! of the location in the tgt system with the unit vectors of the src system,
  ! expressed in tgt coordinates.  These are available, according to the above,
  ! as the columns of the e-matrix.  The coordinates that should go into the
  ! operation are the tgt cordinates relative to the src position; i.e., the
  ! distance from tgt to src, expressed in tgt coordinates, should be subtracted
  ! before the multiplication (if is done after it, the d vector would first
  ! have to be  transformded to src coordinates).
  !-----------------------------------------------------------------------------
  c = 1d0/src%mesh%d
  o = src%mesh%o
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        ps(i1,i2,i3,:) = (matmul(r(i1,i2,i3,:)-d,e))*c + o
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
        available(i1,i2,i3) = (all(int(ps(i1,i2,i3,:)-eps) <  us).and. &
                               all(int(ps(i1,i2,i3,:)+eps) >= ls))
      end do
    end do
  end do
  call os%deallocate (r, 'r')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE indices1


!===============================================================================
!> Vectorized transformation from target index space to source index space.
!> For cell-centered quantities.
!> The r(m,m,m,1:3) array holds the floating index distances from the patch
!> position, in units of cells.  The patch%erot(j,k) array holds the absolute
!> Cartesian displacement in the j-direction for one cell along the local
!> unit vector in the k-direction.
!===============================================================================
SUBROUTINE indices2 (tgt, lt, ut, src, ps, available, d, e)
  class(patch_t)       :: tgt, src
  integer, intent(in)  :: lt(3), ut(3)
  real(8), contiguous  :: ps(:,:,:,:)
  logical, pointer     :: available(:,:,:)
  real(8)              :: d(3), e(3,3)
  !.............................................................................
  real(8)                      :: c(3), o(3), rt(3), rs(3)
  integer                      :: i1, i2, i3, l(3), u(3), ls(3), us(3)
  integer,                save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%indices2', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! To get the coordinates in the src system, we need to take the scalar product
  ! of the location in the tgt system with the unit vectors of the src system,
  ! expressed in tgt coordinates.  These are the rows of the e-matrix.
  !-----------------------------------------------------------------------------
  c = 1d0/src%mesh%d
  o = src%mesh%o
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      !$omp simd
      do i1=lt(1),ut(1)
        rt(1) = tgt%mesh(1)%r(i1) + tgt%centre_nat(1)
        rt(2) = tgt%mesh(2)%r(i2)
        rt(3) = tgt%mesh(3)%r(i3) + tgt%centre_nat(3)
        rs(1) = rt(1)*e(1,1) + rt(2)*e(2,1) + rt(3)*e(3,1) - src%centre_nat(1)
        rs(2) = rt(1)*e(1,2) + rt(2)*e(2,2) + rt(3)*e(3,2)
        rs(3) = rt(1)*e(1,3) + rt(2)*e(2,3) + rt(3)*e(3,3) - src%centre_nat(3)
        ps(i1,i2,i3,1) = c(1)*rs(1) + o(1)
        ps(i1,i2,i3,2) = c(2)*rs(2) + o(2)
        ps(i1,i2,i3,3) = c(3)*rs(3) + o(3)
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
        available(i1,i2,i3) = (all(int(ps(i1,i2,i3,:)-eps) <  us).and. &
                               all(int(ps(i1,i2,i3,:)+eps) >= ls))
      end do
    end do
  end do
  ! call os%deallocate (r, 'r')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE indices2

!===============================================================================
!> Determine the smallest source index cube that covers all points needed in
!> interpolations, allowing a significant reduction in time-interpolation cost
!===============================================================================
SUBROUTINE source_indices (lt, ut, src, p, available, all_src, ls, us)
  class(patch_t):: src
  real(8), contiguous:: p(:,:,:,:)
  logical, contiguous:: available(:,:,:)
  logical:: all_src
  integer:: lt(3), ut(3), ls(3), us(3)
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, j1, j2, j3
  !-----------------------------------------------------------------------------
  ! Determine source index limits
  !-----------------------------------------------------------------------------
  if (all_src) then
    l = src%mesh%lb
    u = src%mesh%ub-1
  else
    l = src%mesh%lo
    u = src%mesh%uo-1
  end if
  ls = u
  us = l
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        if (available(i1,i2,i3)) then
          j1 = p(i1,i2,i3,1)
          j2 = p(i1,i2,i3,2)
          j3 = p(i1,i2,i3,3)
          ls(1) = max (l(1), min (ls(1), j1))
          ls(2) = max (l(2), min (ls(2), j2))
          ls(3) = max (l(3), min (ls(3), j3))
          us(1) = min (u(1), max (us(1), j1))
          us(2) = min (u(2), max (us(2), j2))
          us(3) = min (u(3), max (us(3), j3))
        end if
      end do
    end do
  end do
  us = us + 1
END SUBROUTINE source_indices

!===============================================================================
!> Potentially expand the source index cube for staggered quantities 
!===============================================================================
SUBROUTINE source_indices_stag (lt, ut, src, p, available, all_src, ls, us)
  class(patch_t):: src
  real(8), contiguous:: p(:,:,:,:,:)
  logical, contiguous:: available(:,:,:,:)
  logical:: all_src
  integer:: lt(3), ut(3), ls(3), us(3)
  !.............................................................................
  integer:: l(3), u(3), i1, i2, i3, i4, j1, j2, j3
  !-----------------------------------------------------------------------------
  us = us - 1 ! undo the top limit
  ! ----------------------------------------------------------------------------
  ! set the outermost limits
  ! ----------------------------------------------------------------------------
  if (all_src) then
    l = src%mesh%lb
    u = src%mesh%ub-1
  else
    l = src%mesh%lo
    u = src%mesh%uo-1
  end if
  ! ----------------------------------------------------------------------------
  ! Determine source index limits
  !-----------------------------------------------------------------------------
  do i4=1,3
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          if (available(i1,i2,i3,i4)) then
            j1 = p(i1,i2,i3,i4,1)
            j2 = p(i1,i2,i3,i4,2)
            j3 = p(i1,i2,i3,i4,3)
            ls(1) = max (l(1), min (ls(1), j1-1))
            ls(2) = max (l(2), min (ls(2), j2-1))
            ls(3) = max (l(3), min (ls(3), j3-1))
            us(1) = min (u(1), max (us(1), j1+1))
            us(2) = min (u(2), max (us(2), j2+1))
            us(3) = min (u(3), max (us(3), j3+1))
          end if
        end do
      end do
    end do
  end do
  us = us + 1
END SUBROUTINE source_indices_stag

!===============================================================================
!> Vectorized transformation from target index space to source index space for
!> staggered quantities.
!> The r(m,m,m,1:3) array holds the floating index distances from the patch
!> position, in units of cells.  The patch%erot(j,k) array holds the absolute
!> Cartesian displacement in the j-direction for one cell along the local
!> unit vector in the k-direction.
!> the ps(m,m,m,3,3) holds the indices for all the staggered quantities and
!> available(m,m,m,3) holds corresponding availability.
!===============================================================================
SUBROUTINE indices_stag (self, tgt, lt, ut, src, ps, available, d, e)
  class(interp_t)    :: self
  class(patch_t)     :: tgt, src
  integer            :: lt(3), ut(3)
  real(8), contiguous:: ps(:,:,:,:,:)
  logical, pointer   :: available(:,:,:,:)
  real(8)            :: d(3), e(3,3)
  !.............................................................................
  real(8)                      :: c(3), o(3), v(3)
  integer                      :: i1, i2, i3, i, j, k, ls(3), us(3), iv
  real, contiguous,    pointer :: r(:,:,:,:,:)
  real(8), contiguous, pointer :: dst(:,:), dss(:,:)
  integer,                save :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%indices_stag', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Staggering offsets
  ! vector component | staggering
  !-----------------------------------------------------------------------------
  call os%allocate(dst, [3,3], 'dst', zero=.true.)
  call os%allocate(dss, [3,3], 'dss', zero=.true.)
  ! ----------------------------------------------------------------------------
  ! Find any staggered index. We assume for now, that if there are more than one
  ! staggered quantity, the staggering is the same for all of them.
  ! ----------------------------------------------------------------------------
  if (src%mesh(1)%h(src%idx%px)/=0.0) then
    iv=src%idx%px
  else if (src%idx%bx>0) then
    iv=src%idx%bx
  else
    call io%abort('interp_mod%indices_stag:: staggered index not found.')
  end if
  ! ----------------------------------------------------------------------------
  ! Collect staggering offsets
  ! ----------------------------------------------------------------------------
  do j=1,3
    do i=1,3
      dst(j,i) = tgt%mesh(i)%h(iv+j-1)*tgt%mesh(i)%d
      dss(j,i) = src%mesh(i)%h(iv+j-1)
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Cache the index space cell coordinates -- common to all patches
  !-----------------------------------------------------------------------------
  call os%allocate (r, [tgt%gn,3,3], 'r')
  do iv=1,3
    do i3=1,tgt%gn(3)
      do i2=1,tgt%gn(2)
        do i1=1,tgt%gn(1)
          r(i1,i2,i3,iv,1) = tgt%mesh(1)%r(i1) + dst(iv,1)
          r(i1,i2,i3,iv,2) = tgt%mesh(2)%r(i2) + dst(iv,2)
          r(i1,i2,i3,iv,3) = tgt%mesh(3)%r(i3) + dst(iv,3)
        end do
      end do
    end do
  end do 
  !-----------------------------------------------------------------------------
  ! To get the coordinates in the src system, we need to take the scalar product
  ! of the location in the tgt system with the unit vectors of the src system,
  ! expressed in tgt coordinates.  These are available, according to the above,
  ! as the columns of the e-matrix.  The coordinates that should go into the
  ! operation are the tgt cordinates relative to the src position; i.e., the
  ! distance from tgt to src, expressed in tgt coordinates, should be subtracted
  ! before the multiplication (if is done after it, the d vector would first
  ! have to be  transformded to src coordinates).
  !-----------------------------------------------------------------------------
  c = 1d0/src%mesh%d
  o = src%mesh%o
  do iv=1,3
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          ps(i1,i2,i3,iv,:) = (matmul(r(i1,i2,i3,iv,:)-d,e))*c + o - dss(iv,:)
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! available(:,:,:) is an array, indicating if the ps() is in the valid range
  ! inside the source li and ui
  !-----------------------------------------------------------------------------
  ls = src%mesh%li
  us = src%mesh%ui
  do iv=1,3
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          available(i1,i2,i3,iv) = (all(int(ps(i1,i2,i3,iv,:)-eps) <  us).and. &
                                    all(int(ps(i1,i2,i3,iv,:)+eps) >= ls))
        end do
      end do
    end do
  end do
  call os%deallocate(r, 'r')
  call os%deallocate(dst, 'dst')
  call os%deallocate(dss, 'dss')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE indices_stag

!===============================================================================
!> Interpolate in time over the relevant source region
!===============================================================================
SUBROUTINE in_t (target, source, ls, us, in, out, lt, ut, all, only, extrap, &
                 nbstep, lock)
  class(patch_t)                  :: source, target
  integer, dimension(3)           :: ls, us
  real                            :: in(:,:,:,:,:,:)
  real                            :: out(:,:,:,:)
!  intent(in)                      :: in
!  contiguous                      :: in, out
  integer, dimension(3), optional :: lt, ut
  logical, optional               :: all, extrap
  integer, optional               :: only, nbstep
  type(lock_t), optional          :: lock
  !.............................................................................
  real, pointer, contiguous       :: in1(:,:,:,:), in2(:,:,:,:)
  integer, dimension(3)           :: l, u, o
  integer                         :: i1, i2, i3, i4, jt(2), nv, n, it, id, m(4)
  integer                         :: unit, iv0, iv1, nonpos(2)
  real                            :: pt(2), v_min, d_min
  real, pointer, contiguous       :: d(:,:,:)
  real(8)                         :: wc
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: slock
  integer, save                   :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%in_t', itimer=itimer)
  unit = merge (stdout, stdlog, verbose > 2)
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
  call source%time_interval (target, target%time, jt, pt, extrap=extrap, &
    nbstep=nbstep)
  !-----------------------------------------------------------------------------
  ! Lock access to source while interpolating
  !-----------------------------------------------------------------------------
  m = shape(out)
  call os%allocate (in1, m, 'in1')
  call os%allocate (in2, m, 'in2')
  !
  slock => source%mem_lock(jt(1))
  call                                   slock%set_read ('interp_t%in_t src(1)')
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = target%nv
  end if
  do i4=iv0,iv1                                         ! loop over variables
  do i3=ls(3),us(3)
  do i2=ls(2),us(2)
  do i1=ls(1),us(1)
    in1(i1,i2,i3,i4) = in(i1,i2,i3,i4,jt(1),1)
  end do
  end do
  end do
  end do
  call                                 slock%unset_read ('interp_t%in_t src(1)')
  slock => source%mem_lock(jt(2))
  call                                 slock%set_read   ('interp_t%in_t src(2)')
  do i4=iv0,iv1                                         ! loop over variables
  do i3=ls(3),us(3)
  do i2=ls(2),us(2)
  do i1=ls(1),us(1)
    in2(i1,i2,i3,i4) = in(i1,i2,i3,i4,jt(2),1)
  end do
  end do
  end do
  end do
  call                                 slock%unset_read ('interp_t%in_t src(2)')
  !.............................................................................
  if (target%id == io%id_debug) then
    write (stdlog,'("in_t:",2(i6,":",i2.2,i3),f10.6,2(i4,f10.6,f7.3))') &
      target%id, target%level, target%istep, &
      source%id, source%istep, source%level, &
      target%time, &
      jt(1), source%t(jt(1)), pt(1), &
      jt(2), source%t(jt(2)), pt(2)
  end if
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(a,6i3,3x,4i3)') &
    'interp_t%in_t: shapes=', shape(in), shape(out)
  end if
  !-----------------------------------------------------------------------------
  ! Shifted from source index space
  !-----------------------------------------------------------------------------
  if (present(lock))                                      call lock%set ('in_t')
  if (present(lt)) then
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,9l1)') &
        wallclock(), pad_trace(1:indent)// &
        'interp_t%in_t: tgt,src,jt,pt,ls,us,lt,ut,usg=', &
        target%id, source%id, jt, pt, ls, us, lt, ut, &
        source%unsigned
      if (unit == stdout) then
        write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,2(1x,3i3)),1x,9l1)') &
        'interp_t%in_t(1): tgt,src,jt,pt,ls,us,lt,ut,usg=', &
        target%id, target%level, target%time, &
        source%id, source%level, source%time, &
        jt, pt, ls, us, lt, ut, &
        source%unsigned
      end if
    end if
    !...........................................................................
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
    !...........................................................................
    id = target%idx%d
    if (any(target%pervolume)) then
      call os%allocate (d, target%gn, 'd')
    end if
    nonpos = 0
    o = lt-ls
    do i4=iv0,iv1                                         ! loop over variables
      do i3=ls(3),us(3)
        !-----------------------------------------------------------------------
        ! Use log interpolation only when needed
        !-----------------------------------------------------------------------
        if (target%unsigned(i4)) then
          do i2=ls(2),us(2)
            if (verbose > 0) then
              do i1=ls(1),us(1)
                if (in1(i1,i2,i3,i4) <= 0.0) then
                  if (verbose > 1) then
                    write (stdlog,1) &
                      'in_t: non-positive for in1 id,ii=', &
                      target%time, source%t(jt(1)), &
                      source%id, jt(1), i1, i2, i3, i4, &
                      in1(i1,i2,i3,i4)
                    1 format(a,2f10.6,i7,5i3,1p,g11.2)
                  end if
                  nonpos(1) = nonpos(1) + 1
                  in1(i1,i2,i3,i4) = in2(i1,i2,i3,i4)
                end if
                if (in2(i1,i2,i3,i4) <= 0.0) then
                  if (verbose > 1) then
                    write (stdlog,1) &
                      'in_t: non-positive for in2 id,ii=', &
                      target%time, source%t(jt(2)), &
                      source%id, jt(2), i1, i2, i3, i4, &
                      in2(i1,i2,i3,i4)
                  end if
                  nonpos(2) = nonpos(2) + 1
                  in2(i1,i2,i3,i4) = in1(i1,i2,i3,i4)
                end if
                out(i1+o(1),i2+o(2),i3+o(3),i4) = &
                exp(log(in1(i1,i2,i3,i4))*pt(1) &
                  + log(in2(i1,i2,i3,i4))*pt(2))
              end do
              flush (stdlog)
            else
              !$omp simd
              do i1=ls(1),us(1)
                out(i1+o(1),i2+o(2),i3+o(3),i4) = &
                exp(log(in1(i1,i2,i3,i4))*pt(1) &
                  + log(in2(i1,i2,i3,i4))*pt(2))
              end do
            end if
          end do
          !---------------------------------------------------------------------
          ! If per-volume variables are used, save the density
          !---------------------------------------------------------------------
          if (i4==id .and. any(target%pervolume)) then
            do i2=ls(2),us(2)
              !$omp simd
              do i1=ls(1),us(1)
                  d(i1+o(1),i2+o(2),i3+o(3)) = &
                out(i1+o(1),i2+o(2),i3+o(3),i4)
              end do
            end do
          end if
        !-----------------------------------------------------------------------
        ! Use per-volume interpolation
        !-----------------------------------------------------------------------
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
        !-----------------------------------------------------------------------
        ! Signed variables, or safe positive, shifted in index
        !-----------------------------------------------------------------------
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
    if (any(nonpos > 0)) then
      write (stdlog,2) &
        'in_t: non-positive id,stp,jt,n,times,iit=', target%id, &
        source%id, source%istep, jt, nonpos, target%time, source%t(jt), &
        source%iit
      write (stdout,2) &
        'in_t: non-positive id,stp,jt,n,times,iit=', target%id, &
        source%id, source%istep, jt, nonpos, target%time, source%t(jt), &
        source%iit
      2 format(a,2i7,3i4,2i6,3f11.7,1x,9i3)
    end if
    if (any(target%pervolume)) then
      call os%deallocate (d, 'd')
    end if
    n = product(us-ls+1)
  !-----------------------------------------------------------------------------
  ! Unshifted interpolation in time
  !-----------------------------------------------------------------------------
  else
    o = 0
    !...........................................................................
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2i3,2f7.3,2(1x,3i3),1x,9l1)') &
        wallclock(), &
        pad_trace(1:indent)//'interp_t%in_t: tgt,source,jt,pt,ls,us,usg=', &
        target%id, source%id, jt, pt, ls, us, &
        source%unsigned
      if (unit == stdout) then
        write (unit,'(a,2(i6,":",i2.2,f10.6),1x,2i3,2f7.3,2(1x,3i3),1x,9l1)') &
        'interp_t%in_t(2): tgt,src,jt,pt,ls,us,lt,ut,usg=', &
        target%id, target%level, target%time, &
        source%id, source%level, source%time, &
        jt, pt, ls, us, &
        source%unsigned
      end if
    end if
    !...........................................................................
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
  if (present(lock))                                    call lock%unset ('in_t')
  !.............................................................................
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
  call os%deallocate (in1, 'in1')
  call os%deallocate (in2, 'in2')
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
SUBROUTINE in_x (tgt, lt, ut, src, ls, us, om, in, out, all, only, lock)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  type(lock_t), optional   :: lock
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), nf(2), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: x, wc, wtimer, sm_p, d(3)
  real(8), pointer, contiguous:: dst(:), dss(:)
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
  call os%allocate(dst, tgt%mv, 'dst', zero=.true.)
  call os%allocate(dss, src%mv, 'dss', zero=.true.)
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
  if (present(lock))                                      call lock%set ('in_x')
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
  if (present(lock))                                    call lock%unset ('in_x')
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
  call os%deallocate(dst, 'dst')
  call os%deallocate(dss, 'dss')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_x

!===============================================================================
!> Interpolate in x over the relevant target region
!===============================================================================
SUBROUTINE in_y (tgt, lt, ut, src, ls, us, om, in, out, all, only, lock)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  type(lock_t), optional   :: lock
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, sh(4), iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: y, wc, wtimer, sm_p, d(3)
  real(8), pointer, contiguous:: dst(:), dss(:)
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
  call os%allocate(dst, tgt%mv, 'dst', zero=.true.)
  call os%allocate(dss, src%mv, 'dss', zero=.true.)
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
  if (present(lock))                                      call lock%set ('in_y')
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
  if (present(lock))                                    call lock%unset ('in_y')
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
  call os%deallocate(dst, 'dst')
  call os%deallocate(dss, 'dss')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_y

!===============================================================================
!> Interpolate in x over the relevant target region
!===============================================================================
SUBROUTINE in_z (tgt, lt, ut, src, ls, us, om, in, out, all, only, lock)
  class(patch_t)           :: src, tgt
  integer, dimension(3)    :: lt, ut, ls, us, l1, u1
  type(meshp_t)            :: om(3)
  real, dimension(:,:,:,:) :: in, out
  pointer                  :: in, out
  logical                  :: all
  integer, optional        :: only
  type(lock_t), optional   :: lock
  !.............................................................................
  integer                  :: j, i1, i2, i3, i4, nv, n, iv0, iv1
  integer, dimension(3)    :: l, u, si, st, li, ui
  real                     :: p, q, v_min, d_min
  real(8)                  :: z, wc, wtimer, sm_p, d(3)
  real(8), pointer, contiguous:: dst(:), dss(:)
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
  call os%allocate(dst, tgt%mv, 'dst', zero=.true.)
  call os%allocate(dss, src%mv, 'dss', zero=.true.)
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
  if (present(lock))                                      call lock%set ('in_z')
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
  if (present(lock))                                    call lock%unset ('in_z')
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
  call os%deallocate(dst, 'dst')
  call os%deallocate(dss, 'dss')
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE in_z

!===============================================================================
!> Not interpolating in time, but simply copying over the relevant source region
!===============================================================================
SUBROUTINE no_t (target, source, ls, us, in, out, lt, ut, only, lock)
  class(patch_t)                  :: source, target
  integer, dimension(3)           :: ls, us
  real                            :: in (:,:,:,:,:,:)
  real                            :: out(:,:,:,:,:,:)
!  intent(in)                      :: in
!  contiguous                      :: in, out
  integer, dimension(3), optional :: lt, ut
  integer, optional               :: only
  type(lock_t), optional          :: lock
  !.............................................................................
  integer, dimension(3)           :: l, u, o
  integer                         :: i1, i2, i3, i4, i5, n, it, id, m(5)
  integer                         :: unit, iv0, iv1
  real(8)                         :: wc
  real(8), save                   :: wtimer
  type(lock_t), pointer           :: slock
  integer, save                   :: itimer=0, ntimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('interp_t%no_t', itimer=itimer)
  unit = merge (stdout, stdlog, verbose > 2)
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
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = target%nv
  end if
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
    write (stdlog,'(a,6i5,3x,6i5)') &
    'interp_t%no_t: shapes=', shape(in), shape(out)
  end if
  !---------------------------------------------------------------------------
  ! Shifted from source index space
  !---------------------------------------------------------------------------
  if (present(lock))                                      call lock%set ('no_t')
  if (present(lt)) then
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,2i6,1x,2(1x,2(1x,3i3)),1x,9l1)') &
        wallclock(), pad_trace(1:indent)//'interp_t%no_t: tgt,src,ls,us,lt,ut,usg=', &
        target%id, source%id, ls, us, lt, ut, source%unsigned
        if (unit == stdout) then
          write (unit,'(a,2(i6,":",i2.2,f10.6),2(1x,2(1x,3i3)),1x,9l1)') &
          'interp_t%no_t(1): tgt,src,ls,us,lt,ut,usg=', &
          target%id, target%level, target%time, &
          source%id, source%level, source%time, &
          ls, us, lt, ut, source%unsigned
        end if
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
  if (present(lock))                                    call lock%unset ('no_t')
  !...........................................................................
  if (verbose > 1 .and. present(lt)) then
    wc = wallclock()-wc
    !lt = ls+o
    write (stdlog,'(a,f6.3,1p,6(2x,3e11.2))') &
    'interp_t%no_t(2): mus/p,min,max,min,max=', wc/n, &
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
