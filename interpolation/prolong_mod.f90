!===============================================================================
!> Specialized module for optimzed, conservative prolongation.  The main
!> interface from the DISPATCH code is
!>
!> USE prolong_mod
!> ...
!>   integer:: ipos(3)
!>   ...
!>   ipos = nint((source%position-target%position)/target%size)
!>   call prolong%interpolate (source%mem, target%mem, ipos)
!>
!> Implementation:
!>
!> There are two possible arrangements of source and target points, relevant
!> for cases with 2 (hllc) and 3 (hlld) guard zones:
!>
!> Example with ng=2, n=4:
!>          1       2   |   3       4       5       6   |   7       8
!>      ...-x---+---x---|---o---+---o---+---o---+---o---|---x---+---x-
!>                x---x-|-o---o---o---o---o---o   o---o-|-x---x-
!>                1   2   3   4   5   6   7   8   1   2   3   4
!>               -+---+---o---o-  o---o---o---o---o---o---x---x-
!>                5   6   7   8   1   2   3   4   5   6   7   8
!>
!> At the lower end target index range is [1,2] and the source range is [1,3]
!> At the upper end target index range is [7,8] and the source range is [6,8]
!>
!> Example with ng=3:
!>  1       2       3   |   4       5       6       7   |   8       9      10
!> -o---+---o---+---o---|---o---+---o---+---o---+---o---|---o---+---o---+---o
!>           -+---+---+-|-+---+--                -+---+-|-+---+---+
!>            1   2   3   4   5                   6   7   8   9  10
!>
!> In the first case only one target point is generated for the first pair of
!> source points, while subsequent source intervals generate two target points.
!> In the second case all source intervals generate two points.  We use
!> mod(ng,2) to separate into case 0 and case 1, to cover any alignment
!===============================================================================
MODULE prolong_mod
  USE omp_lib
  USE io_mod
  USE trace_mod
  USE kinds_mod
  implicit none
  private
  type, public:: prolong_t
    integer           :: case=0                         ! mod(ng,2)
    integer           :: ng=2                           ! guard zones
    integer           :: lb, lo, li, ui, uo, ub         ! convention
    integer           :: n=8                            ! internal
    integer           :: gn=12                          ! total
    real, pointer     :: w1d(:)     => null()           ! 1D weights
    real, pointer     :: w3d(:,:,:,:) => null()         ! 3D weights
  contains
    procedure         :: init
    procedure         :: weights_1d
    procedure         :: weights_3d
    procedure         :: dealloc
    procedure         :: finalize
    procedure         :: info
    procedure         :: x_interp1
    procedure         :: x_interp3
    generic           :: x_interp => x_interp1, x_interp3
    procedure         :: y_interp1
    procedure         :: y_interp3
    generic           :: y_interp => y_interp1, y_interp3
    procedure         :: z_interp1
    procedure         :: z_interp3
    generic           :: z_interp => z_interp1, z_interp3
    procedure         :: interp
    procedure         :: interpolate
    procedure         :: test
    procedure         :: patch_test
    procedure         :: dnload_test
    procedure         :: print
  end type
  logical, save              :: array_notation=.false.
  logical, save              :: fix_nml=.true.
  integer, save              :: verbose=0
  logical, save, allocatable :: outside(:,:,:)
  logical, save              :: first_time=.true.
CONTAINS

!===============================================================================
!> Initalize weight factors, based on number of total and guard cells.  The
!> arrays are assumed to be cubical, for simplicity and efficiency.
!===============================================================================
SUBROUTINE init (self, n, ng)
  class(prolong_t):: self
  integer:: n, ng
  !.............................................................................
  integer:: gn, iostat, i1, i2, i3
  namelist /prolong_params/ verbose, array_notation, fix_nml
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%init')
  if (first_time) then
    call io%lock%set ('prolong_t%init')
    if (first_time) then
      call trace%begin ('prolong_t%init')
      rewind (io%input)
      read (io%input, prolong_params, iostat=iostat)
      call io%namelist_warning ('prolong_params', iostat)
      write (io%output, prolong_params)
      !-------------------------------------------------------------------------
      gn = n + ng*2
      self%n  = n
      self%ng = ng
      self%gn = gn
      self%lb = 1
      self%li = ng+1
      self%lo = ng
      self%ub = gn
      self%ui = gn-ng
      self%uo = gn-ng+1
      !-------------------------------------------------------------------------
      ! mask(:,:,:) is true only in the guard zones
      !-------------------------------------------------------------------------
      allocate (outside(gn,gn,gn))
      outside = .true.
      outside(self%li:self%ui,self%li:self%ui,self%li:self%ui) = .false.
      write (io%output,*) 'inside:', self%li, self%ui, &
        count(outside), count(.not.outside)
      flush (io%output)
      self%case = mod(ng,2)
      call self%weights_1d
      call trace%end ()
      first_time = .false.
    end if
    call io%lock%unset ('prolong_t%init')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initalize weight factors, based on number of total and guard cells.  The
!> arrays are assumed to be cubical, for simplicity and efficiency.
!===============================================================================
SUBROUTINE weights_1d (self)
  class(prolong_t):: self
  integer:: i0, i1, gn
  !-----------------------------------------------------------------------------
  ! Case 0 starts with w=0.75, case 1 with 0.25
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%weights_1d')
  call self%dealloc
  gn = self%gn
  allocate (self%w1d(gn))
  i0 = mod(self%ng,2)
  do i1=1,gn
    self%w1d(i1) = 0.25 + mod(i1+i0,2)*0.50       ! starts with 0.75
  end do
  call trace%end ()
END SUBROUTINE weights_1d

!===============================================================================
!> Initialize weight factors, based on number of total and guard cells.  The
!> arrays are assumed to be cubical, for simplicity and efficiency.
!===============================================================================
SUBROUTINE weights_3d (self)
  class(prolong_t):: self
  integer:: i1, i2, i3, i0, gn
  !-----------------------------------------------------------------------------
  ! Case 0 starts with w=0.75, case 1 with 0.25
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%weights_3d')
  call self%dealloc
  gn = self%gn
  allocate (self%w3d(gn,gn,gn,3))
  i0 = mod(self%ng,2)
  do i3=1,gn
  do i2=1,gn
  do i1=1,gn
    self%w3d(i1,i2,i3,1) = 0.25 + mod(i1+i0,2)*0.50
    self%w3d(i1,i2,i3,2) = 0.25 + mod(i2+i0,2)*0.50
    self%w3d(i1,i2,i3,3) = 0.25 + mod(i3+i0,2)*0.50
  end do
  end do
  end do
  call trace%end ()
END SUBROUTINE weights_3d

!===============================================================================
!> Deallocate the wights array (when the task update is finished)
!===============================================================================
SUBROUTINE dealloc (self)
  class(prolong_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%dealloc')
  if (associated(self%w1d)) then
    deallocate (self%w1d)
    nullify (self%w1d)
  end if
  if (associated(self%w3d)) then
    deallocate (self%w3d)
    nullify (self%w3d)
  end if
  call trace%end ()
END SUBROUTINE

!===============================================================================
!> Deallocate the wights array (when the task update is finished)
!===============================================================================
SUBROUTINE finalize (self)
  class(prolong_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%finalize')
  call self%dealloc
  deallocate (outside)
  first_time = .true.
  call trace%end ()
END SUBROUTINE finalize

!===============================================================================
!===============================================================================
SUBROUTINE info (self)
  class(prolong_t):: self
  write (io_unit%log,'(3(3x,a,":",i4))') 'n', self%n, 'ng', self%ng, 'gn', self%gn
END SUBROUTINE

!===============================================================================
!> Prolongation of cell values for a target guard zone.  The location of the
!> source is given by the integer(3) ipos, where -1 means "left", +1 means
!> "right", and 0 means the same position.
!===============================================================================
SUBROUTINE interpolate (self, source, jt, pt, target, ipos, all, l, u)
  class(prolong_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:) :: source
  real(kind=KindScalarVar), dimension(:,:,:)   :: target
  integer                  :: jt(2), ipos(3)
  real                     :: pt(2)
  intent(in)               :: source, jt, pt, ipos
  intent(out)              :: target
  logical, optional        :: all
  integer, optional        :: l(3), u(3)
  !.............................................................................
  integer                  :: ls(3), us(3), lt(3), ut(3), i
  integer                  :: lb, lo, li, ui, uo, ub, n, ng, gn
  integer                  :: itimer=0
  logical                  :: all_l
  real(kind=KindScalarVar), allocatable        :: src_mem(:,:,:)
  !-----------------------------------------------------------------------------
  ! Cf. the sketch at top for index choices
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%interpolate',itimer=itimer)
  if (verbose > 0) &
    write (io_unit%log,*) 'interpolate', ipos
  if (present(all)) then
    all_l = all
  else
    all_l = .false.
  end if
  n  = self%n
  ng = self%ng
  gn = n + 2*ng
  lb = 1
  li = lb + ng
  lo = li - 1
  ui = lo + n
  uo = ui + 1
  ub = ui + ng
  !-----------------------------------------------------------------------------
  ! Loop over directions.  For each direction the relative position is giveń,
  ! in terms of the displacement of the target position relative to the source,
  ! measured in half target size units (integers -3 -1, +1, +3)
  !-----------------------------------------------------------------------------
  do i=1,3
    select case (ipos(i))
  !-----------------------------------------------------------------------------
  ! Target is to the left of the source
  !-----------------------------------------------------------------------------
!>          1       2   |   3       4
!>      ...-x---+---x---|---o---+---o-
!>               -+---+---o---o-
!>                5   6   7   8
    case (+3)
      if (verbose > 2) &
        write(io_unit%log,*) 'source in direction', i,' is to the right'
      lt(i) = uo
      ut(i) = ub
      ls(i) = lo               ! = 2
      us(i) = lo + ng/2 + 1    ! = 2 + 1 + 1 = 4
  !-----------------------------------------------------------------------------
  ! Target is in the left part of the source
  !-----------------------------------------------------------------------------
!>          1       2   |   3       4       5       6
!>      ...-x---+---x---|---o---+---o---+---o---+---o--
!>                x---x-|-o---o---o---o---o---o
!>                1   2   3   4   5   6   7   8   1   2   3   4
    case (+1)
      if (verbose > 2) &
        write(io_unit%log,*) 'source in direction', i,' is inside right'
      lt(i) = li
      ut(i) = ub
      ls(i) = lo                      ! 2
      us(i) = lo + ng/2 + 1 + n/2     ! 2 + 1 + 1 + 2 = 6
  !-----------------------------------------------------------------------------
  ! Target is in the right part of the source
  !-----------------------------------------------------------------------------
!>          1       2   |   3       4       5       6   |   7       8
!>      ...-x---+---x---|---o---+---o---+---o---+---o---|---x---+---x-
!>                                o---o---o---o---o---o---x---x-
!>                                1   2   3   4   5   6   7   8
    case (-1)
      if (verbose > 2) &
        write(io_unit%log,*) 'source in direction', i,' is inside left'
      lt(i) = lb
      ut(i) = ui
      ls(i) = uo - ng/2 - 1 - n/2   ! 7 - 1 - 1 - 2 = 3
      us(i) = uo                    ! 7
  !-----------------------------------------------------------------------------
  ! Target is to the right of the source
  !-----------------------------------------------------------------------------
!>          1       2   |   3       4       5       6   |   7       8
!>      ...-x---+---x---|---o---+---o---+---o---+---o---|---x---+---x-
!>                      |                         o---o-|-x---x-
!>                                                1   2   3   4
    case (-3)
      if (verbose > 2) &
        write(io_unit%log,*) 'source in direction', i,' is to the left'
      lt(i) = lb
      ut(i) = lo
      ls(i) = ui - 1              ! = 5
      us(i) = ui + ng/2           ! = 6 + 1 = 7
    case default
      write (io_unit%log,'(1x,a,3i3)') &
        'prolong_t%interpolat: unexpcted ipos =', ipos
      call trace%end (itimer)
      return
    end select
  end do
  if (verbose > 2) then
    call self%info
    write (io_unit%log,'(4(3x,a,":",3i4))') 'lt', lt, 'ut', ut, 'ls', ls
  end if
  !-----------------------------------------------------------------------------
  ! Interpolate in time, but only in the source region that is used
  !-----------------------------------------------------------------------------
  allocate (src_mem(gn,gn,gn))
  if (any(ls < 1) .or. any (us > ub)) then
    write(io_unit%log,*) ipos, 'ipos'
    write(io_unit%log,*) ls, 'ls'
    write(io_unit%log,*) us, 'us'
    write(io_unit%log,*) n, ng, 'n, ng'
    write(io_unit%log,*) lb, lo, li, 'lb, lo, li'
    write(io_unit%log,*) ui, uo, ub, 'ui, uo, ub'
  end if
  src_mem   (ls(1):us(1),ls(2):us(2),ls(3):us(3)) &
    = source(ls(1):us(1),ls(2):us(2),ls(3):us(3),jt(1))*pt(1) &
    + source(ls(1):us(1),ls(2):us(2),ls(3):us(3),jt(2))*pt(2)
  !-----------------------------------------------------------------------------
  ! Interpolate in 3-D
  !-----------------------------------------------------------------------------
  call io%assert (.not.any(shape(outside)/=shape(src_mem)), &
    'interpolation%interpolate')
  if (all_l) then
    call self%interp (src_mem, ls, us, target, lt, ut)
  else
    call self%interp (src_mem, ls, us, target, lt, ut, mask=outside)
  end if
  deallocate (src_mem)
  if (present(l)) then
    l = lt
  end if
  if (present(u)) then
    u = ut
  end if
  call trace%end (itimer)
END SUBROUTINE interpolate

!===============================================================================
!> Case with 3D weights
!===============================================================================
SUBROUTINE interp (self, source, ls, us, target, lt, ut, mask)
  class(prolong_t):: self
  real(kind=KindScalarVar), intent(in):: source(:,:,:)
  real(kind=KindScalarVar), intent(out):: target(:,:,:)
  integer, intent(in):: lt(3), ut(3), ls(3), us(3)
  logical, optional:: mask(:,:,:)
  !.............................................................................
  integer:: m(3), l(3), u(3)
  real(kind=KindScalarVar), allocatable:: target1(:,:,:), target2(:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('prolong_t%interp')
  m = shape(target)
  allocate (target1(m(1),m(2),m(3)), target2(m(1),m(2),m(3)))
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=ls(2); u(2)=us(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  if (associated(self%w1d)) then
    call self%x_interp (source , self%w1d         , ls(1), target1, l, u)
  else
    call self%x_interp (source , self%w3d(:,:,:,1), ls(1), target1, l, u)
  end if
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  if (associated(self%w1d)) then
    call self%y_interp (target1, self%w1d         , ls(2), target2, l, u)
  else
    call self%y_interp (target1, self%w3d(:,:,:,2), ls(2), target2, l, u)
  end if
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=lt(3); u(3)=ut(3)        ! source range in z
  if (associated(self%w1d)) then
    call self%z_interp (target2, self%w1d         , ls(3), target , l, u, &
      mask=mask)
  else
    call self%z_interp (target2, self%w3d(:,:,:,3), ls(3), target , l, u, &
      mask=mask)
  end if
  deallocate (target1, target2)
  ! --- NML fix ----------------------------------------------------------------
  if (fix_nml) &
    call nml_fix (self, target, l, u)
  call trace%end ()
END SUBROUTINE interp

!===============================================================================
!> Patch up the no-mans-land points, by interpolating btw the next point --
!> which originates from the proper interior of the nbor -- and the first
!> interior point
!===============================================================================
SUBROUTINE nml_fix (self, out, l, u)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  integer, intent(in)           :: l(3), u(3)
  !.............................................................................
  integer                       :: i1, i2, i3
  !-----------------------------------------------------------------------------
  if      (l(1) == self%lo) then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      out(self%lo,i2,i3) = 0.5*(out(self%li,i2,i3)+out(self%lo-1,i2,i3))
    end do
    end do
  else if (u(1) == self%uo) then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      out(self%uo,i2,i3) = 0.5*(out(self%ui,i2,i3)+out(self%uo+1,i2,i3))
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  if (l(2) == self%lo) then
    do i3=l(3),u(3)
    do i1=l(1),u(1)
      out(i1,self%lo,i3) = 0.5*(out(i1,self%li,i3)+out(i1,self%lo-1,i3))
    end do
    end do
  else if (u(2) == self%uo) then
    do i3=l(3),u(3)
    do i1=l(1),u(1)
      out(i1,self%uo,i3) = 0.5*(out(i1,self%ui,i3)+out(i1,self%uo+1,i3))
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  if (l(3) == self%lo) then
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,self%lo) = 0.5*(out(i1,i2,self%li)+out(i1,i2,self%lo-1))
    end do
    end do
  else if (u(3) == self%uo) then
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,self%uo) = 0.5*(out(i1,i2,self%ui)+out(i1,i2,self%uo+1))
    end do
    end do
  end if
END SUBROUTINE nml_fix

!===============================================================================
!> Interpolate in the x-direction
!===============================================================================
SUBROUTINE x_interp1 (self, in, w, ls, out, l, u)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:)
  integer, intent(in)           :: l(3), u(3), ls
  !.............................................................................
  integer                       :: i1, i2, i3, j1, k1
  !-----------------------------------------------------------------------------
  j1 = ls
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        out(i1,i2,i3) = in(j1  ,i2,i3)*(1.-w(i1))&
                      + in(j1+1,i2,i3)*    w(i1)
      end do
    end do
    k1 = j1+1
    if (mod(i1,2)==1) &
      j1 = j1+1
  end do
  if (verbose > 3) &
    write (io_unit%log,'(" X1",5x,4(a,":",2i4,5x),1p,4e11.2)') &
    'xs', ls,k1, 'x', l(1),u(1), 'y', l(2),u(2), 'z',l(3),u(3), &
    minval(in (ls  :k1  ,l(2):u(2),l(3):u(3))), &
    maxval(in (ls  :k1  ,l(2):u(2),l(3):u(3))), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE x_interp1
!===============================================================================
SUBROUTINE x_interp3 (self, in, w, ls, out, l, u)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:,:,:)
  integer, intent(in)           :: l(3), u(3), ls
  !.............................................................................
  integer                       :: i1, i2, i3, j1, k1
  !----------------------------------------------------------------------------.
  j1 = ls
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        out(i1,i2,i3) = in(j1  ,i2,i3)*(1.-w(i1,i2,i3))&
                      + in(j1+1,i2,i3)*    w(i1,i2,i3)
      end do
    end do
    if (verbose > 5) then
      i2=l(2); i3=l(3)
      write (io_unit%log,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j1,in(j1  ,i2,i3),in(j1+1,i2,i3),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k1 = j1+1
    if (mod(i1,2)==1) &
      j1 = j1+1
  end do
  if (verbose > 3) &
    write (io_unit%log,'(" X3",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'xs', ls,k1, 'x ', l(1),u(1), 'y ', l(2),u(2), 'z ',l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE x_interp3

!===============================================================================
!> Interpolate in the y-direction
!===============================================================================
SUBROUTINE y_interp1 (self, in, w, ls, out, l, u)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:)
  integer, intent(in)           :: l(3), u(3), ls
  !.............................................................................
  integer                       :: i1, i2, i3, j2, k2
  !-----------------------------------------------------------------------------
  j2 = ls
  do i2=l(2),u(2)
    do i3=l(3),u(3)
      do i1=l(1),u(1)
        out(i1,i2,i3) = in(i1,j2  ,i3)*(1.-w(i2))&
                      + in(i1,j2+1,i3)*    w(i2)
      end do
    end do
    k2 = j2+1
    if (mod(i2,2)==1) &
      j2 = j2+1
  end do
  if (verbose > 3) &
    write (io_unit%log,'(" Y1",5x,4(a,":",2i4,5x),1p,4e11.2)') &
    'x', l(1),u(1), 'ys', ls,k2, 'y', l(2),u(2), 'z',l(3),u(3), &
    minval(in (l(1):u(1),ls  :k2  ,l(3):u(3))), &
    maxval(in (l(1):u(1),ls  :k2  ,l(3):u(3))), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE y_interp1
!===============================================================================
!===============================================================================
SUBROUTINE y_interp3 (self, in, w, ls, out, l, u)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:,:,:)
  integer, intent(in)           :: l(3), u(3), ls
  !.............................................................................
  integer                       :: i1, i2, i3, j2, k2
  !-----------------------------------------------------------------------------
  j2 = ls
  do i2=l(2),u(2)
    if (array_notation) then
      out(:,i2,:) = in(:,j2  ,:)*(1.-w(:,i2,:))&
                  + in(:,j2+1,:)*    w(:,i2,:)
    else
      do i3=l(3),u(3)
        do i1=l(1),u(1)
          out(i1,i2,i3) = in(i1,j2  ,i3)*(1.-w(i1,i2,i3))&
                        + in(i1,j2+1,i3)*    w(i1,i2,i3)
        end do
      end do
    end if
    if (verbose > 5) then
      i1=l(1); i3=l(3)
      write (io_unit%log,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j2,in(i1,j2,i3),in(i1,j2+1,i3),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k2 = j2+1
    if (mod(i2,2)==1) &
      j2 = j2+1
  end do
  if (verbose > 3) &
    write (io_unit%log,'(" Y3",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'x ', l(1),u(1), 'ys', ls,k2, 'y ', l(2),u(2), 'z ',l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE y_interp3

!===============================================================================
!> Interpolate in the z-direction
!===============================================================================
SUBROUTINE z_interp1 (self, in, w, ls, out, l, u, mask)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:)
  logical, intent(in), optional :: mask(:,:,:)
  integer, intent(in)           :: l(3), u(3), ls
  !............................................................................
  integer                       :: i1, i2, i3, j3, k3
  !----------------------------------------------------------------------------
  j3 = ls
  do i3=l(3),u(3)
    if (present(mask)) then
      if (array_notation) then
        where(mask(:,:,i3))
          out(:,:,i3) = in(:,:,j3  )*(1.-w(i3)) &
                      + in(:,:,j3+1)*    w(i3)
        end where
      else
        do i2=l(2),u(2)
          do i1=l(1),u(1)
             out(i1,i2,i3) = merge(in(i1,i2,j3  )*(1.-w(i3)) &
                                 + in(i1,i2,j3+1)*    w(i3), &
                                  out(i1,i2,i3), mask(i1,i2,i3))
          end do
        end do
      end if
    else if (array_notation) then
      out(:,:,i3) = in(:,:,j3  )*(1.-w(i3)) &
                  + in(:,:,j3+1)*    w(i3)
    else
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          out(i1,i2,i3) = in(i1,i2,j3  )*(1.-w(i3)) &
                        + in(i1,i2,j3+1)*    w(i3)
        end do
      end do
    end if
    if (verbose > 5) then
      i1=l(1); i2=l(2)
      write (io_unit%log,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j3,in(i1,i2,j3),in(i1,i2,j3+1),w(i3),out(i1,i2,i3)
    end if
    k3 = j3+1
    if (mod(i3,2)==1) &
      j3 = j3+1
  end do
  if (verbose > 3) &
    write (io_unit%log, '(" Z1",5x,4(a,":",2i4,4x),1p,4e11.2)') &
    ' x', l(1),u(1), ' y', l(2),u(2), 'zs',ls,k3, ' z', l(3),u(3), &
    minval(in (l(1):u(1),l(2):u(2),ls  :k3  )), &
    maxval(in (l(1):u(1),l(2):u(2),ls  :k3  )), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE z_interp1
!===============================================================================
!===============================================================================
SUBROUTINE z_interp3 (self, in, w, ls, out, l, u, mask)
  class(prolong_t)              :: self
  real(kind=KindScalarVar), contiguous, intent(in)  :: in(:,:,:)
  real(kind=KindScalarVar), contiguous, intent(out) :: out(:,:,:)
  real, contiguous, intent(in)  :: w(:,:,:)
  logical, intent(in), optional :: mask(:,:,:)
  integer, intent(in)           :: l(3), u(3), ls
  !.............................................................................
  integer                       :: i1, i2, i3, j3, k3
  real:: tmp
  !-----------------------------------------------------------------------------
  j3 = ls
  do i3=l(3),u(3)
    if (present(mask)) then
      if (array_notation) then
        where(mask(:,:,i3))
          out(:,:,i3) = in(:,:,j3  )*(1.-w(:,:,i3)) &
                      + in(:,:,j3+1)*    w(:,:,i3)
        end where
      else
        do i2=l(2),u(2)
          do i1=l(1),u(1)
             out(i1,i2,i3) = merge(in(i1,i2,j3  )*(1.-w(i1,i2,i3)) &
                                 + in(i1,i2,j3+1)*    w(i1,i2,i3), &
                                  out(i1,i2,i3), mask(i1,i2,i3))
          end do
        end do
      end if
    else if (array_notation) then
      out(:,:,i3) = in(:,:,j3  )*(1.-w(:,:,i3))&
                  + in(:,:,j3+1)*    w(:,:,i3)
    else
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          out(i1,i2,i3) = in(i1,i2,j3  )*(1.-w(i1,i2,i3))&
                        + in(i1,i2,j3+1)*    w(i1,i2,i3)
        end do
      end do
    end if
    if (verbose > 5) then
      i1=l(1); i2=l(2)
      write (io_unit%log,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j3,in(i1,i2,j3),in(i1,i2,j3+1),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k3 = j3
    if (mod(i3,2)==1) &
      j3 = j3+1
  end do
  if (verbose > 3) &
    write (io_unit%log, '(" Z3",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'x ', l(1),u(1), 'y ', l(2),u(2), 'zs',ls,k3, 'z ', l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE z_interp3

!===============================================================================
SUBROUTINE print (self, a, label, j)
  class(prolong_t):: self
  real:: a(:,:,:)
  character(len=*), optional:: label
  integer:: j
  character(len=32):: fmt
  !............................................................................
  integer:: i, idir
  integer, allocatable:: ival(:,:)
  !----------------------------------------------------------------------------
  ! Print array slices
  !----------------------------------------------------------------------------
  write (fmt,'(a,i2,a)') '(3(4x,',self%gn,'(1x,i3)))'
  if (verbose > 4) then
    if (present(label)) write (io_unit%log,*) label
    allocate (ival(self%gn,3))
    do i=1,self%gn
      do idir=1,3
        select case (idir)
        case (1)
          ! -- yz ----------------------
          ival(:,idir) = nint(a(j,:,j))
        case (2)
          ! -- xz ----------------------
          ival(:,idir) = nint(a(:,j,i))
        case (3)
          ! -- xy ----------------------
          ival(:,idir) = nint(a(:,i,j))
        end select
      end do
      write (stdout,fmt) ival
    end do
    write (stdout,*) ''
    deallocate (ival)
  end if
END SUBROUTINE print

!===============================================================================
!> Patch test, for debugging and validation
!===============================================================================
SUBROUTINE patch_test (self, n, ng)
  class(prolong_t)                    :: self
  real, allocatable, dimension(:,:,:) :: source, target
  integer, allocatable, dimension(:)  :: ival
  integer                             :: n, ng, gn, i, i1, i2, i3, j, k
  !-----------------------------------------------------------------------------
  ! Start up, allocate source and target arrays
  !-----------------------------------------------------------------------------
  gn = n + ng*2
  call self%init (n, ng)
  if (verbose >= 0) &
    write(stdout,'(a)') 'patch_test'
  allocate (source(gn,gn,gn), target(gn,gn,gn), ival(gn))
  !-----------------------------------------------------------------------------
  ! Loop over all 5 possible relative x-positions
  !-----------------------------------------------------------------------------
  target = 1000
  j = self%gn/2                         ! z-index
  do i3=-3,3,2                          ! varying z-position
  do i2=-3,3,2                          ! varying y-position
  do i1=-3,3,2                          ! varying x-position
    source = (i1+4)
    !call self%interpolate (source, target, [i1,i2,i3])
    if (verbose > 5) then
      do k=1,self%gn
        write(stdout,'(1x,20i6)') nint(target(:,k,j))
      end do
    end if
  end do
  end do
  end do
  call self%print (target,'target', j)
  call test_contamination (target)
  deallocate (source, target, ival)
  call self%finalize
  if (verbose >= 0) &
    write (stdout,'(a)') 'end patch_test'
END SUBROUTINE patch_test

!===============================================================================
SUBROUTINE test_contamination (target)
  real:: target(:,:,:)
  if (all(target == 1000. .or. outside)) then
    write (stdout,*) 'INFO: interior not contaminated'
  else
    write (stderr,*) 'ERROR: interior contaminated!'
  end if
END SUBROUTINE

!===============================================================================
!> Download test, for validation and timing
!===============================================================================
SUBROUTINE dnload_test (self, n, ng)
  class(prolong_t)                      :: self
  real, allocatable, dimension(:,:,:,:) :: source
  real, allocatable, dimension(:,:,:)   :: target
  integer, allocatable, dimension(:)    :: ival
  integer                               :: n, ng, gn, i1, i2, i3, n_dn, iv, nv
  real(8)                               :: start, used
  !-----------------------------------------------------------------------------
  ! Start up, allocate source and target arrays
  !-----------------------------------------------------------------------------
  if (verbose >= 0) &
    write (stdout,'(a)') 'dnload_test'
  gn = n + ng*2
  call self%init (n, ng)
  allocate (source(gn,gn,gn,27), target(gn,gn,gn))
  target = 1000.
  source = 1.0
  n_dn = 0
  do i3=-3,3,3
  do i2=-3,3,3
  do i1=-3,3,3
    n_dn = n_dn+1
    if (i1==0.and.i2==0.and.i3==0) cycle
    source (:,:,:,n_dn) = n_dn
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Loop over all 26 possible relative positions
  !-----------------------------------------------------------------------------
  start = omp_get_wtime()
  nv = 5
  do iv=1,nv
    n_dn = 0
    do i3=-3,3,3
    do i2=-3,3,3
    do i1=-3,3,3
      if (i1==0.and.i2==0.and.i3==0) cycle
      n_dn = n_dn+1
      !call self%interpolate (source(:,:,:,n_dn), target, [i1,i2,i3])
    end do
    end do
    end do
  end do
  if (verbose >= 0) then
    used = omp_get_wtime()-start
    write (stdout,'(3(3x,a,":",i3),1x,1p,2(f8.3,1x,a))') &
      'nv',nv, 'n',n, 'n_dn',n_dn, used*1e3,'ms', used*1e6/n**3,'mus'
  end if
  !-----------------------------------------------------------------------------
  ! Again, with the same mem slot
  !-----------------------------------------------------------------------------
  start = omp_get_wtime()
  nv = 5
  n_dn = 1
  do iv=1,nv
    do i3=-3,3,3
    do i2=-3,3,3
    do i1=-3,3,3
      if (i1==0.and.i2==0.and.i3==0) cycle
      !call self%interpolate (source(:,:,:,n_dn), target, [i1,i2,i3])
    end do
    end do
    end do
  end do
  if (verbose >= 0) then
    used = omp_get_wtime()-start
    write (stdout,'(3(3x,a,":",i3),1x,1p,2(f8.3,1x,a))') &
      'nv',nv, 'n',n, 'n_dn',n_dn, used*1e3,'ms', used*1e6/n**3,'mus'
  end if
  call test_contamination (target)
  deallocate (source, target)
  call self%finalize
  if (verbose >= 0) &
    write (stdout,'(a)') 'end dnload_test'
END SUBROUTINE dnload_test

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(prolong_t):: self
  call self%patch_test   (4, 3)
  call self%dnload_test (32, 2)
  call self%dnload_test (16, 2)
END SUBROUTINE

END MODULE prolong_mod
