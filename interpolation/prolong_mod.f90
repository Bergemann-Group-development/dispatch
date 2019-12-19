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
!>      ...-o---+---o---|---o---+---o---+---o---+---o---|---o---+---o-
!>               -+---+-|-+---+---+---+-          +---+-|-+---+-
!>                1   2   3   4   5   6           5   6   7   8 
!>       -+---+-|-+---+---+---+-          +---+-|-+---+-
!>        3   4   5   6   7   8           5   6   7   8 
!>                9  10  11  12
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
  implicit none
  private
  type, public:: prolong_t
    integer           :: case=0                         ! mod(ng,2)
    integer           :: ng=2                           ! guard zones
    integer           :: n=8                            ! internal
    integer           :: gn=12                          ! total
    real, allocatable :: weights(:,:,:)                 ! weight, precomputed
    contains
    procedure         :: init
    procedure         :: dealloc
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
    procedure         :: interp1
    procedure         :: interp3
    generic           :: interp => interp1, interp3
    procedure         :: interpolate
    procedure         :: test
    procedure         :: patch_test
    procedure         :: dnload_test
    procedure         :: print
  end type
  integer:: verbose=0
  type(prolong_t), public:: prolong
CONTAINS

!===============================================================================
!> Initalize weight factors, based on number of total and guard cells.  The
!> arrays are assumed to be cubical, for simplicity and efficiency.
!===============================================================================
SUBROUTINE init (self, n, ng)
  class(prolong_t):: self
  integer:: n, ng
  !.............................................................................
  integer:: gn, i1, i2, i3, iostat
  namelist /prolong_params/ verbose
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, prolong_params, iostat=iostat)
    call io%namelist_warning ('prolong_params', iostat)
    write (stdout, prolong_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  gn = n + ng*2
  self%n  = n
  self%ng = ng
  self%gn = gn
  self%case = mod(ng,2)
  allocate (self%weights(gn,gn,gn))
  select case (self%case)
  !-----------------------------------------------------------------------------
  ! Case 0 starts with w=0.75
  !-----------------------------------------------------------------------------
  case (0)
    do i3=1,gn
    do i2=1,gn
    do i1=1,gn
      self%weights(i1,i2,i3) = 0.25 + mod(i1+i2+i3  ,2)*0.50
    end do
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Case 1 starts with w=0.25
  !-----------------------------------------------------------------------------
  case (1)
    do i3=1,gn
    do i2=1,gn
    do i1=1,gn
      self%weights(i1,i2,i3) = 0.25 + mod(i1+i2+i3+1,2)*0.50
    end do
    end do
    end do
  end select
  call self%print (self%weights*100., 'weights', gn/2)
END SUBROUTINE

!===============================================================================
!> Deallocate the wights array (when the task update is finished)
!===============================================================================
SUBROUTINE dealloc (self)
  class(prolong_t):: self
  !-----------------------------------------------------------------------------
  deallocate (self%weights)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE info (self)
  class(prolong_t):: self
  write (stdout,'(3(3x,a,":",i4))') 'n', self%n, 'ng', self%ng, 'gn', self%gn
END SUBROUTINE

!===============================================================================
!> Prolongation of cell values for a target guard zone.  The location of the
!> source is given by the integer(3) ipos, where -1 means "left", +1 means
!> "right", and 0 means the same position.
!===============================================================================
SUBROUTINE interpolate (self, source, target, ipos)
  class(prolong_t):: self
  real, dimension(:,:,:) :: source, target
  integer, dimension(3)  :: ipos
  intent(in)             :: source
  intent(out)            :: target
  !.............................................................................
  integer                :: ls(3), us(3), lt(3), ut(3), i
  integer                :: lb, lo, li, ui, uo, ub, n, ng
  !-----------------------------------------------------------------------------
  ! Cf. the sketch at top for index choices
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    print *, 'interpolate'
  n  = self%n
  ng = self%ng
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
    ! -- target is to the left of the source
    case (-3)
      if (verbose > 2) &
        write(stdout,*) 'target in direction', i,' is to the left'
      lt(i) = uo
      ut(i) = ub
      ls(i) = li - (ng/2+1)
      us(i) = li
    ! -- target is in the left part of the source
    case (-1)
      if (verbose > 2) &
        write(stdout,*) 'target in direction', i,' is inside left'
      lt(i) = li
      ut(i) = li + (n-1)
      ls(i) = lo
      us(i) = li + (n/2)
    case (0)
      if (verbose > 2) &
        write(stdout,*) 'target in direction', i,' is inside left'
      lt(i) = li
      ut(i) = li + (n-1)
      ls(i) = lo
      us(i) = li + (n/2)
    ! -- target is in the right part of the source
    case (+1)
      if (verbose > 2) &
        write(stdout,*) 'target in direction', i,' is inside right'
      lt(i) = ui - (n-1)
      ut(i) = ui
      ls(i) = ui - (n/2)
      us(i) = uo
    ! -- target is to the right of the source
    case (+3)
      if (verbose > 2) &
        write(stdout,*) 'target in direction', i,' is to the right'
      lt(i) = lb
      ut(i) = lo
      ls(i) = ui
      us(i) = ui + (ng/2+1)
    end select
  end do
  if (verbose > 2) then
    call self%info
    write (stdout,'(4(3x,a,":",3i4))') 'lt', lt, 'ut', ut, 'ls', ls
    write (stdout,*) 'weights:', maxval(self%weights)
  end if
  call self%interp (source, self%weights, ls, us, target, lt, ut)
END SUBROUTINE interpolate

!===============================================================================
!> Compute the first and last source index, and the first target index
!===============================================================================
SUBROUTINE interp1 (self, source, w, ls, us, target, lt, ut)
  class(prolong_t):: self
  real, intent(in):: source(:,:,:), w(:)
  real, intent(out):: target(:,:,:)
  integer, intent(in):: lt(3), ut(3), ls(3), us(3)
  !............................................................................
  integer:: m(3), l(3), u(3)
  real, allocatable:: target1(:,:,:), target2(:,:,:)
  !----------------------------------------------------------------------------
  if (verbose > 3) &
    print *, 'interp1'
  m = shape(source)
  allocate (target1(m(1),m(2),m(3)), target2(m(1),m(2),m(3)))
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=ls(2); u(2)=us(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  call self%x_interp (source , w, ls(1), target1, l, u)
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  call self%y_interp (target1, w, ls(2), target2, l, u)
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=lt(3); u(3)=ut(3)        ! source range in z
  call self%z_interp (target2, w, ls(3), target , l, u)
  deallocate (target1, target2)
END SUBROUTINE interp1

!===============================================================================
!> Case with 3D weights
!===============================================================================
SUBROUTINE interp3 (self, source, w, ls, us, target, lt, ut)
  class(prolong_t):: self
  real, intent(in):: source(:,:,:), w(:,:,:)
  real, intent(out):: target(:,:,:)
  integer, intent(in):: lt(3), ut(3), ls(3), us(3)
  !............................................................................
  integer:: m(3), l(3), u(3)
  real, allocatable:: target1(:,:,:), target2(:,:,:)
  !----------------------------------------------------------------------------
  if (verbose > 3) &
    print *, 'interp3'
  m = shape(source)
  allocate (target1(m(1),m(2),m(3)), target2(m(1),m(2),m(3)))
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=ls(2); u(2)=us(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  call self%x_interp (source , w, ls(1), target1, l, u)
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=ls(3); u(3)=us(3)        ! source range in z
  call self%y_interp (target1, w, ls(2), target2, l, u)
  l(1)=lt(1); u(1)=ut(1)        ! target range in x
  l(2)=lt(2); u(2)=ut(2)        ! source range in y
  l(3)=lt(3); u(3)=ut(3)        ! source range in z
  call self%z_interp (target2, w, ls(3), target , l, u)
  deallocate (target1, target2)
END SUBROUTINE interp3

!===============================================================================
!> Interpolate in the x-direction
!===============================================================================
SUBROUTINE x_interp1 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j1, n(3)
  real:: p, q
  !----------------------------------------------------------------------------
  n = shape(out)
  !--- reset the target output range ------------------------------------------
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  j1 = ls
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        out(i1,i2,i3) = out(i1,i2,i3) + in(j1  ,i2,i3)*(1.-w(i1))&
                                      + in(j1+1,i2,i3)*    w(i1)
      end do
    end do
    if (mod(i1,2)==1) &
      j1 = j1+1
  end do
END SUBROUTINE
!===============================================================================
SUBROUTINE x_interp3 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j1, k1, n(3)
  real:: p, q
  !----------------------------------------------------------------------------
  if (verbose > 3) &
    print *, 'x_interp3'
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  j1 = ls
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        out(i1,i2,i3) = out(i1,i2,i3) + in(j1  ,i2,i3)*(1.-w(i1,i2,i3))&
                                      + in(j1+1,i2,i3)*    w(i1,i2,i3)
      end do
    end do
    if (verbose > 5) then
      i2=l(2); i3=l(3)
      write (stdout,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j1,in(j1  ,i2,i3),in(j1+1,i2,i3),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k1 = j1+1
    if (mod(i1,2)==1) &
      j1 = j1+1
  end do
  if (verbose > 3) &
    write (stdout,'("X",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'xs', ls,k1, 'x', l(1),u(1), 'y', l(2),u(2), 'z',l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE

!===============================================================================
!> Interpolate in the y-direction
!===============================================================================
SUBROUTINE y_interp1 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j2, n(3)
  !----------------------------------------------------------------------------
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  j2 = ls
  do i2=l(2),u(2)
    do i3=l(3),u(3)
      do i1=l(1),u(1)
        out(i1,i2,i3) = out(i1,i2,i3) + in(i1,j2  ,i3)*(1.-w(i1))&
                                      + in(i1,j2+1,i3)*    w(i1)
      end do
    end do
    if (mod(i2,2)==1) &
      j2 = j2+1
  end do
  !print '("Y",5x,3(a,":",2i4,5x),f10.6)', &
    !'x', l(1),u(1), 'y', lo,j2-1, 'z',l(3),u(3), sum(out)/product(shape(out))
END SUBROUTINE
!===============================================================================
!===============================================================================
SUBROUTINE y_interp3 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j2, k2, n(3)
  !----------------------------------------------------------------------------
  if (verbose > 3) &
    print *, 'y_interp3'
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  j2 = ls
  do i2=l(2),u(2)
    do i3=l(3),u(3)
      do i1=l(1),u(1)
        out(i1,i2,i3) = out(i1,i2,i3) + in(i1,j2  ,i3)*(1.-w(i1,i2,i3))&
                                      + in(i1,j2+1,i3)*    w(i1,i2,i3)
      end do
    end do
    if (verbose > 5) then
      i1=l(1); i3=l(3)
      write (stdout,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j2,in(i1,j2,i3),in(i1,j2+1,i3),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k2 = j2+1
    if (mod(i2,2)==1) &
      j2 = j2+1
  end do
  if (verbose > 3) &
    write (stdout,'("Y",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'x', l(1),u(1), 'ys', ls,k2, 'y', l(2),u(2), 'z',l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE

!===============================================================================
!> Interpolate in the z-direction
!===============================================================================
SUBROUTINE z_interp1 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j3, n(3)
  !----------------------------------------------------------------------------
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  j3 = ls
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        out(i1,i2,i3) = out(i1,i2,i3) + in(i1,i2,j3  )*(1.-w(i1))&
                                      + in(i1,i2,j3+1)*    w(i1)
      end do
    end do
    if (mod(i3,2)==1) &
      j3 = j3+1
  end do
  !print '("Z",5x,3(a,":",2i4,5x),3f10.6)', &
    !'x', l(1),u(1), 'y', l(2),u(2), 'z',lo,j3-1, sum(out)/product(shape(out)), &
    !minval(out), maxval(out)
END SUBROUTINE
!===============================================================================
!===============================================================================
SUBROUTINE z_interp3 (self, in, w, ls, out, l, u)
  class(prolong_t):: self
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), ls
  !............................................................................
  integer:: i1, i2, i3, j3, k3, n(3)
  !----------------------------------------------------------------------------
  if (verbose > 3) &
    print *, 'z_interp3'
  j3 = ls
  out(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        out(i1,i2,i3) = out(i1,i2,i3) + in(i1,i2,j3  )*(1.-w(i1,i2,i3))&
                                      + in(i1,i2,j3+1)*    w(i1,i2,i3)
      end do
    end do
    if (verbose > 5) then
      i1=l(1); i2=l(2)
      write (stdout,'(3i4,i6,1p,4e10.2)') &
        i1,i2,i3,j3,in(i1,i2,j3),in(i1,i2,j3+1),w(i1,i2,i3),out(i1,i2,i3)
    end if
    k3 = j3
    if (mod(i3,2)==1) &
      j3 = j3+1
  end do
  if (verbose > 3) &
    write (stdout, '("Z",5x,4(a,":",2i4,5x),1p,3e11.2)') &
    'x', l(1),u(1), 'y', l(2),u(2), 'zs',ls,k3, 'z', l(3),u(3), &
    minval(out(l(1):u(1),l(2):u(2),l(3):u(3))), &
    maxval(out(l(1):u(1),l(2):u(2),l(3):u(3)))
END SUBROUTINE

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(prolong_t):: self
  call self%patch_test (4, 3)
  call self%dnload_test (32, 2)
  call self%dnload_test (16, 2)
END SUBROUTINE

!===============================================================================
!> Patch test, for debugging and validation
!===============================================================================
SUBROUTINE patch_test (self, n, ng)
  class(prolong_t)                    :: self
  real, allocatable, dimension(:,:,:) :: source, target
  integer, allocatable, dimension(:)  :: ival
  integer                             :: n, ng, gn, i, i1, i2, i3, j
  !-----------------------------------------------------------------------------
  ! Start up, allocate source and target arrays
  !-----------------------------------------------------------------------------
  gn = n + ng*2
  call self%init (n, ng)
  if (verbose >= 0) &
    write(stdout,*) 'patch_test'
  allocate (source(gn,gn,gn), target(gn,gn,gn), ival(gn))
  !-----------------------------------------------------------------------------
  ! Loop over all 5 possible relative positions
  !-----------------------------------------------------------------------------
  target = 10000
  i2 = 1                                ! mid-position y
  i3 = 1                                ! mid-position z
  j = self%gn/2                         ! z-index
  do i1=-3,3,2                          ! varying x-position
    source = (i1+4)
    if (verbose > 0) &
      call self%print (source,'source', self%gn/2)
    call self%interpolate (source, target, [i1,i2,i3])
    call self%print (target,'target', j)
  end do
  deallocate (source, target, ival)
  call self%dealloc
  if (verbose >= 0) &
    write (stdout,*) 'end patch_test'
END SUBROUTINE patch_test

!===============================================================================
!> Download test, for validation and timing
!===============================================================================
SUBROUTINE dnload_test (self, n, ng)
  class(prolong_t)                    :: self
  real, allocatable, dimension(:,:,:) :: source, target
  integer, allocatable, dimension(:)  :: ival
  integer                             :: n, ng, gn, i1, i2, i3, n_dn
  real(8)                             :: start, used
  !-----------------------------------------------------------------------------
  ! Start up, allocate source and target arrays
  !-----------------------------------------------------------------------------
  if (verbose >= 0) &
    write (stdout,*) 'dnload_test'
  gn = n + ng*2
  call self%init (n, ng)
  allocate (source(gn,gn,gn), target(gn,gn,gn))
  target = 1.0
  source = 1.0
  !-----------------------------------------------------------------------------
  ! Loop over all 26 possible relative positions
  !-----------------------------------------------------------------------------
  start = omp_get_wtime()
  n_dn = 0
  do i3=-3,3,3
  do i2=-3,3,3
  do i1=-3,3,3
    if (i1==0.and.i2==0.and.i3==0) cycle
    call self%interpolate (source, target, [i1,i2,i3])
    n_dn = n_dn+1
  end do
  end do
  end do
  if (verbose >= 0) then
    used = omp_get_wtime()-start
    write (stdout,'(2(3x,a,i3),1x,1p,2(f8.3,1x,a))') &
      'n:',n, 'n_dn:',n_dn, used*1e3,'ms', used*1e6/n**3,'mus'
    deallocate (source, target)
    call self%dealloc
    write (stdout,*) 'end dnload_test'
  end if
END SUBROUTINE dnload_test

!===============================================================================
SUBROUTINE print (self, a, label, j)
  class(prolong_t):: self
  real:: a(:,:,:)
  character(len=*), optional:: label
  integer:: j
  !............................................................................
  integer:: i
  integer, allocatable:: ival(:)
  !----------------------------------------------------------------------------
  ! Print array slices
  !----------------------------------------------------------------------------
  if (verbose > 4) then
    if (present(label)) write (stdout,*) label
    allocate (ival(self%gn))
    do i=1,self%gn
      select case (i)
      case (1)
        ! -- yz ----------------------
        ival(:) = nint(a(j,:,i))
      case (2)
        ! -- xz ----------------------
        ival(:) = nint(a(:,j,i))
      case (3)
        ! -- xy ----------------------
        ival(:) = nint(a(:,i,j))
      end select
      write (stdout,'(4x,20(1x,i3),$)') ival
    end do
    write (stdout,*) ''
    deallocate (ival)
  end if
END SUBROUTINE

END MODULE prolong_mod
