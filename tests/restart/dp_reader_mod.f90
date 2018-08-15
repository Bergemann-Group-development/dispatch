MODULE dp_reader_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE index_mod
  USE interpolate_mod
  USE link_mod
  USE trace_mod
  USE patch_mod
  USE mesh_mod
  USE task_mod
  USE random_mod
  implicit none
  private
  !---------------------------------------------------------------------
  type, extends(patch_t) :: patch_list_t
    class(patch_list_t), pointer:: head
    class(patch_list_t), pointer:: next
    class(patch_list_t), pointer:: prev
    class(patch_list_t), pointer:: tail
  end type
  !---------------------------------------------------------------------
  type, public:: dp_reader_t
    logical:: active =.false.
    character(len=32):: folder
    integer:: snaptime, ids(2)
    integer:: nw = 1, nread=0
    integer:: iread=0, iwrite=0
    real(8):: time
    class(patch_list_t), pointer :: patches => null()
    procedure(interpolator_interface), pointer, nopass:: interpolate => null()
  contains
    procedure :: cleanup
    procedure :: init
    procedure :: snapshot
    procedure :: read_patches
  end type
  !---------------------------------------------------------------------
  type(dp_reader_t), public:: dp_reader
CONTAINS

!> =====================================================================
!> 
!> =====================================================================
SUBROUTINE init(self)
  class(dp_reader_t) :: self
  logical, save:: first_time=.true.
  integer, save :: snap_number = 0
  integer, save:: ids(2) = 1
  real:: row(21)
  character(len=32):: folder='data/'
  integer:: level, iostat
  integer :: ord_int = -1
  logical,save :: time_derivs = .false.
  namelist /snapshot_params/ folder, snap_number, ids, time_derivs
  !.....................................................................
  integer :: ii, recl
  !---------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    self%interpolate => SelectInterpolator(ord_int)
    rewind (io%input); read (io%input, snapshot_params)
    if (io%master) write (*, snapshot_params)
    self%folder = folder
    self%snaptime = snap_number
    self%ids=ids
    self%active=.true.
    if(time_derivs) self%nw=2
    do ii=self%ids(1), self%ids(2)
      call self%read_patches(ii)
    end do
    if (io%verbose > 0) &
    print '(a,i4.4,a)', 'DP reader initialized. ',self%nread, ' patches read.'
  end if
  
END SUBROUTINE init

!> =====================================================================
!> Read in the data from input files
!> =====================================================================
SUBROUTINE read_patches(self, pid)
  USE iso_fortran_env, only : real32
  class(dp_reader_t) :: self
  integer :: pid
  !.....................................................................
  class(patch_list_t), pointer:: patch
  character(len=64):: inputfile
  character(len=24):: fmt = '(a,i4.4,"_",i4.4)'
  integer :: dum(3), nw, iv, iw, nn(3), l(3), u(3), n(3)
  real    :: rdum(3)
  real, dimension(:,:,:), pointer:: buf
  ! --------------------------------------------------------------------
  write(inputfile, fmt) trim(self%folder),pid, self%snaptime
  print *, inputfile
  open (io_unit%tmp, file=trim(inputfile)//'.txt', form='formatted', status='old', action='read', err=100)
  allocate(patch)
  read(io_unit%tmp,*), patch%time, patch%position, patch%ds, &
  patch%velocity, patch%llc_nat, patch%llc_cart, rdum, &
  patch%level, patch%quality
  read(io_unit%tmp,*), dum, nn, patch%li, patch%ui, patch%n, patch%gn, patch%ng, patch%nv
  read(io_unit%tmp,*), patch%id, patch%iout, patch%istep, patch%dtime, patch%mesh_type, patch%kind
  close (io_unit%tmp)
  if (all(nn==patch%gn)) then
    l = 1
    u = patch%gn
  else
    l = 1 + patch%ng
    u = patch%gn - patch%ng
  end if
  n=u-l+1
  open (io_unit%data, file=trim(inputfile)//'.dat', form='unformatted', status='old', action='read',access='direct', recl=real32*product(n))
  allocate (patch%mem(patch%gn(1),patch%gn(2),patch%gn(3),patch%nv,1,self%nw))
  allocate(buf(n(1),n(2),n(3)))
  do iw=1, self%nw
    do iv=1, patch%nv
      read (io_unit%data,rec=iv+(iw-1)*patch%nv) buf 
      patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,1,iw) = buf(:,:,:)
    end do
  end do
  close (io_unit%data)
  if (associated(self%patches)) then
    patch%prev => self%patches%tail
    patch%head => self%patches
    patch%prev%next => patch
    self%patches%tail => patch
  else
    self%patches => patch
    patch%head => patch
    patch%tail => patch
  end if
  deallocate(buf)
  self%nread = self%nread+1
100 continue  
END SUBROUTINE read_patches

!> =====================================================================
!> Pickup the data from loaded initial data sets
!> =====================================================================
SUBROUTINE snapshot (self, mesh, mem, level, used)
  class(dp_reader_t) :: self
  class(mesh_t), pointer, dimension(:) :: mesh
  real, dimension(:,:,:,:), pointer  :: mem
  integer :: level
  !.....................................................................
  class(task_t), pointer :: task
  class(patch_list_t), pointer :: dpatch, bpatch                                
  logical :: same_found, used
  integer :: i
  !---------------------------------------------------------------------
  
  same_found = .false.
  dpatch => self%patches
  
  do while (associated(dpatch).and..not.(same_found))
    if (the_same(mesh,dpatch, level)) then
      print *, '=== OK, same found ==='
      same_found = .true.
      bpatch => dpatch
    else
      dpatch => dpatch%next
    end if
  end do
  if (.not.same_found) then
    print *, 'same not found...'
    used = .false.
  end if
  
  if (same_found) then
    mem = dpatch%mem(:,:,:,:,1,1)
    used = .true.
  end if



END SUBROUTINE snapshot

LOGICAL FUNCTION the_same(m, p, level)
  class(mesh_t), pointer, dimension(:)  :: m
  class(patch_list_t), pointer :: p
  integer :: level
  !.....................................................................
  integer :: i
  real (kind=8) :: sz, pllc, purc, mllc, murc
  !---------------------------------------------------------------------
  the_same = .true.
  do i = 1,3
    if (m(i)%n > 1) then                                                ! ignore collapsed dimensions
      sz = m(i)%s/m(i)%n
      pllc = p%position(i)-0.5_8*(p%n(i)*p%ds(i))
      purc = p%position(i)+0.5_8*(p%n(i)*p%ds(i))
      mllc = m(i)%p-0.5_8*m(i)%s
      murc = m(i)%p+0.5_8*m(i)%s
    if (abs(m(i)%p - p%position(i)) > 1d-7) the_same = .false.
     ! if (sz /= p%ds(i))           the_same = .false.
     ! if (p1%level /= p2%level)            the_same = .false.
      if (m(i)%n /= p%n(i))        the_same = .false.
      if (abs(pllc-mllc)> 1d-7) the_same =.false.
      if (abs(purc-murc)> 1d-7) the_same =.false.
    end if
  end do
END FUNCTION the_same
  
!> =====================================================================
!> Remove the linked list of temporary patches
!> =====================================================================
SUBROUTINE cleanup (self)
  class(dp_reader_t) :: self
  class(patch_list_t), pointer :: this,  prev
  this => self%patches
  do while (associated(this))
    deallocate(this%mem)
    prev => this
    this => this%next
    deallocate(prev)
  end do
END SUBROUTINE cleanup  

END MODULE dp_reader_mod
