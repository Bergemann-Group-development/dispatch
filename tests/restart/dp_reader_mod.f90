MODULE dp_reader_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE index_mod
  USE link_mod
  USE trace_mod
  USE patch_mod
  USE mesh_mod
  USE task_mod
  USE random_mod
  implicit none
  private
  type, extends(patch_t) :: patch_list_t
    class(patch_list_t), pointer:: head
    class(patch_list_t), pointer:: next
    class(patch_list_t), pointer:: prev
    class(patch_list_t), pointer:: tail
  end type
  type, public:: dp_reader_t
    logical:: active =.false.
    character(len=32):: folder
    integer:: snaptime, ids(2)
    integer:: nw = 1
    integer:: iread=0, iwrite=0
    real(8):: time
    class(patch_list_t), pointer :: patches => null()
  contains
    procedure :: cleanup
    procedure :: init
    procedure :: pick_data
    procedure :: read_patches
  end type
  type(dp_reader_t), public:: dp_reader
CONTAINS

SUBROUTINE init(self)
  class(dp_reader_t) :: self
  logical, save:: first_time=.true.
  integer, save :: snaptime = 0
  integer, save:: ids(2) = 1
  real:: row(21)
  character(len=32):: folder='data/'
  character(len=64):: inputfile
  character(len=24):: fmt
  real(kind=8):: time,position(3),ds(3),vel(3),llc_nat(3),llc_cart(3),centre_nat(3)
  real :: quality
  integer:: level, iostat
  logical,save :: time_derivs = .false.
  namelist /snapshot_params/ folder, snaptime, ids, time_derivs
  !.....................................................................
  integer :: ii, recl
  !---------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, snapshot_params)
    if (io%master) write (*, snapshot_params)
    self%folder = folder
    self%snaptime = snaptime
    self%ids=ids
    self%active=.true.
    if(time_derivs) self%nw=2
    do ii=self%ids(1), self%ids(2)
      call self%read_patches(ii)
    end do
    if (io%verbose > 0) &
    print '(a,i4.4,a)', 'DP reader initialized. ',self%ids(2) - self%ids(1) + 1, ' patches read.'
  end if
  
END SUBROUTINE init


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
  open (io_unit%tmp, file=trim(inputfile)//'.txt', form='formatted', status='old', action='read', err=100)
  allocate(patch)
  read(io_unit%tmp,*), patch%time, patch%position, patch%ds, &
  patch%velocity, patch%llc_nat, patch%llc_cart, rdum, &
  patch%level, patch%quality
  read(io_unit%tmp,*), dum, nn, dum, dum, patch%n, patch%gn, patch%ng, patch%nv
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
100 continue  
END SUBROUTINE read_patches

SUBROUTINE pick_data(self, link, mem)
  class(dp_reader_t) :: self
  class(link_t) :: link
  real, dimension(:,:,:,:), pointer  :: mem
  !.....................................................................
  type(index_t):: idx
  class(task_t), pointer :: task
  class(patch_t), pointer :: patch, dpatch
  !---------------------------------------------------------------------
  
  task => link%task
  select type(task)
  class is(patch_t)
    patch => task
  end select
  
  
  
  


END SUBROUTINE pick_data

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
