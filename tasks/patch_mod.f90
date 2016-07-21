!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mesg_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE task_mod
  USE mesh_mod
  USE buffered_output_mod
  implicit none
  private
  integer, parameter:: ndim = 3, i_d=1
  type, extends(task_t), public:: patch_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    real, dimension(:,:,:,:,:), pointer:: vgas
    real, dimension(:,:,:,:), pointer:: density
    real, dimension(:,:,:,:), pointer:: pressure
    type(mesh_t), dimension(3):: mesh
    real:: courant                      ! courant number
    real:: gamma=1.0                    ! gas gamma
    real:: dt_fixed, u_fixed, u_max     ! for fixed time step Courant condition
    integer:: ndim = ndim               ! number of dimensions
    integer:: nv=5                      ! number of variables
    integer:: nw=2                      ! number of work indices
    integer:: msplit=2                  ! how to split a patch when refining
    integer:: max_files                 ! prevent disk overrun when NaN happens
    logical:: guard_zones, time_derivs  ! switches for I/O
    logical:: root_grid=.false.
    logical:: mem_allocated=.false.
    logical:: do_pebbles=.false.
    !---------------------------------------------------------------------------
    ! HACK below:  FIXME
    !---------------------------------------------------------------------------
    real:: grav                         ! gravitational consant (code units)
    integer:: i_d                       ! density inddex
    real(8):: box(3)=1d0                ! the periodic wrapping size
    real(8):: ds(3), gsize(3)           ! convenience
    integer:: n(3)=0, ng(3)=2           ! convient default
    integer:: li(3), ui(3), gn(3)       ! convenience -- avoid patch%mesg%li etc
    integer:: i,j,k                     ! should go!  FIXME
    real(8):: offset(3)                 ! this is used to compute indices
  contains
    procedure:: check_error
    procedure:: init
    procedure:: rotate
    procedure:: output                  ! legacy output
    procedure:: output_new              ! new output
    procedure:: output_legacy           ! legacy output
    procedure:: pack                    ! pack patch into buffer
    procedure:: unpack                  ! unpack buffer into patch
    procedure:: info                    ! print info to stdout
    procedure:: distance                ! return the distance btw two patches
    procedure:: point_distance          ! distance btw two points
    procedure:: point_distances         ! distances btw two points
    procedure:: overlaps                ! return true if they overlap (incl. ghost zones)
    procedure:: overlaps_point          ! return true if they overlap (incl. ghost zones)
    procedure:: upload_range            ! returns the index ranges btw two points
    procedure:: dnload_range            ! returns the index ranges btw two points
    procedure:: index                   ! compute the index for a single point
    procedure:: nearest                 ! the nearest index that maps to patch point
    procedure:: lowest                  ! the lowest index that maps to patch point
    procedure:: highest                 ! the highest index that maps to patch point
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: myposition
    procedure:: periodic_grid
    procedure:: extrapolate_guards
    procedure:: fmaxval                 ! interior max value
    procedure:: fminval                 ! interior min value
    procedure:: courant_condition
    procedure:: check_density
    procedure:: patch_to_header
    procedure:: header_to_patch
    procedure:: local_iit
    procedure:: time_interpolate
  end type
  !-----------------------------------------------------------------------------
  ! Sequential header, for MPI and I/O
  !-----------------------------------------------------------------------------
  type header_t
    sequence
    integer:: n_header
    integer:: n_data
    integer:: id
    integer:: status
    integer:: rank
    integer:: istep
    integer:: iout
    integer:: level
    integer:: nt
    integer:: nv                        ! number of variables
    integer:: it
    integer:: padding                   ! for alignment
    real(8):: time
    real(8):: end_time
    real(8):: out_next
    real(8):: out_time
    real(8):: print_next
    real(8):: print_time
    real(8):: dtime
    real(8):: position(3)
    real(8):: velocity(3)
    real(8):: size(3)
    real(8), dimension(8):: t
    real(8), dimension(8):: dt
    integer, dimension(8):: iit
    real(8):: box(3)                    ! the periodic wrapping size
    real(8):: ds(3)
    real(8):: gsize(3)                  ! convenience
    integer:: n(3)=16                   ! convenience
    integer:: ng(3)=2                   ! convenience
    integer:: nq=0
    integer:: npad=0
  end type
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type box_t
    real:: size(3)=1d0
    real:: num(3)
    integer:: dim(3)
  end type
CONTAINS

!> =============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!> =============================================================================
SUBROUTINE patch_to_header (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  head%n_header   = sizeof(head)/4
  head%n_data     = product(self%gn)*self%nv
  head%id         = self%id
  head%status     = self%status
  head%rank       = self%rank
  head%istep      = self%istep
  head%iout       = self%iout
  head%level      = self%level
  head%time       = self%time
  head%end_time   = self%end_time
  head%out_next   = self%out_next
  head%out_time   = self%out_time
  head%print_next = self%print_next
  head%print_time = self%print_time
  head%dtime      = self%dtime
  head%position   = self%position
  head%velocity   = self%velocity
  head%size       = self%size
  head%it         = self%it
  head%nt         = self%nt
  nt              = self%nt
  head%t(1:nt)    = self%t
  head%dt(1:nt)   = self%dt
  head%iit(1:nt)  = self%iit
  head%nq         = self%nq
  if (io%verbose > 1) &
    write (io_unit%log,*) 'patch_to_header: id, nq =', head%id, head%nq, self%is_set(bits%swap_request)
END SUBROUTINE patch_to_header

!> =============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!> =============================================================================
SUBROUTINE header_to_patch (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  self%id         = head%id
  self%status     = head%status
  self%rank       = head%rank
  self%istep      = head%istep
  self%iout       = head%iout
  self%level      = head%level
  self%time       = head%time
  self%end_time   = head%end_time
  self%out_next   = head%out_next
  self%out_time   = head%out_time
  self%print_next = head%print_next
  self%print_time = head%print_time
  self%dtime      = head%dtime
  self%position   = head%position
  self%velocity   = head%velocity
  self%size       = head%size
  self%it         = head%it
  self%nt         = head%nt
  nt              = head%nt
  self%t          = head%t(1:nt)
  self%dt         = head%dt(1:nt)
  self%iit        = head%iit(1:nt)
  self%nq         = head%nq
  if (io%verbose > 1) &
    write (io_unit%log,*) 'header_to_patch: id, nq =', self%id, self%nq, self%is_set(bits%swap_request)
END SUBROUTINE header_to_patch

!> =============================================================================
!> Check for y-redundancy
!> =============================================================================
SUBROUTINE check_error (self, a, label, error)
  class(patch_t):: self
  character(len=*), optional:: label
  integer, optional:: error(3)
  integer:: loc(3), it, l(3), u(3), n(3), iy
  real, dimension(:,:,:), pointer:: a, b
  integer, save:: nerror=4
  !-----------------------------------------------------------------------------
  if (.not.io%do_debug) return
  l = self%mesh%lb
  u = self%mesh%ub
  n = u-l+1
  allocate (b(n(1),n(2),n(3)))
  do iy=l(2),u(2)
    b(:,iy,:)=a(:,iy,:)-a(:,1,:)
  end do
  if (any(b /= 0.0)) then
    loc = maxloc(abs(b)) + l-1
    if (present(label)) print*,'patch_t::error ', label
    print*,'patch_t::error loc        :', loc
    print*,'patch_t::error it,new,iit :', self%it, self%new, self%iit
    do iy=l(2),u(2)
      print*,'patch_t::error iy, a:', iy, a(loc(1),iy,loc(3))
    end do
    if (present(error)) then
      error=loc
      return
    else
      if (nerror<=0) stop
      nerror=nerror-1
    end if
  else
    if (present(label)) then
      print*,'patch_t::check_error maxval ', label, maxval(abs(a))
    else
      print*,'patch_t::check_error maxval', maxval(abs(a))
    end if
  end if
  deallocate(b)
END SUBROUTINE

!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  integer, save:: n(3), ng(3), nv, nt, nw, max_files=100
  logical, save:: first_time = .true.
  logical, save:: guard_zones=.false., time_derivs=.false.
  real, save:: u_fixed=-1.0, dt_fixed=-1.0, grace=0.0, courant=0.2
  real, save:: end_time=1.0, out_time=0.1, print_time=0.1
  real(8):: position(3), size(3)
  integer:: iostat, i
  namelist /patch_params/ size, n, nv, nw, ng, nt, courant, u_fixed, dt_fixed, grace
  namelist /out_params/ out_time, end_time, print_time, max_files, &
    guard_zones, time_derivs
  !.............................................................................
  call trace_begin ('patch_mod::init', 2)
  call self%task_t%init                         ! sets unique id
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters; these may have been set by the calling code
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    size = self%size
    n    = self%n
    ng   = self%ng
    nv   = self%nv
    nw   = self%nw
    nt   = self%nt
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
    rewind (io%input)
    read (io%input, out_params)
    if (io%master) write (*, out_params)
  end if
  !print *, 'MK', omp_mythread, self%id, self%nt, nt
  !$omp end critical (input_cr)
  ! reasonable values of self%n should take precedence over a saved value of `n`.
  if (any(self%n == 0)) then
    self%n           = n                        ! patch dimensions
  end if
  self%ng          = ng                         ! number of guard zones
  self%nv          = nv                         ! number of variables
  self%nw          = nw                         ! number of work indices
  self%nt          = nt                         ! number of time slices
  self%courant     = courant                    ! courant number
  self%grace       = grace                      ! grace fraction of dtime
  self%dt_fixed    = dt_fixed
  self%u_fixed     = u_fixed
  self%max_files   = max_files                  ! prevent run-away NaN
  !-----------------------------------------------------------------------------
  ! Output related parameters
  !-----------------------------------------------------------------------------
  self%guard_zones = guard_zones
  self%time_derivs = time_derivs
  self%end_time    = end_time
  self%out_time    = out_time
  self%print_time  = print_time
  self%print_next  = print_time
  self%ds = real(self%size,kind=8) / self%n     ! ds in patch
  if (all(self%size==self%box)) then
    write (io_unit%log,*) mpi%rank, 'task', self%id, ' is root grid'
    self%root_grid = .true.
    call self%set(bits%root_grid)
  end if
  !-----------------------------------------------------------------------------
  ! Special handling of collapsed dimensions
  !-----------------------------------------------------------------------------
  self%mesh%ng = self%ng
  where (self%n==1)                             ! for collapsed dimensions
    self%mesh%ng = 0                            ! no guard zones
    self%box = 1                                ! make sure it is non-zero
    self%ds = minval(self%ds)                   ! then set to reasonable val
    self%size = self%ds                         ! one cell thick
  end where
  !-----------------------------------------------------------------------------
  ! Initialize the meshes, afterwards imposing the external values box size and
  ! position
  !-----------------------------------------------------------------------------
  self%mesh%id = self%id
  self%mesh%n  = self%n
  self%mesh%s  = self%size
  do i=1,3
    call self%mesh(i)%init
    self%mesh(i)%b = self%box(i)
    self%mesh(i)%p = self%position(i)
  end do
  !-----------------------------------------------------------------------------
  ! Hack for backwards compatibility
  !-----------------------------------------------------------------------------
  self%gn = self%mesh%gn
  self%li = self%mesh%li
  self%ui = self%mesh%ui
  self%gsize = (self%size*self%gn/self%n)
  self%offset = (self%mesh%lb+self%mesh%ub)/2.0 ! mid point index offset
  !-----------------------------------------------------------------------------
  ! Measure levels in terms of the reinement factor
  !-----------------------------------------------------------------------------
  self%level = log(maxval(self%box)/minval(self%ds))/log(real(self%msplit)) + 0.1
  if (io%verbose>2) &
    write (io_unit%log,'(1x,a,i7,2(1p,3e12.4,2x))') 'patch_mod::init: id, size, pos =', &
      self%id, self%size, self%position
  if (io%verbose>3) &
    write (io_unit%log,*) 'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not already allocated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
!    if (associated(self%mem)) deallocate (self%mem)
!    if (self%rank == mpi%rank .or. self%is_set (bits%virtual)) then
      allocate (self%mem(self%mesh(1)%gn, &
                         self%mesh(2)%gn, &
                         self%mesh(3)%gn, &
                         self%nv,self%nt,self%nw))
      self%mem_allocated = .true.               ! mark
!    end if
  end if
  self%dt = 0.0
  self%t =  0.0
  self%iit = 1                                  ! initial memory slots
  self%iit(nt) = 2                              ! first output slot
  self%it = self%iit(nt-1)                      ! short-hand (=1 initially)
  self%new = self%iit(nt)                       ! short-hand
  self%t(self%it) = 0.0                         ! the only valid slot (no. 1)
  self%wallclock = wallclock()                  ! starting time
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE rotate (self)
  class(patch_t):: self
  integer:: iv, nv, new
  real(8):: time
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin ('patch_mod::rotate',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,'(a,7i10)')   'patch_mod::rotate iit =', self%iit
    write (io_unit%log,'(a,7f10.6)') 'patch_mod::rotate   t =', self%t(self%iit)
  end if
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation
  !-----------------------------------------------------------------------------
  if (self%root_grid.or.self%is_set(bits%root_grid)) then
    call self%periodic_grid                     ! periodic root grid
  else
    !call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  !-----------------------------------------------------------------------------
  ! One line summary
  !-----------------------------------------------------------------------------
  if (io%verbose>1) then
    write (io_unit%log,'(a,i6,2i3,2e13.6,2x,a,1p,6e10.2)') 'update: id,j,k,l,n,dt,t', &
      self%id, self%it, self%new, self%dtime, self%time, 'minmax:', &
      minval(self%mem(:,:,:,1,self%it ,self%nw)), &
      maxval(self%mem(:,:,:,1,self%it ,self%nw)), &
      minval(self%mem(:,:,:,1,self%new,1)), &
      maxval(self%mem(:,:,:,1,self%new,1))
  end if
  !-----------------------------------------------------------------------------
  ! Time slot rotation.  Note that self%time should be the LAST of these
  ! updates, since we do not want to enforce a critical region, and the
  ! new self%time is what a task is judged on!
  !-----------------------------------------------------------------------------
  call trace_tag ('patch_mod::update mem_cr')
  !call self%lock%set ('rotate')
  !!$omp flush
  time = self%time + self%dtime
  self%dt(self%it) = self%dtime                         ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = time                               ! initial time
  self%iit(self%nt) = self%new                          ! right-most slot
  new = self%new
  !$omp atomic write
  self%it = new                                         ! update time slot, new
  self%new = mod(self%new,self%nt)+1                    ! increment / rotate
  self%iit(1:self%nt-1) = self%iit(2:self%nt)           ! rotate left
  self%iit(self%nt) = self%new                          ! new right-most slot
  self%istep = self%istep + 1
  !$omp atomic write
  self%time = time
  if (io%verbose>2) &
   write (io_unit%log,*) 'patch id:',self%id,'  it:',self%it,'  new:',self%new
  if (io%verbose>3) &
   write (io_unit%log,*) 'patch iit:', self%iit
  if (io%verbose>3) &
   write (io_unit%log,*) 'patch t(iit):', self%t(self%iit)
  !!$omp flush
  !call self%lock%unset ('rotate')
  call trace_end (itimer)
END SUBROUTINE rotate

!> =============================================================================
!> This function gives the physical position of an index triple (patch point)
!> =============================================================================
FUNCTION myposition (self, ii) RESULT(out)
  class(patch_t):: self
  real(8):: out(3)
  integer:: ii(3)
  !.............................................................................
  out = self%position + (ii-self%offset)*self%ds
END FUNCTION

!> =============================================================================
!> Periodic mapping of a grid -- typically the root grid
!> =============================================================================
SUBROUTINE periodic_grid (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), j(3), lb(3), ub(3), li(3), ui(3)
  !.............................................................................
  call trace_begin ('patch_mod::periodic_grid')
  lb = self%mesh%lb
  ub = self%mesh%ub
  li = self%mesh%li
  ui = self%mesh%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    i = [ix,iy,iz]
    if (all(i >= li .and. i <= ui)) cycle
    j = modulo(i-li,self%n)+li
    self%mem(ix,iy,iz,:,self%it,1) = self%mem(j(1),j(2),j(3),:,self%it,1)
  end do
  end do
  end do
  call trace_end
END SUBROUTINE periodic_grid

!> =============================================================================
!> Time extrapolate guard cell values, to allow using also guard zone values to
!> get boundary conditions for higher level patches.  This is done before time
!> slot rotation, so the slot that has just been updated is iit(nt).
!> =============================================================================
SUBROUTINE extrapolate_guards (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), it0, it1, it2, nt, iit(self%nt)
  !.............................................................................
  call trace_begin ('patches_mod:extrapolate_guards', 3)
  nt = self%nt
  iit = self%local_iit()
  it0 = iit(nt)                                         ! current time slot
  it1 = iit(nt-1)                                       ! previous time slot
  it2 = iit(nt-2)                                       ! the one before that
  do iz=self%mesh(3)%lb,self%mesh(3)%ub                 ! loop over all cells
  do iy=self%mesh(2)%lb,self%mesh(2)%ub
  do ix=self%mesh(1)%lb,self%mesh(1)%ub
    i = [ix,iy,iz]
    if (all(i >= self%mesh%li .and. i <= self%mesh%ui)) cycle     ! skip interior
    self%mem(ix,iy,iz,:,it0,1)     = &
    self%mem(ix,iy,iz,:,it1,1)*2.0 - &
    self%mem(ix,iy,iz,:,it2,1)
  end do
  end do
  end do
  !print*,'min,max:',minval(self%mem(:,:,:,1,self%new,1)),maxval(self%mem(:,:,:,1,self%new,1))
  call trace_end
END SUBROUTINE extrapolate_guards

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(patch_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace_begin('patch_mod::output', itimer=itimer)
  do while (self%time > self%out_next .and. self%out_next <= self%end_time)
    if (io%do_legacy) then
      call self%output_legacy
    else
      call self%output_new
    end if
    if (self%out_time > 0.0) then
      self%out_next = (int(self%out_next/self%out_time)+1)*self%out_time
    else
      self%out_next = self%time
    end if
    self%iout = self%iout+1
  end do
  call trace_end (itimer)
END SUBROUTINE output

!> =============================================================================
!> Write results to disk, using the buffered I/O in output_mod
!> =============================================================================
SUBROUTINE output_new (self)
  class(patch_t):: self
  !.............................................................................
  type(header_t):: header
  integer:: n_header, n_buf, n_data, ibuf, it1, it2, iv
  real(4), dimension(:), pointer:: buffer
  real, dimension(:,:,:), pointer:: var
  real:: pt
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::output_new')
  n_header = sizeof(header)/4
  n_buf = product(self%gn)
  n_data = n_header + self%nv*n_buf
  allocate (buffer(n_data), var(self%gn(1),self%gn(2),self%gn(3)))
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer; these are first interpolated in time
  !-----------------------------------------------------------------------------
  it1 = self%iit(self%nt-2)
  it2 = self%iit(self%nt-1)
  pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1e-30)
  do iv=1,self%nv
    var = self%mem(:,:,:,iv,it1,1)*(1.-pt) + self%mem(:,:,:,iv,it2,1)*pt
    call anonymous_copy (n_buf, var, buffer(ibuf))
    ibuf = ibuf + n_buf
  end do
  deallocate (var)
  call buffered_output%buffer (self%id, self%iout, self%out_next, buffer)
  call buffered_output%check
  call trace_end
END SUBROUTINE output_new

!> =============================================================================
!> Pack patch info into a buffer, for MPI communication
!> =============================================================================
SUBROUTINE pack (self, mesg)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  integer:: n_header, n_var, ibuf, iv, it
  real, dimension(:,:,:), pointer:: var
  real:: pt
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::pack', itimer=itimer)
  allocate (mesg)
  n_header = sizeof(header)/4
  n_var = product(self%gn)
  mesg%nbuf = n_header + self%nv*n_var
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    mesg%nbuf = mesg%nbuf+ self%nv*n_var
  end if
  allocate (mesg%buffer(mesg%nbuf))
  allocate (mesg%reqs(self%n_nbors))
  mesg%id = self%id                                             ! MPI message tag
  mesg%nreq = 0                                                 ! reset counter
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,1) 'patch_mod::pack, id, nq, time, status:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request)
    1 format(a,2i9,f12.6,2x,5l1)
    flush (io_unit%log)
  end if
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer.  If the bdry+swap bit are set,
  ! copy as a minimum also the previous time slot
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    call anonymous_copy (n_var, self%mem(:,:,:,iv,self%it,1), mesg%buffer(ibuf))
    ibuf = ibuf + n_var
  end do
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    if (io%verbose>1) then
      write (io_unit%log,*) self%id, 'extra time slot packed', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_var, self%mem(:,:,:,iv,it,1), mesg%buffer(ibuf))
      ibuf = ibuf + n_var
    end do
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after pack', minval(self%mem(:,:,:,1,self%it,1))
  call trace_end (itimer)
END SUBROUTINE pack

!> =============================================================================
!> Pack patch info into a buffer, for MPI communication
!> =============================================================================
SUBROUTINE unpack (self, buffer)
  class(patch_t):: self
  !.............................................................................
  type(header_t):: header
  integer:: n_header, n_buf, ibuf, iv, it
  integer, dimension(:), pointer:: buffer
  real, dimension(:,:,:), pointer:: var
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::unpack', itimer=itimer)
  !write (io_unit%log,*) 'unpack patch id', self%id
  n_header = sizeof(header)/4
  n_buf = product(self%gn)
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info from sequenced header
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  ibuf = ibuf + n_header
  if (io%verbose > 2) write (io_unit%log,1) 'patch_mod::unpk id, nq, time before:', &
    self%id, self%nq, self%time
  call self%header_to_patch (header)
  if (io%verbose > 1) write (io_unit%log,1) 'patch_mod::unpk id, nq, time  after:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request)
    1 format(a,2i9,f12.6,2x,5l1)
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming patches
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%boundary)) then
    call self%set   (bits%virtual)
    call self%clear (bits%boundary)
  else if (self%is_set(bits%virtual)) then
    call self%set   (bits%boundary)
    call self%clear (bits%virtual)
  end if
  !-----------------------------------------------------------------------------
  ! Copy the variables from the output buffer.  If this is swap of patch roles
  ! expect at least one more time slot.
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,self%it,1))
    ibuf = ibuf + n_buf
  end do
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%virtual)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    if (io%verbose>0) then
      write (io_unit%log,*) self%id, 'extra time slot unpacked', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,it,1))
      ibuf = ibuf + n_buf
    end do
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after unpack', minval(self%mem(:,:,:,1,self%it,1))
  call trace_end (itimer)
END SUBROUTINE unpack

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output_legacy (self, experiment_name)
  class(patch_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  !-----------------------------------------------------------------------------
  call trace_begin('patch_mod::output_legacy')
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (self%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (self%time_derivs) then
    no = 2
  else
    no = 1
  end if

  !-----------------------------------------------------------------------------
  ! Write a log file, with files written by this rank
  !-----------------------------------------------------------------------------
  !$omp critical (io_cr)
  write (logname,'(a)') trim(io%outputname)//'log'
  if (first_time) then
    first_time = .false.
    if (io%master) print*, 'output data file log written to '//trim(logname)
    open (io%log_unit,file=trim(logname),status='unknown')
  else
    open (io%log_unit,file=trim(logname),status='unknown',position='append')
  end if
  write (filename,'(a,i4.4,"_",i4.4)') trim(io%outputname), self%id, self%iout
  write (io%log_unit,'(a)') trim(filename)
  ds = self%box/2.0**io%levelmax
  do i=1,3
    if (self%size(i) > 0.0) then
      li(i) = (self%position(i)+0.5*(self%box(i)-self%size(i)))/ds(i)
      ui(i) = (self%position(i)+0.5*(self%box(i)+self%size(i)))/ds(i) - 1.0
    end if
  end do
  write (io%log_unit,'(1p,e14.6,2(2x,3e14.6))') self%time, self%position, self%ds
  write (io%log_unit,'(4(3i7,2x))') int(self%box/self%ds+0.5), u-l+1, li, ui
  close (io%log_unit)

  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  if (io%master .and. io%verbose>=0) &
    print'(a,3x,"time:",f9.5)', 'patch_mod::output: file= '//trim(filename), self%out_next
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  do i=1,3
    if (self%size(i) > 0.0) then
      li(i) = (self%position(i)+0.5*(self%box(i)-self%size(i)))/self%ds(i)
      ui(i) = (self%position(i)+0.5*(self%box(i)+self%size(i)))/self%ds(i) - 1.0
    end if
  end do
  write (io%data_unit,'(1p,e14.6,2(2x,3e14.6))') self%time, self%position, self%ds
  write (io%data_unit,'(4(3i9,2x))') int(self%box/self%ds+0.5), u-l+1, li, ui
  close (io%data_unit)
  !$omp end critical (io_cr)

  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  !$omp critical (io_cr)
  n = u-l+1
  open (io%data_unit, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl=sizeof(self%mem(1,1,1,1,1,1))*product(n))
  do jo=1,no
    do iv=1,self%nv
      allocate (buf(n(1),n(2),n(3)))
      !-------------------------------------------------------------------------
      ! If we write out time derivatives we want a pair of consistent values
      ! and time derivatives, which we have at iit(nt-2).
      !-------------------------------------------------------------------------
      it1 = self%iit(self%nt-2)
      if (self%time_derivs) then
        it = it1
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)
        call check_error(self, buf, "output 1")
      !-------------------------------------------------------------------------
      ! Normally, we write out only variable values, which can then be interpolated
      ! to the exact output time, self%out_next, the update of which must come
      ! after this point
      !-------------------------------------------------------------------------
      else
        it2 = self%iit(self%nt-1)
        it = it2
        pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1e-30)
        if (self%t(it2)==self%t(it1)) pt=0.0
        qt = 1.0-pt
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)*qt &
            + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
        call check_error(self, buf, "output 2")
      end if
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*self%nv) buf
      deallocate (buf)
    end do
  end do
  close (io%data_unit)
  !$omp end critical (io_cr)

  call trace_end
END SUBROUTINE output_legacy

!> =============================================================================
!> Print info to stdout.  If
!>   verbose==0:        print only id=1 to stdout and io_unit%log
!>   verbose==1:        print only id=1 to stdout but all to io_unit%log
!>   verbose==2:        print iv=1 to stdout and io_unit%log for all id
!>   verbose==3:        print iv=1,nv to stdout and io_unit%log for all id
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(patch_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
  integer:: iv, nv
  !-----------------------------------------------------------------------------
  if (io%verbose>=0) then
    if (self%time>=self%print_next) then
      if (self%print_time > 0) &
        self%print_next = (nint(self%time/self%print_time)+1)*self%print_time
      if (io%verbose<1 .and. self%id/=1) return
      nv = 1
      if (io%verbose>2) nv = self%nv
      do iv=1,nv
        write (io_unit%log,1) &
          'upd',self%id, self%level, iv, omp_mythread, self%time, self%dtime, self%u_max, &
          self%fminval(self%mem(:,:,:,iv,self%it ,self%nw)), &
          self%fmaxval(self%mem(:,:,:,iv,self%it ,self%nw)), &
          self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
          self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
          real(n_update), ntask, nq, self%is_set(bits%boundary), self%is_set(bits%virtual), self%status
      end do
      if (io%verbose<2 .and. self%id/=1) return
      do iv=1,nv
        1 format(a,i8,2i4,i3,f10.5,1p,e12.4,0p,f7.2,2x,1p,2e10.2,2x,2e14.6,2x,e10.2,2i5,2l2,o12)
        if (io%master .or. self%id==1) &
          print 1, &
          'upd',self%id, self%level, iv, omp_mythread, self%time, self%dtime, self%u_max, &
          self%fminval(self%mem(:,:,:,iv,self%it ,self%nw)), &
          self%fmaxval(self%mem(:,:,:,iv,self%it ,self%nw)), &
          self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
          self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
          real(n_update), ntask, nq, self%is_set(bits%boundary), self%is_set(bits%virtual), self%status
      end do
    end if
  end if
END SUBROUTINE info

!> =============================================================================
!> Find the p1 target and p2 source index ranges that correspond to the overlap
!> region, for upload from p2 to p1 (p1 has lower level than p2)
!> =============================================================================
SUBROUTINE upload_range (p1, p2, l1, u1, l2, u2)
  class(patch_t):: p1
  class(patch_t):: p2
  integer, dimension(3):: l1, u1, l2, u2
  real, dimension(3):: f, pos1, pos2
  !-----------------------------------------------------------------------------
  ! Lower index limits
  !-----------------------------------------------------------------------------
  !print*, p1%id, ' target position, size:', p1%position, p1%size
  !print*, p2%id, ' source position, size:', p2%position, p2%size
  l2 = max(p2%mesh%li, p2%lowest  (p1, p1%mesh%lb))
  u2 = min(p2%mesh%ui, p2%highest (p1, p1%mesh%ub))
  l1 = p1%nearest (p2, l2)
  u1 = p1%nearest (p2, u2)
  !print*,'l1,l2', l1, l2
  !print*,'u1,u2', u1, u2
END SUBROUTINE upload_range

!> =============================================================================
!> Find the p1 target index ranges that correspond to the overlap region, for
!> download from p2 to p1 (p1 has higher or the same level as p2)
!> =============================================================================
SUBROUTINE dnload_range (target, source, l, u, o)
  class(patch_t):: target
  class(patch_t):: source
  integer, dimension(3):: l, u, o
  real, dimension(3):: f, pos1, pos2, dist
  !-----------------------------------------------------------------------------
  if (target%level == source%level) then
    l = target%mesh%li
    u = target%mesh%ui
    o = 0
    dist = source%distance(target)  ! positive if source to the right of target
    where (dist > 0.0)
      l = target%mesh%uo
      u = target%mesh%ub
      o = -target%mesh%n
    else where (dist < 0)
      l = target%mesh%lb
      u = target%mesh%lo
      o = target%mesh%n
    end where
  else
    print*,'dnload_range: ERROR, not implemented yet'
    stop
  end if
END SUBROUTINE dnload_range
!> =============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%mesh%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%mesh%lb and larger than self%mesh%ub.  For the root
!> grid patch, the index is in the range (li-1,ui) inside the interior.
!> =============================================================================
FUNCTION index (self, other, oindex, f)
  class(patch_t):: self, other
  real, dimension(3):: position, f, pos, pos_lb
  integer, dimension(3):: index, oindex
  !.............................................................................
  pos = other%position + (oindex-other%offset)*other%ds         ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  index = f                                                     ! smaller integer
  f = f - index                                                 ! remainder
END FUNCTION index
!> =============================================================================
!> The nearest index position in the self patch of the point in another patch.
!> The cost of this operation is 4 plus and 2 multiply per direction = 18 flops
!> =============================================================================
FUNCTION nearest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset)*other%ds          ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! nearest integer
END FUNCTION nearest
!> =============================================================================
!> The lowest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION lowest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset-0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! lowest integer
END FUNCTION lowest
!> =============================================================================
!> The highest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION highest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset+0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! highest integer
END FUNCTION highest

!> =============================================================================
!> Courant condition, one way or another
!> =============================================================================
SUBROUTINE courant_condition (self)
  class(patch_t):: self
  real(8):: dtime
  !.............................................................................
  call trace_begin('patch_mod::courant')
  if (io%verbose>2) &
    write (io_unit%log,*) 'courant: u_max, ds', self%u_max, minval(self%ds)
  if (self%dt_fixed > 0.0) then
    dtime = self%dt_fixed
    self%courant = dtime*self%u_max/minval(self%ds)
  else if (self%u_fixed > 0.0) then
    dtime = self%courant*minval(self%ds)/self%u_fixed
  else
    dtime = self%courant*minval(self%ds)/self%u_max
  endif
  !$omp atomic write
  self%dtime = dtime
  !-----------------------------------------------------------------------------
  ! Prevent large (cheap) patches from taking time steps longer than the out_time
  !-----------------------------------------------------------------------------
  if (self%out_time > 0.0) &
    self%dtime = min (self%dtime, self%out_time)
  self%dt(self%it) = self%dtime
  if (io%verbose>3) &
    print'(a,i7,f12.6,f8.4,1p,e10.2)', 'courant:', &
      self%id, self%dtime, self%courant, self%u_max
  call trace_end
END SUBROUTINE courant_condition

!> =============================================================================
!> Signed distance between the centers of two patches in the periodic box
!> =============================================================================
FUNCTION distance (self, patch) RESULT (out)
  class(patch_t):: self
  class(patch_t):: patch
  real(8):: out(3)
  !.............................................................................
  out = self%position-patch%position
  out = modulo(out+0.5_8*self%box,max(self%box,1e-30))-0.5_8*self%box
END FUNCTION distance

!> =============================================================================
!> Max absolute distance between two points in the periodic box
!> =============================================================================
FUNCTION point_distance (self, p1, p2) RESULT (out)
  class(patch_t):: self
  real(8):: p1(3), p2(3), out
  !.............................................................................
  out = maxval(abs(modulo(p1-p2+0.5_8*self%box,max(self%box,1e-30))-0.5_8*self%box))
END FUNCTION point_distance

!> =============================================================================
!> Max absolute distance between two points in the periodic box
!> =============================================================================
FUNCTION point_distances (self, p1, p2) RESULT (out)
  class(patch_t):: self
  real(8):: p1(3), p2(3), out(3)
  !.............................................................................
  out = abs(modulo(p1-p2+0.5_8*self%box,max(self%box,1e-30))-0.5_8*self%box)
END FUNCTION point_distances

!> =============================================================================
!> Check if there is overlap, using the full size (including guard zones), since
!> the guard zones may be used in interpolations for qualified time slices. The
!> definition of overlap also MUST be symmetric, since it is being used to define
!> neighbor relations.
!> =============================================================================
FUNCTION overlaps (self, patch)
  class(patch_t):: self
  class(patch_t):: patch
  logical:: overlaps
  !.............................................................................
  overlaps = all(abs(distance(self,patch)) <= 0.5000000001_8*(self%gsize+patch%gsize))
END FUNCTION overlaps

!> =============================================================================
FUNCTION overlaps_point (self, idir, point)
  class(patch_t):: self
  integer:: idir, adir
  real(8):: point(3)
  logical:: overlaps_point
  !.............................................................................
  overlaps_point = all(abs(point-self%position) < 0.5*(self%gn-1)*self%ds)
END FUNCTION overlaps_point

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_mod::set')
  self%status = ior(self%status, bit)
!  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_mod::clear')
  self%status = iand(self%status, not(bit))
!  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_mod::is_set')
  is_set = (iand(self%status, bit) /= 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_mod::is_clear')
  is_clear = (iand(self%status, bit) == 0)
!  call trace_end
END FUNCTION
!> =============================================================================
FUNCTION fmaxval (self, f)
  real:: fmaxval
  class(patch_t):: self
  real, dimension(:,:,:):: f
  fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION
!> =============================================================================
FUNCTION fminval (self, f)
  real:: fminval
  class(patch_t):: self
  real, dimension(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn):: f
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION

!> =============================================================================
!> Check if the density is negative somewhere, counting points in the interior
!> and guard zones separately.  Attempt to patch the guard zones, if applicable.
!> =============================================================================
SUBROUTINE check_density (self, label)
  class(patch_t):: self
  character(len=*), optional:: label
  real:: mind, posd
  integer:: ix, iy, iz, ii(3), ntot, nins
  logical:: inside
  !-----------------------------------------------------------------------------
  if (io%verbose < 1) return
  mind = minval(self%mem(:,:,:,1,self%it,1))
  if (present(label) .and. io%verbose>1) write (io_unit%log,*) label, self%id, mind
  if (mind <= 0.0) then
    !---------------------------------------------------------------------------
    ! Negative densities appear to develop at a few points in the guard zones,
    ! so this tries to patch up, sampling densitites in the guard zone.
    !---------------------------------------------------------------------------
    posd = self%mem(3,3,2,1,self%it,1)
    ntot = 0
    nins = 0
    do iz=1,self%gn(3)
    do iy=1,self%gn(2)
    do ix=1,self%gn(1)
      ii = [ix,iy,iz]
      inside = all([ix,iy,iz]>self%mesh%lo.and.ii<self%mesh%lo)
      posd = merge (posd, max(posd,self%mem(ix,iy,iz,1,self%it,1)), inside)
      ntot = ntot + merge(1,0,self%mem(ix,iy,iz,1,self%it,1)<=0.0)
      nins = nins + merge(1,0,self%mem(ix,iy,iz,1,self%it,1)<=0.0.and.inside)
      !-------------------------------------------------------------------------
      ! Pick the sampled positive value if a negative one is found in the
      ! guard zone
      !-------------------------------------------------------------------------
      !self%mem(ix,iy,iz,1,self%it,1) = merge (self%mem(ix,iy,iz,1,self%it,1), &
      !  posd, self%mem(ix,iy,iz,1,self%it,1) > 0.0 .or. inside)
    end do
    end do
    end do
    if (present(label)) then
      if (io%verbose>1 .or. nins>0) &
        write (io_unit%log,1) self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside, '// &
         trim(label)
      1 format(i8,3x,a,i8,a,i8,a)
    else
      if (io%verbose>1 .or. nins>0) &
        write (io_unit%log,1) self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside'
    end if
  end if
END SUBROUTINE check_density

!> =============================================================================
!> =============================================================================
FUNCTION local_iit (self) RESULT (iit)
  class(patch_t):: self
  integer:: i, it, iit(self%nt)
  !-----------------------------------------------------------------------------
  ! Make a consistent local copy of self%iit, based on self%it, which is
  ! written and read in atomic operations.  The array iit is such that its last
  ! value, iit(self%nt) is equal to self%new.
  !-----------------------------------------------------------------------------
  !!$omp atomic read
  !it = self%it
  !iit = [(1+mod(it+i,self%nt), i=1,self%nt)]
  iit = self%iit
END FUNCTION local_iit

!> =============================================================================
!> Find the right time interval in source%t(source%iit), which are the source
!> times in increasing order.  To avoid having to lock access to source%iit while
!> reading it, we read instead (with an atomic read) source%new, from which a
!> consisten local copy of itt can be constructed.
!> =============================================================================
SUBROUTINE time_interpolate (self, source, jt, pt)
  class(patch_t):: self, source
  integer:: jt(2)
  integer:: i, it, iit(self%nt)
  real(8):: pt(2), dt
  !.............................................................................
  call trace_begin ('task_t::time_interpolate',3)
  !-----------------------------------------------------------------------------
  ! In case of e.g. an error in initial conditions, choose a valid default
  !-----------------------------------------------------------------------------
  iit = source%local_iit()
  jt = iit(1:2)
  !-----------------------------------------------------------------------------
  ! As long as the first source time is earlier than, or equal to the target
  ! time, renew the interval.  Initially, when all times are zero, we are left
  ! with the interval nt-2:nt-1 and dt=0d0, which means slot nt-1 is used.
  !-----------------------------------------------------------------------------
  do i=2,source%nt-2
    if (source%t(iit(i)) <= self%time) then
      jt = iit(i:i+1)
    end if
  end do
  !-----------------------------------------------------------------------------
  ! A vanishing source time interval happens only at the initial condition ...
  !-----------------------------------------------------------------------------
  dt = source%t(jt(2))-source%t(jt(1))
  if (dt==0d0) then
    pt(2) = 1d0
  else
    pt(2) = (self%time - source%t(jt(1)))/dt
  end if
  pt(1) = 1d0-pt(2)
  if (io%verbose > 2) then
    write (io_unit%log,'(a,2i5,2f7.3,1p,g16.6)') &
      'patch_mod::time_interpolation jt, pt, target', &
      jt, pt, self%time
    do i=1,source%nt
      write (io_unit%log,'(a,3i4,1p,2g16.6)') &
        'patch_mod::time_interpolation i,iit,src%iit,src%t,src&dt:', &
        i, iit(i), source%iit(i), source%t(i), source%dt(i)
    end do
  end if
  !print'(a,5i3,3x,2i3,2x,2f7.3)', 'MK:', iit-source%iit, jt, pt
  !-----------------------------------------------------------------------------
  ! Facilitate exact time evolution when dt_fixed is being used
  !-----------------------------------------------------------------------------
  if (self%dt_fixed > 0.0) then
    if (abs(pt(1))< 1d-5) then
      if (io%verbose>2) write (io_unit%log,*) 'EXACT'
      pt(1)=0d0; pt(2)=1d0
    else if (abs(pt(2)) < 1d-5) then
      if (io%verbose>2) write (io_unit%log,*) 'EXACT'
      pt(1)=1d0; pt(2)=0d0
    end if
  end if
  call trace_end
END SUBROUTINE time_interpolate

END MODULE
