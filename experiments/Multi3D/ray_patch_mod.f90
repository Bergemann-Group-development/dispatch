MODULE ray_patch_mod
  USE mpi_mesg_mod
  USE halo_mod
  USE patch_mod, only: header_t, n_header
  USE gpatch_mod
  USE link_mod
  USE task_mod
  USE io_mod
  USE timer_mod
  USE scaling_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE shared_mod
  USE bits_mod
  USE os_mod
  USE math_mod
  USE handyfuncs_mod
  USE patch_mod, only: patch_t
  USE op_math_mod, only: searchsorted
  USE ray_mod, only: ray_type, one_dimensional
  USE m3d_patch_mod, only: m3d_patch_t, vec, wmu, long_vec, long_wmu, &
  converged, print_send, atom_update_counter, comoving_frame, use_vel
  USE const_mod, only: cc, over_4pi, PI
  implicit none
  private

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_patch_t) :: ray_patch_t
    type(ray_type) :: ray
    real(8) :: vec(3), wmu
    integer :: fs, fe, nf
    integer :: pre_pad_fs, pre_pad_fe
    logical :: padded, level_boundary=.false.
  contains
    procedure :: check_in_bounds
    procedure :: nbor_relations => ray_nbor_relations
    procedure :: set_freq_indices
    procedure :: is_upstream => ray_is_upstream
    procedure :: pack
    procedure :: unpack
    procedure :: allocate_mesg
    procedure :: total_mesg_buffer_size
    procedure :: init_nbors_hash
    procedure :: init_bdries
    procedure :: bndries
  end type ray_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0

contains

!===============================================================================
SUBROUTINE set_freq_indices (self)
  class(ray_patch_t) :: self
  integer :: fs, fe, i, nv
  integer, pointer :: counter
  character(len=16) :: label
  integer, dimension(:), allocatable :: order
  integer, save, target :: counter1=-999, counter2=-999
  !-----------------------------------------------------------------------------
  call unequal_intervals(self%sp%nnu, self%i_nu, self%n_nu, fs, fe)
  self%fs = fs
  self%fe = fe
  self%pre_pad_fs = fs
  self%pre_pad_fe = fe
  self%nf = self%fe - self%fs + 1

  if (self%out_ray) then
    counter => counter2
    label = 'output'
  else
    counter => counter1
    label = 'internal'
  endif
  if (self%i_nu > counter .and. self%verbose > 2) then
    if (counter == -999) write(stdout,'(I2, A10, A18)') self%n_nu, trim(label), 'frequency splits'
    write(stdout,'(I4, I10,A3,I8)') self%i_nu, fs, ' - ', fe
    counter = self%i_nu
  endif
  !if (self%has_atom .and. self%in_ray) then
  !  if (associated(self%atom%line) .and. .not. associated(self%atom%cont)) then
  !    nv = size(self%atom%line)
  !    call unequal_intervals(nv, self%i_nu, self%n_nu, fs, fe)
  !    allocate(order(nv))
  !    call argsort(self%atom%line%nu0, order)
  !    fs = self%atom%line(order(fs))%ired
  !    fe = self%atom%line(order(fe))%iblue
  !  endif
  !endif
END SUBROUTINE set_freq_indices

!===============================================================================
!> Determine if a task is upstream.  
!===============================================================================
logical FUNCTION ray_is_upstream (self, other) result(is_upstream)
  class(ray_patch_t):: self, other
  !-----------------------------------------------------------------------------
  is_upstream = self%i_mu == other%i_mu
  is_upstream = is_upstream .and. self%key(1)-self%dir == other%key(1)
END FUNCTION ray_is_upstream

!===============================================================================
!> Initialize boundary bits
!===============================================================================
SUBROUTINE init_bdries (self)
  class(ray_patch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%init_bdries')
  call self%m3d_patch_t%init_bdries()
  if (self%dir == 0) then
    call self%set(bits%no_download)
  elseif (self%dir > 0) then
    call self%boundaries%clear(bits%zu)
    if (self%boundaries%is_set(bits%zl)) then
      call self%set(bits%no_download)
    endif
  elseif (self%dir < 0) then
    call self%boundaries%clear(bits%zl)
    if (self%boundaries%is_set(bits%zu)) then
      call self%set(bits%no_download)
    endif
  endif
  call trace%end()
END SUBROUTINE init_bdries

!===============================================================================
SUBROUTINE ray_nbor_relations(self, another, needed, needs_me, download)
  class(ray_patch_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_nbor_relations')
  call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    select type(another)
    class is(ray_patch_t)
      if (self%kind==another%kind) then
        if (self%is_upstream(another)) then
          needs_me = .false.
          self%level_boundary = another%level /= self%level
        elseif (another%is_upstream(self)) then
          needed = .false.
          download = .false.
        elseif (self%ipos(3) /= another%ipos(3)) then
          needed = .false.; needs_me = .false.; download = .false.
        endif
      endif
    class is(m3d_patch_t)
      if (another%id /= self%m3d_parent%id) then
        needed = .false.; needs_me = .false.; download = .false.
        if (associated(another%m3d_parent)) then
          if (another%m3d_parent%id == self%m3d_parent%id) then
            ! could be an atom_patch needed by net_opac
            needed = .true.
          endif
        endif
      endif
    class default
      needed = .false.; needs_me = .false.; download = .false.
    end select
  endif

  call trace%end()
END SUBROUTINE ray_nbor_relations

!===============================================================================
!> Pack patch info into a buffer, for MPI communication
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE pack (self)
  class(ray_patch_t):: self
  !.............................................................................
  class(mesg_t), pointer:: mesg
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it, iw, iz
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_pack', itimer=itimer)
  if ((self%verbose > 2 .or. self%print_update_info) .and. print_send) then
    write(stdout, '(I6, A32, I3, A6, I10)') self%id, 'packing message on rank', & 
                                      mpi%rank, 'tag=', self%mesg%tag
  endif
  !-----------------------------------------------------------------------------
  ! Allocate a task%mesg instance, mesg%reqs() and similar array
  !-----------------------------------------------------------------------------
  if (omp%nthreads > 1) then
    allocate (self%mesg%lock)
    call self%mesg%lock%init ('smsg')
  end if
  mesg => self%mesg
  call os%allocate (mesg%reqs,  mpi%size, zero=.true.)
  call os%allocate (mesg%fail,  mpi%size, zero=.true.)
  call os%allocate (mesg%ranks, mpi%size, zero=.true.)
  !-----------------------------------------------------------------------------
  ! Packing meta-data is a task-owner read operation, and does not need locking
  !-----------------------------------------------------------------------------
  call self%patch_to_header (header)
  if (converged .and. self%in_ray) then
    header%level = 0
  endif
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  if (.not. converged .or. self%out_ray) then
    !call self%check_mem()
    !-----------------------------------------------------------------------------
    ! Copy the variables to the output buffer.
    !-----------------------------------------------------------------------------
    iz = merge(self%n(3), 1, self%dir > 0)
    n_buf = self%total_mesg_buffer_size() / (self%nw*self%nv)
    do iw=1,self%nw
      do iv=1,self%nv
        call anonymous_copy (n_buf, self%authority_mem(:,:,iz,iv,self%it,iw), mesg%buffer(ibuf))
        ibuf = ibuf + n_buf
      end do
    end do
  endif
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Unpack patch info from an MPI buffer.  To ensure consistency, the task meta-
!> data should not be updated until after the time slot data are in place.
!> Any task that has this task as virtual nbor and is waiting for the time slot
!> will thus not be scheduled until both the meta-data and the time slot data
!> are up-to-date.  This also reduces the chance of lock contention with respect
!> to the new time slot.
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE unpack (self, buffer)
  class(ray_patch_t):: self
  integer, pointer:: buffer(:)
  !.............................................................................
  type(header_t):: header
  type(lock_t), pointer:: lock
  real(8), dimension(3) :: size, position 
  integer:: n_buf, ibuf, iv, it, i, iw, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_unpack', itimer=itimer)
  call self%set (bits%busy)
  if ((self%verbose > 2 .or. self%print_update_info) .and. print_send) then
    write(stdout, '(I6, A32, I3, A6, I10)') self%id, 'unpacking message on rank', & 
                                      mpi%rank, 'tag=', self%mesg%tag
  end if
  self%unpack_time = wallclock()
  n_buf = self%total_mesg_buffer_size() / (self%nw*self%nv)
  !-----------------------------------------------------------------------------
  ! Copy buffer into sequenced header, to get the new %it slot
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  size = self%size
  position = self%position
  call self%header_to_patch (header)
  self%size = size
  self%position = position
  if (header%level == 0 .and. .not. converged) then
    converged = .true.
    call self%init_level()
    write(stdout, '(A, I3, A, I3)'), 'Converged message has been passed by MPI from rank', header%rank, ' to', mpi%rank
  endif
  ibuf = ibuf + n_header
  if (.not. converged .or. self%out_ray) then
    !call self%check_mem()
    !-----------------------------------------------------------------------------
    ! Copy array data from the buffer to the time slot, under time slot lock
    !-----------------------------------------------------------------------------
    lock => self%mem_lock(header%it)
    iz=1
                                                          call lock%set ('unpack')
    do iw=1,self%nw
      do iv=1,self%nv
        call anonymous_copy (n_buf, buffer(ibuf), self%authority_mem(:,:,iz,iv,header%it,iw))
        ibuf = ibuf + n_buf
      end do
    end do
                                                        call lock%unset ('unpack')
  endif
  !-----------------------------------------------------------------------------
  ! Update the patch meta-data under lock, now that the time slices are in place
  !-----------------------------------------------------------------------------
                                                   call self%lock%set ('unpack')
  self%it         = header%it
  self%time       = header%time
  if (header%nt > 0) then
    self%nt       = header%nt
    self%new      = mod(header%it,header%nt) + 1
    do i=1,self%nt
      self%t(i)   = header%t(i)
      self%dt(i)  = header%dt(i)
      self%iit(i) = header%iit(i)
    end do
  end if
  !.............................................................................
                                               call self%lock%unset ('unpack')
  call self%clear (bits%busy)
  self%wc_last = wallclock()
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE unpack

!===============================================================================
!> Allocate a message structure, for communication
!===============================================================================
SUBROUTINE allocate_mesg (self)
  class(ray_patch_t):: self
  class(mesg_t), pointer:: mesg
  integer:: n_buf
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%allocate_mesg')
  allocate (mesg)
  call os%mem_register (storage_size(mesg), 1)
  if (self%is_set(bits%boundary) .and. associated(self%mesg)) then
    if (self%mesg%id /= self%id .and. self%mesg%id > 0) then
      write (stdlog,*) &
        'allocate_mesg: wrong task%id,mesg%id=', self%id, self%mesg%id
      call trace%end ()
      return
    end if
  end if
  self%mesg => mesg
  call mesg%init
  !-----------------------------------------------------------------------------
  ! Store the standard buffer size (not valid for swap requests!)
  !-----------------------------------------------------------------------------
  n_buf = self%total_mesg_buffer_size()
  mesg%nbuf = n_header + n_buf
  if (associated(mesg%buffer)) then
    call os%deallocate(mesg%buffer)
  endif
  call os%allocate(mesg%buffer, mesg%nbuf)
  call trace%end ()
END SUBROUTINE allocate_mesg

!===============================================================================
FUNCTION total_mesg_buffer_size(self) result(n_buf)
  class(ray_patch_t):: self
  integer :: n_buf
  !-----------------------------------------------------------------------------
  n_buf = product(self%n(:2))
  if (kind(self%mem)==8) n_buf = n_buf*2
  n_buf = n_buf * self%nw * self%nv
END FUNCTION total_mesg_buffer_size

!===============================================================================
FUNCTION check_in_bounds (self) result (cib)
  class(ray_patch_t) :: self
  logical :: cib
  !-----------------------------------------------------------------------------
  cib = self%fcheck>=self%fs .and. self%fcheck<=self%fe
  cib = cib .and. self%m3d_patch_t%check_in_bounds()
END FUNCTION check_in_bounds

!===============================================================================
!> Fast method to generate nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(ray_patch_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  integer :: i, key(4)
  !-----------------------------------------------------------------------------
  call trace%begin('ray_patch_t%init_nbors_hash', itimer=itimer)
  !call self%m3d_patch_t%init_nbors_hash(nbors)

  key(2:4) = self%key(2:4)
  key(1) = self%key(1)-1
  call self%check_nbor(key, nbors)

  key(1) = self%key(1)+1
  call self%check_nbor(key, nbors)
  call trace%end (itimer)
END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE bndries (self)
    class(ray_patch_t) :: self
    integer :: i, j, k, iv, iw, b, it, ng
    real(4) :: xxold(3), xxnew(self%ng(3))
    real(4), dimension(:,:,:,:,:,:), pointer :: mem
    ! ----------------------------------------------------------------------------
    if (self%ng(3) == 0) return
    it = self%it
    ng = self%ng(3)
    mem => self%mem

    if (self%boundaries%is_set(bits%zu)) then
        b = self%ui(3)
        mem(:,:,b+1:b+ng,:,:,:) = 0.
    endif

    if (self%boundaries%is_set(bits%zl)) then
        b = self%li(3)
        xxold = [(i*self%ds(3), i=1, 3)]
        xxnew = [(i*self%ds(3), i=1-ng, 0)]
        do iw=1, self%nw
        do iv=1, self%nv
        do j=1,  self%gn(2)
        do i=1,  self%gn(1)
          call pchip_mono(xxold, mem(i,j,b:b+2,   iv,it,iw), &
                          xxnew, mem(i,j,b-ng:b-1,iv,it,iw), sorted=.true.)
          mem(i,j,b-ng:b-1,iv,it,iw) = max(mem(i,j,b-ng:b-1,iv,it,iw), 0.)
        end do
        end do
        end do
        end do
    endif
  END SUBROUTINE bndries

END MODULE ray_patch_mod
