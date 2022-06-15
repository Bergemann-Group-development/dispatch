module ray_patch_mod

  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE mpi_mesg_mod
  USE halo_mod
  USE gpatch_mod
  USE link_mod
  USE task_mod
  !USE index_mod
  !USE mesh_mod
  USE io_mod
  !USE aux_mod
  USE timer_mod
  USE scaling_mod
  USE trace_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE shared_mod
  !USE math_mod
  USE bits_mod
  USE os_mod
  !USE units_mod
  !USE scratch_mod
  USE ray_mod, only:ray_type
  USE m3d_patch_mod, only: m3d_patch_t, vec, wmu, long_vec, long_wmu
  USE mconst, only: cc
  implicit none
  private

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_patch_t) :: ray_patch_t
    type(ray_type) :: ray
    real(8) :: vec(3), wmu
    integer :: dir=0, fs, fe, i_out=0
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
  contains
    procedure :: pre_init => ray_patch_pre_init
    procedure :: init => ray_patch_init
    procedure :: post_update => ray_post_update
    procedure :: nbor_relations => ray_nbor_relations
    procedure :: set_freq_indices
    procedure :: check_topbot
    procedure :: is_upstream => ray_is_upstream
    procedure :: comp_vproj
    procedure :: deallocate_mem => ray_patch_deallocate_mem
    procedure :: ff_index
    procedure :: pack => ray_pack
    procedure :: unpack => ray_unpack
    procedure :: allocate_mesg => ray_allocate_mesg
    procedure :: ray_header_to_patch
    procedure :: ray_patch_to_header
    procedure :: total_mesg_buffer_size
  end type ray_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0
  !-----------------------------------------------------------------------------
  ! Sequential header, for MPI and I/O
  !-----------------------------------------------------------------------------
  integer, public, save:: n_header = (14*4 + 4*8 + 5*8*3 + 1*8 + 2*8*8 + 1*4*8 + 2*4*3)/4
  type, public:: ray_header_t
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
    integer:: nv                        ! number of frequencies
    integer:: nw                        ! number of variables
    integer:: it
    integer:: nq=0
    integer:: rank_steps=0
    real(8):: time
    real(8):: out_next
    real(8):: print_next
    real(8):: dtime
    real(8):: position(3)
    real(8):: velocity(3)
    real(8):: size(3)
    real(8):: box(3)                    ! the periodic wrapping size
    real(8):: ds(3)
    real(8):: send_time
    real(8), dimension(8):: t
    real(8), dimension(8):: dt
    integer, dimension(8):: iit
    integer:: n(3)=16                   ! convenience
    integer:: ng(3)=2                   ! convenience
  end type
    
contains

!===============================================================================
SUBROUTINE ray_patch_pre_init (self)
  class(ray_patch_t), target :: self
  ! ----------------------------------------------------------------------------
  call self%m3d_patch_t%pre_init()
  
  if (self%out_ray .and. self%i_mu /= 0) then
    self%i_out = self%i_mu - self%n_omega 
    self%vec = long_vec(:,self%i_out) 
    self%wmu = long_wmu(self%i_out)
    self%dir = 1
  elseif(self%i_mu > 0) then
    self%vec = vec(:,self%i_mu)
    self%wmu = wmu(self%i_mu)
    self%dir = merge(1, -1, vec(3,self%i_mu) > 0)
  endif

  call self%check_topbot()

  !if (self%dir < 0) then
  !  self%boundaries%lower = .false.
  !elseif (self%dir > 0) then
  !  self%boundaries%upper = .false.
  !endif

  !select case (self%kind)
  !case ('line_opac', 'cont_opac')
  !  if (self%boundaries%lower .or. self%boundaries%upper) then
  !    self%time = self%time - 0.1
  !  endif
  !end select
  
  call self%set_freq_indices()
END SUBROUTINE ray_patch_pre_init

!===============================================================================
SUBROUTINE set_freq_indices (self)
  class(ray_patch_t) :: self
  integer :: fs, fe, nv, md, i_nu
  integer, save :: counter=0
  !-----------------------------------------------------------------------------
  if (self%kind == 'out_rad') then
    self%fs = 1
    self%fe = self%sp%nnu
    return
  endif
  nv = self%sp%nnu / self%n_nu
  md = mod(self%sp%nnu, self%n_nu)

  i_nu = self%i_nu
  if (i_nu > self%n_nu) then
    i_nu = i_nu - self%n_nu
  endif

  fs = (i_nu-1) * nv + 1 + min(md, i_nu-1)

  if (i_nu <= md) then
    nv = nv + 1
  endif

  fe = fs + nv - 1

  self%fs = fs
  self%fe = fe

  if (i_nu > counter) then
    if (counter == 0) write(stdout,'(I2, A18)') self%n_nu, 'frequency splits'
    write(stdout,'(I4, I10,A3,I8)') i_nu, fs, ' - ', fe
    counter = i_nu
  endif
END SUBROUTINE set_freq_indices

!===============================================================================
SUBROUTINE comp_vproj (self, vproj)
  class(ray_patch_t) :: self
  real(4), dimension(:,:,:), pointer, contiguous :: vproj
  real(4), dimension(:,:,:), pointer :: d, px, py, pz
  integer :: li(3), ui(3)
  !-----------------------------------------------------------------------------
  li = self%patch%li
  ui = self%patch%ui

  d  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%d, self%it,1)
  px => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%px,self%it,1)
  py => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%py,self%it,1)
  pz => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%pz,self%it,1)

  if (self%out_ray) then
    vproj = (-self%vec(1)*px - self%vec(2)*py + self%vec(3)*pz) * scaling%u / (d*cc)
  else
    vproj = (self%vec(1)*px + self%vec(2)*py + self%vec(3)*pz) * scaling%u / (d*cc)
  endif
END SUBROUTINE comp_vproj

!===============================================================================
SUBROUTINE ray_post_update (self)
  class(ray_patch_t) :: self
  logical :: in_bounds
  integer :: f, p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_post_update')
  if (self%verbose > 1) then
    p = self%pcheck
    call self%ff_index(self%fcheck, f)
    in_bounds = f > 0

    if (in_bounds) then
      if (self%kind /= 'out_rad') in_bounds = all(p > 0)
    endif

    if (in_bounds) then
      select case (self%kind)
      case('rate_patch')
        write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), &
        self%mem(p(1),p(2),p(3),self%linecheck,1,:)
      case ('out_rad')
        write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), &
        sum(self%mem(:,:,1,f,1,1)) / product(self%n(1:2))
      case default
        write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), &
        self%mem(p(1),p(2),p(3),f,1,:)
      end select
    elseif (self%verbose > 2) then
      write(stdout,'(A,A)') self%m3d_patch_info()
    endif
  endif 

  call self%m3d_patch_t%post_update()

  call trace%end()
END SUBROUTINE ray_post_update

!===============================================================================
SUBROUTINE ray_patch_init (self)
  class(ray_patch_t) :: self
  real(8) :: offset
  ! ----------------------------------------------------------------------------
  if (self%is_set(bits%virtual)) then
    offset = self%dir * (self%n(3)-1)/2 * self%ds(3)
    self%position(3) = self%position(3) + offset
    self%origin(3) = self%origin(3) + offset
    self%size(3) = self%ds(3)
    self%n(3) = 1
    self%ng(3) = 0
  elseif (self%i_mu > 0) then
    call os%allocate(self%vproj, [self%n(1), self%n(2), self%n(3)])
  endif
  call self%m3d_patch_t%init()
END SUBROUTINE ray_patch_init

!===============================================================================
SUBROUTINE ray_patch_deallocate_mem (self)
  class(ray_patch_t):: self
  integer :: m(6), li(3), ui(3)
  !-----------------------------------------------------------------------------
  call self%m3d_patch_t%deallocate_mem()
  if (associated(self%vproj)) then
    call os%deallocate(self%vproj)
  endif
END SUBROUTINE ray_patch_deallocate_mem

!===============================================================================
!> Determine if a task is upstream.  
!===============================================================================
logical FUNCTION ray_is_upstream (self, other) result(is_upstream)
  class(ray_patch_t):: self, other
  real(8):: ds1(3), ds2(3)
  real:: sinth
  real, parameter:: eps=1e-4
  integer:: is1(3), is2(3), is3(3), is4(3), i, lbit(3), ubit(3)
  logical:: lower, upper, upstream
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and other points inward, or the task has no RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  if (.not.self%on) then
    is_upstream=.false.
    return
  end if
  !
  if (self%i_mu /= other%i_mu) then
    is_upstream=.false.
    return
  end if
  !
  call trace%begin ('ray_patch_t%ray_is_upstream')
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.other%on) then
    is_upstream = .false.
    call trace%end ()
    return
  end if
  !
  is_upstream=.true.
  ds1 = self%position-other%position
  where (self%periodic)
    ds1 = modulo(ds1 + 0.5_8*self%box, self%box) - 0.5_8*self%box
  end where
  ds1 = ds1/self%box                                                            ! relative task offsets
  is1 = 0; is2=0
  do i=1,3
    if (ds1(i)      > +eps) is1(i)=1
    if (self%vec(i) > +eps) is2(i)=1
    if (ds1(i)      < -eps) is1(i)=-1
    if (self%vec(i) < -eps) is2(i)=-1
  end do
  if (all(is1*is2 <= 0)) is_upstream = .false.
  ! ----------------------------------------------------------------------------
  ! If the nbor lies outside the boundaries it is not needed
  ! ----------------------------------------------------------------------------
  lbit = [bits%xl, bits%yl, bits%zl]
  ubit = [bits%xu, bits%yu, bits%zu]
  do i=1,3
    lower = self%boundaries%is_set(lbit(i))
    upper = self%boundaries%is_set(ubit(i))
    if (upper .and. ds1(i)==1) then
      is_upstream = .false.
    end if
    if (lower .and. ds1(i)==-1) then
      is_upstream = .false.
    end if
  end do
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.other%on) then
    is_upstream = .false.
    call trace%end ()
    return
  end if
  call trace%end ()
END FUNCTION ray_is_upstream

!===============================================================================
SUBROUTINE ff_index(self, ff, f)
  class(ray_patch_t):: self
  integer, intent(in)  :: ff
  integer, intent(out) :: f
  logical :: ff_in_range
  !-----------------------------------------------------------------------------
  ff_in_range = ff>=self%fs .and. ff<=self%fe
  if (ff_in_range) then
    f = ff - self%fs + 1
  else
    f = -1
  endif
END SUBROUTINE ff_index

!===============================================================================
SUBROUTINE check_topbot (self)
  class(ray_patch_t):: self
  !.............................................................................
  self%boundaries%upper = .false.
  self%boundaries%lower = .false.
  
  if (self%patch%boundaries%upper .and. self%dir <= 0) then
    self%boundaries%upper = .true.
    call self%boundaries%set(bits%zu)
  elseif (self%ipos(3)==0 .and. self%dir >= 0) then
    self%boundaries%lower = .true.
    call self%boundaries%set(bits%zl)
  endif
END SUBROUTINE check_topbot

!===============================================================================
SUBROUTINE ray_nbor_relations(self, another, needed, needs_me, download)
  class(ray_patch_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  class(task_t), pointer :: other
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_nbor_relations')
  call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

  if (self%kind == 'net_opac' .or. self%kind == 'angle_rad') then
    call self%set(bits%boundary)
  endif

  if (needed .or. needs_me) then
    select type(another)
    class is(ray_patch_t)
      if (self%kind==another%kind) then
        if (self%is_upstream(another)) then
          needs_me = .false.
          if (self%is_set(bits%virtual)) then
            call self%set(bits%external)
          endif
        elseif (another%is_upstream(self)) then
          needed = .false.
          download = .false.
        else
          needed = .false.; needs_me = .false.; download = .false.
        endif
      endif
    class default
      if (another%id /= self%m3d_parent%id) then
        needed = .false.; needs_me = .false.; download = .false.
      endif
    end select
  endif

  call trace%end()
END SUBROUTINE ray_nbor_relations

!===============================================================================
!> Pack patch info into a buffer, for MPI communication
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE ray_pack (self)
  class(ray_patch_t):: self
  !.............................................................................
  class(mesg_t), pointer:: mesg
  type(ray_header_t):: header
  integer:: n_buf, ibuf, iv, it, iw, iz
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_pack', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate a task%mesg instance, mesg%reqs() and similar array
  !-----------------------------------------------------------------------------
  call self%allocate_mesg
  if (omp%nthreads > 1) then
    allocate (self%mesg%lock)
    call self%mesg%lock%init ('smsg')
  end if
  mesg => self%mesg
  allocate (mesg%reqs(mpi%size))
  allocate (mesg%fail(mpi%size))
  allocate (mesg%ranks(mpi%size))
  mesg%fail = 0
  mesg%reqs = 0
  mesg%ranks = 0
  call io%bits_mem (storage_size(mesg%reqs), &
                    product(shape(mesg%reqs)))
  !-----------------------------------------------------------------------------
  ! Packing meta-data is a task-owner read operation, and does not need locking
  !-----------------------------------------------------------------------------
  call self%ray_patch_to_header (header)
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer.
  !-----------------------------------------------------------------------------
  iz = merge(self%ui(3), self%li(3), self%dir > 0)
  n_buf = self%total_mesg_buffer_size() / (self%nw*self%nv)
  do iw=1,self%nw
    do iv=1,self%nv
      call anonymous_copy (n_buf, self%mem(:,:,iz,iv,self%it,iw), mesg%buffer(ibuf))
      ibuf = ibuf + n_buf
    end do
  end do
  !-----------------------------------------------------------------------------
  if (self%verbose > 2) &
    print *, self%id, 'packing message on rank', mpi%rank
  call trace%end (itimer)
END SUBROUTINE ray_pack

!===============================================================================
!> Unpack patch info from an MPI buffer.  To ensure consistency, the task meta-
!> data should not be updated until after the time slot data are in place.
!> Any task that has this task as virtual nbor and is waiting for the time slot
!> will thus not be scheduled until both the meta-data and the time slot data
!> are up-to-date.  This also reduces the chance of lock contention with respect
!> to the new time slot.
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE ray_unpack (self, mesg)
  class(ray_patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(ray_header_t):: header
  type(lock_t), pointer:: lock
  integer:: n_buf, ibuf, iv, it, i, iw, iz
  integer, dimension(:), pointer:: buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_unpack', itimer=itimer)
  call self%set (bits%busy)
  self%unpack_time = wallclock()
  buffer => mesg%buffer
  n_buf = self%total_mesg_buffer_size() / (self%nw*self%nv)
  !-----------------------------------------------------------------------------
  ! Copy buffer into sequenced header, to get the new %it slot
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  ibuf = ibuf + n_header
  call self%ray_header_to_patch (header)
  call self%check_mem()
  !-----------------------------------------------------------------------------
  ! Copy array data from the buffer to the time slot, under time slot lock
  !-----------------------------------------------------------------------------
  lock => self%mem_lock(header%it)
  !iz = merge(self%ui(3), self%li(3), self%dir > 0)
  iz=1
                                                        call lock%set ('unpack')
  do iw=1,self%nw
    do iv=1,self%nv
      call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,iz,iv,header%it,iw))
      ibuf = ibuf + n_buf
    end do
  end do
                                                      call lock%unset ('unpack')
  !-----------------------------------------------------------------------------
  ! Update the patch meta-data under lock, now that the time slices are in place
  !-----------------------------------------------------------------------------
                                                   call self%lock%set ('unpack')
  self%it         = header%it
  self%time       = header%time
  self%nt         = header%nt
  if (header%nt > 0) then
    self%new      = mod(header%it,header%nt) + 1
    do i=1,self%nt
      self%t(i)   = header%t(i)
      self%dt(i)  = header%dt(i)
      self%iit(i) = header%iit(i)
    end do
  else
    write (io_unit%log,*) omp%thread, &
      'header_to_patch WARNING: nt non-positive', header%nt
    flush (io_unit%log)
  end if
  !.............................................................................
                                               call self%lock%unset ('unpack')
  call self%clear (bits%busy)
  self%wc_last = wallclock()
  !-----------------------------------------------------------------------------
  ! Diagnostics
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,'(f12.6,3x,a,i7,2i4,2x,3l1,2i4,f8.4,f11.6)') &
      wallclock(), 'patch_t%unpack: id, nq, BVS, nsnd, nrcv, latency, time =', &
      self%id, mesg%sender, self%nq, &
      self%is_set(bits%boundary), &
      self%is_set(bits%virtual), &
      self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, self%latency, self%time
    write (io_unit%log,1) &
      'patch_t%unpack id, it, nq, time  after:', &
      self%id, self%it, self%nq, self%time, &
      self%is_set (bits%internal), &
      self%is_set (bits%boundary), &
      self%is_set (bits%virtual), &
      self%is_set (bits%external), &
      self%is_set (bits%swap_request), &
      self%wc_last
    1 format(a,3i6,f12.6,2x,5l1,f12.1)
  end if
  if (self%verbose > 2) then
    print *, self%id, 'unpacking message on rank', mpi%rank
  end if
  !-----------------------------------------------------------------------------
    call trace%end (itimer)
END SUBROUTINE ray_unpack

!===============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE ray_patch_to_header (self, head)
  class(ray_patch_t):: self
  type(ray_header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  head%n_header   = n_header
  head%n_data     = product(self%gn)*self%nv
  head%id         = self%id
  head%status     = self%status
  head%rank       = self%rank
  head%istep      = self%istep
  head%iout       = self%iout
  head%level      = self%level
  head%time       = self%time
  head%out_next   = self%out_next
  head%print_next = self%print_next
  head%dtime      = self%dtime
  head%position   = self%position
  head%velocity   = self%velocity
  head%size       = self%size
  head%nt         = self%nt
  head%nv         = self%nv
  head%nw         = self%nw
  head%it         = self%it
  nt              = self%nt
  head%t(1:nt)    = self%t
  head%dt(1:nt)   = self%dt
  head%iit(1:nt)  = self%iit
  head%nq         = shared%n_queue
  head%rank_steps = self%rank_steps
  head%send_time  = wallclock()
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,i4,l3,2i4,f8.4)') &
      wallclock(), 'ray_patch_to_header: id, nq =', self%id, &
      self%nq, self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n
END SUBROUTINE ray_patch_to_header

!===============================================================================
!> Copy he most important fields from a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE ray_header_to_patch (self, head)
  class(ray_patch_t):: self
  type(ray_header_t):: head
  integer:: nt, new, i, save_bits
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming patches.  Also clear the
  ! ready bit, which is an indicator of a patch being in the ready queue, which
  ! an incoming patch is not.  By doing this here, we avoid having to do it in
  ! a critical region, otherwise needed to protect an operation on the patch
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%ray_header_to_patch', itimer=itimer)
  if (iand(head%status, bits%boundary) /= 0) then
    head%status = ior(head%status, bits%virtual)
    head%status = iand(head%status, not(bits%boundary))
    if (io%verbose>2) &
      write (io_unit%mpi, *) head%id, 'swapping to virtual'
  else if (iand(head%status, bits%virtual) /= 0) then
    head%status = ior(head%status, bits%boundary)
    head%status = iand(head%status, not(bits%virtual))
    if (io%verbose>2) &
      write (io_unit%mpi, *) head%id, 'swapping to boundary'
  end if
  head%status = iand(head%status, not(bits%ready+bits%busy))
  !-----------------------------------------------------------------------------
  ! Lock task while modifying metadata
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Except for the first and last appearance, istep should increase be one
  !-----------------------------------------------------------------------------
  if (head%id /= self%id) then
    call io%abort('wrong message ID unpacked')
  end if
  if ((self%istep > 0) .and. (iand(head%status, bits%remove) == 0)) then
    if (head%istep > self%istep+1) then
      write(io_unit%log,1) &
        'WARNING: early arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
    1 format(a,4i6,2f12.6)
      flush (io_unit%log)
    else if (head%istep < self%istep+1) then
      write(io_unit%log,1) &
        'WARNING: late arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
      flush (io_unit%log)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If the task is virtual, preserve bits%init_nbors and bits%update_status
  ! across the task update
  !-----------------------------------------------------------------------------
                                        call self%lock%set ('header_to_patch')
  if (self%is_set(bits%virtual)) then
    save_bits = iand(self%status,bits%init_nbors+bits%update_status)
    self%status = ior(self%status,save_bits)
  else
    self%status     = head%status
  end if
  self%id         = head%id
  self%rank       = head%rank
  self%istep      = head%istep
  self%iout       = head%iout
  self%level      = head%level
  self%out_next   = head%out_next
  self%print_next = head%print_next
  self%dtime      = head%dtime
  !self%position   = head%position
  self%velocity   = head%velocity
  !self%size       = head%size
  self%rank_steps = head%rank_steps
  self%it         = head%it
  self%nv         = head%nv
  self%nw         = head%nw
  self%nt         = head%nt
  self%nq         = head%nq
  self%t          = head%t  (1:self%nt)
  self%dt         = head%dt (1:self%nt)
  self%iit        = head%iit(1:self%nt)
  self%time       = head%time
  !
  if (head%istep == 0) then
    write (io_unit%log,'(1x,a,i6,f12.6,2x,9f12.6)') &
      'patch_t%header_to_patch: id, time, t =', &
      self%id, self%time, self%t
    write (io_unit%log,'(1x,a,i6,2x,3i4,2x,9i4)') &
      'patch_t%header_to_patch: id, istep, it, new, iit =', &
      self%id, self%istep, self%it, self%new, self%iit
  end if
  if (self%it > self%nt .or. self%it < 1) then
    write (io_unit%log,*) 'header_to_patch ERROR: it =', self%it, self%iit(nt-1)
    self%it = self%iit(nt-1)
  end if
  self%nq         = head%nq
  self%latency    = wallclock() - head%send_time
                                        call self%lock%unset ('header_to_patch')
  !-----------------------------------------------------------------------------
  ! Latency stats
  !-----------------------------------------------------------------------------
  !$omp atomic
  timer%latency%max  = max(timer%latency%max,self%latency)
  !$omp atomic
  timer%latency%aver = timer%latency%aver+self%latency
  !$omp atomic
  timer%latency%n    = timer%latency%n+1
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE ray_header_to_patch

!===============================================================================
!> Allocate a message structure, for communication
!===============================================================================
SUBROUTINE ray_allocate_mesg (self)
  class(ray_patch_t):: self
  class(mesg_t), pointer:: mesg
  integer:: n_buf
  !-----------------------------------------------------------------------------
  call trace%begin ('ray_patch_t%allocate_mesg')
  allocate (self%mesg)
  mesg => self%mesg
  mesg%id = self%id
  call mesg%init
  !-----------------------------------------------------------------------------
  ! Store the standard buffer size (not valid for swap requests!)
  !-----------------------------------------------------------------------------
  n_buf = self%total_mesg_buffer_size()
  mesg%nbuf = n_header + n_buf
  allocate (mesg%buffer(mesg%nbuf))
  call io%bits_mem (storage_size(self%mesg%buffer), &
                    product(shape(self%mesg%buffer)))
  call trace%end ()
END SUBROUTINE ray_allocate_mesg

!===============================================================================
FUNCTION total_mesg_buffer_size(self) result(n_buf)
  class(ray_patch_t):: self
  integer :: n_buf
  !-----------------------------------------------------------------------------
  n_buf = product(self%n(:2) + 2*self%ng(:2))
  if (kind(self%mem)==8) n_buf = n_buf*2
  n_buf = n_buf * self%nw * self%nv
END FUNCTION total_mesg_buffer_size

end module ray_patch_mod
