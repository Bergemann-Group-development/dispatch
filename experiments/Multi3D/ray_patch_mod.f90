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
  USE op_math_mod, only: searchsorted
  USE ray_mod, only: ray_type
  USE m3d_patch_mod, only: m3d_patch_t, vec, wmu, long_vec, long_wmu, detailed_prints, &
  converged, print_send, atom_update_counter, hash_string, comoving_frame
  USE mconst, only: cc
  implicit none
  private

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_patch_t) :: ray_patch_t
    type(ray_type) :: ray
    real(8) :: vec(3), wmu
    integer :: fs, fe, i_out=0
    integer :: pre_pad_fs, pre_pad_fe
    logical :: padded
  contains
    procedure :: pre_init => ray_patch_pre_init
    procedure :: init => ray_patch_init
    procedure :: print_post_update_info
    procedure :: check_in_bounds
    procedure :: finish => ray_patch_finish
    procedure :: nbor_relations => ray_nbor_relations
    procedure :: set_freq_indices
    procedure :: check_topbot
    procedure :: is_upstream => ray_is_upstream
    procedure :: comp_vproj
    procedure :: ff_index
    procedure :: compress_mem
    procedure :: decompress_mem
    procedure :: pack
    procedure :: unpack
    procedure :: allocate_mesg
    procedure :: total_mesg_buffer_size
    procedure :: init_nbors_hash
    procedure :: frequency_padding
  end type ray_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0
  character(len=12), public, save :: formal_solver='linear'
  real(4), public, save :: maxvel

contains

!===============================================================================
SUBROUTINE ray_patch_pre_init (self)
  class(ray_patch_t), target :: self
  integer :: iostat
  logical, save :: first_time=.true.
  integer, save :: itimer=0
  namelist /formal_solver_params/ formal_solver
  ! ----------------------------------------------------------------------------
  call trace%begin('ray_patch_t%pre_init', itimer=itimer)

  if (first_time) then
    rewind (io%input)
    read (io%input, formal_solver_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, formal_solver_params)
    end if
    call io%namelist_warning ('formal_solver_params', iostat)
    write (io%namelist, formal_solver_params)
    first_time = .false.
  endif

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
  elseif (self%i_mu == 0 .and. comoving_frame) then
    self%wmu = math%pi4
  endif

  call self%check_topbot()
  call self%set_freq_indices()
  call trace%end(itimer)
END SUBROUTINE ray_patch_pre_init

!===============================================================================
SUBROUTINE ray_patch_init (self)
  class(ray_patch_t) :: self
  real(8) :: offset
  real(4), dimension(:,:,:), pointer, contiguous :: vproj
  integer :: tmp, fs, fe
  ! ----------------------------------------------------------------------------
  self%padded = .false.
  if (comoving_frame .and. .not. self%out_ray) then
    select case (self%kind)
    case ('atom_line', 'atom_rate')
      ! frequency padding for atom_lines is done in atom_line%init
    !case ('back_cont', 'back_line', 'atom_cont', 'net_opac', 'mean_rad')
    case default
      call self%frequency_padding()
      self%padded = .true.
    end select
  endif
  if (self%is_set(bits%virtual)) then
    offset = self%dir * (self%n(3)-1)/2 * self%ds(3)
    self%position(3) = self%position(3) + offset
    self%origin(3) = self%origin(3) + offset
    self%size(3) = self%ds(3)
    self%n(3) = 1
    self%ng = 0
  endif
  
  call self%m3d_patch_t%init()

  if (self%save_output) then
    if (self%padded) then
      fs = self%pre_pad_fs - self%fs + 1
      fe = fs + self%pre_pad_fe - self%pre_pad_fs
      self%out_mem => self%authority_mem(:,:,:,fs:fe,self%it,:)
    else
      self%out_mem => self%authority_mem(:,:,:,:,self%it,:)
    endif
  endif
END SUBROUTINE ray_patch_init

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
  self%pre_pad_fs = fs
  self%pre_pad_fe = fe
  if (i_nu > counter) then
    if (counter == 0) write(stdout,'(I2, A18)') self%n_nu, 'frequency splits'
    write(stdout,'(I4, I10,A3,I8)') i_nu, fs, ' - ', fe
    counter = i_nu
  endif
END SUBROUTINE set_freq_indices

!===============================================================================
SUBROUTINE comp_vproj (self, vproj)
  class(ray_patch_t) :: self
  real(4), dimension(:,:,:) :: vproj
  real(4), dimension(:,:,:), pointer :: d, px, py, pz
  integer :: li(3), ui(3)
  !-----------------------------------------------------------------------------
  if (self%i_mu == 0) then
    vproj = 0
  else
    li = self%patch%li
    ui = self%patch%ui

    d  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%d, self%it,1)
    px => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%px,self%it,1)
    py => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%py,self%it,1)
    pz => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%pz,self%it,1)

    !if (self%out_ray) then
    !  vproj = (-self%vec(1)*px - self%vec(2)*py + self%vec(3)*pz) * scaling%u / (d*cc)
    !else
      vproj = (self%vec(1)*px + self%vec(2)*py + self%vec(3)*pz) * scaling%u / (d*cc)
    !endif
  endif
END SUBROUTINE comp_vproj

!===============================================================================
SUBROUTINE frequency_padding (self)
  class(ray_patch_t) :: self
  real(8) :: max_nu, min_nu, dnu
  integer :: li(3), ui(3), fs, fe, nv
  !-----------------------------------------------------------------------------
  dnu = maxval(self%sp%nu) * maxvel
  min_nu = self%sp%nu(self%fs) - dnu
  max_nu = self%sp%nu(self%fe) + dnu
  fs = searchsorted(self%sp%nu, min_nu, outbounds_ok=.true.)
  fe = searchsorted(self%sp%nu, max_nu, outbounds_ok=.true.)
  self%fs = max(1, fs)
  self%fe = min(self%sp%nnu, fe+1)
  
  nv = self%fe-self%fs+1
  call self%set_nv_and_nw (nv,self%nw)
END SUBROUTINE frequency_padding

!===============================================================================
SUBROUTINE ray_patch_finish (self)
  class(ray_patch_t) :: self
  !-----------------------------------------------------------------------------
  call self%m3d_patch_t%finish()
  self%static_mem = .false.
END SUBROUTINE ray_patch_finish

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
    call self%set(bits%no_download)
  endif
  if (self%patch%boundaries%lower .and. self%dir >= 0) then
    self%boundaries%lower = .true.
    call self%boundaries%set(bits%zl)
    call self%set(bits%no_download)
  endif
END SUBROUTINE check_topbot

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
        elseif (another%is_upstream(self)) then
          needed = .false.
          download = .false.
        else
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
SUBROUTINE compress_mem (self, compressor)
  class(ray_patch_t):: self
  class(m3d_patch_t), intent(in), target :: compressor
  integer :: i, k, m(6)
  real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: dn
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%compressed) then
    call trace%begin ('ray_patch_t%compress_mem', itimer=itimer)
    call self%mem_lock(1)%set ('compress_mem')
    self%static_mem = .true.
    self%compressed = .true.
    m = shape(self%mem)
    m(3) = 1
    k = merge(self%ui(3), self%li(3),self%dir>0)
    call os%allocate(dn, m)
    do i=1, self%nv
      dn(:,:,1:1,i,:,:) = self%mem(:,:,k:k,i,:,:)
    end do
    call self%deallocate_mem(compressor)
    call self%m3d_patch_t%allocate_mem(compressor, m)
    self%allocation_time = self%time
    do i=1, self%nv
      self%mem(:,:,1:1,i,:,:) = dn(:,:,1:1,i,:,:)
    end do
    call os%deallocate(dn)
    call self%mem_lock(1)%unset ('compress_mem')
    if (self%print_update_info .and. detailed_prints) then
      write(stdout,*) self%m3d_patch_info2(), ' has been compressed'
    endif
    call trace%end(itimer)
  endif
END SUBROUTINE compress_mem

!===============================================================================
SUBROUTINE decompress_mem (self, decompressor)
  class(ray_patch_t) :: self
  class(m3d_patch_t), intent(in), target :: decompressor
  integer :: i, k, m(6)
  real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: up
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%compressed) then
    call trace%begin ('ray_patch_t%decompress_mem', itimer=itimer)
    call self%mem_lock(1)%set ('decompress_mem')
    self%static_mem = .false.
    self%compressed = .false.
    m = shape(self%mem)
    call os%allocate(up, m)
    do i=1, self%nv
      up(:,:,1:1,i,:,:) = self%mem(:,:,1:1,i,:,:)
    end do
    call self%m3d_patch_t%deallocate_mem(decompressor)
    call self%allocate_mem(decompressor)
    k = merge(1,self%gn(3),self%dir>0)
    do i=1, self%nv
      self%mem(:,:,k:k,i,:,:) = up(:,:,1:1,i,:,:)
    end do
    call os%deallocate(up)
    call self%mem_lock(1)%unset ('decompress_mem')
    if (self%print_update_info .and. detailed_prints) then
      write(stdout,*) self%m3d_patch_info2(), ' has been decompressed'
    endif
    call trace%end(itimer)
  endif
END SUBROUTINE decompress_mem
  
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
    write(stdout, '(I6, A32, I3, A6, I6)') self%id, 'packing message on rank', & 
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
  if (converged .and. .not. self%out_ray) then
    header%level = 0
  endif
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  if (.not. converged .or. self%out_ray) then
    call self%check_mem()
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
    write(stdout, '(I6, A32, I3, A6, I6)') self%id, 'unpacking message on rank', & 
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
    call self%check_mem()
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
SUBROUTINE print_post_update_info (self)
  class(ray_patch_t) :: self
  integer :: f, p(3)
  !-----------------------------------------------------------------------------
  f = self%fcheck
  p = self%pcheck
  call self%print_to_stdout(arr=self%mem(p(1),p(2),p(3),f,1,:))
END SUBROUTINE print_post_update_info

!===============================================================================
FUNCTION check_in_bounds (self) result (cib)
  class(ray_patch_t) :: self
  logical :: cib
  integer :: f
  !-----------------------------------------------------------------------------
  cib = .false.
  call self%ff_index(self%fcheck, f)
  if (f > 0) then
    if (self%kind == 'out_rad') then
      self%pcheck(3) = 1
    endif
    cib = self%m3d_patch_t%check_in_bounds()
  endif
  self%fcheck = f
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

END MODULE ray_patch_mod
