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
  USE interpol_mod, only: pchip_mono
  USE console_mod
  USE m3d_patch_mod, only: m3d_patch_t, &
    converged, print_send, n_nu_out
  USE spectrum_mod, only: spectrum_type
  USE const_mod, only: cc, over_4pi, PI
  implicit none
  private

  !-----------------------------------------------------------------------------

  public :: output_patch_window, output_patch_range

  type, public, extends(m3d_patch_t) :: ray_patch_t
    real(8) :: vec(3), wmu
    integer :: fs, fe, nf
    integer :: pre_pad_fs, pre_pad_fe
    logical :: padded
  contains
    procedure :: check_in_bounds
    procedure :: set_freq_indices
    procedure :: pack
    procedure :: unpack
    procedure :: allocate_mesg
    procedure :: total_mesg_buffer_size
    procedure :: init_bdries
    procedure :: bndries
  end type ray_patch_t

  !-----------------------------------------------------------------------------

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
    if (self%do_out_sp .and. self%i_nu < 0 .and. self%sp%n_segments > 0) then
      call output_patch_range(self%sp, self%i_nu, n_nu_out, fs, fe)
    else
      call unequal_intervals(self%sp%nnu, self%i_nu, self%n_nu, fs, fe)
    endif

    self%fs = fs
    self%fe = fe
    self%pre_pad_fs = fs
    self%pre_pad_fe = fe
    self%nf = self%fe - self%fs + 1

    if (self%do_out_sp) then
      counter => counter2
      label = 'output'
    else
      counter => counter1
      label = 'internal'
    endif
    if (self%i_nu > counter .and. console_on(ch_m3d, 2)) then
      if (counter == -999) write(console_unit,'(I2, A10, A18)') self%n_nu, trim(label), 'frequency splits'
      write(console_unit,'(I4, I10,A3,I8)') self%i_nu, fs, ' - ', fe
      counter = self%i_nu
    endif
  END SUBROUTINE set_freq_indices

  !===============================================================================
  !> Window (segment) index of output patch i_nu (in -n_out..-1) on a segmented
  !> (window_file) grid.  The n_out patches split every window into the same number
  !> n_out/n_segments of near-equal parts, numbered in wavelength order.
  !===============================================================================
  PURE INTEGER FUNCTION output_patch_window (sp, i_nu, n_out) result(i_seg)
    class(spectrum_type), intent(in) :: sp
    integer, intent(in) :: i_nu, n_out
    !-----------------------------------------------------------------------------
    i_seg = (abs(i_nu) - 1) / (n_out / sp%n_segments) + 1
  END FUNCTION output_patch_window

  !===============================================================================
  !> Frequency range [fs, fe] of output patch i_nu on a segmented grid: part
  !> mod(|i_nu|-1, n_split)+1 of window output_patch_window(), so no patch spans
  !> the gap between two windows.
  !===============================================================================
  SUBROUTINE output_patch_range (sp, i_nu, n_out, fs, fe)
    class(spectrum_type), intent(in) :: sp
    integer, intent(in) :: i_nu, n_out
    integer, intent(out) :: fs, fe
    integer :: n_split, i_seg, first
    !-----------------------------------------------------------------------------
    n_split = n_out / sp%n_segments
    i_seg = output_patch_window(sp, i_nu, n_out)
    first = sp%segment_first(i_seg)
    call unequal_intervals(sp%segment_last(i_seg) - first + 1, &
                           mod(abs(i_nu) - 1, n_split) + 1, n_split, fs, fe)
    fs = fs + first - 1
    fe = fe + first - 1
  END SUBROUTINE output_patch_range

!===============================================================================
!> Keep only the physical boundary bit of the face through which the rays
!> enter; set_upstream_boundary reads it to find the start of a sweep.
!===============================================================================
  SUBROUTINE init_bdries (self)
    class(ray_patch_t):: self
    !-----------------------------------------------------------------------------
    call trace%begin ('ray_patch_t%init_bdries')
    call self%m3d_patch_t%init_bdries()
    if (self%dir > 0) then
      call self%boundaries%clear(bits%zu)
    elseif (self%dir < 0) then
      call self%boundaries%clear(bits%zl)
    endif
    call trace%end()
  END SUBROUTINE init_bdries

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
    if ((console_on(ch_m3d, 3, all_ranks=.true.) .or. self%print_update_info) .and. print_send) then
      write(console_unit, '(I6, A32, I3, A6, I10)') self%id, 'packing message on rank', &
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
    if (converged .and. self%do_full_sp) then
      header%level = 0
    endif
    ibuf = 1
    call anonymous_copy (n_header, header, mesg%buffer(ibuf))
    ibuf = ibuf + n_header
    if (.not. converged .or. self%do_out_sp) then
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
    if ((console_on(ch_m3d, 3, all_ranks=.true.) .or. self%print_update_info) .and. print_send) then
      write(console_unit, '(I6, A32, I3, A6, I10)') self%id, 'unpacking message on rank', &
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
      if (console_on(ch_m3d, 2, all_ranks=.true.)) write(console_unit,'(A, I3, A, I3)') &
        'Converged message has been passed by MPI from rank', header%rank, ' to', mpi%rank
    endif
    ibuf = ibuf + n_header
    if (.not. converged .or. self%do_out_sp) then
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
