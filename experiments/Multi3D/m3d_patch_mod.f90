module m3d_patch_mod

  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE halo_mod
  USE gpatch_mod
  USE solver_mod
  USE link_mod
  USE task_mod
  !USE index_mod
  !USE mesh_mod
  USE io_mod
  !USE aux_mod
  !USE timer_mod
  !USE scaling_mod
  USE trace_mod
  !USE mpi_mod
  !USE omp_mod
  USE omp_lock_mod
  USE shared_mod
  !USE omp_timer_mod
  !USE math_mod
  USE bits_mod
  !USE units_mod
  !USE scratch_mod
  USE atom_var, only: atom_type
  USE composition_mod, only:composition_type
  USE spectrum_mod, only: spectrum_type
  implicit none
  private

  !-----------------------------------------------------------------------------

  type, public, extends(gpatch_t) :: m3d_t
    class(gpatch_t), pointer :: patch            ! pointer to MHD patch
    integer                  :: master_id
    logical                  :: on=.false.       ! allows local turn-off
    logical                  :: in_task_list = .false.

    type(spectrum_type) :: sp

    real(8), dimension(:,:,:), allocatable :: temp
    real(8), dimension(:,:,:), allocatable :: rho 
    real(8), dimension(:,:,:), allocatable :: ne  
    real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
    real(8) :: vec(3), wmu, t_offset
    integer :: i_nu=0, dir=0, i_omega=0, fs, fe, iter=0, maxiter
  contains
    procedure :: pre_init => m3d_pre_init
    procedure :: init => m3d_init
    procedure :: post_update => m3d_post_update
    procedure :: finish => m3d_finish
    procedure :: nbor_relations => m3d_nbor_relations
    procedure :: bndries
    procedure :: m3d_name
    procedure :: init_level
    procedure :: task_info
    procedure :: set_freq_indices
    procedure :: append_to_task_list
    procedure :: extrapolate_guards
    procedure :: check_topbot
    procedure :: is_ahead_of
    procedure :: is_upstream => m3d_is_upstream
    procedure :: periodic_grid => m3d_periodic_grid
  end type m3d_t

  !-----------------------------------------------------------------------------

  type, public, extends(halo_t):: halo_m3d_t
  contains
    procedure:: info
  end type

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0                           ! verbosity
    
contains

!===============================================================================
SUBROUTINE m3d_pre_init (self)
  class(m3d_t), target:: self
  !-----------------------------------------------------------------------------
  logical, save :: first_time=.true.
  logical, save :: on=.true.
  integer :: iostat
  !-----------------------------------------------------------------------------
  namelist /m3d_params/ on, verbose

  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_t%pre_init')

  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, m3d_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, m3d_params)
    end if
    call io%namelist_warning ('m3d_params', iostat)
  end if

  if (self%pre_initialized) print *, 'pre_init issued somewhere'
  call io%lock%unset ('m3d_t%pre_init')

  call self%gpatch_t%pre_init()

  ! ----------------------------------------------------------------------------
  ! Fill namelist variables
  ! ----------------------------------------------------------------------------
  self%on = on
  self%in_task_list = .false.

  ! ----------------------------------------------------------------------------
  ! relevant synchronization items for cartesian setup
  ! ----------------------------------------------------------------------------
  self%n          = self%patch%n
  self%ng         = 0 !self%patch%ng
  self%ipos       = self%patch%ipos
  self%rank       = self%patch%rank
  self%box        = self%patch%box
  self%size       = self%patch%size
  self%position   = self%patch%position
  self%ds         = self%patch%ds
  self%origin     = self%patch%origin
  self%llc_cart   = self%patch%llc_cart
  self%centre_nat = self%patch%centre_nat
  self%llc_nat    = self%patch%llc_nat
  self%status     = self%patch%status
  self%periodic   = self%patch%periodic
  self%level      = 1

  if (self%on) then 
    call self%task_t%init_unique (same=.true.) ! obtain an unique ID
  end if
  call trace%end ()
END SUBROUTINE m3d_pre_init

!===============================================================================
SUBROUTINE m3d_init (self)
  class(m3d_t) :: self
  integer, save:: n_warmup=1
  ! ----------------------------------------------------------------------------
  call self%gpatch_t%init()
  call self%check_topbot()
  self%mem = tiny(1.)
  self%it = 1
  self%new = self%it
  self%authority_mem => self%mem(self%li(1):self%ui(1),self%li(2):self%ui(2),self%li(3):self%ui(3),:,:,:)
  !-----------------------------------------------------------------------------
  ! Synchronize with MHD patch
  ! ----------------------------------------------------------------------------
  self%dt_fixed = self%patch%dtime
  self%dtime    = self%patch%dtime
  self%dt       = self%patch%dtime
  self%out_next = self%patch%out_next

  ! ----------------------------------------------------------------------------
  ! Do not check density or validate with it. This is rather important to be set
  ! ----------------------------------------------------------------------------
  call self%set (bits%no_density)
  call self%set (bits%internal)
  
  ! ----------------------------------------------------------------------------
  ! Allocate and set unsigned and pervolume. Opacity, source function,
  ! intensity and populations are all positive definite
  ! Although quantities are generally per volume, but density is not stored
  ! locally, thus this flag must not be set.
  ! ----------------------------------------------------------------------------
  allocate (self%unsigned(self%nv), self%pervolume(self%nv))
  self%pervolume = .false.
  self%unsigned  = .true.
END SUBROUTINE m3d_init

!===============================================================================
SUBROUTINE m3d_post_update (self)
  class(m3d_t) :: self
  !-----------------------------------------------------------------------------
  if (self%verbose > 1) then
    write(*,'(F5.1, L2, A12, 3(I4), 6(ES15.7))') self%time, self%finished, trim(self%kind), self%patch%id, self%i_nu, self%i_omega, &
    self%authority_mem(1,1,1,1,1,:), self%authority_mem(1,1,self%n(3),1,1,:)
  endif

  if (self%iter == self%maxiter) then
    call self%finish()
  endif

  self%iter = self%iter + 1
END SUBROUTINE m3d_post_update

!===============================================================================
SUBROUTINE m3d_finish (self)
  class(m3d_t):: self
  !.............................................................................
  self%finished = .true.
  self%time = io%end_time + self%t_offset*self%dtime
  call self%set(bits%frozen)
END SUBROUTINE m3d_finish

!===============================================================================
SUBROUTINE set_freq_indices (self, n, nnu, n_nu)
  class(m3d_t) :: self
  integer, intent(in) :: n    ! index  of frequency bin
  integer, intent(in) :: nnu  ! number of total frequencies
  integer, intent(in) :: n_nu ! number of frequency bins
  integer :: fs, fe, nv
  !-----------------------------------------------------------------------------
  nv = nnu / n_nu
  fs = (n-1) * nv + 1
  if (n == n_nu) nv = nv + mod(nnu, n_nu)
  fe = fs + nv - 1

  if (self%kind /= 'rad_rate') then
    self%nv = nv
  endif

  self%i_nu = n
  self%fs = fs
  self%fe = fe
END SUBROUTINE set_freq_indices

!===============================================================================
SUBROUTINE m3d_nbor_relations(self, another, needed, needs_me, download)
  class(m3d_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  class(task_t), pointer :: other
  
  call trace%begin ('m3d_t%m3d_nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  if (another%id == self%patch%id) then
    needed = .true. ; needs_me = .false.; download = .false.
  else
    needed = .false.; needs_me = .false.; download = .false.
  endif

  call trace%end()
END SUBROUTINE m3d_nbor_relations

!===============================================================================
SUBROUTINE append_to_task_list(self)
  class (m3d_t), target :: self
  class(link_t),  pointer :: link
  class (task_t), pointer :: task
  !-----------------------------------------------------------------------------
  allocate (halo_m3d_t:: link)
  link%task => self
  self%link => link
  call self%patch%task_list%append (link)
END SUBROUTINE append_to_task_list

!===============================================================================
SUBROUTINE info (self, label)
  class(halo_m3d_t):: self
  character(*), optional:: label
  class(nbor_t), pointer:: nbor
  real(8):: distance(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_t%info')
  if (verbose < 2) then
    call self%task%task_info()
    nbor => self%nbor
    do while (associated(nbor))
      distance = nbor%task%distance(self%task)
      select type(nt => nbor%task)
      class is(m3d_t)        
        write(io_unit%mpi,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3, a10, i4)') &
          nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
          nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
          nbor%download, nt%is_set(bits%boundary), &
          nt%is_set(bits%virtual), distance/ &
          (0.5d0*(self%task%size+nt%size)), ' i_omega: ', nt%i_omega
      class default
        write(io_unit%mpi,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3)') &
          nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
          nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
          nbor%download, nt%is_set(bits%boundary), &
          nt%is_set(bits%virtual), distance/ &
          (0.5d0*(self%task%size+nt%size))
      end select
      nbor => nbor%next
    end do
  else
    call self%halo_t%info(label)
  end if
  call trace%end ()
END SUBROUTINE info

!===============================================================================
!> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
!===============================================================================
SUBROUTINE task_info (self)
  class(m3d_t):: self
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(1x,a,2x,a10,i6,2i4,7x,2l1,2x,3f9.3,i4)') &
    'm3d_t%task_info: kind, id, rank, level, BV, pos, i_omega =', &
    self%kind, self%id, self%rank, self%level, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), self%position, &
    self%i_omega
END SUBROUTINE task_info

!===============================================================================
!> Determine if a task is is upstream.  
!===============================================================================
logical FUNCTION m3d_is_upstream (self, other) result(is_upstream)
  class(m3d_t):: self, other
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
  if (self%i_omega /= other%i_omega) then
    is_upstream=.false.
    return
  end if
  !
  call trace%begin ('m3d_t%m3d_is_upstream')
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
END FUNCTION m3d_is_upstream

!===============================================================================
SUBROUTINE init_level (self)
  class(m3d_t):: self
  !-----------------------------------------------------------------------------
  if (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  end if
END SUBROUTINE init_level

!===============================================================================
SUBROUTINE bndries (self)
  class(m3d_t):: self
  real(4), dimension(:,:,:,:,:), pointer :: mem
  !----------------------------------------------------------------------------
  call trace_begin('m3d_t%bndries')
  mem => self%mem(:,:,:,:,self%it,:)
  if (size(self%boundaries%topbot) /= 2) then
    print *, 'boundary conditions not initialised for task ', self%id 
  elseif (any(self%boundaries%topbot)) then
    call self%boundaries%boundary(mem)
  endif
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE m3d_periodic_grid (self, only)
  class(m3d_t):: self
  integer, optional:: only
  integer, dimension(3) :: li, ui, lb, ub
  real(4), dimension(:,:,:,:,:), pointer :: mem
  integer :: z, dir, off, is, ie, v, w
  !----------------------------------------------------------------------------
  if (any(self%ng(:2) > 0)) then
    li = self%mesh%li; ui = self%mesh%ui
    lb = self%mesh%lb; ub = self%mesh%ub

    mem => self%mem(:,:,:,:,self%it,:)

    if (all(self%n(:2) == 1)) then
      do w=1, self%nw
        do v=1, self%nv
          do z=1, self%gn(3)
            mem(:,:,z,v,w) = mem(li(1),li(2),z,v,w)
          end do
        end do 
      end do
    else
      do dir=-1,1,2
        is = dir*max(dir*li(1), dir*ui(1))
        ie = dir*min(dir*li(1), dir*ui(1))

        do off=dir, dir*self%ng(1), dir
          mem(is+off,li(2):ui(2),:,:,:) = mem(ie-dir+off,li(2):ui(2),:,:,:)
        end do
      end do

      do dir=-1,1,2
        is = dir*max(dir*li(2), dir*ui(2))
        ie = dir*min(dir*li(2), dir*ui(2))

        do off=dir, dir*self%ng(2), dir
          mem(:,is+off,:,:,:) = mem(:,ie-dir+off,:,:,:)
        end do
      end do
    endif
  endif
END SUBROUTINE m3d_periodic_grid

!===============================================================================
character(len=17) FUNCTION m3d_name (self)
  class(m3d_t):: self
  !----------------------------------------------------------------------------
  write(m3d_name, '(A12,I5)') trim(self%kind), self%id
END FUNCTION

!===============================================================================
SUBROUTINE extrapolate_guards (self, target_mem)
  class(m3d_t):: self
  real(4), pointer, optional:: target_mem(:,:,:,:)
  !.............................................................................
  ! do nothing
END SUBROUTINE

!===============================================================================
SUBROUTINE check_topbot (self)
  class(m3d_t):: self
  !.............................................................................
  if (self%patch%boundaries%topbot(1) .and. self%dir < 0) then
    self%boundaries%topbot(1) = .true.
  elseif (self%patch%boundaries%topbot(2) .and. self%dir > 0) then
    self%boundaries%topbot(2) = .true.
  endif
END SUBROUTINE

!===============================================================================
!> Check if self (which is an nbor task) is ahead of target (which is the one to
!> possibly move to the ready queue), using self%dtime*target%grace as the grace
!> period, since we want to limit the extrapolation in the nbor task to at most
!> target%grace*self%dtime.
!===============================================================================
LOGICAL FUNCTION is_ahead_of (self, target)
  class (m3d_t):: self
  class (task_t):: target
  !.............................................................................
  real:: grace
  real(8):: nbtime, nbdtime, tgtime
  !-----------------------------------------------------------------------------
  grace = self%grace + shared%extra_grace
  if (omp_lock%tasks) then
                                            call self%lock%set ('is_ahead_of')
    !---------------------------------------------------------------------------
    ! dt(self%dt) is always consistent with t(:), while dtime changes in
    ! courant_condition(), before other attributes change in rotate()
    !---------------------------------------------------------------------------
    nbtime = self%time
    nbdtime = self%dt(self%it)
    tgtime = target%time
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    else if (self%istep < 3) then
      is_ahead_of = nbtime > tgtime
    else
      is_ahead_of = nbtime + nbdtime*grace > tgtime
    end if
                                            call self%lock%unset ('is_ahead_of')
  else
    !$omp atomic read
    nbtime = self%time
    !$omp atomic read
    nbdtime = self%dtime
    !$omp atomic read
    tgtime = target%time
    ! -- if a task is frozen (in time) assume it is forever ahead of other tasks
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    ! -- a negative nbor time signals an invalid (e.g. initial RT) state
  !  else if (nbtime < 0.0) then
  !    is_ahead_of = .false.
    ! -- do not use a grace interval for the first few steps
    else if (self%istep < 3) then
      is_ahead_of = nbtime > tgtime
    ! -- use a grace interval that is fraction of the nbor time step
    else
      is_ahead_of = nbtime + nbdtime*grace > tgtime
    end if
  end if
END FUNCTION is_ahead_of
    
end module m3d_patch_mod
