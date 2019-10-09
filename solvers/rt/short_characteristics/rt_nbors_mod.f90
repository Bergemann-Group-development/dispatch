!===============================================================================
!> Module to add nbor lists for RT sub-tasks, and modify the parent MHD tasks
!===============================================================================
MODULE rt_nbors_mod
  USE solver_mod
  USE rt_mod
  USE patch_mod
  USE link_mod
  USE task_mod
  USE trace_mod
  USE io_unit_mod
  USE bits_mod
  implicit none
  private
  !.............................................................................
  type, public:: rt_nbors_t
  contains
    procedure, nopass:: init
  end type
  type(rt_nbors_t), public:: rt_nbors
CONTAINS

!===============================================================================
!> In the syncronized RT/MHD arrangement the steps and relations needed are:
!>
!> 1) The omega0 task needs the MHD nbors to dnload and call EOS for RT
!> 2) Each omega task needs the upstream omega task
!> 3) Each omega task neds the omega0 task
!> 4) The MHD task needs all omega task, while not vice versa.
!>
!> The upstream nbors should have the settings needed=T, needs_me=F, download=T,
!> while the downstream  tasks should have needed=F, needs_me=T, download=F.
!===============================================================================
SUBROUTINE init (mhd)
  class (solver_t):: mhd
  !.............................................................................
  integer:: i_omega
  class(rt_t), pointer:: omega0, omega, upstream
  class(link_t), pointer:: nbor
  class(task_t), pointer:: nb_task
  logical:: needs, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_nbors_t%init')
  omega0 => mhd%rt
  call omega0%clear(bits%init_nbors)
  !-----------------------------------------------------------------------------
  ! 1) The omega0 needs the MHD task for the MHD variables, but no download
  !-----------------------------------------------------------------------------
  needs=.true.; needs_me=.false.; download=.false.
  call init_nbor_pair (omega0, needs, mhd, needs_me, download)
  !-----------------------------------------------------------------------------
  ! Loop over all ray directions, and add nbors and nbor relations. The only
  ! downloads that are needed are from upstream omega tasks
  !-----------------------------------------------------------------------------
  do i_omega=1,omega0%n_omega
    omega => omega0%omega(i_omega)
    !---------------------------------------------------------------------------
    ! Loop over all MHD nbors
    !---------------------------------------------------------------------------
    nbor => mhd%link%nbor
    do while (associated(nbor))
      !-------------------------------------------------------------------------
      ! Cast the nbors to solver_t type
      !-------------------------------------------------------------------------
      nb_task => nbor%task
      select type (nb_task)
      class is (solver_t)
      !-------------------------------------------------------------------------
      ! Select omega sub-tasks
      !-------------------------------------------------------------------------
      if (associated(nb_task%rt)) then
      if (nb_task%rt%n_omega == 0) then
        !-----------------------------------------------------------------------
        ! 2) The omega0 task needs the MHD nbors, and needs download
        !-----------------------------------------------------------------------
        needs=.true.; needs_me=.false.; download=.true.
        call init_nbor_pair (omega0, needs, nb_task, needs_me, download)
      else
        !-----------------------------------------------------------------------
        ! The omega0%needs() function tests for upstream nbor omega
        !-----------------------------------------------------------------------
        upstream => nb_task%rt%omega(i_omega)
        if (omega%needs (omega, upstream) .and. upstream%on) then
          !---------------------------------------------------------------------
          ! 3) The omega tasks need the upstream tasks, and should download from
          ! them, but the upstream tasks do not need the downstream omega tasks
          !---------------------------------------------------------------------
          needs=.true.; needs_me=.false.; download=.true.
          call init_nbor_pair (omega, needs, upstream, needs_me, download)
        end if
      end if
      end if
      end select
      nbor => nbor%next
    end do
    !-----------------------------------------------------------------------
    ! 4) The omega task needs the RT task, which does not need the omega tasks
    !-----------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (omega, needs, omega0, needs_me, download)
    !---------------------------------------------------------------------------
    ! 5) The MHD task needs the omega tasks, which do not need the MHD task
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (mhd, needs, omega, needs_me, download)
  end do
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Add mutual nbor relations to two nbor lists, making sure that the needed and
!> needs_me flags are mutually consistent.
!===============================================================================
SUBROUTINE init_nbor_pair (task1, needed, task2, needs_me, download)
  class(task_t):: task1, task2
  logical:: needed, needs_me, download
  class(link_t), pointer:: link1, link2, nbor1, nbor2
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbor_pairs')
  select type (task1)
  class is (patch_t)
  link1 => task1%link
  end select
  select type (task2)
  class is (patch_t)
  link2 => task2%link
  end select
  !-----------------------------------------------------------------------------
  ! Add the 2nd linked task as an nbor to the 1st linked task
  !-----------------------------------------------------------------------------
  allocate (nbor1)
  nbor1%link => link2
  nbor1%task => nbor1%link%task
  nbor1%needed = needed
  nbor1%needs_me = needs_me
  nbor1%download = download
  !-----------------------------------------------------------------------------
  ! Add the 1st linked task as an nbor to the 2nd linked task
  !-----------------------------------------------------------------------------
  allocate (nbor2)
  nbor2%link => link1
  nbor2%task => nbor2%link%task
  nbor2%needed = needs_me
  nbor2%needs_me = needed
  nbor2%download = .false.
  !-----------------------------------------------------------------------------
  ! If the 1st task is a virtual task, it should not download from the 2nd task
  ! and also does not need it (in the sense that it gets updated by MPI). If
  ! the needs me flag is also false, then none of two tasks should be in the
  ! nbor list of the other (this applies if the 1st task is a virtual downstream
  ! RT task)
  !-----------------------------------------------------------------------------
  if (link1%task%is_set (bits%virtual)) then
    nbor1%download = .false.
    nbor1%needed = .false.
    nbor2%needs_me = .false.
    if (.not.nbor1%needs_me) then
      call trace%end()
      return
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If the 2nd task is a virtual task the same applies, except the meaning of
  ! the flags is reversed (this applies if the 2nd task is a virtual downstream
  ! RT task)
  !-----------------------------------------------------------------------------
  if (link2%task%is_set (bits%virtual)) then
    nbor2%download = .false.
    nbor2%needs_me = .false.
    nbor1%needed = .false.
    if (.not.nbor2%needed) then
      call trace%end()
      return
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Add the two tasks, which now are guaranteed to have matching flags
  !-----------------------------------------------------------------------------
  call link1%add_nbor_by_rank (link1%nbor, nbor1)
  call link2%add_nbor_by_rank (link2%nbor, nbor2)
  call trace%end()
END SUBROUTINE init_nbor_pair

END MODULE rt_nbors_mod
