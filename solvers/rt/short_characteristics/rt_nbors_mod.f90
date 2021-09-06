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
  USE io_mod
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
  class(nbor_t), pointer:: nbor
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
  call init_nbor_pair (omega0%link, needs, mhd%link, needs_me, download)
  omega => omega0%omega
  do while (associated(omega))
    !---------------------------------------------------------------------------
    ! 2) The omega task needs the RT task, which does not need the omega tasks
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (omega%link, needs, omega0%link, needs_me, download)
    !---------------------------------------------------------------------------
    ! 3) The MHD task needs the omega tasks, which do not need the MHD task
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (mhd%link, needs, omega%link, needs_me, download)
    omega => omega%next
  end do    
  !-----------------------------------------------------------------------------
  ! Now loop over all MHD nbors
  !-----------------------------------------------------------------------------
  nbor => mhd%link%nbors_by_level
  do while (associated(nbor))
    !---------------------------------------------------------------------------
    ! Cast the nbors to solver_t type
    !---------------------------------------------------------------------------
    nb_task => nbor%task
    select type (nb_task)
    class is (solver_t)
    !---------------------------------------------------------------------------
    ! 4) The omega0 task needs the MHD nbors, and needs download
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.true.
    call init_nbor_pair (omega0%link, needs, nb_task%link, needs_me, download)
    !---------------------------------------------------------------------------
    ! Loop over all ray directions, and add nbors and nbor relations. The only
    ! downloads that are needed are from upstream omega tasks
    !---------------------------------------------------------------------------
    omega => omega0%omega
    do while (associated(omega))
      !-------------------------------------------------------------------------
      ! Select omega sub-tasks
      !-------------------------------------------------------------------------
      if (associated(nb_task%rt)) then
        !-----------------------------------------------------------------------
        ! The omega0%needs() function tests for upstream nbor omega
        !-----------------------------------------------------------------------
        upstream => nb_task%rt%omega
        do while (associated(upstream) .and. (upstream%i_omega/=omega%i_omega))
          upstream => upstream%next
        end do
        if (.not.associated(upstream)) call io%abort('rt_nbors_mod: upstream err')
        if (omega%needs (omega, upstream) .and. upstream%on) then
          !---------------------------------------------------------------------
          ! 5) The omega tasks need the upstream tasks, and should download from
          ! them, but the upstream tasks do not need the downstream omega tasks
          !---------------------------------------------------------------------
          needs=.true.; needs_me=.false.; download=.true.
          call init_nbor_pair (omega%link, needs, upstream%link, needs_me, download)
        end if
      end if
      omega => omega%next
    end do
    end select
    nbor => nbor%next
  end do
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Add mutual nbor relations to two nbor lists, making sure that the needed and
!> needs_me flags are mutually consistent.
!===============================================================================
SUBROUTINE init_nbor_pair (link1, needed, link2, needs_me, download)
  class(link_t), pointer:: link1, link2
  logical:: needed, needs_me, download
  !.............................................................................
  class(nbor_t), pointer:: nbor1, nbor2
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbor_pairs')
  call link1%add_nbor_by_rank (link2, nbor1)
  call link2%add_nbor_by_rank (link1, nbor2)
  !-----------------------------------------------------------------------------
  ! Add the 2nd linked task as an nbor to the 1st linked task
  !-----------------------------------------------------------------------------
  nbor1%needed = needed
  nbor1%needs_me = needs_me
  nbor1%download = download
  !-----------------------------------------------------------------------------
  ! Add the 1st linked task as an nbor to the 2nd linked task
  !-----------------------------------------------------------------------------
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
  call trace%end()
END SUBROUTINE init_nbor_pair

END MODULE rt_nbors_mod
