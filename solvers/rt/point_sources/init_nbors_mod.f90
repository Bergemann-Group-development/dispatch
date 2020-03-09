MODULE init_nbor_mod
  USE rt_mod
  USE solver_mod
  USE link_mod
  USE task_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
CONTAINS

!===============================================================================
!> In the new, synced RT/MHD arrangement the steps and relations needed are:
!>
!> 1) The omega0 task needs the MHD nbors to dnload and call EOS for RT
!> 2) Each omega task needs the upstream omega task
!> 3) Each omega task neds the omega0 task
!> 4) The MHD task needs all omega task, while not vice versa.
!>
!> The upstream nbors should have the settings needed=T, needs_me=F, download=T,
!> while the downstream  tasks should have needed=F, needs_me=T, download=F.
!===============================================================================
SUBROUTINE init_rt_nbors (self)
  class (solver_t):: self
  !.............................................................................
  integer:: i_omega
  class(rt_t), pointer:: omega0, omega, upstream
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: nb_task
  logical:: needs, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup_t%init_nbors')
  omega0 => self%rt
  call omega0%clear(bits%init_nbors)
  !-----------------------------------------------------------------------------
  ! 1) The omega0 needs the MHD task for the MHD variables, but no download
  !-----------------------------------------------------------------------------
  needs=.true.; needs_me=.false.; download=.false.
  call init_nbor_pair (omega0, needs, self, needs_me, download)
  !-----------------------------------------------------------------------------
  ! Loop over all ray directions, and add nbors and nbor relations. The only
  ! downloads that are needed are from upstream omega tasks
  !-----------------------------------------------------------------------------
  do i_omega=1,omega0%n_omega
    omega => omega0%omega(i_omega)
    !---------------------------------------------------------------------------
    ! Loop over all MHD nbors
    !---------------------------------------------------------------------------
    nbor => self%link%nbor
    do while (associated(nbor))
      !-------------------------------------------------------------------------
      ! Cast the nbors to solver_t type
      !-------------------------------------------------------------------------
      nb_task => nbor%task
      select type (nb_task)
      class is (solver_t)
      !-------------------------------------------------------------------------
      ! Select only omega tasks
      !-------------------------------------------------------------------------
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
        !-----------------------------------------------------------------------
        ! 4) The omega task needs the RT task, which does not need the omega tasks
        !-----------------------------------------------------------------------
        needs=.true.; needs_me=.false.; download=.false.
        call init_nbor_pair (omega, needs, omega0, needs_me, download)
      end if
      end select
      nbor => nbor%next
    end do
    !---------------------------------------------------------------------------
    ! 5) The MHD task needs the omega tasks, which do not need the MHD task
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (self, needs, omega, needs_me, download)
  end do
  call trace%end
END SUBROUTINE init_rt_nbors_sc

END MODULE init_nbor_mod
