!===============================================================================
!> OpenMULTI3D main program.
!>
!> The call hierarchy is
!>
!> dispatch
!>   setup%init                 ! MPI, I/O and log files
!>   dispatcher%init            ! select and initialize the dispatcher
!>   read_external%init         ! read the atmosphere and build the task list
!>   dispatcher%execute         ! run the task list
!>   write_linelist_out         ! line-list diagnostics, if collected
!>   setup%end                  ! shut down
!>
!===============================================================================
PROGRAM dispatch
  USE mpi_mod
  USE os_mod
  USE setup_mod
  USE console_mod
  USE dispatcher_mod
  USE task_list_mod
  USE experiment_mod
  USE read_external_mod
  USE timer_mod
  USE mpi_debug_mod
  USE m3d_patch_mod, only: write_linelist_out, n_nlte_patches, converged, convlim, convmax, &
    nlte_last_change, nlte_last_iter, n_convmax_finished, stop_requested, stop_exit
  USE rt_params_mod, only: rt_iter
  USE rt_column_mod, only: report_memory_estimate, check_memory_limit
  class(task_list_t), pointer:: task_list => null()
  real(8):: memory_estimate
  character(len=160):: text
  character(len=16):: peak, estimate, nit, change, limit
  !.............................................................................
  call setup%init                               ! Generic setup (MPI, I/O, ...)
  call dispatcher%init                          ! Initialize the dispatcher
  call check_mpi_debug()
  call read_external%init (task_list)           ! Initialize task_list
  call report_memory_estimate (task_list, memory_estimate)
  call check_memory_limit (memory_estimate)
  call dispatcher%execute (task_list)           ! Run dispatcher on the task_list
  call write_linelist_out()                     ! Line-list diagnostics, if collected
  if (timer%verbose > 0) call timer%print() 
  call mpi%barrier()
  ! The registered peak at level 0; from level 1 on it is set against the estimate printed
  ! before the run.
  if (console_on(ch_m3d, 0, all_ranks=.true.)) then
    write(peak,'(f10.3)') os%max_mem
    write(text,'(a,a,a,i0,a)') 'Peak memory: ', trim(adjustl(peak)), ' GiB (rank ', mpi%rank, ')'
    if (console_on(ch_m3d, 1, all_ranks=.true.)) then
      write(estimate,'(f10.3)') memory_estimate
      text = trim(text) // ', estimate ' // trim(adjustl(estimate)) // ' GiB'
    endif
    write(console_unit,'(a)') trim(text)
  endif
  ! Success is an LTE run, NLTE iterations that met convlim, or NLTE iterations in
  ! which every depth patch met convmax and finished early (premature_finish); the
  ! last of these never completes the iteration that would measure a final global
  ! change.  A run that used up its iterations says so with the change it ended on.
  write(nit,'(i0)') nlte_last_iter
  if (stop_requested) then
    write(change,'(f10.2)') abs(nlte_last_change)*100
    text = 'OpenMULTI3D stopped on request after ' // trim(nit) // ' NLTE iterations (largest ' // &
      'relative change ' // trim(adjustl(change)) // ' %); departures written'
    if (stop_exit) text = trim(text) // ', no output spectrum'
    call console_warn (ch_m3d, trim(text))
  else if (n_nlte_patches == 0 .or. rt_iter%maxiter == 0 .or. converged .or. &
      (nlte_last_change /= 0d0 .and. abs(nlte_last_change) < convlim)) then
    if (console_on(ch_m3d, 0)) write(console_unit,'(a)') 'OpenMULTI3D terminated successfully.'
  else if (n_convmax_finished == n_nlte_patches) then
    write(limit,'(es10.1)') convmax*100
    if (console_on(ch_m3d, 0)) write(console_unit,'(a)') 'OpenMULTI3D terminated successfully: ' // &
      'every NLTE patch met convmax (' // trim(adjustl(limit)) // ' %) after ' // trim(nit) // ' iterations.'
  else
    write(change,'(f10.2)') abs(nlte_last_change)*100
    write(limit,'(f10.2)') convlim*100
    text = 'OpenMULTI3D terminated, but the NLTE iterations did not converge: after ' // trim(nit) // &
      ' iterations the largest relative change was ' // trim(adjustl(change)) // ' % (convlim ' // &
      trim(adjustl(limit)) // ' %)'
    call console_warn (ch_m3d, trim(text))
  endif
  call setup%end                                ! End setup
!===============================================================================
END PROGRAM dispatch

