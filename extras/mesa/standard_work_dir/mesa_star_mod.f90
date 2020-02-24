! ---------------------------------------------------------------------
!! MODULE FOR MESA WORKFLOW CONTROL
! ---------------------------------------------------------------------
!
! x_logical_ctrl(1) : logical for terminating (if .false.) evolution of star in stellar_evolution.f90
! x_integer_ctrl(30), x_integer_ctrl(31) = star_id, mpi_rank
!
! Plenty of MESA subroutines are included with limited changes to the work flow
! -> these are labelled the original MESA name and '_multi' to distingish from the original MESA routines.
!
! Other rutines such as 'init_mesa_star' and 'continue_star' are build around larger mesa subroutines
!  but cut short to disrupt the classic MESA work flow as needed.

module mesa_star_mod
  use run_star_extras
  use star_lib
  use run_star_support
  use star_def
  use init, only: alloc_star_data
  implicit none

  type (star_info), pointer :: s, s1

  integer :: n_stars = 0
  integer :: n_stars_total = 0


  contains

        character(len=20) function str(k)
            !   "Convert an integer to string. From Stackoverflow"
            integer, intent(in) :: k
            write (str, *) k
            str = adjustl(str)
        end function str

        ! ####################################################################
        ! ####################################################################

        subroutine create_MESA_star(myid, birth_time, fname_in)
          integer :: id, ierr
      	  integer, intent(in) :: myid
      	  real(kind=8), intent(in) :: birth_time
          character(len=255), intent(in), optional :: fname_in
          character(len=255) :: fname
          logical, parameter :: &
             do_alloc_star = .true., &
             do_free_star = .true., &
             okay_to_restart = .false.
          logical :: restart

          if (present(fname_in)) then
            fname = fname_in
          else
      	    fname = 'mesa_inlist'
          endif

          !! KEEPING TRACK OF STARS:
          n_stars = n_stars + 1
          n_stars_total = n_stars_total + 1

          if (n_stars > 10) then
            write(*,*)
            write(*,*) "ERROR: Cannot initiate more than 10 stars per process (MPI rank)"
            write(*,*)
            n_stars = n_stars - 1
            return
          else
            !! READ INLIST:
            ierr = 0
            call do_read_star_job_multi(fname, ierr, id)
            if (ierr /= 0) stop 1
            !! INIT / allocate star in MESA:
            call init_mesa_star(do_alloc_star, do_free_star, okay_to_restart, restart, id, ierr, fname, myid)
            if (ierr == 0) then
              write(*,*) "started star, ID: ", id
      	      write(*,*) "started star by MPI id", myid
            endif
      	    !! Set start time in star pointer
      	    call star_ptr(id, s, ierr)
      	    s% x_ctrl(26) = birth_time
      	    s% x_logical_ctrl(1) = .true.
      	    s% x_integer_ctrl(30) = id ! star id
      	    s% x_integer_ctrl(31) = myid ! rank owner of sink/star
          endif

        end subroutine create_MESA_star


        ! ####################################################################
        ! ####################################################################

	      ! Almost identical to create_MESA_star - but with restart logicals to true.
        subroutine restart_MESA_star(myid, birth_time, photo_path, fname_in)
          integer :: id, ierr
	        integer, intent(in) :: myid
	        real(kind=8), intent(in) :: birth_time
          character(len=255), intent(inout) :: photo_path
          character(len=255), intent(in), optional :: fname_in
          character(len=255) :: fname
          logical, parameter :: &
             do_alloc_star = .true., &
             do_free_star = .true., &
             okay_to_restart = .true.
          logical :: restart = .true.

          if (present(fname_in)) then
            fname = fname_in
          else
	          fname = 'mesa_inlist'
          endif

          !! KEEPING TRACK OF STARS:
          n_stars = n_stars + 1
          n_stars_total = n_stars_total + 1

          if (n_stars > 10) then
            write(*,*)
            write(*,*) "ERROR: Cannot initiate more than 10 stars per process (MPI rank)"
            write(*,*)
            n_stars = n_stars - 1
            return
          else
            !! READ INLIST:
            ierr = 0
            call do_read_star_job_multi(fname, ierr, id)
            if (ierr /= 0) stop 1

            !! INIT / allocate star in MESA:
            call init_mesa_star(do_alloc_star, do_free_star, okay_to_restart, restart, id, ierr, fname, myid, photo_path)
            if (ierr == 0) then
              write(*,*) "started star, ID: ", id
	            write(*,*) "started star by MPI id", myid
            endif
      	    !! Set start time in star pointer
      	    call star_ptr(id, s, ierr)
      	    s% x_ctrl(26) = birth_time
      	    s% x_logical_ctrl(1) = .true.
      	    s% x_integer_ctrl(30) = id
      	    s% x_integer_ctrl(31) = myid
          endif

        end subroutine restart_MESA_star

        ! ####################################################################
        ! ####################################################################

        subroutine init_mesa_star(do_alloc_star, do_free_star, okay_to_restart, restart, id, &
          ierr, inlist_fname_arg, mpi_id, photo_path)
      	  ! init_mesa_star contains the first part of the run1_star routine in run_star_support.f90
      	  ! It cuts of before the evolve-loop commands as we break the MESA work flow here.
          integer :: ierr, mpi_id
          logical, intent(in) :: do_alloc_star, do_free_star, okay_to_restart
          integer, intent(inout) :: id ! input if not do_alloc_star
          logical, intent(inout) :: restart ! input if not do_alloc_star
          character (len=*) :: inlist_fname_arg
          optional inlist_fname_arg
          character(len=255), intent(inout), optional :: photo_path
          character(len=255) :: log_dir
          character(len=255) :: photo_dir

          integer :: id_extra, result, model_number, source, nsteps, j, ci, nz
          logical :: continue_evolve_loop, first_try, just_did_backup
          character (len=strlen) :: inlist_fname

          logical, parameter :: dbg = .false.

      	  inlist_fname = trim(inlist_fname_arg)
          ! dir name for star:
          write(log_dir,*) './MESA/LOGS_mpi_'//trim(str(mpi_id))//'_star_'//trim(str(n_stars_total))
          write(photo_dir,*) './MESA/photos_mpi_'//trim(str(mpi_id))//'_star_'//trim(str(n_stars_total))

          log_dir = log_dir(2:) !! REMOVING WHISPACE IN FRONT OF STRING.
          photo_dir = photo_dir(2:) !! REMOVING WHISPACE IN FRONT OF STRING.

          ! star is initiated here
          if (present(photo_path) .and. okay_to_restart) then                     ! Restart
            call before_evolve_loop_multi(do_alloc_star, okay_to_restart, restart, &
                 null_binary_controls, extras_controls, &
                 id_from_read_star_job, inlist_fname, photo_path, &
                 dbg, 0, id, id_extra, ierr, log_dir, photo_dir)
          else                                                                    ! Create star
            call before_evolve_loop_multi(do_alloc_star, okay_to_restart, restart, &
                 null_binary_controls, extras_controls, &
                 id_from_read_star_job, inlist_fname, "restart_photo", &
                 dbg, 0, id, id_extra, ierr, log_dir, photo_dir)
          endif
          if (failed('before_evolve_loop',ierr)) return

        end subroutine init_mesa_star


        ! ####################################################################
        ! ####################################################################

        subroutine continue_star(id, ierr)
          integer :: id, ierr
          integer :: id_extra, result, model_number, source, nsteps, j, ci, nz
          logical :: continue_evolve_loop, first_try, just_did_backup
	  logical :: custom_stopping_condition
          logical :: dbg
          character (len=255) :: inlist_fname

          inlist_fname = 'mesa_inlist'

      	  print *, 'id: ', id
          call star_ptr(id, s, ierr)
          if (failed('star_ptr',ierr)) return

          continue_evolve_loop = .true.
          s% doing_timing = .false.
          s% job% check_before_step_timing = 0
          s% job% check_step_loop_timing = 0
          s% job% check_after_step_timing = 0
          s% job% time0_initial = 0

          1 format(a35, 99(1pe26.16))
          2 format(a35, i7, 1pe26.16)
          3 format(a15, 2x, f15.6)
          4 format(a15, 2x, e15.6)

          11 format(a35, f20.10)

          dbg = .false.
          if (dbg) write(*,*) 'start evolve_loop'
          evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop

             call before_step_loop(s, ierr)
             if (failed('before_step_loop',ierr)) return

             first_try = .true.
             just_did_backup = .false.

             step_loop: do ! may need to repeat this loop

                if (stop_is_requested(s)) then
                   continue_evolve_loop = .false.
                   result = terminate
                   exit
                end if

                result = star_evolve_step(id, first_try, just_did_backup)
                if (result == keep_going) result = star_check_model(id)
                if (result == keep_going) result = s% extras_check_model(id, id_extra)
                if (result == keep_going) result = star_pick_next_timestep(id)
                if (result == keep_going) exit step_loop

                model_number = get_model_number(id, ierr)
                if (failed('get_model_number',ierr)) return

                if (result == retry .and. s% job% report_retries) then
                   write(*,'(i6,3x,a,/)') model_number, &
                      'retry reason ' // trim(result_reason_str(s% result_reason))
                else if (result == backup .and. s% job% report_backups) then
                   write(*,'(i6,3x,a,/)') model_number, &
                      'backup reason ' // trim(result_reason_str(s% result_reason))
                end if

                if (result == redo) then
                   result = star_prepare_to_redo(id)
                end if
                if (result == retry) then
                   result = star_prepare_to_retry(id)
                end if
                if (result == backup) then
                   result = star_do1_backup(id)
                   just_did_backup = .true.
                else
                   just_did_backup = .false.
                end if
                if (result == terminate) then
                   continue_evolve_loop = .false.
                   exit step_loop
                end if
                first_try = .false.

             end do step_loop

             ! once we get here, the only options are keep_going or terminate.
             ! redo, retry, or backup must be done inside the step_loop

             call after_step_loop(s, inlist_fname, &
                 id_extra, dbg, result, ierr)
             if (failed('after_step_loop',ierr)) return


      	     ! Custom stopping conditions go here:
	     custom_stopping_condition = .false.
	     if (custom_stopping_condition) then
		if (s% star_mass > 0.5) then ! This is for custom stopping conditions for the stars - may need extension!
        		    result = terminate
        		    write(*,*)
        		    write(*,*) 'Terminating star due to custom condition in continue_star'
        		    write(*,*)
        		    call stop_MESA_star(id, ierr)   ! Free MESA handles.
        		    s% x_logical_ctrl(1) = .false.  ! Ignore this star in DISPATCH stellar_evolution module.
        		    continue_evolve_loop = .false.
		endif
	     endif

             if (result /= keep_going) then
                if (result /= terminate) then
                   write(*,2) 'ERROR in result value in run_star_extras: model', &
                      s% model_number
                   write(*,2) 'result', result
                   exit evolve_loop
                end if
                if (s% result_reason == result_reason_normal) then
                   call terminate_normal_evolve_loop(s, &
                      id_extra, dbg, result, ierr)
                   if (failed('terminate_normal_evolve_loop',ierr)) return
                else
                  write(*,*)
                  write(*,*) 'calling stop_MESA_star due to non-normal result:'
                  write(*,*) '----------------------------------'
                  write(*,*) result_reason_str(s% result_reason)
                  write(*,*) '----------------------------------'
                  write(*,*)
                  call stop_MESA_star(id, ierr)
            		  s% x_logical_ctrl(1) = .false.  !! Maybe this should not be here?
                end if
                exit evolve_loop
             end if


             call do_saves( &
                id, id_extra, s)

             if (s% doing_timing) then
                call system_clock(s% job% time1_extra,s% job% clock_rate)
                s% job% after_step_timing = s% job% after_step_timing + &
                   dble(s% job% time1_extra - s% job% time0_extra) / s% job% clock_rate
                s% job% check_time_end = eval_total_times(s% id, ierr)
                s% job% check_after_step_timing = s% job% check_after_step_timing + &
                   (s% job% check_time_end - s% job% check_time_start)
             end if

          end do evolve_loop
        end subroutine continue_star


        ! ####################################################################
        ! ####################################################################

	subroutine write_photos(dir, id)
	  use evolve_support, only: output_to_file
	  integer, intent(inout) :: id
	  integer :: ierr
	  character(len=*), intent(in) :: dir
	  character(len=255) :: cp_command, photo_str, destination_str

	  call star_ptr(id, s, ierr)
	  if (failed('star_ptr',ierr)) return

	  photo_str = 'photo_mpi_'// trim(str(s% x_integer_ctrl(31))) // '_id_' // trim(str(s% x_integer_ctrl(30)))
	  print *,photo_str
	  destination_str = './data/'// trim(dir) // 'MESA/'
	  write(cp_command,*) 'mv '// trim(photo_str) //' '// trim(destination_str) ! command string to move current photo to RAMSES output directory

	  print *, 'copy photo: '
	  print *, cp_command
	  call output_to_file(photo_str, id, ierr) ! Write photo to be moved.
	  call execute_command_line(cp_command) ! Move to photo...


	end subroutine write_photos

    ! ####################################################################
    ! ####################################################################

    subroutine stop_MESA_star(id, ierr)
      ! Stop MESA star model - SN explosion?
      integer :: id, ierr, id_extra
      logical, parameter :: do_free_star = .true.
      call star_ptr(id, s, ierr)
      if (failed('star_ptr',ierr)) return

      call free_star_handle(s)

      call after_evolve_loop(s, id_extra, do_free_star, ierr)
      if (failed('after_evolve_loop',ierr)) return

      n_stars = n_stars - 1

    end subroutine stop_MESA_star

    ! ####################################################################
    ! ####################################################################

    subroutine free_star_handle(s)
  	   ! Free handles if a star is no longer evolved - SN explosion?
       use eos_def
       use kap_def
       type (star_info), pointer :: s
       star_handles(s% id)% in_use = .false.
       eos_handles(s% id)% in_use = .false.
       kap_handles(s% id)% in_use = .false.
    end subroutine free_star_handle

    ! ####################################################################
    ! ####################################################################

    subroutine MESA_shutdown_starlib
      ! Shutdown MESA libraries if ever needed
      call starlib_shutdown()
    end subroutine MESA_shutdown_starlib

    ! ####################################################################
    ! ####################################################################

    subroutine do_read_star_job_multi(filename, ierr, id)
       character(*), intent(in) :: filename
       integer, intent(out) :: ierr
       character(len=strlen) :: inlist_fname
       integer, intent(inout) :: id
       type (star_info), pointer :: s

       include 'formats'
       ierr = 0

      !  if (id_from_read_star_job /= 0) then
      !     write(*,2) 'id_from_read_star_job', id_from_read_star_job
      !     ierr = -1
      !     return
      !  end if

       id = alloc_star(ierr)
       if (ierr /= 0) then
          write(*,*) 'do_read_star_job failed in alloc_star'
          return
       end if

       call star_ptr(id, s, ierr)
       if (ierr /= 0) then
          write(*,*) 'do_read_star_job failed in star_ptr'
          return
       end if

  	  inlist_fname = trim(filename)
       call read_star_job(s, inlist_fname, ierr)
       if (ierr /= 0) then
          write(*,*) 'ierr from read_star_job ' // trim(inlist_fname)
          return
       end if

       id_from_read_star_job = id

       if (s% job% save_star_job_namelist) then
          call write_star_job(s, s% job% star_job_namelist_name, ierr)
          if (ierr /= 0) then
             write(*,*) 'ierr from write_star_job ' // &
                trim(s% job% star_job_namelist_name)
             return
          end if
       end if

    end subroutine do_read_star_job_multi

    ! ####################################################################
    ! ####################################################################

    subroutine step_MESA(mdots, timesteps, istar)
      !real(kind=8), dimension(:), allocatable, intent(in) :: mdots, timesteps
      real(kind=8), parameter :: dt_start = 3600.0*24
  	  real(kind=8), intent(inout) :: mdots, timesteps
      integer :: ierr, i, istar
      type (star_info), pointer :: s1

  	  call star_ptr(istar, s1, ierr)
  	  s1% dt_next = max(dt_start, s1% dt_next_unclipped)
  	  s1% max_age_in_seconds = timesteps + (s1% star_age*secyer)
  	  s1% x_ctrl(50) = mdots
  	  call continue_star(istar, ierr)

    end subroutine step_MESA

    ! ####################################################################
    ! ####################################################################

    subroutine before_evolve_loop_multi(do_alloc_star, okay_to_restart, restart, &
            binary_controls, extras_controls, &
            id_from_read_star_job, inlist_fname, restart_filename, &
            dbg, binary_id, id, id_extra, ierr, log_dir, photo_dir)
       !use se_support, only: se_startup
       use utils_lib
       use ctrls_io
       logical, intent(in) :: do_alloc_star, okay_to_restart
       logical :: restart
       interface
          subroutine binary_controls(id, binary_id, ierr)
             integer, intent(in) :: id, binary_id
             integer, intent(out) :: ierr
          end subroutine binary_controls
          subroutine extras_controls(id, ierr)
             integer, intent(in) :: id
             integer, intent(out) :: ierr
          end subroutine extras_controls
       end interface
       integer :: id_from_read_star_job
       character (len=*) :: inlist_fname, restart_filename
       logical, intent(in) :: dbg
       integer, intent(in) :: binary_id
       integer, intent(out) :: id, id_extra, ierr
       character(len=255) :: log_dir
       character(len=255) :: photo_dir

       type (star_info), pointer :: s

       if (do_alloc_star) then
          if (id_from_read_star_job /= 0) then
             ! already allocated by read_star_job
             id = id_from_read_star_job
             id_from_read_star_job = 0
          else
             id = alloc_star(ierr)
             if (failed('alloc_star',ierr)) return
          end if
          call star_ptr(id, s, ierr)
          if (failed('star_ptr',ierr)) return
       else
          call star_ptr(id, s, ierr)
          if (failed('star_ptr',ierr)) return
          call init_starting_star_data(s, ierr)
          if (failed('init_starting_star_data',ierr)) return
       end if

       if (dbg) write(*,*) 'call starlib_init'
       call starlib_init(s, ierr) ! okay to do extra calls on this
       if (failed('star_init',ierr)) return

       if (dbg) write(*,*) 'call star_set_kap_and_eos_handles'
       call star_set_kap_and_eos_handles(id, ierr)
       if (failed('set_star_kap_and_eos_handles',ierr)) return

       if (dbg) write(*,*) 'call star_setup'
       call star_setup(id, inlist_fname, ierr)
       if (failed('star_setup',ierr)) return

       !! -------------------------------------------
       !! SPECIFY DIRECTORIES:
       s% log_directory = log_dir
       call mkdir(s% log_directory)
       s% photo_directory = photo_dir
       call mkdir(s% photo_directory)
       !! -------------------------------------------

       if (len_trim(s% op_mono_data_path) == 0) &
          call get_environment_variable( &
             "MESA_OP_MONO_DATA_PATH", s% op_mono_data_path)

       if (len_trim(s% op_mono_data_cache_filename) == 0) &
          call get_environment_variable( &
             "MESA_OP_MONO_DATA_CACHE_FILENAME", s% op_mono_data_cache_filename)
      !  if (restart_filename /= "restart_photo") &
      !     restart_filename  = trim(s% photo_directory) // '/' // trim(restart_filename)

       if (okay_to_restart) then
          restart = doing_a_restart(restart_filename)
          print *,'Doing a restart'
          print *, restart_filename
       else
          restart = .false.
       end if

       if (s% job% show_log_description_at_start .and. .not. restart) then
          write(*,*)
          call show_log_description(id, ierr)
          if (failed('show_log_description',ierr)) return
       end if

       s% extras_startup => null_extras_startup
       s% extras_check_model => null_extras_check_model
       s% extras_finish_step => null_extras_finish_step
       s% extras_after_evolve => null_extras_after_evolve
       s% how_many_extra_history_columns => null_how_many_extra_history_columns
       s% data_for_extra_history_columns => null_data_for_extra_history_columns
       s% how_many_extra_profile_columns => null_how_many_extra_profile_columns
       s% data_for_extra_profile_columns => null_data_for_extra_profile_columns

       call extras_controls(id, ierr)
       if (ierr /= 0) return

       call binary_controls(id, binary_id, ierr)
       if (ierr /= 0) return

       if (dbg) write(*,*) 'call do_star_job_controls_before'
       call do_star_job_controls_before(id, s, restart, ierr)
       if (ierr /= 0) return

       if (dbg) write(*,*) 'call do_load1_star'
       call do_load1_star(id, s, restart, restart_filename, ierr)
       if (failed('do_load1_star',ierr)) return

       if (dbg) write(*,*) 'call do_star_job_controls_after'
       call do_star_job_controls_after(id, s, restart, ierr)
       if (failed('do_star_job_controls_after',ierr)) return

       write(*,*)
       write(*,*)

       if (.not. restart) then
          call before_evolve(id, ierr)
          if (failed('before_evolve',ierr)) return
          call start_new_run_for_pgstar(s, ierr)
          if (failed('start_new_run_for_pgstar',ierr)) return
       else
          call show_terminal_header(id, ierr)
          if (failed('show_terminal_header',ierr)) return
          call restart_run_for_pgstar(s, ierr)
          if (failed('restart_run_for_pgstar',ierr)) return
       end if

       id_extra = s% extras_startup(id, restart, ierr)
       if (failed('extras_startup',ierr)) return

       !call se_startup(s, id, restart, s% job% use_se_output, ierr) ! Legacy stuff from r8845
       !if (failed('se_startup',ierr)) return

       if (s% job% profile_starting_model) then
          write(*, '(a, i12)') 'save profile for model number ', s% model_number
          call star_set_vars(id, id_extra, 0d0, ierr)
          if (failed('star_set_vars',ierr)) return
          call save_profile(id, id_extra,3, ierr)
          if (failed('save_profile',ierr)) return
       end if

    end subroutine before_evolve_loop_multi

    ! ####################################################################
    ! ####################################################################

end module mesa_star_mod
