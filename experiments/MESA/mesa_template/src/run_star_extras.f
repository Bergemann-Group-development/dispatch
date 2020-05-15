      module run_star_extras

         use star_lib
         use star_def
         use const_def
         
         implicit none       
         
         type acc_struct 
         sequence
            real(8) :: time, dm_acc
         end type acc_struct
         
         
         TYPE(acc_struct), DIMENSION(:), allocatable :: acc_hist
         TYPE(integer) :: acc_len = 0
         character(64) :: next_model_filename
         real(dp) :: next_model_save = huge(0d0)
         integer :: read_unit, write_unit

         contains        
         
         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s

            call star_ptr(id, s, ierr)

            call init(s)

            s% use_other_energy = .true.
            s% use_other_adjust_mdot = .true.
            
            s% extras_start_step => extras_start_step           ! Get accretion up to a given syncronization point
            s% other_adjust_mdot => other_adjust_mdot ! Calculate mdot from loaded data
            s% other_energy => other_energy          ! put energy onto the star
            s% extras_after_evolve => extras_after_evolve ! Report supernova
            ! todo: other_photo_write and other_photo_read to conserve acc_hist
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
         end subroutine extras_controls
         
         ! opens the pipes and initialize the acc_hist array
         subroutine init(s)
            use utils_lib,only:alloc_iounit
            integer :: ierr
            type(star_info) :: s
            character(64):: start_year
            integer :: num_threads
            real(8) :: age, initial_mass
            namelist /startup/ num_threads,age,initial_mass
            
            ! setup pipes
            read_unit = alloc_iounit(ierr)
            write_unit = alloc_iounit(ierr)
            open (read_unit, file="pipe1", form='formatted', &
               status='unknown', action='read')
            open(write_unit, file="pipe2", form='formatted', &
            status='unknown', action='write')

            ! read initial data
            read(read_unit,startup)
            call    OMP_SET_NUM_THREADS(num_threads)
            s%job%set_initial_age = .true.
            s%job%initial_age = age/secyer
            s%initial_mass = initial_mass/msun
            ! initialize acc_hist 
            allocate(acc_hist(1))
            acc_hist%time=age
            acc_len = 1
         end subroutine init

         integer function extras_start_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            real(dp) :: unbounded_dt = 0;
            logical :: step_forced
            step_forced = .false.
            call star_ptr(id, s, ierr)

            ! handle save point
            if (s%time >= next_model_save) then
               call star_write_model(id, next_model_filename, ierr)
               step_forced = .true.
               next_model_save = huge(0d0)
            end if 
            ! handle no more accretion data
            if (s%time >= acc_hist(acc_len)%time) then
               if (acc_len > 1) & ! write feedback if this is not first invocation
                  call write_feedback(s, .false.)
               call read_next
               step_forced = .true.
            end if
            ! Ensure resetting dt_next after artificially short steps
            if (step_forced) then
               s% dt_next = max(s% dt_next, unbounded_dt * s% max_timestep_factor)
            else
               unbounded_dt = s% dt_next
            end if
            ! Bound to dt_next to hit sync and save times
            s% dt_next = min(acc_hist(acc_len)% time-s%time, &
                             next_model_save-s%time, &
                             s% dt_next)
            extras_start_step = keep_going
         end function extras_start_step

         subroutine other_adjust_mdot(id, ierr)
            integer, intent(in):: id
            integer,intent(out) :: ierr
            real(dp) :: change
            type(star_info), pointer :: s
            call star_ptr(id, s, ierr)
            ! I use time_old here, because of some reason when there is a restart, s% time has been progressed,
            ! And s% time + s% dt is not == the next s% time
            change = get_accretion(s% time_old,  s% time_old + s% dt)
            s% mstar_dot = s% mstar_dot + change
         end subroutine other_adjust_mdot

         function get_accretion(tstart, tend)
            real(dp) :: tstart, tend, dm, lower_t, upper_t, frac, get_accretion
            integer :: i
            ! accretet mass in the given period is the sum of dm's - partial periods.
            ! TODO: this is 0'th order in rate. Possibly upgrade to 1'st order
            dm = 0
            do i = 2,  acc_len
               upper_t = min(acc_hist(i)%time, tend)
               lower_t = max(acc_hist(i-1)%time, tstart)
               frac = (upper_t-lower_t)/(acc_hist(i)%time-acc_hist(i-1)%time)
               if (frac <= 0) cycle
               dm = dm + acc_hist(i)%dm_acc * frac
            end do

            get_accretion = dm / (tend - tstart)
         end function get_accretion


         ! not mine, from sigurd
         subroutine other_energy(id, ierr)
            use star_def


            implicit none

            type (star_info), pointer :: s

            integer, intent(in) :: id
            integer, intent(out) :: ierr
            real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul
            real(dp) :: Lacc, extra_heat_var
            real(dp), parameter :: A = -1.977
            real(dp), parameter :: B = -3.078
            real(dp), parameter :: zp = log(2.792e-6)
            real(dp), parameter :: secyer = 3.1558149984d7 ! seconds per year
            real(dp), parameter :: msun   = 1.9892e33      ! solar mass (g)
            real(dp) :: omega
            integer :: k
            ierr  = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) then ! OOPS
               return
            end if

            mass = s% mstar
            mdot = s% mstar_dot

            if (mdot == 0.) return

            radius = s% r(1)
            ! ----------------------------------------------------------
            ! x_ctrl(1) is the default alpha parameter from inlist.
            !alpha = s% x_ctrl(1)

            ! ----------------------------------------------------------
            ! Time varying alpha, alpha = 0.25 if mdot > 1e-5
            ! crossover for alpha in range mdot [1e-5, 1e-6]
            ! settle at alpha = 0.01 when mdot < 1e-6
            mdot_msun = mdot*(secyer/msun)
            if (mdot_msun >= 1e-5) then
            alpha = 0.25
            else if (mdot_msun <= 1e-6) then
            alpha = 0.01
            else
            alpha = exp(B+A*tanh(zp-log(mdot_msun)))
            endif

            mcumul = 0.
            k = 0
            Lacc = 0.5 * (s% cgrav(1) * (mdot) * (mass)) / radius
            extra_heat_var = Lacc * alpha
            do while ((mcumul < mdot*s% dt .or. extra_heat_var / max(mcumul,1.) > 1e4) .and. mcumul < 0.8 * mass)
            k = k + 1
            mcumul = mcumul + s% dm(k)
            enddo

            extra_heat_var = Lacc * alpha / mcumul
            s% extra_heat(1:k) = extra_heat_var
      end subroutine other_energy

      subroutine extras_after_evolve(id, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         call get_star_ptr(id, s, ierr)
         call write_feedback(s, .true.)
      end subroutine extras_after_evolve
      
      subroutine read_next()
         logical :: save_query
         real(8) :: time, dm_acc
         character(64) :: filename
         namelist /to_mesa/ save_query, time, dm_acc, filename

         read(read_unit,to_mesa)
         if (save_query) then
            next_model_save = time
            next_model_filename = filename
            read(read_unit,to_mesa) ! only 1 save in a row
         end if

         acc_len = acc_len + 1
         ! If no more capacity in accretion_hist, resize using doubling strategy
         if (acc_len > size(acc_hist)) then
            call resize_acc_hist(size(acc_hist)*2)
         endif
         acc_hist(acc_len)%time = time
         acc_hist(acc_len)%dm_acc = dm_acc
      end subroutine read_next
      
      subroutine resize_acc_hist(newSize)
         INTEGER, INTENT(IN) :: newSize
         type(acc_struct), DIMENSION(:), ALLOCATABLE :: temp
         
         allocate(temp(size(acc_hist)))
         temp = acc_hist
         deallocate(acc_hist)
         allocate(acc_hist(newSize))
         acc_hist(1:size(temp)) = temp
      end subroutine resize_acc_hist


      subroutine write_feedback(s, supernova)
         use chem_def, only : ih1, ihe4, nuclide_not_found
         type (star_info), pointer :: s
         real(dp) :: radius, dm_wind, loum, t_eff, Zsurf
         real(dp), save :: lastmass, lasttime = 0
         logical :: supernova
         namelist /to_dispatch/ radius, dm_wind, loum, t_eff, Zsurf, supernova
         if (lasttime == 0) then
            lasttime = acc_hist(1)%time
            lastmass = s% initial_mass * Msun
         end if
         ! mass loss = wind + accretion
         radius = s% r(1)
         dm_wind = (lastmass - s% mstar) + &
            get_accretion(lasttime, s% time) * (s% time - lasttime)
         lasttime = s% time
         lastmass = s% mstar

         loum = s% L(1)
         t_eff = s%teff

         Zsurf = 1
         if (ih1 .NE. nuclide_not_found) Zsurf = Zsurf - s% xa(s% net_iso(ih1),1)
         if (ihe4 .NE. nuclide_not_found) Zsurf = Zsurf - s% xa(s% net_iso(ihe4),1)
         
         write(write_unit, to_dispatch)
      end subroutine write_feedback

      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
      
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      
         names(1) = "mdot_accretion"
         vals(1) = get_accretion(s% time-s%dt, s%time) * secyer / msun 
      
         ! wind = total mass loss + mass added
         names(2) = "mdot_wind"
         vals(2) = -s% mstar_dot * secyer / msun  + vals(1)

         ierr = 0
      end subroutine data_for_extra_history_columns
      
      end module run_star_extras
      
