      module run_star_extras

         use star_lib
         use star_def
         use const_def
         
         implicit none       
         
         !TODO: put this in ../com_structs so that it can be compiled from both MESA and Dispatch, and compatibility is maintained by the compiler
         ! communication between MESA and Dispatch is done in cgs
         type acc_struct 
         sequence
            real(8) :: time, mdot
         end type acc_struct

         type feedback_struct
         sequence
            real :: time, dm_wind, loum
            logical :: supernova = .false. !! booooommmm !!
         end type feedback_struct  
         
         TYPE(acc_struct), DIMENSION(:), allocatable :: acc_hist
         TYPE(integer) :: acc_len = 0
         character(256) :: photo_path
         logical :: save_point = .false.
         integer :: read_unit, write_unit

         contains        
         
         subroutine extras_controls(id, ierr)
            integer :: id, ierr
            type (star_info), pointer :: s

            call star_ptr(id, s, ierr)

            call init(s)

            s% use_other_energy = .true.
            s% use_other_adjust_mdot = .true.
            
            s% extras_finish_step => other_after_step           ! Get accretion up to a given syncronization point
            s% other_adjust_mdot => other_dispatch_adjust_mdot ! Calculate mdot from loaded data
            s% other_energy => accretion_other_energy          ! put energy onto the star
            ! todo: other_photo_write and other_photo_read to conserve acc_hist
         end subroutine extras_controls
         
         ! opens the pipes and initialize the acc_hist array
         subroutine init(s)
            use utils_lib, only : alloc_iounit
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
            call read_next ! read first
         end subroutine init

         integer function other_after_step(id, id_extra)
            integer, intent(in) :: id, id_extra
            integer :: ierr
            type (star_info), pointer :: s

            call star_ptr(id, s, ierr)
            if (s% time >= acc_hist(acc_len)%time) then
               ! When we get here, the timestep was usually decreased in order to hit. This is not neccecary
               s% dt_next = max(s% dt_next,s%dt_old*1.1)

               call write_feedback(0d0, s% L(1), s% time + s% dt_next, .false.)

               if (save_point)  call star_save_for_restart(id, photo_path, ierr)
               
               call read_next ! read next batch of inputs
            else 

            end if
            s% dt_next = min(acc_hist(acc_len)% time-s%time, s% dt_next)
            
            other_after_step = keep_going
         end function other_after_step

         subroutine other_dispatch_adjust_mdot(id, ierr)
            integer, intent(in):: id
            integer,intent(out) :: ierr
            type(star_info), pointer :: s
            call star_ptr(id, s, ierr)

            s% mstar_dot = acc_rate(s% time, s% time + s% dt)
         end subroutine other_dispatch_adjust_mdot
         
         ! calulate acc_rate assuming lineara accretion rate between subsequent timestamps
         function acc_rate(tstart, tend)
            real(dp) :: tstart, tend, dm, acc_rate, lower_t, upper_t
            integer:: i

            ! accretet mass in the given period is the sum of dm's - partial periods.
            ! TODO: this is 0'th order in rate. Possibly upgrade to 1'st order
            dm = 0
            do i = acc_len,  2, -1
               lower_t = max(acc_hist(i-1)%time, tstart)
               upper_t = min(acc_hist(i)%time, tend)
               if (lower_t >= upper_t) cycle
                  dm = dm + acc_hist(i)%mdot*(upper_t-lower_t)
            end do

            acc_rate = dm / (tend - tstart)
         end function acc_rate


         ! blocking read of inputs untill syncronization point
         ! returns time of syncronization point
         subroutine read_next()
            type(acc_struct) :: accretion
            logical :: sync_point
            namelist /query/ accretion, sync_point, save_point, photo_path
            
      100   read (read_unit,query)
            ! if add accretion to list
            acc_len = acc_len + 1
            ! If no more capacity in accretion_hist, resize using doubling strategy
            if (acc_len > size(acc_hist)) then
               call resize_acc_hist(size(acc_hist)*2)
            endif
            acc_hist(acc_len) = accretion
            ! continue untill all queries are read
            if (.not.sync_point) goto 100
         end subroutine read_next

         subroutine write_feedback(dm_wind, loum, time, supernova)
            type(feedback_struct) :: feedback
            real(dp) :: dm_wind, loum, time
            logical :: supernova
            namelist /query/ feedback

            feedback%dm_wind = dm_wind
            feedback%loum = loum
            feedback%time = time
            feedback%supernova = supernova

            write(write_unit, query)
         end subroutine write_feedback


         subroutine resize_acc_hist(newSize)
            INTEGER, INTENT(IN) :: newSize
            type(acc_struct), DIMENSION(:), ALLOCATABLE :: temp
            
            allocate(temp(size(acc_hist)))
            temp = acc_hist
            deallocate(acc_hist)
            allocate(acc_hist(newSize))
            acc_hist(1:size(temp)) = temp
         end subroutine resize_acc_hist


         
         subroutine accretion_other_energy(id, ierr)

            use crlibm_lib, only: safe_log10_cr
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
            !write(*,*) "Mdot (msun)", mdot_msun
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

            !omega = sqrt((s% cgrav(1) * mass)/(2*radius)**3)
            !s% omega(1:k) = omega
            !write(*,*) "Omega: ", omega
            !write(*,*) "k = ", k
            !write(*,*) "Radius(k): " , s% r(k)
            !write(*,*) "Ratio, r(k) / r(1), Heat mass / mass : " , s% r(k) / radius, mcumul / mass
            !write(*,*) "Extra heat: " , extra_heat_var

      end subroutine accretion_other_energy

      
      end module run_star_extras
      
