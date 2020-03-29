
module mesa_mod
    use io_unit_mod, only : io_unit
    use mpi_mod, only: mpi_counter_t
    implicit none
    private 
    public mesa_sink, create_sink, submit_rate

    type :: mesa_sink
      integer :: unit1, unit2 !pipes to and from mesa, respectively
      integer :: sink_id ! unique sink id across whole simulation
      real :: creation_time
    end type

    type struct_from_dispatch
      sequence
      real :: time = 0 ! time in formation of 2.nd larson core (s)
      real :: acc_mass = 0 ! accretedMass (g/s)
    end type struct_from_dispatch
    type struct_from_mesa
      sequence
      real :: time = 0 ! time in formation of 2.nd larson core (s)
      real :: wind_mass = 0 ! accretedMass (g/s)
      logical :: supernova = .false. !! booooommmm !!
    end type struct_from_mesa
    integer :: stars = 0!type(mpi_counter_t):: stars ! todo: implement so that thsi can be used without mpi

    contains
    ! Makes a folder for each mesa instance, where pipes, photos and history files are saved to.
    ! returns the mesa_sink representing the created star
    function create_sink(creation_time)
      real :: creation_time
      type(mesa_sink) :: create_sink
      character(64):: folder
      ! TODO actually implement the MPI thing
      !$omp atomic capture
        stars = stars + 1
        create_sink%sink_id = stars
      !$omp end atomic
      create_sink%creation_time = creation_time

      write (folder, "(A, I4.4, A)") "mesa_processes/star_", create_sink%sink_id, "/"

      ! system calls and newunit() is not threadsafe
      !$omp critical (mesa_init)
        call system("rm -r "  // trim(folder) //  & ! Make sure folder is empty& !";mkdir " // trim(folder) //  & ! create folder
                    ";cp -r ../mesa_template "    // trim(folder) //  & ! cd to it so that this will be the working directory
                    ";cd "    // trim(folder) //  & ! cd to it so that this will be the working directory
                    ";mkfifo pipe1 pipe2" // & ! create the two pipes
                    ";./star &") ! start mesa instance
        create_sink%unit1 = newunit()
        open (create_sink%unit1, file=trim(folder) // 'pipe1', form='formatted', status='unknown', action='write')
        create_sink%unit2 = newunit()
        open (create_sink%unit2, file=trim(folder) // 'pipe2', form='formatted', status='unknown', action='read')
      !$OMP END CRITICAL (mesa_init)
      call submit_rate(create_sink, 1e-5, creation_time)
    end function create_sink

    ! alloc_iounit
    ! submits accretion rates in cgs to MESA
    subroutine submit_rate(sink, accreted_mass, time)
      type(mesa_sink), intent(inout) :: sink
      real, intent(in) :: accreted_mass, time
      type(struct_from_dispatch) :: acc_rate
      namelist /query/ acc_rate
      acc_rate%acc_mass = accreted_mass
      acc_rate%time = time
      write (sink%unit1,*) acc_rate
    end subroutine submit_rate

    ! get a new unit for a mesa pipe
    integer function newunit()
      integer :: LUN_MIN, LUN_MAX = 300 !should there be a io_unit%mesa_max?
      logical :: opened
      
      LUN_MIN = io_unit%mesa
      do newunit=LUN_MIN,LUN_MAX
        inquire(unit=newunit,opened=opened)
        if (.not. opened) exit
      end do
      ! TODO: handle all units occupied
    end function newunit
  END
  
  
  
  