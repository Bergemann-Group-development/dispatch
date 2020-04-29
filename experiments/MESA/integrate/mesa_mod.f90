!===============================================================================
!> MESA module
!===============================================================================
MODULE mesa_mod
  use scaling_mod
  use units_mod
  USE io_mod
  USE trace_mod, only: trace
  USE sinks_mod
  USE gpatch_mod, only: gpatch_t
  implicit none
  private
  
  type, public, extends(sinks_t):: mesa_t
    integer :: read_unit, write_unit   ! Units kept open to allow communication with mesa
    logical :: mesa_simulated=.false. ! The evolution is handled by mesa
    real(8) :: last_mass, last_time      ! The state at previous submittion
    character(64)::folder=""          ! folder of pipe, used for reopening the write pipe
  contains
  procedure:: pre_init
  procedure:: post_update 

  ! private procedures:
  procedure:: start_mesa
  procedure:: send_accrection
  procedure:: recive_feedback
  end type
  real :: mesa_threshhold, mesa_min_mass
  type(mesa_t), public:: mesa
  logical:: first_time = .true.
  CONTAINS

  !===============================================================================
  !> Read settings from namelist
  !===============================================================================
  SUBROUTINE pre_init (self, patch)
    class(mesa_t):: self
    class(gpatch_t), pointer:: patch
    integer:: iostat
    real :: mesa_threshhold_msun, mesa_min_mass_msun
    namelist /mesa_params/ mesa_threshhold_msun, mesa_min_mass_msun
    call trace%begin ('mesa_t%pre_init')
    call self%sinks_t%pre_init(patch)
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, mesa_params, iostat=iostat)
      mesa_threshhold = mesa_threshhold_msun * cgs%m_sun/ scaling%m
      mesa_min_mass   = mesa_min_mass_msun   * cgs%m_sun/ scaling%m
      call io%namelist_warning ('sinks_params', iostat)
      write (io%output, mesa_params)
      write (io%namelist, mesa_params)
      first_time = .false.
    end if
    call io%lock%unset
    !-----------------------------------------------------------------------------
    call trace%end ()
  END SUBROUTINE pre_init

  SUBROUTINE post_update (self, patch)
    USE omp_timer_mod, only : wallclock
    class(mesa_t):: self
    class(gpatch_t), pointer:: patch
    integer, save :: call = 0
    real(8) :: start
    !-----------------------------------------------------------------------------
    call trace%begin ('mesa_t%post_update')
    call self%sinks_t%post_update (patch)
    if (self%mass > mesa_threshhold .and. patch%contains(self%position)) then
      if (.not.self%mesa_simulated) then
        call self%start_mesa
      else
        call self%send_accrection(.true., .false., "")
        call self%recive_feedback()
      end if
    end if 
    call trace%end ()
  end subroutine post_update

  ! create working directory fro mesa instance, start it up and connect
  ! the star is assumed to have constant accretion rate from sink creation to mesa handover.
  ! In practice, mesa cannot handle 0 mass stars. Therefore, we do not let mesa accrete the 
  ! first self%mesa_min_mass
  subroutine start_mesa(self)
    class(mesa_t):: self
    character(12):: start_year
    integer :: num_threads
    real(8) :: age, initial_mass
    namelist /startup/ num_threads,age,initial_mass

    write (self%folder, "(A, I4.4, A)") "mesa_processes/star_", self%id, "/"

  ! system calls and newunit() is not threadsafe
  !$omp critical
    ! todo: Handle case where masa_processes is not created, or where any files are locked
    call system("(rm -r "  // trim(self%folder) //  & ! Make sure folder is empty& !";mkdir " // trim(folder) //  & ! create folder
                ";cp -r ../mesa_template "    // trim(self%folder) //  & ! cd to it so that this will be the working directory
                ";cd "    // trim(self%folder) //  & ! cd to it so that this will be the working directory
                ";mkfifo pipe1 pipe2" // & ! create the two pipes
                ";./star inlist &)") ! start mesa instance
    open (newunit(self%write_unit), file=trim(self%folder) // 'pipe1', form='formatted', status='unknown', action='write')
    open (newunit(self%read_unit), file=trim(self%folder) // 'pipe2', form='formatted', status='unknown', action='read')
  !$OMP END CRITICAL
    
    ! Send startup data to MESA
    num_threads = 8
    age = self%time*scaling%t
    initial_mass = self%mass*scaling%m
    write(self%write_unit,startup)

    self%mesa_simulated=.true.
    self%last_mass = self%mass
    self%last_time = self%time
  end subroutine start_mesa

  subroutine send_accrection(self, sync_point, save_point, path)
    type acc_struct 
    sequence
      real(8) :: time, mdot
    end type acc_struct

    class(mesa_t):: self
    logical :: sync_point, save_point
    type(acc_struct) :: accretion
    character(*) :: path
    character(256) :: photo_path

    namelist /query/ accretion, sync_point, save_point, photo_path

    photo_path = path
    accretion%time = self%time*scaling%t
    accretion%mdot =  (self%mass-self%last_mass)*scaling%m/((self%time-self%last_time)*scaling%t)
    write(self%write_unit, query)

    self%last_mass = self%mass
    self%last_time = self%time
  end subroutine send_accrection

  subroutine recive_feedback(self)
    type feedback_struct
    sequence
      real :: time, dm_wind, loum
      logical :: supernova = .false. !! booooommmm !!
    end type feedback_struct  

    class(mesa_t):: self
    type(feedback_struct) :: feedback
    namelist /query/ feedback
    read(self%read_unit, query)
  end subroutine recive_feedback

  ! get a new unit for a mesa pipe. This ensures that closed units can be reused
  integer function newunit(unit)
    integer, intent(out) :: unit
    integer :: LUN_MIN, LUN_MAX = 300 !should there be a io_unit%mesa_max?
    logical :: opened

    LUN_MIN = io_unit%mesa
    do newunit=LUN_MIN,LUN_MAX
      unit = newunit
      inquire(unit=newunit,opened=opened)
      if (.not. opened) exit
    end do
    ! TODO: handle all units occupied
  end function newunit


END MODULE mesa_mod
