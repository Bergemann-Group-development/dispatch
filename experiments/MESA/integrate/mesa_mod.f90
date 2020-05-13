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
    integer :: read_unit, write_unit      ! Units kept open to allow communication with mesa
    logical :: mesa_simulated=.false.     ! Whether stellar simulation is done by MESA
    real(8) :: last_time                  ! Previously submitted data
    character(64) :: folder               ! Folder where this sink's data is located (should this just be infered?)
    logical :: exploded = .false.         ! This star has exploded
  contains
  procedure:: pre_init
  procedure:: post_update 
  procedure:: output
  
  ! private procedures:
  procedure, private :: start_mesa
  procedure, private :: eval_mesa
  end type
  real :: mesa_threshhold
  
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
    real :: mesa_threshhold_msun
    namelist /mesa_params/ mesa_threshhold_msun
    !-----------------------------------------------------------------------------

    call trace%begin ('mesa_t%pre_init')
    call self%sinks_t%pre_init(patch)
    call io%lock%set
    if (first_time) then
      rewind (io%input)
      read (io%input, mesa_params, iostat=iostat)
      mesa_threshhold = mesa_threshhold_msun * cgs%m_sun/ scaling%m
      call io%namelist_warning ('sinks_params', iostat)
      write (io%namelist, mesa_params)
      first_time = .false.
    end if
    call io%lock%unset
    !-----------------------------------------------------------------------------
    call trace%end ()

    !TODO: initialize MESA model from mesa_threshhold-deltaM with log(mdot)=-6.5
    ! and then save. Thereby the startup cost will will be severely reduced
  END SUBROUTINE pre_init

  !===============================================================================
  !> Send accretion to mesa and respond to feedback
  !===============================================================================
  SUBROUTINE post_update (self, patch)
    class(mesa_t):: self
    class(gpatch_t), pointer:: patch
    !-----------------------------------------------------------------------------
    call trace%begin ('mesa_t%post_update')
    call self%sinks_t%post_update (patch)
    if (self%mass > mesa_threshhold .and. patch% contains(self%position) .and. .not. self%exploded) then
      if (.not.self%mesa_simulated) then
        call self% start_mesa
      else
        !!! insert criteria for when to syncronize !!!
        if (.true.) call self%eval_mesa
      end if
    end if 
    call trace%end()
  end subroutine post_update
  
!===============================================================================
!> Creates a working directory for the mesa instance, 
!> starts up an instance, and gives startup data
!===============================================================================
  subroutine start_mesa(self)
    class(mesa_t):: self
    integer :: num_threads
    real(8) :: age, initial_mass
    namelist /startup/ num_threads,age,initial_mass
    !-----------------------------------------------------------------------------

    write (self%folder, "(A, I4.4, A)") "mesa_processes/star_", self%id, "/"

  ! system calls and newunit() is not threadsafe
  !$omp critical
    ! todo: Handle case where masa_processes is not created, or where any files are locked
    call system("(rm -r "  // trim(self%folder) //  & ! Make sure folder is empty& !";mkdir " // trim(folder) //  & ! create folder
                ";mkdir "  // trim(self%folder) //  & ! cd to it so that this will be the working directory
                ";cd "     // trim(self%folder) //  & ! cd to it so that this will be the working directory
                ";mkfifo pipe1 pipe2" // & ! create the two pipes
                ";../../../mesa_template/star ../../../mesa_template/inlist &)") ! start mesa instance
    open (newunit(self%write_unit), file=trim(self%folder) // 'pipe1', form='formatted', status='unknown', action='write')
    open (newunit(self%read_unit),  file=trim(self%folder) // 'pipe2', form='formatted', status='unknown', action='read')
  !$OMP END CRITICAL
    
    ! Send startup data to MESA
    num_threads = 16
    age = self%time*scaling%t
    initial_mass = self%mass*scaling%m
    write(self%write_unit,startup)

    self%mesa_simulated=.true.
    self%last_time = self%time
  end subroutine start_mesa
  
!===============================================================================
!> Get a new unit for a mesa pipe. This ensures that closed units can be reused
!> Cannot be called concurrently
!===============================================================================
  integer function newunit(unit)
    integer, intent(out) :: unit
    integer :: LUN_MIN, LUN_MAX = 300 !should there be a io_unit%mesa_max?
    logical :: opened
    !-----------------------------------------------------------------------------
    LUN_MIN = io_unit%mesa
    do newunit=LUN_MIN,LUN_MAX
      unit = newunit
      inquire(unit=newunit,opened=opened)
      if (.not. opened) exit
    end do
    ! TODO: handle all units occupied
  end function newunit

!===============================================================================
!> blockingly sends updated accretion data to mesa, and updates wind info
!===============================================================================
  subroutine eval_mesa(self)
    class(mesa_t):: self
    real(8) :: dm_wind, loum, t_eff, Zsurf
    logical :: supernova
    character(64) :: dummy = ""
    namelist /to_dispatch/ dm_wind, loum, t_eff, Zsurf, supernova
    !-----------------------------------------------------------------------------
    call send_query(self%write_unit, .false., self%time*scaling%t, self%dm_accreted*scaling%m, dummy)
    self%dm_accreted = 0

    read(self%read_unit, to_dispatch)
    if (supernova) then 
      self%exploded = .true.
      self%mesa_simulated = .false.
    end if 

    ! set mass loss rate, compensating for 
    ! excess/lack caused by non-syncronos timestepping. If all mass has 
    ! been put out in surrounding, then self% dm_wind = 0
    self%dm_wind = dm_wind - self% dm_wind / scaling%m
    self%mdot_wind = self%dm_wind/(self%time - self%last_time)

    self%last_time = self%time
  end subroutine eval_mesa

!===============================================================================
!> Output mesa task
!===============================================================================
  SUBROUTINE output (self, patch)
    class(mesa_t):: self
    class(gpatch_t), pointer:: patch
    character(64):: filename
    !-----------------------------------------------------------------------------
    call trace%begin ('mesa_t%output')
    if (self%mesa_simulated .and. patch% contains(self%position)) then
      if (self%time + self%dtime > patch% out_next) then
        write (filename,'("../../",a,i5.5,"/sink_",i5.5,".data")') trim(io%outputname), patch%iout, self%id;
        call send_query(self%write_unit, .true., patch% out_next*scaling%t, 0d0, filename)
      end if
    end if
    call self%sinks_t%output(patch)
    call trace%end
  END SUBROUTINE output

  !===============================================================================
  !> Submits the specified data to MESA, namelist can only be defined once. 
  !> IN CGS, AND NOT 2 save_query IN A ROW
  !===============================================================================
    subroutine send_query(unit, save_query, time, dm_acc, filename)
      integer :: unit
      logical :: save_query
      real(8) :: time, dm_acc
      character(64) :: filename
      namelist /to_mesa/ save_query, time, dm_acc, filename
      !-----------------------------------------------------------------------------
      write(unit, to_mesa, DELIM="APOSTROPHE")
    end subroutine  
END MODULE mesa_mod
