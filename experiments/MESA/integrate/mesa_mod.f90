!===============================================================================
!> MESA module
!===============================================================================
MODULE mesa_mod
  use scaling_mod
  use units_mod
  USE io_mod
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod, only: trace
  USE sinks_mod
  USE gpatch_mod, only: gpatch_t
  implicit none
  private
  !
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
  integer:: verbose = 0
  integer:: cadence = 4
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
  namelist /mesa_params/ mesa_threshhold_msun, cadence, verbose
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
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mesa_t%post_update', itimer=itimer)
  call self%sinks_t%post_update (patch)
  if (self%mass > mesa_threshhold .and. patch% contains(self%position) .and. &
      .not. self%exploded) then
    if (.not.self%mesa_simulated) then
      call self% start_mesa
    else
      !!! insert criteria for when to syncronize !!!
      if (.true.) call self%eval_mesa
    end if
  end if 
  call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> Creates a working directory for the mesa instance, 
!> starts up an instance, and gives startup data
!===============================================================================
SUBROUTINE start_mesa (self)
  class(mesa_t):: self
  integer :: num_threads
  real(8) :: age, initial_mass
  namelist /startup/ num_threads,age,initial_mass
  !-----------------------------------------------------------------------------
  call trace%begin ('mesa_t%start_mesa')
  write (self%folder, "(a,i4.4,a)") trim(io%rundir)//"/star_", self%id, "/"
  if (verbose > 0) then
    write (stdout,*) trim(self%folder)//' starting'
  end if
  !
  ! system calls and newunit() are not threadsafe
  !$omp critical (mesa_cr)
  ! todo: Handle case where masa_processes is not created, or where any files are locked
  call system("bash ../mesa_template/startup.sh "//trim(self%folder)) ! start mesa instance
  open (newunit(self%write_unit), file=trim(self%folder) // 'pipe1', form='formatted', status='unknown', action='write')
  open (newunit(self%read_unit),  file=trim(self%folder) // 'pipe2', form='formatted', status='unknown', action='read')
  !
  ! Send startup data to MESA
  num_threads = 16
  age = self%time*scaling%t
  initial_mass = self%mass*scaling%m
  write (self%write_unit, startup)
  !
  self%mesa_simulated=.true.
  self%last_time = self%time
  self%dm_accreted = 0
  !$omp end critical (mesa_cr)
  !.............................................................................
  if (verbose > 1) then
    write (io_unit%mpi, startup)
  end if
  call trace%end ()
END SUBROUTINE start_mesa

!===============================================================================
!> Get a new unit for a mesa pipe. This ensures that closed units can be reused
!> Cannot be called concurrently
!===============================================================================
INTEGER FUNCTION newunit (unit)
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
END FUNCTION newunit

!===============================================================================
!> blockingly sends updated accretion data to mesa, and updates wind info
!===============================================================================
SUBROUTINE eval_mesa (self)
  class(mesa_t):: self
  real(8) :: radius, dm_wind, loum, t_eff, Zsurf, v_esc
  logical :: supernova
  character(64) :: dummy = ""
  integer :: count=0
  integer :: itimer=0
  namelist /to_dispatch/ radius, dm_wind, loum, t_eff, Zsurf, supernova
  !-----------------------------------------------------------------------------
  call trace%begin ('mesa_t%eval_mesa', itimer=itimer)
  call self%patch%lock%set ('mesa')
  !!omp critical (mesa_cr)
  if (mod(count,cadence)==0) then
    if (count == 0) &
      write (stdout,'(f12.6,2x,a)') wallclock(), 'SENDING QUERY'
    if (verbose > 0) &
      write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'SENDING QUERY'
    call send_query (self%write_unit, .false., self%time*scaling%t, &
                    self%dm_accreted*scaling%m, dummy)
    self%dm_accreted = 0
    !...........................................................................
    if (verbose > 0) &
      write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'WAITING FOR REPLY'
    if (count == 0) &
      write (stdout,'(f12.6,2x,a)') wallclock(), 'WAITING FOR REPLY'
    !...........................................................................
    read (self%read_unit, to_dispatch)
    !...........................................................................
    if (verbose > 0) &
      write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'RECEIVED REPLY'
    if (count == 0) &
      write (stdout,'(f12.6,2x,a)') wallclock(), 'RECEIVED REPLY'
    !...........................................................................
    if (supernova) then 
      self%exploded = .true.
      self%mesa_simulated = .false.
    end if 
    v_esc = sqrt(2*scaling%grav*self%mass/(radius/scaling%l))
    self%v_wind = get_vink_alpha(t_eff, Zsurf) * v_esc
    ! set mass loss rate, compensating for 
    ! excess/lack caused by non-syncronos timestepping. If all mass has 
    ! been put out in surrounding, then self% dm_wind = 0
    !$omp atomic update
    self%dm_wind = self%dm_wind - dm_wind / scaling%m
    self%mdot_wind = -self%dm_wind/(self%time - self%last_time)
    !.............................................................................
    if (verbose > 0) then
      write (io_unit%mpi,'(a,1p,9g11.3)') &
        'mesa_t%eval_mesa: t,m_sun,v_e,v,dm_sun.dm_mesa,mdot,mdot=', &
        self%time, self%mass*scaling%m/cgs%m_sun, v_esc*scaling%u/cgs%kms, &
        self%v_wind*scaling%u/cgs%kms, self%dm_wind, dm_wind/scaling%m, &
        self%mdot_wind*scaling%m/cgs%m_sun/(scaling%t/cgs%yr), &
        dm_wind/(self%time-self%last_time)/cgs%m_sun/(scaling%t/cgs%yr)
    end if
    if (verbose > 1) then
      write (io_unit%mpi, to_dispatch)
    end if
    !.............................................................................
    self%last_time = self%time
  end if
  count = count+1
  call self%patch%lock%unset ('mesa')
  !!omp end critical (mesa_cr)
  call trace%end (itimer)
END SUBROUTINE eval_mesa
!
!===============================================================================
!> Code inspired by MESA impleentation, returns v_inf/v_esc
!===============================================================================
REAL(8) FUNCTION get_vink_alpha(T, Z)
  real(8), intent(in) :: T, Z
  real(8) :: alfa, Teff_jump, dT
  real(8), parameter :: Zsolar=0.019d0
  !-----------------------------------------------------------------------------
  ! alfa = 1 for hot side, = 0 for cool side
  if (T > 27500d0) then
      alfa = 1
  else if (T < 22500d0) then
      alfa = 0
  else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
      Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
      dT = 100d0
      if (T > Teff_jump + dT) then
        alfa = 1
      else if (T < Teff_jump - dT) then
        alfa = 0
      else
        alfa = (T - (Teff_jump - dT)) / (2*dT)
      end if
  end if
  !
  ! weighted averagy from v(alfa=0)=1.3 to v(alfa=1)=2.6
  get_vink_alpha = (Z/Zsolar)**0.13d0 * (1+alfa) * 1.3d0
END FUNCTION get_vink_alpha

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
      write (filename,'("../",i5.5,"/sink_",i5.5,".data")') patch%iout, self%id;
      call send_query(self%write_unit, .true., patch% out_next*scaling%t, 0d0, filename)
    end if
  end if
  call self%sinks_t%output(patch)
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Submits the specified data to MESA, namelist can only be defined once. 
!> IN CGS, AND NOT 2 save_query IN A ROW
!===============================================================================
SUBROUTINE send_query (unit, save_query, time, dm_acc, filename)
  integer :: unit
  logical :: save_query
  real(8) :: time, dm_acc
  character(64) :: filename
  namelist /to_mesa/ save_query, time, dm_acc, filename
  !-----------------------------------------------------------------------------
  write (unit, to_mesa, DELIM="APOSTROPHE")
  !.............................................................................
  if (verbose > 1) then
    write (io_unit%mpi, to_mesa)
  end if
END SUBROUTINE send_query

END MODULE mesa_mod
