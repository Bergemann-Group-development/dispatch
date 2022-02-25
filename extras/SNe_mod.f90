!===============================================================================
!> Supernovae module, which keeps track of the lifetime, and time since birth
!> of massive stars that will eventually go supernova.  It needs these procs:
!>
!>     init: Read tables with SN life time as a function of mass
!> lifetime: Interpolate in lifetime table
!>   update: Update the time and decide whether to explode
!>      add: Add a sink in the mass ranke to the SN list
!>  explode: impart the SN energy onto a patch that the SN overlaps with
!>  courant: Impose a timestep condition on overlapping patches
!===============================================================================
MODULE SNe_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE units_mod
  USE scaling_mod
  USE lagrange_mod
  USE list_mod
  USE bits_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! SN table
  !-----------------------------------------------------------------------------
  integer:: n_table
  type sn_table_t
    real:: mass, lifetime
  end type
  type(sn_table_t), allocatable:: sn_table(:)
  !-----------------------------------------------------------------------------
  ! Supernova data type -- extendsion of doubly-linked list node
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: sn_t
    real(8):: position(3)
    real(8):: start, time, lifetime, mass
  contains
  end type
  !-----------------------------------------------------------------------------
  ! Linked list of -- as of yet unexploded -- supernovae
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: sne_t
    class(list_t), pointer:: task_list => null()
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure, nopass:: lifetime
    procedure:: pre_update
    procedure:: update
    procedure:: courant
    procedure, nopass:: SN_test
  end type
  logical:: on=.false.
  integer:: verbose=0
  real:: mass_limit
  real:: remnant_msun=1.2
  character(len=16):: file='Schallerlife.tbl'
  type(sne_t):: sn_list
  type(sne_t), public:: sne
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(sne_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer:: iostat
  logical:: first_time=.true.
  namelist /sne_params/ on, verbose, file
  !-----------------------------------------------------------------------------
  call trace%begin ('sne_t%pre_init')
  self%patch => patch
  call io%lock%set ('sne_t%pre_init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, sne_params, iostat=iostat)
    call io%namelist_warning ('sne_params', iostat)
    call read_table
    call SN_test
  end if
  call io%lock%unset ('sne_t%pre_init')
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self, name)
  class(sne_t):: self
  character(len=*), optional:: name
  !.............................................................................
  integer:: iostat
  logical:: first_time=.true.
  namelist /sne_params/ on, verbose, file
  !-----------------------------------------------------------------------------
  call trace%begin ('sne_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, sne_params, iostat=iostat)
    call read_table
    call SN_test
  end if
  call io%lock%unset
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize the SN lifetime table
!===============================================================================
SUBROUTINE read_table
  integer:: iostat, n, i
  real:: mass, lifetime
  !-----------------------------------------------------------------------------
  call trace%begin ('sne_t%read_table')
  write (stdout,*) 'SNe_t%read_table: file =', 'input_data/'//trim(file)
  open (unit=io_unit%tmp, file=trim(io%top)//'data/SNe/'//trim(file), form='formatted', status='old')
  !-----------------------------------------------------------------------------
  ! Read data from simple text file, with mass(M_Sun), lifetime(Myr)
  !-----------------------------------------------------------------------------
  read (io_unit%tmp,*) n_table
  allocate (sn_table(n_table))
  do i=1,n_table
    read (io_unit%tmp, *) mass, lifetime
    sn_table(i)%mass = mass*cgs%m_sun/scaling%m
    sn_table(i)%lifetime = lifetime*1e6*cgs%yr/scaling%t
  end do
  mass_limit = sn_table(1)%mass
  close (io_unit%tmp)
  call trace%end()
END SUBROUTINE read_table

!===============================================================================
!> Initialize task list pointer
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(sne_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('sne_t%init_task_list')
  self%task_list => task_list
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Interpolate in the table with SN lifetimes
!===============================================================================
REAL FUNCTION lifetime (mass)
  real:: mass
  !.............................................................................
  call trace%begin ('sne_t%lifetime')
  lifetime = lagrange%interpolate1d (mass, sn_table%mass, sn_table%lifetime)
  call trace%end()
END FUNCTION lifetime

!===============================================================================
!> Timestep condition imposed on a patch (selected for overlapping).  Make sure
!> that the patch does not evolve beyond the next SN explosion time
!===============================================================================
SUBROUTINE courant (self)
  class(sne_t):: self
  class(dll_node_t), pointer:: sn
  !.............................................................................
  call trace%begin ('sne_t%courant')
  sn => sn_list%head
  do while (associated(sn))
    select type (sn)
    class is (sn_t)
    if (self%patch%contains (sn%position)) then
      select type (sn)
      class is (sn_t)
      self%patch%dtime = min(self%patch%dtime, &
                             sn%start+sn%lifetime - self%patch%time)
      end select
    end if
    end select
    sn => sn%next
  end do
  call trace%end()
END SUBROUTINE courant

!===============================================================================
!> Update and manage explosion. This gets called for all MHD patches
!===============================================================================
SUBROUTINE pre_update (self)
  class(sne_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sne_t%pre_update')
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> Add a star to the list of future SNe, if the mass is large enough. This gets
!> called every time a sink patch updates, so should also be used to update the
!> SN position and time
!===============================================================================
SUBROUTINE update (self, mass, task, sn)
  class(sne_t):: self
  class(task_t), pointer:: task
  real(4):: mass
  type(sn_t), pointer:: sn
  class(dll_node_t), pointer:: dll_node
  !.............................................................................
  if (.not.on) &
    return
  call trace%begin ('sne_t%update')
  if (verbose > 1) &
    print *, 'SNe_t%update: mass, mass_limit =', mass, mass_limit
  select type (task)
  class is (gpatch_t)
  if (associated(sn)) then
    sn%mass = mass
    sn%time = task%time
    sn%position = task%position
    if (verbose > 0) &
      write (stdout,'(a,i6,3f12.6,2x,3f12.6)') &
        'SNe_t%update: id, time, remains, mass, pos =', &
        sn%id, sn%time, sn%start+sn%lifetime-sn%time, sn%mass, sn%position
    if (sn%time >= sn%start+sn%lifetime) then
      call explode (sn)
      mass = remnant_msun*cgs%m_sun/scaling%m
      call task%set (bits%no_download)
      deallocate (sn)
      nullify (sn)
    end if
  !-----------------------------------------------------------------------------
  ! If the mass exceeds the limiting mass, create a sn_t data type and append
  ! to the sne_t list
  !-----------------------------------------------------------------------------
  else if (mass > mass_limit) then
    allocate (sn)
    sn%id = task%id
    sn%mass = mass
    sn%start = task%time
    sn%position = task%position
    sn%lifetime = lifetime (mass)
    dll_node => sn
    call sne%append (dll_node)
    if (verbose > 0) &
      write(stdout,*) 'SNe_t%update: id, mass, M_sun, lifetime, Myr =', sn%id, &
        sn%mass, scaling%m_sun(sn%mass), sn%lifetime, scaling%myr(sn%lifetime)
  end if
  end select
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Impart the explosion onto a patch
!==============================================================================
SUBROUTINE explode (sn)
  type(sn_t), pointer:: sn
  !.............................................................................
  call trace%begin ('sne_t%update')
  if (verbose > 0) &
    write (stdout,*) &
      'SNe_t%explode: time, mass[M_Sun] =', sn%time, scaling%m_sun(sn%mass)
  write (stderr,*) '***********************************************************'
  write (stderr,*) '************************* BOOOM! **************************'
  write (stderr,*) '***********************************************************'
  call trace%end()
END SUBROUTINE explode

!==============================================================================
!> Test functionality
!==============================================================================
SUBROUTINE SN_test
  real:: mass, life
  !----------------------------------------------------------------------------
  if (io%master) then
    mass = 6.0
    do while (mass < 100.)
      life = lifetime (real(mass*cgs%m_sun/scaling%m,kind=4))
      write(stdout,*) 'SNe_t%test: mass, lifetime =', &
        mass, real(life*scaling%t/(1e6*cgs%yr))
      mass = mass*1.5
    end do
  end if
END SUBROUTINE SN_test

END MODULE SNe_mod
