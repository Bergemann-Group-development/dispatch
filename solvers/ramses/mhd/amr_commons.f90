#define MAXLEVEL 30

module amr_commons
  use amr_parameters

  integer::mpi_comm_use                         ! communicator to nodes in use
  integer::bitwise_test=0                       ! =1 -> write, =2 -> compare
  logical::idle_nodes=.false.                   ! test for bad nodes and idle them
  logical::do_output=.true.                     ! Output just performed
  logical::output_done=.false.                  ! Output just performed
  logical::init=.false.                         ! Set up or run
  logical::balance=.false.                      ! Load balance or run
  logical::shrink=.false.                       ! Shrink mesh or run
  logical::orphan=.false.                       ! Thread with myid > ncpu2
  logical::extra_dump=.false.                   ! trigger extra dump
  logical::extra_load_balance=.false.           ! trigger extra load_balance
  logical::tight_nbor=.false.                   ! require 3 neighbors in 3rd flag step
  integer::nstep=0                              ! Time step
  integer::nstep_coarse=0                       ! Coarse step
  integer::nstep_coarse_old=0                   ! Old coarse step
  integer::nflag                                ! Refinements
  integer::ncoarse                              ! nx.ny.nz
  integer::ngrid_current                        ! Actual number of octs
  integer(8)::n_cell                               ! count cells per main step
  integer(kind=8)::allocate_size=0

  real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
  real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
  real(dp)::eint_tot=0.0D0                      ! Total internal energy
  real(dp)::epot_tot=0.0D0                      ! Total potential energy
  real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
  real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
  real(dp)::const=0.0D0                         ! Energy conservation
  real(dp)::aexp_old=1.0D0                      ! Old expansion factor
  real(dp)::rho_tot=0.0D0                       ! Mean density in the box
  real(dp)::t=0.0D0                             ! Time variable
  real(dp)::tset=0.0D0                          ! Time variable
  real(kind=8)::t8=0.0D0                        ! Time variable used internally. 64 bit.

  logical::local_debug=.false.
  real(kind=8)::dr_debug(3),r_debug(3)
  logical::dbg_grid(nvector)
  integer::dbg_ind(nvector)

  ! MPI variables
  integer::ncpu,ncpu_dump,ndomain,myid,overload=1,ncore
  integer, allocatable, dimension(:)   :: isub_domain, cpu_domain               ! Given domain, return local sub domain and cpu 
  integer, allocatable, dimension(:,:) :: map_domain                            ! Given local sub domain and cpu return domain

  ! Friedman model variables
  integer::n_frw
  real(dp),allocatable,dimension(:)::aexp_frw,hexp_frw,tau_frw,t_frw

  ! Initial conditions parameters from grafic
  integer                  ::nlevelmax_part
  real(dp)                 ::aexp_ini=10.
  real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
  real(dp),dimension(1:MAXLEVEL)::vfact,dtlev
  real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
  integer ,dimension(1:MAXLEVEL)::n1,n2,n3

  ! Level related arrays
  real(dp)                      ::dtlevelmin  ! Largest possible timestep at the coarsest level
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
  real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level

  ! Pointers for each level linked list
  integer,allocatable,dimension(:,:)::headl
  integer,allocatable,dimension(:,:)::taill
  integer,allocatable,dimension(:,:)::numbl
  integer,allocatable,dimension(:,:)::numbtot

  ! Pointers for each level boundary linked list
  integer,allocatable,dimension(:,:)::headb
  integer,allocatable,dimension(:,:)::tailb
  integer,allocatable,dimension(:,:)::numbb

  ! Pointers for free memory grid linked list
  integer::headf,tailf,numbf,used_mem,used_mem_tot

  ! Tree arrays
  real(dp),allocatable,dimension(:,:)::xg      ! grids position
  integer ,allocatable,dimension(:,:)::nbor    ! neighboring father cells
  integer ,allocatable,dimension(:)  ::father  ! father cell
  integer ,allocatable,dimension(:)  ::next    ! next grid in list
  integer ,allocatable,dimension(:)  ::prev    ! previous grid in list
  integer ,allocatable,dimension(:)  ::son     ! sons grids
  integer ,target, allocatable,dimension(:)  ::flag1   ! flag for refine
  integer ,target, allocatable,dimension(:)  ::flag2   ! flag for expansion

  ! Global indexing
  integer ,allocatable,dimension(:)  ::cpu_map  ! domain decomposition
  integer ,allocatable,dimension(:)  ::cpu_map2 ! new domain decomposition

  ! Hilbert key
  real(qdp),allocatable,dimension(:)::hilbert_key
  real(qdp),allocatable,dimension(:)::bound_key,bound_key2
  real(qdp)                         ::order_all_min,order_all_max

  ! Recursive bisection                                                                               
  real(dp),allocatable,dimension(:)    ::bisec_wall         ! bisection wall positions                
  integer ,allocatable,dimension(:,:)  ::bisec_next         ! next 2 child cells in bisection         
  integer::bisec_root                                       ! root of bisection tree                  

  integer,allocatable,dimension(:)     ::bisec_indx         ! map from leaf cell id to cpu id         
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min   ! cpu domains boxes                       
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min2  ! cpu domains boxes for new decomp        
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max2

  integer,allocatable,dimension(:)     ::bisec_cpu_load     ! CPU loads (for stats)                   
  integer,allocatable,dimension(:,:)   ::bisec_hist         ! histograms for load computation         
  integer,allocatable,dimension(:)     ::bisec_hist_bounds  ! histogram splitting boundaries          
  integer,allocatable,dimension(:)     ::new_hist_bounds
  integer,allocatable,dimension(:)     ::bisec_ind_cell     ! histo swap id -> cell id map (big)      
  integer,allocatable,dimension(:)     ::cell_level         ! store the level of the cells (big)      

  real(dp)::bisec_res                                       ! resolution parameters                   
  integer(kind=8)::bisec_nres

  ! Communication structure
  type communicator
     integer                            ::ngrid=0
     integer                            ::npart=0
     integer     ,dimension(:)  ,pointer::igrid => null()
     integer     ,dimension(:,:),pointer::f
     real(kind=8),dimension(:,:),pointer::u
     integer     ,dimension(:,:),pointer::fp
     real(kind=8),dimension(:,:),pointer::up
#ifdef ATON
     real(kind=8),dimension(:,:),pointer::u_radiation
#endif
  end type communicator

  ! Active grid, emission and reception communicators
  type(communicator),allocatable,dimension(:)  ::active
  type(communicator),allocatable,dimension(:,:)::boundary
  type(communicator),allocatable,dimension(:,:)::emission
  type(communicator),allocatable,dimension(:,:)::reception

  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro  ='hydro'
  CHARACTER(LEN=20)::type_accel  ='accel'
  CHARACTER(LEN=20)::type_flag   ='flag'

  ! Units specified by the user in the UNITS_PARAMS namelist for non-cosmo runs.
  ! These values shouldn't be used directly. Instead call units() in amr/units.f90.
  real(dp)::units_density=1.0  ! [g/cm^3]
  real(dp)::units_time=1.0     ! [seconds]
  real(dp)::units_length=1.0   ! [cm]
  real(dp)::units_velocity=0.0 ! [cm/s] (if set, replaces units_time)

  ! Counters for refinement
  integer n_refine_d, n_refine_p, n_refine_u, n_refine_c, n_refine_s, n_refine_o, &
          n_refine_m, n_refine_l, n_refine_q, n_refine_b
  integer n_gfloor_u, n_gfloor_b

  ! Type for supporting restart from arbitrary type of files
  type io_restart_t
     ! Switches detailing the input data file format
#ifndef NPRE
     logical :: double_precision_run = .false.
#else
#if NPRE==4
     logical :: double_precision_run = .false.
#else
     logical :: double_precision_run = .true.
#endif
#endif
#ifdef QUADHILBERT
     logical :: quad_precision_run = .true.
#else
     logical :: quad_precision_run = .false.
#endif
     logical :: big_endian, double_precision, quad_precision
  end type
  type(io_restart_t) :: io_restart
!
  interface load_float
    module procedure load_float_scalar, load_float_1d
  end interface
contains
!==============================================================================!
subroutine print_id (id)
  implicit none
  character(len=*) id
  !real(kind=8) wallclock
!-----------------------------------------------------------------------
  !$omp critical
  if (myid==1 .and. id .ne. ' ') then
    !print '(1x,a,f12.2)', trim(id), wallclock()
    id = ' '
  end if
  !$omp end critical
end subroutine
!==============================================================================!
subroutine load_float_scalar(ilun,x,n)
  implicit none
  integer, intent(in) :: ilun, n
  real(dp) :: x
  real(kind=8) :: x_dp
  real(kind=4) :: x_sp
  if (io_restart%double_precision) then
     read(ilun) x_dp
     x = x_dp
  else 
     read(ilun) x_sp
     x = x_sp
  end if
end subroutine load_float_scalar
!==============================================================================!
subroutine load_float_1d(ilun,x,n)
  implicit none
  integer, intent(in) :: ilun, n
  real(dp), dimension(:) :: x
  real(kind=8), dimension(:), allocatable :: x_dp
  real(kind=4), dimension(:), allocatable :: x_sp
  if (io_restart%double_precision .eqv. io_restart%double_precision_run) then  ! no need for scratch variable
     read(ilun) x(1:n)
  else
     if (io_restart%double_precision) then
        allocate(x_dp(n))
        read(ilun) x_dp
        x(1:n) = x_dp
        deallocate(x_dp)
     else 
        allocate(x_sp(n))
        read(ilun) x_sp
        x(1:n) = x_sp
        deallocate(x_sp)
     end if
  end if
end subroutine load_float_1d
!==============================================================================!
end module amr_commons
!==============================================================================!
module reduce_comm_m
  USE amr_commons
  implicit none
  ! Generic type for encapsulating information for a level-communicator
  type level_communicator_t
    logical :: inuse=.false.                                                   ! is it allocated ?
    integer :: comm                                                            ! communicator handle
    integer :: nodes                                                           ! #nodes in comm
    integer :: rank                                                            ! rank in communicator
    integer, dimension(:), allocatable :: g2l                                  ! mapping from global ranks to communicator
    integer, dimension(:), allocatable :: l2g                                  ! mapping from communicator ranks to global
    integer, dimension(:), allocatable :: nemission                            ! nr of grids to emit in global coordinates
  end type
  
  ! communicator for exchanging ghostzones between nodes on a given level
  type(level_communicator_t), dimension(:), allocatable :: full_level

  ! communicator for exchanging ghostzones between nodes that have active grids on a given level
  type(level_communicator_t), dimension(:), allocatable :: virtual_level
end module reduce_comm_m
!==============================================================================!
module openmp_support
  use amr_commons
  implicit none
  integer :: omp_nthreads, omp_mythread,mpi_omp_threads
  logical :: omp_master
  !$omp threadprivate(omp_mythread,omp_master)
  character(len=80), save:: id='amr_commons.f90 $Id$'
contains
! Function that resets the stacksize of OpenMP threads to avoid segfaults
#if defined (_OPENMP)
subroutine set_openmp_stacksize
#ifdef __INTEL_COMPILER
  use hydro_parameters, only : nvar
  use omp_lib, only : kmp_set_stacksize_s, kmp_get_stacksize_s, kmp_size_t_kind
#endif
  implicit none
#ifdef __INTEL_COMPILER
  integer(kind=kmp_size_t_kind) :: omp_stacksize, new_omp_stacksize

  omp_stacksize = kmp_get_stacksize_s()
  new_omp_stacksize = (4*1024**2 * nvar) / 8 ! Adjust according to number of passive scalars, with baseline being 4M for MHD
  if (dp==kind(1.0E0_8)) new_omp_stacksize = new_omp_stacksize * 2 ! Double for double precision runs

  if (new_omp_stacksize > omp_stacksize) then
     call kmp_set_stacksize_s(new_omp_stacksize) ! Set stacksize to needed value
     if (myid==1) write(*,*) 'WARNING! OpenMP stacksize has been reset to ', new_omp_stacksize/1024**2, 'MB to avoid segfault'
  endif
#endif
end subroutine set_openmp_stacksize
#endif
!
subroutine init_openmp
#if defined (_OPENMP)
  use omp_lib, only : omp_get_thread_num, omp_get_num_threads
#endif
  implicit none
#if defined (_OPENMP)
  ! Make sure stack size is large enough -- if compiler supports setting it on-the-fly
  call set_openmp_stacksize
  !
!$omp parallel
!$omp master
  omp_nthreads = omp_get_num_threads()
!$omp end master
  omp_mythread = omp_get_thread_num()
#else
  omp_mythread = 0
  omp_nthreads = 1
#endif
  omp_master = omp_mythread .eq. 0
!$omp end parallel
  mpi_omp_threads = (ncpu-1)/omp_nthreads+1 ! Ratio of mpi to openmp threads. For hybrid communication.
  ncore = ncpu*omp_nthreads
!
  call print_id(id)
  if (myid==1) print'(1x,a,4i8)','init_openmp: ncpu, ncore, omp_nthreads, mpi_omp_threads =', &
    ncpu, ncore, omp_nthreads, mpi_omp_threads
end subroutine init_openmp
end module openmp_support
!==============================================================================!
