MODULE rt_data_mod
  USE bits_mod
  USE microphysics_mod
  USE rt_mod
  USE trace_mod
  USE io_mod
  implicit none
!----------------------------------------------------------------------- 
  type :: int_data_t
    integer :: ray, point, lwall(3), uwall(3)
    real(kind=8), dimension(:), pointer :: I => null()
    class(int_data_t), pointer :: head => null()
    class(int_data_t), pointer :: next => null()
    class(int_data_t), pointer :: tail => null()
  contains
    procedure :: init => int_data_init
  end type
!----------------------------------------------------------------------- 
  type, extends(rt_t), public :: rt_data_t
    integer :: n_cells(3)
    real(kind=8),     dimension(:,:,:,:), pointer :: opacity => null()  ! (x,y,z,wl_bin)
    real(kind=8),     dimension(:,:,:,:), pointer :: source  => null()  ! (x,y,z,wl_bin)
    type(int_data_t), dimension(:,:,:,:), pointer :: intensity => null()! (x,y,z,direction)
  contains
    procedure :: init   => init_rt_data
    procedure :: update => update_rt_data
    procedure :: check  => check_rt_data
    procedure :: I_store
    procedure :: I_retrv
  end type
!----------------------------------------------------------------------- 
CONTAINS

!> =====================================================================
!> Initialise RT_data tables, which takes required MHD variables and gets
!> the source function and opacity per cell, mimicking the MHD patch cell
!> structure. It also stores intensity per direction per cell to communi-
!> cate it if patches overlap.
!> =====================================================================
SUBROUTINE init_rt_data(self)
  class(rt_data_t) :: self
  !---------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_data')
  call self%task_t%init
  associate (mhd_patch => self%mhd_patch, n => self%n_cells, &
                   nwb => microphysics%n_wl_bins)
  !---------------------------------------------------------------------
  ! Copy relevant data from mhd_task and it's patch to perfectly shadow it
  !---------------------------------------------------------------------
  self%li       = mhd_patch%li
  self%ui       = mhd_patch%ui
  self%n_cells  = mhd_patch%gn
  self%position = mhd_patch%position                ! task position
  self%offset   = mhd_patch%offset
  self%level    = mhd_patch%level                   ! task level
  self%velocity = mhd_patch%velocity                ! task velocity
  self%size     = mhd_patch%size                    ! task size
  self%cell_ds  = mhd_patch%ds                      ! cell size
  self%it       = mhd_patch%it                      ! task time
  self%end_time = self%mhd_patch%end_time
  self%out_time = self%mhd_patch%out_time
  self%time     = -3
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'rt_data'
  self%n_ray_dirs = ray_geometry%ray_dirs
  allocate(self%opacity(n(1), n(2), n(3), nwb)) ! opacity table
  allocate(self%source (n(1), n(2), n(3), nwb)) ! Source function table
  allocate(self%intensity(self%ui(1),self%ui(2),self%ui(3),self%n_ray_dirs))
  self%opacity = 0; self%source = 0
  end associate
  if (io%verbose>2) call self%check(1)      
  call trace_end  
END SUBROUTINE init_rt_data

!> =====================================================================
!> Update RT data task values
!> =====================================================================
SUBROUTINE update_rt_data (self)
  class(rt_data_t) :: self
  character(len=16):: stime
  real, dimension(:,:,:), pointer :: rho, temp
  integer, save:: itimer=0
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_data '//stime,itimer=itimer)
  
  self%it  = self%mhd_patch%it
  
  allocate (temp(size(self%source,1),size(self%source,2),size(self%source,3)))
  associate  (l => self%li, u=>self%ui) 
  rho => self%mhd_patch%mem(:,:,:,1,self%it,1)
  
  temp(l(1):u(1),l(2):u(2),l(3):u(3)) = self%mhd_patch%gas_temperature_ngz() 
       
  self%opacity(l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%opacity%solve (rho(l(1):u(1),l(2):u(2),l(3):u(3)),&
  temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%n_wl_bins)  !FIXME : should pass wl_bins with info
  self%source (l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%source_fn%solve (temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%wl_bins)
  
  ! add cosmic ray heating
  call self%mhd_patch%add_cr_heating(self%time,temp(l(1):u(1),l(2):u(2),l(3):u(3)))
  call self%update_time
  
  end associate
  deallocate (temp)
  
  if (io%verbose>2) call self%check(2)
    
  call trace_end (itimer)
END SUBROUTINE update_rt_data

!> =====================================================================
!> Store intensity
!> =====================================================================
SUBROUTINE I_store (self,I_st,c,pt,ray,dir,lw,uw)
  class (rt_data_t) :: self
  type (int_data_t), pointer :: this
  real (kind=8), dimension(microphysics%n_wl_bins) :: I_st
  integer :: pt, ray, dir, lw(3), uw(3), c(3)
  logical :: found
!-----------------------------------------------------------------------
  call trace_begin('rt_data_mod%I_store')
  associate (I_nt => self%intensity(c(1),c(2),c(3),dir),&
             nwb => microphysics%n_wl_bins)
  
  if (.not.associated(I_nt%head)) then
    call I_nt%init(nwb,ray,pt,lw,uw)
  else
    this  => I_nt%head
    found = .false.
    do while (associated(this%next).and.(found.eqv..false.))
      if (all(this%lwall == lw).and.all(this%uwall == uw).and.&
             (this%point == pt).and.   (this%ray   == ray)) then
        found = .true.
        this%I = I_st
      else
        this => this%next
      end if
    end do
    
    if (.not.found) then
      call I_nt%init(nwb,ray,pt,lw,uw)
      I_nt%tail%I = I_st
    end if
  end if
  end associate
  call trace_end
END SUBROUTINE I_store

!> =====================================================================
!> Retrieve intensity
!> =====================================================================
FUNCTION I_retrv (self,c,dir,lw,uw) result (I_st)
  class (rt_data_t) :: self
  type (int_data_t), pointer :: this
  real (kind=8), dimension(microphysics%n_wl_bins) :: I_st
  integer :: dir, lw(3), uw(3), c(3)
  logical :: found
! ----------------------------------------------------------------------
  call trace_begin('rt_data_mod%I_retrv')
  associate (this => self%intensity(c(1),c(2),c(3),dir)%head)
  found = .false.
  do while(associated(this).and..not.found)
    if (all(this%lwall == lw).and.all(this%uwall == uw)) then
      found = .true.
      I_st = this%I
    end if
  end do
  
  if (.not.found) error stop "ERROR: rt_data_mod%I_retrv did not find data entry"
  
  end associate
  call trace_end
END FUNCTION I_retrv

!> =====================================================================
!> Initialize intensity storage for this cell
!> =====================================================================
SUBROUTINE int_data_init(self,nwb,ray,pt,lw,uw)
  class (int_data_t) :: self
  type  (int_data_t), pointer :: new
  integer :: nwb, ray, pt, lw(3), uw(3)
!-----------------------------------------------------------------------
  call trace_begin('rt_data_mod%int_data_init')
  allocate (new)
  allocate (new%I(nwb))
  new%ray = ray
  new%point = pt
  new%lwall = lw
  new%uwall = uw
  if (.not.associated(self%head)) then
    self%head => new
    self%tail => new
    new%head => new
  else
    self%tail%next => new
    new%head  => self%head
    self%tail => new 
  end if
  call trace_end
END SUBROUTINE int_data_init

SUBROUTINE check_rt_data(self, plc)
  class(rt_data_t) :: self
  integer :: plc
!-----------------------------------------------------------------------
  if (any(isnan(self%opacity))) then
    print *,'=================== opacity ================'
    print *,self%opacity
    if (plc ==1) error stop "rt_data%init: self%opacity has NaN"
    if (plc ==2) error stop "rt_data%update: self%opacity has NaN"
  end if
  
  if (any(isnan(self%source))) then
    print *,'=================== source ================'
    print *,self%source
    if (plc ==1) error stop "rt_data%init: self%source has NaN"
    if (plc ==2) error stop "rt_data%update: self%source has NaN"
  end if
END SUBROUTINE check_rt_data

END MODULE rt_data_mod
