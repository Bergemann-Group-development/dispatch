MODULE rt_data_mod
  USE bits_mod
  USE io_mod
  USE microphysics_mod
  USE ray_geometry_v2_mod
  USE rt_mod
  USE trace_mod
  implicit none
  private
!-----------------------------------------------------------------------
!> Special intensity storage type. All cells can have more than one ray
!> per direction.
  !type :: int_data_t
  !  integer :: ray, point
    !integer ::lwall(3), uwall(3)
   ! real, dimension(:), allocatable :: I
    !class(int_data_t), pointer :: head => null()
    !class(int_data_t), pointer :: prev => null()
    !class(int_data_t), pointer :: next => null()
    !class(int_data_t), pointer :: tail => null()
  !contains
  !  procedure :: init => int_data_init
  !end type
!-----------------------------------------------------------------------
  type, extends(rt_t), public :: rt_data_t
    integer :: n_cells(3)
  !  integer :: n_intensity_dat = 0                                      ! used for mpi comm
    integer :: n_tau_dat = 0                                            ! used for mpi comm
    real, dimension(:,:,:,:), pointer :: opacity => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: source  => null()              ! (x,y,z,wl_bin)
    !type(int_data_t), dimension(:,:,:,:), pointer :: intensity => null()! (x,y,z,direction)
  contains
    procedure :: init   => init_rt_data
    procedure :: update => update_rt_data
    procedure :: check  => check_rt_data
    !procedure :: I_store
    !procedure :: I_retrv_cell
    !procedure :: I_retrv_ray
    procedure :: add_cr_heating
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
  integer :: i, j, k, d
  !---------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_data')
  call self%rt_t%init
  self%type = 'rt_data_t'
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
  self%time     = -3                                ! initial time
  self%dtime    = 1e-6                              ! initial dtime
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'rt_data'
  self%n_ray_dirs = ray_geometry%head%n_ray_dirs
  allocate(self%opacity(n(1), n(2), n(3), nwb)) ! opacity table
  allocate(self%source (n(1), n(2), n(3), nwb)) ! Source function table
  !allocate(self%intensity(self%ui(1),self%ui(2),self%ui(3),self%n_ray_dirs))

 ! do i = 1, self%ui(1)
 !   do j = 1, self%ui(2)
 !     do k = 1, self%ui(3)
 !       do d = 1, self%n_ray_dirs
 !         allocate (self%intensity(i,j,k,d)%I(nwb))
 !       end do
 !     end do
 !   end do
 ! end do


  self%opacity = 0;   self%source = 0
  !self%n_intensity_dat = product(n)*self%n_ray_dirs*microphysics%n_wl_bins
  end associate
  if (io%verbose>1) then
    write (*,'(a,i5,3x,a,i5,3x,a,i3)') &
      '    rt_data_t%id:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id!, &
   !   'n_intensity_dat:', self%n_intensity_dat
  end if
  if (io%verbose>2) call self%check(1)             ! check if rt_data initialized ok
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
  call trace_begin('rt_t%update_rt_data', itimer=itimer)

  self%it  = self%mhd_patch%it

  allocate (temp(size(self%source,1),size(self%source,2),size(self%source,3)))
  associate  (l => self%li, u=>self%ui)
  rho => self%mhd_patch%mem(:,:,:,1,self%it,1)

  temp(l(1):u(1),l(2):u(2),l(3):u(3)) = self%mhd_patch%gas_temperature_ngz() ! estimate temperature

  self%opacity(l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%opacity%solve (rho(l(1):u(1),l(2):u(2),l(3):u(3)),&
  temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%n_wl_bins)  !FIXME : should pass wl_bins with info ! estimate opacity
  self%source (l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%source_fn%solve (temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%wl_bins) ! estimate the source function

  ! add cosmic ray heating
  call self%add_cr_heating(self%mhd_patch%time,temp(l(1):u(1),l(2):u(2),l(3):u(3)))  ! add cr heating
  call self%update_time

  end associate
  deallocate (temp)

  if (io%verbose>2) call self%check(2)

  call trace_end (itimer)
END SUBROUTINE update_rt_data

!> =============================================================================
!> Optionally add cosmic ray heating -- generic if applied before 'apply_heating'
!> =============================================================================
SUBROUTINE add_cr_heating (self, time, temp)
  class(rt_data_t) :: self
  real, dimension(:,:,:) :: temp
  real(kind=8) :: time
  integer :: l(3), u(3)
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%add_cr_heating')
  l = self%mhd_patch%mesh%li
  u = self%mhd_patch%mesh%ui
  associate (rho => self%mhd_patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1), &
               q => self%mhd_patch%rt_q%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it))
  q = q + rho * microphysics%cr_heating
  end associate
  call trace%end
END SUBROUTINE add_cr_heating

!> =====================================================================
!> Store intensity, cell-by-cell.  This is likely to be inefficient.
!>
!> Since all rays in a ray bundle end on the same wall, and all rays
!> in a ray packed have the same length, one should be able to rip off
!> the final intensity or q-value as an array, and store it in the
!> rt_internal or rt_data task.
!>
!> Likewise, if a task contains a child task, which needs boundary
!> values, these can be computed and inserted by the child.
!> =====================================================================
!SUBROUTINE I_store (self, I_st, rays, dir_ID, rt_type)
!  class (rt_data_t) :: self
!  class (ray_t), pointer :: rays
!  type (int_data_t), pointer  :: this
!  real, dimension(:,:,:) :: I_st
!  integer :: dir_ID, rt_type, c(3), pt, ray
!  logical :: found
!  integer, save :: itimer=0
!-----------------------------------------------------------------------
!  call trace_begin('rt_data_mod%I_store',4, itimer=itimer)
!  associate (I => self%intensity(:,:,:,dir_ID),&
!             nwb => microphysics%n_wl_bins)

!  do ray = 1, rays%n_rays
!    do pt = 1, rays%n_points
!      c = rays%cell_coords(:,pt,ray)
!      if (io%verbose > 2) &
!        print '(a,2i5,5i4)', 'I_store: id, rays%id, ray, pt, c =', self%id, rays%id, ray, pt, c

!      I(c(1),c(2),c(3))%I = I_st(ray,:,pt)

!    end do
!  end do
!  end associate
!  nullify (this, rays)
!  call trace_end
!END SUBROUTINE I_store

!> =====================================================================
!> Retrieve intensity by cell
!> =====================================================================
!FUNCTION I_retrv_cell (self,c,dir) result (I_st)
!  class (rt_data_t) :: self
!  real, dimension(microphysics%n_wl_bins) :: I_st
!  integer :: dir, c(3)
! ----------------------------------------------------------------------
!  call trace_begin('rt_data_mod%I_retrv',4)
!  I_st = self%intensity(c(1),c(2),c(3),dir)%I                           ! retrieve the right intensity
!  call trace_end
!END FUNCTION I_retrv_cell

!> =====================================================================
!> Retrieve intensity by rays. Doesn't work.
!> =====================================================================
!FUNCTION I_retrv_ray (self,ray,pt,dir) result (I_st)
!  class (rt_data_t) :: self
!  type (int_data_t), pointer :: this
!  real, dimension(microphysics%n_wl_bins) :: I_st
!  integer :: dir, ray,pt, c(3)
!  logical :: found
! ----------------------------------------------------------------------
!  call trace_begin('rt_data_mod%I_retrv_ray',4)
!  associate (this => self%intensity(c(1),c(2),c(3),dir))
!  found = .false.


!  if (.not.found) error stop "ERROR: rt_data_mod%I_retrv_ray did not find data entry"
!  error stop "ERROR: I_retrv_ray in rt_data_mod doesn't work!!"
!  end associate
!  call trace_end
!END FUNCTION I_retrv_ray

!> =====================================================================
!> Initialize intensity storage for this cell
!> =====================================================================
!SUBROUTINE int_data_init(self,ray,pt)
!  class (int_data_t) :: self
!  type  (int_data_t), pointer :: new
!  integer :: ray, pt
!-----------------------------------------------------------------------
!  call trace_begin('rt_data_mod%int_data_init',4)
!  allocate (new)
!  allocate (new%I(microphysics%n_wl_bins))
!  new%ray = ray                                                         ! store ray number to ID the data
!  new%point = pt                                                        ! store point in ray to ID the data
!  if (.not.associated(self%head)) then                                  ! is it an empty cell
!    self%head => new
!    self%tail => new
!    new%head => new
!  else
!    new%head  => self%head
!    new%prev  => self%tail
!    new%prev%next => new
!    self%tail => new
!  end if
!  call trace_end
!END SUBROUTINE int_data_init

!> =====================================================================
!> Test if rt_data module has any errors
!> =====================================================================
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
