MODULE rt_internal_mod
  USE bits_mod
  USE microphysics_mod
  USE rt_mod
  USE rt_data_mod
  USE ray_geometry_mod
  USE trace_mod
  implicit none
!-----------------------------------------------------------------------
  type,extends(rt_t), public :: rt_internal_t
    type (rt_data_t),               pointer :: my_rt_data => null()
    real,         dimension(:,:),   pointer :: ds         => null()  ! (point,ray)
    real(kind=8), dimension(:,:,:), pointer :: qfw        => null()  ! (wl_bin, point, ray)
    real(kind=8), dimension(:,:,:), pointer :: qrv        => null()  ! (wl_bin, point, ray)
    real(kind=8), dimension(:,:,:), pointer :: ex0        => null()  ! (wl_bin, point, ray), first exponential coeff, for Integral solver
    integer :: solver
    procedure (solve_internal_integral), pointer :: solve => null()
  contains
    procedure :: init   => init_rt_internal
    procedure :: update => update_rt_internal
    procedure :: check  => check_rt_internal
    procedure :: solve_internal_integral
    procedure :: solve_internal_feautrier
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise tables for internal RT solution
!> =====================================================================
SUBROUTINE init_rt_internal (self)
  class(rt_internal_t) :: self
  class(ray_list_t),    pointer :: my_rays
  integer :: point, ray
  logical :: fnd
!-----------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_internal')
  call self%task_t%init
  associate (mhd_patch => self%mhd_patch)
  self%position  = mhd_patch%position
  self%offset    = mhd_patch%offset
  self%level     = mhd_patch%level
  self%velocity  = mhd_patch%velocity
  self%size      = mhd_patch%size
  self%cell_ds   = mhd_patch%ds
  self%it        = mhd_patch%it
  self%end_time  = mhd_patch%end_time
  self%out_time  = mhd_patch%out_time
  self%li        = mhd_patch%li
  self%ui        = mhd_patch%ui
  self%time      = -2
  self%dtime     = 2e-6
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
  !if (mhd_patch%is_set(bits%has_BC))    call self%set (bits%has_BC)
  call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false.)                ! retrieve ray list for this direction
  self%my_rays => my_rays
  allocate(self%name)
  self%name      = 'rt_internal'

  if (self%solver ==  0) then                                        ! set the right solver
    self%solve => solve_internal_feautrier
  else
    self%solve => solve_internal_integral
  end if

  !---------------------------------------------------------------------
  ! Allocate memory slots.
  !---------------------------------------------------------------------
  allocate(self%ds    (my_rays%max_n_points, my_rays%n_rays))
  allocate(self%ex0   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))
  allocate(self%qfw   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))
  allocate(self%qrv   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))

  do ray = 1, self%my_rays%n_rays                                       ! loop over rays to copy coords and ds
    do point = 1, self%my_rays%n_points(ray)
      self%ds(point, ray) = my_rays%ds(point, ray) * self%cell_ds(1)
    end do
  end do

  self%qfw = 0; self%qrv = 0; self%ex0 = 0

  if (io%verbose>2) call self%check
  end associate
  call trace_end
END SUBROUTINE init_rt_internal

!> =====================================================================
!> Update internal RT solution
!> =====================================================================
SUBROUTINE update_rt_internal (self)
  class(rt_internal_t) :: self
  integer :: wl, ray, point, p(3)
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_internal '//stime)

  call self%solve
  call self%update_time

  if (io%verbose>2) call self%check

  call trace_end

END SUBROUTINE update_rt_internal

!> =====================================================================
!> Solution of RTE using the Integral formulation. Forward and reverse
!> incoming radiation is ignored. dtau is stored for every point to later
!> apply the corrections. Iwf and Irv are stored in case of overlapping
!> patches.
!> =====================================================================

SUBROUTINE solve_internal_integral (self)
  class(rt_internal_t) :: self
  real(kind=8), dimension(:,:), allocatable :: ex0, ex1, ex2, ds1, d1s, d2s
  real(kind=8), dimension(:,:), allocatable :: s, dtau, idtau, ctau
  integer :: pt, ipt, it, ray, wl, p(3)
  !---------------------------------------------------------------------
  call trace_begin('rt_t%update_internal%solve_internal_integral')

  associate (npt   => self%my_rays%n_points,nwb => microphysics%n_wl_bins, &
             nrays => self%my_rays%n_rays,  rtd => self%my_rt_data, &
               qfw => self%qfw,             qrv => self%qrv, &
               cell=> self%my_rays%cell_coords, &
               mnp => self%my_rays%max_n_points)

  allocate( dtau(nwb, mnp))
  allocate(idtau(nwb, mnp))
  allocate(ex0  (nwb, mnp))
  allocate(ex1  (nwb, mnp))
  allocate(ex2  (nwb, mnp))
  allocate(ds1  (nwb, mnp))
  allocate(d1s  (nwb, mnp))
  allocate(d2s  (nwb, mnp))
  allocate(s    (nwb, mnp))
  allocate(ctau (nwb, mnp))


  do ray = 1, nrays

    do pt = 1, npt(ray)
      do wl = 1, nwb
        p = cell(:,pt,ray)
        dtau (wl,pt) = rtd%opacity(p(1),p(2),p(3),wl)*self%ds(pt,ray)
           s (wl,pt) = rtd%source (p(1),p(2),p(3),wl)
      end do
    end do

    if (npt(ray) > 1) then
      pt = npt(ray)
      idtau(:,1:pt) = 1_8 / dtau(:,1:pt)
      ex0(:,2:pt) = dexp(-dtau(:,2:pt))                                  ! first exponential coeff
      ex1(:,2:pt) = 1_8-ex0(:,2:pt)                                      ! second exponential coeff
      ex2(:,2:pt) = ex1(:,2:pt)-dtau(:,2:pt) *   ex0(:,2:pt)             ! third exponential coeff
      ds1(:,2:pt) = (s(:,2:pt) -s(:,1:pt-1)) * idtau(:,2:pt)             ! half-point centered dS/dtau
      !-----------------------------------------------------------------
      !  dS/dtau and d2S/dtau2
      !-----------------------------------------------------------------
      ctau(:,2:pt-1) = 1_8 / (dtau(:,3:pt) + dtau(:,2:pt-1))
       d1s(:,2:pt-1) = &
      (ds1(:,3:pt  ) * dtau(:,2:pt-1) + ds1(:,2:pt-1) * &
      dtau(:,3:pt) ) * ctau(:,2:pt-1)                                    ! dS/dtau, 2nd. order accurate
       d2s(:,2:pt-1) = (ds1(:,3:pt  ) - ds1(:,2:pt-1))* &
               2_8*ctau(:,2:pt-1)                                        ! d2S/dtau2, 2nd. order accurate


      d2s(:,1)  = d2s(:,2)                                               ! assume constant d2S/dtau2
      d1s(:,1)  = d1s(:,2)    - dtau(:,2)  * d2s(:,2)                    ! extrapolate to 1
      d2s(:,pt) = d2s(:,pt-1)                                            ! assume constant d2S/dtau2
      d1s(:,pt) = d1s(:,pt-1) + dtau(:,pt) * d2s(:,pt)                   ! extrapolate to n
      !-----------------------------------------------------------------
      !  Internal forward and reverse formal solutions, skipping first and last points
      !-----------------------------------------------------------------

      qfw(:,2:pt  ,ray) =-d1s(:,2:pt)  *ex1(:,2:pt)+d2s(:,2:pt)  *ex2(:,2:pt)
      qrv(:,1:pt-1,ray) = d1s(:,1:pt-1)*ex1(:,2:pt)+d2s(:,1:pt-1)*ex2(:,2:pt)

    else
      ex0(:,npt(ray)) = exp(-dtau(:,npt(ray)))                           ! exp_dtau is still needed
    end if
    self%ex0(:,:,ray) = ex0

  end do

  if (io%verbose >2) then
    if (any(isnan(self%ex0)).or.any(isnan(qfw)).or.any(isnan(qrv))) then
      error stop "internal has NaN"
    end if
  end if

  deallocate(ex0, ex1, ex2, ds1, d1s, d2s, ctau, s, dtau, idtau)
  end associate
  call trace_end
END SUBROUTINE solve_internal_integral

!> =====================================================================
!> Solution of RTE using the Feautrier formulation. Forward and reverse
!> incoming radiation is ignored.
!> =====================================================================
SUBROUTINE solve_internal_feautrier (self)
  class(rt_Internal_t) :: self
  real, dimension(:,:,:), pointer :: a1, a2, a3, d
  integer                         :: j, k, l
  real                            :: a, dsdtau_1, dsdtau_2, d2sdtau2

END SUBROUTINE solve_internal_feautrier

SUBROUTINE check_rt_internal(self)
  class(rt_internal_t) :: self
!-----------------------------------------------------------------------
  if(.not.associated(self%my_rt_data)) then
    error stop "rt_internal%check: my_rt_data is not associated"
  end if

  if(.not.associated(self%solve)) then
    error stop "rt_internal%check: solver is not set"
  end if

  if (any(isnan(self%ds))) then
    print *,'=================== ds ================'
    print *,self%ds
    error stop "rt_internal%check: self%ds has NaN"
  end if

  if (any(isnan(self%qfw))) then
    print *,'=================== ds ================'
    print *,self%qfw
    error stop "rt_internal%check: self%qfw has NaN"
  end if

  if (any(isnan(self%qrv))) then
    print *,'=================== ds ================'
    print *,self%qrv
    error stop "rt_internal%check: self%qrv has NaN"
  end if

  if (any(isnan(self%ex0))) then
    print *,'=================== ds ================'
    print *,self%ex0
    error stop "rt_internal%check: self%ex0 has NaN"
  end if

END SUBROUTINE check_rt_internal

END MODULE rt_internal_mod
