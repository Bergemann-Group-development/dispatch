MODULE microphysics_mod
  USE bits_mod
  USE eos_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mesh_mod
  USE mhd_mod
  USE omp_mod
  USE opacity_mod
  USE patch_mod
  USE ray_geometry_v3_mod
  USE rt_boundary_mod
  USE rt_data_mod
  USE scaling_mod
!  USE source_function_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod  
  implicit none
  private
  type, extends(patch_t), public :: microphysics_t
    class(mhd_t),  pointer :: mhd_patch  => null()
    class(link_t), pointer :: mhd_link   => null()
    class(rt_data_t), pointer :: rt_data => null()                      ! pointer to rt_data object
    logical :: rt = .false.                                             ! do rt?
    logical :: chemistry = .false.                                      ! do chemistry?
    logical :: store_rt  = .false.                                      ! if rt is done, should intensities be stored?
    real, dimension(:,:,:,:,:), pointer :: intensity => null()          ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:),   pointer :: q         => null()          ! (x,y,z,nt)
    real(kind=8) :: cr_heating
  contains
    procedure :: init   => init_mp
    procedure :: update => update_mp
    procedure :: apply_heating_rates
    procedure :: add_cr_heating
    procedure :: store_mean_I
  end type
  
CONTAINS

!> =====================================================================
!> Initialize a single forge_t instance for a single mhd patch
!> =====================================================================
SUBROUTINE mp_setup(self)
  class(link_t) :: self
END SUBROUTINE mp_setup


!> =====================================================================
!> Initialize a forge_t object, whatever it means
!> =====================================================================
SUBROUTINE init_mp(self)
  class(microphysics_t) :: self
  !---------------------------------------------------------------------
  call trace_begin('microphysics_t%init_mp')
  call self%patch_t%init
  self%type = 'microphysics_t'
  associate (mhd_patch => self%mhd_patch, nwb => opacity%nwb, n => self%gn)
  !---------------------------------------------------------------------
  ! Copy relevant data from mhd_task and it's patch to perfectly shadow it
  !---------------------------------------------------------------------
  self%li       = mhd_patch%li
  self%ui       = mhd_patch%ui
  self%gn       = mhd_patch%gn
  self%position = mhd_patch%position                                    ! task position
  self%offset   = mhd_patch%offset
  self%level    = mhd_patch%level                                       ! task level
  self%velocity = mhd_patch%velocity                                    ! task velocity
  self%size     = mhd_patch%size                                        ! task size
  self%ds       = mhd_patch%ds                                          ! cell size
  self%it       = mhd_patch%it                                          ! task time slot
  self%nt       = mhd_patch%nt
  self%end_time = self%mhd_patch%end_time
  self%out_time = self%mhd_patch%out_time
  self%time     = -3                                                    ! initial time
  self%dtime    = 1e-6                                                  ! initial dtime
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'microphysics'
  if (self%rt) then
    allocate(self%q      (n(1), n(2), n(3), self%nt))                     ! net heating
    if (self%store_rt)  then
      allocate(self%intensity (n(1), n(2), n(3), nwb, self%nt)) 
      self%intensity = 0
    end if
  end if
  end associate
END SUBROUTINE init_mp


!> =====================================================================
!> Update a forge_t object, whatever it means
!> =====================================================================
SUBROUTINE update_mp(self)
  class(microphysics_t) :: self
  !---------------------------------------------------------------------
  call self%add_cr_heating                                              ! add cosmic ray heating
  if (associated(self%rt_data)) then
    deallocate(self%rt_data%opacity,self%rt_data%source)                ! free some space
  end if
END SUBROUTINE update_mp

!> =============================================================================
!> Optionally add cosmic ray heating -- generic if applied before 'apply_heating'
!> =============================================================================
SUBROUTINE add_cr_heating (self)
  class(microphysics_t) :: self
  integer :: l(3), u(3)
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%add_cr_heating')
  l = self%mhd_patch%mesh%li
  u = self%mhd_patch%mesh%ui
  associate (rho => self%mhd_patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1), &
               q => self%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it))
  q = q + rho * self%cr_heating
  end associate
  call trace%end
END SUBROUTINE add_cr_heating

!> =====================================================================
!> Apply the heating rates to the MHD patch
!> =====================================================================
SUBROUTINE apply_heating_rates(self)
  class(microphysics_t) :: self
END SUBROUTINE apply_heating_rates

!> =====================================================================
!> Store intensity
!> =====================================================================
SUBROUTINE store_mean_I (self,I_st,rays,dirw)
  class (microphysics_t) :: self
  class (ray_t), pointer :: rays
  real, dimension(:,:,:) :: I_st
  real :: dirw
  integer :: c(3), pt, ray
!-----------------------------------------------------------------------
  call trace_begin('rt_data_mod%I_store',4)
  associate (I => self%intensity(:,:,:,:, self%it),&
             nwb => opacity%nwb)
             
  
  do ray = 1, rays%n_rays
    do pt = 1, rays%n_points
      c = rays%coords(:,pt,ray)
      I(c(1),c(2),c(3),:) = I(c(1),c(2),c(3),:) + I_st(ray,:,pt)*dirw
    end do
  end do
 ! print *,I(:,:,:,:)
  end associate

  call trace_end
END SUBROUTINE store_mean_I

END MODULE microphysics_mod
