!===============================================================================
!> Add radiative heat diffusion to any solver.
!===============================================================================
MODULE rad_diffusion_mod
  USE io_mod
  USE os_mod
  USE link_mod
  USE gpatch_mod
  USE mesh_mod
  USE io_unit_mod
  USE trace_mod
  USE bits_mod
  USE scratch_mod
  USE eos_mod
  implicit none
  private
  type, public:: rad_diffusion_t
    logical:: on=.false.
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
    procedure:: courant_condition
  end type
  integer, save:: verbose=0
  real, save:: dtau_min=3
  logical, save:: check_dtau=.false.
CONTAINS

!===============================================================================
!> Initialize Newton cooling
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(rad_diffusion_t):: self
  class(gpatch_t), target:: patch
  ! ............................................................................
  integer:: iostat, m(3)
  logical, save:: on=.false.
  logical, save:: first_time=.true.
  namelist /rad_diffusion_params/ on, verbose, check_dtau, dtau_min
  !-----------------------------------------------------------------------------
  call trace%begin ('rad_diffusion_t%pre_init')
  self%patch => patch
  call io%lock%set ('rad_diffusion_t%pre_init')
  if (first_time) then
    rewind (io_unit%input)
    read (io%input, rad_diffusion_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io%input, rad_diffusion_params)
    end if
    call io%namelist_warning ('rad_diffusion_params', iostat)
    write (io%namelist, rad_diffusion_params)
    first_time = .false.
  end if
  call io%lock%unset ('rad_diffusion_t%pre_init')
  self%on = on
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize Newton cooling, allocating heating array knowing patch%gn
!===============================================================================
SUBROUTINE init (self)
  class(rad_diffusion_t):: self
  ! ............................................................................
  integer:: m(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('rad_diffusion_t%init')
  ! nothing to do here...
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Apply Newton cooling
!===============================================================================
SUBROUTINE pre_update (self)
  class(rad_diffusion_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  real, dimension(:,:,:), pointer:: d, e, s, kp, kr, Src, Q
  real(8), parameter:: c=1d0/3d0, c2=log(1d4)
  real:: dtau1, dtau2, dtau, qd,min_dtau
  real(8):: ds(3)
  integer:: ix, iy, iz, lb(3), ub(3), n(3)
  ! ----------------------------------------------------------------------------
  if (.not. self%on) return
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  if (trim(self%patch%etype)=='thermal') then
    e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  else if (trim(self%patch%etype)=='entropy') then
    s => self%patch%mem(:,:,:,self%patch%idx%s,self%patch%it,1)
  else
    call io%abort('rad_diffusion_mod:: total energy not set up')
  end if
  lb = self%patch%mesh%lb+1
  ub = self%patch%mesh%ub-1
  ds = self%patch%mesh%d
  ! ----------------------------------------------------------------------------
  ! Allocate scratches
  ! ----------------------------------------------------------------------------  
  n=self%patch%mesh%gn
  call os%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  Q => scratch%heating_per_unit_volume
  allocate(kr(n(1),n(2),n(3)))
  allocate(kp(n(1),n(2),n(3)))
  allocate(Src(n(1),n(2),n(3)))

  ! ----------------------------------------------------------------------------
  ! Lookup opacities and source functions
  ! ----------------------------------------------------------------------------
  call eos%lookup(e=e,d=d,var=kr,name='kr')
  call eos%lookup(e=e,d=d,var=Src,name='Src')
  ! ----------------------------------------------------------------------------
  ! Compute dtau and heating rates
  ! Range is mesh%lb+1 and mesh%ub-1 to account for shifts.
  ! ----------------------------------------------------------------------------
  if (check_dtau) then
    min_dtau=9999.0
    do iz=lb(3),ub(3)
      do iy=lb(2),ub(2)
        do ix=lb(1),ub(1)
          ! --- x ---
          dtau1 = 0.5*(kr(ix-1,iy,iz)+kr(ix,iy,iz))*ds(1)
          dtau2 = 0.5*(kr(ix+1,iy,iz)+kr(ix,iy,iz))*ds(1)
          dtau  = 0.5*(dtau1+dtau2)
          min_dtau=min(dtau,min_dtau)
          qd = &
             + ((Src(ix+1,iy,iz)-Src(ix  ,iy,iz))/dtau2 &
              - (Src(ix  ,iy,iz)-Src(ix-1,iy,iz))/dtau1)/dtau*c
          ! --- y ---
          dtau1 = 0.5*(kr(ix,iy-1,iz)+kr(ix,iy,iz))*ds(2)
          dtau2 = 0.5*(kr(ix,iy+1,iz)+kr(ix,iy,iz))*ds(2)
          dtau  = 0.5*(dtau1+dtau2)
          min_dtau=min(dtau,min_dtau)
          qd = qd &
             + ((Src(ix,iy+1,iz)-Src(ix,iy  ,iz))/dtau2 &
              - (Src(ix,iy  ,iz)-Src(ix,iy-1,iz))/dtau1)/dtau*c
          ! --- z ---
          dtau1 = 0.5*(kr(ix,iy,iz-1)+kr(ix,iy,iz))*ds(3)
          dtau2 = 0.5*(kr(ix,iy,iz+1)+kr(ix,iy,iz))*ds(3)
          dtau  = 0.5*(dtau1+dtau2)
          min_dtau=min(dtau,min_dtau)
          qd = qd &
             + ((Src(ix,iy,iz+1)-Src(ix,iy,iz  ))/dtau2 &
              - (Src(ix,iy,iz  )-Src(ix,iy,iz-1))/dtau1)/dtau*c
          ! --- combine ---
          Q(ix,iy,iz) = Q(ix,iy,iz) + qd
        end do
      end do
    end do
    if (min_dtau<=dtau_min) print *, 'WARNING: dtau leq dtau_min', min_dtau
  else
    do iz=lb(3),ub(3)
      do iy=lb(2),ub(2)
        do ix=lb(1),ub(1)
          ! --- x ---
          dtau1 = 0.5*(kr(ix-1,iy,iz)+kr(ix,iy,iz))*ds(1)
          dtau2 = 0.5*(kr(ix+1,iy,iz)+kr(ix,iy,iz))*ds(1)
          dtau  = 0.5*(dtau1+dtau2)
          qd = &
             + ((Src(ix+1,iy,iz)-Src(ix  ,iy,iz))/dtau2 &
              - (Src(ix  ,iy,iz)-Src(ix-1,iy,iz))/dtau1)/dtau*c
          ! --- y ---
          dtau1 = 0.5*(kr(ix,iy-1,iz)+kr(ix,iy,iz))*ds(2)
          dtau2 = 0.5*(kr(ix,iy+1,iz)+kr(ix,iy,iz))*ds(2)
          dtau  = 0.5*(dtau1+dtau2)
          qd = qd &
             + ((Src(ix,iy+1,iz)-Src(ix,iy  ,iz))/dtau2 &
              - (Src(ix,iy  ,iz)-Src(ix,iy-1,iz))/dtau1)/dtau*c
          ! --- z ---
          dtau1 = 0.5*(kr(ix,iy,iz-1)+kr(ix,iy,iz))*ds(3)
          dtau2 = 0.5*(kr(ix,iy,iz+1)+kr(ix,iy,iz))*ds(3)
          dtau  = 0.5*(dtau1+dtau2)
          qd = qd &
             + ((Src(ix,iy,iz+1)-Src(ix,iy,iz  ))/dtau2 &
              - (Src(ix,iy,iz  )-Src(ix,iy,iz-1))/dtau1)/dtau*c
          ! --- combine ---
          Q(ix,iy,iz) = Q(ix,iy,iz) + qd
        end do
      end do
    end do
  end if
  ! ----------------------------------------------------------------------------
  ! Cleanup
  ! ----------------------------------------------------------------------------
  deallocate(kr, kp, Src)
    
END SUBROUTINE pre_update

SUBROUTINE courant_condition(self)
  class(rad_diffusion_t):: self
END SUBROUTINE courant_condition

END MODULE rad_diffusion_mod
