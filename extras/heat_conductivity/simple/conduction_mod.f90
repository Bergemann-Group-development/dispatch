!=======================================================================
!>
!=======================================================================
MODULE conduction_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE eos_mod
  USE scaling_mod
  USE gpatch_mod
  USE scratch_mod
  USE vector_ops
  USE timer_mod
  implicit none
  private
  type, public:: conduction_t
    real(kind=var_k):: sigma_sc = 0.0
    real(kind=var_k):: u_max = 0.0
    class(gpatch_t), pointer:: patch
    logical:: on = .false.
  contains
    procedure:: pre_init
    procedure:: courant_condition
    procedure:: pre_update
  end type
  real, save:: sigma=1e-3
  logical, save:: on = .false.
  logical, save:: first_time = .true.
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, gpatch)
  class (conduction_t):: self
  class (gpatch_t), target :: gpatch
  ! ----------------------------------------------------------------------------
  integer:: iostat
  namelist /conduction_params/ on, sigma
  !.............................................................................
  call trace%begin ('conduction_t%init')
  call io%lock%set ('conduction_t%init')
  if (first_time) then
    rewind (io_unit%input)
    read (io_unit%input, conduction_params,iostat=iostat)
    if (io%master) write (*, conduction_params)
    first_time = .false.
  end if
  call io%lock%unset('conduction_t%init')
  if (on) then
    self%patch => gpatch
    self%sigma_sc = sigma * scaling%t / scaling%l**2
    self%on = on
  end if
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Compute heat conductivity
!===============================================================================
SUBROUTINE pre_update(self)
  class (conduction_t):: self
  !.............................................................................
  real(8) :: ds(3)
  real(var_k):: vol, area
  integer :: n(3), lb(3), ub(3), i
  real(kind=var_k), dimension(:,:,:),   pointer :: d, e
  real(kind=var_k), dimension(:,:,:),   allocatable:: tg
  real(kind=var_k), dimension(:,:,:,:), allocatable:: tg_grad
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('conduction_t%pre_update', itimer, &
    detailed_timer=timer%detailed)
  !
  lb = self%patch%mesh%lb
  ub = self%patch%mesh%ub
  ds = self%patch%mesh%d
  vol = product(ds)
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  n =  shape(d)
  !-----------------------------------------------------------------------------
  ! The negative heat flux is proportional to the density times the gradient
  ! of the temperature
  !-----------------------------------------------------------------------------
  allocate (tg(n(1),n(2),n(3)))
  allocate (tg_grad(n(1),n(2),n(3),3))
  call eos%lookup_table (n, d=d, e=e, tt=tg)
  tg_grad = grad (self%patch%mesh, tg)
  do i=1,3
    tg_grad(:,:,:,i) = - tg_grad(:,:,:,i) * self%sigma_sc
  end do
  !-----------------------------------------------------------------------------
  ! The heating per unit volume is the divergence of the negative flux
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  scratch%heating_per_unit_volume = &
  scratch%heating_per_unit_volume - div (self%patch%ds, tg_grad)
  self%u_max = self%sigma_sc*(1.0/ds(1)**2 + 1.0/ds(2)**2 + 1.0/ds(3)**2)
  !
  deallocate(tg,tg_grad)
  call trace%end (itimer, detailed_timer=timer%detailed )
END SUBROUTINE pre_update

!===============================================================================
!> Intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(conduction_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('conduction_t%courant_condition')
  !print *, 'u_max(ee after spitzer)', self%patch%id, self%u_max, self%patch%u_max
  self%patch%u_max = max (self%patch%u_max, self%u_max)
  call trace%end ()
END SUBROUTINE courant_condition

END MODULE conduction_mod
