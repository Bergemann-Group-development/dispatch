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
    real(kind=var_k):: kappa = 0.07
    real(kind=var_k):: u_max = 0.0
    class(gpatch_t), pointer:: patch
    logical:: on = .false.
  contains
    procedure:: pre_init
    procedure:: courant_condition
    procedure:: pre_update
  end type
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, gpatch)
  class (conduction_t):: self
  class (gpatch_t), target :: gpatch
  ! ----------------------------------------------------------------------------
  integer:: iostat
  real::    kappa = 0.0
  logical:: on    = .false.
  logical:: first_time = .true.
  namelist /conduction_params/ on, kappa
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
    self%kappa = kappa
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
  real(var_k):: vol, area, c, del
  integer :: i1, i2, i3, n(3), lb(3), ub(3), i, l(3), u(3)
  real(kind=var_k), dimension(:,:,:),   pointer :: d, e
  real(kind=var_k), dimension(:,:,:),   allocatable:: tg
  real(kind=var_k), dimension(:,:,:,:), allocatable:: tg_grad
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('conduction_t%pre_update', itimer=itimer)
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
  !-----------------------------------------------------------------------------
  ! The heating per unit volume is the divergence of the negative flux
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  c = self%kappa * (eos%gamma-1d0)
  l = self%patch%li
  u = self%patch%li
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    del = tg(i1-1,i2,i3) -6.0*tg(i1,i2,i3) + tg(i1,i2,i3+1) &
        + tg(i1,i2-1,i3)                   + tg(i1,i2+1,i3) &
        + tg(i1,i2,i3-1)                   + tg(i1+1,i2,i3)
    scratch%heating_per_unit_volume(i1,i2,i3) = &
    scratch%heating_per_unit_volume(i1,i2,i3) + c*del
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  self%u_max = c/6.0 &
             * (1.0/ds(1)**2 + 1.0/ds(2)**2 + 1.0/ds(3)**2)
  !
  deallocate(tg,tg_grad)
  call trace%end (itimer)
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
