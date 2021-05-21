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
    real:: kappa = 0.07
    real:: u_max = 0.0
    real:: stencil = 1.5
    real:: courant = 1.0
    class(gpatch_t), pointer:: patch
    logical:: kinematic = .false.
    logical:: order2 = .false.
    logical:: on = .false.
  contains
    procedure:: pre_init
    procedure:: courant_condition
    procedure:: pre_update
  end type
  real(8), save:: start_time=0d0
  logical, save:: extrapolate=.false.
  real,save:: mult=1.0
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, gpatch)
  class (conduction_t):: self
  class (gpatch_t), target :: gpatch
  ! ----------------------------------------------------------------------------
  integer:: iostat
  real,save::    kappa = 0.0
  logical,save:: on    = .false.
  logical,save:: first_time = .true.
  namelist /conduction_params/ on, kappa, start_time, extrapolate, mult
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
    ! --------------------------------------------------------------------------
    ! Conductivity needs 2 guard zones to get the good heating rates for 
    ! li:ui range. However, for Riemann solvers, one innermost guard zone
    ! needs good values as well. Thus there are two options:
    ! either extrapolate the temperature gradient, or use 3 guard zones
    ! --------------------------------------------------------------------------
    if (.not.extrapolate) self%patch%ng=max(self%patch%ng,3)
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
  real(var_k):: vol, area, c1, c2, del1, del2
  integer :: i1, i2, i3, n(3), lb(3), ub(3), i, l(3), u(3)
  real(kind=var_k), dimension(:,:,:),   pointer :: d, e, px, py, pz, e_th
  real(kind=var_k), dimension(:,:,:),   allocatable:: tg
  real(kind=var_k), dimension(:,:,:,:), allocatable:: tg_grad
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  !if (self%patch%time<=start_time) return
  call trace%begin ('conduction_t%pre_update', itimer=itimer)
  !
  lb = self%patch%mesh%lb
  ub = self%patch%mesh%ub
  ds = self%patch%mesh%d
  vol = product(ds)
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  px => self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1)
  py => self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1)
  pz => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1)
  n =  shape(d)
  ! ----------------------------------------------------------------------------
  ! get internal energy
  ! ----------------------------------------------------------------------------
  allocate (e_th(n(1),n(2),n(3)))
  l = self%patch%mesh%lb
  u = self%patch%mesh%ub
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        e_th(i1,i2,i3) = e(i1,i2,i3) - 0.5 * (px(i1,i2,i3)**2 + &
                                              py(i1,i2,i3)**2 + &
                                              pz(i1,i2,i3)**2) / d(i1,i2,i3)
      end do
    end do
  end do

  !-----------------------------------------------------------------------------
  ! The negative heat flux is proportional to the density times the gradient
  ! of the temperature
  !-----------------------------------------------------------------------------
  allocate (tg(n(1),n(2),n(3)))
  allocate (tg_grad(n(1),n(2),n(3),3))
  call eos%lookup_table (n, d=d, e=e_th, tt=tg)
  !l = self%patch%li
  !u = self%patch%ui
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1)+1,u(1)-1
        tg_grad(i1,i2,i3,1) = tg(i1+1,i2,i3) - tg(i1-1,i2,i3)
      end do
      tg_grad(l(1),i2,i3,1) = tg_grad(l(1)+1,i2,i3,1)
      tg_grad(u(1),i2,i3,1) = tg_grad(u(1)-1,i2,i3,1)
    end do
  end do
  !
  do i3=l(3),u(3)
    do i2=l(2)+1,u(2)-1
      do i1=l(1),u(1)
        tg_grad(i1,i2,i3,2) = tg(i1,i2+1,i3) - tg(i1,i2-1,i3)
      end do
    end do
    tg_grad(:,l(2),i3,2) = tg_grad(:,l(2)+1,i3,2)
    tg_grad(:,u(2),i3,2) = tg_grad(:,u(2)-1,i3,2)
  end do
  !
  do i3=l(3)+1,u(3)-1
    do i2=l(2)+1,u(2)
      do i1=l(1)+1,u(1)
        tg_grad(i1,i2,i3,3) = tg(i1,i2,i3+1) - tg(i1,i2,i3-1)
      end do
    end do
  end do
  tg_grad(:,:,l(3),3) = tg_grad(:,:,l(3)+1,3)
  tg_grad(:,:,u(3),3) = tg_grad(:,:,u(3)-1,3)

  if (self%kinematic) then
    do i1=1,3
      tg_grad(:,:,:,i1) = tg_grad(:,:,:,i1) * d
    end do
  end if

!-----------------------------------------------------------------------------
  ! The heating per unit volume is the negative divergence of the negative flux.
  ! Here we linearly combine a 1 and 2 delta stencil approximations when using
  ! per unit mass conductivity, while for now we use a purely 2 delta stencil
  ! when using per unit volume conductivity
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  c2 = self%kappa*(eos%gamma-1d0)/ds(1)**2
  c1 = c2*0.25
  if (extrapolate) then
    l = self%patch%mesh%lb+1
    u = self%patch%mesh%ub-1
  else
    l = self%patch%mesh%lb+2
    u = self%patch%mesh%ub-2
  end if
  if (self%kinematic) then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      del1 = tg_grad(i1+1,i2,i3,1) - tg_grad(i1-1,i2,i3,1) &
           + tg_grad(i1,i2+1,i3,2) - tg_grad(i1,i2-1,i3,2) &
           + tg_grad(i1,i2,i3+1,3) - tg_grad(i1,i2,i3-1,3)
      scratch%heating_per_unit_volume(i1,i2,i3) = &
      scratch%heating_per_unit_volume(i1,i2,i3) + c1*del1
    end do
    end do
    end do
  else
    c1 = c1*(2.0-self%stencil)*0.5
    c2 = c2*(self%stencil-1.0)
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      del1 = tg_grad(i1+1,i2,i3,1) - tg_grad(i1-1,i2,i3,1) &
           + tg_grad(i1,i2+1,i3,2) - tg_grad(i1,i2-1,i3,2) &
           + tg_grad(i1,i2,i3+1,3) - tg_grad(i1,i2,i3-1,3)
      del2 = tg(i1-1,i2,i3) -6.0*tg(i1,i2,i3) + tg(i1+1,i2,i3) &
           + tg(i1,i2-1,i3)                   + tg(i1,i2+1,i3) &
           + tg(i1,i2,i3-1)                   + tg(i1,i2,i3+1)
      scratch%heating_per_unit_volume(i1,i2,i3) = &
      scratch%heating_per_unit_volume(i1,i2,i3) + c1*del1 + c2*del2
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  self%u_max = (c1/12. + c2/6.)*mult
  !
  deallocate (tg, tg_grad, e_th)
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> Intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(conduction_t):: self
  logical, optional:: detailed_timer
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('conduction_t%courant_condition')
  if (first_time) then
    call io%lock%set
    if (first_time .and. self%u_max > 0.0) then
      print *, 'cond u_max', self%patch%id, self%u_max, self%patch%u_max
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%patch%u_max = max (self%patch%u_max, self%u_max)
  call trace%end ()
END SUBROUTINE courant_condition

END MODULE conduction_mod
