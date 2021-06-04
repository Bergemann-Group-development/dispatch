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
    real:: kappa
    real:: stencil = 1.8
    logical:: kinematic = .false.
    real:: u_max = 0.0
    real:: courant = 1.0
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: courant_condition
    procedure:: pre_update
  end type
  real,    save:: mult(2) = 1.0
  logical, save:: extrapolate = .false.
  logical, save:: debug = .false.
  logical, save:: on = .false.
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, gpatch)
  class (conduction_t):: self
  class (gpatch_t), target :: gpatch
  ! ----------------------------------------------------------------------------
  integer:: iostat
  real, save:: stencil = 2.0
  logical, save:: first_time = .true.
  namelist /conduction_params/ on, extrapolate, mult, stencil
  !.............................................................................
  call trace%begin ('conduction_t%init')
  call io%lock%set ('conduction_t%init')
  if (first_time) then
    rewind (io_unit%input)
    read (io_unit%input, conduction_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io_unit%input, conduction_params)
    end if
    call io%namelist_warning ('conduction_params', iostat)
    !---------------------------------------------------------------------------
    ! The compact stencil cannot be used in the kinematic case
    !---------------------------------------------------------------------------
    if (self%kinematic) then
      stencil = 2.0
    end if
    write (io%namelist, conduction_params)
    first_time = .false.
  end if
  call io%lock%unset('conduction_t%init')
  if (on) then
    self%patch => gpatch
    self%stencil = stencil
    !---------------------------------------------------------------------------
    ! Conductivity needs 2 guard zones to get the good heating rates for
    ! li:ui range. However, for Riemann solvers, one innermost guard zone
    ! needs good values as well. Thus there are two options:
    ! either extrapolate the temperature gradient, or use 3 guard zones
    !---------------------------------------------------------------------------
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
  real(var_k):: vol, area, del1, del2, c1(3), c2(3)
  integer :: i1, i2, i3, n(3), lb(3), ub(3), i, l(3), u(3)
  real(kind=var_k), dimension(:,:,:),   pointer :: d, e, s, px, py, pz
  real(kind=var_k), dimension(:,:,:),   pointer:: tg
  real(kind=var_k), dimension(:,:,:,:), pointer:: flux
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('conduction_t%pre_update', itimer=itimer)
  !
  lb = self%patch%mesh%lb
  ub = self%patch%mesh%ub
  ds = self%patch%mesh%d
  vol = product(ds)
  d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  if (self%patch%idx%e>0) e => self%patch%mem(:,:,:,self%patch%idx%e,self%patch%it,1)
  if (self%patch%idx%s>0) s => self%patch%mem(:,:,:,self%patch%idx%s,self%patch%it,1)
  px => self%patch%mem(:,:,:,self%patch%idx%px,self%patch%it,1)
  py => self%patch%mem(:,:,:,self%patch%idx%py,self%patch%it,1)
  pz => self%patch%mem(:,:,:,self%patch%idx%pz,self%patch%it,1)
  n =  shape(d)
  ! ----------------------------------------------------------------------------
  ! get temperature
  ! ----------------------------------------------------------------------------
  allocate (tg(n(1),n(2),n(3)))
  l = self%patch%mesh%lb
  u = self%patch%mesh%ub
  if (self%patch%idx%e>0) then
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          tg(i1,i2,i3) = (e(i1,i2,i3) - 0.5*(px(i1,i2,i3)**2 + &
                                            py(i1,i2,i3)**2 + &
                                            pz(i1,i2,i3)**2)/d(i1,i2,i3)) &
                      * (eos%gamma-1d0) / d(i1,i2,i3)
        end do
      end do
    end do
  end if
  if (self%patch%idx%s>0) then
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          tg(i1,i2,i3) = d(i1,i2,i3)**eos%gamma*&
                   exp(s(i1,i2,i3)/d(i1,i2,i3)*(eos%gamma-1d0))/d(i1,i2,i3)
        end do
      end do
    end do
  end if

  !-----------------------------------------------------------------------------
  ! The negative heat flux is proportional to (the density times) the gradient
  ! of the temperature, for dynamic (kinematic) formulations
  !-----------------------------------------------------------------------------
  allocate (flux(n(1),n(2),n(3),3))
  if (self%stencil /= 1.0) then
    l = self%patch%mesh%lb+1
    u = self%patch%mesh%ub-1
    c1 = self%kappa*0.5/ds
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          flux(i1,i2,i3,1) = (tg(i1+1,i2,i3) - tg(i1-1,i2,i3))*c1(1)
          flux(i1,i2,i3,2) = (tg(i1,i2+1,i3) - tg(i1,i2-1,i3))*c1(2)
          flux(i1,i2,i3,3) = (tg(i1,i2,i3+1) - tg(i1,i2,i3-1))*c1(3)
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    ! Note than only the parallel direction is used in the divergence, so there
    ! is no need to extrapolate all components in all directions
    !---------------------------------------------------------------------------
    if (extrapolate) then
      flux(l(1)-1,:,:,1) = flux(l(1),:,:,1)
      flux(u(1)+1,:,:,1) = flux(u(1),:,:,1)
      flux(:,l(2)-1,:,2) = flux(:,l(2),:,2)
      flux(:,u(2)+1,:,2) = flux(:,u(2),:,2)
      flux(:,:,l(3)-1,3) = flux(:,:,l(3),3)
      flux(:,:,u(3)+1,3) = flux(:,:,u(3),3)
    end if
    if (self%kinematic) then
      do i1=1,3
        flux(:,:,:,i1) = flux(:,:,:,i1) * d
      end do
    end if
  end if
!-------------------------------------------------------------------------------
  ! The heating per unit volume is the negative divergence of the negative flux.
  ! Here we linearly combine a 1 and 2 delta stencil approximations when using
  ! per unit mass conductivity, while for now we use a purely 2 delta stencil
  ! when using per unit volume conductivity
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%heating_per_unit_volume, self%patch%mesh%gn)
  if (self%kinematic) then
    if (extrapolate) then
      l = self%patch%mesh%lb+1
      u = self%patch%mesh%ub-1
    else
      l = self%patch%mesh%lb+2
      u = self%patch%mesh%ub-2
    end if
    c1 = 0.5/ds
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      del1 = (flux(i1+1,i2,i3,1) - flux(i1-1,i2,i3,1))*c1(1) &
           + (flux(i1,i2+1,i3,2) - flux(i1,i2-1,i3,2))*c1(2) &
           + (flux(i1,i2,i3+1,3) - flux(i1,i2,i3-1,3))*c1(3)
      scratch%heating_per_unit_volume(i1,i2,i3) = &
      scratch%heating_per_unit_volume(i1,i2,i3) + del1
    end do
    end do
    end do
  else
    c1 = self%kappa/ds**2
    c2 = 0.5/ds
    if (self%stencil /= 2.0) then
      l = self%patch%mesh%lb+1
      u = self%patch%mesh%ub-1
      do i3=l(3),u(3)
      do i2=l(2),u(2)
      do i1=l(1),u(1)
        del1 = (tg(i1-1,i2,i3) -2.0*tg(i1,i2,i3) + tg(i1+1,i2,i3))*c1(1) &
             + (tg(i1,i2-1,i3) -2.0*tg(i1,i2,i3) + tg(i1,i2+1,i3))*c1(2) &
             + (tg(i1,i2,i3-1) -2.0*tg(i1,i2,i3) + tg(i1,i2,i3+1))*c1(3)
        scratch%heating_per_unit_volume(i1,i2,i3) = &
        scratch%heating_per_unit_volume(i1,i2,i3) + (2.0-self%stencil)*del1
      end do
      end do
      end do
    end if
    if (self%stencil /= 1.0) then
      l = self%patch%mesh%lb+2
      u = self%patch%mesh%ub-2
      do i3=l(3),u(3)
      do i2=l(2),u(2)
      do i1=l(1),u(1)
        del2 = (flux(i1+1,i2,i3,1) - flux(i1-1,i2,i3,1))*c2(1) &
             + (flux(i1,i2+1,i3,2) - flux(i1,i2-1,i3,2))*c2(2) &
             + (flux(i1,i2,i3+1,3) - flux(i1,i2,i3-1,3))*c2(3)
        scratch%heating_per_unit_volume(i1,i2,i3) = &
        scratch%heating_per_unit_volume(i1,i2,i3) + (self%stencil-1.0)*del2
      end do
      end do
      end do
    end if
  end if
  if (debug .and. self%patch%id==1) then
    i1 = self%patch%mesh(1)%li
    i2 = self%patch%mesh(2)%li
    do i3=l(3),u(3)
      print *, i3, flux(i1,i2,i3,3), tg(i1,i2,i3), &
        scratch%heating_per_unit_volume(i1,i2,i3)
    end do
  end if
  !-----------------------------------------------------------------------------
  ! The expression below assumes that all ds are the same.  The factor of 3 is
  ! for 2*ds stencils, while the factor 12 is for 1*ds stencils.
  !-----------------------------------------------------------------------------
  l = self%patch%mesh%li
  u = self%patch%mesh%ui
  self%u_max = self%kappa*minval(ds)/self%courant &
             * (mult(1)*4.*(2.0-self%stencil) + mult(2)*(self%stencil-1.0)) &
             * (1./ds(1)**2 + 1./ds(2)**2 + 1./ds(3)**2) &
             / minval(d(l(1):u(1),l(2):u(2),l(3):u(3)))
  deallocate (tg, flux)
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
  if (.not.on) return
  call trace%begin ('conduction_t%courant_condition')
  if (first_time) then
    call io%lock%set
    if (first_time .and. self%u_max > 0.0) then
      print *, 'cond id,kappa,u_max', &
        self%patch%id, self%kappa, self%u_max, self%patch%u_max
      first_time = .false.
    end if
    call io%lock%unset
  end if
  self%patch%u_max = max (self%patch%u_max, self%u_max)
  call trace%end ()
END SUBROUTINE courant_condition

END MODULE conduction_mod
