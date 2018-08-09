!> ============================================================================
!> Module to solve the short characteristics (SC) radiative transfer problem
!> efficiently.  The selv%vars array holds all variables in the 3-D volume,
!> indexed (i1,i2,i3,i4,i5), where
!>
!>   i1-i3 = the horizontal spatial indices
!>      i4 = wavelength bin
!>      i5 = 1: opacity, 2: source function, 3: net specific heating (I-S)
!>
!> The layers (scratch) array has one more index:
!>
!>      i6 = the three layers before, at, and after i5
!>
!> The propagation of EOS information is as follows:
!>
!> eos%lookup_table (d=d, e=e, ...)            rt_patch_t%vars(:,:,:,:,1:2)
!> rt_sc_t%prepare_layers                      rt_sc_t%vars
!>   rt_scr_t%bilinear                         rt_sc_t%layers
!> rt_sc_t%update
!>   rt_sc_t%layer_solver                      rt_sc_t%layers
!> rt_patch_t%update                           rt_patch_t%vars(:,:,:,:,3)
!> experiment_t%update                         patch_t%mem(:,:,:,idx%qr,2)
!> ============================================================================
MODULE rt_sc_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE mesh_mod
  USE link_mod
  USE task_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public:: omega_t
    integer:: id
    real:: mu, phi, w
  end type
  !-----------------------------------------------------------------------------
  type, public:: rt_sc_t
    real(8):: position(3)
    type(omega_t):: omega
    type(link_t), pointer:: fw_nbor => null()
    type(link_t), pointer:: rv_nbor => null()
    integer:: id
    integer:: n_bin
    integer:: o(2)
    integer:: verbose=0
    real:: p(2)
    real(8):: time=0.0
    logical:: done=.true.
    class(mesh_t), pointer:: mesh(:)
    real, dimension(:), pointer:: w_bin
    real, dimension(:,:,:), pointer:: q
    real, dimension(:,:,:,:,:), pointer:: vars
    real, dimension(:,:,:,:,:,:), pointer:: layers
  contains
    procedure:: init
    procedure:: update
    procedure:: layer_solver
    procedure:: prepare_layers
    procedure:: calculate_offsets
    procedure:: bilinear
    procedure:: needs
  end type
CONTAINS

!> ============================================================================
!> Read parameters
!> ============================================================================
SUBROUTINE init (self)
  class (rt_sc_t):: self
  integer:: n_bin
  !-----------------------------------------------------------------------------
  ! Allocate memory for all variables and for three layers of interpolated vars
  ! The variables are: rk, s, q for each wavelength
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%init')
  !                                nx,              ny,              nz,   nlam, ksq
  allocate (self%q   (self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn))
  call io%bits_mem (storage_size(self%q), product(self%mesh%gn))
  allocate (self%w_bin (self%n_bin))
  self%w_bin = 1.0/self%n_bin
  call trace%end
END SUBROUTINE init

!> ============================================================================
!> Integral SC/RT solver
!> ============================================================================
SUBROUTINE update (self, dir)
  class (rt_sc_t):: self
  integer:: dir, i3, i_bin
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: rk
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%update', itimer=itimer)
  if (self%done) then
    !-----------------------------------------------------------------------------
    ! Decide if we have enough upstream info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! Allocate scratch space, and compute three depth layers of everything
    !-----------------------------------------------------------------------------
    allocate (self%layers(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,self%n_bin,3,3))
    call io%bits_mem (storage_size(self%layers), product(self%mesh%gn)*3*3)
    call self%prepare_layers
  end if
  !-----------------------------------------------------------------------------
  ! Forward-index solution
  !-----------------------------------------------------------------------------
  if (dir > 0) then
    do i3=self%mesh(3)%li,self%mesh(3)%ui
      call self%layer_solver (self%layers(:,:,:,:,:,1), &   ! k = opacity
                              self%layers(:,:,:,:,:,2), &   ! s = source
                              self%layers(:,:,:,:,:,3), &   ! q = input
                              self%vars  (:,:,:,:  ,3), &   ! q = solution
                              i3, 1, 3)
    end do
  else
  !-----------------------------------------------------------------------------
  ! Reverse-index solution
  !-----------------------------------------------------------------------------
    do i3=self%mesh(3)%ui,self%mesh(3)%li,-1
      call self%layer_solver (self%layers(:,:,:,:,:,1), &   ! k = opacity
                              self%layers(:,:,:,:,:,2), &   ! s = source
                              self%layers(:,:,:,:,:,3), &   ! q = input
                              self%vars  (:,:,:,:  ,3), &   ! q = solution
                              i3, 3, 1)
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Check if both directions done, and integrate over bins if so
  !-----------------------------------------------------------------------------
  if (self%done) then
    self%done = .false.
  else
    self%q = 0.0
    do i_bin=1,self%n_bin
      self%q = self%q + self%vars(:,:,:,i_bin,3)*self%w_bin(i_bin)
    end do
    deallocate (self%layers)
    call io%bits_mem (-storage_size(self%layers), product(self%mesh%gn)*3*3)
    self%done = .true.
  end if
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Integral solver for one layer.  The exact integal solution Q=I-S for a 2nd
!> order polynomial has coefficient for the upstream Q, and for the 1st and 2nd
!> derivative of S.  However, if one is doing a solution both forwards and 
!> backward, the term in the 1st derivative cancels.
!> =============================================================================
SUBROUTINE layer_solver (self, k, s, qin, qout, i3, l1, l3)
  class (rt_sc_t):: self
  real, dimension(:,:,:,:,:) :: k, s, qin
  real, dimension(:,:,:,:)   :: qout
  integer:: i3, l1, l3, iv
  !.............................................................................
  integer:: l(3), u(3), il, i1, i2
  real:: dz1, dz2, dtau1, dtau2, dsdtau1, dsdtau2, d2sdtau2, ex0, ex2
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%integral_forward', itimer=itimer)
  dz2 = self%mesh(3)%d/2.0
  l = self%mesh(1:3)%li
  u = self%mesh(1:3)%ui
  z => self%mesh(3)%r
  dz1      = abs(z(l1)-z(2))*0.5
  dz2      = abs(z(2)-z(l3))*0.5
  do il=1,self%n_bin
    do i2=l(2),u(2)
      do i1=l(1),u(2)
        dtau1    = dz1*(k(i1,i2,i3,il, 2) + k(i1,i2,i3,il,l1))
        dtau2    = dz2*(k(i1,i2,i3,il,l3) + k(i1,i2,i3,il, 2))
        dsdtau1  =     (s(i1,i2,i3,il, 2) - s(i1,i2,i3,il,l1))/dtau1
        dsdtau2  =     (s(i1,i2,i3,il,l3) - s(i1,i2,i3,il, 2))/dtau2
        d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1+dtau2)
        ex0      = exp(-dtau1)
        ex2      = (1.0-ex0) - dtau1*ex0
        qout(i1,i2,i3,il) = qin(i1,i2,i3,il,l1)*ex0 &
                          + d2sdtau2           *ex2
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE layer_solver

!> ============================================================================
!> Prepare three layers; one centered on the mesh points, one above and below
!> The indices in the layers are (i1,i2,i3,il,iv,1:3), where i1-i3 are 3D coords,
!> ib is bin index, iv is variable index (k,s,q), and the last is layer index
!> ============================================================================
SUBROUTINE prepare_layers (self)
  class (rt_sc_t):: self
  integer:: ib, iv, l(3), u(3), i3
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: rk
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%prepare_layers', itimer=itimer)
  l = self%mesh(1:3)%li
  u = self%mesh(1:3)%ui
  !-----------------------------------------------------------------------------
  rk  => self%vars(:,:,:,:,1)
  do ib=1,self%n_bin
    do iv=1,3
      do i3=2,self%mesh(3)%gn-1
        call self%calculate_offsets (i3, i3-1)
        call self%bilinear (self%vars(:,:,i3,ib,iv), self%layers(:,:,i3,ib,1,iv))
        !-------------------------------------------------------------------------
        self%layers(:,:,i3,ib,2,iv) = self%vars(:,:,i3,ib,iv)
        !-------------------------------------------------------------------------
        call self%calculate_offsets (i3, i3+1)
        call self%bilinear (self%vars(:,:,i3,ib,iv), self%layers(:,:,i3,ib,3,iv))
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE prepare_layers

!> ============================================================================
!> Calculate offsets from one layer to the next
!> ============================================================================
SUBROUTINE calculate_offsets (self, i3a, i3b)
  class (rt_sc_t):: self
  integer:: i3a, i3b
  integer, save:: itimer=0
  real(8), pointer:: z(:)
  real(8):: dx, dy, dz
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%calculate_offsets', 2, itimer=itimer)
  z => self%mesh(3)%r
  dz = z(i3b)-z(i3a)
  dx = self%mesh(1)%d
  dy = self%mesh(2)%d
  self%p = dz*self%omega%mu*[cos(self%omega%phi)/dx,sin(self%omega%phi)/dy]
  self%o = floor(self%p)
  self%p = self%p - self%o
  call trace%end (itimer)
END SUBROUTINE calculate_offsets

!> ============================================================================
!> Bilinear interpolation in a plane
!> ============================================================================
SUBROUTINE bilinear (self, in, out)
  class (rt_sc_t):: self
  real, dimension(:,:), intent(in)  :: in
  real, dimension(:,:), intent(out) :: out
  real:: p(2), q(2)
  integer:: i1, i2, i3, l(3), u(3), o(2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%bilinear', 2, itimer=itimer)
  l = self%mesh%li
  u = self%mesh%ui
  o = self%o
  p = self%p
  q = 1.0-p
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    out(i1,i2) = &
      p(1)*p(2)*in(i1+o(1)  ,i2+o(2)  ) + &
      q(1)*p(2)*in(i1+o(1)+1,i2+o(2)  ) + &
      p(1)*q(2)*in(i1+o(1)+1,i2+o(2)+1) + &
      q(1)*q(2)*in(i1+o(1)+1,i2+o(2)+1)
  end do
  end do
  call trace%end (itimer)
END SUBROUTINE bilinear

!> ============================================================================
!> Determine if a task is needed (is upstream)
!> ============================================================================
FUNCTION needs (self, task, dir)
  logical:: needs
  class(rt_sc_t):: self
  class(task_t), pointer:: task
  integer:: dir
  real:: ds1(3), ds2(3), sinmu
  real, parameter:: eps=1e-4
  integer:: i, is1(3), is2(3)
  !----------------------------------------------------------------------------
  needs = .true.
  ds1 = (self%mesh%p-task%position)/self%mesh%b
  sinmu = sqrt(1d0-self%omega%mu**2)
  ds2 = dir*[sinmu*cos(self%omega%phi),sinmu*sin(self%omega%phi),self%omega%mu]
  is1 = 0
  is2 = 0
  do i=1,3
    if (ds1(i) > +eps) is1(i)=1
    if (ds2(i) > +eps) is2(i)=1
    if (ds1(i) < -eps) is1(i)=-1
    if (ds2(i) < -eps) is2(i)=-1
  end do
  !-----------------------------------------------------------------------------
  ! For purely vertical rays, only patches below or above are needed
  !-----------------------------------------------------------------------------
  if (self%omega%mu==1.0) then
    do i=1,3
      if (is1(i)/=is2(i)) then
        needs = .false.
        exit 
      end if
    end do
  !-----------------------------------------------------------------------------
  ! For general directions, a patch "on the wrong side" is not needed.  This
  ! requires that it is not at the same position in that direction, and is on
  ! the opposite side, relative to the ray direction
  !-----------------------------------------------------------------------------
  else
    do i=1,3
      if (is1(i) /= 0 .and. is1(i) /= is2(i)) then
        needs = .false.
        exit 
      end if
    end do
  end if
  if (self%verbose>0) then
    if (needs) print &
    '("RT task",i5,i4,4x,"needs patch",i5,4x,"rel pos:",3f8.3,4x,"ray dir:",3f8.3,6i3)', &
    self%id, dir*self%omega%id, task%id, ds1, ds2, is1, is2
   end if
END FUNCTION

END MODULE rt_sc_mod
