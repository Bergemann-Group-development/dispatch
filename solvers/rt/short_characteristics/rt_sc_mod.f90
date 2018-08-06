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
!> ============================================================================
MODULE rt_sc_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public:: mesh_t
    integer:: gn, ng, n, li, ui
    real(8):: d=1.
  end type
  !-----------------------------------------------------------------------------
  type, public:: omega_t
    real:: mu, phi, w
  end type
  !-----------------------------------------------------------------------------
  type, public:: rt_sc_t
    type(omega_t):: omega
    integer:: nl, n_bin
    integer:: o(2)
    real:: p(2)
    real(8):: time
    logical:: double_solution=.true.
    type(mesh_t), pointer:: mesh(:)
    real, dimension(:), allocatable:: w_bin
    real, dimension(:,:,:), allocatable:: q
    real, dimension(:,:,:,:,:), allocatable:: vars
    real, dimension(:,:,:,:,:,:), allocatable:: layers
    real(8), dimension(:), allocatable:: upstream_times
  contains
    procedure:: init
    procedure:: test
    procedure:: update
    procedure:: layer_solver
    procedure:: prepare_layers
    procedure:: calculate_offsets
    procedure:: bilinear
    procedure:: eos_volume
  end type
CONTAINS

!> ============================================================================
!> Read parameters
!> ============================================================================
SUBROUTINE init (self)
  class (rt_sc_t):: self
  integer:: n=32, nl=8
  integer:: ng=2, iostat
  namelist /test_params/ n, ng, nl
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, test_params, iostat=iostat)
  if (io%master) write (io%output, test_params)
  self%nl = nl
  !-----------------------------------------------------------------------------
  ! Make a mesh
  !-----------------------------------------------------------------------------
  allocate (self%mesh(3))
  self%mesh%n  = n
  self%mesh%ng = ng
  self%mesh%gn = n+2*ng
  self%mesh%li = 1+ng
  self%mesh%ui = n+ng
  !-----------------------------------------------------------------------------
  ! Allocate memory for all variables and for three layers of interpolated vars
  ! The variables are: rk, s, q for each wavelength
  !-----------------------------------------------------------------------------
  !                                nx,              ny,              nz,   nlam, ksq
  allocate (self%q   (self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn))
  allocate (self%vars(self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn, self%nl,  3))
  self%vars = 1.0
  call self%update
  call timer%print
END SUBROUTINE init

!> ============================================================================
!> Setup a test and run it
!> ============================================================================
SUBROUTINE test (self)
  class (rt_sc_t):: self
  real, dimension(:,:,:,:), allocatable:: in, out
  integer:: i1, i2, i3, il, iv
  !----------------------------------------------------------------------------
  allocate (in (self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,self%nl))
  allocate (out(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,self%nl))
  self%p = 0.5
  call timer%tic
  do il=1,self%nl
  do i2=-1,1
  do i1=-1,1
    self%o = [i1,i2]
    call self%bilinear (in(:,:,:,il), out(:,:,:,il))
  end do
  end do
  end do
  call timer%toc ('layers', self%nl*9*product(self%mesh%n))
  deallocate (in, out)
END SUBROUTINE test

!> ============================================================================
!> Integral SC/RT solver
!> ============================================================================
SUBROUTINE update (self)
  class (rt_sc_t):: self
  integer:: i3, i_bin
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Decide if we have enough upstream info
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Compute the opacity and source function for all bins
  !-----------------------------------------------------------------------------
  call self%eos_volume
  !-----------------------------------------------------------------------------
  ! Allocate scratch space, and compute three depth layers of everything
  !-----------------------------------------------------------------------------
  allocate (self%layers(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn,self%nl,3,3))
  call self%prepare_layers
  !-----------------------------------------------------------------------------
  ! Forward-index solution
  !-----------------------------------------------------------------------------
  do i3=self%mesh(3)%li,self%mesh(3)%ui
    call self%layer_solver (self%layers(:,:,:,:,:,1), &   ! k = opacity
                            self%layers(:,:,:,:,:,2), &   ! s = source
                            self%layers(:,:,:,:,:,3), &   ! q = input
                            self%vars  (:,:,:,:  ,3), &   ! q = solution
                            i3, 1, 3)
  end do
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
  !-----------------------------------------------------------------------------
  ! Integrate over bins
  !-----------------------------------------------------------------------------
  self%q = 0.0
  do i_bin=1,self%n_bin
    self%q = self%q + self%vars(:,:,:,i_bin,3)*self%w_bin(i_bin)
  end do
  deallocate (self%layers)
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Integral solver for one layer.  The exact integal solution Q=I-S for a 2nd
!> order polynomial has coefficient for the upstream Q, and for the 1st and 2nd
!> derivative of S.  However, if one is doing a solution both forwards and 
!> backward, the term in the 1st derivative cancels.  This is left as an option
!> here, in case the solutions forward and backward are made at different times.
!> =============================================================================
SUBROUTINE layer_solver (self, k, s, qin, qout, i3, l1, l3)
  class (rt_sc_t):: self
  real, dimension(:,:,:,:,:) :: k, s, qin
  real, dimension(:,:,:,:)   :: qout
  integer:: i3, l1, l3
  !.............................................................................
  integer:: l(2), u(2), il, i1, i2
  real:: dz2, dtau1, dtau2, dsdtau1, dsdtau2, dsdtau, d2sdtau2, ex0, ex1, ex2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%integral_forward', itimer=itimer)
  dz2 = self%mesh(3)%d/2.0
  l = self%mesh(1:2)%li
  u = self%mesh(1:2)%ui
  if (self%double_solution) then
    do il=1,self%nl
    do i2=l(2),u(2)
    do i1=l(1),u(2)
      dtau1    = dz2*(k(i1,i2,i3,il, 2) + k(i1,i2,i3,il,l1))
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
  else
    do il=1,self%nl
    do i2=l(2),u(2)
    do i1=l(1),u(2)
      dtau1    = dz2*(k(i1,i2,i3,il, 2) + k(i1,i2,i3,il,l1))
      dtau2    = dz2*(k(i1,i2,i3,il,l3) + k(i1,i2,i3,il, 2))
      dsdtau1  =     (s(i1,i2,i3,il, 2) - s(i1,i2,i3,il,l1))/dtau1
      dsdtau2  =     (s(i1,i2,i3,il,l3) - s(i1,i2,i3,il, 2))/dtau2
      dsdtau   = (dtau1*dsdtau2 + dtau2*dsdtau1)/(dtau1+dtau2)
      d2sdtau2 = (dsdtau2 - dsdtau1)*2.0/(dtau1+dtau2)
      ex0      = exp(-dtau1)
      ex1      = 1.0-ex0
      ex2      = ex1 - dtau1*ex0
      qout(i1,i2,i3,il) = qin(i1,i2,i3,il,l1)*ex0 &
                        + dsdtau             *ex1 &
                        + d2sdtau2           *ex2
    end do
    end do
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE layer_solver

!> ============================================================================
!> Integral SC/RT solver
!> ============================================================================
SUBROUTINE prepare_layers (self)
  class (rt_sc_t):: self
  integer:: il, iv
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%prepare_layers', itimer=itimer)
  do il=1,self%nl
    do iv=1,3
      call self%bilinear (self%vars(:,:,:,il,iv), self%layers(:,:,:,il,iv,1))
      self%layers(:,:,:,il,iv,2) = self%vars(:,:,:,il,iv)
      self%o = -self%o
      self%p = -self%p
      call self%bilinear (self%vars(:,:,:,il,iv), self%layers(:,:,:,il,iv,3))
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE prepare_layers

!> ============================================================================
!> Calculate offsets from one layer to the next
!> ============================================================================
SUBROUTINE calculate_offsets (self, i3)
  class (rt_sc_t):: self
  integer:: i3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%calculate_offsets', itimer=itimer)
  self%p = -self%mesh(3)%d*self%omega%mu*[cos(self%omega%phi),sin(self%omega%phi)]
  self%o = floor(self%p)
  self%p = self%p - self%o
  call trace%end (itimer)
END SUBROUTINE calculate_offsets

!> ============================================================================
!> Bilinear interpolation in a plane
!> ============================================================================
SUBROUTINE bilinear (self, in, out)
  class (rt_sc_t):: self
  real, dimension(:,:,:), intent(in)  :: in
  real, dimension(:,:,:), intent(out) :: out
  real:: p(2), q(2)
  integer:: i1, i2, i3, l(3), u(3), o(2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%bilinear', itimer=itimer)
  do i3=self%mesh(3)%li,self%mesh(3)%ui
    call self%calculate_offsets (i3)
    l = self%mesh%li
    u = self%mesh%ui
    o = self%o
    p = self%p
    q = 1.0-p
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      out(i1,i2,i3) = &
        p(1)*p(2)*in(i1+o(1)  ,i2+o(2)  ,i3) + &
        q(1)*p(2)*in(i1+o(1)+1,i2+o(2)  ,i3) + &
        p(1)*q(2)*in(i1+o(1)+1,i2+o(2)+1,i3) + &
        q(1)*q(2)*in(i1+o(1)+1,i2+o(2)+1,i3)
    end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE bilinear

!> ============================================================================
!> Compute the EOS quantities for the volume.  The 1,2 slots of the last index
!> correspond to opacity and source function, respectively.
!> ============================================================================
SUBROUTINE eos_volume (self)
  class (rt_sc_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_sc_t%eos_volume', itimer=itimer)
!  call eos%lookup (d, ee, self%vars(:,:,:,:,1), &
!                          self%vars(:,:,:,:,2))
  self%vars(:,:,:,:,1) = 1.0      ! fake opacity
  self%vars(:,:,:,:,2) = 1.0      ! fake source function
  call trace%end (itimer)
END SUBROUTINE eos_volume

END MODULE rt_sc_mod
