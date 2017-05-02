!> =============================================================================
!> $Id$
!> =============================================================================
MODULE stagger_mod
  USE omp_timer_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_16
  USE stagger_20
  USE stagger_24
  USE stagger_32
  USE stagger_36
  USE mesh_mod
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
    real:: ds(3)
    class(mesh_t), dimension(:), pointer:: m
  contains
    procedure:: ddxdn1
    procedure:: ddydn1
    procedure:: ddzdn1
    procedure:: ddxup1
    procedure:: ddyup1
    procedure:: ddzup1
    procedure, nopass:: xdn1
    procedure, nopass:: ydn1
    procedure, nopass:: zdn1
    procedure, nopass:: xup1
    procedure, nopass:: yup1
    procedure, nopass:: zup1
    procedure:: ddxdn
    procedure:: ddydn
    procedure:: ddzdn
    procedure:: ddxup
    procedure:: ddyup
    procedure:: ddzup
    procedure, nopass:: xdn
    procedure, nopass:: ydn
    procedure, nopass:: zdn
    procedure, nopass:: xup
    procedure, nopass:: yup
    procedure, nopass:: zup
    procedure:: div     => scalar_div_vector
    procedure:: curl_dn => vector_curl_dn_vector
    procedure:: curl_up => vector_curl_up_vector
    procedure:: grad    => vector_grad_scalar
    procedure:: ddown   => scalar_ddown_vector
  end type
  interface norm
    module procedure vector_norm_scalar
  end interface
  interface norm2
    module procedure vector_norm2_scalar
  end interface
  interface edge
    module procedure scalar_edge_vector
  end interface
  interface down
    module procedure scalar_down_vector
  end interface
  interface up
    module procedure vector_up_vector
  end interface
  interface dup
    module procedure vector_dup_vector
  end interface
CONTAINS

!> =============================================================================
FUNCTION ddxdn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddxdn (self, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddydn (self, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddzdn (self, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddxup (self, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddyup (self, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddzup (self, a)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xdn (a)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ydn (a)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zdn (a)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xup (a)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = yup (a)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zup (a)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: c
  real    :: ds
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(1)%d
    if (all(n==[16,16,16])) then
      call ddxdn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddxdn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddxdn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddxdn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddxdn_36 (ds,a,b)
    end if
  else if (n(1) > 1) then
    b(1,:,:) = 0.0
    c = 1./self%m(1)%d
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8) :: c
  real    :: ds
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddydn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(2)%d
    if (all(n==[16,16,16])) then
      call ddydn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddydn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddydn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddydn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddydn_36 (ds,a,b)
    end if
  else if (n(2) > 1) then
    b(:,1,:) = 0.0
    c = 1./self%m(2)%d
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: c
  real    :: ds
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(3)%d
    if (all(n==[16,16,16])) then
      call ddzdn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddzdn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddzdn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddzdn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddzdn_36 (ds,a,b)
    end if
  else if (n(3) > 1) then
    b(:,:,1) = 0.0
    c = 1./self%m(3)%d
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8) :: c
  real    :: ds
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(1)%d
    if (all(n==[16,16,16])) then
      call ddxup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddxup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddxup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddxup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddxup_36 (ds,a,b)
    end if
  else if (n(1) > 1) then
    b(n(1),:,:) = 0.0
    c = 1./self%m(1)%d
    !print*,'mk'
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=1,n(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8) :: c
  real    :: ds
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddyup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(2)%d
    if (all(n==[16,16,16])) then
      call ddyup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddyup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddyup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddyup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddyup_36 (ds,a,b)
    end if
  else if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,n(2),:) = 0.0
    c = 1./self%m(2)%d
    do iz=1,n(3)
    do iy=1,n(2)-1
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8) :: c
  real    :: ds
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    ds = self%m(3)%d
    if (all(n==[16,16,16])) then
      call ddzup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      call ddzup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      call ddzup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      call ddzup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      call ddzup_36 (ds,a,b)
    end if
  else if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,n(3)) = 0.0
    c = 1./self%m(3)%d
    do iz=1,n(3)-1
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call xdn_16 (a, b)
    else if (all(n==[20,20,20])) then
      call xdn_20 (a, b)
    else if (all(n==[24,24,24])) then
      call xdn_24 (a, b)
    else if (all(n==[32,32,32])) then
      call xdn_32 (a, b)
    else if (all(n==[36,36,36])) then
      call xdn_36 (a, b)
    end if
  else if (n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(1,:,:) = a(1,:,:)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ydn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call ydn_16 (a, b)
    else if (all(n==[20,20,20])) then
      call ydn_20 (a, b)
    else if (all(n==[24,24,24])) then
      call ydn_24 (a, b)
    else if (all(n==[32,32,32])) then
      call ydn_32 (a, b)
    else if (all(n==[36,36,36])) then
      call ydn_36 (a, b)
    end if
  else if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,1,:) = a(:,1,:)
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::zdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call zdn_16 (a, b)
    else if (all(n==[20,20,20])) then
      call zdn_20 (a, b)
    else if (all(n==[24,24,24])) then
      call zdn_24 (a, b)
    else if (all(n==[32,32,32])) then
      call zdn_32 (a, b)
    else if (all(n==[36,36,36])) then
      call zdn_36 (a, b)
    end if
  else if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,:,1) = a(:,:,1)
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call xup_16 (a, b)
    else if (all(n==[20,20,20])) then
      call xup_20 (a, b)
    else if (all(n==[24,24,24])) then
      call xup_24 (a, b)
    else if (all(n==[32,32,32])) then
      call xup_32 (a, b)
    else if (all(n==[36,36,36])) then
      call xup_36 (a, b)
    end if
  else if (n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(n(1),:,:) = a(n(1),:,:)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=1,n(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::yup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call yup_16 (a, b)
    else if (all(n==[20,20,20])) then
      call yup_20 (a, b)
    else if (all(n==[24,24,24])) then
      call yup_24 (a, b)
    else if (all(n==[32,32,32])) then
      call yup_32 (a, b)
    else if (all(n==[36,36,36])) then
      call yup_36 (a, b)
    end if
  else if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,n(2),:) = a(:,n(2),:)
    do iz=1,n(3)
    do iy=1,n(2)-1
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::zup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      call zup_16 (a, b)
    else if (all(n==[20,20,20])) then
      call zup_20 (a, b)
    else if (all(n==[24,24,24])) then
      call zup_24 (a, b)
    else if (all(n==[32,32,32])) then
      call zup_32 (a, b)
    else if (all(n==[36,36,36])) then
      call zup_36 (a, b)
    end if
  else if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,:,n(3)) = a(:,:,n(3))
    do iz=1,n(3)-1
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION

!> =============================================================================
FUNCTION vector_scalar_grad (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (self, in)
  out(:,:,:,2) = self%ddydn (in)
  out(:,:,:,3) = ddzdn (self, in)
END FUNCTION vector_scalar_grad

!> =============================================================================
FUNCTION vector_norm2_scalar (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)
END FUNCTION vector_norm2_scalar
!> =============================================================================
FUNCTION vector_norm_scalar (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = sqrt(in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)**2)
END FUNCTION vector_norm_scalar

!> =============================================================================
FUNCTION scalar_edge_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ydn(zdn(in))
  out(:,:,:,2) = zdn(xdn(in))
  out(:,:,:,3) = xdn(ydn(in))
END FUNCTION scalar_edge_vector
!> =============================================================================
FUNCTION scalar_down_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xdn (in)
  out(:,:,:,2) = ydn (in)
  out(:,:,:,3) = zdn (in)
END FUNCTION scalar_down_vector
!> =============================================================================
FUNCTION vector_up_vector (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xup (in(:,:,:,1))
  out(:,:,:,2) = yup (in(:,:,:,2))
  out(:,:,:,3) = zup (in(:,:,:,3))
END FUNCTION vector_up_vector

!> =============================================================================
FUNCTION scalar_ddown_vector (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = self%ddxdn (in)
  out(:,:,:,2) = self%ddydn (in)
  out(:,:,:,3) = self%ddzdn (in)
END FUNCTION scalar_ddown_vector
!> =============================================================================
FUNCTION vector_dup_vector (self, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(stagger_t):: self
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxup (self, in(:,:,:,1))
  out(:,:,:,2) = ddyup (self, in(:,:,:,2))
  out(:,:,:,3) = ddzup (self, in(:,:,:,3))
END FUNCTION vector_dup_vector

!> =============================================================================
FUNCTION vector_grad_scalar (self, in) RESULT(out)
  class(stagger_t):: self
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (self, in)
  out(:,:,:,2) = ddydn (self, in)
  out(:,:,:,3) = ddzdn (self, in)
END FUNCTION
!> =============================================================================
FUNCTION scalar_div_vector (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = ddxup(self, in(:,:,:,1)) + ddyup(self, in(:,:,:,2)) + ddzup(self, in(:,:,:,3))
END FUNCTION
!> =============================================================================
FUNCTION scalar_minusdiv_vector (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = - ddxup(self, in(:,:,:,1)) - ddyup(self, in(:,:,:,2)) - ddzup(self, in(:,:,:,3))
END FUNCTION
!> =============================================================================
FUNCTION vector_curl_dn_vector (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddydn(self, in(:,:,:,3)) - ddzdn(self, in(:,:,:,2))
  out(:,:,:,2) = ddzdn(self, in(:,:,:,1)) - ddxdn(self, in(:,:,:,3))
  out(:,:,:,3) = ddxdn(self, in(:,:,:,2)) - ddydn(self, in(:,:,:,1))
END FUNCTION
!> =============================================================================
FUNCTION vector_curl_up_vector (self, in) RESULT (out)
  class(stagger_t):: self
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddyup(self, in(:,:,:,3)) - ddzup(self, in(:,:,:,2))
  out(:,:,:,2) = ddzup(self, in(:,:,:,1)) - ddxup(self, in(:,:,:,3))
  out(:,:,:,3) = ddxup(self, in(:,:,:,2)) - ddyup(self, in(:,:,:,1))
END FUNCTION
!*******************************************************************************
END MODULE
