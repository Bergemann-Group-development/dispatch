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
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
    real(8):: ds(3)
  contains
    procedure, nopass:: ddxdn1
    procedure, nopass:: ddydn1
    procedure, nopass:: ddzdn1
    procedure, nopass:: ddxup1
    procedure, nopass:: ddyup1
    procedure, nopass:: ddzup1
    procedure, nopass:: xdn1
    procedure, nopass:: ydn1
    procedure, nopass:: zdn1
    procedure, nopass:: xup1
    procedure, nopass:: yup1
    procedure, nopass:: zup1
    procedure, nopass:: ddxdn
    procedure, nopass:: ddydn
    procedure, nopass:: ddzdn
    procedure, nopass:: ddxup
    procedure, nopass:: ddyup
    procedure, nopass:: ddzup
    procedure, nopass:: xdn
    procedure, nopass:: ydn
    procedure, nopass:: zdn
    procedure, nopass:: xup
    procedure, nopass:: yup
    procedure, nopass:: zup
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
  interface ddown
    module procedure scalar_ddown_vector
  end interface
  interface dup
    module procedure vector_dup_vector
  end interface
  interface grad
    module procedure vector_grad_scalar
  end interface
  interface div
    module procedure scalar_div_vector
  end interface
  interface minusdiv
    module procedure scalar_minusdiv_vector
  end interface
  interface curl_dn
    module procedure vector_curl_dn_vector
  end interface
  interface curl_up
    module procedure vector_curl_up_vector
  end interface
CONTAINS

!> =============================================================================
FUNCTION ddxdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddxdn (ds, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddydn (ds, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddzdn (ds, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddxup (ds, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddyup (ds, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3)
    b= ddzup (ds, a)
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
FUNCTION ddxdn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dxdn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dxdn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dxdn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dxdn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dxdn_36 (ds,a,b)
    end if
  else if (n(1) > 1) then
    b(1,:,:) = 0.0
    c = 1./ds(1)
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
FUNCTION ddydn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddydn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dydn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dydn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dydn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dydn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dydn_36 (ds,a,b)
    end if
  else if (n(2) > 1) then
    b(:,1,:) = 0.0
    c = 1./ds(2)
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
FUNCTION ddzdn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzdn', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dzdn_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dzdn_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dzdn_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dzdn_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dzdn_36 (ds,a,b)
    end if
  else if (n(3) > 1) then
    b(:,:,1) = 0.0
    c = 1./ds(3)
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
FUNCTION ddxup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dxup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dxup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dxup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dxup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dxup_36 (ds,a,b)
    end if
  else if (n(1) > 1) then
    b(n(1),:,:) = 0.0
    c = 1./ds(1)
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
FUNCTION ddyup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddyup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dyup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dyup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dyup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dyup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dyup_36 (ds,a,b)
    end if
  else if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,n(2),:) = 0.0
    c = 1./ds(2)
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
FUNCTION ddzup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzup', itimer=itimer)
  n = shape(a)
  if (hardwire) then
    if (all(n==[16,16,16])) then
      !call dzup_16 (ds,a,b)
    else if (all(n==[20,20,20])) then
      !call dzup_20 (ds,a,b)
    else if (all(n==[24,24,24])) then
      !call dzup_24 (ds,a,b)
    else if (all(n==[32,32,32])) then
      !call dzup_32 (ds,a,b)
    else if (all(n==[36,36,36])) then
      !call dzup_36 (ds,a,b)
    end if
  else if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,n(3)) = 0.0
    c = 1./ds(3)
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
FUNCTION vector_scalar_grad (ds, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (ds, in)
  out(:,:,:,2) = ddydn (ds, in)
  out(:,:,:,3) = ddzdn (ds, in)
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
FUNCTION scalar_ddown_vector (ds, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (ds, in)
  out(:,:,:,2) = ddydn (ds, in)
  out(:,:,:,3) = ddzdn (ds, in)
END FUNCTION scalar_ddown_vector
!> =============================================================================
FUNCTION vector_dup_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxup (ds, in(:,:,:,1))
  out(:,:,:,2) = ddyup (ds, in(:,:,:,2))
  out(:,:,:,3) = ddzup (ds, in(:,:,:,3))
END FUNCTION vector_dup_vector

!> =============================================================================
FUNCTION vector_grad_scalar (ds, in) RESULT(out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (ds, in)
  out(:,:,:,2) = ddydn (ds, in)
  out(:,:,:,3) = ddzdn (ds, in)
END FUNCTION
!> =============================================================================
FUNCTION scalar_div_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out = ddxup(ds, in(:,:,:,1)) + ddyup(ds, in(:,:,:,2)) + ddzup(ds, in(:,:,:,3))
END FUNCTION
!> =============================================================================
FUNCTION scalar_minusdiv_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out = - ddxup(ds, in(:,:,:,1)) - ddyup(ds, in(:,:,:,2)) - ddzup(ds, in(:,:,:,3))
END FUNCTION
!> =============================================================================
FUNCTION vector_curl_dn_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddydn(ds, in(:,:,:,3)) - ddzdn(ds, in(:,:,:,2))
  out(:,:,:,2) = ddzdn(ds, in(:,:,:,1)) - ddxdn(ds, in(:,:,:,3))
  out(:,:,:,3) = ddxdn(ds, in(:,:,:,2)) - ddydn(ds, in(:,:,:,1))
END FUNCTION
!> =============================================================================
FUNCTION vector_curl_up_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddyup(ds, in(:,:,:,3)) - ddzup(ds, in(:,:,:,2))
  out(:,:,:,2) = ddzup(ds, in(:,:,:,1)) - ddxup(ds, in(:,:,:,3))
  out(:,:,:,3) = ddxup(ds, in(:,:,:,2)) - ddyup(ds, in(:,:,:,1))
END FUNCTION
!*******************************************************************************
END MODULE
