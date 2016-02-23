!> =============================================================================
!> $Id$
!> =============================================================================
MODULE stagger_mod
  USE scalar_mod
  USE vector_mod
  implicit none
  public
  integer, private:: verbose=0
  type, public:: stagger_t
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
    module procedure scalar_up_vector
  end interface
  interface ddown
    module procedure scalar_ddown_vector
  end interface
  interface dup
    module procedure vector_dup_vector
  end interface
  interface grad
    module procedure vector_real_grad, vector_int_grad, vector_scalar_grad
  end interface
  interface div
    module procedure vector_scalar_div
  end interface
  interface minus_div
    module procedure vector_scalar_minusdiv
  end interface
  interface cross
    module procedure real4_vector_cross, real8_vector_cross, vector_cross, real4_array_cross
  end interface
  interface curl
    module procedure vector_curl
  end interface
CONTAINS

!> =============================================================================
FUNCTION ddxdn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddxdn (a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddydn (a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddzdn (a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddxup (a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddyup (a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ddzup (a)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = xdn (a)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = ydn (a)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = zdn (a)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = xup (a)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = yup (a)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  b = zup (a)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'ddxdn')
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./a%ds(1)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1)+1,a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix  ,iy,iz) - a%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  b%f(a%l(1),:,:) = 0.0
  call a%tempdealloc('ddxdn')
END FUNCTION ddxdn

!*******************************************************************************
SUBROUTINE ddxdn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./a%ds(1)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1)+1,a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                c*(a%f(ix,iy,iz) - a%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddxdn_add

!*******************************************************************************
FUNCTION ddydn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'ddydn')
  if (a%n(2) > 1) then
    c = 1./a%ds(2)
    b%f(:,a%l(2),:) = 0.0
    do iz=a%l(3),a%u(3)
    do iy=a%l(2)+1,a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy  ,iz) - a%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  b%f(:,a%l(2),:) = 0.0
  call a%tempdealloc('ddydn')
END FUNCTION ddydn

!*******************************************************************************
SUBROUTINE ddydn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(2) > 1) then
    c = 1./a%ds(2)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2)+1,a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy  ,iz) - a%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddydn_add

!*******************************************************************************
FUNCTION ddzdn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  print*,'b',b%id,b%temporary,b%allocated,associated(b%f)
  call b%tempalloc (a,'ddzdn')
  print*,'b',b%id,b%temporary,b%allocated,associated(b%f)
  print*,'l',a%l,a%id,a%n
  print*,'u',a%u
  if (a%n(3) > 1) then
    c = 1./a%ds(3)
    b%f(:,:,a%l(3)) = 0.0
    do iz=a%l(3)+1,a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy,iz  ) - a%f(ix,iy,iz-1))
    end do
    end do
    end do
    b%f(:,:,a%l(3)) = 0.0
  else
    b%f = 0.0
  endif
  call a%tempdealloc('ddzdn')
END FUNCTION ddzdn

!*******************************************************************************
SUBROUTINE ddzdn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(3) > 1) then
    c = 1./a%ds(3)
    do iz=a%l(3)+1,a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy,iz  ) - a%f(ix,iy,iz-1))
    end do
    end do
    end do
  endif
END SUBROUTINE ddzdn_add

!*******************************************************************************
FUNCTION ddxup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'ddxup')
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./a%ds(1)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)-1
      b%f(ix,iy,iz) = c*(a%f(ix+1,iy,iz) - a%f(ix  ,iy,iz))
    end do
    end do
    end do
    b%f(a%u(1),:,:) = 0.0
  else
    b%f = 0.0
  endif
  call a%tempdealloc('ddxup')
END FUNCTION ddxup

!*******************************************************************************
SUBROUTINE ddxup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./a%ds(1)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)-1
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix+1,iy,iz) - a%f(ix  ,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddxup_add

!*******************************************************************************
FUNCTION ddyup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'ddyup')
  if (a%n(2) > 1) then
    c = 1./a%ds(2)
    b%f(:,a%u(2),:) = 0.0
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)-1
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy+1,iz) - a%f(ix,iy  ,iz))
    end do
    end do
    end do
    b%f(:,a%u(2),:) = 0.0
  else
    b%f = 0.0
  endif
  call a%tempdealloc('ddyup')
END FUNCTION ddyup

!*******************************************************************************
SUBROUTINE ddyup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(2) > 1) then
    c = 1./a%ds(2)
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)-1
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy+1,iz) - a%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddyup_add

!*******************************************************************************
FUNCTION ddzup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'ddzup')
  print*,'l',a%l,a%id,a%n
  print*,'u',a%u
  if (a%n(3) > 1) then
    c = 1./a%ds(3)
    do iz=a%l(3),a%u(3)-1
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy,iz+1) - a%f(ix,iy,iz  ))
    end do
    end do
    end do
    b%f(:,:,a%u(3)) = 0.0
  else
    b%f = 0.0
  endif
  call a%tempdealloc('ddzup')
END FUNCTION ddzup

!*******************************************************************************
SUBROUTINE ddzup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (a%n(3) > 1) then
    c = 1./a%ds(3)
    do iz=a%l(3),a%u(3)-1
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy,iz+1) - a%f(ix,iy,iz  ))
    end do
    end do
    end do
  endif
END SUBROUTINE ddzup_add

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  call b%tempalloc (a,'xdn')
  if (a%n(1) > 1) then
    c = 0.5
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1)+1,a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix  ,iy,iz) + a%f(ix-1,iy,iz))
    end do
    end do
    end do
    b%f(a%l(1),:,:) = 0.0
  else
    b%f=a%f
  end if
  call a%tempdealloc('xdn')
END FUNCTION xdn

!*******************************************************************************
SUBROUTINE xdn_subr (a, b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=a%l(3),a%u(3)
  do iy=a%l(2),a%u(2)
  do ix=a%l(1)+1,a%u(1)
    b%f(ix,iy,iz) = c*(a%f(ix  ,iy,iz) + a%f(ix-1,iy,iz))
  end do
  end do
  end do
  b%f(a%l(1),:,:) = 0.0
END SUBROUTINE

!*******************************************************************************
SUBROUTINE xdn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  !
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1)+1,a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix  ,iy,iz) + a%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE xdn_add

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  call b%tempalloc (a,'ydn')
  if (a%n(2) > 1) then
    c = 0.5
    do iz=a%l(3),a%u(3)
    do iy=a%l(2)+1,a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy,iz) + a%f(ix,iy-1,iz))
    end do
    end do
    end do
    b%f(:,a%l(2),:) = 0.0
  else
    b%f = a%f
  end if
  call a%tempdealloc('ydn')
END FUNCTION ydn

!*******************************************************************************
SUBROUTINE ydn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (a%n(2) > 1) then
    do iz=a%l(3),a%u(3)
    do iy=a%l(2)+1,a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy  ,iz) + a%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ydn_add

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  call b%tempalloc (a,'zdn')
  if (a%n(3) > 1) then
    c = 0.5
    do iz=a%l(3)+1,a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy,iz) + a%f(ix,iy,iz-1))
    end do
    end do
    end do
    b%f(:,:,a%l(3)) = 0.0
  else
    b%f = 0.0
  end if
  call a%tempdealloc('zdn')
END FUNCTION zdn

!*******************************************************************************
SUBROUTINE zdn_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (a%n(3) > 1) then
    do iz=a%l(3)+1,a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                c*(a%f(ix,iy,iz) + a%f(ix,iy,iz-1))
    end do
    end do
    end do
  endif
END SUBROUTINE zdn_add

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'xup')
  c = 0.5
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)-1
      b%f(ix,iy,iz) = c*(a%f(ix+1,iy,iz) + a%f(ix  ,iy,iz))
    end do
    end do
    end do
    b%f(a%u(1),:,:) = 0.0
  else
    b%f=a%f
  endif
  call a%tempdealloc('xup')
END FUNCTION

!*******************************************************************************
SUBROUTINE xup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (a%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)-1
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix+1,iy,iz) + a%f(ix  ,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE xup_add

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'yup')
  c = 0.5
  if (a%n(2) > 1) then
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)-1
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy+1,iz) + a%f(ix,iy  ,iz))
    end do
    end do
    end do
    b%f(:,a%u(2),:) = 0.0
  else
    b%f = a%f
  endif
  call a%tempdealloc('yup')
END FUNCTION

!*******************************************************************************
SUBROUTINE yup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, by, sign
!...............................................................................
  c = 0.5
  if (a%n(2) > 1) then
    do iz=a%l(3),a%u(3)
    do iy=a%l(2),a%u(2)-1
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy+1,iz) + a%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE yup_add

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (a,'zup')
  c = 0.5
  if (a%n(3) > 1) then
    do iz=a%l(3),a%u(3)-1
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = c*(a%f(ix,iy,iz+1) + a%f(ix,iy,iz  ))
    end do
    end do
    end do
    b%f(:,:,a%u(3)) = 0.0
  else
    b%f = 0.0
  endif
  call a%tempdealloc('zup')
END FUNCTION

!*******************************************************************************
SUBROUTINE zup_add (a, b, sign)
  type(scalar_t), intent(in):: a
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (a%n(3) > 1) then
    do iz=a%l(3),a%u(3)-1
    do iy=a%l(2),a%u(2)
    do ix=a%l(1),a%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(a%f(ix,iy,iz+1) + a%f(ix,iy,iz  ))
    end do
    end do
    end do
  endif
END SUBROUTINE zup_add

!> =============================================================================
FUNCTION vector_curl (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  if (verbose>0) print*,'curl_x'
  out%x = ddyup(in%z) - ddzup(in%y)
  if (verbose>0) print*,'curl_y'
  out%y = ddzup(in%x) - ddxup(in%z)
  if (verbose>0) print*,'curl_z'
  out%z = ddxup(in%y) - ddyup(in%x)
  if (verbose>0) print*,'curl_end'
END FUNCTION vector_curl
!> =============================================================================
FUNCTION vector_real_grad (in) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
END FUNCTION vector_real_grad
!> =============================================================================
FUNCTION vector_int_grad (in) RESULT (out)
  type(vector_t):: out
  integer, intent(in):: in
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
END FUNCTION vector_int_grad
!> =============================================================================
FUNCTION vector_scalar_grad (in) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(inout):: in
  !-----------------------------------------------------------------------------
  out%x = ddxdn (in)
  out%y = ddydn (in)
  out%z = ddzdn (in)
  call in%tempdealloc
END FUNCTION vector_scalar_grad
!> =============================================================================
FUNCTION vector_norm2_scalar (in) RESULT (out)
  type(vector_t), intent(inout):: in
  type(scalar_t):: out
  !-----------------------------------------------------------------------------
  call trace_begin('norm')
  call out%tempalloc (in%x)
  out = xup(in%x)**2 + yup(in%y)**2 + zup(in%z)**2
  call in%tempdealloc
  call trace_end
END FUNCTION vector_norm2_scalar
!> =============================================================================
FUNCTION vector_norm_scalar (in) RESULT (out)
  type(vector_t), intent(inout):: in
  type(scalar_t):: out
  !-----------------------------------------------------------------------------
  call trace_begin('norm')
  call out%tempalloc (in%x)
  out = sqrt(norm2(in))
  call in%tempdealloc
  call trace_end
END FUNCTION vector_norm_scalar
!> =============================================================================
FUNCTION scalar_edge_vector (in) RESULT (out)
  type(scalar_t), intent(inout):: in
  type(vector_t):: out
  !-----------------------------------------------------------------------------
  out%x = ydn(zdn(in))
  out%y = zdn(xdn(in))
  out%z = xdn(ydn(in))
  call in%tempdealloc
END FUNCTION scalar_edge_vector
!> =============================================================================
FUNCTION scalar_down_vector (in) RESULT (out)
  type(scalar_t), intent(inout):: in
  type(vector_t):: out
  !-----------------------------------------------------------------------------
  if (verbose>1) print*, 'down  (in):', in%id, in%n, in%temporary, in%allocated
  out%x = xdn (in)
  out%y = ydn (in)
  out%z = zdn (in)
  if (verbose>1) print*, 'down (out):', out%id, out%n, out%temporary, out%allocated
  call in%tempdealloc
END FUNCTION scalar_down_vector
!> =============================================================================
FUNCTION scalar_up_vector (in) RESULT (out)
  type(scalar_t), intent(inout):: in
  type(vector_t):: out
  !-----------------------------------------------------------------------------
  out%x = xup (in)
  out%y = yup (in)
  out%z = zup (in)
  call in%tempdealloc
END FUNCTION scalar_up_vector
!> =============================================================================
FUNCTION scalar_ddown_vector (in) RESULT (out)
  type(scalar_t), intent(inout):: in
  type(vector_t):: out
  !-----------------------------------------------------------------------------
  out%x = ddxdn (in)
  out%y = ddydn (in)
  out%z = ddzdn (in)
  call in%tempdealloc
END FUNCTION scalar_ddown_vector
!> =============================================================================
FUNCTION vector_dup_vector (in) RESULT (out)
  type(vector_t), intent(inout):: in
  type(vector_t):: out
  !-----------------------------------------------------------------------------
  out%x = ddxup (in%x)
  out%y = ddyup (in%y)
  out%z = ddzup (in%z)
  call in%tempdealloc
END FUNCTION vector_dup_vector
!> =============================================================================
FUNCTION grads (in) RESULT(out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in
  !-----------------------------------------------------------------------------
  out%x = ddxdn (in)
  out%y = ddydn (in)
  out%z = ddzdn (in)
END FUNCTION grads
!> =============================================================================
FUNCTION vector_scalar_div (in) RESULT (out)
  type(vector_t), intent(in):: in
  type(scalar_t):: out
  !-----------------------------------------------------------------------------
  if (verbose>0) print*,in%x%l,in%x%u
  out = ddxup(in%x) + ddyup(in%y) + ddzup(in%z)
  !call in%tempdealloc
END FUNCTION vector_scalar_div
!> =============================================================================
FUNCTION vector_scalar_minusdiv (in) RESULT(out)
  type(vector_t), intent(in):: in
  type(scalar_t):: out
  !-----------------------------------------------------------------------------
  out = -ddxup (in%x) - ddyup (in%y) - ddzup (in%z)
END FUNCTION vector_scalar_minusdiv

!*******************************************************************************
END MODULE
