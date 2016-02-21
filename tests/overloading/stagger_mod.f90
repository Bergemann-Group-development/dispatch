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
  interface grad
    module procedure vector_real_grad, vector_int_grad, vector_scalar_grad
  end interface
  interface div
    module procedure vector_scalar_div
  end interface
  interface minusdiv
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
FUNCTION ddxdn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddxdn (self)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddydn (self)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddzdn (self)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddxup (self)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddyup (self)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ddzup (self)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = xdn (self)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = ydn (self)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = zdn (self)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = xup (self)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = yup (self)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  b = zup (self)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddxdn')
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix  ,iy,iz) - self%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  b%f(self%l(1),:,:) = 0.0
  call self%tempdealloc('ddxdn')
END FUNCTION ddxdn

!*******************************************************************************
SUBROUTINE ddxdn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                c*(self%f(ix,iy,iz) - self%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddxdn_add

!*******************************************************************************
FUNCTION ddydn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddydn')
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    b%f(:,self%l(2),:) = 0.0
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy  ,iz) - self%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  b%f(:,self%l(2),:) = 0.0
  call self%tempdealloc('ddydn')
END FUNCTION ddydn

!*******************************************************************************
SUBROUTINE ddydn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy  ,iz) - self%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddydn_add

!*******************************************************************************
FUNCTION ddzdn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddydn')
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    b%f(:,:,self%l(3)) = 0.0
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy,iz  ) - self%f(ix,iy,iz-1))
    end do
    end do
    end do
    b%f(:,:,self%l(3)) = 0.0
  else
    b%f = 0.0
  endif
  call self%tempdealloc('ddzdn')
END FUNCTION ddzdn

!*******************************************************************************
SUBROUTINE ddzdn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy,iz  ) - self%f(ix,iy,iz-1))
    end do
    end do
    end do
  endif
END SUBROUTINE ddzdn_add

!*******************************************************************************
FUNCTION ddxup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddxup')
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b%f(ix,iy,iz) = c*(self%f(ix+1,iy,iz) - self%f(ix  ,iy,iz))
    end do
    end do
    end do
    b%f(self%u(1),:,:) = 0.0
  else
    b%f = 0.0
  endif
  call self%tempdealloc('ddxup')
END FUNCTION ddxup

!*******************************************************************************
SUBROUTINE ddxup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix+1,iy,iz) - self%f(ix  ,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddxup_add

!*******************************************************************************
FUNCTION ddyup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddyup')
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    b%f(:,self%u(2),:) = 0.0
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy+1,iz) - self%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b%f(:,self%u(2),:) = 0.0
  call self%tempdealloc('ddyup')
END FUNCTION ddyup

!*******************************************************************************
SUBROUTINE ddyup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy+1,iz) - self%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ddyup_add

!*******************************************************************************
FUNCTION ddzup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  call b%tempalloc (self,'ddzup')
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy,iz+1) - self%f(ix,iy,iz  ))
    end do
    end do
    end do
    b%f(:,:,self%u(3)) = 0.0
  else
    b%f = 0.0
  endif
  call self%tempdealloc('ddzup')
END FUNCTION ddzup

!*******************************************************************************
SUBROUTINE ddzup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy,iz+1) - self%f(ix,iy,iz  ))
    end do
    end do
    end do
  endif
END SUBROUTINE ddzup_add

!*******************************************************************************
FUNCTION xdn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2),self%u(2)
  do ix=self%l(1)+1,self%u(1)
    b%f(ix,iy,iz) = c*(self%f(ix  ,iy,iz) + self%f(ix-1,iy,iz))
  end do
  end do
  end do
  b%f(self%l(1),:,:) = 0.0
END FUNCTION xdn

!*******************************************************************************
SUBROUTINE xdn_subr (self, b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2),self%u(2)
  do ix=self%l(1)+1,self%u(1)
    b%f(ix,iy,iz) = c*(self%f(ix  ,iy,iz) + self%f(ix-1,iy,iz))
  end do
  end do
  end do
  b%f(self%l(1),:,:) = 0.0
END SUBROUTINE

!*******************************************************************************
SUBROUTINE xdn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  !
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix  ,iy,iz) + self%f(ix-1,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE xdn_add

!*******************************************************************************
FUNCTION ydn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2)+1,self%u(2)
  do ix=self%l(1),self%u(1)
    b%f(ix,iy,iz) = c*(self%f(ix,iy,iz) + self%f(ix,iy-1,iz))
  end do
  end do
  end do
  b%f(:,self%l(2),:) = 0.0
END FUNCTION ydn

!*******************************************************************************
SUBROUTINE ydn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy  ,iz) + self%f(ix,iy-1,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE ydn_add

!*******************************************************************************
FUNCTION zdn (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  if (self%n(3) > 1) then
    c = 0.5
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy,iz) + self%f(ix,iy,iz-1))
    end do
    end do
    end do
    b%f(:,:,self%l(3)) = 0.0
  else
    b%f = 0.0
  end if
END FUNCTION zdn

!*******************************************************************************
SUBROUTINE zdn_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                c*(self%f(ix,iy,iz) + self%f(ix,iy,iz-1))
    end do
    end do
    end do
  endif
END SUBROUTINE zdn_add

!*******************************************************************************
FUNCTION xup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b%f(ix,iy,iz) = c*(self%f(ix+1,iy,iz) + self%f(ix  ,iy,iz))
    end do
    end do
    end do
  b%f(self%u(1),:,:) = 0.0
  endif
END FUNCTION

!*******************************************************************************
SUBROUTINE xup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix+1,iy,iz) + self%f(ix  ,iy,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE xup_add

!*******************************************************************************
FUNCTION yup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy+1,iz) + self%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b%f(:,self%u(2),:) = 0.0
END FUNCTION

!*******************************************************************************
SUBROUTINE yup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, by, sign
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy+1,iz) + self%f(ix,iy  ,iz))
    end do
    end do
    end do
  endif
END SUBROUTINE yup_add

!*******************************************************************************
FUNCTION zup (self) RESULT (b)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = c*(self%f(ix,iy,iz+1) + self%f(ix,iy,iz  ))
    end do
    end do
    end do
    b%f(:,:,self%u(3)) = 0.0
  else
    b%f = 0.0
  endif
END FUNCTION

!*******************************************************************************
SUBROUTINE zup_add (self, b, sign)
  type(scalar_t), intent(in):: self
  type(scalar_t):: b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b%f(ix,iy,iz) = b%f(ix,iy,iz) + sign* &
                   c*(self%f(ix,iy,iz+1) + self%f(ix,iy,iz  ))
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
