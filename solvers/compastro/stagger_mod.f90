!===============================================================================
!> $Id$
!===============================================================================
MODULE operators_mod
  USE io_mod
  implicit none
  private
  type, public:: operator_t
    integer:: n(3), l(3), u(3)
    real:: ds(3)=1.0
  contains
    procedure:: init
    procedure ddxdn1
    procedure ddydn1
    procedure ddzdn1
    procedure ddxup1
    procedure ddyup1
    procedure ddzup1
    procedure xdn1
    procedure ydn1
    procedure zdn1
    procedure xup1
    procedure yup1
    procedure zup1
    procedure ddxdn
    procedure ddydn
    procedure ddzdn
    procedure ddxup
    procedure ddyup
    procedure ddzup
    procedure xdn
    procedure ydn
    procedure zdn
    procedure xup
    procedure yup
    procedure zup
  end type
CONTAINS
!===============================================================================
SUBROUTINE init (self, n, ds)
  class(operator_t):: self
  integer:: n(3)
  real:: ds(3)
  !.............................................................................
  self%ds = ds
  self%n = n
  self%l = 1
  self%u = n
END SUBROUTINE init
!===============================================================================

!*******************************************************************************
FUNCTION ddxdn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddxdn (self, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddydn (self, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddzdn (self, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddxup (self, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddyup (self, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ddzup (self, a)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = xdn (self, a)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = ydn (self, a)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = zdn (self, a)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = xup (self, a)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = yup (self, a)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  b = zup (self, a)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  b(self%l(1),:,:) = 0.0
END FUNCTION ddxdn

!*******************************************************************************
SUBROUTINE ddxdn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                c*(a(ix,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddx_add:in')
  !call scalar_stats (b, 'ddx_add:out')
END SUBROUTINE ddxdn_add

!*******************************************************************************
FUNCTION ddydn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    b(:,self%l(2),:) = 0.0
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  b(:,self%l(2),:) = 0.0
END FUNCTION ddydn

!*******************************************************************************
SUBROUTINE ddydn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddydn_add:in')
  !call scalar_stats (b, 'ddydn_add:out')
END SUBROUTINE ddydn_add

!*******************************************************************************
FUNCTION ddzdn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    b(:,:,self%l(3)) = 0.0
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
    b(:,:,self%l(3)) = 0.0
  else
    b = 0.0
  endif
END FUNCTION ddzdn

!*******************************************************************************
SUBROUTINE ddzdn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddzdn_add:in')
  !call scalar_stats (b, 'ddzdn_add:out')
END SUBROUTINE ddzdn_add

!*******************************************************************************
FUNCTION ddxup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
    b(self%u(1),:,:) = 0.0
  else
    b = 0.0
  endif
END FUNCTION ddxup

!*******************************************************************************
SUBROUTINE ddxup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./self%ds(1)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddx_add:in')
  !call scalar_stats (b, 'ddx_add:out')
END SUBROUTINE ddxup_add

!*******************************************************************************
FUNCTION ddyup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    b(:,self%u(2),:) = 0.0
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b(:,self%u(2),:) = 0.0
END FUNCTION ddyup

!*******************************************************************************
SUBROUTINE ddyup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(2) > 1) then
    c = 1./self%ds(2)
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddyup_add:in')
  !call scalar_stats (b, 'ddyup_add:out')
END SUBROUTINE ddyup_add

!*******************************************************************************
FUNCTION ddzup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
    b(:,:,self%u(3)) = 0.0
  else
    b = 0.0
  endif
END FUNCTION ddzup

!*******************************************************************************
SUBROUTINE ddzup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  if (self%n(3) > 1) then
    c = 1./self%ds(3)
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ddzup_add:in')
  !call scalar_stats (b, 'ddzup_add:out')
END SUBROUTINE ddzup_add

!*******************************************************************************
FUNCTION xdn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2),self%u(2)
  do ix=self%l(1)+1,self%u(1)
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
  b(self%l(1),:,:) = 0.0
END FUNCTION xdn

!*******************************************************************************
SUBROUTINE xdn_subr (self, a, b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2),self%u(2)
  do ix=self%l(1)+1,self%u(1)
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
  b(self%l(1),:,:) = 0.0
END SUBROUTINE

!*******************************************************************************
SUBROUTINE xdn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  !
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1)+1,self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'x_a:in')
  !call scalar_stats (b, 'x_a:out')
END SUBROUTINE xdn_add

!*******************************************************************************
FUNCTION ydn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  c = 0.5
  do iz=self%l(3),self%u(3)
  do iy=self%l(2)+1,self%u(2)
  do ix=self%l(1),self%u(1)
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
  b(:,self%l(2),:) = 0.0
END FUNCTION ydn

!*******************************************************************************
SUBROUTINE ydn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2)+1,self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy  ,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'ydn_a:in')
  !call scalar_stats (b, 'ydn_a:out')
END SUBROUTINE ydn_add

!*******************************************************************************
FUNCTION zdn (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  real:: c
  integer:: ix, iy, iz
!...............................................................................
  if (self%n(3) > 1) then
    c = 0.5
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
    b(:,:,self%l(3)) = 0.0
  else
    b = 0.0
  end if
END FUNCTION zdn

!*******************************************************************************
SUBROUTINE zdn_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3)+1,self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'zdn_a:in')
  !call scalar_stats (b, 'zdn_a:out')
END SUBROUTINE zdn_add

!*******************************************************************************
FUNCTION xup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
  b(self%u(1),:,:) = 0.0
  endif
END FUNCTION

!*******************************************************************************
SUBROUTINE xup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)-1
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'x_a:in')
  !call scalar_stats (b, 'x_a:out')
END SUBROUTINE xup_add

!*******************************************************************************
FUNCTION yup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b(:,self%u(2),:) = 0.0
END FUNCTION

!*******************************************************************************
SUBROUTINE yup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, by, sign
!...............................................................................
  c = 0.5
  if (self%n(2) > 1) then
    do iz=self%l(3),self%u(3)
    do iy=self%l(2),self%u(2)-1
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'yup_a:in')
  !call scalar_stats (b, 'yup_a:out')
END SUBROUTINE yup_add

!*******************************************************************************
FUNCTION zup (self, a) RESULT (b)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
    end do
    end do
    end do
    b(:,:,self%u(3)) = 0.0
  else
    b = 0.0
  endif
END FUNCTION

!*******************************************************************************
SUBROUTINE zup_add (self, a, b, sign)
  class(operator_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: a, b
  integer :: ix, iy, iz
  real    :: c, sign
!...............................................................................
  c = 0.5
  if (self%n(3) > 1) then
    do iz=self%l(3),self%u(3)-1
    do iy=self%l(2),self%u(2)
    do ix=self%l(1),self%u(1)
      b(ix,iy,iz) = b(ix,iy,iz) + sign* &
                   c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
    end do
    end do
    end do
  endif
  !call scalar_stats (self, a, 'zup_a:in')
  !call scalar_stats (b, 'zup_a:out')
END SUBROUTINE zup_add

!*******************************************************************************
END MODULE
