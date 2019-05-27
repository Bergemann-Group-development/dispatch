!===============================================================================
!> $Id$
!===============================================================================
MODULE vector_ops
  USE stagger_mod
  USE mesh_mod
  implicit none
  public
  interface norm
    module procedure vector_norm_scalar
  end interface
  interface norm2
    module procedure vector_norm2_scalar
  end interface
  interface edge
    module procedure scalar_edge_vector
  end interface
  interface edge_up
    module procedure scalar_edge_vector_up
  end interface
  interface face
    module procedure scalar_face_vector
  end interface
  interface face_up
    module procedure scalar_face_vector_up
  end interface
  interface down
    module procedure scalar_down_vector
    module procedure vector_down_vector
  end interface
  interface up
    module procedure scalar_up_vector
    module procedure vector_up_vector
  end interface
  interface ddown
    module procedure scalar_ddown_vector
    module procedure vector_ddown_vector
  end interface
  interface dup
    module procedure scalar_dup_vector
    module procedure vector_dup_vector
  end interface
  interface grad
    module procedure vector_grad_scalar
  end interface
  interface div
    module procedure scalar_div_vector
  end interface
    interface div_dn
    module procedure scalar_div_dn_vector
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

!===============================================================================
FUNCTION vector_norm2_scalar (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)**2
END FUNCTION vector_norm2_scalar
!===============================================================================
FUNCTION vector_norm_scalar (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = sqrt(in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)**2)
END FUNCTION vector_norm_scalar

!===============================================================================
FUNCTION scalar_edge_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ydn(zdn(in))
  out(:,:,:,2) = zdn(xdn(in))
  out(:,:,:,3) = xdn(ydn(in))
END FUNCTION scalar_edge_vector
!===============================================================================
FUNCTION scalar_edge_vector_up (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = yup(zup(in))
  out(:,:,:,2) = zup(xup(in))
  out(:,:,:,3) = xup(yup(in))
END FUNCTION scalar_edge_vector_up
!===============================================================================
FUNCTION scalar_face_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real:: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xdn(in)
  out(:,:,:,2) = ydn(in)
  out(:,:,:,3) = zdn(in)
END FUNCTION scalar_face_vector
!===============================================================================
FUNCTION scalar_face_vector_up (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real:: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xup(in)
  out(:,:,:,2) = yup(in)
  out(:,:,:,3) = zup(in)
END FUNCTION scalar_face_vector_up
!===============================================================================
FUNCTION scalar_down_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xdn (in)
  out(:,:,:,2) = ydn (in)
  out(:,:,:,3) = zdn (in)
END FUNCTION scalar_down_vector
!===============================================================================
FUNCTION vector_down_vector (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xdn (in(:,:,:,1))
  out(:,:,:,2) = ydn (in(:,:,:,2))
  out(:,:,:,3) = zdn (in(:,:,:,3))
END FUNCTION vector_down_vector
!===============================================================================
FUNCTION scalar_up_vector (in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xup (in)
  out(:,:,:,2) = yup (in)
  out(:,:,:,3) = zup (in)
END FUNCTION scalar_up_vector
!===============================================================================
FUNCTION vector_up_vector (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = xup (in(:,:,:,1))
  out(:,:,:,2) = yup (in(:,:,:,2))
  out(:,:,:,3) = zup (in(:,:,:,3))
END FUNCTION vector_up_vector

!===============================================================================
FUNCTION scalar_ddown_vector (ds, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (ds, in)
  out(:,:,:,2) = ddydn (ds, in)
  out(:,:,:,3) = ddzdn (ds, in)
END FUNCTION scalar_ddown_vector

!===============================================================================
FUNCTION vector_ddown_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxdn (ds, in(:,:,:,1))
  out(:,:,:,2) = ddydn (ds, in(:,:,:,2))
  out(:,:,:,3) = ddzdn (ds, in(:,:,:,3))
END FUNCTION vector_ddown_vector

!===============================================================================
FUNCTION scalar_dup_vector (ds, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxup (ds, in)
  out(:,:,:,2) = ddyup (ds, in)
  out(:,:,:,3) = ddzup (ds, in)
END FUNCTION scalar_dup_vector

!===============================================================================
FUNCTION vector_dup_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddxup (ds, in(:,:,:,1))
  out(:,:,:,2) = ddyup (ds, in(:,:,:,2))
  out(:,:,:,3) = ddzup (ds, in(:,:,:,3))
END FUNCTION vector_dup_vector

!===============================================================================
FUNCTION vector_grad_scalar (mesh, in) RESULT(out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:):: mesh
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  ds = mesh%d
  out(:,:,:,1) = ddxdn (ds, in)
  out(:,:,:,2) = ddydn (ds, in)
  out(:,:,:,3) = ddzdn (ds, in)
END FUNCTION
!===============================================================================
FUNCTION scalar_div_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out = ddxup(ds, in(:,:,:,1)) + ddyup(ds, in(:,:,:,2)) + ddzup(ds, in(:,:,:,3))
END FUNCTION
!===============================================================================
FUNCTION scalar_div_dn_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out = ddxdn(ds, in(:,:,:,1)) + ddydn(ds, in(:,:,:,2)) + ddzdn(ds, in(:,:,:,3))
END FUNCTION
!===============================================================================
FUNCTION scalar_minusdiv_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out = - ddxup(ds, in(:,:,:,1)) - ddyup(ds, in(:,:,:,2)) - ddzup(ds, in(:,:,:,3))
END FUNCTION
!===============================================================================
FUNCTION vector_curl_dn_vector (ds, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real(8):: ds(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = ddydn(ds, in(:,:,:,3)) - ddzdn(ds, in(:,:,:,2))
  out(:,:,:,2) = ddzdn(ds, in(:,:,:,1)) - ddxdn(ds, in(:,:,:,3))
  out(:,:,:,3) = ddxdn(ds, in(:,:,:,2)) - ddydn(ds, in(:,:,:,1))
END FUNCTION
!===============================================================================
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
