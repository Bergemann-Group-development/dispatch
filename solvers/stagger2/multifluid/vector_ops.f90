MODULE vector_ops
  USE omp_timer_mod
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
  interface laplacian
    module procedure scalar_laplacian_scalar
  end interface
  interface laplacian_vector
    module procedure vector_laplacian_vector
  end interface
  interface grad_vector
    module procedure tensor_grad_vector
  end interface
  interface divx
    module procedure vector_divx_tensor
  end interface
  interface divy
    module procedure vector_divy_tensor
  end interface
  interface divz
    module procedure vector_divz_tensor
  end interface
  interface hinv_vector
    module procedure vector_hinv_vector
  end interface
CONTAINS

!> =============================================================================
FUNCTION vector_scalar_grad (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
        d2(i,j,k) = mesh(2)%d
        d3(i,j,k) = mesh(3)%d
      end do
    end do
  end do
  out(:,:,:,1) = ddxdn (d1, in)
  out(:,:,:,2) = ddydn (d2, in)
  out(:,:,:,3) = ddzdn (d3, in)
END FUNCTION vector_scalar_grad

!> =============================================================================
FUNCTION vector_norm2_scalar (in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  !-----------------------------------------------------------------------------
  out = in(:,:,:,1)**2 + in(:,:,:,2)**2 + in(:,:,:,3)**2
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
FUNCTION scalar_ddown_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
        d2(i,j,k) = mesh(2)%d
        d3(i,j,k) = mesh(3)%d
      end do
    end do
  end do
  out(:,:,:,1) = ddxdn (d1, in)
  out(:,:,:,2) = ddydn (d2, in)
  out(:,:,:,3) = ddzdn (d3, in)
END FUNCTION scalar_ddown_vector
!> =============================================================================
FUNCTION vector_dup_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
        d2(i,j,k) = mesh(2)%d
        d3(i,j,k) = mesh(3)%d
      end do
    end do
  end do
  out(:,:,:,1) = ddxup (d1, in(:,:,:,1))
  out(:,:,:,2) = ddyup (d2, in(:,:,:,2))
  out(:,:,:,3) = ddzup (d3, in(:,:,:,3))
END FUNCTION vector_dup_vector

!> =============================================================================
FUNCTION vector_grad_scalar (mesh, in) RESULT(out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
        d2(i,j,k) = mesh(1)%h2c(i) *mesh(2)%d
        d3(i,j,k) = mesh(1)%h31c(i)*mesh(2)%h32c(j)*mesh(3)%d
      end do
    end do
  end do
  out(:,:,:,1) = ddxdn (d1, in)
  out(:,:,:,2) = ddydn (d2, in)
  out(:,:,:,3) = ddzdn (d3, in)
END FUNCTION vector_grad_scalar
!> =============================================================================
FUNCTION scalar_div_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: n1, n2, n3, d1, d2, d3
  class(mesh_t), dimension(:), pointer:: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n1(i,j,k) = mesh(1)%h2f(i)   *mesh(1)%h31f(i)
        n2(i,j,k) = mesh(2)%h32f(j)  *mesh(1)%h31c(i)  *mesh(1)%d
        n3(i,j,k) = mesh(1)%h2c(i)   *mesh(1)%d        *mesh(2)%d
        d1(i,j,k) = mesh(1)%dvol1c(min(i,n(1)-1))
        d2(i,j,k) = mesh(2)%dvol2c(min(j,n(2)-1))*mesh(1)%dvol1c(min(i,n(1)-1))
        d3(i,j,k) = mesh(3)%dvol3c(min(k,n(3)-1))*mesh(1)%dvol1c(min(i,n(1)-1))*mesh(2)%dvol2c(min(j,n(2)-1))
      end do
    end do
  end do
  out =   ddxup( d1, n1*in(:,:,:,1) ) &
        + ddyup( d2, n2*in(:,:,:,2) ) &
        + ddzup( d3, n3*in(:,:,:,3) )
END FUNCTION scalar_div_vector
!> =============================================================================
FUNCTION scalar_minusdiv_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
    real, dimension(size(in,1),size(in,2),size(in,3)):: n1, n2, n3, d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n1(i,j,k) = mesh(1)%h2f(i)   *mesh(1)%h31f(i)
        n2(i,j,k) = mesh(2)%h32f(j)  *mesh(1)%h31c(i)  *mesh(1)%d
        n3(i,j,k) = mesh(1)%h2c(i)   *mesh(1)%d        *mesh(2)%d
        d1(i,j,k) = mesh(1)%dvol1c(min(i,n(1)-1))
        d2(i,j,k) = mesh(2)%dvol2c(min(j,n(2)-1))*mesh(1)%dvol1c(min(i,n(1)-1))
        d3(i,j,k) = mesh(3)%dvol3c(min(k,n(3)-1))*mesh(1)%dvol1c(min(i,n(1)-1))*mesh(2)%dvol2c(min(j,n(2)-1))
      end do
    end do
  end do
  out = - ddxup( d1, n1*in(:,:,:,1) ) &
        - ddyup( d2, n2*in(:,:,:,2) ) &
        - ddzup( d3, n3*in(:,:,:,3) )
END FUNCTION scalar_minusdiv_vector
!> =============================================================================
FUNCTION vector_curl_dn_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: n12, n13,      n21, n31
  real, dimension(size(in,1),size(in,2),size(in,3)):: d12, d13, d23, d21, d31, d32
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n12(i,j,k) = mesh(2)%h32c(j)
        n13(i,j,k) = mesh(2)%d
        n21(i,j,k) = mesh(1)%h31c(i)
        n31(i,j,k) = mesh(1)%h2c(i)
        d12(i,j,k) = mesh(1)%h2c(i) *mesh(2)%dvol2f(j)
        d13(i,j,k) = mesh(1)%h31c(i)*mesh(3)%dvol3f(k)
        d23(i,j,k) = mesh(1)%h31f(i)*mesh(2)%h32c(j)  *mesh(3)%d
        d21(i,j,k) = mesh(1)%h31f(i)*mesh(1)%d
        d31(i,j,k) = mesh(1)%h2f(i) *mesh(1)%d
        d32(i,j,k) = mesh(1)%h2f(i) *mesh(2)%d
      end do
    end do
  end do
  out(:,:,:,1) =   ddydn( d12, n12*in(:,:,:,3) ) &
                 - ddzdn( d13, n13*in(:,:,:,2) )
  out(:,:,:,2) =   ddzdn( d23,     in(:,:,:,1) ) &
                 - ddxdn( d21, n21*in(:,:,:,3) )
  out(:,:,:,3) =   ddxdn( d31, n31*in(:,:,:,2) ) &
                 - ddydn( d32,     in(:,:,:,1) )
END FUNCTION vector_curl_dn_vector
!> =============================================================================
FUNCTION vector_curl_up_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: n12, n13,      n21, n31
  real, dimension(size(in,1),size(in,2),size(in,3)):: d12, d13, d23, d21, d31, d32
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n12(i,j,k) = mesh(2)%h32f(j)
        n13(i,j,k) = mesh(2)%d
        n21(i,j,k) = mesh(1)%h31f(i)
        n31(i,j,k) = mesh(1)%h2f(i)
        d12(i,j,k) = mesh(1)%h2f(i) *mesh(2)%dvol2c(j)
        d13(i,j,k) = mesh(1)%h31f(i)*mesh(3)%dvol3c(k)
        d23(i,j,k) = mesh(1)%h31c(i)*mesh(2)%h32f(j)  *mesh(3)%d
        d21(i,j,k) = mesh(1)%h31c(i)*mesh(1)%d
        d31(i,j,k) = mesh(1)%h2c(i) *mesh(1)%d
        d32(i,j,k) = mesh(1)%h2c(i) *mesh(2)%d
      end do
    end do
  end do
  out(:,:,:,1) =   ddyup( d12, n12*in(:,:,:,3) ) &
                 - ddzup( d13, n13*in(:,:,:,2) )
  out(:,:,:,2) =   ddzup( d23,     in(:,:,:,1) ) &
                 - ddxup( d21, n21*in(:,:,:,3) )
  out(:,:,:,3) =   ddxup( d31, n31*in(:,:,:,2) ) &
                 - ddyup( d32,     in(:,:,:,1) )
END FUNCTION vector_curl_up_vector
!> =============================================================================
FUNCTION scalar_laplacian_scalar (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:), pointer:: mesh
  !-----------------------------------------------------------------------------
  out = div(mesh,grad(mesh,in))
END FUNCTION scalar_laplacian_scalar
!> =============================================================================
FUNCTION vector_laplacian_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:), pointer:: mesh
 !-----------------------------------------------------------------------------
  out = grad(mesh,div(mesh,in)) - curl_up(mesh,curl_dn(mesh,in))
END FUNCTION vector_laplacian_vector
!> =============================================================================
FUNCTION tensor_grad_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3,3):: out
  real, dimension(size(in,1),size(in,2),size(in,3)):: h31f, h31c, h2f, h2c, d1, &
                                                      h32f, h32c, d2 , d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1  (i,j,k) = mesh(1)%d
        d2  (i,j,k) = mesh(2)%d
        d3  (i,j,k) = mesh(3)%d
        h2f (i,j,k) = mesh(1)%h2f(i)
        h2c (i,j,k) = mesh(1)%h2c(i)
        h31f(i,j,k) = mesh(1)%h31f(i)
        h31c(i,j,k) = mesh(1)%h31c(i)
        h32f(i,j,k) = mesh(2)%h32f(j)
        h32c(i,j,k) = mesh(2)%h32c(j)
      end do
    end do
  end do
  out(:,:,:,1,1) = ddxup(d1          ,in(:,:,:,1))
  out(:,:,:,1,2) = ddxdn(d1          ,in(:,:,:,2))
  out(:,:,:,1,3) = ddxdn(d1          ,in(:,:,:,3))
  out(:,:,:,2,1) = ddydn(h2f      *d2,in(:,:,:,1)) - xdn(in(:,:,:,2))*ddxdn(h2f*d1,h2c)
  out(:,:,:,2,2) = ddyup(h2c      *d2,in(:,:,:,2)) + xup(in(:,:,:,1))*ddxup(h2c*d1,h2f)
  out(:,:,:,2,3) = ddydn(h2c      *d2,in(:,:,:,3))
  out(:,:,:,3,1) = ddzdn(h31f*h32c*d3,in(:,:,:,1)) - xdn(in(:,:,:,3))*ddxdn(h31f*    d1,h31c)
  out(:,:,:,3,2) = ddzdn(h31c*h32f*d3,in(:,:,:,2)) - ydn(in(:,:,:,3))*ddydn(h32f*h2c*d2,h32c)
  out(:,:,:,3,3) = ddzup(h31c*h32c*d3,in(:,:,:,3)) + xup(in(:,:,:,1))*ddxup(h31c*    d1,h31f) &
                                                   + yup(in(:,:,:,2))*ddyup(h32c*h2c*d2,h32f)
END FUNCTION tensor_grad_vector
!> =============================================================================
FUNCTION vector_divx_tensor (mesh, in11, in21, in31, in22, in33) RESULT (out)
  real, dimension(:,:,:), intent(in):: in11, in21, in31, in22, in33
  real, dimension(size(in11,1),size(in11,2),size(in11,3)):: out
  real, dimension(size(in11,1),size(in11,2),size(in11,3)):: n1, n2, n3, n4, n5
  real, dimension(size(in11,1),size(in11,2),size(in11,3)):: d1, d2, d3, d4, d5
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in11,1),size(in11,2),size(in11,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n1(i,j,k) = mesh(1)%h2c(i)   *mesh(1)%h31c(i)
        n2(i,j,k) = mesh(2)%h32f(j)  *mesh(1)%h31f(i)  *mesh(1)%d
        n3(i,j,k) = mesh(1)%h2f(i)   *mesh(1)%d        *mesh(2)%d
        n4(i,j,k) = mesh(1)%h2c(i)
        n5(i,j,k) = mesh(1)%h31c(i)
        d1(i,j,k) = mesh(1)%dvol1f(max(i,2))
        d2(i,j,k) = mesh(2)%dvol2c(min(j,n(2)-1))*mesh(1)%dvol1f(max(i,2))
        d3(i,j,k) = mesh(3)%dvol3c(min(k,n(3)-1))*mesh(1)%dvol1f(max(i,2))*mesh(2)%dvol2c(min(j,n(2)-1))
        d4(i,j,k) = mesh(1)%h2f(i)   *mesh(1)%d
        d5(i,j,k) = mesh(1)%h31f(i)  *mesh(1)%d
      end do
    end do
  end do
  out =   ddxdn( d1, n1*in11 ) &
        + ddyup( d2, n2*in21 ) &
        + ddzup( d3, n3*in31 ) &
        - ddxdn( d4, n4      ) * in22 & !xdn(in22) & !exp(xdn(log(in22))) &
        - ddxdn( d5, n5      ) * in33   !xdn(in33)   !exp(xdn(log(in33)))
END FUNCTION vector_divx_tensor
!> =============================================================================
FUNCTION vector_divy_tensor (mesh, in12, in22, in32, in33) RESULT (out)
  real, dimension(:,:,:), intent(in):: in12, in22, in32, in33
  real, dimension(size(in12,1),size(in12,2),size(in12,3)):: out
  real, dimension(size(in12,1),size(in12,2),size(in12,3)):: n1, n2, n3, n4
  real, dimension(size(in12,1),size(in12,2),size(in12,3)):: d1, d2, d3, d4
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in12,1),size(in12,2),size(in12,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n1(i,j,k) = (mesh(1)%h2f(i)**2)*mesh(1)%h31f(i)
        n2(i,j,k) = mesh(2)%h32c(j)  *mesh(1)%h31c(i)  *mesh(1)%d*mesh(1)%h2c(i)
        n3(i,j,k) = mesh(1)%h2c(i)   *mesh(1)%d        *mesh(2)%d*mesh(1)%h2c(i)
        n4(i,j,k) = mesh(2)%h32c(j)
        d1(i,j,k) = mesh(1)%dvol1c(min(i,n(1)-1))
        d2(i,j,k) = mesh(2)%dvol2f(max(j,2))*mesh(1)%dvol1c(min(i,n(1)-1))
        d3(i,j,k) = mesh(3)%dvol3c(min(k,n(3)-1))*mesh(1)%dvol1c(min(i,n(1)-1))*mesh(2)%dvol2f(max(j,2))
        d4(i,j,k) = mesh(2)%h32f(j)  *mesh(2)%d
      end do
    end do
  end do
  out =   ddxup( d1, n1*in12 ) &
        + ddydn( d2, n2*in22 ) &
        + ddzup( d3, n3*in32 ) &
        - ddydn( d4, n4      ) * in33 !ydn(in33) !exp(ydn(log(in33)))
END FUNCTION vector_divy_tensor
!> =============================================================================
FUNCTION vector_divz_tensor (mesh, in13, in23, in33) RESULT (out)
  real, dimension(:,:,:), intent(in):: in13, in23, in33
  real, dimension(size(in13,1),size(in13,2),size(in13,3)):: out
  real, dimension(size(in13,1),size(in13,2),size(in13,3)):: n1, n2, n3
  real, dimension(size(in13,1),size(in13,2),size(in13,3)):: d1, d2, d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in13,1),size(in13,2),size(in13,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        n1(i,j,k) = mesh(1)%h2f(i)   *mesh(2)%h32c(j)  *mesh(1)%h31f(i)**2
        n2(i,j,k) = mesh(2)%h32f(j)  *mesh(1)%h31c(i)  *mesh(1)%d      *mesh(1)%h31c(i)*mesh(2)%h32f(j)
        n3(i,j,k) = mesh(1)%h2c(i)   *mesh(1)%d        *mesh(2)%d      *mesh(1)%h31c(i)*mesh(2)%h32c(j)
        d1(i,j,k) = mesh(1)%dvol1c(min(i,n(1)-1))
        d2(i,j,k) = mesh(2)%dvol2c(min(j,n(2)-1))*mesh(1)%dvol1c(min(i,n(1)-1))
        d3(i,j,k) = mesh(3)%dvol3f(max(k,2))*mesh(1)%dvol1c(min(i,n(1)-1))*mesh(2)%dvol2c(min(j,n(2)-1))
      end do
    end do
  end do
  out =   ddxup( d1, n1*in13 ) &
        + ddyup( d2, n2*in23 ) &
        + ddzdn( d3, n3*in33 )
END FUNCTION vector_divz_tensor
!> =============================================================================
FUNCTION vector_hinv_vector (mesh, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k,1) = in(i,j,k,1)
        out(i,j,k,2) = in(i,j,k,2)/mesh(1)%h2c(i)
        out(i,j,k,3) = in(i,j,k,3)/mesh(1)%h31c(i)/mesh(2)%h32c(j)
      end do
    end do
  end do
END FUNCTION vector_hinv_vector
!*******************************************************************************
END MODULE
