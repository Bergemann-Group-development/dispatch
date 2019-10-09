MODULE vector_ops
  USE trace_mod
  USE mesh_mod
  USE kinds_mod
  implicit none
  public

  interface curl_dn
    module procedure vector_curl_dn_vector
  end interface
  interface curl_up
    module procedure vector_curl_up_vector
  end interface
  interface curl_emfs
    module procedure dbdt_curl_up_emfs
  end interface
  interface grad
    module procedure grad_dn_scalar
  end interface
  interface diverg
    module procedure divergence_up_vector
  end interface

CONTAINS
!===============================================================================
!> Calculate the gradient of a zone-centred quantity. The result is face-centred.
!> e.g. pressure gradient, gravitational acceleration
!===============================================================================
FUNCTION grad_dn_scalar (mesh, in) RESULT (out)
  class(mesh_t), dimension(:):: mesh
  real(kind=KindScalarVar), dimension(:,:,:), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3),3):: out
  integer:: n(3), ione, jone, kone, i, j, k, im1, jm1, km1
  !-----------------------------------------------------------------------------
  call trace%begin('vector_ops::grad_dn_scalar')

  n = shape(in)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)

  associate (m1 =>mesh(1), m2=>mesh(2), m3=>mesh(3))
  do k=2,n(3)-kone
    km1 = k - kone
    do j=2,n(2)-jone
      jm1 = j - jone
      do i=2,n(1)-ione
        im1 = i - ione
        out(i,j,k,1) = (in(i,j,k) - in(im1,j,k)) / m1%d
        out(i,j,k,2) = (in(i,j,k) - in(i,jm1,k)) / (m1%h2c(i) * m2%d)
        out(i,j,k,3) = (in(i,j,k) - in(i,j,km1)) / (m1%h31c(i) * m2%h32c(j) * m3%d)
      end do
    end do
  end do

  end associate
  call trace%end()
END FUNCTION grad_dn_scalar

!===============================================================================
!> Calculate the curl of a face-centred quantity. Result is edge-centred.
!> e.g. The current density from the magnetic field.
!===============================================================================
FUNCTION vector_curl_dn_vector (mesh, in) RESULT (out)
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:):: mesh
  
  integer:: i, j, k, ione, jone, kone, im1, jm1, km1, n(4)
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin('vector_ops::vector_curl_dn_vector')

  n = shape(in)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)

  associate (m1 => mesh(1), m2 => mesh(2), m3 => mesh(3))

  ! 1-component
  do k=2,n(3)
    km1 = k - kone
    do j=2,n(2)
      jm1 = j - jone
      do i=1,n(1)
        out(i,j,k,1) = (m2%h32c(j) * in(i,j,k,3) - m2%h32c(jm1) * in(i,jm1,k,3)) &
                     / (m1%h2c(i) * m2%dvol2c(j)) &
                     - (in(i,j,k,2) - in(i,j,km1,2)) &
                     * m2%dar32c(j) / (m1%h31c(i) * m3%d)
      end do
    end do
  end do

  ! 2-component
  do k=2,n(3)
    km1 = k - kone
    do j=1,n(2)
      do i=2,n(1)
        im1 = i - ione
        out(i,j,k,2) = (in(i,j,k,1) - in(i,j,km1,1)) &
                     / (m1%h31f(i) * m2%h32c(j) * m3%d) &
                     - (m1%h31c(i) * in(i,j,k,3) - m1%h31c(im1) * in(i,j,k,3)) &
                     / (m1%h31f(i) * m1%d)
      end do
    end do
  end do

  ! 3-component
  do k=1,n(3)
    do j=2,n(2)
      jm1 = j - jone
      do i=2,n(1)
        im1 = i - ione
        out(i,j,k,3) = (m1%h2c(i) * in(i,j,k,2) - m1%h2c(im1) * in(im1,j,k,2)) &
                     / (m1%h2f(i) * m1%d) &
                     - (in(i,j,k,1) - in(i,jm1,k,1)) &
                     / (m1%h2f(i) * m2%d)
      end do
    end do
  end do

  end associate
  call trace%end()
END FUNCTION vector_curl_dn_vector

!===============================================================================
!> Calculate the curl of a zone-centred quantity. Result is face-centred.
!> e.g. Time rate change of magnetic field from the EMFs; the magnetic field
!> from the vector potential.
!> The difference between this and `vector_curl_up_emfs` is that the definition
!> of EMFs in the latter include metric and length factors, e.g. `E_i * g_i * dx_i`.
!===============================================================================
FUNCTION vector_curl_up_vector (mesh, in) RESULT (out)
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:):: mesh
  
  integer:: i, j, k, ione, jone, kone, ip1, jp1, kp1, n(4)
  !-----------------------------------------------------------------------------
  call trace%begin('vector_ops::vector_curl_up_vector')

  n = shape(in)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)

  associate (m1 =>mesh(1), m2=>mesh(2), m3=>mesh(3))

  ! 1-component
  do k=1,n(3)-1
    kp1 = k + kone
    do j=1,n(2)-1
      jp1 = j + jone
      do i=1,n(1)
        out(i,j,k,1) = (m2%h32f(jp1) * in(i,jp1,k,3) - m2%h32f(j) * in(i,j,k,3)) &
                     / (m1%h2f(i) * m2%dvol2f(j)) &
                     - (in(i,j,kp1,2) - in(i,j,k,2)) &
                     * m2%dar32f(j) / (m1%h31f(i) * m3%dvol3f(k))
      end do
    end do
  end do
  
  ! 2-component
  do k=1,n(3)-1
    kp1 = k + kone
    do j=1,n(2)
      do i=1,n(1)-1
        ip1 = i + ione
        out(i,j,k,2) = (in(i,j,kp1,1) - in(i,j,k,1)) &
                     / (m1%h31c(i) * m2%h32f(j) * m3%d) &
                     - (m1%h31f(ip1) * in(ip1,j,k,3) - m1%h31f(i) * in(i,j,k,3)) &
                     / (m1%h31c(i) * m1%d)
      end do
    end do
  end do

  ! 3-component
  do k=1,n(3)
    do j=1,n(2)-1
      jp1 = j + jone
      do i=1,n(1)-1
        ip1 = i + ione
        out(i,j,k,3) = (m1%h2f(ip1) * in(ip1,j,k,2) - m1%h2f(i) * in(i,j,k,2)) &
                     / (m1%h2c(i) * m1%d) &
                     - (in(i,jp1,k,1) - in(i,j,k,1)) &
                     / (m1%h2c(i) * m2%d)
      end do
    end do
  end do

  end associate
  call trace%end()
END FUNCTION vector_curl_up_vector

!===============================================================================
!> Calculate the curl of a zone-centred quantity. Result is face-centred.
!> e.g. Time rate change of magnetic field from the EMFs; the magnetic field
!> from the vector potential.
!> The difference between this and `vector_curl_up_vector` is that the latter is
!> more generic and does *not* assume that the EMFs alrady include metric and
!> length factors.
!===============================================================================
FUNCTION dbdt_curl_up_emfs (mesh, in) RESULT (out)
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3),3):: out
  class(mesh_t), dimension(:):: mesh
  
  integer:: i, j, k, ione, jone, kone, ip1, jp1, kp1, n(4)
  !-----------------------------------------------------------------------------
  call trace%begin('vector_ops::vector_curl_up_vector')

  n = shape(in)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)

  associate (m1 =>mesh(1), m2=>mesh(2), m3=>mesh(3))

  ! REMEMBER: As stored, the EMFs include metric factors and grid spacings,
  ! e.g., `E_i * g_i * dx_i`, which is why the curl looks a bit funny here.
  ! This form mimics the constraint transport implementation (CT) in ZEUS.
  ! 1-component
  do k=1,max(1,n(3)-1)
    kp1 = k + kone
    do j=1,max(1,n(2)-1)
      jp1 = j + jone
      do i=1,n(1)
        out(i,j,k,1) = (in(i,jp1,k,3) - in(i,j,k,3) - in(i,j,kp1,2) + in(i,j,k,2)) &
                     / (m1%h2f(i) * m1%h31f(i) * m2%dvol2c(j) * m3%d)
      end do
    end do
  end do
  
  ! 2-component
  do k=1,max(1,n(3)-1)
    kp1 = k + kone
    do j=1,n(2)
      do i=1,max(1,n(1)-1)
        ip1 = i + ione
        out(i,j,k,2) = (in(i,j,kp1,1) - in(i,j,k,1) - in(ip1,j,k,3) + in(i,j,k,3)) &
                     / (m1%h31c(i) * m1%d * m2%h32f(j) * m3%d)
      end do
    end do
  end do

  ! 3-component
  do k=1,n(3)
    do j=1,max(1,n(2)-1)
      jp1 = j + jone
      do i=1,max(1,n(1)-1)
        ip1 = i + ione
        out(i,j,k,3) = (in(ip1,j,k,2) - in(i,j,k,2) - in(i,jp1,k,1) + in(i,j,k,1)) &
                     / (m1%h2c(i) * m1%d * m2%d)
      end do
    end do
  end do

  end associate
  call trace%end()
END FUNCTION dbdt_curl_up_emfs

!===============================================================================
!> Calculate the divergence of a face-centred quantity (e.g. velocities or
!> magnetic field. This function is derived heavily from `diverg` in ZEUS-3D,
!> except that it does not perform any kind of normalization of the result.
!===============================================================================
FUNCTION divergence_up_vector (mesh, in) RESULT (out)
  real(kind=KindScalarVar), dimension(:,:,:,:), intent(in):: in
  real(kind=KindScalarVar), dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh

  integer:: i, j, k, ione, jone, kone, ip1, jp1, kp1, n(4)!-----------------------------------------------------------------------------
  call trace%begin('vector_ops::divergence_up_vector')

  n = shape(in)
  ione = merge(1, 0, n(1) > 1)
  jone = merge(1, 0, n(2) > 1)
  kone = merge(1, 0, n(3) > 1)

  associate (m1 =>mesh(1), m2=>mesh(2), m3=>mesh(3))
  do k=1,n(3)-kone
    kp1 = k + kone
    do j=1,n(2)-jone
      jp1 = j + jone
      do i=1,n(1)-ione
        ip1 = i + ione
        out(i,j,k) = (m1%dar1c(ip1) * in(ip1,j,k,1) - m1%dar1c(i) * in(i,j,k,1)) / m1%dvol1f(i) &
                   + (m2%h32f (jp1) * in(i,jp1,k,2) - m2%h32f (j) * in(i,j,k,2)) &
                     * m1%dar2f(i) / m2%dvol2f(j) &
                   + (in(i,j,kp1,3) - in(i,j,k,3)) * m1%dar31f(i) * m2%dar32f(j) / m3%dvol3f(k)
      end do
    end do
  end do

  end associate
  call trace%end()
END FUNCTION divergence_up_vector

END MODULE vector_ops
