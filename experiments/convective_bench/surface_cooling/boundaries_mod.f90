!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!>
!> For velocity components boundary conditions are:
!> 1) antisymmetric for vertical velocity
!> 2) symmetric for horizontal velocities
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE patch_mod
  USE eos_mod
  USE benchmark_mod
  USE scratch_mod
  !
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    !---------------------------------------------------------------------------
    ! Experiment - specific variables
    !---------------------------------------------------------------------------.
    real:: tt_top = 1.0
    real:: g = 20.
    real:: theta = 10.0
    real:: mpoly = 1.0
  contains
    procedure:: init
    procedure:: primary
  end type
  real,save:: bc_rate=0.5
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t):: self
  integer:: iostat
  real:: theta = 10.0
  real:: mpoly = 1.0
  !-----------------------------------------------------------------------------
  call trace%begin ('boundaries_t%init')
  self%set_bc_bits=.false.
  call self%standard_boundaries_t%init ()
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Boundary conditions as specified in the benchmark
!===============================================================================
SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !.............................................................................
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: d, e, s, eos_d, eos_e, eos_s
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: pg, tt, q
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: px, py, pz
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  integer:: ix, iy, iz, iz1, n(3), id, lb, ub, l, u
  real(8), pointer:: z(:)
  real:: h, f, dz, ds(3), d0, pg0, tt0, d1, pg1, tt1, dtt, dss
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('boundaries_t%primary', itimer=itimer)
  id    =  self%patch%id
  ds    =  self%patch%ds
  patch => self%patch
  m     => self%mesh
  allocate (z(patch%gn(3)))
  z = patch%mesh(3)%r + patch%mesh(3)%p

  !-----------------------------------------------------------------------------
  ! Allocate scratches
  !-----------------------------------------------------------------------------
  n = m%gn
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (e (n(1),n(2),n(3)))
  !-----------------------------------------------------------------------------
  ! Point pointers
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,patch%idx%d)
  px => mem(:,:,:,patch%idx%px)
  py => mem(:,:,:,patch%idx%py)
  pz => mem(:,:,:,patch%idx%pz)
  !-----------------------------------------------------------------------------
  ! Get velocities
  !-----------------------------------------------------------------------------
  l=m(1)%li
  u=m(1)%ui
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = px(ix,iy,iz)/d(ix,iy,iz)
        uy(ix,iy,iz) = py(ix,iy,iz)/d(ix,iy,iz)
        uz(ix,iy,iz) = pz(ix,iy,iz)/d(ix,iy,iz)
      end do
    end do
    ! if (any(px(l:u,l:u,iz) /= px(l,l,iz))) print*,patch%id,iz,'px'
    ! if (any(py(l:u,l:u,iz) /= py(l,l,iz))) print*,patch%id,iz,'py'
    ! if (any(pz(l:u,l:u,iz) /= pz(l,l,iz))) print*,patch%id,iz,'pz'
    ! if (any(ux(l:u,l:u,iz) /= ux(l,l,iz))) print*,patch%id,iz,'ux'
    ! if (any(uy(l:u,l:u,iz) /= uy(l,l,iz))) print*,patch%id,iz,'uy'
    ! if (any(uz(l:u,l:u,iz) /= uz(l,l,iz))) print*,patch%id,iz,'uz'
  end do
  !-----------------------------------------------------------------------------
  ! Get thermal energy, if applicable
  !-----------------------------------------------------------------------------
  if (patch%idx%e > 0) then
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,iz,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                         ux(ix,iy,iz)**2 &
                       + uy(ix,iy,iz)**2 &
                       + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (patch%idx%s > 0) then
    s => mem(:,:,:,patch%idx%s)
    e = d**eos%gamma*exp(s*(eos%gamma-1d0)/d)/(eos%gamma-1d0)
    ! do iz=m(3)%lb,m(3)%ub
    !   if (any(d(l:u,l:u,iz) /= d(l,l,iz))) print*,patch%id,iz,'d'
    !   if (any(s(l:u,l:u,iz) /= s(l,l,iz))) print*,patch%id,iz,'s'
    !   if (any(e(l:u,l:u,iz) /= e(l,l,iz))) print*,patch%id,iz,'e'
    ! end do
  end if
  !-----------------------------------------------------------------------------
  ! The lb and ub indices are the index limits for BC modified values.  The BCs
  ! are placed at cell centers -- using face centered BCs does not work well
  !-----------------------------------------------------------------------------
  lb = m(3)%li
  if (m(3)%lower_boundary) then
    call benchmark%top (z(:lb), &
                     d(:,:,:lb), &
                     e(:,:,:lb), &
                    ux(:,:,:lb), &
                    uy(:,:,:lb), &
                    uz(:,:,:lb))
  end if
  !-----------------------------------------------------------------------------
  ! For pressure and temperature at the bottom boundary:  Leave the initial
  ! values in the T-case, and impose the given gradient in the F-case
  !-----------------------------------------------------------------------------
  ub = m(3)%ui
  if (m(3)%upper_boundary) then
    call benchmark%bottom (z(ub:), &
                        d(:,:,ub:), &
                        e(:,:,ub:), &
                       ux(:,:,ub:), &
                       uy(:,:,ub:), &
                       uz(:,:,ub:))
  end if
  !-----------------------------------------------------------------------------
  ! Put the total energy back into mem, if applicable
  !-----------------------------------------------------------------------------
  if (patch%idx%e>0) then
    if (m(3)%lower_boundary) then
      do iz=m(3)%lb,lb
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            px (ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
            py (ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
            pz (ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
            mem(ix,iy,iz,patch%idx%e) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                        ux(ix,iy,iz)**2 &
                                      + uy(ix,iy,iz)**2 &
                                      + uz(ix,iy,iz)**2)
          end do
        end do
      end do
    end if
    if (m(3)%upper_boundary) then
      do iz=ub,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            px (ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
            py (ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
            pz (ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
            mem(ix,iy,iz,patch%idx%e) = e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                       ux(ix,iy,iz)**2 &
                                     + uy(ix,iy,iz)**2 &
                                     + uz(ix,iy,iz)**2)
          end do
        end do
      end do
    end if
  else if (patch%idx%s > 0) then
    if (m(3)%lower_boundary .or. m(3)%upper_boundary) then
      px = d*ux
      py = d*uy
      pz = d*uz
      s = log((eos%gamma-1d0)*e/d**eos%gamma)*d/(eos%gamma-1d0)
      ! write (stdout,'(a,i6,1p,4e11.3)') &
      ! 'BC: id,dmin,dmax,smin,smax=', patch%id, &
      ! minval(d), maxval(d), minval(s), maxval(s)
    else
      !px = d*ux
      !py = d*uy
      !pz = d*uz
      s = log((eos%gamma-1d0)*e/d**eos%gamma)*d/(eos%gamma-1d0)
      ! l=m(1)%lb
      ! u=m(1)%ub
      ! do iz=m(3)%li,m(3)%ui
      !   if (any(px(l:u,l:u,iz)/=px(l,l,iz))) then
      !     print *,patch%id,iz,patch%position(3),'px3'
      !     print *, minval(px(l:u,l:u,iz)-px(l,l,iz)), maxval(px(l:u,l:u,iz)-px(l,l,iz))
      !     print *, minloc(px(l:u,l:u,iz)-px(l,l,iz)), maxloc(px(l:u,l:u,iz)-px(l,l,iz))
      !   end if
      !   if (any(py(l:u,l:u,iz)/=py(l,l,iz))) print *,patch%id,iz,patch%position(3),'py3'
      !   if (any(pz(l:u,l:u,iz)/=pz(l,l,iz))) print *,patch%id,iz,patch%position(3),'pz3'
      !   if (any( d(l:u,l:u,iz)/= d(l,l,iz))) print *,patch%id,iz,patch%position(3),'d3'
      !   !if (any( s(l:u,l:u,iz)/= s(l,l,iz))) print *,patch%id,iz,patch%position(3),'s3'
      !   if (any( s(l:u,l:u,iz)/= s(l,l,iz))) then
      !     print *,patch%id,iz,patch%position(3),'s3'
      !     print *, minval(s(l:u,l:u,iz)-s(l,l,iz)), maxval(s(l:u,l:u,iz)-s(l,l,iz))
      !     print *, minloc(s(l:u,l:u,iz)-s(l,l,iz)), maxloc(s(l:u,l:u,iz)-s(l,l,iz))
      !   end if
      ! end do
    end if
  end if
  deallocate (z,e,ux,uy,uz)

  !  define filling

  if (associated(patch%filled)) then
    do iz=m(3)%uo,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          patch%filled(ix,iy,iz) = patch%level
        end do
      end do
    end do
  end if
  
  if (associated(patch%filled)) then
    do iz=m(3)%lb,m(3)%lo
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          patch%filled(ix,iy,iz) = patch%level
        end do
      end do
    end do
  end if

  call trace%end (itimer)
END SUBROUTINE primary

END MODULE boundaries_mod
