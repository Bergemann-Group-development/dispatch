!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
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
    logical:: flux_BC = .true.
  contains
    procedure:: init
    procedure:: primary
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t):: self
  integer:: iostat
  logical:: flux_BC = .true.
  real:: theta = 10.0
  real:: mpoly = 1.0
  namelist /boundaries_params/ theta, mpoly, flux_BC
  !-----------------------------------------------------------------------------
  call trace%begin ('boundaries_t%init')
  call self%standard_boundaries_t%init ()
  call io%lock%set ('boundaries_t%init')
  rewind (io%input)
  read (io%input, boundaries_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, boundaries_params)
  end if
  call io%namelist_warning ('boundaries_params', iostat)
  write (io%namelist,boundaries_params)
  call io%lock%unset ('boundaries_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Boundary conditions as specified in the benchmark
!===============================================================================
SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !.............................................................................
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: d, e
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  integer:: ix, iy, iz, n(3), id
  real(8), pointer:: z(:)
  real:: h, f, tt0, tt1, dz, ds(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('boundaries_t%primary', itimer=itimer)
  id    =  self%patch%id
  ds    =  self%patch%ds
  patch => self%patch
  m     => self%mesh
  z     => m(3)%r
  !-----------------------------------------------------------------------------
  ! Allocate scratches
  !-----------------------------------------------------------------------------
  n = m%gn
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (e (n(1),n(2),n(3)))
  !-----------------------------------------------------------------------------
  ! Compute velocities
  !-----------------------------------------------------------------------------
  d => mem(:,:,:,patch%idx%d)
  !-----------------------------------------------------------------------------
  ! Get thermal energy
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = mem(ix,iy,iz,patch%idx%px)/d(ix,iy,iz)
        uy(ix,iy,iz) = mem(ix,iy,iz,patch%idx%py)/d(ix,iy,iz)
        uz(ix,iy,iz) = mem(ix,iy,iz,patch%idx%pz)/d(ix,iy,iz)
        e (ix,iy,iz) = mem(ix,iy,iz,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                        ux(ix,iy,iz)**2 &
                      + uy(ix,iy,iz)**2 &
                      + uz(ix,iy,iz)**2)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Compute the scale heights at the boundaries and fill in the guard zone
  ! values of pressure and density consistent with constant temperature HSE
  !-----------------------------------------------------------------------------
  if (m(3)%lower_boundary) then
    do iy=m(2)%lb,m(2)%ub
      do iz=m(3)%lb,m(3)%lo
        do ix=m(1)%lb,m(1)%ub
          h = self%tt_top/self%g
          dz = (iz-m(3)%li)*ds(3)
          f = exp(dz/h)
          d(ix,iy,iz) = d(ix,iy,m(3)%li)*f
        end do
      end do
      do iz=m(3)%lb,m(3)%li
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = d(ix,iy,iz)*self%tt_top/(eos%gamma-1d0)
        end do
      end do
    end do
  end if
  if (m(3)%upper_boundary) then
    if (self%flux_BC) then
      do iz=m(3)%uo,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            tt0 = e(ix,iy,m(3)%ui)/d(ix,iy,m(3)%ui)*(eos%gamma-1d0)
            dz = (iz-m(3)%ui)*ds(3)
            tt1         = tt0              + self%theta*dz
            d(ix,iy,iz) = d(ix,iy,m(3)%ui) + self%theta*dz
            e(ix,iy,iz) = d(ix,iy,iz)*tt1/(eos%gamma-1d0)
          end do
        end do
      end do
    end if
  end if
  !----------------------------------------------------------------------------
  ! Both top and bottom boundary are identical and are an "inpenetrable wall"
  ! Boundary conditions should be:
  ! 1) antisymmetric for vertical velocity
  ! 2) symmetric for horizontal velocities
  !----------------------------------------------------------------------------
  call self%symmetric     (ux,idx%px)
  call self%symmetric     (uy,idx%py)
  call self%antisymmetric (uz,idx%pz)
  !----------------------------------------------------------------------------
  ! Put the quantities back
  !----------------------------------------------------------------------------
  if (m(3)%lower_boundary) then
    do iz=m(3)%lb,m(3)%li
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,iz,patch%idx%px) =  d(ix,iy,iz)*ux(ix,iy,iz)
          mem(ix,iy,iz,patch%idx%py) =  d(ix,iy,iz)*uy(ix,iy,iz)
          mem(ix,iy,iz,patch%idx%pz) =  d(ix,iy,iz)*uz(ix,iy,iz)
          mem(ix,iy,iz,patch%idx%e ) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                       ux(ix,iy,iz)**2 &
                                     + uy(ix,iy,iz)**2 &
                                     + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  end if
  !
  if (m(3)%upper_boundary) then
    do iz=m(3)%ui,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,iz,patch%idx%px) =  d(ix,iy,iz)*ux(ix,iy,iz)
          mem(ix,iy,iz,patch%idx%py) =  d(ix,iy,iz)*uy(ix,iy,iz)
          mem(ix,iy,iz,patch%idx%pz) =  d(ix,iy,iz)*uz(ix,iy,iz)
        end do
      end do
    end do
    if (self%flux_BC) then
      do iz=m(3)%uo,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            mem(ix,iy,iz,patch%idx%e ) = e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                        ux(ix,iy,iz)**2 &
                                      + uy(ix,iy,iz)**2 &
                                      + uz(ix,iy,iz)**2)
          end do
        end do
      end do
    end if
  end if
  deallocate (ux,uy,uz,e)
  !print *, self%patch%id, minval(d), maxval(d)
  call trace%end (itimer)
END SUBROUTINE primary

END MODULE boundaries_mod
