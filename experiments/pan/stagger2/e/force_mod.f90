!*******************************************************************************
!> $Id$
!*******************************************************************************
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  implicit none
  private
  type, public:: force_t
    integer :: id = 0
    integer :: seed
    integer :: verbose
    real :: t_turb = 0.0
    real :: t_damp = 0.0
    real, allocatable, dimension(:,:,:,:,:)   :: fran
    procedure(turbulence), pointer :: selected=>null()
    class(mesh_t), dimension(:), pointer  :: mesh
    type(random_t):: random
  contains
    procedure :: init
    procedure :: dealloc
  end type
  real :: ampl_turb = 0.1
  real :: a_helmh=1.
  real :: k1=1.
  real :: k2=1.5
  real :: pk=7./6.
  real :: t_turn = 0.1
  real :: t_turb = 0.0
  real :: t_damp = 3e-4
  real :: pi=3.14159365
  logical :: do_force = .true.
  logical :: do_helmh = .true.
  integer :: verbose = 0
  integer :: seed = -3
CONTAINS

!===============================================================================
!> Initialize the force module
!===============================================================================
SUBROUTINE init (self, solver, id, mesh)
  class(force_t):: self
  integer:: id
  character(len=64):: solver
  class(mesh_t), dimension(:), pointer:: mesh
  !.............................................................................
  integer :: mx, my, mz
  real    :: rhoav, xav, yav, zav
  namelist /turbulence_params/ do_force, do_helmh, a_helmh, k1, k2, pk, &
                               ampl_turb, t_turn, t_turb, seed, t_damp, verbose
  logical, save:: first_time = .true.
  !---------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input,turbulence_params)
    if (io%master) write (*,turbulence_params)
  end if
  !$omp end critical (init_cr)
  if (t_turb == 0.) t_turb = -1.5*t_turn
  self%id     = id
  self%seed   = seed
  self%t_turb = t_turb
  self%t_damp = t_damp
  self%verbose= verbose
  self%mesh   => mesh
  mx = self%mesh(1)%gn
  my = self%mesh(2)%gn
  mz = self%mesh(3)%gn
  allocate (self%fran(mx,my,mz,3,2))
  call io%gb_mem (storage_size(self%fran)*(product(shape(self%fran))/(8.*1024.**3)))
  self%fran = 0.0
  self%selected => turbulence
  call self%random%init (seed)
END SUBROUTINE init

!===============================================================================
!> Deallocate force update arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(force_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('force_t%dealloc')
  if (allocated(self%fran)) then
    call io%bits_mem (-storage_size(self%fran), product(shape(self%fran)), '-fran')
    deallocate (self%fran)
  end if
  call trace%end
END SUBROUTINE dealloc

!===============================================================================
!> Calculate a random force, confined to a shell in k-space, at regular intervals
!> in time.  We assume rho, P, and the sound speed to be near unity.  The velocity
!> amplitude should be of the order ampl_turb.  The size (scale) of the drivings
!> motions is of the order of 1.0/k1.  Thus the turnover time is given by
!> t_turb*ampl_turb = 1.0/k1; and the acceleration is of the order ampl_turb/t_turb.
!>
!> Note:  In order to restart properly, one needs to retrace the force updates.
!>
!>  06-aug-93/aake:  bug corrected; seed was not negative at start
!>  09-aug-93/aake:  changed to allow restart w same random coeffs
!>  21-feb-95/paolo: changed in order to use segments of constant
!>                   time der. of acceleration (so continuous accel.)
!>  23-jan-96/aake:  bug corrected; ran1() -> 2.*ran1()-1. to span [-1,1]
!>
!>            Restart must be done AFTER t=t_turn.
!>            ------------------------------------
!===============================================================================
FUNCTION turbulence (self, time, d, p, Ux, Uy, Uz, mesh) RESULT (ff)
  class(force_t)                             :: self
  real(8)                                    :: time
  class(mesh_t), pointer, dimension(:)       :: mesh
  real, dimension(mesh(1)%gn,mesh(2)%gn,mesh(3)%gn,3) :: ff
  real, dimension(:,:,:)                     :: Ux, Uy, Uz, d
  real, dimension(:,:,:,:)                   :: p
  !.............................................................................
  integer                :: i, j, k, jx, jy, jz, mx, my, mz, nrand, kmax
  complex                :: fxx, fyy, fzz, kx, ky, kz, corr
  complex                :: expikrx, expikry, expikrz
  real                   :: x, y, z, w, fact, fpow, accel
  real                   :: dx, dy, dz, eps, fk
  class(mesh_t), pointer :: m
  integer, save          :: itimer=0, nprint=5
  !-----------------------------------------------------------------------------
  if (.not.do_force) then
    ff(:,:,:,:) = 0.0
    return
  end if
  call trace_begin ('force_t%turbulence', itimer=itimer)
  eps = 1e-5
  mx = mesh(1)%gn
  my = mesh(2)%gn
  mz = mesh(3)%gn
  !-----------------------------------------------------------------------------
  !  Loop until self%t_turb > t
  !-----------------------------------------------------------------------------
  do while (self%t_turb <= time)
    self%t_turb = self%t_turb + t_turn
    !---------------------------------------------------------------------------
    ! Maintain two force snapshots
    !---------------------------------------------------------------------------
    self%fran(:,:,:,:,1) = self%fran(:,:,:,:,2)
    self%fran(:,:,:,:,2) = 0.0
    !---------------------------------------------------------------------------
    ! Count the number of wavenumbers
    !---------------------------------------------------------------------------
    nrand = 0
    kmax=ifix(k2+1.)
    do jx=-kmax,kmax
      do jy=-kmax,kmax
        do jz=-kmax,kmax
          fk = sqrt(float(jx**2+jy**2+jz**2))
          if (fk >= k1-eps .and. fk <= k2+eps) nrand=nrand+1
        enddo
      enddo
    enddo
    !-----------------------------------------------------------------------------
    ! Normalized acceleration factor
    !-----------------------------------------------------------------------------
    accel  = ampl_turb/t_turn/sqrt(float(nrand)/8.)    ! rms=1./8. per comp.
    !-----------------------------------------------------------------------------
    ! To obtain a Kolmogorov slope of the driving alone, one should have amplitudes
    ! a(k) that drop with k^(-11./6.), to have a(k)^2*k^2 = k^(-5./3.).  This ASSUMES
    ! that the amplitudes are proportional to the driving.  But the energy drain may
    ! be rather inversely proportional to the turnover time, which goes as k^(2./3.)
    !-----------------------------------------------------------------------------
    fpow = 0.
    do jx=-kmax,kmax
      kx = cmplx (0., jx*2.*pi/mesh(1)%b)
      do jy=-kmax,kmax
        ky = cmplx (0., jy*2.*pi/mesh(2)%b)
        do jz=-kmax,kmax
          kz = cmplx (0., jz*2.*pi/mesh(3)%b)
          fk = sqrt(float(jx**2+jy**2+jz**2))
          if (fk >= k1-eps .and. fk <= k2+eps) then
            fxx = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            fyy = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            fzz = cexp(cmplx(0., 2.*pi*self%random%ran1()))/fk**pk
            !------------------
            ! solenoidal field:
            !------------------
            if (do_helmh) then
               corr=(kx*fxx+ky*fyy+kz*fzz)/(kx*kx+ky*ky+kz*kz+1e-20) 
               if (jx /= 0) fxx = fxx - a_helmh*corr*kx
               if (jy /= 0) fyy = fyy - a_helmh*corr*ky
               if (jz /= 0) fzz = fzz - a_helmh*corr*kz
            endif
            !------------------
            fact=1.
            if (jx /= 0) fact=fact*0.5
            if (jy /= 0) fact=fact*0.5
            if (jz /= 0) fact=fact*0.5
            fpow = fpow+fact*(cabs(fxx)**2+cabs(fyy)**2+cabs(fzz)**2)
            dx = mesh(1)%d
            dy = mesh(2)%d
            dz = mesh(3)%d
            do k=1,mz
              m => mesh(3)
              z = m%p + m%r(k)
              do j=1,my
                m => mesh(2)
                y = m%p + m%r(j)
                m => mesh(1)
                do i=1,mx
                  x = m%p + m%r(i)
                  expikrx = accel*cexp(kx*(x-0.5*dx)+ky*y+kz*z)
                  expikry = accel*cexp(kx*x+ky*(y-0.5*dy)+kz*z)
                  expikrz = accel*cexp(kx*x+ky*y+kz*(z-0.5*dz))
                  self%fran(i,j,k,1,2) = self%fran(i,j,k,1,2) + real(fxx*expikrx)
                  self%fran(i,j,k,2,2) = self%fran(i,j,k,2,2) + real(fyy*expikry)
                  self%fran(i,j,k,3,2) = self%fran(i,j,k,3,2) + real(fzz*expikrz)
                end do
              end do
            end do
          endif
        enddo
      enddo
    enddo
    if (nprint>0 .and. self%id==1) then
      !$omp atomic
      nprint = nprint-1
      print '(a,i6,f10.4,i4,f10.4)', 'turbulence: id, t_turb, nrand, fpow =', &
        self%id, self%t_turb, nrand, fpow
    end if
  enddo
  !-----------------------------------------------------------------------------
  ! Time interpolation of the force
  !-----------------------------------------------------------------------------
  w = (t_turn-(self%t_turb-time))/t_turn
  w = 0.5*(1.-cos(w*pi))
  ff(:,:,:,:) = self%fran(:,:,:,:,1) + (self%fran(:,:,:,:,2)-self%fran(:,:,:,:,1))*w
  call trace_end (itimer)
END FUNCTION turbulence

END MODULE force_mod
