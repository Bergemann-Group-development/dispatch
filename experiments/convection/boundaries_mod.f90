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
  real,save:: bc_rate=0.5
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t):: self
  integer:: iostat
  logical:: flux_BC = .true.
  real:: theta = 10.0
  real:: mpoly = 1.0
  namelist /boundaries_params/ theta, mpoly, flux_BC, bc_rate
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
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: d, e, s, eos_d, eos_e, eos_s
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: pg, tt
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: px, py, pz
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  integer:: ix, iy, iz, iz1, n(3), id, lb, ub
  real(8), pointer:: z(:)
  real:: h, f, dz, ds(3), d0, pg0, tt0, d1, pg1, tt1, dtt, dss
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
  if (patch%idx%e>0) allocate (e (n(1),n(2),n(3)))
  !-----------------------------------------------------------------------------
  ! Point pointers
  !-----------------------------------------------------------------------------
  d  => mem(:,:,:,patch%idx%d)
  px => mem(:,:,:,patch%idx%px)
  py => mem(:,:,:,patch%idx%py)
  pz => mem(:,:,:,patch%idx%pz)
  if (patch%idx%s>0) s=> mem(:,:,:,patch%idx%s)
  !-----------------------------------------------------------------------------
  ! Get velocities
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = px(ix,iy,iz)/d(ix,iy,iz)
        uy(ix,iy,iz) = py(ix,iy,iz)/d(ix,iy,iz)
        uz(ix,iy,iz) = pz(ix,iy,iz)/d(ix,iy,iz)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Get thermal energy, if applicable
  !-----------------------------------------------------------------------------
  if (patch%idx%e>0) then
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e (ix,iy,iz) = mem(ix,iy,iz,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                          ux(ix,iy,iz)**2 &
                        + uy(ix,iy,iz)**2 &
                        + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! The lb and ub indices are the index limits for BC modified values.  The BCs
  ! are placed at cell centers -- using face centered BCs does not work well
  !-----------------------------------------------------------------------------
  lb = m(3)%li
  if (m(3)%lower_boundary) then
    !---------------------------------------------------------------------
    ! Compute the scale height and fill in the guard zone with values of
    ! pressure and density consistent with constant temperature HSE
    !---------------------------------------------------------------------
    tt0 = self%tt_top
    h = tt0/self%g
    do iy=m(2)%lb,m(2)%ub
      do iz=m(3)%lb,lb-1
        iz1 = lb + (lb - iz)
        do ix=m(1)%lb,m(1)%ub
          !---------------------------------------------------------------------
          ! Compute the scale height and fill in the guard zone with values of
          ! pressure and density consistent with constant temperature HSE
          !---------------------------------------------------------------------
          dz = (iz-m(3)%li)*ds(3)
          f = exp(dz/h)
          d(ix,iy,iz)  = d(ix,iy,m(3)%li)*f
          if (patch%idx%e>0) then
            ! Constant temperature BC
            !e(ix,iy,iz)  = d(ix,iy,iz)*tt0/(eos%gamma-1d0)
            tt1 =  self%tt_top+max(0.0,self%theta*(z(iz)+m(3)%p))
            e(ix,iy,iz) = tt1*d(ix,iy,iz)/(eos%gamma-1d0)
          end if
          if (patch%idx%s>0) then
            tt1 =  self%tt_top+max(0.0,self%theta*(z(iz)+m(3)%p))
            pg1 = tt1*d(ix,iy,iz)
            s(ix,iy,iz) = log(pg1*d(ix,iy,iz)**-eos%gamma) &
                        * d(ix,iy,iz) / (eos%gamma-1d0)
            ! below - constant entropy gradient
            ! dss = (s(ix,iy,lb+1)-s(ix,iy,lb))/ds(3)
            ! s(ix,iy,iz) = s(ix,iy,lb)+dss*dz
          end if
          ux(ix,iy,iz) = +ux(ix,iy,iz1)
          uy(ix,iy,iz) = +uy(ix,iy,iz1)
          uz(ix,iy,iz) = -uz(ix,iy,iz1)
          px(ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
          py(ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
          pz(ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
        end do
      end do
      !-------------------------------------------------------------------------
      ! Impose constant temperature and vanishing mass flux at the boundary
      !-------------------------------------------------------------------------
      !e (:,iy,lb) = d(:,iy,lb)*self%tt_top/(eos%gamma-1d0)
      if (patch%idx%e>0) &
        e (:,iy,lb) = (1.+self%theta*(z(iz)+m(3)%p))*d(:,iy,iz)/(eos%gamma-1d0)
      pz(:,iy,lb) = 0.0
      uz(:,iy,lb) = 0.0
    end do
  end if
  !-----------------------------------------------------------------------------
  ! For pressure and temperature at the bottom boundary:  Leave the initial
  ! values in the T-case, and impose the given gradient in the F-case
  !-----------------------------------------------------------------------------
  ub = m(3)%ui
  if (m(3)%upper_boundary) then
    if (self%flux_BC) then
      ! ------------------------------------------------------------------------
      ! Look up pressure and temperature at ub-1
      ! ------------------------------------------------------------------------
      allocate(pg(n(1),n(2),1))
      allocate(tt(n(1),n(2),1))
      eos_d => d(:,:,ub-1:ub-1)
      if (patch%idx%e>0) then
        eos_e => e(:,:,ub-1:ub-1)
        call eos%lookup_table(shape(pg),d=eos_d,e=eos_e,pg=pg,tt=tt)
      end if
      if (patch%idx%s>0) then
        eos_s => s(:,:,ub-1:ub-1)
        call eos%lookup_table(shape(pg),d=eos_d,s=eos_s,pg=pg,tt=tt)
      end if
      do iz=ub+1,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            !-------------------------------------------------------------------
            ! Compute the pressure scale height from the temperature, and
            ! extrapolate the pressure in the guard zones from that, then get
            ! the density from the ratio of gas pressure to temperature.
            !-------------------------------------------------------------------

            !!d0  = d(ix,iy,ub-1)
            !!pg0 = e(ix,iy,ub-1)*(eos%gamma-1d0)
            !!tt0 = pg0/d0                        ! T(ui)
            dz = (iz-(ub-1))*ds(3)              ! depth interval
            !tt1 = e(ix,iy,ub-2)*(eos%gamma-1d0)/d(ix,iy,ub-1)
            !dtt=(tt0-tt1)/ds(3)
            !tt1 = tt0+dtt*dz
            tt1 = tt(ix,iy,1) + self%theta*dz           ! T(iz)
            h = 0.5*(tt(ix,iy,1)+tt1)/self%g            ! average scale height
            pg1 = pg(ix,iy,1)*exp(dz/h)                 ! Pg(iz)
            d(ix,iy,iz) = pg1/tt1               ! d(iz)
            if (patch%idx%e>0) e(ix,iy,iz) = pg1/(eos%gamma-1d0)   ! e(iz)
            if (patch%idx%s>0) then
              s(ix,iy,iz) = log(pg1*d(ix,iy,iz)**-eos%gamma) &
                        * d(ix,iy,iz) / (eos%gamma-1d0)
            ! below - constant entropy gradient
            ! dss = (s(ix,iy,ub)-s(ix,iy,ub-1))/ds(3)
            ! s(ix,iy,iz) = s(ix,iy,lb)+dss*dz
            end if
          end do
        end do
      end do
      deallocate(pg,tt)
    end if
    !---------------------------------------------------------------------------
    ! Antisymmetric vertical and symmetric horizontal velocity componentsd
    !---------------------------------------------------------------------------
    do iy=m(2)%lb,m(2)%ub
      do iz=ub+1,m(3)%ub
        iz1 = ub - (iz-ub)
        do ix=m(1)%lb,m(1)%ub
          ux(ix,iy,iz) = +ux(ix,iy,iz1)
          uy(ix,iy,iz) = +uy(ix,iy,iz1)
          uz(ix,iy,iz) = -uz(ix,iy,iz1)
          px(ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
          py(ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
          pz(ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
        end do
      end do
      pz(:,iy,ub) = 0.0
      uz(:,iy,ub) = 0.0
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Put the total energy back into mem, if applicable
  !-----------------------------------------------------------------------------
  if (patch%idx%e>0) then
    if (m(3)%lower_boundary) then
      do iz=m(3)%lb,lb
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            mem(ix,iy,iz,patch%idx%e) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                        ux(ix,iy,iz)**2 &
                                      + uy(ix,iy,iz)**2 &
                                      + uz(ix,iy,iz)**2)
          end do
        end do
      end do
    end if
    if (m(3)%upper_boundary) then
      if (self%flux_BC) then
        do iz=ub,m(3)%ub
          do iy=m(2)%lb,m(2)%ub
            do ix=m(1)%lb,m(1)%ub
              mem(ix,iy,iz,patch%idx%e) = e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                        ux(ix,iy,iz)**2 &
                                      + uy(ix,iy,iz)**2 &
                                      + uz(ix,iy,iz)**2)
            end do
          end do
        end do
      end if
    end if
    deallocate(e)
  end if
  deallocate (ux,uy,uz)
  call trace%end (itimer)
END SUBROUTINE primary

END MODULE boundaries_mod
