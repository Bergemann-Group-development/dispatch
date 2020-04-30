!===============================================================================
!===============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE mpi_mod
  USE trace_mod
  USE kinds_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE vector_ops
  USE stagger_mod
  USE gpatch_mod
  USE force_mod
  USE index_mod
  USE timestep_mod
  USE non_ideal_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, e, dedt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: nu(6), csound
    logical:: mhd=.true.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: pde
    procedure:: update
    procedure:: output
    procedure:: gas_pressure
  end type
  real, save:: smag_s=0.2, smag_d=0.1
CONTAINS

!===============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!===============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  integer:: i, iostat
  logical, save:: hardwire=.false.
  logical, save:: first_time=.true.
  character(len=16), save:: eos
  namelist /mhd_params/ nu, gamma, csound, cdtd, hardwire, eos
  namelist /smag_params/ smag_s, smag_d
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  call self%gpatch_t%init
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, mhd_params)
    rewind (io%input); read (io%input, smag_params, iostat=iostat)
    if (io%master) write (*, mhd_params)
    if (io%master) write (*, smag_params)
    call non_ideal%init
    call timestep%init
    ! fool-proofing
    if (timestep%time_order <= 0) call mpi%abort("The Stagger solvers require `time_order` > 0. Abort!")
  end if
  call io%lock%unset
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    self%mesh(i)%h(5+i) = -0.5
  end do
  call self%force%init (self%kind)
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, emf, J, fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: du, lnd, ee, cs, pg, pa, pp, pb, Q
  real, dimension(:,:,:),   pointer:: divu, s2, sxx, syy, szz, sxy, sxz, syz
  !.............................................................................
  real(8):: ds(3), dx
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('mhd_t%pde',itimer=itimer)
  call self%output
  associate (i => self%idx)
  !-----------------------------------------------------------------------------
  ! Pointer assignment and scratch allocation
  !-----------------------------------------------------------------------------
  self%d => self%mem(:,:,:,i%d,      self%it,1); self%dddt => self%mem(:,:,:,i%d,      self%it,2)
  self%e => self%mem(:,:,:,i%e,      self%it,1); self%dedt => self%mem(:,:,:,i%e,      self%it,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%it,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%it,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
  end if
  end associate
  call allocate_vectors (self%gn, ld, dd, emf, J, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pa, du, ee, cs, pp, pb, Q)
  call allocate_scalars (self%gn, divu, s2, sxx, syy, szz, sxy, sxz, syz)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,i%d,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
associate (d=>self%d, dddt=>self%dddt, &
           e=>self%e, dedt=>self%dedt, &
           p=>self%p, dpdt=>self%dpdt, &
           B=>self%B, dBdt=>self%dBdt, &
           px=>self%p(:,:,:,1), py=>self%p(:,:,:,2), pz=>self%p(:,:,:,3), &
           Bx=>self%B(:,:,:,1), By=>self%B(:,:,:,2), Bz=>self%B(:,:,:,3), &
           dpdtx=>self%dpdt(:,:,:,1), dpdty=>self%dpdt(:,:,:,2), dpdtz=>self%dpdt(:,:,:,3), &
           dBdtx=>self%dBdt(:,:,:,1), dBdty=>self%dBdt(:,:,:,2), dBdtz=>self%dBdt(:,:,:,3), &
           Jx =>  J(:,:,:,1), Jy =>  J(:,:,:,2), Jz =>  J(:,:,:,3), &
           Ux =>  U(:,:,:,1), Uy =>  U(:,:,:,2), Uz =>  U(:,:,:,3), &
           fxx=> fx(:,:,:,1), fxy=> fx(:,:,:,2), fxz=> fx(:,:,:,3), &
           fyx=> fy(:,:,:,1), fyy=> fy(:,:,:,2), fyz=> fy(:,:,:,3), &
           fzx=> fz(:,:,:,1), fzy=> fz(:,:,:,2), fzz=> fz(:,:,:,3), &
           fdx=> fd(:,:,:,1), fdy=> fd(:,:,:,2), fdz=> fd(:,:,:,3), &
           ldx=> ld(:,:,:,1), ldy=> ld(:,:,:,2), ldz=> ld(:,:,:,3), &
           ddx=> dd(:,:,:,1), ddy=> dd(:,:,:,2), ddz=> dd(:,:,:,3), &
           Ex =>emf(:,:,:,1), Ey =>emf(:,:,:,2), Ez =>emf(:,:,:,3))
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = down(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure
  !-----------------------------------------------------------------------------
  pg = self%gas_pressure()
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure
  !-----------------------------------------------------------------------------
  ds = self%mesh%d
  dx = minval(ds,self%n > 1)                            ! smallest {dx,dy,dz}
  sxx = ddxup(ds,Ux)
  syy = ddyup(ds,Uy)
  szz = ddzup(ds,Uz)
  divu = sxx+syy+szz
  if (self%gamma /= 1.0) then                           ! unless isothermal
    dedt = -pg*divu                                     ! -Pg*div(U) term
  else
    dedt = 0.0
  end if
  du = -dx*divu                                         ! convergence rate
  du = max(du,0.0)                                      ! positive part
  self%u_max = 0.0
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  !-----------------------------------------------------------------------------
  ! Diagonal Reynolds stress and artificial pressure
  fxx = xup(Ux)**2
  fyy = yup(Uy)**2
  fzz = zup(Uz)**2
  pp = (fxx+fyy+fzz)                                    ! velocity squared
  pa = self%nu(2)*d*du**2                               ! dimension: pressure
  if (self%mhd) then
    ! Diagonal (negative) Maxwell stress (borrowing emf for scratch)
    Ex = xup(Bx)**2
    Ey = yup(By)**2
    Ez = zup(Bz)**2
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)
    ! Fast mode speed + velocity
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>2) print*,'u_max',self%u_max
    ! Isotropic pressure part
    pp = pg + pa + pb
    ! Total diagonal stress
    du = dx*self%nu(1)*cs
    fxx = pp + d*(fxx - du*sxx) - Ex
    fyy = pp + d*(fyy - du*syy) - Ey
    fzz = pp + d*(fzz - du*szz) - Ez
    Q = d*du*(sxx**2+syy**2+szz**2)
  else
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>2) print*,'u_max',self%u_max
    pp = pg + pa
    ! Total diagondal stress, with density factor
    du = dx*self%nu(1)*cs
    fxx = pp + d*(fxx - du*sxx)
    fyy = pp + d*(fyy - du*syy)
    fzz = pp + d*(fzz - du*szz)
    Q = d*du*(ddxup(ds,Ux)**2+ddyup(ds,Uy)**2+ddzup(ds,Uz))
  end if
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -div(ds,fm)                                    ! density derivative
  self%u_max = max(self%u_max,-0.4*dx*self%fminval(dddt/d))
  if (io%verbose>2) print*,'u_max',self%u_max
  !-----------------------------------------------------------------------------
  ! Conservation of energy
  !-----------------------------------------------------------------------------
  if (self%gamma /= 1.0) then
    ee = e/d
    fm = fm*down(ee) - (self%nu(5)*dd*down(du))*ddown(ds,ee)
    dedt = dedt - div(ds,fm) + Q
  end if
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)
  fxz = zdn(Ux)*xdn(Uz)
  fyz = zdn(Uy)*ydn(Uz)
  !---------------------------------------------------------------------------
  ! Add viscous part -- Smagorinsky style
  !---------------------------------------------------------------------------
  sxy = 0.5*(ddxdn(ds,Uy)+ddydn(ds,Ux))
  sxz = 0.5*(ddxdn(ds,Uz)+ddzdn(ds,Ux))
  syz = 0.5*(ddydn(ds,Uz)+ddzdn(ds,Uy))
  s2 = 2.0*(xup(yup(sxy))**2 + xup(zup(sxz))**2 + yup(zup(syz))**2) &
     + ddxup(ds,Ux)**2 + ddyup(ds,Uy)**2 + ddzup(ds,Uz)**2
  s2 = dx**2*(smag_s*sqrt(2.0*s2) + smag_d*abs(divu))
  fxy = fxy - xdn1(ydn1(s2))*sxy
  fxz = fxz - ydn1(zdn1(s2))*sxz
  fyz = fyz - zdn1(xdn1(s2))*syz
  fxx = fxx - d*s2*sxx
  fyy = fyy - d*s2*syy
  fzz = fzz - d*s2*szz
  !---------------------------------------------------------------------------
  ! Multiply non-diagnoal part with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy
  fxz = exp(zdn(ldx))*fxz
  fyz = exp(zdn(ldy))*fyz
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - ydn(Bx)*xdn(By)
    fxz = fxz - zdn(Bx)*xdn(Bz)
    fyz = fyz - zdn(By)*ydn(Bz)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)
    fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)
    fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)
    fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
    fxz = fxz + (xdn1(fdz)*zdn(Ux)+zdn1(fdx)*xdn(Uz))
    fyz = fyz + (zdn1(fdy)*ydn(Uz)+ydn1(fdz)*zdn(Uy))
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations
  !-----------------------------------------------------------------------------
  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fxz)
  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)
  dpdtz = - ddxup1(ds,fxz) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    emf = self%nu(6)*edge(du)*J
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))
    if (self%gamma /= 1.0) &
      dedt = dedt + Q
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)
    call non_ideal%update (self%gn, minval(ds), self%d, Jx, Jy, Jz, Bx, By, Bz, Q, emf, self%u_max)
    dBdt = -curl_up(ds,emf)
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, emf, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ee, cs, pp, pb, Q)
  call deallocate_scalars (divu, s2, sxx, syy, szz, sxy, sxz, syz)
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace_end (itimer)
end associate
END SUBROUTINE pde

!===============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod.
!===============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  call self%pde
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
END SUBROUTINE update

!===============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!===============================================================================
SUBROUTINE output (self, experiment_name)
  class(mhd_t):: self
  character(len=64):: experiment_name, filename
  optional:: experiment_name
  logical:: exists
  !-----------------------------------------------------------------------------
  if (self%time < self%out_next .or. .not.io%do_output .or. .not.io%do_legacy) return
  call self%gpatch_t%output (experiment_name)
  !$omp critical (output_cr)
  write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
    self%id, self%iout-1
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (output_cr)
END SUBROUTINE output

!===============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%d*self%csound**2
  else
    pg = self%e*(self%gamma-1.0)
  end if
END FUNCTION gas_pressure

FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(mhd_t) :: self
  integer :: c(3)
  real :: pg(c(1),c(2),c(3))
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3), 1,self%it,1), &
             e => self%mem(c(1),c(2),c(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d*self%csound**2
  else
    pg = self%e*(self%gamma-1.0)
  end if
  end associate
END FUNCTION cell_gas_pressure

END MODULE mhd_mod
