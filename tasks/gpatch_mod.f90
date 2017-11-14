!> =============================================================================
!> The gpatch_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!> =============================================================================
MODULE gpatch_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE initial_mod
  USE kinds_mod
  USE data_hub_mod
  USE mesh_mod
  USE units_mod
  USE scaling_mod
  USE poisson_mod
  USE omp_timer_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t)  :: initial
    type(data_hub_t) :: data_hub
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real(8):: d0=0.0
  contains
    procedure:: init
    procedure:: dnload
    procedure:: gravity
    procedure:: update
  end type
  real(8), save :: total_mass=0d0, total_volume=0d0
  real(8), save :: restart_time = 0d0
  logical, save :: selfgravity=.false., prediction=.true., detailed_timer=.false.
  integer, save :: max_iter=100, verbose=0
  real,    save :: repsilon=1e-2, aepsilon=1e-3, fourpiG=1.0
  integer, save :: order=1
CONTAINS

!> =============================================================================
!> Add restart functionality to patch initialization
!> =============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  integer                           :: iostat
  logical, save                     :: first_time=.true.
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  namelist /selfgravity_params/ selfgravity, prediction, fourpiG, repsilon, aepsilon, &
                                max_iter, order, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%init')
  !$omp critical (poisson_cr)
  if (first_time) then
    fourpiG = scaling%grav
    rewind (io%input)
    read (io%input, selfgravity_params, iostat=iostat)
    if (io%master) write (io%output, selfgravity_params)
  end if
  !$omp end critical (poisson_cr)
  if (selfgravity) then
    self%selfgravity = selfgravity
    if (trim(self%kind) /= 'zeus_mhd_patch') self%nv = self%nv + 1
    self%idx%phi = self%nv
    if (self%nw == 1) then
      self%nv = self%nv + 1
      self%idx%dphi = self%nv
    end if
    call self%patch_t%init
    call self%poisson%init (self, selfgravity)
  else
    call self%patch_t%init
  end if
  call self%initial%init (self%kind)
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    if (self%use_data_hub) then
      !call data_hub%test
    else
      call download%test
    end if
  end if
  !$omp end critical (poisson_cr)
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
  !-----------------------------------------------------------------------------
  ! If the restart parameter is not set, or the file does not exist, initialize
  ! the patch from scratch, resetting other necessary counters only if positive
  ! restart_time indicates that restart files were read for previous patches.
  !-----------------------------------------------------------------------------
  else
    self%mem = 0.0
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff, self%idx)
  end if
  !-----------------------------------------------------------------------------
  ! Sum up the total mass and total volume, to get the average density d0.
  ! FIXME: This needs to be modified for MPI
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  !$omp atomic
  total_mass = total_mass + self%fsum(d)*product(self%ds)
  !$omp atomic
  total_volume = total_volume + product(self%size)
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = (int(self%time/io%out_time)+1)*io%out_time
    self%iout = self%restart+1
    !$omp master
    write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
    if (io%iodir/=self%iout) call io%mkdir (trim(filename))
    io%iodir = self%iout
    !$omp end master
  end if
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(gpatch_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  if (self%use_data_hub) then
    call self%data_hub%update (self%link, only=only)
  else
    call download%download_link (self%link, only=only)
  end if
END SUBROUTINE dnload

!> =============================================================================
!> Prepare the external forces for the pde call, in particular selfgravity. At
!> this point the nbors are declared up-to-date with respect to everything,
!> including their potentials, which we need for guard zones.  They may in fact
!> be providing only their preliminary potentials in their (nt-1) time slots,
!> while their older time slots should by now contain the converged potential.
!> =============================================================================
SUBROUTINE gravity (self)
  class(gpatch_t):: self
  !.............................................................................
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, phi, source, resid
  real:: error, error0
  real(8):: sec
  integer:: iter, ix ,iy, iz
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('gpatch_t%gravity', itimer=itimer)
  d   => self%mem(:,:,:,self%idx%d  ,self%it,1)
  phi => self%mem(:,:,:,self%idx%phi,self%it,1)
  allocate (source(size(d,1),size(d,2),size(d,3)), resid(size(d,1),size(d,2),size(d,3)))
  source = 0.0
  if (self%d0==0.0) then
    self%d0 = total_mass/total_volume
    print *,'setting self%d0 =', self%d0, self%faver(d)
  end if
  do iz=self%mesh(3)%li,self%mesh(3)%ui
    do iy=self%mesh(2)%li,self%mesh(2)%ui
      do ix=self%mesh(1)%li,self%mesh(1)%ui
        source(ix,iy,iz) = fourpiG*(d(ix,iy,iz)-self%d0)
      end do
    end do
  end do
  if (self%is_periodic()) call self%make_periodic (source)
  sec = wallclock()
  call self%poisson%laplace (self%ds, phi, source, resid)
!  print *,'src1', self%fminval(source), self%fmaxval(source)
!  print *,'phi1', self%fminval(phi), self%fmaxval(phi)
!  print *,'res1', self%fminval(resid), self%fmaxval(resid)
!  write (10) shape(phi); write (10) source, phi, resid
  if (self%poisson%kind == 'zeus_sor') self%poisson%omega_sor = 1.0 / ( 1.0 - 0.5 * self%poisson%rjac2 )
  if (self%poisson%kind=='cg') then
    call self%poisson%cg (self, phi, source, error0, error, iter, detailed_timer)
  else if (self%poisson%kind=='sor_iter') then
    call self%poisson%sor_iter (self, phi, source, error0, error, iter, detailed_timer)
  else
    !---------------------------------------------------------------------------
    ! call the Poisson solver iteratively
    !---------------------------------------------------------------------------
    if (self%poisson%kind == 'zeus_sor') self%poisson%omega_sor = 1.0 / ( 1.0 - 0.5 * self%poisson%rjac2 )
    do iter=1,max_iter
    call self%poisson%update (self, phi, resid, source, d, error, detailed_timer)
    if (self%is_set(bits%root_grid)) call self%periodic_grid (only=self%idx%phi)
      !---------------------------------------------------------------------------
      ! The error coming back from poisson should be (Laplace(phi)-source)/d, which
      ! we then here normalize with grho, to get the relative error in density.
      !---------------------------------------------------------------------------
    error = error / fourpiG
      if (iter == 1) error0 = error
      if (verbose > 1) print "(i6,2x,a,': iter = ',i4.1,3x,'error =',1p,g12.4)", &
        self%id, trim(self%poisson%kind), iter, error
      !---------------------------------------------------------------------------
      ! check for convergence
      !---------------------------------------------------------------------------
      if (error < repsilon * error0 .or. error < aepsilon) then
        exit
      end if
    end do
    !---------------------------------------------------------------------------
    ! too many iterations; throw warning before continuing
    !---------------------------------------------------------------------------
    if (iter > max_iter) then
      print "(i6,' WARNING: poisson solution not converged in', i4, ' iterations')", &
        self%id, max_iter
      print "('error0 = ',g14.5,'; error = ',g14.5)",error0,error
    end if
  end if
  call self%poisson%laplace (self%ds, phi, source, resid)
!  print *,'src2', self%fminval(source), self%fmaxval(source)
!  print *,'phi2', self%fminval(phi), self%fmaxval(phi)
!  print *,'res2', self%fminval(resid), self%fmaxval(resid)
!  write (10) shape(phi); write (10) source, phi, resid

  ! too many iterations; throw warning before continuing
  if (iter > max_iter) then
    print "(i6,' WARNING: poisson solution not converged in', i4, ' iterations')", &
      self%id, max_iter
    print "('error0 = ',g14.5,'; error = ',g14.5)",error0,error
  end if
  if (verbose>0) &
    print "(i6,2x,a,' converged in ',i4,' iterations; initial, final error, time taken:',1p,3e10.3,2x,2e12.3)", &
      self%id, trim(self%poisson%kind), iter, error0, error, wallclock()-sec, &
      self%fminval(phi), self%fmaxval(phi)
  deallocate (source, resid)
  !-----------------------------------------------------------------------------
  ! compute estimate of dphi/dt
  !-----------------------------------------------------------------------------
  call time_derivs (self)
  call trace%end (itimer)
END SUBROUTINE gravity

!> =============================================================================
!> Having advanced the MHD state to the new time, we make a preliminary
!> call to the Poisson solver, based on extrapolation of the nbor potentials.
!> We then provide this preliminary potential to nbors, until we have a
!> better offer.
!> =============================================================================
SUBROUTINE update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  if (.not. selfgravity) return
  call trace%begin ('gpatch_t%update')
  call self%rotate
  call self%dnload (only=self%idx%d)
  call self%dnload (only=self%idx%phi)
  call self%gravity
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Compute an approximate, conservative estimate of dphi/dt, which will cause
!> the timestep routine to make a good guess at the next phi, at the same time
!> as it is updating the other variables.  The 1st order prediction is in
!> practice better at reducing the number of iterations than the 2nd order.
!> =============================================================================
SUBROUTINE time_derivs (self)
  class(gpatch_t):: self
  integer:: it1, it2, it3, nt, l(3), u(3)
  real(8):: t1, t2, t3
  integer, save:: itimer=0
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: dphidt
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('gpatch_t%time_derivs', itimer=itimer)
  if (self%nw == 2) then
    dphidt => self%mem(:,:,:,self%idx%phi,self%it,2)
  else
    dphidt => self%mem(:,:,:,self%idx%dphi,self%it,1)
  end if

  if (prediction) then
    nt = self%nt
    it1 = self%iit(nt-1)
    it2 = self%iit(nt-2)
    it3 = self%iit(nt-3)
    t1 = self%t(it1)
    t2 = self%t(it2)
    t3 = self%t(it3)
    dphidt = 0.0
    if (it2/=it1 .and. t1/=t2) then
      l = self%mesh%li
      u = self%mesh%ui
      if (order==2 .and. t2/=t3) then
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)*1.5 - &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it3,1))/(t2-t3)*0.5
      else
        dphidt(l(1):u(1),l(2):u(2),l(3):u(3)) = &
                 (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it1,1) - &
                  self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%phi,it2,1))/(t1-t2)
      end if
      if (verbose>1) print *,self%id,'d(phi)/dt min,max,nv =', &
          self%fminval(dphidt), self%fmaxval(dphidt), size(self%mem,4)
    end if
  else
    dphidt = 0.0
  end if
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE time_derivs

END MODULE gpatch_mod
