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
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t)  :: initial
    type(data_hub_t) :: data_hub
    type(poisson_t)  :: poisson
    logical:: selfgravity=.false.
    real:: d0=1.0
  contains
    procedure:: init
    procedure:: dnload
    procedure:: gravity
  end type
  real(8), save :: restart_time = 0d0
  logical, save :: selfgravity=.false., prediction=.true.
  logical, save :: detailed_timer=.false.
  integer, save :: max_iter=100, verbose=0
  real,    save :: epsilon=1e-3, fourpiG=1.0
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
  namelist /selfgravity_params/ selfgravity, prediction, fourpiG, epsilon, &
                                max_iter, verbose
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
!> USE a Poisson solver module to solve for the potential phi in the current
!> patch.  The task of the Poisson solver is to make one correction of phi,
!> which is alredy from the start based on the time extrapolation of the
!> previous solution.
!> =============================================================================
SUBROUTINE gravity (self, error, only)
  class(gpatch_t):: self
  integer, optional:: only !unused
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
  do iz=self%mesh(3)%li,self%mesh(3)%ui
    do iy=self%mesh(2)%li,self%mesh(2)%ui
      do ix=self%mesh(1)%li,self%mesh(1)%ui
        source(ix,iy,iz) = fourpiG*(d(ix,iy,iz)-self%d0)
      end do
    end do
  end do

  ! call the Poisson solver iteratively
  sec = wallclock()
  do iter=1,max_iter
    call self%poisson%update (self, phi, resid, source, d, error, detailed_timer)
    !---------------------------------------------------------------------------
    ! The error coming back from poisson should be (Laplace(phi)-source)/d, which
    ! we then here normalize with grho, to get the relative error in density.
    !---------------------------------------------------------------------------
    error = error / fourpiG
    if (iter == 1) error0 = error
    if (verbose > 1) print "(i6,' self_gravity: iter =',i3.1,' error =',1p,g14.5)", self%id, iter, error

    ! check for convergence
    if (error < epsilon * error0) then
      if (verbose>0) &
        print "(i6,' poisson converged in',i4,' iterations; initial, final error, time taken:',1p,3e10.3)", &
          self%id, iter, error0, error, wallclock()-sec
      exit
    end if
  end do

  ! too many iterations; throw warning before continuing
  if (iter > max_iter) then
    print "(i6,' WARNING: poisson solution not converged in', i4, ' iterations')", &
      self%id, max_iter
    print "('error0 = ',g14.5,'; error = ',g14.5)",error0,error
  end if

  deallocate (source, resid)
  ! compute estimate of dphi/dt
  call time_derivs (self)

  call trace%end (itimer)
END SUBROUTINE gravity

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
