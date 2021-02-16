!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the C-shock problem to test ambipolar diffusion
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
!  USE scalar_mod
  USE boundary_mod
!  USE timestep_mod
  USE index_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Setup a patch with a sharp shock front as present in a C-shock
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  real:: x, y, z, g1
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: b, p
  real, dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, px, py, pz, s, e1, et
  real,save:: P_l, P_r, rho_i, rho_l, vx_l, vy_l, vz_l=0., &
              Bx_l, By_l, Bz_l=0., rho_r, vx_r, vy_r, vz_r=0., Bx_r, By_r, Bz_r=0., &
              x0=0.0
  namelist /experiment_params/ P_l, P_r, rho_l, vx_l, vy_l, Bx_l, By_l, rho_r, &
                               vx_r, vy_r, Bx_r, By_r, x0

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')

  self%nv = 8
  self%mhd = .true.
  call self%solver_t%init

  call io%lock%set ('experiment_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  call io%lock%unset ('experiment_t%init')

  if (self%kind(1:6) == 'ramses') then 
    et => self%mem(:,:,:,self%idx%e,self%it,1)
    allocate (e1(self%gn(1),self%gn(2),self%gn(3)))
  endif

  d => self%mem(:,:,:,self%idx%d,self%it,1)
  if (self%kind(1:7) == 'stagger') s => self%mem(:,:,:,self%idx%s,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  B => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);

  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  print*, 'gn is:', self%gn(1)
  print*, 'P and rho and gamma:', P_l, P_r, rho_l, rho_r, eos%gamma

  if (eos%gamma == 1d0) then
    g1 = 1.0
  else
    g1 = eos%gamma-1d0
  end if
  do i=self%mesh(1)%lb,self%mesh(1)%ub
    x = self%mesh(1)%p + self%mesh(1)%r(i)
    !y = self%mesh(2)%p + self%mesh(2)%r(i)
    !z = self%mesh(3)%p + self%mesh(3)%r(i)
    if (x < x0) then
      d (i,:,:)  = rho_l
      px(i,:,:) = vx_l*rho_l
      py(i,:,:) = vy_l*rho_l
      pz(i,:,:) = 0d0*rho_l
      Bx(i,:,:) = Bx_l
      By(i,:,:) = By_l
      Bz(i,:,:) = Bz_l
      if (self%kind(1:6) == 'ramses') then
          et(i,:,:)= P_l/g1 + 0.5*rho_l * (vx_l**2+vy_l**2+vz_l**2) + 0.5 * (Bx(i,:,:)**2 + By(i,:,:)**2 + Bz(i,:,:)**2)  
      end if
      if (self%kind(1:7) == 'stagger') then
        s (i,:,:) = (log(P_l)-log(rho_l)*eos%gamma)/g1 * rho_l
      end if
    else
      d (i,:,:) = rho_r
      px(i,:,:) = vx_r*rho_r
      py(i,:,:) = vy_r*rho_r
      pz(i,:,:) = 0d0*rho_r
      Bx(i,:,:) = Bx_r 
      By(i,:,:) = By_r
      Bz(i,:,:) = Bz_r 
      if (self%kind(1:6) == 'ramses') then
        et(i,:,:)= P_r/g1 + 0.5*rho_r * (vx_r**2+vy_r**2+vz_r**2) + 0.5 * (Bx(i,:,:)**2 + By(i,:,:)**2 + Bz(i,:,:)**2)            
      end if
      if (self%kind(1:7) == 'stagger') then
        s (i,:,:) = (log(P_r)-log(rho_r)*eos%gamma)/g1 * rho_r
      end if
    endif
  end do

  self%periodic = .true.
  self%mhd = .true.
  print*, 'done!'

  ! assign constant boundary conditions.
  boundary%P_l   = P_l
  boundary%rho_l = rho_l
  if (self%kind(1:7) == 'stagger') then
    if (eos%gamma .eq. 1.) then
      boundary%s_l = -log(rho_l)*rho_l
    else
      boundary%s_l = (log(P_l)-log(rho_l)*eos%gamma)/(eos%gamma-1.0) * rho_l
    end if
  end if
  if (self%kind(1:6) == 'ramses') then
    if (eos%gamma .eq. 1.) then
      boundary%e_l = P_l
    else
      boundary%e_l = P_l/(eos%gamma-1d0) + 0.5*rho_l * (vx_l**2+vy_l**2+vz_l**2) + 0.5 * (Bx_l**2 + By_l**2 + Bz_l**2)
    end if
  end if
  boundary%vx_l  = vx_l
  boundary%vy_l  = vy_l
  boundary%vz_l  = vz_l
  boundary%Bx_l  = Bx_l
  boundary%By_l  = By_l
  boundary%Bz_l  = Bz_l
  boundary%P_r   = P_r
  boundary%rho_r = rho_r
  if (self%kind(1:7) == 'stagger') then
    if (eos%gamma .eq. 1.) then
      boundary%s_r = -log(rho_r)*rho_r
    else
      boundary%s_r = (log(P_r)-log(rho_r)*eos%gamma)/(eos%gamma-1.0) * rho_r
    end if
  end if
  if (self%kind(1:6) == 'ramses') then
    if (eos%gamma .eq. 1.) then
      boundary%e_r = P_r
    else
      boundary%e_r = P_r/(eos%gamma-1d0) + 0.5*rho_r * (vx_r**2+vy_r**2+vz_r**2) + 0.5 * (Bx_r**2 + By_r**2 + Bz_r**2)
    end if
  end if
  boundary%vx_r  = vx_r
  boundary%vy_r  = vy_r
  boundary%vz_r  = vz_r
  boundary%Bx_r  = Bx_r
  boundary%By_r  = By_r
  boundary%Bz_r  = Bz_r

  call trace_end
END SUBROUTINE init

!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  if (self%llc_nat(1)-self%ds(1) .lt. self%mesh(1)%origin) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,1,self%no_mans_land)
    end do 
  endif
  if (self%llc_nat(1)+self%size(1)+self%ds(1) .gt. self%mesh(1)%origin + self%box(1)) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,1,self%no_mans_land)
    end do 
  endif

!  if (self%llc_nat(2)-self%ds(2) .lt. self%mesh(2)%origin) then 
!    do iv=1,self%nv
!      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,2,self%no_mans_land)
!    end do 
!  endif
!  if (self%llc_nat(2)+self%size(2)+self%ds(2) .gt. self%mesh(2)%origin + self%box(2)) then 
!    do iv=1,self%nv
!      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,2,self%no_mans_land)
!    end do 
!  endif
!
!  if (self%llc_nat(3)-self%ds(3) .lt. self%mesh(3)%origin) then 
!    do iv=1,self%nv
!      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,3,self%no_mans_land)
!    end do 
!  endif
!  if (self%llc_nat(3)+self%size(3)+self%ds(3) .gt. self%mesh(3)%origin + self%box(3)) then
!    do iv=1,self%nv
!      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,3,self%no_mans_land)
!    end do 
!  endif

  call self%solver_t%update

  !call self%courant_condition
  call self%check_density('before timestep')

!  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density('after timestep')
  call trace_end
END SUBROUTINE update
!===============================================================================
END MODULE experiment_mod
