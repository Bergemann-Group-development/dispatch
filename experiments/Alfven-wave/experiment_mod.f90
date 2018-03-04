!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the C-shock problem to test ambipolar diffusion
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
!  USE scalar_mod
  USE boundary_mod
  USE timestep_mod
  USE index_mod
!  USE solver_mod
  USE mhd_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type

  type(index_t):: idx

CONTAINS

!> =============================================================================
!> Setup a patch with a sharp shock front as present in a C-shock
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real, dimension(:,:,:,:), pointer:: b, p
  real, dimension(:,:,:), pointer:: d, vx, vy, vz, Bx, By, Bz, px, py, pz, s, e1, et
  real,save:: P_l, P_r, rho_i, rho_l, vx_l, vy_l, vz_l=0., &
              Bx_l, By_l, Bz_l=0., rho_r, vx_r, vy_r, vz_r=0., Bx_r, By_r, Bz_r=0., &
              x0=0.0, rhov2, B2,et_l,et_r,lambda,k_w
  namelist /experiment_params/ P_l, P_r, rho_l, vx_l, vy_l, Bx_l, By_l, rho_r, &
                               vx_r, vy_r, Bx_r, By_r, x0

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  self%nv = 8
  self%mhd = .true.
  call self%mhd_t%init
 
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  if (self%kind(1:6) == 'ramses') then 
    et => self%mem(:,:,:,idx%e,self%it,1)
    allocate (e1(self%gn(1),self%gn(2),self%gn(3)))
  endif

  d => self%mem(:,:,:,idx%d,self%it,1)
  if (self%kind(1:7) == 'stagger') s => self%mem(:,:,:,idx%s,self%it,1)
  p => self%mem(:,:,:,idx%px:idx%pz,self%it,1)
  B => self%mem(:,:,:,idx%bx:idx%bz,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  print*, 'gn is:', self%gn(1)

  print*,'self%mesh%lo!!!',self%mesh%lo

  lambda=1.
  k_w=2.*acos(-1.0d0)/lambda

  do i=self%mesh(1)%lb,self%mesh(1)%ub
    do j=self%mesh(2)%lb,self%mesh(2)%ub
      do k=self%mesh(3)%lb,self%mesh(3)%ub
        x = self%mesh(1)%p + self%mesh(1)%r(i)
        y = self%mesh(2)%p + self%mesh(2)%r(j)
        z = self%mesh(3)%p + self%mesh(3)%r(k)
        d (i,j,k)= 1.
        px(i,j,k)= d(i,j,k) * 1.*cos(k_w*z) !  0.9992286*cos(k*zz)
        py(i,j,k)=-d(i,j,k) * 1.*sin(k_w*z) ! -0.9992286*sin(k*zz)
        pz(i,j,k)= 0.
        Bx(i,j,k)= 1.*cos(k_w*z) 
        By(i,j,k)=-1.*sin(k_w*z)
        Bz(i,j,k)= 1.
        if (self%kind(1:6) == 'ramses') then
          et(i,j,k)= 1./(self%gamma-1d0) + 0.5*d(i,j,k)*(cos(k_w*z)**2 + sin(k_w*z)**2) + 0.5 * (Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)  
        end if
        if (self%kind(1:7) == 'stagger') then
          s(i,j,k)= (log(1.)-log(d(i,j,k))*self%gamma)/(self%gamma-1.0) * d(i,j,k)
        end if
      end do
    end do
    print*,'x',x
  end do

  call self%clear(bits%root_grid)
  self%periodic = .true.
  self%mhd = .true.
  print*, 'done!'

  print*,'self%mem(1,1,1,:,self%it,1) in experiment setup',self%mem(1,1,1,:,self%it,1)
!  if (self%kind(1:6) == 'ramses') deallocate (e1)
!    nullify (d, e1, et, px, py, pz, Bx, By, Bz)
  if (self%initial%time > 0d0) then
    self%time = self%initial%time
    self%t(self%it) = self%initial%time
  end if

  call trace_end
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv,i1,i2
  logical :: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  if (self%llc_nat(1)-self%ds(1) .lt. self%mesh(1)%origin) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,1)
!      do i1=self%mesh(1)%lb,self%mesh(1)%lo
!        self%mem(i1,:,:,iv,self%it,1)=self%mem(self%mesh(1)%li,:,:,iv,self%it,1)
!      end do
    end do 
  endif
  if (self%llc_nat(1)+self%size(1)+self%ds(1) .gt. self%mesh(1)%origin + self%box(1)) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,1)
!      do i1=self%mesh(1)%uo,self%mesh(1)%ub
!        self%mem(i1,:,:,iv,self%it,1)=self%mem(self%mesh(1)%ui,:,:,iv,self%it,1)
!      end do
    end do 
  endif

  if (self%llc_nat(2)-self%ds(2) .lt. self%mesh(2)%origin) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,2)
!      do i2=self%mesh(2)%lb,self%mesh(2)%lo
!        self%mem(:,i2,:,iv,self%it,1)=self%mem(:,self%mesh(2)%li,:,iv,self%it,1)
!      end do
    end do 
  endif
  if (self%llc_nat(2)+self%size(2)+self%ds(2) .gt. self%mesh(2)%origin + self%box(2)) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,2)
!      do i2=self%mesh(2)%uo,self%mesh(2)%ub
!        self%mem(:,i2,:,iv,self%it,1)=self%mem(:,self%mesh(2)%ui,:,iv,self%it,1)
!      end do
    end do 
  endif


  if (self%llc_nat(3)-self%ds(3) .lt. self%mesh(3)%origin) then 
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .true.,3)
    end do 
  endif
  if (self%llc_nat(3)+self%size(3)+self%ds(3) .gt. self%mesh(3)%origin + self%box(3)) then
    do iv=1,self%nv
      call boundary%condition (iv, self%mem(:,:,:,iv,self%it,1), self%mesh, .false.,3)
    end do 
  endif

  call self%mhd_t%update
  
  !call self%courant_condition
!  call self%check_density('before timestep')

!  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
!  call self%check_density('after timestep')
  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE experiment_mod
