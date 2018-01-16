!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Module for an immersed boundary patch: a set of cells interior to a domain
!> which act as a physical boundary.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE immersed_boundary_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE patch_mod
  USE mesh_mod
  USE mhd_mod
  implicit none
  private

  type, extends(patch_t), public :: immersed_boundary_t
    ! possibilities: reflecting, absorbing, constant
    character(len=64) :: boundary_type
  contains
    procedure :: set_type => set_immersed_boundary_type
    procedure :: init => init_immersed_boundary
    procedure :: update => update_immersed_boundary
    procedure :: conditions => slip_boundary_conditions
    procedure :: output
  end type

  type, extends(immersed_boundary_t), public :: immersed_cylinder
    real(4) :: radius, height, mass, cyl_density
  contains
  end type
  
  type, extends(immersed_boundary_t), public :: immersed_sphere
    real(4) :: radius, mass, sph_density
  contains
  end type

  ! possible boundary types
  character(len=64), dimension(3), parameter :: type_options = ["reflecting", &
                                                                "absorbing ", &
                                                                "constant  "]

CONTAINS

SUBROUTINE void (self, patch, ii, it)
  class (immersed_boundary_t):: self
  class (patch_t), target :: patch
  integer, intent(in) :: it, ii(3)
  !-----------------------------------------------------------------------------
  call mpi%abort('ERROR: immersed_boundary_t%void called')
END SUBROUTINE

SUBROUTINE init_immersed_boundary (self)
  class(immersed_boundary_t) :: self
  !-----------------------------------------------------------------------------
  call trace_begin('immersed_boundary_t%init')

  call self%task_t%init()

  self%kind = 'immersed_boundary_t'

  call trace_end
END SUBROUTINE init_immersed_boundary

!> =============================================================================
!> Set the type for an immersed boundary patch (incl. the boundary type)
!> =============================================================================
SUBROUTINE set_immersed_boundary_type (self, bndy_type)
  class(immersed_boundary_t) :: self
  character(len=64), optional :: bndy_type
  !-----------------------------------------------------------------------------
  call trace_begin ('immersed_boundary%init')
  if (present(bndy_type)) then
    self%boundary_type = bndy_type
  else
    self%boundary_type = 'reflective'  ! default value
  endif

  call trace_end
END SUBROUTINE set_immersed_boundary_type

!> =============================================================================
!> =============================================================================
SUBROUTINE update_immersed_boundary (self)
  class(immersed_boundary_t) :: self
  !-----------------------------------------------------------------------------
  call trace_begin('immersed_boundary_t%update')

  ! FIXME: Do nothing.
  self%dtime = self%dt(self%it)
  call self%task_t%update()

  call trace_end
END SUBROUTINE update_immersed_boundary

!> =============================================================================
!> Apply the appropriate boundary conditions to the input patch.
!> =============================================================================
SUBROUTINE slip_boundary_conditions (self, patch, ii, it)
  class(immersed_boundary_t) :: self
  class(patch_t), target :: patch
  class(mesh_t), pointer :: mx=>null(), my=>null(), mz=>null()
  integer, intent(in) :: it, ii(3)
  integer :: ix, iy, iz, offi
  real(4) :: x, y, radius
  !-----------------------------------------------------------------------------
  call trace_begin('immersed_boundary_t%conditions')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  ix = ii(1)
  iy = ii(2)
  iz = ii(3)
  
  select type (self)
  type is (immersed_cylinder)
    radius = self%radius
  type is (immersed_sphere)
    radius = self%radius
  class default
    print*,'WARNING! Slip boundary conditions have not been implemented for this class!'
  end select

  select type (mx)
  type is (Cartesian_mesh)
    y = my%centre_nat + patch%mesh(2)%r(iy) - self%position(2)
    x = mx%centre_nat + patch%mesh(1)%r(ix) - self%position(1)
    if (sqrt(x**2+y**2) < radius) then
      patch%mem(ix,iy,iz,patch%idx%d ,it,1) = 1.4 ! FIXME
      patch%mem(ix,iy,iz,patch%idx%e ,it,1) = 2.5
      patch%mem(ix,iy,iz,patch%idx%px,it,1) = 0.0
      patch%mem(ix,iy,iz,patch%idx%py,it,1) = 0.0
      patch%mem(ix,iy,iz,patch%idx%pz,it,1) = 0.0
    end if
  type is (cylindrical_mesh)
    y = my%centre_nat + patch%mesh(2)%r(iy) - self%position(2)
    if (y <= radius) then
      offi = 2 * my%lo + 1
      patch%mem(ix,iy,iz,patch%idx%d ,it,1) =  patch%mem(ix,offi-iy,iz,patch%idx%d ,it,1)
      patch%mem(ix,iy,iz,patch%idx%e ,it,1) =  patch%mem(ix,offi-iy,iz,patch%idx%e ,it,1)
      patch%mem(ix,iy,iz,patch%idx%px,it,1) =  patch%mem(ix,offi-iy,iz,patch%idx%px,it,1)
      patch%mem(ix,iy,iz,patch%idx%py,it,1) = -patch%mem(ix,offi-iy,iz,patch%idx%py,it,1)
      patch%mem(ix,iy,iz,patch%idx%pz,it,1) =  patch%mem(ix,offi-iy,iz,patch%idx%pz,it,1)
    end if
  type is (spherical_mesh)
    x = mx%centre_nat + patch%mesh(1)%r(ix) - self%position(1)
    if (x <= radius) then
      offi = 2 * mx%lo + 1
      patch%mem(ix,iy,iz,patch%idx%d ,it,1) =  patch%mem(offi-ix,iy,iz,patch%idx%d ,it,1)
      patch%mem(ix,iy,iz,patch%idx%e ,it,1) =  patch%mem(offi-ix,iy,iz,patch%idx%e ,it,1)
      patch%mem(ix,iy,iz,patch%idx%px,it,1) =  patch%mem(offi-ix,iy,iz,patch%idx%px,it,1)
      patch%mem(ix,iy,iz,patch%idx%py,it,1) = -patch%mem(offi-ix,iy,iz,patch%idx%py,it,1)
      patch%mem(ix,iy,iz,patch%idx%pz,it,1) =  patch%mem(offi-ix,iy,iz,patch%idx%pz,it,1)
    end if
  class default
    print*,"WARNING!: No method to interface the current patch with the immersed boundary!"
  end select

  nullify(mx, my, mz)
  call trace_end
END SUBROUTINE slip_boundary_conditions

!> =============================================================================
!> Write results to disk; preserve the format relative to `patch_mod::output_legacy`
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(immersed_boundary_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: i, l(3), u(3), no, jo, it
  real:: ds(3)
  logical:: first_time=.true.
  real(8):: time
  !-----------------------------------------------------------------------------
  call trace_begin('immersed_boundary_t%output')
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (io%time_derivs) then
    no = min(2,self%nw)
  else
    no = 1
  end if
  if (io%out_time==0.0) then
    time = self%time
  else if (io%time_derivs) then
    time = self%t(self%nt-2)
  else
    time = self%out_next
  end if
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  self%ioformat = 1
  write (filename,'(a,i4.4,"_",i4.4)') trim(io%outputname), self%id, self%iout
  if (io%verbose>=0) &
    print'(a,3x,"time:",2f13.6)', 'immersed_boundary_t%output: file= '//trim(filename), time, self%time
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  write (io%data_unit,'(i2.2)') self%ioformat
  write (io%data_unit,'(a)') trim(self%kind)
  write (io%data_unit,'(a)') trim(self%eos)
  write (io%data_unit,'(a)') trim(self%opacity)
  write (io%data_unit,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') time, self%dtime, &
        self%position, self%ds, self%velocity, self%llc_nat, self%llc_cart, &
        self%mesh%centre_nat, self%quality, self%gamma
  write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(self%box/self%ds+0.5), &
        u-l+1, self%li, self%ui, self%n, self%gn, self%ng, self%nv, self%level
  write (io%data_unit,'(3(1i8.1,1x),i2.1)') self%id, self%iout, self%istep, self%mesh_type
  close (io%data_unit)

  call trace_end
END SUBROUTINE output

END MODULE immersed_boundary_mod
