!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE scaling_mod
  USE units_mod
  USE link_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: dnload
  end type
  integer, save:: verbose=0
  real, save:: density_cgs=1e-20, cs_kms=0.18
CONTAINS

!===============================================================================
!> Setup an initial state.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real:: d, e
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /experiment_params/ density_cgs, cs_kms, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    first_time = .false.
  end if
  call io%lock%unset
  d = density_cgs/scaling%d
  e = (cs_kms*cgs%kms/scaling%u)**2*d/(eos%gamma*(eos%gamma-1.0))
  self%mem(:,:,:,self%idx%d,:,1) = d
  self%mem(:,:,:,self%idx%e,:,1) = e
  self%unsigned(self%idx%d) = .true.
  self%unsigned(self%idx%e) = .true.
  self%pervolume(self%idx%px) = .true.
  self%pervolume(self%idx%py) = .false.
  self%pervolume(self%idx%pz) = .true.
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  integer:: iv, i1, i2, m, m1, m2, m3, m4
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%solver_t%update
  i1 = self%task_list%na/2
  i2 = i1+1
  m = self%gn(1)/2
  m1 = self%mesh(1)%ui-1
  m2 = self%mesh(1)%ub
  m3 = 1
  m4 = self%mesh(1)%li+1
  if (verbose > 0) then
    do iv=1,2
      if (self%id==i1) write (stdout,'(i5,2i2,1p,3e12.3,05x,4e12.3)') self%istep, self%new, iv, self%time, self%dtime, self%u_max, self%mem(m1:m2,m,m,iv,self%it ,1)
      if (self%id==i1) write (stdout,'(i5,2i2,1p,3e12.3,05x,4e12.3)') self%istep, self%new, iv, self%time, self%dtime, self%u_max, self%mem(m1:m2,m,m,iv,self%new,1)
      if (self%id==i2) write (stdout,'(i5,2i2,1p,3e12.3,55x,4e12.3)') self%istep, self%new, iv, self%time, self%dtime, self%u_max, self%mem(m3:m4,m,m,iv,self%it ,1)
      if (self%id==i2) write (stdout,'(i5,2i2,1p,3e12.3,55x,4e12.3)') self%istep, self%new, iv, self%time, self%dtime, self%u_max, self%mem(m3:m4,m,m,iv,self%new,1)
    end do
  end if
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  call self%gpatch_t%dnload (only)
  call trace%end ()
END SUBROUTINE dnload

END MODULE
