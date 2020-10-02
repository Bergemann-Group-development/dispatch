!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Cylinder-organized Cartesian patches
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE cylinder_mod
  USE io_mod
  USE io_unit_mod
  USE task_mod
  USE experiment_mod
  USE scaling_mod
  USE units_mod
  USE task_list_mod
  USE bits_mod
  USE trace_mod
  implicit none
  private
  type, public:: cylinder_t
    integer:: n_r
    integer:: n_phi
    integer:: n_z
    real(8):: r
    real(8):: z
    real(8):: phi, dphi
    real(8):: lnr_min, dlnr
    class(task_list_t), pointer:: task_list
  contains
    procedure:: init
  end type
CONTAINS

!===============================================================================
!> Initialize a number of Cartesian patches in cylinder arrangement
!===============================================================================
SUBROUTINE init (self)
  class(cylinder_t):: self
  integer, save:: n_r=10, n_phi=100, n_z=10
  real(8), save:: r_min=1.0
  namelist /cylinder_params/ n_r, n_phi, n_z, r_min
  !.............................................................................
  integer:: i_r, i_z, i_phi, iostat
  class(experiment_t), pointer:: cpatch
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('cylinder_t%init')
  call io%lock%set ('cylinder_t%init')
  rewind (io%input)
  read (io%input, cylinder_params, iostat=iostat)
  call io%namelist_warning ('cylinder_params', iostat)
  write (io%output, cylinder_params)
  call io%lock%unset ('cylinder_t%init')
  self%n_r   = n_r
  self%n_z   = n_z
  self%n_phi = n_phi
  self%lnr_min = log(r_min)
  !-----------------------------------------------------------------------------
  !> Loop over positions, creating a task list with cpatches
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  call self%task_list%init
  self%dlnr = 1d0/n_phi
  self%dphi = 2d0*cgs%pi*self%dlnr
  do i_r=1,n_r
    do i_z=1,n_z
      do i_phi=1,n_phi
        allocate (cpatch)
        cpatch%r = exp(self%lnr_min + (i_r-1)*self%dlnr)
        cpatch%size = cpatch%r*self%dphi/(1d0-self%dphi/2d0)
        cpatch%phi = (i_phi-1)*self%dphi
        cpatch%z = (i_z - (1+n_z)*0.5d0)*cpatch%size(3)
        cpatch%position = &
          [cpatch%r*cos(cpatch%phi), cpatch%r*sin(cpatch%phi), cpatch%z]
        task => cpatch
        allocate (cpatch%link)
        cpatch%link%task => cpatch
        call cpatch%link%init
        call cpatch%init
        cpatch%size = cpatch%r*self%dphi/(1d0-self%dphi/2d0)
        cpatch%ds = cpatch%size/cpatch%n
        cpatch%box = 1.0
        cpatch%position = &
          [cpatch%r*cos(cpatch%phi), cpatch%r*sin(cpatch%phi), cpatch%z]
        call self%task_list%append_link (cpatch%link)
        !-----------------------------------------------------------------------
        ! Set boundary bits
        !-----------------------------------------------------------------------
        if (i_r == 1) then
          call cpatch%boundaries%set (bits%xl)
        endif
        if (i_r == n_r) then
          call cpatch%boundaries%set (bits%xu)
        end if
        if (i_z == 1) then
          call cpatch%boundaries%set (bits%zl)
        endif
        if (i_z == n_z) then
          call cpatch%boundaries%set (bits%zu)
        end if
      end do
    end do
  end do
  call trace%end()
END SUBROUTINE init

END MODULE cylinder_mod
