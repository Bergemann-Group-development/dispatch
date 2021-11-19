!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  USE units_mod
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    logical,     dimension(:), pointer, contiguous :: bottop
    real(var_k), dimension(:), pointer, contiguous :: temp, rho, pe
    procedure(boundary_mhd), pointer     :: boundary
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: boundary_mhd
    procedure:: boundary_sc0
    procedure:: boundary_scX
    procedure:: radiation
  end type
  integer, save:: verbose=1
CONTAINS

!===============================================================================
!> Initialize the boundary data type.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t) :: self
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /boundary_params/ verbose
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  patch => self%patch
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, boundary_params)
    end if
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  if (.not. io%master) verbose=0

  select case(patch%kind)
  case ('hlld_mhd_patch')
    allocate(self%bottop(2), self%rho(2), self%temp(2), self%pe(2))
    self%bottop = .false. ! default value
    self%boundary => boundary_mhd
  case ('m3d_rt_sc0')
    self%boundary => boundary_sc0
  case ('m3d_rt_sc')
    self%boundary => boundary_scX
  case default
  end select

  call self%standard_boundaries_t%init ()
END SUBROUTINE init


SUBROUTINE boundary_mhd(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer :: i, b

  do i=1,2 ! loop over bottom - top
    if (self%bottop(i)) then
      if (i == 1) b = lbound(mem,3)
      if (i == 2) b = ubound(mem,3)

      mem(:,:,b,idx%d ,1) = self%rho(i)
      mem(:,:,b,idx%e ,1) = self%temp(i)
      mem(:,:,b,6     ,1) = self%pe(i)
      mem(:,:,b,idx%px,1) = 0
      mem(:,:,b,idx%py,1) = 0
      mem(:,:,b,idx%pz,1) = 0
    endif
  end do
END SUBROUTINE boundary_mhd


SUBROUTINE boundary_sc0(self, mem)
  use background_mod, only: uppsala_get_totn
  use scs_mod, only: lte_pop, partition_function
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real(8), dimension(:,:,:), allocatable :: temp, rho, ne, totn, partition
  real(8), dimension(:,:,:,:), allocatable :: nh
  real(8) :: ggpgh
  integer :: i, lvl, b, nx, ny, nlevel

  nx = size(mem,1)
  ny = size(mem,2)
  nlevel = size(mem, 5)

  allocate(temp (nx, ny, 1))
  allocate(nh   (nx, ny, 1, 6))
  allocate(rho, totn, partition, mold=temp)

  do i=1,2 ! loop over bottom - top
    if (self%bottop(i)) then
      if (i == 1) b = lbound(mem,3)
      if (i == 2) b = ubound(mem,3)

      temp = dble(self%temp(i))
      rho  = dble(self%rho(i) * scaling%d)
      ne   = dble(self%pe(i) * scaling%p / (cgs%k_b * temp))

      call uppsala_get_totn(temp, rho, ne, nh, totn, ggpgh)
      partition = partition_function(temp, ne, totn)

      do lvl=1, nlevel
        mem(:,:,b,2,lvl) = lte_pop(lvl, temp(:,:,1), ne(:,:,1), partition(:,:,1))
        mem(:,:,b,1,lvl) = mem(:,:,b,2,lvl)
      end do
    endif
  end do
END SUBROUTINE boundary_sc0


SUBROUTINE boundary_scX(self, mem)
  use background_mod, only: uppsala_background_opac
  use parallel_mod, only: fs
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real(8), dimension(:,:,:), allocatable :: temp, rho, ne, chi_c, eta_c, Jnu
  integer :: i, b, nx, ny, f, nfreq

  nx = size(mem,1)
  ny = size(mem,2)
  nfreq = size(mem,5)

  allocate(temp (nx, ny, 1))
  allocate(rho, ne, chi_c, eta_c, Jnu, mold=temp)

  Jnu = 0

  do i=1,2 ! loop over bottom - top
    if (self%bottop(i)) then
      if (i == 1) b = lbound(mem,3)
      if (i == 2) b = ubound(mem,3)

      temp = dble(self%temp(i))
      rho  = dble(self%rho(i) * scaling%d)
      ne   = dble(self%pe(i) * scaling%p / (cgs%k_b * temp))

      do f=1,nfreq
        call uppsala_background_opac(f+fs-1, temp, rho, ne, Jnu, chi_c, eta_c)

        mem(:,:,b,1,f) = chi_c(:,:,1)
        mem(:,:,b,2,f) = eta_c(:,:,1) / chi_c(:,:,1)

        if (i == 1) mem(:,:,b,3,f) = mem(:,:,b,2,f) ! lower boundary
        if (i == 2) mem(:,:,b,3,f) = 0              ! upper boundary
      end do
    endif
  end do
END SUBROUTINE boundary_scX


SUBROUTINE radiation (self, mem, src, rk, tt, dir, mu)
  class(boundaries_t) :: self
  real, dimension(:,:,:,:), pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt,mem2
  real   :: mu(3)
  integer:: dir
  !.....................................................................
  ! keep this just for rt_mod to compile
END SUBROUTINE radiation

END MODULE boundaries_mod
