!===============================================================================
!> EM-fields module -- for now explicit and simplest possible
!===============================================================================
MODULE fields_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE vector_ops
  USE kinds_mod
  implicit none
  private
  type, public:: fields_t
    class(patch_t), pointer:: patch  => null()
    integer:: lb(3), ub(3)
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
  end type
  integer :: verbose = 1
  real    :: time_center = 0.5
  logical :: first_time = .true.
CONTAINS

!===============================================================================
!> Read field parameters
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(fields_t):: self
  class(patch_t), pointer:: patch
  integer:: iostat
  namelist /fields_params/ verbose, time_center
  !-----------------------------------------------------------------------------
  self%patch => patch
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read   (io%input, fields_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read   (io%input, fields_params)
    end if
    call io%namelist_warning ('fields_params', iostat)
    write (io%namelist, fields_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
END SUBROUTINE pre_init

!===============================================================================
!> Allocate field arrays
!===============================================================================
SUBROUTINE init (self, patch)
  class(fields_t):: self
  class(patch_t), pointer:: patch
  !-----------------------------------------------------------------------------
END SUBROUTINE init

!===============================================================================
!> Update the fields, using Maxwell's equations in explicit leap-frog form
!===============================================================================
SUBROUTINE update (self, J)
  class(fields_t):: self
  class(patch_t), pointer:: patch
  real(var_k), pointer, dimension(:,:,:,:):: J, B0, B1, E0, E1, dBdt, dEdt
  integer:: m(6)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('fields_t%update', itimer=itimer)
  patch => self%patch
  m = shape(patch%mem)
  B0 => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%it ,1)
  E0 => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,patch%it ,1)
  B1 => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%new,1)
  E1 => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,patch%new,1)
  allocate (dEdt(m(1),m(2),m(3),3),dBdt(m(1),m(2),m(3),3))
  dBdt = curl_up (patch%ds, E0)                               ! negative curl
  B1 = B0 - patch%dtime * dBdT
  dEdt = curl_dn (patch%ds, B1)                               ! positive curl
  dEdt = dEdt - J                                             ! Ampere
  E1 = E0 + patch%dtime * dEdT
  deallocate (dEdt, dBdt)
  call trace%end (itimer)
END SUBROUTINE update

END MODULE fields_mod
