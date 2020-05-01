!===============================================================================
!> Construct a spherical layer, with volleyball distribution of patches
!===============================================================================
MODULE volleyball_layers_mod
  USE io_mod
  USE sphere_points_mod
  USE dll_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
  !.............................................................................
  type, public:: volleyball_layer_t
    integer:: ntheta=10
    real(8):: theta0, dtheta, radius
    type(dll_t):: list
  end type
  !.............................................................................
  type, public:: volleyball_layers_t
    integer:: nlayer
    type(dll_t):: list
  contains
    procedure:: init
  end type
  !.............................................................................
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Given the number of intervals in the theta direction, generate the patch
!> positions, saving the point information in a dll_t data type
!===============================================================================
SUBROUTINE init (self)
  class(volleyball_layers_t):: self
  !.............................................................................
  type(sphere_points_t), pointer:: point, save, tmp
  class(dll_node_t), pointer:: dll_node, layer_node
  type(volleyball_layer_t), pointer:: layer
  integer:: nlayer=2, iostat, i
  real(8):: radius=1d0
  namelist /volleyball_layers_params/ verbose, nlayer, radius
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layers_t%init')
  call io%lock%set
  rewind (io_unit%input)
  read (io_unit%input, volleyball_layers_params, iostat=iostat)
  write (io_unit%output, volleyball_layers_params)
  self%nlayer = nlayer
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Construct nlayer layers and append to self%list
  !-----------------------------------------------------------------------------
  do i=1,nlayer
    allocate (layer, dll_node)
    call init_layer (self, layer, radius)
  end do
  write (io_unit%output,*) 'layers:', self%list%n
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Given the number of intervals in the theta direction, generate the patch
!> positions, saving the point information in a dll_t data type
!===============================================================================
SUBROUTINE init_layer (self, layer, radius)
  class(volleyball_layers_t):: self
  class(volleyball_layer_t):: layer
  real(8):: radius
  !.............................................................................
  type(sphere_points_t), pointer:: point, save
  class(dll_node_t), pointer:: dll_node
  integer:: ntheta=10, iostat
  namelist /volleyball_layer_params/ ntheta
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layer_t%init_layer')
  call io%lock%set
  read (io_unit%input, volleyball_layer_params, iostat=iostat)
  write (io_unit%output, volleyball_layer_params)
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Construct a set of patch positions that define a volleyball front face
  !-----------------------------------------------------------------------------
  layer%radius = radius
  layer%ntheta = ntheta
  layer%theta0 = atan(sqrt(0.5d0))
  layer%dtheta = layer%theta0/layer%ntheta
  allocate (point)
  call point%init
  if (verbose > 0) &
    call point%info ('init')
  allocate (save, source=point)
  !-----------------------------------------------------------------------------
  ! March in the theta direction
  !-----------------------------------------------------------------------------
  do while (point%is_inside (layer%theta0, layer%dtheta))
    point%r = radius
    point%size = radius*layer%dtheta
    write (stdout,*) point%size, 'size'
    !---------------------------------------------------------------------------
    ! Add half a patch size north-east offset
    !---------------------------------------------------------------------------
    call point%go_north (0.5d0*layer%dtheta)
    call point%go_east (0.5d0*layer%dtheta)
    if (verbose > 0) &
      call point%info ('row start')
    !---------------------------------------------------------------------------
    ! March in the phi direction while patches are inside the front section,
    ! while appending a point and its three mirror images to the list of points
    !---------------------------------------------------------------------------
    do while (point%is_inside (layer%theta0, layer%dtheta))
      call append
      point%y = -point%y
      call append
      point%z = -point%z
      call append
      point%y = -point%y
      call append
      point%z = -point%z
      call point%go_east (layer%dtheta)
      if (verbose > 1) &
        call point%info ('east')
    end do
    point = save
    call point%go_north (layer%dtheta)
    if (verbose > 0) &
      call point%info ('north')
    save = point
  end do
  deallocate (save)
  radius = radius*(1d0-layer%dtheta)
  call trace%end
contains
  !-----------------------------------------------------------------------------
  ! Append a (copy of a) point to the list of points
  !-----------------------------------------------------------------------------
  subroutine append
    class(sphere_points_t), pointer:: tmp
    allocate (tmp, source=point)
    allocate (dll_node)
    dll_node%car => tmp
    call self%list%append (dll_node)
  end subroutine
END SUBROUTINE init_layer

END MODULE volleyball_layers_mod
