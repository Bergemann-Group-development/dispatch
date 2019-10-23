!===============================================================================
!> Construct a spherical layer, with volleyball distribution of patches
!===============================================================================
MODULE volleyball_layers_mod
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
  contains
    procedure:: init_layer
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
  class(dll_node_t), pointer:: dll_node
  type(volleyball_layer_t), pointer:: layer
  integer:: nlayer=2, iostat, i
  real(8):: radius=1d0
  namelist /volleyball_layers_params/ verbose, nlayer, radius
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layers_t%init')
  !$omp critical (input_cr)
  rewind (io_unit%input)
  read (io_unit%input, volleyball_layers_params, iostat=iostat)
  write (io_unit%output, volleyball_layers_params)
  self%nlayer = nlayer
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Construct nlayer layers and append to self%list
  !-----------------------------------------------------------------------------
  do i=1,nlayer
    allocate (layer, dll_node)
    call layer%init_layer (radius)
    dll_node%car => layer
    call self%list%append (dll_node)
    write (io_unit%output,*) 'layer points, radius:', layer%list%n, radius
  end do
  write (io_unit%output,*) 'layers:', self%list%n
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Given the number of intervals in the theta direction, generate the patch
!> positions, saving the point information in a dll_t data type
!===============================================================================
SUBROUTINE init_layer (self, radius)
  class(volleyball_layer_t):: self
  real(8):: radius
  !.............................................................................
  type(sphere_points_t), pointer:: point, save
  class(dll_node_t), pointer:: dll_node
  integer:: ntheta=10, iostat
  namelist /volleyball_layer_params/ ntheta
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layer_t%init_layer')
  !$omp critical (input_cr)
  read (io_unit%input, volleyball_layer_params, iostat=iostat)
  write (io_unit%output, volleyball_layer_params)
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Construct a set of patch positions that define a volleyball front face
  !-----------------------------------------------------------------------------
  self%radius = radius
  self%ntheta = ntheta
  self%theta0 = atan(sqrt(0.5d0))
  self%dtheta = self%theta0/self%ntheta
  allocate (point)
  call point%init
  if (verbose > 0) &
    call point%info ('init')
  allocate (save, source=point)
  !-----------------------------------------------------------------------------
  ! March in the theta direction
  !-----------------------------------------------------------------------------
  do while (point%is_inside(self%theta0))
    !---------------------------------------------------------------------------
    ! Add half a patch size north-east offset
    !---------------------------------------------------------------------------
    call point%go_north (0.5d0*self%dtheta)
    call point%go_east (0.5d0*self%dtheta)
    if (verbose > 0) &
      call point%info ('row start')
    !---------------------------------------------------------------------------
    ! March in the phi direction while patches are inside the front section,
    ! while appending a point and its three mirror images to the list of points
    !---------------------------------------------------------------------------
    do while (point%is_inside(self%theta0))
      call append
      point%y = -point%y
      call append
      point%z = -point%z
      call append
      point%y = -point%y
      call append
      point%z = -point%z
      call point%go_east (self%dtheta)
      if (verbose > 1) &
        call point%info ('east')
    end do
    point = save
    call point%go_north (self%dtheta)
    if (verbose > 0) &
      call point%info ('north')
    save = point
  end do
  deallocate (save)
  radius = radius*(1d0-self%dtheta)
  call trace%end
contains
  !-----------------------------------------------------------------------------
  ! Append a (copy of a) point to the list of points
  !-----------------------------------------------------------------------------
  subroutine append
    type(sphere_points_t), pointer:: tmp
    allocate (tmp, source=point)
    allocate (dll_node)
    dll_node%car => tmp
    call self%list%append (dll_node)
  end subroutine
END SUBROUTINE init_layer

END MODULE volleyball_layers_mod
