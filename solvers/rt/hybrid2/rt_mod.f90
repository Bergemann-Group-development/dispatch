MODULE rt_mod
  USE bits_mod
  USE io_mod
  USE trace_mod
  USE link_mod
  USE omp_mod
  USE ray_geometry_v3_mod
  USE task_mod
  USE mesh_mod
  USE mhd_mod
  USE rt_boundary_mod
  implicit none
  private
!-----------------------------------------------------------------------
  type, extends(task_t), public :: rt_t
    class(mhd_t),       pointer :: mhd_patch  => null()                 ! pointer to MHD task
    class(link_t),      pointer :: mhd_link   => null()                 ! pointer to MHD link to access it's nbors
    class(link_t),      pointer :: my_nb_link => null()                 ! pointer for corrections towards itself as nbor of internal rt task
    class(link_t),      pointer :: link       => null()                 ! pointer to self link in the task list
    class (ray_bndl_t), pointer :: bundle     => null() 
    ! common items
    integer :: n_ray_dirs = 3
    integer :: geom = 1                                                 !cartesian = 1, spherical = 2, cylindrical = 3
    integer :: n_var                                                    ! number of variables
    integer :: rt_freq = 1                                              ! how often rt tasks are done wrt experiment
    integer :: rt_type, li(3), ui(3)
    integer :: ibc = 0, obc = 0, ibct = 3, obct = 0                     ! inner/outer boundary condition, treatment
    real(kind=8) :: cell_ds(3)
    real(8):: offset(3)                                                 ! this is used to compute indices, locations
    !---------------------------------------------------------------------------
    ! The pointers below, and the procedure 'init_mesh' that sets them, are
    ! currently not used, but are to enable future implementations of RT solvers
    ! that only know about density, temperature, ! and the mesh.  Such solvers
    ! would not need a pointer to an experiment_t data type.
    !---------------------------------------------------------------------------
    class(mesh_t),            pointer :: mesh        => null()          ! pointer to mesh
    real, dimension(:,:,:),   pointer :: density     => null()          ! pointer to density
    real, dimension(:,:,:),   pointer :: temperature => null()          ! pointer to temperature
    real, dimension(:,:,:,:), pointer :: q           => null()          ! (x,y,z,nt)
    logical :: do_output = .false.    
  contains
    procedure :: absolute_point_position
    procedure :: cell2pos
    procedure :: find_upstream
    procedure :: init => init_rt
    procedure :: init_mesh
    procedure :: init_obc
    procedure :: point_in_cell
    procedure :: pos2cell
    procedure :: ray2cell
    procedure :: ray2pos
    procedure :: the_right_wall
    procedure :: update => update_rt
    procedure :: update_time
  end type
!-----------------------------------------------------------------------------
! Sequential header, for MPI and I/O
!-----------------------------------------------------------------------------
  integer, save:: rt_n_header = (14*4 + 4*8 + 5*8*3 + 2*8*8 + 1*4*8 + 2*4*3)/4
  type rt_header_t
    sequence
    integer:: version=1
    integer:: n_header                                                  ! number of header entries
    integer:: n_data                                                    ! number of data entries
    integer:: id
    integer:: status
    integer:: rank
    integer:: istep
    integer:: iout
    integer:: level
    integer:: nt
    integer:: nv                        ! number of variables
    integer:: it
    integer:: nq=0
    integer:: npad=0
    real(8):: time
    real(8):: out_next
    real(8):: print_next
    real(8):: dtime
    real(8):: position(3)
    real(8):: velocity(3)
    real(8):: size(3)
    real(8):: box(3)                    ! the periodic wrapping size
    real(8):: ds(3)
    real(8), dimension(8):: t
    real(8), dimension(8):: dt
    integer, dimension(8):: iit
    integer:: n(3)=16                   ! convenience
    integer:: ng(3)=2                   ! convenience
  end type
  CONTAINS

!=======================================================================
!> Initialize rt_t tasks, by first getting a uniq id, and then allocating
!> the extra arrays needed to handle time slices.
!=======================================================================
SUBROUTINE init_rt (self)
  class(rt_t) :: self
  !---------------------------------------------------------------------
  call self%task_t%init
  self%nt = 5
  self%it = 1
  self%new = 2
  self%grace = 0.0
  allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
END SUBROUTINE init_rt

!=======================================================================
!> Initialize pointers to mesh, density and temperature
!=======================================================================
SUBROUTINE init_mesh (self, mesh, density, temperature)
  class(rt_t) :: self
  class(mesh_t), pointer :: mesh
  real, dimension(:,:,:), pointer :: density, temperature
  self%mesh => mesh
  self%density => density
  self%temperature => temperature
END SUBROUTINE init_mesh

SUBROUTINE update_rt (self)
  class(rt_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: rt_t%update_rt was called, but it is an empty subroutine"
END SUBROUTINE update_rt

SUBROUTINE set_reflective_BC(self)
  class(rt_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: rt_t: set_reflective_BC not implemented yet!"
END SUBROUTINE set_reflective_BC

SUBROUTINE set_periodic_BC(self)
  class(rt_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: rt_t: set_periodic_BC not implemented yet!"
END SUBROUTINE set_periodic_BC

!=======================================================================
!> Determine of patches 'me' and 'nbor' share the wall(s) indicated by
!> the triple index 'wall', where e.g. the right and left x-walls are
!> indicated by wall(1) = +1/-1.
!=======================================================================
LOGICAL FUNCTION the_right_wall (self, nbor, wall)
  class(rt_t) :: self
  class (link_t),        pointer :: nbor
  class (mhd_t), pointer :: mhd, nb_mhd
  integer :: wall(3), x
  real (kind=8) :: d(3)
  logical :: answer, allw(3)
!-----------------------------------------------------------------------
  the_right_wall = .false.
  associate (my_type => self%MHD_link%task)
  SELECT type (my_type)
  class is (mhd_t)
    mhd => my_type
  class is (rt_t)
    mhd => my_type%mhd_patch
  class default
    error stop '(1) in rt_mod::the_right_wall'
  END SELECT
  end associate
  !
  associate (nbor_type => nbor%task)
  SELECT type (nbor_type)
  class is (mhd_t)
    nb_mhd => nbor_type
  class is (rt_t)
    nb_mhd => nbor_type%mhd_patch                                       ! if I got an RT task to test, switch to corresponding MHD patch
  class default
    error stop '(2) in rt_mod::the_right_wall'
  end SELECT
  end associate
  !
  allw = .false.
  if (abs(nb_mhd%level-mhd%level)<2) then                               ! connect only +/- 1 or same level

    do x = 1,3

    if (wall(x) > 0) then
      if (mhd%position(x) < nb_mhd%position(x)) then
        if (mhd%position(x)+mhd%size(x)*0.5d0 > &
        nb_mhd%position(x) - nb_mhd%size(x)*0.5d0 - nb_mhd%ds(x)) allw(x) = .true.
      end if
    else if (wall(x) < 0) then
      if (mhd%position(x) > nb_mhd%position(x)) then
        if (mhd%position(x)-mhd%size(x)*0.5d0 < &
        nb_mhd%position(x) + nb_mhd%size(x)*0.5d0 + nb_mhd%ds(x)) allw(x) = .true.
      end if
    else
      allw(x) = .true.
    end if

    end do
  end if
  if (all(allw)) the_right_wall = .true.
  nullify (mhd, nb_mhd)

END FUNCTION the_right_wall

!=======================================================================
!> Return a linked list starting with nbor_l with upstream walls of the
!> self rt_t task.
!=======================================================================
SUBROUTINE find_upstream (self, nbor_l, found)
  class(rt_t)  :: self
  class(link_t), pointer :: nbor_l, nul
  class(link_t), pointer :: nbor, my_mhd_l, new
  class(mhd_t), pointer :: my_task, nbor_task
  integer :: found, wall(3)
! ----------------------------------------------------------------------
  call trace_begin('rt_t%find_upstream')
  
  select case(self%rt_type)                                             ! select the right upstream wall
  case(2)
    wall = self%bundle%lwall
  case(3)
    wall = self%bundle%uwall
  end select
  
  found = 0
  allocate (nbor_l)
  my_task => self%mhd_patch
  nbor => self%mhd_link%nbor
  do while (associated(nbor))
    associate (nbtask => nbor%task)
    SELECT type (nbtask)
    class is (mhd_t)
      nbor_task => nbtask
      if (.not.nbor_is_my_child(my_task, nbor_task)) then
        ! if the nbor has a wall common with our upstream wall
        ! allocate a new link, pointing to the nbor task, and
        ! prepend it to the list starting with nbor_l
        if (self%the_right_wall (nbor, wall)) then
          allocate(new)
          new%task => nbor_task
          new%link => nbor%link
          if (associated(nbor_l%next)) then
            new%next => nbor_l%next
            nbor_l%next%prev => new
            nbor_l%next => new
          else
            nbor_l%next => new
            new%prev => nbor_l
          end if
          found = found + 1
        end if
      end if
    end SELECT
    end associate
    nbor => nbor%next
  end do
  nul => nbor_l
  nbor_l => nbor_l%next
  deallocate(nul)
  call trace_end
END SUBROUTINE find_upstream

!=======================================================================
!=======================================================================
SUBROUTINE init_obc (self)
  class(rt_t) :: self
!-----------------------------------------------------------------------
  call self%set(bits%has_BC)
  self%obc = 1
  self%my_nb_link%check = .true.
END SUBROUTINE init_obc

!=======================================================================
!> Check if 'nbor' is fully contained inside 'me'
!=======================================================================
FUNCTION nbor_is_my_child (me, nbor) result (out)
  class (mhd_t), pointer :: me, nbor
  real(kind=8) :: me_li, me_ui, nbor_li, nbor_ui
  logical :: out
!-----------------------------------------------------------------------
  out = .false.
  if ((all(me%position + me%size*0.5_8 > nbor%position + nbor%size*0.5_8)) &
.and. (all(me%position - me%size*0.5_8 < nbor%position - nbor%size*0.5_8))) &
  out = .true.
END FUNCTION nbor_is_my_child

!=======================================================================
!> Find the absolute position of index 'ray' in ray packet 'rays'.
!> If the optional bn parameter is present, the answer depends on
!> the direction, otherwise 'pt' chooses the point along the ray.
!=======================================================================
FUNCTION absolute_point_position (self, rays, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  class (ray_t), pointer :: rays
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  call trace_begin ('rt_t%absolute_point_position', 4)
  if (present(bn)) then
    select case(self%rt_type)
    case(2)
      ! self%bundle
      cp = rays%fwd_coords(:,ray)
    case(3)
      cp = rays%rev_coords(:,ray)
    end select
  else
    cp = rays%coords(:,pt,ray)
  end if
  cell_pos = self%position-0.5*self%size+(cp-self%li)*self%cell_ds
  call trace_end
END FUNCTION absolute_point_position

!=======================================================================
!> Return the cell integer coordinates of point 'pt' of index 'ray' in
!> ray packet 'rays'.
!=======================================================================
FUNCTION ray2cell(self, rays, ray, pt) result (cell)
  class (rt_t) :: self
  class (ray_t), pointer :: rays
  integer :: ray, pt, cell(3)
!-----------------------------------------------------------------------
  cell = rays%coords(:,pt,ray)
END FUNCTION ray2cell

!=======================================================================
!> Return the absolute position of integer index triple 'cell'
!=======================================================================
FUNCTION cell2pos (self, cell) result(p)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  p = self%position - 0.5_8*self%size + (cell-self%li)*self%cell_ds
END FUNCTION cell2pos

!=======================================================================
!> Convert absolute position to integer index triple 'cell'
!=======================================================================
FUNCTION pos2cell (self, p) result (cell)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  cell = nint((p-self%position+0.5_8*self%size)/self%cell_ds)+self%li
END FUNCTION pos2cell

!=======================================================================
!> Convert absolute position to integer index triple 'cell'
!=======================================================================
FUNCTION ray2pos (self, rays, ray, pt) result (p)
  class(rt_t)  :: self
  class(ray_t), pointer :: rays
  real(kind=8) :: p(3)
  integer :: ray, pt
!-----------------------------------------------------------------------
  p = self%cell2pos(self%ray2cell(rays,ray,pt))
END FUNCTION ray2pos

!=======================================================================
!> Check if an absolute coordinate is inside the index triplet 'cell'
!=======================================================================
FUNCTION point_in_cell (self, pt, cell, ds) result (out)
  class(rt_t)  :: self
  real(kind=8) :: pt(3), cell(3), li(3), ui(3), ds(3)
  logical :: out
!-----------------------------------------------------------------------
  out = .false.
  li = cell - 0.5*ds
  ui = cell + 0.5*ds
  if ((all(li < pt)).and.(all(ui > pt))) out = .true.
END FUNCTION point_in_cell

!=======================================================================
!> Update the time of an rt_t task.  Set self%time to the earliest point
!> in time when the task should next be updated, while setting the time
!> slot in self%t() to the actual time of the MHD task, so inter- and
!> extrapolations are correct.
!>
!> NOTE: This procedure should be called exactly when the physical data
!> on which the task is based is collected
!=======================================================================
SUBROUTINE update_time(self)
  class(rt_t) :: self
!-----------------------------------------------------------------------
  if (self%time < 0) then
    self%t(self%new) = self%mhd_patch%time
    self%time  = self%mhd_patch%time + self%dtime*(self%rt_freq)
  else
    self%dtime = self%mhd_patch%dtime * self%rt_freq
    call self%rotate
  end if
END SUBROUTINE update_time



END MODULE rt_mod
