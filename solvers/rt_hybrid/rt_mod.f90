MODULE rt_mod
  USE experiment_mod
  USE bits_mod
  USE io_mod
  USE link_mod
  USE omp_mod
  USE ray_geometry_mod
  USE task_mod
  USE experiment_mod
  USE rt_boundary_mod
  USE microphysics_mod
!----------------------------------------------------------------------- 
  type, extends(task_t), public :: rt_t
    type (experiment_t), pointer :: mhd_patch  => null()  ! pointer to MHD task
    type (link_t)      , pointer :: mhd_link   => null()  ! pointer to MHD link to access it's nbors
    class(link_t),       pointer :: my_nb_link => null()  ! pointer for corrections towards itself as nbor of internal rt task 
    class(link_t),       pointer :: link       => null()  ! pointer to self link in the task list 
    type (ray_list_t),   pointer :: my_rays
    ! common items
    integer :: n_ray_dirs = 3
    integer :: geom = 1                                                 !cartesian = 1, spherical = 2, cylindrical = 3
    integer :: rt_freq = 1                                              ! how often rt tasks are done wrt experiment
    integer :: bundle_ID, dir_ID, rt_type, li(3), ui(3)
    integer :: ibc = 0, obc = 0, ibct = 3, obct = 0                     ! inner/outer boundary condition, treatment 
    real(kind=8) :: cell_ds(3)
    real(8):: offset(3)                 ! this is used to compute indices, locations
  contains
    procedure :: find_upstream                                             
    procedure :: init => init_rt
    procedure :: update => update_rt
    procedure :: init_obc
    procedure :: absolute_point_position
    procedure :: ray2cell
    procedure :: cell2pos
    procedure :: cell2ray
    procedure :: pos2cell
    procedure :: pos2ray
    procedure :: ray2pos
    procedure :: the_right_wall
    procedure :: update_time
  end type
!-----------------------------------------------------------------------
  CONTAINS

SUBROUTINE init_rt (self)
  class(rt_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: rt_t%init_rt was called, but it is an empty subroutine"
END SUBROUTINE init_rt

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

FUNCTION the_right_wall(self,me, nbor, wall) result (answer)
  class(rt_t) :: self
  type (link_t),        pointer :: me, nbor
  class (experiment_t), pointer :: mhd, nb_mhd
  integer :: wall(3)
  real (kind=8) :: d(3)
  logical :: answer
!-----------------------------------------------------------------------
  answer = .false.
  associate (my_type => me%task)
  SELECT type (my_type)
  type is (experiment_t)
    mhd => my_type
  class is (rt_t)
    mhd => my_type%mhd_patch
  END SELECT
  end associate
  
  associate (nbor_type => nbor%task)
  SELECT type (nbor_type)
  type is (experiment_t)
    nb_mhd => nbor_type
  class is (rt_t)
    nb_mhd => nbor_type%mhd_patch                                       ! if I got an RT task to test, switch to corresponding MHD patch
  end SELECT
  end associate
    
  if (abs(nb_mhd%level-mhd%level)<2) then                               ! connect only +/- 1 or same level
    if (wall(1) > 0) then
      if (mhd%position(1) < nb_mhd%position(1)) then
        if (mhd%position(1)+mhd%size(1)*0.5d0 > &
        nb_mhd%position(1) - nb_mhd%size(1)*0.5d0 - nb_mhd%ds(1)) answer = .true.
      end if
    else if (wall(1) < 0) then
      if (mhd%position(1) > nb_mhd%position(1)) then
        if (mhd%position(1)-mhd%size(1)*0.5d0 < &
        nb_mhd%position(1) + nb_mhd%size(1)*0.5d0 + nb_mhd%ds(1)) answer = .true.
      end if
    end if
    
    if (wall(2) > 0) then
      if (mhd%position(2) < nb_mhd%position(2)) then
        if (mhd%position(2)+mhd%size(2)*0.5d0 > &
        nb_mhd%position(2) - nb_mhd%size(2)*0.5d0 - nb_mhd%ds(2)) answer = .true.
      end if
    else if (wall(2) < 0) then
      if (mhd%position(2) > nb_mhd%position(2)) then
        if (mhd%position(2)-mhd%size(2)*0.5d0 < &
        nb_mhd%position(2) + nb_mhd%size(2)*0.5d0 + nb_mhd%ds(2)) answer = .true.
      end if
    end if
    
    if (wall(3) > 0) then
      if (mhd%position(3) < nb_mhd%position(3)) then
        if (mhd%position(3)+mhd%size(3)*0.5d0 > &
        nb_mhd%position(3) - nb_mhd%size(3)*0.5d0 - nb_mhd%ds(3)) answer = .true.
      end if
    else if (wall(3) < 0) then
      if (mhd%position(3) > nb_mhd%position(3)) then
        if (mhd%position(3)-mhd%size(3)*0.5d0 < &
        nb_mhd%position(3) + nb_mhd%size(3)*0.5d0 + nb_mhd%ds(3)) answer = .true.
      end if
    end if
    
  end if 
END FUNCTION the_right_wall

SUBROUTINE find_upstream (self, nbor_l, found)
  class (rt_t)  :: self
  class(link_t), pointer :: my_mhd_l
  class(link_t), pointer :: nbor, nbor_l, new
  class(experiment_t), pointer :: my_task, nbor_task
  integer :: found
! ----------------------------------------------------------------------
  call trace_begin('rt_t%find_upstream')
  found = 0
  allocate (nbor_l)
  nbor => self%mhd_link%nbor
  my_task => self%mhd_patch
  
  do while (associated(nbor))
    associate (nbtask => nbor%task)
    SELECT type (nbtask)
    class is (experiment_t)
      nbor_task => nbtask
      if (.not.nbor_is_my_child(my_task, nbor_task)) then
        if (self%the_right_wall(self%mhd_link, nbor, self%my_rays%lwall)) then
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
  
  nbor_l => nbor_l%next
  call trace_end
END SUBROUTINE find_upstream

SUBROUTINE init_obc (self)
  class(rt_t) :: self
!-----------------------------------------------------------------------
  call self%set(bits%has_BC)
  self%obc = 1
  self%my_nb_link%check = .true.
END SUBROUTINE init_obc

SUBROUTINE check_common_surface(me, nbor, wall, right, it)
  type (experiment_t), pointer :: me, nbor
  integer :: wall, w, s1, s2, it
  logical :: right, ok
!-----------------------------------------------------------------------
  SELECT case(it)
    case(1)
      w = 1; s1 = 2; s2 = 3
    case(2)
      w = 2; s1 = 1; s2 = 3
    case(3)
      w = 3; s1 = 1; s2 = 2
  end SELECT
  ok = .false.
  if (wall < 0) then
    if ((nbor%position(w)-0.5_8*nbor%size(w)) > &           
          (me%position(w)+0.5_8*  me%size(w)+me%ds(w))) ok = .true.
  else if (wall > 0) then
    if ((nbor%position(w)+0.5_8*nbor%size(w)) > &           
          (me%position(w)-0.5_8*  me%size(w)-me%ds(w))) ok = .true.
  else
    right = .true.
  end if
        
  if (ok) then
    if ((((nbor%position(s1)+0.5_8*nbor%size(s1)) >= &
            (me%position(s1)-0.5_8*  me%size(s1)))  .and. &
         ((nbor%position(s1)-0.5_8*nbor%size(s1)) <= &
            (me%position(s1)+0.5_8*  me%size(s1)))) .and. &
                
        (((nbor%position(s2)+0.5_8*nbor%size(s2)) >= &
            (me%position(s2)-0.5_8*  me%size(s2)))  .and. &
         ((nbor%position(s2)-0.5_8*nbor%size(s2)) <= &
            (me%position(s2)+0.5_8*  me%size(s2))))) right = .true.
  end if

END SUBROUTINE check_common_surface

FUNCTION nbor_is_my_child(me, nbor) result (out)
  type (experiment_t), pointer :: me, nbor
  !real(kind=8) :: me_li, me_ui, nbor_li, nbor_ui
  logical :: out
!-----------------------------------------------------------------------
  out = .false.
  if ((all(me%position + me%size*0.5_8 > nbor%position + nbor%size*0.5_8)) &
.and. (all(me%position - me%size*0.5_8 < nbor%position - nbor%size*0.5_8))) &
  out = .true.

END FUNCTION nbor_is_my_child


FUNCTION absolute_point_position(self, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  if (present(bn)) then
    cp = self%my_rays%fwd_coords(:,ray)
  else
    cp = self%my_rays%cell_coords(:,pt,ray)
  end if
  cell_pos = self%position-0.5*self%size+(cp-self%li)*self%cell_ds
END FUNCTION absolute_point_position

FUNCTION ray2cell(self, ray, pt) result (cell)
  class (rt_t) :: self
  integer :: ray, pt, cell(3)
!-----------------------------------------------------------------------
  cell = self%my_rays%cell_coords(:,pt,ray)
  
END FUNCTION ray2cell

FUNCTION cell2pos (self, cell) result(p)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  p = self%position - 0.5_8*self%size + (cell-self%li)*self%cell_ds

END FUNCTION cell2pos

FUNCTION pos2cell (self, p) result (cell)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  cell = nint((p-self%position+0.5_8*self%size)/self%cell_ds)+self%li
  
END FUNCTION pos2cell

! orient = 1 - lwall, orient = 2 - uwall
SUBROUTINE cell2ray (self, c, w, ray, pt, found)
  class(rt_t) :: self
  type(rt_cell_t), pointer :: cell
  integer :: c(3), ray, pt, w(3), orient, wall(3)
  logical :: found
!-----------------------------------------------------------------------
  cell => self%my_rays%cell_domain(c(1),c(2),c(3))%head
  wall = cell%lwall
  found = .false.
  do while(associated(cell))
    if(.not.found) then
      if (all(wall == w)) then
        ray = cell%ray
        pt = cell%point
        found = .true.
        return
      else
        cell => cell%next
      end if
    end if
  end do
END SUBROUTINE cell2ray

SUBROUTINE pos2ray (self, p, w, ray, pt, found)
  class(rt_t) :: self
  real(kind=8) :: p(3)
  integer :: w(3), ray, pt
  logical :: found
!-----------------------------------------------------------------------
  call self%cell2ray(self%pos2cell(p), w, ray, pt, found)
  
END SUBROUTINE pos2ray

FUNCTION ray2pos (self, ray, pt) result (p)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer :: ray, pt
!-----------------------------------------------------------------------
  p = self%cell2pos(self%ray2cell(ray,pt))
  
END FUNCTION ray2pos



FUNCTION point_in_cell (pt, cell, ds) result (out)
  real(kind=8) :: pt(3), cell(3), li(3), ui(3), ds(3)
  logical :: out
!-----------------------------------------------------------------------
  out = .false.
  li = cell - 0.5*ds
  ui = cell + 0.5*ds
  if ((all(li < pt)).and.(all(ui > pt))) out = .true.
END FUNCTION point_in_cell

SUBROUTINE update_time(self)
  class(rt_t) :: self
!-----------------------------------------------------------------------  
  if (self%time < 0) then
    self%time  = self%mhd_patch%time + self%dtime
    self%dtime = self%mhd_patch%dtime * self%rt_freq
  else   
    self%dtime = self%mhd_patch%dtime * self%rt_freq
    self%time  = self%time + self%dtime
  end if

END SUBROUTINE update_time


END MODULE rt_mod
