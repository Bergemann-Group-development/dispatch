!===============================================================================
!> Create a nested set of patches, in a hierarchy of Rubik's cubes, where the
!> innermost patch has a given size, and there are n_levels-1 sets of larger
!> Rubik's cubes, each three times larger than the previous, with the previous
!> sets of patches inside the center patch.
!===============================================================================
MODULE rubiks_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE timer_mod
  USE task_mod
  USE patch_mod
  USE task_list_mod
  USE experiment_mod
  USE units_mod
  USE scaling_mod
  USE shared_mod
  USE link_mod
  implicit none
  private
  type, public:: rubiks_t
    real(8):: size(3) = [1.0, 1.0, 1.0]
    real(8):: position(3) = [0.0, 0.0, 0.0]
    integer:: n_levels, n_repeat(3)
    integer:: n_patch(3), m_zoom(3)
    class(task_list_t), pointer:: task_list
  contains
    procedure init_hierarchy
    procedure init_level
  end type
  integer, save:: verbose=0
  real(8), save:: box(3), origin(3)
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!===============================================================================
SUBROUTINE init_hierarchy (self, label)
  class(rubiks_t):: self
  character(len=*), optional:: label
  !.............................................................................
  class(experiment_t), pointer:: patch, parent, new_parent
  class(task_t), pointer:: task
  integer, save:: n_patch(3)=0, m_zoom(3)=3, n_levels=5
  integer, save:: n_repeat(3)=0
  integer:: iostat
  real(8), save:: position(3), size(3), f_zoom=1d0
  integer:: i, j, k, irep, jrep, krep, level, levels
  namelist /rubiks_params/ size, position, n_levels, n_patch, m_zoom, n_repeat, &
                           f_zoom, verbose
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_t%init')
  size = self%size
  position = self%position
  n_levels = self%n_levels
  rewind (io%input)
  read (io%input, rubiks_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, rubiks_params)
  end if
  write (io%namelist, rubiks_params)
  !-----------------------------------------------------------------------------
  ! Avoid having to edit dozens of existing .nml files, by interpreting a single
  ! n_repeat value as being intended for the y-direction.
  !-----------------------------------------------------------------------------
  if (n_repeat(1)>0 .and. n_repeat(2)==0 .and. n_repeat(3)==0) then
    n_repeat(2) = n_repeat(1)
    n_repeat(1) = 1
    n_repeat(3) = 1
  end if
  !-----------------------------------------------------------------------------
  ! Similarly for n_patch, but assume that a single value applies to all
  ! directions.
  !-----------------------------------------------------------------------------
  if (n_patch(1)>0 .and. n_patch(2)==0 .and. n_patch(3)==0) then
    n_patch(2) = n_patch(1)
    n_patch(3) = n_patch(1)
  end if
  !-----------------------------------------------------------------------------
  ! Similarly for m_zoom, assuming that a single value applies to all
  ! directions.
  !-----------------------------------------------------------------------------
  if (m_zoom(1)>0 .and. m_zoom(2)==0 .and. m_zoom(3)==0) then
    m_zoom(2) = m_zoom(1)
    m_zoom(3) = m_zoom(1)
  end if
  box = n_patch*size*n_repeat*m_zoom**(n_levels-1)
  origin = -box/2d0
  shared%refine_ratio = m_zoom(1)
  shared%n_patch_total = product(n_patch)*(product(n_repeat)+n_levels-1)
  self%n_levels = n_levels
  self%position = position
  self%size = size
  self%n_patch = n_patch
  self%m_zoom = m_zoom
  shared%f_zoom = m_zoom*f_zoom
  !-----------------------------------------------------------------------------
  ! Save the outer size of the hierarchy as the task list size, for distance
  ! measurements.  The self%task_list%size is the full size of the objects in
  ! the task list, and the size of the period in the periodic wrapping.  The
  ! self%size, on the other hand, is the size of the smallest, central patch.
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  call self%task_list%init (label)
  self%task_list%size = self%size*self%n_patch*shared%f_zoom**(n_levels-1)*n_repeat
  n_repeat = (n_repeat/2)*2 + 1
  shared%n_repeat = n_repeat
  level = 1
  ! repeated grids at the coarsest resolution
  do krep=0,n_repeat(3)-1
    do jrep=0,n_repeat(2)-1
      do irep=0,n_repeat(1)-1
        position = (([irep,jrep,krep]-n_repeat/2)*self%task_list%size)/n_repeat &
                 + self%position
        if (verbose > 0) &
          write (stdout,'(a,2(i1,"x"),i1,a,1p,3g14.5)') &
            ' generating a cube of ', n_patch, &
            ' patches centred at position', position
        call self%init_level (level, position)
        if (all(position == self%position)) then
        end if
      end do
    end do
  end do
  ! succesively refined hierarchies
  position = self%position
  if (verbose > 0) &
    write (stdout,'(a,2(i1,"x"),i1,a,1p,3g14.5)') &
     ' generating a hiearchy of ', n_patch, &
     ' patches centered at position', position
  do level=2,n_levels
    if (verbose > 0) &
      write (stdout,'(2x,a,i3)') 'level', level
    call self%init_level (level, position)
    !call timer%print ()
  end do
  !-----------------------------------------------------------------------------
  ! Initialize the task list, setting bits for internal, boundary, virtual and
  ! external patches, counting these, setting physical boundary bits, and 
  ! initialize nbor lists based on this.
  !-----------------------------------------------------------------------------
  call self%task_list%reset_status
  call self%task_list%count_status
  io%nwrite = self%task_list%na
  io%ntotal = shared%n_patch_total
  call self%task_list%init_levels
  call self%task_list%find_parents
  call self%task_list%init_bdries
  call self%task_list%init_all_nbors
  call self%task_list%print_tasks
  call trace%end()
END SUBROUTINE init_hierarchy

!===============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!===============================================================================
SUBROUTINE init_level (self, level, position)
  class(rubiks_t):: self
  integer:: level
  real(8):: position(3)
  integer, save:: itimer=0
  !.............................................................................
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k, ndrop, task_level
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_revenge_t%init_level', itimer=itimer)
  ndrop = 0
  do k=1,self%n_patch(3)
  do j=1,self%n_patch(2)
  do i=1,self%n_patch(1)
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    allocate (experiment_t:: task)
    select type (task)
    class is (experiment_t)
    patch => task
    end select
    !---------------------------------------------------------------------------
    ! The patch size should equal the self%size when level=n_levels, since that
    ! is the level of the smallest, central patch.  The size of the largest
    ! patches is f_zoom**(n_levels-1) larger.
    !---------------------------------------------------------------------------
    patch%size = self%size*real(shared%f_zoom,kind=8)**(self%n_levels-level)
    patch%position = ([i,j,k]-(1.+self%n_patch)/2.0)*patch%size + position
    !---------------------------------------------------------------------------
    ! The patch%box value is used in periodic wrap-around, so this implements
    ! periodic conditions along the x-direction (= phi-direction in cylindrical
    ! coordinates).
    !---------------------------------------------------------------------------
    patch%level = level
    if (patch%is_relevant()) then
      patch%box = box
      patch%origin = origin
      call patch%init
      call patch%init_level
      !-------------------------------------------------------------------------
      ! Set bits
      !-------------------------------------------------------------------------
      call patch%set (bits%static)
      if (k==1.or.k==self%n_patch(3)) call patch%set (bits%top_bottom)
      if (verbose > 2) &
        write (stdout,'(a,i6,i4,2(2x,3g14.5),3x,6l1)') &
          'rubiks_mod%init_level: id, level, size, pos =', &
          patch%id, patch%level, patch%size, patch%position, &
          patch%boundaries%is_set (bits%xl), &
          patch%boundaries%is_set (bits%yl), &
          patch%boundaries%is_set (bits%zl), &
          patch%boundaries%is_set (bits%xu), &
          patch%boundaries%is_set (bits%yu), &
          patch%boundaries%is_set (bits%zu)
      !-------------------------------------------------------------------------
      ! Add to a list of links, with pointers from patch to patch
      !-------------------------------------------------------------------------
      call self%task_list%add (task)
      !-------------------------------------------------------------------------
      ! If at edge, set boundary bits
      !-------------------------------------------------------------------------
      if (level==1) then
        if (.not.patch%periodic(1).and.i==1)               call patch%boundaries%set (bits%xl)
        if (.not.patch%periodic(2).and.j==1)               call patch%boundaries%set (bits%yl)
        if (.not.patch%periodic(3).and.k==1)               call patch%boundaries%set (bits%zl)
        if (.not.patch%periodic(1).and.i==self%n_patch(1)) call patch%boundaries%set (bits%xu)
        if (.not.patch%periodic(2).and.j==self%n_patch(2)) call patch%boundaries%set (bits%yu)
        if (.not.patch%periodic(3).and.k==self%n_patch(3)) call patch%boundaries%set (bits%zu)
        if (patch%boundaries%bits /= 0 .and. verbose == 2) &
          write (stdout,'(a,i6,5x,a,3x,6l1)') &
            ' id:', patch%id, 'boundary bits:', &
            patch%boundaries%is_set (bits%xl), &
            patch%boundaries%is_set (bits%yl), &
            patch%boundaries%is_set (bits%zl), &
            patch%boundaries%is_set (bits%xu), &
            patch%boundaries%is_set (bits%yu), &
            patch%boundaries%is_set (bits%zu)
      end if
    else
      ndrop = ndrop+1
      if (verbose > 0) write (stdout,'(a,1p,3e11.3,a)') &
        'patch at position',patch%position,' not relevant -- dropping'
      deallocate (task)
    end if
  end do
  end do
  end do
  if (ndrop > 0) &
    write (stdout,'("dropped",i4," not relevant patches on level",i3)') &
      ndrop, level
  call trace_end (itimer)
END SUBROUTINE init_level

END MODULE rubiks_mod
