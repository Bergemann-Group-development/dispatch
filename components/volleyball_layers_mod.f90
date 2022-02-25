!===============================================================================
!> Construct a spherical layer, with volleyball distribution of patches
!===============================================================================
MODULE volleyball_layers_mod
  USE sphere_point_mod
  USE experiment_mod
  USE task_mod
  USE mesh_mod
  USE dll_mod
  USE io_mod
  USE trace_mod
  USE link_mod
  USE halo_mod
  USE task_list_mod
  USE timer_mod
  USE omp_mod
  USE scaling_mod   ! FIXME: remove later
  USE units_mod
  implicit none
  private
  !.............................................................................
  type, public:: volleyball_layer_t
    integer                            :: number=0
    integer                            :: ntheta=10
    real(8)                            :: theta0, dtheta, radius
    real(8)                            :: r_min, r_max
    type(dll_t)                        :: list
    class(link_t),             pointer :: link => null()
    class(volleyball_layer_t), pointer :: prev => null()
    class(volleyball_layer_t), pointer :: next => null()
  contains
    procedure:: init => init_layer
  end type
  !.............................................................................
  type, public:: volleyball_layers_t
    integer:: n_wide
    real(8):: r_min, r_max
    type(dll_t):: list
    type(volleyball_layer_t), pointer :: head => null()  
  contains
    procedure:: init
    procedure:: init_all_nbors
    procedure:: init_nbors
  end type
  !.............................................................................
  integer:: verbose=0
  integer:: n_wide=8

  logical:: detailed_timer=.false.
  real(8):: r_min=0.6
  real(8):: r_max=0.95
CONTAINS

!===============================================================================
!> Given the number of intervals in the theta direction, generate the patch
!> positions, saving the point information in a dll_t data type
!===============================================================================
SUBROUTINE init (self)
  class(volleyball_layers_t):: self
  !.............................................................................
  type(sphere_point_t), pointer:: point, save, tmp
  class(dll_node_t), pointer:: dll_node, layer_node
  type(volleyball_layer_t), pointer:: layer
  type(volleyball_layer_t), pointer:: prev, head
  integer:: iostat, i, n_layer
  real(8):: f, radius
  namelist /vb_layers_params/ verbose, n_wide, r_min, r_max, detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layers_t%init')
  !$omp critical (input_cr)
  rewind (io_unit%input)
  read (io_unit%input, vb_layers_params, iostat=iostat)
  write (io%namelist, vb_layers_params)
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Construct nlayer layers and append to self%list, arranging also for a
  ! singly-linked list of layers
  !-----------------------------------------------------------------------------
  call self%list%init
  self%n_wide = n_wide
  f = 1d0-atan(sqrt(0.5d0))/n_wide/2d0
  radius = r_max*f
  self%r_max = radius
  n_layer=0
  nullify (head, prev)
  do while (radius*f > r_min)
    allocate (layer, dll_node)
    if (.not. associated(self%head)) then
      self%head => layer
    end if
    n_layer=n_layer+1
    layer%radius = radius
    layer%number = n_layer
    call layer%init (radius)
    dll_node%car => layer
    call self%list%append (dll_node)
    if (associated(prev)) then
      prev%next => layer
      layer%prev => prev
    end if
    prev => layer
  end do
  self%r_min = radius*f
  write (stdlog,*) 'volleyball_layers: number of patches =', self%list%n
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Given the number of intervals in the theta direction, generate the patch
!> positions, saving the point information in a dll_t data type
!===============================================================================
SUBROUTINE init_layer (self, radius)
  class(volleyball_layer_t):: self
  real(8):: radius, size, step
  !.............................................................................
  type(sphere_point_t), pointer:: point, save
  class(dll_node_t), pointer:: dll_node
  integer:: i, j
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layer_t%init_layer')
  call self%list%init
  !-----------------------------------------------------------------------------
  ! Construct a set of patch positions that define a volleyball front face
  !-----------------------------------------------------------------------------
  self%radius = radius
  self%ntheta = n_wide
  self%theta0 = atan(sqrt(0.5d0))
  self%dtheta = self%theta0/self%ntheta
  write (stdlog,*) 'volleyball_layer: theta0, dtheta =', self%theta0, self%dtheta
  !-----------------------------------------------------------------------------
  ! The point%init sets the coordinates to the mid-point of the front face
  !-----------------------------------------------------------------------------
  allocate (point)
  call point%init
  point%r       = radius
  point%size    = radius*self%dtheta/(1d0-self%dtheta)
  point%verbose = verbose
  point%r       = point%r - point%size/2d0
  if (verbose > 0) &
    call point%info ('init')
  allocate (save, source=point)
  point%verbose = verbose
  !---------------------------------------------------------------------------
  ! Add half a patch size north-east offset
  !---------------------------------------------------------------------------
  call point%go_north (0.5d0*self%dtheta)
  call point%go_east (0.5d0*self%dtheta)
  call point%info ('row start')
  !-----------------------------------------------------------------------------
  ! March in the theta direction
  !-----------------------------------------------------------------------------
  do while (point%is_inside (self%theta0, self%dtheta))
    write (stdlog,*) 'volleyball_layer: size, r =', point%size, point%r
    !---------------------------------------------------------------------------
    ! March in the phi direction while patches are inside the front section,
    ! while appending a point and its three mirror images to the list of points
    !---------------------------------------------------------------------------
    save = point
    do while (point%is_inside (self%theta0, self%dtheta))
      call symmetries
      call point%east_west              ! swap East-West
      call symmetries
      call point%north_south            ! swap Upper-Lower
      call symmetries
      call point%east_west
      call symmetries
      call point%north_south
      !--- Back to start -----
      call point%go_east (self%dtheta)
      call point%info ('east')
    end do
    !--- Finished with longitude sweep
    call point%info ('east exit')
    !--- Back to starting point
    save%id = point%id
    point = save
    !--- Next row
    call point%go_north (self%dtheta)
    call point%info ('north')
    save = point
  end do
  call point%info ('north exit')
  deallocate (save)
  size = point%size
  print '(a,i4,2f8.3)', &
    ' layer,radius,size = ', self%number, radius, size
  !-----------------------------------------------------------------------------
  ! The step in radius needs to be the average of the current size and the
  ! size in the next layer -- which is smaller by a factor (1-step/R); i..e
  ! step = old_step*(1.+(1.-old_step/R))/2=old_step*(1.-old_step/R/2)
  !-----------------------------------------------------------------------------
  step = size*(1d0-self%dtheta)*(1d0-self%dtheta**2)
  if (verbose>0) write (stdout,'(a,3f8.2,i6)') 'radius, size, step, n_task =', &
    radius*cgs%r_sun/scaling%l, size*cgs%r_sun/scaling%l, &
    step*cgs%r_sun/scaling%l, self%list%n
  radius = radius - step
  call trace%end
contains
  !=============================================================================
  ! Append a (copy of a) point to the list of points
  !=============================================================================
  subroutine append
    class(sphere_point_t), pointer:: tmp
    integer, save:: id=0
    allocate (tmp, source=point)
    point%id = point%id+1
    call tmp%unit_vectors()
    call tmp%info('append',verbose=1)
    if (verbose>0) write (stdout,'(i6,2(2x,3f8.4))') &
      point%id, point%phi, point%theta, point%r, point%p
    allocate (dll_node)
    dll_node%car => tmp
    call self%list%append (dll_node)
  end subroutine
  !=============================================================================
  subroutine symmetries
    !--- loop over the front and back
    do i=1,2
      !--- loop over the three front faces: East, North, Up
      do j=1,3
        call append
        call point%permute ()             ! East -> Up -> North
      end do
      call point%reflect ()               ! front -> back
    end do
  end subroutine symmetries
END SUBROUTINE init_layer

!===============================================================================
! Construct all nbor lists
!===============================================================================
SUBROUTINE init_all_nbors (self, task_list)
  class(volleyball_layers_t)         :: self
  class(task_list_t),         pointer :: task_list
  !.............................................................................
  class(volleyball_layer_t), pointer :: layer, prev, next
  class(link_t),             pointer :: link, first, nextl
  class(dll_node_t),         pointer :: node
  integer                            :: i, i0, i1, i2, n, m
  integer, save                      :: itimer=0
  !.............................................................................
  call trace%begin('vb_layers_t%init_all_nbors', itimer=itimer)
  call timer%tic ()
  !call init_nonleaf (task_list)                           ! set non_leaf bit?
  call task_list%reset_status                             ! update task bits
  call task_list%count_status                             ! update task counts
  i2 = task_list%n/omp%nthreads + 1                       ! tasks per thread
  write (stdout,*) task_list%n, 'tasks in task_list'
  n = 0
  !$omp parallel private(i,i0,i1,link,first,layer,node, nextl,m) &
  !$omp          shared(self,i2,n,stdout) default(none)
  i0 = omp%thread*i2
  i1 = i0 + i2
  i = 0                                                   ! task counter
  layer => self%head                                      ! start checking links
  do while (associated(layer))
    node => layer%list%head                               ! loop over layer
    do while (associated(node))                           ! process layer
      link => to_link (node)                              ! to task link
      i = i+1                                             ! increment counter
      if (i >= i0 .and. i<i1) then                        ! thread interval
        call self%init_nbors (layer, link)
        !$omp atomic capture
        n = n+1
        m = n
        !$omp end atomic
        if (mod(m,2000) == 0) &
          write (stdout,*) m,'nbor lists'
        if (associated(layer%next)) then
          first => layer%next%link
          call self%init_nbors (layer%next, link)
        end if
        if (associated(layer%prev)) then
          first => layer%prev%link
          call self%init_nbors (layer%prev, link)
        end if
      end if
      node => node%next
    end do
    first => layer%link                                   ! previous layer
    layer => layer%next                                   ! continue list
  end do
  !$omp end parallel
  !call self%reset_status                                ! set bits and count
  call timer%toc ('time to initialize nbor lists')
  call trace%end (itimer)
contains
END SUBROUTINE init_all_nbors

!===============================================================================
function to_link (node) result (link)
  class(dll_node_t), pointer:: node
  class(link_t),     pointer:: link
  class(*),          pointer:: ptr
  !-----------------------------------------------------------------------------
  ptr => node%car
  select type (ptr)
  class is (link_t)
    link => ptr
  end select
end function to_link

!===============================================================================
! Construct the nbor lists for one task
!===============================================================================
SUBROUTINE init_nbors (self, layer, link)
  class(volleyball_layers_t)         :: self
  class(volleyball_layer_t), pointer :: layer
  class(link_t),             pointer :: link
  !.............................................................................
  class(task_t),             pointer :: task, nbtask
  class(experiment_t),       pointer :: task1, task2
  real(8)                            :: limit
  class(dll_node_t),         pointer :: node
  class(link_t),             pointer :: nblink
  integer                            :: id, n
  integer, save                      :: itimer=0, nprint=10
  !-----------------------------------------------------------------------------
  ! Loop over tasks in that and subsequent layers, stopping below radius_limit
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_layers_t%init_nbors', itimer=itimer)
  id = link%task%id
  task => link%task
  select type (task)
  class is (experiment_t)
    task1 => task
    if (verbose > 0) &
      write (stdout,'("init_nbors:",2i6,5(2x,3f8.4))') &
        task%id, task%point%id, task%point%phi, task%point%theta, &
        task%point%r, task%point%p, task%point%e
    !---------------------------------------------------------------------------
    ! Loop over one layer
    !---------------------------------------------------------------------------
    n = 0
    node => layer%list%head                               ! loop over layer
    do while (associated(node))                           ! process layer
      nblink => to_link (node)                            ! to task link
      nbtask => nblink%task
      select type (nbtask)
      class is (experiment_t)
        task2 => nbtask
        limit = nbtask%size(1)*1.5                        ! 1.5 patch size
        if (abs(nbtask%radius-task%radius) <= limit) then
          if (nbtask%id /= task%id) then
            if (task1%overlaps2 (task2) .or. task2%overlaps2 (task1)) then
              n = n+1
              call link%add_nbor (nblink)
              if (verbose > 1) then
                write (stdout,'("init_nbors:",6x,2i6,5(2x,3f8.4))') nbtask%id, &
                  nbtask%point%id, nbtask%point%phi, nbtask%point%theta, &
                  nbtask%point%r, nbtask%point%p, nbtask%point%e
              end if
            end if
          end if
        end if
      end select
      node => node%next
    end do
    if (verbose > 0 .and. nprint > 0) then
      nprint = nprint-1
      write (stdout,'(a,2i6)') &
        'volleyball_layers_t%init_nbors: id,n=', task%id, n
  end if
  end select
  call trace%end (itimer)
contains

!===============================================================================
function overlaps0 (task1, task2) result (out)
  class(experiment_t) :: task1, task2
  logical             :: out
  real(8)             :: p1(3), p2(3), d(3), limit(3), s
  logical             :: debug=.false.
  integer,       save :: itimer=0
  !.............................................................................
  if (detailed_timer) &
    call trace%begin ('overlaps', itimer=itimer)
  call timer%tic ()
  limit = (task1%size + task2%size + task1%ds + task2%ds)*0.75_8
  select type (task1)
  class is (experiment_t)
    d = task1%position - task2%position
    s = sqrt(sum(d**2))
    out = all (abs(d) < limit)
    if (verbose > 1) then
      if (all(abs(d) < 2.5*limit)) then
        d = d/limit
        write (stdout,1) task1%id, task1%size, task1%position, &
                        task2%id, task2%size, task2%position, &
                        d, limit, out
      1 format("overlaps: ",2(i6,2(2x,3f8.4)),2(2x,3f8.4),l3)
        flush (stdout)
      end if
    end if
  end select
  call timer%toc ('overlaps1')
  if (detailed_timer) &
    call trace%end (itimer)
end function overlaps0

!===============================================================================
!> Exact overlap determination
!===============================================================================
function overlaps1 (task1, task2) result (out)
  class(experiment_t), pointer :: task1, task2
  class(mesh_t),  pointer :: m
  logical                 :: out
  real(8), dimension(3)   :: p1, s1, p2, s2, llc, urc, corner
  real(8)                 :: erot(3,3)
  logical                 :: debug=.false.
  integer                 :: i1, i2, i3, i
  !-----------------------------------------------------------------------------
  ! Rotation matrix transforming from task2 to task1 local coordinates
  !-----------------------------------------------------------------------------
  erot = matmul(transpose(task1%erot), task2%erot)
  !-----------------------------------------------------------------------------
  ! Task1 and task2 positions and half sizes in local coordinates
  !-----------------------------------------------------------------------------
  p1 = matmul(task1%position, task1%erot)
  p2 = matmul(task2%position, task2%erot)
  s1 = task1%size*(task1%n+2*task1%ng)/(2d0*task1%n)
  s2 = task2%size/2d0
  !-----------------------------------------------------------------------------
  ! Task1 LLC and URC corner
  !-----------------------------------------------------------------------------
  llc = p1 - s1
  urc = p1 + s1
  out = .false.
  !-----------------------------------------------------------------------------
  ! Loop over and check corners
  !-----------------------------------------------------------------------------
  do i3=-1,1,2
  do i2=-1,1,2
  do i1=-1,1,2
    !---------------------------------------------------------------------------
    ! Transform task2 coordinates to task1 local, and check corners
    !---------------------------------------------------------------------------
    corner = p2 + [i1,i2,i3]*s2
    corner = matmul (task2%erot, corner)
    corner = matmul (corner, task1%erot)
    out = out .or. all(corner > llc .and. corner < urc)
    if (out) then
      if (task1%id == io%id_debug) then
        p2 = matmul(task2%position, task1%erot)
        PRINT '("overlaps1",2i6,6(1x,3f8.4))', &
          task1%id, task2%id, task1%position, task2%position, p1, p2, s1, s2
      end if
      return
    end if
  end do
  end do
  end do
end function overlaps1

END SUBROUTINE init_nbors

!===============================================================================
FUNCTION to_radius (link) result (out)
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  real(8):: out
  !.............................................................................
  if (associated(link)) then
    task => link%task
    select type (task)
    class is (experiment_t)
      out = task%radius
    end select
  else
    out = huge(1.)
  end if
END FUNCTION to_radius

END MODULE volleyball_layers_mod
