!===============================================================================
!> A communications hub between patches. Three methods are used: download_same,
!> download_higher and download_lower
!===============================================================================
MODULE data_hub_mod
  USE io_mod
  USE omp_mod
  USE bits_mod
  USE kinds_mod
  USE mesh_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public:: data_hub_t
    class(patch_t), pointer:: target => null()                          ! pointer to a target patch
    class(patch_t), pointer:: source => null()                          ! pointer to a source patch
    class(data_hub_t), pointer :: head => null()
    class(data_hub_t), pointer :: next => null()
    class(data_hub_t), pointer :: prev => null()
    integer, dimension(:,:), pointer :: t_l => null(), t_u => null()
    integer, dimension(:,:), pointer :: s_l => null(), s_u => null()
    logical :: first_time = .true.                                      ! if first time, hib has to be initialized.
    logical :: rt_target  = .false.
    logical :: all_cells  = .false.
    procedure(dnload_from_same), pointer :: dnload => null()
    integer :: int_order = -1
    integer :: level
    integer :: ID
  contains
    procedure:: add_hub
    procedure:: check_in_nbors
    procedure:: check_in_target
    procedure:: dnload_from_higher
    procedure:: dnload_from_lower
    procedure:: dnload_from_same
    procedure:: init
    procedure:: init_higher
    procedure:: init_lower
    procedure:: init_same
    procedure:: interpolate
    procedure:: interpolate_unsigned
    procedure:: interpolate_one
    procedure:: interpolate_unsigned_one
    procedure:: remove_hub
    procedure:: select_dnloader
    procedure:: update
  end type
CONTAINS

!=======================================================================
!> Initialise the data communication hub.
!=======================================================================
SUBROUTINE init (self, link, all_cells)
  class(data_hub_t):: self
  class(link_t), pointer :: link
  logical, optional :: all_cells
  ! ....................................................................
  class(patch_t), pointer :: patch
  ! --------------------------------------------------------------------
  call trace%begin ('data_hub_t%init')
  if (self%first_time) then
    patch => task2patch(link%task)                                      ! convert to patch
    self%target => patch                                                ! save a pointer to target patch
    call self%select_dnloader(link, all_cells)                          ! map the dh-cell shell
    self%first_time = .false.                                           ! should not be repeated again.
  end if
  call trace%end()
END SUBROUTINE init

!=======================================================================
!> Select the download method for each nbor.
!=======================================================================
SUBROUTINE select_dnloader (self, link, all_cells)
  class(data_hub_t)      :: self
  class(link_t), pointer :: link
  logical, optional      :: all_cells
  ! ....................................................................
  class(link_t),     pointer :: nbor
  class(patch_t),    pointer :: source, patch
  class(data_hub_t), pointer :: hub
  ! --------------------------------------------------------------------
  patch => self%target
  nbor  => link%nbor
  do while (associated(nbor))
    source => task2patch(nbor%task)
    if (nbor%download) then
      if (trim(patch%kind)/='rt_solver') then
        if ((source%is_clear(bits%no_download))) then
          if (abs(source%level-patch%level) < 2) then
            call self%add_hub (patch, source, hub, all_cells)
            if (source%level == patch%level) then
              hub%dnload => dnload_from_same                                ! target and source are of the same level
            call hub%init_same
            else if (source%level > patch%level) then
              hub%dnload => dnload_from_higher                              ! source is of higher level than the target
              call hub%init_higher
            else
              hub%dnload => dnload_from_lower                               ! source is of lower level than the target
              call hub%init_lower
            end if
          end if
        end if
      else if (trim(source%kind)=='rt_solver') then
        if (abs(source%level-patch%level) < 2) then
          call self%add_hub (patch, source, hub, all_cells)
          hub%rt_target = .true.
          if (source%level == patch%level) then
            hub%dnload => dnload_from_same                                ! target and source are of the same level
            call hub%init_same
          else if (source%level > patch%level) then
            hub%dnload => dnload_from_higher                              ! source is of higher level than the target
            call hub%init_higher
          else
            hub%dnload => dnload_from_lower                               ! source is of lower level than the target
            call hub%init_lower
          end if
        end if
      end if
    end if
    nbor => nbor%next                                                   ! move to next nbor
  end do
END SUBROUTINE select_dnloader

!=======================================================================
!> Add a new data hub, responsible for the source in question. Data hubs
!> are arranged in the source level of refinement from highest to lowest,
!> as for levels, lower than the target, additional cells are needed for
!> interpolation.
!=======================================================================
SUBROUTINE add_hub (self, target, source, hub, all_cells)
  class(data_hub_t)          :: self
  class(patch_t),    pointer :: target, source
  class(data_hub_t), pointer :: hub
  logical, optional:: all_cells
  ! ....................................................................
  class (data_hub_t), pointer :: next, prev
  logical :: found
  ! --------------------------------------------------------------------
  allocate (hub)
  hub%target => target
  hub%source => source
  hub%level  = source%level
  hub%ID     = source%id
  if (present(all_cells)) hub%all_cells = all_cells
  if (.not.associated(self%head)) then
    self%head => hub
  else
    nullify (prev)
    found = .false.
    next => self%head
    do while (associated(next).and.(.not.found))
      if (.not.associated(next%source,hub%source)) then
        if (next%source%level > hub%source%level) then
          found = .true.
          hub%next => next
          if (associated(prev)) then
            prev%next => hub
          else
            self%head => hub
          end if
        end if
      end if
      prev => next
      next => next%next
    end do
    if (.not.found) then
      if (associated(prev)) then
        prev%next => hub
      else
        call io%abort("data_hub_t::add_hub:: location not found?")
      end if
    end if
  end if
END SUBROUTINE add_hub

!=======================================================================
!> Remove a particular data hub. If clear is set to be .true., all hubs
!> in the list are removed.
!=======================================================================
SUBROUTINE remove_hub (self, id, clear)
  class(data_hub_t) :: self
  integer, optional :: id
  logical, optional :: clear
  ! ....................................................................
  class(data_hub_t), pointer :: this, next
  logical:: cleared
  ! --------------------------------------------------------------------
  if (present(id)) then
    cleared = .false.
    if (associated(self%head)) then
      this => self%head
      do while(associated(this).and.(.not.cleared))
        if (this%id == id) then
          if (associated(this%prev)) then
            this%prev%next => this%next
            this%next%prev => this%prev
            deallocate(this%s_l, this%s_u)
            deallocate(this%t_l, this%t_u)
            deallocate(this)
            cleared = .true.
          else
            if (associated(this%next)) then
              self%head => this%next
            else
              self%first_time = .true.
            end if
            deallocate(this%s_l, this%s_u)
            deallocate(this%t_l, this%t_u)
            deallocate(this)
            cleared = .true.
          end if
        end if
        this => this%next
      end do
    end if
  end if
  if (present(clear)) then
    this => self%head
    do while(associated(this))
      next => this%next
      deallocate(this%s_l, this%s_u)
      deallocate(this%t_l, this%t_u)
      deallocate(this)
      this => next
    end do
  end if
END SUBROUTINE remove_hub

!=======================================================================
!> Initialize the method when both the target and the source are of the
!> same level of refinement. Approach is identical to the dnload_same in
!> download_mod.
!> If "no no-mans land" mode is used, m%n has to be changed into (m%n-1)
!=======================================================================
SUBROUTINE init_same (self)
  class(data_hub_t) :: self
  ! ....................................................................
  class(mesh_t), pointer ::m
  real(8) :: dist(3)
  integer :: i
  ! --------------------------------------------------------------------
  dist = self%source%distance(self%target)
  allocate (self%t_l(3,1))
  allocate (self%t_u(3,1))
  allocate (self%s_l(3,1))
  allocate (self%s_u(3,1))
  do i = 1,3
    m => self%target%mesh(i)
    if (dist(i) > self%source%ds(i)) then
      self%t_l(i,1) = m%uo
      self%t_u(i,1) = m%ub
      self%s_l(i,1) = self%t_l(i,1) - m%n
      self%s_u(i,1) = self%t_u(i,1) - m%n
    else if (dist(i) < -self%source%ds(i)) then
      self%t_l(i,1) = m%lb
      self%t_u(i,1) = m%lo
      self%s_l(i,1) = self%t_l(i,1) + m%n
      self%s_u(i,1) = self%t_u(i,1) + m%n
    else
      self%t_l(i,1) = m%li
      self%t_u(i,1) = m%ui
      self%s_l(i,1) = self%t_l(i,1)
      self%s_u(i,1) = self%t_u(i,1)
    end if
  end do
END SUBROUTINE init_same

!=======================================================================
!> Download guard cell values from the source, which has the same level
!> of refinement, as the target patch. Download range is pre-determined
!> by the initialization subroutine above.
!=======================================================================
SUBROUTINE dnload_from_same (self, only)
  class(data_hub_t) :: self
  integer, optional :: only
  ! ....................................................................
  integer :: jt(2), iv, sl(3), su(3), tl(3), tu(3), iv1, iv2
  real    :: pt(2)
  integer, save:: itimer=0
  ! --------------------------------------------------------------------
  call trace%begin ('data_hub_t%dnload_same',itimer=itimer)
  call self%source%time_interval (self%target%time, jt, pt)
  sl = self%s_l(:,1); su = self%s_u(:,1)
  tl = self%t_l(:,1); tu = self%t_u(:,1)
  if (self%rt_target) then
    self%target%mem(tl(1):tu(1),tl(2):tu(2),tl(3):tu(3),:,self%target%it,1)= &
         pt(1) *self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),:,jt(1),1) + &
         pt(2) *self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),:,jt(2),1)
  else
    if (present(only)) then
      iv1 = only
      iv2 = only
    else
      iv1 = 1
      iv2 = self%source%nv
    end if
    do iv=iv1, iv2
      if (self%source%unsigned(iv)) then
        self%target%mem(tl(1):tu(1),tl(2):tu(2),tl(3):tu(3),iv,self%target%it,1) = exp &
        (pt(1) *log(self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),iv,jt(1),1)) + &
         pt(2) *log(self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),iv,jt(2),1)))
      else
        self%target%mem(tl(1):tu(1),tl(2):tu(2),tl(3):tu(3),iv,self%target%it,1)= &
         pt(1) *self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),iv,jt(1),1) + &
         pt(2) *self%source%mem(sl(1):su(1),sl(2):su(2),sl(3):su(3),iv,jt(2),1)
      end if
    end do
  end if
  call trace%end(itimer)
END SUBROUTINE dnload_from_same

!=======================================================================
!> Initialize the method when the source is of a higher level of refinement
!> than the target patch.
!=======================================================================
SUBROUTINE init_higher (self)
  class(data_hub_t) :: self
  ! ....................................................................
  class(mesh_t),  pointer:: m, m2
  class(patch_t), pointer:: source, target
  real(8) :: pos(3)
  integer :: iv, i, j(3), nv
  real(8) :: eps, pl(3), pu(3), dp
  ! --------------------------------------------------------------------
  source => self%source
  target => self%target
  if (self%rt_target) then
    allocate (self%t_l(3,1))
    allocate (self%t_u(3,1))
    nv = 1
  else
    allocate (self%t_l(3,target%nv))
    allocate (self%t_u(3,target%nv))
    nv = target%nv
  end if
  eps = 0.001_8
  do iv = 1, nv
  ! --------------------------------------------------------------------
  ! First and last points inside the source domain
  ! --------------------------------------------------------------------
    do i = 1,3
      m => source%mesh(i)
      pl(i) = m%r(m%li) - eps
      pu(i) = m%r(m%ui) + eps
    end do
  !---------------------------------------------------------------------
  ! Map to the first and last point inside the source domain
  !---------------------------------------------------------------------
    self%t_l(:,iv) = target%index_only2 (source, pl, iv, roundup=.true.)
    self%t_u(:,iv) = target%index_only2 (source, pu, iv)
  ! --------------------------------------------------------------------
  ! Make sure the target range is legal
  ! --------------------------------------------------------------------
  if (.not.self%all_cells) then
    do i=1,3
      m => target%mesh(i)
      m2=> source%mesh(i)
      self%t_l(i,iv) = max(m%lb,min(m%uo,self%t_l(i,iv)))
      self%t_u(i,iv) = max(m%lo,min(m%ub,self%t_u(i,iv)))
      dp = round(m%p-m2%p,12)+m%h(iv)*(m%d-m2%d)
      if (self%t_u(i,iv)>m%lo) then
        if (m%r(m%li)+dp > m2%r(m2%ui)) &
        self%t_u(i,iv) = min(self%t_u(i,iv),m%lo)
      end if
      if (self%t_l(i,iv)<m%uo) then
        if (m%r(m%ui)+dp < m2%r(m2%li)) &
        self%t_l(i,iv) = max(self%t_l(i,iv),m%uo)
      end if
      if (self%t_u(i,iv)>m%ui) then
        if (m%r(m%uo)+dp > m2%uf) then
          self%t_u(i,iv) = min(self%t_u(i,iv),m%ui)
        else
          self%t_u(i,iv) = max(self%t_l(i,iv),m%ub)
        end if
      end if
      if (self%t_l(i,iv)<m%li) then
        if (m%r(m%lo)+dp < m2%lf) then
          self%t_l(i,iv) = max(self%t_l(i,iv),m%li)
        else
          self%t_l(i,iv) = min(self%t_l(i,iv),m%lb)
        end if
      end if
    end do
  end if
  end do
END SUBROUTINE init_higher

!=======================================================================
!> Download guard cell values from the source, which has a larger level
!> of refinement, than the target patch. A temporary scratch cube is
!> created and values the right entries in the source memory slots are
!> copied. Then an interpolation is done. Can be easily extended to
!> accommodate the van Leer interpolation scheme.
!=======================================================================
SUBROUTINE dnload_from_higher (self, only)
  class(data_hub_t) :: self
  integer, optional :: only
  ! ....................................................................
  integer :: jt(2), iv, tl(3), tu(3), i, j, k, iv1, iv2
  integer :: d, l(3), u(3), ii
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: cube
  real(kind=KindScalarVar) :: tcell
  integer, save:: itimer=0
  real    :: pt(2), p(3)
  class(patch_t), pointer :: source, target
  class(mesh_t), pointer :: m
  real(8) :: pos(3)
  ! --------------------------------------------------------------------
  call trace%begin ('data_hub_t%dload_higher',itimer=itimer)
  source => self%source
  target => self%target
  allocate(cube(2,2,2,2))
  call source%time_interval (target%time, jt, pt)
  if (self%rt_target) then
    tl = self%t_l(:,1); tu = self%t_u(:,1)
    do k = tl(3), tu(3)
      do j = tl(2), tu(2)
        do i = tl(1), tu(1)
          l=[i,j,k]
          do ii = 1,3
          m => target%mesh(ii)
          pos(ii) = m%r(l(ii))
          end do
          !
          if ((self%all_cells) .or. &
                any((l-target%mesh%li) < 0 .or. (l-target%mesh%ui) > 0)) then
            call source%index_space2(target,pos, 1, l, p)
            u=l+1
            mem => source%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,:,1)
            do iv = 1, target%nv
              cube(:,:,:,1) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
              cube(:,:,:,2) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
              call self%interpolate(cube, pt, p,tcell)
              target%mem(i,j,k,iv,target%it,1) = tcell
            end do
          end if
        end do
      end do
    end do
  else
    if (present(only)) then
      iv1 = only
      iv2 = only
    else
      iv1 = 1
      iv2 = target%nv
    end if
    d = self%source%idx%d
    do iv = iv1, iv2
      tl = self%t_l(:,iv); tu = self%t_u(:,iv)
      do k = tl(3), tu(3)
        do j = tl(2), tu(2)
          do i = tl(1), tu(1)
            l=[i,j,k]
            do ii = 1,3
            m => target%mesh(ii)
            pos(ii) = m%r(l(ii))
            end do
            if ((self%all_cells) .or. &
                any((l-target%mesh%li) < 0 .or. (l-target%mesh%ui) > 0)) then
              call source%index_space2(target,pos, iv, l, p)
              u=l+1
              mem => source%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,:,1)
              if (source%pervolume(iv)) then
                cube(:,:,:,1) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1) / &
                                source%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(1),1)
                cube(:,:,:,2) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1) / &
                                source%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(2),1)
              else
                cube(:,:,:,1) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
                cube(:,:,:,2) = source%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
              end if
              if (source%unsigned(iv)) then
                call self%interpolate_unsigned(cube, pt, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell
              else if (source%pervolume(iv)) then
                call self%interpolate(cube, pt, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell * target%mem(i,j,k,d,target%it,1)
              else
                call self%interpolate(cube, pt, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell
              end if
            end if
          end do
        end do
      end do
    end do
  end if
  deallocate(cube)
  call trace%end(itimer)
END SUBROUTINE dnload_from_higher

!=======================================================================
!> Initialize the method when the source is of a lower level of refinement
!> than the target patch.
!=======================================================================
SUBROUTINE init_lower (self)
  class(data_hub_t) :: self
  ! ....................................................................
  class(mesh_t),  pointer:: m, m2
  class(patch_t), pointer:: source, target
  integer :: iv, i, j(3), nv
  real(8) :: eps, pl(3), pu(3), dp
  ! --------------------------------------------------------------------
  source => self%source
  target => self%target
  eps = 0.001_8
  if (self%rt_target) then
    allocate (self%t_l(3,1))
    allocate (self%t_u(3,1))
    allocate (self%s_l(3,1))
    allocate (self%s_u(3,1))
    nv = 1
  else
    allocate (self%t_l(3,target%nv))
    allocate (self%t_u(3,target%nv))
    allocate (self%s_l(3,target%nv))
    allocate (self%s_u(3,target%nv))
    nv = target%nv
  end if
  ! --------------------------------------------------------------------
  ! Corners of the source domain
  ! --------------------------------------------------------------------
  do iv = 1, nv
    do i = 1,3
      m => source%mesh(i)
      pl(i)  = m%lf - eps
      pu(i)  = m%uf + eps
    end do
  !---------------------------------------------------------------------
  ! Map to the range inside the target domain
  !---------------------------------------------------------------------
    self%t_l(:,iv) = target%index_only2 (source, pl, iv, roundup = .true.)
    self%t_u(:,iv) = target%index_only2 (source, pu, iv)
  ! --------------------------------------------------------------------
  ! Enforce proper limits
  ! --------------------------------------------------------------------
  if (.not.self%all_cells) then
    do i=1,3
      m => target%mesh(i)
      m2=> source%mesh(i)
      self%t_l(i,iv) = max(m%lb,min(m%uo,self%t_l(i,iv)))
      self%t_u(i,iv) = max(m%lo,min(m%ub,self%t_u(i,iv)))
      dp = round(m%p-m2%p,12)+m%h(iv)*(m%d-m2%d)
      if (self%t_u(i,iv)>m%lo) then
        if (m%r(m%li)+dp > m2%r(m2%ui)) &
        self%t_u(i,iv) = min(self%t_u(i,iv),m%lo)
      end if
      if (self%t_l(i,iv)<m%uo) then
        if (m%r(m%ui)+dp < m2%r(m2%li)) &
        self%t_l(i,iv) = max(self%t_l(i,iv),m%uo)
      end if
      if (self%t_u(i,iv)>m%ui) then
        if (m%r(m%uo)+dp > m2%uf) then
          self%t_u(i,iv) = min(self%t_u(i,iv),m%ui)
        else
          self%t_u(i,iv) = max(self%t_l(i,iv),m%ub)
        end if
      end if
      if (self%t_l(i,iv)<m%li) then
        if (m%r(m%lo)+dp < m2%lf) then
          self%t_l(i,iv) = max(self%t_l(i,iv),m%li)
        else
          self%t_l(i,iv) = min(self%t_l(i,iv),m%lb)
        end if
      end if
    end do
  end if
  ! --------------------------------------------------------------------
  ! Map the acquired indices back to the source
  ! --------------------------------------------------------------------
    do i = 1,3
      m => target%mesh(i)
      pl(i) = m%r(self%t_l(i,iv))
      pu(i) = m%r(self%t_u(i,iv))
    end do
    self%s_l(:,iv) = source%index_only2 (target,pl, iv)
    self%s_u(:,iv) = source%index_only2 (target,pu, iv, roundup=.true.)
    if (.not.self%all_cells) then
      do i = 1,3
        m => source%mesh(i)
        self%s_l(i,iv) = max(m%lo,min(m%ui,self%s_l(i,iv)))
        self%s_u(i,iv) = min(m%uo,max(m%li,self%s_u(i,iv)))
      end do
    end if
  end do
END SUBROUTINE init_lower

!=======================================================================
!> Download guard cell values from the source, which has a lower level
!> of refinement, than the target patch. A temporary scratch array is
!> created and are taken from the source patch. If the needed values lay
!> outside of the source%li/%ui, then a search in the target is done. If
!> the right location for the data is not in the target, a search in the
!> neighbors is done. All values are initially interpolated in time to the
!> target%time and then interpolated in space to get the required quantities.
!> Can be extended to accommodate the van Leer interpolation scheme.
!=======================================================================
SUBROUTINE dnload_from_lower (self, only)
  class(data_hub_t) :: self
  integer, optional :: only
  ! ....................................................................
  integer :: jt(2),jt2(2), iv, sl(3), su(3), tl(3), tu(3), i, j, k
  integer :: d, cl(3), ii, l(3), u(3), id, n(3), iv1, iv2
  class(mesh_t), pointer :: m
  class(patch_t), pointer:: source, target
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: scr
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: rt_scr, rt_cube
  real(kind=KindScalarVar), dimension(:,:,:),   pointer :: cube
  real(kind=KindScalarVar), pointer :: scr3
  real(kind=KindScalarVar), dimension(:),       pointer :: rt_scr3
  real(kind=KindScalarVar) :: tcell
  real    :: pt(2), pt2(2), p(3)
  real(8) :: pos(3), time, pos2(3)
  logical :: outside, found
  integer, save:: itimer=0
  ! --------------------------------------------------------------------
  call trace%begin ('data_hub_t%dnload_lower',itimer=itimer)
  source => self%source
  target => self%target
  d = source%idx%d
  mem => self%source%mem(:,:,:,:,:,1)
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = source%nv
  end if
  ! --------------------------------------------------------------------
  ! interpolate in time and fill in a scratch cube. When downloading from
  ! a larger patch, not all needed cells are usually available inside the
  ! source patch, and have to be retrieved either from the target or from
  ! other nbors.
  ! --------------------------------------------------------------------
  call source%time_interval (target%time, jt2, pt2)
  time = target%time
  id = -1
  if (self%rt_target) then
    sl = self%s_l(:,1); su = self%s_u(:,1)
    tl = self%t_l(:,1); tu = self%t_u(:,1)
    n = su - sl +2
    allocate (rt_scr(n(1), n(2), n(3),target%nv))
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          cl = [sl(1)+i-1,sl(2)+j-1,sl(3)+k-1]
          outside = .false.
          do ii = 1,3
            if ((cl(ii)<source%li(ii)).or.(cl(ii)>source%ui(ii))) outside = .true.
          end do
          if (outside) then
            found = .false.
            do ii = 1,3
              m => source%mesh(ii)
              pos(ii) = m%r(cl(ii))
              pos2(ii) = m%r(cl(ii))+m%p+m%h(1)*m%d
            end do
            rt_scr3 => rt_scr(i,j,k,:)
            call self%check_in_target(pos, 1, time, jt, pt, id, found, rt_scr=rt_scr3)    ! many cells are needed from the target for good interpolation
            if (.not.found) &
              call self%check_in_nbors (pos, 1, time, jt, pt, id, found, rt_scr=rt_scr3)  ! worst case.. at least necessary nbors will be same size or smaller.
            if (.not.found) then
                print *, iv, cl, pos
                print *, pos2
                call io%abort("data_hub::dnload_lower:: cell not found")    ! very unusual, but a necessary check and a stop.
            end if
          else
            rt_scr(i,j,k,:) = pt2(1) * mem(cl(1),cl(2),cl(3),:, jt2(1)) + &
                              pt2(2) * mem(cl(1),cl(2),cl(3),:, jt2(2))
          end if
        end do
      end do
    end do
   ! -------------------------------------------------------------------
   ! Now that the scratch cube is ready, we can spatially interpolate
   ! back to target
   ! -------------------------------------------------------------------
    do k=tl(3), tu(3)
      do j=tl(2), tu(2)
        do i=tl(1), tu(1)
          cl = [i,j,k]
          if ((self%all_cells) .or. &
              any((cl-target%mesh%li) < 0 .or. (cl-target%mesh%ui) > 0)) then
            do ii = 1,3
              m => target%mesh(ii)
              pos(ii) = m%r(cl(ii))
            end do
            call source%index_space2(target,pos, 1, cl, p)
            l  = cl - sl + 1
            u  = l +1
            do iv = 1, target%nv
              cube => rt_scr(l(1):u(1),l(2):u(2),l(3):u(3), iv)
              call self%interpolate_one(cube, p,tcell)
              target%mem(i,j,k,iv,target%it,1) = tcell
            end do
          end if
        end do
      end do
    end do
    deallocate(rt_scr)
  else
    do iv = iv1, iv2
      sl = self%s_l(:,iv); su = self%s_u(:,iv)
      tl = self%t_l(:,iv); tu = self%t_u(:,iv)
      n = su - sl +2
      allocate (scr(n(1), n(2), n(3)))
      do k = 1, n(3)
        do j = 1, n(2)
          do i = 1, n(1)
            cl = [sl(1)+i-1,sl(2)+j-1,sl(3)+k-1]
            outside = .false.
            do ii = 1,3
              if ((cl(ii)<source%li(ii)).or.(cl(ii)>source%ui(ii))) outside = .true.
            end do
            if (outside) then
              found = .false.
              do ii = 1,3
                m => source%mesh(ii)
                pos(ii) = m%r(cl(ii))
                pos2(ii) = m%r(cl(ii))+m%p+m%h(iv)*m%d
              end do
              scr3 => scr(i,j,k)
              call self%check_in_target(pos, iv, time, jt, pt, id, found, scr=scr3)    ! many cells are needed from the target for good interpolation
              if (.not.found) &
                call self%check_in_nbors (pos, iv, time, jt, pt, id, found, scr=scr3)  ! worst case.. at least necessary nbors will be same size or smaller.
              if (.not.found) then
                  print *, iv, cl, pos
                  print *, pos2
                  call io%abort("data_hub::dnload_lower:: cell not found")    ! very unusual, but a necessary check and a stop.
              end if
            else
              if (source%unsigned(iv)) then
                scr(i,j,k) = exp(pt2(1) * log(mem(cl(1),cl(2),cl(3),iv, jt2(1))) + &
                                 pt2(2) * log(mem(cl(1),cl(2),cl(3),iv, jt2(2))))
              else if (source%pervolume(iv)) then
                scr(i,j,k) = pt2(1) * (mem(cl(1),cl(2),cl(3),iv, jt2(1))  / &
                                       mem(cl(1),cl(2),cl(3), d, jt2(1))) + &
                             pt2(2) * (mem(cl(1),cl(2),cl(3),iv, jt2(2))  / &
                                       mem(cl(1),cl(2),cl(3), d, jt2(2)))
              else
                scr(i,j,k) = pt2(1) * mem(cl(1),cl(2),cl(3),iv, jt2(1)) + &
                             pt2(2) * mem(cl(1),cl(2),cl(3),iv, jt2(2))
              end if
            end if
          end do
        end do
      end do
     ! -------------------------------------------------------------------
     ! Now that the scratch cube is ready, we can spatially interpolate
     ! back to target
     ! -------------------------------------------------------------------
      do k=tl(3), tu(3)
        do j=tl(2), tu(2)
          do i=tl(1), tu(1)
            cl = [i,j,k]
            if ((self%all_cells) .or. &
                any((cl-target%mesh%li) < 0 .or. (cl-target%mesh%ui) > 0)) then
              do ii = 1,3
                m => target%mesh(ii)
                pos(ii) = m%r(cl(ii))
              end do
              call source%index_space2(target,pos, iv, cl, p)
              l  = cl - sl + 1
              u  = l +1
              cube => scr(l(1):u(1),l(2):u(2),l(3):u(3))
              if (source%unsigned(iv)) then
                call self%interpolate_unsigned_one(cube, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell
              else if (self%source%pervolume(iv)) then
                call self%interpolate_one(cube, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell * target%mem(i,j,k,d,target%it,1)
              else
                call self%interpolate_one(cube, p,tcell)
                target%mem(i,j,k,iv,target%it,1) = tcell
              end if
            end if
          end do
        end do
      end do
      deallocate(scr)
    end do
 end if
 call trace%end(itimer)
END SUBROUTINE dnload_from_lower

!=======================================================================
!> Some support data is needed, when download_lower is done. Usually it
!> resides in the target patch. This data is spatially interpolated to
!> the source resolution and added to the scratch array.
!=======================================================================
SUBROUTINE check_in_target (self, pos, iv, time, jt, pt, id, found, scr, rt_scr)
  class(data_hub_t)                 :: self
  real(8)                           :: pos(3), time
  real(kind=KindScalarVar), pointer, optional :: scr
  real(kind=KindScalarVar), dimension(:),pointer, optional :: rt_scr
  integer                           :: iv, jt(2), id
  real                              :: pt(2)
  logical                           :: found
  ! ....................................................................
  integer :: sl(3), su(3), tl(3), tu(3), i, j, k, aa, bb, cc
  integer :: d, cl(3), ii, l(3), u(3), v
  real    :: p(3)
  real(8) :: eps = 0.001_8
  class(mesh_t), pointer :: m
  class(patch_t), pointer:: source, target
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: cube
  logical :: inside
  ! --------------------------------------------------------------------
  source => self%source
  target => self%target
  d = source%idx%d
  inside = .true.
  call target%index_space2(source,pos,iv, cl, p)
  inside = all(cl>=target%mesh%li .and. cl<=target%mesh%ui)
  if (inside) then                                                      ! do the interpolation from smaller
  ! --------------------------------------------------------------------
  ! do spatial interpolation.
  ! --------------------------------------------------------------------
    if (id /=target%id) then
      call target%time_interval (time, jt, pt)
      id = target%id
    end if
    allocate(cube(2,2,2,2))
    l = cl
    u = cl+1
    if (self%rt_target) then
      do v=1,target%nv
        cube(:,:,:,1) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),v,jt(1),1)
        cube(:,:,:,2) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),v,jt(2),1)
        call self%interpolate(cube, pt, p,rt_scr(v))
      end do
    else
      if (source%unsigned(iv)) then
        cube(:,:,:,1) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
        cube(:,:,:,2) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
        call self%interpolate_unsigned(cube, pt, p,scr)
      else if (source%pervolume(iv)) then
        cube(:,:,:,1) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1) / &
                        target%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(1),1)
        cube(:,:,:,2) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1) / &
                        target%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(2),1)
        call self%interpolate(cube, pt, p,scr)
      else
        cube(:,:,:,1) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
        cube(:,:,:,2) = target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
        call self%interpolate(cube, pt, p,scr)
      end if
    end if
    deallocate(cube)
    found = .true.
  end if
END SUBROUTINE check_in_target

!=======================================================================
!> In some cases, e.g. corners, ends of edges, some support data is located
!> not in the target patch, but rather in adjacent nbors. They have the
!> same level of refinement as the target or the source. This data is
!> spatially interpolated to the source resolution, if needed as well as
!> interpolated in time to the needed target time. Then it is added
!> to the scratch array.
!=======================================================================
SUBROUTINE check_in_nbors (self, pos, iv, time, jt, pt, id, found, scr, rt_scr)
  class(data_hub_t)                 :: self
  real(8)                           :: pos(3), time
  real(kind=KindScalarVar), pointer, optional :: scr
  real(kind=KindScalarVar), dimension(:),pointer, optional :: rt_scr
  integer                           :: jt(2), id
  real                              :: pt(2)
  logical                           :: found
  ! ....................................................................
  integer :: iv, sl(3), su(3), tl(3), tu(3), i, j, k
  integer :: d, cl(3), ii, l(3), u(3), v
  real    :: p(3)
  real(8) :: eps = 0.001_8
  class(mesh_t), pointer :: m
  class(link_t), pointer :: nbor
  class(patch_t), pointer:: source, target, nbpatch
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: cube
  logical :: inside
  ! --------------------------------------------------------------------
  source => self%source
  target => self%target
  d = source%idx%d
  nbor => source%link%nbor
  do while((associated(nbor)).and..not.found)
    if ((nbor%task%id /= source%id).and.(nbor%task%id /= target%id)) then ! source is already checked, target should not be there
      nbpatch => task2patch(nbor%task)
      if (abs(nbpatch%level-source%level)<2) then
        if (nbpatch%level == source%level) then                         ! same level
          cl = nbpatch%index_only2 (source,pos, iv, closest=.true.)
          inside = all(cl>=nbpatch%mesh%li .and. cl<=nbpatch%mesh%ui)
          if (inside) then                                              ! same size, just take the value by temporal interpolation
            if (id /=nbpatch%id) then
              call nbpatch%time_interval (time, jt, pt)
              id = nbpatch%id
            end if
            if(self%rt_target) then
              rt_scr =  &
              pt(1) * nbpatch%mem(cl(1),cl(2),cl(3),:,jt(1),1) + &
              pt(2) * nbpatch%mem(cl(1),cl(2),cl(3),:,jt(2),1)
              found = .true.
            else
              if (source%unsigned(iv)) then
                scr = exp ( &
                pt(1) * log(nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(1),1)) + &
                pt(2) * log(nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(2),1)))
              else if (source%pervolume(iv)) then
                scr =  &
                pt(1) * (nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(1),1) / &
                           nbpatch%mem(cl(1),cl(2),cl(3), d,jt(1),1))+ &
                pt(2) * (nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(2),1) / &
                           nbpatch%mem(cl(1),cl(2),cl(3), d,jt(2),1))
              else
                scr =  &
                pt(1) * nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(1),1) + &
                pt(2) * nbpatch%mem(cl(1),cl(2),cl(3),iv,jt(2),1)
              end if
              found = .true.
            end if
          end if
        else
          call nbpatch%index_space2(source,pos,iv, cl, p)
          inside = all(cl>=nbpatch%mesh%li .and. cl<=nbpatch%mesh%ui)
          if (inside) then                                              ! different size, spatial and temporal interpolation
            if (id /=nbpatch%id) then
              call nbpatch%time_interval (time, jt, pt)
              id = nbpatch%id
            end if
            allocate(cube(2,2,2,2))
            l = cl
            u = cl+1
            if (self%rt_target) then
              do v = 1, target%nv
                cube(:,:,:,1) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),v,jt(1),1)
                cube(:,:,:,2) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),v,jt(2),1)
                call self%interpolate(cube, pt, p,rt_scr(v))
              end do
            else
              if (source%unsigned(iv)) then
                cube(:,:,:,1) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
                cube(:,:,:,2) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
                call self%interpolate_unsigned(cube, pt, p,scr)
              else if (source%pervolume(iv)) then
                cube(:,:,:,1) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1) / &
                                nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(1),1)
                cube(:,:,:,2) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1) / &
                                nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3), d,jt(2),1)
                call self%interpolate(cube, pt, p,scr)
              else
                cube(:,:,:,1) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(1),1)
                cube(:,:,:,2) = nbpatch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,jt(2),1)
                call self%interpolate(cube, pt, p,scr)
              end if
            end if
            deallocate(cube)
            found = .true.
          end if
        end if
      end if
    end if
    nbor => nbor%next
  end do
END SUBROUTINE check_in_nbors

!=======================================================================
!> Update the guard zone values
!=======================================================================
SUBROUTINE update (self, link, all_cells, only)
  class(data_hub_t):: self
  class(link_t), pointer :: link
  logical, optional :: all_cells
  integer, optional :: only
  ! ....................................................................
  class  (data_hub_t), pointer :: hub
  integer, save:: itimer=0
  ! --------------------------------------------------------------------
  !call trace%begin ('data_hub_t%update',itimer=itimer)
  if (self%first_time) call self%init(link, all_cells)                  ! map data shell if it is the first time
  !
  hub => self%head
  do while(associated(hub))
    call hub%dnload(only=only)
    hub => hub%next
  end do
  !call trace%end(itimer)
END SUBROUTINE update

SUBROUTINE interpolate (self, src, pt, p, res)
  class(data_hub_t)        :: self
  real(kind=KindScalarVar), dimension(:,:,:,:):: src
  real                     :: pt(2), p(3)
  real(kind=KindScalarVar) :: res
  !---------------------------------------------------------------------
    if (any(p < 0.0).or.any(p>1.0)) then
    print *, p
    call io%abort("p negative")
    end if
    res = &
     pt(1)*((1.-p(3))*((1.-p(2))*((1.-p(1))*src(1,1,1,1)    + &
                                  (   p(1))*src(2,1,1,1))   + &
                       (   p(2))*((1.-p(1))*src(1,2,1,1)    + &
                                  (   p(1))*src(2,2,1,1)))  + &
            (   p(3))*((1.-p(2))*((1.-p(1))*src(1,1,2,1)    + &
                                  (   p(1))*src(2,1,2,1))   + &
                       (   p(2))*((1.-p(1))*src(1,2,2,1)    + &
                                  (   p(1))*src(2,2,2,1)))) + &
     pt(2)*((1.-p(3))*((1.-p(2))*((1.-p(1))*src(1,1,1,2)    + &
                                  (   p(1))*src(2,1,1,2))   + &
                       (   p(2))*((1.-p(1))*src(1,2,1,2)    + &
                                  (   p(1))*src(2,2,1,2)))  + &
            (   p(3))*((1.-p(2))*((1.-p(1))*src(1,1,2,2)    + &
                                  (   p(1))*src(2,1,2,2))   + &
                       (   p(2))*((1.-p(1))*src(1,2,2,2)    + &
                                  (   p(1))*src(2,2,2,2))))
END SUBROUTINE interpolate

SUBROUTINE interpolate_unsigned (self, src, pt, p, res)
  class(data_hub_t)        :: self
  real(kind=KindScalarVar), dimension(:,:,:,:):: src
  real                     :: pt(2), p(3)
  real(kind=KindScalarVar) :: res
  !---------------------------------------------------------------------
  if (any(p < 0.0).or.any(p>1.0)) then
    print *, p
    call io%abort("p negative")
    end if
    res = exp ( &
      pt(1)*((1.-p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,1,1))    + &
                                   (   p(1))*log(src(2,1,1,1)))   + &
                        (   p(2))*((1.-p(1))*log(src(1,2,1,1))    + &
                                   (   p(1))*log(src(2,2,1,1))))  + &
             (   p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,2,1))    + &
                                   (   p(1))*log(src(2,1,2,1)))   + &
                        (   p(2))*((1.-p(1))*log(src(1,2,2,1))    + &
                                   (   p(1))*log(src(2,2,2,1))))) + &
      pt(2)*((1.-p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,1,2))    + &
                                   (   p(1))*log(src(2,1,1,2)))   + &
                        (   p(2))*((1.-p(1))*log(src(1,2,1,2))    + &
                                   (   p(1))*log(src(2,2,1,2))))  + &
             (   p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,2,2))    + &
                                   (   p(1))*log(src(2,1,2,2)))   + &
                        (   p(2))*((1.-p(1))*log(src(1,2,2,2))    + &
                                   (   p(1))*log(src(2,2,2,2))))))
END SUBROUTINE interpolate_unsigned

SUBROUTINE interpolate_one (self, src, p, res)
  class(data_hub_t)        :: self
  real(kind=KindScalarVar), dimension(:,:,:):: src
  real                     :: p(3)
  real(kind=KindScalarVar) :: res
  !---------------------------------------------------------------------
    if (any(p < 0.0).or.any(p>1.0)) then
    print *, p
    call io%abort("p negative")
    end if
    res = &
    ((1.-p(3))*((1.-p(2))*((1.-p(1))*src(1,1,1)    + &
                           (   p(1))*src(2,1,1))   + &
                (   p(2))*((1.-p(1))*src(1,2,1)    + &
                           (   p(1))*src(2,2,1)))  + &
     (   p(3))*((1.-p(2))*((1.-p(1))*src(1,1,2)    + &
                           (   p(1))*src(2,1,2))   + &
                (   p(2))*((1.-p(1))*src(1,2,2)    + &
                           (   p(1))*src(2,2,2))))
END SUBROUTINE interpolate_one

SUBROUTINE interpolate_unsigned_one (self, src, p, res)
  class(data_hub_t)        :: self
  real(kind=KindScalarVar), dimension(:,:,:):: src
  real                     :: p(3)
  real(kind=KindScalarVar) :: res
  !---------------------------------------------------------------------
  if (any(p < 0.0).or.any(p>1.0)) then
    print *, p
    call io%abort("p negative")
    end if
    res = exp &
    ((1.-p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,1))    + &
                           (   p(1))*log(src(2,1,1)))   + &
                (   p(2))*((1.-p(1))*log(src(1,2,1))    + &
                           (   p(1))*log(src(2,2,1))))  + &
     (   p(3))*((1.-p(2))*((1.-p(1))*log(src(1,1,2))    + &
                           (   p(1))*log(src(2,1,2)))   + &
                (   p(2))*((1.-p(1))*log(src(1,2,2))    + &
                           (   p(1))*log(src(2,2,2)))))
END SUBROUTINE interpolate_unsigned_one

END MODULE data_hub_mod
