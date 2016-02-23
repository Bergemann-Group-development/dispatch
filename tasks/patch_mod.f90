!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE trace_mod
  USE mpi_mod
  USE task_mod
  USE mesh_mod
  USE io_mod
  implicit none
  private
  integer, parameter:: ndim = 3
  type, extends(tasks_t), public:: patch_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    type(mesh_t), dimension(3):: mesh
    real:: courant
    real:: gamma=1.0
    integer:: ndim = ndim           ! Number of dimensions
    !integer:: nvar                  ! Number of variables
    integer:: nv                    ! FIXME
    integer:: n(3)
    integer:: msplit=2
    !---------------------------------------------------------------------------
    ! HACK below
    !---------------------------------------------------------------------------
    integer:: li(3), ui(3), gn(3)
    integer:: i,j,k
    real:: box(3), ds(3), gsize(3)
    real(8):: offset(3)
    real:: dt_fixed, u_fixed, u_max
    logical:: mem_allocated=.false.
  contains
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: distance        ! return the distance btw two patches
    procedure:: overlaps        ! return true if they overlap (incl. ghost zones)
    procedure:: upload_range    ! returns the index ranges btw two points
    procedure:: dnload_range    ! returns the index ranges btw two points
    procedure:: index           ! compute the index for a single point
    procedure:: nearest         ! the nearest index that maps to patch point
    procedure:: lowest          ! the lowest index that maps to patch point
    procedure:: highest         ! the highest index that maps to patch point
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: ahead_of
    procedure:: download_to
    procedure:: myposition
    procedure:: periodic_grid
    procedure:: extrapolate_guards
    procedure:: fmaxval         ! interior max value
    procedure:: fminval         ! interior min value
    procedure:: courant_condition
  end type
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type box_t
    real:: size(3)
    real:: num(3)
    integer:: dim(3)
  end type
  integer:: id=0
CONTAINS

!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  integer:: n(3)=16, ng=2, nt=7, i
  logical:: first_time = .true.
  real:: position(3)
  real, save:: courant=0.2, u_fixed=-1.0, dt_fixed=-1.0
  real, save:: end_time=1.0, out_time=0.1
  namelist /patch_params/ n, ng, nt, courant, u_fixed, dt_fixed, end_time, out_time
  !.............................................................................
  call trace_begin ('patch_t%init', 3)
  !$omp atomic
  id = id + 1                                   ! increment global
  if (io%verbose>2) &
    print*,'patch_t%init: id=', id
  self%id = id                                  ! new patch id
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
  end if
  self%courant = courant
  self%dt_fixed = dt_fixed
  self%u_fixed = u_fixed
  self%status=0
  self%end_time=end_time
  self%out_time=out_time
  self%n  = n                                   ! patch dimensions
  self%ds = real(self%size,kind=8) / self%n     ! ds in patch
  self%mesh%n  = n
  self%mesh%ng = ng
  self%mesh%s  = self%size
  where (self%n==1)                             ! for collapsed dimensions
    self%mesh%ng = 0                            ! no guard zones
    self%size = 0.0                             ! zero size
    self%ds = 0.0                               ! first zero
    self%box = 1                                ! make sure it is non-zero
    self%ds = maxval(self%ds)                   ! then set to reasonable val
  end where
  !-----------------------------------------------------------------------------
  ! Initialize the meshes
  !-----------------------------------------------------------------------------
  do i=1,3
    call self%mesh(i)%init
  end do
  !-----------------------------------------------------------------------------
  ! Hack for backwards compatibility
  !-----------------------------------------------------------------------------
  self%gn = self%mesh%gn
  self%li = self%mesh%li
  self%ui = self%mesh%ui
  self%offset = (self%mesh%lb+self%mesh%ub)/2.0 ! mid point index offset
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  self%nt = nt                                  ! number of time slices
  self%level = log(maxval(self%box)/minval(self%ds))/log(real(self%msplit)) + 0.1
  if (io%verbose>3) &
    print*,'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not allocated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
!    if (associated(self%mem)) deallocate (self%mem)
    allocate (self%mem(self%mesh(1)%gn, &
                       self%mesh(2)%gn, &
                       self%mesh(3)%gn, &
                       self%nv,self%nt,2))
    self%mem_allocated = .true.                 ! mark
  end if
  self%t = 0.0                                  ! initial time
  self%iit = 1                                  ! initial memory slots
  self%iit(nt) = 2                              ! first output slot
  self%it = self%iit(nt-1)                      ! short-hand
  self%new = self%iit(nt)                       ! short-hand
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE update (self)
  class(patch_t):: self
  !.............................................................................
  call trace_begin ('patch_t%update')
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    print'(a,7i10)',   'experiment_t%update: iit =', self%iit
    print'(a,7f10.6)', 'experiment_t%update:   t =', self%t(self%iit)
  end if
  ! FIXME
  !call self%courant_condition
  !call self%timestep%update (self%id, self%iit, self%dt, self%mem)
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation
  !-----------------------------------------------------------------------------
  !if (self%id==1) then
  !  call self%periodic_grid                     ! periodic root grid
  !else
  !  call self%extrapolate_guards                ! time-extrapolate guard cells
  !end if
  !-----------------------------------------------------------------------------
  ! One line summary
  !-----------------------------------------------------------------------------
  self%time = self%time + self%dtime                    ! update self time
  if (io%verbose>0) then
    print'(a,i4,2i3,2f11.6,2x,a,1p,6e10.2)', 'update: id,j,k,l,n,dt,t', &
      self%id, self%it, self%new, self%dtime, self%time, 'minmax:', &
      minval(self%mem(:,:,:,1,self%it ,2)), &
      maxval(self%mem(:,:,:,1,self%it ,2)), &
      minval(self%mem(:,:,:,1,self%new,1)), &
      maxval(self%mem(:,:,:,1,self%new,1))
  end if
  print '(i8,f10.6,f9.6,f7.2,2x,1p,2e10.2,2x,2e14.6,i6,2x,7i3,2x,2i3)', self%id, &
    self%time, self%dtime, self%courant, &
    minval(self%mem(:,:,:,1,self%it ,2)), &
    maxval(self%mem(:,:,:,1,self%it ,2)), &
    minval(self%mem(:,:,:,1,self%new,1)), &
    maxval(self%mem(:,:,:,1,self%new,1)), &
    self%status, self%iit, self%it, self%new
  !-----------------------------------------------------------------------------
  ! Time slot rotation
  !-----------------------------------------------------------------------------
  self%dt(self%it) = self%dtime                         ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = self%time                          ! initial time
  self%iit(self%nt) = self%new                          ! right-most slot
  self%it = self%new                                    ! update time slot
  self%new = mod(self%new,self%nt)+1                    ! increment / rotate
  self%iit(1:self%nt-1) = self%iit(2:self%nt)           ! rotate left
  self%iit(self%nt) = self%new                          ! new right-most slot
  self%istep = self%istep + 1
  if (io%verbose>2) &
   print*,'patch id:',self%id,'  it:',self%it,'  new:',self%new
  if (io%verbose>3) &
   print*,'patch iit:', self%iit
  if (io%verbose>3) &
   print*,'patch t(iit):', self%t(self%iit)
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Download values from the source (self) to the target (patch), trying to get
!> all cells in the target if "all_target" is true, and only guard cells if not.
!> We use myposition() to map an index position in the target patch to an
!> absolute position, and then index() to get the integer and floating point
!> position in the source patch.  With that, we make a tri-linear interpolation.
!> This is useful for both guard cell fill-in, and for "prolongation" = creating
!> patches with higher resolution.
!> =============================================================================
SUBROUTINE download_to (self, task, special, all_cells)
  class(patch_t):: self
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  logical, optional:: all_cells
  integer:: special
  logical:: all_target
  integer:: ii(3), ix, iy, iz, it, iv, j0(3), j1(3), ntot, nuse
  integer, dimension(3):: l1, u1, l2, u2
  real:: f(3), pos(3), weight
  !-----------------------------------------------------------------------------
  ! Upgrade to a patch link
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
  patch => task
  end select
!  if (self%id==patch%id .and. all(patch%size /= box%size)) return
  if (io%verbose>2) print'(a,2(i3,3f8.3),2(3x,3i4))','download_to:self,patch=', &
    self%id,self%position,patch%id,patch%position, self%mesh%li, self%mesh%ui
associate (source=>self%mem, target=>patch%mem)
  it = patch%it
  !-----------------------------------------------------------------------------
  ! Download rules:  1) Always download guard cells.  2) Download all cells if
  ! a) all_cells is set, and the source level is at most one level higher. In
  !    this case (first fill-in), the source level may be lower than target level
  ! b) the level is one level higher
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    all_target = .true.                                   ! default: use all points
    if (self%level <= patch%level) all_target = .false.   ! only upload from deeper levels
    if (all_cells) all_target = .true.                    ! unless told to do so
    if ((self%level-patch%level) > 1) return              ! upload at most 1 level
    if (iv==special) then                               ! if the variable is phi
      if (self%level >= patch%level) cycle                ! only download phi from lower levels
      !all_target = .true.                                 ! interior + guard cells
      !if (io%verbose>0) &
      !  print'(10(1x,a,i6,3x))', 'patch', patch%id, 'at level', patch%level, &
      !  'downloading variable', special, 'from patch', self%id, 'level', self%level
    end if

    !---------------------------------------------------------------------------
    ! When the source resolution is higher than the target resolution, all cells
    ! in the source need to be distributed over the nearest cells in the target.
    ! First, we need to find the index ranges in the source that map entirely
    ! into the target. For simplicity, we first ignore the periodic mapping of
    ! the box, and assume all patches are inside the box.  Then the range of an
    ! index is the cross section of the whole source range with the range of
    ! target mapped into the source. Since this is to be used for uploading,
    ! there is no need to go the very limit, so to keep things simple we require
    ! a target point to exist both above and below a source point, for it to be
    ! included in the range.
    !---------------------------------------------------------------------------
    if (self%level > patch%level) then                    ! source has finer cells
      !-------------------------------------------------------------------------
      ! upload_range (p1, p2, li, ui, l1, u1, l2, u2) gives the index ranges in
      ! patches p1 and p2, respectively, corresponding to p2%li and p2%mesh%ui, of a
      ! region inside p2%li-p2%mesh%ui, which is also inside p1%mesh%lb-p1%mesh%ub. In this 
      ! context p1 (patch) is the target patch, and p2 (self) is the source patch.
      !-------------------------------------------------------------------------
      call patch%upload_range (self, l1, u1, l2, u2)
      !print '(2(a,i4,2(2x,3i4)))', &
        !'uploading to patch', patch%id, l1, u1, ' from patch', self%id, l2, u2
      target(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),iv,it,1) = 0.0
      weight = product(self%ds)/product(patch%ds)
      !-------------------------------------------------------------------------
      ! If we are using all target points, this is a straight loop, in which
      ! the calculation of the target index is the most costly part, so it may
      ! be a good idea to use an inline expression here.
      !-------------------------------------------------------------------------
      if (all_target) then
        do iz=l2(3),u2(3)
        do iy=l2(2),u2(2)
        do ix=l2(1),u2(1)
          ii = [ix,iy,iz]
          j0 = patch%nearest (self, ii)
         !j0 = nint(patch%offset+(self%position+(index-self%offset)*self%ds-patch%position)/patch%ds)
          target(j0(1),j0(2),j0(3),iv,it,1) = &
          target(j0(1),j0(2),j0(3),iv,it,1) + weight*source(ix,iy,iz,iv,it,1)
        end do
        end do
        end do
      !-------------------------------------------------------------------------
      ! If we are only using guard cells we just need to exclude the interior of
      ! the target.  This is much simpler and cheaper than treating all 26 guard
      ! zone regions separately. For 15-16 cell interiors with 2+2 guard cells
      ! the number of excluded points is roughly half of all points considered.
      !-------------------------------------------------------------------------
      else
        do iz=l2(3),u2(3)
        do iy=l2(2),u2(2)
        do ix=l2(1),u2(1)
          ii = [ix,iy,iz]
          j0 = patch%nearest (self, ii)
          if (any(j0 >= patch%mesh%li .or. j0 <= patch%mesh%ui)) cycle
          target(j0(1),j0(2),j0(3),iv,it,1) = &
          target(j0(1),j0(2),j0(3),iv,it,1) + weight*source(ix,iy,iz,iv,it,1)
        end do
        end do
        end do
      end if
    else
      call patch%dnload_range (self, l1, u1)
      !if (iv==1) then
      !  print 1, self%id , ' source position, size:', self%position , self%size , self%ub , self%it , minval(self%mem(:,:,:,1,self%it,1))
      !  print 1, patch%id, ' target position, size:', patch%position, patch%size, patch%ub, patch%it, minval(self%mem(:,:,:,1,patch%it,1))
      !  1 format(i4,a,2(3f8.3,2x),3i4,i5,f8.3)
      !  print*,'l1, u1', l1, u1
      !end if
      ntot = 0
      nuse = 0
      do iz=l1(3),u1(3)
      do iy=l1(2),u1(2)
      do ix=l1(1),u1(1)
        ntot = ntot+1
        ii = [ix,iy,iz]
        !---------------------------------------------------------------------------
        ! Check if the target point is in the guard zone region (= not inside)
        !---------------------------------------------------------------------------
        if (all_target .or. .not.all((ii >= patch%mesh%li) .and. (ii <= patch%mesh%ui))) then
          j0 = self%index (patch, ii, f)
          !---------------------------------------------------------------------------
          ! Check that the source is in the valid region
          !---------------------------------------------------------------------------
          if (all((j0+f >= self%mesh%lb-0.01 .and. j0+f <= self%mesh%ub+0.01) .or. self%n==1)) then
            j1 = j0+1
            where (j1 > self%mesh%ub)
              j0 = j0-1
              j1 = j1-1
              f  = f+1.0
            end where
            where (j0 < self%mesh%lb)
              j0 = j0+1
              j1 = j1+1
              f  = f-1.0
            end where
            where (self%n==1)
              j0 = 1
              j1 = 1
            end where
            target(ix,iy,iz,iv,it,1) = &
              (1.-f(3))*(1.-f(2))*((1.-f(1))*source(j0(1),j0(2),j0(3),iv,it,1)    &
                                     + f(1) *source(j1(1),j0(2),j0(3),iv,it,1)) + &
              (1.-f(3))*(   f(2))*((1.-f(1))*source(j0(1),j1(2),j0(3),iv,it,1)    &
                                     + f(1) *source(j1(1),j1(2),j0(3),iv,it,1)) + &
              (   f(3))*(1.-f(2))*((1.-f(1))*source(j0(1),j0(2),j1(3),iv,it,1)    &
                                     + f(1) *source(j1(1),j0(2),j1(3),iv,it,1)) + &
              (   f(3))*(   f(2))*((1.-f(1))*source(j0(1),j1(2),j1(3),iv,it,1)    &
                                     + f(1) *source(j1(1),j1(2),j1(3),iv,it,1))
            !if (patch%id==3 .and. iv==1) &
            !  print'(a,2(i4,2x,3i3),3f8.2,2x,8f8.2,3x,f8.2)','download:',self%id,ix,iy,iz,it,j0,f, &
            !                                 source(j0(1),j0(2),j0(3),iv,it,1), &
            !                                 source(j1(1),j0(2),j0(3),iv,it,1), &
            !                                 source(j0(1),j1(2),j0(3),iv,it,1), &
            !                                 source(j1(1),j1(2),j0(3),iv,it,1), &
            !                                 source(j0(1),j0(2),j1(3),iv,it,1), &
            !                                 source(j1(1),j0(2),j1(3),iv,it,1), &
            !                                 source(j0(1),j1(2),j1(3),iv,it,1), &
            !                                 source(j1(1),j1(2),j1(3),iv,it,1), &
            !                                 target(ix,iy,iz,iv,it,1)
            nuse = nuse+1
          else
            if (iv==1) print'("outside",i4,2x,3i4,2x,i4,3f8.3,3(2x,3i4))',patch%id,ix,iy,iz,self%id,j0+f,j0,self%mesh%lb,self%mesh%ub
          end if
        end if
      end do
      end do
      end do
    end if
    !if (patch%id==1) print'(7(1x,a,i7,1x))', 'patch 1, var', iv, ' at level', patch%level, &
    !  ' picked up', nuse, 'of', ntot, ' points from patch', self%id, ' at level', self%level
    if (iv==1.and.io%verbose>1) print'(7(1x,a,i7,1x))', 'patch', patch%id, ' level', patch%level, &
      ' picked up', nuse, 'of', ntot, ' points from patch', self%id, ' at level', self%level
    if (io%verbose>2) print'(7(1x,a,i7,1x),2l4)', 'patch', patch%id, ' level', patch%level, 'var', iv+10*self%nv, &
      ' picked up', nuse, 'of', ntot, ' points from patch', self%id, ' at level', self%level, &
      all_target, all_cells
  end do
  if (patch%t(patch%it)/=self%t(patch%it)) &
    print'(2(a,i4,f9.4,4x),2i3,2x,5i3,7f8.3)','ERROR: patch:', patch%id, patch%t(patch%it), &
      'self:', self%id, self%t(self%it), patch%it, patch%new, patch%iit, patch%t
end associate
END SUBROUTINE download_to

!> =============================================================================
!> This function gives the physical position of an index triple (patch point)
!> =============================================================================
FUNCTION myposition (self, ii) RESULT(out)
  class(patch_t):: self
  real:: out(3)
  integer:: ii(3)
  !.............................................................................
  out = self%position + (ii-self%offset)*self%ds
END FUNCTION
!> =============================================================================
!> Periodic mapping of a grid -- typically the root grid
!> =============================================================================
SUBROUTINE periodic_grid (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), j(3), lb(3), ub(3), li(3), ui(3)
  !.............................................................................
  lb = self%mesh%lb
  ub = self%mesh%ub
  li = self%mesh%li
  ui = self%mesh%ui
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    i = [ix,iy,iz]
    if (all(i >= li .and. i <= ui)) cycle
    j = modulo(i-lb,self%n)+lb
    self%mem(ix,iy,iz,:,self%new,1) = self%mem(j(1),j(2),j(3),:,self%new,1)
  end do
  end do
  end do
END SUBROUTINE periodic_grid
!> =============================================================================
!> Time extrapolate guard cell values, to allow using also guard zone values to
!> get boundary conditions for higher level patches
!> =============================================================================
SUBROUTINE extrapolate_guards (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), it0, it1, it2
  !.............................................................................
  call trace_begin ('patches_mod:extrapolate_guards', 3)
  !print*,'it,new,iit(nt):',self%id,self%it,self%new,self%iit(self%nt)
  !print*,'min,max:',minval(self%mem(:,:,:,1,self%it,1)),maxval(self%mem(:,:,:,1,self%it,1))
  it0 = self%new                                        ! current time slot
  it1 = modulo(it0-2,self%nt)+1                         ! previous time slot
  it2 = modulo(it0-3,self%nt)+1                         ! the one before that
  do iz=self%mesh(3)%lb,self%mesh(3)%ub                 ! loop over all cells
  do iy=self%mesh(2)%lb,self%mesh(2)%ub
  do ix=self%mesh(1)%lb,self%mesh(1)%ub
    i = [ix,iy,iz]
    if (all(i >= self%mesh%li .and. i <= self%mesh%ui)) cycle     ! skip interior
    self%mem(ix,iy,iz,:,it0,1)     = &
    self%mem(ix,iy,iz,:,it1,1)*2.0 - &
    self%mem(ix,iy,iz,:,it2,1)
  end do
  end do
  end do
  !print*,'min,max:',minval(self%mem(:,:,:,1,self%new,1)),maxval(self%mem(:,:,:,1,self%new,1))
  call trace_end
END SUBROUTINE extrapolate_guards
!> =============================================================================

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(patch_t):: self
  character(len=32) filename, logname, experiment_name
  optional:: experiment_name
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3)
  real:: ds(3)
  logical:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! If the out_next time has passed, set a new time = the integer number of
  ! out_time's that have passed, plus one more
  !-----------------------------------------------------------------------------
  if (self%time+0.5*self%dtime < self%out_next) return
  call trace_begin('patch:t%output')
  if (self%time==0.0) then
    self%out_next = self%out_time
  else
    self%out_next = self%out_time*int((self%time+0.5*self%dtime)/self%out_time) + self%out_time
  end if
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace ('1')
  call system ('mkdir -p data/')
  write (filename,'(a,i3.3,"_",i4.4)') &
    'data/'//trim(experiment_name)//'_', self%id, self%iout
  if (io%verbose>=0) &
    print'(a,3x,"time:",f9.5)', 'patches_t%output: file= '//trim(filename), self%time
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  where (self%size > 0.0)
    li = (self%position+0.5*(self%box-self%size))/self%ds
    ui = (self%position+0.5*(self%box+self%size))/self%ds - 1.0
  end where
  write (io%data_unit,'(1p,e14.6,2(2x,3e14.6))') self%time, self%position, self%ds
  write (io%data_unit,'(4(3i7,2x))') int(self%box/self%ds+0.5), self%n, li, ui
  close (io%data_unit)
  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube
  !-----------------------------------------------------------------------------
  open (io%data_unit, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl =sizeof(self%mem(1,1,1,1,1,1))*product(self%n))
  l = self%mesh%li
  u = self%mesh%ui
  do iv=1,self%nv
    if (io%verbose>1) &
      print*,trim(filename),iv, &
              minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)), &
              maxval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1))
    write (io%data_unit,rec=iv) self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)
  end do
  close (io%data_unit)
  !-----------------------------------------------------------------------------
  ! Write a log file, with files written by this rank
  !-----------------------------------------------------------------------------
  write (logname,'(a)') 'data/'//trim(experiment_name)//'.log'
  if (first_time) then
    first_time = .false.
    if (io%master) print*, 'output data file log written to '//trim(logname)
    open (io%log_unit,file=trim(logname),status='unknown')
  else
    open (io%log_unit,file=trim(logname),status='unknown',position='append')
  end if
  write (io%log_unit,'(a)') trim(filename)
  ds = self%box/2.0**io%levelmax
  where (self%size > 0.0)
    li = (self%position+0.5*(self%box-self%size))/ds
    ui = (self%position+0.5*(self%box+self%size))/ds - 1.0
  end where
  write (io%log_unit,'(1p,e14.6,2(2x,3e14.6))') self%time, self%position, self%ds
  write (io%log_unit,'(4(3i7,2x))') int(self%box/self%ds+0.5), self%n, li, ui
  close (io%log_unit)
  self%iout = self%iout+1
  call trace_end
END SUBROUTINE output
!> =============================================================================
!> Find the p1 target and p2 source index ranges that correspond to the overlap
!> region, for upload from p2 to p1 (p1 has lower level than p2)
!> =============================================================================
SUBROUTINE upload_range (p1, p2, l1, u1, l2, u2)
  class(patch_t):: p1
  class(patch_t):: p2
  integer, dimension(3):: l1, u1, l2, u2
  real, dimension(3):: f, pos1, pos2
  !-----------------------------------------------------------------------------
  ! Lower index limits
  !-----------------------------------------------------------------------------
  !print*, p1%id, ' target position, size:', p1%position, p1%size
  !print*, p2%id, ' source position, size:', p2%position, p2%size
  l2 = max(p2%mesh%li, p2%lowest  (p1, p1%mesh%lb))
  u2 = min(p2%mesh%ui, p2%highest (p1, p1%mesh%ub))
  l1 = p1%nearest (p2, l2)
  u1 = p1%nearest (p2, u2)
  !print*,'l1,l2', l1, l2
  !print*,'u1,u2', u1, u2
END SUBROUTINE upload_range
!> =============================================================================
!> Find the p1 target index ranges that correspond to the overlap region, for
!> download from p2 to p1 (p1 has higher or the same level as p2)
!> =============================================================================
SUBROUTINE dnload_range (p1, p2, l1, u1)
  class(patch_t):: p1
  class(patch_t):: p2
  integer, dimension(3):: l1, u1
  real, dimension(3):: f, pos1, pos2, dist
  !-----------------------------------------------------------------------------
  ! Index limits
  !-----------------------------------------------------------------------------
  l1 = p1%mesh%lb
  u1 = p1%mesh%ub
  if (p2%id>1) then
    ! FIXME
    !dist = p2%distance (p1)     ! positive if p2 (source) to the right of p1 (target)
    where (dist >=0)
      l1 = max(l1,p1%nearest(p2,p2%mesh%lb))
      u1 = p1%mesh%ub
    else where
      l1 = p1%mesh%lb
      u1 = min(u1,p1%nearest(p2,p2%mesh%ub))
    end where
  end if
END SUBROUTINE dnload_range
!> =============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%mesh%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%mesh%lb and larger than self%mesh%ub.  For the root
!> grid patch, the index is in the range (li-1,ui) inside the interior.
!> =============================================================================
FUNCTION index (self, other, oindex, f)
  class(patch_t):: self, other
  real, dimension(3):: position, f, pos, pos_lb
  integer, dimension(3):: index, oindex
  !.............................................................................
  pos = other%position + (oindex-other%offset)*other%ds         ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer 
  index = f                                                     ! smaller integer
  f = f - index                                                 ! remainder
END FUNCTION index
!> =============================================================================
!> The nearest index position in the self patch of the point in another patch.
!> The cost of this operation is 4 plus and 2 multiply per direction = 18 flops
!> =============================================================================
FUNCTION nearest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset)*other%ds          ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer 
  out = nint(f)                                                 ! nearest integer
END FUNCTION nearest
!> =============================================================================
!> The lowest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION lowest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset-0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer 
  out = nint(f)                                                 ! lowest integer
END FUNCTION lowest
!> =============================================================================
!> The highest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION highest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset+0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer 
  out = nint(f)                                                 ! highest integer
END FUNCTION highest

!> =============================================================================
!> Signed distance between the centers of two patches in the periodic box
!> =============================================================================
SUBROUTINE courant_condition (self)
  class(patch_t):: self
  !.............................................................................
  if (self%dt_fixed > 0.0) then
    self%dtime = self%dt_fixed
    self%courant = self%dtime*self%u_max/minval(self%ds)
  else if (self%u_fixed > 0.0) then
    self%dtime = self%courant*minval(self%ds)/self%u_fixed
  else
    self%dtime = self%courant*minval(self%ds)/self%u_max
  endif
  if (io%verbose>4) &
    print'(a,i7,f12.6,f8.4)', 'courant:', self%id, self%dtime, self%courant
END SUBROUTINE courant_condition

!> =============================================================================
!> Signed distance between the centers of two patches in the periodic box
!> =============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(patch_t):: self
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  real(8):: out(3)
  !.............................................................................
  select type (task)
  class is (patch_t)
  patch => task
  end select
  out = self%position-patch%position
  out = modulo(out+0.5_8*self%box,max(self%box,1e-30))-0.5_8*self%box
END FUNCTION distance

!> =============================================================================
!> Check if there is overlap, using the full size (including guard zones), since
!> the guard zones may be used in interpolations for qualified time slices. The
!> definition of overlap also MUST be symmetric, since it is being used to define
!> neighbor relations.
!> =============================================================================
FUNCTION overlaps (self, task)
  class(patch_t):: self
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  logical:: overlaps
  !.............................................................................
  select type (task)
  class is (patch_t)
  patch => task
  end select
  ! FIXME
  overlaps = all(abs(distance(self,task)) <= 0.5000000001_8*(self%gsize+patch%gsize))
END FUNCTION overlaps

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%set')
  self%status = ior(self%status, bit)
!  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%clear')
  self%status = iand(self%status, not(bit))
!  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%is_set')
  is_set = (iand(self%status, bit) /= 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (patch_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('patch_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION ahead_of (self, task)
  class (patch_t):: self
  class (task_t):: task
  !.............................................................................
  ahead_of = self%time > task%time + task%grace
!  print*,'ahead_of: ', self%time, patch%time, patch%grace, ahead_of
END FUNCTION
!> =============================================================================
FUNCTION fmaxval (self, f)
  real:: fmaxval
  class(patch_t):: self
  real, dimension(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn):: f
  fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui,self%mesh(2)%li:self%mesh(2)%ui,self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION
!> =============================================================================
FUNCTION fminval (self, f)
  real:: fminval
  class(patch_t):: self
  real, dimension(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn):: f
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui,self%mesh(2)%li:self%mesh(2)%ui,self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION

END MODULE
