!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Generic module for MHD patch.  To be overloaded by an MHD layer and a timestep
!> layer.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patches_mod
  USE trace_mod
  USE io_mod
  USE omp_timer_mod
  USE bits_mod
  implicit none
  integer, parameter:: ntime=7
  private
  type, public:: patch_t
    integer:: id = 0
    integer:: nt = ntime, nv = 5
    integer:: n(3), ng(3)=2, gn(3), li(3), ui(3), lb(3), ub(3)
    integer:: it=1, new=2
    integer:: i, j, k
    integer:: status
    integer, pointer:: iit(:)
    real, dimension(3):: position, size, gsize, box=1, ds, origin, offset
    real, pointer:: t(:)
    real, pointer:: dt(:)
    real:: time=0.0, dtime=0.1, end_time=1.0, grace=0.0
    real:: out_time=0.1, out_next=0.0           ! controls output
    real:: courant=0.2                          ! courant constant
    integer:: istep=0                           ! time step counter
    integer:: iout=0                            ! output counter
    integer:: level
    real, dimension(:,:,:,:,:,:), pointer:: mem
    integer:: verbose = 0
    logical:: is_leaf = .true.
    logical:: mem_allocated = .false.
    character(len=64), pointer:: experiment_name
  contains
    procedure:: init
    procedure:: update => update_patch
    procedure:: solve           ! static solver backup
    procedure:: output
    procedure:: distance        ! return the distance btw two patches
    procedure:: overlaps        ! return true if they overlap (incl. ghost zones)
    procedure:: upload_range    ! returns the index ranges btw two points
    procedure:: dnload_range    ! returns the index ranges btw two points
    procedure:: indices         ! computes indices and fractional indices
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
  end type
  real, public:: end_time=1.0, out_time=0.1
  integer, public:: msplit=2    ! shared with refine_mod and patch_list_mod
  integer:: id = 0
CONTAINS
!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  integer:: dims(5) = 6
  logical:: first_time = .true.
  real:: position(3)
  namelist /patch_params/ dims
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
  dims(4) = self%nt                             ! set by timestep_t:init
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
  end if
  self%experiment_name => experimentname
  self%status=0
  self%end_time=end_time
  self%out_time=out_time
  self%origin = self%position - 0.5*self%size   ! physical origin -- already defined
  self%n  = dims(1:3)                           ! patch dimensions
  self%ds = real(self%size,kind=8) / self%n     ! ds in patch
  where (self%n==1)                             ! for collapsed dimensions
    self%ng = 0                                 ! no guard zones
    self%size = 0.0                             ! zero size
    self%ds = 0.0                               ! first zero
    self%box = 1                                ! make sure it is non-zero
    self%ds = maxval(self%ds)                   ! then set to reasonable val
  end where
  self%gn = self%n + 2*self%ng                  ! including guard cells
  self%gsize = self%size*real(self%gn,kind=8)/self%n ! size incl guard cells
  self%lb = 1                                   ! lower guard zone index
  self%li = 1       + self%ng                   ! lower interior index
  self%ui = self%gn - self%ng                   ! upper interior index
  self%ub = self%gn                             ! upper guard zone index
  self%offset = (self%lb+self%ub)/2.0           ! mid point index offset
  self%nt = dims(4)                             ! number of time slices
  self%nv = dims(5)                             ! number of variables
  self%level = log(maxval(self%box)/minval(self%ds))/log(real(msplit)) + 0.1
  if (io%verbose>3) &
    print*,'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not allocated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
!    if (associated(self%mem)) deallocate (self%mem)
    allocate (self%mem(self%gn(1),self%gn(2),self%gn(3),self%nv,self%nt,2))
    self%mem_allocated = .true.                 ! mark
  end if
  self%t = 0.0                                  ! initial time
  self%iit = 1                                  ! initial memory slots
  if (self%id<=3) call test_patch (self)
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE test_patch (self)
  class(patch_t):: self
  real, dimension(3):: pos, f
  integer, dimension(3):: ii, nn
  integer:: i
  !.............................................................................
  print *, 'test_patch:', self%id
  do i=self%lb(1),self%ub(1)
    pos = self%myposition([i,i,i])
    ii = self%index(self, [i,i,i], f)
    nn = self%nearest(self, [i,i,i])
    print '(i4,2(3f7.3,3i4))', i, pos, ii, f, nn
  end do
END SUBROUTINE test_patch

!> =============================================================================
SUBROUTINE solve (self)
  class(patch_t):: self
  !.............................................................................
END SUBROUTINE solve

!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE update_patch (self)
  class(patch_t):: self
  !.............................................................................
  call trace_begin ('patch_t%update', 3)
  self%dt(self%it) = self%dtime                         ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = self%time                          ! initial time
  self%it = self%new                                    ! update time slot
  self%new = mod(self%new,self%nt)+1                      ! increment / rotate
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
END SUBROUTINE update_patch

!> =============================================================================
!> Download values from the source (self) to the target (patch), trying to get
!> all cells in the target if "all_target" is true, and only guard cells if not.
!> We use myposition() to map an index position in the target patch to an
!> absolute position, and then index() to get the integer and floating point
!> position in the source patch.  With that, we make a tri-linear interpolation.
!> This is useful for both guard cell fill-in, and for "prolongation" = creating
!> patches with higher resolution.
!> =============================================================================
SUBROUTINE download_to (self, patch, only_down, all_cells)
  class(patch_t):: self
  class(patch_t), pointer:: patch
  logical:: all_cells
  integer:: only_down
  logical:: all_target
  integer:: ii(3), ix, iy, iz, it, iv, j0(3), j1(3), ntot, nuse, only_down_
  integer, dimension(3):: l1, u1, l2, u2
  real:: f(3), pos(3), weight
  !.............................................................................
!  if (self%id==patch%id .and. all(patch%size /= box%size)) return
  if (io%verbose>2) print'(a,2(i3,3f8.3),2(3x,3i4))','download_to:self,patch=', &
    self%id,self%position,patch%id,patch%position, self%li, self%ui
associate (source=>self%mem, target=>patch%mem, jm=>self%gn)
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
    if (iv==only_down) then                               ! if the variable is phi
      if (self%level >= patch%level) cycle                ! only download phi from lower levels
      !all_target = .true.                                 ! interior + guard cells
      !if (io%verbose>0) &
      !  print'(10(1x,a,i6,3x))', 'patch', patch%id, 'at level', patch%level, &
      !  'downloading variable', only_down, 'from patch', self%id, 'level', self%level
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
      ! patches p1 and p2, respectively, corresponding to p2%li and p2%ui, of a
      ! region inside p2%li-p2%ui, which is also inside p1%lb-p1%ub. In this 
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
          if (any(j0 >= patch%li .or. j0 <= patch%ui)) cycle
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
        if (all_target .or. .not.all((ii >= patch%li) .and. (ii <= patch%ui))) then
          j0 = self%index (patch, ii, f)
          !---------------------------------------------------------------------------
          ! Check that the source is in the valid region
          !---------------------------------------------------------------------------
          if (all((j0+f >= self%lb-0.01 .and. j0+f <= self%ub+0.01) .or. self%n==1)) then
            j1 = j0+1
            where (j1 > self%ub)
              j0 = j0-1
              j1 = j1-1
              f  = f+1.0
            end where
            where (j0 < self%lb)
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
            if (iv==1) print'("outside",i4,2x,3i4,2x,i4,3f8.3,3(2x,3i4))',patch%id,ix,iy,iz,self%id,j0+f,j0,self%lb,self%ub
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
  integer:: ix, iy, iz, i(3), j(3)
  !.............................................................................
  do iz=self%lb(3),self%ub(3)
  do iy=self%lb(2),self%ub(2)
  do ix=self%lb(1),self%ub(1)
    i = [ix,iy,iz]
    if (all(i >= self%li .and. i <= self%ui)) cycle
    j = modulo(i-self%lb,self%n)+self%lb
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
  do iz=self%lb(3),self%ub(3)                           ! loop over all cells
  do iy=self%lb(2),self%ub(2)
  do ix=self%lb(1),self%ub(1)
    i = [ix,iy,iz]
    if (all(i >= self%li .and. i <= self%ui)) cycle     ! skip interior
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
!> Get guard cells for "patch" using active cells in "self"
SUBROUTINE troels_get_guardcells_for (self, patch)
  class(patch_t):: self
  class(patch_t), pointer:: patch
  !
  integer :: ls(3), us(3), lp(3), up(3), dir, it
  real(kind=8) :: rlb(3), rub(3), origin(3), hs
  integer :: nprint=100
  !
  ! Check if the two patches have same size
  if (maxval(abs(patch%ds-self%ds)/max(patch%ds,self%ds)) < 1e-5) then  ! we are on the same level
    ! Get real lower and upper boundary with possible wrap-around
    rlb = self%position-(self%n/2+1e-5_8)*self%ds ! go slightly below lower boundary
    rub = self%position+(self%n/2+1e-5_8)*self%ds ! go slightly above upper boundary

    ! Troels: if you need a global origin you have to put it in a different variable -- presumably
    ! a globale one (e.g. in patch_list_mod) -- this one is taken.  But it is also obvious, that
    ! the only information needed is the relative position btw self and patch, so this "origin"
    ! is redundant.
    !origin = self%origin                            ! global lower boundary in box
    rlb = modulo(rlb - origin, max(self%box,1e-30)) + origin   ! wrap around to stay inside box
    rub = modulo(rub - origin, max(self%box,1e-30)) + origin
    where(self%n==1)
      rlb=0.
      rub=0.
    endwhere
    if (io%verbose > 4) then
      print *, 'self, patch id :', self%id, patch%id
      print *, 'rlb for self:', rlb
      print *, 'rub for self:', rub
      print *, 'patch positi:', patch%position
      print *, 'patch size  :', patch%size
    endif
    !
    ! check if we can get guard cells in a given direction
    do dir=1,3
      if (patch%n(dir)==1) cycle                                         ! exclude inactive dimensions
      hs = 0.5_8 * patch%size(dir)
      ! check lower boundary
      if (abs(patch%position(dir) - rub(dir)) < hs) then
        ! We are overlapping in this direction, exchange data
        ls = 1; us=self%gn
        ls(dir) = self%ui(dir) - patch%ng(dir) + 1
        us(dir) = self%ui(dir)

        lp = 1; up=patch%gn
        lp(dir) = 1
        up(dir) = patch%li(dir) - 1

        ! copy data
        do it=self%nt,1,-1
          if (self%t(self%iit(it))==patch%time) then
            patch%mem(lp(1):up(1),lp(2):up(2),lp(3):up(3),:,patch%it,:) = &
             self%mem(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,self%iit(it),:)
            exit
          endif
        enddo

        if (io%verbose > 4) then
          print *, 'Overlapping in lower boundary:', dir
          print *, 'Patch (receiver):', lp, up, up-lp+1
          print *, 'Self  (donor)   :', ls, us, us-ls+1
          print *, 'Patch time      :', patch%it, patch%time
          print *, 'Self  time      :', it, self%t(self%iit(it))
          nprint = nprint - 1
          if (nprint==0) stop
        endif
      endif
      ! check upper boundary
      if (abs(rlb(dir)-patch%position(dir)) < hs) then
        ! We are overlapping in this direction, exchange data
        ls = 1; us=self%gn
        ls(dir) = patch%li(dir)
        us(dir) = patch%li(dir) + patch%ng(dir) - 1

        lp = 1; up=patch%gn
        lp(dir) = patch%ui(dir) + 1
        up(dir) = patch%gn(dir)

        ! copy data
        do it=self%nt,1,-1
          if (self%t(self%iit(it))==patch%time) then
            patch%mem(lp(1):up(1),lp(2):up(2),lp(3):up(3),:,patch%it,:) = &
             self%mem(ls(1):us(1),ls(2):us(2),ls(3):us(3),:,self%iit(it),:)
            exit
          endif
        enddo

        if (io%verbose > 4) then
          print *, 'Overlapping in upper boundary:', dir
          print *, 'Patch (receiver):', lp, up, up-lp+1
          print *, 'Self  (donor)   :', ls, us, us-ls+1
          print *, 'Patch time      :', patch%it, patch%time
          print *, 'Self  time      :', it, self%t(self%iit(it))
          nprint = nprint - 1
          if (nprint==0) then
            print *, 'Max prints reached in guard cells. Stopping.'
            stop
          endif
        endif
      endif
    enddo
  endif
END SUBROUTINE troels_get_guardcells_for

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self)
  class(patch_t):: self
  character(len=32) filename, logname
  integer:: iv, li(3)=0, ui(3)=0
  real:: ds(3)
  logical:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! If the out_next time has passed, set a new time = the integer number of
  ! out_time's that have passed, plus one more
  !-----------------------------------------------------------------------------
  if (self%time+0.5*self%dtime < self%out_next) return
  if (self%time==0.0) then
    self%out_next = self%out_time
  else
    self%out_next = self%out_time*int((self%time+0.5*self%dtime)/self%out_time) + self%out_time
  end if
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_t%output')
  call system ('mkdir -p data/')
  write (filename,'(a,i3.3,"_",i4.4)') &
    'data/'//trim(self%experiment_name)//'_', self%id, self%iout
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
associate (l=>self%li, u=>self%ui)
  do iv=1,self%nv
    if (io%verbose>1) &
      print*,trim(filename),iv, &
              minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)), &
              maxval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1))
    write (io%data_unit,rec=iv) self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)
  end do
end associate
  close (io%data_unit)
  !-----------------------------------------------------------------------------
  ! Write a log file, with files written by this rank
  !-----------------------------------------------------------------------------
  write (logname,'(a)') 'data/'//trim(self%experiment_name)//'.log'
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
  l2 = max(p2%li, p2%lowest  (p1, p1%lb))
  u2 = min(p2%ui, p2%highest (p1, p1%ub))
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
  l1 = p1%lb
  u1 = p1%ub
  if (p2%id>1) then
    dist = p2%distance (p1)     ! positive if p2 (source) to the right of p1 (target)
    where (dist >=0)
      l1 = max(l1,p1%nearest(p2,p2%lb))
      u1 = p1%ub
    else where
      l1 = p1%lb
      u1 = min(u1,p1%nearest(p2,p2%ub))
    end where
  end if
END SUBROUTINE dnload_range
!> =============================================================================
SUBROUTINE indices (self, np, p, ii, f)
  class(patch_t):: self
  real:: p(np), f(np,3)
  integer:: np, dir, ii(np,3), i
  !.............................................................................
  do dir=1,3
    f(:,dir) = modulo(p-self%origin(dir), self%box(dir))/self%ds(dir) + 1.5
    ii(:,dir) = max(1, min(self%gn(dir), int(f(:,dir))))
      end do
  f = f-ii
  do i=1,np
    print'(i4,3x,a,3i4,3x,a,3f8.3)', i, 'index:', ii(i,:), 'fraction:', f(i,:)
  end do
END SUBROUTINE indices
!> =============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%lb and larger than self%ub.  For the root
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
FUNCTION distance (self, patch) RESULT (out)
  class(patch_t):: self, patch
  real(8):: out(3)
  !.............................................................................
  out = self%position-patch%position
  out = modulo(out+0.5_8*self%box,max(self%box,1e-30))-0.5_8*self%box
END FUNCTION distance
!> =============================================================================
!> Check if there is overlap, using the full size (including guard zones), since
!> the guard zones may be used in interpolations for qualified time slices. The
!> definition of overlap also MUST be symmetric, since it is being used to define
!> neighbor relations.
!> =============================================================================
FUNCTION overlaps (self, patch)
  class(patch_t):: self, patch
  logical:: overlaps
  !.............................................................................
  overlaps = all(abs(distance(self,patch)) <= 0.5000000001_8*(self%gsize+patch%gsize))
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
LOGICAL FUNCTION ahead_of (self, patch)
  class (patch_t):: self, patch
  !.............................................................................
  ahead_of = self%time > patch%time + patch%grace
!  print*,'ahead_of: ', self%time, patch%time, patch%grace, ahead_of
END FUNCTION
!> =============================================================================
FUNCTION fmaxval (self, f)
  real:: fmaxval
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: f
  fmaxval = maxval(f(self%li(1):self%ui(1),self%li(2):self%ui(2),self%li(3):self%ui(3)))
END FUNCTION
!> =============================================================================
FUNCTION fminval (self, f)
  real:: fminval
  class(patch_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: f
  fminval = minval(f(self%li(1):self%ui(1),self%li(2):self%ui(2),self%li(3):self%ui(3)))
END FUNCTION
END MODULE
