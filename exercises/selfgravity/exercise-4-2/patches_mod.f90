!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Generic module for MHD patch.  To be overloaded by an MHD layer and a timestep
!> layer.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patches_mod
  USE trace_mod
  USE io_mod
  USE bits_mod
  implicit none
  integer, parameter:: ntime=7
  private
  type, public:: patch_t
    integer:: id = 0
    integer:: nt = ntime, nv = 5
    integer:: n(3), ng(3)=2, gn(3), li(3), ui(3)
    integer:: it=1, new=2
    integer:: i, j, k
    integer:: status
    integer, pointer:: iit(:)
    real, dimension(3):: position, size, gsize, box=1, ds, origin
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
    procedure:: index_range     ! returns the index ranges btw two points
    procedure:: indices         ! computes indices and fractional indices
    procedure:: index           ! compute the index for a single point
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: ahead_of
    procedure:: is_relevant_for
    procedure:: download_to
    procedure:: myposition
    procedure:: fmaxval         ! interior max value
    procedure:: fminval         ! interior min value
  end type
  real, public:: end_time=1.0, out_time=0.1
  integer:: id = 0
CONTAINS
!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self, position, size)
  class(patch_t):: self
  integer:: dims(5) = 6
  real:: position(3), size(3)
  logical:: first_time = .true.
  namelist /patch_params/ dims, size
  !.............................................................................
  call trace_begin ('patch_t%init', 3)
  !$omp atomic
  id = id + 1                                   ! increment global
  print*,'patch_t%init: id=', id
  self%id = id                                  ! new patch id
  self%position = position                      ! as per call
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
  self%size = size                              ! patch size
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
  self%li = 1       + self%ng                   ! lower interior index
  self%ui = self%gn - self%ng                   ! upper interior index
  self%nt = dims(4)                             ! number of time slices
  self%nv = dims(5)                             ! number of variables
  self%level = log(maxval(self%box)/minval(self%ds))/log(2.) + 0.1
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
  call trace_end
END SUBROUTINE init

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
!> all cells in the tager if "all_target" is true, and only guard cells if not.
!> We use myposition() to map an index position in the target patch to an
!> absolute position, and then index() to get the integer and floating point
!> position in the source patch.  With that, we make a tri-linear interpolation.
!> This is useful for both guard cell fill-in, and for "prolongation" = creating
!> patches with higher resolution.
!> =============================================================================
SUBROUTINE download_to (self, patch, only_down)
  class(patch_t):: self
  class(patch_t), pointer:: patch
  integer:: only_down
  logical:: all_target
  integer:: ii(3), ix, iy, iz, it, iv, j0(3), j1(3), n, only_down_
  real:: f(3), pos(3)
  !.............................................................................
!  if (self%id==patch%id .and. all(patch%size /= box%size)) return
  if (io%verbose>2) print'(a,2(i3,3f8.3),2(3x,3i4))','index:self,patch:', &
    self%id,self%position,patch%id,patch%position, self%li, self%ui
associate (source=>self%mem, target=>patch%mem, jm=>self%gn)
  it = patch%it
  n = 0
  !-----------------------------------------------------------------------------
  ! Download rules:  1) Always download guard cells.  2) Download all cells if
  ! a) all_cells is set, and the source level is at most one level higher. In
  !    this case (first fill-in), the source level may be lower than target level
  ! b) the level is one level higher
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    all_target = .true.
    if (self%level <= patch%level) all_target = .false.   ! do not upload from below
    if ((self%level-patch%level) > 1) return              ! upload at most 1 level
    if (iv==only_down) then
      if (self%level >= patch%level) cycle                ! only download
      all_target = .true.                                 ! interior + guard cells
      if (io%verbose>0) &
        print'(10(1x,a,i6,3x))', 'patch', patch%id, 'at level', patch%level, &
        'downloading variable', only_down_, 'from patch', self%id, 'level', self%level
    end if
    do iz=1,patch%gn(3)
    do iy=1,patch%gn(2)
    do ix=1,patch%gn(1)
      ii = [ix,iy,iz]
      !---------------------------------------------------------------------------
      ! Check if the target point is the guard zone region (= not inside)
      !---------------------------------------------------------------------------
      if (all_target .or. .not.all((ii >= patch%li) .and. (ii <= patch%ui))) then
        pos = patch%myposition(ii)
        j0 = self%index (pos, f)
        !---------------------------------------------------------------------------
        ! Check that the source in in the internal region (= not guard zone)
        !---------------------------------------------------------------------------
        if (all((j0+f >= self%li-0.1) .and. (j0+f <= self%ui+0.1) .or. (self%n==1))) then
          if (io%verbose>3.and.patch%id==10) print'(a,2(i4,2x,3i4,2x),1p,2(3f9.4,2x),2x,2f11.6," ****")', &
            'index:ii,j0,f:', patch%id, ii,  self%id, j0, f, pos, self%t(it)-patch%t(it)
          j1 = j0+1
          where (j1 > self%ui)
            j0 = j0-1
            j1 = j1-1
            f  = f+1.0
          end where
          where (j0 < self%li)
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
          n = n+1
        end if
      end if
    end do
    end do
    end do
  end do
  if (io%verbose > 1) print'(6(1x,a,i5))', 'patch', patch%id, ' level', patch%level, &
    ' picked up', n, ' points from patch', self%id, ' at level', self%level
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
  out = self%origin + (ii-self%li + 0.5)*self%ds
END FUNCTION

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
  if (self%time < self%out_next) return
  self%out_next = self%out_time*int(self%time/self%out_time) + self%out_time
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
FUNCTION index_range (self, position)
  class(patch_t):: self
  real:: position(3)
  integer:: index_range(3,2)
  integer:: dir
  !.............................................................................
  print '(1x,a,2(3x,3f8.3))',   'positions:', position
  do dir=1,3
    index_range(:,dir) = int(modulo(position(dir)-self%origin, self%box)/self%ds + 1.5)
    index_range(:,dir) = max(1, min(self%gn, index_range(:,dir)))
    print '(1x,a,2(3x,i3,a,i3))', 'index range:', &
      index_range(1,dir), ' - ', index_range(2,dir)
  end do
END FUNCTION index_range
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
FUNCTION index (self, position, f)
  class(patch_t):: self
  real, dimension(3):: position, f
  integer:: index(3)
  !.............................................................................
  f = modulo(position - self%origin + self%li*self%ds, self%box)/self%ds - 0.5 + 1e-5
  index = f
  f = f - index
END FUNCTION index

!> =============================================================================
FUNCTION distance (self, patch)
  class(patch_t):: self, patch
  real:: distance(3)
  distance = self%position-patch%position
  distance = abs(modulo(distance+0.5*self%box,max(self%box,1e-30))-0.5*self%box)
END FUNCTION distance
!> =============================================================================
FUNCTION overlaps (self, patch)
  class(patch_t):: self, patch
  logical:: overlaps
  integer:: idir
  !.............................................................................
  if (minval(self%ds) <= minval(patch%ds)) then                                 ! make sure smallest guard zone is used
    overlaps = all(distance(self,patch) <= 0.5000000001_8*(self%gsize+patch%size)) ! prologation needs at least an additional zone
  else
    overlaps = all(distance(self,patch) <= 0.5000000001_8*(self%size+patch%gsize))
  endif
END FUNCTION
!> =============================================================================
!> This is the relevance criterion if we don't need corner nbors
!> We are checking if self is relevant for patch for guard cells
!> This does not work if patch is embedded (a descendant) of self.
!> =============================================================================
LOGICAL FUNCTION is_relevant_for (self, patch)
  class (patch_t):: self
  class (patch_t), pointer:: patch
  !
  integer      :: dir, overlap
  real(kind=8) :: dist(3), right_distance
  logical      :: face, edge, corner
  !............................................................................
  dist = distance(self,patch)
  overlap = 0
  do dir=1,3
    if (patch%n(dir)==1) cycle ! cycle if direction is passive
    ! this is the right distance
    right_distance = 0.5_8*(self%ds(dir)*self%n(dir) + patch%ds(dir)*patch%n(dir))
    ! inc overlap counter if we are at the right distance. Compare with ds to avoid round-off errors
    if (abs(dist(dir) - right_distance) < 0.1_8*min(self%ds(dir),patch%ds(dir))) &
      overlap = overlap + 1
  enddo

  face = (overlap==1)
  edge = (overlap==2)
  corner = (overlap==3)

  ! For now only check for face-overlap
  is_relevant_for = face
END FUNCTION

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
