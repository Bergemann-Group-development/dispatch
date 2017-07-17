!> =============================================================================
!>  download_link: takes care of downloads to link%task
!>           same: called for patches on the same level
!>         differ: called for patches on different levels
!> =============================================================================
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
  use interpolate_mod
  implicit none
  private
  type :: dh_t
    real(kind=KindScalarVar), dimension(:), pointer :: target => null()
    real(kind=KindScalarVar), dimension(:), pointer :: source => null()
    real(kind=KindScalarVar), dimension(:), pointer :: my_d   => null()
    real(kind=KindScalarVar), dimension(:), pointer :: nb_d   => null()
    class(patch_t), pointer :: src_patch
    integer :: iv
    real    :: wgt(3)
    logical :: set = .false.
    integer :: lvl = 0
    type (dh_t), pointer :: head => null()
    type (dh_t), pointer :: next => null()
    type (dh_t), pointer :: tail => null()
  end type
  ! --------------------------------------------------------------------
  type, public:: data_hub_t
    type (dh_t), dimension(:), pointer :: hub => null()
    class(link_t), pointer :: my_link => null()
    logical :: first_time = .true.
  contains
    procedure:: add_dh
    procedure:: clear_hub
    procedure:: init
    procedure:: map_shell
    procedure:: map_nbors
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Initialise downloading and interpolating options.
!> =============================================================================
SUBROUTINE init (self, link)
  class(data_hub_t):: self
  class(link_t), pointer :: link
  ! ....................................................................
  class(patch_t), pointer :: patch
  type (dh_t),    pointer :: dh
  !---------------------------------------------------------------------
  call trace%begin ('data_hub_t%init')

  if (self%first_time) then
    self%my_link => link
    patch => task2patch(link%task)
    allocate(self%hub((product(patch%gn)-product(patch%n))*patch%nv))
    call self%map_shell
    self%first_time = .false.
  end if

  call trace%end
END SUBROUTINE init

SUBROUTINE update (self, link)
  class(data_hub_t):: self
  class(link_t), pointer :: link
  ! ....................................................................
  class (patch_t), pointer :: patch
  type  (dh_t)   , pointer :: dh
  integer :: i, ij, jt(2), id
  real    :: pt(2), p(3)
  real(kind=KindScalarVar) :: j(16), d(16)
  integer, save :: itimer = 0
  ! --------------------------------------------------------------------
  call trace%begin ('data_hub_t%update', itimer=itimer)
  if (self%first_time) call self%init(link)                             ! map data shell if it is the first time
  patch => task2patch(self%my_link%task)
  id = -1
  do i=1,size(self%hub)
    if (self%hub(i)%set) then
      if ((id<0).or.self%hub(i)%src_patch%id/=id) then
        call self%hub(i)%src_patch%time_interval (patch%time, jt, pt)
      end if
      if (associated(self%hub(i)%next)) then                            ! multiple entries
        dh => self%hub(i)%head
        p = self%hub(i)%wgt
        if (patch%pervolume(self%hub(i)%iv)) then
          do ij = 1, 8
            j(ij)   = dh%source(jt(1))/dh%nb_d(jt(1))
            j(ij+8) = dh%source(jt(2))/dh%nb_d(jt(2))
            dh => dh%next
          end do
          self%hub(i)%target(patch%it) = (&
          pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*j(1) + &
                                (    p(1))*j(2))   + &
                     (   p(2))*((1.0-p(1))*j(3)    + &
                                (    p(1))*j(4)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*j(5)    + &
                                (    p(1))*j(6))   + &
                     (   p(2))*((1.0-p(1))*j(7)    + &
                                (    p(1))*j(8)))) + &
   pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*j(9)    + &
                                (    p(1))*j(10))  + &
                     (   p(2))*((1.0-p(1))*j(11)   + &
                                (    p(1))*j(12))) + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*j(13)   + &
                                (    p(1))*j(14))  + &
                     (   p(2))*((1.0-p(1))*j(15)   + &
                                (    p(1))*j(16)))))*self%hub(i)%my_d(patch%it)
        else
          do ij = 1, 8
            j(ij)   = dh%source(jt(1))
            j(ij+8) = dh%source(jt(2))
            dh => dh%next
          end do
          if (patch%unsigned(self%hub(i)%iv)) then
            self%hub(i)%target(patch%it) = exp( &
     pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*log(j(1 ))    + &
                                  (    p(1))*log(j(2 )))   + &
                       (   p(2))*((1.0-p(1))*log(j(3 ))    + &
                                  (    p(1))*log(j(4 ))))  + &
            (   p(3))*((1.-p(2))*((1.0-p(1))*log(j(5 ))    + &
                                  (    p(1))*log(j(6 )))   + &
                       (   p(2))*((1.0-p(1))*log(j(7 ))    + &
                                  (    p(1))*log(j(8 ))))) + &
     pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*log(j(9 ))    + &
                                  (    p(1))*log(j(10)))   + &
                       (   p(2))*((1.0-p(1))*log(j(11))    + &
                                  (    p(1))*log(j(12))))  + &
            (   p(3))*((1.-p(2))*((1.0-p(1))*log(j(13))    + &
                                  (    p(1))*log(j(14)))   + &
                       (   p(2))*((1.0-p(1))*log(j(15))    + &
                                  (    p(1))*log(j(16))))))
          else
            self%hub(i)%target(patch%it) =  &
     pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*j(1 )    + &
                                  (    p(1))*j(2 ))   + &
                       (   p(2))*((1.0-p(1))*j(3 )    + &
                                  (    p(1))*j(4 )))  + &
            (   p(3))*((1.-p(2))*((1.0-p(1))*j(5 )    + &
                                  (    p(1))*j(6 ))   + &
                       (   p(2))*((1.0-p(1))*j(7 )    + &
                                  (    p(1))*j(8 )))) + &
     pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*j(9 )    + &
                                  (    p(1))*j(10))   + &
                       (   p(2))*((1.0-p(1))*j(11)    + &
                                  (    p(1))*j(12)))  + &
            (   p(3))*((1.-p(2))*((1.0-p(1))*j(13)    + &
                                  (    p(1))*j(14))   + &
                       (   p(2))*((1.0-p(1))*j(15)    + &
                                  (    p(1))*j(16))))
          end if
        end if
      else                                                              ! single entry
        if (patch%unsigned(self%hub(i)%iv)) then
          self%hub(i)%target(patch%it) = exp(    &
          pt(1)*log(self%hub(i)%source(jt(1))) + &
          pt(2)*log(self%hub(i)%source(jt(2))))
        else
          self%hub(i)%target(patch%it) =    &
          pt(1)*self%hub(i)%source(jt(1)) + &
          pt(2)*self%hub(i)%source(jt(2))
        end if
      end if
    end if
    if (isnan(self%hub(i)%target(patch%it))) then
    print *, i, self%hub(i)%iv, associated(self%hub(i)%next)
    error stop "nan in target"
    end if
  end do

  call trace%end
END SUBROUTINE update

SUBROUTINE map_shell(self)
  class(data_hub_t) :: self
  ! ....................................................................
  class(patch_t), pointer :: patch, source
  class(link_t), pointer :: link,nbor
  integer :: l(3), u(3), it, iv
  integer :: ix, iy, iz
  ! --------------------------------------------------------------------

  patch => task2patch(self%my_link%task)
  it = 1
  ! lower x boundary
  l(1) = patch%mesh(1)%lb
  u(1) = patch%mesh(1)%lo
  l(2) = patch%mesh(2)%lb
  u(2) = patch%mesh(2)%ub
  l(3) = patch%mesh(3)%lb
  u(3) = patch%mesh(3)%ub
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do
  ! upper x boundary
  l(1) = patch%mesh(1)%uo
  u(1) = patch%mesh(1)%ub
  l(2) = patch%mesh(2)%lb
  u(2) = patch%mesh(2)%ub
  l(3) = patch%mesh(3)%lb
  u(3) = patch%mesh(3)%ub
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do
  ! lower y boundary
  l(1) = patch%mesh(1)%li
  u(1) = patch%mesh(1)%ui
  l(2) = patch%mesh(2)%lb
  u(2) = patch%mesh(2)%lo
  l(3) = patch%mesh(3)%lb
  u(3) = patch%mesh(3)%ub
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do
  ! upper y boundary
  l(1) = patch%mesh(1)%li
  u(1) = patch%mesh(1)%ui
  l(2) = patch%mesh(2)%uo
  u(2) = patch%mesh(2)%ub
  l(3) = patch%mesh(3)%lb
  u(3) = patch%mesh(3)%ub
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do
    ! lower z boundary
  l(1) = patch%mesh(1)%li
  u(1) = patch%mesh(1)%ui
  l(2) = patch%mesh(2)%li
  u(2) = patch%mesh(2)%ui
  l(3) = patch%mesh(3)%lb
  u(3) = patch%mesh(3)%lo
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do
  ! upper z boundary
  l(1) = patch%mesh(1)%li
  u(1) = patch%mesh(1)%ui
  l(2) = patch%mesh(2)%li
  u(2) = patch%mesh(2)%ui
  l(3) = patch%mesh(3)%uo
  u(3) = patch%mesh(3)%ub
  do iz=l(3), u(3)
    do iy = l(2), u(2)
      do ix = l(1), u(1)
        do iv = 1, patch%nv
          call self%map_nbors(patch, it, iv, ix, iy, iz)
          if (patch%pervolume(iv)) self%hub(it)%my_d => patch%mem(ix,iy,iz,patch%i_d,:,1)
          it = it + 1
        end do
      end do
    end do
  end do

END SUBROUTINE map_shell

SUBROUTINE add_dh (self, it, source, iv, x, y, z, overwrite)
  class (data_hub_t)      :: self
  class(patch_t), pointer :: source
  integer                 :: it, iv, x, y, z
  logical                 :: overwrite
  ! ....................................................................
  type(dh_t), pointer :: dh, this
  ! --------------------------------------------------------------------
  if (overwrite) then
    call self%clear_hub(it)
  end if
  if (associated(self%hub(it)%head)) then
    allocate(dh)
    dh%head => self%hub(it)
    self%hub(it)%tail%next => dh
    self%hub(it)%tail => dh
    dh%source => source%mem(x,y,z,iv,:,1)
    if (source%pervolume(iv)) dh%nb_d => source%mem(x,y,z,source%i_d,:,1)
  else
    self%hub(it)%head => self%hub(it)
    self%hub(it)%tail => self%hub(it)
    self%hub(it)%source => source%mem(x,y,z,iv,:,1)
    if (source%pervolume(iv)) self%hub(it)%nb_d => source%mem(x,y,z,source%i_d,:,1)
    self%hub(it)%src_patch => source
  end if
END SUBROUTINE add_dh

SUBROUTINE clear_hub (self, it)
  class(data_hub_t)    :: self
  integer :: it
  ! ....................................................................
  type(dh_t), pointer :: this, next
  if (associated(self%hub(it)%head%next)) then
    this => self%hub(it)%head%next
    do while(associated(this))
      next => this%next
      deallocate(this)
      this => next
      next => this%next
    end do
    nullify(this, next)
  else
    nullify(self%hub(it)%head)
  end if
  nullify (self%hub(it)%src_patch)
  self%hub(it)%set = .false.
END SUBROUTINE clear_hub

 
FUNCTION pos2cell (self, pos, iv) result (out)
  class(patch_t):: self
  real(8)           :: pos(3)
  integer           :: iv, out(3)
  ! ....................................................................
  integer           ::  ii
  real              :: p(3)
  real              :: h(3)
  real, parameter   :: eps=0.001
  !---------------------------------------------------------------------
  h = self%index_stagger(iv)
  p = pos - self%position                               ! relative position
  p = modulo(p+0.5_8*self%mesh%b,self%mesh%b) &         ! wrap
              -0.5_8*self%mesh%b
  p = p/self%mesh%d + self%mesh%o - h                   ! float index
  p = p + eps                                           ! ensure same points ok
  out = p                                                 ! integer part
END FUNCTION pos2cell

LOGICAL FUNCTION inside(self, cell)
  class(patch_t):: self
  integer       :: cell(3)
  !.....................................................................
  integer :: ii
  logical :: ok(3)
  ! --------------------------------------------------------------------
  inside = .false.
  ok(:) = .false.
  do ii=1,3
    if ((self%li(ii)<=cell(ii)).and.(self%ui(ii)>=cell(ii))) ok(ii) = .true.
  end do
  if (all(ok)) inside = .true.
END FUNCTION inside

SUBROUTINE map_nbors(self,patch, it, iv, ix, iy, iz)
  class(data_hub_t) :: self
  class(patch_t), pointer :: patch, source
  integer :: it, iv, ix, iy, iz
  ! ....................................................................
  class(link_t), pointer :: link,nbor
  integer :: nbpos(3), cpos(3)
  real(8) :: pos(3)
  integer :: j(3), itr
  real    :: p(3)
  ! --------------------------------------------------------------------
  cpos = [ix,iy,iz]
  pos = patch%myposition(cpos, iv)
  self%hub(it)%target => patch%mem(ix,iy,iz,iv,:,1)
  nbor => self%my_link%nbor
  do while (associated(nbor))
    source => task2patch(nbor%task)
    if (source%is_clear(bits%no_download)) then
      if (abs(source%level-patch%level) < 2) then
        nbpos = pos2cell(source, pos, iv)
        if (inside(source,nbpos)) then
          if (patch%level == source%level) then                 ! same level = no interpolation
            if (.not.self%hub(it)%set) then
              self%hub(it)%src_patch => source
              self%hub(it)%source => source%mem &
              (nbpos(1),nbpos(2),nbpos(3),iv,:,1)
              self%hub(it)%lvl = source%level
              self%hub(it)%iv  = iv
              self%hub(it)%head => self%hub(it)
              self%hub(it)%tail => self%hub(it)
              self%hub(it)%set = .true.
            else if (source%level > self%hub(it)%lvl) then
              call self%clear_hub(it)
              self%hub(it)%src_patch => source
              self%hub(it)%source => source%mem &
              (nbpos(1),nbpos(2),nbpos(3),iv,:,1)
              self%hub(it)%lvl = source%level
              self%hub(it)%iv  = iv
              self%hub(it)%head => self%hub(it)
              self%hub(it)%tail => self%hub(it)
              self%hub(it)%set = .true.
            end if
          else                                                  ! different level = interpolation
            if (.not.self%hub(it)%set) then
              call source%index_space(pos, iv, j, p)
              self%hub(it)%wgt = p
              call self%add_dh(it,source,iv,j(1)  ,j(2)  ,j(3)  ,.false.)
              call self%add_dh(it,source,iv,j(1)+1,j(2)  ,j(3)  ,.false.)
              call self%add_dh(it,source,iv,j(1)  ,j(2)+1,j(3)  ,.false.)
              call self%add_dh(it,source,iv,j(1)+1,j(2)+1,j(3)  ,.false.)
              call self%add_dh(it,source,iv,j(1)  ,j(2)  ,j(3)+1,.false.)
              call self%add_dh(it,source,iv,j(1)+1,j(2)  ,j(3)+1,.false.)
              call self%add_dh(it,source,iv,j(1)  ,j(2)+1,j(3)+1,.false.)
              call self%add_dh(it,source,iv,j(1)+1,j(2)+1,j(3)+1,.false.)
              self%hub(it)%lvl = source%level
              self%hub(it)%iv  = iv
              self%hub(it)%set = .true.
            else
              if (source%level > self%hub(it)%lvl) then
                call source%index_space(pos, iv, j, p)
                self%hub(it)%wgt = p
                call self%add_dh(it,source,iv,j(1)  ,j(2)  ,j(3)  ,.true.)
                call self%add_dh(it,source,iv,j(1)+1,j(2)  ,j(3)  ,.false.)
                call self%add_dh(it,source,iv,j(1)  ,j(2)+1,j(3)  ,.false.)
                call self%add_dh(it,source,iv,j(1)+1,j(2)+1,j(3)  ,.false.)
                call self%add_dh(it,source,iv,j(1)  ,j(2)  ,j(3)+1,.false.)
                call self%add_dh(it,source,iv,j(1)+1,j(2)  ,j(3)+1,.false.)
                call self%add_dh(it,source,iv,j(1)  ,j(2)+1,j(3)+1,.false.)
                call self%add_dh(it,source,iv,j(1)+1,j(2)+1,j(3)+1,.false.)
                self%hub(it)%lvl = source%level
                self%hub(it)%set = .true.
              end if
            end if
          end if
        end if
      end if
    end if
    nbor => nbor%next
  end do
  
END SUBROUTINE map_nbors

END MODULE data_hub_mod
