module ray_patch_mod

  !USE rt_integral_mod
  !USE radau_mod
  !USE boundaries_mod
  !USE eos_mod
  !USE list_mod
  USE halo_mod
  USE gpatch_mod
  USE link_mod
  USE task_mod
  !USE index_mod
  !USE mesh_mod
  USE io_mod
  !USE aux_mod
  !USE timer_mod
  USE scaling_mod
  USE trace_mod
  !USE mpi_mod
  !USE omp_mod
  USE omp_lock_mod
  USE shared_mod
  !USE omp_timer_mod
  !USE math_mod
  USE bits_mod
  USE os_mod
  !USE units_mod
  !USE scratch_mod
  USE ray_mod, only:ray_type
  USE m3d_patch_mod, only: m3d_patch_t, vec, wmu, long_vec, long_wmu
  implicit none
  private

  !-----------------------------------------------------------------------------

  type, public, extends(m3d_patch_t) :: ray_patch_t
    type(ray_type) :: ray
    real(8) :: vec(3), wmu
    integer :: dir=0, fs, fe, i_out=0
    real(4), dimension(:,:,:), pointer, contiguous :: vproj
  contains
    procedure :: pre_init => ray_patch_pre_init
    procedure :: init => ray_patch_init
    procedure :: set_freq_indices
    procedure :: check_topbot
    procedure :: is_upstream => ray_is_upstream
    procedure :: comp_vproj
  end type ray_patch_t

  !-----------------------------------------------------------------------------
  integer, save:: verbose=0
    
contains

!===============================================================================
SUBROUTINE ray_patch_pre_init (self)
  class(ray_patch_t), target :: self
  ! ----------------------------------------------------------------------------
  call self%m3d_patch_t%pre_init()
  
  if (self%i_mu > self%n_omega) then
    self%out_ray = .true.
    self%i_out = self%i_mu - self%n_omega 
    self%vec = long_vec(:,self%i_out) 
    self%wmu = long_wmu(self%i_out)
    self%dir = 1
  elseif(self%i_mu > 0) then
    self%vec = vec(:,self%i_mu)
    self%wmu = wmu(self%i_mu)
    self%dir = merge(1, -1, vec(3,self%i_mu) > 0)
  endif

  if (self%dir < 0) then
    self%boundaries%lower = .false.
  elseif (self%dir > 0) then
    self%boundaries%upper = .false.
  endif

  call self%set_freq_indices()
END SUBROUTINE ray_patch_pre_init

!===============================================================================
SUBROUTINE set_freq_indices (self)
  class(ray_patch_t) :: self
  integer :: fs, fe, nv
  !-----------------------------------------------------------------------------
  nv = self%sp%nnu / self%n_nu
  fs = (self%i_nu-1) * nv + 1
  if (self%i_nu == self%n_nu) nv = nv + mod(self%sp%nnu, self%n_nu)
  fe = fs + nv - 1

  self%fs = fs
  self%fe = fe
END SUBROUTINE set_freq_indices

!===============================================================================
SUBROUTINE comp_vproj (self, vproj)
  class(ray_patch_t) :: self
  real(4), dimension(:,:,:), pointer, contiguous :: vproj
  real(4), dimension(:,:,:), pointer :: d, px, py, pz
  integer :: li(3), ui(3)
  !-----------------------------------------------------------------------------
  li = self%patch%li
  ui = self%patch%ui

  d  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%d, self%it,1)
  px => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%px,self%it,1)
  py => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%py,self%it,1)
  pz => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%pz,self%it,1)

  if (self%out_ray) then
    vproj = (-self%vec(1)*px - self%vec(2)*py + self%vec(3)*pz) * scaling%u / d
  else
    vproj = (self%vec(1)*px + self%vec(2)*py + self%vec(3)*pz) * scaling%u / d
  endif
END SUBROUTINE comp_vproj

!===============================================================================
SUBROUTINE ray_patch_init (self)
  class(ray_patch_t) :: self
  ! ----------------------------------------------------------------------------
  call self%m3d_patch_t%init()
  call self%check_topbot()

  if (self%i_mu > 0) then
    call os%allocate(self%vproj, [self%n(1), self%n(2), self%n(3)])
  endif
END SUBROUTINE ray_patch_init

!===============================================================================
!> Determine if a task is is upstream.  
!===============================================================================
logical FUNCTION ray_is_upstream (self, other) result(is_upstream)
  class(ray_patch_t):: self, other
  real(8):: ds1(3), ds2(3)
  real:: sinth
  real, parameter:: eps=1e-4
  integer:: is1(3), is2(3), is3(3), is4(3), i, lbit(3), ubit(3)
  logical:: lower, upper, upstream
  !-----------------------------------------------------------------------------
  ! If self is at a boundary and other points inward, or the task has no RT,
  ! then the task is outside and is not needed
  !-----------------------------------------------------------------------------
  if (.not.self%on) then
    is_upstream=.false.
    return
  end if
  !
  if (self%i_mu /= other%i_mu) then
    is_upstream=.false.
    return
  end if
  !
  call trace%begin ('ray_patch_t%ray_is_upstream')
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.other%on) then
    is_upstream = .false.
    call trace%end ()
    return
  end if
  !
  is_upstream=.true.
  ds1 = self%position-other%position
  where (self%periodic)
    ds1 = modulo(ds1 + 0.5_8*self%box, self%box) - 0.5_8*self%box
  end where
  ds1 = ds1/self%box                                                            ! relative task offsets
  is1 = 0; is2=0
  do i=1,3
    if (ds1(i)      > +eps) is1(i)=1
    if (self%vec(i) > +eps) is2(i)=1
    if (ds1(i)      < -eps) is1(i)=-1
    if (self%vec(i) < -eps) is2(i)=-1
  end do
  if (all(is1*is2 <= 0)) is_upstream = .false.
  ! ----------------------------------------------------------------------------
  ! If the nbor lies outside the boundaries it is not needed
  ! ----------------------------------------------------------------------------
  lbit = [bits%xl, bits%yl, bits%zl]
  ubit = [bits%xu, bits%yu, bits%zu]
  do i=1,3
    lower = self%boundaries%is_set(lbit(i))
    upper = self%boundaries%is_set(ubit(i))
    if (upper .and. ds1(i)==1) then
      is_upstream = .false.
    end if
    if (lower .and. ds1(i)==-1) then
      is_upstream = .false.
    end if
  end do
  !-----------------------------------------------------------------------------
  ! If the task does not have RT it is not needed
  !-----------------------------------------------------------------------------
  if (.not.other%on) then
    is_upstream = .false.
    call trace%end ()
    return
  end if
  call trace%end ()
END FUNCTION ray_is_upstream

!===============================================================================
SUBROUTINE check_topbot (self)
  class(ray_patch_t):: self
  !.............................................................................
  if (self%patch%boundaries%upper .and. self%dir < 0) then
    self%boundaries%upper = .true.
  elseif (self%patch%boundaries%lower .and. self%dir > 0) then
    self%boundaries%lower = .true.
  endif
END SUBROUTINE

end module ray_patch_mod
