!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Default variables
    ! --------------------------------------------------------------------------
    logical:: first_lower=.true.
    logical:: first_upper=.true.
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real:: cool, hot
    logical:: wall(6)
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: radiation
  end type
  integer, save:: verbose=1
  logical, save:: smudge_boundary = .false.
  logical, save:: scale = .true.
  logical, save:: deriv_bnd = .false.
  real, save:: d_bot=0.0, e_bot=0.0, ee_top=0.0
  real, save:: bottom_bc_rate=0.01
  real, save:: damp_rate=0.01
  real, save:: r_sponge = 0.5
  real, save:: pb_fact = 1.0
  integer, save:: divb_iter = 4
  ! --- lower_damp parameters
  real, save:: beta_lim   = 1e+2
  real, save:: spot_lower = 0.0
  real, save:: pb_spot    = 0.0
  real, save:: t_bdry     = 3e+0
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t)              :: self
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer :: d, e
  real, dimension(:,:,:,:), pointer :: src, rk
  real :: t, d0
  integer:: m(4)
  logical, save:: wall(6)= [.false.,.false.,.false.,.false.,.true.,.false.]
  real, save :: hot_temp = 1000, cool_temp = 10
  real, save :: cool = 0.0, hot = 0.0
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /boundary_params/ cool_temp, hot_temp, wall
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (io%master) write (*, boundary_params)
    !d_bot = d_bot/scaling%d
    !e_bot = e_bot/scaling%e
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  if (.not. io%master) verbose=0

    allocate(d(1,1,2))
    allocate(e(1,1,2))
    allocate(src(1,1,2, eos%n_lambda))
    allocate(rk(1,1,2, eos%n_lambda))
    d = 1.0/scaling%d
    e(:,:,1) = (cool_temp/scaling%temp)/(eos%gamma-1.0)
    e(:,:,2) = (hot_temp/scaling%temp)/(eos%gamma-1.0)
    m = shape(rk)
    call eos%lookup_table (m, d=d, e=e, rk=rk, src=src)
    cool = 0.0      
    hot  = src(1,1,2,1) - src(1,1,1,1)          
    deallocate(d,e,src, rk)
  self%cool = cool
  self%hot  = hot
  self%wall = wall
  call self%standard_boundaries_t%init()
  call trace%end ()
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem2
  real(8):: dtime
  integer:: iv, n
  ! ----------------------------------------------------------------------------
  n = size(mem,4)
  do iv=1, n
    mem2 => mem(:,:,:,iv)
    call self%zero_gradient(mem2, iv)
  end do

END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
END SUBROUTINE derivative

SUBROUTINE radiation (self, mem, src, rk, tt, dir, mu)
  class(boundaries_t) :: self
  real, dimension(:,:,:,:), pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt,mem2
  real   :: mu(3)
  !.....................................................................
  real   :: w, tau, dr
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i3, i3b, i2, i1, n(4), iv, dir, li, ui
  integer:: lbit(3), ubit(3)
  !-----------------------------------------------------------------------------
  n=shape(mem)
  lbit = [bits%xl, bits%yl, bits%zl]        ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]        ! upper boundary bits -- convenience
  dr = self%mesh(dir)%d/abs(mu(dir))
  do iv=1,n(4)
    mem2 => src(:,:,:,iv)
    call self%zero_gradient(mem2, 1)
    mem2 => rk(:,:,:,iv)
    call self%zero_gradient(mem2, 1)
  end do
  if (self%is_set(lbit(3)) ) then
    dr  = self%mesh(3)%d/abs(mu(3))
    li  = self%mesh(3)%li
    i3 =  self%mesh(3)%lo
    i3b = self%mesh(3)%lb
    if (self%wall(5)) then
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
        do i1=self%mesh(1)%lb,self%mesh(1)%ub
          do iv=1, n(4)
            tau = dr*0.5*(rk(i1,i2,li,iv)+rk(i1,i2,i3,iv))
            mem(i1,i2,i3b:i3,iv) = self%hot*exp(-tau)
          end do
        end do
      end do
    else
      dr  = self%mesh(3)%d/abs(mu(3))
      ui  = self%mesh(3)%ui
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
        do i1=self%mesh(1)%lb,self%mesh(1)%ub
          do iv=1, n(4)
            tau = dr*0.5*(rk(i1,i2,ui,iv)+rk(i1,i2,i3,iv))
            mem(i1,i2,i3b:i3,iv) = self%cool*exp(-tau)
          end do
        end do
      end do
    end if
  end if
  if (self%is_set(ubit(3)) ) then
    i3  = self%mesh(3)%uo
    i3b = self%mesh(3)%ub
    if (self%wall(6)) then
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i2,i3:i3b,iv) = self%hot*rk(i1,i2,i3:i3b,iv)*dr
      end do
      end do
      end do
    else
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i2,i3:i3b,iv) = self%cool*rk(i1,i2,i3:i3b,iv)
      end do
      end do
      end do
    end if
  end if
  !
  if (self%is_set(lbit(2)) ) then
    i3  = self%mesh(2)%lo
    i3b = self%mesh(2)%lb
    if (self%wall(3)) then
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3b:i3,i2,iv) = self%hot*rk(i1,i3b:i3,i2,iv)*dr
      end do
      end do
      end do
    else
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3b:i3,i2,iv) = self%cool*rk(i1,i3b:i3,i2,iv)
      end do
      end do
      end do
    end if
  end if
  if (self%is_set(ubit(2)) ) then
    i3  = self%mesh(2)%uo
    i3b = self%mesh(2)%ub
    if (self%wall(4)) then
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3:i3b,i2,iv) = self%hot*rk(i1,i3:i3b,i2,iv)*dr
      end do
      end do
      end do
    else
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3:i3b,i2,iv) = self%cool*rk(i1,i3:i3b,i2,iv)
      end do
      end do
      end do
    end if
  end if
  !
  if (self%is_set(lbit(1)) ) then
    i3  = self%mesh(1)%lo
    i3b = self%mesh(1)%lb
    if (self%wall(1)) then
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3b:i3,i1,i2,iv) = self%hot*rk(i3b:i3,i1,i2,iv)*dr
      end do
      end do
      end do
    else
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3b:i3,i1,i2,iv) = self%cool*rk(i3b:i3,i1,i2,iv)
      end do
      end do
      end do
    end if
  end if  
  if (self%is_set(ubit(1)) ) then
    i3  = self%mesh(1)%uo
    i3b = self%mesh(1)%ub
    if (self%wall(2)) then
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3:i3b,i1,i2,iv) = self%hot*rk(i3:i3b,i1,i2,iv)*dr
      end do
      end do
      end do
    else
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3:i3b,i1,i2,iv) = self%cool*rk(i3:i3b,i1,i2,iv)
      end do
      end do
      end do
    end if
  end if  
END SUBROUTINE radiation

END MODULE boundaries_mod
