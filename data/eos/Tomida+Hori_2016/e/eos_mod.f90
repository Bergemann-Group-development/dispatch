!> =============================================================================
!> Module using a table prepared by a Python procedure
!> =============================================================================
MODULE eos_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE scaling_mod
  USE units_mod
  private
  type, public:: eos_t
    character(len=64) :: file
    integer           :: ioformat, n1, n2, n3
    integer           :: i_pg, i_tt, i_ss, i_rk, i_src
    real              :: x_min, x_step, y_min, y_step
    real, pointer     :: table(:,:,:)
    logical           :: first_time=.true.
  contains
    procedure         :: init
    procedure         :: init_table
    procedure         :: lookup
  end type
  type(eos_t), public:: eos, eose
CONTAINS

!> =============================================================================
!> Initialize equation-of-state handling
!> =============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  call self%init_table ('/data/eos/Tomida+Hori_2016/e/eos.dat')
  self%x_min = self%x_min - log(scaling%d)
  self%y_min = self%y_min - log(scaling%ee)
END SUBROUTINE init

!> =============================================================================
!> Initialize table
!> =============================================================================
SUBROUTINE init_table (self, file)
  class(eos_t):: self
  character(len=*):: file
  !-----------------------------------------------------------------------------
  if (self%first_time) then
    !$omp critical (io_cr)
    if (self%first_time) then
      self%first_time = .false.
      !-------------------------------------------------------------------------
      ! Read a custom binary file
      !-------------------------------------------------------------------------
      print *, 'init_table: reading '//trim(io%top)//trim(file)
      open (io_unit%data, file=trim(io%top)//trim(file), form='unformatted', status='old')
      read (io_unit%data) self%ioformat, self%n1, self%n2, self%n3
      read (io_unit%data) self%x_min, self%x_step
      read (io_unit%data) self%y_min, self%y_step
      read (io_unit%data) self%i_pg, self%i_tt, self%i_ss, self%i_rk, self%i_src
      if (io%master) &
        print '(a,7i4)', ' init_table: i_pg, i_tt, i_ss, i_rk, i_src =', &
        self%i_pg, self%i_tt, self%i_ss, self%i_rk, self%i_src
      allocate (self%table(self%n1,self%n2,self%n3))
      read (io_unit%data) self%table
      close (io_unit%data)
      !-------------------------------------------------------------------------
      ! Rescale to the current code units
      !-------------------------------------------------------------------------
      self%table(:,:,self%i_pg ) = self%table(:,:,self%i_pg )/scaling%p
      self%table(:,:,self%i_tt ) = self%table(:,:,self%i_tt )/scaling%temp
      if (self%i_rk>0) &
      self%table(:,:,self%i_rk ) = self%table(:,:,self%i_rk )/scaling%rk
      if (self%i_src>0) &
      self%table(:,:,self%i_src) = self%table(:,:,self%i_src)/scaling%src
    end if
    !$omp end critical (io_cr)
  end if
END SUBROUTINE init_table

!> =============================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =============================================================================
SUBROUTINE lookup (self, dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma
  !.............................................................................
  real, dimension(:,:,:), pointer               :: x_loc, y_loc
  real, dimension(dim(1))                       :: p1, p2, q1, q2
  real, dimension(dim(1),4,size(self%table,3))  :: f
  real, dimension(dim(1),4)                     :: g
  integer, dimension(dim(1))                    :: i1, i2
  integer                                       :: ix, iy, iz, iv, m1, m2, m3
  real                                          :: clnd, css
  logical                                       :: skip
  logical, save                                 :: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays, if necessary to take the log of an input
  !-----------------------------------------------------------------------------
  m3=  self%n3
  m1 = self%n1
  m2 = self%n2
  if (present(pg).or.present(tt).or.present(rk).or.present(src)) then
    skip = .false.
  else
    skip = .true.
  end if
  if (present(lnx)) then
    x_loc => lnx
  else if (present(x)) then
    allocate (x_loc(dim(1),dim(2),dim(3)))
    x_loc = log(x)
  else
    call mpi%abort ('eos_t: missing 1st argument')
  end if
  if (present(lny)) then
    y_loc => lny
  else if (present(y)) then
    allocate (y_loc(dim(1),dim(2),dim(3)))
    y_loc = log(y)
  else
    call mpi%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  ! Loop over points
  !-----------------------------------------------------------------------------
  cx = 1./self%x_step
  cy = 1./self%y_step
  do iz=1,dim(3)
  do iy=1,dim(2)
    !-------------------------------------------------------------------------
    ! Addresses and weights, vectorized
    !-------------------------------------------------------------------------
    p1 = 1.+(x_loc(:,iy,iz)-self%x_min)*cx
    p2 = 1.+(y_loc(:,iy,iz)-self%y_min)*cy
    i1 = p1
    i2 = p2
    i1 = max(1,min(i1,m1))
    i2 = max(1,min(i2,m2))
    p1 = p1-i1
    p2 = p2-i2
    q1 = 1.0-p1
    q2 = 1.0-p2
    g(:,1) = q1*q2
    g(:,2) = p1*q2
    g(:,3) = q1*p2
    g(:,4) = p1*p2
    !---------------------------------------------------------------------------
    !  Table values, indirect addressing
    !---------------------------------------------------------------------------
    do ix=1,dim(1)
    do iv=1,m3
      f(ix,1,iv) = self%table(i1(ix)  ,i2(ix)  ,iv)
      f(ix,2,iv) = self%table(i1(ix)+1,i2(ix)  ,iv)
      f(ix,3,iv) = self%table(i1(ix)  ,i2(ix)+1,iv)
      f(ix,4,iv) = self%table(i1(ix)+1,i2(ix)+1,iv)
    end do
    end do
    if (skip) cycle
    !---------------------------------------------------------------------------
    !  Loop over table entries
    !---------------------------------------------------------------------------
    if (present(pg)) then
      pg(:,iy,iz) = &
         f(:,1,self%i_pg)*g(:,1) &
       + f(:,2,self%i_pg)*g(:,2) &
       + f(:,3,self%i_pg)*g(:,3) &
       + f(:,4,self%i_pg)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(tt)) then
      tt(:,iy,iz) = &
         f(:,1,self%i_tt)*g(:,1) &
       + f(:,2,self%i_tt)*g(:,2) &
       + f(:,3,self%i_tt)*g(:,3) &
       + f(:,4,self%i_tt)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(rk)) then
      do iv=0,size(rk,4)-1
        rk(:,iy,iz, iv+1) = &
           f(:,1,self%i_rk+iv)*g(:,1) &
         + f(:,2,self%i_rk+iv)*g(:,2) &
         + f(:,3,self%i_rk+iv)*g(:,3) &
         + f(:,4,self%i_rk+iv)*g(:,4)
      end do
    end if
    !---------------------------------------------------------------------------
    if (present(src)) then
     do iv=0,size(src,4)-1
      src(:,iy,iz,iv+1) = &
         f(:,1,self%i_src+iv)*g(:,1) &
       + f(:,2,self%i_src+iv)*g(:,2) &
       + f(:,3,self%i_src+iv)*g(:,3) &
       + f(:,4,self%i_src+iv)*g(:,4)
     end do
    end if
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Test printout
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (io_cr)
    if (first_time) then
      first_time = .false.
      print *, 'lookup: x_loc =', x_loc(1,1,1)
      print *, 'lookup: y_loc =', y_loc(1,1,1)
      if (present(pg)) &
      print *, 'lookup: pg    =', pg   (1,1,1)
    end if
    !$omp end critical (io_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate if allocated
  !-----------------------------------------------------------------------------
  if (present(x))  deallocate (x_loc)
  if (present(y))  deallocate (y_loc)
END SUBROUTINE lookup

END MODULE eos_mod
