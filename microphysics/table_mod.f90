!===============================================================================
!> Module using a table prepared by a Python procedure
!===============================================================================
MODULE table_mod
  USE trace_mod
  USE io_mod
  USE io_unit_mod
  USE scaling_mod
  private
  type, public:: table_t
    integer:: ioformat
    integer:: i_pg, i_tt, i_ss, i_src, i_ct, i_cs, i_gamma, i_cv, i_kr, i_kp
    integer:: n1, n2, n3
    logical:: first_time=.true.
    logical:: first_print=.true.
    real:: x_min, x_step
    real:: y_min, y_step
    real, dimension(:,:,:), pointer:: data
  contains
    procedure:: init
    procedure:: lookup
  end type
CONTAINS

!===============================================================================
!> Initialize equation-of-state handling
!===============================================================================
SUBROUTINE init (self, file)
  class(table_t):: self
  character(len=*):: file
  !-----------------------------------------------------------------------------
  call trace%begin ('table_t%init')
  if (self%first_time) then
    call io%lock%set ('table_t%init')
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
      read (io_unit%data) self%i_pg, self%i_tt, self%i_gamma, self%i_ss, self%i_kr, self%i_kp, self%i_src
      if (io%master) then
        print '(a,4i4)', ' init_table: ioformat, n1, n2, nvar =', &
        self%ioformat, self%n1, self%n2, self%n3
        print *,' init_table: dmin, dstep =', self%x_min, self%x_step
        print *,' init_table: emin, estep =', self%y_min, self%y_step
        print '(a,9i4)', ' items: i_pg, i_tt, i_gamma, i_ss, i_kr, i_kp, i_src =', &
        self%i_pg, self%i_tt, self%i_gamma,self%i_ss, self%i_kr, self%i_kp, self%i_src
      end if
      allocate (self%data(self%n1,self%n2,self%n3))
      read (io_unit%data) self%data
      close (io_unit%data)
      !-------------------------------------------------------------------------
      ! Rescale to the current code units
      !-------------------------------------------------------------------------
      self%data(:,:,self%i_pg ) = self%data(:,:,self%i_pg )/scaling%p
      self%data(:,:,self%i_tt ) = self%data(:,:,self%i_tt )/scaling%temp
      if (self%i_kr>0) &
      self%data(:,:,self%i_kr ) = self%data(:,:,self%i_kr )*scaling%l
      if (self%i_kp>0) &
      self%data(:,:,self%i_kp ) = self%data(:,:,self%i_kp )*scaling%l
      if (self%i_src>0) &
      self%data(:,:,self%i_src) = self%data(:,:,self%i_src)/(scaling%e/scaling%l**2/scaling%t)
    end if
    call io%lock%unset ('table_t%init')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!>
!> Optimized to vectorize where possible.  In gfortran taking a log takes 
!> about 8 ns/pt, and the address computation takes about 10 ns.  Each 
!> variable to compute after that takes about 6 ns, and thus the average
!> for 4 opacities and 4 source functions is about 10/8 + 6 = 7 ns/pt
!===============================================================================
SUBROUTINE lookup (self, dim, lnd, lny, pg, tt, ss, rk, src, gamma)
  class (table_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, lny
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma(:,:,:)
  !.............................................................................
  real                                          :: p1, p2, q1, q2
  real, dimension(dim(1),4,size(self%data,3))  :: f
  real, dimension(dim(1),4)                     :: g
  integer, dimension(dim(1))                    :: i1, i2
  integer                                       :: ix, iy, iz, iv, m0, m1, m2, m3
  real                                          :: cx, cy
  logical, save                                 :: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Figure out first and last table variable index needed
  !-----------------------------------------------------------------------------
  m1 = self%n1
  m2 = self%n2
  m3 = 0
  if (present(pg))     m3 = max(m3,self%i_pg)
  if (present(tt))     m3 = max(m3,self%i_tt)
  if (present(gamma))  m3 = max(m3,self%i_gamma)
  if (present(ss))     m3 = max(m3,self%i_ss)
  if (present(rk))     m3 = max(m3,self%i_kr +size(rk ,4)-1)
  if (present(src))    m3 = max(m3,self%i_src+size(src,4)-1)
  m0 = self%n3
  if (present(pg))     m0 = min(m0,self%i_pg)
  if (present(tt))     m0 = min(m0,self%i_tt)
  if (present(gamma))  m0 = min(m0,self%i_gamma)
  if (present(ss))     m0 = min(m0,self%i_ss)
  if (present(rk))     m0 = min(m0,self%i_kr +size(rk ,4)-1)
  if (present(src))    m0 = min(m0,self%i_src+size(src,4)-1)
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays, if necessary to take the log of an input
  !-----------------------------------------------------------------------------
  if (self%first_print) then
    print *,'m0,m3',m0,m3
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
    do ix=1,dim(1)
      p1 = 1.+(lnd(ix,iy,iz)-self%x_min)*cx
      p2 = 1.+(lny(ix,iy,iz)-self%y_min)*cy
      i1(ix) = p1
      i2(ix) = p2
      if (i1(ix)>m1) print *, 'i1,m1', i1(ix), m1
      if (i2(ix)>m2) print *, 'i2,m2', i2(ix), m2
      i1(ix) = max(1,min(i1(ix),m1))
      i2(ix) = max(1,min(i2(ix),m2))
      p1 = p1-i1(ix)
      p2 = p2-i2(ix)
      q1 = 1.0-p1
      q2 = 1.0-p2
      g(ix,1) = q1*q2
      g(ix,2) = p1*q2
      g(ix,3) = q1*p2
      g(ix,4) = p1*p2
    end do
    !---------------------------------------------------------------------------
    ! Table values via indirect addressing, possibly vectorizes (?) over iv, or
    ! at least does not require indirect addressing for every value.
    !---------------------------------------------------------------------------
    do ix=1,dim(1)
    do iv=m0,m3
      f(ix,1,iv) = self%data(i1(ix)  ,i2(ix)  ,iv)
      f(ix,2,iv) = self%data(i1(ix)+1,i2(ix)  ,iv)
      f(ix,3,iv) = self%data(i1(ix)  ,i2(ix)+1,iv)
      f(ix,4,iv) = self%data(i1(ix)+1,i2(ix)+1,iv)
    end do
    end do
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
    if (present(ss)) then
      ss(:,iy,iz) = &
         f(:,1,self%i_ss)*g(:,1) &
       + f(:,2,self%i_ss)*g(:,2) &
       + f(:,3,self%i_ss)*g(:,3) &
       + f(:,4,self%i_ss)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(rk)) then
      rk(:,iy,iz, 1) = &
         f(:,1,self%i_kr)*g(:,1) &
       + f(:,2,self%i_kr)*g(:,2) &
       + f(:,3,self%i_kr)*g(:,3) &
       + f(:,4,self%i_kr)*g(:,4)
       rk(:,iy,iz, 2) = &
         f(:,1,self%i_kp)*g(:,1) &
       + f(:,2,self%i_kp)*g(:,2) &
       + f(:,3,self%i_kp)*g(:,3) &
       + f(:,4,self%i_kp)*g(:,4)
    end if
    ! --------------------------------------------------------------------------
    if (present(gamma)) then
      gamma(:,iy,iz) = &
         f(:,1,self%i_gamma)*g(:,1) &
       + f(:,2,self%i_gamma)*g(:,2) &
       + f(:,3,self%i_gamma)*g(:,3) &
       + f(:,4,self%i_gamma)*g(:,4)
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
  if (self%first_print) then
    call io%lock%set ('table_t%test')
    if (self%first_print) then
      self%first_print = .false.
      print *, 'lookup: x_loc =', lnd(1,1,1)
      print *, 'lookup: y_loc =', lny(1,1,1)
      if (present(pg)) &
      print *, 'lookup: pg    =', pg(1,1,1)
      if (present(ss)) &
      print *, 'lookup: ss    =', ss(1,1,1)
    end if
    call io%lock%unset ('table_t%test')
  end if

END SUBROUTINE lookup
END MODULE table_mod
