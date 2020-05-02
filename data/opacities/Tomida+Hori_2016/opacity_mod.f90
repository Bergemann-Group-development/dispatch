!===============================================================================
!> Read the kappa data from Tomida & Hori 2016, which comes in two tables, with
!> Rosseland and Planck means
!===============================================================================
MODULE opacity_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE eos_mod
  USE lagrange_mod
  USE scaling_mod
  USE timer_mod
  implicit none
  private
  type, public:: opacity_t
    integer:: nT=450, nd=450, nv=2, nwb =1
    real, dimension(:), pointer:: log_T => null()
    real, dimension(:), pointer:: log_d => null()
  contains
    procedure:: init
  end type
  real, dimension(:,:,:), pointer:: table1, table2, table
  type (opacity_t), public:: opacity
CONTAINS

!===============================================================================
!> Read in opacity data from a table file, the size of which we might as well
!> hardwire, since it is specific
!===============================================================================
SUBROUTINE init (self)
  class (opacity_t):: self
  real:: row(3), log_d, log_T
  integer:: iT, id, iv, i, j, k
  real, dimension(:), pointer:: x, y
  real, dimension(:,:,:), pointer:: table1, table2
  !-----------------------------------------------------------------------------
  call trace%begin ('opacity_t%init (Tomida+Hori_2016)')
  call timer%tic
  if (io_unit%master) print *, &
 '---------------- reading opacity data by Tomida & Hori 2016 ------------------'
  allocate (self%log_T(self%nT), self%log_d(self%nd), table1(self%nT,self%nd,self%nv))
  open (io_unit%tmp+1, form='formatted', status='old', &
        file=trim(io_unit%top)//'/data/opacities/Tomida+Hori_2016/kr.dat')
  open (io_unit%tmp+2, form='formatted', status='old', &
        file=trim(io_unit%top)//'/data/opacities/Tomida+Hori_2016/kp.dat')
  do k=1,self%nv
    do i=1,self%nT
      do j=1,self%nd
        read (io_unit%tmp+k,*) row
        table1(i,j,k) = row(3)
        self%log_d(j) = row(1)
      end do
      self%log_T(i) = row(2)
    end do
    if (io%verbose > 0 .and. io%master) then
      print 1, k, 'log_kappa_min, log_kappa_max [cgs] =', minval(table1(:,:,k)), &
                                                          maxval(table1(:,:,k))
      1 format(i3,2x,a,1p,2e12.4)
    end if
    close (io_unit%tmp+k)
  end do
  if (io%verbose > 0 .and. io%master) then
    print 2, 'T_min, T_max [K]    =', 10.**self%log_T(1), 10.**self%log_T(self%nT)
    print 2, 'd_min, d_max [cgs0] =', 10.**self%log_d(1), 10.**self%log_d(self%nd)
  end if
  2 format(1x,a,1p,2e15.6)
!===============================================================================
!> We now have data as a function of (log_T,log_d), and we need to be interpolate
!> and append to the eos%table
!===============================================================================
  if (io%master) print *, &
  '------------------------- append to fast EOS table --------------------------'
  !-----------------------------------------------------------------------------
  ! For each temperature, interpolate in density to the EOS table values
  !-----------------------------------------------------------------------------
  allocate (table2(self%nT,eos%nd,self%nv))
  x => self%log_d
  do k=1,self%nv
    do i=1,self%nT
      y => table1(i,:,k)
      do j=1,eos%nd
        log_d = eos%table%y_min + (j-1)*eos%table%y_step + log10(scaling%d)
        table2(i,j,k) = interpolate(x, y, log_d)
        if (io%verbose > 2 .and. io%master) then
          if ((i==1.or.i==252.or.i==253.or.i==self%nT) .and. &
              (j==1.or.j==eos%nd)) &
            print *, i, j, k, self%log_T(i), &
              eos%lnd(j)/log(10.)+log10(scaling%d), table2(i,j,k)
        end if
      end do
    end do
    if (io%verbose > 0 .and. io%master) then
      print 3, k, 'log_kappa_min, log_kappa_max [cgs] =', &
        minval(table2(:,:,k)), minloc(table2(:,:,k)), &
        maxval(table2(:,:,k)), maxloc(table2(:,:,k))
      3 format(i3,2x,a,1p,2(e12.4,2i4))
    end if
  end do
  !-----------------------------------------------------------------------------
  ! For each density, interpolate in temperature to the EOS table values for T,
  ! which are eos%table%data(:,:,eos%table%i_tt)
  !-----------------------------------------------------------------------------
  allocate (table(eos%ne,eos%nd,eos%nv+self%nv))
  table(:,:,1:eos%nv) = eos%table%data
  x => self%log_T                       ! opacity table temperatures
  do k=1,self%nv
    do j=1,eos%nd
      y => table2(:,j,k)                   ! interpolated table values
      log_d = eos%table%y_min + (j-1)*eos%table%y_step + log10(scaling%d)
      do i=1,eos%ne
        log_T =  log10(eos%table%data(i,j,eos%table%i_tt))
        table(i,j,eos%nv+k) = 10.**interpolate(x, y, log_T)
        if (io%verbose > 2 .and. io%master) then
          if ((i==1.or.i==eos%ne).and.(j==1.or.j==eos%nd)) &
            print *, i, j, k, log_T, log_d, table(i,j,eos%nv+k)
        end if
      end do
    end do
    !---------------------------------------------------------------------------
    ! Rescaling table to code units, area per unit mass
    !---------------------------------------------------------------------------
    if (io%verbose > 0 .and. io%master) then
      print 3, k, 'kappa_min, kappa_max [cgs]  =', &
        minval(table(:,:,k)), minloc(table(:,:,k)), &
        maxval(table(:,:,k)), maxloc(table(:,:,k))
      table(:,:,eos%nv+k) = table(:,:,eos%nv+k)/(scaling%l**2/scaling%m)
      print 3, k, 'kappa_min, kappa_max [code] =', &
        minval(table(:,:,k)), minloc(table(:,:,k)), &
        maxval(table(:,:,k)), maxloc(table(:,:,k))
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Deallocate the reader table and the old EOS table, insert the new EOS table
  ! and register the new table size
  !-----------------------------------------------------------------------------
  deallocate (eos%table%data, table1, table2)
  eos%table%data => table
  eos%table%i_rk = eos%nv+1
  eos%table%i_src = eos%nv+2
  eos%nv = eos%nv+self%nv
  call timer%toc ('opacity_t%init')
  call trace%end()
END SUBROUTINE init

!-------------------------------------------------------------------------------
! Interpolation using Lagrange interpolation in short sequences
!-------------------------------------------------------------------------------
FUNCTION interpolate (x, y, x0) RESULT (out)
  !USE lagrange_mod
  real, dimension(:):: x, y
  real:: x0, out, p
  integer:: i
  !-----------------------------------------------------------------------------
  out = lagrange%sequence (x0, x, y)
END FUNCTION interpolate

END MODULE opacity_mod
