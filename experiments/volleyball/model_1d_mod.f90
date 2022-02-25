!===============================================================================
!> Generic 1-D table handling
!===============================================================================
MODULE model_1d_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public:: model_1d_t
    integer              :: verbose=0
    integer              :: n, m
    integer, allocatable :: idx(:)
    real(8), allocatable :: r(:), d(:), e(:), g(:), tg(:)
  contains
    procedure:: init
    procedure:: read
    procedure:: interpolate
  end type
  real(8):: margin=1d0
  character(64):: top='../../data/', dir='solar_models/', file='JCD-model-S'
  type(model_1d_t), public:: model_1d
CONTAINS

!===============================================================================
!> Initialize the table, prepare for interpolation
!===============================================================================
SUBROUTINE init (self)
  class(model_1d_t):: self
  integer:: iostat
  namelist /model_1d_params/ top, dir, file, margin
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------¨
  call trace%begin ('model_1d_t%init')
  call io%lock%set ('model_1d_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, model_1d_params, iostat=iostat)
    write (io%namelist, model_1d_params)
    call model_1d%read (trim(top)//trim(dir)//trim(file))
  end if
  call io%lock%unset ('model_1d_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Read a table of (n,m) values, from n lines with m values.  The 2nd line in
!> the file gives the variable indices that the table covers, eg corresponding
!> to density and internal energy
!===============================================================================
SUBROUTINE read (self, file)
  class(model_1d_t):: self
  character(*):: file
  integer:: unit, i
  !-----------------------------------------------------------------------------
  call trace%begin ('model_1d_t%read')
  open (newunit=unit, file=file, form='unformatted', status='old')
  read (unit) self%n, self%m
  allocate (self%r(self%n), self%g(self%n))
  allocate (self%d(self%n), self%e(self%n), self%tg(self%n))
  do i=1,self%n
    read (unit) self%r(i), self%g(i), self%d(i), self%e(i), self%tg(i)
  end do
  ! ---------------------------------
  ! Convert specific internal energy to internal energy
  ! ---------------------------------
  self%e=self%e*self%d
  ! ---------------------------------
  ! Rescale to code units
  ! ---------------------------------
  call info ('model_1d%read: before scaling (CGS)')
  self%r = self%r*cgs%r_sun/scaling%l
  self%g = self%g*scaling%t**2/scaling%l
  self%d = self%d/scaling%d
  self%e = self%e/scaling%e
  call info ('model_1d%read: after scaling (CODE)')
  close (unit)
  call trace%end ()
contains
subroutine info (label)
  character(*):: label
  write (stdout,*) '---------------------- '//label//' ------------------------'
  write (stdout,1) &
    'model_1d:         min/max mass density:', minval(self%d), maxval(self%d)
  write (stdout,1) &
    'model_1d: min/max energy per unit mass:', minval(self%e), maxval(self%e)
  write (stdout,1) &
    'model_1d: min/max radius:', minval(self%r), maxval(self%r)
  1 format(a,1p,2e14.3)
end subroutine
END SUBROUTINE read

!===============================================================================
!> Interpolate in the model; this is only done at upstart, so no need for speed
!===============================================================================
SUBROUTINE interpolate (self, r, d, e)
  class(model_1d_t):: self
  real(8):: r, p, d, e
  integer:: i
  !-----------------------------------------------------------------------------
  ! Find the table index k that corresponds to i
  !-----------------------------------------------------------------------------
  call trace%begin ('model_1d_t%interpolate')
  do i=2,self%n-2
    if (self%r(i) < r) then
      exit
    end if
  end do
  p = (r-self%r(i-1))/(self%r(i)-self%r(i-1))
  p = max(min(p,1d0+margin),-margin)
  d = self%d(i-1)*(1d0-p) + self%d(i)*p
  e = self%e(i-1)*(1d0-p) + self%e(i)*p
  call trace%end ()
END SUBROUTINE interpolate

END MODULE model_1d_mod
