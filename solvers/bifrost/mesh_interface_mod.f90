!===============================================================================
!> Interface module, with a data type that can store all required BIFROST mesh
!> details as part of the solver_t data type, for get/put operations.
!===============================================================================
MODULE mesh_interface_mod
  ! Dispatch modules
  USE io_unit_mod
  USE io_mod
  USE mpi_mod
  USE mesh_mod
  USE trace_mod
  ! Bifrost modules
  USE params_mod
  USE parallel_mod
  USE stagger_mod
  implicit none
  private
  type, public:: bifrost_mesh_t
    class(mesh_t), dimension(:), pointer:: mesh
    real,dimension(:),   allocatable :: xm,xmdn,dxidxup,dxidxdn
    real,dimension(:),   allocatable :: dx1d
    real,dimension(:),   allocatable :: ym,ymdn,dyidyup,dyidydn
    real,dimension(:),   allocatable :: dy1d
    real,dimension(:),   allocatable :: zm,zmdn,dzidzup,dzidzdn
    real,dimension(:),   allocatable :: dz1d
    real,dimension(:,:), allocatable :: zlincoeff,dzlincoeff
    real,dimension(:,:), allocatable :: zuincoeff,dzuincoeff
    real :: dx,dy,dz
    integer(4):: xs,xe,xsb,xeb
    integer(4):: ys,ye,ysb,yeb
    integer(4):: zs,ze,zsb,zeb
    character(len=132) :: meshfile
    integer(4)::myid
  contains
    procedure:: check
    procedure:: get
    procedure:: init
    procedure:: init_bifrost_vars
    procedure:: one_sided_coeff                                         ! FIXME: should be moved to stagger_t in phase2 of dispatch-bifrost interface
    procedure:: put
    procedure:: read_bifrost_mesh
    procedure:: setup_mesh
    procedure:: stagger_inv                                             ! FIXME: should be moved to stagger_t in phase2 of dispatch-bifrost interface
  end type
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!===============================================================================
!> Allocate patch space for Bifrost mesh data
!===============================================================================
SUBROUTINE init (self, mesh)
  class(bifrost_mesh_t):: self
  class(mesh_t),dimension(:), pointer :: mesh
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_mesh_t%init')
  self%mesh => mesh
  call self%read_bifrost_mesh
  call self%init_bifrost_vars
  call trace%end()
END SUBROUTINE init

!=======================================================================
!> Allocate global variables, used by Bifrost
!=======================================================================
SUBROUTINE init_bifrost_vars (self)
  class(bifrost_mesh_t):: self
  !---------------------------------------------------------------------
  class(mesh_t), pointer :: m
  logical, save :: first_time = .true.
  !.....................................................................
  call trace%begin ('bifrost_mesh_t%init_bifrost_vars')
  !first_time=.false.
  if (first_time) then
    call io%lock%set
    first_time = .false.
    m => self%mesh(1)   ! X direction
    allocate(xm     (m%b_sb:m%b_eb))
    allocate(xmdn   (m%b_sb:m%b_eb))
    allocate(dxidxup(m%b_sb:m%b_eb))
    allocate(dxidxdn(m%b_sb:m%b_eb))
    allocate(dx1d   (m%b_sb:m%b_eb))
    m=> self%mesh(2)    ! Y direction
    allocate(ym     (m%b_sb:m%b_eb))
    allocate(ymdn   (m%b_sb:m%b_eb))
    allocate(dyidyup(m%b_sb:m%b_eb))
    allocate(dyidydn(m%b_sb:m%b_eb))
    allocate(dy1d   (m%b_sb:m%b_eb))
    m => self%mesh(3)   ! Z direction
    allocate(zm     (m%b_sb:m%b_eb))
    allocate(zmdn   (m%b_sb:m%b_eb))
    allocate(dzidzup(m%b_sb:m%b_eb))
    allocate(dzidzdn(m%b_sb:m%b_eb))
    allocate(dz1d   (m%b_sb:m%b_eb))
    allocate(zlincoeff (m%b_sb:m%b_eb, 0:5))
    allocate(dzlincoeff(m%b_sb:m%b_eb, 0:5))
    allocate(zuincoeff (m%b_sb:m%b_eb, 0:5))
    allocate(dzuincoeff(m%b_sb:m%b_eb, 0:5))
    c = 3./256.
    b = -25./256.
    a = .5-b-c
    a1st = 0.5
    call io%lock%unset
  end if
  call trace%end()

END SUBROUTINE init_bifrost_vars

!===============================================================================
!> Put the patch values back into the Bifrost variables
!===============================================================================
SUBROUTINE put (self)
  class(bifrost_mesh_t):: self
  !---------------------------------------------------------------------
  class(mesh_t), pointer :: m
  call trace%begin ('bifrost_mesh_t%put')
  meshfile   = self%meshfile
  m => self%mesh(1)   ! X direction
  xs         = m%b_s
  xsb        = m%b_sb
  xe         = m%b_e
  xeb        = m%b_eb
  adx        = m%b_ad
  bdx        = m%b_bd
  cdx        = m%b_cd
  dx         = m%b_d
  xm         = m%b_m
  xmdn       = m%b_mdn
  dxidxup    = m%b_didup
  dxidxdn    = m%b_diddn
  dx1d       = m%b_1d
  m=> self%mesh(2)    ! Y direction
  ys         = m%b_s
  ysb        = m%b_sb
  ye         = m%b_e
  yeb        = m%b_eb
  ady        = m%b_ad
  bdy        = m%b_bd
  cdy        = m%b_cd
  dy         = m%b_d
  ym         = m%b_m
  ymdn       = m%b_mdn
  dyidyup    = m%b_didup
  dyidydn    = m%b_diddn
  dy1d       = m%b_1d
  m => self%mesh(3)
  zs         = m%b_s
  zsb        = m%b_sb
  ze         = m%b_e
  zeb        = m%b_eb
  adz        = m%b_ad
  bdz        = m%b_bd
  cdz        = m%b_cd
  dz         = m%b_d
  zm         = m%b_m
  zmdn       = m%b_mdn
  dzidzup    = m%b_didup
  dzidzdn    = m%b_diddn
  dz1d       = m%b_1d
  zlincoeff  = m%zlincoeff
  dzlincoeff = m%dzlincoeff
  zuincoeff  = m%zuincoeff
  dzuincoeff = m%dzuincoeff
  call trace%end()

END SUBROUTINE put

SUBROUTINE check (self)
  class(bifrost_mesh_t):: self
  !---------------------------------------------------------------------
  class(mesh_t), pointer :: m
  integer:: i
  call trace%begin ('bifrost_mesh_t%put')
  meshfile   = self%meshfile
  m => self%mesh(1)   ! X direction
  if (xs/=m%b_s) print *, 'ERROR: xs neq m%b_s' 
  if (xsb/=m%b_sb) print *, 'ERROR: xsb neq m%b_sb' 
  if (xe/=m%b_e) print *, 'ERROR: xe neq m%b_e'
  if (xeb/=m%b_eb) print *, 'ERROR: xeb neq m%b_eb'  
  if (adx/=m%b_ad) print *, 'ERROR: adx neq m%b_ad'
  if (bdx/=m%b_bd) print *, 'ERROR: bdx neq m%b_bd'
  if (cdx/=m%b_cd) print *, 'ERROR: cdx neq m%b_cd'
  if (dx/=m%b_d) print *, 'ERROR: dx neq m%b_d', dx, m%b_d 
  if (any(xm/=m%b_m))  print *, 'ERROR: xm neq m%b_m'
  if (any(xmdn/=m%b_mdn)) print *, 'ERROR: xmdn neq m%b_mdn'
  if (any(dxidxup/=m%b_didup)) print *, 'ERROR: dxidxup/=m%b_didup'
  if (any(dxidxdn/=m%b_diddn)) then 
  print *, 'ERROR: dxidxdn/=m%b_diddn'
      do i =xsb,xeb
      print *, i, dxidxdn(i),m%b_diddn(i)
    end do
    error stop "end"
  end if
  if (any(dx1d/=m%b_1d)) print *, 'ERROR: dx1d/=m%b_1d'
  
  m=> self%mesh(2)    ! Y direction
  if (ys/=m%b_s) print *, 'ERROR: ys neq m%b_s' 
  if (ysb/=m%b_sb) print *, 'ERROR: ysb neq m%b_sb' 
  if (ye/=m%b_e) print *, 'ERROR: ye neq m%b_e'
  if (yeb/=m%b_eb) print *, 'ERROR: yeb neq m%b_eb'  
  if (ady/=m%b_ad) print *, 'ERROR: ady neq m%b_ad'
  if (bdy/=m%b_bd) print *, 'ERROR: bdy neq m%b_bd'
  if (cdy/=m%b_cd) print *, 'ERROR: cdy neq m%b_cd'
  if (dy/=m%b_d) print *, 'ERROR: dy neq m%b_d', dy, m%b_d
  if (any(ym/=m%b_m)) print *, 'ERROR: ym neq m%b_m'
  if (any(ymdn/=m%b_mdn)) print *, 'ERROR: ymdn neq m%b_mdn'
  if (any(dyidyup/=m%b_didup)) print *, 'ERROR: dyidyup/=m%b_didup'
  if (any(dyidydn/=m%b_diddn)) print *, 'ERROR: dyidydn/=m%b_diddn'
  if (any(dy1d/=m%b_1d)) print *, 'ERROR: dy1d/=m%b_1d'
  m => self%mesh(3)
  if (zs/=m%b_s) print *, 'ERROR: zs neq m%b_s' 
  if (zsb/=m%b_sb) print *, 'ERROR: zsb neq m%b_sb' 
  if (ze/=m%b_e) print *, 'ERROR: ze neq m%b_e'
  if (zeb/=m%b_eb) print *, 'ERROR: zeb neq m%b_eb'  
  if (adz/=m%b_ad) print *, 'ERROR: adz neq m%b_ad'
  if (bdZ/=m%b_bd) print *, 'ERROR: bdz neq m%b_bd'
  if (cdZ/=m%b_cd) print *, 'ERROR: cdz neq m%b_cd'
  if (dz/=m%b_d) print *, 'ERROR: dz neq m%b_d', dz, m%b_d
  if (any(Zm/=m%b_m)) print *, 'ERROR: zm neq m%b_m'
  if (any(Zmdn/=m%b_mdn)) print *, 'ERROR: zmdn neq m%b_mdn'
  if (any(dzidzup/=m%b_didup)) print *, 'ERROR: dzidzup/=m%b_didup'
  if (any(dzidzdn/=m%b_diddn)) then
   print *, 'ERROR: dzidzdn/=m%b_diddn'
  do i =zsb,zeb
      print *, i, dzidzdn(i),m%b_diddn(i)
    end do
    error stop "end"
  end if
  if (any(dz1d/=m%b_1d)) print *, 'ERROR: dz1d/=m%b_1d'
  if (any(zlincoeff/=m%zlincoeff)) print *, 'zlincoeff/=m%zlincoeff'
  if (any(dzlincoeff/=m%dzlincoeff)) print *, 'dzlincoeff/=m%dzlincoeff'
  if (any(zuincoeff/=m%zuincoeff)) print *, 'zuincoeff/=m%zuincoeff'
  if (any(dzuincoeff/=m%dzuincoeff)) print *, 'dzuincoeff/=m%dzuincoeff'

  call trace%end()

END SUBROUTINE check

!===============================================================================
!> Copy the Bifrost mesh data to patch memory
!===============================================================================
SUBROUTINE get (self)
  class(bifrost_mesh_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_mesh_t%get')
  error stop "error: bifrost_mesh_t%get should not be called."
  self%meshfile   = meshfile
  self%dx         = dx
  self%xm         = xm
  self%xmdn       = xmdn
  self%dxidxup    = dxidxup
  self%dxidxdn    = dxidxdn
  self%dx1d       = dx1d
  self%dy         = dy
  self%ym         = ym
  self%ymdn       = ymdn
  self%dyidyup    = dyidyup
  self%dyidydn    = dyidydn
  self%dy1d       = dy1d
  self%dz         = dz
  self%zm         = zm
  self%zmdn       = zmdn
  self%dzidzup    = dzidzup
  self%dzidzdn    = dzidzdn
  self%dz1d       = dz1d
  self%zlincoeff  = zlincoeff
  self%dzlincoeff = dzlincoeff
  self%zuincoeff  = zuincoeff
  self%dzuincoeff = dzuincoeff
  self%xs         = xs
  self%xe         = xe
  self%xsb        = xsb
  self%xeb        = xeb
  self%ys         = ys
  self%ye         = ye
  self%ysb        = ysb
  self%yeb        = yeb
  self%zs         = zs
  self%ze         = ze
  self%zsb        = zsb
  self%zeb        = zeb
  call trace%end()
END SUBROUTINE get

!=======================================================================
!> Read external mesh file
!=======================================================================
SUBROUTINE read_bifrost_mesh (self)
  class(bifrost_mesh_t)     :: self
  ! --------------------------------------------------------------------
  class(mesh_t), pointer :: m
  logical, save:: first_time=.true.
  character(len=132), save :: meshfile='qt_128x128x128.mesh'
  real, save:: dx=1.875e-1
  real, save:: dy=1.875e-1
  real, save:: dz=1.969e-1
  integer, save:: dims(3)
  integer, save:: mx = 128
  integer, save:: my = 128
  integer, save:: mz = 128
  real, save, dimension(:), pointer:: xm, xmdn, dxidxup, dxidxdn, dx1d
  real, save, dimension(:), pointer:: ym, ymdn, dyidyup, dyidydn, dy1d
  real, save, dimension(:), pointer:: zm, zmdn, dzidzup, dzidzdn, dz1d
  integer:: i(3), iostat, im
  real:: dd(3)
  namelist /bifrost_mesh_params/ meshfile, dims, dx, dy, dz
  ! ....................................................................
  call trace%begin ('bifrost_t%read_mesh')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, bifrost_mesh_params, iostat=iostat)
    if (io%master) write (*, bifrost_mesh_params)
    call read_params("stagger")                                         !FIXME: should be a local parameter set, already present in the namelist.
  end if
  !---------------------------------------------------------------------------
  ! Read the mesh from the .mesh file
  !---------------------------------------------------------------------------
  if (io%master) print *,'reading mesh data from ',trim(meshfile),' ', dims
  open (io_unit%data, file=trim(meshfile),status='old', form='formatted')
  ! read dimensions
  ! --- X dimension ---
  read(io_unit%data, *, iostat=iostat) mx
  allocate(xm(mx))
  allocate(xmdn(mx))
  allocate(dxidxup(mx))
  allocate(dxidxdn(mx))
  allocate(dx1d(mx))
  read(io_unit%data, *) xm, xmdn, dxidxup, dxidxdn
  dx1d(:)=abs(cshift(xm(:),1)-xm(:))
  if (mx > 1) then
    dx1d(mx)=abs(xm(mx)-xm(mx-1))
  else
    dx1d(mx)=1.
  endif
  ! --- Y dimension ---
  read(io_unit%data, *, iostat=iostat) my
  allocate(ym(my))
  allocate(ymdn(my))
  allocate(dyidyup(my))
  allocate(dyidydn(my))
  allocate(dy1d(my))
  read(io_unit%data, *, iostat=iostat) ym, ymdn, dyidyup, dyidydn
  dy1d(:)=abs(cshift(ym(:),1)-ym(:))
  if (my > 1) then
    dy1d(my)=abs(ym(my)-ym(my-1))
  else
    dx1d(my)=1.
  endif
  ! --- Z dimension
  read(io_unit%data, *, iostat=iostat) mz
  allocate(zm(mz))
  allocate(zmdn(mz))
  allocate(dzidzup(mz))
  allocate(dzidzdn(mz))
  allocate(dz1d(mz))
  read(io_unit%data, *, iostat=iostat) zm, zmdn, dzidzup, dzidzdn
  dz1d(:)=abs(cshift(zm(:),1)-zm(:))
  if (mz > 1) then
    dz1d(mz)=abs(zm(mz)-zm(mz-1))
  else
    dz1d(mz)=1.
  endif
  !
  if (.not. all(dims==[mx,my,mz])) then
    if (io%master) print *, 'ERROR: file dims =', dims
    if (io%master) print *, 'ERROR: mesh dims =', [mx,my,mz]
    call mpi%abort
  end if
  close (io_unit%data)
  call io%lock%unset
  ! --------------------------------------------------------------------
  ! Set up index space
  ! --------------------------------------------------------------------
  self%meshfile = trim(meshfile)
  dd=[dx,dy,dz]
  do im=1,3
    m => self%mesh(im)
    i = 1 + nint(m%p/m%d) - (m%n+1)/2             ! index offsets
    m%b_s =i(im)
    m%b_d = dd(im)
    m%b_e =i(im)+m%n-1
    if ((dims(im) == m%n).and.periodic(im)) then
      m%b_sb=m%b_s
      m%b_eb=m%b_e 
    else
      m%b_sb=i(im)-m%ng
      m%b_eb=m%b_e+m%ng 
    end if
    m%b_r=m%b_e-m%b_s+1
    m%b_r=m%b_eb-m%b_sb+1
  end do
  ! x index space
  m => self%mesh(1)
  call self%setup_mesh(m,1,mx,xm, xmdn,dxidxup,dxidxdn,dx1d)
  m%b_d = sum(abs(m%b_1d(m%b_s:m%b_e)))/size(m%b_1d(m%b_s:m%b_e))
  ! y index space
  m => self%mesh(2)
  call self%setup_mesh(m,2,my,ym, ymdn,dyidyup,dyidydn,dy1d)
  m%b_d = sum(abs(m%b_1d(m%b_s:m%b_e)))/size(m%b_1d(m%b_s:m%b_e))
  ! z index space
  m => self%mesh(3)
  call self%setup_mesh(m,3,mz,zm, zmdn,dzidzup,dzidzdn,dz1d)
  m%b_d = sum(abs(m%b_1d(m%b_s:m%b_e)))/size(m%b_1d(m%b_s:m%b_e))
  !
  allocate(m%zlincoeff (m%b_sb:m%b_eb,0:5))
  allocate(m%dzlincoeff(m%b_sb:m%b_eb,0:5))
  allocate(m%zuincoeff (m%b_sb:m%b_eb,0:5))
  allocate(m%dzuincoeff(m%b_sb:m%b_eb,0:5))
  call self%one_sided_coeff(m,mz)
  call trace%end()
END SUBROUTINE read_bifrost_mesh

SUBROUTINE setup_mesh(self,m,i, mx, xm, xmdn,dxidxup,dxidxdn,dx1d)
  class(bifrost_mesh_t)      :: self
  class(mesh_t),pointer :: m
  real, dimension(:), pointer :: xm, xmdn, dxidxup, dxidxdn, dx1d
  integer:: i, mx
  ! --------------------------------------------------------------------
  real,dimension(:),allocatable :: tmp1,tmp2,tmp3,tmp4,tmp5
  real :: tmp
  integer:: j
  ! ....................................................................
  
  allocate(tmp1(m%b_sb:m%b_eb))
  allocate(tmp2(m%b_sb:m%b_eb))
  allocate(tmp3(m%b_sb:m%b_eb))
  allocate(tmp4(m%b_sb:m%b_eb))
  allocate(tmp5(m%b_sb:m%b_eb))
  ! lower guard zones
  if (m%b_sb<1) then
    if (periodic(i)) then
      tmp1(m%b_sb:m%b_s-1)= xm     (mx-m%ng+1:mx)
      tmp2(m%b_sb:m%b_s-1)= xmdn   (mx-m%ng+1:mx)
      tmp3(m%b_sb:m%b_s-1)= dxidxup(mx-m%ng+1:mx)
      tmp4(m%b_sb:m%b_s-1)= dxidxdn(mx-m%ng+1:mx)
      tmp5(m%b_sb:m%b_s-1)= dx1d   (mx-m%ng+1:mx)
    else
      do j=1,m%ng
        tmp1(1-j)     = 2*xm      (1) - xm      (1+j)
        tmp         = 2./8.*(15.*xmdn(1  )-10.*xmdn(2  )+3.*xmdn(3  ))
        tmp2(1-j  ) = tmp - xmdn(  j)
        tmp         = 2./8.*(15.*dxidxdn(1  )-10.*dxidxdn(2  )+3.*dxidxdn(3))
        tmp4(1-j)   = tmp - dxidxdn(  j)
        tmp3(1-j)   = 2*dxidxup (1) - dxidxup (1+j)
        tmp5(1-j)   = 2*dx1d    (1) - dx1d    (1+j)
      end do
    end if
  else
    if (dims(i) /= 1) then
      tmp1(m%b_sb:m%b_s-1)=xm     (m%b_sb:m%b_s-1)
      tmp2(m%b_sb:m%b_s-1)=xmdn   (m%b_sb:m%b_s-1)
      tmp3(m%b_sb:m%b_s-1)=dxidxup(m%b_sb:m%b_s-1)
      tmp4(m%b_sb:m%b_s-1)=dxidxdn(m%b_sb:m%b_s-1)
      tmp5(m%b_sb:m%b_s-1)=dx1d   (m%b_sb:m%b_s-1)
    end if
  endif
  ! internal part
  tmp1(m%b_s:m%b_e)=xm     (m%b_s:m%b_e)
  tmp2(m%b_s:m%b_e)=xmdn   (m%b_s:m%b_e)
  tmp3(m%b_s:m%b_e)=dxidxup(m%b_s:m%b_e)
  tmp4(m%b_s:m%b_e)=dxidxdn(m%b_s:m%b_e)
  tmp5(m%b_s:m%b_e)=dx1d   (m%b_s:m%b_e)
  ! upper guard zones
  if (m%b_eb>mx) then
    if (periodic(i)) then
      tmp1(m%b_e+1:m%b_eb)=xm     (1:m%ng)
      tmp2(m%b_e+1:m%b_eb)=xmdn   (1:m%ng)
      tmp3(m%b_e+1:m%b_eb)=dxidxup(1:m%ng)
      tmp4(m%b_e+1:m%b_eb)=dxidxdn(1:m%ng)
      tmp5(m%b_e+1:m%b_eb)=dx1d   (1:m%ng)
    else
      do j=1,m%ng
        tmp1(mx+j) = 2*xm     (mx) - xm  (mx-j)
        tmp        = 2./8.*(15.*xmdn(mx)-10.*xmdn(mx-1)+3.*xmdn(mx-2))
        tmp2(mx+j) =   tmp         - xmdn(mx-j+1)
        tmp3(mx+j) = 2*dxidxup(mx) - dxidxup(mx-j)
        tmp        = 2./8.*(15.*dxidxdn(mx)-10.*dxidxdn(mx-1)+3.*dxidxdn(mx-2))
        tmp4(mx+j) =   tmp         - dxidxdn(mx-j+1)
        tmp5(mx+j) = 2*dx1d   (mx) - dx1d (mx-j)
      end do
    end if
  else
    if (dims(i) /= 1) then
      tmp1(m%b_e+1:m%b_eb)=xm     (m%b_e+1:m%b_eb)
      tmp2(m%b_e+1:m%b_eb)=xmdn   (m%b_e+1:m%b_eb)
      tmp3(m%b_e+1:m%b_eb)=dxidxup(m%b_e+1:m%b_eb)
      tmp4(m%b_e+1:m%b_eb)=dxidxdn(m%b_e+1:m%b_eb)
      tmp5(m%b_e+1:m%b_eb)=dx1d   (m%b_e+1:m%b_eb)
    end if
  end if
  
  allocate(m%b_m    (m%b_sb:m%b_eb))
  allocate(m%b_mdn  (m%b_sb:m%b_eb))
  allocate(m%b_1d   (m%b_sb:m%b_eb))
  allocate(m%b_didup(m%b_sb:m%b_eb))
  allocate(m%b_diddn(m%b_sb:m%b_eb))
  m%b_m(:)     = tmp1(:)
  m%b_mdn(:)   = tmp2(:)
  m%b_didup(:) = tmp3(:)
  m%b_diddn(:) = tmp4(:)
  m%b_1d(:)    = tmp5(:)
  deallocate(tmp1,tmp2,tmp3,tmp4,tmp5)
  !
  m%b_cd = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
  m%b_bd = (-1.-120.*m%b_cd)/24.
  m%b_ad = (1.-3.*m%b_bd-5.*m%b_cd)
  !
END SUBROUTINE setup_mesh


!=======================================================================
!> compute coefficients for one sided derivatives
!=======================================================================
subroutine one_sided_coeff(self,m,mz)
  class(bifrost_mesh_t) :: self
  class(mesh_t), pointer :: m
  integer :: mz
  ! --------------------------------------------------------------------
  integer :: j,k,ordsh
  integer, dimension(5) :: iordl=(/1,1,5,5,5/),iordu=(/1,1,5,5,5/)
  integer, dimension(5) :: dordl=(/1,1,5,5,5/),dordu=(/1,1,5,5,5/)
  real, dimension(0:5) :: gridcoeff
  ! ....................................................................

  m%zlincoeff(:,:)=0
  m%dzlincoeff(:,:)=0
  m%zuincoeff(:,:)=0
  m%dzuincoeff(:,:)=0
  k=m%b_sb
  do j=0,5
    gridcoeff(j)=float(j)+0.5
  enddo
  m%zlincoeff (k,:) = self%stagger_inv(gridcoeff,iordl(k-m%b_sb+1),0,0)
  m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,dordl(k-m%b_sb+1),1,0)
  do j=0,5
    gridcoeff(j)=float(j)-0.5
  enddo
  m%zuincoeff (k,:) = self%stagger_inv(gridcoeff,iordl(k-m%b_sb+1),0,0)
  m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,dordl(k-m%b_sb+1),1,0)
  if (mz >= 2) then
    k=m%b_sb+1
    do j=0,5
      gridcoeff(j)=float(j)-0.5
    enddo
    m%zlincoeff (k,:) = self%stagger_inv(gridcoeff,iordl(k-m%b_sb+1),0,0)
    m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,dordl(k-m%b_sb+1),1,0)
    ordsh=0
    if(iordl(k-m%b_sb+1)==1) ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%zuincoeff(k,:) = self%stagger_inv(gridcoeff,iordl(k-m%b_sb+1),0,ordsh)
    ordsh=0
    if(dordl(k-m%b_sb+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,dordl(k-m%b_sb+1),1,ordsh)
    k=m%b_sb+2
    ordsh=0
    if(iordl(k-m%b_sb+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%zlincoeff(k,:) = self%stagger_inv(gridcoeff,iordl(k-m%b_sb+1),0,ordsh)
    ordsh=0
    if(dordl(k-m%b_sb+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,dordl(k-m%b_sb+1),1,ordsh)
    do j=0,5
      gridcoeff(j)=float(j)-2.5
    enddo
    m%zuincoeff (k,:) = self%stagger_inv(gridcoeff,5,0,0)
    m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,5,1,0)
!
    do k=m%b_sb+3,m%b_eb-3
      do j=0,5
        gridcoeff(j)=float(j)-2.5
      enddo
      m%zlincoeff (k,:) = self%stagger_inv(gridcoeff,5,0,0)
      m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,5,1,0)
      do j=0,5
        gridcoeff(j)=float(j)-2.5
      enddo
      m%zuincoeff (k,:) = self%stagger_inv(gridcoeff,5,0,0)
      m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,5,1,0)
    enddo
!
    k=m%b_eb-2
    ordsh=0
    if(iordu(m%b_eb-k+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%zuincoeff(k,:) = self%stagger_inv(gridcoeff,iordu(m%b_eb-k+1),0,ordsh)
    ordsh=0
    if(dordu(m%b_eb-k+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,dordu(m%b_eb-k+1),1,ordsh)
    do j=0,5
      gridcoeff(j)=float(j)-2.5
    enddo
    m%zlincoeff (k,:) = self%stagger_inv(gridcoeff,5,0,0)
    m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,5,1,0)
    k=m%b_eb-1
    do j=0,5
      gridcoeff(j)=float(j)-0.5
    enddo
    m%zuincoeff (k,:) = self%stagger_inv(gridcoeff,iordu(m%b_eb-k+1),0,0)
    m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,dordu(m%b_eb-k+1),1,0)
    ordsh=0
    if(iordu(m%b_eb-k+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%zlincoeff(k,:) = self%stagger_inv(gridcoeff,iordu(m%b_eb-k+1),0,ordsh)
    ordsh=0
    if(dordl(m%b_eb-k+1)==1)ordsh=1
    do j=0,5
      gridcoeff(j)=float(j)-1.5+ordsh
    enddo
    m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,dordu(m%b_eb-k+1),1,ordsh)
    k=m%b_eb
    do j=0,5
      gridcoeff(j)=float(j)+0.5
    enddo
    m%zuincoeff (k,:) = self%stagger_inv(gridcoeff,iordu(m%b_eb-k+1),0,0)
    m%dzuincoeff(k,:) = self%stagger_inv(gridcoeff,dordu(m%b_eb-k+1),1,0)
    do j=0,5
      gridcoeff(j)=float(j)-0.5
    enddo
    m%zlincoeff (k,:) = self%stagger_inv(gridcoeff,iordu(m%b_eb-k+1),0,0)
    m%dzlincoeff(k,:) = self%stagger_inv(gridcoeff,dordu(m%b_eb-k+1),1,0)

  endif

end subroutine one_sided_coeff

!=======================================================================
!> Calculate the inverse of the Vandermonde matrix
!> Taken from stagger_omp in the original Oslo Stagger Code
!=======================================================================
function stagger_inv (self,x, n, y, o) result(out)
  class(bifrost_mesh_t):: self
  real, dimension(0:5):: x
  real, dimension(0:5):: out
  integer,intent(in) :: n, y, o
  ! --------------------------------------------------------------------
  real, dimension(0:5):: c, b
  real:: t
  integer:: i,j

  do i=0,n
     c(i) = 0.0
  end do
  c(n) = -x(0)
  do i=1,n
    do j=n-i,n-1
       c(j) = c(j) - x(i)*c(j+1)
    end do
    c(n) = c(n) - x(i)
  end do
  do i=0,5
    out(i)=0.0
  end do
  do i=0,n
    t = 1.0
    b(n) = 1.0
    do j=n,1,-1
      b(j-1) = c(j) + x(i)*b(j)
      t = x(i)*t + b(j-1)
    end do
    out(i+o) = b(y)/t
  end do
end function stagger_inv


END MODULE mesh_interface_mod
