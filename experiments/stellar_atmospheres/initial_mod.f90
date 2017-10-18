!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE index_mod
  USE lagrange_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    logical:: mhd
    integer:: id=0
    integer:: nv
    integer:: iread=0, iwrite=0
    real(8):: time
    character(len=64):: solver
    procedure(single_solenoidal), pointer:: condition=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  integer, save:: iread=0, iwrite=1
CONTAINS

!> =============================================================================
!> Setup an initial state -- from reading a STAGGER code snapshot, a test case,
!> or (FIXME) by reading a snapshot from a production run
!> =============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  integer:: iostat
  logical, save:: first_time=.true.
  !.............................................................................
  character(len=20), save:: type='stagger'
  logical, save:: mhd=.true.
  namelist /initial_params/ type, mhd
  !-----------------------------------------------------------------------------
  call trace_begin('initial_mod:init')
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    mhd = self%mhd
    rewind (io%input); read (io%input, initial_params, iostat=iostat)
    if (io%master) write (*, initial_params)
  end if
  !$omp end critical (init_cr)
  self%solver = trim(solver)
  !-----------------------------------------------------------------------------
  select case (trim(type))
  case ('snapshot')
    self%condition => snapshot                   ! snapshot reader
    self%nv = 8
  case ('stagger')
    self%condition => stagger                    ! legacy STAGGER code reader
    self%nv = 9
  case ('stagger_hd')
    self%condition => stagger                    ! legacy STAGGER code reader
    self%nv = 6
  case ('raw')
    self%condition => raw                        ! generic raw data reader
    self%nv = 5
  case ('single_solenoidal')
    self%condition => single_solenoidal          ! for testing
    self%nv = 8
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE'
    call mpi%abort
  end select
  self%mhd    = mhd .or. self%nv>5
  self%time   = 0d0
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!> =============================================================================
!> Initialize from a snapshot file
!> =============================================================================
SUBROUTINE snapshot (self, m, ff)
  class(initial_t):: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:):: ff
  logical, save:: first_time=.true.
  !.............................................................................
  integer:: dims(4)=0
  namelist /snapshot_params/ iread, iwrite, dims
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%snapshot')
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, snapshot_params)
    if (io%master) write (*, snapshot_params)
  end if
  !$omp end critical (init_cr)
  self%iread = iread
  self%iwrite = iwrite
  self%nv = dims(4)
  call trace_end
END SUBROUTINE snapshot

!> =============================================================================
!> Initialize from a legacy raw data file, with auxciliary grid file
!> =============================================================================
SUBROUTINE stagger (self, m, ff)
  class(initial_t):: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:):: ff
  integer:: n(3), o(3), iv, recl, iy
  logical, save:: first_time=.true.
  integer:: mx, my, mz, ii(9)
  real, save, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn, dxidxup, dxidxdn
  real, save, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn, dyidyup, dyidydn
  real, save, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn, dzidzup, dzidzdn
  real, dimension(:), pointer:: xi, yi, zi, xo, yo, zo
  integer:: nx, ny, nz
  !.............................................................................
  integer, save:: dims(3)=63
  integer, save:: order=2
  type(index_t):: i
  character(len=64):: file='snapshot'
  namelist /stagger_params/ file, dims, iread, iwrite, order
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%stagger')
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, stagger_params)
    if (io%master) write (*, stagger_params)
    n = nint(m%b/m%d)
    n(2) = n(2)+4
    if (.not.all(n==dims)) then
      if (io%master) print *, 'file  dims =', dims, self%nv
      if (io%master) print *, 'setup dims =', n
    end if
    !---------------------------------------------------------------------------
    ! At first call, allocate a buffer and read the entire file -- this is OK since
    ! we only use one thread per rank to initialize, computing also how many calls
    ! will be made, so the buffer can be deallocated when all patches are done
    !---------------------------------------------------------------------------
    allocate (buffer(dims(1),dims(2),dims(3),self%nv))
    recl = 4*product(dims)                                           ! FIXME
    if (io%master) print *,'reading MHD data from ' //trim(file)//'.dat'
    open (io_unit%data, file=trim(file)//'.dat', form='unformatted', access='direct', recl=recl)
    i%tt = 9
    ii=[i%d,i%px,i%py,i%pz,i%e,i%tt,i%bx,i%by,i%bz]
    do iv=1,self%nv
      read (io_unit%data, rec=iv) buffer(:,:,:,ii(iv))
    end do
    close (io_unit%data)
    !---------------------------------------------------------------------------
    ! Read the vertical meshs from the .msh file
    !---------------------------------------------------------------------------
    open (io_unit%data, file=trim(file)//'.msh', form='unformatted')
    ! read dimensions
      read(io_unit%data) mx, my, mz
      if (.not. all(dims==[mx,my,mz])) then
        if (io%master) print *, ' file dims =', dims
        if (io%master) print *, ' mesh dims =', m%n, m%gn
        if (io%master) print *, 'setup dims =', nint(m%b/m%d), m%b, m%d
      end if
      if (io%master) print *,'reading mesh data from ' //trim(file)//'.msh', mx, my, mz
      ! allocate single-precision variables for reading
      allocate( dxm(mx), dxmdn(mx), xm(mx), xmdn(mx), dxidxup(mx), dxidxdn(mx) )
      allocate( dym(my), dymdn(my), ym(my), ymdn(my), dyidyup(my), dyidydn(my) )
      allocate( dzm(mz), dzmdn(mz), zm(mz), zmdn(mz), dzidzup(mz), dzidzdn(mz) )
      ! read single-precision variables
      read(io_unit%data) dxm, dxmdn, xm, xmdn, dxidxup, dxidxdn
      read(io_unit%data) dym, dymdn, ym, ymdn, dyidyup, dyidydn
      read(io_unit%data) dzm, dzmdn, zm, zmdn, dzidzup, dzidzdn
      if (io%verbose>1) then
        do iy=1,size(dym)
          print *, 'initial_t%stagger: min(dx),dy,min(dz)', &
            iy, minval(dxm), dym(iy), minval(dzm)
        end do
      end if
      print *, 'initial_t%stagger: min(dx),min(dy),min(dz)', &
            minval(dxm), minval(dym), minval(dzm)
      print *,'initial_t%stagger: ny for max res:', int((ym(my)-ym(1))/minval(dym))+1
    close (io_unit%data)
  end if
  !$omp end critical (init_cr)
  self%nv = self%nv
  !-----------------------------------------------------------------------------
  ! Pick out the piece that belongs to this patch, based on its position
  !-----------------------------------------------------------------------------
  nx = m(1)%gn
  ny = m(2)%gn
  nz = m(3)%gn
  allocate (xo(nx))
  allocate (yo(ny))
  allocate (zo(nz))
  do iv=1,self%nv
    xo = m(1)%p + m(1)%r + m(1)%h(iv)
    yo = m(2)%p + m(2)%r + m(2)%h(iv)
    zo = m(3)%p + m(3)%r + m(3)%h(iv)
    xi => xm
    yi => ym
    zi => zm
    if (iv==i%px.or.iv==i%bx) xi => xmdn
    if (iv==i%py.or.iv==i%by) yi => ymdn
    if (iv==i%pz.or.iv==i%bz) zi => zmdn
    !print '(i4,2x,2(2x,3(2x,2f7.2)))', iv, &
      !minval(xi), maxval(xi), minval(yi), maxval(yi), minval(zi), maxval(zi), &
      !minval(xo), maxval(xo), minval(yo), maxval(yo), minval(zo), maxval(zo)
    if (order==1) then
      call lagrange%trilinear3d (xi, yi, zi, buffer(:,:,:,iv), &
                                 xo, yo, zo,     ff(:,:,:,iv), order)
    else
      call lagrange%interpolate3d (xi, yi, zi, buffer(:,:,:,iv), &
                                   xo, yo, zo,     ff(:,:,:,iv), order)
    end if
    !print '(i4,i3,2x,3f8.3,1p,2e12.4)', self%id, iv, m%p, minval(ff(:,:,:,iv)), maxval(ff(:,:,:,iv))
  end do
  deallocate (xo, yo, zo)
  if (io%verbose>0) print *, 'initial_t%stagger: id, pos =', self%id, m%p
  !-----------------------------------------------------------------------------
  ! Diagnostic output
  !-----------------------------------------------------------------------------
  if (io%verbose > 0) then
    do iv=1,self%nv
      write (io_unit%log,'(a,3i4,i5,1p,2g12.3)') 'ix,iy,iz,fmin,fmax =', o, iv, &
        minval(ff(m(1)%li:m(1)%ui,m(2)%lb:m(2)%ub,m(3)%li:m(3)%ui,iv)), &
        maxval(ff(m(1)%li:m(1)%ui,m(2)%lb:m(2)%ub,m(3)%li:m(3)%ui,iv))
    end do
  end if
  call trace_end
END SUBROUTINE stagger

!> =============================================================================
!> Initialize from a raw data file
!> =============================================================================
SUBROUTINE raw (self, m, ff)
  class(initial_t):: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:):: ff
  integer:: n(3), i(3), iv, recl
  logical, save:: first_time=.true.
  !.............................................................................
  character(len=64):: file='../../data/exa256/raw0.02.dat'
  integer, save:: dims(4)=[256,256,256,4]
  real, save:: t0=0.02
  namelist /raw_params/ file, dims, t0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%raw')
  !-----------------------------------------------------------------------------
  ! At first call, allocate a buffer and read the entire file -- this is OK since
  ! we only use one thread per rank to initialize, computing also how many calls
  ! will be made, so the buffer can be deallocated when all patches are done
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, raw_params)
    if (io%master) write (*, raw_params)
    n = nint(m%b/m%d)
    if (.not.all(n==dims(1:3))) then
      if (io%master) print *, 'file  dims =', dims
      if (io%master) print *, 'setup dims =', n, m%b, m%d
    end if
    allocate (buffer(dims(1),dims(2),dims(3),dims(4)))
    recl = 4*product(n)                                         ! FIXME
    open (io_unit%data, file=file, form='unformatted', access='direct', recl=recl)
    do iv=1,dims(4)
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  !$omp end critical (init_cr)
  self%nv = dims(4)+1
  i = 1 + nint(m%p/m%d) - (m%n+1)/2                             ! index offsets
  ff = 0.0
  do iv=1,dims(4)
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,iv) = buffer(i(1):i(1)+m(1)%n-1, &
                                    i(2):i(2)+m(2)%n-1, &
                                    i(3):i(3)+m(3)%n-1,iv)
  end do
  ff(:,:,:,5) = 0.0
  self%time = t0
  call trace_end
END SUBROUTINE raw

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
SUBROUTINE single_solenoidal (self, m, ff)
  class(initial_t)      :: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t) :: random
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, phi(3)=0d0
  logical, save  :: first_time=.true.
  !.............................................................................
  integer, save:: seed=-33
  real, save:: k(3)=1.0, a0(3)=0.1, d0=1.0, e0=1.0, u0(3)=0.0, b0(3)=0.0
  namelist /solenoidal_params/ seed, k, a0, d0, e0, u0, b0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_mod:single_solenoidal')
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, solenoidal_params)
    if (io%master) write (*, solenoidal_params)
  end if
  !$omp end critical (init_cr)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,1) = d0
        ff(ix,iy,iz,2) = e0
        ff(ix,iy,iz,3) = a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,4) = a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,5) = a0(3)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,6) = b0(1)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,7) = b0(2)*(sin(kx)+sin(ky))
        ff(ix,iy,iz,8) = b0(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE single_solenoidal

END MODULE
