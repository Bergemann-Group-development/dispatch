!> =============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object mhd_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE timer_mod
  USE lagrange_mod
  USE mesh_mod
  USE gpatch_mod
  USE link_mod
  USE rt_solver_mod
  implicit none
  private
  type, public, extends(rt_solver_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure boundary_conditions
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  integer, save:: verbose=1
  integer, save:: dims(3)=63
  integer, save:: order=2
  real, save:: grav=2.75
  character(len=64):: file='snapshot'
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time
  integer, save:: itimer=0
  namelist /experiment_params/ file, grav, dims, order, verbose
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (init_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, experiment_params, iostat=iostat)
      if (io%master) write (*, experiment_params)
    end if
    !$omp end critical (init_cr)
  end if
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%rt_solver_t%init
  call self%set(bits%static)
  self%force_per_unit_mass => force
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  self%periodic(3) = .false.                     ! non-periodic vertically
  call self%initial_condition
  call self%init_bdries
  call self%boundary_conditions
  call trace_end (itimer=itimer)
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if (first_time) then
    first_time = .false.
    call self%initial%cleanup
  end if
  call self%rt_solver_t%update
  call self%boundary_conditions
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  integer:: iy, li(3), ui(3)
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  li = self%mesh%li
  ui = self%mesh%ui
  if (self%mesh(2)%lower_boundary) then
    do iy=self%mesh(2)%lb,self%mesh(2)%lo
      self%mem(:,iy,:,:,self%it,1) = self%mem(:,self%mesh(2)%li,:,:,self%it,1)
    end do
  end if
  if (self%mesh(2)%upper_boundary) then
    do iy=self%mesh(2)%uo,self%mesh(2)%ub
      self%mem(:,iy,:,:,self%it,1) = self%mem(:,self%mesh(2)%ui,:,:,self%it,1)
    end do
  end if
  call trace_end
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Initialize from a legacy raw data file, with auxciliary grid file
!> =============================================================================
SUBROUTINE initial_condition (self)
  class(experiment_t):: self
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), o(3), iv, recl, iy
  logical, save:: first_time=.true.
  integer:: mx, my, mz, ii(9)
  real, save, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn, dxidxup, dxidxdn
  real, save, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn, dyidyup, dyidydn
  real, save, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn, dzidzup, dzidzdn
  real, dimension(:), pointer:: xi, yi, zi, xo, yo, zo
  integer:: nx, ny, nz
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%initial_condition')
  !-----------------------------------------------------------------------------
  !$omp critical (init_cr)
  m => self%mesh
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
  ii=[self%idx%d,self%idx%px,self%idx%py,self%idx%pz, &
      self%idx%e,self%idx%tt,self%idx%bx,self%idx%by,self%idx%bz]
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
  associate (ff=>self%mem(:,:,:,:,self%it,1))
  do iv=1,self%nv
    xo = m(1)%p + m(1)%r + m(1)%h(iv)
    yo = m(2)%p + m(2)%r + m(2)%h(iv)
    zo = m(3)%p + m(3)%r + m(3)%h(iv)
    xi => xm
    yi => ym
    zi => zm
    if (iv==self%idx%px.or.iv==self%idx%bx) xi => xmdn
    if (iv==self%idx%py.or.iv==self%idx%by) yi => ymdn
    if (iv==self%idx%pz.or.iv==self%idx%bz) zi => zmdn
    if (verbose>1) &
      print '(i4,2x,2(2x,3(2x,2f7.2)))', iv, &
      minval(xi), maxval(xi), minval(yi), maxval(yi), minval(zi), maxval(zi), &
      minval(xo), maxval(xo), minval(yo), maxval(yo), minval(zo), maxval(zo)
    if (order==1) then
      call lagrange%trilinear3d (xi, yi, zi, buffer(:,:,:,iv), &
                                 xo, yo, zo,     ff(:,:,:,iv), order)
    else
      call lagrange%interpolate3d (xi, yi, zi, buffer(:,:,:,iv), &
                                   xo, yo, zo,     ff(:,:,:,iv), order)
    end if
    if (verbose>0) &
      print '(i4,i3,2x,3f8.3,1p,2e12.4)', self%id, iv, m%p, minval(ff(:,:,:,iv)), maxval(ff(:,:,:,iv))
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
  end associate
  call trace_end
END SUBROUTINE initial_condition

!> =============================================================================
!> Stellar atmoshpere with gravity pointing in the positive y-direction
!> =============================================================================
FUNCTION force (self) RESULT (ff)
  class(gpatch_t)   :: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3) :: ff
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%force')
  ff(:,:,:,1) = 0.0
  ff(:,:,:,2) = 0.0
  ff(:,:,:,3) = grav
  call trace_end
END FUNCTION force

END MODULE
