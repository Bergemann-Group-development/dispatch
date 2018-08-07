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
  USE solver_mod
  implicit none
  private
  type, public, extends(rt_solver_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure initial_condition
    procedure boundary_conditions
  end type
  real, save, allocatable, dimension(:,:,:,:):: buffer
  integer, save:: verbose=1
  integer, save:: dims(3)=63
  integer, save:: order=2
  real, save:: grav=2.75
  character(len=64):: file='mhd63'
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
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
  if (.not. io%master) verbose=0
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
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if (first_time) then
    first_time = .false.
    if (allocated(buffer)) &
      deallocate (buffer)
  end if
  call self%solver_t%update
  call self%rt_solver_t%update
  call self%boundary_conditions
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> For now, uses initiallu frozen BCs, which will work for at least a while --
!> enough to develope and validate the RT
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  integer:: iz, m(6)
  real, save, dimension(:,:,:,:), pointer:: lower => null()
  real, save, dimension(:,:,:,:), pointer:: upper => null()
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  m = shape(self%mem)
  m(3) = self%ng(3)
  if (self%mesh(3)%lower_boundary) then
    if (.not.associated(lower)) then
      allocate (lower(m(1),m(2),m(3),m(4)))
      lower = self%mem(:,:,lb(3):lo(3),:,self%it,1)
      if (verbose > 0) &
        print *,'LBC',self%id,self%position, &
        minval(lower(:,:,:,1)),maxval(lower(:,:,:,1))
    else
      self%mem(:,:,lb(3):lo(3),:,self%it,1) = lower
    end if
  end if
  if (self%mesh(3)%upper_boundary) then
    if (.not.associated(upper)) then
      allocate (upper(m(1),m(2),m(3),m(4)))
      upper = self%mem(:,:,uo(3):ub(3),:,self%it,1)
      if (verbose > 0) &
        print *,'UBC',self%id,self%position, &
        minval(upper(:,:,:,1)),maxval(upper(:,:,:,1))
    else
      self%mem(:,:,uo(3):ub(3),:,self%it,1) = upper
    end if
  end if
  call trace_end
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Initialize from a legacy raw data file, with auxciliary grid file
!> =============================================================================
SUBROUTINE initial_condition (self)
  class(experiment_t):: self
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), o(3), iv, recl, ix, iy, iz
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
  m => self%mesh
  n = nint(m%b/m%d)
  n(3) = n(3)+4
  !$omp critical (init_cr)
  if (first_time) then
    first_time = .false.
    if (io%master) &
      print *,'reading MHD data from ' //trim(file)//'.*'
    if (.not.all(n==dims)) then
      if (io%master) print *, 'file  dimensions =', dims, self%nv
      if (io%master) print *, 'setup dimensions =', n
    end if
    !---------------------------------------------------------------------------
    ! Read the vertical meshs from the .msh file
    !---------------------------------------------------------------------------
    open (io_unit%data, file=trim(file)//'.msh', form='unformatted')
    ! read dimensions
    read(io_unit%data) mx, mz, my
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
    read(io_unit%data) dym, dymdn, ym, ymdn, dyidyup, dyidydn
    read(io_unit%data) dzm, dzmdn, zm, zmdn, dzidzup, dzidzdn
    read(io_unit%data) dxm, dxmdn, xm, xmdn, dxidxup, dxidxdn
    if (verbose>1) then
      do iz=1,size(dzm)
        print *, 'initial_t%stagger: iz, min(dx),min(dy), dz', &
          iy, minval(dxm), minval(dym), dzm(iz)
      end do
    end if
    print *, 'experiment_t%ICs: min(dx),min(dy),min(dz)', &
          minval(dxm), minval(dym), minval(dzm)
    print *,'experiment_t%ICs: ny for max res:', int((zm(my)-zm(1))/minval(dzm))+1
    close (io_unit%data)
    !---------------------------------------------------------------------------
    ! At first call, allocate a buffer and read the entire file -- this is OK since
    ! we only use one thread per rank to initialize, computing also how many calls
    ! will be made, so the buffer can be deallocated when all patches are done
    !---------------------------------------------------------------------------
    allocate (buffer(dims(1),dims(2),dims(3),self%nv))
    recl = 4*product(dims)                                           ! FIXME
    if (io%master) &
      print *,'reading MHD data from ' //trim(file)//'.dat'
    open (io_unit%data, file=trim(file)//'.dat', form='unformatted', access='direct', recl=recl)
    ii=[self%idx%d,self%idx%py,self%idx%pz,self%idx%px, &
        self%idx%e,self%idx%tt,self%idx%by,self%idx%bz,self%idx%bx]
    !---------------------------------------------------------------------------
    ! Read data, and swap indices, so z becomes vertical
    !---------------------------------------------------------------------------
    do iv=1,min(self%nv,9)
      read (io_unit%data, rec=iv) (((buffer(ix,iy,iz,ii(iv)), iy=1,my), iz=1,mz), ix=1,mx)
      if (verbose > 0) &
        print *, iv, minval(buffer(:,:,:,ii(iv))), maxval(buffer(:,:,:,ii(iv)))
    end do
    close (io_unit%data)
    if (verbose > 0) then
      print *,'input model:'
      do iz=1,dims(3)
        print '(i4,1p,6e12.3)',iz, (buffer(1,1,iz,ii(iv)), iv=1,6)
      end do
    end if
  end if
  !$omp end critical (init_cr)
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
  if (verbose>0) print *, 'experiment_t%IC: id, pos =', self%id, m%p
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
    if (verbose>2) &
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
  end do
  deallocate (xo, yo, zo)
  !-----------------------------------------------------------------------------
  ! Diagnostic output
  !-----------------------------------------------------------------------------
  if (verbose>0) then
    do iv=1,self%nv
      o = maxloc(ff(:,:,:,iv))
      print '(a,i3,3i4,2x,3f8.3,1p,2(2x,2e12.4))', &
        'iv,ix,iy,iz,pos,bufmin,bufmax,fmin,fmax =', iv, o, m%p, &
        minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv)), &
        minval(    ff(:,:,:,iv)), maxval(    ff(:,:,:,iv))
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
