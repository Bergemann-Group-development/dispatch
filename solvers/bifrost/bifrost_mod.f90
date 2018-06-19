!> =============================================================================
!> Interface module, with date type procedures that are called from solver_t
!> procedures and which mainly USE Bifrost side module files
!> =============================================================================
MODULE bifrost_mod
  ! Dispatch modules
  USE trace_mod
  USE mesh_interface_mod
  ! Bifrost modules
  USE params_mod
  USE parallel_mod
  USE stagger_mod
  USE bifrost_timestep_mod
  USE mhd_mod
  USE math_mod
  USE adaption_mod
  USE index_mod
  USE work_mod
  implicit none
  private
  type, public:: bifrost_t
    integer :: it
    integer :: rereadpars
    logical, dimension(:), allocatable:: scratches
  contains
    procedure:: alloc_scratches
    procedure:: check_scratches
    procedure:: init
    procedure:: put
    procedure:: get
    !procedure:: dnload
    procedure:: update => update_bifrost
    procedure:: pack
    procedure:: unpack
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
  end type
  type(bifrost_t), public:: bifrost
  type(index_t):: idx
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!> =============================================================================
!> Initialize mem storage that fits the Bifrost MHD array sizes.
!> =============================================================================
SUBROUTINE init (self, n, mesh)
  class(bifrost_t):: self
  class(bifrost_mesh_t):: mesh
  integer:: n(3)
  ! ............................................................................
  call trace%begin ('bifrost_t%init')
  n = shape(e)
  self%rereadpars=rereadpars
  self%it = nstepstart+1
  !allocate(self%scratches(12))
  !self%scratches = .false.
  call mesh%init
  call mesh%get
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Put patch mem data back to Bifrost MHD variables
!> =============================================================================
SUBROUTINE put (self, mem, m, time)
  class(bifrost_t):: self
  real, dimension(:,:,:,:):: mem
  type(bifrost_mesh_t):: m
  real(8):: time
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_t%put')
  r(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)  = mem(:,:,:,idx%d)
  e(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)  = mem(:,:,:,idx%e)
  px(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%px)
  py(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%py)
  pz(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%pz)
  bx(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%bx)
  by(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%by)
  bz(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb) = mem(:,:,:,idx%bz)
  !call self%alloc_scratches
  call trace%end
END SUBROUTINE put

!> =============================================================================
!> Get patch mem data from Bifrost MHD variables
!> =============================================================================
SUBROUTINE get (self, mem, m, time, dtime)
  class(bifrost_t):: self
  type(bifrost_mesh_t):: m
  real, dimension(:,:,:,:):: mem
  real(8):: time, dtime
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_t%get')
  mem(:,:,:,idx%d)  = r(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%e)  = e(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%px) = px(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%py) = py(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%pz) = pz(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%bx) = bx(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%by) = by(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  mem(:,:,:,idx%bz) = bz(m%xsb:m%xeb,m%ysb:m%yeb,m%zsb:m%zeb)
  !call self%check_scratches
  call trace%end
END SUBROUTINE get

!> =============================================================================
!> Download guard zone values
!> =============================================================================
SUBROUTINE dnload (self)
  class(bifrost_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_t%dnload')
  call trace%end
END SUBROUTINE dnload

!> =============================================================================
!> Take one time step, whatever that means.  The patch memory, time and time
!> step should be updated.  Whether the code does RK updates internally, or
!> leaves it to bifrost_mod.f90 calling DISPATCH timestep_m.f90 is optional.
!> =============================================================================
SUBROUTINE update_bifrost (self, time, dtime)
  class(bifrost_t):: self
  real(kind=8):: time, dtime
  logical:: ss
  integer :: counter_start,i,tl_tp=0,tl_tpold=0,tl_cr,tl_cm,tl_passed
  !.............................................................................
  call trace%begin ('bifrost_t%update')
  dtold = dt
  print *, "BIFROST MOD, IT, READPARS", self%it, self%rereadpars
  if (mod(self%it,self%rereadpars) == 0) call re_read_params_file
  if (do_adapt) call adaption_start
  if (do_mhd .and. clean <= it .and. mhdclean > 0) then
    if (.not.periodic(3)) then
      if (mhdclean_ub .and. it.gt.(nstepstart+1)) then
        if (mhdclean_lb) then
          call div_clean(bx,by,bz,1-mb,mz+mb)
        else
          call div_clean(bx,by,bz,1,mz+mb)
        endif
      else
        if (mhdclean_lb .and. it.gt.(nstepstart+1)) then
          call div_clean(bx,by,bz,1-mb,mz)
        else
          call div_clean(bx,by,bz,1,mz)
        endif
      endif
    else
      call div_clean(bx,by,bz,1,mz)
    endif
  endif
!
  call timeit('timestep')
  call timestep
  imax_time=imax
  call timeit('timestep')
!
  if (do_mhd .and. clean <= it .and. mhdclean > 0) then
     if (clean < it) clean = it
     clean=clean+mhdclean
  endif
  !  call system_clock(counter_end)
  if (do_adapt) call adaption_end
!
  call io_timestep
!
  if (mod(it,nscr) == 0) then
    isnap=-1*isnap
    call write_params
    call write_snapshot
    call write_aux
    isnap=-1*isnap
  endif
!
  if (mod(it,nsnap) == 0) then
    isnap=isnap+1
    call write_params
    call write_snapshot
    call write_aux
  endif
!
  if ( t >= tsnap+dtsnap ) then
    tsnap=tsnap+dtsnap
    isnap=isnap+1
    call write_params
    call write_snapshot
    call write_aux
  end if
!
  if (t >= tscr+dtscr ) then
    tscr=tscr+dtscr
    isnap=-1*isnap
    call write_params
    call write_snapshot
    call write_aux
    isnap=-1*isnap
  end if
!
!  call check_lower_bound('ts_end')
  ss=.false.
  if (myid == R_MASTER) then
    if (time_lim > 0.0) then
      tl_tpold=tl_tp
      call system_clock(COUNT=tl_tp)
      if (tl_tp < tl_tpold) then
        tl_passed=tl_passed+(tl_tp+tl_cm-tl_tpold)
      else
        tl_passed=tl_passed+(tl_tp-tl_tpold)
      end if
      if ((1.0*tl_passed)/(1.0*tl_cr) > time_lim) then
        ss=.true.
        write (*,*)
        write (*,*) "Time limit reached - exiting"
        write (*,*)
      end if
    endif
  end if
  self%it = self%it + 1
  time = t
  dtime = dt
  !call MPI_BCAST(ss,1,MPI_LOGICAL,R_MASTER,C_CART,ierr)
  !if (ss) exit
  call trace%end
END SUBROUTINE update_bifrost

!> =============================================================================
!> Pack information about a patch into a buffer.  Let the patch_t%pack take
!> care of most of it, and then append whatever more is needed
!> =============================================================================
SUBROUTINE pack (self)
  class(bifrost_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_t%pack')
  call trace%end
END SUBROUTINE pack

!> =============================================================================
!> Unpack an MPI buffer.  The patch_t%unpack takes care of locating the patch,
!> and unpacking the patch related header, after which remains only to unpack
!> whatever else we appended
!> =============================================================================
SUBROUTINE unpack (self)
  class(bifrost_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('bifrost_t%unpack')
  call trace%end
END SUBROUTINE unpack

SUBROUTINE check_scratches(self)
  class(bifrost_t):: self
  ! ----------------------------------------------------------------------------
  call trace%begin('bifrost_t%check_scratches')
    if (allocated(scr1)) self%scratches(1)=.true.
    if (allocated(scr2)) self%scratches(2)=.true.
    if (allocated(scr3)) self%scratches(3)=.true.
    if (allocated(scr4)) self%scratches(4)=.true.
    if (allocated(scr5)) self%scratches(5)=.true.
    if (allocated(scr6)) self%scratches(6)=.true.
    if (allocated(scr7)) self%scratches(7)=.true.
    if (allocated(scr8)) self%scratches(8)=.true.
    if (allocated(scr9)) self%scratches(9)=.true.
    if (allocated(scr10)) self%scratches(10)=.true.
    if (allocated(scr11)) self%scratches(11)=.true.
    if (allocated(scr12)) self%scratches(12)=.true.
  call trace%end
END SUBROUTINE check_scratches

SUBROUTINE alloc_scratches(self)
  class(bifrost_t):: self
  ! ----------------------------------------------------------------------------
  call trace%begin('bifrost_t%alloc_scratches')
    if (self%scratches( 1).and..not.(allocated( scr1))) allocate( scr1(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 2).and..not.(allocated( scr2))) allocate( scr2(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 3).and..not.(allocated( scr3))) allocate( scr3(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 4).and..not.(allocated( scr4))) allocate( scr4(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 5).and..not.(allocated( scr5))) allocate( scr5(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 6).and..not.(allocated( scr6))) allocate( scr6(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 7).and..not.(allocated( scr7))) allocate( scr7(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 8).and..not.(allocated( scr8))) allocate( scr8(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches( 9).and..not.(allocated( scr9))) allocate( scr9(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches(10).and..not.(allocated(scr10))) allocate(scr10(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches(11).and..not.(allocated(scr11))) allocate(scr11(xsb:xeb,ysb:yeb,zsb:zeb))
    if (self%scratches(12).and..not.(allocated(scr12))) allocate(scr12(xsb:xeb,ysb:yeb,zsb:zeb))
  call trace%end
END SUBROUTINE alloc_scratches

!> =============================================================================
SUBROUTINE void (self, out)
  class (bifrost_t):: self
  real, dimension(:,:,:):: out
END SUBROUTINE void

END MODULE bifrost_mod
