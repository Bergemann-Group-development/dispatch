!> =============================================================================
!> A generic interface to an external code, which needs to supply the procedures
!> exemplified in code/dispatch_mod.f90 -- which is the ONLY new file needed in
!> the external code. The data to be exchanged is defined in interface_mod.f90
!>
!> There are two cases / possible modes of operation:
!>
!> 1) The external code handles the entire update itself.  In this case it is 
!>    also responsible for choosing the time step (must have an internal Courant
!>    condition). This is the simplest alternative.
!>
!> 2) The external code computes time derivatives, while DISPATCH takes care of
!>    updating. This is the preferred mode, since DISPATCH in any case needs the
!>    time slice information. This mode therefore needs less memory, but if one
!>    of the memory conserving RK schemes is used then mode 1) might have
!>    essentially the same low memory requirement.
!>
!> In both modes DISPATCH needs much less memory than the original code, since
!> scratch arrays (which should be allocated and deallocaed in each time step)
!> are only needed for active threads.  The number of patches per thread should
!> normally be of the order 50, and effectively speaking the scratch memory need
!> is thus reduced by this factor.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE trace_mod
  USE gpatch_mod
  USE interface_mod
  USE mpi_mesg_mod
  USE patch_mod
  USE timestep_mod
  USE dispatch_mod
  USE download_mod
  implicit none
  private
  type, public, extends(gpatch_t):: experiment_t
  contains
    procedure:: begin
    procedure:: init
    procedure:: connect
    procedure:: disconnect
    procedure:: dnload
    procedure:: update => update_exp
    procedure:: output
    procedure:: info
    procedure:: pack
    procedure:: unpack
    procedure:: end
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
  end type
  type(experiment_t), public:: experiment
CONTAINS

!> =============================================================================
SUBROUTINE begin (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call bifrost_begin
  call trace%end
END SUBROUTINE begin

!> =============================================================================
!> Interface to let the external code initialize a patch
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  self%nv = 8                                   ! set number of variables
  self%nt = 7                                   ! set number of time slots
  self%ng = 4                                   ! set number of guard zones
  self%nw = 2                                   ! time derivatives allocated
  call self%patch_t%init                        ! allocate memory and set size
  interface%mem => self%mem(:,:,:,:,self%it,:)
  call self%connect
  call bifrost_init                                ! ICs, from scratch or snapshot
  !call timestep%init
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Setup connection to the patch values via interface pointers, and copy any
!> dynamic scalar variables to the interface module, where pointers point to
!> targets in the external code.
!>
!> Note that the external code never needs to know about the multiple time slots
!> used in DISPATCH.  It gets pointers to the variables, and is responsible 
!> either for updating them, or for producing time derivatices, allowing in that
!> case DISPATCH to handle the updating.
!>
!> The required order of execution is that dispatch_mod::connect_code must first
!> set up the pointer connection from the interface to the code values, before
!> the conndect procedure below can copy values from DISPATCH to the code.
!> =============================================================================
SUBROUTINE connect (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%connect')
  interface%mem => self%mem(:,:,:,:,self%it,:)
  call bifrost_connect
  interface%position   =  self%position
!  interface%size       =  self%size
  interface%time       =  self%time
  interface%dtime      =  self%dtime
  interface%id         =  self%id
  !interface%mpi_rank   =  mpi%rank
  interface%log_unit   =  io_unit%log
  !interface%input_unit =  io%input
  call trace%end
END SUBROUTINE connect

!> =============================================================================
!> Copy back those dynamic variables that have changed
!> =============================================================================
SUBROUTINE disconnect (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%disconnect')
  self%time     = interface%time
  self%dtime    = interface%dtime
  self%u_max    = interface%u_max
  call trace%end
END SUBROUTINE disconnect

!> =============================================================================
!> Download guard zone values
!> =============================================================================
SUBROUTINE dnload (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  call download%download_link (self%link)
  call trace%end
END SUBROUTINE dnload

!> =============================================================================
!> Update the external code with one time step.  Note that this requires the
!> bifrost_update to set the interface_t%u_max value, which is used by the Courant
!> condition (patch_t%courant_condition); cf. code/dispatch_mod.f90
!> =============================================================================
SUBROUTINE update_exp (self)
  class(experiment_t):: self
  real, dimension(:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%connect
  call bifrost_update
  call self%disconnect
  call self%patch_t%courant_condition
  !call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call trace%end
END SUBROUTINE update_exp

!> =============================================================================
!> Give the external code a chance to do I/O.  Note that this may not be 
!> necessary, if calling self%patch_t%output is enough.
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(experiment_t):: self
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%output')
  call bifrost_output (io%out_time)
  call self%patch_t%output
  call trace%end
END SUBROUTINE output

!> =============================================================================
!> Let the external code print status info -- typically a line in a log file
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(experiment_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%info')
  call bifrost_info
  call trace%end
END SUBROUTINE info

!> =============================================================================
!> Pack information about a patch into a buffer.  Let the patch_t%pack take 
!> care of most of it, and then append whatever more is needed
!> =============================================================================
SUBROUTINE pack (self, mesg)
  class(experiment_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%pack')
  call self%patch_t%pack (mesg)
  !call bifrost_pack (mesg%buffer, mesg%nbuf, mesg%ibuf)
  call trace%end
END SUBROUTINE pack

!> =============================================================================
!> Unpack an MPI buffer.  The patch_t%unpack takes care of locating the patch,
!> and unpacking the patch related header, after which remains only to unpack
!> whatever else we appended
!> =============================================================================
SUBROUTINE unpack (self, mesg)
  class(experiment_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%unpack')
  call self%patch_t%unpack (mesg)
  !call bifrost_unpack (mesg%buffer, mesg%nbuf, mesg%ibuf)
  call trace%end
END SUBROUTINE unpack


!> =============================================================================
SUBROUTINE end (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%end')
  call bifrost_end
  call trace%end
END SUBROUTINE end

!> =============================================================================
SUBROUTINE void (self, out)
  class (experiment_t):: self
  real, dimension(:,:,:):: out
END SUBROUTINE

END MODULE experiment_mod
