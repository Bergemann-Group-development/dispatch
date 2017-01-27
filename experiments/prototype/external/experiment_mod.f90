!> =============================================================================
!> A generic interface to an external code, which needs to supply the procedures
!> exemplified in code/dispatch_mod.f90 -- which is the ONLY new file needed in
!> the external code. The data to be exchanged is defined in interface_mod.f90
!>
!> There are two cases / possible choices:
!>
!> 1) The external code handles the entire update itself.  In this case it is 
!>    also responsible for choosing the time step (must have an internal Courant
!>    condition)
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE interface_mod
  USE mpi_mesg_mod
  USE patch_mod
  USE dispatch_mod
  implicit none
  private
  type, public, extends(patch_t):: experiment_t
  contains
    procedure:: begin
    procedure:: init
    procedure:: connect
    procedure:: disconnect
    procedure:: update
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
  call code_begin
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
  self%nw = 1                                   ! no time derivatives allocated
  call self%patch_t%init                        ! allocate memory and set size
  call code_init                                ! ICs, from scratch or snapshot
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
!> =============================================================================
SUBROUTINE connect (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%connect')
  interface%mem      => self%mem(:,:,:,:,self%it,1)
  interface%position =  self%position
  interface%size     =  self%size
  interface%time     =  self%time
  interface%dtime    =  self%dtime
  interface%id       =  self%id
  interface%log_unit =  io_unit%log
  call trace%end
END SUBROUTINE connect

!> =============================================================================
!> Copy back those dynamic variables that have changed
!> =============================================================================
SUBROUTINE disconnect (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%connect')
  self%time     = interface%time
  self%dtime    = interface%dtime
  self%u_max    = interface%u_max
  call trace%end
END SUBROUTINE disconnect

!> =============================================================================
!> Update the external code with one time step.  Note that this requires the
!> code_update to set the interface_t%u_max value, which is used by the Courant
!> condition (patch_t%courant_condition); cf. code/dispatch_mod.f90
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real, dimension(:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%connect
  call code_update
  call self%disconnect
  call self%patch_t%update
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Give the external code a chance to do I/O.  Note that this may not be 
!> necessary, if calling self%patch_t%output is enough.
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(experiment_t):: self
  character(len=*), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%output')
  call code_output (self%out_time)
  call self%patch_t%output
  call trace%end
END SUBROUTINE output

!> =============================================================================
!> Let the external code print status info -- typically a line in a log file
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(experiment_t):: self
  integer, optional:: nq, ntask
  character(len=*), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%info')
  call code_info
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
  call code_pack (mesg%buffer, mesg%nbuf, mesg%ibuf)
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
  call code_unpack (mesg%buffer, mesg%nbuf, mesg%ibuf)
  call trace%end
END SUBROUTINE unpack


!> =============================================================================
SUBROUTINE end (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%end')
  call code_end
  call trace%end
END SUBROUTINE end

!> =============================================================================
SUBROUTINE void (self, out)
  class (experiment_t):: self
  real, dimension(:,:,:):: out
END SUBROUTINE

END MODULE experiment_mod
