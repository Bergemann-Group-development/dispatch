MODULE background_io_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE omp_timer_mod
  USE mpi_mod
  USE omp_lib
  implicit none
  private
  type, public:: background_io_t
  contains
    procedure, nopass:: output
    procedure, nopass:: input
  end type
  integer, save:: snapshot = -1
  integer, save:: record = 0
  type(background_io_t), public:: background_io
CONTAINS

!===============================================================================
!> Interface to prepare background writing into a specific record of the
!> snapshot file for one rank and snapshot
!===============================================================================
SUBROUTINE output (self)
  class(patch_t):: self
  integer:: m(6), l(3), u(3), snap, rec, it(1)
  real(var_k), pointer:: buffer(:,:,:,:)
  character(len=64) filename
  !-----------------------------------------------------------------------------
  ! In brief critical region, check if starting a new snapshot
  !-----------------------------------------------------------------------------
  call trace%begin ('background_t%output')
  call io%lock%set ('background_t%output')
  if (snapshot < self%iout) then
    snapshot = self%iout
    record = 0
    write (filename,'(a,i5.5,"/snapshot_*.dat")') trim(io%run_dir), snapshot
  end if
  io%nml_version = max(4,io%nml_version)
  record = record + 1
  rec = record
  snap = snapshot
  self%record = record
  it = minloc(abs(self%t-self%out_next))
  call io%lock%unset ('background_t%output')
  write (stdlog,*) &
    'background_io_t%output: snapshot,record,it=', snapshot, record, it
  !-----------------------------------------------------------------------------
  ! Take a copy of the output data, and start a background task that writes it
  !-----------------------------------------------------------------------------
  m = shape(self%mem)
  l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
  u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
  m(1:3) = u-l+1
  allocate (buffer(m(1),m(2),m(3),m(4)))
  buffer = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,it(1),1)
  !-----------------------------------------------------------------------------
  ! Write out with background task
  !-----------------------------------------------------------------------------
  call log (rec, 1)
  !call omp_set_max_active_levels (3)
  !$omp task firstprivate (buffer, snap, rec)
    call write_out (buffer, snap, rec)
    deallocate (buffer)
  !$omp end task
  call log (rec, 2)
  call trace%end ()
END SUBROUTINE output

!===============================================================================
SUBROUTINE log (rec, step)
  integer:: step, rec
  write (io_unit%mpi,'(a,2i5,5x,"time(",i1,"):",f12.6)') &
    'output:', snapshot, rec, step, wallclock()
END SUBROUTINE log

!===============================================================================
!> Write results to disk, placing the files belonging to each snapshot in a
!> separate directory, to avoid the slowdown sometimes associated with having
!> many files in the same directory.
!===============================================================================
SUBROUTINE write_out (buffer, snap, rec)
  real:: buffer(:,:,:,:)
  integer:: snap, rec
  !.............................................................................
  character(len=64) filename
  integer:: unit
  !-----------------------------------------------------------------------------
  call trace%begin('background_io%write_out')
  call log (rec, 3)
  !-----------------------------------------------------------------------------
  ! Write a record with the data cube to a direct access file
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/snapshot_",i5.5,".dat")') &
    trim(io%run_dir), snap, mpi%rank
  call io%lock%set ('background_io_t%write_out')
  open (newunit=unit, file=trim(filename), form='unformatted', &
        recl=kind(buffer)*product(shape(buffer)), &
        access='direct', status='unknown')
  write (unit, rec=rec) buffer
  close (unit)
  call io%lock%unset ('background_io_t%write_out')
  !-----------------------------------------------------------------------------
  call log (rec, 4)
  call trace%end()
END SUBROUTINE write_out

!===============================================================================
!> Read a record with a data cube from a direct access file
!===============================================================================
SUBROUTINE input (self, ok)
  class(patch_t):: self
  logical:: ok
  !.............................................................................
  character(len=64) filename
  real, pointer:: buffer(:,:,:,:)
  integer:: unit, l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace%begin('background_io_t%input')
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/snapshot_",i5.5,".dat")') &
    trim(io%inputdir), io%restart, self%rank
  l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
  u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
  buffer => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,self%it,1)
  open (newunit=unit, file=trim(filename), form='unformatted', &
        recl=kind(buffer)*product(shape(buffer)), &
        access='direct', status='unknown')
  read (unit, rec=self%record) buffer
  close (unit)
  !-----------------------------------------------------------------------------
  call log (self%record, 5)
  ok = .true.
  call trace%end ()
END SUBROUTINE input

END MODULE background_io_mod
