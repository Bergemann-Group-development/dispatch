MODULE background_io_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE task_mod
  USE kinds_mod
  USE omp_timer_mod
  USE mpi_mod
  USE omp_lib
  USE snapshots_mod
  implicit none
  private
  type, public:: background_io_t
  contains
    procedure, nopass:: output
    procedure, nopass:: input
  end type
  integer, save:: verbose = 0
  type(background_io_t), public:: background_io
CONTAINS

!===============================================================================
!> Interface to prepare background writing into a specific record of the
!> snapshot file for one rank and snapshot
!===============================================================================
SUBROUTINE output (self)
  class(patch_t):: self
  class(task_t), pointer:: task
  integer:: l(3), u(3), snap, rec, it(1), recl, unit
  real(var_k), pointer:: buffer(:,:,:,:)
  real(8):: start
  character(len=64) filename
  class(snapshot_t), pointer:: snapshot
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! In brief critical region, check if starting a new snapshot
  !-----------------------------------------------------------------------------
  call trace%begin ('background_t%output', itimer=itimer)
  call io%lock%set ('background_t%output')
  io%nml_version = max(4,io%nml_version)
  snap = self%iout
  it = minloc(abs(self%t-self%out_next))
  !-----------------------------------------------------------------------------
  ! Take a copy of the output data, and start a background task that writes it
  !-----------------------------------------------------------------------------
  l = merge (self%mesh%lb, self%mesh%li, io%guard_zones)
  u = merge (self%mesh%ub, self%mesh%ui, io%guard_zones)
  allocate (buffer(l(1):u(1),l(2):u(2),l(3):u(3),self%mv))
  if (io%do_generic) then
    task => self%link%task
    call task%allocate_generic()
    call task%to_generic ()
    buffer = self%generic
    call task%deallocate_generic()
  else
    buffer = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1:self%mv,self%it,1)
  end if
  if (verbose > 2) then
    call self%stats (buffer, 'background_io')
  end if
  !-----------------------------------------------------------------------------
  ! Open a direct access file if not already open
  !-----------------------------------------------------------------------------
  snapshot => snapshots%get (self%iout)
  if (snapshot%data_unit == 0) then
    snapshot%record = 0
    write (filename,'(a,i5.5,"/snapshot_",i5.5,".dat")') &
      trim(io%run_dir), snap, mpi%rank
    recl = kind(buffer)*product(shape(buffer))
    open (newunit=snapshot%data_unit, file=trim(filename), form='unformatted', &
        recl=recl,  access='direct', status='unknown')
    if (verbose > 0) &
      write (stdout,'(1x,a,3i5,2x,a)') 'background_io: unit,iout,rec,file=', &
        snapshot%data_unit, self%iout, snapshot%record, trim(filename)
  end if
  unit = snapshot%data_unit
  snapshot%record = snapshot%record + 1
  self%record = snapshot%record
  rec = snapshot%record
  if (verbose > 1) &
    write (stdout,*) 'background_io_t%output: unit,iout,rec=', &
      snapshot%data_unit, snapshot%iout, snapshot%record
  call io%lock%unset ('background_t%output')
  !.............................................................................
  if (verbose > 1) &
    write (stdlog,*) &
      'background_io_t%output: snapshot,record,it=', self%iout, rec, it
  !-----------------------------------------------------------------------------
  ! Write out with background task
  !-----------------------------------------------------------------------------
  !$omp task firstprivate (buffer, snap, rec, start, unit)
    start = wallclock()
    call write_out (buffer, snap, rec, unit, start)
    deallocate (buffer)
  !$omp end task
  call trace%end (itimer)
END SUBROUTINE output

!===============================================================================
SUBROUTINE log (rec, step, snapshot)
  integer:: step, rec, snapshot
  if (verbose > 0) then
    write (io_unit%mpi,'(a,2i5,5x,"time(",i1,"):",f12.6)') &
      'output:', snapshot, rec, step, wallclock()
  end if
END SUBROUTINE log

!===============================================================================
!> Write results to disk, placing the files belonging to each snapshot in a
!> separate directory, to avoid the slowdown sometimes associated with having
!> many files in the same directory.
!===============================================================================
SUBROUTINE write_out (buffer, snap, rec, unit, start)
  real:: buffer(:,:,:,:)
  integer:: snap, rec, unit
  real(8):: start
  !.............................................................................
  character(len=64) filename
  integer:: recl
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('background_io%write_out', itimer=itimer)
  if (verbose > 0) then
    write (io_unit%mpi,'(a,2i5,f10.6)') &
      'write_out begin:',snap, rec, wallclock()-start
  end if
  if (verbose > 1) &
    write (io_unit%mpi,*) &
      'write_out: record,recl=', rec, recl
  !-----------------------------------------------------------------------------
  call io%lock%set ('background_io_t%write_out')
  write (unit, rec=rec) buffer
  call io%lock%unset ('background_io_t%write_out')
  !-----------------------------------------------------------------------------
  !call log (rec, 4)
  if (verbose > 0) then
    write (io_unit%mpi,'(a,2i5,f10.6)') &
      'write_out   end:', snap, rec, wallclock()-start
  end if
  call trace%end (itimer)
END SUBROUTINE write_out

!===============================================================================
!> Read a record with a data cube from a direct access file
!===============================================================================
SUBROUTINE input (self, ok)
  class(patch_t):: self
  class(task_t), pointer:: task
  logical:: ok
  !.............................................................................
  character(len=64) filename
  real, pointer:: buffer(:,:,:,:)
  integer:: unit, l(3), u(3), recl
  !-----------------------------------------------------------------------------
  call trace%begin('background_io_t%input')
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/snapshot_",i5.5,".dat")') &
    trim(io%inputdir), io%restart, self%in%rank
  !-----------------------------------------------------------------------------
  ! Use the guard zones value from the input metadata to set the buffer size
  !-----------------------------------------------------------------------------
  l = merge (self%mesh%lb, self%mesh%li, io%in%guard_zones)
  u = merge (self%mesh%ub, self%mesh%ui, io%in%guard_zones)
  buffer => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1:self%mv,self%it,1)
  recl = kind(buffer)*product(shape(buffer))
  open (newunit=unit, file=trim(filename), form='unformatted', &
        recl=recl, access='direct', status='unknown')
  write (io_unit%mpi,*) 'background_io_t%input: recl,rec=', recl, self%record
  flush (io_unit%mpi)
  read (unit, rec=self%record) buffer
  close (unit)
  !-----------------------------------------------------------------------------
  if (io%in%do_generic) then
    call task%allocate_generic()
    self%generic = buffer
    call task%from_generic ()
    call task%deallocate_generic()
  end if 
  !-----------------------------------------------------------------------------
  !call log (self%record, 5)
  ok = .true.
  call trace%end ()
END SUBROUTINE input

END MODULE background_io_mod
