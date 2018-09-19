!> =============================================================================
!> Use MPI parallel I/O to write everythin to a single file.  No critical regions
!> are needed here; all thread protection is one inside mpi_io_file_mod and
!> mpi_io_mod.
!> =============================================================================
MODULE parallel_io_mod
  USE patch_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_io_mod
  USE mpi_io_file_mod
  USE trace_mod
  USE central_mod
  implicit none
  private
  type, public:: parallel_io_t
  contains
    procedure, nopass:: init
    procedure, nopass:: input
    procedure, nopass:: output
    procedure, nopass:: close
  end type
  integer, parameter:: verbose=0
  type(central_t):: central
  type(mpi_io_file_t), save:: file_in, file_out
  type(parallel_io_t), public:: parallel_io
CONTAINS

!> =============================================================================
!> Initialize for MPI paraallel I/O.  This is an MPI collective operation, and
!> needs to be done at the start of the run.
!> =============================================================================
SUBROUTINE init (patch)
  class(patch_t):: patch
  !.............................................................................
  character(len=5):: snapname
  character(len=120):: filename
  logical:: exists
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%init')
  call mpi_io%init
  if (io%restart > 0) then
    filename=trim(io%inputdir)//'snapshots.dat'
    inquire (file=filename, exist=exists)
    if (.not.exists) then
      filename=trim(io%inputdir)//'00000/snapshot.dat'
      inquire (file=filename, exist=exists)
      if (.not.exists) then
        if (io%master) print *, 'neither '//trim(io%inputdir) &
          //'snapshots.dat, nor '//trim(filename)//' exist'
        call mpi%abort()
      end if
    end if
    if (mpi%master) print *,' input file: ',filename
    call file_in%openr (filename)
  end if
  if (io%format==5) then
    filename=trim(io%outputname)//'00000/snapshot.dat'
  else
    filename=trim(io%outputname)//'snapshots.dat'
  end if
  if (mpi%master) print *,'output file: ',filename
  call file_out%openw (filename)
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Close the data I/O method
!> =============================================================================
SUBROUTINE close
  !.............................................................................
  call trace%begin ('parallel_io_t%close')
  call file_out%close
  call file_in%close
  call trace%end
END SUBROUTINE close

!> =============================================================================
!> Write output in increments of one patch data cube
!> =============================================================================
SUBROUTINE output (patch)
  class(patch_t):: patch
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), nv, iw, nw
  real, dimension(:,:,:,:), pointer:: buffer, mem
  real:: pt(2)
  logical:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allow the output to have fewer variables than the solver (e.g. for isothermal)
  !-----------------------------------------------------------------------------
  nv = merge (io%nv, patch%nv, io%nv>0)
  !-----------------------------------------------------------------------------
  ! Set size parameters; a buffer is all variables in the patch, and the
  ! snspshot offset is that times the number of patches per snapshot.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  !-----------------------------------------------------------------------------
  ! A call to mpi_io%iwrite (iout, id) result in slots in the file ordered
  ! by (iv,id,iout), where iv=1,nv; id=1,io%ntotal, iout=1,nsnap. To output
  ! time derivatives, while still using a buffer with nv slots, we pretend that
  ! io%ntotal is twice as large, and output the time derivs as patch 2+2*(id-1),
  ! with variables as patch 1+2*(id-1).
  !-----------------------------------------------------------------------------
  buffer_words = product(n)*nv
  if (io%time_derivs) then
    snapshot_words = buffer_words*io%ntotal*2
  else
    snapshot_words = buffer_words*io%ntotal
  end if
  call mpi_io%use (file_out)
  call mpi_io%set (snapshot_words, buffer_words)
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  if (io%verbose==1) print '(a,2i6)', ' parallel_io_t%output: iout, id =', &
    patch%iout, patch%id
  if (first_time) then
    first_time = .false.
    if (io%master) print '(a,f8.3,a,i10,a)', &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
  end if
  if (io%verbose>1) print '(a,4i6,2i12)', &
    ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  !-----------------------------------------------------------------------------
  ! Interpolate to the given output time
  !-----------------------------------------------------------------------------
  allocate (buffer(n(1),n(2),n(3),patch%nv))
  call patch%time_interval (patch%out_next, jt, pt)
  nw = merge (patch%nw, 1, io%time_derivs)
  do iw=1,nw
    do iv=1,nv
      if (io%guard_zones) then
        buffer(:,:,:,iv) = patch%mem(:,:,:,iv,jt(1),1)*pt(1) &
                         + patch%mem(:,:,:,iv,jt(2),1)*pt(2)
      else
        buffer(:,:,:,iv) = patch%mem(patch%li(1):patch%ui(1), &
                                     patch%li(2):patch%ui(2), &
                                     patch%li(3):patch%ui(3),iv,jt(1),iw)*pt(1) &
                         + patch%mem(patch%li(1):patch%ui(1), &
                                     patch%li(2):patch%ui(2), &
                                     patch%li(3):patch%ui(3),iv,jt(2),iw)*pt(2)
      end if
      if (io%verbose>2 .or. (io%verbose>1 .and. iv==patch%idx%d)) &
        write (io_unit%output,'(2a,i7,i4,1p,e14.6,i4,2e12.4)') &
         ' parallel_io_t%output: ', &
         ' id,iout,time,iv,vmin,vmax =', patch%id, patch%iout, patch%time, iv, &
          minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv))
      if (verbose>0) &
        write (io_unit%mpi,'(2a,i7,i4,1p,e14.6,i3,6e12.4,2i3,0p,2f7.3)') &
         ' parallel_io_t%output: ', &
         'id,iout,time,iv,vmin,vmax =', patch%id, patch%iout, patch%time, iv, &
          minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv)), &
          patch%fminval(patch%mem(:,:,:,iv,jt(1),1)), &
          patch%fmaxval(patch%mem(:,:,:,iv,jt(1),1)), &
          patch%fminval(patch%mem(:,:,:,iv,jt(2),1)), &
          patch%fmaxval(patch%mem(:,:,:,iv,jt(2),1)), &
          jt, pt
    end do
    !---------------------------------------------------------------------------
    ! For RAMSES consistency, convert momenta to velocities for ramses solver
    !---------------------------------------------------------------------------
    if (patch%solver_is('ramses')) then
      do iv=patch%idx%px,patch%idx%pz
        buffer(:,:,:,iv) = buffer(:,:,:,iv)/buffer(:,:,:,patch%idx%d)
      end do
    end if
    if (verbose > 0) then
      pos = 4_8*(patch%iout*snapshot_words + (patch%id-1_8)*buffer_words)
      write (io_unit%mpi,*) 'WR6', patch%id, patch%iout, pos
    end if
    !---------------------------------------------------------------------------
    ! If native time derivatives are available (nw=2), output these to buffer,
    ! using offsets that effectively double the number of variables per snapshot
    !---------------------------------------------------------------------------
    if (io%time_derivs) then
      if (patch%nw==1) then
        call mpi_io%iwrite (buffer, 1+patch%iout, 1+2*(patch%id-1))
      else
        call mpi_io%iwrite (buffer, 1+patch%iout, 2+2*(patch%id-1))
      end if
    !---------------------------------------------------------------------------
    ! If no time derivs, increment with one per snapshot
    !---------------------------------------------------------------------------
    else
      call mpi_io%iwrite (buffer, 1+patch%iout, patch%id)
    end if
  end do
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivative output
  !-----------------------------------------------------------------------------
  if (io%time_derivs .and. patch%nw==1) then
    do iv=1,nv
      if (patch%solver_is('ramses')) then
        allocate (mem(patch%gn(1),patch%gn(2),patch%gn(3),patch%nt))
        if ((iv >= patch%idx%px) .and. (iv <= patch%idx%pz)) then
          mem = patch%mem(:,:,:,iv,:,1)/patch%mem(:,:,:,patch%idx%d,:,1)
        else
          mem = patch%mem(:,:,:,iv,:,1)
        end if
        call central%time_derivative (patch, mem, buffer(:,:,:,iv))
        deallocate (mem)
      else
        call central%time_derivative (patch, patch%mem(:,:,:,iv,:,1), buffer(:,:,:,iv))
      end if
    end do
    call mpi_io%iwrite (buffer, 1+patch%iout, 2+2*(patch%id-1))
  end if
  call trace%end (itimer)
END SUBROUTINE output

!> =============================================================================
!> Read input in increments of one variable data cube
!> =============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), nv
  real, dimension(:,:,:,:), pointer:: buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Set size parameters
  !-----------------------------------------------------------------------------
  if (file_in%handle==0) then
    call init (patch)
  end if
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  !-----------------------------------------------------------------------------
  ! Allow the data to have fewer variables than the solver.  Note that this
  ! requires using namelist parameters when reading that are consistent with the
  ! parameters used when writing.
  !-----------------------------------------------------------------------------
  nv = merge (io%nv, patch%nv, io%nv>0)
  !-----------------------------------------------------------------------------
  ! Compute offset size constants for the mpi_io reader
  !-----------------------------------------------------------------------------
  buffer_words = product(n)*nv
  if (io%time_derivs) then
    snapshot_words = buffer_words*io%ntotal*2
  else
    snapshot_words = buffer_words*io%ntotal
  end if
  call mpi_io%use (file_in)
  call mpi_io%set (snapshot_words, buffer_words)
  !-----------------------------------------------------------------------------
  ! Read via buffer, in case guard zones not included. FIXME: read snapshot
  ! guard zone parameter, to allow it to differ
  !-----------------------------------------------------------------------------
  allocate (buffer(n(1),n(2),n(3),nv))
  if (io%time_derivs) then
    call mpi_io%read (buffer, 1+io%restart, 1 + 2*(patch%id-1))
  else
    call mpi_io%read (buffer, 1+io%restart, 1 +   (patch%id-1))
  end if
  if (mpi_io%err /= 0) then
    ok = .false.
    call trace%end
    return
  end if
  if (io%guard_zones) then
    patch%mem(:,:,:,:,patch%it,1) = buffer
  else
    patch%mem(patch%li(1):patch%ui(1), &
              patch%li(2):patch%ui(2), &
              patch%li(3):patch%ui(3),:,patch%it,1) = buffer
  end if
  if (io%verbose>1) then
    do iv=1,patch%nv
      write (io_unit%output,*) 'parallel_io_t%input: id, iv, iout, min, max =', &
        patch%id, iv, patch%iout, patch%fminval(iv), patch%fmaxval(iv)
    end do
  end if
  deallocate (buffer)
  !-----------------------------------------------------------------------------
  ! Restore to momentum variables for RAMSES snapshots
  !-----------------------------------------------------------------------------
  if (patch%solver_is('ramses')) then
    do iv=patch%idx%px,patch%idx%pz
      patch%mem(:,:,:,iv         ,patch%it,1) = &
      patch%mem(:,:,:,iv         ,patch%it,1) * &
      patch%mem(:,:,:,patch%idx%d,patch%it,1)
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Signal success
  !-----------------------------------------------------------------------------
  ok = .true.
  patch%time = io%restart*io%out_time
  if (io%verbose==1) print '(a,2i6,1p,e16.6)', ' parallel_io_t%input: iout, id, time =', &
    patch%iout, patch%id, patch%time
  if (io%verbose>1) print '(a,4i6,2i7)', &
    ' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call trace%end (itimer)
END SUBROUTINE input

END MODULE parallel_io_mod
