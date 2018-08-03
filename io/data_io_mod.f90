!> =============================================================================
!> Interface from gpatch_mod to a choice of binary data I/O methods, controlled
!> by the io%method text string -- optionally (but deprecated) by io%do_legacy
!> or io%do_direct.
!>
!> Concerning thread safety: Simultaneous calls from many threads are typical,
!> and it is the responsibility of each sub-module here (e.g. parallel_io_mod)
!> to introduce critical regions to protect agains simultaneus calls to non-
!> threadsafe procedures, or if/when accessing modified shared data.
!> =============================================================================
MODULE data_io_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE omp_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE direct_io_mod
  USE legacy_io_mod
  USE parallel_io_mod
  USE buffered_io_mod
  USE snapshot_list_mod
  USE index_mod
  implicit none
  private
  type, public:: data_io_t
    logical:: initialized=.false.
  contains
    procedure:: init
    procedure:: output
    procedure, nopass:: output_txt
    procedure, nopass:: output_nml
    procedure:: input
    procedure, nopass:: close
  end type
  type(data_io_t), public:: data_io
CONTAINS

!> =============================================================================
!> Initialize the I/O method, if not already done
!> =============================================================================
SUBROUTINE init (self, task)
  class(data_io_t):: self
  class(patch_t):: task
  !.............................................................................
  if (self%initialized) return
  !$omp critical (data_io_cr)
  if (.not.self%initialized) then
    call trace%begin('data_io_t%init')
    if (io%do_output) then
      select type (task)
      class is (patch_t)
        if      (io%method=='direct') then
          call direct_io%init (task)
        else if (io%method=='parallel') then
          call parallel_io%init (task)
        else if (io%method=='paralle5') then
          io%format = 5
          call parallel_io%init (task)
        end if
      end select
    end if
    self%initialized = .true.
    call trace%end
  end if
  !$omp end critical (data_io_cr)
END SUBROUTINE init

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self, patch, experiment_name)
  class(data_io_t):: self
  class(patch_t):: patch
  !.............................................................................
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  logical, save:: checked=.false.
  type(snapshot_t), pointer:: snapshot
  integer, save:: itimer=0
  logical:: append
  integer:: it
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin('data_io_t%output', itimer=itimer)
  call self%init (patch)
  if (io%time_derivs) then
    it = patch%iit(patch%nt-2)
  else
    it = patch%iit(patch%nt-1)
  end if
  do while (patch%t(it) >= patch%out_next)
    !$omp critical (outputname_cr)
    if (.not.checked) then
      checked = .true.
      call os%mkdir (trim(io%outputname))
    end if
    !$omp end critical (outputname_cr)
    call os%mkdir_no (patch%iout)
    call snapshot_list%count (patch, snapshot, .false.)
    append = (snapshot%ibuf < snapshot%nbuf-1)
    if      (io%method=='direct') then
      call direct_io%output (patch)
    else if (io%method=='legacy') then
      call legacy_io%output (patch)
      call output_txt (patch, append)
    else if (io%method=='buffered') then
      call buffered_io%output (patch)
    else if (io%method=='paralle5') then
      io%format = 5
      call parallel_io%output (patch)
    else if (io%method=='parallel') then
      io%format = 6
      call parallel_io%output (patch)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
    call output_nml (patch, append)
    if (snapshot%ibuf==0) &
      call snapshot_list%remove (snapshot)
    patch%iout = patch%iout+1
    if (io%out_time==0.0) exit
    patch%out_next = (nint(patch%out_next/io%out_time)+1)*io%out_time
    call os%mkdir_no (patch%iout)
  end do
  call trace%end (itimer)
END SUBROUTINE output

!> =============================================================================
!> Read snapshot from disk
!> =============================================================================
SUBROUTINE input (self, patch, ok)
  class(data_io_t):: self
  class(patch_t):: patch
  logical:: ok
  type(snapshot_t), pointer:: snapshot
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%input')
  call self%init (patch)
  ok = .false.
  if (patch%restart > 0) then
    if      (io%method=='direct') then
      call direct_io%input (patch, ok)
    else if (io%method=='legacy') then
      call legacy_io%input (patch, ok)
    else if (io%method=='parallel') then
      call parallel_io%input (patch, ok)
    else
      call mpi%abort ('unknown io%method = '//io%method)
    end if
  end if
  if (ok) &
    call snapshot_list%count (patch, snapshot, .true.)
  patch%t(patch%it) = patch%time
  if (ok .and. io%verbose>2) print *, patch%id, ' LEGACY INPUT'
  call trace%end
END SUBROUTINE input

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
SUBROUTINE output_txt (patch, append)
  class(patch_t):: patch
  logical:: append
  !.............................................................................
  real(8):: time
  integer:: l(3), u(3)
  character(len=128) filename
  integer:: nlines=6
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_txt')
  !$omp critical (output_cr)
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  if (io%format > 2) then
     write (filename,'(a,i5.5,"/rank_",i5.5,".txt")') trim(io%outputname), patch%iout, mpi%rank
  else
     write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%iout, patch%id
  end if
  if (io%debug(2)) &
    print'(a,3x,"time:",2f13.6)', ' data_io_t%output_txt: '//trim(filename), &
    patch%out_next, patch%time
  !-----------------------------------------------------------------------------
  ! Open the file, appending or not
  !-----------------------------------------------------------------------------
  open  (io_unit%data, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
  write (io_unit%data,'(i2.2)') io%format
  if (io%format==3 .or. io%format==4) write (io_unit%data,*) patch%id, nlines
  !-----------------------------------------------------------------------------
  ! Write the text info
  !-----------------------------------------------------------------------------
  write (io_unit%data,'(a)') trim(patch%kind)
  write (io_unit%data,'(a)') trim(patch%eos)
  write (io_unit%data,'(a)') trim(patch%opacity)
  write (io_unit%data,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') patch%out_next, &
        patch%dtime, patch%position, patch%ds, patch%velocity, patch%llc_nat, &
        patch%llc_cart, patch%mesh%centre_nat, patch%quality, patch%gamma
  write (io_unit%data,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(patch%box/patch%ds+0.5), &
        u-l+1, patch%li, patch%ui, patch%n, patch%gn, patch%ng, patch%nv, patch%level
  write (io_unit%data,'(3(1i8.1,1x),i2.1)') patch%id, patch%iout, patch%istep, patch%mesh_type
  close (io_unit%data)
  !$omp end critical (output_cr)
  call trace%end
END SUBROUTINE output_txt

!> =============================================================================
!> Output text with patch info, either to individual files (io%format=1,2) or
!> to one file per rank (io%format>2).  This is not threadsafe, both because of
!> the text file I/O, and because snapshot%ibuf is consantly being modified by
!> other threads (must be protected by the same critical region!).
!> =============================================================================
SUBROUTINE output_nml (patch, append)
  class(patch_t):: patch
  logical :: append
  !.............................................................................
  character(len=128)   :: filename
  integer              :: ioformat, id, iout, istep, mesh_type, level,nv, nt, &
                          nw, ntotal, format
  integer, dimension(3):: ncell, li, ui, n, ng, gn, l, u
  real                 :: gamma, quality
  real(8)              :: time, dtime, out_next, out_time
  real(8), dimension(3):: size, position, ds, box, velocity, &
                          llc_nat, llc_cart, centre_nat
  logical              :: guard_zones, time_derivs, no_mans_land, periodic(3)
  character(len=16)    :: kind, eos, opacity
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs
  namelist /snapshot_nml/ iout, time, dtime, ntotal
  namelist /patch_nml/ ioformat, iout, time, dtime, id, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, &
    ncell, li, ui, n, ng, gn, nv, nt, nw, ntotal, &
    llc_nat, llc_cart, centre_nat, kind, eos, opacity, &
    periodic, guard_zones, time_derivs, no_mans_land
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  call trace%begin ('data_io_t%output_nml')
  !$omp critical (output_cr)
  !-----------------------------------------------------------------------------
  ! Write snapshot_nml to the file "data/run/NNNNN/snapshot.nml"
  !-----------------------------------------------------------------------------
  format      = io%format
  ioformat    = io%format
  ntotal      = io%ntotal
  out_time    = io%out_time
  guard_zones = io%guard_zones
  time_derivs = io%time_derivs
  iout        = patch%iout
  time        = merge (patch%time, patch%out_next, io%out_time==0)
  dtime       = patch%dtime
  if (mpi%master .and. .not.append) then
    write (filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), patch%iout
    open  (io_unit%data, file=trim(filename), form='formatted', status='unknown')
    write (io_unit%data, io_nml)
    write (io_unit%data, snapshot_nml)
    call patch%idx%output
    call file_append (io_unit%nml, io_unit%data)
    close (io_unit%data)
    if (io%debug(2)) print'(a,3x,"time:",2f13.6)', &
      ' data_io_t%output_nml: '//trim(filename), patch%out_next, patch%time
  end if
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml", appending or not
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
         trim(io%outputname), patch%iout, mpi%rank
  open  (io_unit%data, file=trim(filename), form='formatted', status='unknown', &
         access=trim(merge('append    ','sequential',append)))
  !-----------------------------------------------------------------------------
  ! Write the namelist info, including a leading idx_nml
  !-----------------------------------------------------------------------------
  if (mpi%master .and. .not.append) then
    call patch%idx%output
  end if
  l = merge(patch%mesh%lb,patch%mesh%li,io%guard_zones)
  u = merge(patch%mesh%ub,patch%mesh%ui,io%guard_zones)
  ncell        = u-l+1
  id           = patch%id
  istep        = patch%istep
  mesh_type    = patch%mesh_type
  position     = patch%position
  size         = patch%size
  ds           = patch%ds
  box          = patch%box
  velocity     = patch%velocity
  level        = patch%level
  quality      = patch%quality
  gamma        = patch%gamma
  llc_nat      = patch%llc_nat
  llc_cart     = patch%llc_cart
  centre_nat   = patch%centre_nat
  li           = patch%li
  ui           = patch%ui
  n            = patch%n
  gn           = patch%gn
  ng           = patch%ng
  nv           = patch%nv
  nv           = merge (io%nv, patch%nv, io%nv>0)
  nv           = merge (2*nv, nv, io%time_derivs)
  nt           = patch%nt
  nw           = patch%nw
  kind         = patch%kind
  eos          = patch%eos
  opacity      = patch%opacity
  no_mans_land = patch%no_mans_land
  periodic     = patch%periodic
  write (io_unit%data, patch_nml)
  close (io_unit%data)
  !$omp end critical (output_cr)
  call trace%end
END SUBROUTINE output_nml

!> =============================================================================
!> Copy one text file to another
!> =============================================================================
SUBROUTINE file_append (unit1, unit2)
  integer:: unit1, unit2
  character(len=128):: line
  integer:: iostat
  !-----------------------------------------------------------------------------
  rewind (unit1)
  do while (.true.)
    read (unit1,'(a)',iostat=iostat) line
    if (iostat /= 0) exit
    write (unit2,'(a)') trim(line)
  end do
END SUBROUTINE file_append

!> =============================================================================
!> Close the data I/O
!> =============================================================================
SUBROUTINE close
  !.............................................................................
  if (.not.io%do_output) return
  call trace%begin('patch_t%close')
  if      (io%method=='legacy') then
    !call legacy_io%close
  else if (io%method=='direct') then
    call direct_io%close
  else if (io%method=='buffered') then
    call buffered_io%close
  else if (io%method=='parallel') then
    call parallel_io%close
  end if
  call trace%end
END SUBROUTINE close

END MODULE data_io_mod
