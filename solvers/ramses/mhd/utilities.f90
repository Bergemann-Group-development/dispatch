! Detection of NaN
!=======================================================================
LOGICAL FUNCTION is_nan(x)
  use amr_parameters, only : dp
#if !defined(__GFORTRAN__)
  use, intrinsic :: ieee_arithmetic
#endif
  implicit none
  real(dp) :: x

#if defined(__GFORTRAN__)
  is_nan = isnan(x)
#else
  is_nan = ieee_is_nan(x)
  is_nan = is_nan .or. abs(x) > 1d100
#endif
END FUNCTION is_nan

! Routine to get the hostname of the process - useful for debug and thread scheduling
!=======================================================================
SUBROUTINE gethostname(name)
  use amr_commons, only : myid
#ifdef __INTEL_COMPILER
  use ifport, only : MAX_HOSTNAM_LENGTH, HOSTNAM
#endif
#ifdef __xlc__
  use xlfutility, only : hostnm_
#endif
  implicit none
  character(len=*)  :: name
  logical, save     :: first_call=.true.
  integer           :: lname, istat, lmax
  integer(kind=1)   :: istat1
  logical           :: compiler_detected
#ifdef __INTEL_COMPILER
  character(MAX_HOSTNAM_LENGTH + 1) :: hostname
#endif
#ifdef __PGI
  integer, external :: hostnm
#endif
#ifdef __xlc__
  character(len=63) :: hostname
#endif
  !
  compiler_detected=.false.
  lname = len(name)
#ifdef __INTEL_COMPILER
  istat = HOSTNAM(hostname)
  lmax  = max(lname,len_trim(hostname))
  name  = ''
  name(1:lmax)  = hostname(1:lmax)
  compiler_detected=.true.
#endif
#ifdef __GFORTRAN__
  istat1 = hostnm(name)
  compiler_detected=.true.
#endif
#ifdef __PGI
  istat  = hostnm(name)
  compiler_detected=.true.
#endif
#ifdef __xlc__
  istat = hostnm_(hostname)
  lmax  = max(lname,len_trim(hostname))
  name  = ''
  name(1:lmax)  = hostname(1:lmax)
  compiler_detected=.true.
#endif
  if (.not. compiler_detected .and. first_call .and. myid==1) then
     print *, 'WARNING from gethostname. Compiler does not support hostname function. Please add compiler specific call'
     name = ''
  end if
  first_call = .false.
END SUBROUTINE gethostname
!=======================================================================
SUBROUTINE mkdir (dir)
  USE amr_commons
  implicit none
  character(len=*) dir
  character(len=80) cmd
  logical exists
  integer ntry, mpi_err
!#ifndef WITHOUTMPI
!  include 'mpif.h'
!#endif  
!-----------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  inquire(directory=trim(dir),exist=exists)
#else
  inquire(file=trim(dir),exist=exists)
#endif
  ntry=0
  do while (.not. exists)
    cmd='mkdir -p '//TRIM(dir)
#ifdef NOSYSTEM
    if(myid==1)call PXFMKDIR(TRIM(dir),LEN(TRIM(dir)),O'755',info)
#else
    if(myid==1)call system(cmd)
#endif
#ifdef __INTEL_COMPILER
    inquire(directory=trim(dir),exist=exists)
#else
    inquire(file=trim(dir),exist=exists)
#endif
    ntry=ntry+1
    if (ntry > 20) then
      if(myid==1) print*,'failed to make directory '//trim(dir)
      dir='./'
      if(myid==1) print*,'using directory '//trim(dir)
      exists=.true.
    endif
    if (.not. exists) call sleep(1)
  enddo
#ifndef WITHOUTMPI
  call MPI_BCAST(dir,len(dir),MPI_CHARACTER,0,mpi_comm_use,mpi_err)
#endif
END
!=======================================================================
INTEGER(KIND=8) FUNCTION file_position(u)
#if __INTEL_COMPILER
  USE IFPORT, only : ftelli8
#endif
  USE amr_commons, only : myid
  implicit none
#ifdef __PGI
#ifdef __i386__
  integer, external   :: ftell
#else
  integer, external   :: ftell64
#endif
#endif
#ifdef __xlc__
  integer(kind=8), external   :: ftell64_
#endif
  integer, intent(in) :: u
  logical, save       :: first_call
#if __INTEL_COMPILER
  file_position = ftelli8(u)
#endif
#if defined(__GFORTRAN__) || defined(__G95__)
  if (first_call .and. myid==1) print *, 'file_position','You have compiled with g95'// &
                                         ' or gfortran, and they only contain a 32 bit'// &
                                         ' version of ftell, the underlying C-function.'// &
                                         ' file_position could overflow for files > 2GB.'
  file_position = ftell(u)
  first_call =.false.
#endif
#ifdef __xlc__
  file_position = ftell64_(u)
#endif
#ifdef __PGI
#ifdef __i386__
  file_position = ftell(u)
#else
  file_position = ftell64(u)
#endif
#endif

END FUNCTION file_position
!=======================================================================
SUBROUTINE trace_mpi(routine,label,level)
  USE amr_commons
  implicit none
  character(len=*) routine,label
  integer level
  character(len=120) file,dir
  character(len=16) fmt
  real(kind=8), save:: t0=0., t1
  logical, save:: first_time=.true., exists=.false.
  integer, save:: n=0, indent=0
!#ifndef WITHOUTMPI
!  include 'mpif.h'
!#endif  
!-----------------------------------------------------------------------
  if (level > trace_level) return

  !$omp master
  if (t0==0.0) t0=wallclock()
  t1=wallclock()
  n=n+1

  write(file,'(a,i5.5,a)') trim(datadir)//'/trace_',myid,'.log'
  if (trim(label)=='init') then
    open(2,file=file,form='formatted',status='unknown')
  else
    open(2,file=file,form='formatted',status='unknown',position='append')
  endif 

  if (trim(label)=='exit' ) indent=indent-1
  write(fmt,'(a,i2,a)') '(i7,f10.3,',3+2*indent,'x,a)'
  write(2,fmt) n,t1-t0,trim(routine)//': '//trim(label)
  if (trim(label)=='enter') indent=indent+1
  close(2)

  t0=wallclock()
  !$omp end master
  end
!=======================================================================
MODULE timer_m
  implicit none
  integer,            parameter         :: mtimer=100                    ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE findit (label)
  implicit none
  character(len=*) label
  do itimer=1,ntimer
    if (trim(label) == trim(labels(itimer))) return
  end do
  ntimer = ntimer+1
  itimer = ntimer
  labels(itimer) = label
  time(itimer) = 0.
END SUBROUTINE
END MODULE
!=======================================================================
SUBROUTINE timer (label, cmd)
  USE amr_commons
  USE timer_m 
  implicit none
  character(len=*) label, cmd
  real(kind=8) current 
  integer ierr
!-----------------------------------------------------------------------
  if (a_verbose>5) print*,'timer:',myid,TRIM(label)//':'//TRIM(cmd)
#ifndef WITHOUTMPI
  if (a_verbose>6) call MPI_BARRIER(mpi_comm_use,ierr)
#endif
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
    time(itimer) = time(itimer) + current - start(itimer)               ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
    start(itimer) = current                                             ! register start time
  else if (cmd == 'stop') then                                          ! stop command
    itimer = 0                                                          ! turn off timer
  end if
END SUBROUTINE
!=======================================================================
SUBROUTINE finalize_timer
  USE amr_parameters
  USE amr_commons
  USE timer_m
  implicit none
!#ifndef WITHOUTMPI
!  include 'mpif.h'
!#endif  
  real(kind=8) :: total, gtotal, gtime
  integer      :: mpi_err
  logical      :: o
!-----------------------------------------------------------------------
  o = myid == 1
  total = 1e-9
  if (o) write (*,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do itimer = 1,ntimer
    total = total + time(itimer)
  end do
  do itimer = 1,ntimer
    if (o .and. time(itimer)/total > 0.001) write (*,'(f12.3,4x,f6.1,4x,a24)') time(itimer), 100.*time(itimer)/total,labels(itimer)
  end do
  if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  if (ncpu > 1) then

#ifndef WITHOUTMPI
    call MPI_BARRIER(mpi_comm_use,mpi_err)
    call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,mpi_comm_use,mpi_err)
    gtotal = gtotal / ncpu
#else
    gtotal = total
#endif
    if (o) write (*,*) hl
    if (o) write (*,'(/a)') '     seconds         %    STEP'
    do itimer = 1,ntimer
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(real(time(itimer),kind=8),gtime,1,MPI_REAL8,MPI_SUM,mpi_comm_use,mpi_err)
      gtime = gtime / ncpu
#else
      gtime = time(itimer)
#endif
      if (o .and. gtime/gtotal > 0.001) write (*,'(f12.3,4x,f6.1,4x,a24)') gtime, 100.*gtime/gtotal,labels(itimer)
    end do
    if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  endif
END SUBROUTINE
!=======================================================================
SUBROUTINE local_timer
  USE amr_parameters
  USE amr_commons
  USE timer_m
  implicit none
!#ifndef WITHOUTMPI
!  include 'mpif.h'
!#endif  
  real(kind=8) :: total, gtotal, gtime, muspt
  integer      :: mpi_err, n_cell_total, igod
  integer, save:: nstep_prv=0
  character(len=80):: filename, outdir
!-----------------------------------------------------------------------
  if (timer_type==0) return
  write(outdir,'(a,i5.5)') TRIM(datadir)//'/output_',ifout
  call mkdir(outdir)
  write(filename,'(a,i5.5,a)') TRIM(outdir)//'/timer_',myid,'.log'
  if (timer_type == 1) then
     open (tmp_unit,file=filename,status='unknown',form='formatted',position='append')
  else
     open (tmp_unit,file=filename,status='unknown',form='formatted',position='rewind')
  end if
  total = 1e-9
  write (tmp_unit,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do itimer = 1,ntimer
    total = total + time(itimer)
  end do
  igod = 0
  do itimer = 1,ntimer
    if (time(itimer)/total > 0.001) &
      write (tmp_unit,'(f12.3,4x,f6.1,4x,a24)') &
      time(itimer), 100.*time(itimer)/total,labels(itimer)
    if (TRIM(labels(itimer)) == 'hydro - godunov') igod = itimer
  end do
  write (tmp_unit,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  n_cell_total = n_cell
  call sum_int_mpi (n_cell_total)
  write (tmp_unit,'(a, 2i9, f7.3)') 'load balance: cells, total, factor =', &
    n_cell, n_cell_total, n_cell*ncpu/real(n_cell_total)
  if (igod > 0) then
     muspt = time(igod)*1e6/(n_cell*(nstep_coarse-nstep_prv)*twotondim*ncpu)*ncore
     write (tmp_unit,'(a,f10.2,2i5)') 'mus/pt in godunov =', muspt, nstep_coarse-nstep_prv
  endif
  time(1:ntimer) = 0.
  close (tmp_unit)
  nstep_prv = nstep_coarse
END SUBROUTINE
!=======================================================================
SUBROUTINE reset_timer
  USE amr_parameters
  USE amr_commons
  USE timer_m
  implicit none
!#ifndef WITHOUTMPI
!  include 'mpif.h'
!#endif  
  real(kind=8) :: total, gtotal, gtime
  integer      :: mpi_err
  logical      :: o
!-----------------------------------------------------------------------
  do itimer = 1,ntimer
    time(itimer)=0.0
  end do
END SUBROUTINE

!=======================================================================
SUBROUTINE bcast_int_mpi (i)
  use amr_commons
  implicit none
  integer i
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer mpi_err
  call MPI_BCAST(i,1,MPI_INTEGER,0,mpi_comm_use,mpi_err)
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE bcast_ints_mpi (i,n)
  use amr_commons
  implicit none
  integer i,n
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer mpi_err
  call MPI_BCAST(i,n,MPI_INTEGER,0,mpi_comm_use,mpi_err)
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE bcast_real_mpi (r,n)
  use amr_commons
  implicit none
  integer n
  real r(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer mpi_err
  call MPI_BCAST(r,n,MPI_REAL,0,mpi_comm_use,mpi_err)
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE bcast_realdp_mpi (r,n)
  use amr_commons
  implicit none
  integer n
  real(kind=dp) r(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer mpi_err
  call MPI_BCAST(r,n,DPMPI_REAL,0,mpi_comm_use,mpi_err)
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE min_int_mpi (i)
  use amr_commons
  implicit none
  integer i
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer gi, mpi_err
  call MPI_ALLREDUCE(i,gi,1,MPI_INTEGER,MPI_MIN,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE max_int_mpi (i)
  use amr_commons
  implicit none
  integer i
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer gi, mpi_err
  call MPI_ALLREDUCE(i,gi,1,MPI_INTEGER,MPI_MAX,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE min_real_mpi (i)
  use amr_commons
  implicit none
  real(kind=dp) i
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=dp) gi
  integer mpi_err
  call MPI_ALLREDUCE(i,gi,1,DPMPI_REAL,MPI_MIN,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE max_real_mpi (r)
  use amr_commons
  implicit none
  real(kind=dp) r
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=dp) gr
  integer mpi_err
  call MPI_ALLREDUCE(r,gr,1,DPMPI_REAL,MPI_MAX,mpi_comm_use,mpi_err)
  r = gr
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE max_reals_mpi (r,n)
  use amr_commons
  implicit none
  integer n
  real(kind=dp) r(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=dp) gr(n)
  integer mpi_err
  call MPI_ALLREDUCE(r,gr,n,DPMPI_REAL,MPI_MAX,mpi_comm_use,mpi_err)
  r = gr
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE max_real8s_mpi (r,n)
  use amr_commons
  implicit none
  integer n
  real(kind=8) r(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=8) gr(n)
  integer mpi_err
  call MPI_ALLREDUCE(r,gr,n,MPI_REAL8,MPI_MAX,mpi_comm_use,mpi_err)
  r = gr
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE sum_int_mpi (i)
  use amr_commons
  implicit none
  integer i
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer gi, mpi_err
  call MPI_ALLREDUCE(i,gi,1,MPI_INTEGER,MPI_SUM,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE sum_int8_mpi (i)
  use amr_commons
  implicit none
  integer(kind=8) i, gi
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer mpi_err
  call MPI_ALLREDUCE(i,gi,1,MPI_INTEGER8,MPI_SUM,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE sum_ints_mpi (i,n)
  use amr_commons
  implicit none
  integer n,i(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer gi(n), mpi_err
  call MPI_ALLREDUCE(i,gi,n,MPI_INTEGER,MPI_SUM,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE sum_realdp_mpi (x)
  use amr_commons
  implicit none
  integer n
  real(kind=dp) x
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=dp) gx
  integer mpi_err
  call MPI_ALLREDUCE(x,gx,n,DPMPI_REAL,MPI_SUM,mpi_comm_use,mpi_err)
  x = gx
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE sum_real8s_mpi (x,n)
  use amr_commons
  implicit none
  integer n
  real(kind=8) x(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  real(kind=8) gx(n)
  integer mpi_err
  call MPI_ALLREDUCE(x,gx,n,MPI_REAL8,MPI_SUM,mpi_comm_use,mpi_err)
  x = gx
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE max_ints_mpi (i,n)
  use amr_commons
  implicit none
  integer n,i(n)
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer gi(n), mpi_err
  call MPI_ALLREDUCE(i,gi,n,MPI_INTEGER,MPI_MAX,mpi_comm_use,mpi_err)
  i = gi
#endif  
END SUBROUTINE

!=======================================================================
SUBROUTINE abort_mpi (msg)
  use amr_commons, only : myid, mpi_comm_use 
  implicit none
  character(len=*) msg
#ifndef WITHOUTMPI
  !include 'mpif.h'
  integer code, mpi_err
  print*,'abort_mpi: rank, msg =', myid, msg
  code = 1
  call MPI_ABORT(mpi_comm_use, code, mpi_err)
#else
  print*,'abort_mpi: msg =', msg
#endif  
  stop
END SUBROUTINE

!=======================================================================
SUBROUTINE read_real (unit, buf, n)
  use amr_parameters, only: real_size, dp
  implicit none
  integer unit, n
  real(dp) buf(n)
  real(kind=4), allocatable:: buf4(:)
  real(kind=8), allocatable:: buf8(:)
  if (real_size == 4) then
    allocate (buf4(n))
    read (unit) buf4(1:n)
    buf(1:n) = buf4(1:n)
    deallocate (buf4)
  else if (real_size == 8) then
    allocate (buf8(n))
    read (unit) buf8(1:n)
    buf(1:n) = buf8(1:n)
    deallocate (buf8)
  else
    read(unit) buf(1:n)
  endif
END SUBROUTINE

!=======================================================================
SUBROUTINE bitwise(buffer,bytes,label)
  USE amr_parameters, only: bitwise_unit,datadir
  USE amr_commons, only: myid,bitwise_test
  implicit none
  integer bytes
  character buffer(bytes)
#ifdef __GFORTRAN__
  character(len=*) label
#else
  character(len=8) label
#endif
  integer ncompare,nlab,ndiff,i
  character,allocatable, dimension(:):: compare
  character, allocatable, dimension(:):: lab
  character(len=4) number
  character(len=80) filename
  logical:: neq,first_time=.true.,unit_open=.false.
  integer, save:: nread, nwrite, ncomp, ndiffer, nagree

  if (bitwise_test <= 0) then
    if (unit_open) then
      close(bitwise_unit)
      unit_open=.false.
      if (myid==1) then
        if (nread  > 0) print'(a,4i8)', &
          'bitwise test: buffers read, compared, differ, agree =', &
          nread,ncomp,ndiffer,nagree
        if (nwrite > 0) print*,'bitwise test: buffers written =',nwrite
      endif
    endif
    return
  endif

  if (first_time) then
    first_time=.false.
    write(number,'(i4.4)') myid
    filename=TRIM(datadir)//'/bitwise.'//number
    open(unit=bitwise_unit,file=filename,status='unknown',form='unformatted')
    unit_open=.true.
    nread=0
    nwrite=0
    ncomp=0
    ndiffer=0
    nagree=0
    if (myid==1) then
      if (bitwise_test == 1) print*,'bitwise test: started writing'
      if (bitwise_test == 2) print*,'bitwise test: started reading'
    endif
  endif

  if (bitwise_test==1) then
    nlab=sizeof(label)
    write(bitwise_unit) nlab,bytes
    write(bitwise_unit) label
    write(bitwise_unit) buffer
    nwrite=nwrite+1
  else
    nread=nread+1
    read (bitwise_unit) nlab,ncompare

    allocate (lab(nlab), compare(ncompare))
    read(bitwise_unit) lab
    neq=.false.
    do i=1,nlab
      if (lab(i) /= label(i:i)) neq=.true.
      ncomp=ncomp+1
    end do
    if (neq) print*,'bitwise: labels differ for myid, nread, labels =',myid,nread,label,lab

    read(bitwise_unit) compare
    ndiff=ncompare
    do i=1,ncompare
      if (buffer(i)==compare(i)) ndiff=ndiff-1
      ncomp=ncomp+1
    end do
    if (ndiff > 0) then
      ndiffer=ndiffer+1
      print'(a,a8,i6,i12,2i8)','bitwise: buffers differ for label, myid, nread, ndiff = ',label,myid,nread,ncompare,ndiff
    else
      nagree=nagree+1
      print'(a,a8,i6,i12,i8)','bitwise: buffers agree for label, myid, nread, ncompare = ',label,myid,nread,ncompare
    end if

    deallocate (lab,compare)
  end if 
END

!=======================================================================
LOGICAL FUNCTION in_debug_region (ind_grid)
  USE amr_parameters
  USE amr_commons
  implicit none
  integer ind_grid, idim
  real(kind=8) scale, r2
!-----------------------------------------------------------------------
  if (local_debug) then
     in_debug_region=.true.
     scale=boxlen/dble(icoarse_max-icoarse_min+1)
     do idim=1,ndim
        in_debug_region = in_debug_region .and. &
          (abs(xg(ind_grid,idim)*scale - r_debug(idim)) < dr_debug(idim))
     end do
  else
     in_debug_region=.false.
  end iF
END

subroutine barrier_mpi
  use amr_parameters
  use amr_commons
  implicit none
  integer info
#ifndef WITHOUTMPI
  !include 'mpif.h'
  call MPI_BARRIER(mpi_comm_use,info)
#endif  
end

!=======================================================================
SUBROUTINE dump_debug_cube (routine,label,cube,nx,ny,nz,ngrid,nvar,ndim,ind_grid)
  USE amr_parameters, only: datadir, dp, tmp_unit
  USE amr_commons, only: myid, r_debug, dr_debug, xg
  implicit none
  integer,parameter:: fmt=0
  character(len=*):: routine,label
  integer:: nx,ny,nz,ngrid,nvar,ndim,ind_grid(ngrid)
  integer:: i,idim
  real(KIND=dp):: cube(ngrid,nx,ny,nz,nvar)
  character(len=5):: rank
  logical in_debug_region
!-----------------------------------------------------------------------
  do i=1,ngrid
     if (in_debug_region(ind_grid(i))) then
        print*,'MKX',myid,ind_grid(i)
        !$omp critical
        write (rank,'(i5.5)') myid
        open (tmp_unit,file=trim(datadir)//'/'//'debug_'//rank//'.dat', &
           form='unformatted', status='unknown',position='append')
        write (tmp_unit) fmt,myid
        write (tmp_unit) nx,ny,nz,nvar,ndim,len(routine),len(label)
        write (tmp_unit) real(cube(i,:,:,:,:),kind=4),real(xg(ind_grid(i),:),kind=8),routine,label
        close (tmp_unit)
        !$omp end critical
     end if
  end do
END SUBROUTINE
