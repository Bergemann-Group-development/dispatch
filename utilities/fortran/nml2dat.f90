  !===============================================================================
!> Utility to convert .nml metadata files to unformatted Fortran files
!===============================================================================
PROGRAM nm2dat
  USE omp_lib
  implicit none
  logical:: last
  character(len=80):: base
  integer:: id, level, task_offset, iostat, n, io, rank, key=1
  real(8):: time, dtime, times(2), dtimes(2), position(3), size(3), start
  namelist /patch_nml/id,time,position,size,level,dtime,times,dtimes,task_offset
  !-----------------------------------------------------------------------------
  ! Loop over snapshots -- terminate when fails to find rank 0 file
  !-----------------------------------------------------------------------------
  start=omp_get_wtime()
  do io=1,99
    !---------------------------------------------------------------------------
    ! Loop over ranks -- terminate when fails to find rank file
    !---------------------------------------------------------------------------
    last = .false.
    do rank=0,99
      write (base,'(a,i5.5,a,i5.5)') 'data/512/',io,'/rank_',rank
      open (1,file=trim(base)//'_patches.nml', form='formatted', &
        status='old', iostat=iostat)
      if (iostat /= 0) then
        if (rank==0) then
          last = .true.
        else
          print '(2(1x,a,i4,2x))', 'snapshot:',io,'    last rank:',rank-1
        end if
        exit
      end if
      open (2,file=trim(base)//'_patches.dat', form='unformatted', &
        status='unknown')
      write (2) key
      n=0
      iostat=0
      !-------------------------------------------------------------------------
      ! Loop over patches -- terminate at EOF
      !-------------------------------------------------------------------------
      do while (iostat==0)
        read (1,patch_nml,iostat=iostat)
        write(2) id,time,position,size,level,dtime,times,dtimes,task_offset
        n=n+1
      end do
      close(1)
      close(2)
    end do
    if (last) then
      exit
    end if
  end do
  print '(f6.2," sec")',omp_get_wtime()-start
END PROGRAM
