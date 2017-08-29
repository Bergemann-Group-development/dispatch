subroutine backup_poisson(filename)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,icpu,mcpu,info
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_poisson'

  !-----------------------------------
  ! Output grav vars in file. Max 128 threads output simultaneously.
  !-----------------------------------  
  mcpu = (ncpu_dump-1)/128+1                       ! number of times we have to dump
  do icpu=0,mcpu-1                            ! Nr of times we dump
     if (myid==1 .and. mcpu > 3) print '(a,i3,a,i3)', 'Dumping GRAV batch ',icpu+1,' out of ',mcpu
     if (modulo(myid-1,mcpu)==icpu .and. ngrid_current>0) then
        ilun=10
           
        call title(myid,nchar)
        fileloc=TRIM(filename)//TRIM(nchar)
        open(unit=ilun,file=fileloc,form='unformatted')
        write(ilun)ncpu_dump
        write(ilun)ndim+1
        write(ilun)nlevelmax
        write(ilun)nboundary
        do ilevel=1,nlevelmax
           do ibound=1,nboundary+ncpu_dump
              if(ibound<=ncpu_dump)then
                 ncache=numbl(ibound,ilevel)
                 istart=headl(ibound,ilevel)
              else
                 ncache=numbb(ibound-ncpu_dump,ilevel)
                 istart=headb(ibound-ncpu_dump,ilevel)
              end if
              write(ilun)ilevel
              write(ilun)ncache
              if(ncache>0)then
                 allocate(xdp(1:ncache))
                 ! Loop over cells
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    do ivar=1,ndim
                       igrid=istart
                       do i=1,ncache
                          xdp(i)=f(igrid+iskip,ivar)
                          igrid=next(igrid)
                       end do
                       write(ilun)xdp
                    end do
                    igrid=istart
                    do i=1,ncache
                       xdp(i)=phi_old(igrid+iskip)
                       igrid=next(igrid)
                    end do
                    write(ilun)xdp
                 end do
                 deallocate(xdp)
              end if
           end do
        end do
        close(ilun)
     end if
#ifndef WITHOUTMPI
     call MPI_BARRIER(mpi_comm_use,info)
#endif
  end do
     
end subroutine backup_poisson
