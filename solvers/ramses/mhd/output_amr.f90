!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_all(advance)
  use amr_commons
  use pm_commons
  use hydro_commons
  !use cooling_module
  implicit none
  logical advance
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  character(LEN=5)::nchar
  character(LEN=80)::filename
  character(LEN=120)::filedir,filecmd
  integer::i,itest,info,ntry,ierr
  logical::exists
  real(kind=8)::wallclock,iotime

  if (ncpu_dump == ncpu .and. .not. do_particle_dump) then
    if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
    if(nstep_coarse==0.and.nrestart>0)return
  endif
  if(verbose)write(*,*)'Entering dump_all'

  if (do_particle_dump) ifout=nrestart                                  ! If updating tarcer particle properties, make sure we dump to the right snapshot

  call write_screen
  call title(ifout,nchar)
  output_done=.true.

  if (advance) then                                                     ! advance to next dump number?
     ifout=ifout+1                                                      ! next dump number
     if(myid==1) print'(1x,a,2i7,1p,2e16.7)', &
       'next snapshot: iout,ifout,tout,aout = ',iout,ifout,tout_iout,aout_iout
  endif

  if(ndim>1)then
     iotime=wallclock()
     filedir=TRIM(datadir)//'/'//'output_'//TRIM(nchar)//'/'
     filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
     if(myid==1)call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
     if(myid==1)call system(filecmd)
#endif
#ifndef WITHOUTMPI
     call MPI_BARRIER(mpi_comm_use,info)
#endif
     exists=.false.
     ntry=0
     do while (.not. exists)
        ntry=ntry+1
#ifdef __INTEL_COMPILER
        inquire(directory=filedir,exist=exists)
#else
        inquire(file=filedir,exist=exists)
#endif
        if (ntry > 20) then
           if(myid==1) print*,'failed to make directory '//trim(filedir)
           filedir='./'
           if(myid==1) print*,'using directory '//trim(filedir)
           exists=.true.
        endif
#ifndef __xlc__
        if (.not. exists) call sleep(1)
#endif
     enddo
#ifndef WITHOUTMPI
     call MPI_BCAST(filedir,80,MPI_CHARACTER,0,mpi_comm_use,ierr)
#endif

     if (do_particle_dump .and. pic)then
        filename=TRIM(filedir)//'part_'//TRIM(nchar)//'.out'
        call backup_part(filename)
        return
     end if

     filename=TRIM(filedir)//'header_'//TRIM(nchar)//'.txt'
     call output_header(filename)
     if(myid==1)then
        filename=TRIM(filedir)//'info_'//TRIM(nchar)//'.txt'
        call output_info(filename)
     endif
     filename=TRIM(filedir)//'amr_'//TRIM(nchar)//'.out'
     call backup_amr(filename)
     if(hydro)then
        filename=TRIM(filedir)//'hydro_'//TRIM(nchar)//'.out'
        call backup_hydro(filename)
     end if
     if(pic)then
        filename=TRIM(filedir)//'part_'//TRIM(nchar)//'.out'
        call backup_part(filename)
     end if
     if(poisson)then
        filename=TRIM(filedir)//'grav_'//TRIM(nchar)//'.out'
        call backup_poisson(filename)
     end if
#ifdef ATON
     if(aton)then
        filename=TRIM(filedir)//'rad_'//TRIM(nchar)//'.out'
        call backup_radiation(filename)
        filename=TRIM(filedir)//'radgpu_'//TRIM(nchar)//'.out'
        call store_radiation(filename)
     end if
#endif
     iotime=wallclock()-iotime
     if(myid==1)print'(1x,a,1p,g17.9,i8,0p,f7.1,a)', &
       'dump_all: '//trim(filedir)//'   t =',t8,iout,iotime,' sec'
  end if

end subroutine dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine backup_amr(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  

  integer::nx_loc,ny_loc,nz_loc,ilun
  integer::ilevel,ibound,ncache,istart,i,igrid,idim,ind,iskip,icpu,mcpu,info
  integer,allocatable,dimension(:)::iig
  real(dp),allocatable,dimension(:)::xdp
  real(dp),dimension(1:3)::skip_loc
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale
  integer ngrid_written

  if(verbose)write(*,*)'Entering backup_amr'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !-----------------------------------
  ! Output amr grid in file. Max 128 threads output simultaneously.
  !-----------------------------------  
  ndomain = ncpu_dump*overload
  mcpu = (ncpu_dump-1)/128+1                       ! number of times we have to dump
  do icpu=0,mcpu-1                            ! Nr of times we dump
     if (myid==1 .and. mcpu > 3) print '(a,i3,a,i3)', 'Dumping AMR batch ',icpu+1,' out of ',mcpu
     if (modulo(myid-1,mcpu)==icpu .and. ngrid_current>1) then
        ilun=10
        call title(myid,nchar)
        fileloc=TRIM(filename)//TRIM(nchar)
        open(unit=ilun,file=fileloc,form='unformatted')
        ! Write grid variables
        write(ilun)ncpu_dump
        write(ilun)ndim
        write(ilun)nx,ny,nz
        write(ilun)nlevelmax
        write(ilun)ngridmax
        write(ilun)nboundary
        if (a_verbose>2) print*,myid,': ngrid_current =',ngrid_current
        write(ilun)ngrid_current
        write(ilun)boxlen
        ! Write time variables
        write(ilun)noutput,iout,ifout
        write(ilun)tout(1:noutput)
        write(ilun)aout(1:noutput)
        write(ilun)t
        write(ilun)dtold(1:nlevelmax)
        write(ilun)dtnew(1:nlevelmax)
        write(ilun)nstep,nstep_coarse
        write(ilun)const,mass_tot_0,rho_tot
        write(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
        write(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
        write(ilun)mass_sph
        ! Write levels variables
        write(ilun)headl(1:ncpu_dump,1:nlevelmax)
        write(ilun)taill(1:ncpu_dump,1:nlevelmax)
        write(ilun)numbl(1:ncpu_dump,1:nlevelmax)
        write(ilun)numbtot(1:10,1:nlevelmax)
        ! Read boundary linked list
        if(simple_boundary)then
           write(ilun)headb(1:nboundary,1:nlevelmax)
           write(ilun)tailb(1:nboundary,1:nlevelmax)
           write(ilun)numbb(1:nboundary,1:nlevelmax)
        end if
        ! Write free memory
        write(ilun)headf,tailf,numbf,used_mem,used_mem_tot
        ! Write cpu boundaries
        write(ilun)ordering
        if(ordering=='bisection') then
           write(ilun)bisec_wall(1:nbinodes)
           write(ilun)bisec_next(1:nbinodes,1:2)
           write(ilun)bisec_indx(1:nbinodes)
           write(ilun)bisec_cpubox_min(1:ncpu_dump,1:ndim)
           write(ilun)bisec_cpubox_max(1:ncpu_dump,1:ndim)
        else
           write(ilun)overload,ndomain,bound_key(0:ndomain),isub_domain(1:ndomain), &
                      cpu_domain(1:ndomain),map_domain(1:overload,1:ncpu_dump)
        endif
       
        ! Write coarse level
        write(ilun)son(1:ncoarse)
        write(ilun)flag1(1:ncoarse)
        write(ilun)cpu_map(1:ncoarse)
        ngrid_written = 0
        ! Write fine levels
        do ilevel=1,nlevelmax
           do ibound=1,nboundary+ncpu_dump
              if(ibound<=ncpu_dump)then
                 ncache=numbl(ibound,ilevel)
                 istart=headl(ibound,ilevel)
              else
                 ncache=numbb(ibound-ncpu_dump,ilevel)
                 istart=headb(ibound-ncpu_dump,ilevel)
              end if
              ngrid_written = ngrid_written + ncache
              if(ncache>0)then
                 allocate(iig(1:ncache))
                 ! Write grid index
                 igrid=istart
                 do i=1,ncache
                    iig(i)=igrid
                    igrid=next(igrid)
                 end do
                 write(ilun)iig
                 ! Write next index
                 igrid=istart
                 do i=1,ncache
                    iig(i)=next(igrid)
                    igrid=next(igrid)
                 end do
                 write(ilun)iig
                 ! Write prev index
                 igrid=istart
                 do i=1,ncache
                    iig(i)=prev(igrid)
                    igrid=next(igrid)
                 end do
                 write(ilun)iig
                 deallocate(iig)
                 ! Write grid center
                 allocate(xdp(1:ncache))
                 do idim=1,ndim
                    igrid=istart
                    do i=1,ncache
                       xdp(i)=xg(igrid,idim)
                       igrid=next(igrid)
                    end do
                    write(ilun)xdp
                 end do
                 deallocate(xdp)
                 allocate(iig(1:ncache))
                 ! Write father index
                 igrid=istart
                 do i=1,ncache
                    iig(i)=father(igrid)
                    igrid=next(igrid)
                 end do
                 write(ilun)iig
                 ! Write nbor index
                 do ind=1,twondim
                    igrid=istart
                    do i=1,ncache
                       iig(i)=nbor(igrid,ind)
                       igrid=next(igrid)
                    end do
                    write(ilun)iig
                 end do
                 ! Write son index
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    igrid=istart
                    do i=1,ncache
                       iig(i)=son(igrid+iskip)
                       igrid=next(igrid)
                    end do
                    write(ilun)iig
                 end do
                 ! Write cpu map
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    igrid=istart
                    do i=1,ncache
                       iig(i)=cpu_map(igrid+iskip)
                       igrid=next(igrid)
                    end do
                    write(ilun)iig
                 end do
                 ! Write refinement map
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    igrid=istart
                    do i=1,ncache
                       iig(i)=flag1(igrid+iskip)
                       igrid=next(igrid)
                    end do
                    write(ilun)iig
                 end do
                 deallocate(iig)
              end if
           end do
        end do
        if (a_verbose>2) print*,myid,' ngrid_written =',ngrid_written
        close(ilun)
     end if
#ifndef WITHOUTMPI
     call MPI_BARRIER(mpi_comm_use,info)
#endif
  end do
  
end subroutine backup_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_info'

  ilun=11

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu_dump
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')scale
  write(ilun,'("time        =",E23.15)')t8
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("H0          =",E23.15)')h0
  write(ilun,'("omega_m     =",E23.15)')omega_m
  write(ilun,'("omega_l     =",E23.15)')omega_l
  write(ilun,'("omega_k     =",E23.15)')omega_k
  write(ilun,'("omega_b     =",E23.15)')omega_b
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
  write(ilun,*)
  
  ! Write ordering information
  write(ilun,'("ordering type=",A80)')ordering
  if(ordering=='bisection') then
     do icpu=1,ncpu_dump
        ! write 2*ndim floats for cpu bound box                                                         
        write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
        ! write 1 float for cpu load                                                                    
        write(ilun,'(E23.15)')dble(bisec_cpu_load(icpu))
     end do
  else
     write(ilun,'("   DOMAIN   ind_min                 ind_max")')
     do idom=1,ncpu_dump*overload
        write(ilun,'(I8,1X,F32.0,1X,F32.0)')idom,bound_key(idom-1),bound_key(idom)
     end do
  endif

  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::info,npart_tot,ilun
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_header'

  ! Compute total number of particles
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,mpi_comm_use,info)
#endif
#ifdef WITHOUTMPI
  npart_tot=npart
#endif

  if(myid==1)then

     ilun=11

     ! Open file
     fileloc=TRIM(filename)
     open(unit=ilun,file=fileloc,form='formatted')
     
     ! Write header information
     write(ilun,*)'Total number of particles'
     write(ilun,*)npart_tot
     write(ilun,*)'Total number of dark matter particles'
     write(ilun,*)npart_tot-nstar_tot
     write(ilun,*)'Total number of star particles'
     write(ilun,*)nstar_tot
     write(ilun,*)'Total number of sink particles'
     write(ilun,*)nsink
     
     close(ilun)

  endif

end subroutine output_header

