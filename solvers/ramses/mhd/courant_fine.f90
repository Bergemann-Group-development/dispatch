subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ix,iy,iz,ig,ivar,idim,ind,ncache,igrid,iskip,kmin
  integer::info,nleaf,ngrid,nx_loc,icpu(1)
  integer,dimension(1:nvector) :: ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,urms_loc,vol_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_tmp,eint_locp
  real(kind=8), dimension(:), allocatable :: dt_all
  real(kind=8),dimension(11)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar+3) :: uu
  real(dp),dimension(1:nvector,1:ndim)   :: gg
  real(dp),dimension(1:twotondim,1:ndim) :: xc
  integer, save :: last_ilevel=0
  logical       :: compute_urms, nan_found, nan_tmp
  logical, save :: first_time=.true.
  real(dp),save :: dln_factor(1:MAXLEVEL)=1.0_dp

  logical, external :: is_nan

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nan_found = .false.

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
                  urms_loc=0.0d0
                   vol_loc=0.0d0
                    dt_loc=dtnew(ilevel)

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
#if NDIM>1
     xc(ind,2)=(dble(iy)-0.5D0)*dx
#endif
#if NDIM>2
     xc(ind,3)=(dble(iz)-0.5D0)*dx
#endif
  end do

  if (ilevel==levelmin) then
    courant_number%cb = 0.0_dp
    courant_number%cv = 0.0_dp
    courant_number%cg = 0.0_dp
  endif

  ! Only add to urms if we are at the base level, or if this is first time in a substep
  compute_urms = ilevel==levelmin .or. ilevel.gt.last_ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  !$omp parallel default(none) firstprivate(ncache,vol,xc,dx) &
  !$omp   reduction(+:urms_loc,vol_loc,mass_loc,ekin_loc,emag_loc) &
  !$omp   private(igrid,ngrid,ig,i,ind_grid,ind,iskip,ind_cell,nleaf, &
  !$omp           ind_leaf,uu,ivar,gg,idim,kmin,dt_lev,eint_locp,nan_tmp) &
  !$omp   shared(active,ilevel,ncoarse,ngridmax,son,uold,f,compute_urms,poisson,xg, &
  !$omp          dt_loc,eint_loc,courant_number,nan_found,myid)
  eint_locp=0.0_8
  !$omp do
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do
        
        ! Gather hydro variables
#ifdef CHECK_NAN
        do i=1,nleaf
           nan_tmp = .false.
           do ivar=1,nvar+3
              uu(i,ivar)=uold(ind_leaf(i),ivar)
              nan_tmp = nan_tmp .or. is_nan(uu(i,ivar))
           end do
           if (nan_tmp) then
              !$omp critical
              if (.not. nan_found) then
                 ig = ind_leaf(i) - iskip
                 write (*,'(a,i5.5,a,i3,a,3(1x,f10.7))') &
                   'R=',myid,' lvl=', ilevel,' Hit a NaN cell. pos=:', xg(ig,:)+xc(ind,:)
                 do idim=1,nvar
                    if (twotondim==2) write (*,'(a,i5.5,a,i3,a,i2,a,2(1x,e13.6))') &
                       'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1/)*ngridmax,idim)
                    if (twotondim==4) write (*,'(a,i5.5,a,i3,a,i2,a,4(1x,e13.6))') &
                       'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1,2,3/)*ngridmax,idim)
                    if (twotondim==8) write (*,'(a,i5.5,a,i3,a,i2,a,8(1x,e13.6))') &
                       'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1,2,3,4,5,6,7/)*ngridmax,idim)
                 enddo
                 nan_found = .true.
              end if
              !$omp end critical
           end if
        end do
#else
        do ivar=1,nvar+3
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
#ifndef _OPENMP
           call bitwise(uu,sizeof(uu),'courant ')
#endif
        end do

#endif
        
        ! Compute Urms
        if (compute_urms) then
           do i=1,nleaf
              gg(i,1)=1.0_dp/uu(i,1)
           end do
           do ivar=1,ndim
           do i=1,nleaf
              urms_loc= urms_loc + (uu(i,1+ivar)*gg(i,1))**2*vol
           end do
           end do
#ifdef CHECK_NAN
           if (is_nan(urms_loc)) then ! We hit a bad cell
              do ivar=1,ndim
              do i=1,nleaf
                 if (nan_found) exit
                 if (is_nan(uu(i,1+ivar)*gg(i,1))) then
                    !$omp critical
                    if (.not. nan_found) then ! Only print for first oct
                       ig = ind_leaf(i) - iskip
                       write (*,'(a,i5.5,a,i3,a,3(1x,f10.7))') 'R=',myid,' lvl=', ilevel,' Hit a NaN cell. pos=:', xg(ig,:)+xc(ind,:)
                       do idim=1,nvar
                          if (twotondim==2) write (*,'(a,i5.5,a,i3,a,i2,a,2(1x,e13.6))') &
                             'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1/)*ngridmax,idim)
                          if (twotondim==4) write (*,'(a,i5.5,a,i3,a,i2,a,4(1x,e13.6))') &
                             'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1,2,3/)*ngridmax,idim)
                          if (twotondim==8) write (*,'(a,i5.5,a,i3,a,i2,a,8(1x,e13.6))') &
                             'R=',myid,' lvl=', ilevel,' ivar=',idim,' uold=', uold(ig + ncoarse+(/0,1,2,3,4,5,6,7/)*ngridmax,idim)
                       enddo
                       nan_found = .true.
                    endif
                    !$omp end critical
                 endif
              end do
              end do
           endif
#endif
           vol_loc = vol_loc + nleaf*vol
        end if
        
        ! Gather gravitational acceleration
        gg=0.0_dp
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        do i=1,nleaf
           ! Compute total mass
           mass_loc=mass_loc+uu(i,1)*vol
           ! Compute total energy
           ekin_loc=ekin_loc+uu(i,5)*vol
           ! Compute total magnetic energy
           do ivar=1,3
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
           ! Compute total internal energy
           eint_locp=eint_locp+uu(i,5)*vol
           do ivar=1,3
              eint_locp=eint_locp-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
        
        ! Compute CFL time-step
        if(nleaf>0)then
           !$omp critical
           dt_lev = dt_loc
           !$omp end critical
           call cmpdt(uu,gg,dx,dt_lev,nleaf,kmin,dt_loc,ilevel)
           if (kmin > 0) then
              !$omp critical
              if (dt_lev <= dt_loc) then
                 courant_number%x = xg(ind_leaf(kmin)-iskip,:)+xc(ind,:)
                 courant_number%rho = uold(ind_leaf(kmin),1)!uu(kmin,1)
                 courant_number%etot = uold(ind_leaf(kmin),5)!uu(kmin,5)
                 courant_number%level = ilevel
                 dt_loc=min(dt_loc,dt_lev)
              end if
              !$omp end critical
           end if
           !$omp critical
           dt_loc=min(dt_loc,dt_lev)
           !$omp end critical
        end if
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids
  !$omp enddo nowait
  !$omp atomic
  eint_loc = eint_loc + eint_locp
  !$omp end parallel

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  comm_buffin(4)=emag_loc
  comm_buffin(5)=urms_loc
  comm_buffin(6)=vol_loc

  call MPI_ALLREDUCE(comm_buffin,comm_buffout,6,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &mpi_comm_use,info)

  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
  courant_number%urms   = courant_number%urms + comm_buffout(5)
  courant_number%volume = courant_number%volume + comm_buffout(6)

  ! Find size of min dt, and which thread has it
  allocate(dt_all(ncpu))

  call MPI_ALLGATHER(dt_loc, 1, MPI_DOUBLE_PRECISION, &
                     dt_all, 1, MPI_DOUBLE_PRECISION, mpi_comm_use, info)

  icpu = minloc(dt_all)
  dt_loc = dt_all(icpu(1))
  deallocate(dt_all)

  ! Reduce dt with dln_factor, based on how much too large dlnr or dlnp was in the previous step,
  ! which might already have had a dln_factor larger than unity. While the factor is above one,
  ! reduce it with at most a factor 0.9, to avoid rapidly changing dt, and to retain memory of
  ! "bad times"
  if (dlnr_max==0.0) dlnr_max=-1.0_dp
  if (dlnp_max==0.0) dlnp_max=-1.0_dp
  dln_factor(ilevel) = max(1.0_dp,dln_factor(ilevel)*max(0.9_dp,dlnr(ilevel)/dlnr_max,dlnp(ilevel)/dlnp_max))
  dt_loc = dt_loc/dln_factor(ilevel)

  ! Broadcast courant numbers from thread with min dt
  if (myid==icpu(1)) then
    comm_buffin(1)=courant_number%cb
    comm_buffin(2)=courant_number%cv
    comm_buffin(3)=courant_number%cg
    comm_buffin(4:3+ndim)=courant_number%x
    comm_buffin(7)=courant_number%rho
    comm_buffin(8)=courant_number%umag
    comm_buffin(9)=courant_number%bmag
    comm_buffin(10)=courant_number%etot
    comm_buffin(11)=courant_number%level
  endif
  
  call MPI_BCAST(comm_buffin,11,MPI_DOUBLE_PRECISION,icpu(1)-1,mpi_comm_use,info)

  courant_number%cb=comm_buffin(1)
  courant_number%cv=comm_buffin(2)
  courant_number%cg=comm_buffin(3)
  courant_number%x =comm_buffin(4:3+ndim)
  courant_number%rho=comm_buffin(7)
  courant_number%umag=comm_buffin(8)
  courant_number%bmag=comm_buffin(9)
  courant_number%etot=comm_buffin(10)
  courant_number%level=comm_buffin(11)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
  courant_number%urms = urms_loc
  courant_number%volume = vol_loc
#endif

  courant_number%cr = dlnr(ilevel)
  courant_number%cp = dlnp(ilevel)

  if (compute_urms) then
     mass_tot=mass_tot+mass_all
     ekin_tot=ekin_tot+ekin_all
     eint_tot=eint_tot+eint_all
     emag_tot=emag_tot+emag_all
     last_ilevel = ilevel
  end if

  dtlev(ilevel)=dt_loc
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_loc)
  dt_tmp=dt_loc
  do i=ilevel,nlevelmax-1
     dt_tmp=dt_tmp/nsubcycle(i)
  end do

  if (myid==1 .and. h_verbose>0) then
    if (first_time) then
      first_time=.false.
      open(unit=courant_unit,file=TRIM(datadir)//'/courant.log',form='formatted',status='unknown')
    end if
    if (h_verbose>1) &
      write(courant_unit,'(a,i7,a,i3,1p,2(a,e13.6),4(a,e12.5),a,0p,3f11.7)') &
        ' step=' , nstep                , &
        ' level=', ilevel               , &
        ' t='    , t                    , &
        ' dt='   , dt_tmp               , &
        ' rho='  , courant_number%rho   , &
        ' umag=' , courant_number%umag  , &
        ' bmag=' , courant_number%bmag  , &
        ' etot=' , courant_number%etot  , &
        ' pos='  , courant_number%x

    if (ilevel==nlevelmax) & 
      write(courant_unit,'(a,i7,a,i3,1p,2(a,e13.6),4(a,e12.5),a,0p,3f11.7)') &
        ' step=' , nstep                , &
        ' level=', courant_number%level , &
        ' t='    , t                    , &
        ' dt='   , dt_tmp               , &
        ' rho='  , courant_number%rho   , &
        ' umag=' , courant_number%umag  , &
        ' bmag=' , courant_number%bmag  , &
        ' etot=' , courant_number%etot  , &
        ' pos='  , courant_number%x
  end if

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
 
  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        
        ! Impose analytical velocity field
        !call velana(xx,vv,dx_loc,t,ngrid)
        
        ! Impose induction variables
        do i=1,ngrid
           uold(ind_cell(i),1)=1.0
        end do
        do idim=1,3
           do i=1,ngrid
              uold(ind_cell(i),idim+1)=vv(i,idim)
           end do
        end do
        ! Update total energy
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           uold(ind_cell(i),neul)=1.0+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
        end do
        
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
