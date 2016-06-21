!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uf_debug (ind_grid,ncache,ilevel,label)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer ncache, ind_grid(ncache)
  integer ilevel, i, ind, ind_cell
  !logical in_debug_region
  character(len=*) label

!  call trace_mpi('uf_debug','enter',3)

!  if (local_debug) then
!     do i=1,ncache
!        if (in_debug_region(ind_grid(i))) then
!           do ind=1,twotondim
!              ind_cell = ind_grid(i)+ncoarse+(ind-1)*ngridmax
!              print'(a,i8,i5,i3,2i8,i3,1p,e12.5,5(1x,3e12.5))', TRIM(label), &
!                nstep_coarse,myid,ilevel,ind_grid(i),ind_cell,ind,unew(ind_cell,1:nvar),f(ind_cell,:)
!           end do
!        endif
!     end do
!  end if

!  call trace_mpi('uf_debug','exit',3)
end subroutine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector)::ind_grid
  real(kind=8) :: speed
  character(len=80):: id='godunov_fine.f90 $Id$'

  call print_id(id)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
!  call trace_mpi('godunov_fine(ilevel)','enter',2)

  courant_number%speed = 0.0_dp

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  n_cell=n_cell+ncache
  !$omp parallel do private(igrid,ngrid,i,ind_grid) schedule(guided)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

  if (h_verbose>1 .or. (h_verbose>0 .and. ilevel==levelmin)) then
     call sum_ints_mpi (n_entropy,2)
     if (myid==1) print*,'n_entropy_test,fix =', n_entropy
     n_entropy=0
  end if

111 format('   Entering godunov_fine for level ',i2)

!  call trace_mpi('godunov_fine(ilevel)','exit',2)
end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel,icount)
  use amr_parameters
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,icount
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,j,ivar,ind,icpu,iskip,not_iso
  real(dp)::d,u,v,w,e,A,B,C,fact,smalle,e_kin,e_mag,ca2,cs2,cc2,ctot,smallt2
  logical in_debug_region,sml_e
  logical, save :: first_time=.true.
  real(dp)      :: scale_l,scale_t,scale_d,scale_v,scale_nH
  real(dp),save :: scale_T2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  not_iso=0
  n_camax=0
  n_cbmax=0
  smallt2=1.5 * 2.73 / 2.   ! Limit for T / mu - set it to 50% above T_CMB to allow for variations of gamma and mu
  if (first_time) then
     first_time = .false.
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)                 ! gets unit conversion factors from RAMSES units to cgs
  end if

!  call trace_mpi('set_unew(ilevel,icount)','enter',2)
  call check_internal_energy(uold,ilevel,'set_unew1',0)

  ! Set unew to uold for myid cells
  !$omp parallel private(ind,iskip,ivar,i,j,d,u,v,w,e,A,B,C,e_kin,e_mag,sml_e, &
  !$omp&                 smalle,fact,ca2,cs2,cc2,ctot)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     !$omp do reduction(min:dmin) reduction(max:dmax)
     do i=1,active(ilevel)%ngrid
        j=active(ilevel)%igrid(i)+iskip
        if (son(j)>0) cycle
        d=uold(j,1)
        dmin=min(dmin,d)
        dmax=max(dmax,d)
        uold(j,1)=max(uold(j,1),smallr)
     end do
     !$omp enddo nowait
     if(camax > 0.0)then
        !$omp do reduction(+:n_camax)
        do i=1,active(ilevel)%ngrid
           j=active(ilevel)%igrid(i)+iskip
           if (son(j)>0) cycle
           d=uold(j,1)
           u=uold(j,2)/d
           v=uold(j,3)/d
           w=uold(j,4)/d
           A=uold(j,6)+uold(j,nvar+1)
           B=uold(j,7)+uold(j,nvar+2)
           C=uold(j,8)+uold(j,nvar+3)
           e_mag=0.125*(A*A+B*B+C*C)
           e_kin=0.5*d*(u*u+v*v+w*w)
           e=uold(j,5)-e_kin-e_mag
           !cs2 = gamma*(gamma-1_dp)*e/d
           !cc2 = 0.5_dp*((A*A+B*B+C*C)/d+cs2)
           !ctot =        abs(u) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*A*A/d,0.0_dp)))
           !ctot = ctot + abs(v) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*B*B/d,0.0_dp)))
           !ctot = ctot + abs(w) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*C*C/d,0.0_dp)))
           ctot = 3_dp*(A*A+B*B+C*C)/d
           ca2=(ctot/camax)**2
           if (ca2 > 1.0 .and. uold(j,1) < daver) then
             uold(j,1)=uold(j,1)*ca2
#ifdef KROME
             uold(j,ichem+1:nvar)=uold(j,ichem+1:nvar)*ca2
#endif
             e_kin=e_kin/ca2
             uold(j,5)=e+e_kin+e_mag
             if (debug .and. debug_level>0) then
                n_camax=n_camax+1
             end if
           end if
        end do
        !$omp enddo nowait
     end if
     if (subcycle_extrap) then
        if (nstep_coarse==0 .or. ilevel==nlevelmax .or. icount==1) then
          fact = 0.0
        else
          fact = min(dtold(ilevel+1)/dtold(ilevel),0.6_dp)   ! extrapolate to the time of the next substep
        end if
        do ivar=1,nvar+3
           !$omp do
           do i=1,active(ilevel)%ngrid
              unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
              uold(active(ilevel)%igrid(i)+iskip,ivar) = &
                (1.0d0+fact)*uold(active(ilevel)%igrid(i)+iskip,ivar) - &
                       fact *usav(active(ilevel)%igrid(i)+iskip,ivar)
              usav(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
           end do
           !$omp enddo nowait
        end do
     else
        do ivar=1,nvar+3
           !$omp do
#if defined(__MIC__)
           !DEC$ NOVECTOR
#endif
           do i=1,active(ilevel)%ngrid
              unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
           end do
           !$omp enddo nowait
        end do
     end if
     if (t_damp > 0.0) then                                     ! option for artificial damping / friction
        !$omp do
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,2) = unew(active(ilevel)%igrid(i)+iskip,2)*(1.-dtnew(ilevel)/t_damp)
           unew(active(ilevel)%igrid(i)+iskip,3) = unew(active(ilevel)%igrid(i)+iskip,3)*(1.-dtnew(ilevel)/t_damp)
           unew(active(ilevel)%igrid(i)+iskip,4) = unew(active(ilevel)%igrid(i)+iskip,4)*(1.-dtnew(ilevel)/t_damp)
        end do
        !$omp enddo nowait
     end if
     if(pressure_fix)then
        !$omp do
        do i=1,active(ilevel)%ngrid
           j=active(ilevel)%igrid(i)+iskip
           if (son(j)>0) cycle
           divu(j) = 0.0_dp
           d=unew(j,1)
           u=unew(j,2)/d
           v=unew(j,3)/d
           w=unew(j,4)/d
           A=unew(j,6)+unew(j,nvar+1)
           B=unew(j,7)+unew(j,nvar+2)
           C=unew(j,8)+unew(j,nvar+3)
           e=unew(j,5)-0.5_dp*d*(u**2+v**2+w**2)-0.125_dp*(A**2+B**2+C**2)
           enew(j)=e
        end do
        !$omp enddo nowait
     end if
     if (do_isothermal .and. debug) then
        !$omp do reduction(+:not_iso)
        do i=1,active(ilevel)%ngrid
           j=active(ilevel)%igrid(i)+iskip
           if (son(j)>0) cycle
           d=unew(j,1)
           u=unew(j,2)
           v=unew(j,3)
           w=unew(j,4)
           A=unew(j,6)+unew(j,nvar+1)
           B=unew(j,7)+unew(j,nvar+2)
           C=unew(j,8)+unew(j,nvar+3)
           e=unew(j,5)-0.5_dp*(u**2+v**2+w**2)/d-0.125_dp*(A**2+B**2+C**2)
           not_iso = not_iso + merge(1, 0, (abs((e/d)/temp_iso-1.0_dp) > 0.05))
        end do
        !$omp enddo nowait
     else if (debug) then
        fact = smallt2 / (scale_T2 * (gamma-1.0_dp))
        !$omp do reduction(+:not_iso)
        do i=1,active(ilevel)%ngrid
           j=active(ilevel)%igrid(i)+iskip
           if (son(j)>0) cycle
           d=unew(j,1)
           u=unew(j,2)
           v=unew(j,3)
           w=unew(j,4)
           A=unew(j,6)+unew(j,nvar+1)
           B=unew(j,7)+unew(j,nvar+2)
           C=unew(j,8)+unew(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d
           e_mag=0.125_dp*(A*A+B*B+C*C)
#ifdef KROME
           if (chemistry) fact = smallt2 / (scale_T2 * (unew(j,ichem)-1.0_dp))   ! Use cell dependent effective adiabatic index
#endif
           smalle = fact * d
           sml_e=(unew(j,5)-e_kin-e_mag) < smalle
           not_iso=not_iso+merge(1,0,sml_e)
           unew(j,5)=merge(smalle+e_kin+e_mag,unew(j,5),sml_e)
        end do
     endif
  end do
  !$omp end parallel

  if (debug) then
     if (debug_level>4) then                  ! 5: output for each relevant rank and level
        if ((not_iso+n_camax+n_cbmax)>0) then
           if (do_isothermal) then
              write(*,'(a,i5,a,i3,a,3i11)') 'set_unew: rank=',myid,' lvl=',ilevel,' not_iso, n_camax, n_cbmax =',not_iso,n_camax,n_cbmax
           else
              write(*,'(a,i5,a,i3,a,3i11)') 'set_unew: rank=',myid,' lvl=',ilevel,' small e, n_camax, n_cbmax =',not_iso,n_camax,n_cbmax
           end if
        end if
     else if (debug_level>2) then             ! 4: output for each level
        call sum_int_mpi(not_iso)
        call sum_int_mpi(n_camax)
        call sum_int_mpi(n_cbmax)
        if ((not_iso+n_camax+n_cbmax)>0 .and. myid==1) then
           if (do_isothermal) then
              write(*,'(a,i3,a,3i11)') 'set_unew: lvl=', ilevel, ' not_iso, n_camax, n_cbmax =', not_iso, n_camax, n_cbmax
           else
              write(*,'(a,i3,a,3i11)') 'set_unew: lvl=', ilevel, ' small e, n_camax, n_cbmax =', not_iso, n_camax, n_cbmax
           end if
        end if
     else if (debug_level==2) then            ! 2: output for each Fine step
        not_iso_tot = not_iso_tot + not_iso
        n_camax_tot = n_camax_tot + n_camax
        n_cbmax_tot = n_cbmax_tot + n_cbmax
        if (ilevel==nlevelmax) then
           call sum_int8_mpi (not_iso_tot)
           call sum_int8_mpi (n_camax_tot)
           call sum_int8_mpi (n_cbmax_tot)
           if ((not_iso_tot+n_camax_tot+n_cbmax_tot)>0 .and. myid==1) then
              if (do_isothermal) then
                 write(*,'(a,3i11)') 'set_unew: not_iso, n_camax, n_cbmax =', not_iso_tot, n_camax_tot, n_cbmax_tot
              else
                 write(*,'(a,3i11)') 'set_unew: small e, n_camax, n_cbmax =', not_iso_tot, n_camax_tot, n_cbmax_tot
              end if
           end if
           not_iso_tot = 0
           n_camax_tot = 0
           n_cbmax_tot = 0
        end if
     else if (debug_level==1) then            ! 1: output for each coarse step (in adaptive_loop)
        not_iso_tot = not_iso_tot + not_iso
        n_camax_tot = n_camax_tot + n_camax
        n_cbmax_tot = n_cbmax_tot + n_cbmax
     end if
  endif

  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,reception(icpu,ilevel)%ngrid
        unew(reception(icpu,ilevel)%igrid(i)+iskip,1:nvar+3)=0.0
     end do
     if(pressure_fix)then
        do i=1,reception(icpu,ilevel)%ngrid
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end do
     end if
  end do
  end do

  call min_real_mpi (dmin)
  call max_real_mpi (dmax)

  call check_internal_energy(unew,ilevel,'set_unew2',1)
111 format('   Entering set_unew for level ',i2)

!  call trace_mpi('set_unew(ilevel,icount)','exit',2)
end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine check_internal_energy(uu,ilevel,name,level)
  use amr_commons
  use hydro_commons
  implicit none
  real(kind=dp):: uu(ncoarse+twotondim*ngridmax,nvar+3)
  integer :: ilevel,level
  character(len=*) :: name
  integer :: ind,iskip,i,j,ix,not_iso,gnot_iso,nr_rank,nr_rankt,info,allred(5),gallred(5),nsmallt,ncell
  real(kind=dp) :: d,u,v,w,A,B,C,e_kin,e_mag,p,smalle,smallt2,fact
  integer, dimension(:), pointer :: g
  real(kind=dp)      :: eps=0.001
  integer, save:: nprint=20, nprintd=20, nprintx=100, nprintx2=100
  integer, save, dimension(3):: n_small=0
#ifdef KROME
  integer, save, dimension(ichem_param+1:nvar):: x_small=0
  real(kind=dp)      :: totX
#endif
  logical, save :: first_time=.true.
  real(dp)      :: scale_l,scale_t,scale_d,scale_v,scale_nH
  real(dp),save :: scale_T2
  logical       :: in_debug_region

  ! If parameter not set, return
  if (.not. debug.or.debug_level<level) return
!  call trace_mpi('check_internal_energy','enter',2)
  not_iso=0
  nsmallt=0
  ncell=0
  smallt2=(1.0_dp - 1e-6)*1.5*2.73/2. ! Limit for T / mu - set it to 50% above T_CMB to allow for variations of gamma and mu

  if (first_time) then
     first_time = .false.
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)                 ! gets unit conversion factors from RAMSES units to cgs
  end if

  ! Check density
  !$omp parallel do private(ind,iskip,g,i,j,d) shared(n_small,ilevel) reduction(+:ncell)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
     do i=1,active(ilevel)%ngrid
        j = g(i)+iskip
        if (son(j)>0) cycle
        ncell = ncell + 1
        d = uu(j,1)
        if (d < 0_dp) then
           !$omp atomic
           n_small(1) = n_small(1) + 1
           if (debug_level>5 .and. nprintd>0) then
              print 1,'Negative density at xg =', xg(g(i),:), d, ilevel
            1 format(a,3f12.9,1p,e12.2,i5)
              if (debug_level<7) nprintd=nprintd-1
           end if
        else if (d == 0_dp) then
           !$omp atomic
           n_small(2) = n_small(2) + 1
           if (debug_level>5 .and. nprintd>0) then
              print 1,'Zero density at xg =', xg(g(i),:), d, ilevel
              if (debug_level<7) nprintd=nprintd-1
           end if
        else if (d < smallr) then
           !$omp atomic
           n_small(3) = n_small(3) + 1
           if (debug_level>5 .and. nprintd>0) then
              print 1,'Small density at xg =', xg(g(i),:), d, ilevel
              if (debug_level<7) nprintd=nprintd-1
           end if
        end if
     end do
  end do

  ! Print diagnostic -- if debug_level==1 accumulate until last ilevel
  if (debug_level>2) then
     if (sum(n_small)>0) then
        print'(a,i7,a,i3,a,3i8)','Small densities in rank',myid,' at level',ilevel,': neg, zero, small =',n_small 
        n_small = 0
     end if
  else if (debug_level==2) then
     call sum_ints_mpi (n_small,3)
     if (sum(n_small)>0) then
       if (myid==1) print'(a,i3,a,3i8)','Small densities at level',ilevel,': neg, zero, small =',n_small 
        n_small = 0
     end if
  else if (debug_level==1 .and. ilevel==levelmin) then
     call sum_ints_mpi (n_small,3)
     if (sum(n_small)>0) then
        if (myid==1) print'(a,3i8)','Small densities: neg, zero, small =',n_small 
        n_small = 0
     end if
  end if

#ifdef KROME
  ! Check density of chemical species
  !$omp parallel do default(none) shared(x_small,ilevel,ichem,nprintx,nprintx2,uu, &
  !$omp                                  ngridmax,ncoarse,xg,active,debug_level,name) &
  !$omp   private(ind,iskip,g,i,j,d,ix,totX)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
     do i=1,active(ilevel)%ngrid
        j = g(i)+iskip
        !if (son(j)>0) cycle
        totX = 0.0_dp
        do ix=ichem+1,nvar
           d = uu(j,ix)
           totX = totX + d
           if (d < 0_dp) then
              !$omp atomic
              x_small(ix) = x_small(ix) + 1
              if (debug_level>5 .and. nprintx>0) then
                 print 2,'Negative density for species ',ix-ichem,' at xg =', xg(g(i),:), d, ilevel
               2 format(a,i3,a,3f12.9,1p,e12.2,i5)
                 if (debug_level<7) then
                    !$omp atomic
                    nprintx=nprintx-1
                 end if
              end if
           end if
        end do
        if (abs(totX/uu(j,1)-1.0_dp) > 1e-5_dp) then
           if (nprintx2>0) then
              print '(a,3f12.9,a,i3,1p,3e12.2)', &
                      'Abundances do not sum at '//trim(name)//'. xg =', xg(g(i),:), &
                      ' lvl,d,totX=', ilevel, uu(j,1), totX, totX/uu(j,1)-1.
              if (totX / uu(j,1) <0.125 .or. totX / uu(j,1) > 8.) &
                print '(a,99e10.2)', 'Abundances :', uu(j,ichem+1:nvar), sum(uu(j,ichem+1:nvar)), totX, uu(j,1)
              if (debug_level<7) then
                 !$omp atomic
                 nprintx2=nprintx2-1
              end if
           end if
           uu(j,ichem+1:nvar) = uu(j,1) / totX * uu(j,ichem+1:nvar) ! Enforce mass abundances sum to 1
        end if
     end do
  end do

  ! Print diagnostic -- if debug_level==1 accumulate until last ilevel
  if (debug_level>2) then
     do ix=ichem+1,nvar
        if (x_small(ix)>0) then
           print'(2(a,i3),a,i7,a,i8)', &
              'Negative densities at '//trim(name)//' for species ', ix-ichem, &
              ' in rank',myid,' at level',ilevel,': ',x_small(ix) 
           x_small(ix) = 0
        end if
     end do
  else if (debug_level==2) then
     call sum_ints_mpi (x_small,nvar-ichem)
     do ix=ichem+1,nvar
        if (x_small(ix)>0) then
          if (myid==1) print'(2(a,i3),a,i8)',&
             'Negative densities at '//trim(name)//' for species', ix-ichem, &
             ' at level',ilevel,': ',x_small(ix)
           x_small(ix) = 0
        end if
     end do
  else if (debug_level==1 .and. ilevel==levelmin) then
     call sum_ints_mpi (x_small,nvar-ichem)
     do ix=ichem+1,nvar
        if (x_small(ix)>0) then
           if (myid==1) print'(a,i3,a,i8)','Negative densities at '//trim(name)//' for species', ix-ichem,': ',x_small(ix)
           x_small(ix) = 0
        end if
     end do
  end if
#endif

  ! Check if energy is isothermal
  if (do_isothermal) then
     !$omp parallel do private(ind,iskip,g,i,j,d,u,v,w,A,B,C,e_kin,e_mag) reduction(+:not_iso)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
        do i=1,active(ilevel)%ngrid
           j = g(i)+iskip
           if (son(j)>0) cycle
           d=uu(j,1)
           u=uu(j,2)
           v=uu(j,3)
           w=uu(j,4)
           A=uu(j,6)+uu(j,nvar+1)
           B=uu(j,7)+uu(j,nvar+2)
           C=uu(j,8)+uu(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d
           e_mag=0.125_dp*(A*A+B*B+C*C)
           if (abs((uu(j,5)-e_kin-e_mag)/d-temp_iso)/temp_iso > eps) not_iso = not_iso + 1
        end do
     end do
  else
     ! Check if energy is positive and temperature not too low
     !$omp parallel private(ind,iskip,g,i,j,d,u,v,w,A,B,C,e_kin,e_mag,smalle,fact)
     fact = smallt2 / (scale_T2 * (gamma-1.0_dp))
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
        !$omp do reduction(+:not_iso,nsmallt)
        do i=1,active(ilevel)%ngrid
           j = g(i)+iskip
           if (son(j)>0) cycle
           d=uu(j,1)
           u=uu(j,2)
           v=uu(j,3)
           w=uu(j,4)
           A=uu(j,6)+uu(j,nvar+1)
           B=uu(j,7)+uu(j,nvar+2)
           C=uu(j,8)+uu(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d
           e_mag=0.125_dp*(A*A+B*B+C*C)
#ifdef KROME
           if (chemistry) fact = smallt2 / (scale_T2 * (uu(j,ichem)-1.0_dp)) ! Use cell dependent effective adiabatic index
#endif
           smalle = fact * d
           nsmallt= nsmallt + merge(1,0,uu(j,5) < e_kin + e_mag + smalle)
           if (uu(j,5)-e_kin-e_mag < 0.0_dp) then
              not_iso = not_iso + 1
              if (debug_level>5 .and. nprint>0) then
                 print'(a,i9,i5,3f11.7,1p,e11.2,3(1x,3e10.2))', &
                   'pressure error; j,lev,xg,d,U,B,ek,em,ei:', &
                    j,ilevel,xg(g(i),:),d,u,v,w,0.5*A,0.5*B,0.5*C, &
                    e_kin/uu(j,5),e_mag/uu(j,5),(uu(j,5)-e_kin-e_mag)/uu(j,5)
                 if (debug_level<7) then
                    !$omp atomic
                    nprint=nprint-1
                 end if
              end if
           end if
        end do
     end do
     !$omp end parallel
  end if
#ifndef WITHOUTMPI
  allred = (/ncell,not_iso,nsmallt,merge(1_4,0_4,not_iso>0),merge(1_4,0_4,nsmallt>0)/)
  call MPI_ALLREDUCE(allred,gallred,5,MPI_INTEGER,MPI_SUM,mpi_comm_use,info)
  ncell   = gallred(1)
  not_iso = gallred(2)
  nsmallt = gallred(3)
  nr_rank = gallred(4)
  nr_rankt= gallred(5)
#else
  nr_rank =1
  nr_rankt=1
#endif
  if (debug_level>3 .and. myid==1) then
     if (do_isothermal) then
        if (not_iso > 0) &
          write(*,'(a,i5,1x,a,i3,1x,a,i11,a12,i12)') &
            'Not isothermal at '//trim(name)//' #ranks=',nr_rank,'lvl=',ilevel,' #cells=',not_iso, ' #tot cells=', ncell
     else
        if (not_iso > 0) &
          write(*,'(a,i5,1x,a,i5,1x,a,i11,a12,i12)') &
            'Negative pressure at '//trim(name)//' #ranks=',nr_rank,'lvl=',ilevel,' #cells=',not_iso, ' #tot cells=', ncell
        if (nsmallt > 0) &
          write(*,'(a,i5,1x,a,i5,1x,a,i11,a12,i12)') &
            'Small pressure at '//trim(name)//' #ranks=',nr_rankt,'lvl=',ilevel,' #cells=',nsmallt, ' #tot cells=', ncell
     endif
  endif

!  call trace_mpi('check_internal_energy','exit',2)
end subroutine check_internal_energy
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_parameters
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the
  ! hydro step.
  !--------------------------------------------------------------------------
  integer::i,j,ivar,ind,iskip,nx_loc,not_iso
  real(dp)::scale,id,d,u,v,w,A,B,C,smallp,dlnr1,dlnp1,e_old,e_new,dlnp0,e_min,smallt2
  real(dp)::e_mag,e_kin,e_prim,e_trunc,div,dx,fact,d_old,d_new
  real, allocatable, dimension(:) :: ubw
  integer, dimension(:), pointer :: g
  logical, save :: first_time=.true.
  real(dp)      :: scale_l,scale_t,scale_d,scale_v,scale_nH
  real(dp),save :: scale_T2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

!  call trace_mpi('set_uold(ilevel)','enter',2)
  call check_internal_energy(unew,ilevel,'set_uold1',0)

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale
  smallp = smallr*smallc**2/gamma
  not_iso=0

  smallt2=1.5 * 2.73 / 2.   ! Limit for T / mu - set it to 50% above T_CMB to allow for variations of gamma and mu
  if (first_time) then
     first_time = .false.
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)                 ! gets unit conversion factors from RAMSES units to cgs
  end if

  dlnr(ilevel)=0.0_dp                                                           ! reset measurements of density rate of change
  dlnp(ilevel)=0.0_dp                                                           ! reset measurements of pressure rate of change

  if (bitwise_test > 0) then
     allocate (ubw(active(ilevel)%ngrid))
     if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,nvar+3
           do i=1,active(ilevel)%ngrid
             j = g(i)+iskip
             ubw(i) = unew(j,ivar)
           end do
           call bitwise (ubw, active(ilevel)%ngrid*sizeof(ubw(1)), 'godunovf')
        end do
     end do
     deallocate (ubw)
  endif

  ! Add gravity source term at time t with half time step 

  !$omp parallel default(none) private(ind,iskip,i,ivar,j,id,d,u,v,w,A,B,C,g, &
  !$omp&  dlnr1,dlnp0,dlnp1,e_mag,e_kin,e_prim,e_trunc,div,fact,d_old,d_new,e_old,e_new,e_min) &
  !$omp&  shared(son,dlnr,dlnp,gamma,ichem,chemistry,Grho,f,do_isothermal,smallr,smallp,not_iso, &
  !$omp&  hexp,beta_fix,ncoarse,ngridmax,uold,unew,poisson,active,temp_iso, &
  !$omp&  ilevel,pressure_fix,dtnew,enew,dx,divu,dlnr_max,smallc,scale_T2,smallt2)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     !$omp do
     do i=1,active(ilevel)%ngrid
        j=active(ilevel)%igrid(i)+iskip
        if (son(j)>0) cycle
        unew(j,1)=max(unew(j,1),smallr)
     end do
  end do
  if(poisson)then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        !$omp do
        do i=1,active(ilevel)%ngrid
           d=unew(active(ilevel)%igrid(i)+iskip,1)
           id = 1./d
           u=unew(active(ilevel)%igrid(i)+iskip,2)*id
           v=unew(active(ilevel)%igrid(i)+iskip,3)*id
           w=unew(active(ilevel)%igrid(i)+iskip,4)*id
           e_kin=0.5*d*(u**2+v**2+w**2)
           e_prim=unew(active(ilevel)%igrid(i)+iskip,5)-e_kin
           d_old=uold(active(ilevel)%igrid(i)+iskip,1)
           fact=d_old*id*0.5*dtnew(ilevel)*Grho
           u=u+f(active(ilevel)%igrid(i)+iskip,1)*fact
           unew(active(ilevel)%igrid(i)+iskip,2)=d*u
#if NDIM>1
           v=v+f(active(ilevel)%igrid(i)+iskip,2)*fact
           unew(active(ilevel)%igrid(i)+iskip,3)=d*v
#endif
#if NDIM>2
           w=w+f(active(ilevel)%igrid(i)+iskip,3)*fact
           unew(active(ilevel)%igrid(i)+iskip,4)=d*w
#endif
           e_kin=0.5*d*(u**2+v**2+w**2)
           unew(active(ilevel)%igrid(i)+iskip,5)=e_prim+e_kin
        end do
        !$omp enddo nowait
     end do
  end if

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     if (associated(active(ilevel)%igrid)) g => active(ilevel)%igrid
     if (do_isothermal) then        
        !$omp do
        do i=1,active(ilevel)%ngrid
           j = g(i)+iskip
           d=unew(j,1)
           u=unew(j,2)
           v=unew(j,3)
           w=unew(j,4)
           A=unew(j,6)+unew(j,nvar+1)
           B=unew(j,7)+unew(j,nvar+2)
           C=unew(j,8)+unew(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d
           e_mag=0.125_dp*(A*A+B*B+C*C)
           uold(j,5)=temp_iso*d+e_kin+e_mag                                     ! new total energy
        end do
        !$omp enddo nowait
        do ivar=1,nvar+3
           if (ivar==5) cycle                                                   ! total energy already updated
           !$omp do
           do i=1,active(ilevel)%ngrid
              uold(g(i)+iskip,ivar) = unew(g(i)+iskip,ivar)
           end do
           !$omp enddo nowait
        end do
     else
        dlnr1=dlnr(ilevel)
        dlnp1=dlnp(ilevel)
        if(pressure_fix)then
           !$omp do
           do i=1,active(ilevel)%ngrid
              j = g(i)+iskip
              if (son(j)>0) cycle
#ifdef KROME
              fact = (merge(unew(j,ichem),gamma,chemistry)-1.0_dp) ! Use cell dependent effective adiabatic index
#else
              fact = (gamma-1.0_dp)
#endif
              d=unew(j,1)
              u=unew(j,2)
              v=unew(j,3)
              w=unew(j,4)
              A=unew(j,6)+unew(j,nvar+1)
              B=unew(j,7)+unew(j,nvar+2)
              C=unew(j,8)+unew(j,nvar+3)
              e_kin=0.5_dp*(u*u+v*v+w*w)/d
              e_mag=0.125_dp*(A*A+B*B+C*C)
              e_new=unew(j,5)-e_kin-e_mag
              e_prim=enew(j)*(1.0d0+fact*divu(j))! Note: here divu=-div.u*dt
              div=abs(divu(j))*dx/dtnew(ilevel)
              ! Estimate of the local truncation errors
              e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
              if(e_new<e_trunc .and. divu(j)<0.0)then
                 unew(j,5)=e_prim+e_kin+e_mag
              end if
           end do
           !$omp enddo nowait
        end if
        fact = smallt2 / (scale_T2 * (gamma-1.0_dp))
        !$omp do reduction(+:not_iso)
        do i=1,active(ilevel)%ngrid
           j = g(i)+iskip
           if (son(j)>0) cycle
           d_old=uold(j,1)
           u=uold(j,2)
           v=uold(j,3)
           w=uold(j,4)
           A=uold(j,6)+uold(j,nvar+1)
           B=uold(j,7)+uold(j,nvar+2)
           C=uold(j,8)+uold(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d_old
           e_mag=0.125_dp*(A*A+B*B+C*C)
           e_old=uold(j,5)-e_kin-e_mag
           d_new=unew(j,1)
           u=unew(j,2)
           v=unew(j,3)
           w=unew(j,4)
           A=unew(j,6)+unew(j,nvar+1)
           B=unew(j,7)+unew(j,nvar+2)
           C=unew(j,8)+unew(j,nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d_new
           e_mag=0.125_dp*(A*A+B*B+C*C)
           e_new=unew(j,5)-e_kin-e_mag
           !e_min=d_new*smallc**2/(fact*(fact+1_dp))                             ! clamp at rho*smallc**2/(gamma*(gamma-1))
#ifdef KROME
           if (chemistry) fact = smallt2 / (scale_T2 * (unew(j,ichem)-1.0_dp))   ! Use cell dependent effective adiabatic index
#endif
           e_min=d_new * fact                                                    ! clamp at minimum temperature (T_CMB)
           if (e_new < e_min) then
              unew(j,5)=e_min+e_kin+e_mag
              not_iso=not_iso+1
           end if
           e_new=unew(j,5)-e_kin-e_mag
           dlnp0=abs(log(max(smallp/fact,e_new)/max(smallp/fact,e_old)))        ! e_new, e_old should never be zero, but just in case ...
           dlnp1=max(dlnp1,0.5_dp*dlnp0,-dlnp0)                                 ! allow twice as large increases compared to decreases
           dlnr1=max(dlnr1,abs(log(max(smallr,d_new)/max(smallr,d_old))))
        end do
        !$omp enddo nowait
        !$omp atomic
        dlnr(ilevel)=max(dlnr(ilevel),dlnr1)
        !$omp atomic
        dlnp(ilevel)=max(dlnp(ilevel),dlnp1)
        do ivar=1,nvar+3
           !$omp do
           do i=1,active(ilevel)%ngrid
              uold(g(i)+iskip,ivar) = unew(g(i)+iskip,ivar)
           end do
           !$omp enddo nowait
        end do     
     end if
  end do
  !$omp end parallel

  call max_real_mpi(dlnr(ilevel))                                               ! max over all MPI processes
  call max_real_mpi(dlnp(ilevel))                                               ! max over all MPI processes

  if (debug) then
     if (debug_level>4) then                  ! 4: output for each relevant rank and level
        if (not_iso>0) then
           if (do_isothermal) then
              write(*,'(a,i5,a,i3,a,i11)') 'set_uold: rank=',myid,' lvl=',ilevel,' not_iso =',not_iso
           else
              write(*,'(a,i5,a,i3,a,i11)') 'set_uold: rank=',myid,' lvl=',ilevel,' small e =',not_iso
           end if
        end if
     else if (debug_level>2) then             ! 3: output for each level
        call sum_int_mpi(not_iso)
        if (not_iso>0 .and. myid==1) then
           if (do_isothermal) then
              write(*,'(a,i3,a,i11)') 'set_uold: lvl=', ilevel, ' not_iso =', not_iso
           else
              write(*,'(a,i3,a,i11)') 'set_uold: lvl=', ilevel, ' small e =', not_iso
           end if
        end if
     else if (debug_level==2) then            ! 2: output for each Fine step
        not_iso_tot = not_iso_tot + not_iso
        if (ilevel==nlevelmax) then
           call sum_int8_mpi (not_iso_tot)
           if (not_iso_tot>0) then
              if (myid==1) then
              if (do_isothermal) then
                 write(*,'(a,i11)') 'set_uold: not_iso =', not_iso_tot
              else
                 write(*,'(a,i11)') 'set_uold: small e =', not_iso_tot
              end if
              end if
              not_iso_tot = 0
           end if
        end if
     else if (debug_level==1) then            ! 1: output for each coarse step (in adaptive_loop)
        not_iso_tot = not_iso_tot + not_iso
     end if
  endif

  call check_internal_energy(uold,ilevel,'set_uold2',1)

111 format('   Entering set_uold for level ',i2)

!  call trace_mpi('set_uold(ilevel)','exit',2)
end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------

  integer, parameter ::neul=5
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dflux_x,dflux_y,dflux_z
  real(dp)::dx,scale,oneontwotondim
  real(dp)::dflux,weight,tm
  real(dp),parameter::zero=0.0_dp
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         )::ind1,jcpu,indc
  logical ,dimension(1:nvector,0:twondim         )::problems
  real(dp),dimension(1:nvector,1:ndim            )::xgscr
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim)::g1
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::g2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gloc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,1:3,1:3,1:3)::emfx,emfy,emfz
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim)::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ok
  logical ,dimension(1:nvector,twotondim,twondim) :: mark_low_order      ! This is for marking cell interfaces
  logical ,dimension(1:nvector,1:3,1:3,1:3,3)     :: mark_low_order_edge ! This is for marking edges

  integer,dimension(1:nvector)::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  logical in_debug_region
  integer ind_gr
  real(dp):: d,u,v,w,A,B,C,e,e_kin,e_mag,e_limit

!  call trace_mpi('godfine1(ind_grid,ncache,ilevel)','enter',2)

  oneontwotondim = 1.d0/dble(twotondim)
  
  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  g1 = 0.0_dp; g2 = 0.0_dp; gloc=0.0_dp

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  call uf_debug (ind_grid,ncache,ilevel, 'god1')

  !dbg_grid = (ind_grid==29753 .or. ind_grid==32210) .and. local_debug
  !dbg_ind = ind_grid
  !if (any(dbg_grid)) print'(a,16i7)','call mag_unsplit with dbg_ind  =',dbg_ind

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  !call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          xgscr(nbuffer,:)=xg(ind_grid(i),:)
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     if (local_debug) then
        if (any(dbg_grid)) print'(a,3i3,2i4)', 'fathers: ',i1,j1,k1,nbuffer,nexist
     end if
     
     ! If not, interpolate variables from parent cells
     !call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
     do j=0,twondim
        do ivar=1,nvar+3
           do i=1,nbuffer
              u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
           end do
        end do
        do i=1,nbuffer
           ind1(i,j)=son(ibuffer_father(i,j))
           indc(i,j)=ibuffer_father(i,j)-ncoarse-((ibuffer_father(i,j)-ncoarse-1)/ngridmax)*ngridmax ! find grid id
           jcpu(i,j)=cpu_map(father(indc(i,j))) ! Find the owner cpu of the grid of the father cell
        end do
        if(poisson)then
           do idim=1,ndim
              do i=1,nbuffer
                 g1(i,j,idim)=f(ibuffer_father(i,j),idim)
              end do
           end do
        end if
     end do

     if (local_debug)then
        if (any(dbg_grid)) then
           do i=1,nbuffer
             print'(a,i3,1p,8e14.6)', 'u1',i,u1(i,:,2)
           end do
        end if
     end if

     problems(1:nbuffer,:) = .false.
     !call interpol_hydro(u1,g1,ind1,indc,u2,g2,nbuffer,xgscr,jcpu,'godfine1',problems)

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        if (local_debug) then
        if (i1==i1min+1 .and. j1==j1min+1 .and. k1==k1min+1) then
         do i=1,nexist
           ind_gr=igrid_nbor(ind_exist(i))
           if (in_debug_region(ind_gr)) then
              print 1, 'uold', &
                ind_cell(i),i2,j2,k2,0, &
                xg(ind_gr,1)+i2*dx, &
                xg(ind_gr,2)+j2*dx, &
                xg(ind_gr,3)+k2*dx, &
                uold(ind_cell(i),1), uold(ind_cell(i),3), uold(ind_cell(i),5)
            1 format(a,i8,4i3,3f10.6,1p,2(1x,3g14.6))
           endif
         end do
        endif
        endif
        
        ! Gather hydro variables
        if (nexist==16 .and. ncache==16) then ! we do not have to do indirect writing
           do ivar=1,nvar+3
              uloc(1,i3,j3,k3,ivar)=uold(ind_cell(1),ivar)
#if NVECTOR>15
              uloc(2,i3,j3,k3,ivar)=uold(ind_cell(2),ivar)
              uloc(3,i3,j3,k3,ivar)=uold(ind_cell(3),ivar)
              uloc(4,i3,j3,k3,ivar)=uold(ind_cell(4),ivar)
              uloc(5,i3,j3,k3,ivar)=uold(ind_cell(5),ivar)
              uloc(6,i3,j3,k3,ivar)=uold(ind_cell(6),ivar)
              uloc(7,i3,j3,k3,ivar)=uold(ind_cell(7),ivar)
              uloc(8,i3,j3,k3,ivar)=uold(ind_cell(8),ivar)
              uloc(9,i3,j3,k3,ivar)=uold(ind_cell(9),ivar)
              uloc(10,i3,j3,k3,ivar)=uold(ind_cell(10),ivar)
              uloc(11,i3,j3,k3,ivar)=uold(ind_cell(11),ivar)
              uloc(12,i3,j3,k3,ivar)=uold(ind_cell(12),ivar)
              uloc(13,i3,j3,k3,ivar)=uold(ind_cell(13),ivar)
              uloc(14,i3,j3,k3,ivar)=uold(ind_cell(14),ivar)
              uloc(15,i3,j3,k3,ivar)=uold(ind_cell(15),ivar)
              uloc(16,i3,j3,k3,ivar)=uold(ind_cell(16),ivar)
#endif
           end do
           if (local_debug) then
              do i=1,ncache
                 if (dbg_grid(i)) then
                    print'(a,i8,i4,2i3,i4,1p,e14.6)', 'gath16:', &
                     dbg_ind(i),i3,j3,k3,i,uloc(i,i3,j3,k3,2)
                end if
              end do
           end if
        else
           do ivar=1,nvar+3
              do i=1,nexist
                 uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              end do
              do i=1,nbuffer
                 uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              end do
              if (local_debug) then
                 do i=1,nexist
                    if (ivar==2 .and. dbg_grid(ind_exist(i))) then
                       print'(a,i8,i4,2i3,i4,1p,e14.6)', 'gather:', &
                        dbg_ind(i),i3,j3,k3,ind_exist(i),uloc(ind_exist(i),i3,j3,k3,ivar)
                    endif
                 end do
                 do i=1,nbuffer
                    if (ivar==2 .and. dbg_grid(ind_nexist(i))) then
                       print'(a,i8,i4,2i3,i4,2i3,1p,e14.6)', 'interp:', &
                        dbg_ind(i),i3,j3,k3,ind_nexist(i),i,ind_son,uloc(ind_nexist(i),i3,j3,k3,ivar)
                       print'(1p,8e14.6)',u2(i,:,ivar)
                    endif
                 end do
              end if
           end do
        endif        
        ! Gather gravitational acceleration
        if(poisson)then
           if (nexist==16 .and. ncache==16) then ! we do not have to do indirect writing
              do idim=1,ndim
                 gloc( 1,i3,j3,k3,idim)=f(ind_cell( 1),idim)
#if NVECTOR>15
                 gloc( 2,i3,j3,k3,idim)=f(ind_cell( 2),idim)
                 gloc( 3,i3,j3,k3,idim)=f(ind_cell( 3),idim)
                 gloc( 4,i3,j3,k3,idim)=f(ind_cell( 4),idim)
                 gloc( 5,i3,j3,k3,idim)=f(ind_cell( 5),idim)
                 gloc( 6,i3,j3,k3,idim)=f(ind_cell( 6),idim)
                 gloc( 7,i3,j3,k3,idim)=f(ind_cell( 7),idim)
                 gloc( 8,i3,j3,k3,idim)=f(ind_cell( 8),idim)
                 gloc( 9,i3,j3,k3,idim)=f(ind_cell( 9),idim)
                 gloc(10,i3,j3,k3,idim)=f(ind_cell(10),idim)
                 gloc(11,i3,j3,k3,idim)=f(ind_cell(11),idim)
                 gloc(12,i3,j3,k3,idim)=f(ind_cell(12),idim)
                 gloc(13,i3,j3,k3,idim)=f(ind_cell(13),idim)
                 gloc(14,i3,j3,k3,idim)=f(ind_cell(14),idim)
                 gloc(15,i3,j3,k3,idim)=f(ind_cell(15),idim)
                 gloc(16,i3,j3,k3,idim)=f(ind_cell(16),idim)
#endif
              end do
           else
              do idim=1,ndim
                 do i=1,nexist
                    gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
                 end do
                 do i=1,nbuffer
                    gloc(ind_nexist(i),i3,j3,k3,idim)=g2(i,ind_son,idim)
                 end do
              end do
           end if
        end if
        
        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  call uf_debug (ind_grid,ncache,ilevel, 'god2')

  ! Check if we have to default to LLF
  if (cbmax > 0) call cmp_fast_mode(mark_low_order,mark_low_order_edge,ind_grid,uloc,ncache,ilevel)

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call mag_unsplit(uloc,gloc,flux,emfx,emfy,emfz,tmp,mark_low_order,mark_low_order_edge,dx,dx,dx,dtnew(ilevel),ncache,ind_grid)

#ifdef CHECK_NAN
  idim=1
  do ivar=1,nvar
     do k3=MIN(1,ku1+2),MAX(1,ku2-2)
     do j3=MIN(1,ju1+2),MAX(1,ju2-2)
     do i3=if1,if2
        do i=1,ncache
           if (isnan(flux(i,i3,j3,k3,ivar,idim))) then
              print '(7i4,a,f10.5)', ilevel, i, i3, j3, k3, ivar, idim, 'flux ', flux(i,i3,j3,k3,ivar,idim)
           end if
        end do
     end do
     end do
     end do
  end do

  idim=2
  do ivar=1,nvar
     do k3=MIN(1,ku1+2),MAX(1,ku2-2)
     do j3=jf1,jf2
     do i3=MIN(1,iu1+2),MAX(1,iu2-2)
        do i=1,ncache
           if (isnan(flux(i,i3,j3,k3,ivar,idim))) then
              print '(7i4,a,f10.5)', ilevel, i, i3, j3, k3, ivar, idim, 'flux ', flux(i,i3,j3,k3,ivar,idim)
           end if
        end do
     end do
     end do
     end do
  end do


   do k3=kf1,ndim
   do j3=jf1,jf2
   do i3=if1,if2
      do i=1,ncache
         if (isnan(emfx(i,i3,j3,k3))) &
            print '(5i4,a,f10.5)', ilevel, i, i3, j3, k3, 'emfx ', emfx(i,i3,j3,k3)
         if (isnan(emfy(i,i3,j3,k3))) &
            print '(5i4,a,f10.5)', ilevel, i, i3, j3, k3, 'emfy ', emfy(i,i3,j3,k3)
         if (isnan(emfz(i,i3,j3,k3))) &
            print '(5i4,a,f10.5)', ilevel, i, i3, j3, k3, 'emfz ', emfz(i,i3,j3,k3)
      end do
   end do
   end do
   end do
#endif
  
  if(scheme.eq.'induction')then
  !---------------------------------
  ! Reset all Euler variables fluxes
  !---------------------------------
  do idim=1,ndim
  do ivar=1,neul
     do k3=k3min,k3max+1
     do j3=j3min,j3max+1
     do i3=i3min,i3max+1
        do i=1,ncache
           flux(i,i3,j3,k3,ivar,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
  end do

  endif

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              flux(i,i3,j3,k3,ivar,idim)=merge(flux(i,i3,j3,k3,ivar,idim),zero, &
                                  .not. (ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3)))
           end do
           if (local_debug.and.ivar==neul) then
              do i=1,ncache
                 if (in_debug_region(ind_grid(i))) &
                   print'(a,6i5,1p,e12.4,2l4)','reset1',i,i3,j3,k3,ivar,idim,flux(i,i3,j3,k3,ivar,idim), &
                   ok(i,i3-i0,j3-j0,k3-k0), ok(i,i3,j3,k3)
              end do
           end if
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              tmp (i,i3,j3,k3,ivar,idim)=merge(tmp(i,i3,j3,k3,ivar,idim),zero, &
                                  .not. (ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3)))
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bx
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,6,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#if NDIM>1
  !---------------------------------------------------------------
  ! Reset electromotive force along direction z at refined edges
  !---------------------------------------------------------------
  do k3=k3min,k3max
  do j3=1,3
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3  ,k3) .or. ok(i,i3  ,j3-1,k3) .or.  &
         & ok(i,i3-1,j3  ,k3) .or. ok(i,i3-1,j3-1,k3))then
           emfz(i,i3,j3,k3)=0.0d0
#if NDIM==2
           emfz(i,i3,j3,k3+1)=0.0d0
#endif
        end if
     end do
     if (local_debug) then
        do i=1,ncache
           if (in_debug_region(ind_grid(i))) &
             print'(a,4i5,1p,e12.4)','reset emfz',i,i3,j3,k3,emfz(i,i3,j3,k3)
        end do
     end if
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for By
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,7,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif
#if NDIM>2
  !---------------------------------------------------------------
  ! Reset electromotive force along direction y at refined edges
  !---------------------------------------------------------------
  do k3=1,3
  do j3=1,2
  do i3=1,3
     do i=1,ncache
        if(ok(i,i3  ,j3,k3  ) .or. ok(i,i3  ,j3,k3-1) .or.  &
         & ok(i,i3-1,j3,k3  ) .or. ok(i,i3-1,j3,k3-1))then
           emfy(i,i3,j3,k3)=0.0d0
        end if
     end do
     if (local_debug) then
        do i=1,ncache
           if (in_debug_region(ind_grid(i))) &
             print'(a,4i5,1p,e12.4)','reset emfy',i,i3,j3,k3,emfy(i,i3,j3,k3)
        end do
     end if
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset electromotive force along direction x at refined edges
  !---------------------------------------------------------------
  do k3=1,3
  do j3=1,3
  do i3=1,2
     do i=1,ncache
        if(ok(i,i3,j3  ,k3  ) .or. ok(i,i3,j3  ,k3-1) .or.  &
         & ok(i,i3,j3-1,k3  ) .or. ok(i,i3,j3-1,k3-1))then
           emfx(i,i3,j3,k3)=0.0d0
        end if
     end do
     if (local_debug) then
        do i=1,ncache
           if (in_debug_region(ind_grid(i))) &
             print'(a,4i5,1p,e12.4)','reset emfx',i,i3,j3,k3,emfx(i,i3,j3,k3)
        end do
     end if
  end do
  end do
  end do
  !---------------------------------------------------------------
  ! Reset Euler fluxes for Bz
  !---------------------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do i=1,ncache
           flux(i,i3,j3,k3,8,idim)=0.0d0
        end do
     end do
     end do
     end do
  end do
#endif

  call uf_debug (ind_grid,ncache,ilevel, 'god3')

  !-----------------------------------------------------
  ! Conservative update at level ilevel for Euler system
  !-----------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
!$omp atomic
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
           if (local_debug.and.ivar==neul) then
              do i=1,ncache
                 if (in_debug_region(ind_grid(i))) then
                   d=unew(ind_cell(i),1)
                   u=unew(ind_cell(i),2)
                   v=unew(ind_cell(i),3)
                   w=unew(ind_cell(i),4)
                   A=unew(ind_cell(i),6)+unew(ind_cell(i),nvar+1)
                   B=unew(ind_cell(i),7)+unew(ind_cell(i),nvar+2)
                   C=unew(ind_cell(i),8)+unew(ind_cell(i),nvar+3)
                   e_kin=0.5_dp*(u*u+v*v+w*w)/d
                   e_mag=0.125_dp*(A*A+B*B+C*C)
                   e=unew(ind_cell(i),5)-e_kin-e_mag
                   print'(a,i8,6i3,1p,7e14.6,1x,a)', 'god4', &
                   ind_cell(i),i,i3,j3,k3,ivar,idim, &
                   flux(i,i3,j3,k3,ivar,idim),flux(i,i3+i0,j3+j0,k3+k0,ivar,idim), &
                   flux(i,i3,j3,k3,ivar,idim)-flux(i,i3+i0,j3+j0,k3+k0,ivar,idim), &
                   unew(ind_cell(i),ivar),e_kin,e_mag,e,merge('NEG','   ',e<0)
                 end if
              end do
           end if
        end do
        do ivar=1,3
           do i=1,ncache
!$omp atomic
              unew(ind_cell(i),nvar+ivar)=unew(ind_cell(i),nvar+ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,neul+ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,neul+ivar,idim))
           end do
        end do

        if (local_debug) then
           do i=1,ncache
              if (in_debug_region(ind_grid(i))) then
                 print 1, 'flux', &
                   ind_cell(i),i3,j3,k3,idim, &
                   xg(ind_grid(i),1)+i2*dx, &
                   xg(ind_grid(i),2)+j2*dx, &
                   xg(ind_grid(i),3)+k2*dx, &
                   flux(i,i3   ,j3   ,k3   ,1,idim), &
                   flux(i,i3   ,j3   ,k3   ,3,idim), &
                   flux(i,i3   ,j3   ,k3   ,5,idim), &
                   flux(i,i3+i0,j3+j0,k3+k0,1,idim), &
                   flux(i,i3+i0,j3+j0,k3+k0,3,idim), &
                   flux(i,i3+i0,j3+j0,k3+k0,5,idim)
                 print 1, 'unew', &
                   ind_cell(i),i3,j3,k3,idim, &
                   xg(ind_grid(i),1)+i2*dx, &
                   xg(ind_grid(i),2)+j2*dx, &
                   xg(ind_grid(i),3)+k2*dx, &
                   unew(ind_cell(i),1) , unew(ind_cell(i),3), unew(ind_cell(i),5)
              endif
           end do
        endif

        if(pressure_fix)then
        ! Update velocity divergence
        do i=1,ncache
!$omp atomic
           divu(ind_cell(i))=divu(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
        end do
        ! Update internal energy
        do i=1,ncache
!$omp atomic
           enew(ind_cell(i))=enew(ind_cell(i))+ &
                & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
        end do
        end if
     end do
     end do
     end do
  end do

  call uf_debug (ind_grid,ncache,ilevel, 'god5')

  !---------------------------------------------------------
  ! Conservative update at level ilevel for induction system
  !---------------------------------------------------------
  do k3=k3min,k3max
  do j3=j3min,j3max
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bx using constraint transport
     do i=1,ncache
        dflux_x=( emfy(i,i3,j3,k3)-emfy(i,i3,j3,k3+1) ) &
          &    -( emfz(i,i3,j3,k3)-emfz(i,i3,j3+1,k3) )
        unew(ind_cell(i),neul+1)=unew(ind_cell(i),neul+1)+dflux_x
        dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
          &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )  
        unew(ind_cell(i),nvar+1)=unew(ind_cell(i),nvar+1)+dflux_x
     end do
  end do
  end do
  end do
#if NDIM>1  
  do k3=k3min,k3max
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update By using constraint transport
     do i=1,ncache
        dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
          &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
        unew(ind_cell(i),neul+2)=unew(ind_cell(i),neul+2)+dflux_y
        dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
          &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
        unew(ind_cell(i),nvar+2)=unew(ind_cell(i),nvar+2)+dflux_y
     end do
  end do
  end do
  end do
#endif
#if NDIM>2
  do k3=1,2
  do j3=1,2
  do i3=1,2
     ind_son=i3+2*(j3-1)+4*(k3-1)
     iskip=ncoarse+(ind_son-1)*ngridmax
     do i=1,ncache
        ind_cell(i)=iskip+ind_grid(i)
     end do
     ! Update Bz using constraint transport
     do i=1,ncache
        dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
          &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
        unew(ind_cell(i),neul+3)=unew(ind_cell(i),neul+3)+dflux_z
        dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
          &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
        unew(ind_cell(i),nvar+3)=unew(ind_cell(i),nvar+3)+dflux_z
     end do
  end do
  end do
  end do
#endif

  if(ilevel>levelmin)then

  !-----------------------------------------------------------
  ! Conservative update at level ilevel-1 for the Euler system
  !-----------------------------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
!$omp atomic
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
!$omp atomic
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,neul+ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
!$omp atomic
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
!$omp atomic
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
!$omp atomic
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     do ivar=1,3
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
!$omp atomic
              unew(ind_buffer(i),nvar+ivar)=unew(ind_buffer(i),nvar+ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,neul+ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
     ! Update velocity divergence
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
!$omp atomic
           divu(ind_buffer(i))=divu(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
!$omp atomic
           enew(ind_buffer(i))=enew(ind_buffer(i)) &
                & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
        end do
     end do
     end do
     end do
     end if

  end do
  ! End loop over dimensions

  call uf_debug (ind_grid,ncache,ilevel, 'god6')

#if NDIM>1  
  !---------------------------------------------------------------
  ! Conservative update at level ilevel-1 for the induction system
  !---------------------------------------------------------------
  i1=1; j1=0; k1=0
  if(ndim>1)j1=1
  if(ndim>2)k1=1

  !--------------------------------------
  ! Deal with 4 EMFz edges
  !--------------------------------------

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=0 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25*weight
!$omp atomic
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux
!$omp atomic
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
!$omp atomic
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux*0.5
!$omp atomic
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=0 and Y=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25*weight
!$omp atomic
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
!$omp atomic
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
!$omp atomic
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
!$omp atomic
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=1 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1+1)+9*(k1  )
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25*weight
!$omp atomic
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
!$omp atomic
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
!$omp atomic
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux*0.5
!$omp atomic
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and By using fine EMFz on X=1 and Y=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1-1)+9*(k1  )
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25*weight
!$omp atomic
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
!$omp atomic
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
!$omp atomic
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
!$omp atomic
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux*0.5
     endif
  end do
#endif
  !--------------------------------------
  ! Deal with 4 EMFx edges
  !--------------------------------------
#if NDIM>2
  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25*weight
!$omp atomic
     unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)+dflux
!$omp atomic
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
!$omp atomic
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux*0.5
!$omp atomic
        unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=0 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1-1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25*weight
!$omp atomic
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
!$omp atomic
     unew(ind_buffer2,2+nvar)=unew(ind_buffer2,2+nvar)-dflux
!$omp atomic
     unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
!$omp atomic
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25*weight
!$omp atomic
     unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)-dflux
!$omp atomic
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
!$omp atomic
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
!$omp atomic
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux*0.5
     endif
  end do

  ! Update coarse By and Bz using fine EMFx on Y=1 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
  ind_father2=1+(i1  )+3*(j1+1)+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25*weight
!$omp atomic
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux
!$omp atomic
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,2+neul)=unew(ind_buffer2,2+neul)+dflux
!$omp atomic
     unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux*0.5
!$omp atomic
        unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
     endif
  end do

  !--------------------------------------
  ! Deal with 4 EMFy edges
  !--------------------------------------

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=0 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25*weight
!$omp atomic
     unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux
!$omp atomic
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)+dflux
!$omp atomic
     unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux*0.5
!$omp atomic
        unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=1 grid edge
  ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1-1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25*weight
!$omp atomic
     unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)+dflux
!$omp atomic
     unew(ind_buffer2,1+nvar)=unew(ind_buffer2,1+nvar)+dflux
!$omp atomic
     unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux*0.5
!$omp atomic
        unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=1 grid edge
  ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1+1)
  ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25*weight
!$omp atomic
     unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux
!$omp atomic
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)+dflux
!$omp atomic
     unew(ind_buffer2,3+neul)=unew(ind_buffer2,3+neul)-dflux
!$omp atomic
     unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux*0.5
!$omp atomic
        unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux*0.5
     endif
  end do

  ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=0 grid edge
  ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
  ind_father2=1+(i1+1)+3*(j1  )+9*(k1-1)
  ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
  do i=1,ncache
     ind_buffer1=nbors_father_cells(i,ind_father1)
     ind_buffer2=nbors_father_cells(i,ind_father2)
     ind_buffer3=nbors_father_cells(i,ind_father3)
     weight=1.0
     if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
     if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
     dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25*weight
!$omp atomic
     unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux
!$omp atomic
     unew(ind_buffer2,3+nvar)=unew(ind_buffer2,3+nvar)-dflux
!$omp atomic
     unew(ind_buffer2,1+neul)=unew(ind_buffer2,1+neul)-dflux
!$omp atomic
     unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux
     if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
!$omp atomic
        unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux*0.5
!$omp atomic
        unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux*0.5
     endif
  end do
#endif

  end if  !  (ilevel>levelmin)

  if (entropy_fix) then
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        do i=1,ncache
           d=unew(ind_cell(i),1)
           u=unew(ind_cell(i),2)
           v=unew(ind_cell(i),3)
           w=unew(ind_cell(i),4)
           A=unew(ind_cell(i),6)+unew(ind_cell(i),nvar+1)
           B=unew(ind_cell(i),7)+unew(ind_cell(i),nvar+2)
           C=unew(ind_cell(i),8)+unew(ind_cell(i),nvar+3)
           e_kin=0.5_dp*(u*u+v*v+w*w)/d
           e_mag=0.125_dp*(A*A+B*B+C*C)
           e=unew(ind_cell(i),5)-e_kin-e_mag
           if (e < entropy_elimit*unew(ind_cell(i),5)) then
              e_limit=exp(uloc(i,i3,j3,k3,5)+gamma*log(d))/(gamma-1.0)
              if (in_debug_region(ind_grid(i))) &
                 print'(a,i8,i3,1p,10e11.3)','MKF: S,etot,eint,elimit =', &
                   ind_grid(i),i2+j2*2+k2*4, &
                   uloc(i,i3,j3,k3,5), &
                   uloc(i,i2,j2,k2,5), &
                   unew(ind_cell(i),5),e,e_limit
              n_entropy(1)=n_entropy(1)+1
              if (e < e_limit) then
                 unew(ind_cell(i),5)=e_limit+e_kin+e_mag
                 n_entropy(2)=n_entropy(2)+1
              end if
           end if
        end do
     end do
     end do
     end do
  end if

  call uf_debug (ind_grid,ncache,ilevel, 'god7')
!  call trace_mpi('godfine1(ind_grid,ncache,ilevel)','exit',2)
end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_fast_mode(mark_low_order,mark_low_order_edge,ind_grid,uloc,ngrid,ilevel)
  use amr_parameters,   only : dp,nvector,ngridmax,ndim,twotondim,twondim,verbose,chemistry
  use amr_commons,      only : ncoarse,xg,myid
  use hydro_parameters, only : gamma,nvar,iu1,iu2,ju1,ju2,ku1,ku2,h_verbose,ichem
  use hydro_commons,    only : uold,smallr,smallc,cbmax,n_switch,n_cbmax
  use const,            only : zero,half,one
  implicit none
  integer :: ngrid
  logical ,dimension(1:nvector,twotondim,twondim)                :: mark_low_order      ! This is for marking cell interfaces
  logical ,dimension(1:nvector,1:3,1:3,1:3,3)                    :: mark_low_order_edge ! This is for marking edges
  integer, dimension(1:nvector)                                  :: ind_grid
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uloc
  real(dp), save :: ctot_max=0.0
  integer                                                        :: ilevel
  !
  integer :: ind, i, iskip, icell, idim, ix, iy, iz, lcx, ucx, lx, ux, ly, uy, lz, uz
  real(dp):: rho, B2, P, cs2, cc2, ctot,smallp
  real(dp), dimension(ndim) :: v, Bc
  character(len=1), dimension(3) :: dir=(/'x','y','z'/)
#ifdef KROME
  real(dp):: eff_gamma
#endif
  integer, save :: nmessage=100
  !
!  call trace_mpi('cmp_fast_mode','enter',2)

  smallp = smallr*smallc**2/gamma
  mark_low_order = .false.
  mark_low_order_edge = .false.
  !
  ! Check if the fast speed in the cell is above the limit. Mark the corresponding
  ! cell interfaces. If it is an interior cell all 6 interfaces are marked
  !
  ! For edges we directly mark an edge if any of the four adjacent cells have cfast > cbmax
  !
  do iz=max(0,ku1),min(3,ku2)
     do iy= max(0,ju1), min(3,ju2)
        lx=max(0,iu1); ux=min(3,iu2)
        do ix=max(0,iu1), min(3,iu2)
           if (corner(ix,iy,iz)) cycle ! exclude corner cells; they do not mark any edge or cell interface
           do i=1,ngrid
              rho = max(uloc(i,ix,iy,iz,1),smallr)
              B2 = zero
              P = uloc(i,ix,iy,iz,5)
              do idim=1,ndim
                 v(idim) = uloc(i,ix,iy,iz,1+idim) / rho
                 Bc(idim) = half*(uloc(i,ix,iy,iz,5+idim)+uloc(i,ix,iy,iz,nvar+idim))
                 B2 = B2 + Bc(idim)*Bc(idim)
                 P = P - half*rho*v(idim)*v(idim)
              enddo
              P = P - half*B2
#ifdef KROME
              eff_gamma = merge(uloc(i,ix,iy,iz,ichem),gamma,chemistry)
              P = max((eff_gamma-one)*P,smallp*gamma/eff_gamma)
              cs2 = eff_gamma*P/rho
#else
              P = max((gamma-one)*P,smallp)
              cs2 = gamma*P/rho
#endif
              cc2 = half*(B2/rho+cs2)
              ctot =        abs(v(1)) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*Bc(1)*Bc(1)/rho,zero)))
#if NDIM>1
              ctot = ctot + abs(v(2)) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*Bc(2)*Bc(2)/rho,zero)))
#endif
#if NDIM>2
              ctot = ctot + abs(v(3)) + sqrt(cc2+sqrt(max(cc2*cc2-cs2*Bc(3)*Bc(3)/rho,zero)))
#endif
              if (h_verbose > 2) then
                if (ctot > ctot_max) print *, 'ctot :', myid, ctot
                ctot_max = max(ctot,ctot_max)
              endif
              if (ctot >= cbmax) then ! we are above, mark in arrays
                 !$omp atomic
                 n_cbmax = n_cbmax+1
                 ! Mark edges. Between 1 and 12 edges are marked
                 lx = max(ix,1); ux = min(ix+1,3)
                 ly = max(iy,1); uy = min(iy+1,3)
                 lz = max(iz,1); uz = min(iz+1,3)
                 if ((ix-1)*(2-ix)>=0) mark_low_order_edge(i,   ix,ly:uy,lz:uz,1) = .true. ! mark edges corresponding to emfx
                 if ((iy-1)*(2-iy)>=0) mark_low_order_edge(i,lx:ux,iy,   lz:uz,2) = .true. ! mark edges corresponding to emfy
                 if ((iz-1)*(2-iz)>=0) mark_low_order_edge(i,lx:ux,ly:uy,iz,   3) = .true. ! mark edges corresponding to emfz

                 ! cell interfaces; only if not on an edge of the 4x4x4 cube
                 if (.not. corner_edge(ix,iy,iz)) then ! mark cell interfaces 
                    ind = ix   + 2*(iy-1 + 2*(iz-1)); if (inside(ix,iy,iz)) mark_low_order(i,ind,:) = .true.  ! mark all interfaces of cell
                    if (verbose) print '(a,3i2,a,i10,a,6L2)', 'Marking cell with coords :', ix, iy, iz, ' oct coord :', ind_grid(i), &
                                   ' Inside: ', inside(ix-1,iy,iz), inside(ix+1,iy,iz), inside(ix,iy-1,iz), &
                                                inside(ix,iy+1,iz), inside(ix,iy,iz-1), inside(ix,iy,iz+1)
                    !
                    ind = ix-1 + 2*(iy-1 + 2*(iz-1)); if (inside(ix-1,iy,iz)) mark_low_order(i,ind,2) = .true. ! mark interface of cell -1 in x-dir
                    ind = ix+1 + 2*(iy-1 + 2*(iz-1)); if (inside(ix+1,iy,iz)) mark_low_order(i,ind,1) = .true. ! mark interface of cell +1 in x-dir
#if NDIM>1
                    ind = ix + 2*(iy-1-1 + 2*(iz-1)); if (inside(ix,iy-1,iz)) mark_low_order(i,ind,4) = .true. ! mark interface of cell -1 in y-dir
                    ind = ix + 2*(iy+1-1 + 2*(iz-1)); if (inside(ix,iy+1,iz)) mark_low_order(i,ind,3) = .true. ! mark interface of cell +1 in y-dir
#endif
#if NDIM>2
                    ind = ix + 2*(iy-1 + 2*(iz-1-1)); if (inside(ix,iy,iz-1)) mark_low_order(i,ind,6) = .true. ! mark interface of cell -1 in z-dir
                    ind = ix + 2*(iy-1 + 2*(iz+1-1)); if (inside(ix,iy,iz+1)) mark_low_order(i,ind,5) = .true. ! mark interface of cell +1 in z-dir
#endif
                 endif
              endif 
           enddo
        enddo
     enddo
  enddo
  !
  if (h_verbose >  2) then
     do i=1,ngrid
        do ind=1,twotondim
           if (any(mark_low_order(i,ind,:))) write(*,'(a,6L2,a,i2,a,i10,a,i5,a,3f10.5,a,i3)') &
              ' Switched to low order llf solver for ',mark_low_order(i,ind,:), &
              ' interfaces in cell ',ind,' in the oct ', ind_grid(i), ' on thread ', myid, &
              ' with position ', xg(ind_grid(i),1:ndim), ' on level ',ilevel
        enddo
        do ix=1,3 ! direction
           if (any(mark_low_order_edge(i,:,:,:,ix))) write(*,'(a,27L2,a,i10,a,i5,a,3f10.5,a,i3)') &
              ' Switched to low order llf solver for ',mark_low_order_edge(i,:,:,:,ix), &
              ' edges along the '//dir(ix)//'-direction in oct ', ind_grid(i), ' on thread ', myid, &
              ' with position ', xg(ind_grid(i),1:ndim), ' on level ',ilevel
        enddo
     enddo
  else if (h_verbose > 1) then
     do i=1,ngrid
        if ((any(mark_low_order(i,:,:)) .or. any(mark_low_order_edge(i,:,:,:,:))) .and. nmessage > 0) then
           write(*,'(2(a,i3),a,i10,a,i5,a,3f10.5,a,i3)') &
             ' Switched to llf solver for ', sum(merge(1,0,mark_low_order(i,:,:))), &
             ' interfaces and ',sum(merge(1,0,mark_low_order_edge(i,:,:,:,:))), &
             ' edges in the oct ', ind_grid(i), ' on thread ', myid, &
             ' with position ', xg(ind_grid(i),1:ndim), ' on level ',ilevel
           nmessage = nmessage - 1
        endif
     enddo
  else if (h_verbose > 0) then
     do i=1,ngrid
        if (any(mark_low_order(i,:,:))) then
           !$omp atomic
           n_switch=n_switch+1
        endif
     enddo
  endif
!  call trace_mpi('cmp_fast_mode','exit',2)
contains
  logical function corner_edge(ix, iy, iz) ! Are we on any of the eight 2-cell edges in the 4x4x4 cube
     integer :: ix, iy, iz
     corner_edge = ((ix-1)*(2-ix) < 0 .and. (iy-1)*(2-iy) < 0) .or. &
                   ((ix-1)*(2-ix) < 0 .and. (iz-1)*(2-iz) < 0) .or. &
                   ((iy-1)*(2-iy) < 0 .and. (iz-1)*(2-iz) < 0)
  end function corner_edge
  !
  logical function corner(ix, iy, iz) ! Are we in any of the 8 corner cells in the 4x4x4 cube
     integer :: ix, iy, iz
     corner = ((ix-1)*(2-ix) < 0 .and. (iy-1)*(2-iy) < 0 .and. (iz-1)*(2-iz) < 0)
  end function corner
  !
  logical function inside(ix, iy, iz) ! Are we inside the central 2x2x2 cube ?
     integer :: ix, iy, iz
     inside = (ix-1)*(2-ix) >= 0 .and. (iy-1)*(2-iy) >= 0 .and. (iz-1)*(2-iz) >= 0
  end function inside
end subroutine cmp_fast_mode
