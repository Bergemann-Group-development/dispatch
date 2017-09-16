PROGRAM tests
  USE omp_lib
  implicit real*8 (a-h,o-z)
  implicit integer*8 (i-n)
  real(8), dimension(2):: r, r0, r1, r2, r3, v, a
  real(8):: dt, t
  real, pointer:: rr(:,:)
  real, parameter:: pi=4.*atan(1.)
  character(len=8):: s1, s2
  !-----------------------------------------------------------------------------
  ! The eccentricity e is defined by r_min/r_max=(1-e)/(1+e).  In terms of speed
  ! this translates into v_min/v_max=(1-e)/(1+e).  The speed reduction at apogee
  ! f can be deduced from energy conservation.
  !-----------------------------------------------------------------------------
  e=0.5
  ratio=(1.0-e)/(1.0+e)
  period=2.*pi*((1+ratio)/2.)**1.5
  f=sqrt(2.*ratio/(1+ratio))
  print '("e =",f4.1,"  v =",f6.3,"  period =",f6.3)', e, f,period
  nstepmax=2e11
  noutput=5000
  n_orbit=1e8
  do k=1,1
    if (k==1) then
      cdt = 0.01
    else
      cdt = 0.005
    end if
    n_period=ceiling(period/cdt)
    n_steps=n_orbit*n_period
    n_output=n_steps/noutput
    allocate (rr(n_period+1,2))
    write (30+k*10+m) n_period+1
    do m=2,2
      r=[0,1]
      v=[1,0]*f
      e0=0.5*sum(v**2)-1./sqrt(sum(r**2))
      t=0.0
      nstep=0
      dtmin=1e9
      dtmax=0.
      rmin=10.
      rmax=0.
      emin=0.
      emax=-10.
      wc = omp_get_wtime()
      iorbit=0
      i=0
      r3 = r; r2 = r; r1 = r
      if (k==0.and.m==3) exit
      do while (t < period*n_orbit)
        if (k==1) then
          r0=1.875*r1-1.25*r2+0.375*r3
          ra=sqrt(sum(r0**2))
          rmin=min(rmin,ra)
          rmax=max(rmax,ra)
          dt=cdt*sqrt(ra)**3
          dtmin=min(dtmin,dt)
          dtmax=max(dtmax,dt)
        else
          !dt=cdt*ratio**2/f
          dt=cdt
        end if
        if      (m==1 .or. (m==3.and.mod(nstep,2)==1)) then
          ! DKD
          r=r+0.5*dt*v
          a=-r/sqrt(sum(r**2))**3
          v=v+dt*a
          r=r+0.5*dt*v
        else if (m==2 .or. (m==3.and.mod(nstep,2)==0)) then
          ! KDK
          a=-r/sqrt(sum(r**2))**3
          v=v+0.5*dt*a
          r=r+dt*v
          a=-r/sqrt(sum(r**2))**3
          v=v+0.5*dt*a
        end if
        r3=r2; r2=r1; r1=r
        t=t+dt
        i=i+1
        if (i<=n_period+1) rr(i,:)=r
        nstep=nstep+1
        e1=0.5*sum(v**2)-1./sqrt(sum(r**2))
        emin=min(emin,e1)
        emax=max(emax,e1)
        if (mod(nstep,n_output)==1) write (10+k*10+m) t/period,dt,e1,emin,emax
        if (nstep>nstepmax) exit
        if (r1(2)>0.0 .and. r2(2)<0.0) then
          if (mod(iorbit,1000000)==0) then
            print '("orbit",1pe9.1)', real(iorbit)
            write (30+k*10+m) rr
          end if
          iorbit=iorbit+1
          !print *,iorbit
          i=0
        end if
        !print 2,i,t,dt,ra,cdt,r0,r,v,a,e1
        !2 format(i5,7(2f11.6,1x))
      end do
      !print *,iorbit
      e1=0.5*sum(v**2)-1./sqrt(sum(r**2))
      n_orbit=t/period
      if (k==1) then
        s2='adaptive'
      else
        s2=''
      end if
      if (m==1) then
        s1='    DKD:'
      else if (m==2) then
        s1='    KDK:'
      else
        s1='DKD+KDK:'
      end if
      !print *,dtmin,dtmax,ratio
      !print *,(rmin+rmax)*0.5,pi*(rmin+rmax),period,period*iorbit/t
      de=max(abs(emax-e0),abs(emin-e0))
      print 1, s2,s1,real(nstep)/n_orbit,'steps/orb',real(n_orbit),'orbits',real(nstep),'steps', &
        omp_get_wtime()-wc,'s',e1,'energy',de/e0,'error',(1.-e1/e0)/t,'trend',t,'time'
    1 format(a8,1x,a,1p,10(e11.2,1x,a,1x))
    end do
    deallocate (rr)
  end do
END
