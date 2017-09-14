PROGRAM tests
  USE omp_lib
  implicit real*8 (a-h,o-z)
  real(8), dimension(2):: r, r0, r1, r2, r3, v, a
  real(8):: dt, t
  real, parameter:: pi=4.*atan(1.)
  integer:: per_orbits(5)
  integer(8):: nstep, nstepmax

  !-----------------------------------------------------------------------------
  ! The eccentricity e is defined by r_min/r_max=(1-e)/(1+e).  In terms of speed
  ! this translates into v_min/v_max=(1-e)/(1+e).  When e=0 the two speeds are
  ! equal -- here =1 -- with energy = -1/2.   If the speed is reduced to zero
  ! the energy is reduced to -1.  The energy when the speed is changed with a
  ! factor is -1 + f^2/2
  !-----------------------------------------------------------------------------
  e=0.9
  v_ratio=(1.0-e)/(1.0+e)
  f=sqrt(2.*v_ratio/(1+v_ratio))
  print '("e =",f4.1,"  v =",f6.3)', e, f
  per_orbits=[50,100,200,300,500]
  nstepmax=2e8
  do j=1,5
    per_orbit=per_orbits(j)
    cdt=8.*atan(1.)/per_orbit
    n=n_orbit*per_orbit
    do k=0,1
      n_orbit=10000
      r=[0,1]
      r1 = r; r2 = r; r3 = r
      v=[1,0]*f
      e0=0.5*sum(v**2)-1./sqrt(sum(r**2))
      t=0.0
      nstep=0
      wc = omp_get_wtime()
      do while (t < 2.*pi*n_orbit)
        if (k==1) then
          r0=1.875*r1-1.25*r2+0.375*r3
          dt=cdt*sum(r0**2)**(3./4.)
        else
          dt=cdt*v_ratio**2/f
        end if
        ! DKD
        r=r+0.5*dt*v
        a=-r/sqrt(sum(r**2))**3
        v=v+dt*a
        r=r+0.5*dt*v
        r3=r2; r2=r1; r1=r
        t=t+dt
        nstep=nstep+1
        if (nstep>nstepmax) exit
      end do
      e1=0.5*sum(v**2)-1./sqrt(sum(r**2))
      print *,'    DKD:',e1,(1.-e1/e0),nstep

      r=[0,1]
      v=[1,0]*f
      e0=0.5*sum(v**2)-1./sqrt(sum(r**2))
      t=0.0
      dt=cdt
      nstep=0
      do while (t < 2.*pi*n_orbit)
        if (k==1) then
          r0=1.875*r1-1.25*r2+0.375*r3
          dt=cdt*sum(r0**2)**(3./4.)
        else
          dt=cdt*v_ratio**2/f
        end if
        ! KDK
        a=-r/sqrt(sum(r**2))**3
        v=v+0.5*dt*a
        r=r+dt*v
        a=-r/sqrt(sum(r**2))**3
        v=v+0.5*dt*a
        r3=r2; r2=r1; r1=r
        t=t+dt
        nstep=nstep+1
        if (nstep>nstepmax) exit
      end do
      e1=0.5*sum(v**2)-1./sqrt(sum(r**2))
      print *,'    KDK:',e1,(1.-e1/e0),nstep

      r=[0,1]
      v=[1,0]*f
      e0=0.5*sum(v**2)-1./sqrt(sum(r**2))
      dt=cdt
      t=0.0
      nstep=0
      do while (t < 2.*pi*n_orbit)
        if (k==1) then
          r0=1.875*r1-1.25*r2+0.375*r3
          dt=cdt*sum(r0**2)**(3./4.)
        else
          dt=cdt*v_ratio**2/f
        end if
        if (mod(nstep,2)==1) then
          ! DKD
          r=r+0.5*dt*v
          a=-r/sqrt(sum(r**2))**3
          v=v+dt*a
          r=r+0.5*dt*v
        else
          ! KDK
          a=-r/sqrt(sum(r**2))**3
          v=v+0.5*dt*a
          r=r+dt*v
          a=-r/sqrt(sum(r**2))**3
          v=v+0.5*dt*a
        end if
        r3=r2; r2=r1; r1=r
        t=t+dt
        nstep=nstep+1
        if (nstep>nstepmax) exit
      end do
      n_orbit=t/(2.*pi)
      e1=0.5*sum(v**2)-1./sqrt(sum(r**2))
      print *,'DKD+KDK:',e1,(1.-e1/e0),nstep,nstepmax
      print '(1p,2(e11.2,a,5x),e11.2,a)', real(nstep)/n_orbit,' steps/orb',real(n_orbit),' orbits', omp_get_wtime()-wc,' seconds'
    end do
  end do
END
