!> =============================================================================
!> Setup the driver conditions for the BP experiment
!> =============================================================================
MODULE driver_mod
  USE io_mod
  USE mesh_mod
  USE scalar_mod
  USE mhd_mod
  USE bndry_mod
  USE scaling_mod
  USE stagger_mod

CONTAINS

!> =============================================================================
!> The cylinder driver on the lower x boundary.
!> =============================================================================
SUBROUTINE cylinder_driver (self, m)
      implicit none
      class(mhd_t)   :: self
      class(mesh_t), pointer:: m(:)
      
      real, dimension(:,:,:), pointer :: Ux, Uy, Uz
      integer        :: i, nxgrid, nygrid, nzgrid, ix, iy, iz, iostat
      integer, save  :: ifirst=1
      integer, save  :: nsmooth=10
      real, save     :: ampl=1., tramp=1., t0=0.
      real, save     :: vv, dist, angle
      real, save     :: centery=0., centerz=0.
      real, save     :: direction=0., radius=.1
      real           :: y,z
      real, allocatable, dimension(:,:)   :: t1, vr, vt, th
      real, allocatable, dimension(:,:)   :: y2, z2
      real, allocatable, dimension(:,:,:) :: t2, t3
      real, pointer    , dimension(:,:)   :: uy_xl, uz_xl
      real, save :: pi=2.*asin(1.)
      
      namelist /cylinder_params/ ampl, centery, centerz, direction, radius, nsmooth, tramp, t0
!
!-----------------------------------------------------------------------
!
!  driver where a simple spot is moving along a straight line for a period of 
!  time allowing only to move the strong flux source while everything else 
!  is swept around it as in an incompressible cylinder flow.

     Ux => self%U(:,:,:,1)
     Uy => self%U(:,:,:,2)
     Uz => self%U(:,:,:,3)

!$OMP CRITICAL (io_crit)
      if (ifirst.gt.0) then
         rewind (io%input)
         read (io%input,cylinder_params, iostat=iostat)
         if (io%master) write (*,cylinder_params)
         ifirst=0
      endif
!$OMP END CRITICAL (io_crit)

! grid size of variables
      nxgrid=size(self%mem, 1)
      nygrid=size(self%mem, 2)
      nzgrid=size(self%mem, 3)
      
! point boundary variables (slices) well into the 9th memory block
      uy_xl => self%mem(:,:,5,9,self%it,1)
      uz_xl => self%mem(:,:,6,9,self%it,1)

! move the velocity pattern in time if lower boundary
      if (int(m(1)%d) - int(m(1)%lb) == 0) then
         
         allocate(t1(nygrid,nzgrid), vr(nygrid,nzgrid), vt(nygrid,nzgrid), th(nygrid,nzgrid))
         allocate(y2(nygrid,nzgrid), z2(nygrid,nzgrid))
         
         vv = ampl*0.5*(1 + tanh((self%time - t0)/tramp - 3.) )
         dist=ampl*0.5*(-t0 + self%time - tramp*dlog(dcosh(3.0d+0)/dcosh((t0 - self%time + 3.0d+0*tramp)/tramp))) 
         angle=direction/180.*pi
         
         do iy=1,nygrid
            y = m(2)%p + m(2)%r(iy)
           do iz=1,nzgrid
             z = m(3)%p + m(3)%r(iz)
             y2(iy,iz)=y-centery-dist*cos(angle)
             z2(iy,iz)=z-centerz-dist*sin(angle)
           enddo
         enddo

!!!  need to know the full size of the domain to see that the pattern is all the time centered on a virtual box!!
         where (y2 < -self%box(2))
            y2=y2+2*self%box(2)
         end where
         where (y2 > self%box(2))
            y2=y2-2*self%box(2)
         end where
         where (z2 < -self%box(3))
            z2=z2+2*self%box(3)
         end where
         where (z2 > self%box(3))
            z2=z2-2*self%box(3)
         end where
         t1=sqrt(y2**2+z2**2)

! angle of all points relative the new center of driver point.
         th = atan2(z2,y2)

! flow profile, returned to coordinate system where the cylinder is moving through fluid in rest
         vr=0
         vt=0
         where (t1 > radius)
            vr= vv*(1-radius**2/t1**2)*cos(pi+th-angle)
            vt=-vv*(1+radius**2/t1**2)*sin(pi+th-angle)
         endwhere

         uy_xl=cos(th)*vr-sin(th)*vt+vv*cos(angle)
         uz_xl=sin(th)*vr+cos(th)*vt+vv*sin(angle)
         
! smooth the profile in both y and z plane allowing a smooth transition
! to zero speed around the flow pattern and a smooth change in flow speed
! close to the cyliner area.
         do i=1,nsmooth
           call smooth2dyz(uy_xl, nygrid, nzgrid,m)
           call smooth2dyz(uz_xl, nygrid, nzgrid,m)
           ! UPDATES GHOSTZONES in all directions -- overkill but simple
           call self%dnload
         enddo

         allocate(t2(nxgrid, nygrid, nzgrid), t3(nxgrid, nygrid, nzgrid))
         t2=0
         t2(1,:,:)=uy_xl
         t3=ydn(t2)
         uy_xl=t3(1,:,:)
         t2(1,:,:)=uz_xl
         t2=zdn(t2)
         uz_xl=t3(1,:,:)
         deallocate(vr, vt, th, y2, z2, t1, t2, t3)

! assign the data to the real velocity arrays
         Ux(m(1)%li,:,:) = -Ux(m(1)%li+1,:,:)
         Uy(m(1)%li,:,:) = uy_xl
         Uz(m(1)%li,:,:) = uz_xl
         Uy(m(1)%li+1,:,:) = (4./9.)*Uy(m(1)%li,:,:) + (6./9.)*Uy(m(1)%li+2,:,:) - (1./9.)*Uy(m(1)%li+2,:,:)
         Uz(m(1)%li+1,:,:) = (4./9.)*Uz(m(1)%li,:,:) + (6./9.)*Uz(m(1)%li+2,:,:) - (1./9.)*Uz(m(1)%li+2,:,:)
         deallocate (uy_xl, uz_xl)
      endif

! assumes zero flow speed at the top boundary
      if (int(m(1)%d) - int(m(1)%ub) == 0) then
         Ux(m(1)%ui,:,:) = -Ux(m(1)%ui-1,:,:)
         Uy(m(1)%ui-1,:,:) = 0.
         Uz(m(1)%ui-1,:,:) = 0.
         Uy(m(1)%ui-1,:,:) = (4./9.)*Uy(m(1)%ui,:,:) + (6./9.)*Uy(m(1)%ui-2,:,:) - (1./9.)*Uy(m(1)%ui-3,:,:)
         Uz(m(1)%ui-1,:,:) = (4./9.)*Uz(m(1)%ui,:,:) + (6./9.)*Uz(m(1)%ui-2,:,:) - (1./9.)*Uz(m(1)%ui-3,:,:)
      endif

END SUBROUTINE cylinder_driver

!>----------------------------------------------------

SUBROUTINE smooth2dyz(f, nygrid, nzgrid, m)

!>----------------------------------------------------
!> do smoothing in the two horizontal directions. Assume the domain is periodic with in the patch!!!!
!> only for the interal points, assume bdry are handled else where, which is not good for the case
!> many patches... need to update in the loop above...

      implicit none
      class(mesh_t), pointer             :: m(:)
      integer                            :: nygrid, nzgrid
      real, dimension(:,:)               :: f
      real, allocatable, dimension(:,:)  :: t1, t2

      allocate (t1(nygrid, nzgrid), t2(nygrid, nzgrid))

      t1(m(2)%li:m(2)%ui,:) = f(m(2)%li+1:m(2)%ui+1,:)
      t2(m(2)%li:m(2)%ui,:) = f(m(2)%li-1:m(2)%ui-1,:)
      f = .5*f + .25*(t1 + t2)

      t1(:,m(3)%li:m(3)%ui) = f(:,m(3)%li+1:m(3)%ui+1)
      t2(:,m(3)%li:m(3)%ui) = f(:,m(3)%li-1:m(3)%ui-1)
      f = .5*f + .25*(t1 + t2)

      deallocate (t1, t2)

END SUBROUTINE smooth2dyz

!>----------------------------------------------------

SUBROUTINE damping(self, dfdt, f, fbound, scale, m)
      implicit none
      class(mhd_t)                           :: self
      class(mesh_t), pointer                 :: m(:)
      real, dimension(:,:,:)                 :: dfdt, f, fbound
      real, allocatable, dimension(:,:,:)    :: dissipation
      integer                                :: scale, iostat
      real, dimension(:), allocatable, save  ::  chix
      integer, save                          :: ifirst=1, low_number=-1, high_number=-1
      real, save                             :: lowdt=-1, highdt=-1, l_scale=6., p_scale=1.
      real, save                             :: sc=1.
      real                                   :: x
      real, save                             :: x_low=-1e30, x_length_low=-1e30
      real, save                             :: x_length_high=1e30, x_high=1e30
      integer                                :: ix, nxgrid, nygrid, nzgrid
      namelist /damping_params/  x_low, x_length_low, lowdt, x_length_high, x_high, highdt, l_scale, p_scale

! uses a newton coolin approach to fix the boundary values of density and thermal 
! energy to the intial values at the top and bottom boundary. Indicates the
! time scale for the process. Notice that it gets longer as the distance from the
! boundaries is incerased.

      nxgrid=size(self%mem,1)
      nygrid=size(self%mem,2)
      nzgrid=size(self%mem,3)

      if (ifirst .eq. 1) then
        ifirst=0
        allocate(chix(nxgrid))
        chix=1.

!  default values for input parameters
        rewind (io%input); read (io%input, damping_params, iostat=iostat)
        write(*,damping_params)

        allocate(dissipation(nxgrid, nygrid, nzgrid))
        dissipation=1.

! lower boundary
        if (lowdt .gt. 0) then
          low_number=-1
          do ix=1,nxgrid
            x=m(1)%p+m(1)%r(ix)
            if (x < x_length_low) then
              low_number=ix
              dissipation(ix,:,:) = dissipation(ix,:,:)*exp(-l_scale*(x_length_low-x)/(x_length_low-x_low))
            endif
          enddo
          if (low_number .ne. -1) then
             do ix=1,low_number
                x=m(1)%p+m(1)%r(ix)
                chix(ix)=lowdt*exp(l_scale*(x-x_low)/(x_length_low-x_low))
             end do
          endif
        endif

! upper boundary
        if (highdt .gt. 0) then
          high_number=1
          do ix=1,nxgrid
            x=m(1)%p+m(1)%r(ix)
            if (x <= x_length_high) then
               high_number=ix
            else 
               dissipation(ix,:,:) = dissipation(ix,:,:)*exp(-l_scale*(x-x_length_high)/(x_high-x_length_high))
            endif
          enddo
          if (high_number .ne. nxgrid) then
            do ix=high_number,nxgrid
               x=m(1)%p+m(1)%r(ix)
               chix(ix)=highdt*exp(l_scale*(x_high-x)/(x_high-x_length_high))
            end do
          endif
        endif
      endif

! special scaling for momentum relative e and r
      if (scale .eq. 0) then
          sc=1.
      else
          sc=p_scale
      endif

! lower boundary
      if (lowdt .gt. 0. .and. low_number .ne. -1 ) then
        do ix=1,low_number
          dfdt(ix,:,:) = dfdt(ix,:,:) - (f(ix,:,:)-fbound(ix,:,:))/(chix(ix)*sc)
        enddo
      endif

! upper boundary
      if (highdt .gt. 0. .and. high_number .ne. nxgrid) then
        do ix=high_number,nxgrid
          dfdt(ix,:,:) = dfdt(ix,:,:) - (f(ix,:,:)-fbound(ix,:,:))/(chix(ix)*sc)
        enddo
      endif

END SUBROUTINE damping

!>------------------------------------------------------------------------------------

END MODULE
