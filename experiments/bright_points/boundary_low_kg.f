      MODULE BOUNDARY_CALL

      USE CPARAM
      USE CDATA
      USE CWORK
      USE CSTAGGER
      USE CPROCESSOR
      USE CCOMINUCATION
      USE STAGGER_X_CALL
      USE STAGGER_Y_CALL
      USE STAGGER_Z_CALL
      USE SEND_CALL
      USE CHEAT

      CONTAINS
!*****************************************************************************
!
!  New boundary correcton routines that only works on half point centered
!  variables, e.i. Bx, Px, Ey, Ez, Iy, Iz... resistivity parameters...
!
!  05-sep-00: KG first version
!  03-04-02 : KG first MPI version taken from the HPF version of
!                boundary_low_stretch.f, split in both the y and the z
!                direction.
!
!*****************************************************************************
      subroutine Bx_divB_correction()

      real, DIMENSION(:,:), allocatable :: db, dby, dbz
!----------------------------------------
!  Uses the criteria that div(B)=0 on the boundary to determine the values of Bx
!  in the one point ghostzone outside the numerical domain.
!  bx is half centered in x.
!----------------------------------------

! lower boundary
      allocate (db(nygrid,nzgrid), dby(nygrid,nzgrid))
      allocate (dbz(nygrid,nzgrid))
      if (coords(1) .eq. 0) then
        db=by(1,:,:)
        call ddyup_2d(db, dby)
        db=bz(1,:,:)
        call ddzup_2d(db, dbz)
        bx(1,:,:) = bx(2,:,:) + (dby 
     &   +dbz
     &   )*dx
#if x_stretch==1
     &   /fsh_x(1,:,:)
#endif
      endif

! upper boundary

      if (coords(1) .eq. dims(1)-1) then
        db=by(nxgrid-1,:,:)
        call ddyup_2d(db, dby)
        db=bz(nxgrid-1,:,:)
        call ddzup_2d(db, dbz)

        bx(nxgrid,:,:) = bx(nxgrid-1,:,:)  - (dby
     &    +dbz
     &   )*dx
#if x_stretch==1
     &   /fsh_x(nxgrid,:,:)
#endif
      endif
      
      deallocate(db, dby)
      deallocate(dbz)

      end subroutine
!*****************************************************************************
      subroutine ddyup_2D (f, dby)

      real, DIMENSION(:,:)  :: f, dby
      real, DIMENSION(:,:), allocatable :: fup, flow, rb2d, lb2d

      if (dims(2) .gt. 1) then
         allocate (flow(3,nzgrid), fup(3,nzgrid), rb2d(3,nzgrid), lb2d(3,nzgrid))
         flow = f(1:3,:)
         fup  = f(nygrid-2:nygrid,:)
         call send_y_ghost_2D (fup, flow, lb2d, rb2d)
         deallocate(fup, flow)
      endif

      dby(3:nygrid-3,:) =  
#if y_stretch==1
     &     fsf_y(1,3:nygrid-3,:)*
#endif
     &   (ady*(- f(3:nygrid-3,:) + f(4:nygrid-2,:)) +
     &   bdy*(- f(2:nygrid-4,:) + f(5:nygrid-1,:)) +
     &   cdy*(- f(1:nygrid-5,:) + f(6:nygrid  ,:)))

! only one processor in the y-direction
      if (dims(2) .eq. 1) then
        dby(nygrid-2,:) = 
#if y_stretch==1
     &     fsf_y(1,nygrid-2,:)*
#endif
     &     (ady*(- f(nygrid-2,:) + f(nygrid-1,:)) +
     &     bdy*(- f(nygrid-3,:) + f(nygrid  ,:)) +
     &     cdy*(- f(nygrid-4,:) + f(1   ,:)))
        dby(nygrid-1,:) = 
#if y_stretch==1
     &     fsf_y(1,nygrid-1,:)*
#endif
     &     (ady*(- f(nygrid-1,:) + f(nygrid  ,:)) +
     &     bdy*(- f(nygrid-2,:) + f(1   ,:)) +
     &     cdy*(- f(nygrid-3,:) + f(2   ,:)))
        dby(nygrid  ,:) =
#if y_stretch==1
     &     fsf_y(1,nygrid,:)*
#endif
     &     (ady*(- f(nygrid  ,:) + f(1   ,:)) +
     &     bdy*(- f(nygrid-1,:) + f(2   ,:)) +
     &     cdy*(- f(nygrid-2,:) + f(3   ,:)))
        dby(1   ,:) =
#if y_stretch==1
     &     fsf_y(1,1,:)*
#endif
     &     (ady*(- f(1   ,:) + f(2   ,:)) +
     &     bdy*(- f(nygrid  ,:) + f(3   ,:)) +
     &     cdy*(- f(nygrid-1,:) + f(4   ,:)))
        dby(2   ,:) = 
#if y_stretch==1
     &     fsf_y(1,2,:)*
#endif
     &     (ady*(- f(2   ,:) + f(3   ,:)) +
     &     bdy*(- f(1   ,:) + f(4   ,:)) +
     &     cdy*(- f(nygrid  ,:) + f(5   ,:)))

! more than one processor in the y direction
      else
!        call wait_2D
        dby(nygrid-2,:) = 
#if y_stretch==1
     &     fsf_y(1,nygrid-2,:)*
#endif
     &     (ady*(- f(nygrid-2,:) + f(nygrid-1,:)) +
     &     bdy*(- f(nygrid-3,:) + f(nygrid  ,:)) +
     &     cdy*(- f(nygrid-4,:) + rb2d(1   ,:)))
        dby(nygrid-1,:) = 
#if y_stretch==1
     &     fsf_y(1,nygrid-1,:)*
#endif
     &     (ady*(- f(nygrid-1,:) + f(nygrid  ,:)) +
     &     bdy*(- f(nygrid-2,:) + rb2d(1   ,:)) +
     &     cdy*(- f(nygrid-3,:) + rb2d(2   ,:)))
        dby(nygrid  ,:) =
#if y_stretch==1
     &     fsf_y(1,nygrid,:)*
#endif
     &     (ady*(- f(nygrid  ,:) + rb2d(1   ,:)) +
     &     bdy*(- f(nygrid-1,:) + rb2d(2   ,:)) +
     &     cdy*(- f(nygrid-2,:) + rb2d(3   ,:)))
        dby(1   ,:) =
#if y_stretch==1
     &     fsf_y(1,1,:)*
#endif
     &     (ady*(- f(1   ,:) + f(2   ,:)) +
     &     bdy*(- lb2d(3  ,:) + f(3   ,:)) +
     &     cdy*(- lb2d(2  ,:) + f(4   ,:)))
        dby(2   ,:) = 
#if y_stretch==1
     &     fsf_y(1,2,:)*
#endif
     &     (ady*(- f(2   ,:) + f(3   ,:)) +
     &     bdy*(- f(1   ,:) + f(4   ,:)) +
     &     cdy*(- lb2d(3  ,:) + f(5   ,:)))
         deallocate(lb2d, rb2d)
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine ddzup_2D(f, dbz)
      real, DIMENSION(:,:) :: f, dbz
      real, DIMENSION(:,:), allocatable :: lb2d, rb2d, t1, t2

      if (nzgrid .lt. 6) then
        dbz = 0.
        return
      end if

      if (dims(3) .gt. 1) then
        allocate (lb2d(nygrid, 3), rb2d(nygrid, 3))
        allocate (t1(nygrid, 3), t2(nygrid, 3))
        t2=f(:,1:3)
        t1=f(:,nzgrid-2:nzgrid)
        call send_z_ghost_2D (t1, t2, lb2d, rb2d)
        deallocate (t1,t2)
      endif

      dbz(:,3:nzgrid-3) = 
#if z_stretch==1
     &   fsf_z(1,:,3:nzgrid-3)*
#endif
     &   (adz*(- f(:,3:nzgrid-3) + f(:,4:nzgrid-2)) +
     &   bdz*(- f(:,2:nzgrid-4) + f(:,5:nzgrid-1)) +
     &   cdz*(- f(:,1:nzgrid-5) + f(:,6:nzgrid  )))

!only one processor in the z direction
      if (dims(3) .eq. 1) then
        dbz(:,nzgrid-2) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid-2)*
#endif
     &     (adz*(- f(:,nzgrid-2) + f(:,nzgrid-1)) +
     &     bdz*(- f(:,nzgrid-3) + f(:,nzgrid  )) +
     &     cdz*(- f(:,nzgrid-4) + f(:,1       )))
        dbz(:,nzgrid-1) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid-1)*
#endif
     &     (adz*(- f(:,nzgrid-1) + f(:,nzgrid)) +
     &     bdz*(- f(:,nzgrid-2) + f(:,1     )) +
     &     cdz*(- f(:,nzgrid-3) + f(:,2     )))
        dbz(:,nzgrid  ) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid)*
#endif
     &     (adz*(- f(:,nzgrid  ) + f(:,1)) +
     &     bdz*(- f(:,nzgrid-1) + f(:,2)) +
     &     cdz*(- f(:,nzgrid-2) + f(:,3)))
        dbz(:,1   ) =
#if z_stretch==1
     &   fsf_z(1,:,1)*
#endif
     &     (adz*(- f(:,1       ) + f(:,2)) +
     &     bdz*(- f(:,nzgrid  ) + f(:,3)) +
     &     cdz*(- f(:,nzgrid-1) + f(:,4)))
        dbz(:,2   ) = 
#if z_stretch==1
     &   fsf_z(1,:,2)*
#endif
     &     (adz*(- f(:,2       ) + f(:,3)) +
     &     bdz*(- f(:,1       ) + f(:,4)) +
     &     cdz*(- f(:,nzgrid  ) + f(:,5)))

! more than one processor in the z direction
      else
!        call wait_2D
        dbz(:,nzgrid-2) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid-2)*
#endif
     &     (adz*(- f(:,nzgrid-2) + f(:,nzgrid-1)) +
     &     bdz*(- f(:,nzgrid-3) + f(:,nzgrid  )) +
     &     cdz*(- f(:,nzgrid-4) + rb2d(:,1   )))
        dbz(:,nzgrid-1) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid-1)*
#endif
     &     (adz*(- f(:,nzgrid-1) + f(:,nzgrid  )) +
     &     bdz*(- f(:,nzgrid-2) + rb2d(:,1   )) +
     &     cdz*(- f(:,nzgrid-3) + rb2d(:,2   )))
        dbz(:,nzgrid  ) =
#if z_stretch==1
     &   fsf_z(1,:,nzgrid)*
#endif
     &     (adz*(- f(:,nzgrid  ) + rb2d(:,1   )) +
     &     bdz*(- f(:,nzgrid-1) + rb2d(:,2   )) +
     &     cdz*(- f(:,nzgrid-2) + rb2d(:,3   )))
        dbz(:,1   ) =
#if z_stretch==1
     &   fsf_z(1,:,1)*
#endif
     &     (adz*(- f(:,1   ) + f(:,2   )) +
     &     bdz*(- lb2d(:,3) + f(:,3   )) +
     &     cdz*(- lb2d(:,2) + f(:,4   )))
        dbz(:,2   ) = 
#if z_stretch==1
     &   fsf_z(1,:,2)*
#endif
     &     (adz*(- f(:,2   ) + f(:,3   )) +
     &     bdz*(- f(:,1   ) + f(:,4   )) +
     &     cdz*(- lb2d(:,3) + f(:,5   )))
        deallocate (lb2d, rb2d)
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine Px_correction
!-----------------------------------------
!  Uses closed boundaries. This implies that there is no mass motion across
!  the x boundary at any time. Setting Vx_b=0 => Px_b = 0. Gives the criteria
!  to determine Px in the ghostzone.
!  This routine may also be used for the vx velosity if required.
!-----------------------------------------

! lower boundary
      if(coords(1) .eq. 0) then
#if x_open_lower==0
        px(1 ,:,:) = -px(2   ,:,:)
#else
        px(1,:,:) = px(2,:,:)+(px(2,:,:)-px(3,:,:))
#if x_stretch==1
     &     *ddx(1,:,:)/ddx(2,:,:)
#endif
        where (px(1,:,:) .lt. 0. .and. px(2,:,:) .gt. 0.)
          px(1,:,:) = 0.
        else where (px(1,:,:) .gt. 0. .and. px(2,:,:) .lt. 0.)
          px(1,:,:) = 0.
        end where
#endif
      endif

! upper boundary
      if(coords(1) .eq. dims(1)-1) then
#if x_open_upper==0
        px(nxgrid,:,:) = -px(nxgrid-1,:,:)
#else
#if x_stretch==0
        px(nxgrid,:,:) = 2*px(nxgrid-1,:,:)-px(nxgrid-2,:,:)
#else
        px(nxgrid,:,:) = px(nxgrid-1,:,:)+(px(nxgrid-1,:,:)-px(nxgrid-2,:,:))*ddx(nxgrid-1,:,:)/ddx(nxgrid-2,:,:)
#endif
        where (px(nxgrid,:,:) .lt. 0. .and. px(nxgrid-1,:,:) .gt. 0.)
          px(nxgrid,:,:) = 0.
        else where (px(nxgrid,:,:) .gt. 0. .and. px(nxgrid-1,:,:) .lt. 0.)
          px(nxgrid,:,:) = 0.
        endwhere
#endif
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine E_field_correction(El,Elb)
      real, DIMENSION(:,:,:)   :: El, Elb
!-----------------------------------------
!  The value of the electric field components, Ey, Ez are know on the boundary
!  from the knowledge of Vy, Vz and Bx on the boundary. This can be used to
!  derive the value of E outside the domain requiring it to give the correct
!  result in the induction equation
!-----------------------------------------
! lower boundary
      if(coords(1) .eq. 0) then
#if x_stretch==0
      El(1,:,:) = 2.*Elb(1,:,:) -  El(2,:,:)
#else
      El(1,:,:) = Elb(1,:,:) + (Elb(1,:,:)-El(2,:,:))*ddx(1,:,:)/ddx(2,:,:)
!*fsh_x(1,:,:)/fsh_x(2,:,:)
#endif
      endif

! upper boundary
      if(coords(1) .eq. dims(1)-1) then
#if x_stretch==0
      El(nxgrid,:,:) = 2.*Elb(nxgrid-1,:,:) -  El(nxgrid-1,:,:)
#else
      El(nxgrid,:,:) = Elb(nxgrid-1,:,:) + (Elb(nxgrid-1,:,:)-El(nxgrid-1,:,:))*ddx(nxgrid-1,:,:)/ddx(nxgrid-2,:,:)
#endif
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine I_field_correction (Ic)
      real, DIMENSION(:,:,:)   :: Ic
!-----------------------------------------
!  The current is only known at half point inside the domain. Due to the use
!  of boundary driving the electric field is known at the boundary. The
!  physical correct value of the current is therefore not needed  on the
!  boundary => simpe symmetric extrapolation prevents runaway in the time
!  stepping.
!-----------------------------------------
! lower boundary
      if(coords(1) .eq. 0) then
!       Ic(1,:,:)  = Ic(2,:,:)

#if x_stretch==0
      Ic(1,:,:)  = 2*Ic(2,:,:) - Ic(3,:,:)
#else
      Ic(1,:,:)  = Ic(2,:,:) + (Ic(2,:,:)-Ic(3,:,:))*ddx(1,:,:)/ddx(2,:,:)
#endif

      endif
! upper boundary      
      if(coords(1) .eq. dims(1)-1) then
!        Ic(nxgrid,:,:) = Ic(nxgrid-1,:,:)
#if x_stretch==0
      Ic(nxgrid,:,:) = 2*Ic(nxgrid-1,:,:)-Ic(nxgrid-2,:,:)
#else
      Ic(nxgrid,:,:) = Ic(nxgrid-1,:,:) + (Ic(nxgrid-1,:,:)-Ic(nxgrid-2,:,:))*ddx(nxgrid-1,:,:)/ddx(nxgrid-2,:,:)
#endif

      endif
!-
      end subroutine
!*****************************************************************************
      subroutine Rho_half_correction()
!-----------------------------------------
! Values of rho are used in half staggered locations in the x direction.
! In doing a xup we want to retrive the same boundary values as was already
! there from the beginning. This can easily be achived by using this 
! information to setting the boundary values outside the domain. Just as in the
! same way as for the electric field just above.
!-----------------------------------------
! lower boundary
      if(coords(1) .eq. 0) then
        xdnr(1 ,:,:) = 2.*r(1   ,:,:) - xdnr(2   ,:,:)
      endif
! upper boundary
      if(coords(1) .eq. dims(1)-1) then
        xdnr(nxgrid,:,:) = 2.*r(nxgrid-1,:,:) - xdnr(nxgrid-1,:,:)
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine Px_correction_2
!-----------------------------------------
!  In some cases it is better to set the gradient of Px to zero, by making
!  a symmetric extrapolation into the outside ghostpoint.
!-----------------------------------------
! lower boundary
      if(coords(1) .eq. 0) then
        px(1,:,:) = px(2,:,:)
      endif

! upper boundary
      if(coords(1) .eq. dims(1)-1) then
        px(nxgrid,:,:) = px(nxgrid-1,:,:)
      endif
!-
      end subroutine
!*****************************************************************************
      subroutine boundary_fill()
!
! update the ghostcell at the upper boundary to prevent runaway
!
#if x_periodic==0
      if(coords(1) .eq. dims(1)-1) then
        r (nxgrid,:,:) = r (nxgrid-1,:,:)
        py(nxgrid,:,:) = py(nxgrid-1,:,:)
        pz(nxgrid,:,:) = pz(nxgrid-1,:,:)
        e (nxgrid,:,:) = e (nxgrid-1,:,:)
        by(nxgrid,:,:) = by(nxgrid-1,:,:)
        bz(nxgrid,:,:) = bz(nxgrid-1,:,:)
      endif
#endif

      end subroutine
!*****************************************************************************
      subroutine Byz_divB_correction()

      real, DIMENSION(:,:), allocatable :: db, dby, dbz
!----------------------------------------
!  Uses the criteria that div(B)=0 on the boundary to determine the values of Bx
!  in the one point ghostzone outside the numerical domain.
!  bx is half centered in x.
!----------------------------------------

! lower boundary
      allocate (db(nygrid,nzgrid), dby(nygrid,nzgrid))
      allocate (dbz(nygrid,nzgrid))
      
      rate=.25
      if (coords(1) .eq. 0) then
        do i=1,200 
          db=by(2,:,:)
          call ddyup_2d(db, dby)
          db=bz(2,:,:)
          call ddzup_2d(db, dbz)
          db=dby+dbz+(adx3*(bx(3,:,:)-bx(2,:,:))+bdx3*(bx(4,:,:)-bx(1,:,:)))
#if x_stretch==1
     &     *fsh_x(2,:,:)
#endif
           by(2,:,:)=by(2,:,:)+db*(dy*rate)
#if y_stretch==1
     &     /fsh_y(2,:,:)
#endif
           bz(2,:,:)=bz(2,:,:)+db*(dz*rate)
#if z_stretch==1
     &     /fsh_z(2,:,:)
#endif
        enddo

        do i=1,100 
          do j=0,2
            db=by(3+j,:,:)
            call ddyup_2d(db, dby)
            db=bz(3+j,:,:)
            call ddzup_2d(db, dbz)
            db=dby+dbz+(adx*(bx(4+j,:,:)-bx(3+j,:,:))+bdx*(bx(5+j,:,:)-bx(2+j,:,:))+cdx*(bx(6+j,:,:)-bx(1+j,:,:)))
#if x_stretch==1
     &       *fsh_x(3+j,:,:)
#endif
             by(3+j,:,:)=by(3+j,:,:)+db*(dy*rate)
#if y_stretch==1
     &       /fsh_y(3+j,:,:)
#endif
             bz(3+j,:,:)=bz(3+j,:,:)+db*(dz*rate)
#if z_stretch==1
     &       /fsh_z(3+j,:,:)
#endif
          enddo
        enddo
      endif

! upper boundary

      if (coords(1) .eq. dims(1)-1) then
      endif
      
      deallocate(db, dby)
      deallocate(dbz)

      end subroutine
!*****************************************************************************
      subroutine Anti_field_correction (Ic)
      real, DIMENSION(:,:,:)   :: Ic
!-----------------------------------------
!  The current is only known at half point inside the domain. Due to the use
!  of boundary driving the electric field is known at the boundary. The
!  physical correct value of the current is therefore not needed  on the
!  boundary => simpe symmetric extrapolation prevents runaway in the time
!  stepping.
!-----------------------------------------
! lower boundary
      if(coords(1) .eq. 0) then
#if x_stretch==0
        Ic(1,:,:)  = -Ic(2,:,:)
#else
        Ic(1,:,:)  = -Ic(2,:,:) !+ (Ic(2,:,:)-Ic(3,:,:))*ddx(1,:,:)/ddx(2,:,:)
#endif
      endif

! upper boundary      
      if(coords(1) .eq. dims(1)-1) then
#if x_stretch==0
        Ic(nxgrid,:,:) = -Ic(nxgrid-1,:,:)
#else
        Ic(nxgrid,:,:) = -Ic(nxgrid-1,:,:) !+ (Ic(nxgrid-1,:,:)-Ic(nxgrid-2,:,:))*ddx(nxgrid-1,:,:)/ddx(nxgrid-2,:,:)
#endif
      endif
!-
      end subroutine
!*****************************************************************************
      END MODULE BOUNDARY_CALL
