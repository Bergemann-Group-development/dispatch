! vim: nowrap
!=====================================================================!
SUBROUTINE e_boundaries(ex,ey,ez)
!----------------------------------------------------------------------!
!=> GENERALLY BOUNDARY CONDITIONS FOR E-FIELD MUST SUPPLY:             !
!   !Ex for x={flb-1},{fub  }, y={flb-1},{fub+1}, z={flb-1},{fub+1}    !
!   !Ey for x={flb-1},{fub+1}  y={flb-1},{fub  }, z={flb-1},{fub+1}    !
!   !Ez for x={flb-1},{fub+1}, y={flb-1},{fub+1}, z={flb-1},{fub  }    !
!Such that:                                                            !
!   !Ex is defined on x={flb-1:fub  }, y={flb-1:fub+1}, z={flb-1:fub+1}!
!   !Ey is defined on x={flb-1:fub+1}, y={flb-1:fub  }, z={flb-1:fub+1}!
!   !Ez is defined on x={flb-1:fub+1}, y={flb-1:fub+1}, z={flb-1:fub  }!
!----------------------------------------------------------------------!
  USE grid_m, ONLY: g, e_bnd
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)) :: ex,ey,ez
  integer:: ix, iy, iz

!SET XYZ-PERIODIC BOUNDARY CONDITIONS ON Ex,Ey and Ez.
  ! ----------------------------------------------------------------------------
  ! Set boundary condiitons on Ex, Ey, Ez, by taking the older time slots with
  ! up-to-date value from nboring patches or whatever boundary conditions are
  ! applied externally.
  ! ----------------------------------------------------------------------------
  do iz=1,g%n(3)
    do iy=1,g%n(2)
      do ix=1, g%lb(1)-1
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      enddo
      do ix=g%ub(1),g%n(1)
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      enddo
    end do
  end do

  do iz=1,g%n(3)
    do iy=1, g%lb(2)-1
      do ix=1,g%n(1)
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      end do
      end do
      do iy=g%ub(2),g%n(2)
      do ix=1,g%n(1)
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      end do
    end do
  end do

  do iz=1, g%lb(3)-1
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      end do
    end do
  end do
  do iz=g%ub(3),g%n(3)
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        Ex(ix,iy,iz) = E_bnd(ix,iy,iz,1)
        Ey(ix,iy,iz) = E_bnd(ix,iy,iz,2)
        Ez(ix,iy,iz) = E_bnd(ix,iy,iz,3)
      end do
    end do
  end do
     !=> Ex is now defined on {flb-1:fub  }{flb-1:fub+1}{flb-1:fub+1}
     !=> Ey is now defined on {flb-1:fub+1}{flb-1:fub  }{flb-1:fub+1}
     !=> Ez is now defined on {flb-1:fub+1}{flb-1:fub+1}{flb-1:fub  }
END SUBROUTINE e_boundaries

!=======================================================================
SUBROUTINE j_boundaries_noadd(jx,jy,jz)
  USE grid_m, ONLY: g
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)):: jx, jy, jz
!-----------------------------------------------------------------------
! The physical BCs for the electric current.
!-----------------------------------------------------------------------
  call overlap(jx)
  call overlap(jy)
  call overlap(jz)
END SUBROUTINE j_boundaries_noadd

!=======================================================================
SUBROUTINE j_boundaries(jx,jy,jz)
  USE grid_m, ONLY: g, J_bnd
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)):: jx, jy, jz
!-----------------------------------------------------------------------
! The physical BCs for the electric current.
!-----------------------------------------------------------------------
  integer:: ix, iy, iz
  do iz=1,g%n(3)
    do iy=1,g%n(2)
      do ix=1, g%lb(1)-1
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      enddo
      do ix=g%ub(1),g%n(1)
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      enddo
    end do
  end do

  do iz=1,g%n(3)
    do iy=1, g%lb(2)-1
      do ix=1,g%n(1)
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      end do
      end do
      do iy=g%ub(2),g%n(2)
      do ix=1,g%n(1)
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      end do
    end do
  end do

  do iz=1, g%lb(3)-1
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      end do
    end do
  end do
  do iz=g%ub(3),g%n(3)
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        Jx(ix,iy,iz) = J_bnd(ix,iy,iz,1)
        Jy(ix,iy,iz) = J_bnd(ix,iy,iz,2)
        Jz(ix,iy,iz) = J_bnd(ix,iy,iz,3)
      end do
    end do
  end do
END SUBROUTINE j_boundaries

!=======================================================================
SUBROUTINE charge_boundaries_noadd(charge)
  USE grid_m, ONLY: g
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)) :: charge

  call overlap(charge)
END SUBROUTINE charge_boundaries_noadd

!=======================================================================
SUBROUTINE charge_boundaries(charge)
  USE grid_m, ONLY: g
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)):: charge

  call overlap_add(charge)
END SUBROUTINE charge_boundaries

!=====================================================================!
SUBROUTINE b_boundaries(bx,by,bz)

!---------------------------------------------------------------------!
!=> GENERALLY BOUNDARY CONDITIONS FOR B-FIELD MUST SUPPLY:            !
!     !Bx for x={flb-1},{fub+1}, y={flb-1},{fub  }, z={flb-1},{fub  } !
!     !By for x={flb-1},{fub  }  y={flb-1},{fub+1}, z={flb-1},{fub  } !
!     !Bz for x={flb-1},{fub  }, y={flb-1},{fub  }, z={flb-1},{fub+1} !
!Such that:                                                           !
!     !Bx defined on x={flb-1:fub+1}, y={flb-1:fub},   z={flb-1:fub}  !
!     !By defined on x={flb-1:fub},   y={flb-1:fub+1}, z={flb-1:fub}  !
!     !Bz defined on x={flb-1:fub},   y={flb-1:fub},   z={flb-1:fub+1}!
!---------------------------------------------------------------------!

  USE grid_m, ONLY: g, b_bnd
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)) :: bx
  real, dimension(g%n(1),g%n(2),g%n(3)) :: by
  real, dimension(g%n(1),g%n(2),g%n(3)) :: bz
  integer:: ix, iy, iz

  ! ----------------------------------------------------------------------------
  ! Set boundary condiitons on Bx, By, Bz, by taking the older time slots with
  ! up-to-date value from nboring patches or whatever boundary conditions are
  ! applied externally.
  ! ----------------------------------------------------------------------------
  do iz=1,g%n(3)
    do iy=1,g%n(2)
      do ix=1, g%lb(1)-1
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      enddo
      do ix=g%ub(1),g%n(1)
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      enddo
    end do
  end do

  do iz=1,g%n(3)
    do iy=1, g%lb(2)-1
      do ix=1,g%n(1)
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      end do
      end do
      do iy=g%ub(2),g%n(2)
      do ix=1,g%n(1)
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      end do
    end do
  end do

  do iz=1, g%lb(3)-1
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      end do
    end do
  end do
  do iz=g%ub(3),g%n(3)
    do iy=1,g%n(2)
      do ix=1,g%n(1)
        bx(ix,iy,iz) = b_bnd(ix,iy,iz,1)
        by(ix,iy,iz) = b_bnd(ix,iy,iz,2)
        bz(ix,iy,iz) = b_bnd(ix,iy,iz,3)
      end do
    end do
  end do

END SUBROUTINE b_boundaries
!=====================================================================!
SUBROUTINE phi_boundaries(phi)
  USE grid_m, ONLY: g
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)) :: phi

!SET BOUNDARY CONDITIONS ON 'PHI' -- phi is needed on [flb,fub+1], and
!is known on [flb+1,fub]. 

  call overlap(phi)

END SUBROUTINE phi_boundaries
!=====================================================================!
SUBROUTINE phib_boundaries(phi)
  USE grid_m, ONLY: g
  implicit none
  real, dimension(g%n(1),g%n(2),g%n(3)) :: phi

!SET BOUNDARY CONDITIONS ON 'PHI' -- phi is needed on [flb,fub+1], and
!is known on [flb+1,fub]. 

  call overlap(phi)

END SUBROUTINE phib_boundaries

!--------------------------------------------------------------------------!
! THE BOUNDARY CONDITIONS ARE TRICKY AT BEST
!
!       z
!       ^
!       |
!       Ez            Bx
!       |
!       |
!    By |
!       |
!       |
!       |
!       *-------------Ey---> y
!      /
!     /
!    Ex             Bz
!   /
!  x
!
!  we need, so to speak, one extra grid point at the
!    min location. There only Bt and En are of
!    relevance and need to be supplied via bc.
!  at the max location, again bc need to be supplied only for Bt and En
!    these have now the same grid index as the outermost cell
!    Bn can be calculated from the existing info, and is
!    already calculated at the min boundaries (the ghost cell is irrelevant)
!
!  sx and sy are half point centered, as bx and by
!----------------------------------------------------------------------!
!=======================================================================
SUBROUTINE impose_particle_boundaries
  implicit none
  ! No special boundaries. SendParticles takes care of sending between mpi-threads
END SUBROUTINE impose_particle_boundaries
!=======================================================================
SUBROUTINE impose_source_boundaries
!-----------------------------------------------------------------------
! Nothing to do. Overlap is done on EM fields and at I/O not on fields%...
!-----------------------------------------------------------------------
  implicit none
END SUBROUTINE impose_source_boundaries
!=======================================================================
SUBROUTINE ddj_boundaries(n,u,v,q,dd,dir,isp)
  implicit none
  integer, intent(in) :: n, dir, isp
  real,    intent(in) :: dd(3)
  real                :: v(n-2,2), q(2,2)
  real, target        :: u(n-2,2)
!-----------------------------------------------------------------------
! Peridic case
!-----------------------------------------------------------------------
  u(1,1) = dd(3);   u(1,2) = dd(2);   u(2,2) = dd(3)
  v(1,1) = dd(3);   v(1,2) = dd(2);   v(2,2) = dd(3)
END SUBROUTINE ddj_boundaries
