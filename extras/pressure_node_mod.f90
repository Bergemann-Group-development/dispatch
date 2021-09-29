!===============================================================================
!> Add a pressure node to any solver
!===============================================================================
MODULE pressure_node_mod
  USE io_mod
  USE link_mod
  USE gpatch_mod
  USE mesh_mod
  USE io_unit_mod
  USE trace_mod
  USE bits_mod
  USE scratch_mod
  USE eos_mod
  USE kinds_mod
  USE trace_mod
  implicit none
  private
  type, public:: pressure_node_t
    logical:: on=.false.
    real(var_k)       :: d_bot, e_bot
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
   ! procedure:: pre_update
    procedure:: post_update
  end type
  integer:: verbose=0
  ! -- the default values are reasonable for the solar photosphere
  real, save:: t_bdry     = 3e+0
  real, save:: bottom_bc_rate=0.5
  real, save:: pb_fact = 1.0
CONTAINS

!===============================================================================
!> Initialize pressure node
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(pressure_node_t):: self
  class(gpatch_t), target:: patch
  ! ............................................................................
  integer:: iostat, m(3)
  logical, save:: on=.false.
  logical, save:: first_time=.true.
  namelist /pressure_node_params/ on, t_bdry, bottom_bc_rate, pb_fact
  !-----------------------------------------------------------------------------
  call trace%begin ('pressure_node_t%pre_init')
  self%patch => patch
  call io%lock%set ('pressure_node_t%pre_init')
  if (first_time) then
    rewind (io_unit%input)
    read (io%input, pressure_node_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io_unit%input)
      read (io%input, pressure_node_params)
    end if
    call io%namelist_warning ('pressure_node_params', iostat)
    write (io%namelist, pressure_node_params)
    first_time = .false.
  end if
  call io%lock%unset ('pressure_node_t%pre_init')
  self%on = on
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize Newton cooling, allocating heating array knowing patch%gn
!===============================================================================
SUBROUTINE init (self)
  class(pressure_node_t):: self
  ! ............................................................................
  integer:: m(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('pressure_node_t%init')
  ! nothing to do here...
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Apply Newton cooling
!===============================================================================
SUBROUTINE post_update (self)
  class(pressure_node_t) :: self
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  class(gpatch_t), pointer:: p
  real(var_k), dimension(:,:,:), contiguous, pointer:: d, e, px, py, pz, bx, by, bz                ! it memory
  real(var_k), dimension(:,:,:), contiguous, pointer:: nd, ne, npx, npy, npz, nbx, nby, nbz        ! new memory
  real(var_k), dimension(:,:),   contiguous, pointer:: u_bnd, p_bot
  real(var_k), dimension(:,:),   contiguous, pointer:: em, nem
  real(var_k), dimension(:,:),   contiguous, pointer:: eth,neth
  real(var_k), dimension(:,:,:), contiguous, pointer:: pg, bx1, by1, bz1
  real(var_k), dimension(:,:,:), contiguous, pointer:: d1, p1, e1
  real(var_k), dimension(:,:,:),             pointer:: p_loc, d_loc, e_loc
  real(var_k), dimension(:,:),   contiguous, pointer:: ee_bnd, d_bnd
  real(var_k), dimension(:,:),   contiguous, pointer:: dlnpdlnr_E, dlnpdE_r
  real(var_k) :: F1,F2,F3,F4,F5,hh,dddt_add, dddt_b, dedt_b, dddt, dedt, ekin
  real(var_k) :: ux, uy, uz
  real:: ct, u_sm, nrm
  real(8):: ds(3), ids(3)
  integer:: it, new, ix, iy, iz, zi, ii, jj, iter
  integer, dimension(3)::lb,ub,lo,uo,li,ui,gn
  integer, save:: itimer=0
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! Check if pressure node is applicable
  ! ----------------------------------------------------------------------------
  if (.not.(self%on)) return
  if (.not.(self%patch%mesh(3)%lower_boundary)) return
  call trace%begin('pressure_node_t%post_update', itimer=itimer)
  self%d_bot = self%patch%boundaries%d_bot
  self%e_bot = self%patch%boundaries%e_bot
  ! ----------------------------------------------------------------------------
  ! Point to memory slots
  ! ----------------------------------------------------------------------------
  p => self%patch
  lb = p%mesh%lb; lo = p%mesh%lo; li = p%mesh%li
  ub = p%mesh%ub; uo = p%mesh%uo; ui = p%mesh%ui
  gn = p%mesh%gn
  ds = p%mesh%d
  ids= 1d0/ds
  it = p%it
  new= p%new
  nd  => p%mem(:,:,:,p%idx%d ,new,1); d  => p%mem(:,:,:,p%idx%d ,it,1)
  ne  => p%mem(:,:,:,p%idx%e ,new,1); e  => p%mem(:,:,:,p%idx%e ,it,1)
  npx => p%mem(:,:,:,p%idx%px,new,1); px => p%mem(:,:,:,p%idx%px,it,1)
  npy => p%mem(:,:,:,p%idx%py,new,1); py => p%mem(:,:,:,p%idx%py,it,1)
  npz => p%mem(:,:,:,p%idx%pz,new,1); pz => p%mem(:,:,:,p%idx%pz,it,1)
  nbx => p%mem(:,:,:,p%idx%bx,new,1); bx => p%mem(:,:,:,p%idx%bx,it,1)
  nby => p%mem(:,:,:,p%idx%by,new,1); by => p%mem(:,:,:,p%idx%by,it,1)
  nbz => p%mem(:,:,:,p%idx%bz,new,1); bz => p%mem(:,:,:,p%idx%bz,it,1)
  ! ----------------------------------------------------------------------------
  ! Compute magnetic energy.
  ! We assume, that %it memory has good magnetic fields, but %new memory needs
  ! them fixed to get the correct thermal energy
  ! ----------------------------------------------------------------------------
  allocate(bx1(gn(1),gn(2),3))
  allocate(by1(gn(1),gn(2),3))
  allocate(bz1(gn(1),gn(2),3))

  allocate(em (gn(1),gn(2)))
  allocate(nem(gn(1),gn(2)))
  !
  do iy=li(2),ui(2)
    do ix=li(1),ui(1)
      em(ix,iy) = ((0.5*(bx(ix,iy,li(3))+bx(ix+1,iy,li(3))))**2 + &
                   (0.5*(by(ix,iy,li(3))+by(ix,iy+1,li(3))))**2 + &
                   (0.5*(bz(ix,iy,li(3))+bz(ix,iy,li(3)+1)))**2)* 0.5
    end do
  end do
  zi=1
  bx1=0.0
  by1=0.0
  bz1=0.0
  do iz=li(3),li(3)+2
    do iy=li(2),ui(2)+1
      do ix=li(1),ui(1)+1
        bx1(ix,iy,zi)=nbx(ix,iy,iz)
        by1(ix,iy,zi)=nby(ix,iy,iz)
        bz1(ix,iy,zi)=nbz(ix,iy,iz)
      end do
    end do
    zi = zi+1
  end do
  ! assume constant upper XY corner
  bx1(uo(1),:,:)    =bx1(ui(1),:,:)
  bx1(:,uo(2),:)    =bx1(:,ui(2),:)
  by1(uo(1),:,:)    =by1(ui(1),:,:)
  by1(:,uo(2),:)    =by1(:,ui(2),:)
  bz1(uo(1),:,:)    =bz1(ui(1),:,:)
  bz1(:,uo(2),:)    =bz1(:,ui(2),:)
  ! Now fix those upper X and Y zones
  do iter=1,4
    iy = ui(2)
    do iz=1,2
      do ix=li(1),ui(1)
        by1(ix,iy+1,iz) =   + by1(ix,iy,iz) &
                            - (ids(1) * (bx1(ix+1,iy,iz) - bx1(ix,iy,iz)) &
                            +  ids(3) * (bz1(ix,iy,iz+1) - bz1(ix,iy,iz))) * ds(2)
      end do
    end do
    ix = ui(2)
    do iz=1,2
      do iy=li(2),ui(2)
        bx1(ix+1,iy,iz) =   + bx1(ix,iy,iz) &
                            - (ids(2) * (by1(ix,iy+1,iz) - by1(ix,iy,iz)) &
                            +  ids(3) * (bz1(ix,iy,iz+1) - bz1(ix,iy,iz))) * ds(1)
      end do
    end do
  end do
  ! Compute %new magnetic energy
  do iy=li(2),ui(2)
    do ix=li(1),ui(1)
      nem(ix,iy) = ((0.5*(bx1(ix,iy,1)+bx1(ix+1,iy,2)))**2 + &
                    (0.5*(by1(ix,iy,1)+by1(ix,iy+1,2)))**2 + &
                    (0.5*(bz1(ix,iy,1)+bz1(ix,iy,  2)))**2)* 0.5
    end do
  end do
  deallocate(bx1,by1,bz1)
  ! ----------------------------------------------------------------------------
  ! Compute thermal energy
  ! ----------------------------------------------------------------------------
  allocate(eth (gn(1),gn(2)))
  allocate(neth(gn(1),gn(2)))
  if (p%kind(1:6)=='ramses') then
    do iy=1,gn(2)
      do ix=1,gn(1)
        eth(ix,iy) =  e(ix,iy,li(3)) - &
              0.5* (px(ix,iy,li(3))**2 + &
                    py(ix,iy,li(3))**2 + &
                    pz(ix,iy,li(3))**2)/ &
                      d(ix,iy,li(3))
      end do
    end do
    do iy=1,gn(2)
      do ix=1,gn(1)
        eth(ix,iy) = eth(ix,iy) - em(ix,iy)
      end do
    end do
    !
    do iy=1,gn(2)
      do ix=1,gn(1)
        neth(ix,iy) =  ne(ix,iy,li(3)) - &
                0.5*(npx(ix,iy,li(3))**2 + &
                      npy(ix,iy,li(3))**2 + &
                      npz(ix,iy,li(3))**2)/ &
                      nd(ix,iy,li(3))
      end do
    end do
    do iy=1,gn(2)
      do ix=1,gn(1)
        neth(ix,iy) = neth(ix,iy) - nem(ix,iy)
      end do
    end do
  else
     eth =>  e(:,:,li(3))
    neth => ne(:,:,li(3))
  end if
  
  ct = 0.5/t_bdry                                                               ! characteristic damping constant
  !ct = bottom_bc_rate!1.0/t_bdry                                                               ! characteristic damping constant
  
  !  --- allocate scratch arrays ---  
  allocate(ee_bnd    (size(d,1),size(d,2)))
  allocate( d_bnd    (size(d,1),size(d,2)))
  allocate(dlnpdE_r  (size(d,1),size(d,2)))
  allocate(dlnpdlnr_E(size(d,1),size(d,2)))
  allocate (p_bot(size(d,1),size(d,2)))
  allocate (u_bnd(size(d,1),size(d,2)))


  d_bnd  = nd(:,:,li(3))
  ee_bnd = neth(:,:)/d_bnd
  !
  ! ----------------------------------------------------------------------------
  ! Compute the pressure pg at the bottom BC, and what we would like it to be
  ! ----------------------------------------------------------------------------
  allocate (d1(size(d,1),size(d,2),2))
  allocate (e1(size(d,1),size(d,2),2))
  allocate (p1(size(d,1),size(d,2),2))
  d1(:,:,1) = self%d_bot                                                        ! mean constant density in the simulation box
  d1(:,:,2) = d_bnd
  e1(:,:,1) = self%e_bot                                                        ! mean constant energy density in the simulation box
  e1(:,:,2) = neth(:,:)                                                         ! new energy density
  !print *, minval(neth), minval(eth), self%e_bot
  p_loc=>p1(li(1):ui(1),li(2):ui(2),:)
  d_loc=>d1(li(1):ui(1),li(2):ui(2),:)
  e_loc=>e1(li(1):ui(1),li(2):ui(2),:)
  call eos%lookup_table (shape(p_loc), pg=p_loc, d=d_loc, e=e_loc, label='pressure node', id=self%patch%id)
  !
  p_bot=bottom_bc_rate*p1(:,:,1) + (1.0-bottom_bc_rate)*p1(:,:,2)
  deallocate (d1, e1, p1)
  ! ----------------------------------------------------------------------------
  ! Get gas flow patterns
  ! ----------------------------------------------------------------------------
  if (p%kind(1:6)=='ramses'.or.p%kind(1:4)=='hlld') then
    u_bnd=npz(:,:,li(3))/d(:,:,li(3))
  else
    u_bnd=0.5*(pz(:,:,li(3))+pz(:,:,li(3)+1))/d(:,:,li(3))
  end if
  ! ----------------------------------------------------------------------------
  ! Get pressure and thermodynamic derivatives for pressure variation, scaling 
  ! factors for the damping of pressure fluctuations.
  ! dP/dee_r with constant density
  ! dP/dr_ee with constant ee=e/d
  ! ----------------------------------------------------------------------------
  call eos%ddeos(ee=ee_bnd, d=d_bnd, dlnpdE_r=dlnpdE_r,dlnpdlnr_E=dlnpdlnr_E)
  
  do iy=li(2), ui(2)
    do ix = li(1), ui(1)
      F1 = dlnpdlnr_E(ix,iy)                                                    ! EOS table dlnpdlnr_E
      F2 = dlnpdE_r(ix,iy)                                                      ! EOS table dlnpdE_r
      hh = ee_bnd(ix,iy)+p_bot(ix,iy)/d_bnd(ix,iy)                              ! E in notes
      F5 = F1 + F2*(hh-ee_bnd(ix,iy))                                           ! H in notes
      F3 = (F1 - F2*ee_bnd(ix,iy))/F5                                           ! denominator of F3 and F4
      F4 = F2/F5
      ! ------------------------------------------------------------------------
      ! Smooth the velocities out to minimize random fluctuations if velocities
      ! are low
      ! ------------------------------------------------------------------------
      !nrm=0.0
      !u_sm = 0.0
      !do jj=iy-1,iy+1
      !  do ii=ix-1,ix+1
      !    u_sm = u_sm + u_bnd(ii,jj)
      !    nrm  = nrm + 1.0
      !  end do
      !end do
      !u_sm = u_sm/nrm
      u_sm = u_bnd(ix,iy)
      ! ------------------------------------------------------------------------
      ! Get time derivatives
      ! ------------------------------------------------------------------------
      dddt = (nd(ix,iy,li(3))-d(ix,iy,li(3)))/self%patch%dtime
      dedt = (neth(ix,iy)-eth(ix,iy))/self%patch%dtime
      ! ------------------------------------------------------------------------
      ! For outflows
      !-------------------------------------------------------------------------
      if (u_sm<=0.) then
        dddt_add = F3*((self%d_bot-d(ix,iy,li(3)))*ct - dddt) + &
                   F4*((self%e_bot-e(ix,iy,li(3)))*ct - dedt)
        dddt = dddt + dddt_add
        dedt = dedt + dddt_add * hh
      ! ------------------------------------------------------------------------
      ! For inflows
      ! ------------------------------------------------------------------------
      else
        dddt = (self%d_bot-  d(ix,iy,li(3))) * ct
        dedt = (self%e_bot-eth(ix,iy))       * ct
      end if
      ! ------------------------------------------------------------------------
      ! Density and energy perturbations from magnetic pressure
      ! ------------------------------------------------------------------------
      dddt_B = -ct*d(ix,iy,li(3))*pb_fact*em(ix,iy)/p_bot(ix,iy)/F5
      dedt_B =  hh*dddt_B
      dddt = dddt + dddt_B
      dedt = dedt + dedt_B
      neth(ix,iy)       = eth(ix,iy)       + dedt*self%patch%dtime
      ! ------------------------------------------------------------------------
      ! Also, modify momentum via velocities
      ! ------------------------------------------------------------------------
      ux = npx(ix,iy,li(3))/nd(ix,iy,li(3))
      uy = npy(ix,iy,li(3))/nd(ix,iy,li(3))
      uz = npz(ix,iy,li(3))/nd(ix,iy,li(3))
      ! modify density
      nd  (ix,iy,li(3)) =   d(ix,iy,li(3)) + dddt*self%patch%dtime
      ! Put the momentum back
      npx(ix,iy,li(3)) = ux*nd(ix,iy,li(3))
      npy(ix,iy,li(3)) = uy*nd(ix,iy,li(3))
      npz(ix,iy,li(3)) = uz*nd(ix,iy,li(3))

      ! ------------------------------------------------------------------------
      ! Put internal energy back to total energy 
      ! ------------------------------------------------------------------------
      if (p%kind(1:6)=='ramses') then
        ne(ix,iy,li(3)) = neth(ix,iy) + nem(ix,iy) + &
                          0.5*(ux**2 + uy**2 + uz**2)*nd(ix,iy,li(3))
      end if
    end do
  end do
  if (p%kind(1:6)=='ramses') then
    deallocate (eth, neth)
  end if

  deallocate(ee_bnd, d_bnd, dlnpdE_r, dlnpdlnr_E, p_bot, u_bnd)
  deallocate(em,nem)
    
call trace%end (itimer)
END SUBROUTINE post_update

END MODULE pressure_node_mod
