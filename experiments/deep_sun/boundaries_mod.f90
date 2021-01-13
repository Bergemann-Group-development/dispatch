!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Default variables
    ! --------------------------------------------------------------------------
    logical:: first_lower=.true.
    logical:: first_upper=.true.
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real(var_k)       :: gravity = -2.74
    real(var_k)       :: d_bot, e_bot, ee_top
    procedure(lower1), pointer     :: lower
    procedure(upper1), pointer     :: upper
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: lower1
    procedure:: lower2
    procedure:: lower_trap
    procedure:: lower_simple
    procedure:: upper1
    procedure:: upper2
    procedure:: radiation => radiation2
    procedure:: radiation2
    procedure:: div_free_lower
    procedure:: div_free_upper
    !
    procedure:: lower_damp
    procedure:: ddt_boundary_bot
    procedure:: ddt_boundary_top
    procedure:: sponge_boundaries
  end type
  integer, save:: verbose=1
  integer, save:: lower_bc=3
  integer, save:: upper_bc=2
  logical, save:: smudge_boundary = .false.
  logical, save:: scale = .true.
  logical, save:: deriv_bnd = .false.
  logical, save:: test = .false.
  real, save:: d_bot=0.0, e_bot=0.0, ee_top=0.0
  real, save:: bottom_bc_rate=0.01
  real, save:: top_bc_rate=0.5
  real, save:: damp_rate=0.01
  real, save:: r_sponge = 0.5
  real, save:: pb_fact = 1.0
  integer, save:: divb_iter = 4
  ! --- lower_damp parameters
  real, save:: beta_lim   = 1e+2
  real, save:: spot_lower = 0.0
  real, save:: pb_spot    = 0.0
  real, save:: t_bdry     = 3e+0
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t) :: self
  class(patch_t), target:: patch
  ! ----------------------------------------------------------------------------
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /boundary_params/ d_bot, e_bot, ee_top, bottom_bc_rate, verbose, divb_iter, &
                             scale, t_bdry, r_sponge, damp_rate, deriv_bnd, pb_fact, &
                             test, lower_bc, upper_bc
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (io%master) write (*, boundary_params)
    !d_bot = d_bot/scaling%d
    !e_bot = e_bot/scaling%e
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  if (.not. io%master) verbose=0

  select case(lower_bc)
  case(1)
    self%lower => lower1
  case(2)
    self%lower => lower2
  case(3)
    self%lower => lower_trap
  case(4)
    self%lower => lower_simple
  case default
    self%lower => lower_trap
  end select

  select case(upper_bc)
  case(1)
    self%upper => upper1
  case(2)
    self%upper => upper2
  case default
    self%upper => upper2
  end select

  self%gravity = -scaling%grav*(cgs%m_sun/scaling%m)/(cgs%r_sun/scaling%l)**2


  call self%standard_boundaries_t%init(patch)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:), pointer:: d, e, px, py, pz
  real, dimension(:,:,:), pointer:: bx, by, bz
  if (test) then
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    if (self%mesh(3)%lower_boundary.or.self%mesh(3)%upper_boundary) then
    call self%log_extrapolate(d,idx%d)
    call self%log_extrapolate(e,idx%e)
    call self%zero_gradient(px,idx%px)
    call self%zero_gradient(py,idx%py)
    call self%zero_gradient(pz,idx%pz)
    if (self%mhd) then
      call self%symmetric(bx,idx%bx)
      call self%symmetric(by,idx%by)
      call self%zero_gradient(bz,idx%bz)
      call self%div_free_lower (bx, by, bz)
    end if
    end if
  else
    if (self%mesh(3)%lower_boundary) call self%lower(mem)
    if (self%mesh(3)%upper_boundary) call self%upper(mem)
  end if

END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  if (deriv_bnd) then
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)
    !if (self%mesh(3)%lower_boundary) call self%sponge_boundaries(mem, lower=.true.)
    !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)
    if (self%mesh(3)%lower_boundary) call self%ddt_boundary_bot(mem)
    if (self%mesh(3)%upper_boundary) call self%ddt_boundary_top(mem)
  end if

END SUBROUTINE derivative

!===============================================================================
!> Top boundary conditions. It is the top of the
!> atmosphere.
!===============================================================================
SUBROUTINE upper1 (self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n, gn, ng
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, eth, emag
  real, dimension(:,:,:), pointer:: bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), allocatable:: pg
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1, i, zi, bn
  real:: d_ratio, e_ratio, h_scale, p, ee, dlnddz, dlnedz, uz
  !-----------------------------------------------------------------------------
  !if (.not.self%first_upper) return

  call trace%begin('boundaries_t%bc_upper')
  if (self%first_upper) then
    self%first_upper=.false.
    if (ee_top/=0.0) self%ee_top = ee_top
  end if
  !
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  gn = self%mesh%gn
  ng = self%mesh%ng
  m3 => self%mesh(3)
  bn = ui(3)
    z  => m3%r
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    ! --------------------------------------------------------------------------
    ! First, fix the magnetic fields, as they are needed to get thermal energy
    ! --------------------------------------------------------------------------
    if (self%mhd) then
    do iz=uo(3),ub(3)
        bx(:,:,iz) = bx(:,:,ui(3))
        by(:,:,iz) = by(:,:,ui(3))
        bz(:,:,iz) = bz(:,:,ui(3))
      end do
    !call self%zero_gradient(bx,idx%bx)
    !call self%zero_gradient(by,idx%by)
    !call self%zero_gradient(bz,idx%bz)
    !call self%symmetric(bx,idx%bx)
    !  call self%symmetric(by,idx%by)
    !  call self%zero_gradient(bz,idx%bz)
    call self%div_free_upper(bx,by,bz)
  end if
    !---------------------------------------------------------------------------
    ! Boundary conditions need the gas pressure
    !---------------------------------------------------------------------------
    d1 => d(:,:,ui(3):ui(3))
    if (self%patch%kind(1:6)/='ramses') then
      e1  => e(:,:,bn:bn)
      eth => e(:,:,uo(3):)
    else
      allocate(e1 (gn(1),gn(2),1))
      allocate(eth(gn(1),gn(2),ng(3)))
      allocate(emag(gn(1),gn(2),ng(3)))
      do iy=1,gn(2)
        do ix=1,gn(1)
          e1(ix,iy,1) =  e(ix,iy,bn) - &
                   0.5*(px(ix,iy,bn)**2 + &
                        py(ix,iy,bn)**2 + &
                        pz(ix,iy,bn)**2)/ &
                         d(ix,iy,bn)
        end do
      end do
      !
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          e1(ix,iy,1) = &
          e1(ix,iy,1) - 0.125* ((bx(ix,iy,bn)+bx(ix+1,iy,bn))**2 + &
                                (by(ix,iy,bn)+by(iy,iy+1,bn))**2 + &
                                (bz(ix,iy,bn)+bz(iz,iy,bn+1))**2)
        end do
      end do
      e1(gn(1),:,1)=e1(gn(1)-1,:,1)
      e1(:,gn(2),1)=e1(:,gn(2)-1,1)
    end if
  !---------------------------------------------------------------------------
  ! Boundary conditions need the gas pressure in one layer
  !---------------------------------------------------------------------------
  allocate (pg(size(d,1),size(d,2),1))

  call eos%lookup_table (shape(pg), pg=pg, d=d1, e=e1, label='upper_bnd', id=self%id)
  zi = 1
  do iz=uo(3),ub(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      !-----------------------------------------------------------------------
      ! 1) The BC for energy per unit mass is that it is = the value at ui,
      ! while the density is assumed to drop with a scale height that is to
      ! lowest order proportional to temperature.
      !-----------------------------------------------------------------------
      iz1 = ui(3)
      h_scale = pg(ix,iy,1)/(d(ix,iy,iz1)*self%gravity)
      ee = e1(ix,iy,1)/d(ix,iy,iz1)
      d_ratio = exp((z(iz)-z(iz1))/h_scale)
      d(ix,iy,iz) = d(ix,iy,iz1)*d_ratio
      eth(ix,iy,zi) = d(ix,iy,iz)*ee
      !-----------------------------------------------------------------------
      ! 2) The horizontal components of the momenta are assumed to taper off
      ! with density, consistent with assuming a no-shear boundary condition
      ! for the horizontal velocity components
      !-----------------------------------------------------------------------
      px(ix,iy,iz) = px(ix,iy,iz1)*d_ratio
      py(ix,iy,iz) = py(ix,iy,iz1)*d_ratio
      !-----------------------------------------------------------------------
      ! 3) The vertical velocity is tapered of with the d_ratio factor.
      !-----------------------------------------------------------------------
      pz(ix,iy,iz) = pz(ix,iy,iz1)*d_ratio**2
    end do
    end do
    zi = zi+1
  end do
  deallocate (pg)
    ! --------------------------------------------------------------------------
    ! For Riemann solver, convert back to total energy
    ! --------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    e(:,:,uo(3):) = eth + &
         0.5*(px(:,:,uo(3):)**2 + &
              py(:,:,uo(3):)**2 + &
              pz(:,:,uo(3):)**2)/ d(:,:,uo(3):)
    zi=1
    do iz=uo(3),ub(3)-1
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,zi) = &
                          ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                           (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                           (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)* 0.125
        end do
      end do
      zi = zi+1
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    emag(:,:,ng(3))=emag(:,:,ng(3)-1)
    e(:,:,uo(3):) = e(:,:,uo(3):) + emag
    deallocate(eth,emag)
  end if 
call trace%end()
END SUBROUTINE upper1

!===============================================================================
!> Top boundary conditions. It is the top of the
!> atmosphere. Same as upper1, but with ddt boundaries incorporated.
!===============================================================================
SUBROUTINE upper2 (self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n, gn, ng
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, eth, emag
  real, dimension(:,:,:), pointer:: bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:), pointer:: em, ep
  real, dimension(:,:,:), allocatable:: pg
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1, i, zi, bn
  real:: d_ratio, e_ratio, h_scale, p, ee, dlnddz, dlnedz, uz, dddt
  !-----------------------------------------------------------------------------

  call trace%begin('boundaries_t%bc_upper')
  if (self%first_upper) then
    self%first_upper=.false.
    if (ee_top/=0.0) self%ee_top = ee_top
  end if
  !
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  gn = self%mesh%gn
  ng = self%mesh%ng
  m3 => self%mesh(3)
  bn = ui(3)
    z  => m3%r
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    ! --------------------------------------------------------------------------
    ! First, fix the magnetic fields, as they are needed to get thermal energy
    ! --------------------------------------------------------------------------
    if (self%mhd) then
      do iz=uo(3),ub(3)
        bx(:,:,iz) = bx(:,:,bn)
        by(:,:,iz) = by(:,:,bn)
        bz(:,:,iz) = bz(:,:,bn)
      end do
    !call self%extrapolate(bx,idx%bx)
    !call self%extrapolate(by,idx%by)
    !call self%zero_gradient(bz,idx%bz)
    !call self%symmetric(bx,idx%bx)
    !  call self%symmetric(by,idx%by)
    !  call self%zero_gradient(bz,idx%bz)
    call self%div_free_upper(bx,by,bz)
  end if
    !---------------------------------------------------------------------------
    ! Boundary conditions need the gas pressure
    !---------------------------------------------------------------------------
    d1 => d(:,:,ui(3):ui(3))
    allocate(e1 (gn(1),gn(2),1))
    allocate(em (gn(1),gn(2)))
    do iy=1,gn(2)-1
      do ix=1,gn(1)-1
        em(ix,iy) = ((0.5*(bx(ix,iy,bn)+bx(ix+1,iy,bn)))**2 + &
                     (0.5*(by(ix,iy,bn)+by(ix,iy+1,bn)))**2 + &
                     (0.5*(bz(ix,iy,bn)+bz(ix,iy,bn+1)))**2)* 0.5
      end do
    end do
    em(gn(1),:)=em(gn(1)-1,:)
    em(:,gn(2))=em(:,gn(2)-1)
    if (self%patch%kind(1:6)/='ramses') then
      e1(:,:,1)  = e(:,:,bn)/d(:,:,bn)
      eth => e(:,:,uo(3):ub(3))
    else
      allocate(ep (gn(1),gn(2)))
      allocate(eth(gn(1),gn(2),ng(3)))
      allocate(emag(gn(1),gn(2),ng(3)))
      do iy=1,gn(2)
        do ix=1,gn(1)
          ep(ix,iy) = 0.5* (px(ix,iy,bn)*px(ix,iy,bn) + &
                            py(ix,iy,bn)*py(ix,iy,bn) + &
                            pz(ix,iy,bn)*pz(ix,iy,bn))/d1(ix,iy,1)
        end do
      end do
      do iy=1,gn(2)
        do ix=1,gn(1)
          e1(ix,iy,1) = (e(ix,iy,bn) - em(ix,iy) - ep(ix,iy)) / d1(ix,iy,1)
        end do
      end do     
    end if
  !---------------------------------------------------------------------------
  ! Boundary conditions need the gas pressure in one layer
  !---------------------------------------------------------------------------
  allocate (pg(size(d,1),size(d,2),1))
!if (self%id==757) then
if (any(e1<eos%ee_min)) then
  print *, self%id
  print *, minval (e1), maxval(e1),eos%ee_min, minloc(e1), self%ee_top
  print *, em(22,17), ep(22,17), e(22,17,bn)
  print *, 'b', bx(22,17,bn), by(22,17,bn), bz(22,17,bn)
  print *, 'pz', pz(ix,iy,bn)
  if (any(e1<eos%ee_min))error stop "e1 leq ee_min"
end if
deallocate(ep)
  call eos%lookup_table (shape(pg), pg=pg, d=d1, ee=e1, label='upper_bnd2', id=self%id)
  zi = 1
  do iz=uo(3),ub(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      !-----------------------------------------------------------------------
      ! 1) The BC for energy per unit mass is that it is = the value at ui,
      ! while the density is assumed to drop with a scale height that is to
      ! lowest order proportional to temperature.
      !-----------------------------------------------------------------------
      h_scale = pg(ix,iy,1)/(d(ix,iy,bn)*self%gravity)
      if (pz(ix,iy,bn) < 0.0) then
        ee = top_bc_rate*self%ee_top + (1.0-top_bc_rate)*e1(ix,iy,1)
      else
        ee = e1(ix,iy,1)
      end if
      
      d_ratio = exp((z(iz)-z(bn))/h_scale)
      d(ix,iy,iz) = d(ix,iy,bn)*d_ratio

      eth(ix,iy,zi) = d(ix,iy,iz)*ee
      !-----------------------------------------------------------------------
      ! 2) The horizontal components of the momenta are assumed to taper off
      ! with density, consistent with assuming a no-shear boundary condition
      ! for the horizontal velocity components
      !-----------------------------------------------------------------------
      px(ix,iy,iz) = px(ix,iy,bn)*d_ratio
      py(ix,iy,iz) = py(ix,iy,bn)*d_ratio
      !-----------------------------------------------------------------------
      ! 3) The vertical velocity is tapered of with the d_ratio factor.
      !-----------------------------------------------------------------------
      pz(ix,iy,iz) = pz(ix,iy,bn)*d_ratio**2
    end do
    end do
    zi = zi+1
  end do
  deallocate (pg)
    ! --------------------------------------------------------------------------
    ! For Riemann solver, convert back to total energy
    ! --------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    zi=1
    do iz=uo(3),ub(3)
      do iy=lb(2),ub(2)
        do ix=lb(1),ub(1)
          e(ix,iy,iz) = eth(ix,iy,zi) + &
              0.5*(px(ix,iy,iz)**2 + &
                   py(ix,iy,iz)**2 + &
                   pz(ix,iy,iz)**2)/ d(ix,iy,iz)
        end do
      end do
      zi = zi+1
    end do
    zi=1
    do iz=uo(3),ub(3)-1
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,zi) = &
                          ((0.5*(bx(ix,iy,iz)+bx(ix+1,iy,iz)))**2 + &
                           (0.5*(by(ix,iy,iz)+by(iy,iy+1,iz)))**2 + &
                           (0.5*(bz(ix,iy,iz)+bz(iz,iy,iz+1)))**2)* 0.5*1.01
        end do
      end do
      zi = zi+1
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    emag(:,:,ng(3))=emag(:,:,ng(3)-1)
    zi=1
    do iz=uo(3),ub(3)
      do iy=lb(2),ub(2)
        do ix=lb(1),ub(1)
          e(ix,iy,iz) = e(ix,iy,iz) + emag(ix,iy,zi)
        end do
      end do
      zi = zi+1
    end do
    deallocate(eth,emag, e1,em)
  end if 
call trace%end()
END SUBROUTINE upper2

!===============================================================================
!> Top boundary conditions.  Since the Z is flipped, it is the bottom of the
!> atmosphere.
!===============================================================================
SUBROUTINE lower1 (self, mem)
  USE vector_ops
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, gn
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), pointer::  e2, eth, emag, ux, uy, uz
  real, dimension(:,:,:), allocatable:: pg
  real, dimension(:,:,:,:), pointer:: dd
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1
  real:: d_ratio, e_ratio, ee_ratio, h_scale, p, ee, dlnddz, dlnedz, dlneedz, uux, uuy, uuz, w
  !-----------------------------------------------------------------------------
  if (.not.self%first_lower) then
    call trace%begin('boundries_t%bc_lower')
    !
    if (self%first_lower) then
      self%first_lower = .false.
      if (d_bot /= 0.0) self%d_bot = d_bot
      if (e_bot /= 0.0) self%e_bot = e_bot
    end if
    lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
    ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
    gn = self%mesh%gn
    m3 => self%mesh(3)
    z  => m3%r
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    !---------------------------------------------------------------------------
    ! Compute the pressure pg at the bottom BC, and what we would like it to be
    !---------------------------------------------------------------------------
    d1 => d(:,:,li(3):li(3))
    if (self%patch%kind(1:6)/='ramses') then
      e1  => e(:,:,li(3):li(3))
      eth => e
    else
      allocate(e1 (gn(1),gn(2),1))
      allocate(eth(gn(1),gn(2),gn(3)))
      allocate(emag(gn(1),gn(2),gn(3)))
      e1 = e(:,:,li(3):li(3)) - &
         0.5*(px(:,:,li(3):li(3))**2 + &
              py(:,:,li(3):li(3))**2 + &
              pz(:,:,li(3):li(3))**2)/ d(:,:,li(3):li(3))
      eth = e - &
           0.5*(px**2 + &
                py**2 + &
                pz**2)/ d
      iz = li(3)
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          e1(ix,iy,1) = &
          e1(ix,iy,1) - 0.125*((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                                (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                                (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
        end do
      end do
      e1(gn(1),:,:)=e1(gn(1)-1,:,:)
      e1(:,gn(2),:)=e1(:,gn(2)-1,:)
      !
      do iz=1,gn(3)-1
        do iy=1,gn(2)-1
          do ix=1,gn(1)-1
            emag(ix,iy,iz) = &
                             ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                             (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                             (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
          end do
        end do
      end do
      emag(gn(1),:,:)=emag(gn(1)-1,:,:)
      emag(:,gn(2),:)=emag(:,gn(2)-1,:)
      emag(:,:,gn(3))=emag(:,:,gn(3)-1)
      eth = eth - 0.125*emag
    end if
    !
    allocate (pg(size(d,1),size(d,2),1))
    call eos%lookup_table (shape(pg), pg=pg, d=d1, e=e1, label='lower1 bnd', id=self%id)
    !
    allocate (d1(size(d,1),size(d,2),1))
    allocate (e2(size(d,1),size(d,2),1))
    allocate (p1(size(d,1),size(d,2),1))
    d1 = self%d_bot
    e2 = self%e_bot
    call eos%lookup_table (shape(p1), pg=p1, d=d1, e=e2, label='lower1 bnd', id=self%id)
    !---------------------------------------------------------------------------
    if (scale) then
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      dlnddz = log(d(ix,iy,li(3)+3)/d(ix,iy,li(3)+1))/(z(li(3)+1)-z(li(3)+3))
      dlnedz = log(eth(ix,iy,li(3)+3)/eth(ix,iy,li(3)+1))/(z(li(3)+1)-z(li(3)+3))
      dlneedz = log((eth(ix,iy,li(3)+3)/d(ix,iy,li(3)+3)) &
                    /eth(ix,iy,li(3)+1)/d(ix,iy,li(3)+1))/(z(li(3)+1)-z(li(3)+3))
      !-------------------------------------------------------------------------
      ! We assume that the iz index of the last free values is li(3), except
      ! for pz, for which it us li(3)+1.  Hence, if there are 2 guard zones, we
      ! set 3 values for pz, and 2 values for the other variables
      !-------------------------------------------------------------------------
      w = bottom_bc_rate
      do iz=li(3),lb(3),-1
        if (iz == li(3)) then
            d1(ix,iy,1) = d(ix,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
            e2(ix,iy,1) = eth(ix,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
            d(ix,iy,iz) = w*d1(ix,iy,1) + (1.0-w)*d(ix,iy,iz)
            eth(ix,iy,iz) = w*e2(ix,iy,1) + (1.0-w)*eth(ix,iy,iz)
          !---------------------------------------------------------------------
          ! Note that pz(ix,iy,iz) is centered half a cell above pressure, so
          ! by imposing pz at ui(3)+1=uo(3), we are in principle imposing the
          ! density evolution at uo(3), but that is overruled in any case here,
          ! hence it is not necessary to stagger the BC for pz
          !---------------------------------------------------------------------
        else
            d_ratio = exp((z(li(3))-z(iz))*dlnddz)
            e_ratio = exp((z(li(3))-z(iz))*dlnedz)
            ee_ratio = exp((z(li(3))-z(iz))*dlneedz)
            d(ix,iy,iz)  =  d(ix,iy,li(3))*d_ratio
            eth(ix,iy,iz)  =  eth(ix,iy,li(3))*e_ratio
        end if
        !-----------------------------------------------------------------------
        ! For inflow, set d,e to d_bot,e_bot at li, and set the momenta constant,
        ! corresponding to velocities that taper off with depth
        !-----------------------------------------------------------------------
        if (pz(ix,iy,li(3)) > 0.0) then
          px(ix,iy,iz) = px(ix,iy,li(3))
          py(ix,iy,iz) = py(ix,iy,li(3))
          pz (ix,iy,iz) = pz(ix,iy,li(3))!/d(ix,iy,li(3)) * (iz-1)/float(m3%ng) * &
                                         ! d(ix,iy,iz)
        !-----------------------------------------------------------------------
        ! For outflow, extrapolate d and e, and assume velocity is constant.
        ! Leave values untouched at iz=ui(3), to affect outflow as little as
        ! possible
        !-----------------------------------------------------------------------
        else
          !---------------------------------------------------------------------
          ! At iz=li(3), nudge the density towards what would give a constant
          ! pressure = pg(d_bot,e_bot), while keeping ee=e/d unchanged
          !---------------------------------------------------------------------
            px(ix,iy,iz) = px(ix,iy,li(3))
            py(ix,iy,iz) = py(ix,iy,li(3))
            pz(ix,iy,iz) = pz(ix,iy,li(3))
        
        end if
      end do
    end do
    end do
    else
      allocate (ux(size(d,1),size(d,2),size(d,3)))
      allocate (uy(size(d,1),size(d,2),size(d,3)))
      allocate (uz(size(d,1),size(d,2),size(d,3)))
      allocate (dd(size(d,1),size(d,2),size(d,3),3))
      dd=down1(d)
      ux = px/dd(:,:,:,1)
      uy = py/dd(:,:,:,2)
      uz = pz/dd(:,:,:,3)
      w = bottom_bc_rate
      do iy=lb(2),ub(2)
        do ix=lb(1),ub(1)
          iz = li(3)
          d1(ix,iy,1) = d(ix,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
          e1(ix,iy,1) = eth(ix,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
          d(ix,iy,iz) = w*d1(ix,iy,1) + (1.0-w)*d(ix,iy,iz)
          eth(ix,iy,iz) = w*e1(ix,iy,1) + (1.0-w)*eth(ix,iy,iz)
        end do
      end do
      call self%log_extrapolate(d,idx%d)
      call self%log_extrapolate(eth,idx%e)
      call self%zero_gradient(ux,idx%px)
      call self%zero_gradient(uy,idx%py)
      call self%zero_gradient(uz,idx%pz)
      dd=down1(d)
      px = ux * dd(:,:,:,1)
      py = uy * dd(:,:,:,2)
      pz = uz * dd(:,:,:,3)
      deallocate(dd, ux, uy, uz)
    end if
    if (self%mhd) then
      call self%symmetric(bx,idx%bx)
      call self%symmetric(by,idx%by)
      call self%zero_gradient(bz,idx%bz)
      call self%div_free_lower (bx, by, bz)
    end if
    deallocate (d1, e1, p1, pg,e2)
    call trace%end()
  else
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    !
    if (self%patch%kind(1:6)/='ramses') then
      eth => e
    else
      allocate(eth(gn(1),gn(2),gn(3)))
      allocate(emag(gn(1),gn(2),gn(3)))
      eth = e - &
           0.5*(px**2 + &
                py**2 + &
                pz**2)/ d
      do iz=1,gn(3)-1
        do iy=1,gn(2)-1
          do ix=1,gn(1)-1
            emag(ix,iy,iz) = &
                             ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                             (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                             (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
          end do
        end do
      end do
      emag(gn(1),:,:)=emag(gn(1)-1,:,:)
      emag(:,gn(2),:)=emag(:,gn(2)-1,:)
      emag(:,:,gn(3))=emag(:,:,gn(3)-1)
      eth = eth - 0.125*emag
    end if

    call self%log_extrapolate(d,idx%d)
    call self%log_extrapolate(eth,idx%e)
    call self%zero_gradient(px,idx%px)
    call self%zero_gradient(py,idx%py)
    call self%zero_gradient(pz,idx%pz)
    if (self%mhd) then
      call self%symmetric(bx,idx%bx)
      call self%symmetric(by,idx%by)
      call self%zero_gradient(bz,idx%bz)
      call self%div_free_lower (bx, by, bz)
    end if
  end if
  
  ! --------------------------------------------------------------------------
  ! For Riemann solver, convert back to total energy
  ! --------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    e = eth + &
         0.5*(px**2 + &
              py**2 + &
              pz**2)/ d
    do iz=1,gn(3)-1
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,iz) = &
                          ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                           (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                           (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
        end do
      end do
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    emag(:,:,gn(3))=emag(:,:,gn(3)-1)
    e = e + 0.125*emag
    deallocate(eth,emag)
  end if 

END SUBROUTINE lower1

!===============================================================================
!> Top boundary conditions.  Since the Z is flipped, it is the bottom of the
!> atmosphere.
!===============================================================================
SUBROUTINE lower2 (self, mem)
  USE vector_ops
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, gn
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), pointer::  ux, uy, uz, eth, emag
  real, dimension(:,:,:), allocatable:: pg
  real, dimension(:,:,:,:), pointer:: dd
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1
  real:: d_ratio, e_ratio, ee_ratio, h_scale, p, ee, dlnddz, dlnedz, dlneedz, uux, uuy, uuz
  !-----------------------------------------------------------------------------
  call trace%begin('boundries_t%bc_lower')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
    ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
    gn = self%mesh%gn
    m3 => self%mesh(3)
    z  => m3%r
    d  => mem(:,:,:,idx%d )
    e  => mem(:,:,:,idx%e )
    px => mem(:,:,:,idx%px)
    py => mem(:,:,:,idx%py)
    pz => mem(:,:,:,idx%pz)
    bx => mem(:,:,:,idx%bx)
    by => mem(:,:,:,idx%by)
    bz => mem(:,:,:,idx%bz)
    if (self%patch%kind(1:6)/='ramses') then
      eth => e
    else
      allocate(eth(gn(1),gn(2),gn(3)))
      allocate(emag(gn(1),gn(2),gn(3)))
      eth = e - &
           0.5*(px**2 + &
                py**2 + &
                pz**2)/ d
      !
      do iz=1,gn(3)-1
        do iy=1,gn(2)-1
          do ix=1,gn(1)-1
            emag(ix,iy,iz) = &
                             ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                             (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                             (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
          end do
        end do
      end do
      emag(gn(1),:,:)=emag(gn(1)-1,:,:)
      emag(:,gn(2),:)=emag(:,gn(2)-1,:)
      emag(:,:,gn(3))=emag(:,:,gn(3)-1)
      eth = eth - 0.125*emag
    end if

  if (self%first_lower.or..not.deriv_bnd) then
    !
    if (self%first_lower) then
      self%first_lower = .false.
      if (d_bot /= 0.0) self%d_bot = d_bot
      if (e_bot /= 0.0) self%e_bot = e_bot
    end if
    !---------------------------------------------------------------------------
    if (scale) then
      do iy=lb(2),ub(2)
      do ix=lb(1),ub(1)
        dlnddz = log(d(ix,iy,li(3)+3)/d(ix,iy,li(3)+2))/(z(li(3)+2)-z(li(3)+3))
        dlnedz = log(eth(ix,iy,li(3)+3)/eth(ix,iy,li(3)+2))/(z(li(3)+2)-z(li(3)+3))
        dlneedz = log((eth(ix,iy,li(3)+3)/d(ix,iy,li(3)+3)) &
                      *d(ix,iy,li(3)+2)/eth(ix,iy,li(3)+2))/(z(li(3)+2)-z(li(3)+3))
        !-------------------------------------------------------------------------
        ! We assume that the iz index of the last free values is li(3), except
        ! for pz, for which it us li(3)+1.  Hence, if there are 2 guard zones, we
        ! set 3 values for pz, and 2 values for the other variables
        !-------------------------------------------------------------------------
        do iz=lo(3),lb(3),-1
              d_ratio = exp((z(li(3))-z(iz))*dlnddz)
              e_ratio = exp((z(li(3))-z(iz))*dlnedz)
              ee_ratio = exp((z(li(3))-z(iz))*dlneedz)
              !print *, 'd ratio', d_ratio, 1/d_ratio,  d(ix,iy,li(3))*d_ratio
              !print *, 'e ratio', e_ratio, 1/e_ratio, e(ix,iy,li(3))*e_ratio
              d(ix,iy,iz)  =  d(ix,iy,li(3))*d_ratio
              
              ee  =  (eth(ix,iy,li(3))/d(ix,iy,li(3)))*ee_ratio*d(ix,iy,iz)
              !
              eth(ix,iy,iz)  = ee! e(ix,iy,li(3))*e_ratio
              px(ix,iy,iz) = px(ix,iy,li(3))*d_ratio
              py(ix,iy,iz) = py(ix,iy,li(3))*d_ratio
              pz(ix,iy,iz) = pz(ix,iy,li(3))*d_ratio
              if (ee<eos%ee_min) &
              print *, 'ee < ee_min', self%id, ix, iy, iz, ee, eos%ee_min, ee_ratio, d(ix,iy,iz), pz(ix,iy,li(3))
            if ((ix==65).and.(iy==65)) then
          print *, 'z', iz
          print *, 'density', d_ratio, d(ix,iy,iz)
          print *, 'e', e_ratio, eth(ix,iy,iz)
          print *, 'ee', ee_ratio, ee*d(ix,iy,iz), ee
          print *, 'px', px(ix,iy,iz)
          print *, 'py', py(ix,iy,iz)
          print *, 'pz', pz(ix,iy,iz)
            end if
        end do
      end do
      end do
    else
      allocate (ux(size(d,1),size(d,2),size(d,3)))
      allocate (uy(size(d,1),size(d,2),size(d,3)))
      allocate (uz(size(d,1),size(d,2),size(d,3)))
      if (self%patch%kind(1:6)/='ramses') then
        allocate (dd(size(d,1),size(d,2),size(d,3),3))
        dd=down1(d)
        ux = px/dd(:,:,:,1)
        uy = py/dd(:,:,:,2)
        uz = pz/dd(:,:,:,3)
      else
        ux = px/d(:,:,:)
        uy = py/d(:,:,:)
        uz = pz/d(:,:,:)
      end if
      call self%log_extrapolate(d,idx%d)
      call self%log_extrapolate(eth,idx%e)
      call self%zero_gradient(ux,idx%px)
      call self%zero_gradient(uy,idx%py)
      call self%zero_gradient(uz,idx%pz)
      if (self%patch%kind(1:6)/='ramses') then
        dd=down1(d)
        px = ux * dd(:,:,:,1)
        py = uy * dd(:,:,:,2)
        pz = uz * dd(:,:,:,3)
      else
        px = ux * d(:,:,:)
        py = uy * d(:,:,:)
        pz = uz * d(:,:,:)
      end if
      
      deallocate(ux, uy, uz)
      if (self%patch%kind(1:6)/='ramses') deallocate(dd)
    end if
    !
    if (self%mhd) then
      !call self%zero_gradient(bx,idx%bx)
      !call self%zero_gradient(by,idx%by)
      call self%symmetric(bx,idx%bx)
      call self%symmetric(by,idx%by)
      call self%zero_gradient(bz,idx%bz)
      call self%div_free_lower (bx, by, bz)
    end if
    !deallocate (d1, e1, p1, pg)
  else
    !
    call self%log_extrapolate(d,idx%d)
    call self%log_extrapolate(eth,idx%e)
    call self%symmetric(px,idx%px)
    call self%symmetric(py,idx%py)
    call self%symmetric(pz,idx%pz)
    if (self%mhd) then
      call self%symmetric(bx,idx%bx)
      call self%symmetric(by,idx%by)
      call self%zero_gradient(bz,idx%bz)
     call self%div_free_lower (bx, by, bz)
    end if
  end if
  ! --------------------------------------------------------------------------
  ! For Riemann solver, convert back to total energy
  ! --------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    e = eth + &
         0.5*(px**2 + &
              py**2 + &
              pz**2)/ d
    do iz=1,gn(3)-1
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,iz) = &
                          ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                           (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                           (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2)
        end do
      end do
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    emag(:,:,gn(3))=emag(:,:,gn(3)-1)
    e = e + 0.125*emag
    deallocate(eth,emag)
  end if 
  call trace%end()
END SUBROUTINE lower2

!===============================================================================
!> Top boundary conditions.  Since the Z is flipped, it is the bottom of the
!> atmosphere.
!===============================================================================
SUBROUTINE lower_trap (self, mem)
  USE vector_ops
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, gn, n, ng
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd, u_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dlnpdlnr_E, dlnpdE_r
  real, dimension(:,:),   pointer:: p_bot, pb
  real(kind=KindScalarVar) :: ct, beta_z, ee_bot, d_bot_loc, e_bot_loc
  real(kind=KindScalarVar) :: F1,F2,F3,F4,F5,hh,dddt_add, dddt_b, dedt_b
  real(kind=KindScalarVar) :: u_sm, dd, deth
  real, dimension(:,:,:), pointer::  ux, uy, uz, eth, emag
  real, dimension(:,:,:), allocatable:: pg
  !real, dimension(:,:,:,:), pointer:: dd
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1
  real:: d_ratio, e_ratio, ee_ratio, h_scale, p, ee, dlnddz, dlnedz, dlneedz, uux, uuy, uuz
  !-----------------------------------------------------------------------------
  call trace%begin('boundries_t%bc_lower')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  gn = self%mesh%gn
  ng = self%mesh%ng
  m3 => self%mesh(3)
  z  => m3%r
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  
  if (self%patch%kind(1:6)/='ramses') then
    eth => e
  else
    allocate(eth(gn(1),gn(2),li(3)+2))
    allocate(emag(gn(1),gn(2),li(3)+2))
    eth = e(:,:,:li(3)+2) - &
         0.5*(px(:,:,:li(3)+2)**2 + &
              py(:,:,:li(3)+2)**2 + &
              pz(:,:,:li(3)+2)**2)/ d(:,:,:li(3)+2)
    !
    do iz=1,li(3)+2
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,iz) = &
                          ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                           (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                           (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2) * 0.125
        end do
      end do
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    eth = eth - emag
  end if
   !
  if (self%first_lower) then
    self%first_lower = .false.
    if (d_bot /= 0.0) self%d_bot = d_bot
    if (e_bot /= 0.0) self%e_bot = e_bot
  !  print *, self%d_bot, self%e_bot
  end if
  !-----------------------------------------------------------------------------
  ! First, do just for lo(3)
  ! ----------------------------------------------------------------------------
  do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      dlnddz = log(d(ix,iy,li(3)+2)/d(ix,iy,li(3)))/(z(li(3))-z(li(3)+2))
      dlnedz = log(eth(ix,iy,li(3)+2)/eth(ix,iy,li(3)))/(z(li(3))-z(li(3)+2))
      dlneedz = log((eth(ix,iy,li(3)+2)/d(ix,iy,li(3)+2)) &
                    *d(ix,iy,li(3))/eth(ix,iy,li(3)))/(z(li(3))-z(li(3)+2))
      !print *, 'dlnddz dlnedz', dlnddz, dlnedz, dlneedz
      !-------------------------------------------------------------------------
      ! We assume that the iz index of the last free values is li(3), except
      ! for pz, for which it us li(3)+1.  Hence, if there are 2 guard zones, we
      ! set 3 values for pz, and 2 values for the other variables
      !-------------------------------------------------------------------------
    !do iz=lo(3),lb(3),-1
    iz = lo(3)
        d_ratio = exp((z(li(3))-z(iz))*dlnddz)
        e_ratio = exp((z(li(3))-z(iz))*dlnedz)
        ee_ratio = exp((z(li(3))-z(iz))*dlneedz)
        !print *, 'd ratio', d_ratio, 1/d_ratio,  d(ix,iy,li(3))*d_ratio
        !print *, 'e ratio', e_ratio, 1/e_ratio, e(ix,iy,li(3))*e_ratio
        d(ix,iy,iz)  =  d(ix,iy,li(3))*d_ratio
        
        !ee  =  (eth(ix,iy,li(3))/d(ix,iy,li(3)))*ee_ratio*d(ix,iy,iz)
        !
        eth(ix,iy,iz)  =  e(ix,iy,li(3))*e_ratio
        px(ix,iy,iz) = px(ix,iy,li(3))*d_ratio
        py(ix,iy,iz) = py(ix,iy,li(3))*d_ratio
        pz(ix,iy,iz) = pz(ix,iy,li(3))*d_ratio
      !end do
    end do
  end do
  ! ----------------------------------------------------------------------------
  ! Magnetic field boundaries -- might need extra work!
  ! ----------------------------------------------------------------------------
  if (self%mhd) then
    call self%symmetric(bx,idx%bx)
    call self%symmetric(by,idx%by)
    call self%zero_gradient(bz,idx%bz)
    call self%div_free_lower (bx, by, bz)
  end if
  ! ----------------------------------------------------------------------------
  ! Add a pressure trap at lo(3)
  ! ----------------------------------------------------------------------------
  n = shape(d)
  ct = 1.0/t_bdry                                                               ! characteristic damping constant
  !ct = bottom_bc_rate!1.0/t_bdry                                                               ! characteristic damping constant
  
  !  --- allocate scratch arrays ---  
  allocate(ee_bnd    (n(1),n(2)))
  allocate( d_bnd    (n(1),n(2)))
  allocate( pb       (n(1),n(2)))                                               ! magnetic pressure
  allocate(dlnpdE_r  (n(1),n(2)))
  allocate(dlnpdlnr_E(n(1),n(2)))
  allocate (p_bot(size(d,1),size(d,2)))
  allocate (u_bnd(size(d,1),size(d,2)))

  d_bnd  = d(:,:,lo(3))
  ee_bnd = eth(:,:,lo(3))/d_bnd
  !
  ! ----------------------------------------------------------------------------
  ! Compute the pressure pg at the bottom BC, and what we would like it to be
  ! ----------------------------------------------------------------------------
  allocate (d1(size(d,1),size(d,2),2))
  allocate (e1(size(d,1),size(d,2),2))
  allocate (p1(size(d,1),size(d,2),2))
  d1(:,:,1) = self%d_bot                                                        ! mean constant density in the simulation box
  d1(:,:,2) = d(:,:,lo(3))                                                      ! current density
  e1(:,:,1) = self%e_bot                                                        ! mean constant energy density in the simulation box
  e1(:,:,2) = eth(:,:,lo(3))                                                    ! current energy density

  call eos%lookup_table (shape(p1), pg=p1, d=d1, e=e1, label='lower trap', id=self%id)
  !
  p_bot=bottom_bc_rate*p1(:,:,1) + (1.0-bottom_bc_rate)*p1(:,:,2)
  deallocate (d1, e1, p1)
  ! ----------------------------------------------------------------------------
  ! Get gas flow patterns
  ! ----------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    u_bnd=pz(:,:,li(3))/d(:,:,li(3))
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

  iz=lo(3)
  do iy=1,gn(2)-1
    do ix=1,gn(1)-1
      pb(ix,iy) = &
                ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                 (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                 (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2) * 0.125
    end do
  end do
  pb(gn(1),:)=pb(gn(1)-1,:)
  pb(:,gn(2))=pb(:,gn(2)-1)
  
  do iy=lb(2), ub(2)
    do ix = lb(1), ub(1)
      F1 = dlnpdlnr_E(ix,iy)                                                    ! EOS table dlnpdlnr_E
      F2 = dlnpdE_r(ix,iy)                                                      ! EOS table dlnpdE_r
      hh = ee_bnd(ix,iy)+p_bot(ix,iy)/d_bnd(ix,iy)                              ! E in notes
      F5 = F1 + F2*(hh-ee_bnd(ix,iy))                                           ! H in notes
      F3 = (F1 - F2*ee_bnd(ix,iy))/F5                                           ! denominator of F3 and F4
      F4 = F2/F5
      
      !u_sm = (u_bnd(ix-1,iy-1) + u_bnd(ix,iy-1) + u_bnd(ix+1,iy-1)  + &
      !        u_bnd(ix-1,iy  ) + u_bnd(ix,iy  ) + u_bnd(ix+1,iy+1)  + &
      !        u_bnd(ix-1,iy+1) + u_bnd(ix,iy  ) + u_bnd(ix+1,iy+1)) / 9.0
      u_sm = u_bnd(ix,iy)
      ! ------------------------------------------------------------------------
      ! For outflows
      !-------------------------------------------------------------------------
      if (u_sm<=0.) then
        dddt_add = F3*((self%d_bot-d(ix,iy,lo(3)))*ct) + &
                   F4*((self%e_bot-e(ix,iy,lo(3)))*ct)
        d  (ix,iy,lo(3)) = d  (ix,iy,lo(3)) +  dddt_add      ! * self%patch%dtime
        eth(ix,iy,lo(3)) = eth(ix,iy,lo(3)) + (dddt_add * hh)! * self%patch%dtime
      ! ------------------------------------------------------------------------
      ! For inflows
      ! ------------------------------------------------------------------------
      else
        dd = (self%d_bot-d  (ix,iy,lo(3))) * ct! * self%patch%dtime
        deth = (self%e_bot-eth(ix,iy,lo(3))) * ct !* self%patch%dtime
        d   (ix,iy,lo(3)) = d   (ix,iy,lo(3)) + dd
        eth (ix,iy,lo(3)) = eth (ix,iy,lo(3)) + deth
      end if
      ! ------------------------------------------------------------------------
      ! Density and energy perturbations from magnetic pressure
      ! ------------------------------------------------------------------------
      dddt_B = -ct*d(ix,iy,lo(3))*pb_fact*pb(ix,iy)/p_bot(ix,iy)/F5
      dedt_B =  hh*dddt_B
      d  (ix,iy,lo(3)) = d  (ix,iy,lo(3)) + dddt_B! * self%patch%dtime
      eth(ix,iy,lo(3)) = eth(ix,iy,lo(3)) + dedt_B! * self%patch%dtime
    end do
  end do

  deallocate(ee_bnd, d_bnd, dlnpdE_r, dlnpdlnr_E, p_bot, u_bnd, pb)

  !-----------------------------------------------------------------------------
  ! Now, extrapolate to lb(3)
  ! ----------------------------------------------------------------------------
  do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      dlnddz = log(d(ix,iy,lo(3)+1)/d(ix,iy,lo(3)))/(z(lo(3))-z(lo(3)+1))
      dlnedz = log(eth(ix,iy,lo(3)+1)/eth(ix,iy,lo(3)))/(z(lo(3))-z(lo(3)+1))
      dlneedz = log((eth(ix,iy,lo(3)+1)/d(ix,iy,lo(3)+1)) &
                    *d(ix,iy,lo(3))/eth(ix,iy,lo(3)))/(z(lo(3))-z(lo(3)+1))
      !-------------------------------------------------------------------------
      ! We assume that the iz index of the last free values is li(3), except
      ! for pz, for which it us li(3)+1.  Hence, if there are 2 guard zones, we
      ! set 3 values for pz, and 2 values for the other variables
      !-------------------------------------------------------------------------
    do iz=lo(3)-1,lb(3),-1
        d_ratio = exp((z(lo(3))-z(iz))*dlnddz)
        e_ratio = exp((z(lo(3))-z(iz))*dlnedz)
        ee_ratio = exp((z(lo(3))-z(iz))*dlneedz)
        !print *, 'd ratio', d_ratio, 1/d_ratio,  d(ix,iy,li(3))*d_ratio
        !print *, 'e ratio', e_ratio, 1/e_ratio, e(ix,iy,li(3))*e_ratio
        d(ix,iy,iz)  =  d(ix,iy,lo(3))*d_ratio
        
        ee  =  (eth(ix,iy,lo(3))/d(ix,iy,lo(3)))*ee_ratio*d(ix,iy,iz)
        !
        eth(ix,iy,iz)  = ee! e(ix,iy,li(3))*e_ratio
        px(ix,iy,iz) = px(ix,iy,lo(3))*d_ratio
        py(ix,iy,iz) = py(ix,iy,lo(3))*d_ratio
        pz(ix,iy,iz) = pz(ix,iy,lo(3))*d_ratio
      end do
    end do
  end do


  ! --------------------------------------------------------------------------
  ! For Riemann solver, convert back to total energy
  ! --------------------------------------------------------------------------
  if (self%patch%kind(1:6)=='ramses') then
    e(:,:,:lo(3)) = eth(:,:,:lo(3)) + &
         0.5*(px(:,:,:lo(3))**2 + &
              py(:,:,:lo(3))**2 + &
              pz(:,:,:lo(3))**2)/ d(:,:,:lo(3))
    do iz=1,lo(3)
      do iy=1,gn(2)-1
        do ix=1,gn(1)-1
          emag(ix,iy,iz) = &
                          ((bx(ix,iy,iz)+bx(ix+1,iy,iz))**2 + &
                           (by(ix,iy,iz)+by(iy,iy+1,iz))**2 + &
                           (bz(ix,iy,iz)+bz(iz,iy,iz+1))**2) * 0.125
        end do
      end do
    end do
    emag(gn(1),:,:)=emag(gn(1)-1,:,:)
    emag(:,gn(2),:)=emag(:,gn(2)-1,:)
    e(:,:,:lo(3)) = e(:,:,:lo(3)) + emag(:,:,:lo(3))
    deallocate(eth,emag)
  end if 
  call trace%end()
END SUBROUTINE lower_trap

!===============================================================================
!> Top boundary conditions.  Since the Z is flipped, it is the bottom of the
!> atmosphere.
!===============================================================================
SUBROUTINE lower_simple (self, mem)
  USE vector_ops
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, gn, n, ng
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd, u_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dlnpdlnr_E, dlnpdE_r
  real, dimension(:,:),   pointer:: p_bot, pb
  real(kind=KindScalarVar) :: ct, beta_z, ee_bot, d_bot_loc, e_bot_loc
  real(kind=KindScalarVar) :: F1,F2,F3,F4,F5,hh,dddt_add, dddt_b, dedt_b
  real(kind=KindScalarVar) :: u_sm, dd, deth
  real, dimension(:,:,:), pointer::  ux, uy, uz, eth, emag
  real, dimension(:,:,:), allocatable:: pg
  !real, dimension(:,:,:,:), pointer:: dd
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1
  real:: d_ratio, e_ratio, ee_ratio, h_scale, p, ee, dlnddz, dlnedz, dlneedz, uux, uuy, uuz
  !-----------------------------------------------------------------------------
  call trace%begin('boundries_t%bc_lower')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  gn = self%mesh%gn
  ng = self%mesh%ng
  m3 => self%mesh(3)
  z  => m3%r
  d  => mem(:,:,:,idx%d )
  e  => mem(:,:,:,idx%e )
  px => mem(:,:,:,idx%px)
  py => mem(:,:,:,idx%py)
  pz => mem(:,:,:,idx%pz)
  bx => mem(:,:,:,idx%bx)
  by => mem(:,:,:,idx%by)
  bz => mem(:,:,:,idx%bz)
  call self%zero_gradient(d,idx%d)
  call self%zero_gradient(e,idx%e)
  call self%zero_gradient(px,idx%px)
  call self%zero_gradient(py,idx%py)
  call self%zero_gradient(pz,idx%pz)
  call self%zero_gradient(bx,idx%bx)
  call self%zero_gradient(by,idx%by)
  call self%zero_gradient(bz,idx%bz)

  call trace%end()
END SUBROUTINE lower_simple

! ------------------------------------------------------------------------------
! Lower damp module develhere
! ------------------------------------------------------------------------------
SUBROUTINE lower_damp(self, mem)
  class(boundaries_t):: self
  ! ----------------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), pointer:: dddt, dedt, dpdtx, dpdty, dpdtz, dbdtx, dbdty
  real, dimension(:,:),   pointer:: pg, dp, bxu, byu, bzu,p_bot, pe_ratio
  real, dimension(:,:,:), pointer:: ux, uy, uz
  real, dimension(:,:,:,:), pointer:: dd
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dpdeecd, dpddcee, dpddcs
  real(kind=KindScalarVar) :: ct, beta_z, ee_bot, d_bot_loc, e_bot_loc
  real(kind=KindScalarVar) :: d_inv, ee, hh, dpddce, dpdecd, dpfluct
  integer :: ix, iy, iz
  integer, save::itimer=0
  ! ............................................................................
  call trace%begin('boundaries_t%lower_damp', itimer=itimer)
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  m3 => self%mesh(3)
  z  => m3%r
  d  => mem(:,:,:,idx%d ,1)
  e  => mem(:,:,:,idx%e ,1)
  px => mem(:,:,:,idx%px,1)
  py => mem(:,:,:,idx%py,1)
  pz => mem(:,:,:,idx%pz,1)
  bx => mem(:,:,:,idx%bx,1)
  by => mem(:,:,:,idx%by,1)
  bz => mem(:,:,:,idx%bz,1)
  !
  dddt  => mem(:,:,:,idx%d ,2)
  dedt  => mem(:,:,:,idx%e ,2)
  dpdtx => mem(:,:,:,idx%px,2)
  dpdty => mem(:,:,:,idx%py,2)
  dpdtz => mem(:,:,:,idx%pz,2)
  dbdtx => mem(:,:,:,idx%bx,2)
  dbdty => mem(:,:,:,idx%by,2)
  !
  n = shape(d)
  ct = 1.0/t_bdry                                                               ! characteristic damping constant
  !ct = bottom_bc_rate!1.0/t_bdry                                                               ! characteristic damping constant
  allocate(ee_bnd(n(1),n(2)))
  allocate( d_bnd(n(1),n(2)))
  allocate(    pg(n(1),n(2)))
  allocate(dpdeecd(n(1),n(2)))
  allocate(dpddcee(n(1),n(2)))
  allocate(dpddcs(n(1),n(2)))
  !
  !allocate(bxu(n(1),n(2)))
  !allocate(byu(n(1),n(2)))
  allocate(bzu(n(1),n(2)))
  !
  allocate (d1(size(d,1),size(d,2),1))
  allocate (e1(size(d,1),size(d,2),1))
  allocate (p1(size(d,1),size(d,2),1))
  allocate (p_bot(size(d,1),size(d,2)))
  allocate (pe_ratio(size(d,1),size(d,2)))
  d1 = self%d_bot
  e1 = self%e_bot
  call eos%lookup_table (shape(p1), pg=p1, d=d1, e=e1, label='lower damp bnd', id=self%id)
  p_bot=p1(:,:,1)
  pe_ratio = p_bot/e_bot
  deallocate (d1, e1, p1)

  ! ----------------------------------------------------------------------------
  ! Get pressure and thermodynamic derivatives for pressure variation, scaling 
  ! factors for the damping of pressure fluctuations.
  ! dP/dee_r with constant density
  ! dP/dr_ee with constant ee=e/d
  ! dP/dr_s  with constant entropy
  ! ----------------------------------------------------------------------------
  d_bnd  = d(:,:,li(3))
  ee_bnd = e(:,:,li(3))
  !call eos%ddeos(ee_bnd,d_bnd,pg,dpdeecd,dpddcee,dpddcs)
  ! ----------------------------------------------------------------------------
  ! Subtract magnetic pressure part 
  ! ----------------------------------------------------------------------------
  !do ix = lb(1), ub(1)-1
  !  bxu(ix,:) = 0.5*(bx(ix,:,li(3))+bx(ix+1,:,li(3)))
  !end do
  !bxu(ub(1),:) = bx(ub(1),:,li(3))
  !
  !do iy = lb(2), ub(2)-1
  !  byu(:,iy) = 0.5*(by(:,iy,li(3))+by(:,iy+1,li(3)))
  !end do
  !byu(:,ub(2))=by(:,ub(2),li(3))
  !!
  bzu = 0.5*(bz(:,:,li(3)+1)+bz(:,:,li(3)))
  !
  !dp = pg - 0.5*(bxu**2+byu**2+bzu**2)
  !
  do iy=lb(2), ub(2)
    do ix=lb(1), ub(1)
      beta_z = 2.*pg(ix,iy)/max(bzu(ix,iy)**2,1e-20)
      d_bot_loc = self%d_bot
      e_bot_loc = self%e_bot
      if (beta_z < beta_lim.and. spot_lower /= 0.0) then
        ee_bot  = e_bot_loc/d_bot_loc * &
                 (1.0-spot_lower*beta_lim/max(beta_z,beta_lim))
        d_bot_loc = (p_bot(ix,iy)-pb_spot)/pe_ratio(ix,iy)/ee_bot
        e_bot_loc = ee_bot*d_bot_loc
      end if
      ! ------------------------------------------------------------------------
      ! For outgoing flows, damp out pressure fluctuations to create a
      ! pressure node
      ! ------------------------------------------------------------------------
      if (pz(ix,iy,li(3)+2)<0.0) then
        d_inv = 1.0/d(ix,iy,li(3))
        ! Internal energy and enthalpy
        ee = e (ix,iy,li(3)) * d_inv
        hh = pg(ix,iy)       * d_inv + ee
        ! ----------------------------------------------------------------------
        ! The table uses ee=e/d, but the MHD solver computes e,
        ! so convert the derivatives dP/dr_ee and
        ! dP/dee_d into dP/dd_e and dP/de_d.
        ! ----------------------------------------------------------------------
        dpddce=(dpddcee(ix,iy) - dpdeecd(ix,iy)*ee*d_inv)
        dpdecd= dpdeecd(ix,iy) * d_inv
        ! ----------------------------------------------------------------------
        ! Compute the (negative) pressure fluctation using
        ! dP = (dP/dd)_e * dd + (dP/de)_d * de
        ! from the density and energy fluctuations dd and de
        ! decrasing rbot and ebot according to change in p needed to
        ! maintain TOTAL pressure equivalent to ebot and dbot.
        ! ----------------------------------------------------------------------
        dpfluct = dpddce*((d_bot_loc-d(ix,iy,li(3)))*ct-dddt(ix,iy,li(3))) + &
                  dpdecd*((e_bot_loc-e(ix,iy,li(3)))*ct-dedt(ix,iy,li(3)))
        ! ----------------------------------------------------------------------
        ! Now damp this pressure fluctuation dP by
        ! following an adiabatic curve: compute
        ! the equivalent adiabatic density fluctuation
        ! using the relation
        ! dP_s = (dP/dd)_s * dd_s
        ! ----------------------------------------------------------------------
        dpfluct = dpfluct/dpddcs(ix,iy)
        ! ----------------------------------------------------------------------
        ! Apply the damping to d and e; de and dd are related by de=(P+e)/d*dd
        ! ----------------------------------------------------------------------
        dddt(ix,iy,li(3)) = dddt(ix,iy,li(3))+dpfluct
        dedt(ix,iy,li(3)) = dedt(ix,iy,li(3))+dpfluct*hh
      ! ------------------------------------------------------------------------
      ! For incoming flows, damp internal energy and density deviations and
      ! damp horizontal motion
      ! ------------------------------------------------------------------------
      
      else
        dedt(ix,iy,li(3))  = (e_bot_loc-e(ix,iy,li(3)))*ct
        dddt(ix,iy,li(3))  = (d_bot_loc-d(ix,iy,li(3)))*ct
        !
        dpdtx(ix,iy,li(3)) = dpdtx(ix,iy,li(3))-px(ix,iy,li(3))*ct
        dpdty(ix,iy,li(3)) = dpdty(ix,iy,li(3))-py(ix,iy,li(3))*ct
      end if
    end do
  end do
  ! ----------------------------------------------------------------------------
  ! Make sure that nothing bad happens in guard zones
  ! ----------------------------------------------------------------------------
  dddt (:,:,lb(3):lo(3))=0.0
  dedt (:,:,lb(3):lo(3))=0.0
  dpdtx(:,:,lb(3):lo(3))=0.0
  dpdty(:,:,lb(3):lo(3))=0.0
  dpdtz(:,:,lb(3):li(3))=0.0
  !
  deallocate(ee_bnd, d_bnd, pg, dpdeecd, dpddcee,dpddcs, bzu, p_bot, pe_ratio)
  !deallocate(bxu,byu,bzu)

  call trace%end(itimer)
END SUBROUTINE lower_damp

! ------------------------------------------------------------------------------
! Lower damp module develhere
! ------------------------------------------------------------------------------
SUBROUTINE ddt_boundary_bot(self, mem)
  class(boundaries_t):: self
  ! ----------------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), pointer:: dddt, dedt, dpdtx, dpdty, dpdtz, dbdtx, dbdty, dbdtz
  real, dimension(:,:),   pointer:: p_bot
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd, u_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dlnpdlnr_E, dlnpdE_r
  real(kind=KindScalarVar) :: ct, beta_z, ee_bot, d_bot_loc, e_bot_loc
  real(kind=KindScalarVar) :: F1,F2,F3,F4,F5,hh,dddt_add, pb, dddt_b, dedt_b
  real(kind=KindScalarVar) :: u_sm
  integer :: ix, iy, iz, aa(2)
  integer, save::itimer=0
  ! ............................................................................
  call trace%begin('boundaries_t%ddt_boundary_bot', itimer=itimer)
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  d  => mem(:,:,:,idx%d ,1)
  e  => mem(:,:,:,idx%e ,1)
  px => mem(:,:,:,idx%px,1)
  py => mem(:,:,:,idx%py,1)
  pz => mem(:,:,:,idx%pz,1)
  bx => mem(:,:,:,idx%bx,1)
  by => mem(:,:,:,idx%by,1)
  bz => mem(:,:,:,idx%bz,1)
  error stop "ddt bnd top"
  !
  dddt  => mem(:,:,:,idx%d ,2)
  dedt  => mem(:,:,:,idx%e ,2)
  dpdtx => mem(:,:,:,idx%px,2)
  dpdty => mem(:,:,:,idx%py,2)
  dpdtz => mem(:,:,:,idx%pz,2)
  dbdtx => mem(:,:,:,idx%bx,2)
  dbdty => mem(:,:,:,idx%by,2)
  dbdtz => mem(:,:,:,idx%bz,2)
  !
  n = shape(d)
  ct = 1.0/t_bdry                                                               ! characteristic damping constant
  !ct = bottom_bc_rate!1.0/t_bdry                                                               ! characteristic damping constant
  
  !  --- allocate scratch arrays ---  
  allocate(ee_bnd    (n(1),n(2)))
  allocate( d_bnd    (n(1),n(2)))
  allocate(dlnpdE_r  (n(1),n(2)))
  allocate(dlnpdlnr_E(n(1),n(2)))
  allocate (p_bot(size(d,1),size(d,2)))
  allocate (u_bnd(size(d,1),size(d,2)))
  
   d_bnd = d(:,:,li(3))
  ee_bnd = e(:,:,li(3))/d_bnd
  
  ! ----------------------------------------------------------------------------
  ! Compute the pressure pg at the bottom BC, and what we would like it to be
  ! ----------------------------------------------------------------------------
  allocate (d1(size(d,1),size(d,2),2))
  allocate (e1(size(d,1),size(d,2),2))
  allocate (p1(size(d,1),size(d,2),2))
  d1(:,:,1) = self%d_bot                                                        ! mean constant density in the simulation box
  d1(:,:,2) = d(:,:,li(3))                                                      ! current density
  e1(:,:,1) = self%e_bot                                                        ! mean constant energy density in the simulation box
  e1(:,:,2) = e(:,:,li(3))                                                      ! current energy density

  call eos%lookup_table (shape(p1), pg=p1, d=d1, e=e1, label='ddt bot', id=self%id)
  aa=maxloc(p1(:,:,2))
  associate (w => bottom_bc_rate)
  p_bot=w*p1(:,:,1) + (1.0-w)*p1(:,:,2)

  end associate
  deallocate (d1, e1, p1)
  !
  u_bnd=0.5*(pz(:,:,li(3)+1)+pz(:,:,li(3)+2))/d(:,:,li(3))
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
      hh = ee_bnd(ix,iy)+p_bot(ix,iy)/d_bnd(ix,iy)                                         ! E in notes
      F5 = F1 + F2*(hh-ee_bnd(ix,iy))                                           ! H in notes
      F3 = (F1 - F2*ee_bnd(ix,iy))/F5                                           ! denominator of F3 and F4
      F4 = F2/F5
      
      u_sm = (u_bnd(ix-1,iy-1) + u_bnd(ix,iy-1) + u_bnd(ix+1,iy-1)  + &
              u_bnd(ix-1,iy  ) + u_bnd(ix,iy  ) + u_bnd(ix+1,iy+1)  + &
              u_bnd(ix-1,iy+1) + u_bnd(ix,iy  ) + u_bnd(ix+1,iy+1)) / 9.0
      ! ------------------------------------------------------------------------
      ! For outflows
      !-------------------------------------------------------------------------
      if (u_sm<=0.) then
      !if (pz(ix,iy,li(3)+1)<=0.) then
        dddt_add = F3*((self%d_bot-d(ix,iy,li(3)))*ct-dddt(ix,iy,li(3))) + &
                   F4*((self%e_bot-e(ix,iy,li(3)))*ct-dedt(ix,iy,li(3)))
        dddt(ix,iy,li(3)) = dddt(ix,iy,li(3)) + dddt_add
        dedt(ix,iy,li(3)) = dedt(ix,iy,li(3)) + dddt_add * hh
      ! ------------------------------------------------------------------------
      ! For inflows
      ! ------------------------------------------------------------------------
      else
        dddt(ix,iy,li(3)) = (self%d_bot-d(ix,iy,li(3))) * ct
        dedt(ix,iy,li(3)) = (self%e_bot-e(ix,iy,li(3))) * ct
        !dddt_add = F3*((self%d_bot-d(ix,iy,li(3)))*ct-dddt(ix,iy,li(3))) + &
        !           F4*((self%e_bot-e(ix,iy,li(3)))*ct-dedt(ix,iy,li(3)))
        !dddt(ix,iy,li(3)) = dddt(ix,iy,li(3)) + dddt_add
        !dedt(ix,iy,li(3)) = dedt(ix,iy,li(3)) + dddt_add * hh
      end if
      ! ------------------------------------------------------------------------
      ! Density and energy perturbations from magnetic pressure
      ! ------------------------------------------------------------------------
      pb = ((0.5*(bx(ix,iy,li(3))+bx(ix+1,iy  ,li(3)  )))**2  + &
                (0.5*(by(ix,iy,li(3))+bx(ix  ,iy+1,li(3)  )))**2  + &
                (0.5*(bz(ix,iy,li(3))+bx(ix  ,iy  ,li(3)+1)))**2) * 0.5
      dddt_B = -ct*d(ix,iy,li(3))*pb_fact*pb/p_bot(ix,iy)/F5
      dedt_B =  hh*dddt_B
      dddt(ix,iy,li(3)) = dddt(ix,iy,li(3)) + dddt_B
      dedt(ix,iy,li(3)) = dedt(ix,iy,li(3)) + dedt_B
    end do
  end do

  ! ----------------------------------------------------------------------------
  ! Make sure that nothing bad happens in guard zones
  ! ----------------------------------------------------------------------------
  dddt (:,:,lb(3):lo(3))=0.0
  dedt (:,:,lb(3):lo(3))=0.0
  dpdtx(:,:,lb(3):lo(3))=0.0
  dpdty(:,:,lb(3):lo(3))=0.0
  dpdtz(:,:,lb(3):li(3))=0.0

  if (self%mhd) then
    dbdtx(:,:,lb(3):lo(3))=0.0
    dbdty(:,:,lb(3):lo(3))=0.0
    dbdtz(:,:,lb(3):li(3))=0.0
  end if

  deallocate(ee_bnd, d_bnd, dlnpdE_r, dlnpdlnr_E, p_bot, u_bnd)

  
  call trace%end(itimer)
END SUBROUTINE ddt_boundary_bot

SUBROUTINE ddt_boundary_top (self, mem)
  class(boundaries_t):: self
  ! ----------------------------------------------------------------------------
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, pg
  real, dimension(:,:,:), pointer:: dddt, dedt, dpdtx, dpdty, dpdtz, dbdtx, dbdty, dbdtz
  real, dimension(:,:),   pointer:: p_bot
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd, u_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dlnpdlnr_E, dlnpdE_r
  real(kind=KindScalarVar) :: h_scale
  real(kind=KindScalarVar) :: F1,F2,F3,F4,F5,hh,dddt_add, pb, dddt_b, dedt_b
  real(kind=KindScalarVar) :: u_sm, ds(3)
  integer :: ix, iy, iz, aa(2)
  integer, save::itimer=0
  ! ............................................................................
  call trace%begin('boundaries_t%ddt_boundary_top', itimer=itimer)
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  ds = self%mesh%d
  error stop "ddt top bnd called"
  d  => mem(:,:,:,idx%d ,1)
  e  => mem(:,:,:,idx%e ,1)
  px => mem(:,:,:,idx%px,1)
  py => mem(:,:,:,idx%py,1)
  pz => mem(:,:,:,idx%pz,1)
  bx => mem(:,:,:,idx%bx,1)
  by => mem(:,:,:,idx%by,1)
  bz => mem(:,:,:,idx%bz,1)
  !
  dddt  => mem(:,:,:,idx%d ,2)
  dedt  => mem(:,:,:,idx%e ,2)
  dpdtx => mem(:,:,:,idx%px,2)
  dpdty => mem(:,:,:,idx%py,2)
  dpdtz => mem(:,:,:,idx%pz,2)
  dbdtx => mem(:,:,:,idx%bx,2)
  dbdty => mem(:,:,:,idx%by,2)
  dbdtz => mem(:,:,:,idx%bz,2)

  !---------------------------------------------------------------------------
  ! Boundary conditions need the gas pressure
  !---------------------------------------------------------------------------
  allocate (pg(size(d,1),size(d,2),1))
  d1 => d(:,:,ui(3):ui(3))
  e1 => e(:,:,ui(3):ui(3))
  iz=ui(3)
  call eos%lookup_table (shape(pg), pg=pg, d=d1, e=e1, label='ddt_top', id=self%id)
  do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      h_scale = pg(ix,iy,1)/(d(ix,iy,iz)*self%gravity)                          ! Pressure scale height
      if (py(ix,iy,iz-1) < 0. .or. py(ix,iy,iz) < 0.) then                      ! For inflows
        dddt(ix,iy,iz  ) = dddt(ix,iy,iz  ) &                                   ! let the density
       + (pz(ix,iy,iz  ) -   pz(ix,iy,iz+1)) / ds(3) - &                        ! remove ddyup(pz) part
       - h_scale*pz(ix,iy,iz-1)                                                 ! replace with -pz/H term
      end if
      dddt (ix,iy,uo(3):ub(3)) = 0.
      dpdtx(ix,iy,uo(3):ub(3)) = 0.
      dpdty(ix,iy,uo(3):ub(3)) = 0.
      dpdtz(ix,iy,uo(3):ub(3)) = 0.
      dedt (ix,iy,uo(3):ub(3)) = 0.
      !
      dedt (ix,iy,ui(3)      ) = dddt(ix,iy,iz)*self%ee_top
      if (self%mhd) then
        dbdtx(ix,iy,uo(3):ub(3)) = 0.
        dbdty(ix,iy,uo(3):ub(3)) = 0.
        dbdtz(ix,iy,uo(3):ub(3)) = 0.
      end if
    end do
  end do
  deallocate(pg)
  call trace%end(itimer)
END SUBROUTINE ddt_boundary_top


SUBROUTINE sponge_boundaries(self, mem, lower)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  logical:: lower
  ! ----------------------------------------------------------------------------
  integer, dimension(3):: lb, lo, li, ui, uo, ub, n
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), pointer:: dddt, dedt, dpdtx, dpdty, dpdtz, dbdtx, dbdty
  real, dimension(:,:),   pointer:: pg, dp, bxu, byu, bzu,p_bot, pe_ratio
  real, dimension(:,:,:), pointer:: ux, uy, uz
  real, dimension(:,:,:,:), pointer:: dd
  class(mesh_t), pointer:: m1, m2, m3
  real(kind=KindScalarVar), dimension(:,:), pointer:: ee_bnd, d_bnd
  real(kind=KindScalarVar), dimension(:,:), pointer:: dpdeecd, dpddcee, dpddcs
  real(kind=KindScalarVar) :: asymptote
  integer :: i, j, k
  real(kind=KindScalarVar):: x, y, z, zs
  integer, save::itimer=0
  ! ----------------------------------------------------------------------------
  real(kind=KindScalarVar):: zkill_lo, zkill_hi, zbndy_lo, zbndy_hi, tau_kill
  ! ............................................................................
  if (damp_rate<=0.0) return


  call trace%begin('boundaries_t%sponge_boundaries', itimer=itimer)
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  d  => mem(:,:,:,idx%d ,1)
  e  => mem(:,:,:,idx%e ,1)
  px => mem(:,:,:,idx%px,1)
  py => mem(:,:,:,idx%py,1)
  pz => mem(:,:,:,idx%pz,1)
  bx => mem(:,:,:,idx%bx,1)
  by => mem(:,:,:,idx%by,1)
  bz => mem(:,:,:,idx%bz,1)
  !
  dddt  => mem(:,:,:,idx%d ,2)
  dedt  => mem(:,:,:,idx%e ,2)
  dpdtx => mem(:,:,:,idx%px,2)
  dpdty => mem(:,:,:,idx%py,2)
  dpdtz => mem(:,:,:,idx%pz,2)
  dbdtx => mem(:,:,:,idx%bx,2)
  dbdty => mem(:,:,:,idx%by,2)
  !
  n = shape(d)
  ! ---
  tau_kill = math%pi2*damp_rate

  zbndy_lo = m3%p + m3%r(m3%li) - 0.5 * m3%d
  zkill_lo = zbndy_lo + r_sponge
  zbndy_hi = m3%p + m3%r(m3%ui) + 0.5 * m3%d
  zkill_hi = zbndy_hi - r_sponge
  
  do k=m3%lb,m3%ub
    z = m3%p + m3%r(k)
    zs = m3%p + m3%r(k) - 0.5_8*m3%d
    if (z <= zkill_lo) then
      do j=m2%lb,m2%ub
        do i=m1%lb,m1%ub
          dddt(i,j,k) = self%sponge(dddt(i,j,k), 0.0, tau_kill, zkill_lo, zbndy_lo, z)
          dedt(i,j,k) = self%sponge(dedt(i,j,k), 0.0, tau_kill, zkill_lo, zbndy_lo, z)
          dpdtx(i,j,k) = self%sponge(dpdtx(i,j,k), 0.0, tau_kill, zkill_lo, zbndy_lo, z)
          dpdty(i,j,k) = self%sponge(dpdty(i,j,k), 0.0, tau_kill, zkill_lo, zbndy_lo, z)
        end do
      end do
    end if
    if (zs <= zkill_lo) then
      do j=m2%lb,m2%ub
        do i=m1%lb,m1%ub
          dpdtz(i,j,k) = self%sponge(dpdtz(i,j,k), 0.0, tau_kill, zkill_lo, zbndy_lo, zs)
        end do
      end do
    end if
  end do
  call trace%end(itimer)
END SUBROUTINE sponge_boundaries


SUBROUTINE radiation (self, mem, src, rk, tt, dir, mu)
  class(boundaries_t) :: self
  real, dimension(:,:,:,:),pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt
  real :: mu(3)
  integer:: dir
  logical:: lower
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3
  integer:: ib, i3, i3b, i2, i1, n(4), ui
  real:: h_scale, tau, f_scale, dr
  !-----------------------------------------------------------------------------
  n=shape(mem)
  if (self%is_set(bits%zu))  then
    dr = self%mesh(3)%d/abs(mu(3))
    i3 = self%mesh(3)%uo
    i3b = self%mesh(3)%ub
    ui = self%mesh(3)%ui
    do ib=1,n(4)
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
        ! -- Estimate the boundary Q from an estimate of the total optical depth
        !h_scale = tt(i1,i2,ui)/(scaling%temp*self%gravity)
        !tau = h_scale*rk(i1,i2,ui,ib)/(mu(3))
        tau = dr*0.5*(rk(i1,i2,ui,ib)+rk(i1,i2,i3,ib))
        mem(i1,i2,i3:i3b,ib) = -src(i1,i2,i3,ib)*exp(-tau)
        !mem(i1,i2,i3:i3b,ib) = 0.0
      end do
      end do
    end do
  end if
  if (self%is_set(bits%zl))  then
    i3 = self%mesh(3)%lo
    i3b = self%mesh(3)%lb
    mem(:,:,i3:i3b,:) = 0.0
  end if 
END SUBROUTINE radiation

SUBROUTINE radiation2 (self, mem, src, rk, tt, dir, mu)
  class(boundaries_t) :: self
  real, dimension(:,:,:,:), pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt,mem2
  real   :: mu(3)
  !.....................................................................
  real   :: w, tau, dr
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i3, i3b, i2, i1, n(4), iv, dir, ui
  integer:: lbit(3), ubit(3)
  !-----------------------------------------------------------------------------
  n=shape(mem)
  lbit = [bits%xl, bits%yl, bits%zl]        ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]        ! upper boundary bits -- convenience
  dr = self%mesh(dir)%d/abs(mu(dir))
  !do iv=1,n(4)
  !  mem2 => src(:,:,:,iv)
  !  call self%log_extrapolate(mem2, 1)
  !  mem2 => rk(:,:,:,iv)
  !  call self%log_extrapolate(mem2, 1)
  !end do
  if (self%is_set(lbit(3)) ) then
    i3 =  self%mesh(3)%lo
    i3b = self%mesh(3)%lb
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
        do i1=self%mesh(1)%lb,self%mesh(1)%ub
          do iv=1, n(4)
            mem(i1,i2,i3b:i3,iv) = 0.0
          end do
        end do
      end do
  end if
  if (self%is_set(ubit(3)) ) then
    dr  = self%mesh(3)%d/abs(mu(3))
    i3  = self%mesh(3)%uo
    i3b = self%mesh(3)%ub
    ui  = self%mesh(3)%ui
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        tau = dr*0.5*(rk(i1,i2,ui,iv)+rk(i1,i2,i3,iv))
        mem(i1,i2,i3:i3b,iv) = -src(i1,i2,i3,iv)*exp(-tau)
        !mem(i1,i2,i3:i3b,iv) = 0.0
      end do
      end do
      end do
  end if
  !
  if (self%is_set(lbit(2)) ) then
    i3  = self%mesh(2)%lo
    i3b = self%mesh(2)%lb
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3b:i3,i2,iv) = 0.0
      end do
      end do
      end do
  end if
  if (self%is_set(ubit(2)) ) then
    i3  = self%mesh(2)%uo
    i3b = self%mesh(2)%ub
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3:i3b,i2,iv) = 0.0
      end do
      end do
      end do
  end if
  !
  if (self%is_set(lbit(1)) ) then
    i3  = self%mesh(1)%lo
    i3b = self%mesh(1)%lb
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3b:i3,i1,i2,iv) = 0.0
      end do
      end do
      end do
  end if  
  if (self%is_set(ubit(1)) ) then
    i3  = self%mesh(1)%uo
    i3b = self%mesh(1)%ub
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3:i3b,i1,i2,iv) = 0.0
      end do
      end do
      end do
  end if  
END SUBROUTINE radiation2


!===============================================================================
!> Divergence-free boundaries. We do div_free operations for Z faces and then
!> for X/Y parts in the boundary.
!===============================================================================
SUBROUTINE div_free_upper(self, b1, b2, b3)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: b1, b2, b3
  ! ----------------------------------------------------------------------------
  integer:: iter, i1, i2, i3, i
  integer, dimension(3):: li, ui, lb, ub, lo, uo
  real(8):: ds(3)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (divb_iter<1) return
  call trace%begin('boundaries_t%div_free_upper', itimer=itimer)
  li = self%mesh%li
  ui = self%mesh%ui
  lb = self%mesh%lb
  ub = self%mesh%ub
  lo = self%mesh%lo
  uo = self%mesh%uo
  ds = self%mesh%d
  do iter=1,divb_iter
    call iter_z
    call iter_x
    call iter_y
    call iter_z
  end do
  call trace%end (itimer)
contains
  subroutine iter_x
    !---------------------------------------------------------------------------
    ! Right x-face.  Initially unconstrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=lb(1),ub(1)-1
    do i3=uo(3),ub(3)-1
    do i2=lb(2),ub(2)-1
      b1(i1+1,i2,i3) =                           + b1(i1,i2,i3)  &
                     - (1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3)) &
                     +  1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3))) * ds(1)
      !b1(i1+1,i2,i3) =                                + b1(i1,i2,i3)  &
      !               - b2(i1,i2+1,i3) + b2(i1,i2,i3) &
      !               - b3(i1,i2,i3+1) + b3(i1,i2,i3)
    end do
      b1(i1+1,ub(2),i3) = b1(i1+1,ub(2)-1,i3)
    end do
      b1(i1+1,:,ub(3)) = b1(i1+1,:,ub(3)-1)
    end do
    b1(ub(1),:,:)=b1(ub(1)-1,:,:)
  end subroutine iter_x
  subroutine iter_y
    !---------------------------------------------------------------------------
    ! Right y-face.
    !---------------------------------------------------------------------------
    do i2=lb(2),ub(2)-1
    do i3=uo(3),ub(3)-1
    do i1=lb(1),ub(1)-1
      b2(i1,i2+1,i3) =                           + b2(i1,i2,i3)  &
                     - (1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3)) &
                     +  1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3))) * ds(2)
      !b2(i1,i2+1,i3) =                                + b2(i1,i2,i3)  &
      !               - b3(i1,i2,i3+1) + b3(i1,i2,i3) &
      !               - b1(i1+1,i2,i3) + b1(i1,i2,i3)
    end do
      b2(ub(1),i2+1,i3) = b2(ub(1)-1,i2+1,i3)
    end do
      b2(:,i2+1,ub(3))  = b2(:,i2+1,ub(3)-1)
    end do
      b2(:,ub(2),:) = b2(:,ub(2)-1,:)
  end subroutine iter_y
  subroutine iter_z
    !---------------------------------------------------------------------------
    ! Right z-face: Not known, b3(i1,i2,i3+1)
    !---------------------------------------------------------------------------
    do i3=ui(3),ub(3)-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3+1) =                           + b3(i1,i2,i3)  &
                     - (1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3)) &
                     +  1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3))) * ds(3)
      !b3(i1,i2,i3+1) =                                + b3(i1,i2,i3)  &
      !               - b1(i1+1,i2,i3) + b1(i1,i2,i3) &
      !               - b2(i1,i2+1,i3) + b2(i1,i2,i3)
    end do
      b3(ub(1),i2,i3+1) = b3(ub(1)-1,i2,i3+1)
    end do
      b3(:,ub(2),i3+1)  = b3(:,ub(2)-1,i3+1)
    end do
    b3(:,:,ub(3)) = b3(:,:,ub(3)-1)
  end subroutine iter_z
END subroutine div_free_upper

SUBROUTINE div_free_lower(self, b1, b2, b3)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: b1, b2, b3
  ! ----------------------------------------------------------------------------
  integer:: iter, i1, i2, i3
  integer, dimension(3):: li, ui, lb, ub, lo, uo
  integer:: itimer=0
  real(8):: ds(3)
  !.............................................................................
  if (divb_iter<1) return
  call trace%begin('boundaries_t%div_free_lower', itimer=itimer)
  li = self%mesh%li
  li = self%mesh%li
  ui = self%mesh%ui
  lb = self%mesh%lb
  ub = self%mesh%ub
  lo = self%mesh%lo
  uo = self%mesh%uo
  ds = self%mesh%d
  do iter=1,divb_iter
    call iter_z
    call iter_x
    call iter_y
    call iter_z
  end do
  call trace%end (itimer)
contains
  subroutine iter_x
    !---------------------------------------------------------------------------
    ! Left x-face.  Initially un-constrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=ub(1)-1,lb(1),-1
    do i3=li(3)  ,lb(3),-1
    do i2=lb(2)  ,ub(2)-1
      b1(i1,i2,i3) =            b1(i1+1,i2,i3)                 &
                   + (1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3)) &
                   +  1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3))) * ds(1)
      !b1(i1,i2,i3) =                 b1(i1+1,i2,i3)                  &
      !             + b2(i1,i2+1,i3) - b2(i1,i2,i3) &
      !             + b3(i1,i2,i3+1) - b3(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_x
  subroutine iter_y
    !---------------------------------------------------------------------------
    ! Left y-face.  .  Initially unconstrained: b3(y) along edges
    !---------------------------------------------------------------------------
    do i2=ub(2)-1,lb(2),-1
    do i3=li(3)  ,lb(3),-1
    do i1=lb(1)  ,ub(1)-1
      b2(i1,i2,i3) =            b2(i1,i2+1,i3) +               &
                   + (1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3)) &
                   +  1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3)))  * ds(2)
      !b2(i1,i2,i3) =                 b2(i1,i2+1,i3) +               &
      !             + b3(i1,i2,i3+1) - b3(i1,i2,i3) &
      !             + b1(i1+1,i2,i3) - b1(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_y
  subroutine iter_z
    !---------------------------------------------------------------------------
    ! Left z-face.  Not known: b3(:,:,l3)
    !---------------------------------------------------------------------------
    do i3=li(3),lb(3),-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3) =                 b3(i1,i2,i3+1)                 &
                   + (1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3)) &
                   +  1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3))) * ds(3)
      !b3(i1,i2,i3) =                 b3(i1,i2,i3+1)                 &
      !             + b1(i1+1,i2,i3) - b1(i1,i2,i3) &
      !             + b2(i1,i2+1,i3) - b2(i1,i2,i3)
    end do
    end do
    end do
  end subroutine iter_z
END subroutine div_free_lower

END MODULE boundaries_mod
