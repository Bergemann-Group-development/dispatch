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
  USE model_1d_mod
  USE patch_mod
  USE eos_mod
  USE interpolate_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
  ! --------------------------------------------------------------------------
    real, pointer, contiguous:: grav_bot(:,:,:)
    real, pointer, contiguous:: grav_top(:,:,:)
    real, pointer, contiguous:: p_bb(:,:), d_bb(:,:)
    real, pointer, contiguous:: p_ub(:,:), d_ub(:,:)
    integer:: l_z, u_z
    logical:: surface=.false.
    logical:: initialized=.false.
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
    procedure:: upper
    procedure:: upper_surface
    procedure:: lower
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t)     :: self
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  integer:: i0, ix, iy, iz, n(3), zi
  real(var_k), dimension(:,:,:), pointer:: d, e, pg
  real(8):: pos(3), rp, xp, yp, zp
  integer:: iostat
  logical, save:: set_bc_bits=.false.
  logical, save:: first_time=.true.
  namelist /boundary_params/ set_bc_bits
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  patch => self%patch
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, boundary_params)
    end if
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  self%set_bc_bits=set_bc_bits
  call self%standard_boundaries_t%init ()
  !
  !
  m     => patch%mesh
  n = m%gn
  n(3)=m(3)%ng+1
  self%l_z=m(3)%li
  self%u_z=m(3)%ui
  
  if (self%is_set(bits%zu)) then
    allocate(self%grav_top(n(1),n(2),n(3)))
    pos = patch%position
    do iz=1,n(3)
      zi=m(3)%ui+iz-1
      do iy=m(2)%lb, m(2)%ub
        yp = m(2)%r(iy)
        do ix=m(1)%lb, m(1)%ub
          xp = m(1)%r(ix)
          rp=sqrt(pos(3)**2)
          zp = m(3)%r(zi) + rp
          rp = sqrt(zp**2)
          ! ----------------------------------------------------------------------
          ! interpolate table values and store them
          ! ----------------------------------------------------------------------
          self%grav_top(ix,iy,iz) = interpolate%one_d (rp, model_1d%r, model_1d%g)
        end do
      end do
    end do
    allocate (d(n(1),n(2),1))
    allocate (e(n(1),n(2),1))
    allocate (pg(n(1),n(2),1))
    iz=self%u_z
    do iy=m(2)%lb, m(2)%ub
      yp = m(2)%r(iy)
      do ix=m(1)%lb, m(1)%ub
        xp = m(1)%r(ix)
        rp=sqrt(pos(3)**2)
        zp = m(3)%r(iz) + rp
        rp = sqrt(zp**2)
        d(ix,iy,1) = interpolate%one_d (rp, model_1d%r, model_1d%d)
        e(ix,iy,1) = interpolate%one_d (rp, model_1d%r, model_1d%e)
      end do
    end do
    call eos%lookup(d=d,e=e,pg=pg)
    allocate(self%d_ub(n(1),n(2)))
    allocate(self%p_ub(n(1),n(2)))
    self%d_ub=d(:,:,1)
    self%p_ub=pg(:,:,1)
    deallocate(d,e,pg)
  end if
  !
  if (self%is_set(bits%zl)) then
    allocate(self%grav_bot(n(1),n(2),n(3)))
    pos = patch%position
    do iz=m(3)%lb, m(3)%li
      do iy=m(2)%lb, m(2)%ub
        yp = m(2)%r(iy)
        do ix=m(1)%lb, m(1)%ub
          xp = m(1)%r(ix)
          rp=sqrt(pos(3)**2)
          zp = m(3)%r(iz) + rp
          rp = sqrt(zp**2)
          ! ----------------------------------------------------------------------
          ! interpolate table values and store them
          ! ----------------------------------------------------------------------
          self%grav_bot(ix,iy,iz) = interpolate%one_d (rp, model_1d%r, model_1d%g)
        end do
      end do
    end do
    allocate (d(n(1),n(2),1))
    allocate (e(n(1),n(2),1))
    allocate (pg(n(1),n(2),1))
    iz=self%l_z
    do iy=m(2)%lb, m(2)%ub
      yp = m(2)%r(iy)
      do ix=m(1)%lb, m(1)%ub
        xp = m(1)%r(ix)
        rp=sqrt(pos(3)**2)
        zp = m(3)%r(iz) + rp
        rp = sqrt(zp**2)
        d(ix,iy,1) = interpolate%one_d (rp, model_1d%r, model_1d%d)
        e(ix,iy,1) = interpolate%one_d (rp, model_1d%r, model_1d%e)
      end do
    end do
    call eos%lookup(d=d,e=e,pg=pg)
    allocate(self%d_bb(n(1),n(2)))
    allocate(self%p_bb(n(1),n(2)))
    self%d_bb=d(:,:,1)
    self%p_bb=pg(:,:,1)

    deallocate(d,e,pg)
  end if
  self%initialized=.true.
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%is_set(bits%zl)) then
    call self%lower(mem)
  end if
  if (self%is_set(bits%zu)) then
   ! if (self%surface) then
   !   print *, 'surf'
      call self%upper_surface(mem)
   ! else
   !   call self%upper(mem)
   ! end if
  end if
END SUBROUTINE primary

!===============================================================================
!> Lower BC: Constant gas pressure p_bb=1, constant entropy s_bb=0 for
!> ascending gas, while descending gas keeps its entropy
!===============================================================================
SUBROUTINE lower(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  real(kind=var_k), dimension(:,:,:), pointer:: d,s,p1, t1, dd
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz, e, pp,ee
  real(kind=var_k), dimension(:,:,:), pointer:: px, py, pz
  real(kind=var_k), dimension(:,:,:), pointer:: bx, by, bz
  integer:: ix, iy, iz, iz1, n(3), id, lb, ub, l, u
  real(8):: ds(3), z_bb
  real:: z1, h,p, d1
  real(8), pointer:: z(:)
  ! ............................................................................
  patch => self%patch
  m     => self%mesh
  allocate (z(patch%gn(3)))
  z = patch%mesh(3)%r + patch%mesh(3)%p
  z_bb=z(self%l_z)
  !-----------------------------------------------------------------------------
  ! Allocate scratches
  !-----------------------------------------------------------------------------
  n = m%gn
  n(3)=m(3)%li
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (e (n(1),n(2),n(3)))
  allocate (ee(n(1),n(2),n(3)-1))
  allocate (pp(n(1),n(2),n(3)-1))
  !-----------------------------------------------------------------------------
  ! Point pointers
  !-----------------------------------------------------------------------------
  d  => mem(:,:,m(3)%lb:m(3)%li,patch%idx%d)
  px => mem(:,:,m(3)%lb:m(3)%li,patch%idx%px)
  py => mem(:,:,m(3)%lb:m(3)%li,patch%idx%py)
  pz => mem(:,:,m(3)%lb:m(3)%li,patch%idx%pz)
  if (patch%idx%bx>0) then
    bx => mem(:,:,m(3)%lb:m(3)%li,patch%idx%px)
    by => mem(:,:,m(3)%lb:m(3)%li,patch%idx%py)
    bz => mem(:,:,m(3)%lb:m(3)%li,patch%idx%pz)
  end if

  !-----------------------------------------------------------------------------
  ! Get velocities
  !-----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%li
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = px(ix,iy,iz)/d(ix,iy,iz)
        uy(ix,iy,iz) = py(ix,iy,iz)/d(ix,iy,iz)
        uz(ix,iy,iz) = pz(ix,iy,iz)/d(ix,iy,iz)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Get thermal energy, if applicable
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='total') then
    do iz=m(3)%lb,m(3)%li
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,iz,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                         ux(ix,iy,iz)**2 &
                       + uy(ix,iy,iz)**2 &
                       + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='thermal') then
    do iz=m(3)%lb,m(3)%li
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,iz,patch%idx%e)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    s => mem(:,:,:,patch%idx%s)
    do iz=m(3)%lb,m(3)%li
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = d(ix,iy,iz)**eos%gamma* &
          exp(s(ix,iy,iz)*(eos%gamma-1d0)/d(ix,iy,iz))/(eos%gamma-1d0)
        end do
      end do
    end do
  else
    call io%abort("boundaries_t%lower:: unknown etype")
  end if

  allocate(p1(n(1),n(2),1))
  allocate(t1(n(1),n(2),1))
  
  do iz=n(3)-1,1,-1
    call eos%lookup(d=d(:,:,iz+1:iz+1),e=e(:,:,iz+1:iz+1),pg=p1, tt=t1)
    do iy=1,n(2)
      do ix=1,n(1)
    !---------------------------------------------------------------------------
    ! Compute the local pressure scale height and extrapolate pressure with that
    !---------------------------------------------------------------------------
    d1 = d(ix,iy,iz+1)
    h = t1(ix,iy,1)/self%grav_bot(ix,iy,iz)
    !---------------------------------------------------------------------------
    ! Set the entropy per unit mass of ascending gas to the initial s_bb=0,
    ! while leaving entropy of descending gas untouched.
    ! We know, that constant \gamma is OK at the bottom.
    !---------------------------------------------------------------------------
    if (uz(ix,iy,n(3)) >= 0.0) then
      p = self%p_bb(ix,iy)*exp((z_bb-z(iz))/h)
      d(ix,iy,iz) = self%d_bb(ix,iy)*(p/self%p_bb(ix,iy))**(1d0/eos%gamma)
    else
      p = p1(ix,iy,1)*exp((z(iz+1)-z(iz))/h)
      d(ix,iy,iz) = d1*(p/p1(ix,iy,1))**(1d0/eos%gamma)
    end if
    pp(ix,iy,iz) = p 
    !---------------------------------------------------------------------------
    ! Vanishing vertical derivative of velocity components
    !---------------------------------------------------------------------------
    ux(ix,iy,iz) = ux(ix,iy,n(3))
    uy(ix,iy,iz) = uy(ix,iy,n(3))
    uz(ix,iy,iz) = uz(ix,iy,n(3))
  end do
  end do
  end do
  dd => d(:,:,:n(3)-1)
  call eos%inverse_lookup(d=dd,v=pp,name='P',ee=ee)
  do iz=1,n(3)-1
  do iy=1,n(2)
  do ix=1,n(1)
    e(ix,iy,iz)=ee(ix,iy,iz)*d(ix,iy,iz)
  end do
  end do
  end do


  
  ! ----------------------------------------------------------------------------
  ! Put momenta back to mem
  ! ----------------------------------------------------------------------------
  do iz=m(3)%lb,m(3)%lo
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        px (ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
        py (ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
        pz (ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
      end do
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Put the total energy or entropy back into mem, if applicable
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='total') then
    do iz=m(3)%lb,m(3)%lo
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,iz,patch%idx%e) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                      ux(ix,iy,iz)**2 &
                                    + uy(ix,iy,iz)**2 &
                                    + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='thermal') then
    do iz=m(3)%lb,m(3)%lo
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,iz,patch%idx%e) =  e(ix,iy,iz)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    do iz=m(3)%lb,m(3)%lo
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,iz,patch%idx%s) = log((eos%gamma-1d0)* &
          e(ix,iy,iz)/d(ix,iy,iz)**eos%gamma)*d(ix,iy,iz)/(eos%gamma-1d0)
        end do
      end do
    end do
  end if

  if (associated(patch%filled)) then
    do iz=m(3)%lb,m(3)%lo
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          patch%filled(ix,iy,iz) = patch%level
        end do
      end do
    end do
  end if

  deallocate(z, ux, uy, uz, e, p1, t1, pp, ee)

END SUBROUTINE lower

!===============================================================================
!> Upper BC at high depth: 
!> Constant gas pressure p_ub=1, constant entropy s_bb=0 for
!> ascending gas, while descending gas keeps its entropy
!===============================================================================
SUBROUTINE upper(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  real(kind=var_k), dimension(:,:,:), pointer:: d, s,p1, t1,dd
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz,e, ee, pp
  real(kind=var_k), dimension(:,:,:), pointer:: px, py, pz
  real(kind=var_k), dimension(:,:,:), pointer:: bx, by, bz
  integer:: ix, iy, iz, iz1, n(3), id, lb, ub, l, u, zi
  real(8):: ds(3), z_ub
  real:: z1, h,p, d1, tt
  real(8), pointer:: z(:)
  ! ............................................................................
  id    =  self%patch%id
  ds    =  self%patch%ds
  patch => self%patch
  m     => self%mesh
  allocate (z(patch%gn(3)))
  z = patch%mesh(3)%r !+ patch%mesh(3)%p
  z_ub=z(self%u_z)
  !-----------------------------------------------------------------------------
  ! Allocate scratches
  !-----------------------------------------------------------------------------
  n = m%gn
  n(3)=m(3)%ng+1
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (e (n(1),n(2),n(3)))
  allocate (ee(n(1),n(2),n(3)-1))
  allocate (pp(n(1),n(2),n(3)-1))

  !-----------------------------------------------------------------------------
  ! Point pointers
  !-----------------------------------------------------------------------------
  d  => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%d)
  px => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%px)
  py => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%py)
  pz => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%pz)
  if (patch%idx%bx>0) then
    bx => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%bx)
    by => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%by)
    bz => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%bz)
  end if

  !-----------------------------------------------------------------------------
  ! Get velocities
  !-----------------------------------------------------------------------------
  do iz=1,n(3)
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = px(ix,iy,iz)/d(ix,iy,iz)
        uy(ix,iy,iz) = py(ix,iy,iz)/d(ix,iy,iz)
        uz(ix,iy,iz) = pz(ix,iy,iz)/d(ix,iy,iz)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Get thermal energy, if applicable
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='total') then
    do iz=1,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,zi,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                         ux(ix,iy,iz)**2 &
                       + uy(ix,iy,iz)**2 &
                       + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='thermal') then
    do iz=1,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,zi,patch%idx%e)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    s => mem(:,:,:,patch%idx%s)
    do iz=1, n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = d(ix,iy,iz)**eos%gamma* &
          exp(s(ix,iy,zi)*(eos%gamma-1d0)/d(ix,iy,iz))/(eos%gamma-1d0)
        end do
      end do
    end do
  else
    call io%abort("boundaries_t%upper:: unknown etype")
  end if

  allocate(p1(n(1),n(2),1))
  allocate(t1(n(1),n(2),1))
  !print *, 'upper bnd minmax e, d', minval(e(:,:,1:1)),maxval(e(:,:,1:1)), minval(d(:,:,1:1)),maxval(d(:,:,1:1))
  call eos%lookup(d=d(:,:,1:1),e=e(:,:,1:1),pg=p1, tt=t1)
ee=e/d
  do iy=1,n(2)
    do ix=1,n(1)
      d1 = d(ix,iy,1)
      do iz=2,n(3)
        zi = m(3)%ui+iz-1
        !-----------------------------------------------------------------------
        ! Compute the local pressure scale height and extrapolate pressure with that
        !-----------------------------------------------------------------------
        !tt=p1(ix,iy,1)/d1
        !h = t1(ix,iy,1)/grav
        h = tt/self%grav_top(ix,iy,1)
        !-----------------------------------------------------------------------
        ! Set the entropy per unit mass of descending gas to the initial s_bb=0,
        ! while leaving entropy of ascending gas untouched
        !-----------------------------------------------------------------------
        
        if (uz(ix,iy,1) < 0.0) then
          !pp(ix,iy,iz-1) = self%p_ub(ix,iy)*exp((z_ub-z(zi))/h)
          !d(ix,iy,iz) = self%d_ub(ix,iy)*(pp(ix,iy,iz-1)/self%p_ub(ix,iy))**(1d0/eos%gamma)
          d(ix,iy,iz)=d(ix,iy,1)*exp((z_ub-z(zi))/h)
        else
          d(ix,iy,iz)=d(ix,iy,1)*exp((z_ub-z(zi))/h)
          !pp(ix,iy,iz-1) = p1(ix,iy,1)*exp((z(m(3)%ui)-z(zi))/h)
          !d(ix,iy,iz) = d1*(pp(ix,iy,iz-1)/p1(ix,iy,1))**(1d0/eos%gamma)
        end if
        ee(ix,iy,iz-1)=ee(ix,iy,iz-1)*d(ix,iy,iz)
       ! if (ix==12 .and.iy==12) then
       !   print *, iz, p, p1(ix,iy,1), self%p_ub(ix,iy), h,exp((z_ub-z(zi))/h),exp((z(m(3)%ui)-z(zi))/h) , z_ub, z(zi), z(m(3)%ui)
       ! end if
        !e(ix,iy,iz) = p/(eos%gamma-1d0)
        !e(ix,iy,iz) = e(ix,iy,1)*exp((z_ub-z(zi))/h)
        !-----------------------------------------------------------------------
        ! Vanishing vertical derivative of velocity components
        !-----------------------------------------------------------------------
        ux(ix,iy,iz) = ux(ix,iy,1)
        uy(ix,iy,iz) = uy(ix,iy,1)
        uz(ix,iy,iz) = uz(ix,iy,1)
      end do
    end do
  end do
 ! dd=> d(:,:,2:)
 ! call eos%inverse_lookup(d=dd,v=pp,name='P',ee=ee)
  do iy=1,n(2)
    do ix=1,n(1)
      d1 = d(ix,iy,1)
      do iz=2,n(3)
      e(ix,iy,iz) = ee(ix,iy,iz-1)!*d(ix,iy,iz)
      end do
    end do
  end do


  
  ! ----------------------------------------------------------------------------
  ! Put momenta back to mem
  ! ----------------------------------------------------------------------------
  do iz=2,n(3)
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        px (ix,iy,iz) = d(ix,iy,iz)*ux(ix,iy,iz)
        py (ix,iy,iz) = d(ix,iy,iz)*uy(ix,iy,iz)
        pz (ix,iy,iz) = d(ix,iy,iz)*uz(ix,iy,iz)
      end do
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Put the energy or entropy back into mem,
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='thermal') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%e) =  e(ix,iy,iz)
        end do
      end do
    end do
  else if (trim(patch%etype)=='total') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%e) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                      ux(ix,iy,iz)**2 &
                                    + uy(ix,iy,iz)**2 &
                                    + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%s) = log((eos%gamma-1d0)* &
          e(ix,iy,iz)/d(ix,iy,iz)**eos%gamma)*d(ix,iy,iz)/(eos%gamma-1d0)
        end do
      end do
    end do
  end if

  if (associated(patch%filled)) then
    do iz=m(3)%uo,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          patch%filled(ix,iy,iz) = patch%level
        end do
      end do
    end do
  end if

  deallocate(z, ux, uy, uz, e, ee, pp, p1, t1)

END SUBROUTINE upper


!===============================================================================
!> Upper BC at Solar surface: 
!>vanishing vertical derivatives of velocity components, constant
!> local temperature and scale heights in the guard zones
!===============================================================================
SUBROUTINE upper_surface(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: m(:)
  real(kind=var_k), dimension(:,:,:), pointer:: d,s, t1, ee
  real(kind=var_k), dimension(:,:,:), pointer, contiguous:: ux, uy, uz,e
  real(kind=var_k), dimension(:,:,:), pointer:: px, py, pz
  real(kind=var_k), dimension(:,:,:), pointer:: bx, by, bz
  integer:: ix, iy, iz, iz1, n(3), id, lb, ub, l, u, zi
  real(8):: ds(3)
  real::  z1, h,p, d_ratio
  real(8), pointer:: z(:)
  ! ............................................................................
  patch => self%patch
  m     => self%mesh
  allocate (z(patch%gn(3)))
  z = patch%mesh(3)%r + patch%mesh(3)%p
  !-----------------------------------------------------------------------------
  ! Allocate scratches
  !-----------------------------------------------------------------------------
  n = m%gn
  n(3)=m(3)%ng+1
  allocate (ux(n(1),n(2),n(3)))
  allocate (uy(n(1),n(2),n(3)))
  allocate (uz(n(1),n(2),n(3)))
  allocate (e (n(1),n(2),n(3)))
  allocate (ee(n(1),n(2),n(3)))
  !-----------------------------------------------------------------------------
  ! Point pointers
  !-----------------------------------------------------------------------------
  d  => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%d)
  px => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%px)
  py => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%py)
  pz => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%pz)
  if (patch%idx%bx>0) then
    bx => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%px)
    by => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%py)
    bz => mem(:,:,m(3)%ui:m(3)%ub,patch%idx%pz)
  end if

  !-----------------------------------------------------------------------------
  ! Get velocities
  !-----------------------------------------------------------------------------
  do iz=1,n(3)
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ux(ix,iy,iz) = px(ix,iy,iz)/d(ix,iy,iz)
        uy(ix,iy,iz) = py(ix,iy,iz)/d(ix,iy,iz)
        uz(ix,iy,iz) = pz(ix,iy,iz)/d(ix,iy,iz)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Get thermal energy, if applicable
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='total') then
    do iz=1,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,zi,patch%idx%e) - 0.5*d(ix,iy,iz)*( &
                         ux(ix,iy,iz)**2 &
                       + uy(ix,iy,iz)**2 &
                       + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='thermal') then
    do iz=1,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = mem(ix,iy,zi,patch%idx%e)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    s => mem(:,:,:,patch%idx%s)
    do iz=1, n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          e(ix,iy,iz) = d(ix,iy,iz)**eos%gamma* &
          exp(s(ix,iy,zi)*(eos%gamma-1d0)/d(ix,iy,iz))/(eos%gamma-1d0)
        end do
      end do
    end do
  else
    call io%abort("boundaries_t%upper_surface:: unknown etype")
  end if

  ! ----------------------------------------------------------------------------
  ! Get specific energy
  ! -----------------------------------------------------------------------------
  ee = e/d
  ! ----------------------------------------------------------------------------
  ! Lookup gas temperature
  ! ----------------------------------------------------------------------------
  allocate(t1(n(1),n(2),1))
  
  z1 = z(self%mesh(3)%ui)
  do iz=2,n(3)
    call eos%lookup(d=d(:,:,iz-1:iz-1),e=e(:,:,iz-1:iz-1),tt=t1)
    do iy=1,n(2)
      do ix=1,n(1)
        h = t1(ix,iy,1)/self%grav_top(ix,iy,iz-1)
        d_ratio = exp((z(iz-1)-z(iz))/h)
        d(ix,iy,iz) = d(ix,iy,1)*d_ratio
        !-----------------------------------------------------------------------
        ! The horizontal components of the momenta are assumed to taper off
        ! with density, consistent with assuming a no-shear boundary condition
        ! for the horizontal velocity components
        !-----------------------------------------------------------------------
        px(ix,iy,iz) = px(ix,iy,1)*d_ratio
        py(ix,iy,iz) = py(ix,iy,1)*d_ratio
        !-----------------------------------------------------------------------
        ! The vertical velocity is tapered of with the d_ratio factor.
        !-----------------------------------------------------------------------
        pz(ix,iy,iz) = pz(ix,iy,1)*d_ratio**2
      end do
    end do
  end do
  ! ----------------------------------------------------------------------------
  ! Convert back to thermal energy, assuming constant specific energy in the
  ! boundaries.
  ! ----------------------------------------------------------------------------
  e = ee*d
  !-----------------------------------------------------------------------------
  ! Put the energy or entropy back into mem,
  !-----------------------------------------------------------------------------
  if (trim(patch%etype)=='thermal') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%e) =  e(ix,iy,iz)
        end do
      end do
    end do
  else if (trim(patch%etype)=='total') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%e) =  e(ix,iy,iz) + 0.5*d(ix,iy,iz)*( &
                                      ux(ix,iy,iz)**2 &
                                    + uy(ix,iy,iz)**2 &
                                    + uz(ix,iy,iz)**2)
        end do
      end do
    end do
  else if (trim(patch%etype)=='entropy') then
    do iz=2,n(3)
      zi = m(3)%ui+iz-1
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          mem(ix,iy,zi,patch%idx%s) = log((eos%gamma-1d0)* &
          e(ix,iy,iz)/d(ix,iy,iz)**eos%gamma)*d(ix,iy,iz)/(eos%gamma-1d0)
        end do
      end do
    end do
  end if

  deallocate(z, ux, uy, uz, e, t1, ee)

  if (associated(patch%filled)) then
    do iz=m(3)%uo,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          patch%filled(ix,iy,iz) = patch%level
        end do
      end do
    end do
  end if

END SUBROUTINE upper_surface

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)

END SUBROUTINE derivative


END MODULE boundaries_mod
