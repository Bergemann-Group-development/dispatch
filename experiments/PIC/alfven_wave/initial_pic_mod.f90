!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_pic_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE params,             ONLY: master, hl, stdout, out_namelists, params_unit, stdin
  USE grid_m,             ONLY: Bx,By,Bz,Ex,Ey,Ez,g
  USE units,              ONLY: c, elm
  USE species_m,          ONLY: sp, fields, nspecies, field
  USE debug,              ONLY: debugit, dbg_init
  USE particleoperations, ONLY: make_particle_pair, make_particle
  USE initial_mod
  implicit none
  private
  
  type, public, extends(initial_t):: initial_pic_t
     real    :: d0_pic                                                         ! background density
     real    :: v0                                                             ! veloccity amplitude 
     integer :: ksign
     logical :: do_pairs, periodic
     real    :: a0_pic, b0_pic
  contains
    procedure :: init_pic => init_pic_experiment
  end type

  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!===============================================================================

SUBROUTINE init_pic_experiment (self)
 implicit none
  class(initial_pic_t) :: self
  logical, save :: mhd=.true.
  real, save    :: d0_pic                                                         ! background density
  real, save    :: v0                                                         ! veloccity amplitude 
  logical, save :: do_pairs
  integer, save :: ksign
  logical, save :: first_time=.true.
  logical, save :: periodic
  real,    save :: x0=0.0, y0=0.0, z0=0.0, r0=0.1, a0_pic=0.1, b0_pic, b0, a0, temp, gamma, beta, expe
  integer, save :: isp, i, ix, iy, iz, isign
  real, save    :: x, y, z, sign, vz, v, w, wpe, m_electron, q_electron, vg(3)
  namelist /init/ d0_pic, v0, do_pairs, ksign, a0_pic
  namelist /IC_params/ x0, y0, z0, r0, b0, a0, temp, gamma, beta, expe
  
  
  !.............................................................................
  call trace_begin ('initial_t%init_experiment')

  ! default values for namelist parameters
   periodic = .true.                                                     ! MUST define for periodic
  
   d0_pic = 0.1                                                              ! background density
   v0 = 0.0                                                              ! streaming velocity w.r.t. rest frame is zero
   do_pairs = .true.                                                     ! initialize in pairs to remove initial electric field
  
  !$omp critical (initial_t)
  if (first_time) then
     rewind(stdin); read (stdin,init)                                       ! rewind file descriptor and read input from file
     rewind(stdin); read(stdin, initial_params)
     if (master) write(params_unit,init)
     if (out_namelists) write(stdout,init)
     first_time = .false.
  end if
  !$omp end critical (initial_t)

  v0 = 0.0                                                               ! streaming velocity w.r.t. rest frame is STILL zero
  
  
  self%d0_pic = d0_pic
  self%v0 = v0
  self%do_pairs = do_pairs
  self%periodic = periodic
  self%x0 = x0
  self%y0 = y0
  self%z0 = z0
  self%r0 = r0
  self%a0_pic = a0_pic
  self%b0_pic = b0_pic
  self%temp = temp

    Bx=0.                                                                 ! No fields initially
    By=0.
    Bz=0.
    Ex=0.
    Ey=0.
    Ez=0.
  
    do isp=1,nspecies                                                     ! loop over species
      if (trim(sp(isp)%name) == 'electron') then
        m_electron = sp(isp)%mass                                         ! mass
        q_electron = abs(sp(isp)%charge)                                  ! size of charge
      end if
    end do

  wpe = sqrt(self%d0_pic*c%fourpi*elm%ke/m_electron)*q_electron                  ! electron plasma freq
   


   
    if (master) then
    write(stdout,*) hl
    write(stdout,*) 'Electron plasma frequency =',wpe
    write(stdout,*) 'Electron skin depth [dx]  =',(c%c / wpe) / g%ds(1)
    write(stdout,*) 'Ion      skin depth [dx]  =',(c%c / wpe) / g%ds(1)
  endif
  
  
  do isp=1,nspecies                                                   ! loop over species
    do iz=g%lb(3),g%ub(3)-1                                           ! loop over grid cells
      !w = self%d0_pic/(sp(isp)%ntarget)                                        ! density weight
      !z = (g%z(iz)-g%grlb(3)) - self%z0                                    ! 0 < z < 1
      do iy=g%lb(2),g%ub(2)-1                                         ! y-loop
        !y = (g%y(iy)-g%grlb(2)) - self%y0                                  ! 0 < y < 1
        do ix=g%lb(1),g%ub(1)-1                                       ! x-loop
          !block
          !real:: r2, r, f, ex, ey, d, j_phi, jx, jy, vx, vy, mtoti, d_pg
          !x = (g%x(ix)-g%grlb(1)) - self%x0                                ! 0 < x < 1
          !r2 = (x**2+y**2)/self%r0**2
          !r = max(sqrt(r2),1e-30)
          !f = exp(-r2)
          !ex = x/r
          !ey = y/r
          !d = self%d0_pic*(1.-self%a0_pic*f**2)
          !temp = sp(1)%temperature
          !d_pg = self%d0_pic*self%a0_pic*self%temp
          !self%b0_pic = sqrt(2.*d_pg)
          !bx(ix,iy,iz) = 0.0
          !by(ix,iy,iz) = 0.0
          !bz(ix,iy,iz) = self%b0_pic*f
          !j_phi = -2.*r*bz(ix,iy,iz)/self%r0**2
          !jy =  j_phi*ex
          !jx = -j_phi*ey
          !mtoti = sum(1./sp%mass)
          !vx = jx/(sp(isp)%charge*d)*mtoti/sp(isp)%mass                 ! species bulk speed
          !vy = jy/(sp(isp)%charge*d)*mtoti/sp(isp)%mass                 ! species bulk speed
          w = fields(ix,iy,iz,isp)%d / sp(isp)%ntarget                   ! numdens/npercell
          vg = fields(ix,iy,iz,isp)%v
          !vg = (/vx,vy,0./)
          !fields(ix,iy,iz,isp)%d = d                                    ! number density
          if (self%do_pairs) then 
            if (mod(float(isp), float(2)) /= 0) then
              do i=1,sp(isp)%ntarget                                      ! particles/cell
                call make_particle_pair(isp,isp+1,ix,iy,iz,vg,w,.false.)  ! add particle pair (nb vg==0.)
                !ff(ix,iy,iz,isp)%d = self%d0                               ! number density, this cell
              end do
            end if
          else
            do i=1,sp(isp)%ntarget                                      ! particles/cell
              call make_particle(isp, ix, iy, iz, vg, w)                ! add particle
              !ff(ix,iy,iz,isp)%d = self%d0                               ! number density, this cell
            end do
          end if
          !end block
        end do
      end do
    end do
  end do

  call trace_end
END SUBROUTINE init_pic_experiment

END MODULE initial_pic_mod
