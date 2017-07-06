! $Id$
! vim: nowrap
!-----------------------------------------------------------------------!
! Initialize a plasmaoscillation, with no E- or Bfield.
! We assume c=1.  The default is a harmonic perturbation wave in
! the z-direction.
!
! SCALINGS: When changing the total volume and nothing else, wpe should
! remain the same. The total energies should scale with the volume. The
! particle weights remain the same, and therefore the charge densities 
! and current densities should also remain the same.  That, again, should
! mean that the field values remain the same, and that the total EM
! energy therefore indeed should scale with the volume.
!
! The coherent initial current density is the same when the wavelength 
! is shorter, but the divergence of it is larger, so the electric charge
! grows faster.  The coherent charge densities scale inversely with the
! size, while the coherent electric potential scales linearly with the
! size.  The noise charge density, on the other hand, starts out the 
! same, so the initial noise electric field scales with the size, and is
! smaller for a smaller size.

!=======================================================================
! Leave boundaries as they are, so no need to even read a namelist
!=======================================================================
SUBROUTINE init_experiment_boundaries
implicit none
END SUBROUTINE init_experiment_boundaries

!=======================================================================
!SUBROUTINE input (inputfile)
!  USE params, ONLY: mfile
!  USE dumps,  ONLY: dump, dump_set
!  implicit none
!  character(len=mfile) inputfile
!.......................................................................
!  inputfile = 'input_pic.nml'                                              ! Where to find input
!END SUBROUTINE

!=======================================================================
!SUBROUTINE init_pic_experiment
!  USE params
!  USE grid_m,             ONLY: Bx, By, Bz, Ex, Ey, Ez, g
!  USE units,              ONLY: c, elm
!  USE species_m,          ONLY: sp, fields
!  USE dumps,              ONLY: dump, dump_set
!  USE particleoperations, ONLY: make_particle
!  implicit none
!  real, save:: d0                                                       ! background density
!  real, save:: v0                                                       ! veloccity amplitude
!  logical, save :: do_pairs
!  integer, save:: kz                                                    ! wave number
!  integer isp, i, ix, iy, iz
!  real :: x, y, z, sign, vz, vg(3), w, wpe, m_electron, q_electron
!  real, parameter:: pi=3.14159265
!  namelist /init/ d0, v0, do_pairs, kz                                            ! alternative perturb
!  logical, save:: first_call = .true.
!.......................................................................
!  call dump_set ('experiment')
!  periodic = .true.                                                     ! MUST define for periodic

!  Bx=0; By=0; Bz=0; Ex=0; Ey=0; Ez=0                                    ! only plasmaoscillations

  !---------------------------------------------------------------------
  ! In the DISPATCH context, this routine is called once for every patch,
  ! in order to initialize its population of particle species.  To avoid
  ! that the namelist is read and written to the log file multiple times
  ! we use the first_call switch.  To ensure that only one thread gets
  ! the privilege of reading the namelist, and that it does so only once,
  ! we enclose the first_call logic inside a critical region. To ensure
  ! that the values read are available to the next thread they must have
  ! the "save" attribute.
  !---------------------------------------------------------------------
 ! !$omp critical (init_cr)
!  if (first_call) then
!    first_call = .false.
!    kz = 1                                                              ! wavenumber in z-dir
!    d0 = 1.                                                             ! background density
!    v0 = 0.1                                                            ! perturb strength
!    rewind(stdin); read (stdin,init)
!    if (master) write(params_unit,init)
!    if (out_namelists) write(stdout,init)
!  end if
!  !$omp end critical (init_cr)

!  do isp=1,nspecies                                                     ! loop over species
!    if (trim(sp(isp)%name) == 'electron') then
!      m_electron = sp(isp)%mass
!      q_electron = sp(isp)%charge
!    end if
!  end do
!  wpe = sqrt(d0*c%fourpi*elm%ke/m_electron)*abs(q_electron)             ! electron plasma freq
!  if (master) write(stdout,*) 'Electron plasma frequency =',wpe

! do isp=1,nspecies                                                     ! loop over species
!    sign = 2*modulo(isp,2)-1.                                           ! alternating signs
!    do iz=g%lb(3),g%ub(3)-1                                             ! loop over grid cells
!      z = (g%z(iz)-g%grlb(3))/g%gs(3)                                   ! 0 < z < 1
!      vz = v0*sign*sin(2.*pi*kz*z)*m_electron/sp(isp)%mass              ! amplitude
!      vg = (/0.,0.,vz/)
!      w = d0/(sp(isp)%ntarget)                                          ! density weight
!      do iy=g%lb(2),g%ub(2)-1                                           ! y-loop
!        do ix=g%lb(1),g%ub(1)-1                                         ! x-loop
!          do i=1,sp(isp)%ntarget                                        ! particles/cell
!            fields(ix,iy,iz,isp)%d = d0                                 ! number density
!            call make_particle(isp, ix, iy, iz, vg, w)                  ! add particle
!          end do
!        end do
!      end do
!    end do
!                                 call dump (fields(:,:,:,isp)%d   , 'd' )
!  end do
!  print*,'Expected kinetic energy:',0.25*v0**2*sp(1)%mass*d0*g%gV
!
!END SUBROUTINE init_pic_experiment
