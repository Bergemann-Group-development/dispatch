MODULE rt_params_mod
  USE io_mod
  USE ray_mod, only: set_carlson_quadratures, set_radau_quadratures, &
                     set_lobatto_quadratures, set_gauss_quadratures, &
                     set_custom_quadratures, set_grid_quadratures
  USE handyfuncs_mod, only: count_words, lowercase
  USE console_mod
  USE composition_mod, only: comp_sample_linelist => sample_linelist
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Configuration of one radiative-transfer pass.  Two instances exist:
  !   rt_iter -- the NLTE iteration pass (short characteristics in 3D, column
  !              long-characteristics solves in 1D), read from &rt_iteration_params
  !   rt_out  -- the final output-spectrum pass (long characteristics), read
  !              from &rt_output_params
  ! Both groups share the same member names (solver, angle_quad, ntheta, nphi,
  ! maxiter, sample_back_opac, sample_linelist, include_linelist, include_absmet);
  ! &rt_output_params additionally has custom_mu, ray_offset and pivot_z, which
  ! only apply to the inclined output rays, and save_linelist, which controls the
  ! linelist_out.txt dump written at the end of the output pass.
  ! For the output pass, maxiter is the number of fixed-population ("Lambda")
  ! iterations appended after the rt_iter%maxiter NLTE iterations (was ilambd).
  ! Defaults differ per pass: include_linelist is F/T and include_absmet T/F
  ! for rt_iter/rt_out (see rt_params_nml).
  !-----------------------------------------------------------------------------
  type, public :: rt_params_t
    character(len=128) :: solver = 'linear'
    character(len=32)  :: angle_quad = ''
    integer :: ntheta=0, nphi=0
    integer :: maxiter=0
    logical :: sample_back_opac=.false.
    logical :: sample_linelist=.false.  ! line list via sampled LTE table instead of detailed lines
    logical :: include_linelist=.false. ! detailed line-list opacity in this pass
    logical :: include_absmet=.false.   ! sampled metal background (absmet) opacity in this pass
    real(8), allocatable :: vec(:,:)  ! ray directions (3,nmu), set by setup_quadrature
    real(8), allocatable :: wmu(:)    ! integration weights,   set by setup_quadrature
  end type rt_params_t

  type(rt_params_t), public, save, target :: rt_iter, rt_out

  ! Output-ray geometry (members of &rt_output_params only)
  character(len=512), public, save :: custom_mu=''
  real(8), public :: ray_offset(2)=[0.,0.], pivot_z=-1
  ! linelist_out.txt lists every transition the run kept, with its strength at the
  ! run's own composition -- a diagnostic, never an input. A curve of growth repeats
  ! one atmosphere at many abundances and only needs the dump once, so it can be
  ! switched off for the rungs that would otherwise rewrite it.
  logical, public, save :: save_linelist=.true.

  public :: read_rt_params, echo_rt_params, setup_quadrature

CONTAINS

  !===============================================================================
  SUBROUTINE read_rt_params ()
    call rt_params_nml (do_echo=.false.)
  END SUBROUTINE read_rt_params

  !===============================================================================
  SUBROUTINE echo_rt_params ()
    call rt_params_nml (do_echo=.true.)
  END SUBROUTINE echo_rt_params

  !===============================================================================
  !> Read (or echo) both parameter groups.  Namelist member names are the local
  !> variable names, shared by the two groups; values are copied to/from the
  !> rt_iter and rt_out instances around the transfer.
  !===============================================================================
  SUBROUTINE rt_params_nml (do_echo)
    logical, intent(in) :: do_echo
    character(len=128) :: solver
    character(len=32)  :: angle_quad
    integer :: ntheta, nphi, maxiter
    logical :: sample_back_opac, sample_linelist, include_linelist, include_absmet
    integer :: iostat
    namelist /rt_iteration_params/ solver, angle_quad, ntheta, nphi, maxiter, sample_back_opac, &
      sample_linelist, include_linelist, include_absmet
    namelist /rt_output_params/ solver, angle_quad, ntheta, nphi, maxiter, sample_back_opac, &
      sample_linelist, include_linelist, include_absmet, &
      custom_mu, ray_offset, pivot_z, save_linelist
    !-----------------------------------------------------------------------------
    if (do_echo) then
      call to_locals (rt_iter)
      write (io%namelist, rt_iteration_params)
      call to_locals (rt_out)
      write (io%namelist, rt_output_params)
    else
      ! iteration-pass defaults
      call to_locals (rt_params_t(angle_quad='radau', ntheta=2, nphi=4, &
        solver='linear', include_linelist=.false., include_absmet=.true.))
      rewind (io%input)
      read (io%input, rt_iteration_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, rt_iteration_params)
      end if
      call io%namelist_warning ('rt_iteration_params', iostat)
      call from_locals (rt_iter)
      ! output-pass defaults
      call to_locals (rt_params_t(angle_quad='lobatto', ntheta=4, nphi=4, &
        solver='quadratic', include_linelist=.true., include_absmet=.false.))
      rewind (io%input)
      read (io%input, rt_output_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, rt_output_params)
      end if
      call io%namelist_warning ('rt_output_params', iostat)
      call from_locals (rt_out)
      ! the composition module builds its sampled LTE line table (and masks
      ! departure species out of it) whenever either pass samples the line list
      comp_sample_linelist = rt_iter%sample_linelist .or. rt_out%sample_linelist
    end if
  contains
    subroutine to_locals (p)
      type(rt_params_t), intent(in) :: p
      solver = p%solver
      angle_quad = p%angle_quad
      ntheta = p%ntheta
      nphi = p%nphi
      maxiter = p%maxiter
      sample_back_opac = p%sample_back_opac
      sample_linelist = p%sample_linelist
      include_linelist = p%include_linelist
      include_absmet = p%include_absmet
    end subroutine to_locals
    subroutine from_locals (p)
      type(rt_params_t), intent(inout) :: p
      p%solver = lowercase(solver)
      p%angle_quad = lowercase(angle_quad)
      p%ntheta = ntheta
      p%nphi = nphi
      p%maxiter = maxiter
      p%sample_back_opac = sample_back_opac
      p%sample_linelist = sample_linelist
      p%include_linelist = include_linelist
      p%include_absmet = include_absmet
    end subroutine from_locals
  END SUBROUTINE rt_params_nml

  !===============================================================================
  !> Assign the pass's ray directions and integration weights from its scheme.
  !> The custom/grid/none schemes only exist for the output pass.
  !===============================================================================
  SUBROUTINE setup_quadrature (p, is_output, ds)
    type(rt_params_t), intent(inout), target :: p
    logical, intent(in) :: is_output
    real(8), intent(in) :: ds(3)
    real(8), dimension(:), allocatable :: cust_mu
    character(len=32) :: quad
    logical :: Jnu_weights
    !-----------------------------------------------------------------------------
    ! solver='none' disables the whole output pass: route it through the empty
    ! 'none' quadrature below (regardless of angle_quad, which is left at its
    ! default) so size(rt_out%vec,2)==0 and every downstream "n_mu > 0" guard
    ! short-circuits -- no output RT patches are created.
    quad = p%angle_quad
    if (p%solver == 'none') quad = 'none'
    Jnu_weights = .not. is_output
    select case (quad)
    case ('disk_center', 'set_a2', 'set_a4', 'set_a6', 'set_a8', 'set_b4', 'set_b6', 'set_b8')
      call set_carlson_quadratures (p%angle_quad, p%vec, p%wmu, Jnu_weights=Jnu_weights)
    case ('radau')
      call set_radau_quadratures (p%ntheta, p%nphi, p%vec, p%wmu, Jnu_weights=Jnu_weights)
    case ('lobatto')
      call set_lobatto_quadratures (p%ntheta, p%nphi, p%vec, p%wmu, Jnu_weights=Jnu_weights)
    case ('gauss')
      call set_gauss_quadratures (p%ntheta, p%nphi, p%vec, p%wmu, Jnu_weights=Jnu_weights)
    case ('custom')
      call assert_output ()
      call count_words (custom_mu, p%ntheta)
      allocate (cust_mu(p%ntheta))
      read (custom_mu, *) cust_mu
      call set_custom_quadratures (cust_mu, p%nphi, p%vec, p%wmu)
    case ('grid+', 'gridx', 'grid*')
      call assert_output ()
      call set_grid_quadratures (p%angle_quad, ds, p%vec, p%wmu)
    case ('none')
      call assert_output ()
      allocate (p%vec(3,0))
      allocate (p%wmu(0))
    case default
      call console_error ('quadrature keyword "'// trim(p%angle_quad) // '" not recognized for ' // &
        trim(merge('rt_output_params   ', 'rt_iteration_params', is_output)))
      stop
    end select
  contains
    subroutine assert_output ()
      if (.not. is_output) then
        call console_error ('quadrature keyword "'// trim(p%angle_quad) // &
          '" is only available for rt_output_params')
        stop
      end if
    end subroutine assert_output
  END SUBROUTINE setup_quadrature

END MODULE rt_params_mod
