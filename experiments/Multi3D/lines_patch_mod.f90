module lines_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE omp_timer_mod
    !USE index_mod
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    USE scaling_mod
    USE trace_mod
    USE mpi_mod
    !USE omp_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE mconst
    USE opac_patch_mod
    USE rad_patch_mod
    USE rate_patch_mod
    USE atom_var, only: atom_type
    USE trans_mod, only: line_type
    USE spectrum_mod, only:spectrum_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, ilambd
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: lines_patch_t
      class(atom_type), pointer :: atom
    contains
      procedure :: pre_init => lines_pre_init
      !procedure :: update => linelist_update
      !procedure :: post_update => linelist_post_update
      !procedure :: finish => linelist_finish
      !procedure :: nbor_relations => linelist_nbor_relations
    end type lines_patch_t

    !-----------------------------------------------------------------------------
  
  contains

  !===============================================================================
  SUBROUTINE lines_pre_init (self)
    class(atom_patch_t) :: self
    class(line_type), pointer :: line
    integer :: i, iostat
    logical, save :: first_time=.true.
    integer, save :: nq=9, qmax=10, q0=0
    namelist /line_params/ nq, qmax, q0
    !---------------------------------------------------------------------------
    call self%m3d_patch_t%pre_init

    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, line_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, line_params)
      end if
      call io%namelist_warning ('line_params', iostat)
      write (io%namelist, line_params)
    endif

    do i=1, self%atom%nline
      line => self%atom%line(i)

      if (allocated(line%nu)) deallocate(line%nu)
      line%nnu = nq
      allocate(line%nu(line%nnu))

      allocate(q(nq))
      allocate(wq(nq))
  
      if (qmax.le.q0) then
        ! linear spacing
        dq = qmax / (nnu-1)
        q(1)= 0
        do j=2,nq
            q(j) = q(j-1) + dq
        end do
  
      else if(qmax.ge.0.0 .and. q0.ge.0.0) then
  
        ! constants
        a=10.**(q0+.5)
        xmax=log10(a*max(0.5,qmax-q0-0.5))
    
        dx=xmax/(nnu-1)
        do j=1,nnu
            x=(j-1)*dx
            x10=10.**x
            q(j)=x+(x10-1./x10)/a
        end do
    
      end if
  
      ! convert to [Hz]
      line%nu = line%nu0*(1.0+q*qnorm*KM_TO_CM/cc)
    end do

  END SUBROUTINE lines_pre_init

  !SUBROUTINE linelist_update ()
  !END SUBROUTINE

  end module
