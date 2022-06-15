module back_line_patch_mod
    USE os_mod
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
    USE m3d_patch_mod, only: m3d_patch_t
    USE ray_patch_mod, only: ray_patch_t
    USE composition_mod
    USE element_mod, only: element_type
    USE ion_mod, only: ion_type
    USE bbtrans_mod, only: bbtrans_type
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(ray_patch_t) :: back_line_patch_t
      class(composition_type), pointer :: composition
      integer :: nline, ired, iblue
      real(8) :: dlam
    contains
      procedure :: pre_init => back_line_patch_pre_init
      procedure :: update => back_line_patch_update
      !procedure :: post_update => back_line_post_update
      !procedure :: finish => back_line_finish
      procedure :: nbor_relations => back_line_patch_nbor_relations
    end type back_line_patch_t

    !-----------------------------------------------------------------------------
  
  contains

  !===============================================================================
  SUBROUTINE back_line_patch_pre_init (self)
    class(back_line_patch_t), target :: self
    integer :: i, nv, nw, iostat, lunit
    logical, save :: first_time=.true.
    real(8), save :: dlam=0.5
    namelist /linelist_params/ dlam
    !---------------------------------------------------------------------------
    call self%ray_patch_t%pre_init()
    call self%append_to_task_list()

    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, linelist_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, linelist_params)
      end if
      call io%namelist_warning ('linelist_params', iostat)
      write (io%namelist, linelist_params)
    endif

    self%dlam = dlam

    nw = 2
    nv = self%fe-self%fs+1
    call self%set_nv_and_nw (nv,nw)
    self%ng = 0
    call self%post_pre_init()
  END SUBROUTINE back_line_patch_pre_init

  !===============================================================================
  SUBROUTINE back_line_patch_update (self)
    class(back_line_patch_t) :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(bbtrans_type), pointer :: line
    real(4), dimension(:,:,:), pointer, contiguous :: tothi, ne, temp, adamp, inv_df
    real(4), dimension(:,:,:), pointer, contiguous :: vel_dop, nh6, inv_part, unsold_temp, stark_pe
    real(8), dimension(:,:,:), pointer, contiguous :: EK_over_temp
    real(4), dimension(:,:,:,:), pointer, contiguous :: nh, phi, eta, chi
    real(8), dimension(:), pointer, contiguous :: nu, h4pi_nu
    integer :: i, j, k, l, f, ired, iblue, nactive, iactive
    real(8) :: cnt_ev, gf, awgt
    real(4) :: scalar
    integer, save:: itimer=0, itimer2=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    call trace%begin('back_line_patch_t%back_line_patch_update', itimer=itimer)

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
    nh6  => nh(:,:,:,6)
    nu   => self%sp%nu(self%fs:self%fe)
    EK_over_temp => self%composition%EK_over_temp

    call os%allocate(inv_df, self%n)
    call os%allocate(tothi, self%n)
    call os%allocate(adamp, self%n)
    call os%allocate(vel_dop, self%n)
    call os%allocate(inv_part, self%n)
    call os%allocate(unsold_temp, self%n)
    call os%allocate(stark_pe, self%n)
    call os%allocate(h4pi_nu, shape(nu))

    call self%comp_vproj(self%vproj)
    tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
    stark_pe = ne * bk * temp**one_sixth
    h4pi_nu = h_4pi * nu

    nactive = 0
    do i=1, self%composition%nel
      element => self%composition%elements(i)
      if (all(element%ions%nline == 0)) cycle

      do j=lbound(element%ions,1), ubound(element%ions,1) - 1
        ion => element%ions(j)
        if (ion%nline == 0) cycle

        do l=1, ion%nline
          line => ion%line(l)
          if (.not. line%active) cycle
          call line%get_range(nu, self%dlam, ired, iblue)
          if (iblue == 0) cycle
          nactive = nactive + 1
        end do
      end do
    end do

    !print *, nactive, ' active lines in this region'

    iactive = 0
    do i=1, self%composition%nel
        element => self%composition%elements(i)
        if (all(element%ions%nline == 0)) cycle

        awgt = element%awgt * UU

        scalar = CC / sqrt(2*BK / awgt)
        inv_df = scalar / sqrt(temp)

        vel_dop = (1 + self%vproj) * inv_df

        scalar = 8.0*BK / PI * (1.0/MH + 1.0/awgt)
        unsold_temp =  scalar * temp

        do j=lbound(element%ions,1), ubound(element%ions,1) - 1
            ion => element%ions(j)
            if (ion%nline == 0) cycle

            ! ion%part is actually then inverse partition
            inv_part = ion%f * self%composition%toth * ion%part

            do l=1, ion%nline
                line => ion%line(l)
                if (.not. line%active) cycle
                call line%get_range(nu, self%dlam, ired, iblue)
                if (iblue == 0) cycle

                iactive = iactive + 1
                !print *, 'working on line ', iactive, ' of ', nactive 

                eta  => self%mem(:,:,:,ired:iblue,1,self%i_eta)
                chi  => self%mem(:,:,:,ired:iblue,1,self%i_chi)

                call trace%begin('bbtrans_type%damping', itimer=itimer2)
                adamp = 0
                call line%damping(adamp, temp, ne, tothi, unsold_temp, stark_pe)
                if (element%id == 'H') then
                  call line%dampingH(adamp, temp, ne, nh6)
                endif
                call trace%end(itimer2)
                call line%add_lte_opacities(nu(ired:iblue), h4pi_nu(ired:iblue), vel_dop, inv_df, EK_over_temp, inv_part, adamp, eta, chi)
            end do
        end do
    end do

    call os%deallocate(inv_df)
    call os%deallocate(tothi)
    call os%deallocate(adamp)
    call os%deallocate(vel_dop)
    call os%deallocate(inv_part)
    call os%deallocate(unsold_temp)
    call os%deallocate(stark_pe)
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE back_line_patch_update

  !===============================================================================
  SUBROUTINE back_line_patch_nbor_relations(self, another, needed, needs_me, download)
    class(back_line_patch_t) :: self
    class(task_t), pointer :: another        ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                            ! does it need me, should I download?
    class(task_t), pointer :: other
    logical :: default
    ! ----------------------------------------------------------------------------
    call trace%begin ('lines_patch_t%nbor_relations')
    call self%ray_patch_t%nbor_relations(another, needed, needs_me, download)

    if (needed .or. needs_me) then
      select type(another)
      class is (m3d_patch_t)
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        elseif (another%kind == 'net_opac') then
          needed = .false.
        else
          needed = .false.; needs_me = .false.
        endif
      class default
        needed = .false.; needs_me = .false.
      end select
    endif
    call trace%end()
  END SUBROUTINE back_line_patch_nbor_relations

  end module
