MODULE back_line_patch_mod
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
    USE m3d_patch_mod, only: m3d_patch_t, hash_string, count_words
    USE ray_patch_mod, only: ray_patch_t
    USE composition_mod
    USE element_mod, only: element_type, isotope_type
    USE ion_mod, only: ion_type
    USE molecule_mod, only: molec_type
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
      procedure :: nbor_relations => back_line_patch_nbor_relations
      procedure :: mask_linelist
      procedure :: init_nbors_hash
    end type back_line_patch_t

    !-----------------------------------------------------------------------------

    character(128), save :: exclude_elements=''
  
  contains

  !===============================================================================
  SUBROUTINE back_line_patch_pre_init (self)
    class(back_line_patch_t), target :: self
    integer :: i, nv, nw, iostat, lunit
    logical, save :: first_time=.true.
    real(8), save :: dlam=0.5
    integer, save :: itimer=0
    namelist /linelist_params/ dlam, exclude_elements
    !---------------------------------------------------------------------------
    call trace%begin('back_line_t%pre_init', itimer=itimer)
    self%kind = 'back_line'
    call self%ray_patch_t%pre_init()

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
    call trace%end(itimer)
  END SUBROUTINE back_line_patch_pre_init

  !===============================================================================
  SUBROUTINE back_line_patch_update (self)
    class(back_line_patch_t) :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(molec_type), pointer :: mol
    class(isotope_type), pointer :: iso
    class(bbtrans_type), pointer :: line
    real(4), dimension(:,:,:), pointer, contiguous :: tothi, ne, temp, adamp, inv_df, vproj
    real(4), dimension(:,:,:), pointer, contiguous :: nh6, inv_part, unsold_temp, stark_pe
    real(8), dimension(:,:,:), pointer, contiguous :: EK_over_temp
    real(4), dimension(:,:,:), pointer, contiguous :: iso_unsold_temp, iso_inv_df, iso_inv_part
    real(4), dimension(:,:,:), pointer, contiguous :: unsoldT, idopf, ipart, vdop
    real(4), dimension(:,:,:,:), pointer, contiguous :: phi
    real(4), dimension(:,:,:,:), pointer :: nh, eta, chi
    real(8), dimension(:), pointer, contiguous :: nu, hnu
    integer :: i, j, k, l, f, ired, iblue, nactive, prev_iso
    real(8) :: cnt_ev, gf, awgt, iso_awgt
    real(4) :: scalar
    integer, save:: itimer=0, itimer2=0
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('back_line_patch_t%back_line_patch_update', itimer=itimer)

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
    nh6  => nh(:,:,:,6)
    nu   => self%sp%nu(self%fs:self%fe)
    EK_over_temp => self%composition%EK_over_temp

    call os%allocate(inv_df, self%n)
    call os%allocate(iso_inv_df, self%n)
    call os%allocate(tothi, self%n)
    call os%allocate(adamp, self%n)
    call os%allocate(inv_part, self%n)
    call os%allocate(iso_inv_part, self%n)
    call os%allocate(unsold_temp, self%n)
    call os%allocate(iso_unsold_temp, self%n)
    call os%allocate(stark_pe, self%n)
    call os%allocate(vproj, self%n)

    call self%comp_vproj(vproj)
    tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
    stark_pe = ne * bk * temp**one_sixth

    call self%mask_linelist()

    nactive = 0
    do i=1, self%composition%nel
      element => self%composition%elements(i)
      if (.not. any(element%ions%in_linelist)) cycle

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

    do i=1, self%composition%nel
        element => self%composition%elements(i)
        if (.not. any(element%ions%in_linelist)) cycle

        awgt = element%awgt * UU

        scalar = CC / sqrt(2*BK / awgt)
        inv_df = scalar / sqrt(temp)

        scalar = 8.0*BK / PI * (1.0/MH + 1.0/awgt)
        unsold_temp =  scalar * temp
        prev_iso = 0

        do j=lbound(element%ions,1), ubound(element%ions,1) - 1
            ion => element%ions(j)
            if (ion%nline == 0) cycle

            ! ion%part is actually the inverse partition
            inv_part = ion%f * self%composition%toth * ion%part

            do l=1, ion%nline
                line => ion%line(l)
                if (.not. line%active) cycle
                call line%get_range(nu, self%dlam, ired, iblue)
                if (iblue == 0) cycle

                eta => self%mem(:,:,:,ired:iblue,1,self%i_eta)
                chi => self%mem(:,:,:,ired:iblue,1,self%i_chi)

                if (line%iso_idx /= prev_iso) then
                  prev_iso = line%iso_idx
                  iso => element%isotope(line%iso_idx)
                  iso_awgt = iso%mass * UU
                  scalar = CC / sqrt(2*BK / iso_awgt)
                  iso_inv_df = scalar / sqrt(temp)
                  scalar = 8.0*BK / PI * (1.0/MH + 1.0/iso_awgt)
                  iso_unsold_temp =  scalar * temp
                  iso_inv_part = inv_part * iso%f
                  unsoldT => iso_unsold_temp
                  idopf   => iso_inv_df
                  ipart   => iso_inv_part
                elseif (line%iso_idx == 0) then
                  unsoldT => unsold_temp
                  idopf   => inv_df
                  ipart   => inv_part
                endif

                call trace%begin('bbtrans_type%damping', itimer=itimer2)
                adamp = 0
                call line%damping(adamp, temp, ne, tothi, unsoldT, stark_pe)
                if (element%id == 'H') then
                  call line%dampingH(adamp, temp, ne, nh6)
                endif
                call trace%end(itimer2)
                call line%add_lte_opacities(nu(ired:iblue), &
                      vproj, idopf, EK_over_temp, ipart, adamp, eta, chi)
            end do
        end do
    end do

    do i=1, self%composition%nmolec
      mol => self%composition%molecules(i)
      if (.not. mol%in_linelist) cycle

      awgt = 0
      do j=1, size(mol%ions)
        awgt = awgt + mol%ions(j)%p%awgt
      end do
      awgt = awgt * UU

      scalar = CC / sqrt(2*BK / awgt)
      inv_df = scalar / sqrt(temp)

      scalar = 8.0*BK / PI * (1.0/MH + 1.0/awgt)
      unsold_temp =  scalar * temp

      inv_part = mol%f * self%composition%toth / mol%irwin_partf(dble(temp))

      do l=1, mol%nline
          line => mol%line(l)
          if (.not. line%active) cycle
          call line%get_range(nu, self%dlam, ired, iblue)
          if (iblue == 0) cycle

          eta  => self%mem(:,:,:,ired:iblue,1,self%i_eta)
          chi  => self%mem(:,:,:,ired:iblue,1,self%i_chi)

          call trace%begin('bbtrans_type%damping', itimer=itimer2)
          adamp = 0
          call line%damping(adamp, temp, ne, tothi, unsold_temp, stark_pe)
          call trace%end(itimer2)
          call line%add_lte_opacities(nu(ired:iblue), &
                vproj, inv_df, EK_over_temp, inv_part, adamp, eta, chi)
      end do
    end do

    call os%deallocate(inv_df)
    call os%deallocate(iso_inv_df)
    call os%deallocate(tothi)
    call os%deallocate(adamp)
    call os%deallocate(inv_part)
    call os%deallocate(iso_inv_part)
    call os%deallocate(unsold_temp)
    call os%deallocate(iso_unsold_temp)
    call os%deallocate(stark_pe)
    call os%deallocate(vproj)
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE back_line_patch_update

  !===============================================================================
  SUBROUTINE mask_linelist (self)
    class(back_line_patch_t) :: self
    class(element_type), pointer :: element
    class(molec_type), pointer :: mol
    integer :: n_excluded, i, ii, p
    character(len=3) :: id
    logical :: new_word
    ! ----------------------------------------------------------------------------
    call count_words(exclude_elements, n_excluded)
    p = 1
    do i=1, n_excluded
      id = ''
      do ii=1, 3
        new_word = .false.
        do while (exclude_elements(p:p) == ' ' .and. p < len(exclude_elements))
          new_word = .true.
          p = p + 1
        end do
        if (new_word) exit
        id(ii:ii) = exclude_elements(p:p)
        p = p + 1
      end do
      call self%composition%element_by_id(id, element)
      if (associated(element)) then
        element%ions%in_linelist = .false.
      else
        call self%composition%molecule_by_id(id, mol)
        if (associated(mol)) then
          mol%in_linelist = .false.
        endif
      endif
    end do
  END SUBROUTINE mask_linelist

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
  
  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(back_line_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: i, i_nu, i_mu, key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('back_line_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)
  
    key(1:3) = self%key(1:3)
    key(4) = hash_string('net_opac')
    call self%check_nbor(key, nbors)
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE
