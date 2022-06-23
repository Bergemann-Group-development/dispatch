module atom_patch_mod
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
    USE opac_patch_mod
    USE rad_patch_mod
    USE rate_patch_mod
    USE atom_mod, only: atom_type, atomic_spectrum_type
    USE trans_mod, only: line_type
    USE mconst, only: CC
    USE ion_mod, only: ion_type
    USE element_mod, only: element_type
    USE bbtrans_mod, only: bbtrans_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, ilambd
    USE bftrans_mod, only: bftrans_type
    USE lines_patch_mod, only:lines_patch_t
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      !type(atom_type), pointer :: atom, output_atom
      type(atom_type) :: atom, output_atom
      class(composition_type), pointer :: composition
      class(atomic_spectrum_type), pointer :: output_sp
      real(8) :: maxerr=huge(1.0)
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: update => atom_patch_update
      procedure :: post_update => atom_patch_post_update
      procedure :: finish => atom_patch_finish
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure :: pre_init_atom
      procedure :: init_atom
      procedure :: normalize_phi
      procedure :: sum_rates
      procedure :: update_populations
      procedure :: exclude_double_continua
      procedure :: exclude_double_lines
    end type atom_patch_t

    real(8), public:: global_maxerr=0
    real(8), public:: convlim=1d-3

    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_patch_pre_init (self)
    class(atom_patch_t), target:: self
    class(rate_patch_t), pointer :: rate
    class(opacity_t), pointer :: opac
    class(lines_patch_t), pointer :: lines
    integer :: m, n, m2, n2
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_pre_init', itimer=itimer)
    call self%m3d_patch_t%pre_init()
    call self%append_to_task_list()
    call self%pre_init_atom()
    call self%set_nv_and_nw(self%atom%nlevel, 2)

    do n=1, self%n_nu
      do m=0, self%n_omega
        if (m==0) then
          allocate(opac)
          opac%atom => self%atom
          opac%kind = 'cont_opac'
          call opac%pre_init_and_link(self, i_mu=m, i_nu=n, sp=self%sp)
        else
          allocate(lines)
          lines%kind = 'lines_patch'
          lines%atom => self%atom
          call lines%pre_init_and_link(self, i_mu=m, i_nu=n, sp=self%sp)

          allocate(rate)
          rate%kind = 'rate_patch'
          rate%atom => self%atom
          call rate%pre_init_and_link(self, i_mu=m, i_nu=n, sp=self%sp)
        endif
      end do
    end do

    do m=0, self%n_out
      m2 = m+self%n_omega
      do n=1, self%n_nu
        n2 = n+self%n_nu
        if (m==0) then
          allocate(opac)
          opac%atom => self%output_atom
          opac%kind = 'cont_opac'
          call opac%pre_init_and_link(self, i_mu=0, i_nu=n2, sp=self%output_sp)
        else
          allocate(lines)
          lines%kind = 'lines_patch'
          lines%atom => self%output_atom
          call lines%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=self%output_sp)
        endif
      end do
    end do
    call self%post_pre_init()
    self%atom%fcheck = self%fcheck
    self%atom%pcheck = self%pcheck
    self%atom%linecheck = self%linecheck
    call trace%end(itimer)
  END SUBROUTINE atom_patch_pre_init

  !===============================================================================
  SUBROUTINE pre_init_atom (self)
    class(atom_patch_t) :: self
    type(atom_type), save, target :: atom, output_atom
    type(atomic_spectrum_type), save, target :: sp
    logical, save :: first_time=.true., atom_read=.false.
    character(len=128), save :: atom_file='./input_multi3d/atoms/atom.ba06'
    character(len=256) :: command
    real(8), save :: qnorm=10.0, logf_cut = -99.9, red_cut=huge(1.0), blue_cut=tiny(1.0)
    integer :: iostat, lunit
    integer, save:: verbose=0
    integer, save :: itimer=0, itimer1=0, n_in_lines=0
    logical, save :: exclude_trace_cont=.true., use_atom_abnd=.true., second_sp=.false.
    integer, save, allocatable :: in_lines(:)
    namelist /atom_params/ atom_file, qnorm, convlim, verbose, exclude_trace_cont, &
    use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines
    namelist /line_mask/ in_lines
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%pre_init_atom', itimer=itimer)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, atom_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, atom_params)
      end if
      call io%namelist_warning ('atom_params', iostat)
      write (io%namelist, atom_params)

      open(newunit=lunit, file='./out/atom_params.nml', delim='apostrophe', form='formatted')
      write (lunit, atom_params)

      if (n_in_lines > 0) then
        allocate(in_lines(n_in_lines))
        rewind (io%input)
        read (io%input, line_mask, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, line_mask)
        end if
        write (lunit, line_mask)
      endif

      close(lunit)

      sp%verbose = verbose
      if (verbose) write (*,*) ' - Reading atom'
      call atom%read(atom_file, qnorm, logf_cut=logf_cut, red_cut=red_cut, blue_cut=blue_cut, in_lines=in_lines, teff=self%teff)
      output_atom = atom
      call atom%create_spectrum(sp)
      call atom%set_active_transitions(sp)
      atom%sp => sp

      second_sp = self%sp%nnu > 0 
      if (second_sp) then
        call output_atom%set_active_transitions(self%sp) 
      endif
      write (command, *) 'cp ', trim(atom_file), ' ./out/atom.txt' 
      call execute_command_line (trim(command), exitstat=iostat)
    endif

    self%atom = atom
    if (second_sp) then
      self%output_atom = output_atom
    else
      self%output_atom = atom
    endif

    if (second_sp) then
      self%output_sp => self%sp
    else
      self%output_sp => sp
    endif

    self%sp => sp

    if (exclude_trace_cont) then
      call self%exclude_double_continua()
    endif

    call self%exclude_double_lines()

    if (use_atom_abnd) then
      call self%composition%overwrite_abundance(atom%id, atom%abnd, atom%awgt)
    endif

    call trace%end(itimer)
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE init_atom (self)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    real(4), dimension(:,:,:), pointer :: ne, temp
    real(4), dimension(:,:,:,:), pointer :: nh
    integer, save:: itimer=0
    integer :: i
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_atom', itimer=itimer)
    atom => self%atom

    allocate(atom%C(atom%nlevel, atom%nlevel))
    allocate(atom%W(atom%nlevel, atom%nlevel))
    allocate(atom%e(atom%nlevel))

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)

    atom%nstar => self%mem(:,:,:,:,self%it,1)
    atom%n     => self%mem(:,:,:,:,self%it,2)

    call os%allocate(atom%totn, shape(temp))
    call os%allocate(atom%doppler_factor, shape(temp))

    atom%doppler_factor = atom%compute_doppler_factor(temp)
    !call atom%compute_line_damping(ne, temp, nh)
    call self%composition%get_totn(atom%id, atom%abnd, atom%totn)
    call atom%compute_lte_populations(temp, ne)
    call atom%compute_cont_factors()
    self%output_atom%n => atom%n
    self%output_atom%nstar => atom%nstar
    self%output_atom%doppler_factor => atom%doppler_factor
    self%output_atom%totn => atom%totn
    do i=1,atom%ncont
      self%output_atom%cont(i)%ni_over_nj_lte => atom%cont(i)%ni_over_nj_lte
    enddo
    call trace%end(itimer)
  END SUBROUTINE init_atom

  !===============================================================================
  SUBROUTINE atom_patch_update (self)
    class(atom_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    real(8), dimension(:,:,:,:), allocatable :: gamma_ij, gamma_ji
    character(len=64) :: text
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call self%pre_update()
    if (self%finished) return
    call trace%begin('atom_patch_t%atom_patch_update', itimer=itimer)

    if (self%iter == 0) then
      write(text,'(A,I2,A,I2,A)') 'LTE Iteration .'
      call io%header (text)

      call self%init_atom()
      self%dtime = self%dtime + self%dtime*ilambd

      nbor => self%link%nbor
      do while (associated(nbor))
        nbtask => nbor%task
        select type(nbtask)
        class is (rate_patch_t)
          nbtask%time = nbtask%time + nbtask%dtime*ilambd
        class is (opacity_t)
          nbtask%dtime = nbtask%dtime
        end select
        nbor => nbor%next
      end do
    else
      write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter, '.'
      call io%header (text)

      if (self%iter == 1) then
        self%dtime = self%dt_fixed
        nbor => self%link%nbor
        do while (associated(nbor))
          nbtask => nbor%task
          select type(nbtask)
          class is (opacity_t)
            nbtask%dtime = nbtask%dt_fixed
          end select
          nbor => nbor%next
        end do
        call self%normalize_phi()
      endif

      allocate(gamma_ij(self%n(1), self%n(2), self%n(3), self%atom%nrad))
      allocate(gamma_ji(self%n(1), self%n(2), self%n(3), self%atom%nrad))

      call self%sum_rates(gamma_ij, gamma_ji)
      call self%update_populations(gamma_ij, gamma_ji)
    endif
  
    call self%post_update()
    call trace%end(itimer)
  END SUBROUTINE atom_patch_update

  !=============================================================================
  SUBROUTINE atom_patch_post_update (self)
    class(atom_patch_t) :: self
    integer :: p(3), i, z
    integer, save :: count=0
    ! --------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_post_update')
    if (self%verbose > 0) then
      p = self%pcheck
      if (all(p > 0)) then
        i = self%lvlcheck
        call io%lock%set
        write(stdout,'(A, 8(ES15.7))') self%m3d_patch_info(), self%atom%abnd, &
          self%atom%n(p(1),p(2),p(3),i), self%atom%totn(p(1),p(2),p(3)), &
          10.0**(self%atom%abnd-12.0)*self%composition%toth(p(1),p(2),p(3))
        write(stdout,'(A19, ES19.7, A10, 3ES15.7)') 'Maxerr =', self%maxerr, &
          ', global =', global_maxerr
        write(stdout,'(1X, A80)') &
        '................................................................................'
        call io%lock%unset
      endif
    end if

    if (self%iter > 0) count = count + 1 
    
    if (count >= self%task_list%dims(3)) then
      if ((global_maxerr > 0) .and. global_maxerr < convlim) then
        print *, self%ipos(3), 'Converged!', count
        call self%finish()
      endif
    endif

    call self%m3d_patch_t%post_update()
    call trace%end()
  END SUBROUTINE atom_patch_post_update

  !===============================================================================
  SUBROUTINE atom_patch_finish (self, rotated)
    use m3d_patch_mod, only: converged
    class(atom_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    logical, optional :: rotated
    integer :: i, j
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_finish')
    self%finished = .true.
    self%time = io%end_time + self%t_offset*self%dtime
    converged = .true.

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needs_me) then
        nbtask => nbor%task
        select type(nbtask)
        class is (m3d_patch_t)
          if (nbtask%i_mu > 0) then
            if (.not. nbtask%out_ray) then
              call nbtask%finish()
            endif
          else
            nbtask%time = io%end_time + (self%t_offset - 1)*nbtask%dtime
            nbtask%iter = nbtask%maxiter
            nbtask%out_ray = .true.
          endif
        end select
      endif
      nbor => nbor%next
    end do
    call trace%end()
  END SUBROUTINE atom_patch_finish

  !===============================================================================
  SUBROUTINE update_populations (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    integer, dimension(3, self%atom%nlevel) :: pos3
    real, dimension(self%atom%nlevel) :: maxerr
    real(4), pointer :: temp, ne, totn
    real(4), dimension(:), pointer :: nh
    real(8), dimension(:), pointer :: Gij, Gji
    real(4), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    integer :: i, j, k, ii
    integer :: cswitch=1
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    atom => self%atom
    maxerr = 0.0d0

    do k=1, self%n(3)
      do j=1, self%n(2)
        do i=1, self%n(1)
          ne    => self%m3d_parent%mem(i,j,k,self%i_ne,            self%it,1)
          temp  => self%m3d_parent%mem(i,j,k,self%i_temp,          self%it,1)
          nh    => self%m3d_parent%mem(i,j,k,self%i_nh1:self%i_nh6,self%it,1)
          totn  => atom%totn  (i,j,k)
          n     => atom%n     (i,j,k,:)
          nstar => atom%nstar (i,j,k,:)
          Gij   => gamma_ij   (i,j,k,:)
          Gji   => gamma_ji   (i,j,k,:)

          call atom%solve_stat_equilibrium(temp, ne, totn, n, nstar, nh, Gij, Gji, cswitch)

          ! keep track of maximum population change
          do ii=1,atom%nlevel
            if (abs(atom%e(ii)) > maxerr(ii)) then
                maxerr(ii) = atom%e(ii)
                pos3(1,ii) = i  
                pos3(2,ii) = j 
                pos3(3,ii) = k
            end if
          end do
        end do
      end do
    end do

    if (global_maxerr == self%maxerr) then
      self%maxerr = maxval(abs(maxerr))
      global_maxerr = self%maxerr
    else
      self%maxerr = maxval(abs(maxerr))
      global_maxerr = max(global_maxerr, self%maxerr)
    endif
  END SUBROUTINE update_populations

  !===============================================================================
  SUBROUTINE normalize_phi (self)
    class(atom_patch_t), target :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(line_type), pointer :: line
    integer :: i, j, i_mu, i_nu
    ! ----------------------------------------------------------------------------
    do i=1, self%atom%nline
      line => self%atom%line(i)
      call os%allocate(line%wphi, self%n, zero=.true.)

      nbor => self%link%nbor
      do while (associated(nbor))
        nbtask => nbor%task
        select type(opac => nbtask)
        class is (lines_patch_t)
          if (.not. opac%out_ray) then
            j = opac%idx_map(i)
            if (j>0) then
              line%wphi = line%wphi + opac%wmu * opac%line(j)%wphi
              call os%deallocate(opac%line(j)%wphi)
            endif
          endif
        class is (opacity_t)
          if (opac%i_mu > 0 .and. .not. opac%out_ray) then
            line%wphi = line%wphi + opac%wmu * opac%wphi(:,:,:,i)
          endif
        end select
        nbor => nbor%next
      end do

      line%wphi = math%pi4 / line%wphi
    end do

    nbor => self%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type(opac => nbtask)
      class is (lines_patch_t)
        if (opac%finished) then
          call os%deallocate(opac%idx_map)
          deallocate(opac%line)
        endif
      end select
      nbor => nbor%next
    end do
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE sum_rates (self, gamma_ij, gamma_ji)
    class(atom_patch_t), target :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    integer :: i, i_mu, i_nu, m(4), p(3)
    real(8), dimension(:,:,:,:), allocatable:: gamma_ij, gamma_ji
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    gamma_ij = 0.0d0
    gamma_ji = 0.0d0

    nbor => self%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type(rate => nbtask)
      class is (rate_patch_t)
        gamma_ij = gamma_ij + rate%wmu * rate%mem(:,:,:,:,rate%it,rate%i_Gij)
        gamma_ji = gamma_ji + rate%wmu * rate%mem(:,:,:,:,rate%it,rate%i_Gji)
      end select
      nbor => nbor%next
    end do

    do i=1, self%atom%nline
      gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * self%atom%line(i)%wphi
      gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * self%atom%line(i)%wphi
    end do
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE exclude_double_continua (self)
    class(atom_patch_t), target :: self
    class(bftrans_type), pointer :: trans
    integer :: i, lower, kc
    logical :: same
    !-----------------------------------------------------------------------------
    do i=1, size(self%composition%trans)
      trans => self%composition%trans(i)
      if (trans%level%g2 > 0) then
        if (uppercase(trans%ion%id, 2) == uppercase(self%atom%id, 2)) then
          do kc=1, self%atom%ncont
            lower = self%atom%cont(kc)%i

            if (self%atom%ev(lower)==0) then
              same = trans%level%xion == self%atom%ev(lower)
            else
              same = abs(trans%level%xion / self%atom%ev(lower) - 1) < 1d-2
            endif

            if (same) then
              trans%exclude = .true.
              print *, 'excluding transition', trans%label
            endif
          end do
        endif
      endif
    end do
  END SUBROUTINE exclude_double_continua

  !===============================================================================
  SUBROUTINE exclude_double_lines (self)
    class(atom_patch_t), target :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(bbtrans_type), pointer :: ll_line, line
    integer :: i, j, l
    logical :: same
    !-----------------------------------------------------------------------------
    call self%composition%element_by_id(trim(self%atom%id), element)
    if (.not. associated(element)) return
    if (all(element%ions%nline == 0)) return

    do j=lbound(element%ions,1), ubound(element%ions,1) - 1
      ion => element%ions(j)
      if (ion%nline == 0) cycle

      do l=1, ion%nline
        ll_line => ion%line(l)
        
        do i=1, self%atom%nline
          line => self%atom%line(i)
          if (abs(line%nu0 / ll_line%nu0 -1) < 1e-6 &
              .and. abs(line%Aji / ll_line%Aji - 1) < 1e-4) then
            ll_line%active = .false.
            print '(A, F10.2)', 'excluding linelist line at ', CC / ll_line%nu0 *1e8
          endif
        end do
      end do
    end do
  END SUBROUTINE exclude_double_lines
  
  !===============================================================================
  SUBROUTINE atom_patch_nbor_relations(self, another, needed, needs_me, download)
    class(atom_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_patch_t%nbor_relations')
    call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)
  
    if (needed .or. needs_me) then
      select type (another)
      class is (lines_patch_t)
        needed = .false.
        if (associated(another%atom)) then
          if (another%atom%id /= self%atom%id) then
            needs_me = .false.
          endif
        else
          needs_me = .false.
        endif
      class is (opacity_t)
        needed = .false.
        if (associated(another%atom)) then
          if (another%atom%id /= self%atom%id) then
            needs_me = .false.
          endif
        else
          needs_me = .false.
        endif
      class is (rate_patch_t)
        if (another%atom%id /= self%atom%id) then
          needed = .false.; needs_me = .false.
        endif
      class default
        if (another%id == self%m3d_parent%id) then
          needs_me = .false.
        else
          needed = .false. ; needs_me = .false.
        endif
      end select
    endif
  
    call trace%end()
  END SUBROUTINE atom_patch_nbor_relations

  !===============================================================================
  pure function uppercase(string,ss)
    ! converts the first ss characters of string to uppercase and returns them as a new string.
    character(len=*),intent(in) :: string
    integer,intent(in) :: ss
    character(len=ss) :: uppercase
    integer :: i,k
    character(len=26),parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    character(len=26),parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
    !-----------------------------------------------------------------------------
    uppercase(:)=' '
    uppercase(1:ss)=string(1:ss)
    do i=1,ss
    k=index(lcase,string(i:i))
    if (k > 0) then
        uppercase(i:i)=ucase(k:k)
    else
        uppercase(i:i)=string(i:i)
    end if
    end do
  end function uppercase

  end module
