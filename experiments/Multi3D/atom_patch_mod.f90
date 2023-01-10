MODULE atom_patch_mod
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
    USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE atom_cont_patch_mod, only: atom_cont_patch_t
    USE angle_rad_patch_mod, only: rad_patch_t
    USE atom_rate_patch_mod, only: atom_rate_t
    USE atom_mod, only: atom_type, atomic_spectrum_type
    USE trans_mod, only: line_type
    USE mconst, only: CC
    USE ion_mod, only: ion_type, particle_type
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
    USE bbtrans_mod, only: bbtrans_type
    USE composition_mod, only: composition_type
    USE m3d_patch_mod, only: m3d_patch_t, maxiter, check_NaN, recl, cswitch, &
    global_maxerr, convlim, atom_update_counter, hash_string, comoving_frame, global_pos
    USE bftrans_mod, only: bftrans_type
    USE atom_line_patch_mod, only:atom_line_t, second_sp
    USE net_opac_patch_mod, only:net_opac_t
    implicit none
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(m3d_patch_t) :: atom_patch_t
      !type(atom_type), pointer :: atom, output_atom
      type(atom_type) :: atom, output_atom
      class(composition_type), pointer :: composition
      class(atomic_spectrum_type), pointer :: output_sp
      real(8) :: maxerr=huge(1.0)
      real(8), dimension(:), allocatable :: lvlerr
      logical :: update_skipped=.false., atom_allocated=.false.
    contains
      procedure :: pre_init => atom_patch_pre_init
      procedure :: update => atom_patch_update
      procedure :: finish => atom_patch_finish
      procedure :: nbor_relations => atom_patch_nbor_relations
      procedure :: pre_init_atom
      procedure :: init_atom
      procedure :: update_populations
      procedure :: check_skip
      procedure :: exclude_double_continua
      procedure :: exclude_double_lines
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: save_abundance_to_file
      procedure :: init_nbors_hash
    end type atom_patch_t

    integer, save :: epos(3), elvl
    logical, save :: level_partf=.true.
    logical, save :: exclude_from_line_list=.true.

    !-----------------------------------------------------------------------------
  
  contains
    
  !===============================================================================
  SUBROUTINE atom_patch_pre_init (self)
    class(atom_patch_t), target:: self
    class(atom_rate_t), pointer :: rate
    class(atom_cont_patch_t), pointer :: cont
    class(atom_line_t), pointer :: lines
    integer :: m, n, m2, n2
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_pre_init', itimer=itimer)
    self%kind = 'atom_patch'
    call self%m3d_patch_t%pre_init()
    call self%pre_init_atom()
    if (.not. self%atom_allocated) then
      self%maxiter = -1 ! prevents addition to task list
      call trace%end(itimer)
      return
    endif
    call self%set_nv_and_nw(self%atom%nlevel, 2)
    self%static_mem=.true.

    do n=1, self%n_nu
      allocate(cont)
      cont%atom => self%atom
      call cont%pre_init_and_link(self, i_mu=0, i_nu=n, sp=self%sp)

      if (comoving_frame) then
        allocate(lines)
        lines%atom => self%atom
        call lines%pre_init_and_link(self, i_mu=0, i_nu=n, sp=self%sp)
        allocate(rate)
        rate%atom => self%atom
        call rate%pre_init_and_link(self, i_mu=0, i_nu=n, sp=self%sp)
      else
        do m=1, self%n_omega
          allocate(lines)
          lines%atom => self%atom
          call lines%pre_init_and_link(self, i_mu=m, i_nu=n, sp=self%sp)
        end do

        do m=1, self%n_omega
            allocate(rate)
            rate%atom => self%atom
            call rate%pre_init_and_link(self, i_mu=m, i_nu=n, sp=self%sp)
        end do
      endif
    end do

    do n=1, self%n_nu
      n2 = n+self%n_nu
      do m=0, self%n_out
        m2 = m+self%n_omega
        if (m==0) then
          allocate(cont)
          cont%atom => self%output_atom
          call cont%pre_init_and_link(self, i_mu=0, i_nu=n2, sp=self%output_sp)
        else
          allocate(lines)
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
    character(len=128), save :: atom_file=''
    character(len=256) :: command
    real(8), save :: qnorm=10.0, logf_cut = -99.9, red_cut=huge(1.0), blue_cut=tiny(1.0)
    real(8), save :: abundance=-99.9
    integer :: iostat
    integer, save :: itimer=0, itimer1=0, n_in_lines=0
    logical, save :: exclude_trace_cont=.true., use_atom_abnd=.true.
    integer, save, allocatable :: in_lines(:)
    namelist /atom_params/ atom_file, qnorm, convlim, exclude_trace_cont, level_partf, &
    use_atom_abnd, logf_cut, red_cut, blue_cut, n_in_lines, abundance, exclude_from_line_list
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
      if (atom_file == '') then
        print *, 'WARNING: no atom specified'
        call trace%end(itimer)
        return
      endif

      if (n_in_lines > 0) then
        allocate(in_lines(n_in_lines))
        rewind (io%input)
        read (io%input, line_mask, iostat=iostat)
        if (iostat > 0) then
          rewind (io%input)
          read (io%input, line_mask)
        end if
        write (io%namelist, line_mask)
      endif

      write (*,*) ' - Reading atom'
      call atom%read(atom_file, qnorm, logf_cut=logf_cut, red_cut=red_cut, &
      blue_cut=blue_cut, in_lines=in_lines, read_coll=self%maxiter>0)
      if (abundance > -99) atom%abnd = abundance
      call atom%create_spectrum(sp)
      !output_atom = atom
      call atom%copy_to(output_atom)
      call atom%set_active_transitions(sp)
      atom%sp => sp

      second_sp = self%sp%nnu > 0 
      if (second_sp) then
        call output_atom%set_active_transitions(self%sp) 
      endif
      call os%copy_file(trim(atom_file), trim(io%run_dir)//'atom.txt')
    endif

    if (atom_file == '') then
      call trace%end(itimer)
      return
    endif

    ! self%atom = atom
    call atom%copy_to(self%atom)
    self%atom_allocated = .true.

    if (second_sp) then
      !self%output_atom = output_atom
      call output_atom%copy_to(self%output_atom)
    else
      !self%output_atom = atom
      call atom%copy_to(self%output_atom)
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

    if (use_atom_abnd .or. abundance > -99) then
      call self%composition%overwrite_abundance(atom%id, atom%abnd, atom%awgt)
    else
      self%atom%abnd = self%atom%abnd + self%composition%FeH
    endif

    call trace%end(itimer)
  END SUBROUTINE pre_init_atom

  !===============================================================================
  SUBROUTINE init_atom (self)
    class(atom_patch_t), target :: self
    class(atom_type), pointer :: atom
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    real(4), dimension(:,:,:), pointer :: temp, ne
    real(4), dimension(:,:,:,:), pointer :: nh
    integer, save:: itimer=0
    integer :: i
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_atom', itimer=itimer)
    atom => self%atom

    allocate(atom%C(atom%nlevel, atom%nlevel))
    allocate(atom%W(atom%nlevel, atom%nlevel))
    allocate(atom%e(atom%nlevel))
    allocate(self%lvlerr(atom%nlevel))

    ne   => self%m3d_parent%mem(:,:,:,self%i_ne,            self%it,1)
    temp => self%m3d_parent%mem(:,:,:,self%i_temp,          self%it,1)
    nh   => self%m3d_parent%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)

    atom%nstar => self%mem(:,:,:,:,self%it,1)
    atom%n     => self%mem(:,:,:,:,self%it,2)

    call os%allocate(atom%totn, shape(temp), label='totn')
    call os%allocate(atom%doppler_factor, shape(temp), label='doppler_factor')

    atom%doppler_factor = atom%compute_doppler_factor(temp)
    call self%composition%get_totn(atom%id, atom%abnd, atom%totn, element, molec)
    if (.not. associated(element) .and. .not.associated(molec)) then
      call self%save_abundance_to_file()
    endif
    if (level_partf) nullify(element)
    call atom%compute_lte_populations(temp, ne, real(self%composition%dxi), element, molec)
    call atom%compute_cont_factors(temp, real(self%composition%dxi), ne)
    self%output_atom%n => atom%n
    self%output_atom%nstar => atom%nstar
    self%output_atom%doppler_factor => atom%doppler_factor
    self%output_atom%totn => atom%totn
    do i=1,atom%ncont
      self%output_atom%cont(i)%ni_over_nj_lte => atom%cont(i)%ni_over_nj_lte
    enddo
    
    recl = max(recl, kind(self%mem)*product(shape(self%authority_mem)))
    call trace%end(itimer)
  END SUBROUTINE init_atom

  !===============================================================================
  SUBROUTINE atom_patch_update (self)
    class(atom_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(line_type), pointer :: line
    real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij, gamma_ji
    integer :: i
    integer, save:: itimer=0, count=0
    real(8) :: global
    !-----------------------------------------------------------------------------
    if (self%is_finished()) return
    call trace%begin('atom_patch_t%atom_patch_update', itimer=itimer)

    if (self%iter == 0) then
      call self%init_atom()
    else
      if (self%iter == 1) then
        self%dtime = self%dt_fixed

        ! Normalize wphi
        do i=1, self%atom%nline
          self%atom%line(i)%wphi = math%pi4 / self%atom%line(i)%wphi
        end do
      endif

      if (.not. self%atom%stiff) then
        gamma_ij => self%atom%gamma_ij
        gamma_ji => self%atom%gamma_ji

        do i=1, self%atom%nline
          gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * self%atom%line(i)%wphi
          gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * self%atom%line(i)%wphi
        end do

        call self%update_populations(gamma_ij, gamma_ji)
        self%update_skipped = .false.
      else
        self%update_skipped = .true.
      endif

      call os%deallocate(self%atom%gamma_ij, label='atom%gamma_ij')
      call os%deallocate(self%atom%gamma_ji, label='atom%gamma_ji')
    endif
    
    if (self%iter > 0) then
      count = atom_update_counter%update(+1)
    endif
    if (count >= self%task_list%dims(3)) then
      global = global_maxerr%update(0.d0)
      if ((global > 0) .and. global < convlim .and. self%maxerr < convlim) then
        write(stdout, *) self%ipos(3)+1, 'Converged!', count
        call self%finish()
      endif
    endif
  
    call self%post_update()

    if (self%atom%stiff .and. .not. self%is_finished()) then
      call self%check_skip()
    endif
    call trace%end(itimer)
  END SUBROUTINE atom_patch_update

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(atom_patch_t):: self
    character(len=64) :: text
    integer :: i, p(3), ipos
    real(4), dimension(4) :: arr
    real(8) :: global
    !-----------------------------------------------------------------------------
    if (self%iter == 0) then
      write(text,'(A,I2,A,I2,A)') 'LTE Iteration .'
    else
      write(text,'(A,I2,A,I2,A)') 'Iteration ', self%iter, '/', self%maxiter, '.'
    endif
    call io%header(text)
    if (self%update_skipped) then
      call self%print_to_stdout(text='update skipped')
    elseif (self%m3d_patch_t%check_in_bounds()) then
      i = self%lvlcheck
      p = self%pcheck
      arr(1) = self%lvlerr(i)
      arr(2) = self%atom%n(p(1),p(2),p(3),i)
      arr(3) = self%atom%totn(p(1),p(2),p(3))
      arr(4) = 10.0**(self%atom%abnd-12.0)*self%composition%toth(p(1),p(2),p(3))
      call self%print_to_stdout(arr=arr)
    else
      call self%print_to_stdout()
    endif
    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    write(stdout,'(A19, ES19.7, A10, ES15.7, A9, I2, 2X, 4I4)') 'local err =', self%maxerr, &
      ', global =', global, ', ipos = ', ipos+1, elvl, epos
    write(stdout,'(1X, A80)') &
    '................................................................................'
  END SUBROUTINE print_post_update_info

  !===============================================================================
  SUBROUTINE check_skip (self)
    class(atom_patch_t) :: self
    class(nbor_t), pointer :: nbor
    class(task_t), pointer :: nbtask
    integer :: iter, atom_iter
    real(8) :: time
    ! ----------------------------------------------------------------------------
    time=huge(1.d0)
    nbor => self%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type(nbtask)
      class is (atom_rate_t)
        if (.not. nbtask%out_ray) then
          time = min(time, nbtask%rad%time)
        endif
      end select
      nbor => nbor%next
    end do
    iter = max(0, floor(time - self%time))
    atom_iter = max(0, min(self%maxiter - self%iter - 1, iter))
    call self%skip_iterations(atom_iter)
    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needs_me) then
        nbtask => nbor%task
        select type(nbtask)
        class is (m3d_patch_t)
          if (.not. nbtask%out_ray .and. .not. nbtask%is_finished()) then
            call nbtask%skip_iterations(iter)
          endif
        end select
      endif
      nbor => nbor%next
    end do
  END SUBROUTINE check_skip

  !===============================================================================
  SUBROUTINE atom_patch_finish (self)
    USE m3d_patch_mod, only: converged
    class(atom_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    integer :: i, j
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%atom_patch_finish')
    call self%m3d_patch_t%finish()
    if (self%iter < self%maxiter) converged = .true.
    !call self%output()

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needs_me) then
        nbtask => nbor%task
        select type(nbtask)
        class is (m3d_patch_t)
          if (.not. nbtask%out_ray .and. associated(nbtask%m3d_parent)) then
            if (nbtask%m3d_parent%id == self%id) then
              if (.not. nbtask%finished) then
                call nbtask%finish()
                call nbtask%free_unused_memory()
              elseif (.not. nbtask%static_mem) then
                call nbtask%mem_lock(1)%set ('atom_patch_finish')
                call nbtask%deallocate_mem(self, print=.true.)
                call nbtask%mem_lock(1)%unset ('atom_patch_finish')
              endif
            endif
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
    integer, dimension(3) :: ep
    real(4), pointer :: temp, ne, totn
    real(4), dimension(:), pointer :: nh
    real(8), dimension(:), pointer :: Gij, Gji
    real(4), dimension(:), pointer :: n, nstar
    real(8), dimension(:,:,:,:), target :: gamma_ij, gamma_ji
    integer :: i, j, k, ii, ipos, il
    !integer :: cswitch=1
    integer, save :: itimer=0
    real(8) :: global
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%update_populations', itimer=itimer)
    atom => self%atom
    self%lvlerr = 0.0d0
    self%maxerr = 0.0d0

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

          if (totn > 0) then
            call atom%solve_stat_equilibrium(temp, ne, totn, n, nstar, nh, Gij, Gji, cswitch, check_NaN)

            ! keep track of maximum population change
            do ii=1,atom%nlevel
              if (abs(atom%e(ii)) > self%lvlerr(ii)) then
                  self%lvlerr(ii) = abs(atom%e(ii))
                  if (self%lvlerr(ii) > self%maxerr) then
                    self%maxerr = self%lvlerr(ii)
                    ep = [i, j, k + self%zs]
                    il = ii
                  endif
              end if
            end do
          endif
        end do
      end do
    end do

    global = global_maxerr%update(0.d0)
    ipos = global_pos%update(+0)
    if (ipos == self%ipos(3) .or. self%maxerr > global) then
      if (ipos /= self%ipos(3)) call global_pos%reset(self%ipos(3))
      call global_maxerr%reset(self%maxerr)
      epos = ep
      elvl = il
    endif
    call trace%end(itimer)
  END SUBROUTINE update_populations

  !===============================================================================
  SUBROUTINE exclude_double_continua (self)
    class(atom_patch_t), target :: self
    class(bftrans_type), pointer :: trans
    integer :: i, lower, kc
    logical :: same
    logical, save :: first_time=.true.
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
              if (first_time) then
                print *, 'excluding transition ', trans%label
                first_time = .false.
              endif
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
    class(molec_type), pointer :: molec
    class(ion_type), pointer :: ion
    class(bbtrans_type), pointer :: ll_line, line
    integer :: i, j, l
    logical, save :: first_time=.true.
    !-----------------------------------------------------------------------------
    call self%composition%element_by_id(trim(self%atom%id), element)
    if (.not. associated(element)) then
      call self%composition%molecule_by_id(trim(self%atom%id), molec)
      if (.not. associated(molec)) then
        return
      else
        if (molec%nline == 0) return
      endif
    else
      if (all(element%ions%nline == 0)) return
    endif

    if (exclude_from_line_list) then
      if (associated(element)) then
        element%ions%in_linelist = .false.
      else
        molec%in_linelist = .false.
      endif
      if (first_time) then
        first_time = .false.
        print *, 'excluding all ', trim(self%atom%id), ' lines from line list'
      endif
    else
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
    endif
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
      class is (atom_line_t)
        needed = .false.
        if (associated(another%atom)) then
          if (another%atom%id /= self%atom%id) then
            needs_me = .false.
          endif
        else
          needs_me = .false.
        endif
      class is (atom_cont_patch_t)
        needed = .false.
        if (another%atom%id /= self%atom%id) then
          needs_me = .false.
        endif
      class is (atom_rate_t)
        if (another%atom%id /= self%atom%id) then
          needed = .false.; needs_me = .false.
        endif
      class is (net_opac_t)
        needed = .false.
      class default
        if (another%id == self%m3d_parent%id) then
          ! the composition patch needs this to be able to
          ! update again before final spectrum is printed
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

  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(atom_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = .false.
    
    if (self%m3d_patch_t%check_in_bounds()) then
      cib = self%lvlcheck > 0 .and. self%lvlcheck <= self%atom%nlevel
    elseif (self%verbose > 2) then
      cib = .true.
    endif
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE save_abundance_to_file (self)
    class(atom_patch_t), target :: self
    integer :: lunit, i
    !-----------------------------------------------------------------------------
    if (self%ipos(3)==0) then
      open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED', access='append')
      write(lunit, '(A3, F8.4)') uppercase(self%atom%id, len(self%atom%id)), self%atom%abnd
      close(lunit)
    endif
  END SUBROUTINE save_abundance_to_file

  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(atom_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    integer, save:: itimer=0
    integer :: i, i_nu, i_mu, key(4)
    !-----------------------------------------------------------------------------
    call trace%begin('atom_patch_t%init_nbors_hash', itimer=itimer)
    call self%m3d_patch_t%init_nbors_hash(nbors)
  
    key(1) = self%key(1)
    do i_nu=1, self%n_nu
      key(2) = i_nu
      do i_mu=0, self%n_omega
        key(3) = i_mu
        if (i_mu == 0) then
          key(4) = hash_string('atom_cont')
          call self%check_nbor(key, nbors)
        endif
        if (comoving_frame .or. i_mu > 0) then
          key(4) = hash_string('atom_line')
          call self%check_nbor(key, nbors)

          key(4) = hash_string('net_opac')
          call self%check_nbor(key, nbors)

          key(4) = hash_string('atom_rate')
          call self%check_nbor(key, nbors)
        endif
      end do
    end do
    do i_nu=self%n_nu+1, 2*self%n_nu
      key(2) = i_nu
      do i=0, self%n_out
        if (i == 0) then
          key(3) = 0
          key(4) = hash_string('atom_cont')
          call self%check_nbor(key, nbors)
        else
          key(3) = i + self%n_omega
          key(4) = hash_string('atom_line')
          call self%check_nbor(key, nbors)

          key(4) = hash_string('net_opac')
          call self%check_nbor(key, nbors)
        endif
      end do
    end do
    call trace%end (itimer)
  END SUBROUTINE init_nbors_hash

  END MODULE
