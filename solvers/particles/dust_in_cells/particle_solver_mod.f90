!===============================================================================
!> Particle solver for handling dust particles represented by array data,
!> allowing vectorization of the KDK part and using only 36 bytes per particle.
!> Updating the particles is done by looping over cells, and looping over bins
!> and particles in bins within each cells, with these steps:
!>
!> 1) compute cell gas properties needed for the update
!> 2) update particle positions with a kick-drift-kick step
!> 3) trade particles with neighboring cells
!> 4) import particles from neighboring patches
!> 5) remove holes in the arrays caused by the trading
!> 6) sumup the particle weights to get dust densities
!> 7) add or remove particles and renormalize weights
!>
!> Step 3 is the most expensive one, using about 35% of the time
!===============================================================================
MODULE particle_solver_mod
  USE omp_lib
  USE io_mod
  USE trace_mod
  USE random_mod
  USE task_mod
  USE gpatch_mod
  USE scaling_mod
  USE units_mod
  USE list_mod
  USE link_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE particles_mod
  USE mesh_mod
  USE timer_mod
  implicit none
  private
  type, public, extends(task_t):: particle_solver_t
    integer                    :: n=32, n_bin=5, mp=60, me, np, lb(3), ub(3)
    real, pointer, contiguous  :: s_cgs(:), m_cgs(:)
    real, pointer, contiguous  :: a(:,:,:,:), v(:,:,:,:)
    real, pointer, contiguous  :: d(:,:,:), dd(:,:,:,:)
    type(particles_t), pointer :: particles(:,:,:,:)
    type(particles_t), pointer :: previous(:,:,:,:)
    type(particles_t), pointer :: exports(:) => null()
    class(gpatch_t), pointer   :: patch
    class(list_t), pointer     :: task_list
    class(nbor_t), pointer     :: nbor
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: initial_condition
    procedure:: update
    procedure, private:: setup
    procedure:: pack
    procedure, private:: import
    procedure, private:: kdk
    procedure, private:: trade
    procedure, private:: sumup
    procedure, private:: renorm
    procedure:: output
  end type
  logical:: on = .true.
  logical:: particle_io = .false.
  integer:: verbose = 0
  integer:: n_bin = 1
  integer:: particles_per_bin = 30
  integer:: output_every = 1000
  integer:: seed = 11                         ! random seed
  real:: dust_to_gas = 0.01                   ! total dust-to-gas ratio
  real:: min_fraction = 0.7                   ! particles per bin (max decrease)
  real:: max_fraction = 1.4                   ! particles per bin (max increase)
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: solid_density_cgs = 2.0              ! g/cm**3
  real:: density_yr_cgs = 3e-12               ! volume density giving 1 yr stop
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  self%patch => patch
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize dust parameters
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t):: self
  integer:: iostat
  namelist /particle_solver_params/ on, verbose, particles_per_bin, n_bin, seed, &
    dust_to_gas, min_size_cgs, max_size_cgs, min_fraction, max_fraction, &
    solid_density_cgs, density_yr_cgs, particle_io
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  !.............................................................................
  call io%lock%set ('particle_solver_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
    first_time = .false.
  end if
  call io%lock%unset ('particle_solver_t%init')
  !.............................................................................
  self%n_bin = n_bin
  self%mp = 2*max_fraction*particles_per_bin
  self%me = particles_per_bin*self%patch%n(1)*self%patch%n(2)
  call self%random%init (seed)
  call self%aux%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize dust bins
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t):: self
  integer:: mp, np, nmax, iostat
  integer:: i1, i2, i3, ib, j, m(3)
  real:: w
  type(particles_t), pointer:: p, o
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('particle_solver_t%initial_condition')
  m = self%patch%gn
  allocate (self%exports(n_bin))
  allocate (self%particles(m(1),m(2),m(3),n_bin))
  allocate (self%previous (m(1),m(2),m(3),n_bin))
  allocate (self%v(      3,m(1),m(2),m(3)))
  allocate (self%a(      3,m(1),m(2),m(3)))
  allocate (self%dd(       m(1),m(2),m(3),n_bin))
  allocate (self%s_cgs(n_bin), self%m_cgs(n_bin))
  call io%bits_mem (storage_size(self%exports), n_bin, 'mem')
  call io%bits_mem (storage_size(self%particles), &
    product(shape(self%particles)), 'mem')
  call io%bits_mem (storage_size(self%previous), &
    product(shape(self%previous)), 'mem')
  call io%bits_mem (storage_size(self%v ), product(shape(self%v )), 'v')
  call io%bits_mem (storage_size(self%a ), product(shape(self%a )), 'a')
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  self%d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  do ib=1,n_bin
    self%s_cgs(ib) = exp(log(min_size_cgs) &
                     + (ib-0.5)*log(max_size_cgs/min_size_cgs)/n_bin)
    self%m_cgs(ib) = solid_density_cgs*self%s_cgs(ib)**3
    write(io_unit%log,'(a,i4,1p,2e11.2)') &
      'dust_t%init: bin, s_cgs, m_cgs =', ib, self%s_cgs(ib), self%m_cgs(ib)
    allocate (self%exports(ib)%lock)
    call self%exports(ib)%lock%init ('exp')
    p => self%exports(ib)
    call p%lock%set ('inii')
    call p%init (self%me, verbose)
    call p%lock%unset ('inii')
    do i3=1,m(3)
    do i2=1,m(2)
    do i1=1,m(1)
      self%dd(i1,i2,i3,ib) = dust_to_gas*self%d(i1,i2,i3)
      self%a(:,i1,i2,i3) = 0.0              ! cell acceleration
      self%v(:,i1,i2,i3) = 0.0              ! cell velocity
      p => self%particles(i1,i2,i3,ib)
      call p%init (self%mp, verbose)
      if (verbose > 0) then
        o => self%previous(i1,i2,i3,ib)
        call o%init (self%mp, verbose)
      end if
      p%id = [self%patch%id,i1,i2,i3,ib]
      p%np = particles_per_bin
      p%lp = p%np
      w = self%dd(i1,i2,i3,ib)/(n_bin*p%np)
      do j=1,p%np
        p%r(j,:) = self%random%uniform(3)
        p%v(j,:) = 0.0
        p%q(j,:) = [i1,i2,i3]
        p%w(j) = w
        p%t(j) = 0d0
        p%flag(j) = 0
      end do
    end do
    end do
    end do
  end do
  call trace%end ()
END SUBROUTINE initial_condition

!===============================================================================
!> Initialize particle solver, after nbor lists are created, making it possible
!> to create a copy of the patch nbor list, with tasks pointing to self
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(particle_solver_t):: self
  class(list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('particle_solver_t%init_task_list')
  self%task_list => task_list
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update particles, using a series of operations.  Constraints: setup must be 
!> done before kdk, trade must be done after kdk, import must be done after 
!> trade and before pack, sumup is cheapest after pack, and renorm needs to be 
!> done after sumup.
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  logical:: dbg
  real(8):: wc
  class(particles_t), pointer:: p
  if (on) then
    call trace%begin ('particle_solver_t%update')
    wc = wallclock()
    call self%output
    call self%setup               ! set up cell properties
    call self%kdk                 ! move particles
    call self%trade               ! trade btw cells and export
    call self%import              ! import particles from nbors
    call self%pack                ! remove array holes
    call self%sumup               ! compute densities
    call self%renorm              ! renormalize particles
    associate (t_particles => timer%t_particles)
    !$omp atomic update
    t_particles = t_particles + (wallclock()-wc)
    end associate
    call trace%end ()
  end if
contains
subroutine ppp
  integer:: j
  if (dbg) then
    do j=1,p%lp
      print 3,j,p%r(j,:),p%w(j)
      3 format(i4,3f7.3,1p,e10.2)
    end do
  end if
end subroutine
END SUBROUTINE update

!===============================================================================
!> Initialize cell data
!===============================================================================
SUBROUTINE setup (self)
  class(particle_solver_t)                      :: self
  !.............................................................................
  integer                                       :: i1, i2, i3, ib, np
  real                                          :: c1, c2, c3
  real(8)                                       :: wc
  class(gpatch_t), pointer                      :: patch
  type(particles_t), pointer                    :: p
  real, pointer, contiguous, dimension(:,:,:)   :: p1, p2, p3
  integer, save                                 :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%setup', itimer=itimer)
  wc = omp_get_wtime()
  patch => self%patch
  self%dtime = patch%dtime
  self%time = patch%time + patch%dtime
  p1 => patch%mem(:,:,:,patch%idx%px,patch%it,1)
  p2 => patch%mem(:,:,:,patch%idx%py,patch%it,1)
  p3 => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
  self%lb = patch%mesh%li
  self%ub = patch%mesh%ui
  np = 0
  c1 = 1.0/self%patch%ds(1)
  c2 = 1.0/self%patch%ds(2)
  c3 = 1.0/self%patch%ds(3)
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    self%v(1,i1,i2,i3) = p1(i1,i2,i3)/self%d(i1,i2,i3)*c1
    self%v(2,i1,i2,i3) = p2(i1,i2,i3)/self%d(i1,i2,i3)*c2
    self%v(3,i1,i2,i3) = p3(i1,i2,i3)/self%d(i1,i2,i3)*c3
    self%a(1,i1,i2,i3) = 0.0
    self%a(2,i1,i2,i3) = 0.0
    self%a(3,i1,i2,i3) = 0.0
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      np = np + p%np
    end do
  end do
  end do
  end do
  self%np = np
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,i6)') &
      wc, wc/np, self%patch%id, np, 'setup', self%patch%id
  end if
  call trace%end (itimer)
END SUBROUTINE setup

!===============================================================================
!> Kick-drift-kick update, with no hole in any cells
!===============================================================================
SUBROUTINE kdk (self)
  class(particle_solver_t)   :: self
  integer                    :: i1, i2, i3, ib, i, j, q, np
  real                       :: dt, drag_c, cs, c, a, v, vmin, vmax, cmin, cmax
  real(8)                    :: wc
  type(particles_t), pointer :: p, o
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%kdk', itimer=itimer)
  wc = omp_get_wtime()
  !-----------------------------------------------------------------------------
  dt = self%dtime*0.5d0
  drag_c = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  if (verbose==2) then
    write (io_unit%log,'("kdk:",i6)') self%patch%id
    flush (io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  vmin = 0.0
  vmax = 0.0
  cmin = huge(1.)
  cmax = 0.0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    !if (verbose==3) then
    !  write (io_unit%log,'("kdk:",i6,3i4)') self%patch%id, i1, i2, i3
    !end if
    cs = drag_c*self%d(i1,i2,i3)
    do ib=1,n_bin
      if (verbose==4) then
        write (io_unit%log,'("kdk:",i6,4i4)') self%patch%id, i1, i2, i3, ib
      end if
      c = cs/self%s_cgs(ib)
      !if (self%patch%id==1.and.i3==self%lb(3).and.i2==self%lb(2).and.i1==self%lb(1)) then
      !  write (stdout,*) 'ib,c:', ib,c
      !end if
      p => self%particles(i1,i2,i3,ib)
      associate (n_particles => timer%n_particles)
      !$omp atomic update
      n_particles = n_particles + p%np
      end associate
      if (verbose > 0) then
        o => self%previous(i1,i2,i3,ib)
        o%lp = p%lp
        o%np = p%np
        do j=1,p%lp
          dt = self%dtime*0.5d0
          cmin = min(c*dt,cmin)
          cmax = max(c*dt,cmax)
          o%t(j) = p%t(j)
          o%w(j) = p%w(j)
          do i=1,3
            o%r(j,i) = p%r(j,i)
            o%v(j,i) = p%v(j,i)
            o%q(j,i) = p%q(j,i)
          end do
        end do
      end if
      np = p%np
      if (verbose==5) then
        write (io_unit%log,'("kdk:",i6,4i4,1p,2(3e10.2,2x))') &
          self%patch%id, i1, i2, i3, ib, self%v(:,i1,i2,i3), self%a(:,i1,i2,i3)
      end if
      do i=1,3
        a = self%a(i,i1,i2,i3)
        v = self%v(i,i1,i2,i3)
        if (verbose==6) then
          do j=1,np
            dt = max(0.0,0.5*(self%time - p%t(j)))
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            p%r(j,i) = p%r(j,i) + dt*p%v(j,i)*2.0
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            q = floor(p%r(j,i))
            p%r(j,i) = p%r(j,i) - q
            p%q(j,i) = p%q(j,i) + q
!            vmin = min(vmin,p%v(j,i))
!            vmax = max(vmax,p%v(j,i))
            write (io_unit%log,'("kdk:",i6,8i4,1p,3e10.2)') self%patch%id, &
              i1, i2, i3, ib, j, i, q, p%q(j,i), p%r(j,i), p%v(j,i)
          end do
        else
          !$omp simd
          do j=1,np
            dt = max(0.0,0.5*(self%time - p%t(j)))
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            p%r(j,i) = p%r(j,i) + dt*p%v(j,i)*2.0
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            q = floor(p%r(j,i))
            p%r(j,i) = p%r(j,i) - q
            p%q(j,i) = p%q(j,i) + q
!            vmin = min(vmin,p%v(j,i))
!            vmax = max(vmax,p%v(j,i))
          end do
        end if
      end do
      do j=1,np
        p%t(j) = self%time
      end do
    end do
  end do
  end do
  end do
  !................................................-+.............................
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,1p,4e9.1)') &
      wc, wc/self%np, self%patch%id, self%np, 'KDK', &
      vmin, vmax, cmin, cmax
  end if
  call trace%end (itimer)
END SUBROUTINE kdk

!===============================================================================
!> Re-assign particles to other cells.  New particles will be added at the
!> end, and holes will only appear for indices up to the initial array lenght
!===============================================================================
SUBROUTINE trade (self)
  class(particle_solver_t)                    :: self
  integer                                     :: i1, i2, i3, ip, ib, i, j, ii(3)
  integer                                     :: q1, q2, q3, pq(3)
  integer                                     :: j2, ne, nm, ne1, nm1
  integer                                     :: n1, n2, n3, n4, n5, me, mp
  real(8)                                     :: wc
  real                                        :: dt
  type(particles_t), pointer                  :: p, pn, pe, o
  real, pointer, contiguous, dimension(:,:)   :: pv, pr
  real, pointer, contiguous, dimension(:,:,:) :: p1, p2, p3
  logical                                     :: export, move
  logical, save                               :: warn(2) = .true.
  integer, save                               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%trade', itimer=itimer)
  wc = omp_get_wtime()
  n1 = 0
  n2 = 0
  n3 = 0
  n4 = 0
  n5 = 0
  nm = 0
  ne = 0
  do ib=1,n_bin
    pe => self%exports(ib)
    call pe%lock%set ('exp')
    if (verbose > 1) then
      write (io_unit%log,'(i6,2x,a,3i10)') &
      self%patch%id, 'exports before:', pe%np, pe%lp, pe%mp
    end if
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      ii = [i1,i2,i3]
      p => self%particles(i1,i2,i3,ib)
      if (verbose > 0) then
        o => self%previous(i1,i2,i3,ib)
      end if
      if (verbose > 1) then
        call  p%validate ('p  before')
        call pe%validate ('pe before')
      end if
      if (p%np == 0) write(io_unit%log,'(a,6i5)') '0p bef trade', p%id, p%lp
      ne1 = 0
      nm1 = 0
      if (verbose > 2) then
        write (io_unit%log,'(a,i6,4i4,6i6)') 'trade: id,i1,i2,i3,ib,nm,ne,nlp,nle=', &
        self%patch%id,i1,i2,i3,ib,nm1,ne1,p%np,p%lp,pe%np,pe%lp
      end if
      do j=1,p%lp
        if (p%w(j) > 0.0) then
          !---------------------------------------------------------------------
          ! If the position index is outside, move the particle to exports
          !---------------------------------------------------------------------
          pq(1) = p%q(j,1)
          pq(2) = p%q(j,2)
          pq(3) = p%q(j,3)
          export = .false.
          move = .false.
          do i=1,3
            if (pq(i) /= ii(i)) then
              if (pq(i) < self%lb(i) .or. pq(i) > self%ub(i)) then
                export = .true.
                exit
              end if
              move = .true.
            end if
          end do
          if (export) then
            ne1 = ne1 + 1
            if (pe%lp >= pe%mp) then
              me = 1.2*pe%mp
              call pe%resize(me)
              self%me = me
              if (verbose > 0) then
                write (io_unit%log,*) &
                  'export resize: mo, mp =', me, self%me
              end if
              n3 = n3+1
            end if
            pe%np = pe%np + 1           ! increment export count
            pe%lp = pe%lp + 1           ! increment export length
            j2 = pe%lp
            do i=1,3
              pe%r(j2,i) = p%r(j,i)
              pe%v(j2,i) = p%v(j,i)
              pe%q(j2,i) = p%q(j,i)
            end do
            pe%w(j2) = p%w(j)
            pe%t(j2) = p%t(j)
            p%w(j) = 0.0                ! mark position empty
            p%np = p%np - 1             ! decrement actual count
            if (j == p%lp) then         ! if last slot
              p%lp = p%lp - 1           ! also decrement length
            end if
            if (verbose > 1) then
              call  p%validate ('p0 after ')
              call pe%validate ('pe after ')
              if (verbose > 2) then
                write (io_unit%log,'(i6,":",3i4,2i6,2x,a,i6,3i4,2i6,1p,e10.2)') &
                j, i1, i2, i3, p%np, p%lp, 'export', j2, pq, pe%np, pe%lp, pe%w(j2)
              end if
            end if
            !-------------------------------------------------------------------
            ! If the position index has changed, add to the target cell and
            ! remove from the source cell
            !-------------------------------------------------------------------
          else if (move) then
            nm1 = nm1 + 1
            !q1 = pq(1)
            !q2 = pq(2)
            !q3 = pq(3)
            pn => self%particles(pq(1),pq(2),pq(3),ib)
            if (pn%lp >= p%mp) then
              mp = 1.2*pn%mp
              call pn%resize(mp)
              pn%mp = mp
              if (verbose > 0) then
                write (io_unit%log,*) &
                  'particles resize: mo, mp =', mp, self%mp
              end if
              n4 = n4+1
            end if
            if (pn%lp < pn%mp) then
              pn%np = pn%np + 1           ! increment target count
              j2 = pn%lp + 1
              pn%lp = j2                  ! increment target length
              do i=1,3
                pn%r(j2,i) = p%r(j,i)
                pn%v(j2,i) = p%v(j,i)
                pn%q(j2,i) = p%q(j,i)
              end do
              pn%flag(j2) = p%flag(j)
              pn%t(j2) = p%t(j)
              pn%w(j2) = p%w(j)
              p%w(j) = 0.0                ! mark position empty
              p%np = p%np - 1             ! decrement source count
              !if (j == p%lp) then         ! if last slot
              !  p%lp = p%lp - 1           ! also decrement length
              !end if
            else
              call io%abort ('should not happen')
            end if
            if (verbose > 1) then
              call  p%validate ('p1 after ')
              call pn%validate ('pn after ')
            end if
          end if
        end if
      end do
      if (p%np <= 0.5*min_fraction*particles_per_bin) then
        n5 = n5+1
      end if
      ne = ne + ne1
      nm = nm + nm1
      if (p%np == 0) then
        if (verbose > 0) then
          write(io_unit%log,'(a,i6,4i4,3f6.1)') &
            '0p after trade', p%id, self%v(:,i1,i2,i3)*self%dtime
          do j=1,o%np
            dt = p%t(j) - o%t(j)
            write (io_unit%log,'(a,1p,3e10.2,0p,2(3i3,1x,3f6.2,1x,3f6.2))') &
              '0p', dt, o%w(j), p%w(j), &
              o%q(j,:), o%r(j,:), o%v(j,:)*dt, &
              p%q(j,:), p%r(j,:), p%v(j,:)*dt
          end do
        end if
      end if
      if (verbose > 2) then
        write (io_unit%log,'(a,i6,4i4,6i6)') 'trade: id,i1,i2,i3,ib,nm,ne,nlp.nle=', &
        self%patch%id,i1,i2,i3,ib,nm1,ne1,p%np,p%lp,pe%np,pe%lp
      end if
      if (verbose > 1) then
        call pe%validate ('pe')
      end if
    end do
    end do
    end do
    if (verbose > 1) then
      write (io_unit%log,'(i6,2x,a,3i10)') &
      self%patch%id, 'exports  after:', pe%np, pe%lp, pe%mp
    end if
    call pe%lock%unset ('exp')
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,7i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'trade', nm, ne, n1, n2, n3, n4, n5
  end if
  call trace%end (itimer)
END SUBROUTINE trade

!===============================================================================
!> Import particles from other patch's exports
!===============================================================================
SUBROUTINE import (self)
  class(particle_solver_t)   :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particles_t), pointer:: pe, p
  real(8)                    :: wc, r(3), dist(3)
  integer                    :: ni, ne, ib, n, j, jp, i, q(3), np
  integer                    :: i1, i2, i3
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%import' ,itimer=itimer)
  if (verbose > 1) then
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      do ib=1,n_bin
        p => self%particles(i1,i2,i3,ib)
        write(io_unit%log,'(a,i6,4i4,3i6)') &
        'after import: i1,i2,i3,ib,np,lp,mp0', &
        self%patch%id,i1,i2,i3,ib,p%np,p%lp,p%mp
        call pe%validate ('before import')
      end do
    end do
    end do
    end do
  end if
  ni = 0
  ne = 0
  wc = omp_get_wtime()
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      dist = nbtask%patch%task_t%distance(self%patch)/nbtask%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      n = 0
      do ib=1,n_bin
        pe => nbtask%exports(ib)
        call pe%lock%set ('import')
        j = 1
        do while (j < pe%lp)
          if (pe%w(j) > 0.0) then
            r = dist + (pe%r(j,:)+pe%q(j,:))
            q = floor(r)
            if (all(q >= self%lb .and. q <= self%ub)) then
              p => self%particles(q(1),q(2),q(3),ib)
              if (p%lp >= p%mp) then
                call p%pack (np)
                if (verbose > 1) call p%validate ('p forced')
              end if
              if (p%lp < p%mp) then
                n = n  + 1
                jp = p%lp + 1
                do i=1,3
                  p%r(jp,i) = pe%r(j,i)
                  p%v(jp,i) = pe%v(j,i)
                  p%q(jp,i) = q(i)
                end do
                p%w(jp) = pe%w(j)
                p%t(jp) = pe%t(j)
                p%lp = jp
                p%np = p%np+1
                pe%w(j) = 0.0
                pe%np = pe%np - 1
              else
                ne = ne + 1
              end if
            end if
          end if
          j = j + 1
        end do
        call pe%lock%unset ('import')
      end do
      ni = ni + n
    end select
    nbor => nbor%next
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 1) then
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      do ib=1,n_bin
        p => self%particles(i1,i2,i3,ib)
        write(io_unit%log,'(a,i6,4i4,3i6)') &
        'after import: i1,i2,i3,ib,np,lp,mp0', &
        self%patch%id,i1,i2,i3,ib,p%np,p%lp,p%mp
        call p%validate ('after import')
      end do
    end do
    end do
    end do
  end if
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,2i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'import', ni, ne
  end if
  call trace%end (itimer)
END SUBROUTINE import

!===============================================================================
!> Remove the holes in the arrays
!===============================================================================
SUBROUTINE pack (self)
  class(particle_solver_t)   :: self
  integer                    :: i1, i2, i3, ib, n, lp(2), np(2), mp(2)
  real(8)                    :: wc
  type(particles_t), pointer :: p
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pack', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Cell arrays
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  n = 0
  lp = 0
  np = 0
  !-----------------------------------------------------------------------------
  ! Pack particles
  !-----------------------------------------------------------------------------
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      p%id = [self%patch%id,i1,i2,i3,ib]
      np(1) = max(np(1),p%np)
      lp(1) = max(lp(1),p%lp)
      if (p%np == 0) write(io_unit%log,'(a,7i5,3f6.1)') &
        '0p bef pack ',p%id,p%np,p%lp, self%v(:,i1,i2,i3)
      call p%pack (n)
      !if (p%np == 0) write(io_unit%log,'(a,7i5,3f6.1)') &
        !'0p aft pack ',p%id,p%np,p%lp, self%v(:,i1,i2,i3)
      if (verbose > 1) call p%validate ('normal p')
    end do
  end do
  end do
  end do
  mp(1) = p%mp
  !-----------------------------------------------------------------------------
  ! Pack exports
  !-----------------------------------------------------------------------------
  do ib=1,n_bin
    p => self%exports(ib)
    call p%lock%set ('pack')
    p%id = [self%patch%id,0,0,0,ib]
    np(2) = max(np(2),p%np)
    lp(2) = max(lp(2),p%lp)
    call p%pack (n)
    call p%lock%unset ('pack')
    if (verbose > 1) call p%validate ('normal pe')
  end do
  mp(2) = p%mp
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,7i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'pack', n, np, lp, mp
    flush (io_unit%log)
  end if
  call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Sum up density
!===============================================================================
SUBROUTINE sumup (self)
  class(particle_solver_t)    :: self
  integer                     :: i1, i2, i3, ib, j, nmax
  real                        :: dd
  real(8)                     :: wc
  type(particles_t), pointer  :: p
  integer, save               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%sumup', itimer=itimer)
  wc = omp_get_wtime()
  nmax = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      nmax = max(nmax,p%lp)
      dd = 0.0
      do j=1,p%lp
        dd = dd + p%w(j)
      end do
      self%dd(i1,i2,i3,ib) = dd
      if (dd == 0.0) then
        write (io_unit%log,'(a,i6,4i4,2i6,1p,4e10.2)') &
        'sumup 0d: id,np,lp,v=', &
        p%id, p%np, p%lp, self%v(:,i1,i2,i3)*self%dtime
      end if
    end do
  end do
  end do
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'density', nmax
  end if
  call trace%end (itimer)
END SUBROUTINE sumup

!===============================================================================
!> Add and remove particles where needed
!===============================================================================
  SUBROUTINE renorm (self)
  class(particle_solver_t)    :: self
  integer                     :: i1, i2, i3, ib, i, j, n, na, nr, np, lp
  real                        :: n_min, n_max, w
  real(8)                     :: t, wc
  type(particles_t), pointer  :: p
  integer, save               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%renorm', itimer=itimer)
  wc = omp_get_wtime()
  n_min = particles_per_bin*min_fraction
  n_max = particles_per_bin*max_fraction
  na = 0
  nr = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      9 format(3f8.3,i4)
      call io%assert (p%np==p%lp, 'The state must be packed before renorm')
      if (self%dd(i1,i2,i3,ib) <= 0.0) then
        write (io_unit%log,'(a,i6,4i3,2i4,1p,5e9.1,i10)') &
        'zero density warning: id,np,lp,v=', &
        p%id,p%np,p%lp, self%v(:,i1,i2,i3), self%d(i1,i2,i3), &
        self%patch%time, self%patch%istep
        cycle
      end if
      !-------------------------------------------------------------------------
      ! If too few particles, add new ones and set new weights to all
      !-------------------------------------------------------------------------
      if (p%np < n_min) then
        w = self%dd(i1,i2,i3,ib)/particles_per_bin
        !call p%pack (n)
        t = 0d0
        do j=1,p%np
          p%w(j) = w
          t = t + p%t(j)
        end do
        t = t/p%np
        do j=p%np+1,particles_per_bin
          do i=1,3
            p%r(j,i) = self%random%ran1()
            p%v(j,i) = self%v(i,i1,i2,i3)
          end do
          p%q(j,:) = [i1,i2,i3]
          p%w(j) = w
          p%t(j) = t
        end do
        na = na + particles_per_bin - p%np
        p%np = particles_per_bin
        p%lp = particles_per_bin
        if (verbose > 1) call p%validate ('renorm add')
      !-------------------------------------------------------------------------
      ! If too many, renormalize and remove extraneous
      !-------------------------------------------------------------------------
      else if (p%np > n_max) then
        nr = nr + (p%np - particles_per_bin)
        j = p%lp
        np = p%np
        lp = p%lp
        do while (p%np > particles_per_bin)
          if (j <= 0) then
            write (io_unit%log,'(a,i6,4i4,2x,4i4)') &
              'error:', self%patch%id, i1, i2, i3, ib, np, lp, p%np, p%lp
            write (stderr,*) 'ERROR'
            stop
          end if
          if (verbose > 2) then
            write (io_unit%log,'("rem:",3i4,1p,e10.2 )') j, p%np, p%lp, p%w(j)
          end if
          if (p%w(j) == 0.0) then                   ! if a hole
            p%lp = p%lp - 1                         ! shorten length
          else                                      ! if a particle
            p%np = p%np - 1                         ! decrement count
            p%lp = p%lp - 1                         ! shorten length
          end if
          j = j-1
        end do
        w = self%dd(i1,i2,i3,ib)/particles_per_bin  ! new weight
        !call p%pack (n)                             ! pack and set weight
        do j=1,p%np                                 ! assumes packed state!
          p%w(j) = w
        end do
        p%np = particles_per_bin
        p%lp = particles_per_bin
        if (verbose > 1) call p%validate ('renorm rem')
      end if
    end do
  end do
  end do
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,2i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'renorm', na, nr
  end if
  call trace%end (itimer)
END SUBROUTINE renorm

!===============================================================================
!> I/O
!===============================================================================
SUBROUTINE output (self)
  class(particle_solver_t):: self
  !.............................................................................
  class(particles_t), pointer:: p
  class(mesh_t), pointer:: mesh(:)
  real, pointer, contiguous:: tmp_r(:,:)
  integer, pointer, contiguous:: tmp_i(:)
  integer:: i, j, n, i1, i2, i3, ib
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%output')
  if (verbose > 1) then
    write (io_unit%log, '(a,i6,2f8.4)') &
      'dust_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_d', self%dd)
    !...........................................................................
    if (verbose > 1) then
      do i=1,n_bin
        write (io_unit%log, '(a,i6, i3,2x,1p,10e11.2)') &
          'dust_t%output: id, bin, n =', self%patch%id, i, self%dd(3,3,1:10,i)
      end do
    end if
    !...........................................................................
    if (particle_io) then
      mesh => self%patch%mesh
      do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(1)
      do i1=self%lb(1),self%ub(1)
        do ib=1,n_bin
          n = 0
          p => self%particles(i1,i2,i3,ib)
          do j=1,p%np
            if (p%flag(j) > 0) then
             n = n+1
            end if
          end do
          allocate (tmp_r(5,n))
          i = 0
          do j=1,p%np
            if (p%flag(j) > 0) then
              i = i+1
              tmp_r(1:3,i) = mesh%p + (p%q(j,:) + p%r(j,:) - mesh%o)*mesh%d
              tmp_r(  4,i) = p%w(j)
              tmp_r(  5,i) = self%s_cgs(ib)
            end if
          end do
        end do
      end do
      end do
      end do
      call self%aux%register ('dust_p' , tmp_r)
      call self%aux%output_sequence (self%patch%id, self%time)
      if (verbose > 0) then
        write (io_unit%log, '(a,i6,f8.4,i3,2i8)') &
         'dust_t%output: id, particle array =', &
          self%patch%id, self%time, shape(tmp_r), i
      end if
      deallocate (tmp_r)
    else
      call self%aux%output_sequence (self%patch%id, self%time)
    end if
    if (verbose > 1) then
      write (io_unit%log, '(a,i6,1p,10e11.2)') &
       'dust_t%output: self%aux%aux_next =', &
        self%patch%id, self%aux%aux_next
    end if
  end if
  call trace%end ()
END SUBROUTINE output

END MODULE particle_solver_mod
