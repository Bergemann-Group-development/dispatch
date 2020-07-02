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
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_lock_mod
  implicit none
  private
  type, public:: particles_t
    real           , pointer, contiguous:: r(:,:), v(:,:), w(:)
    integer(kind=2), pointer, contiguous:: q(:,:)
    type(lock_t), pointer               :: lock => null()
    integer                             :: id(5)
    integer                             :: np = 0         ! actual count
    integer                             :: lp = 0         ! array length
    integer                             :: mp = 0         ! max count
  contains
    procedure:: init
    procedure:: pack => particle_pack
    procedure:: validate
    procedure:: resize
  end type
CONTAINS

!===============================================================================
!> Initialize particle arrays
!===============================================================================
SUBROUTINE init (self, mp)
  class(particles_t):: self
  integer:: mp
  !-----------------------------------------------------------------------------
  allocate (self%r(mp,3), self%v(mp,3), self%q(mp,3), self%w(mp))
  call io%bits_mem (storage_size(self%w), mp*(1+3*3), 'p')
  self%mp = mp
  self%np = 0
  self%lp = 0
END SUBROUTINE init

!===============================================================================
!> Resize
!===============================================================================
SUBROUTINE resize (self, mp)
  class(particles_t):: self
  integer:: mp, mo
  real, pointer:: r2(:,:), r1(:)
  integer(kind=2), pointer:: i2(:,:)
  !-----------------------------------------------------------------------------
  mo = self%mp
  call io%bits_mem (storage_size(self%w), -mo*(1+3*3), 'p')
  call io%bits_mem (storage_size(self%w), +mp*(1+3*3), 'p')
  r2 => self%r; allocate (self%r(mp,3)); self%r(1:mo,:) = r2(:,:)
  r2 => self%v; allocate (self%v(mp,3)); self%v(1:mo,:) = r2(:,:)
  i2 => self%q; allocate (self%q(mp,3)); self%q(1:mo,:) = i2(:,:)
  r1 => self%w; allocate (self%w(mp  )); self%w(1:mo  ) = r1(:)
  self%mp = mp
END SUBROUTINE resize

!===============================================================================
!> Remove holes in particle arrays
!===============================================================================
SUBROUTINE particle_pack (self, n)
  class(particles_t):: self
  integer:: n
  integer:: j1, j2, i
  !-----------------------------------------------------------------------------
  if (self%lp > self%np) then
    !---------------------------------------------------------------------------
    ! Search for additional holes, if any
    !---------------------------------------------------------------------------
    j2 = self%lp
    j1 = 1
    n = self%lp-self%np 
    do while (self%lp > self%np)
      if (self%w(j1) == 0.0)   then
        !-----------------------------------------------------------------------
        ! Shrink from the end, looking for a valid entry beyond j1
        !-----------------------------------------------------------------------
        do while (self%w(j2) == 0.0 .and. j2 > 0)
          j2 = j2-1
        end do
        !-----------------------------------------------------------------------
        ! Unless j2 has shrunk below j1, copy the valid entry to the hole
        !-----------------------------------------------------------------------
        if (j2 > j1) then
          do i=1,3
            self%q(j1,i) = self%q(j2,i)
            self%r(j1,i) = self%r(j2,i)
            self%v(j1,i) = self%v(j2,i)
          end do
          self%w(j1) = self%w(j2)
          j2 = j2-1                                   ! previous end slot
          end if
        self%lp = j2
      end if
      j1 = j1+1
    end do
  end if
  if (self%np /= self%lp) then
    write (io_unit%log,'(a,i6,8i4)') &
    'pack ERROR: id,j1,j2,np,lp=', self%id, j1, j2, self%np, self%lp
    call io%abort ('particles_t%validate')
  end if
END SUBROUTINE particle_pack

!===============================================================================
!> Remove holes in particle arrays
!===============================================================================
SUBROUTINE validate (self, label)
  class(particles_t):: self
  character(len=*):: label
  integer:: j, n
  !-----------------------------------------------------------------------------
  n = self%lp - self%np
  do j=1,self%lp
    if (self%w(j) == 0.0) then
      n = n-1
      if (n < 0) then
        write (stderr,'(a,i4,i6,7i4)') &
          'validate ERROR: zero-w particle remains at '//label, &
          j, self%id, self%np, self%lp
        write (io_unit%log,'(a,i4,i6,7i4)') &
          'validate ERROR: zero-w particle remains at '//label, &
          j, self%id, self%np, self%lp
        call io%abort ('particles_t%validate')
      end if
    end if
  end do
END SUBROUTINE validate
END MODULE particles_mod

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
  USE particles_mod
  implicit none
  private
  type, public, extends(task_t):: particle_solver_t
    integer                    :: n=32, n_bin=5, mp=60, me, np, lb(3), ub(3)
    real, pointer, contiguous  :: s(:)
    real, pointer, contiguous  :: a(:,:,:,:), v(:,:,:,:)
    real, pointer, contiguous  :: d(:,:,:), dd(:,:,:,:)
    type(particles_t), pointer :: particles(:,:,:,:)
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
  end type
  logical:: on=.true.
  integer:: verbose = 0
  integer:: n_bin = 1
  integer:: particles_per_bin = 30
  integer:: seed = 11                         ! random seed
  real:: dust_to_gas = 0.01                   ! total dust-to-gas ratio
  real:: min_fraction = 0.7                   ! particles per bin (max decrease)
  real:: max_fraction = 1.4                   ! particles per bin (max increase)
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: solid_density_cgs = 2.0              ! g/cm**3
  real:: density_yr_cgs = 1e-12               ! volume density giving 1 yr stop
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
  namelist /dust_cell_params/ on, verbose, particles_per_bin, n_bin, seed, &
  dust_to_gas, min_size_cgs, max_size_cgs, min_fraction, max_fraction
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  !.............................................................................
  call io%lock%set ('particle_solver_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_cell_params, iostat=iostat)
    call io%namelist_warning ('dust_cell_params', iostat)
    write (io%namelist, dust_cell_params)
    first_time = .false.
  end if
  call io%lock%unset ('particle_solver_t%init')
  !.............................................................................
  self%n_bin = n_bin
  self%mp = 6*particles_per_bin
  self%me = 2*self%mp*self%patch%n(1)*self%patch%n(2)
  call self%random%init (seed)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize dust bins
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t):: self
  integer:: mp, np, nmax, iostat
  integer:: i1, i2, i3, ib, j, m(3)
  type(particles_t), pointer:: i, p
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('particle_solver_t%initial_condition')
  m = self%patch%gn
  allocate (self%exports(n_bin))
  allocate (self%particles(m(1),m(2),m(3),n_bin))
  allocate (self%v(      3,m(1),m(2),m(3)))
  allocate (self%a(      3,m(1),m(2),m(3)))
  allocate (self%dd(       m(1),m(2),m(3),n_bin))
  allocate (self%s(n_bin))
  call io%bits_mem (storage_size(self%exports), n_bin, 'mem')
  call io%bits_mem (storage_size(self%particles), &
    product(shape(self%particles)), 'mem')
  call io%bits_mem (storage_size(self%v ), product(shape(self%v )), 'v')
  call io%bits_mem (storage_size(self%a ), product(shape(self%a )), 'a')
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  self%d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%new,1)
  do ib=1,n_bin
    self%s(ib) = exp(log(min_size_cgs) &
                     + (ib-0.5)*log(max_size_cgs/min_size_cgs)/n_bin)
    allocate (self%exports(ib)%lock)
    p => self%exports(ib)
    call p%init (self%me)
    do i3=1,m(3)
    do i2=1,m(2)
    do i1=1,m(1)
      self%dd(i1,i2,i3,ib) = dust_to_gas*self%d(i1,i2,i3)
      self%a(:,i1,i2,i3) = 0.0              ! cell acceleration
      self%v(:,i1,i2,i3) = 0.0              ! cell velocity
      p => self%particles(i1,i2,i3,ib)
      call p%init (self%mp)
      p%id = [self%patch%id,i1,i2,i3,ib]
      p%np = particles_per_bin
      p%lp = p%np
      do j=1,p%np
        p%r(j,:) = self%random%uniform(3)
        p%v(j,:) = 0.0
        p%q(j,:) = [i1,i2,i3]
        p%w(j) = self%dd(i1,i2,i3,ib)/p%np
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
  class(particles_t), pointer:: p
  if (on) then
    call trace%begin ('particle_solver_t%update')
    call self%setup               ! set up cell properties
    call self%kdk                 ! move particles
    call self%trade               ! trade btw cells and export
    call self%import              ! import particles from nbors
    call self%pack                ! remove array holes
    call self%sumup               ! compute densities
    call self%renorm              ! renormalize particles
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
  p1 => patch%mem(:,:,:,patch%idx%px,patch%new,1)
  p2 => patch%mem(:,:,:,patch%idx%py,patch%new,1)
  p3 => patch%mem(:,:,:,patch%idx%pz,patch%new,1)
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
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10)') &
      wc, wc/np, self%patch%id, np, 'setup'
  end if
  call trace%end (itimer)
END SUBROUTINE setup

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
  integer                    :: i1, i2, i3, ib, n
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
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      p%id = [self%patch%id,i1,i2,i3,ib]
      if (p%np == 0) write(io_unit%log,'(a,7i5)') '0p bef pack',p%id,p%np,p%lp
      call p%pack (n)
      if (p%np == 0) write(io_unit%log,'(a,7i5)') '0p aft pack',p%id,p%np,p%lp
      if (verbose > 1) call p%validate ('normal p')
    end do
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Also in exports
  !-----------------------------------------------------------------------------
  do ib=1,n_bin
    p => self%exports(ib)
    p%id = [self%patch%id,0,0,0,ib]
    call p%pack (n)
    if (verbose > 1) call p%validate ('normal pe')
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,4i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'pack', n, p%np, p%lp, p%mp
    flush (io_unit%log)
  end if
  call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Kick-drift-kick update, with no hole in any cells
!===============================================================================
SUBROUTINE kdk (self)
  class(particle_solver_t)   :: self
  integer                    :: i1, i2, i3, ib, i, j, q, np
  real                       :: dt, drag_c, cs, c, a, v, vmin, vmax
  real(8)                    :: wc
  type(particles_t), pointer :: p
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
      c = cs/self%s(ib)
      p => self%particles(i1,i2,i3,ib)
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
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            p%r(j,i) = p%r(j,i) + dt*p%v(j,i)*2.0
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            q = floor(p%r(j,i))
            p%r(j,i) = p%r(j,i) - q
            p%q(j,i) = p%q(j,i) + q
            write (io_unit%log,'("kdk:",i6,8i4,1p,3e10.2)') self%patch%id, &
              i1, i2, i3, ib, j, i, q, p%q(j,i), p%r(j,i), p%v(j,i)
          end do
        else
          !$omp simd
          do j=1,np
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            p%r(j,i) = p%r(j,i) + dt*p%v(j,i)*2.0
            p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
            q = floor(p%r(j,i))
            p%r(j,i) = p%r(j,i) - q
            p%q(j,i) = p%q(j,i) + q
            vmin = min(vmin,p%v(j,i))
            vmax = max(vmax,p%v(j,i))
          end do
        end if
      end do
    end do
  end do
  end do
  end do
  !................................................-+.............................
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,0p,2f10.1)') &
      wc, wc/self%np, self%patch%id, self%np, 'KDK', vmin, vmax
  end if
  call trace%end (itimer)
END SUBROUTINE kdk

!===============================================================================
!> Re-assign particles to other cells.  New particles will be added at the
!> end, and holes will only appear for indices up to the initial array lenght
!===============================================================================
SUBROUTINE trade (self)
  class(particle_solver_t)                    :: self
  integer                                     :: i1, i2, i3, ip, ib, i, j
  integer                                     :: q1, q2, q3, pq(3), n3, me
  integer                                     :: j2, ne, nm, n1, n2, ne1, nm1
  real(8)                                     :: wc
  type(particles_t), pointer                  :: p, pn, pe
  real, pointer, contiguous, dimension(:,:)   :: pv, pr
  real, pointer, contiguous, dimension(:,:,:) :: p1, p2, p3
  logical, save                               :: warn(2) = .true.
  integer, save                               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%trade', itimer=itimer)
  wc = omp_get_wtime()
  n1 = 0
  n2 = 0
  n3 = 0
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
      p => self%particles(i1,i2,i3,ib)
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
          pq = p%q(j,:)
          if (any(pq < self%lb .or. pq > self%ub)) then
            ne1 = ne1 + 1
            if (pe%lp >= pe%mp) then
              me = 1.2*pe%mp
              !call pe%resize(me)
              if (verbose > 0) then
                write (io_unit%log,*) &
                  'export resize: mo, mp =', me, self%me
              end if
!              self%me = me
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
          else if (any(pq /= [i1,i2,i3])) then
            nm1 = nm1 + 1
            !q1 = pq(1)
            !q2 = pq(2)
            !q3 = pq(3)
            pn => self%particles(pq(1),pq(2),pq(3),ib)
            if (pn%lp < pn%mp) then
              pn%np = pn%np + 1           ! increment target count
              j2 = pn%lp + 1
              pn%lp = j2                  ! increment target length
              do i=1,3
                pn%r(j2,i) = p%r(j,i)
                pn%v(j2,i) = p%v(j,i)
                pn%q(j2,i) = p%q(j,i)
              end do
              pn%w(j2) = p%w(j)
              p%w(j) = 0.0                ! mark position empty
              p%np = p%np - 1             ! decrement source count
              !if (j == p%lp) then         ! if last slot
              !  p%lp = p%lp - 1           ! also decrement length
              !end if
            else
              if (warn(2)) then
                write (stderr,     *) 'WARNING: overflow in particle buffer', pn%mp
                write (io_unit%log,'(a,i6,4i4,3i7)') &
                  'WARNING: overflow in particle buffer', &
                  self%patch%id, q1, q2, q3, ib, pn%np, pn%lp, pn%mp
                warn(2) = .false.
                stop
              end if
              p%q(j,:) = [i1,i2,i3]       ! reset position
              n2 = n2 + 1
            end if
            if (verbose > 1) then
              call  p%validate ('p1 after ')
              call pn%validate ('pn after ')
            end if
          end if
        end if
      end do
      ne = ne + ne1
      nm = nm + nm1
      if (p%np == 0) write(io_unit%log,'(a,6i5)') '0p aft trade', p%id, p%lp
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
    write (io_unit%log,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,5i10)') &
      wc, wc/self%np, self%patch%id, self%np, 'trade', nm, ne, n1, n2, n3
  end if
  call trace%end (itimer)
END SUBROUTINE trade

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
  real(8)                     :: wc
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
        write (io_unit%log,'(a,i6,4i4,2i6,1p,e10.2)') &
        'ERROR: id,np,lp,dd=', &
        p%id,p%np,p%lp, self%dd(i1,i2,i3,ib)
        call io%abort ('Dust density non-positive')
      end if
      call io%assert (self%dd(i1,i2,i3,ib) > 0.0, 'Dust density non-positive')
      !-------------------------------------------------------------------------
      ! If too few particles, add new ones and set new weights to all
      !-------------------------------------------------------------------------
      if (p%np < n_min) then
        w = self%dd(i1,i2,i3,ib)/particles_per_bin
        !call p%pack (n)
        do j=1,p%np
          p%w(j) = w
        end do
        do j=p%np+1,particles_per_bin
          do i=1,3
            p%r(j,i) = self%random%ran1()
            p%v(j,i) = self%v(i,i1,i2,i3)
          end do
          p%q(j,:) = [i1,i2,i3]
          p%w(j) = w
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

END MODULE particle_solver_mod
