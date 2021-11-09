!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
MODULE read_external_mod
  USE io_mod
  USE units_mod
  USE scaling_mod
  USE trace_mod
  USE task_mod
  USE link_mod
  USE mesh_mod
  USE bits_mod
  USE halo_mod
  USE task_list_mod
  USE experiment_mod

  type, public:: read_external_t
    character(len=126) :: atmos_file
    class(task_list_t), pointer:: task_list
    integer :: nz
    real(8), dimension(:), allocatable :: zi, temp, pe, rho, vturb
    real(8) :: size
    character(len=126) :: atmos_name
  contains
    procedure:: init
    procedure:: read_snapshot
    procedure:: create_patches
    procedure:: load_variables_to
    procedure:: log_nbor_lists
  end type read_external_t

  type(read_external_t), public:: read_external
  integer :: verbose = 1

CONTAINS

SUBROUTINE init(self, task_list, name)
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(task_list_t),  pointer:: task_list
  character(len=*),   optional:: name
  !-----------------------------------------------------------------------------
  class(task_t),       pointer:: task
  class(experiment_t), pointer:: patch
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  logical :: searching
  integer :: iostat, npatches, patch_nz, iz, i
  character(len=126) :: atmos_file
  namelist /atmos_params/ atmos_file, npatches
  !-----------------------------------------------------------------------------
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%init ('read_external_t%init')
  call io%header ('read_external_t%init')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and saves them
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    ! default values
    atmos_file = './input_multi3d/atmos/atmos.sun_MARCS'

    rewind (io%input)
    read (io%input, nml=atmos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, atmos_params)
    end if
    call io%namelist_warning ('atmos_params', iostat)
    write (io%namelist, atmos_params)

    self%atmos_file = atmos_file
  endif
  call io%lock%unset

  allocate (task_list)
  self%task_list => task_list
  call self%task_list%init ('multi3d_mod')

  call scaling%init()
  call self%read_snapshot()

  self%size = self%zi(self%nz) - self%zi(1)

  patch_nz = self%nz / npatches

  do i=1, npatches
    iz = (i - 1) * patch_nz + 1

    if (i==npatches) patch_nz = self%nz - iz + 1

    call self%create_patches(patch,iz,patch_nz) ! patch is allocated here
    call self%load_variables_to(patch, iz)

    if (i==1) write(stdout,'(I2, A28)') npatches, ' subdomains in z-direction:'
    write(stdout,'(I4, I6,A3,I4)') i, iz, ' - ', iz + patch_nz - 1
  end do

  !-----------------------------------------------------------------------------
  ! Add all tasks to the ready queue, preventing a slow default add_by_time
  !-----------------------------------------------------------------------------
  call self%task_list%queue_all
  call self%task_list%print_queue

  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------

  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_levels'
  call self%task_list%init_levels
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_all_nbors'
  call self%task_list%init_all_nbors
  if (verbose > 0) write(stdout,*) 'read_external_t%init: reset_status'
  call self%task_list%reset_status

  call task_list%info
  call self%log_nbor_lists ()

  i = 1
  link => self%task_list%head
  do while (associated(link))
    !print *, 'info on task number ', i
    !print '(A15, 3F12.6, A4)', 'postition = ', link%task%position, 'Mm'
    !print '(A15, 3F12.6, A4)', 'size = ', link%task%size, 'Mm'
    !print '(A15, 3F12.6, A4)', 'lower = ', link%task%position + link%task%authority%l, 'Mm'
    !print '(A15, 3F12.6, A4)', 'upper = ', link%task%position + link%task%authority%u, 'Mm'

    !print *, link%task%mesh%lower_boundary
    !print '(A15, 3F12.6, A4)', 'box = ', link%task%box, 'Mm'
    link => link%next
    i = i + 1
  end do

  !print *, shape(patch%mem)
  !print *, 'density', patch%mem(3,3,:,patch%idx%d,1,1)
  !print *, 'temperature', patch%mem(3,3,:,patch%idx%e,1,1)
  !print *, 'pressure', patch%mem(3,3,:,6,1,1)

  call trace%end()
END SUBROUTINE init

!===============================================================================
subroutine log_nbor_lists(self)
  use m3d_rt_mod, only: m3d_rt_t
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    
    select type(lt => link%task) 
    class is (m3d_rt_t)
      write (stdlog,'(a,i6,":",i2.2,i5.5)') 'log_nbor_lists: link=', &
        lt%id, lt%level, lt%i_omega
      nbor => link%nbor
      do while (associated(nbor))
        
        select type (nb => nbor%task)
        class is (m3d_rt_t)
          write (stdlog,'(2x,a,i6,":",i2.2,i5.5)') 'log_nbor_lists: nbor=', &
            nb%id, nb%level ,nb%i_omega
        end select
        
        nbor => nbor%next
      end do
    class default
      write (stdlog,*) 'non-m3d_rt patch.',link%task%id,link%task%kind

    end select
    link => link%next

end do

end subroutine log_nbor_lists

!===============================================================================
!> Initialize a minimum amount of patch information in a Cartesian arrangement
!> This assumes that positions lie btw 0.0 and self%size
!===============================================================================
SUBROUTINE create_patches (self, patch, iz, patch_nz)
  class (read_external_t):: self
  class (experiment_t), pointer:: patch
  class(link_t), pointer:: link
  !.............................................................................
  integer :: iz, patch_nz
  class (task_t), pointer:: task
  class(*), pointer:: ptr
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%create_patches')
  allocate (patch)

  patch%task_list => self%task_list

  patch%mhd = .false.
  patch%kind = 'hlld_mhd_patch'
  call patch%idx%init (energy=5, mhd=patch%mhd)
  call patch%extras_t%pre_init
  patch%pre_initialized = .true.

  patch%n = [1, 1, patch_nz]
  patch%ng = [5,5,1]
  patch%nv = 6
  patch%nt = 1
  patch%nw = 1

  call patch%init
  patch%periodic = [.true., .true., .false.]

  patch%level  = 1
  patch%box(1) = 1.0d0
  patch%box(2) = 1.0d0
  patch%box(3) = self%size

  patch%size(1) = patch%box(1)
  patch%size(2) = patch%box(2)
  patch%size(3) = self%zi(iz + patch_nz - 1) - self%zi(iz)

  if (iz > 1) then
    patch%size(3) = patch%size(3) + 0.5d0*(self%zi(iz) - self%zi(iz-1))
  endif

  if ((iz + patch_nz) < self%nz) then
    patch%size(3) = patch%size(3) + 0.5d0*(self%zi(iz + patch_nz) - self%zi(iz + patch_nz - 1))
  endif

  patch%ds = patch%size / (patch_nz - 1)

  patch%position(1) = 0.5d0*patch%size(1)
  patch%position(2) = 0.5d0*patch%size(2)
  patch%position(3) = self%zi(iz) + 0.5d0*patch%size(3)

  !---------------------------------------------------------------------------
  ! Make sure interpolations are done optimally
  !---------------------------------------------------------------------------
  patch%unsigned(patch%idx%d) = .true.
  patch%unsigned(patch%idx%e) = .true.
  patch%pervolume(patch%idx%px) = .true.
  patch%pervolume(patch%idx%py) = .true.
  patch%pervolume(patch%idx%pz) = .true.

  call patch%set (bits%has_output)
  !---------------------------------------------------------------------------
  allocate (patch%filled(patch%gn(1),patch%gn(2),patch%gn(3)))
  call io%bits_mem (storage_size(patch%filled), &
                    product(shape(patch%filled)), 'fill')
  patch%filled = -1
  !---------------------------------------------------------------------------
  ! Add the patch to the task list -- this creates the patch%link
  !---------------------------------------------------------------------------
  allocate (halo_t:: link)
  link%task => patch
  patch%link => link
  call self%task_list%prepend_link (link)
  call trace%end ()
END SUBROUTINE create_patches


!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
SUBROUTINE read_snapshot(self)
  !-----------------------------------------------------------------------------
  USE m3io_mod, only: NextEntry
  class(read_external_t):: self
  !-----------------------------------------------------------------------------
  integer :: lunit, i, iostat
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  open(file=self%atmos_file, newunit=lunit)
  call nextentry(lunit)
  read(lunit, *) self%atmos_name
  call nextentry(lunit)
  read(lunit, *) self%nz

  allocate(self%zi(self%nz))
  allocate(self%temp, self%pe, self%rho, self%vturb, mold=self%zi)

  do i=1, self%nz
    call nextentry(lunit)
    read(lunit, *) self%zi(i), self%temp(i), self%pe(i), self%rho(i), self%vturb(i)
  end do

  self%zi  = self%zi  / scaling%l
  self%rho = self%rho / scaling%d
  self%pe  = self%pe  / scaling%p
END SUBROUTINE read_snapshot

!===============================================================================
!> Takes variables read from text file and loads them into patch memory
!===============================================================================
SUBROUTINE load_variables_to(self, patch, iz)
  !-----------------------------------------------------------------------------
  class (read_external_t):: self
  class (experiment_t), pointer:: patch
  integer :: i,j,k, it, iz, ii

  it = 1 ! time index, it=1 is the initial condition

  do k=1, patch%gn(3)
    ii = iz + (k - patch%li(3))

    if (ii < 1) ii = 1
    if (ii > self%nz) ii = self%nz

    do j=1, patch%gn(2)
      do i=1, patch%gn(1)
        patch%mem(i,j,k,patch%idx%d ,it,1) = self%rho(ii)
        patch%mem(i,j,k,patch%idx%px,it,1) = 0
        patch%mem(i,j,k,patch%idx%py,it,1) = 0
        patch%mem(i,j,k,patch%idx%pz,it,1) = 0
        patch%mem(i,j,k,patch%idx%e ,it,1) = self%temp(ii)
        patch%mem(i,j,k,          6 ,it,1) = self%pe(ii)
      end do
    end do
  end do

  patch%vgas = 0

  !-----------------------------------------------------------------------------
  ! Convert variables into solver's preference
  !-----------------------------------------------------------------------------
  !call patch%from_generic()
END SUBROUTINE

END MODULE