!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE bits_mod
  USE kinds_mod
  USE link_mod
  USE list_mod
  USE patch_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private

  type:: refine_t
    integer:: check_interval                    ! How often do we check refinement criteria
                                                ! (given in number of steps/updates between checks)
    integer:: refinement_ratio                  ! How does a patch split when it is refined
    real:: min_dx, max_dx
    real:: refine_from_t, refine_until_t        ! Time limits for occurence of refinement
    real:: min_jeans                            ! Jeans number below which a patch is flagged for refinement
    real:: tol_shock, tol_contact               ! Tolerance values for shock/discontinuity detectors
    integer:: ngradvar=10
    real, dimension(:), pointer:: tol_grad      ! Tolerance values for gradient detector
    character(len=32), dimension(:), pointer:: grad_var
    real, dimension(3):: region_llc, region_urc ! Sometimes you want to refine a region (manually)
    logical:: do_refine                         ! Flag to turn the entire process on or off
  contains
    procedure:: init
    procedure:: check_current
    procedure:: check_jeans
    procedure:: shock_and_cd_detector
    procedure:: refine_region
    procedure:: split_patch
    procedure:: gradient_detector
  end type

  type(refine_t), public:: refine

CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer:: check_interval=50, refinement_ratio=2
  real:: refine_from_t=0.0, refine_until_t=1.0e32, min_jeans=0.0, min_dx=0.0, max_dx=1.0e32, &
         tol_shock=0.0, tol_contact=0.0
  real, dimension(3):: region_llc=0.0, region_urc=0.0
  real, dimension(self%ngradvar):: tol_grad
  character(len=32), dimension(self%ngradvar):: grad_var
  logical:: do_refine=.false.
  logical, save:: first_time = .true.
  namelist /refine_params/ check_interval, refinement_ratio, refine_from_t, refine_until_t, &
                           min_dx, max_dx, min_jeans, tol_shock, tol_contact, tol_grad, &
                           grad_var, region_llc, region_urc, do_refine
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%init')

  if (first_time) then
    first_time=.false.
    ! allocate and initialise gradient detector variables
    allocate (self%tol_grad(self%ngradvar), self%grad_var(self%ngradvar))
    tol_grad(:) = 0.0
    grad_var(:) = ""
    rewind (io%input)
    read (io%input, refine_params, iostat=iostat)
    if (io%master) write (*, refine_params)
  end if
  self%check_interval = check_interval
  self%refinement_ratio = refinement_ratio
  self%refine_from_t = refine_from_t
  self%refine_until_t = refine_until_t
  self%min_dx = min_dx
  self%max_dx = max_dx
  self%min_jeans = min_jeans
  self%tol_shock = tol_shock
  self%tol_contact = tol_contact
  self%tol_grad = tol_grad
  self%grad_var = grad_var
  self%region_llc = region_llc
  self%region_urc = region_urc
  self%do_refine = do_refine

  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Check if refinement is desired and needed for the current task; if so, push
!> new tasks onto the (ready) queue with the same time as the parent, i.e., to
!> the head of the queue.
!> =============================================================================
SUBROUTINE check_current (self, tasklist, link, was_refined)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  logical:: needed
  logical, intent(out):: was_refined
  real:: dxmin, dxmax
  integer, save:: itimer=0
  !.............................................................................
  was_refined = .false.
  if (.not. self%do_refine) return
  call trace%begin('refine_t%check_current', itimer=itimer)

  patch => task2patch(link%task)
  needed = .false.

  ! filter out patches which do not match the necessary criteria
  if (patch%is_set(bits%virtual)) then
    call trace%end (itimer)
    return
  end if
  if (patch%istep < patch%check_refine_next) then
    call trace%end (itimer)
    return
  end if
  if (patch%is_set(bits%not_leaf)) then
    call trace%end (itimer)
    return
  end if
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (dxmin < self%min_dx .or. dxmax >= self%max_dx) then
    call trace%end (itimer)
    return
  end if
  if (patch%time < self%refine_from_t .or. patch%time > self%refine_until_t) then
    call trace%end (itimer)
    return
  end if

  !------------------------------------------------------------------------
  ! Refinement Criteria
  !------------------------------------------------------------------------
  ! Jeans criterion
  if (self%min_jeans > 0.0) needed = needed .or. self%check_Jeans(patch)

  ! shock and/or contact discontinuity detector
  if (self%tol_shock > 0.0 .or. self%tol_contact > 0.0) &
    needed = needed .or. self%shock_and_cd_detector(patch)

  ! gradient detector
  if (any(self%tol_grad > 0.0) ) needed = needed .or. self%gradient_detector(patch)

  ! refine on (a) user-specified region(s)
  if (any(self%region_llc /= 0.0) .and. any(self%region_urc /= 0.0)) &
    needed = needed .or. self%refine_region(patch)

  patch%check_refine_next = patch%istep + self%check_interval

  !------------------------------------------------------------------------
  ! If refinement has been deemed necessary, do it!
  !------------------------------------------------------------------------
  if (needed) then
    if (io%verbose > 1) print*,'refine_t%needed!; patch:', patch%id
    call self%split_patch (tasklist, link)
    was_refined = .true.
  end if

  call trace%end (itimer)
END SUBROUTINE check_current

!> =============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!> =============================================================================
SUBROUTINE make_child_patch (tlist, parent_link, parent, pos, child_size)
  USE download_mod

  class(list_t) :: tlist
  class(link_t), pointer :: parent_link, child_link
  class(patch_t), pointer :: parent, child
  real, dimension(3), intent(in) :: pos, child_size
  !------------------------------------------------------------------------
  call trace%begin('refine_mod::make_child_patch')
  !------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !
  ! Using sourced allocation below effectively clones the parent.
  !------------------------------------------------------------------------
  allocate (child, source=parent)
  child%mem_allocated = .false.
  nullify (child%mem)
  nullify (child%mesh)
  child%id = 0                                       ! to ensure it gets its own ID
  child%size = child_size                            ! new size
  child%position = pos                               ! patch position
  child%llc_cart = child%position - merge(0.5,0.0,child%n>1) * child%size  ! lower left corner
  child%llc_nat = child%llc_cart                     ! FIXME!: this form only applies for Cartesian coords.
  call child%init()                                  ! patch, for a new id
  child%level = parent%level + 1                     ! override the calculation in `patch_t%init`.
  child%dt(child%it) = parent%dt(parent%it)          ! reset the time-step
  child%t(child%it) = child%time
  child%out_next = child%time                        ! force output now
  child%iout = 0                                     ! reset output counter
  !------------------------------------------------------------------------
  ! Reset flags (FIXME: Would it be better to just reset `status` to zero?)
  !------------------------------------------------------------------------
  call child%clear(bits%root_grid)
  call child%clear(bits%ready)
  call child%clear(bits%busy)
  !------------------------------------------------------------------------
  ! Append links to the new patch in the patch list, sorted by quality
  ! Neighbour lists are constructed at the end of `split_patch`.
  !------------------------------------------------------------------------
  allocate (child_link)
  child_link%task => child
  child_link%parent => parent_link                          ! link is the link to the new patch
  call tlist%append_link (child_link)                       ! which is added to the patch list
  !------------------------------------------------------------------------
  ! Download values from parent and neighbors, except at time = 0, to allow
  ! recursive refinement there.
  ! If time = 0, patch should have been filled with values during `child%init()`
  ! call.
  !------------------------------------------------------------------------
  if (parent%time > 0.0) call download%download_link (child_link, all_cells=.true.)
  if (io%verbose>0) then
    1 format(a,3i6,f10.6,2(2x,3f8.4))
    print 1,'refine_mod::make_child_patch:  child, level, iout, time, size =', &
      child%id, child%level, child%iout, child%time, child%size
    print 1,'                              parent, level, iout, time, size =', &
      parent%id, parent%level, parent%iout, parent%time, parent%size
  end if

  call tlist%check_ready(child_link)

  call trace%end
END SUBROUTINE make_child_patch

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split_patch (self, tlist, this)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: this, link
  class(patch_t), pointer:: this_patch
  integer:: ix, iy, iz
  integer:: nsplit(3)
  real, dimension(3):: pos, rpos, split_size
  !.............................................................................
  if (this%task%is_set(bits%not_leaf)) return
  call trace%begin('refine_t%split_patch', 1)
  nsplit = self%refinement_ratio

  this_patch => task2patch(this%task)

  where (this_patch%size==0 .or. this_patch%n==1)
    nsplit = 1
  end where

  split_size = this_patch%size / nsplit
  do ix=1,nsplit(1)
    rpos(1) = (ix-(nsplit(1)+1)*0.5) * split_size(1)
    do iy=1,nsplit(2)
      rpos(2) = (iy-(nsplit(2)+1)*0.5) * split_size(2)
      do iz=1,nsplit(3)
        rpos(3) = (iz-(nsplit(3)+1)*0.5) * split_size(3)
        pos = rpos + this_patch%position
        call make_child_patch(tlist, this, this_patch, pos, split_size)
      end do
    end do
  end do
  call this_patch%set(bits%not_leaf)
  
  call tlist%remove_parents
  call tlist%init_all_nbors
  call tlist%reset_status
  link => tlist%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call tlist%check_ready (link)
    link => link%next
  end do

  if (io%verbose > 0) call tlist%print_tasks()

  call trace%end
END SUBROUTINE split_patch

!> =============================================================================
!> Calculate the Jeans length and compare it to the patch resolution.
!> =============================================================================
FUNCTION check_jeans (self, patch) RESULT(out)
  USE scaling_mod
  class(refine_t):: self
  class(patch_t), pointer:: patch
  logical:: out
  real:: jeans_length, ratio
  real, dimension(:,:,:), allocatable:: c_sound2
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  !.............................................................................
  call trace%begin ('refine_t%check_check_jeans', 2)

  allocate(c_sound2(patch%mesh(1)%gn,patch%mesh(2)%gn,patch%mesh(3)%gn))

  d => patch%mem(:,:,:,patch%i_d,patch%it,1)
  e => patch%mem(:,:,:,patch%i_e,patch%it,1)
  c_sound2 = patch%gamma * eos%pressure(d,e) / d
  jeans_length = minval(sqrt(c_sound2 / (d*scaling%grav)))
  ratio = jeans_length / minval(patch%mesh%d)
  out = ratio < self%min_jeans

  if (io%verbose > 2) print*,'Jeans:', patch%id, ratio, self%min_jeans, out

  deallocate (c_sound2)
  nullify (d,e)
  call trace%end
END FUNCTION check_jeans

!> =============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!> =============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine_it)
  USE mesh_mod

  class(refine_t):: self
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: mx, my, mz
  integer:: ix, iy, iz, ione=0, jone=0, kone=0
  real(kind=KindScalarVar):: dp1, dp2, dp3, dv1, dv2, dv3, dd1, dd2, dd3
  real, parameter:: small = 1.0e-10
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: pg
  real(kind=KindScalarVar), dimension(:,:,:,:), allocatable:: v
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  logical:: refine_it
  !.............................................................................
  call trace%begin('refine_t%shocks_and_cd_detector')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  if (mx%n > 1) ione = 1
  if (my%n > 1) jone = 1
  if (mz%n > 1) kone = 1

  allocate(pg(patch%mesh(1)%gn,patch%mesh(2)%gn,patch%mesh(3)%gn), &
           v(patch%mesh(1)%gn,patch%mesh(2)%gn,patch%mesh(3)%gn,3))
  ! we need the gas pressure, density and velocity
  select type (patch)
  class is (solver_t)
    d => patch%mem(:,:,:,patch%i_d,patch%it,1)
    e => patch%mem(:,:,:,patch%i_e,patch%it,1)
    pg = eos%pressure(d, e)
    v  = patch%gas_velocity()
  end select

  refine_it = .false.
  do iz=mz%li,mz%ui
    do iy=my%li,my%ui
      do ix=mx%li,mx%ui
        dp1 = abs(pg(ix+ione,iy,iz) - pg(ix-ione,iy,iz)) &
            / (min(pg(ix+ione,iy,iz), pg(ix-ione,iy,iz)) + small )
        dp2 = abs(pg(ix,iy+jone,iz) - pg(ix,iy-jone,iz)) &
            / (min(pg(ix,iy+jone,iz), pg(ix,iy-jone,iz)) + small )
        dp3 = abs(pg(ix,iy,iz+kone) - pg(ix,iy,iz-kone)) &
            / (min(pg(ix,iy,iz+kone), pg(ix,iy,iz-kone)) + small )

        if (self%tol_shock > 0.0) then
          dv1 = v(ix  ,iy,iz,1) - v(ix+ione,iy,iz,1)
          dv2 = v(ix,iy  ,iz,2) - v(ix,iy+jone,iz,2)
          dv3 = v(ix,iy,iz  ,3) - v(ix,iy,iz+kone,3)

          if ( dp1 > self%tol_shock .and. dv1 > small ) refine_it = .true.
          if ( dp2 > self%tol_shock .and. dv2 > small ) refine_it = .true.
          if ( dp3 > self%tol_shock .and. dv3 > small ) refine_it = .true.
        end if

        if (self%tol_contact > 0.0) then
          dd1 = abs(d(ix+ione,iy,iz) - d(ix-ione,iy,iz)) &
              / (min(d(ix+ione,iy,iz), d(ix-ione,iy,iz)) + small )
          dd2 = abs(d(ix,iy+jone,iz) - d(ix,iy-jone,iz)) &
              / (min(d(ix,iy+jone,iz), d(ix,iy-jone,iz)) + small )
          dd3 = abs(d(ix,iy,iz+kone) - d(ix,iy,iz-kone)) &
              / (min(d(ix,iy,iz+kone), d(ix,iy,iz-kone)) + small )

          if ( dd1 > self%tol_contact .and. dp1 < 0.01*self%tol_contact ) refine_it = .true.
          if ( dd2 > self%tol_contact .and. dp2 < 0.01*self%tol_contact ) refine_it = .true.
          if ( dd3 > self%tol_contact .and. dp3 < 0.01*self%tol_contact ) refine_it = .true.
        end if
      end do
    end do
  end do

  nullify (mx, my, mz, d)
  deallocate (pg, v)

  call trace%end
END FUNCTION shock_and_cd_detector

!> =============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!> =============================================================================
FUNCTION gradient_detector (self, patch) RESULT (refine_it)
  USE mesh_mod
  USE solver_mod

  class(refine_t):: self
  class(patch_t), pointer:: patch
  class(mesh_t), pointer:: mx, my, mz
  integer:: i, j, k, ione=0, jone=0, kone=0, im1, jm1, km1, ip1, jp1, kp1, n
  real, parameter:: small=1e-10
  real(kind=KindScalarVar):: dqmax
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: q, dq1, dq2, dq3
  real, parameter:: rsmall = 1.0e-10
  logical:: refine_it
  !.............................................................................
  call trace%begin('refine_t%gradient_detector')

  allocate (dq1(patch%gn(1),patch%gn(2),patch%gn(3)), &
            dq2(patch%gn(1),patch%gn(2),patch%gn(3)), &
            dq3(patch%gn(1),patch%gn(2),patch%gn(3)), &
            q  (patch%gn(1),patch%gn(2),patch%gn(3)))

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  if (mx%n > 1) ione = 1
  if (my%n > 1) jone = 1
  if (mz%n > 1) kone = 1

  refine_it = .false.
  do n=1,size(self%grad_var)
    if (trim(self%grad_var(n)) == '') exit
    if (trim(self%grad_var(n)) == 'd')  q = patch%mem(:,:,:,patch%i_d,patch%it,1)
    if (trim(self%grad_var(n)) == 'e')  q = patch%mem(:,:,:,patch%i_e,patch%it,1)
    if (trim(self%grad_var(n)) == 'p1') q = patch%mem(:,:,:,patch%i_px,patch%it,1)
    if (trim(self%grad_var(n)) == 'p2') q = patch%mem(:,:,:,patch%i_py,patch%it,1)
    if (trim(self%grad_var(n)) == 'p3') q = patch%mem(:,:,:,patch%i_pz,patch%it,1)
    if (trim(self%grad_var(n)) == 'b1') q = patch%mem(:,:,:,patch%i_bx,patch%it,1)
    if (trim(self%grad_var(n)) == 'b2') q = patch%mem(:,:,:,patch%i_by,patch%it,1)
    if (trim(self%grad_var(n)) == 'b3') q = patch%mem(:,:,:,patch%i_bz,patch%it,1)
    select type(patch)
    class is (solver_t)
      if (trim(self%grad_var(n)) == 'p')  q = patch%gas_pressure()
      if (trim(self%grad_var(n)) == 'v1') q = patch%gas_velocity(1)
      if (trim(self%grad_var(n)) == 'v2') q = patch%gas_velocity(2)
      if (trim(self%grad_var(n)) == 'v3') q = patch%gas_velocity(3)
    end select

    ! calculate the gradient in each direction
    do k=mz%lo,mz%uo
      km1 = k - kone
      kp1 = k + kone
      do j=my%lo,my%uo
        jm1 = j - jone
        jp1 = j + jone
        do i=mx%lo,mx%uo
          im1 = i - ione
          ip1 = i + ione
          dq1(i,j,k) = abs ( q(ip1,j,k) - q(im1,j,k) ) &
                           / ( max ( q(ip1,j,k), q(im1,j,k) ) &
                             + small )
          dq2(i,j,k) = abs ( q(i,jp1,k) - q(i,jm1,k) ) &
                           / ( max ( q(i,jp1,k), q(i,jm1,k) ) &
                             + small )
          dq3(i,j,k) = abs ( q(i,j,kp1) - q(i,j,km1) ) &
                           / ( max ( q(i,j,kp1), q(i,j,km1) ) &
                             + small )
        end do
      end do
    end do

    ! choose max. gradient; if the max. is greater than the threshold, flag
    ! for refinement
    do k=mz%li+kone,mz%ui-kone
      km1 = k - kone
      kp1 = k + kone
      do j=my%li+jone,my%ui-jone
        jm1 = j - jone
        jp1 = j + jone
        do i=mx%li+ione,mx%ui-ione
          im1 = i - ione
          ip1 = i + ione
          dqmax = max ( dq1(i,j,k), dq1(ip1,j,k), dq1(im1,j,k) &
                      , dq2(i,j,k), dq2(i,jp1,k), dq2(i,jm1,k) &
                      , dq3(i,j,k), dq3(i,j,kp1), dq3(i,j,km1) )
          if (dqmax > self%tol_grad(n)) refine_it = .true.
        end do
      end do
    end do

  end do

  deallocate (dq1, dq2, dq3, q)
  call trace%end
END FUNCTION gradient_detector

!> =============================================================================
!> Given a box defined by its lower and upper corners, flag a patch for refinement
!> if it lies inside this box, up to a maximum level.
!> This particular criteria is primarily used for testing and debugging.
!> =============================================================================
FUNCTION refine_region (self, patch) RESULT (refine_it)
  USE mesh_mod

  class(refine_t) :: self
  class(patch_t), pointer :: patch
  class(mesh_t), dimension(:), pointer:: m
  logical :: intersect
  logical :: refine_it
  integer :: ix, iy, iz
  real :: x, y, z
  real, dimension(3) :: p_llc, p_urc, b_llc, b_urc
  !.............................................................................
  call trace%begin('refine_t%refine_on_box')

  refine_it = .false.
  m => patch%mesh

  ! check if the patch (as a whole) intersects the box-shaped region of interest
  do ix=1,3
    associate (mr => m(ix)%r)
    p_llc(ix) = m(ix)%p + mr(m(ix)%li) - 0.5 * m(ix)%d
    p_urc(ix) = m(ix)%p + mr(m(ix)%ui) + 0.5 * m(ix)%d
    end associate
  end do
  b_llc = self%region_llc
  b_urc = self%region_urc

  intersect = all(b_llc > p_llc .and. b_llc < p_urc) .or. &
              all(b_urc > p_llc .and. b_urc < p_urc)
  if (.not. intersect) then
    call trace_end
    return
  endif

  do iz=m(3)%li,m(3)%ui
    z = m(3)%p + m(3)%r(iz)
    do iy=m(2)%li,m(2)%ui
      y = m(2)%p + m(2)%r(iy)
      do ix=m(1)%li,m(1)%ui
        x = m(1)%p + m(1)%r(ix)
        refine_it = refine_it .or. (x >= b_llc(1) .and. x <= b_urc(1) .and. &
                                    y >= b_llc(2) .and. y <= b_urc(2) .and. &
                                    z >= b_llc(3) .and. z <= b_urc(3) )
      end do
    end do
  end do

  call trace%end
END FUNCTION refine_region

END MODULE refine_mod
