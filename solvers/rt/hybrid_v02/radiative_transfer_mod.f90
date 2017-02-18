!> =============================================================================
!> Module for setting up the radiative transfer, and the relations between the
!> RT tasks and ray systems.
!> =============================================================================
MODULE radiative_transfer_mod
  USE bits_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mesh_mod
  USE mhd_mod
  USE omp_mod
  USE patch_mod
  USE ray_geometry_v3_mod
  USE rt_mod
  USE rt_boundary_mod
  USE scaling_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod
  implicit none
  private
  PUBLIC rt_setup, connect_nbors
CONTAINS

!> =============================================================================
!> The only procedure callable from the outside
!> =============================================================================
SUBROUTINE rt_setup (self)
  class(task_list_t)          :: self
  class(task_list_t), pointer :: rt_task_list                           ! local task list
  class(list_t),      pointer :: list
  class(link_t),      pointer :: t_list, link, next, rt_link
  class(task_t),      pointer :: task
  class(mesh_t),      pointer :: mesh
  class(mhd_t),       pointer :: mhd_patch
  class(rt_t),        pointer :: rt_task
  logical, save               :: first_time=.true.
  character(len=16)           :: RT_solver ='integral'      ! choice of RT solver
  integer :: n_rt
  real,    save :: cr_heating=4.4e-27
  integer, save :: solver
  integer, save :: rt_freq=1
  logical, save :: store_rt = .false.
  logical, save :: do_chains = .true.
  namelist /rt_params/ RT_solver, cr_heating, rt_freq, store_rt, do_chains
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer::rt_setup')
  associate (mhd => self%head%task)
  SELECT TYPE(mhd)
  class is (mhd_t)
    mhd_patch => mhd
  class default
    error stop '(1) in rt_setup'
  end SELECT
  end associate
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, rt_params)
    if (io%master) write (*, rt_params)
    SELECT case(RT_solver)
    case('feautrier')
      solver = 0
    case('integral')
      solver = 1
    end SELECT
  end if

  mesh => mhd_patch%mesh(1)
  select type(mesh)
    type is(Cartesian_mesh)
      call ray_geometry%init(geometry=1, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds, do_chains=do_chains)  
    type is(spherical_mesh)
      call ray_geometry%init(geometry=2, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds, do_chains=do_chains) 
    type is(cylindrical_mesh)
      call ray_geometry%init(geometry=3, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds, do_chains=do_chains)   
  end select
  nullify(mesh)

  write (io_unit%log,*) 'Ray geometry: n_ray_dirs:', ray_geometry%head%n_ray_dirs
  !call microphysics%init
  call rt_boundary%init
  rt_boundary%ubox = mhd_patch%box                                      ! FIXME : add bndry patches
  rt_boundary%lbox = 0.0
  allocate(rt_task_list)
  call rt_task_list%list_t%init ('rt task list')
  t_list => self%head
  do while (associated(t_list))
    !-------------------------------------------------------------------
    ! Get the right type of the task to pick up other properties (mhd_t)
    !-------------------------------------------------------------------
    associate (mhd => t_list%task)
    SELECT TYPE(mhd)
    class is (mhd_t)
      allocate (rt_task, link)                                          ! the internal RT task and its link
      mhd_patch => mhd
      n_rt = 0
      !-----------------------------------------------------------------
      ! Construct the RT task and add it to the task list, 
      ! giving it a nbor list consisting of only the MHD task.
      !-----------------------------------------------------------------
      link%task         => rt_task                                      ! connect to RT task
      rt_task%link => link
      call connect_nbors(link,t_list)
      if (store_rt) rt_task%do_output = .true.
      rt_task%rt_freq    =  rt_freq                                     ! integer frequency, relative to MHD task
      rt_task%position = mhd_patch%position
      rt_task%box      = mhd_patch%box
      rt_task%size     = mhd_patch%size
      rt_task%dim      = mhd_patch%mesh%gn                              ! RT dimensions 
      rt_task%level    = mhd_patch%level 
      call rt_task%init                                                 ! initialise
      rt_task%size = mhd_patch%size
      rt_task%ds   = mhd_patch%ds
      n_rt = n_rt+1
      call rt_task_list%append_link(link)                               ! insert rt_task into the rt task list
    end SELECT
    end associate
    t_list => t_list%next                                               ! go to the next mhd task
  end do
  
  call make_chains (rt_task_list)                                       ! inter-connect all the chains
  link => rt_task_list%head
  do while (associated(link))
    task => link%task
    next => link%next
    !---------------------------------------------------------------------------
    ! Append the new task list link to the task_list and queue it
    !---------------------------------------------------------------------------
    select type (task)
    class is (rt_t)
      call self%append_link (link)
      link%task%time = -999d0
    end select
    link => next
  end do
  call self%print
  call trace_end
  if (io%master) then
    print *, '=================================='
    print *, 'INIT RT DONE'
    print *, '=================================='
  end if
END SUBROUTINE rt_setup


!> =====================================================================
!> Finalize the bundle chains by connecting them to the appropriate
!> upstream and downstream.
!> INPUT: rt_task_list
!> OUTPUT: none
!> RESULT: rt chains are connected to each other. Master chains are
!> identified and all chain inter-connections are done.
!> SUBROUTINES, CALLED FROM WITHIN: 
!>   rt_mod::rt_chain_t%finish_chain
!>   rt_mod::rt_chain_t%connect_to_upstream_diff_level
!>   rt_mod::rt_chain_t%ring_connects_to
!> FLOW:
!> 1) Subroutine loops over the rt_task list, takes one patch and loops
!>   over all it's chain rings.
!> 2) The appropriate chain ring is tested against being at a boundary.
!>   If that is the case, it is marked appropriately.
!> 3) Initially the ring is marked as being the "master" ring, meaning 
!>   the chain reaction starts from it.
!> 4) The search for an upstream ring begins. A loop over all patches in
!>   the rt_task list is done until the right chain is found. If an
!>   upstream is found, the "master" flag is removed and chains are connected.
!> 5) This new connection is reported through the chain (if there are
!>   more rings connected downstream), so every ring would know, that a
!>   new master has taken the throne. The loop (4) ends.
!> 6) If the ring encounters an upstream with different refinement level,
!>   it doesn't lose it's master title, but connects the upstream to
!>   a separate link (upstream, which is an array of pointers). In such
!>   case, a subroutine to connect the appropriate intensity positions 
!>   is called
!> 6) When the loop through rt tasks is done, and the ring still hasn't
!>   lost it's title as a master, it is put to the master list. For sim-
!>   plicity, chain_ring_t structure is used here, since rings are not
!>   tasks. In order not to mess up the ring structure, the connections
!>   are done backwards (ring%prev) - master rings do not have %prev 
!>   companions, so this is safe. Upon finalizing the chains, the %prev
!>   connection is broken.
!> 7) When the loop (1) ends, another loop, over the master chains (6) is
!>   done. Starting from each master ring, all rings are finalized to point
!>   their incoming intensity locations to where they are supposed to
!>   point, see rt_mod::rt_chain_t%finish_chain for more details.
!> 8) After (7) is complete, this subroutine has done it's job.
!> =====================================================================
SUBROUTINE make_chains (self)
  class(list_t) :: self
  !.....................................................................
  class(link_t),       pointer :: rt_list
  class(link_t),       pointer :: rt, try
  class(chain_ring_t), pointer :: master, origin, temp
  class(task_t),       pointer :: this, that
  class(rt_t),         pointer :: ithis, ithat
  logical :: connects = .false.
  logical :: is_upstream = .false.
  logical :: connected = .false.
  integer:: ibnd, irn
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  allocate(master)                                                      ! temporary master chain list
  origin => master
  rt_list => self%head
  do while (associated(rt_list))                                        ! loop through rt_task_list
    this => rt_list%task
    SELECT TYPE (this)
    class is (rt_t)                                                     ! select the right task
      ithis => this
    end SELECT
    do ibnd = 1, ithis%n_bundle                                         ! loop over chain rings
      connected = .false.                                               ! no connections have been made
      ithis%my_rings(ibnd)%master = .true.                              ! flag the ring as being the master
      if (ithis%has_OB(ithis%my_rings(ibnd))) then                      ! if this ring is at the outer boundary
        ithis%my_rings(ibnd)%BC(1) = .true.                             ! then flag it to have an outer boundary condition
        master%prev => ithis%my_rings(ibnd)
        master => master%prev                                           ! move to this ring
        goto 9                                                          ! since it is a boundary ring, it can't have upstreams or other masters
      end if
        ! TODO: inner boundaries.
        !if (ithis%has_IB(ithis%my_rings(ibnd))) then                    ! if this ring is at the inner boundary
        !  error stop "radiative_transfer_mod:: make_chains (1)"         ! not implemented yet!
        !  ithis%my_rings(ibnd)%BC(2) = .true.                           ! then flag it to have an inner boundary condition
        !  master%prev => ithis%my_rings(ibnd)
        !  goto 9                                                        ! since it is a boundary ring, it can't have upstreams or other masters
        !end if

        
      try => self%head
      do while (associated(try))                                        ! loop through rt_task_list (4)
        that => try%task
        SELECT TYPE (that)
        class is (rt_t)
          ithat => that
        end SELECT
        
        if (ithis%id /= ithat%id) then                                  ! not the same patch
          if (ithis%overlaps(ithat)) then                               ! patches overlap
            ! loop over ithat rings
            do irn=1, ithat%n_bundle                                    ! loop over nbor rings
              connects = .false.                                        ! can rings connect?
              is_upstream = .false.                                     ! nbor is different lvl upstreamer?
              call ithis%my_rings(ibnd)%ring_connects_to(ithat%my_rings(irn),connects, is_upstream)
              if (connects) then                                        ! the rings connect
                ithis%my_rings(ibnd)%master = .false.                   ! i am no longer the master
                ithis%my_rings(ibnd)%prev => ithat%my_rings(irn)        ! this is my superior
                ithat%my_rings(irn)%next => ithis%my_rings(ibnd)        ! this is my inferior
                
                if(.not.associated(ithis%my_rings(ibnd)%tail)) then
                  ithat%my_rings(irn)%tail => ithis%my_rings(ibnd)      ! connect the tails
                else
                  ithat%my_rings(irn)%tail => ithis%my_rings(ibnd)%tail 
                end if
                
                temp => ithis%my_rings(ibnd)                            ! inform the rest of the chain, that the new master has arrived
                do while(associated(temp))
                  temp%head => ithat%my_rings(irn)
                  temp => temp%next
                end do
                connected = .true.
                go to 9
              end if
              if (is_upstream) call &
              ithis%my_rings(ibnd)%connect_to_upstream_diff_level(ithat%my_rings(irn))
            end do
          end if
        end if
        try => try%next
      end do
      if (.not.connected) then
        master%prev => ithis%my_rings(ibnd)
        master => master%prev  
      end if
9     CONTINUE
    end do
    rt_list => rt_list%next                                             ! move to next rt task
  end do
  call origin%finish_chain
  deallocate(master)
  nullify(origin, temp, this, that, ithis, ithat)
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE make_chains

!> =====================================================================
!> My separate nbor connectors. The original one in link_mod does not do
!> the job right for rt. Second one is for connecting rt_internal and
!> corresponding corrections. Both checks are deactivated by default. If
!> task becomes a BC one, it activates itself's link.
!> =====================================================================
SUBROUTINE connect_nbors(self, nbor)
  class(link_t), pointer :: self
  class(link_t), pointer :: nbor, link1, link2
! -------------------------------------------------------------------
  allocate(link1, link2)
  link1%task   => self%task
  link1%link   => self
  link1%parent => self%parent
  link1%needed   = .true.
  link1%needs_me = .true.
  call nbor%add_nbor_by_rank(link1)
  link2%task   => nbor%task
  link2%link   => nbor
  link2%parent => nbor%parent
  link2%needed   = .true.
  link2%needs_me = .true.
  call self%add_nbor_by_rank(link2)
END SUBROUTINE connect_nbors

END MODULE radiative_transfer_mod
