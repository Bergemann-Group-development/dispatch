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
  USE rmhd_mod
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
  class(link_t),      pointer :: t_list, link, rt_link
  class(task_t),      pointer :: task
  class(mesh_t),      pointer :: mesh
  class(rmhd_t),      pointer :: mhd_patch
  class(rt_t),        pointer :: rt_task
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer::rt_setup')
  associate (mhd => self%head%task)
  SELECT TYPE(mhd)
  class is (rmhd_t)
    mhd_patch => mhd
  class default
    error stop '(1) in rt_setup'
  end SELECT
  end associate

  mesh => mhd_patch%mesh(1)
  select type(mesh)
    type is(Cartesian_mesh)
      call ray_geometry%init(geometry=1, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds)  
    type is(spherical_mesh)
      call ray_geometry%init(geometry=2, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds) 
    type is(cylindrical_mesh)
      call ray_geometry%init(geometry=3, li=mhd_patch%li,ui=mhd_patch%ui, cell_dims=mhd_patch%ds)   
  end select
  nullify(mesh)

  write (io_unit%log,*) 'Ray geometry: n_ray_dirs:', ray_geometry%head%n_ray_dirs

  call rt_boundary%init(lbox=mhd_patch%box-self%size,ubox=mhd_patch%box)   !!!!!!! FIXME !!!!!!

  allocate(rt_task_list)
  call rt_task_list%list_t%init ('rt task list')
  t_list => self%head
  do while (associated(t_list))
    !-------------------------------------------------------------------
    ! Get the right type of the task to pick up other properties (mhd_t)
    !-------------------------------------------------------------------
    associate (mhd => t_list%task)
    SELECT TYPE(mhd)
    class is (rmhd_t)
      allocate (rt_task, link)                                          ! the internal RT task and its link
      mhd_patch => mhd
      !-----------------------------------------------------------------
      ! Construct the RT task and add it to the task list, 
      ! giving it a nbor list consisting of only the MHD task.
      !-----------------------------------------------------------------
      link%task         => rt_task                                      ! connect to RT task
      rt_task%link => link
      call connect_nbors(link,t_list)
      mhd%rt => rt_task
      rt_task%mhd_link => t_list
      rt_task%position = mhd_patch%position
      rt_task%dim      = mhd_patch%mesh%gn                              ! RT dimensions 
      rt_task%n        = mhd_patch%n
      rt_task%gn       = mhd_patch%gn
      rt_task%size     = mhd_patch%size
      call rt_task%init                                                 ! initialise
      rt_task%box      = mhd_patch%box
      rt_task%level    = mhd_patch%level                                                
      rt_task%ds       = mhd_patch%ds
      rt_task%li       = mhd_patch%li
      rt_task%ui       = mhd_patch%ui
      call rt_task_list%append_link(link)                               ! insert rt_task into the rt task list
    end SELECT
    end associate
    t_list => t_list%next                                               ! go to the next mhd task
  end do
  call find_rt_nbors(rt_task_list)
  if (ray_geometry%head%do_chains) call make_chains (rt_task_list)      ! inter-connect all the chains
  link => rt_task_list%head
  do while (associated(link))
    call self%append_link (link)
    link => link%next
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
!>   identified and all chain inter-connections are done.
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
  logical :: is_bndry
  integer:: ibnd, irn, n_rings, n_masters, n_slaves, n_conn,n_bn
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  allocate(master)                                                      ! temporary master chain list
  n_rings=0
  n_masters=0
  n_slaves=0
  n_conn=0
  n_bn=0
  origin => master
  rt_list => self%head
  do while (associated(rt_list))                                        ! loop through rt_task_list
    this => rt_list%task
    SELECT TYPE (this)
    class is (rt_t)                                                     ! select the right task
      ithis => this
    end SELECT
    do ibnd = 1, ithis%n_bundle                                         ! loop over chain rings
      is_bndry=.false.
      n_rings=n_rings+1
      connected = .false.                                               ! no connections have been made
      ithis%my_rings(ibnd)%master = .true.                              ! flag the ring as being the master
      
      if (ithis%has_OB(ithis%my_rings(ibnd))) then                      ! if this ring is at the outer boundary
        ithis%my_rings(ibnd)%BC(1) = .true.                             ! then flag it to have an outer boundary condition
        is_bndry=.true.                                                 ! since it is a boundary ring, it can't have upstreams or other masters 
        n_bn=n_bn+1                                        
      end if
        ! TODO: inner boundaries.
        !if (ithis%has_IB(ithis%my_rings(ibnd))) then                    ! if this ring is at the inner boundary
        !  error stop "radiative_transfer_mod:: make_chains (1)"         ! not implemented yet!
        !  ithis%my_rings(ibnd)%BC(2) = .true.                           ! then flag it to have an inner boundary condition
        !  master%prev => ithis%my_rings(ibnd)
        !  goto 9                                                        ! since it is a boundary ring, it can't have upstreams or other masters
        !end if

      if(.not.is_bndry) then                                            ! if it is not a boundary, found a master on an upstream
        try => rt_list%nbor
      tryloop: do while (associated(try))                               ! loop through rt_task_list (4)
          that => try%task
          SELECT TYPE (that)
          class is (rt_t)
            ithat => that
          end SELECT
        
          if (ithis%id /= ithat%id) then                                  ! not the same patch
            if (ithis%overlaps(ithat)) then                               ! patches overlap
            
              ! loop over ithat rings
              nbring_loop: do irn=1, ithat%n_bundle                       ! loop over nbor rings
                connects = .false.                                        ! can rings connect?
                is_upstream = .false.                                     ! nbor is different lvl upstreamer?
                call ithis%ring_connects_to(ithat,ithis%my_rings(ibnd),ithat%my_rings(irn),connects, is_upstream)
                !..................
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
                  n_conn=n_conn+1
                  exit tryloop
                end if
                !...................
                if (is_upstream) call &
                  ithis%my_rings(ibnd)%connect_to_upstream_diff_level(ithat%my_rings(irn))
              end do nbring_loop
            end if
          end if
          try => try%next
        end do tryloop
      end if
      if (.not.connected) then
        n_masters = n_masters+1
        master%prev => ithis%my_rings(ibnd)
        master => master%prev  
      else
        n_slaves = n_slaves+1
      end if
    end do
    rt_list => rt_list%next                                             ! move to next rt task
  end do
  if (io%verbose>1) print *, 'n_rings',n_rings,'n_masters',n_masters, &
        'n_slaves',n_slaves,'n_connected',n_conn,'n_boundary',n_bn

  call origin%finish_chain                                              ! finalize the chain by giving them right upstream positions
  deallocate(origin)
  nullify(temp, this, that, ithis, ithat)
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE make_chains

!> =====================================================================
!> My separate nbor connectors. The original one in link_mod does not do
!> the job right for rt.
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

SUBROUTINE find_rt_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link1, link2
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: p1, p2
  link1 => self%head                                                    ! start on the self list
  do while (associated(link1))
    link2 => link1                                                  ! start searching throuch the chain
    do while (associated(link2))
      if (.not.associated(link1%task,link2%task)) then
        p1 => task2patch(link1%task)
        p2 => task2patch(link2%task)
        if (p1%overlaps(p2)) then
          call connect_nbors(link1,link2)
        end if
      end if
      link2 => link2%next                                               ! advance in the inner chain
    end do
    link1 => link1%next                                                 ! advance in the outer chain
  end do
  call trace_end
END SUBROUTINE find_rt_nbors

END MODULE radiative_transfer_mod
