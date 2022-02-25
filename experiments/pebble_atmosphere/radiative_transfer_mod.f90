!===============================================================================
!> Module for setting up the radiative transfer, and the relations between the
!> RT tasks and ray systems.
!===============================================================================
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
  USE rt_boundaries_mod
  USE rt_components_mod
  USE scaling_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod
  implicit none
  private
  PUBLIC rt_setup, connect_nbors
CONTAINS

!===============================================================================
!> The only procedure callable from the outside
!===============================================================================
SUBROUTINE rt_setup (self, level_min)
  class(task_list_t)     :: self
  integer, optional:: level_min                                         ! lowest level, at which RT is started
  !.....................................................................
  class(list_t), pointer :: rt_task_list                                ! local task list
  class(list_t), pointer :: list
  class(link_t), pointer :: t_list, link, rt_link
  class(task_t), pointer :: task
#ifdef __GFORTRAN__
  class(mesh_t), pointer:: mesh
#endif
  class(experiment_t), pointer :: mhd_patch
  class(rt_t),   pointer :: rt_task
  real(8) :: eps = 1d-3
  real(8):: lc(3) = +huge(0d0)        ! lower corner
  real(8):: uc(3) = -huge(0d0)        ! upper corner
  real(8):: llc(3) = +huge(0d0)       ! lower left corner
  real(8):: urc(3) = -huge(0d0)       ! upper right corner
  integer :: nd, lmin, ii
  integer :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin('radiative_transfer::rt_setup', itimer=itimer)
  associate (mhd => self%head%task)
  SELECT TYPE(mhd)
  class is (experiment_t)
    mhd_patch => mhd
  class default
    error stop '(1) in rt_setup: unknown task type'
  end SELECT
  end associate

  !> -------------------------------------------------------------------
  !> Init the ray_geometry objet.
  !> -------------------------------------------------------------------
#ifndef __GFORTRAN__
  associate (mesh => mhd_patch%mesh(1))
#else
  mesh => mhd_patch%mesh(1)
#endif
  select type(mesh)
    type is(Cartesian_mesh)
      !call ray_geometry%init(geometry=1, &
      !n=mhd_patch%n, ng=mhd_patch%ng, cell_dims=mhd_patch%ds)
    type is(spherical_mesh)
      !call ray_geometry%init(geometry=2, &
      !n=mhd_patch%n,  ng=mhd_patch%ng, cell_dims=mhd_patch%ds)
    type is(cylindrical_mesh)
      !call ray_geometry%init(geometry=3, &
      !n=mhd_patch%n,  ng=mhd_patch%ng, cell_dims=mhd_patch%ds)
    class default
      error stop '(2) in rt_setup: unknown mesh type'
  end select
#ifndef __GFORTRAN__
  end associate
#else
  nullify(mesh)
#endif
  if (ray_geometry%head%ray_list%directionality(1) == 1) then
    nd = 2
  else
    nd = 1
  end if

  write (io_unit%log,*) 'Ray geometry: n_ray_dirs:', ray_geometry%head%n_ray_dirs
  !> -------------------------------------------------------------------
  !> Initialize the RT_boundary object
  !> -------------------------------------------------------------------
  call self%init_bdries
  if (present(level_min)) then
    t_list => self%head
    lmin = level_min
    do while (associated(t_list))
      if (t_list%task%level >= level_min) then
        lc = min(lc, t_list%task%position)
        uc = max(uc, t_list%task%position)
        llc = min(llc, t_list%task%position-0.5000000001_8*t_list%task%size)
        urc = max(urc, t_list%task%position+0.5000000001_8*t_list%task%size)
      end if
      t_list => t_list%next
    end do
  else
    lmin = 0
    lc = self%lc; uc = self%uc
    llc = self%llc; urc = self%urc
  end if
  call rt_boundaries%init(llc, urc)
  !> -------------------------------------------------------------------
  !> Go through the task list, select the rmhd_t tasks and initialize the
  !> radiative transfer patches for them.
  !> -------------------------------------------------------------------
  allocate(rt_task_list)
  call rt_task_list%init ('rt task list')
  t_list => self%head
  do while (associated(t_list))
    !> -----------------------------------------------------------------
    !> Get the right type of the task to pick up other properties (rmhd_t)
    !> -----------------------------------------------------------------
    if (t_list%task%level >=lmin) then
      associate (mhd => t_list%task)
      SELECT TYPE(mhd)
      class is (experiment_t)
        allocate (rt_task, link)                                        ! the internal RT task and its link
        mhd_patch => mhd
        !> -------------------------------------------------------------
        !> Construct the RT task and add it to the task list,
        !> giving it a nbor list consisting of only the MHD task.
        !> -------------------------------------------------------------
        link%task         => rt_task                                    ! connect to RT task
        rt_task%link => link                                            ! save it's link for convenience
        !call connect_nbors(link,t_list)                                ! connect rt task and rmhd task
        mhd%rt => rt_task                                               ! save link to rt task
        rt_task%mhd_link => t_list                                      ! save link to mhd task
        rt_task%position = mhd_patch%position                           ! shadow the position
        rt_task%n        = mhd_patch%n                                  ! active innre cells
        rt_task%gn       = mhd_patch%gn                                 ! guard zones
        rt_task%size     = mhd_patch%size                               ! shadow the size
        rt_task%box      = mhd_patch%box                                ! same box
        rt_task%level    = mhd_patch%level                              ! same level
        rt_task%ds       = mhd_patch%ds                                 ! same cell size
        rt_task%time     = mhd_patch%time                               ! same cell size
        rt_task%iout     = mhd_patch%iout                               ! same cell size
        rt_task%out_next = mhd_patch%out_next                           ! same cell size
        call rt_task%init                                               ! initialise
        !> -------------------------------------------------------------
        !> Set up the boundary bits
        !> -------------------------------------------------------------
        if (.not.rt_task%periodic(1)) then
          if (rt_task%position(1) <= lc(1)+eps) call rt_task%rt_boundaries%set (bits%xl)
          if (rt_task%position(1) >= uc(1)-eps) call rt_task%rt_boundaries%set (bits%xu)
        end if
        if (.not.rt_task%periodic(2)) then
          if (rt_task%position(2) <= lc(2)+eps) call rt_task%rt_boundaries%set (bits%yl)
          if (rt_task%position(2) >= uc(2)-eps) call rt_task%rt_boundaries%set (bits%yu)
        end if
        if (.not.rt_task%periodic(3)) then
          if (rt_task%position(3) <= lc(3)+eps) call rt_task%rt_boundaries%set (bits%zl)
          if (rt_task%position(3) >= uc(3)-eps) call rt_task%rt_boundaries%set (bits%zu)
        end if
        call rt_task%rt_boundaries%init (rt_task%gn, ray_geometry%head%n_ray_dirs*nd, rt_task%mesh)
        call rt_task_list%append_link(link)                             ! insert rt_task into the rt task list
      end SELECT
      end associate
    end if
    t_list => t_list%next                                               ! go to the next mhd task
  end do

  call find_rt_nbors(rt_task_list)                                      ! connect nboring rt tasks to be nbors
  !> -------------------------------------------------------------------
  !> Depending whether the RT chains are used (default, recommended) or
  !> not, call the right subroutine to initialize them.
  !> -------------------------------------------------------------------
  if (ray_geometry%head%do_chains) call make_chains (rt_task_list)      ! inter-connect all the chains

  !> -------------------------------------------------------------------
  !> TODO:  for mesh refinement and mpi, add a call to find upstreams.
  !> -------------------------------------------------------------------

  !> -------------------------------------------------------------------
  !> append the rt tasks to the task list.
  !> -------------------------------------------------------------------
  do ii = 1,ray_geometry%head%n_warm_up
    link => rt_task_list%head
    print *, 'Warming up the RT communication...', ii
    do while (associated(link))
      associate (mhd => link%task)
      SELECT TYPE(mhd)
      class is (rt_t)
        rt_task => mhd
      class default
        error stop '(2) in rt_setup'
      end SELECT
      end associate
      rt_link => link%next
      rt_task%dry = .true.
      call rt_task%update 
      rt_task%dry = .false.
      call rt_task%finish()
      !print *, ii,'dry-update for rt task', rt_task%id, rt_task%time
      link => rt_link
    end do
  end do
  !call self%print
  call trace%end (itimer)
  if (io%master) then
    print *, '=================================='
    print *, 'RT SETUP IS COMPLETE'
    print *, '=================================='
  end if
END SUBROUTINE rt_setup


!=======================================================================
!> Finalize the bundle chains by connecting them to the appropriate
!> upstream and downstream.
!> INPUT: rt_task_list
!> OUTPUT: none
!> RESULT: rt chains are connected to each other. Master chains are
!>   identified and all chain inter-connections are done.
!> SUBROUTINES, CALLED FROM WITHIN:
!>   rt_mod::rt_chain_t%finish_chain
!>   rt_mod::rt_chain_t%connect_to_upstream_diff_level
!>   rt_mod::rt_chain_t%link_connects_to
!> FLOW:
!> 1) Subroutine loops over the rt_task list, takes one patch and loops
!>   over all it's chain links.
!> 2) The appropriate chain link is tested against being at a boundary.
!>   If that is the case, it is marked appropriately.
!> 3) Initially the link is marked as being the "master" link, meaning
!>   the chain reaction starts from it.
!> 4) The search for an upstream link begins. A loop over all patches in
!>   the rt_task list is done until the right chain is found. If an
!>   upstream is found, the "master" flag is removed and chains are connected.
!> 5) This new connection is reported through the chain (if there are
!>   more links connected downstream), so every link would know, that a
!>   new master has taken the throne. The loop (4) ends.
!> 6) If the link encounters an upstream with different refinement level,
!>   it doesn't lose it's master title, but connects the upstream to
!>   a separate link (upstream, which is an array of pointers). In such
!>   case, a subroutine to connect the appropriate intensity positions
!>   is called
!> 6) When the loop through rt tasks is done, and the link still hasn't
!>   lost it's title as a master, it is put to the master list. For sim-
!>   plicity, chain_link_t structure is used here, since links are not
!>   tasks. In order not to mess up the link structure, the connections
!>   are done backwards (link%prev) - master links do not have %prev
!>   companions, so this is safe. Upon finalizing the chains, the %prev
!>   connection is broken.
!> 7) When the loop (1) ends, another loop, over the master chains (6) is
!>   done. Starting from each master link, all links are finalized to point
!>   their incoming intensity locations to where they are supposed to
!>   point, see rt_mod::rt_chain_t%finish_chain for more details.
!> 8) After (7) is complete, this subroutine has done it's job.
!=======================================================================
SUBROUTINE make_chains (self)
  class(list_t) :: self
  !.....................................................................
  class(link_t),       pointer :: rt_list
  class(link_t),       pointer :: rt, try
  class(chain_link_t), pointer :: master, origin, temp, my_link, nb_link
  class(task_t),       pointer :: this, that
  class(rt_t),         pointer :: ithis, ithat
  logical :: connects = .false.
  logical :: upstream = .false.
  logical :: connected = .false.
  logical :: is_bndry
  integer:: ip, ipn, n_links, n_masters, n_slaves, n_conn,n_bn
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  print *, 'making chains...'
  !$omp critical (list_cr)
  allocate(master)                                                      ! temporary master chain list
  n_links=0
  n_masters=0
  n_slaves=0
  n_conn=0
  n_bn=0
  origin => master                                                      ! this will be the head of the temporary master chain
  rt_list => self%head
  do while (associated(rt_list))                                        ! loop through rt_task_list

    this => rt_list%task
    SELECT TYPE (this)
    class is (rt_t)                                                     ! select the right task
      ithis => this
    end SELECT

    do ip = 1, ithis%n_packet                                           ! loop over chain links
      is_bndry=.false.                                                  ! reset the boundary flag
      connected = .false.                                               ! no connections have been made
      n_links=n_links+1
      my_link => ithis%my_links(ip)
      my_link%master = .true.                                           ! flag the link as being the master for now

      !> ---------------------------------------------------------------
      !> Check if the link is adjacent to an outer on an inner boundary
      !> ---------------------------------------------------------------
      if (ithis%has_OB(my_link)) then                                   ! if this link is at the outer boundary
        my_link%BC(1) = .true.                                          ! then flag it to have an outer boundary condition
        call ithis%map_bndry(my_link)
        is_bndry=.true.                                                 ! since it is a boundary link, it can't have upstreams or other masters
        n_bn=n_bn+1
      end if

      !> ---------------------------------------------------------------
      !> If it is not adjacent to a boundary, find a master link for it
      !> ---------------------------------------------------------------
      if(.not.is_bndry) then                                            ! if it is not a boundary, found a master on an upstream
        !try => rt_list%nbor
      tryloop: do while (associated(try))                               ! loop through rt_task_list (4)
          that => try%task
          SELECT TYPE (that)
          class is (rt_t)
            ithat => that
          end SELECT

          if (ithis%id /= ithat%id) then                                ! not the same patch
            if (ithis%overlaps(ithat)) then                             ! patches overlap
              ! loop over ithat links
              nblink_loop: do ipn=1, ithat%n_packet                     ! loop over nbor links
                connects = .false.                                      ! can links connect?
                upstream = .false.                                      ! nbor is different lvl upstreamer?
                nb_link => ithat%my_links(ipn)
                call ithis%link_connects_to(ithat,my_link,nb_link, &
                connects, upstream)
                !> -----------------------------------------------------
                !> if an upstream link is found, do the following
                !> -----------------------------------------------------
                if (connects.and..not.upstream) then                    ! the links connect
                  my_link%master = .false.                              ! i am no longer the master
                  my_link%prev => nb_link                               ! this is my superior
                  nb_link%next => my_link                               ! this is my inferior

                  if(.not.associated(my_link%tail)) then
                    nb_link%tail => my_link                             ! connect the tails
                  else
                    nb_link%tail => my_link%tail
                  end if
                  !> ---------------------------------------------------
                  !> go through the chain (if there is any) of the newly
                  !> connected link and repoint their head pointers
                  !> ---------------------------------------------------
                  temp => my_link                                       ! inform the rest of the chain, that the new master has arrived
                  do while(associated(temp))
                    temp%head => nb_link
                    temp => temp%next
                  end do
                  connected = .true.                                    ! a connection has been made
                  n_conn=n_conn+1
                  exit tryloop                                          ! the loop over nbor links can be stopped
                end if

                if (upstream) then
                print *, 'found upstream, my id, up id', ithis%id, ithat%id
                 my_link%BC(3) = .true.                                 ! if the upstream link is of different level or is a virtual patch
                end if
              end do nblink_loop
            end if
          end if
          try => try%next
        end do tryloop
      end if
      !> ---------------------------------------------------------------
      !> if no connections have been done, that means this link is still
      !> a  master link, therefore it has to be added to the master ting
      !> chain.
      !> ---------------------------------------------------------------
      if (.not.connected) then
          n_masters = n_masters+1
          master%prev => my_link
          master => master%prev
      else
        n_slaves = n_slaves+1
      end if
    end do
    rt_list => rt_list%next                                             ! move to next rt task
  end do

  if (io%verbose>-1) print *, 'n_links',n_links,'n_masters',n_masters, &
        'n_slaves',n_slaves,'n_connected',n_conn,'n_boundary',n_bn
  print *, 'made chains'
  !> -------------------------------------------------------------------
  !> finalize the chains by mapping their rays to the correct upstream
  !> positions.
  !> -------------------------------------------------------------------
  call origin%finish_chain
  print *, 'finished chains'
  deallocate(origin)
  nullify(temp, this, that, ithis, ithat, my_link, nb_link)
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE make_chains

!=======================================================================
!> My separate nbor connectors. The original one in link_mod does not do
!> the job right for rt.
!=======================================================================
SUBROUTINE connect_nbors(self, nbor)
  class(nbor_t), pointer :: self
  class(nbor_t), pointer :: nbor
  class(link_t), pointer :: link1, link2
! -------------------------------------------------------------------
  !allocate(link1, link2)
  !link1%task   => self%task
  !link1%link   => self
  !link1%parent => self%parent
  !nbor1%needed   = .true.
  !nbor1%needs_me = .true.
  !call nbor%add_nbor_by_rank(nbor%nbor,link1)
  !link2%task   => nbor%task
  !link2%link   => nbor
  !link2%parent => nbor%parent
  !nbor2%needed   = .true.
  !nbor2%needs_me = .true.
  !call self%add_nbor_by_rank(self%nbor,link2)
END SUBROUTINE connect_nbors

!=======================================================================
!> Go through the rt_task list, find the nboring rt tasks ant make them
!> nbors. Original init_nbors subroutine is not suitable, because all
!> rt_patches already have a nbor - the rmhd_patch.
!=======================================================================
SUBROUTINE find_rt_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link1, link2
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: p1, p2
  integer, save:: itimer=0
  call trace%begin('radiative_transfer_mod%find_rt_nbors',itimer=itimer)
  link1 => self%head                                                    ! start on the self list
  do while (associated(link1))
    link2 => link1                                                      ! start searching throuch the chain
    do while (associated(link2))
      if (.not.associated(link1%task,link2%task)) then
        p1 => task2patch(link1%task)
        p2 => task2patch(link2%task)
        if (p1%overlaps(p2)) then
          !call connect_nbors(link1,link2)
        end if
      end if
      link2 => link2%next                                               ! advance in the inner chain
    end do
    link1 => link1%next                                                 ! advance in the outer chain
  end do
  call trace_end(itimer)
END SUBROUTINE find_rt_nbors

END MODULE radiative_transfer_mod
