MODULE interface_mod
  USE kinds_mod
  USE trace_mod
  USE io_mod
  USE dll_mod
  USE link_mod
  USE patch_mod
  USE task_mod
  USE bdry_interface_mod
  implicit none
  private
  !=============================================================================
  !> interface_t is a separate module. It's object should be a part patch_t
  !>
  !> NOTES: 
  !> 1) in AMR a child patch may have more than one interface with parent!!!
  !> 2) prevent booking same node multiple times!
  !>
  !=============================================================================

  type, public:: fields_halo_t
    class(dll_t), dimension(:,:,:), pointer:: field
  contains
    procedure:: collect_from_nbor
    procedure, private:: deposit_for_nbors4
    procedure, private:: deposit_for_nbors5
    generic, public:: deposit_for_nbors => deposit_for_nbors4, deposit_for_nbors5
    procedure:: init => init_fields_halo
    procedure:: get_halo_coord
    procedure:: get_fields_hub
  end type


  type, public:: interface_list_t
    class(patch_t), pointer:: patch                                 ! link to master patch
  contains
    procedure:: add_quantity
    procedure:: create_nbor_node
    procedure:: get_halo
    procedure:: init
    procedure:: remove_quantity
    procedure, private:: update4
    procedure, private:: update5
    generic, public:: update => update4, update5
    !procedure:: check_accounts -- check if any accounts should be opened/closed
    !procedure:: open_account  ! add_nbor
    !procedure:: close_account ! remove_nbor
  end type
  logical, save:: debug=.true.
  logical, save:: detailed_timer=.false.
  logical, save:: flux_correction = .false.
  integer, save:: verbose=0
  logical, save:: on=.false.
  CONTAINS

  !=============================================================================
  !> Definitions:
  !> p1 - master patch
  !> p2..px - nbor patches
  !> 
  !> We first allocate an area cube for p1. The actual area nodes are put in place
  !> later, when the interface quantities are requested.
  !> Create a dll of interface areas between nbors. Later, this can be maintained
  !> separately.
  !> for the 
  !> For each nbor it identifies, what area (face,edge,corner) the patch has an
  !> interface with. The according item in the 3x3 area cube is allocated
  !=============================================================================
  SUBROUTINE init (self, link)
    class(interface_list_t):: self
    class(link_t), pointer :: link
    ! --------------------------------------------------------------------------
    class(patch_t), pointer:: patch
    logical, save:: first_time=.true.
    integer:: iostat
    namelist /interface_params/ on, flux_correction
    ! ..........................................................................
    call io%lock%set ('interface_mod%init')
    if (first_time) then
      rewind (io%input)
      read (io%input, interface_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, interface_params)
      end if
      call io%namelist_warning ('interface_params', iostat)
      write (io%namelist, interface_params)
      first_time = .false.
    end if
    call io%lock%unset ('interface_mod%init')
    patch => task2patch(link%task)
    self%patch => patch
    if (on) call patch%bdry_interface%init('bndry_interface_list')
    verbose = io%verbose
  END SUBROUTINE init

  !=============================================================================
  !> Add flux or field areas to interface nodes.
  !> Called from mhd_t%init
  !> This is done in two steps:
  !> 1) for a given variable iv, it first creates a halo
  !> for p1 and then 
  !> 2) continues to open node accounts for all the nbors. It 
  !> goes through the nbors list, check nbor location in the halo and creates the
  !> account, if it is needed.
  !=============================================================================
  SUBROUTINE add_quantity (self, label)
    class(interface_list_t):: self
    character(len=*):: label
    ! --------------------------------------------------------------------------
    integer:: id
    logical, dimension(:,:,:), pointer:: halo
    class(nbor_dll_t), pointer:: nbor_dll
    class(dll_node_t), pointer:: node, nb_node
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(patch_t), pointer:: patch,nbpatch
    class(*), pointer:: car
    ! ..........................................................................
    integer:: ix, iy, iz
    integer:: rng(3)          ! interval range
    logical,dimension(3,3,3)::halo_loc
    ! --------------------------------------------------------------------------
    if (.not.on) return
    call trace%begin ('interface_list_t::add_quantity')
    ! --------------------------------------------------------------------------
    ! obtain the halo configuration and corresponding id.
    ! --------------------------------------------------------------------------
    call self%get_halo(label, id, halo)
    ! --------------------------------------------------------------------------
    ! Add a new halo node
    ! --------------------------------------------------------------------------
    patch => self%patch
    call patch%bdry_interface%add_halo_node(id,label,node)
    deallocate(halo)
  END SUBROUTINE add_quantity

  SUBROUTINE update4(self, field, label)
    class(interface_list_t):: self
    real(kind=var_k), dimension(:,:,:,:), pointer:: field
    character(len=*):: label
    ! --------------------------------------------------------------------------
    class(patch_t), pointer:: patch, nbpatch
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(dll_node_t), pointer:: node, nb_node, node2, my_node
    class(dll_t), pointer:: field_dll, nb_field_dll
    class(dll_node_t), pointer:: field_node, nb_field_node
    class(fields_hub_t), pointer:: fields_hub
    class(nbor_dll_t), pointer:: nbor_dll, nb_nbor_dll
    integer:: id, crd(3), inode, l(3),u(3), iv, dir, it, nb_it
    class(fields_halo_t), pointer:: fields_halo, nb_fields_halo
    logical, dimension(:,:,:), pointer:: halo
    class(*), pointer:: car
    integer, save:: itimer=0
    logical:: nbor_operations
    ! ..........................................................................
    ! ..........................................................................
    if (.not.on) return
    call trace%begin('interface_list_t%update4', itimer=itimer)
    call self%get_halo(label, id, halo)
    !
    patch => self%patch
    !it = patch%it
    it = 1
    nb_it=1
    ! --------------------------------------------------------------------------
    ! Divide by dtime
    ! --------------------------------------------------------------------------
    !field = field/patch%dtime
    ! --------------------------------------------------------------------------
    ! get the node for this field
    ! --------------------------------------------------------------------------
    nullify(node)
    call patch%bdry_interface%find_id(id,node)
    call patch%bdry_interface%lock%set('interface_update')
    ! --------------------------------------------------------------------------
    ! get the nbor_dll
    ! --------------------------------------------------------------------------
    call patch%bdry_interface%get_nbor_dll(node,nbor_dll,label='update4')
    ! --------------------------------------------------------------------------
    !> Go through nbors list and update the fields. Nbor might be not present
    !> due to several reasons:
    !> 1) it is the first time
    !> 2) new nbor has arrived, e.g. in AMR or moving patches
    !> if that is the case, create a new nbor node
    ! --------------------------------------------------------------------------
    nbor => patch%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      nbpatch=> task2patch(nbtask)
      nullify(nb_node)
      call nbor_dll%find_id(nbor%task%id,nb_node)
      if (.not.associated(nb_node)) then
        ! ----------------------------------------------------------------------
        ! Create new nbor node
        ! ----------------------------------------------------------------------
        if (nbor%download.and.needed_for_field(patch,nbpatch,halo)) then
          call self%create_nbor_node(nbor_dll,nbor, halo, fields_halo,patch, nbpatch, label)
        else
          nbor => nbor%next
          cycle
        end if
      else
        ! ----------------------------------------------------------------------
        ! Select nbor node
        ! ----------------------------------------------------------------------
        car => nb_node%car
        select type(car)
        class is (fields_halo_t)
          fields_halo => car
        end select
      end if
      !
      ! ------------------------------------------------------------------------
      ! Do the operations with nbors only if nbors are available
      ! ------------------------------------------------------------------------
      nbor_operations=.true.
      ! ------------------------------------------------------------------------
      ! Find my node at nbor
      ! ------------------------------------------------------------------------
      nullify(node2)
      call nbpatch%bdry_interface%find_id(id,node2)
      if (associated(node2)) then
        ! ----------------------------------------------------------------------
        ! get the nbor nbor_dll
        ! ----------------------------------------------------------------------
        call nbpatch%bdry_interface%get_nbor_dll(node2,nb_nbor_dll,label='update4-3')
        nullify(my_node)
        if (nb_nbor_dll%updated) then
          call nb_nbor_dll%find_id(patch%id,my_node)
          ! ----------------------------------------------------------------------
          ! Select nbor node
          ! ----------------------------------------------------------------------
          if (associated(my_node)) then
            car => my_node%car
            select type(car)
            class is (fields_halo_t)
              nb_fields_halo => car
              call nb_fields_halo%get_halo_coord(nbpatch,patch,crd)
              nb_field_dll => nb_fields_halo%field(crd(1),crd(2),crd(3))
            end select
          else
            nbor_operations = .false.
          end if
        else
          nbor_operations=.false.
        end if

      else
        nbor_operations=.false.
      end if
      !
      ! ------------------------------------------------------------------------
      ! Update fields
      ! ------------------------------------------------------------------------
      call fields_halo%get_halo_coord(patch,nbpatch,crd)
      field_dll =>fields_halo%field(crd(1),crd(2),crd(3))
      !            !
      ! ------------------------------------------------------------------------
      ! We first collect the field data from nbors, if that is applicable
      ! ------------------------------------------------------------------------
      if (nbor_operations) then
        call nb_field_dll%lock%set('nb_field_update')
        call fields_halo%collect_from_nbor(field_dll,nb_field_dll, it, nb_it)
        call nb_field_dll%lock%unset('nb_field_update')
      end if
      ! ------------------------------------------------------------------------
      ! And now we update the local fields
      ! ------------------------------------------------------------------------
      call field_dll%lock%set('interface_update')
      field_node => field_dll%head
      do while (associated(field_node))
        call fields_halo%get_fields_hub(field_node,fields_hub,label='update4-2')
        l  = fields_hub%my_area(it)%l
        u  = fields_hub%my_area(it)%u
        iv = fields_hub%iv
        fields_hub%f1 => field(l(1):u(1),l(2):u(2),l(3):u(3),iv)
        call fields_hub%update(patch%time, patch%dtime, it, iv==1)
        field_node => field_node%next
      end do
      ! ------------------------------------------------------------------------
      ! Lastly, we update the data for nbors to pick up
      ! ------------------------------------------------------------------------
      call fields_halo%deposit_for_nbors(field_dll,field, it, patch%time, patch%dtime)
      ! ------------------------------------------------------------------------
      ! Unlock the dll
      ! ------------------------------------------------------------------------
      call field_dll%lock%unset('interface_update')
      ! ------------------------------------------------------------------------
      ! Go to the next nbor
      ! ------------------------------------------------------------------------
      nbor => nbor%next
    end do
    nbor_dll%updated=.true.

    call patch%bdry_interface%lock%unset('interface_update')

    ! --------------------------------------------------------------------------
    ! Multiply by dtime
    ! --------------------------------------------------------------------------
    !field = field*patch%dtime
    ! --------------------------------------------------------------------------
    ! clean up
    ! --------------------------------------------------------------------------
    deallocate(halo)
    call trace%end(itimer)
    
  END SUBROUTINE update4

  SUBROUTINE update5(self, field, label, all)
    class(interface_list_t):: self
    real(kind=var_k), dimension(:,:,:,:,:), pointer:: field
    character(len=*):: label
    logical, optional:: all
    ! --------------------------------------------------------------------------
    class(patch_t), pointer:: patch, nbpatch
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    class(dll_node_t), pointer:: node, nb_node, node2, my_node
    class(dll_t), pointer:: field_dll, nb_field_dll
    class(dll_node_t), pointer:: field_node, nb_field_node
    class(fields_hub_t), pointer:: fields_hub
    class(nbor_dll_t), pointer:: nbor_dll, nb_nbor_dll
    integer:: id, crd(3), inode, l(3),u(3), iv, dir, it, nb_it
    class(fields_halo_t), pointer:: fields_halo, nb_fields_halo
    logical, dimension(:,:,:), pointer:: halo
    class(*), pointer:: car
    integer, save:: itimer=0
    logical:: nbor_operations
    ! ..........................................................................
    if (.not.on) return
    call trace%begin('interface_list_t%update5', itimer=itimer)
    !
    call self%get_halo(label,id, halo)
    !
    patch => self%patch
    !it = patch%it
    it = 1
    nb_it=1
    ! --------------------------------------------------------------------------
    ! Divide by dtime
    ! --------------------------------------------------------------------------
    !field = field/patch%dtime
    ! --------------------------------------------------------------------------
    ! get the node for this field
    ! --------------------------------------------------------------------------
    nullify(node)
    call patch%bdry_interface%find_id(id,node)
    call patch%bdry_interface%lock%set('interface_update')
    ! --------------------------------------------------------------------------
    ! get the nbor_dll from this node
    ! --------------------------------------------------------------------------
    call patch%bdry_interface%get_nbor_dll(node,nbor_dll,label='update5')
    ! --------------------------------------------------------------------------
    !> Go through nbors list and update the fields. Nbor might be not present
    !> due to several reasons:
    !> 1) it is the first time
    !> 2) new nbor has arrived, e.g. in AMR or moving patches
    !> if that is the case, create a new nbor node
    ! --------------------------------------------------------------------------
    nbor => patch%link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      nbpatch=> task2patch(nbtask)
      nullify(nb_node)
      call nbor_dll%find_id(nbor%task%id,nb_node)
      if (.not.associated(nb_node)) then
        ! ----------------------------------------------------------------------
        ! Create new nbor node
        ! ----------------------------------------------------------------------
        if (nbor%download.and.needed_for_field(patch,nbpatch,halo)) then
          call self%create_nbor_node(nbor_dll,nbor, halo, fields_halo, patch, nbpatch, label)
        else
          nbor => nbor%next
          cycle
        end if
      else
        ! ----------------------------------------------------------------------
        ! Select nbor node
        ! ----------------------------------------------------------------------
        car => nb_node%car
        select type(car)
        class is (fields_halo_t)
          fields_halo => car
        end select
      end if
      !
      ! ------------------------------------------------------------------------
      ! Do the operations with nbors only if nbors are available
      ! ------------------------------------------------------------------------
      nbor_operations=.true.
      ! ------------------------------------------------------------------------
      ! Find my node at nbor
      ! ------------------------------------------------------------------------
      nullify(node2)
      call nbpatch%bdry_interface%find_id(id,node2)
      if (associated(node2)) then
        ! ----------------------------------------------------------------------
        ! get the nbor nbor_dll
        ! ----------------------------------------------------------------------
        call nbpatch%bdry_interface%get_nbor_dll(node2,nb_nbor_dll,label='update5-3')
        nullify(my_node)
        if (nb_nbor_dll%updated) then
          call nb_nbor_dll%find_id(patch%id,my_node)
          ! ----------------------------------------------------------------------
          ! Select nbor node
          ! ----------------------------------------------------------------------
          if (associated(my_node)) then
            car => my_node%car
            select type(car)
            class is (fields_halo_t)
              nb_fields_halo => car
              call nb_fields_halo%get_halo_coord(nbpatch,patch,crd)
              nb_field_dll => nb_fields_halo%field(crd(1),crd(2),crd(3))
            end select
          else
            nbor_operations = .false.
          end if
        else
          nbor_operations=.false.
        end if

      else
        nbor_operations=.false.
      end if
      !
      ! ------------------------------------------------------------------------
      ! Update fields
      ! ------------------------------------------------------------------------
      call fields_halo%get_halo_coord(patch,nbpatch,crd)
      field_dll =>fields_halo%field(crd(1),crd(2),crd(3))
      !            !
      ! ------------------------------------------------------------------------
      ! We first collect the field data from nbors, if that is applicable
      ! ------------------------------------------------------------------------
      if (nbor_operations) then
        call nb_field_dll%lock%set('nb_field_update')
        call fields_halo%collect_from_nbor(field_dll,nb_field_dll, it, nb_it)
        call nb_field_dll%lock%unset('nb_field_update')
      end if

      ! ------------------------------------------------------------------------
      ! And now we update the local fields
      ! ------------------------------------------------------------------------
      call field_dll%lock%set('interface_update')
      field_node => field_dll%head
      do while (associated(field_node))
        call fields_halo%get_fields_hub(field_node,fields_hub,label='update4-2')
        l  = fields_hub%my_area(it)%l
        u  = fields_hub%my_area(it)%u
        iv = fields_hub%iv
        dir= fields_hub%dir
        fields_hub%f1 => field(l(1):u(1),l(2):u(2),l(3):u(3),iv,dir)
        call fields_hub%update(patch%time, patch%dtime, it, iv==1)
        field_node => field_node%next
      end do
      ! ------------------------------------------------------------------------
      ! Lastly, we update the data for nbors to pick up
      ! ------------------------------------------------------------------------

      call fields_halo%deposit_for_nbors(field_dll,field, it, patch%time, patch%dtime)
      ! ------------------------------------------------------------------------
      ! Unlock the dll
      ! ------------------------------------------------------------------------
      call field_dll%lock%unset('interface_update')
      ! ------------------------------------------------------------------------
      ! Go to the next nbor
      ! ------------------------------------------------------------------------
      nbor => nbor%next
    end do
    nbor_dll%updated=.true.

    call patch%bdry_interface%lock%unset('interface_update')

    ! --------------------------------------------------------------------------
    ! Multiply by dtime
    ! --------------------------------------------------------------------------
    !field = field*patch%dtime
    
    ! --------------------------------------------------------------------------
    ! clean up
    ! --------------------------------------------------------------------------
    deallocate(halo)
    call trace%end (itimer)
  END SUBROUTINE update5

  SUBROUTINE create_nbor_node(self, nbor_dll,nbor, halo, fields_halo, patch, nbpatch, label)
  class(interface_list_t):: self
  class(nbor_dll_t), pointer:: nbor_dll
  class(nbor_t), pointer:: nbor
  logical, dimension(:,:,:), pointer:: halo
  class(fields_halo_t), pointer:: fields_halo
  class(patch_t), pointer:: patch, nbpatch
  character(len=*):: label
  ! ----------------------------------------------------------------------------
  class(dll_node_t), pointer:: nb_node
  ! ............................................................................
  allocate (nb_node)
  nb_node%id = nbor%task%id
  call nbor_dll%append(nb_node)
  allocate (fields_halo)
  nb_node%car => fields_halo
  call fields_halo%init(halo, patch, nbpatch, label)
END SUBROUTINE create_nbor_node

  SUBROUTINE deposit_for_nbors5(self, field_dll, field, it, t, dt)
    class(fields_halo_t):: self
    class(dll_t), pointer:: field_dll
    real(kind=var_k), dimension(:,:,:,:,:), pointer:: field
    integer:: it
    real(8):: t, dt
    ! --------------------------------------------------------------------------
    class(dll_node_t), pointer:: node
    class(fields_hub_t), pointer:: hub
    integer:: l(3), u(3)
    ! ..........................................................................
    node => field_dll%head
    do while (associated(node))
      ! pick up local hub
      call self%get_fields_hub(node,hub,label='deposit_for_nbors5')
      call hub%deposit(t,dt,it,field)
      ! ------------------------------------------------------------------------
      ! move to next field
      ! ------------------------------------------------------------------------
      node => node%next
    end do

  END SUBROUTINE deposit_for_nbors5

  SUBROUTINE deposit_for_nbors4(self, field_dll, field, it, t, dt)
    class(fields_halo_t):: self
    class(dll_t), pointer:: field_dll
    real(kind=var_k), dimension(:,:,:,:), pointer:: field
    integer:: it
    real(8):: t, dt
    ! --------------------------------------------------------------------------
    class(dll_node_t), pointer:: node
    class(fields_hub_t), pointer:: hub
    integer:: l(3), u(3)
    ! ..........................................................................
    node => field_dll%head
    do while (associated(node))
      ! pick up local hub
      call self%get_fields_hub(node,hub,label='deposit_for_nbors4')
      call hub%deposit(t,dt,it,field)
      ! ------------------------------------------------------------------------
      ! move to next field
      ! ------------------------------------------------------------------------
      node => node%next
    end do
  END SUBROUTINE deposit_for_nbors4

  SUBROUTINE collect_from_nbor(self, field_dll, nb_field_dll, it, nb_it)
    class(fields_halo_t):: self
    class(dll_t), pointer:: field_dll, nb_field_dll
    integer:: it, nb_it
    ! --------------------------------------------------------------------------
    class(dll_node_t), pointer:: field, nb_field
    class(fields_hub_t), pointer:: hub, nb_hub
    class(*), pointer:: car
    ! ..........................................................................
    if (.not.on) return
    ! --------------------------------------------------------------------------
    ! the field_dll and nd_field_dll are dlls, containing nodes, which each 
    ! carry a fields_hub object, which has two entries: my_area and nb_area.
    ! To collect from nbors, we:
    ! 1) go through field nodes on local an nbor dlls
    ! 2) pick up corresponding hubs
    ! 3) update field%my_area%f0 with nb_field%nb_area%f0 if 
    ! --------------------------------------------------------------------------
    field => field_dll%head
    nb_field => nb_field_dll%head
    do while (associated(field))
      ! pick up local hub
      call self%get_fields_hub(field,hub,label='collect_from_nbor1')
      ! pick up nb hub
      call self%get_fields_hub(nb_field,nb_hub,label='collect_from_nbor1')
      if (hub%my_area(it)%t<=nb_hub%nb_area(nb_it)%t) then
        ! shift data to old slot
        hub%my_area(2)%t  = hub%my_area(it)%t
        hub%my_area(2)%dt = hub%my_area(it)%dt
        hub%my_area(2)%f0 = hub%my_area(it)%f0
        ! pick up new slot
        hub%my_area(it)%t  = nb_hub%nb_area(nb_it)%t
        hub%my_area(it)%dt = nb_hub%nb_area(nb_it)%dt
        hub%my_area(it)%f0 = nb_hub%nb_area(nb_it)%f0
        hub%is_ahead=.false.
        ! perhaps old slot is too outdated
        if (hub%my_area(2)%t<=nb_hub%nb_area(2)%t) then
          hub%my_area(2)%t  = nb_hub%nb_area(2)%t
          hub%my_area(2)%dt = nb_hub%nb_area(2)%dt
          hub%my_area(2)%f0 = nb_hub%nb_area(2)%f0
        end if
      else if (hub%my_area(2)%t<=nb_hub%nb_area(nb_it)%t) then
        hub%my_area(2)%t = nb_hub%nb_area(nb_it)%t
        hub%my_area(2)%dt = nb_hub%nb_area(nb_it)%dt
        hub%my_area(2)%f0= nb_hub%nb_area(nb_it)%f0
        hub%is_ahead=.false. 
      end if
      field => field%next
      nb_field => nb_field%next
    end do
  END SUBROUTINE collect_from_nbor

    SUBROUTINE get_fields_hub(self, node, field, label)
    class(fields_halo_t):: self
    class(dll_node_t),   pointer:: node
    class(fields_hub_t), pointer:: field
    character(len=*), optional:: label
    ! --------------------------------------------------------------------------
    class(*), pointer:: car
    ! ..........................................................................
    car => node%car
    select type(car)
    class is (fields_hub_t)
      field => car
    class default
      call io%abort('garbage in car at get_fields_dll, caller: '//label)
    end select
  END SUBROUTINE get_fields_hub


  SUBROUTINE get_halo(self, label, id, halo)
    class(interface_list_t):: self
    character(len=*):: label
    integer:: id
    logical, dimension(:,:,:), pointer:: halo
    ! --------------------------------------------------------------------------
    if (.not.on) return
    select case(trim(label))
    ! --------------------------------------------------------------------------
    ! hydrodynamic fluxes
    ! --------------------------------------------------------------------------
    case('hd')
      id=1
      allocate(halo(3,3,3))
      halo = .false.
      ! ------------------------------------------------------------------------
      ! only faces are needed
      ! ------------------------------------------------------------------------
      halo(2,2,1)=.true.
      halo(2,2,3)=.true.
      halo(1,2,2)=.true.
      halo(3,2,2)=.true.
      halo(2,1,2)=.true.
      halo(2,3,2)=.true.
    ! ------------------------------------------------------------------------
    ! Electromotive force
    ! --------------------------------------------------------------------------
    case('emf')
      id=2
      allocate(halo(3,3,3))
      ! ------------------------------------------------------------------------
      ! faces, edges and perhaps corners are needed
      ! ------------------------------------------------------------------------
      halo = .true.
      halo(2,2,2)=.false.
    case default
      call io%abort('interface_mod%add_quantity:: interface quantity not known')
    end select
  END SUBROUTINE get_halo

  !=============================================================================
  !> Remove a specific quantity from our accounting
  !=============================================================================
  SUBROUTINE remove_quantity (self, iv)
    class(interface_list_t):: self
    integer:: iv
    ! --------------------------------------------------------------------------
    class(dll_node_t), pointer:: node
    class(area_t), dimension(:,:,:), pointer :: area, nb_area
    class(fields_halo_t), pointer:: fields_halo
    class(*), pointer:: car
    integer:: x, y, z
    ! ..........................................................................
    if (.not.on) return
    nullify(node)
    call self%patch%bdry_interface%find_id(iv,node)
    if (associated(node)) then
      if (associated(node%car)) then
        car => node%car
        select type(car)
        class is (fields_halo_t)
          fields_halo=> car
          !area => fields_halo%my_area
          !nb_area => area_halo%nb_area
        end select
        do z=1,3
          do y=1,3
            do x=1,3
            !  call area(x,y,z)%dealloc
            !  call nb_area(x,y,z)%dealloc
            end do
          end do
        end do
      !  deallocate (area)
      !  deallocate (nb_area)
      !  deallocate(area_halo)
      end if
      call self%patch%bdry_interface%remove(node)
      deallocate(node)
    end if
  END SUBROUTINE remove_quantity

  SUBROUTINE get_halo_coord(self,patch,nbpatch,coord)
    class(fields_halo_t):: self
    class(patch_t), pointer:: patch, nbpatch
    integer:: coord(3)
    ! --------------------------------------------------------------------------
    integer:: offset(3)
    real(8):: dist(3), ds(3)
    ! ..........................................................................
    if (.not.on) return
    !
    dist = nbpatch%position - patch%position
    where (nbpatch%periodic)
      dist = modulo(dist+0.5_8*nbpatch%box,nbpatch%box) - 0.5_8*nbpatch%box
    end where
    !---------------------------------------------------------------------------
    ! For same level patch, the integer offset is what should be added to the
    ! task index to get the nbtask index.  For aligned meshes, offset=0.
    !---------------------------------------------------------------------------
    if (patch%level==nbpatch%level) then
      ds = patch%size/patch%n
      offset = nint(dist/ds)
    else
      call io%abort("interface_mod::get_halo_coord:: different level tasks not set up yet")
    end if
    where (offset==-patch%n) coord=1
    where (offset==patch%n) coord=3
    where (offset==0) coord=2
  END SUBROUTINE get_halo_coord


  SUBROUTINE init_fields_halo (self, halo, patch, nbpatch, label, all)
    class(fields_halo_t):: self
    logical, dimension(:,:,:), pointer:: halo
    class(patch_t), pointer:: patch, nbpatch
    character(len=*):: label
    logical, optional:: all
    ! ..........................................................................
    class(dll_t), pointer:: field
    class(dll_node_t), pointer:: field_node
    class(fields_hub_t), pointer:: fields_hub
    integer:: x, y, z, iv, idir, i, it, nt, nv
    integer:: offset(3), hal(3), l(3), u(3), n(3)
    real(8):: dist(3), ds(3)
    ! --------------------------------------------------------------------------
    if (.not.on) return
    nt = 2
    nullify(self%field)
    allocate(self%field(3,3,3))
    !---------------------------------------------------------------------------
    ! First take care of mapping the distance between the nbtask and task
    ! centers, periodically, so this can be ignored in the next steps.
    !---------------------------------------------------------------------------
    dist = nbpatch%position - patch%position
    where (nbpatch%periodic)
      dist = modulo(dist+0.5_8*nbpatch%box,nbpatch%box) - 0.5_8*nbpatch%box
    end where
    !---------------------------------------------------------------------------
    ! For same level patch, the integer offset is what should be added to the
    ! task index to get the nbtask index.  For aligned meshes, offset=0.
    !---------------------------------------------------------------------------
    if (patch%level==nbpatch%level) then
      ds = patch%size/patch%n
      offset = nint(dist/ds)
    else
      call io%abort("interface_mod::range:: different level tasks not set up yet")
    end if
    where (offset==-patch%n) hal=1
    where (offset==patch%n) hal=3
    where (offset==0) hal=2
    field =>self%field(hal(1),hal(2),hal(3))
    call field%init('field')
    ! --------------------------------------------------------------------------
    ! Since we know how fields should be arranged for specific solvers, it is 
    ! simplest to work explicitly here.
    ! --------------------------------------------------------------------------
    select case(trim(label))
    case('hd')
      ! ------------------------------------------------------------------------
      ! in case of HD, all fluxes are through faces, so only one item in hal
      ! will have a not-centred value
      ! ------------------------------------------------------------------------
      idir=0
      do i=1,3
        if (hal(i)/=2) idir=i
      end do
      ! ------------------------------------------------------------------------
      ! set personal index ranges; 
      ! ------------------------------------------------------------------------
      l = patch%li
      u = patch%ui+1
      n = patch%n+1
      where (hal==1) u=patch%li
      where (hal==3) l=patch%ui+1
      where (hal/=2) n = 1
      if (present(all)) then
        nv = 8
      else
        nv= 5
      end if
      ! ------------------------------------------------------------------------
      ! HD has fluxes for 5 (or 8) quantities. Only fluxes in perpendicular 
      ! direction to the face (idir direction) are
      ! relevant, thus one item per IV
      ! ------------------------------------------------------------------------
      do iv = 1,nv
        allocate(field_node)
        allocate(fields_hub)
        field_node%id = iv
        field_node%car => fields_hub
        fields_hub%iv = iv
        fields_hub%flux_correction = flux_correction
        fields_hub%dir= idir
        fields_hub%id = patch%id
        allocate(fields_hub%my_area(2))
        allocate(fields_hub%nb_area(2))
        call field%append(field_node)
        !
        ! ----------------------------------------------------------------------
        ! Loop over time slots
        ! ----------------------------------------------------------------------
        do it=1,nt
          ! --------------------------------------------------------------------
          ! Set info for personal patch
          ! --------------------------------------------------------------------
          fields_hub%my_area(it)%l = l
          fields_hub%my_area(it)%u = u
          allocate(fields_hub%my_area(it)%f0(n(1),n(2),n(3)))
          fields_hub%my_area(it)%f0 = 0.0
          fields_hub%my_area(it)%initialized = .true.
          ! --------------------------------------------------------------------
          ! Set info for nbor patch.
          ! The ranges for nbor are the range in local flux block, that nbor is 
          ! interested in. For same level patches they coincide, for different
          ! level patches, they might be in different locations. 
          ! --------------------------------------------------------------------
          fields_hub%nb_area(it)%l = l
          fields_hub%nb_area(it)%u = u
          allocate(fields_hub%nb_area(it)%f0(n(1),n(2),n(3)))
          fields_hub%nb_area(it)%f0 = 0.0
          fields_hub%nb_area(it)%initialized = .true.
        end do
      end do

    case('emf')
      ! ------------------------------------------------------------------------
      ! in case of EMF, we need two orthogonal fields for the direction
      ! ------------------------------------------------------------------------
      !idir=0
      !do i=1,3
      !  if (hal(i)/=2) idir=i
      !end do
      ! ------------------------------------------------------------------------
      ! set personal index ranges; 
      ! ------------------------------------------------------------------------
      l = patch%li
      u = patch%ui+1
      n = patch%n+1
      where (hal==1) u=patch%li
      where (hal==3) l=patch%ui+1
      where (hal/=2) n = 1
      ! ------------------------------------------------------------------------
      ! MHD has fields in 3 directions. Only orthogonal fields to particular
      ! directions are relevant, thus two items per IV
      ! ------------------------------------------------------------------------
      do iv = 1,3
       ! if (iv==idir) cycle ! skip the idir direction
        allocate(field_node)
        allocate(fields_hub)
        field_node%id = iv
        field_node%car => fields_hub
        fields_hub%iv = iv
        fields_hub%flux_correction = flux_correction
        fields_hub%dir= iv
        fields_hub%id = patch%id
        allocate(fields_hub%my_area(2))
        allocate(fields_hub%nb_area(2))
        call field%append(field_node)
        !
        ! ----------------------------------------------------------------------
        ! Loop over time slots
        ! ----------------------------------------------------------------------
        do it=1,nt
          ! --------------------------------------------------------------------
          ! Set info for personal patch
          ! --------------------------------------------------------------------
          fields_hub%my_area(it)%l = l
          fields_hub%my_area(it)%u = u
          allocate(fields_hub%my_area(it)%f0(n(1),n(2),n(3)))
          fields_hub%my_area(it)%f0 = 0.0
          fields_hub%my_area(it)%initialized = .true.
          ! --------------------------------------------------------------------
          ! Set info for nbor patch.
          ! The ranges for nbor are the range in local flux block, that nbor is 
          ! interested in. For same level patches they coincide, for different
          ! level patches, they might be in different locations. 
          ! --------------------------------------------------------------------
          fields_hub%nb_area(it)%l = l
          fields_hub%nb_area(it)%u = u
          allocate(fields_hub%nb_area(it)%f0(n(1),n(2),n(3)))
          fields_hub%nb_area(it)%f0 = 0.0
          fields_hub%nb_area(it)%initialized = .true.
        end do
      end do

    case default
      call io%abort('interface_mod::init_fields_halo:: fields halo not set for this kind')
    end select
  END SUBROUTINE init_fields_halo

    !=============================================================================
  !> Determine the axis and face bordering a neighbor task
  !> TO DO: should return a halo instead
  !=============================================================================
  LOGICAL FUNCTION needed_for_field (task, nbtask, halo)
    class(patch_t)          :: task, nbtask
    logical, dimension(:,:,:), pointer:: halo
    !...........................................................................
    integer                 :: offset(3)
    real(8), dimension(3)   :: dist, ds(3)
    integer                 :: hal(3)
    !---------------------------------------------------------------------------
    if (.not.on) return
    needed_for_field=.false.
    !---------------------------------------------------------------------------
    ! First take care of mapping the distance between the nbtask and task
    ! centers, periodically, so this can be ignored in the next steps.
    !---------------------------------------------------------------------------
    dist = nbtask%position - task%position
    where (nbtask%periodic)
      dist = modulo(dist+0.5_8*nbtask%box,nbtask%box) - 0.5_8*nbtask%box
    end where
    !---------------------------------------------------------------------------
    ! For same level patch, the integer offset is what should be added to the
    ! task index to get the nbtask index.  For aligned meshes, offset=0.
    !---------------------------------------------------------------------------
    !if (task%level==nbtask%level) then
      ds = task%size/task%n
      offset = nint(dist/ds)
    !else
    !  call io%abort("interface_mod::range:: different level tasks not set up yet")
    !end if
    where (offset==-task%n) hal=1
    where (offset==task%n) hal=3
    where (offset==0) hal=2
    if (halo(hal(1),hal(2),hal(3)).eqv..true.) needed_for_field=.true.
    
    END FUNCTION needed_for_field


END MODULE interface_mod
