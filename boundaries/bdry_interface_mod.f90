MODULE bdry_interface_mod
  USE kinds_mod
  USE trace_mod
  USE io_mod
  USE dll_mod
  USE link_mod
  implicit none
  private
! ------------------------------------------------------------------------------
  type, extends(dll_t), public:: nbor_dll_t
    character(len=8):: label
    logical:: updated=.false.
  contains
    procedure:: add_halo_node
    procedure:: get_nbor_dll
  end type
  !=============================================================================
  !> area_t is a data type holding interface information for all fluxes and fields
  !=============================================================================
  type, public:: area_t
    integer:: l(3), u(3)                                          ! range to point f1 to
    logical:: initialized = .false.
    real(8):: t=0d0, dt=0d0, dt0=0d0, dt1=0d0                     ! previous t,dt
    real(kind=var_k), dimension(:,:,:), allocatable:: f0          ! previous flux
  contains
    procedure:: init => init_area
    procedure:: dealloc => dealloc_area
  end type

  !=============================================================================
  !> Areas in the halo are split into two parts:
  !> 1) my_area, where I pick up halo information, collected from nbors
  !> 2) nbor_area, where I prepare data for nbor to be picked up
  !> while the two areas coincide in unigrid case, but if patch%level differs,
  !> they might not, e.g. in the case of a parent/child
  !>============================================================================ 
  type, public:: fields_hub_t
    type(area_t), dimension(:), pointer:: my_area
    type(area_t), dimension(:), pointer:: nb_area
    real(kind=var_k), dimension(:,:,:), pointer:: f1 => null()    ! current flux
    logical:: is_ahead=.false.
    logical:: flux_correction=.false.
    integer:: dir
    integer:: iv
    integer:: id
  contains
    procedure:: update => update_fields_hub
    procedure,private:: deposit4
    procedure,private:: deposit5
    generic,public:: deposit => deposit4, deposit5
  end type
  ! ----------------------------------------------------------------------------
  integer:: verbose=0

  CONTAINS
  !=============================================================================
  !> Adds a halo node and returns
  !> a pointer to the new node. It can be hydro (id==1), emf(id==2), etc.
  !=============================================================================
  SUBROUTINE add_halo_node(self, id, label, node)
    class(nbor_dll_t):: self
    integer:: id
    character(len=*):: label
    class(dll_node_t), pointer:: node
    ! --------------------------------------------------------------------------
    class(nbor_dll_t), pointer:: nbor_dll
    ! ..........................................................................
    nullify(node)
    call self%find_id(id,node)
    if (.not.associated(node)) then
      allocate(node)
      node%id=id
      call self%append(node)
      allocate(nbor_dll)
      node%car => nbor_dll
      call nbor_dll%init('nbor_dll')
      nbor_dll%label = trim(label)
    end if
  END SUBROUTINE add_halo_node

  SUBROUTINE get_nbor_dll(self, node, dll, label)
    class(nbor_dll_t):: self
    class(dll_node_t), pointer:: node
    class(nbor_dll_t), pointer:: dll
    character(len=*), optional:: label
    ! --------------------------------------------------------------------------
    class(*), pointer:: car
    ! ..........................................................................
    car => node%car
    select type(car)
    class is (nbor_dll_t)
      dll => car
    class default
      call io%abort('garbage in car at get_nbor_dll, caller: '//label)
    end select
  END SUBROUTINE get_nbor_dll

  SUBROUTINE deposit5(self, t, dt, it, field)
    class(fields_hub_t):: self
    real(kind=var_k), dimension(:,:,:,:,:), pointer:: field
    integer:: it
    real(8):: t, dt
    ! --------------------------------------------------------------------------
    integer:: l(3), u(3)
    ! ..........................................................................

    ! --------------------------------------------------------------------------
    ! Here, some interpolation in space might be needed for different level
    ! patches.
    ! --------------------------------------------------------------------------
    l = self%nb_area(it)%l
    u = self%nb_area(it)%u
    self%nb_area(2)%t = self%nb_area(it)%t+self%nb_area(it)%dt
    self%nb_area(2)%dt= self%nb_area(it)%dt
    if (self%flux_correction) then
      self%nb_area(2)%f0 = self%nb_area(it)%f0
    end if
    self%nb_area(it)%t = t+dt
    self%nb_area(it)%dt= dt
    if (self%flux_correction) then
      self%nb_area(it)%f0 = field(l(1):u(1),l(2):u(2), l(3):u(3),self%iv,self%dir)
    end if
  END SUBROUTINE deposit5

  SUBROUTINE deposit4(self, t, dt, it, field)
    class(fields_hub_t):: self
    real(kind=var_k), dimension(:,:,:,:), pointer:: field
    integer:: it
    real(8):: t, dt
    ! --------------------------------------------------------------------------
    integer:: l(3), u(3)
    ! ..........................................................................

    ! --------------------------------------------------------------------------
    ! Here, some interpolation in space might be needed for different level
    ! patches.
    ! --------------------------------------------------------------------------
    l = self%nb_area(it)%l
    u = self%nb_area(it)%u
    self%nb_area(2)%t = self%nb_area(it)%t+self%nb_area(it)%dt
    self%nb_area(2)%dt= self%nb_area(it)%dt
    if (self%flux_correction) then
      self%nb_area(2)%f0 = self%nb_area(it)%f0
    end if
    self%nb_area(it)%t = t+dt
    self%nb_area(it)%dt= dt
    if (self%flux_correction) then
      self%nb_area(it)%f0 = field(l(1):u(1),l(2):u(2), l(3):u(3),self%dir)
    end if
  END SUBROUTINE deposit4


  SUBROUTINE update_fields_hub(self, t, dt, it, debug)
    class(fields_hub_t):: self
    real(8):: t, dt
    integer:: it
    logical::debug
    ! --------------------------------------------------------------------------
    real(8):: t0, t1, t2, t3, dt0, dt1, dtold, new_t, new_dt
    real(kind=var_k):: f1(size(self%my_area(it)%f0,1), &
                          size(self%my_area(it)%f0,2), &
                          size(self%my_area(it)%f0,3))
    real(kind=var_k):: f0(size(self%my_area(it)%f0,1), &
                          size(self%my_area(it)%f0,2), &
                          size(self%my_area(it)%f0,3))
    real(kind=var_k):: fold(size(self%my_area(it)%f0,1), &
                            size(self%my_area(it)%f0,2), &
                            size(self%my_area(it)%f0,3))                      
    integer, save:: n(4)=0
    ! ..........................................................................
    !-----------------------------------------------------------------------------
    ! If the task on the two sides of a boundary update alternatingly, then one
    ! task pushing its update information to the other task is sufficient.  But
    ! if the same task updates twice, it also needs to have the corresponding
    ! information.  Therefore, after updating the flux at the boundary, the task
    ! should push the information both to the other task, and to itself
    !-----------------------------------------------------------------------------
    call trace%begin ('face_t%update')
    t0 = self%my_area(it)%t - self%my_area(it)%dt                               ! previous start of interval
    t1 = t                                                                      ! start of time step for this task
    t2 = self%my_area(it)%t                                                     ! current face time
    t3 = t+dt                                                                   ! end of time step for this task
    if (verbose > 1) &
      write (stdout,'(a,i6,6f10.6)') 'id, times:', &
        self%id, t0, t1, t2, t3, dt
    ! -------------------------------------------------------------------------- 
    if (t1 >= t2) then               
    !---------------------------------------------------------------------------
    ! the initial state was extrapolated
    ! This task advances the face time and is fully authoritative
    ! time: ----t2---------t1---------------------t3
    ! flux: f0f0f0f0f0f0f0f0f1f1f1f1f1f1f1f1f1f1f1f1
    !---------------------------------------------------------------------------
      dt0 = 0.0                     ! nbor f0 is not the authoritative flux
      dt1 = dt                      ! self%f1 is the authoritative flux
      dtold=0.0
      new_t = t3                    ! new face time
      new_dt = t3-t1                  ! interval where the current %f1 is valid
      n(1) = n(1)+1
      f0 = self%my_area(it)%f0
      f1 = self%f1
      fold = self%my_area(2)%f0

    else if (t3 >= t2 .and.t0<t1) then        
    !---------------------------------------------------------------------------
    ! the time steps overlap
    ! This task advances the face time, and is partially authoritative
    ! time: ----t1--------t2----------------------t3
    ! flux: f1f1f1f1f1f1f1f1f0f0f0f0f0f0f0f0f0f0f0f0
    !---------------------------------------------------------------------------
      dt0 = t2-t1                   ! nbor f0 is 1st authoritative flux
      dt1 = t3-t2                   ! self%f1 is 2nd authoritative flux
      dtold=0.0
      new_t = t3           ! this tasks advanced the face time
      new_dt = t3-t1         ! interval where the current %f1 is valid
      n(2) = n(2)+1
      f0 = self%my_area(it)%f0
      f1 = self%f1
      fold=self%my_area(2)%f0
    else if (t3 >= t2 .and.t0>=t1) then         
    !---------------------------------------------------------------------------
    ! the time steps overlap
    ! This task advances the face time, and is partially authoritative
    ! time: --t1----t0--------t2--------------------t3
    ! flux: f0f1f1f1f1f0f0f0f0f0f1f1f1f1f1f1f1f1f1f1f1
    !---------------------------------------------------------------------------
      if (verbose>0) print *, 'using old t'
      dt0 = t2-t0                   ! nbor f0 is 1st authoritative flux
      dt1 = t3-t2                   ! self%f1 is 2nd authoritative flux
      dtold=t0-t1
      new_t = t3           ! this tasks advanced the face time
      new_dt = t3-t1         ! interval where the current %f1 is valid
      n(2) = n(2)+1
      f0 = self%my_area(it)%f0
      f1 = self%f1
      fold=self%my_area(2)%f0         ! f1 has authority over the first part of the interval 
     else if (t0>=t1 .and. t2>t3) then
      !---------------------------------------------------------------------------
      ! Time steps overlap.
      ! This task advances the face time and is partially authoritative
      ! time: t1------t0-------------------t3------t2
      ! flux: f1f1f1f1f1f0f00f0f0f0f0f0f0f0f0f1f1f1f1
      !---------------------------------------------------------------------------
      !dt0 = t0-t1      ! self%f1 is 1st authoritative flux
      !dt1 = t3-t0      ! nbor f0 is 2nd authoritative flux
      !dtold=t0-t1
      !f0 = self%f1
      !f1 = self%my_area(it)%f0
      !fold = self%my_area(2)%f0
      !new_dt = dt0+dt1
      !new_t = t3
      !n(4) = n(4)+1

      dt0 = t0-t1      ! self%f1 is 1st authoritative flux
      dt1 = t3-t0      ! nbor f0 is 2nd authoritative flux
      dtold=0.0
      f0 = self%my_area(2)%f0
      f1 = self%my_area(it)%f0
      fold = self%f1
      new_dt = t3-t1
      new_t = t3
      n(4) = n(4)+1
    else if (t1 > t0 .and.t3<=t2) then         ! t3 is still inside the [t0,t2]
      !print *,  'face_t%update(): task is behind, 3'
      !---------------------------------------------------------------------------
      ! This task does not advance the face time, and the stored flux is used,
      ! which is OK as long as the current task is fully inside the same timestep:
      ! time: ----t0----t1-------------------t3-------t2
      ! flux: fxfxfxf0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0
      !---------------------------------------------------------------------------
      dt0 = dt                      ! nbor f0 is 100% authoritative still
      dt1 = 0.0                     ! self%f1 has zere weight
      dtold = 0.0
      !self%t = t1                  ! self%t remains unchanged
      self%my_area(it)%dt = dt0
      n(3) = n(3)+1
      f0 = self%my_area(it)%f0
      f1 = self%f1
      fold = self%my_area(2)%f0
      new_t = t3
      new_dt = t3-t1
    else if (t3<t0) then
      call io%abort('t3<to')
      !-------------------------------------------------------------------------
      ! The current task is more than one step behind. Old nbor data is evaluated
      ! ------------------------------------------------------------------------
      ! This task advances the face time and is fully authoritative
      ! time: told------t1------------t3---------t0------------t2
      ! flux:   fofofofof1f1f1f1f1f1f1f1f0f0f0f0f0f0f0f0f0f0f0f0f
      !-------------------------------------------------------------------------
      dt0 = 0.0                     ! nbor f0 is not the authoritative flux
      dt1 = 0.0                     ! self%f1 is not the authoritative flux
      dtold=t1-self%my_area(2)%t    ! fold is the authoritative flux
      new_t = t3           ! new face time
      new_dt = dt1         ! interval where the current %f1 is valid
      n(1) = n(1)+1
      f0 = self%my_area(it)%f0
      f1 = self%f1
      fold = self%my_area(2)%f0
    else
      print *, t0,t1,t2,t3
      print *, t0>t1, t1>t2, t2>t3
      call io%abort('else.... so more cases needed.')
      !---------------------------------------------------------------------------
      ! FIXME: If the current task is more than one step behind, more needed,
      ! since part of the time interval had a different flux, no longer there
      ! time: ----t1------t0-------------------t3------t2
      ! flux: fxfxfxfxfxfxfxf0f00f0f0f0f0f0f0f0f0f0f0f0f0
      !---------------------------------------------------------------------------
    !  print *,  'face_t%update(): task is behind, 5, need more time slots..'
!
    !  dt0 = dt                      ! self%f0 is 100% authoritative still
    !  dt1 = 0.0                     ! self%f1 has zere weight
      !self%t = t1                  ! self%t remains unchanged
    !  new_dt = dt0
    !  n(4) = n(4)+1
    !  f1 = self%f1
    !  f0 = self%my_area(2)%f0
    end if
    !.............................................................................
    if (verbose > 1) then
      if (debug) then
      write (stdout,'(1x,a,4f9.6,1 x,3f9.6,1x,3e11.3,4i6)') &
        'face_t%update: t0,t1,t2,t3,dt0,dt1,f0,f1 =', &
        t0, t1, t2, t3, dt, dt0, dt1, sum(self%my_area(it)%f0), sum(self%f1), &
        sum((dt0*f0 + dt1*f1 + dtold*fold)/(dt0+dt1+dtold+tiny(1d0))), n
    end if
    end if
    if (verbose > 3) then
      write (stdout,'(1x,a,4(2x,2i4))') &
        'face_t%update: shapes =', shape(self%my_area(it)%f0), shape(self%f1)
    end if
    !-----------------------------------------------------------------------------
    ! self%f1 points to the current flux variable in the face plane.  It
    ! needs to be replace with a weighted average of itself and the previous
    ! flux, as delivered by the opposite side task
    !-----------------------------------------------------------------------------
    !f1 = self%f1                                        ! save current flux
    if (self%flux_correction) then
      self%f1 = (dt0*f0 + dt1*f1 + dtold*fold)/(dt0+dt1+dtold+tiny(1d0))   ! flux to use
      if (debug) then
      if (verbose > 0) &
        write (stdout,'(a,1p,6e12.4)') &
        'correction:', sum(self%my_area(it)%f0), sum(f1), sum(self%f1), dt0, dt1, dt0+dt1
      end if
    end if
    !-----------------------------------------------------------------------------
    ! If (and only if) the current times extends into new territory, the flux
    ! registered for the face should be the flux computed in the current task.
    ! This flux will be used only if the current task updates again, before the
    ! nbor task updates and replaces the flux value and the time value.
    !-----------------------------------------------------------------------------
    !if (t3 >= t2) then                                  ! if ahead of face time

      self%my_area(2)%f0 = self%my_area(it)%f0                ! f0 is old f0 
      self%my_area(2)%dt = self%my_area(it)%dt                  ! save time interval
      self%my_area(2)%t = self%my_area(it)%t                    ! save time interval
      self%my_area(it)%f0 = self%f1                              ! current f1 is new f0
      self%my_area(it)%dt = new_dt
      self%my_area(it)%t = new_t
      self%my_area(it)%dt0 = self%my_area(it)%dt               ! save time interval

      self%is_ahead = .true.                    ! indicate update
      if (verbose > 1) &
        write (stdout,*) 'is_ahead: axis, iv =',self%dir,self%iv
    !else
    !end if
    call trace%end ()
    !print *,n 
!print *, t0,t1,t2,t3

  END SUBROUTINE update_fields_hub

  SUBROUTINE init_area (self,rng,l,u,pos)
    class(area_t):: self
    integer:: rng(3), pos(3), l(3), u(3)
    ! ..........................................................................
    integer:: i
    ! --------------------------------------------------------------------------
    allocate(self%f0(rng(1),rng(2),rng(3)))
    self%f0 = 0.0
    do i=1,3
      if (pos(i)==1) then
        self%l(i)=l(i)
        self%u(i)=l(i)
      else if (pos(i)==2) then
        self%l(i)=l(i)
        self%u(i)=u(i)
      else
        self%l(i)=u(i)
        self%u(i)=u(i)
      end if
    end do
  END SUBROUTINE init_area

  SUBROUTINE dealloc_area(self)
    class(area_t):: self
    ! --------------------------------------------------------------------------
    if (allocated(self%f0)) deallocate(self%f0)
    self%initialized = .false.
  END SUBROUTINE dealloc_area




END MODULE bdry_interface_mod
