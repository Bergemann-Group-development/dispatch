MODULE rt_internal_mod
  USE io_mod
  USE bits_mod
  USE link_mod
  USE mesh_mod
  USE microphysics_mod
  USE rt_mod
  USE rt_data_mod
  USE rt_boundary_mod
  USE ray_geometry_v3_mod
  USE trace_mod
  implicit none
  private
!-----------------------------------------------------------------------
  type :: internal_data_t
    real, dimension(:,:), pointer :: I_in  => null()                    ! (ray, wl_bin) incoming intensity.
    real, dimension(:,:), pointer :: I_out => null()                    ! leaving intensity
    real, dimension(:),   pointer :: time => null()                     ! (time)
  end type
!-----------------------------------------------------------------------
  type,extends(rt_t), public :: rt_internal_t
    type (link_t), pointer :: upstream => null()                        ! quick access for intensities. list, if there is more than one
    type (rt_data_t),  pointer :: my_rt_data => null()                  ! self explanatory            >>  upstream (e.g. in case of AMR)
    class(internal_data_t), dimension(:), pointer :: packet_data => null() ! incoming intensities
    real, dimension(:,:), pointer :: I_in, I_out
    integer :: solver
    integer :: bundle_ID, dir_ID
    integer :: lwall(3)
    integer :: n_upstream                                               ! Number of upstream sources
    real :: dir_weight, ray_ds
    procedure (solve_internal_integral_2nd), pointer :: solve => null()
  contains
!> Task management
    procedure :: init   => init_rt_internal
    procedure :: update => update_rt_internal
    procedure :: print
    procedure :: check  => check_rt_internal
    procedure :: check_intensities => check_rt_get_intensities 
!> Boundary handling
    procedure :: add_upstream
    procedure :: bundle_connects_to
    procedure :: get_dynamic_legacy
    procedure :: get_same_static
    procedure :: get_static_legacy
    procedure :: get_correction_from_BC
    procedure :: get_upstream
    procedure :: has_OB
    procedure :: I_retrv
    procedure :: remove_upstream
!> Internal procedures
    procedure :: dtau
    procedure :: ray_vol_map
    procedure :: solve_internal_integral => solve_internal_integral_2nd
    procedure :: solve_internal_feautrier
    procedure :: update_q
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise tables for internal RT solution
!> =====================================================================
SUBROUTINE init_rt_internal (self)
  class(rt_internal_t) :: self
  class(ray_bndl_t),    pointer :: bundle
  class(ray_t), pointer :: packet
  integer :: pt, ray, pck, ir, irm, r
!-----------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_internal')
  call self%rt_t%init
  self%type = 'rt_internal_t'
  associate (mhd_patch => self%mhd_patch)
  self%position  = mhd_patch%position
  self%offset    = mhd_patch%offset
  self%level     = mhd_patch%level
  self%velocity  = mhd_patch%velocity
  self%size      = mhd_patch%size
  self%cell_ds   = mhd_patch%ds
  self%it        = mhd_patch%it
  self%end_time  = mhd_patch%end_time
  self%out_time  = mhd_patch%out_time
  self%li        = mhd_patch%li
  self%ui        = mhd_patch%ui
  self%time      = -2
  self%dtime     = 2e-6
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
  end associate

  call ray_geometry%head%get_my_bundle(self%bundle, self%dir_ID, self%bundle_ID)

  allocate(self%name)
  SELECT case(self%rt_type)
    case(2)
      self%name = 'rt_forward'
      self%lwall = self%bundle%lwall
    case(3)
      self%name = 'rt_reverse'
      self%lwall = self%bundle%uwall
  end SELECT

  if (self%solver ==  0) then                                           ! set the right solver
    self%solve => solve_internal_feautrier
  else
    self%solve => solve_internal_integral_2nd
  end if

  self%dir_weight = 2.*pi/self%n_ray_dirs                               ! 4pi/(2*n_ray_dirs)

  !---------------------------------------------------------------------
  ! Allocate memory slots.
  !---------------------------------------------------------------------
  allocate(self%packet_data(self%bundle%n_packets))
  allocate(self%I_in(self%bundle%n_rays, microphysics%n_wl_bins))
  allocate(self%I_out(self%bundle%n_rays, microphysics%n_wl_bins))
  self%I_in = 0; self%I_out = 0
  packet => self%bundle%rays
  self%ray_ds = packet%ds * self%cell_ds(self%bundle%axis)
  ir = 1
  do ray = 1, self%bundle%n_packets
    irm = ir + packet%n_rays - 1
    self%packet_data(ray)%I_in  => self%I_in(ir:irm,:)
    self%packet_data(ray)%I_out => self%I_out(ir:irm,:)
    allocate(self%packet_data(ray)%time(packet%n_rays))
    self%packet_data(ray)%time = self%time
    packet => packet%next
    ir = irm + 1
  end do
  nullify (packet)
  if (io%verbose>2) call self%check
  call trace_end
END SUBROUTINE init_rt_internal

SUBROUTINE print (self)
  class (rt_internal_t):: self
  !.............................................................................
  if (io%verbose>1) then
    write (*,'(a,i5,3x,a,i5,3x,a,i4,3x,a,i3,3x,a,i2,2(3x,a,f5.0),2(3x,a,3i2),2(3x,a,i4))') &
      'rt_internal_t%id:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id, &
      'dir_ID:',self%bundle%dir_ID, &
      'n_packets',self%bundle%n_packets, &
      'axis:', self%bundle%axis, &
      'theta:', self%bundle%theta*radian, &
      'phi:', self%bundle%phi*radian, &
      'lwall', self%bundle%lwall, &
      'uwall', self%bundle%uwall, &
      'nbors', self%n_nbors, &
      'n_upstream', self%n_upstream
  end if
END SUBROUTINE print

!> =====================================================================
!> Update internal RT solution
!> =====================================================================
SUBROUTINE update_rt_internal (self)
  class(rt_internal_t) :: self
  character(len=16) :: stime
  integer, save:: itimer=0
  real::ot
!-----------------------------------------------------------------------
  !write (stime,'(1p,g16.6)') self%time
  !call trace_begin('rt_t%update_rt_internal '//stime,itimer=itimer)
  call self%get_upstream
  call self%solve
  ot=self%time
  call self%update_time
  !print *,'updated internal',self%id,'old time, new time',ot,self%time, self%is_set(bits%has_BC)

  if (io%verbose > 2) then
    write (*,'(a,i5,3x,a,i5,3x,a,i4,3x,a,i3,3x,a,i2,2(3x,a,f5.0))') &
      'rt_internal_t%update:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id, &
      'dir_ID:',self%bundle%dir_ID, &
      'n_packets',self%bundle%n_packets, &
      'axis:', self%bundle%axis, &
      'theta:', self%bundle%theta*radian, &
      'phi:', self%bundle%phi*radian
  end if
  if (io%verbose > 1) call self%check
  !call trace_end (itimer)
END SUBROUTINE update_rt_internal

!> =====================================================================
!> Select appropriate method to get the upstream values
!> =====================================================================
SUBROUTINE get_upstream (self)
  class(rt_internal_t) :: self
  logical:: static, local, same
  class(link_t), pointer:: nbor
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%get_upstream', itimer=itimer)

  if (self%is_set(bits%has_BC)) then                                    ! FIXME - make inner/outer BCs separate.
    if (self%obc == 0) then                                             ! Currently only outer BC
      error stop "huh?"
    end if
    call self%get_correction_from_BC                                    ! get boundary value
  else
    static = .true.; local = .true.; same = .true.                      ! not moving? on this rank? same level?
    if (self%is_clear(bits%static)) static = .false.
    nbor => self%upstream
    do while (associated(nbor))                                         ! check nbors for these logicals
      associate (this => nbor%task)
      SELECT type(this)
        class is (rt_internal_t)
          if (this%is_clear(bits%static)) static = .false.
          if (this%is_set(bits%virtual)) local = .false.
          if (self%level /= this%level) same = .false.
      end SELECT
      end associate
      nbor => nbor%next
    end do

    if (same) then
      call self%get_same_static
    else
      if (static) then
        call self%get_static_legacy
      else
        call self%get_dynamic_legacy
      end if
    end if
  end if

  nullify (nbor)
  call trace_end (itimer)
END SUBROUTINE get_upstream

!> =====================================================================
!> Upstream data is on same MPI rank and tasks are not moving.
!> =====================================================================
SUBROUTINE get_static_legacy (self)
  class(rt_internal_t) :: self
  class(rt_internal_t), pointer :: rt_nb
  class(link_t), pointer :: nbor
  class(ray_t), pointer :: packet, nb_packet
  integer :: ray, nb_cell(3), nb_ray, nb_pt, try
  integer :: fw_c(3), pck, nb_pck
  logical, dimension(:), allocatable :: updated, all_upd
  real, dimension(microphysics%n_wl_bins) :: mean, Iin
  real :: n_mean
  real(8) :: my_c_pos(3)
  !---------------------------------------------------------------------
  call trace_begin('rt_corrections_mod::get_static')

  allocate(all_upd(self%bundle%n_packets))
  all_upd = .false.
  nbor => self%link%nbor
  do while(associated(nbor))
    nullify(rt_nb)

    associate (rt => nbor%task)
    SELECT type (rt)
    class is (rt_internal_t)
      if (self%the_right_wall(nbor, self%lwall)) then                   ! make sure nbor is on the right wall
        if (rt%rt_type == self%rt_type) then                            ! make sure nbor is of the right RT type
          if (io%verbose > 3) write (*,'(a,i5,3x,a,i5)') &
            'rt_internal task',self%id,'loading data from',rt%id

          rt_nb => rt
          packet => self%bundle%rays

          do pck = 1, self%bundle%n_packets                            ! loop over my ray packets
            allocate (updated(packet%n_rays))                              ! make sure all intensities are updated

            if (io%verbose > 3) write (*,*) 'packet', pck
            do ray = 1, packet%n_rays                                   ! loop  over my rays

              if (io%verbose > 4) write (*,*) 'ray', ray
              select case(self%rt_type)
              case(2)
                fw_c = packet%fwd_coords(:,ray)                         ! select wall coordinates
              case(3)
                fw_c =  packet%rev_coords(:,ray)                        ! select wall coordinates
              end select

              if (self%level >= rt_nb%level) then                       ! if i am not larger than my nbor
                nb_cell = rt_nb%pos2cell(self%cell2pos(fw_c))           ! corresponding cell at the nbor mesh
                if (all(nb_cell>=rt_nb%li).and.&
                    all(nb_cell<=rt_nb%ui)) then                        ! double check that we pointed to the right spot
                   self%packet_data(pck)%I_in(ray,:) = &
                     self%I_retrv(rt_nb, nb_cell)

                  self%packet_data(pck)%time(ray) = rt_nb%time          ! store the time of the update
                  updated(ray) = .true.                                 ! the ray was updated
                end if
              else                                                      ! i am larger than my nbor

                my_c_pos = self%cell2pos(fw_c)                          ! get the absolute cell position in the comp. domain
                n_mean = 0                                              ! number of rays to be averaged
                mean = 0                                                ! mean intensity

                nb_packet => rt_nb%bundle%rays
                do nb_pck = 1, rt_nb%bundle%n_packets                  ! loop over nbor ray packets
                  do nb_ray = 1, nb_packet%n_rays                       ! loop over nbor rays
                    try = -1                                            ! indicator whether this ray meets the criterion

                    select case(self%rt_type)
                    case(2)
                      do nb_pt = 1, nb_packet%n_points                  ! find the last point in the ray before the wall in question
                        if (self%point_in_cell(rt_nb%ray2pos(nb_packet,nb_ray,nb_pt),&
                            my_c_pos, self%cell_ds)) try = nb_pt
                      end do
                    case(3)
                      do nb_pt = nb_packet%n_points, 1, -1              ! find the last point in the ray before the wall in question
                        if (self%point_in_cell(rt_nb%ray2pos(nb_packet,nb_ray,nb_pt),&
                            my_c_pos, self%cell_ds)) try = nb_pt
                      end do
                    end select

                    if (try > 0) then
                      Iin = self%I_retrv(rt_nb, nb_packet%coords(:,nb_ray, try))
                      mean = mean + Iin
                      n_mean = n_mean + 1.
                    end if
                  end do
                  nb_packet => nb_packet%next
                end do

                if ((n_mean > 0).and.(self%time <= rt_nb%time)) then
                  self%packet_data(pck)%I_in(ray,:) = mean/n_mean
                  self%packet_data(pck)%time(ray) = rt_nb%time          ! store the time of the update
                  updated(ray) = .true.                                 ! the ray was updated
                end if
              end if
            end do
            if (all(updated)) all_upd(pck) = .true.
            deallocate(updated)
          packet => packet%next
          end do
        end if
      end if

    end SELECT
    end associate

  nbor => nbor%next
  end do

  if(.not.all(all_upd)) then                                                              ! TO DO: old test, to be moved to self%check
    print *, 'task ID, level, time, n_nbors', self%id, self%level, self%time, self%n_nbors
    print *, updated
    error stop "not all rays updated for in local static rt"
  end if
  deallocate(all_upd)
  nullify(rt_nb, nbor, packet, nb_packet)

  if (io%verbose>2) call self%check_intensities(1)
  call trace_end
END SUBROUTINE get_static_legacy

!> =====================================================================
!> Upstream data is on same refinement level. Direct lookup
!> =====================================================================
SUBROUTINE get_same_static (self)
  class(rt_internal_t) :: self
  class(rt_internal_t), pointer :: rt_nb => null()
  class(link_t), pointer :: nbor => null()
  integer :: mbid, ray
  integer, dimension(:), pointer :: nb_ray_pos
  !---------------------------------------------------------------------
  call trace_begin('rt_corrections_mod::get_same_static')
  if (self%n_upstream > 1) then                                         ! this should not happen.
    error stop "rt_internal::get_same_static:: (1)"
  end if
  
  select case(self%rt_type)
  case(2)
  mbid = self%bundle%next_chain_bndl_rev
  nb_ray_pos => self%bundle%ray_chain_rev
  case(3)
  mbid = self%bundle%next_chain_bndl_fw
  nb_ray_pos => self%bundle%ray_chain_fw
  end select
  nbor => self%upstream
  associate(this => nbor%task)
  SELECT type(this)
  type is (rt_internal_t)
    rt_nb => this
  end SELECT
  end associate
  
  select case(self%bundle%l_chain)
  case(:2)
    self%I_in = rt_nb%I_out
  case default
    do ray = 1, self%bundle%n_rays
      self%I_in(ray,:) = rt_nb%I_out(nb_ray_pos(ray),:)
    end do
  end select
  nullify(rt_nb)
  call trace_end
END SUBROUTINE get_same_static

!> Slow for now. Will be optimised.
FUNCTION I_retrv(self, nbor, c) result (I_out)
  class(rt_internal_t) :: self
  class(rt_internal_t), pointer :: nbor
  class(ray_t), pointer :: packet
  integer :: c(3),  my_ray, nb_ray, rid
  real, dimension(microphysics%n_wl_bins) :: I_out
  logical :: fnd
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_mod%I_retrv',4, itimer=itimer)
  fnd = .false.
  packet => nbor%bundle%rays
  do while (associated(packet).and..not.fnd)                            ! loop over nbor ray packets
    do nb_ray = 1, packet%n_rays                                        ! loop over rays in the nbor ray packet
    
      select case(nbor%rt_type)                                         ! is it fwd or reverse ray direction
        case(2)                                                         ! forward
          rid = 1                                                       ! looking for first point
        case(3)                                                         ! reverse
          rid = packet%n_points                                         ! looking for last point
      end select
      
      if (all(packet%coords(:,rid,nb_ray) == c)) then                   ! if cell coords match with what we are looking for
        fnd = .true.                                                    ! we have found the cell
        I_out = nbor%packet_data(packet%id)%I_out(nb_ray,:)
        return
      end if
      
    end do
    packet => packet%next
  end do

  call trace_end
END FUNCTION I_retrv

SUBROUTINE get_dynamic_legacy (self)
  class(rt_internal_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: local different dynamic not implemented"
END SUBROUTINE get_dynamic_legacy

!> =====================================================================
!> Upstream data is from boundary conditions.
!> =====================================================================
SUBROUTINE get_correction_from_BC (self)
  class(rt_internal_t) :: self
  integer :: ray, pck
  class(ray_t), pointer :: packet
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%get_correction_from_BC')
  SELECT case(self%rt_type)
  case(0)
  case(1)
  case(2,3)
    packet => self%bundle%rays
    do pck =1, self%bundle%n_packets
      do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
        ppos = self%absolute_point_position (packet, ray, 1, self%rt_type)
        if (self%obct == 0) then
          call rt_boundary%set_obc (ppos, self%cell_ds, self%packet_data(pck)%I_in(ray,:))
        end if
      end do
      packet => packet%next
    end do
  nullify (packet)
  end SELECT
  call trace_end
END SUBROUTINE get_correction_from_BC

!> =====================================================================
!> compute dtau given opacity and step ds along the ray
!> =====================================================================
SUBROUTINE dtau (self, opacity, dt)
  class(rt_internal_t) :: self
  real, dimension(:,:,:), pointer :: opacity, dt
  integer :: i, j
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%dtau_staggered')
  associate (n=>size(opacity,3))
  dt = opacity
  dt(:,:,2:) = dt(:,:,1:n-1) + opacity(:,:,1:n-1)
  dt = 0.5*self%ray_ds*dt
  end associate
  call trace_end
END SUBROUTINE dtau

!> =====================================================================
!> Solution of RTE using the 1st order Integral formulation. Iwf and Irv
!> are stored in case of overlapping patches.
!> =====================================================================
SUBROUTINE solve_internal_integral_1st (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet
  real, dimension(:,:,:), pointer :: ex0, ds, dtau
  real, dimension(:,:,:), pointer :: s, q, Ifw, ropac
  real, dimension(:,:), pointer :: I_in, I_out
  integer :: pt, p(3), pck, nwb, ipk, ray, id
  integer, save :: itimer(2)=0
  logical, parameter :: detailed_timer=.false.
  !---------------------------------------------------------------------
  if (.not.detailed_timer) &
  call trace_begin('rt_internal_t%integral_1st', itimer=itimer(1))
  nwb = microphysics%n_wl_bins
  packet => self%bundle%rays
  do pck = 1, self%bundle%n_packets
    if (detailed_timer) &
      call trace_begin('rt_internal_t%integral_1st(1)', itimer=itimer(1))
    associate (npt => packet%n_points, nrays => packet%n_rays, &
               rtd => self%my_rt_data, cell => packet%coords)
    allocate( dtau(nrays, nwb, npt))
    allocate(ex0  (nrays, nwb, npt))
    allocate(ds   (nrays, nwb, npt))
    allocate(s    (nrays, nwb, npt))
    allocate(q    (nrays, nwb, npt))
    allocate(Ifw  (nrays, nwb, npt))
    allocate(I_in (nrays, nwb))
    allocate(ropac(nrays, nwb,npt))
    do ray = 1, nrays
      do pt = 1, npt
        p = cell(:,pt,ray)
            s(ray,:,pt) = rtd%source (p(1),p(2),p(3),:)
        ropac(ray,:,pt) = rtd%opacity(p(1),p(2),p(3),:)
      end do
    end do
    call self%dtau(ropac, dtau)
    I_in = self%packet_data(pck)%I_in

    do pt = 2, npt    
#if defined __INTEL_COMPILER
    ex0(:,:,pt) = 2**(-dtau(:,:,pt)*1.44269504)                                 ! first exponential coeff
#else
    ex0(:,:,pt) = exp(-dtau(:,:,pt))                                            ! first exponential coeff
#endif
    ds(:,:,pt)  = (1.-ex0(:,:,pt))*s(:,:,pt-1)                                      ! Source factor
    end do
    
    if (detailed_timer) then
      call trace_end (itimer(1))
      call trace_begin('rt_internal_t%integral_1st(2)', itimer=itimer(2))
    end if
    
      !-----------------------------------------------------------------
      ! Forward / reverse solutions
      !-----------------------------------------------------------------
    pt = npt
    select case(self%rt_type)
    case(2)                                                             ! forward
      Ifw(:,:,1) = I_in * ex0(:,:,1) + ds(:,:,1)                        ! Intensity after the first point
      q  (:,:,1) = Ifw(:,:,1) - I_in                                    ! initial forward Q
      if (npt > 1) then
        Ifw(:,:,2:pt) = Ifw(:,:,1:pt-1) * ex0(:,:,2:pt) + ds(:,:,2:pt)  ! Intensity along the ray
        q  (:,:,2:pt) = Ifw(:,:,2:pt  ) - Ifw(:,:,1:pt-1)
      end if
      self%packet_data(pck)%I_out = Ifw(:,:,npt)
    case(3)                                                             ! reverse
      Ifw(:,:,pt) = I_in * ex0(:,:,pt) + ds(:,:,pt)
      q  (:,:,pt) = Ifw(:,:,pt) - I_in                                  ! initial reverse Q
      if (npt > 1) then
        do pt = npt-1, 1, -1
          Ifw(:,:,pt) = Ifw(:,:,pt+1) * ex0(:,:,pt) + ds(:,:,pt)
          q  (:,:,pt) = Ifw(:,:,pt  ) - Ifw(:,:,pt+1)
        end do
      end if
      self%packet_data(pck)%I_out = Ifw(:,:,1)
    end select
    call self%update_q (packet, q, ropac)
    deallocate(ex0, ds, s, dtau, q, Ifw, I_in, ropac)
    end associate
    packet => packet%next
    if (detailed_timer) &
      call trace_end (itimer(2))
  end do
  nullify (packet)
  if (.not.detailed_timer) &
  call trace_end (itimer(1))
END SUBROUTINE solve_internal_integral_1st

!> =====================================================================
!> Solution of RTE using the 2nd order Integral formulation.
!> Iwf and Irv are stored in case of overlapping patches.
!> =====================================================================
SUBROUTINE solve_internal_integral_2nd (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet => null()
  class(rt_data_t), pointer :: rtd
  integer, dimension (:,:,:), pointer :: cell => null()
  real, dimension(:,:,:,:), pointer :: c_source, c_opac
  real, dimension(:,:,:), pointer :: ex0, ds1, d1s, d2s
  real, dimension(:,:,:), pointer :: s, dtau, q, opac
  real, dimension(:,:), pointer :: I_in, ctau
  integer :: pt, npt, it, ray, wl, p(3), pck, nwb, ipk, iwb, nrays
  real    :: temp, dir_wgt
  integer, save :: itimer=0, itimer2=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%integral_2nd', itimer=itimer)

  dir_wgt = self%dir_weight

  nwb = microphysics%n_wl_bins
  packet => self%bundle%rays
  rtd    => self%my_rt_data
  
  do pck = 1, self%bundle%n_packets
    cell     => packet%coords
    nrays    =  packet%n_rays
    npt      =  packet%n_points
    c_opac   => rtd%opacity                                             ! cell - based opacity
    c_source => rtd%source                                              ! cell - based source function

    allocate(dtau(nrays, nwb, npt+2))
    allocate(ex0 (nrays, nwb, npt+2))
    allocate(ds1 (nrays, nwb, npt+2))
    allocate(d1s (nrays, nwb, npt))
    allocate(d2s (nrays, nwb, npt))
    allocate(s   (nrays, nwb, npt+2))
    allocate(q   (nrays, nwb, npt))
    allocate(opac(nrays, nwb, npt+2))
    allocate(ctau(nrays, nwb))
    allocate(I_in(nrays, nwb))
   ! call trace_begin('rt_t%update_rt_internal%ray_vol_map ',itimer=itimer2)
    call self%ray_vol_map (packet, nwb, c_source, s)                    ! Map source function onto rays
    call self%ray_vol_map (packet, nwb, c_opac, opac)                   ! Map opacities onto rays
   ! call trace_end (itimer2)
    call self%dtau(opac, dtau)                                          ! Calculate dtau
    I_in = self%packet_data(pck)%I_in
    
#if defined __INTEL_COMPILER
    do pt = 1, npt+1
      ex0(:,:,pt) = 2**(-dtau(:,:,pt+1)*1.44269504)                     ! first exponential coeff
      ds1(:,:,pt) = (s(:,:,pt+1) - s(:,:,pt)) / dtau(:,:,pt+1)          ! half-point centered dS/dtau
    end do   
#else
    do pt = 1, npt+1
      ex0(:,:,pt) = exp(-dtau(:,:,pt+1))                                  ! first exponential coeff
      if (any(dtau(:,:,pt+1) ==0)) print *,'dtau zero'
      ds1(:,:,pt) = (s(:,:,pt+1) - s(:,:,pt)) / dtau(:,:,pt+1)            ! half-point centered dS/dtau
    end do                                                 
#endif

    !-------------------------------------------------------------------
    !  dS/dtau and d2S/dtau2
    !-------------------------------------------------------------------
    do pt=1, npt  
      ctau = 1. / (dtau(:,:,pt+2) + dtau(:,:,pt+1))
      d1s(:,:,pt) =(ds1(:,:,pt+1) * dtau(:,:,pt+1) &
                  + ds1(:,:,pt  ) * dtau(:,:,pt+2))      * ctau         ! dS/dtau, 2nd. order accurate
      d2s(:,:,pt)= (ds1(:,:,pt+1) - ds1 (:,:,pt  )) * 2. * ctau         ! d2S/dtau2, 2nd. order accurate
    end do

    !-----------------------------------------------------------------
    !  Internal forward or reverse solutions
    !-----------------------------------------------------------------
    select case(self%rt_type)
    case(2)                                                             ! forward
      q  (:,:,1) = I_in - s(:,:,1)                                      ! initial forward Q
      do pt = 2, npt
        q(:,:,pt) = q(:,:,pt-1) * ex0(:,:,pt) - d1s(:,:,pt) &
                  * (1. - ex0(:,:,pt)) + d2s(:,:,pt) &
                  * (1. - ex0(:,:,pt) - dtau(:,:,pt)*ex0(:,:,pt))
      end do      
      self%packet_data(pck)%I_out = q(:,:,npt) + s(:,:,npt)
    case(3)                                                             ! reverse
      q  (:,:,npt) = I_in - s(:,:,npt)                                  ! initial reverse Q
      do pt = npt-1, 1, -1
        q(:,:,pt) =   q(:,:,pt+1) * ex0(:,:,pt+1) + d1s(:,:,pt) &
                  * (1. - ex0(:,:,pt+1)) + d2s(:,:,pt  ) &
                  * (1. - ex0(:,:,pt+1) - dtau(:,:,pt+1)*ex0(:,:,pt+1))
      end do
      self%packet_data(pck)%I_out = q(:,:,1) + s(:,:,1)
    case default
      error stop "rt_internal_mod::solve_internal_integral_2nd:: unknown rt type"     
    end select

    call self%update_q (packet, q, opac)
    deallocate(dtau, ex0, ds1, d1s, d2s, s, q, opac, ctau, I_in)
    packet => packet%next
  end do
  nullify (packet,rtd,cell)
  call trace_end (itimer)
END SUBROUTINE solve_internal_integral_2nd

!> =====================================================================
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map (self, rays, nwb, in, out)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt, nwb, wl
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  cell => rays%coords
  !if (nwb > 1) then
    do ray = 1, size(out,1)
      do pt = 2, size(out,3)-1                                          ! First and last points ommitted
        p = cell(:,pt-1,ray)
        !DIR$ SIMD
        do wl=1,nwb
          out(ray,wl,pt) = in(p(1),p(2),p(3),wl)
        end do
      end do
      !DIR$ SIMD
      do wl=1,nwb
        out(ray,wl,1) = in(rays%fwd_coords(1,ray), &
        rays%fwd_coords(2,ray), rays%fwd_coords(3,ray), wl)             ! Get the first point
        out(ray,wl,size(out,3)) = in(rays%rev_coords(1,ray), &
        rays%rev_coords(2,ray), rays%rev_coords(3,ray), wl)             ! Get the last point
      end do
    end do
    
  if (any(in==0)) then
  error stop "ray vol map in=0"
!  else
!  print *,"ray_vol_map says OK for IN"
  end if
    if (any(out==0)) then
    do pt = 1, size(out,3)
    print *,"pt",pt, out(:,:,pt)
    end do
  error stop "ray vol map out=0"
!  else
!  print *,"ray_vol_map says OK for OUT"
  end if
END SUBROUTINE ray_vol_map

!> =====================================================================
!> Solution of RTE using the Feautrier formulation. Forward and reverse
!> incoming radiation is ignored.
!> =====================================================================
SUBROUTINE solve_internal_feautrier (self)
  class(rt_Internal_t) :: self
  real, dimension(:,:,:), pointer :: a1, a2, a3, d
  integer                         :: j, k, l
  real                            :: a, dsdtau_1, dsdtau_2, d2sdtau2

END SUBROUTINE solve_internal_feautrier

SUBROUTINE update_q (self, rays, q, op)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), wl, ray, pt
  real, dimension(:,:,:) :: op, q
  real, dimension(:,:,:,:), pointer :: qq
  integer, dimension(:,:,:), pointer :: cell
  real :: dirw
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace_begin ('rt_internal_t_%update_q', 4, itimer=itimer)
  dirw = self%dir_weight
  qq => self%my_rt_data%q
  cell => rays%coords
  do ray = 1, size(q,1)
    do pt = 1, size(q,3)
      p = cell(:,pt,ray)
      do wl =1, size(q,2)
        qq(p(1),p(2),p(3),self%it) = &
        qq(p(1),p(2),p(3),self%it) + dirw*op(ray,wl,pt)*q(ray,wl,pt)
      end do
    end do
  end do
  call trace_end (itimer)
END SUBROUTINE update_q

!> =====================================================================
!> Check if the rays in self match up in the upstream with the rays in other
!> =====================================================================
LOGICAL FUNCTION bundle_connects_to (self, other)
  class (rt_internal_t)          :: self
  class (rt_internal_t), pointer :: other
  integer :: pp, fwd(3), rev(3)
  !.............................................................................
  bundle_connects_to = .false.
  if (.not.self%mhd_patch%overlaps(other%mhd_patch)) return             ! no overlap
  if (self%is_set(bits%has_BC+bits%boundary))        return             ! boundary patch
  
  select case(self%rt_type)
  case(2)
    fwd =  self%bundle%rays%fwd_coords(:,1) + self%bundle%dstep         ! testing only one ray in the bundle for now... For AMR: needs to run
    if (self%bundle%l_chain <2) then                                    ! through all the rays and if ANY has a match - add to the list...
      pp = 1
    else
      pp = self%bundle%ray_chain_fw(1)
    end if      
    rev = other%bundle%rays%rev_coords(:,pp)                             
  case(3)
    fwd =  self%bundle%rays%rev_coords(:,1) - self%bundle%dstep   
    if (self%bundle%l_chain <2) then
      pp = 1
    else
      pp = self%bundle%ray_chain_rev(1)      
    end if 
    rev = other%bundle%rays%fwd_coords(:,pp)                             
  end select
  
  if (all(other%pos2cell(self%cell2pos(fwd))==rev)) bundle_connects_to = .true.        
   
END function bundle_connects_to

!> =====================================================================
!> =====================================================================
LOGICAL FUNCTION has_OB (self)
  class (rt_internal_t) :: self
  class (mesh_t), pointer :: m
  integer :: i(3), a, d(3)
  !.....................................................................
  has_OB = .false.
  d = self%bundle%dstep
  
  if (self%rt_type == 2) then
    i = self%bundle%rays%fwd_coords(:,1)
    do a=1,3
      m => self%mhd_patch%mesh(a)
      if (d(a)<0) then
        if (m%p + m%r(i(a)) > m%b) has_OB = .true.
      else 
        if (m%p + m%r(i(a)) < 0.0) has_OB = .true.
      end if
    end do
  else if (self%rt_type == 3) then
    i = self%bundle%rays%rev_coords(:,1)
    do a=1,3
      m => self%mhd_patch%mesh(a)
      if (d(a) < 0) then
        if (m%p + m%r(i(a)) < 0.0) has_OB = .true.
      else
        if (m%p + m%r(i(a)) > m%b) has_OB = .true.
      end if
    end do
  end if
END FUNCTION has_OB

!> =====================================================================
!> Add a new nbor for a quick lookup of incoming intensities.
!> =====================================================================
SUBROUTINE add_upstream(self, ups)
  class(rt_internal_t) :: self
  class(link_t), pointer :: ups, new
  !---------------------------------------------------------------------
  allocate(new)
  new%task => ups%task                                                  ! simple task connection
  if(associated(self%upstream)) then
    self%upstream%prev => new
    new%next => self%upstream
    self%upstream => new
  else
    self%upstream => new
  end if
  self%n_upstream = self%n_upstream + 1
END SUBROUTINE add_upstream

!> =====================================================================
!> Remove  nbor from the quick lookup list.
!> =====================================================================
SUBROUTINE remove_upstream(self, ups)
  class(rt_internal_t) :: self
  class(link_t), pointer :: ups, this, prev
  logical :: done = .false.
  !--------------------------------------------------------------------- 
  this => self%upstream
  do while(associated(this).and..not.done)
    if (associated(this%task,ups%task)) then
      if (associated(prev)) then
        prev%next => ups%next
      else
        self%upstream => ups%next
      end if
      self%n_upstream = self%n_upstream - 1
      done = .true.
    end if
    prev => this
    this => this%next
  end do  
END SUBROUTINE remove_upstream

!> =====================================================================
!> Test if rt_internal module has any errors
!> =====================================================================
SUBROUTINE check_rt_internal(self)
  class(rt_internal_t) :: self
  integer :: i
!-----------------------------------------------------------------------
  if(.not.associated(self%my_rt_data)) then
    error stop "rt_internal%check: my_rt_data is not associated"
  end if

  if(.not.associated(self%solve)) then
    error stop "rt_internal%check: solver is not set"
  end if

  if(.not.associated(self%packet_data)) then
    error stop "rt_internal%check: packet data is missing"
  end if

  do i=1, size(self%packet_data)
    if(.not.associated(self%packet_data(i)%I_in)) then
      print *,'iter',i
      error stop "rt_internal%check: packet data bndry is missing"
    end if
    if(.not.associated(self%packet_data(i)%time)) then
      print *,'iter',i
      error stop "rt_internal%check: packet data time is missing"
    end if
    !if(.not.associated(self%packet_data(i)%ds)) then
    !  print *,'iter',i
    !  error stop "rt_internal%check: packet data ds is missing"
    !end if

    !if (any(isnan(self%packet_data(i)%ds))) then
    !!  print *,'=================== ds ================'
    !  print *,self%packet_data(i)%ds
    !  error stop "rt_internal%check: self%ds has NaN"
    !end if

  end do

END SUBROUTINE check_rt_internal

SUBROUTINE check_rt_get_intensities(self,method)
  class(rt_internal_t) :: self
  integer :: method, i
  !---------------------------------------------------------------------
  do i=1, size(self%packet_data)
    if (any(isnan(self%packet_data(i)%I_in))) then
      select case(method)
      case(1)
        print *,'ERROR: rt_internal, get intensities: local static has NaN in I_in'
      case(2)
        print *,'ERROR: rt_internal, get intensities: local dynamic has NaN in I_in'
      end select

      print *,'rt_type, task ID, BC bit:', self%rt_type, self%id, self%is_set(bits%has_bc)
      error stop "STOPPING"
    end if
  end do

END SUBROUTINE check_rt_get_intensities

SUBROUTINE check_rt_corrections(self)
  class(rt_internal_t) :: self
!-----------------------------------------------------------------------

  if(.not.associated(self%my_rt_data)) then
    error stop "rt_corrections%check: my_rt_data is not associated"
  end if


  if(.not.associated(self%mhd_patch)) then
    error stop "rt_corrections%check: rt_t%mhd_patch is not associated"
  end if

  if(.not.associated(self%mhd_link)) then
    error stop "rt_corrections%check: rt_t%mhd_link is not associated"
  end if

  if(.not.associated(self%my_nb_link)) then
    error stop "rt_corrections%check: rt_t%my_nb_link is not associated"
  end if

  if(.not.associated(self%bundle)) then
    error stop "rt_corrections%check: rt_t%my_rays is not associated"
  end if
END SUBROUTINE check_rt_corrections

END MODULE rt_internal_mod
