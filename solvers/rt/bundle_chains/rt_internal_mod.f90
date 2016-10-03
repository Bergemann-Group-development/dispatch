MODULE rt_internal_mod
  USE io_mod
  USE bits_mod
  USE link_mod
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
    real, dimension(:,:), pointer :: bndry => null()                    ! (ray, wl_bin) incoming intensity.
    real, dimension(:),   pointer :: time => null()                     ! (time)
    !real, dimension(:,:), pointer :: ds   => null()                     ! (ray,point) NOTE: not needed for MARCONI
    real, dimension(:,:), pointer :: I_out => null()                    ! leaving intensity
    integer, dimension(:,:), pointer :: cell => null()                  ! convenience, to identify which cell the ray point corresponds to
  end type
!-----------------------------------------------------------------------
  type,extends(rt_t), public :: rt_internal_t
    type (rt_data_t),     pointer :: my_rt_data => null()
    class(internal_data_t), dimension(:), pointer :: packet_data => null() ! incoming intensities
    integer :: solver
    integer :: BC_direction=0
    integer :: bundle_ID, dir_ID
    integer :: lwall(3)
    real :: dir_weight, ray_ds
    procedure (solve_internal_integral_1st), pointer :: solve => null()
  contains
    procedure :: dtau
    procedure :: init   => init_rt_internal
    procedure :: update => update_rt_internal
    procedure :: check  => check_rt_internal
    procedure :: check_intensities => check_rt_get_intensities
    procedure :: needs_BC
    procedure :: get_upstream
    procedure :: get_static
    procedure :: get_dynamic
    procedure :: get_correction_from_BC
    procedure :: solve_internal_integral => solve_internal_integral_2nd
    procedure :: solve_internal_feautrier
    procedure :: update_q
    procedure :: I_retrv
    procedure :: print
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise tables for internal RT solution
!> =====================================================================
SUBROUTINE init_rt_internal (self)
  class(rt_internal_t) :: self
  class(ray_bndl_t),    pointer :: my_rays
  class(ray_t), pointer :: packet
  integer :: pt, ray, pck, r
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

  call ray_geometry%head%get_my_bundle(self%my_rays, self%dir_ID, self%bundle_ID)

  allocate(self%name)
  SELECT case(self%rt_type)
    case(2)
      self%name = 'rt_forward'
      self%lwall = self%my_rays%lwall
    case(3)
      self%name = 'rt_reverse'
      self%lwall = self%my_rays%uwall
  end SELECT

  if (self%solver ==  0) then                                           ! set the right solver
    self%solve => solve_internal_feautrier
  else
    self%solve => solve_internal_integral_1st
  end if

  self%dir_weight = 4.*pi/self%n_ray_dirs

  !---------------------------------------------------------------------
  ! Allocate memory slots.
  !---------------------------------------------------------------------
  allocate(self%packet_data(self%my_rays%n_packets))
  packet => self%my_rays%rays
  self%ray_ds = packet%ds * self%cell_ds(self%my_rays%axis)
  do ray = 1, self%my_rays%n_packets
    allocate(self%packet_data(ray)%bndry (packet%n_rays,microphysics%n_wl_bins))
    allocate(self%packet_data(ray)%I_out (packet%n_rays,microphysics%n_wl_bins))
    allocate(self%packet_data(ray)%cell (3,packet%n_rays))
    allocate(self%packet_data(ray)%time(packet%n_rays))
    !allocate(self%packet_data(ray)%ds  (packet%n_rays,packet%n_points))  - not needed for MARCONI
    self%packet_data(ray)%bndry(:,:) = 0
    self%packet_data(ray)%I_out(:,:) = 0
    self%packet_data(ray)%time = self%time
    !self%packet_data(ray)%ds = packet%ds * self%cell_ds(self%my_rays%axis) ! FIXME: should not be like this for non-cubic cells
    select case(self%rt_type)
    case(2)
      pt = packet%n_points
    case(3)
      pt = 1
    end select
    do r = 1, packet%n_rays
      self%packet_data(ray)%cell(:,r) = packet%coords(:,pt,r)
    end do
    packet => packet%next
  end do
  nullify (packet)
  if (io%verbose>2) call self%check
  call trace_end
END SUBROUTINE init_rt_internal

SUBROUTINE print (self)
  class (rt_internal_t):: self
  !.............................................................................
  if (io%verbose>1) then
    write (*,'(a,i5,3x,a,i5,3x,a,i4,3x,a,i3,3x,a,i2,2(3x,a,f5.0),2(3x,a,3i2),3x,a,i4)') &
      'rt_internal_t%id:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id, &
      'dir_ID:',self%my_rays%dir_ID, &
      'n_packets',self%my_rays%n_packets, &
      'axis:', self%my_rays%axis, &
      'theta:', self%my_rays%theta*radian, &
      'phi:', self%my_rays%phi*radian, &
      'lwall', self%my_rays%lwall, &
      'uwall', self%my_rays%uwall, &
      'nbors', self%n_nbors
  end if
END SUBROUTINE print

!> =====================================================================
!> Update internal RT solution
!> =====================================================================
SUBROUTINE update_rt_internal (self)
  class(rt_internal_t) :: self
  character(len=16) :: stime
  integer, save:: itimer=0
!-----------------------------------------------------------------------
  !write (stime,'(1p,g16.6)') self%time
  !call trace_begin('rt_t%update_rt_internal '//stime,itimer=itimer)
  call self%get_upstream
  call self%solve
  call self%update_time
  if (io%verbose > 2) then
    write (*,'(a,i5,3x,a,i5,3x,a,i4,3x,a,i3,3x,a,i2,2(3x,a,f5.0))') &
      'rt_internal_t%update:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id, &
      'dir_ID:',self%my_rays%dir_ID, &
      'n_packets',self%my_rays%n_packets, &
      'axis:', self%my_rays%axis, &
      'theta:', self%my_rays%theta*radian, &
      'phi:', self%my_rays%phi*radian
  end if
  if (io%verbose > 1) call self%check
  !call trace_end (itimer)
END SUBROUTINE update_rt_internal

!> =============================================================================
!> =============================================================================
LOGICAL FUNCTION needs_BC (self)
  class (rt_internal_t)          :: self
  !.............................................................................
  needs_BC = .false.
  if (self%mhd_patch%is_set(bits%has_BC+bits%boundary)) then            ! boundary
    if (self%BC_direction > 0) then                                     ! valid BC_direction
      if (self%bundle%dstep(self%BC_direction) > 0) needs_BC = .true.   ! step in BC_direction
    end if
  end if
END FUNCTION needs_BC

!> =====================================================================
!> Select appropriate method to get the upstream values
!> =====================================================================
SUBROUTINE get_upstream (self)
  class(rt_internal_t) :: self
  logical:: static, local
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
    static = .true.; local = .true.                                     ! not moving? on this rank?
    if (self%is_clear(bits%static)) static = .false.
    nbor => self%link%nbor
    do while (associated(nbor))                                         ! check nbors for these logicals
      associate (this => nbor%task)
      SELECT type(this)
        class is (rt_internal_t)
          if (this%is_clear(bits%static)) static = .false.
          if (this%is_set(bits%virtual)) local = .false.
      end SELECT
      end associate
      nbor => nbor%next
    end do

    if (static) then
      call self%get_static
    else
      call self%get_dynamic
    end if
  end if

  nullify (nbor)
  call trace_end (itimer)
END SUBROUTINE get_upstream

!> =====================================================================
!> Upstream data is on same MPI rank and tasks are not moving.
!> =====================================================================
SUBROUTINE get_static (self)
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

  allocate(all_upd(self%my_rays%n_packets))
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
          packet => self%my_rays%rays

          do pck = 1, self%my_rays%n_packets                            ! loop over my ray packets
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
                   self%packet_data(pck)%bndry(ray,:) = &
                     self%I_retrv(rt_nb, nb_cell)

                  self%packet_data(pck)%time(ray) = rt_nb%time          ! store the time of the update
                  updated(ray) = .true.                                 ! the ray was updated
                end if
              else                                                      ! i am larger than my nbor

                my_c_pos = self%cell2pos(fw_c)                          ! get the absolute cell position in the comp. domain
                n_mean = 0                                              ! number of rays to be averaged
                mean = 0                                                ! mean intensity

                nb_packet => rt_nb%my_rays%rays
                do nb_pck = 1, rt_nb%my_rays%n_packets                  ! loop over nbor ray packets
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
                  self%packet_data(pck)%bndry(ray,:) = mean/n_mean
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
END SUBROUTINE get_static

!> Slow for now. Will be optimised.
FUNCTION I_retrv(self, nbor, c) result (I_out)
  class(rt_internal_t) :: self
  class(rt_internal_t), pointer :: nbor
  class(ray_t), pointer :: packet
  integer :: c(3),  my_ray, nb_ray
  real, dimension(microphysics%n_wl_bins) :: I_out
  logical :: fnd
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_mod%I_retrv',4, itimer=itimer)
  fnd = .false.
  packet => nbor%my_rays%rays
  do while (associated(packet).and..not.fnd)
    do nb_ray = 1, packet%n_rays
      if (all(nbor%packet_data(packet%id)%cell(:,nb_ray) == c)) then
        fnd = .true.
        I_out = nbor%packet_data(packet%id)%I_out(nb_ray,:)
      end if
    end do
    packet => packet%next
  end do

  call trace_end
END FUNCTION I_retrv

SUBROUTINE get_dynamic (self)
  class(rt_internal_t) :: self
  !---------------------------------------------------------------------
  error stop "ERROR: local different dynamic not implemented"
END SUBROUTINE get_dynamic

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
    packet => self%my_rays%rays
    do pck =1, self%my_rays%n_packets
      do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
        ppos = self%absolute_point_position (packet, ray, 1, self%rt_type)
        if (self%obct == 0) then
          call rt_boundary%set_obc (ppos, self%cell_ds, self%packet_data(pck)%bndry(ray,:))
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
FUNCTION dtau (self, opacity) result (dt)
  class(rt_internal_t) :: self
  real, dimension(:,:,:)  :: opacity
  real, dimension(size(opacity,1), size(opacity,2), size(opacity,3)):: dt
  integer :: i, j
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%dtau')
  associate (n=>size(opacity,3))
  dt = opacity
  dt(:,:,2:) = dt(:,:,1:n-1) + opacity(:,:,1:n-1)
  dt = 0.5*self%ray_ds*dt
  end associate
  call trace_end
END FUNCTION dtau

!> =====================================================================
!> Solution of RTE using the 1st order Integral formulation. Iwf and Irv
!> are stored in case of overlapping patches.
!> =====================================================================
SUBROUTINE solve_internal_integral_1st (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet
  real, dimension(:,:,:), allocatable :: ex0, ds, dtau
  real, dimension(:,:,:), allocatable :: s, q, Ifw, ropac
  real, dimension(:,:), allocatable :: I_in, I_out
  integer :: pt, p(3), pck, nwb, ipk, ray
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%integral_1st', itimer=itimer)

  nwb = microphysics%n_wl_bins
  packet => self%my_rays%rays

  do pck = 1, self%my_rays%n_packets

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
    dtau = self%dtau(ropac)
    I_in = self%packet_data(pck)%bndry
    pt = npt
#if defined __INTEL_COMPILER
    ex0 = 2**(-dtau*1.44269504)                                         ! first exponential coeff
#else
    ex0 = exp(-dtau)                                                    ! first exponential coeff
#endif
    ds  = (1.  - ex0) * s                                               ! Source factor
      !-----------------------------------------------------------------
      ! Forward / reverse solutions
      !-----------------------------------------------------------------
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
  end do
  nullify (packet)
  call trace_end (itimer)
END SUBROUTINE solve_internal_integral_1st

!> =====================================================================
!> Solution of RTE using the 2nd order Integral formulation.
!> Iwf and Irv are stored in case of overlapping patches.
!> =====================================================================

SUBROUTINE solve_internal_integral_2nd (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet
  real, dimension(:,:,:), allocatable :: ex0, ex1, ex2, ds1, d1s, d2s
  real, dimension(:,:,:), allocatable :: s, dtau, idtau, ctau, q, ropac
  real, dimension(:,:), allocatable :: I_in
  integer :: pt, ipt, it, ray, wl, p(3), pck, nwb, ipk
  real    :: temp, dir_wgt
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%integral_2nd', itimer=itimer)

  dir_wgt = self%dir_weight

  nwb = microphysics%n_wl_bins
  packet => self%my_rays%rays

  do pck = 1, self%my_rays%n_packets

    associate (npt => packet%n_points, nrays => packet%n_rays, &
               rtd => self%my_rt_data, cell => packet%coords)

    allocate( dtau(nrays, nwb, npt))
    allocate(idtau(nrays, nwb, npt))
    allocate(ex0  (nrays, nwb, npt))
    allocate(ex1  (nrays, nwb, npt))
    allocate(ex2  (nrays, nwb, npt))
    allocate(ds1  (nrays, nwb, npt))
    allocate(d1s  (nrays, nwb, npt))
    allocate(d2s  (nrays, nwb, npt))
    allocate(s    (nrays, nwb, npt))
    allocate(q    (nrays, nwb, npt))
    allocate(ctau (nrays, nwb, npt))
    allocate(I_in (nrays, nwb))
    allocate(ropac(nrays, nwb,npt))

    do ray = 1, nrays
      do pt = 1, npt
        p = cell(:,pt,ray)
            s(ray,:,pt) = rtd%source (p(1),p(2),p(3),:)
        ropac(ray,:,pt) = rtd%opacity(p(1),p(2),p(3),:)
      end do
    end do

    dtau = self%dtau(ropac)
    I_in = self%packet_data(pck)%bndry
    pt = npt
    idtau = 1. / dtau
#if defined __INTEL_COMPILER
      ex0 = 2**(-dtau*1.44269504)                                       ! first exponential coeff
#else
      ex0 = exp(-dtau)                                                  ! first exponential coeff
#endif
      ex1 = 1.  - ex0                                                   ! second exponential coeff
      ex2 = ex1 - dtau * ex0                                            ! third exponential coeff
      ds1(:,:,2:pt) = (s(:,:,2:pt) -s(:,:,1:pt-1)) * idtau(:,:,1:pt-1)  ! half-point centered dS/dtau
      !-----------------------------------------------------------------
      !  dS/dtau and d2S/dtau2
      !-----------------------------------------------------------------
     ctau(:,:,2:pt-1) = 1. / (dtau(:,:,2:pt-1) + dtau(:,:,3:pt))
      d1s(:,:,2:pt-1) = &
     (ds1(:,:,3:pt  ) * dtau(:,:,2:pt-1) + ds1(:,:,2:pt-1) * &
     dtau(:,:,3:pt) ) * ctau(:,:,2:pt-1)                                ! dS/dtau, 2nd. order accurate
      d2s(:,:,2:pt-1) = (ds1(:,:,3:pt  ) - ds1(:,:,2:pt-1))* &
                   2. * ctau(:,:,2:pt-1)                                ! d2S/dtau2, 2nd. order accurate


      d2s(:,:,1)  = d2s(:,:,2)                                          ! assume constant d2S/dtau2
      d1s(:,:,1)  = d1s(:,:,2)    - dtau(:,:,2)  * d2s(:,:,2)           ! extrapolate to 1
      d2s(:,:,pt) = d2s(:,:,pt-1)                                       ! assume constant d2S/dtau2
      d1s(:,:,pt) = d1s(:,:,pt-1) + dtau(:,:,pt) * d2s(:,:,pt)          ! extrapolate to n
      !-----------------------------------------------------------------
      !  Internal forward or reverse solutions
      !-----------------------------------------------------------------

    select case(self%rt_type)
    case(2)                                                             ! forward
      q  (:,:,1) = I_in -   s(:,:,1)                                    ! initial forward Q
      do pt = 2, npt
        q(:,:,pt) = q(:,:,pt-1) * ex0(:,:,pt) - d1s(:,:,pt) &
                  * ex1(:,:,pt) + d2s(:,:,pt) * ex2(:,:,pt)
      end do
      self%packet_data(pck)%I_out = q(:,:,npt) + s(:,:,npt)

    case(3)                                                             ! reverse
      q  (:,:,pt) = I_in -   s(:,:,pt)                                  ! initial reverse Q
      do pt = npt-1, 1, -1
        q(:,:,pt) =   q(:,:,pt+1) * ex0(:,:,pt+1) + d1s(:,:,pt) &
                  * ex1(:,:,pt+1) + d2s(:,:,pt  ) * ex2(:,:,pt+1)
      end do
      self%packet_data(pck)%I_out = q(:,:,1) + s(:,:,1)
    end select

    call self%update_q (packet, q, ropac)
    deallocate(ex0, ex1, ex2, ds1, d1s, d2s, ctau, s, dtau, &
               idtau, q, I_in, ropac)
    end associate
    packet => packet%next
  end do
  nullify (packet)
  call trace_end (itimer)
END SUBROUTINE solve_internal_integral_2nd

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
    if(.not.associated(self%packet_data(i)%bndry)) then
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
    if (any(isnan(self%packet_data(i)%bndry))) then
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

  if(.not.associated(self%my_rays)) then
    error stop "rt_corrections%check: rt_t%my_rays is not associated"
  end if


END SUBROUTINE check_rt_corrections

END MODULE rt_internal_mod
