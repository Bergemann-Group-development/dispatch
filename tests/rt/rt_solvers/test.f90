!===============================================================================
!===============================================================================
MODULE radiation_mod
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod
  USE timer_mod
  USE ray_geometry_v3_mod
  implicit none
  type, public :: rt_data_t
    real, dimension(:,:,:,:), pointer :: opacity => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: source  => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: Q0      => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: Q       => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: dtau    => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: tau     => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: P       => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: error   => null()              ! (x,y,z,wl_bin)
    integer :: n(3), nrun
  contains
    procedure :: init => init_data
    procedure :: set => set_data
  end type
!-----------------------------------------------------------------------
  type :: internal_data_t
    real, dimension(:,:), pointer :: I_in  => null()                    ! (ray, wl_bin) incoming intensity.
    real, dimension(:,:), pointer :: I_out => null()                    ! leaving intensity
  end type
  !-----------------------------------------------------------------------
  type, public :: rt_internal_t
    class (ray_bndl_t),   pointer :: bundle
    class (rt_data_t),    pointer :: my_rt_data => null()
    class(internal_data_t), dimension(:), pointer :: packet_data => null() ! incoming intensities
    real, dimension(:,:), pointer :: I_in, I_out
    integer :: rt_type, nrun
    real :: dir_weight
  contains
    procedure :: init => init_internal
    !procedure :: solve => solve_1st
    procedure :: solve => solve_2nd1
    procedure :: copy_p
    procedure :: copy_p2
    procedure :: update_q
  end type
  !logical, parameter :: detailed_timer=.false.
  logical, parameter :: detailed_timer=.true.
  integer, save :: itimer_p = 0
  integer, save :: itimer_q = 0
CONTAINS

SUBROUTINE init_data(self, ng, nrun)
  class(rt_data_t) :: self
  integer :: ng(3), nrun
  !---------------------------------------------------------------------
  allocate (self%source(ng(1),ng(2),ng(3),nrun))
  allocate (self%q0(ng(1),ng(2),ng(3),nrun))
  allocate (self%q(ng(1),ng(2),ng(3),nrun))
  allocate (self%dtau(ng(1),ng(2),ng(3),nrun))
  allocate (self%tau(ng(1),ng(2),ng(3),nrun))
  allocate (self%P(ng(1),ng(2),ng(3),nrun))
  allocate (self%error(ng(1),ng(2),ng(3),nrun))
  self%n = ng
  self%q = 0.
  self%q0 = 0.
  self%nrun = nrun
END SUBROUTINE init_data

SUBROUTINE set_data(self)
  class(rt_data_t) :: self

  real, dimension(:,:), pointer :: ifw, irv
  ! --------------------------------------------------------------------
  real :: logtaumin, logtaumax, dlogtau,sigma, logtau, tmin,tmax,dtau1
  real, dimension(:,:,:), pointer :: tau
  integer:: n(3), x,run, k

  k=1
  n = self%n
  logtaumin = -6.
  logtaumax = +3.
  dlogtau = (logtaumax-logtaumin)/float(n(1)-1)
  tmin=-pi
  tmax=+pi
  dtau1=(tmax-tmin)/float(n(1)-1)
  sigma = 1.4
  do run =1, self%nrun
    do x = 1, n(1)
      logtau  = logtaumin+dlogtau*(x-1)
      self%tau(x,:,:,run) = 10.**logtau
      if (x > 1) self%dtau(x,:,:,run) = self%tau(x,:,:,run)-self%tau(x-1,:,:,run)
      select case (run)
      case (1)
        self%source(x,:,:,run)   = self%tau(x,:,:,run)**run
        self%q0(x,:,:,run)  = 0.5*exp(-self%tau(x,:,:,run))
      case (2)
        self%source(x,:,:,run)   = self%tau(x,:,:,run)**(run)
        self%q0(x,:,:,run)  = 2.-exp(-self%tau(x,:,:,run))
      case (3)
        self%source(x,:,:,run)   = self%tau(x,:,:,run)**run
        self%q0(x,:,:,run)  = 6.*self%tau(x,:,:,run)+3.*exp(-self%tau(x,:,:,run))
      case (5)
        self%source(x,:,:,run)   = (erf(logtau/sigma)+1)*0.5
        self%q0(x,:,:,run)  = 0.
        self%tau(1,:,:,run) = 1e-10
        self%source(1,:,:,run)   = (erf(alog10(self%tau(1,:,:,run))/sigma)+1)*0.5
        self%dtau(2,:,:,run) = self%tau(2,:,:,run)-self%tau(1,:,:,run)
      case (4)
        self%tau(x,:,:,run) = tmin+dtau1*(x-1)
        if (x > 1) self%dtau(x,:,:,run) = self%tau(x,:,:,run)-self%tau(x-1,:,:,run)
        self%source(x,:,:,run) = 1. + sin(k*self%tau(x,:,:,run))
        self%q0(x,:,:,run)     = -k**2*sin(k*self%tau(x,:,:,run))/(1.+k**2)
      end select
    end do
    !self%dtau(1,:,:,run) = self%dtau(2,:,:,run)*0.5
  end do
END SUBROUTINE set_data

SUBROUTINE init_internal(self, nrun)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet
  integer :: nrun, ir, ray, irm, run, c(3), pck
  real :: tau, k
  !---------------------------------------------------------------------
  k=1
  allocate(self%packet_data(self%bundle%n_packets))
  allocate(self%I_in(self%bundle%n_rays, nrun))
  allocate(self%I_out(self%bundle%n_rays, nrun))
  self%I_in = 0; self%I_out = 0
  self%nrun = nrun
  self%dir_weight = 0.5
  packet => self%bundle%rays
  ir = 1
  do ray = 1, self%bundle%n_packets
    irm = ir + packet%n_rays -1
    Self%packet_data(ray)%I_in  => Self%I_in(ir:irm,:)
    Self%packet_data(ray)%I_out => Self%I_out(ir:irm,:)
    packet => packet%next
    ir = irm + 1
  end do

  packet => self%bundle%rays
  do pck = 1, self%bundle%n_packets
    do run = 1, nrun
      do ray = 1, packet%n_rays
        select case(run)
        case(1)
          if (self%rt_type == 2) then
            c = packet%coords(:,1,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = -1. + tau + exp(-tau)
          else if (self%rt_type == 3) then
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = +1. + tau
          end if
        case(2)
          if (self%rt_type == 2) then
            c = packet%coords(:,1,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = 2. - 2.*tau + tau**2 - 2.*exp(-tau)
          else if (self%rt_type == 3) then
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = 2. + 2.*tau + tau**2
          end if
        case(3)
          if (self%rt_type == 2) then
            c = packet%coords(:,1,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = -6. + 6.*tau - 3.*tau**2 + tau**3 + 6.*exp(-tau)
          else if (self%rt_type == 3) then
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = +6. + 6.*tau + 3.*tau**2 + tau**3
          end if
        case(5)
          if (self%rt_type == 2) then
            self%packet_data(pck)%I_in(ray,run) = 0.
          else if (self%rt_type == 3) then
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%source(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = tau
          end if
        case(4)
          if (self%rt_type == 2) then
            c = packet%coords(:,1,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = 1. + (-k*cos(k*tau) + &
              sin(k*tau))/(1. + k**2)
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_out(ray,run) = 1. + (-k*cos(k*tau) + &
              sin(k*tau))/(1. + k**2)
          else if (self%rt_type == 3) then
            c = packet%coords(:,packet%n_points,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_in(ray,run) = 1. + (+k*cos(k*tau) + &
              sin(k*tau))/(1. + k**2)
            c = packet%coords(:,1,ray)
            tau = self%my_rt_data%tau(c(1),c(2),c(3),run)
            self%packet_data(pck)%I_out(ray,run) = 1. + (+k*cos(k*tau) + &
              sin(k*tau))/(1. + k**2)
          end if
        end select
      end do
    end do
  end do

END SUBROUTINE init_internal

SUBROUTINE update_q (self, rays, q, dummy)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt, run
  real, dimension(:,:,:), intent(in) :: q
  real, dimension(:,:,:,:), pointer :: qq
  real, dimension(:,:,:), pointer, optional :: dummy
  integer, dimension(:,:,:), pointer :: cell
  real :: dirw
  !---------------------------------------------------------------------
  call trace_begin ('rt_internal_t_%update_q', 4, itimer=itimer_q)
  dirw = self%dir_weight
  qq => self%my_rt_data%q
  cell => rays%coords
  do ray = 1, size(q,1)
    do pt = 1, size(q,3)
      p = cell(:,pt,ray)
      do run=1,self%nrun
        qq(p(1),p(2),p(3),run) = &
        qq(p(1),p(2),p(3),run) + dirw*q(ray,run,pt)
      end do
    end do
  end do
  call trace_end (itimer_q)
END SUBROUTINE update_q

SUBROUTINE copy_p (self, rays, in, out)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt, run
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  call trace_begin ('copy_p', 4, itimer=itimer_p)
  cell => rays%coords
  do ray = 1, size(out,1)
    do pt = 1, size(out,3)
      p = cell(:,pt,ray)
      !DIR$ SIMD
      do run=1,self%nrun
        out(ray,run,pt) = in(p(1),p(2),p(3),run)
      end do
    end do
  end do
  call trace_end (itimer_p)
END SUBROUTINE copy_p

SUBROUTINE copy_p2 (self, rays, in1, out1, in2, out2)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt, run
  real, dimension(:,:,:)   , pointer, intent(out) :: out1, out2
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in1, in2
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  call trace_begin ('copy_p2', 4, itimer=itimer_p)
  cell => rays%coords
  do ray = 1, size(out1,1)
    do pt = 1, size(out1,3)
      p = cell(:,pt,ray)
      !DIR$ SIMD
      do run=1,self%nrun
        out1(ray,run,pt) = in1(p(1),p(2),p(3),run)
        out2(ray,run,pt) = in2(p(1),p(2),p(3),run)
      end do
    end do
  end do
  call trace_end (itimer_p)
END SUBROUTINE copy_p2

SUBROUTINE solve_1st (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet
  real, dimension(:,:,:), pointer :: ex0, ds, dtau
  real, dimension(:,:,:), pointer :: s, q, Ifw, ropac
  real, dimension(:,:), pointer :: I_in, I_out
  integer :: pt, p(3), pck, nwb, ipk, ray, id
  integer, save :: itimer(4)=0
  !---------------------------------------------------------------------
  if (detailed_timer) then
    itimer = 0
    itimer_q = 0
    call timer%init
  end if
  call tic
  nwb = self%nrun
  packet => self%bundle%rays
  do pck = 1, self%bundle%n_packets
    if (detailed_timer) &
      call trace_begin('solve_1st', itimer=itimer(1))
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
        dtau(ray,:,pt) = rtd%dtau(p(1),p(2),p(3),:)
      end do
    end do
    if (detailed_timer) then
      call trace_end (itimer(1))
      call trace_begin('solve_1st', itimer=itimer(2))
    end if

    I_in = self%packet_data(pck)%I_in

    do pt = 2, npt
      ex0(:,:,pt) = exp(-dtau(:,:,pt))                                  ! first exponential coeff
      ds(:,:,pt)  = (1.-ex0(:,:,pt))*s(:,:,pt-1)                          ! Source factor
    end do
      ds(:,:,1)  = (1.-ex0(:,:,1))*s(:,:,1)                          ! Source factor

    if (detailed_timer) then
      call trace_end (itimer(2))
      call trace_begin('solve_1st', itimer=itimer(3))
    end if

      !-----------------------------------------------------------------
      ! Forward / reverse solutions
      !-----------------------------------------------------------------
    select case(self%rt_type)
    case(2)                                                             ! forward
      Ifw(:,:,1) = I_in                                                 ! Intensity after the first point
      do pt=2,npt
        Ifw(:,:,pt) = Ifw(:,:,pt-1) * ex0(:,:,pt) + ds(:,:,pt)
      end do
      q = Ifw - s
      self%packet_data(pck)%I_out = Ifw(:,:,npt) !* ex0(:,:,npt) + ds(:,:,npt)

    case(3)                                                             ! reverse
      Ifw(:,:,npt) = I_in
      do pt = npt-1, 1, -1
        Ifw(:,:,pt) = Ifw(:,:,pt+1) * ex0(:,:,pt) + ds(:,:,pt)
      end do
      q  = Ifw - s
      self%packet_data(pck)%I_out = Ifw(:,:,1) !* ex0(:,:,1) + ds(:,:,1)
    end select
    if (detailed_timer) then
      call trace_end (itimer(3))
    end if
    call self%update_q (packet, q, ropac)
    deallocate(ex0, ds, s, dtau, q, Ifw, I_in, ropac)
    end associate
    packet => packet%next
  end do
  nullify (packet)
  call timer%toc('solve 1st', product(self%my_rt_data%n)*self%nrun)
  if (detailed_timer) then
    call timer%print
  end if
END SUBROUTINE solve_1st

!===============================================================================
!> This solver formulation uses loops over points, and precomputes  ex0 and ex1
!> factors, but not ex2.
!===============================================================================
SUBROUTINE solve_2nd1 (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet => null()
  integer, dimension(:,:,:), pointer   :: cell => null()
  class(rt_data_t), pointer :: rtd
  real, dimension(:,:,:), pointer :: ex0, ds1, d1s, d2s
  real, dimension(:,:,:), pointer :: s, dtau, q, ropac
  real, dimension(:,:,:,:), pointer :: ss, ddtau
  real, dimension(:,:), pointer :: I_in, ctau
  integer :: pt, npt, nrays, it, ray, wl, p(3), pck, nwb, ipk, iwb
  integer, save :: itimer(4)=0
  !---------------------------------------------------------------------
  if (detailed_timer) then
    itimer = 0
    itimer_p = 0
    itimer_q = 0
    call timer%init
  end if
  call tic
  nwb = self%nrun                                                       ! wl bins represent different tests here
  packet => self%bundle%rays

  rtd => self%my_rt_data
  do pck = 1, self%bundle%n_packets
    if (detailed_timer) then
      call trace_begin('solve_2nd1: begin', itimer=itimer(1))
    end if
    cell => packet%coords
    nrays = packet%n_rays
    npt =  packet%n_points
    ddtau => rtd%dtau
    ss    => rtd%source

    allocate( dtau(nrays, nwb, npt))
    allocate(ex0  (nrays, nwb, npt))
    allocate(ds1  (nrays, nwb, npt))
    allocate(d1s  (nrays, nwb, npt))
    allocate(d2s  (nrays, nwb, npt))
    allocate(s    (nrays, nwb, npt))
    allocate(q    (nrays, nwb, npt))
    allocate(ctau (nrays, nwb))
    allocate(I_in (nrays, nwb))

    if (detailed_timer) then
      call trace_end (itimer(1))
    end if

    call self%copy_p (packet, ss, s)
    call self%copy_p (packet, ddtau, dtau)

    if (detailed_timer) then
      call trace_begin('solve_2nd1: ex0,ex1', itimer=itimer(3))
    end if

    I_in = self%packet_data(pck)%I_in

!#if defined __INTEL_COMPILER
!      ex0 = 2**(-dtau*1.44269504)                                       ! first exponential coeff
!#else
    do pt = 2, npt
      ex0(:,:,pt) = exp(-dtau(:,:,pt))                                   ! first exponential coeff
      ds1(:,:,pt) = (s(:,:,pt) - s(:,:,pt-1)) / dtau(:,:,pt)             ! half-point centered dS/dtau
    end do
!#endif
    if (detailed_timer) then
      call trace_end (itimer(3))
      call trace_begin('solve_2nd1: transfer', itimer=itimer(4))
    end if

    !-------------------------------------------------------------------
    !  dS/dtau and d2S/dtau2
    !-------------------------------------------------------------------
    do pt=2, npt-1
      ctau = 1. / (dtau(:,:,pt+1) + dtau(:,:,pt))
      d1s(:,:,pt) = (ds1(:,:,pt+1) * dtau(:,:,pt) &
                  + ds1(:,:,pt) * dtau(:,:,pt+1) ) * ctau               ! dS/dtau, 2nd. order accurate
      d2s(:,:,pt) = (ds1(:,:,pt+1) - ds1(:,:,pt))* 2. * ctau            ! d2S/dtau2, 2nd. order accurate
    end do
    d2s(:,:,1)  = d2s(:,:,2)                                            ! assume constant d2S/dtau2
    d1s(:,:,1)  = d1s(:,:,2)    - dtau(:,:,2)  * d2s(:,:,2)             ! extrapolate to 1
    d2s(:,:,npt) = d2s(:,:,npt-1)                                       ! assume constant d2S/dtau2
    d1s(:,:,npt) = d1s(:,:,npt-1) + dtau(:,:,npt) * d2s(:,:,npt)        ! extrapolate to n

    !-------------------------------------------------------------------
    !  Internal forward or reverse solutions
    !-------------------------------------------------------------------
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
      error stop "unknown rt type"
    end select
    if (detailed_timer) then
      call trace_end (itimer(4))
    end if

    call self%update_q (packet, q, ropac)
    deallocate(ex0, ds1, d1s, d2s, ctau, s, dtau, q, I_in)
    packet => packet%next
  end do
  nullify (packet,rtd,cell)
  call timer%toc('solve 2nd1', product(self%my_rt_data%n)*self%nrun)
  if (detailed_timer) then
    call timer%print
  end if
END SUBROUTINE solve_2nd1

SUBROUTINE solve_2nd12 (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet => null()
  integer, dimension(:,:,:), pointer   :: cell => null()
  class(rt_data_t), pointer :: rtd
  real, dimension(:,:,:), pointer :: ex0, ds1, d1s, d2s
  real, dimension(:,:,:), pointer :: s, dtau, q, ropac, ctau
  real, dimension(:,:,:,:), pointer :: ss, ddtau
  real, dimension(:,:), pointer :: I_in
  integer :: pt, npt, nrays, it, ray, wl, p(3), pck, nwb, ipk, iwb
  integer, save :: itimer(4)=0
  !---------------------------------------------------------------------
  if (detailed_timer) then
    itimer = 0
    itimer_p = 0
    itimer_q = 0
    call timer%init
  end if
  call tic
  nwb = self%nrun                                                       ! wl bins represent different tests here
  packet => self%bundle%rays

  rtd => self%my_rt_data
  do pck = 1, self%bundle%n_packets
    if (detailed_timer) then
      call trace_begin('solve_2nd1-2: begin', itimer=itimer(1))
    end if
    cell => packet%coords
    nrays = packet%n_rays
    npt =  packet%n_points
    ddtau => rtd%dtau
    ss    => rtd%source

    allocate( dtau(nrays, nwb, npt))
    allocate(ex0  (nrays, nwb, npt))
    allocate(ds1  (nrays, nwb, npt))
    allocate(d1s  (nrays, nwb, npt))
    allocate(d2s  (nrays, nwb, npt))
    allocate(s    (nrays, nwb, npt))
    allocate(q    (nrays, nwb, npt))
    allocate(ctau (nrays, nwb,npt-2))
    allocate(I_in (nrays, nwb))

    if (detailed_timer) then
      call trace_end (itimer(1))
    end if

    call self%copy_p (packet, ss, s)
    call self%copy_p (packet, ddtau, dtau)

    if (detailed_timer) then
      call trace_begin('solve_2nd1-2: ex0,ex1', itimer=itimer(3))
    end if

    I_in = self%packet_data(pck)%I_in

!#if defined __INTEL_COMPILER
!      ex0 = 2**(-dtau*1.44269504)                                       ! first exponential coeff
!#else
    !do pt = 2, npt
    !  ex0(:,:,pt) = exp(-dtau(:,:,pt))                                   ! first exponential coeff
    !  ds1(:,:,pt) = (s(:,:,pt) - s(:,:,pt-1)) / dtau(:,:,pt)             ! half-point centered dS/dtau
    !end do
    ex0(:,:,2:npt) = exp(-dtau(:,:,2:npt))                              ! first exponential coeff
    ds1(:,:,2:npt) = (s(:,:,2:npt) - s(:,:,1:npt-1)) / dtau(:,:,2:npt)  ! half-point centered dS/dtau
!#endif
    if (detailed_timer) then
      call trace_end (itimer(3))
      call trace_begin('solve_2nd1-2: transfer', itimer=itimer(4))
    end if

    !-------------------------------------------------------------------
    !  dS/dtau and d2S/dtau2
    !-------------------------------------------------------------------

    ctau = 1. / (dtau(:,:,3:npt) + dtau(:,:,2:npt-1))
    d1s(:,:,2:npt-1) = (ds1(:,:,3:npt) * dtau(:,:,2:npt-1) &
                     + ds1(:,:,2:npt-1) * dtau(:,:,3:npt) ) * ctau      ! dS/dtau, 2nd. order accurate
    d2s(:,:,2:npt-1) = (ds1(:,:,3:npt) - ds1(:,:,2:npt-1))* 2. * ctau   ! d2S/dtau2, 2nd. order accurate
    d2s(:,:,1)  = d2s(:,:,2)                                            ! assume constant d2S/dtau2
    d1s(:,:,1)  = d1s(:,:,2)    - dtau(:,:,2)  * d2s(:,:,2)             ! extrapolate to 1
    d2s(:,:,npt) = d2s(:,:,npt-1)                                       ! assume constant d2S/dtau2
    d1s(:,:,npt) = d1s(:,:,npt-1) + dtau(:,:,npt) * d2s(:,:,npt)        ! extrapolate to n

    !-------------------------------------------------------------------
    !  Internal forward or reverse solutions
    !-------------------------------------------------------------------
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
      error stop "unknown rt type"
    end select
    if (detailed_timer) then
      call trace_end (itimer(4))
    end if

    call self%update_q (packet, q, ropac)
    deallocate(ex0, ds1, d1s, d2s, ctau, s, dtau, q, I_in)
    packet => packet%next
  end do
  nullify (packet,rtd,cell)
  call timer%toc('solve 2nd1-2', product(self%my_rt_data%n)*self%nrun)
  if (detailed_timer) then
    call timer%print
  end if
END SUBROUTINE solve_2nd12

!===============================================================================
!> This solver formulation uses loops over points, and precomputes  ex0, ex1,
!> and ex2 factors.
!===============================================================================
SUBROUTINE solve_2nd2(self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet => null()
  integer, dimension(:,:,:), pointer   :: cell => null()
  class(rt_data_t), pointer :: rtd
  real, dimension(:,:,:), pointer :: ex0, ex1, ex2, ds1dtau, dsdtau
  real, dimension(:,:,:), pointer :: d2sdtau2, s, dtau, q
  real, dimension(:,:,:,:), pointer :: ss, ddtau

  real, dimension(:,:), pointer :: I_in, ctau
  integer :: pt, npt, run, nrun, nrays, p(3), nr, pck, ray
  integer, save :: itimer(6)=0
  !---------------------------------------------------------------------
  if (detailed_timer) then
    itimer = 0
    itimer_p = 0
    itimer_q = 0
    call timer%init
  end if
  call tic
  nrun = self%nrun                                                       ! wl bins represent different tests here
  packet => self%bundle%rays

  rtd => self%my_rt_data

  do pck = 1, self%bundle%n_packets
    if (detailed_timer) then
      call trace_begin('solve_2nd2: begin', itimer=itimer(1))
    end if
    cell => packet%coords
    nrays = packet%n_rays
    npt =  packet%n_points
    ddtau => rtd%dtau
    ss    => rtd%source

    allocate( dtau(nrays, nrun, npt))  !
    allocate(ex0  (nrays, nrun, npt))
    allocate(ex1  (nrays, nrun, npt))
    allocate(ex2  (nrays, nrun, npt))
    allocate(ds1dtau(nrays, nrun, npt))
    allocate(dsdtau(nrays, nrun, npt))
    allocate(d2sdtau2(nrays, nrun, npt))
    allocate(s(nrays, nrun, npt))  !
    allocate(q(nrays, nrun, npt))
    allocate(I_in (nrays, nrun))
    allocate(ctau (nrays, nrun))

    if (detailed_timer) then
      call trace_end (itimer(1))
      call trace_begin('solve_2nd2: s, dtau', itimer=itimer(2))
    end if
    do ray = 1, nrays
      do pt = 1, npt
        p = cell(:,pt,ray)
        !DIR$ SIMD
        do run=1,nrun
             s(ray,run,pt) =    ss(p(1),p(2),p(3),run)
        end do
      end do
    end do
    do ray = 1, nrays
      do pt = 1, npt
        p = cell(:,pt,ray)
        !DIR$ SIMD
        do run=1,nrun
          dtau(ray,run,pt) = ddtau(p(1),p(2),p(3),run)
        end do
      end do
    end do

    if (detailed_timer) then
      call trace_end (itimer(2))
      call trace_begin('solve_2nd2: ex0,ex1', itimer=itimer(3))
    end if

    do pt = 2, npt
      ex0(:,:,pt) = exp(-dtau(:,:,pt))
      ex1(:,:,pt) = 1.-ex0(:,:,pt)
      ex2(:,:,pt) = ex1(:,:,pt)-dtau(:,:,pt)*ex0(:,:,pt)
      ds1dtau(:,:,pt) = (s(:,:,pt)-s(:,:,pt-1))/dtau(:,:,pt)
    end do
    if (detailed_timer) then
      call trace_end (itimer(3))
      call trace_begin('solve_2nd2: transfer', itimer=itimer(4))
    end if

    do pt = 2, npt-1
      ctau = 1./(dtau(:,:,pt+1)+dtau(:,:,pt))
      dsdtau(:,:,pt) = (ds1dtau(:,:,pt+1)*dtau(:,:,pt)+ds1dtau(:,:,pt)*dtau(:,:,pt+1))*ctau
      d2sdtau2(:,:,pt) = (ds1dtau(:,:,pt+1)-ds1dtau(:,:,pt))*2.*ctau
    end do
    d2sdtau2(:,:,1) = d2sdtau2(:,:,2)
    dsdtau(:,:,1) = dsdtau(:,:,2) - dtau(:,:,2)*d2sdtau2(:,:,2)
    d2sdtau2(:,:,npt) = d2sdtau2(:,:,npt-1)
    dsdtau(:,:,npt) = dsdtau(:,:,npt-1)+dtau(:,:,npt)*d2sdtau2(:,:,npt)

    I_in = self%packet_data(pck)%I_in
    select case(self%rt_type)
    case(2) !forward
      q(:,:,1) = I_in - s(:,:,1)
      do pt=2,npt
      q(:,:,pt)=q(:,:,pt-1)*ex0(:,:,pt)-dsdtau(:,:,pt)*ex1(:,:,pt)+d2sdtau2(:,:,pt)*ex2(:,:,pt)
      end do
      self%packet_data(pck)%I_out = q(:,:,npt) + s(:,:,npt)
    case(3) !reverse
      q(:,:,npt)=I_in-s(:,:,npt)
      do pt=npt-1,1,-1
      q(:,:,pt)=q(:,:,pt+1)*ex0(:,:,pt+1)+dsdtau(:,:,pt)*ex1(:,:,pt+1)+d2sdtau2(:,:,pt)*ex2(:,:,pt+1)
      end do
      self%packet_data(pck)%I_out = q(:,:,1) + s(:,:,1)
    end select
    if (detailed_timer) then
      call trace_end (itimer(4))
    end if

    call self%update_q (packet, q)
    deallocate(dtau,ex0,ex1,ex2,ds1dtau,dsdtau,d2sdtau2,s,q,I_in,ctau)
    packet => packet%next
  end do
  nullify (packet,rtd,cell)
  call timer%toc('solve 2nd2', product(self%my_rt_data%n)*self%nrun)
  if (detailed_timer) then
    call timer%print
  end if
END SUBROUTINE solve_2nd2

!===============================================================================
!> This solver formulation uses index notation
!===============================================================================
SUBROUTINE solve_2nd3 (self)
  class(rt_internal_t) :: self
  class(ray_t), pointer :: packet => null()
  integer, dimension(:,:,:), pointer   :: cell => null()
  class(rt_data_t), pointer :: rtd
  real, dimension(:,:,:), pointer :: ex0, ds1, d1s, d2s
  real, dimension(:,:,:), pointer :: s, dtau, q, ropac, ctau
  real, dimension(:,:), pointer :: I_in
  real, dimension(:,:,:,:), pointer :: ss, ddtau
  integer :: pt, npt, nrays, it, ray, wl, p(3), pck, iwb, nwb, ipk
  integer, save :: itimer(4)=0
  !---------------------------------------------------------------------
  if (detailed_timer) then
    itimer = 0
    itimer_p = 0
    itimer_q = 0
    call timer%init
  end if
  call tic

  nwb = self%nrun                                                       ! wl bins represent different tests here
  packet => self%bundle%rays

  rtd => self%my_rt_data
  do pck = 1, self%bundle%n_packets
    if (detailed_timer) then
      call trace_begin('solve_2nd3: begin', itimer=itimer(1))
    end if
    cell => packet%coords
    nrays = packet%n_rays
    npt =  packet%n_points
    ddtau => rtd%dtau
    ss    => rtd%source

    allocate( dtau(nrays, nwb, npt))
    allocate(ex0  (nrays, nwb, npt))
    allocate(ds1  (nrays, nwb, npt))
    allocate(d1s  (nrays, nwb, npt))
    allocate(d2s  (nrays, nwb, npt))
    allocate(s    (nrays, nwb, npt))
    allocate(q    (nrays, nwb, npt))
    allocate(ctau (nrays, nwb,npt-2))
    allocate(I_in (nrays, nwb))

    if (detailed_timer) then
      call trace_end (itimer(1))
    end if
    call self%copy_p2 (packet, ss, s, ddtau, dtau)

    if (detailed_timer) then
      call trace_begin('solve_2nd3: s, dtau', itimer=itimer(2))
    end if
    do ray = 1, nrays
      do pt = 1, npt
        p = cell(:,pt,ray)
        !DIR$ SIMD
        do iwb=1,nwb
             s(ray,iwb,pt) = ss(p(1),p(2),p(3),iwb)
          dtau(ray,iwb,pt) = ddtau(p(1),p(2),p(3),iwb)
        end do
      end do
    end do
    if (detailed_timer) then
      call trace_end (itimer(2))
      call trace_begin('solve_2nd3: ex0,ex1', itimer=itimer(3))
    end if

!#if defined __INTEL_COMPILER
!      ex0 = 2**(-dtau*1.44269504)                                       ! first exponential coeff
!#else
    ex0(:,:,2:npt) = exp(-dtau(:,:,2:npt))                              ! first exponential coeff
    ds1(:,:,2:npt) = (s(:,:,2:npt) - s(:,:,1:npt-1)) / dtau(:,:,2:npt)  ! half-point centered dS/dtau
!#endif
    if (detailed_timer) then
      call trace_end (itimer(3))
      call trace_begin('solve_2nd3: transfer', itimer=itimer(4))
    end if

    !-------------------------------------------------------------------
    !  dS/dtau and d2S/dtau2
    !-------------------------------------------------------------------
    ctau = 1. / (dtau(:,:,3:npt) + dtau(:,:,2:npt-1))
    d1s(:,:,2:npt-1) = (ds1(:,:,3:npt) * dtau(:,:,2:npt-1) &
                     + ds1(:,:,2:npt-1) * dtau(:,:,3:npt) ) * ctau      ! dS/dtau, 2nd. order accurate
    d2s(:,:,2:npt-1) = (ds1(:,:,3:npt) - ds1(:,:,2:npt-1))* 2. * ctau   ! d2S/dtau2, 2nd. order accurate
    d2s(:,:,1)  = d2s(:,:,2)                                            ! assume constant d2S/dtau2
    d1s(:,:,1)  = d1s(:,:,2)    - dtau(:,:,2)  * d2s(:,:,2)             ! extrapolate to 1
    d2s(:,:,npt) = d2s(:,:,npt-1)                                       ! assume constant d2S/dtau2
    d1s(:,:,npt) = d1s(:,:,npt-1) + dtau(:,:,npt) * d2s(:,:,npt)        ! extrapolate to n

    !-------------------------------------------------------------------
    !  Internal forward or reverse solutions
    !-------------------------------------------------------------------
    I_in = self%packet_data(pck)%I_in
    select case(self%rt_type)
    case(2)                                                             ! forward
      q  (:,:,1) = I_in - s(:,:,1)                                      ! initial forward Q
       ! q(:,:,2:npt) = q(:,:,1:npt-1) * ex0(:,:,2:npt) - d1s(:,:,2:npt) &
       !           * (1. - ex0(:,:,2:npt)) + d2s(:,:,2:npt) * &
       !           (1. - ex0(:,:,2:npt) - dtau(:,:,2:npt)*ex0(:,:,2:npt))
      do pt=2,npt
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
      error stop "unknown rt type"
    end select
    if (detailed_timer) then
      call trace_end (itimer(4))
    end if

    call self%update_q (packet, q, ropac)
    deallocate(ex0, ds1, d1s, d2s, ctau, s, dtau, q, I_in)
    packet => packet%next
  end do
  nullify (packet,rtd,cell)
  call timer%toc('solve 2nd3', product(self%my_rt_data%n)*self%nrun)
  if (detailed_timer) then
    call timer%print
  end if
END SUBROUTINE solve_2nd3

END MODULE radiation_mod


PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod
  USE timer_mod
  USE radiation_mod
  USE ray_geometry_v3_mod
  USE io_mod
  implicit none
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  class(rt_internal_t), pointer :: rt_obj_fw, rt_obj_rv
  class(rt_data_t), pointer :: rt_data
  integer :: n(3), l(3), u(3), ng(3), nwb
  integer :: dir, bnd, ray, ir, irm, x, y, z, run, nrun,i, c(3)
  class(ray_t), pointer :: packet
  real :: cput, cpu(2), dtime, maxerr

  !call tic
  dir = 1; bnd = 1  !test just on dir 1, bundle 1 now
  call omp%init
  call io%init
  call timer%init
  n = 100   ! active inner cells
  l = 1    ! li
  u = l+n-1! ui
  ng = u+2
  nwb = 1  ! n wl bins
  nrun = 4

  call tic
  call ray_geometry%init(1,l,n)
  call timer%toc('ray_geometry%init')

  allocate (rt_data)
  call tic
  call rt_data%init(n,nrun)
  call timer%toc('rt_data%init')
  call tic
  call rt_data%set
  call timer%toc('rt_data%set')

  allocate (rt_obj_fw, rt_obj_rv)
  rt_obj_fw%rt_type = 2
  rt_obj_rv%rt_type = 3
  rt_obj_fw%my_rt_data => rt_data
  rt_obj_rv%my_rt_data => rt_data
  dir = ray_geometry%head%ray_list%id
  bnd = ray_geometry%head%ray_list%bundle%id
  call ray_geometry%head%get_my_bundle(rt_obj_fw%bundle, dir,bnd)
  call ray_geometry%head%get_my_bundle(rt_obj_rv%bundle, dir,bnd)

  call tic
  call rt_obj_fw%init(nrun)
  call rt_obj_rv%init(nrun)
  call timer%toc('rt_obj_xx%init')

  rt_data%q = 0.0
  call rt_obj_fw%solve
  rt_data%q = 0.0
  call rt_obj_fw%solve
  call rt_obj_rv%solve
  rt_data%P = rt_data%source + rt_data%q
  rt_data%error = (rt_data%q-rt_data%q0)/rt_data%P
  do run =1,nrun
    print *,'q0 run',run,maxval(abs(rt_data%error(:,:,:,run)))
  end do

    rt_data%q = 0.0
  call solve_2nd12(rt_obj_fw)
  call solve_2nd12(rt_obj_Rv)
  rt_data%P = rt_data%source + rt_data%q
  rt_data%error = (rt_data%q-rt_data%q0)/rt_data%P
  do run =1,nrun
    print *,'q0 run',run,maxval(abs(rt_data%error(:,:,:,run)))
  end do

  rt_data%q = 0.0
  call solve_2nd2(rt_obj_fw)
  call solve_2nd2(rt_obj_Rv)
  rt_data%P = rt_data%source + rt_data%q
  rt_data%error = (rt_data%q-rt_data%q0)/rt_data%P
  do run =1,nrun
    print *,'q0 run',run,maxval(abs(rt_data%error(:,:,:,run)))
  end do

  rt_data%q = 0.0
  call solve_2nd3(rt_obj_fw)
  call solve_2nd3(rt_obj_Rv)
  rt_data%P = rt_data%source + rt_data%q
  rt_data%error = (rt_data%q-rt_data%q0)/rt_data%P
  do run =1,nrun
    print *,'q0 run',run,maxval(abs(rt_data%error(:,:,:,run)))
  end do

!do run =1, nrun
!print *, 'forward I_in: ', rt_obj_fw%I_in(1,run)
!print *,'  i run   X   Y   Z     tau         S           P           Q          Q0          err'
  !do i=1,n(1)
    !c = rt_obj_fw%bundle%rays%coords(:,i,1)
    !print '(5i4,5e12.4,f12.6)',i,run, c(1), c(2), c(3), rt_data%tau(i,1,1,run), &
    !rt_data%source(i,1,1,run), rt_data%P(i,1,1,run), rt_data%q(i,1,1,run), rt_data%q0(i,1,1,run), rt_data%error(i,1,1,run)
  !end do
!end do

!cput=dtime(cpu)
! do iter=1,niter
!   call transfer (n,dtau,s,ifw1,irvN,q,irv1,ifwN)
! end do
! cput=dtime(cpu)
! print *,'cpu =',cput
! niter = niter*2./cput                                                  ! number of iterations for 2 seconds CPU

!  run = 1

!  call test1(run, rt_obj%my_rt_data%source, rt_obj%my_rt_data%q0, &
!             rt_obj%my_rt_data%dtau, rt_obj%i_in)

!SUBROUTINE test1(run,source,q0,dtau, ifw, irv)
!  integer :: run


!END SUBROUTINE test1

 !call timer%toc('all', product(n)*nrun)

END
