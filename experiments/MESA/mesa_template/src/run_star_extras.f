! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

! For now, just useing file handles 96 and 97. TODO: fix.
      module DISPATCH_COMMUNICATION
         use const_def
         implicit none ! TODO: import these form the same compilation unit as mesa_mod.f90 to ensure compiler-enforced compatibility
         type struct_from_dispatch
            sequence
            real :: time = 0 ! time since formation of 2.nd larson core (s)
            real :: acc_mass = 0 ! accretedMass (g/s)
         end type struct_from_dispatch
         type struct_from_mesa
            sequence
            real :: time = 0 ! time since formation of 2.nd larson core (s)
            real :: wind_mass = 0 ! accretedMass (g/s)
            logical :: supernova = .false. !! booooommmm !!
         end type struct_from_mesa

         
         TYPE(struct_from_dispatch), DIMENSION(:), allocatable :: acc_hist
         TYPE(integer) :: acc_len = 0
         integer :: unit1, unit2 !pipes to and from mesa, respectively

         
         contains
         
         subroutine init()
            use utils_lib, only: alloc_iounit
            integer :: ierr
            unit1 = alloc_iounit(ierr)
            open (unit1, file="pipe1", form='formatted', &
               status='unknown', action='read')
            unit2 = alloc_iounit(ierr)
            open (unit2, file="pipe2", form='formatted', &
               status='unknown', action='write')
            allocate(acc_hist(1))
            acc_len = 1
         end subroutine init


         function acc_rate(tstart, tend)
            real(dp) :: tstart, tend, acc_rate, lb, mp, ub
            integer:: i, istart, iend

            ! Make sure our accretion history covers thre requested timespan [tstart:tend]
            do while (acc_hist(acc_len)%time < tend)
               acc_len = acc_len + 1
               ! If no more capacity in accretion_hist, resize using doubling strategy
               if (acc_len > size(acc_hist)) then
                  call resize_acc_hist(size(acc_hist)*2)
               endif
               acc_hist(acc_len) = request_rate()
            end do

            ! Find the start and end of the interval in the history array
            ! Due to retries, we have no guarantee that this is the end
            i = acc_len
            do while (acc_hist(i)%time > tstart)  
               if (acc_hist(i)%time > tend) iend = i
               i = i - 1
               if (i == 1) exit
            end do
            istart = i
            ! now istart is last element which happended before tstart, 
            ! and iend the first element after tend, 
            ! ie acc_hist(istart)%time < tstart < tend < acc_hist(iend)%time


            ! Estimate rate as linear sum (Thanks old code)
            acc_rate = 0.0
            do i=istart, iend-1
               lb = max(acc_hist(i)%time,tstart)
               mp = max(min(0.5*(acc_hist(i)%time+acc_hist(i+1)%time),tend),tstart)
               ub = min(acc_hist(i+1)%time,tend)
               acc_rate = acc_rate + (mp-lb)*acc_hist(i)%acc_mass+(ub-mp)*acc_hist(i+1)%acc_mass
            enddo

            acc_rate = acc_rate / (tend - tstart)
         end function acc_rate

         type(struct_from_dispatch) function request_rate()
            ! namelist /query/ request_rate
            read (unit1,*) request_rate
         end function request_rate


         ! simple resize function using a temp array. This copies to a temp array and then back
         ! which is not neccecary. I think SAVE or pointers could be used.
         subroutine resize_acc_hist(newSize)
            INTEGER, INTENT(IN) :: newSize
            type(struct_from_dispatch), DIMENSION(:), ALLOCATABLE :: temp
            
            allocate(temp(size(acc_hist)))
            
            temp = acc_hist
            
            deallocate(acc_hist)
            allocate(acc_hist(newSize))
            
            acc_hist(1:size(temp)) = temp
            !call allocate(tmp(size_old+n_enlarge))
            !tmp(1:size_old) = array(1:size_old)
            !call move_alloc(tmp, array)

         end subroutine resize_acc_hist
            
            
         


      end module DISPATCH_COMMUNICATION
 


      module run_star_extras

         use star_lib
         use star_def
         use const_def
         
         implicit none
         
         integer :: time0, time1, clock_rate

         
         contains
         
         
         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            !s% use_extras_startup = .false.
            s% use_other_energy = .true.
            s% use_other_adjust_mdot = .true.
            s% extras_startup => extras_startup
            s% other_energy => accretion_other_energy
            s% other_adjust_mdot => other_dispatch_adjust_mdot
         end subroutine extras_controls

         integer function extras_startup(id, restart, ierr)
            use DISPATCH_COMMUNICATION
         
            implicit none

            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr

            call init
         end function extras_startup

         subroutine other_dispatch_adjust_mdot(id, ierr)
            use DISPATCH_COMMUNICATION
            use star_def

            implicit none

            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type(star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) then 
                return
            end if

            s% mstar_dot = acc_rate(s% time_old, s% time_old + s% dt)
            write(*,*) s% time_old, s% time_old + s% dt, s% mstar_dot

         end subroutine other_dispatch_adjust_mdot


         
         subroutine accretion_other_energy(id, ierr)

            use crlibm_lib, only: safe_log10_cr
            use star_def


            implicit none

            type (star_info), pointer :: s

            integer, intent(in) :: id
            integer, intent(out) :: ierr
            real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul
            real(dp) :: Lacc, extra_heat_var
            real(dp), parameter :: A = -1.977
            real(dp), parameter :: B = -3.078
            real(dp), parameter :: zp = log(2.792e-6)
            real(dp), parameter :: secyer = 3.1558149984d7 ! seconds per year
            real(dp), parameter :: msun   = 1.9892e33      ! solar mass (g)
            real(dp) :: omega
            integer :: k
            ierr  = 0

            call star_ptr(id, s, ierr)
            if (ierr /= 0) then ! OOPS
               return
            end if

            mass = s% mstar
            mdot = s% mstar_dot

            if (mdot == 0.) return

            radius = s% r(1)
            ! ----------------------------------------------------------
            ! x_ctrl(1) is the default alpha parameter from inlist.
            !alpha = s% x_ctrl(1)

            ! ----------------------------------------------------------
            ! Time varying alpha, alpha = 0.25 if mdot > 1e-5
            ! crossover for alpha in range mdot [1e-5, 1e-6]
            ! settle at alpha = 0.01 when mdot < 1e-6
            mdot_msun = mdot*(secyer/msun)
            !write(*,*) "Mdot (msun)", mdot_msun
            if (mdot_msun >= 1e-5) then
            alpha = 0.25
            else if (mdot_msun <= 1e-6) then
            alpha = 0.01
            else
            alpha = exp(B+A*tanh(zp-log(mdot_msun)))
            endif

            mcumul = 0.
            k = 0
            Lacc = 0.5 * (s% cgrav(1) * (mdot) * (mass)) / radius
            extra_heat_var = Lacc * alpha
            do while ((mcumul < mdot*s% dt .or. extra_heat_var / max(mcumul,1.) > 1e4) .and. mcumul < 0.8 * mass)
            k = k + 1
            mcumul = mcumul + s% dm(k)
            enddo

            extra_heat_var = Lacc * alpha / mcumul
            s% extra_heat(1:k) = extra_heat_var

            !omega = sqrt((s% cgrav(1) * mass)/(2*radius)**3)
            !s% omega(1:k) = omega
            !write(*,*) "Omega: ", omega
            !write(*,*) "k = ", k
            !write(*,*) "Radius(k): " , s% r(k)
            write(*,*) "Ratio, r(k) / r(1), Heat mass / mass : " , s% r(k) / radius, mcumul / mass
            write(*,*) "Extra heat: " , extra_heat_var

      end subroutine accretion_other_energy
      end module run_star_extras
      
