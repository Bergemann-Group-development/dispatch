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

! ***********************************************************************
!     x_ctrl's:
!     x_ctrl(1): alpha
!     x_ctrl(19): initial mass
!     x_ctrl(5): Lacc_out
!     x_ctrl(6): Lacc_in
!     x_ctrl(50): mdot from MESA_driver
!
! ***********************************************************************


!! ---------------------------------------------------------------------------
!! ---------------------------------------------------------------------------

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none

      ! these routines are called by the standard run_star check_model
      contains

      !include 'standard_run_star_extras.inc'

        subroutine extras_controls(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           s% extras_startup => extras_startup
           !s% extras_check_model => extras_check_model
           !s% extras_finish_step => extras_finish_step
           !s% extras_after_evolve => extras_after_evolve
           s% how_many_extra_history_columns => how_many_extra_history_columns
           s% data_for_extra_history_columns => data_for_extra_history_columns
           !s% how_many_extra_profile_columns => how_many_extra_profile_columns
           !s% data_for_extra_profile_columns => data_for_extra_profile_columns
           !s% job% warn_run_star_extras =.false.


        	 !New routines to add
        	 s% other_energy => accretion_other_energy
           s% other_adjust_mdot => other_ramses_adjust_mdot
           s% job% warn_run_star_extras =.false.

        end subroutine extras_controls

        integer function extras_startup(id, restart, ierr)
           integer, intent(in) :: id
           logical, intent(in) :: restart
           integer, intent(out) :: ierr

           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_startup = 0

           !!! ------------------------------
           s% x_ctrl(19) = 0.0014!s% mstar
           !!! ------------------------------

        end function extras_startup

        subroutine other_ramses_adjust_mdot(id, ierr)
            use star_def
            implicit none
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type(star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr) ! retrieve star id

            s% mstar_dot = s% x_ctrl(50)
            !write(*,*) "Adjust_mdot: ", s% mstar_dot

        end subroutine other_ramses_adjust_mdot

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
            alpha = s% x_ctrl(1)

            ! ----------------------------------------------------------
            ! Time varying alpha, alpha = 0.25 if mdot > 1e-5
            ! crossover for alpha in range mdot [1e-5, 1e-6]
            ! settle at alpha = 0.01 when mdot < 1e-6
          !   mdot_msun = mdot*(secyer/msun)
          !  !write(*,*) "Mdot (msun)", mdot_msun
          !   if (mdot_msun >= 1e-5) then
          !     alpha = 0.25
          !   else if (mdot_msun <= 1e-6) then
          !     alpha = 0.01
          !   else
          !     alpha = exp(B+A*tanh(zp-log(mdot_msun)))
          !   endif

            mcumul = 0.
            k = 0
            Lacc = 0.5 * (s% cgrav(1) * (mdot) * (mass)) / radius
            extra_heat_var = Lacc * alpha
            s% x_ctrl(6) = extra_heat_var
            s% x_ctrl(5) = Lacc * (1-alpha)
            do while ((mcumul < mdot*s% dt .or. extra_heat_var / max(mcumul,1.) > 1e4) .and. mcumul < 0.8 * mass)
              k = k + 1
              mcumul = mcumul + s% dm(k)
            enddo

            extra_heat_var = Lacc * alpha / mcumul
            s% extra_heat(1:k) = extra_heat_var
            !write(*,*) "k = ", k
            !write(*,*) "Radius(k): " , s% r(k)
            !write(*,*) "Ratio, r(k) / r(1), Heat mass / mass : " , s% r(k) / radius, mcumul / mass
            !write(*,*) "Extra heat: " , extra_heat_var

        end subroutine accretion_other_energy

        integer function how_many_extra_history_columns(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        how_many_extra_history_columns = 3
        end function how_many_extra_history_columns


        subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
        type (star_info), pointer :: s
        integer, intent(in) :: id, id_extra, n
        character (len=maxlen_history_column_name) :: names(n)
        real(dp) :: vals(n)
        integer, intent(out) :: ierr
        real(dp) :: mdot
        ierr  = 0

        call star_ptr(id, s, ierr)
        if (ierr /= 0) then ! OOPS
            return
        end if

        mdot = s% mstar_dot

        if (mdot == 0.) then
          s% x_ctrl(5) = 0.0
          s% x_ctrl(6) = 0.0
        endif

        names(1) = 'extra_lum'
        vals(1) = s% x_ctrl(5)

        names(2) = 'extra_heat'
        vals(2) = s% x_ctrl(6)

        names(3) = 'mdot'
        vals(3) = mdot

        ierr = 0

        !note: do NOT add the extras names to history_columns.list
        ! the history_columns.list is only for the built-in log column options.
        ! it must not include the new column names you are adding here.


        end subroutine data_for_extra_history_columns



      end module run_star_extras
