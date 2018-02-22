! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use crlibm_lib
      
      implicit none
      
      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 120 ! minutes

      integer, parameter :: restart_info_alloc = 1
      integer, parameter :: restart_info_get = 2
      integer, parameter :: restart_info_put = 3
      
      
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% job% warn_run_star_extras=.false.

      end subroutine extras_controls
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: restart_time, prev_time_used
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_restart_info(s)
            s% lxtra1 = .true.
            s% lxtra2 = .true.
         else
            call unpack_restart_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
         end if
         extras_startup = keep_going
      end function extras_startup
      

      integer function extras_check_model(id, id_extra)
         use chem_def
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         integer :: ierr,i
         logical :: is_ne_biggest
         extras_check_model = keep_going

      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         names(1) = 'runtime_minutes'
         vals(1) = dt
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: center_h1, center_h1_old, center_he4, center_he4_old
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         write(*,*) "omega/omega_crit is", s% w_div_w_crit_avg_surf, s% w_div_w_crit_roche(1)

         if(.not. s% lxtra1 .and. abs(s% xtra1_old - s% center_h1) > 0.005) then
             s% dt_next = min(s% dt_next, s% dt * s% min_timestep_factor)
             write(*,*) "reducing dt due to large change in central hydrogen"
         end if
         s% xtra1 = s% center_h1

         !save zams profile one model after re-relax
         if (.not. s% lxtra1 .and. s% lxtra2) then
            call star_write_profile_info(s% id, "LOGS/prof_1ZAMS.data", s% id, ierr)
            if (ierr /= 0) return ! failure in profile
            s% lxtra2 = .false.
         end if

         if (s% lxtra1 .and. &
            abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) < 0.005 .and. &
            s% star_age > 1d2) then
            ! if here, star reached thermal equilibrium (reached ZAMS), so activate mass loss
            s% lxtra1 = .false.
            ! save ZAMS profiles
            call star_relax_uniform_omega(s% id, 2, s% job% new_surface_rotation_v,&
               s% job% num_steps_to_relax_rotation, s% job% relax_omega_max_yrs_dt, ierr)
         end if

         if (s% center_h1 < 1d-3) then
            extras_finish_step = terminate
            write(*,*) "Terminate due to hydrogen depletion"
         end if

         if (extras_finish_step == terminate) then
            call star_write_profile_info(s% id, "LOGS/prof_9FINAL.data", s% id, ierr)
            if (ierr /= 0) return ! failure in profile
         else
            !additional profiles to be saved
            center_h1 = s% xa(s% net_iso(ih1),s% nz)
            center_h1_old = s% xa_old(s% net_iso(ih1),s% nz_old)
            center_he4 = s% xa(s% net_iso(ihe4),s% nz)
            center_he4_old = s% xa_old(s% net_iso(ihe4), s% nz_old)
            if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
               call star_write_profile_info(s% id, "LOGS/prof_2H50.data", s% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 0.3 .and. center_h1_old > 0.3) then
               call star_write_profile_info(s% id, "LOGS/prof_XH30.data", s% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 0.25 .and. center_h1_old > 0.25) then
               call star_write_profile_info(s% id, "LOGS/prof_3H25.data", s% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
               call star_write_profile_info(s% id, "LOGS/prof_4H00.data", s% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6) then
               if (center_he4 < 0.75 .and. center_he4_old > 0.75) then
                  call star_write_profile_info(s% id, "LOGS/prof_5He75.data", s% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
                  call star_write_profile_info(s% id, "LOGS/prof_6He50.data", s% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.25 .and. center_he4_old > 0.25) then
                  call star_write_profile_info(s% id, "LOGS/prof_7He25.data", s% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 1d-6 .and. center_he4_old > 1d-6) then
                  call star_write_profile_info(s% id, "LOGS/prof_8He00.data", s% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               end if
            end if
         end if

         call system_clock(time1,clock_rate)
         call store_restart_info(s)
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         write(*,'(/,a50,2f18.6,99i10/)') 'runtime, retries, backups, steps', &
            dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring data so can do restarts

      
      subroutine alloc_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_alloc)
      end subroutine alloc_restart_info
      
      
      subroutine unpack_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_get)
      end subroutine unpack_restart_info
      
      
      subroutine store_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_put)
      end subroutine store_restart_info
      
      
      subroutine move_restart_info(s,op)
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         call move_int(time0)
         call move_int(time1)
         
         num_ints = i
         
         i = 0
         ! call move_dbl 
         
         num_dbls = i
         
         if (op /= restart_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (restart_info_get)
               dbl = s% extra_work(i)
            case (restart_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            include 'formats'
            i = i+1
            select case (op)
            case (restart_info_get)
               !write(*,3) 'restore int', i, s% extra_iwork(i)
               int = s% extra_iwork(i)
            case (restart_info_put)
               !write(*,3) 'save int', i, int
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (restart_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (restart_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_restart_info
      
      


      end module run_star_extras
      
