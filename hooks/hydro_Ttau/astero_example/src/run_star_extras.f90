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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none

      include 'hydro_Ttau/hydro_Ttau_def.inc'
      
      contains

      include 'hydro_Ttau/624.dek'
      include 'hydro_Ttau/hydro_Ttau_proc.inc'
      
      subroutine extras_controls(id, ierr)
         use astero_def, only: star_astero_procs
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         call hydro_Ttau_setup(id, ierr)
         if (ierr /= 0) stop 'failed to setup hydro_Ttau in extras_controls'

         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% other_gradr_factor => hydro_Ttau_gradr_factor
         s% other_surface_PT => hydro_Ttau_surface_PT
         
         s% job% warn_run_star_extras =.false.       
            
         include 'set_star_astero_procs.inc'
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call hydro_Ttau_update(id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

      end subroutine extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         use astero_def, only: my_var1, my_var2, my_var3
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         
         include 'formats'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call hydro_Ttau_update(id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

         extras_check_model = keep_going
         
         my_var1 = s% delta_Pg
         !write(*,2) 'delta_Pg', s% model_number, my_var1


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depenending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination conditon'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      subroutine set_my_vars(id, ierr) ! called from star_astero code
         !use astero_search_data, only: include_my_var1_in_chi2, my_var1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ! my_var's are predefined in the simplex_search_data.
         ! this routine's job is to assign those variables to current value in the model.
         ! it is called whenever a new value of chi2 is calculated.
         ! only necessary to set the my_var's you are actually using.
         ierr = 0
         !if (include_my_var1_in_chi2) then
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            !my_var1 = s% Teff
         !end if
      end subroutine set_my_vars


      subroutine will_set_my_param(id, i, new_value, ierr) ! called from star_astero code
         !use astero_search_data, only: vary_my_param1
         integer, intent(in) :: id
         integer, intent(in) :: i ! which of my_param's will be set
         real(dp), intent(in) :: new_value
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0

         if (i == 1) then
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            s% x_ctrl(1) = new_value
         end if

      end subroutine will_set_my_param


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         how_many_extra_history_columns = 2 ! 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         
         
         names(1) = 'alpha'
         vals(1) = s% mixing_length_alpha

         names(2) = 'tau_base'
         vals(2) = s% tau_base

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         
         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 2
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         names(1) = 'alpha'
         vals(1) = s% mixing_length_alpha

         names(2) = 'tau_base'
         vals(2) = s% tau_base

      end subroutine data_for_extra_profile_header_items
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         use astero_def
         use utils_lib, only: mv
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character (len=256) :: format_string, num_string, basename
         ierr = 0

         write(format_string,'( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_string,format_string) sample_number+1 ! sample number hasn't been incremented yet
         basename = trim(sample_results_prefix) // trim(num_string)
         call mv(best_model_fgong_filename, trim(basename) // trim('.fgong'), skip_errors=.true.)
         call mv(best_model_profile_filename, trim(basename) // trim('.profile'), skip_errors=.true.)

      end subroutine extras_after_evolve

      end module run_star_extras
      
