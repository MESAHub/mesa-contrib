! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      
      implicit none
      
      include 'wd_builder/wd_builder_def.inc'

      ! data for composition from file
      integer :: num_pts
      real(dp), allocatable :: xq_data(:), xa_data(:,:)

      contains

      include 'wd_builder/wd_builder.inc'


      subroutine get_xa_from_code(s, q, xa)

        use chem_def

        type (star_info), pointer :: s
        real(dp), intent(in) :: q
        real(dp) :: xa(:)

        xa = 0

        if (q < 0.50d0) then
           xa(s% net_iso(ihe4))  = 0.0d0
           xa(s% net_iso(ic12))  = 0.23d0
           xa(s% net_iso(in14))  = 0.0d0
           xa(s% net_iso(io16))  = 0.75d0
           xa(s% net_iso(ine22)) = 0.02d0

        else if (0.50 <= q .and. q < 0.90) then
           xa(s% net_iso(ihe4))  = 0.0d0
           xa(s% net_iso(ic12))  = 0.48d0
           xa(s% net_iso(in14))  = 0.0d0
           xa(s% net_iso(io16))  = 0.50d0
           xa(s% net_iso(ine22)) = 0.02d0

        else if (0.90 <= q .and. q < 0.999) then
           xa(s% net_iso(ihe4))  = 0.0d0
           xa(s% net_iso(ic12))  = 0.65d0
           xa(s% net_iso(in14))  = 0.02d0
           xa(s% net_iso(io16))  = 0.33d0
           xa(s% net_iso(ine22)) = 0.00d0

        else if (0.999 < q) then
           xa(s% net_iso(ihe4))  = 0.98d0
           xa(s% net_iso(ic12))  = 0.0d0
           xa(s% net_iso(in14))  = 0.02d0
           xa(s% net_iso(io16))  = 0.0d0
           xa(s% net_iso(ine22)) = 0.0d0
        end if


      end subroutine get_xa_from_code


      subroutine read_composition_from_file(s, ierr)

        type (star_info), pointer :: s

        integer :: num_species
        integer :: i, iounit, ierr

        open(newunit=iounit, file=trim(s% job% relax_composition_filename), &
             status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
           write(*,*) 'open failed', ierr, iounit
           write(*, '(a)') 'failed to open ' // trim(s% job% relax_composition_filename)
           close(iounit)
           return
        end if

        read(iounit, *, iostat=ierr) num_pts, num_species
        if (ierr /= 0) then
           close(iounit)
           write(*, '(a)') 'failed while trying to read 1st line of ' // &
                trim(s% job% relax_composition_filename)
           return
        end if

        ! if (num_species .ne. s% species) then
        !    write(*,*) 'Error in ',trim(s% job% relax_composition_filename)
        !    write(*,'(a,I4,a)') 'got ',num_species,' species'
        !    write(*,'(a,I4,a)') 'expected ', s% species,' species'
        !    write(*,*)
        !    ierr=-1
        !    return
        ! end if

        allocate(xq_data(num_pts), xa_data(num_species,num_pts))
        do i = 1, num_pts
           read(iounit,*,iostat=ierr) xq_data(i), xa_data(1:num_species,i)
           if (ierr /= 0) then
              close(iounit)
              write(*, '(a)') &
                   'failed while trying to read ' // trim(s% job% relax_composition_filename)
              write(*,*) 'line', i+1
              write(*,*) 'perhaps wrong info in 1st line?'
              write(*,*) '1st line must have num_pts and num_species in that order'
              deallocate(xq_data, xa_data)
              return
           end if
        end do
        close(iounit)

      end subroutine read_composition_from_file



      subroutine get_xa_from_file(s, q, xa)

        use interp_1d_def, only: pm_work_size
        use interp_1d_lib, only: interpolate_vector, interp_pm

        type (star_info), pointer :: s
        real(dp), intent(in) :: q
        real(dp) :: xa(:)

        real(dp), pointer :: work(:)

        real(dp) :: x_new(1), v_new(1)

        integer :: j, ierr

        allocate(work(pm_work_size*num_pts))

        ! file coordinate is xq
        x_new(1) = 1d0 - q

        do j = 1, s% species
           call interpolate_vector( &
                num_pts, xq_data, 1, x_new, xa_data(j,:), v_new, &
                interp_pm, pm_work_size, work, 'get_xa_from_target', ierr)
           xa(j) = v_new(1)
        end do

        deallocate(work)

      end subroutine get_xa_from_file

      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_build_initial_model => wd_builder

         ! get xa from file
         ! this implementation reads the file specified by
         ! the star_job option `relax_composition_filename`
         ! and expects the same file format
         wd_builder_get_xa => get_xa_from_file
         call read_composition_from_file(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to read composition'
            return
         end if

         ! get xq from subroutine
         ! wd_builder_get_xa => get_xa_from_code

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
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

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
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
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
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
         how_many_extra_profile_header_items = 0
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

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

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
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      end module run_star_extras
      
