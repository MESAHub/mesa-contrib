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
      use chem_def

      implicit none

      contains

      include 'dbconvpen/dbconvpen.inc'

      subroutine extras_controls(id, ierr)
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

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

          ! s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
          s% lxtra(1) = .true. ! turns it on!
          s% other_overshooting_scheme => extended_convective_penetration
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
          logical :: do_retry
          integer k
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_check_model = keep_going

          ! Flag PZ as anonymous_mixing
          if (s% lxtra(1)) then
            do k=s% ixtra(2), s% ixtra(1), -1
                s%mixing_type(k) = anonymous_mixing
            end do
          endif

          do_retry = .false.

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
          how_many_extra_history_columns = 7
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
          integer, intent(in) :: id, n
          character (len=maxlen_history_column_name) :: names(n)
          real(dp) :: vals(n)
          real(dp) :: r_cb
          integer :: k
          integer, intent(out) :: ierr
          type (star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          call star_eval_conv_bdy_r(s, 1, r_cb, ierr)

          ! s% xtra(2)=0
          ! s% xtra(3)=0
          ! s% xtra(4)=0

          ! if (s% mixing_type(s% nz) /= convective_mixing) then
          !    s% xtra(5) = 0
          !    s% xtra(1) = 0
          !    s% xtra(6) = 0
          ! else
          !    call star_eval_conv_bdy_k(s, 1, k, ierr)
          !    s% xtra(5) = r_cb
          !    s% xtra(1) = s%m(k)
          !    s% xtra(6) = s%rho(k)
          ! endif

          names(1) = 'm_core'
          names(2) = 'mass_pen_zone'
          names(3) = 'delta_r_pen_zone'
          names(4) = 'alpha_pen_zone'
          names(5) = 's% xtra(5)'
          names(6) = 's% xtra(6)_pen'
          names(7) = 'r_cb'

          vals(1) = s% xtra(1)/Msun
          vals(2) = s% xtra(2)/Msun
          vals(3) = s% xtra(3)/Rsun
          vals(4) = s% xtra(4)
          vals(5) = s% xtra(5)/Rsun
          vals(6) = s% xtra(6)
          vals(7) = r_cb/Rsun

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
          ! integer :: vals_nr
          type (star_info), pointer :: s
          integer :: k
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

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
          character (len=200) :: fname
          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_finish_step = keep_going

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
