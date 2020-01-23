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
      use crlibm_lib
      
      implicit none


      ! these routines are called by the standard run_star check_model
      contains
      
      include 'enhanced_Tayler_Spruit/enhanced_Tayler_Spruit_procs.inc'
      !include 'standard_run_star_extras.inc'

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
       
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% job% warn_run_star_extras = .false.     

         s% other_am_mixing => TSF

      end subroutine extras_controls

      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 6
      end function how_many_extra_profile_columns
       
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)

         real(dp) :: vals(nz,n)
         real(dp) :: brunts,alpha,nu_tsf,nu_tsf_t,diffm,difft,omegac,omegaa,omegat,lognuomega,bruntsn2,shearsmooth
         real(dp) :: xmagfmu, xmagft, xmagfdif, xmagfnu, &
            xkap, xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig, &
            xeta

         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k,j,op_err,nsmooth,nsmootham
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         !if (n /= 1) stop 'data_for_extra_profile_columns'
         names(1) = 'log_N'
         names(2) = 'q_smooth'
         names(3) = 'nu_TSF'
         names(4) = 'log_nu_TSF'
         names(5) = 'diff_mag'
         names(6) = 'log_diff_mag'

         vals(:,1) = 1d-50
         vals(:,2) = 1d-30
         vals(:,3) = 1d-50
         vals(:,4) = -50
         vals(:,5) = 1d-50
         vals(:,6) = -50

         op_err = 0
         alpha = 1d0
         nsmooth=5
         nsmootham=nsmooth-3     
         do k=nsmooth+1,s% nz-(nsmooth+1)

            vals(k,1) = safe_log10_cr(sqrt_cr(abs(s% brunt_N2(k))))
            vals(k,2) = s% omega_shear(k)/(2d0*nsmooth+1d0)
            nu_tsf=1d-30
            nu_tsf_t=1d-30

            do j=1,nsmooth
                vals(k,2) = vals(k,2) + (1d0/(2d0*nsmooth+1d0))*( s% omega_shear(k-j) + s% omega_shear(k+j) )
            end do
            shearsmooth = vals(k,2)

                    diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
            vals(k,5) = diffm
            vals(k,6) =  safe_log10_cr(vals(k,5))

            omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt_cr(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt_cr(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0.) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))*pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10_cr(vals(k,3))
            end if

         end do

          do k=s% nz-nsmooth,s% nz

            vals(k,1) = safe_log10_cr(sqrt_cr(abs(s% brunt_N2(k))))
            vals(k,2) = vals(k-1,2)
            nu_tsf=1d-30
            nu_tsf_t=1d-30

            diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
            vals(k,5) = diffm
            vals(k,6) = safe_log10_cr(vals(k,5))

            omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt_cr(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt_cr(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0.) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))*pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10_cr(vals(k,3))
            end if

          end do

          do k=nsmooth,1

            vals(k,1) = safe_log10_cr(sqrt_cr(abs(s% brunt_N2(k))))
            vals(k,2) = vals(k+1,2)
            nu_tsf=1d-30
            nu_tsf_t=1d-30

            diffm = diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
            vals(k,5) = diffm
            vals(k,6) =  safe_log10_cr(vals(k,5))

            omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt_cr(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt_cr(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0d0) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))*pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10_cr(vals(k,3))
            end if
          end do


         !Smooth nu_omega
         do k=nsmootham+1,s% nz-(nsmootham+1)
            if (s% mixing_type(k)==1) then
                vals(k,4) = vals(k,4)
            else
                vals(k,4) = vals(k,4)/(2d0*nsmootham+1d0)
            end if          
            do j=1,nsmootham
                !Don't smooth convective diffusivity into non-convective zones
                if (s% mixing_type(k-j)<3.5d0) then
                    vals(k,4) = vals(k,4)
                !Smooth zones if not including a convective zone
                else 
                    vals(k,4) = vals(k,4) + (1d0/(2d0*nsmootham+1d0))*vals(k-j,4) 
                end if
            end do
            do j=1,nsmootham
                !Don't smooth convective diffusivity into non-convective zones
                if (s% mixing_type(k+j)<3.5d0) then
                    vals(k,4) = vals(k,4)
                !Smooth zones if not including a convective zone
                else 
                    vals(k,4) = vals(k,4) + (1d0/(2d0*nsmootham+1d0))*vals(k+j,4) 
                end if
            end do
            vals(k,3) = pow_cr(10d0,vals(k,4))
         end do

          !Values near inner boundary
         do k=s% nz-nsmootham,s% nz
            vals(k,4) = vals(k-1,4)
            vals(k,3) = pow_cr(10d0,vals(k,4))
         end do

          !Values near outer boundary
         do k=nsmootham,1
            vals(k,4) = vals(k-1,4)
            vals(k,3) = pow_cr(10d0,vals(k,3))
         end do


      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 6
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         use const_def, only: pi
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         real(dp) :: ominteg,ninteg,brunts,dr,J1,J2,J3,JHe,nu_pulse
         integer :: j,k,op_err
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         op_err = 0

         names(1) = 'nu_max' !Aseteroseismic nu_max
         names(2) = 'om_g' !Core rotation as sensed by g-modes  
         names(3) = 'AM_1' !Angular momentum in inner 0.1 Msun
         names(4) = 'AM_2' !Angular momentum in inner 0.2 Msun
         names(5) = 'AM_3' !Angular momentum in inner 0.3 Msun
         names(6) = 'AM_He' !Angular momentum in helium core

         vals(1) = 1d-50
         vals(2) = 1d-50
         vals(3) = 1d-50
         vals(4) = 1d-50
         vals(5) = 1d-50
         vals(6) = 1d-50

         vals(1) = s% nu_max
         nu_pulse = s% nu_max !Set pulsation period via scaling relations
 
         if (s% log_surface_radius < -1d0) then
            nu_pulse = 1000d0 !If white dwarf, set pulsation period = 1000 microhertz
         end if
 
         ominteg=1d-99
         ninteg=1d-50
         J1=1d-50
         J2=1d-50
         J3=1d-50
         JHe=1d-50

         do k = 2, s% nz-1

            if (s% brunt_N2(k) < 1d-14) brunts = 1d-14 
            if (s% brunt_N2(k) > 1d-14) brunts = s% brunt_N2(k)

            if (2d0*pi*nu_pulse/1d6<sqrt_cr(brunts)) then
                if (2d0*pi*nu_pulse/1d6<sqrt_cr(2d0)*s% csound(k)/s% r(k)) then
                   dr = (s% r(k+1)-s% r(k-1))/2d0
                   ominteg = ominteg + sqrt_cr(brunts)*(dr/s% r(k))*s% omega(k)
                    ninteg = ninteg + sqrt_cr(brunts)*(dr/s% r(k))
                end if
            end if

            if (s% m(k)/1.99d33 < 0.1d0) then
                J1 = J1 + (2d0/3d0)*pow2(s% r(k))*s% dm(k)*s% omega(k)
            end if 

            if (s% m(k)/1.99d33 < 0.2d0) then
                J2 = J2 + (2d0/3d0)*pow2(s% r(k))*s% dm(k)*s% omega(k)
            end if

            if (s% m(k)/1.99d33 < 0.3d0) then
                J3 = J3 + (2d0/3d0)*pow2(s% r(k))*s% dm(k)*s% omega(k)
            end if  

            if (s% m(k)/1.99d33 < s% he_core_mass) then
                JHe = JHe + (2d0/3d0)*pow2(s% r(k))*s% dm(k)*s% omega(k)
            end if 

         end do

         vals(2) = ominteg/ninteg
         vals(3) = J1
         vals(4) = J2
         vals(5) = J3
         vals(6) = JHe

      end subroutine data_for_extra_history_columns



      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
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


      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
      ! put these defs at the top and delete from the following routines
      !integer, parameter :: extra_info_alloc = 1
      !integer, parameter :: extra_info_get = 2
      !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
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
            case (extra_info_get)
               dbl = s% extra_work(i)
               !dbl = s% xtra1_array(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
           !s% xtra1_array(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info


      end module run_star_extras
      
