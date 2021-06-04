      subroutine low_mass_torques(id, ierr)
         use math_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, k, k_end, wind_distribution_option, Jdot_option
         real(dp) :: jdot, net_I, solar_photosphere_p, solar_tcz
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         Jdot_option = s% x_integer_ctrl(2)


         ! If there's a disk, lock to it
         if (s% star_age < s% x_ctrl(1)) then
            do j=1,s% nz
               s% omega(j) = s% x_ctrl(2)
            end do
         ! Otherwise, find the net torque
         else if (Jdot_option == 1) then
            jdot = s% x_ctrl(3) * 1.1265d47 * s%omega(1) * s%omega(1) * s%omega(1) ! Constant from Kawaler 1988
         else if (Jdot_option == 2) then
            jdot = s% x_ctrl(3) * 1.1265d47 * s%omega(1) * pow2(min(s% omega(1), s% x_ctrl(4)))
         else if (Jdot_option == 3) then
            solar_photosphere_p = 1.12612d5
            solar_tcz = 1.09d6
            jdot = s% x_ctrl(3) * 5.557d46 * pow(s%r(1) / Rsun, 3d0) * pow(s%m(1) / Msun, -0.22d0) ! True solar torque of 1.3e30 divided by omega_sun^3
            jdot = jdot * pow(s%L(1) / Lsun, 0.56d0)
            jdot = jdot * pow(photosphere_pressure(id) / solar_photosphere_p, 0.44d0)
            jdot = jdot * s%omega(1) * pow2(min(s%omega(1) * s% x_ctrl(5) / solar_tcz, s% x_ctrl(4)))
         end if

         ! Distribute the torque
         if (s% star_age > s% x_ctrl(1)) then
 
            wind_distribution_option = s% x_integer_ctrl(1)

            if (wind_distribution_option == 1) then
               ! Distribute jdot proportional to i_rot in the top convection zone

               k_end = s% x_integer_ctrl(5)

               ! Catches if the CZ disappeared
               if (k_end < 10) then
                  k_end = 10
               else if (k_end > s% nz) then
                  k_end = s% nz
               end if

               ! Figure out how massive the top convection zone is
               net_I = 0d0
               do k=1,k_end   
                  net_I = net_I + s% i_rot(k) * s% dm(k) 
               end do

               ! Distribute Jdot
               do k=1,k_end
                  s% extra_jdot(k) = -jdot * s% i_rot(k) * s% dm(k)  / net_I
               end do
            else if (wind_distribution_option == 2) then
               ! Distribute jdot proportional to i_rot in the whole star
               net_I = 0d0
               do k=1,s% nz
                  net_I = net_I + s% i_rot(k) * s% dm(k) 
               end do

               do k=1,s% nz   
                  s% extra_jdot(k) = -jdot * s% i_rot(k) * s% dm(k) / net_I
               end do
            end if
   
            do k=1,s% nz
                  s% extra_jdot(k) = s% extra_jdot(k) / s% dm(k)
            end do

         end if

   

         ! Shoulduse one or the other but not both
         s% extra_omegadot(:) = 0
      end subroutine low_mass_torques

      ! None of the following functions are called unless you set their
      ! function point in extras_control.

      real(dp) function convective_tau(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         real(dp) :: conv_tau
         integer j, j_cz
         integer ierr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         j_cz = 1
         do j=2,s%nz
            if (s%brunt_N2(j-1) < 1d-8 .and. s% brunt_N2(j) > 1d-8) then
               ! Found a lower CZ boundary
               j_cz = j
               exit ! Leave the loop
            end if
            if (j == s% nz) then
               j_cz = s%nz
            end if
         end do

         s%x_integer_ctrl(5) = j_cz

         do j=1,s%nz
            if (s%r(j)-s%scale_height(j) < s%r(j_cz)) then
               conv_tau = s%scale_height(j) / s%conv_vel(j)
               exit ! Leave the loop
            end if
         end do

         convective_tau = conv_tau
      end function convective_tau      
      
      real(dp) function photosphere_pressure(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         real(dp) :: p_tau, tau23
         integer j
         integer :: ierr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         tau23 = 2d0 / 3d0
         do j=1,s%nz
            if (s% tau(j) > tau23) then
               if (j == 1) then
                  ! Extrapolate
                  p_tau = s% peos(1) * tau23 / s%tau(1)
               else
                  ! Interpolate
                  p_tau = (s% peos(j) * (s% tau(j) - tau23) + s% peos(j-1) *(tau23 - s% tau(j-1))) &
                           / (s%tau(j) -  s% tau(j-1))
               end if
               exit ! Leave the loop
            end if
         end do

         photosphere_pressure = p_tau
      end function photosphere_pressure   
