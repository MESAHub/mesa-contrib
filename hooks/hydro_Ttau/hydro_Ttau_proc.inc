      subroutine hydro_Ttau_setup(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp) :: Teff, grv, alpha
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call init_hydro_Ttau_data(ierr)
         if (ierr /= 0) stop 'failed to initialize hydro_Ttau data in hydro_Ttau_setup'
         ! write(*,*) 'initialized hydro_Ttau in hydro_Ttau_setup'

         Teff = 5772d0
         grv = 27437d0
         alpha = -1d0
         
         call get_qhpf(Teff, grv, alpha, ierr)
         ! write(*,*) 'called get_qhpf in extras_controls'
         s% mixing_length_alpha = alpha*s% x_ctrl(1)
         s% tau_base = tau_base_hydro_Ttau

      end subroutine hydro_Ttau_setup


      subroutine hydro_Ttau_update(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp) :: grv, Teff, alpha

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         Teff = s% Teff
         grv = exp10(s% log_surface_gravity) ! exp10(4.438d0) 
         call get_qhpf(s% Teff, grv, alpha, ierr)
         s% mixing_length_alpha = alpha*s% x_ctrl(1)
         s% tau_base = tau_base_hydro_Ttau
      end subroutine hydro_Ttau_update


      subroutine hydro_Ttau_gradr_factor(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         do k = 1, s% nz
            s% gradr_factor(k) = 1d0 + dq_dtau(s% tau(k))
         end do         
      end subroutine hydro_Ttau_gradr_factor


      real(dp) function q(tau)
         real(dp), intent(in) :: tau
         real(dp) :: dq_dtau
         integer :: ierr

         ierr = 0
         call eval1_qhpf_tau(tau, q, dq_dtau, ierr)
         if (ierr /= 0) stop 'eval1_qhpf_tau failed in run_star_extras:q'
         ! write(*,*) q, dq_dtau
      end function q


      real(dp) function dq_dtau(tau)
         real(dp), intent(in) :: tau
         real(dp) :: q
         integer :: ierr

         ierr = 0
         call eval1_qhpf_tau(tau, q, dq_dtau, ierr)
         if (ierr /= 0) stop 'eval1_qhpf_tau failed in run_star_extras:dq_dtau'
         ! write(*,*) q, dq_dtau
       end function dq_dtau


       subroutine hydro_Ttau_surface_PT( &
            id, skip_partials, Teff, &
            lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

         use atm_lib, only : atm_Teff

         integer, intent(in) :: id ! star id if available; 0 otherwise
         logical, intent(in) :: skip_partials

         real(dp), intent(out) :: Teff ! temperature at photosphere
         real(dp), intent(out) :: lnT ! natural log of temperature at base of atmosphere
         real(dp), intent(out) :: lnP ! natural log of pressure at base of atmosphere (Pgas + Prad)
         
         ! partial derivatives of lnT and lnP
         real(dp), intent(out) :: dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap
         real(dp), intent(out) :: dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap

         integer, intent(out) :: ierr  ! == 0 means AOK

         real(dp) :: L, R, M, kap
         real(dp) :: phot_r, phot_m, phot_v, phot_L, phot_cs, &
              phot_kap, phot_ysum
         integer :: k_phot

         type (star_info), pointer :: s
         real(dp) :: tau, Pextra_factor, g
         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L = s% L(1)
         R = s% R(1)
         M = s% M(1)

         kap = s% opacity(1)

         call get_phot_info(s, phot_r, phot_m, phot_v, &
              phot_L, phot_cs, phot_kap, phot_ysum, k_phot)
         
         tau = s% tau_base*s% tau_factor
         Teff = atm_Teff(L, phot_r)
         Pextra_factor = 1d0 ! 0d0
         g = s% cgrav(1)*M/(R*R)

         call simple_photosphere(tau, &
              g, Pextra_factor, kap, L, s% cgrav(1), M, Teff, skip_partials, &
              lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, & 
              lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
              ierr)
         
      end subroutine hydro_Ttau_surface_PT


      subroutine simple_photosphere( &
            tau, g, Pextra_factor, kap, L, cgrav, M, Teff, skip_partials, &
            lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, & 
            lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
         real(dp), intent(in) :: tau, g, Pextra_factor, kap, L, cgrav, M, Teff
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap
         real(dp), intent(out) :: lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

         real(dp) :: &
            P0, Pextra, Pextra2, Pfactor, P, d_logg_dlnR, d_logg_dlnM, &
            dP0_dlnR, dP0_dlnkap, dP0_dL, dP0_dlnM, &
            dPfactor_dlnR, dPfactor_dlnkap, dPfactor_dL, dPfactor_dlnM, &
            dlnTeff_dL, dlnTeff_dlnR, Teff4, T4, T, dlnT_dlnTeff
         
         include 'formats'
         
         ierr = 0
         
         ! P = (tau*g/kap)*(1 + 1.6d-4*kap*(L/Lsun)/(M/Msun))
         ! the factor comes from including nonzero Prad at tau=0
         ! see, e.g., Cox & Giuli, Section 20.1
         P0 = tau*g/kap   
         if (Pextra_factor < 0) then ! old form
            Pextra = 1.6d-4*kap*(L/Lsun)/(M/Msun)
         else      
            ! Pextra = Pextra_factor*(kap/tau)*(L/M)/(6d0*pi*clight*cgrav)
            Pextra = Pextra_factor*(kap/tau)*(L/M)/(4d0*pi*clight*cgrav)*q(tau)
         end if
         ! thanks to Konstantin Pavlovskii for correcting an earlier version of this
         ! for tau=2/3 the new Pextra is smaller by about a factor of 4
         ! since massive stars can have L/M > 10^4,
         ! this can make a big difference in P surface.

         Pfactor = 1 + Pextra
         P = P0*Pfactor
         lnP = log(P)

         Teff4 = Teff*Teff*Teff*Teff
         T4 = 0.75d0*Teff4*(tau + q(tau)) ! q(tau) was 2/3, eddington
         lnT = log(T4)*0.25d0
         T = exp(lnT)

         if (.not. skip_partials) then

            d_logg_dlnR = -2
            d_logg_dlnM = 1

            dP0_dlnR = d_logg_dlnR*P0
            dP0_dlnkap = -P0
            dP0_dL = 0
            dP0_dlnM = d_logg_dlnM*P0
         
            dPfactor_dlnR = 0
            dPfactor_dlnkap = Pextra
            dPfactor_dL = Pextra/L
            dPfactor_dlnM = -Pextra

            dlnP_dL = (dP0_dL*Pfactor + P0*dPfactor_dL)/P
            dlnP_dlnR  = (dP0_dlnR*Pfactor + P0*dPfactor_dlnR)/P
            dlnP_dlnM = (dP0_dlnM*Pfactor + P0*dPfactor_dlnM)/P
            dlnP_dlnkap = (dP0_dlnkap*Pfactor + P0*dPfactor_dlnkap)/P
         
            dlnTeff_dL = 1/(4*L)
            dlnTeff_dlnR = -0.5d0
         
            dlnT_dlnTeff = 1
         
            dlnT_dL = dlnTeff_dL*dlnT_dlnTeff
            dlnT_dlnR  = dlnTeff_dlnR*dlnT_dlnTeff
            dlnT_dlnM = 0
            dlnT_dlnkap = 0

       else

          dlnP_dL = 0._dp
          dlnP_dlnR  = 0._dp
          dlnP_dlnM = 0._dp
          dlnP_dlnkap = 0._dp

          dlnT_dL = 0._dp
          dlnT_dlnR = 0._dp
          dlnT_dlnM = 0._dp
          dlnT_dlnkap = 0._dp

       endif
         
      end subroutine simple_photosphere


      subroutine init_hydro_Ttau_data(ierr)

        implicit none

        integer :: Mqtau, Msims, iounit, i, j
        integer, intent(inout) :: ierr
        character(len=256) :: mesa_contrib_dir, hydro_Ttau_data

        ist = 1 ! might be faster if this can be re-used
        logtau_base_hydro_Ttau = -0.3d0 ! this is just a reasonable initial guess
        tau_base_hydro_Ttau = exp10(logtau_base_hydro_Ttau)
        lntau_base_hydro_Ttau = log(tau_base_hydro_Ttau)

        call get_environment_variable("MESA_CONTRIB_DIR", mesa_contrib_dir)
        hydro_Ttau_data = trim(mesa_contrib_dir)//'/hooks/hydro_Ttau/TtauFeH0.dat'

        write(*,*) 'Loading '//trim(hydro_Ttau_data)//'...'

        open(newunit=iounit, file=trim(hydro_Ttau_data), &
             action='read', status='old', iostat=ierr)

        read(iounit, '(12x,2I12)') Mqtau, Msims
        read(iounit, '(12x,'//sNsims//'f12.3)') Teffsl
        read(iounit, '(12x,'//sNsims//'f12.3)') gravsl
        read(iounit, '(12x,'//sNsims//'f12.3)') FeH
        read(iounit, '(12x,'//sNsims//'f12.3)') alph
        read(iounit, '(12x,'//sNsims//'f12.3)') sigalph
        read(iounit, '('//sNsims1//'f12.8)') aa1

        close(iounit)

        ! write(*,*) Teffsl
        ! write(*,*) gravsl
        ! write(*,*) FeH

        ! if ((Mqtau /= Nqtau) .or. (Msims /= Nsims)) ierr = -1

        ! Split off the log10(tau) array from aa1 and store in qltau
        do j=1, Nqtau
           qltau(j) = aa1(1,j)
           do i=1, Nsims
              aa(i,j) = aa1(i+1,j)
           end do
        end do

        ! data in table is actually Teff but we need log10(Teff)
        do i=1, Nsims
           Teffsl(i) = log10(Teffsl(i))
        end do

        ! Use LTeff3 as basis for triangulation since Teff varies a
        ! lot less than g - this ensures less acute angles.
        do i=1, Nsims
           LTeff3(i) = (Teffsl(i)-3d0)*8d0
           Logg3(i) = gravsl(i)
        end do

        ! Remove these 11 lines if used with tables other than
        ! TtauFeH0.dat from Trampedach et al. (2014).
        LTeff3(19) = LTeff3(19) - 0.15d0
        Logg3 (29) = Logg3 (29) + 0.08d0
        Logg3 ( 9) = Logg3 ( 9) + 0.3d0
        Logg3 ( 7) = Logg3 ( 7) - 0.2d0
        Logg3 (11) = Logg3 (11) + 0.2d0
        Logg3 (33) = Logg3 (33) - 0.07d0
        Logg3 (28) = Logg3 (28) - 0.4d0
        LTeff3(28) = LTeff3(28) - 0.02d0
        LTeff3(32) = LTeff3(32) - 0.1d0
        LTeff3(37) = LTeff3(37) + 0.03d0
        Logg3 (36) = Logg3 (36) - 0.15d0

        ! Sort Teffsl, gravsl and alpha [=aa(*,1)] w.r.t. Teffsl and store
        ! indexes in isrt
        ! do i =1, Nsims
        !    write(*,*) gravsl(i), Teffsl(i), alph(i)
        ! end do
        
        call reordr(Nsims,3,Teffsl,gravsl,alph,isrt)
        call permut(Nsims,isrt,LTeff3)
        call permut(Nsims,isrt,Logg3)

        ! Set-up triangulation in the reordered log10(Teff) and log10(grav)
        call trmesh(Nsims, LTeff3, Logg3, iadj, iend, ierr)
        if (ierr /= 0) then
           stop 'failed to create triangulation mesh in init_hydro_Ttau_data'
        end if

        do j=1, Nqtau
           call permut(Nsims,isrt,aa(1,j))
        enddo

        allocate(f1(4*Nqtau))
        allocate(work1(nwork*Nqtau))

        f1 = 0d0
        work1 = 0d0

      end subroutine init_hydro_Ttau_data

      subroutine get_qhpf(Teff, grv, alpha, ierr)
        use num_lib, only : safe_root_with_guess
        use interp_1d_lib, only : interp_pm, interp_value_and_slope
        real(dp), intent(in) :: Teff, grv
        real(dp), intent(out) :: alpha
        real(dp) :: q_base, ignore
        integer :: ierr, j !, ist

        real(dp) :: Teffl, gravl

        ! working arrays for safe root
        integer, pointer :: ipar(:)
        integer, target :: ipara(1)
        real(dp), pointer :: rpar(:)
        real(dp), target :: rpara(1)
        integer :: newt_imax, imax, lipar, lrpar
        real(dp) :: x_guess, dx, x1, x3, y1, y3, epsx, epsy

        Teffl = log10(Teff)
        gravl = log10(grv)

        ! ist = 1
        ! Perform linear interpolation of alpha, in log10(Teff) and
        ! log10(grav)
        call intrc0(Nsims, Teffl, gravl, Teffsl, gravsl, alph, &
             iadj, iend, ist, alpha, ierr)
        if (ierr < 0) then
           write(*,*) 'STOP after interpolation of alpha: ierr =', ierr
           stop 'ierr /= 0 after interpolation of alpha'
        end if
        ! if (ierr > 0) write(*,*) 'WARNING extrapolated to get alpha'

        ! Perform linear interpolation of Hopf function, in
        ! log10(Teff) and log10(grav)
        do j=1, Nqtau
           call intrc0(Nsims, Teffl, gravl, Teffsl, gravsl, &
                aa(1,j), iadj, iend, ist, qhpf(j), ierr)
        end do
        if (ierr < 0) then
           write(*,*) 'STOP after interpolation of qhpf: ierr =', ierr
           stop 'ierr /= 0 after interpolation of qhpf'
        end if
        ! if (ierr > 0) write(*,*) 'WARNING extrapolated to get qhpf'

        ! now create an interpolant
        f1 = 0d0
        do j=1, Nqtau
           f1(4*(j-1)+1) = qhpf(j)
        end do
        call interp_pm(qltau, Nqtau, f1, nwork, work1, interp_dbg_str, ierr)
        if (ierr /= 0) stop 'failed to create interpolant for qhpf'

        ! get tau_base
        nullify(rpar, ipar)
        rpara = 0d0
        ipara = 0
        rpar => rpara
        ipar => ipara
        lrpar = 1
        lipar = 1
        newt_imax = 100 ! was 10
        imax = 100
        dx = 0.1d0
        x_guess = logtau_base_hydro_Ttau
        ! bounds possibly generous to allow more extrapolation?
        x1 = -0.6d0 ! was -0.5d0
        x3 =  0.0d0 ! was -0.2d0
        call interp_value_and_slope(qltau, Nqtau, f1, x1, y1, ignore, ierr)
        call interp_value_and_slope(qltau, Nqtau, f1, x3, y3, ignore, ierr)
        ! write(*,*) y1, y3
        ! y1 = y1 + exp10(x1) - 4d0/3d0
        ! y3 = y3 + exp10(x3) - 4d0/3d0
        ! write(*,*) y1, y3
        epsx = 1d-9
        epsy = 1d-9

        ! call signature is
     !  real(dp) function safe_root_with_guess(
     ! >      f, x_guess, dx, x1, x3, y1, y3, newt_imax, imax, 
     ! >      epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
        logtau_base_hydro_Ttau = safe_root_with_guess( &
             f, x_guess, dx, x1, x3, y1, y3, newt_imax, imax, epsx, epsy, &
             lrpar, rpar, lipar, ipar, ierr)
        tau_base_hydro_Ttau = exp10(logtau_base_hydro_Ttau)
        lntau_base_hydro_Ttau = log(tau_base_hydro_Ttau)
        ! write(*,*) logtau_base_hydro_Ttau, tau_base_hydro_Ttau, lntau_base_hydro_Ttau

        ! y3 = f(logtau_base_hydro_Ttau, ignore, lrpar, rpar, lipar, ipar, ierr)
        ! write(*,*) y3
        ! call interp_value_and_slope(qltau, Nqtau, f1, logtau_base_hydro_Ttau, &
        !      y3, ignore, ierr)
        ! write(*,*) y3 + tau_base_hydro_Ttau - 4d0/3d0

        contains

        real(dp) function f(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
          real(dp), intent(in) :: x
          real(dp), intent(out) :: dfdx
          integer, intent(in) :: lrpar, lipar
          integer, intent(inout), pointer :: ipar(:) ! (lipar)
          real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
          integer, intent(out) :: ierr
          real(dp) :: qhpf_out, dqhpf_out
          call interp_value_and_slope(qltau, Nqtau, f1, x, qhpf_out, dqhpf_out, ierr)

          dfdx = dqhpf_out + log(10d0)*exp10(x)
          f = qhpf_out + exp10(x) - 4d0/3d0
        end function f

      end subroutine get_qhpf
      
      subroutine eval1_qhpf(logtau, qhpf_out, dqhpf_out, ierr)
        use interp_1d_lib, only : interp_value_and_slope

        real(dp), intent(in) :: logtau
        real(dp), intent(out) :: qhpf_out, dqhpf_out
        integer :: ierr

        if (logtau > qltau(Nqtau)) then
           qhpf_out = qhpf(Nqtau)
           dqhpf_out = 0d0
           return
        end if

        call interp_value_and_slope(qltau, Nqtau, f1, logtau, &
             qhpf_out, dqhpf_out, ierr)

      end subroutine eval1_qhpf

      ! WARNING: not sure what will happen if logtau < -4.5 or > 2
      subroutine eval_qhpf(nlogtau, logtau, qhpf_out, dqhpf_out, ierr)
        use interp_1d_lib, only : interp_values_and_slopes
        integer, intent(in) :: nlogtau
        real(dp), intent(in), dimension(nlogtau) :: logtau
        real(dp), intent(out), dimension(nlogtau) :: qhpf_out, dqhpf_out
        integer :: ierr

        call interp_values_and_slopes(qltau, Nqtau, f1, nlogtau, logtau, &
             qhpf_out, dqhpf_out, ierr)

      end subroutine eval_qhpf

      ! you probably want to interpolate in tau, and get derivative w.r.t tau
      subroutine eval1_qhpf_tau(tau, qhpf_out, dqhpf_out, ierr)
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: qhpf_out, dqhpf_out
        real(dp) :: logtau
        integer :: ierr

        logtau = log10(tau)
        call eval1_qhpf(logtau, qhpf_out, dqhpf_out, ierr)
        dqhpf_out = dqhpf_out/log(10d0)/tau
      end subroutine eval1_qhpf_tau

      subroutine eval_qhpf_tau(ntau, tau, qhpf_out, dqhpf_out, ierr)
        integer, intent(in) :: ntau
        real(dp), intent(in), dimension(ntau) :: tau
        real(dp), intent(out), dimension(ntau) :: qhpf_out, dqhpf_out
        real(dp), dimension(ntau) :: logtau
        integer :: ierr, i

        do i=1,ntau
           logtau(i) = log10(tau(i))
        end do
        call eval_qhpf(ntau, logtau, qhpf_out, dqhpf_out, ierr)
        dqhpf_out = dqhpf_out/log(10d0)/tau
      end subroutine eval_qhpf_tau

      
      ! cloned from $MESA_DIR/star/private.star_utils.f90:868-930
      subroutine get_phot_info(s, r, m, v, L, cs, kap, ysum, k_phot)
         type (star_info), pointer :: s
         real(dp), intent(out) :: r, m, v, L, cs, kap, ysum
         integer, intent(out) :: k_phot

         integer :: k
         real(dp) :: tau00, taup1, dtau, r003, rp13, r3, tau_phot

         include 'formats'

         tau00 = 0
         taup1 = 0
         ysum = 0
         r = s% r(1)
         m = s% m(1)

         L = s% L(1)
         if (s% u_flag) then
            v = s% u(1)
         else if (s% v_flag) then
            v = s% v(1)
         else
            v = 0d0
         end if
         cs = s% csound(1)
         kap = s% opacity(1)
         
         k_phot = 1
         tau_phot = s% tau_base
         tau00 = s% tau_factor*tau_phot
         if (tau00 >= tau_phot) return
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            ysum = ysum + s% rho(k)*(s% r(k) - s% r(k+1))
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               r003 = s% r(k)*s% r(k)*s% r(k)
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
               r3 = r003 + (rp13 - r003)*(tau_phot - tau00)/dtau
               r = pow(r3,1d0/3d0)
               m = s% m(k) - s% dm(k)*(tau_phot - tau00)/dtau
               if (s% u_flag) then
                  v = s% u_face(k) + &
                     (s% u_face(k+1) - s% u_face(k))*(tau_phot - tau00)/dtau
               else if (s% v_flag) then
                  v = s% v(k) + (s% v(k+1) - s% v(k))*(tau_phot - tau00)/dtau
               end if
               L = s% L(k) + (s% L(k+1) - s% L(k))*(tau_phot - tau00)/dtau
               k_phot = k
               cs = s% csound(k_phot)
               kap = s% opacity(k_phot)
               return
            end if
            tau00 = taup1
         end do
         k_phot = s% nz
         r = s% R_center
         m = s% m_center
         v = s% v_center
         L = s% L_center
         cs = s% csound(k_phot)
         kap = s% opacity(k_phot)
      end subroutine get_phot_info
