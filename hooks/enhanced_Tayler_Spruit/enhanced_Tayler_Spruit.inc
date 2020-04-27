

      subroutine TSF(id, ierr)

        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        integer :: k,j,op_err,nsmooth,nsmootham
        real(dp) :: alpha,shearsmooth,nu_tsf,nu_tsf_t,omegac,omegag,omegaa,&
        omegat,difft,diffm,brunts,bruntsn2,logamnuomega,alphaq

        call star_ptr(id,s,ierr)
        if (ierr /= 0) return

        alpha=1d0
        nsmooth=5
        nsmootham=nsmooth-3
        shearsmooth=1d-30
        op_err = 0

        !Calculate shear at each zone, then calculate TSF torque
        do k=nsmooth+1,s% nz-(nsmooth+1)

            nu_tsf=1d-30
            nu_tsf_t=1d-30
            !Calculate smoothed shear, q= dlnOmega/dlnr
            shearsmooth = s% omega_shear(k)/(2d0*nsmooth+1d0)
            do j=1,nsmooth
                shearsmooth = shearsmooth + (1d0/(2d0*nsmooth+1d0))*( s% omega_shear(k-j) + s% omega_shear(k+j) )
            end do

            diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) !Magnetic diffusivity
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k)) !Thermal diffusivity
            omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(omegaa*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
                    +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
                (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0d0) then
                if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
                end if
                s% am_nu_omega(k) = s% am_nu_omega(k) + max(nu_tsf,nu_tsf_t) + 1d-1
            end if

         end do


          !Values near inner boundary
          do k=s% nz-nsmooth,s% nz
            nu_tsf=1d-30
            nu_tsf_t=1d-30
            shearsmooth = shearsmooth

            diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) !Magnetic diffusivity
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k)) !Thermal diffusivity
            omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
                (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0d0) then
                if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
                end if
                s% am_nu_omega(k) = s% am_nu_omega(k) + max(nu_tsf,nu_tsf_t) + 1d-1
            end if
          end do

      !Values near outer boundary
      do k=nsmooth,1
        nu_tsf=1d-30
        nu_tsf_t=1d-30
        shearsmooth = shearsmooth

        diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) !Magnetic diffusivity
        difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k)) !Thermal diffusivity
        omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
        omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
        brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
            +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
        bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
        brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
        brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
        omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

        ! Calculate nu_TSF
        if (s% brunt_N2(k) > 0d0) then
            if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
            end if
            ! Add TSF enabled by thermal diffusion
            if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
                nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
            end if
            s% am_nu_omega(k) = s% am_nu_omega(k) + max(nu_tsf,nu_tsf_t) + 1d-1
        end if
      end do

      !Smooth nu_omega
      logamnuomega=-3d1
      do k=nsmootham+1,s% nz-(nsmootham+1)
        !Don't smooth convective diffusivity into non-convective zones
        if (s% mixing_type(k)==1) then
            s% am_nu_omega(k) = s% am_nu_omega(k)
        !Smooth zones if not including a convective zone
        else
            logamnuomega = log10(s% am_nu_omega(k))/(2d0*nsmootham+1d0)
        end if 
        do j=1,nsmootham
            !Don't smooth convective diffusivity into non-convective zones
            if (s% mixing_type(k-j)<3.5d0) then
                logamnuomega = log10(s% am_nu_omega(k))
            !Smooth zones if not including a convective zone
            else 
                logamnuomega = logamnuomega + (1d0/(2d0*nsmootham+1d0))*log10(s% am_nu_omega(k-j)) 
            end if
        end do
        do j=1,nsmootham
            !Don't smooth convective diffusivity into non-convective zones
            if (s% mixing_type(k+j)<3.5d0) then
                logamnuomega = logamnuomega
            !Smooth zones if not including a convective zone
            else 
                logamnuomega = logamnuomega + (1d0/(2d0*nsmootham+1d0))*log10(s% am_nu_omega(k+j))
            end if
        end do
        s% am_nu_omega(k) = pow(10d0,logamnuomega)
      end do

      !Values near inner boundary
      do k=s% nz-nsmootham,s% nz
            s% am_nu_omega(k) = s% am_nu_omega(k-1)
      end do

      !Values near outer boundary
      do k=nsmootham,1
            s% am_nu_omega(k) = s% am_nu_omega(k-1)
      end do

      end subroutine TSF

       
      subroutine TS_data_for_extra_history_columns(id, n, names, vals, ierr)
         use const_def, only: pi
         integer, intent(in) :: id, n
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

            if (2d0*pi*nu_pulse/1d6<sqrt(brunts)) then
                if (2d0*pi*nu_pulse/1d6<sqrt(2d0)*s% csound(k)/s% r(k)) then
                   dr = (s% r(k+1)-s% r(k-1))/2d0
                   ominteg = ominteg + sqrt(brunts)*(dr/s% r(k))*s% omega(k)
                    ninteg = ninteg + sqrt(brunts)*(dr/s% r(k))
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

      end subroutine TS_data_for_extra_history_columns


      subroutine TS_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
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

            vals(k,1) = safe_log10(sqrt(abs(s% brunt_N2(k))))
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
            vals(k,6) =  safe_log10(vals(k,5))

            omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0.) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10(vals(k,3))
            end if

         end do

          do k=s% nz-nsmooth,s% nz

            vals(k,1) = safe_log10(sqrt(abs(s% brunt_N2(k))))
            vals(k,2) = vals(k-1,2)
            nu_tsf=1d-30
            nu_tsf_t=1d-30

            diffm =  diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
            vals(k,5) = diffm
            vals(k,6) = safe_log10(vals(k,5))

            omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0.) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10(vals(k,3))
            end if

          end do

          do k=nsmooth,1

            vals(k,1) = safe_log10(sqrt(abs(s% brunt_N2(k))))
            vals(k,2) = vals(k+1,2)
            nu_tsf=1d-30
            nu_tsf_t=1d-30

            diffm = diffmag(s% rho(k),s% T(k),s% abar(k),s% zbar(k),op_err) 
            difft = 16d0*5.67d-5*pow3(s% T(k))/(3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
            vals(k,5) = diffm
            vals(k,6) =  safe_log10(vals(k,5))

            omegaa = s% omega(k)*pow(shearsmooth*s% omega(k)/sqrt(abs(s% brunt_N2(k))),1d0/3d0) !Alfven frequency at saturation, assuming adiabatic instability
            omegat = difft*pow2(sqrt(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k))) !Thermal damping rate assuming adiabatic instability
            brunts = sqrt(abs( s% brunt_N2_composition_term(k)&
                +(s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0 + omegat/omegaa) )) !Suppress thermal part of brunt
            bruntsn2 = sqrt(abs( s% brunt_N2_composition_term(k)+&
            (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft) )) !Effective brunt for isothermal instability
            brunts = max(brunts,bruntsn2) !Choose max between suppressed brunt and isothermal brunt
            brunts = max(s% omega(k),brunts) !Don't let Brunt be smaller than omega
            omegaa = s% omega(k)*pow(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0) !Recalculate omegaa

            ! Calculate nu_TSF
            if (s% brunt_N2(k) > 0d0) then
                if (pow2(brunts) > 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    omegac = 1d0*s% omega(k)*sqrt(brunts/s% omega(k))*pow(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)  !Critical field strength
                    nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac)) !Suppress AM transport if omega_a<omega_c
                    nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts) !nu_omega for revised Tayler instability
                end if
                ! Add TSF enabled by thermal diffusion
                if (pow2(brunts) < 2d0*pow2(vals(k,2))*pow2(s% omega(k))) then
                    nu_tsf_t = alpha*abs(vals(k,2))*s% omega(k)*pow2(s% r(k))
                end if
                vals(k,3) = max(nu_tsf,nu_tsf_t) + 1d-1
                vals(k,4) = safe_log10(vals(k,3))
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
            vals(k,3) = pow(10d0,vals(k,4))
         end do

          !Values near inner boundary
         do k=s% nz-nsmootham,s% nz
            vals(k,4) = vals(k-1,4)
            vals(k,3) = pow(10d0,vals(k,4))
         end do

          !Values near outer boundary
         do k=nsmootham,1,-1
            vals(k,4) = vals(k+1,4)
            vals(k,3) = pow(10d0,vals(k,3))
         end do


      end subroutine TS_data_for_extra_profile_columns



      real(dp) function diffmag(rho,T,abar,zbar,ierr)

         ! Written by S.-C. Yoon, Oct. 10, 2003
         ! Electrical conductivity according to Spitzer 1962
         ! See also Wendell et al. 1987, ApJ 313:284
         real(dp), intent(in) :: rho, T, abar, zbar
         integer, intent(out) :: ierr
         real(dp) :: xmagfmu, xmagft, xmagfdif, xmagfnu, &
            xkap, xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig, &
            xeta
    
        if (ierr /= 0) return

        xgamma = 0.2275d0*zbar*zbar*pow(rho*1d-6/abar,1d0/3d0)*1d8/T
        xlg = log10(xgamma)
        if (xlg < -1.5d0) then
            xsig1 = sige1(zbar,T,xgamma)
            xsig = xsig1
        else if (xlg >= -1.5d0 .and. xlg <= 0d0) then
            xxx = (xlg + 0.75d0)*4d0/3d0
            ffff = 0.25d0*(2d0-3d0*xxx + xxx*xxx*xxx)
            xsig1 = sige1(zbar,T,xgamma)
            
            xsig2 = sige2(T,rho,zbar,ierr)
            if (ierr /= 0) return
      
            xsig = (1d0-ffff)*xsig2 + ffff*xsig1
        else if (xlg > 0d0 .and. xlg < 0.5d0) then
            xsig2 = sige2(T,rho,zbar,ierr)
            if (ierr /= 0) return

            xsig = xsig2
        else if (xlg >= 0.5d0 .and. xlg < 1d0) then
            xxx = (xlg-0.75d0)*4d0
            ffff = 0.25d0*(2d0-3d0*xxx + xxx*xxx*xxx)
            xsig2 = sige2(T,rho,zbar,ierr)
            if (ierr /= 0) return

            xsig3 = sige3(zbar,T,xgamma)
            xsig = (1d0-ffff)*xsig3 + ffff*xsig2
        else
            xsig3 = sige3(zbar,T,xgamma)
            xsig = xsig3
        endif

        diffmag = 7.1520663d19/xsig ! magnetic diffusivity

      end function diffmag


      ! Helper functions

      real(dp) function sige1(z,t,xgamma)
         ! Written by S.-C. Yoon, Oct. 10, 2003
         ! Electrical conductivity according to Spitzer 1962
         ! See also Wendell et al. 1987, ApJ 313:284
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: etan, xlambda,f
         if (t >= 4.2d5) then
            f = sqrt(4.2d5/t)
         else
            f = 1d0
         end if
         xlambda = sqrt(3d0*z*z*z)*pow(xgamma,-1.5d0)*f + 1d0
         etan = 3d11*z*log(xlambda)*pow(t,-1.5d0)             ! magnetic diffusivity
         etan = etan/(1d0-1.20487d0*exp(-1.0576d0*pow(z,0.347044d0))) ! correction: gammae
         sige1 = clight*clight/(4d0*pi*etan)                    ! sigma = c^2/(4pi*eta)
      end function sige1


      real(dp) function sige2(T,rho,zbar,ierr)
         ! writen by S.-C. YOON Oct. 10, 2003
         ! electrical conductivity using conductive opacity
         ! see Wendell et al. 1987 ApJ 313:284
         use kap_lib, only: kap_get_elect_cond_opacity
         real(dp), intent(in) :: t,rho,zbar
         integer, intent(out) :: ierr
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         call kap_get_elect_cond_opacity( &
            zbar, log10(rho), log10(T),  &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         sige2 = 1.11d9*T*T/(rho*kap)
      end function sige2

      real(dp) function sige3(z,t,xgamma)
         ! writen by S.-C. YOON Oct. 10, 2003
         ! electrical conductivity in degenerate matter,
         ! according to Nandkumar & Pethick (1984)
         real(dp), intent(in) :: z, t, xgamma
         real(dp) :: rme, rm23, ctmp, xi
         rme = 8.5646d-23*t*t*t*xgamma*xgamma*xgamma/pow5(z)  ! rme = rho6/mue
         rm23 = pow(rme,2d0/3d0)
         ctmp = 1d0 + 1.018d0*rm23
         xi= sqrt(3.14159d0/3d0)*log(z)/3d0 + 2d0*log(1.32d0+2.33d0/sqrt(xgamma))/3d0-0.484d0*rm23/ctmp
         sige3 = 8.630d21*rme/(z*ctmp*xi)
      end function sige3
