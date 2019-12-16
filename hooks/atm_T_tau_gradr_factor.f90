! Modifies the radiative gradient so that the "interior" model
! correctly reproduces the T(τ) relation selected by the
! `atm_T_tau_relation` control.
!
! If the T(τ) relation is expressed as
!
!     T⁴(τ) = 3Teff⁴/4 (τ + q(τ))
!
! then the radiative gradient factor is 1+q'(τ).  For an Eddington
! atmosphere, q(τ)=2/3 so q'(τ)=0 and the factor isn't necessary.
!
! Added by @warrickball (2019-12-16).


      subroutine atm_T_tau_gradr_factor(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         do k = 1, s% nz
            s% gradr_factor(k) = 1d0 + dq_dtau(s% atm_T_tau_relation, s% tau(k))
         end do

      end subroutine atm_T_tau_gradr_factor


      real(dp) function dq_dtau(atm_T_tau_relation, tau)
         character (len=strlen), intent(in) :: atm_T_tau_relation
         real(dp), intent(in) :: tau
         real(dp) :: e1, e2, de1_dtau, de2_dtau, h1
         real(dp), parameter :: &
            q1 = 1.0361d0, &
            q2 = -0.3134d0, &
            q3 = 2.44799995d0, &
            q4 = -0.29589999d0, &
            q5 = 30.0d0

         select case (atm_T_tau_relation)
            case ('Eddington')
               dq_dtau = 0d0
            case ('Krishna_Swamy')
               if (tau == 0) then
                  e1 = 1; de1_dtau = -2.54d0
                  e2 = 1; de2_dtau = -30d0
               else
                  e1 = exp(-2.54d0*tau); de1_dtau = -2.54d0*e1
                  e2 = exp(-30d0*tau); de2_dtau = -30d0*e2
               end if
               dq_dtau = - 0.815d0*de1_dtau - 0.025d0*de2_dtau
            case ('solar_Hopf')
               if (tau == 0) then
                  e1 = 1; de1_dtau = -q3
                  e2 = 1; de2_dtau = -q5
               else
                  e1 = exp(-q3*tau); de1_dtau = -q3*e1
                  e2 = exp(-q5*tau); de2_dtau = -q5*e2
               end if
               dq_dtau = q2*de1_dtau + q4*de2_dtau
            case default
               dq_dtau = -1
         end select

      end function dq_dtau
