  ! Renamed variable w.r.t. original run_star_extras.f90
  ! lxtra(1) is doing_DBP
  ! xtra(1) is m_core
  ! xtra(2) is mass_PZ
  ! xtra(3) is delta_r_PZ
  ! xtra(4) is alpha_PZ
  ! xtra(5) is r_core
  ! xtra(6) is rho_core_top
  ! ixtra(1) is k_PZ_top
  ! ixtra(2) is k_PZ_bottom


  subroutine other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    real(dp) :: fraction, Peclet_number, diffusivity    ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
    integer :: k
    logical, parameter :: DEBUG = .false.

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (s%D_mix(1) /= s%D_mix(1)) return  ! To ignore iterations where Dmix and gradT are NaNs

    if (s%num_conv_boundaries <= 1) then  ! Is zero at initialisation of the run
       if (DEBUG) then
          write(*,*) 'runstarex_gradT: skip since there are no convective boundaries'
       end if
       return
    endif

    do k= s%nz, 1, -1
       if (s%D_mix(k) <= s% min_D_mix) exit

       diffusivity = 16.0_dp * boltz_sigma * pow3(s% T(k)) / ( 3.0_dp * s% opacity(k) * pow2(s% rho(k)) * s% cp(k) )
       Peclet_number = s% conv_vel(k) * s% scale_height(k) * s% mixing_length_alpha / diffusivity

       if (Peclet_number >= 100.0_dp) then
          fraction = 1.0_dp
       else if (Peclet_number <= 0.01_dp) then
          fraction = 0.0_dp
       else
          fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
       end if

       s% adjust_mlt_gradT_fraction(k) = fraction
    end do

  end subroutine other_adjust_mlt_gradT_fraction_Peclet


  subroutine extended_convective_penetration(id, i, j, k_a, k_b, D, vc, ierr)
    integer, intent(in) :: id, i, j
    integer, intent(out) :: k_a, k_b
    real(dp), intent(out), dimension(:) :: D, vc
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    logical, parameter :: DEBUG = .false.
    real(dp) :: f, f2, f0
    real(dp) :: D0, Delta0
    real(dp) :: w
    real(dp) :: factor
    real(dp) :: r_cb, Hp_cb
    real(dp) :: r_ob, D_ob, vc_ob
    logical  :: outward
    integer  :: dk, k, k_ob
    real(dp) :: r, dr, r_step

    ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
    ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
    ! using the j'th set of overshoot parameters. The overshoot
    ! follows the extended convective penetration scheme description by Mathias
    ! Michielsen, "Probing the shape of the mixing profile and of the thermal
    ! structure at the convective core boundary through asteroseismology",
    ! A&A, 628, 76 (2019)

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if ((i /= 1) .or. (s%mixing_type(s%nz) /= convective_mixing)) then
       write(*,'(A,i2,A,i2)') 'ERROR: dissipation_balanced_penetration can only be used for core convection, &
            &so the first convective boundary. The routine got called for convective boundary number ',i, &
            &', and the mixing type in the core was', s%mixing_type(s%nz)
       ierr = -1
       return
    end if

    call dissipation_balanced_penetration(s, id) !, s%xtra(1), s% xtra(2), s% xtra(3), s% xtra(4), s% xtra(5), s% xtra(6))
    ! s% xtra(4) is distance from core boundary outward, so add f0 to it to make PZ zone reach that region
    s% xtra(4) = s% xtra(4) + s%overshoot_f0(j)
    ! Extract parameters
    f = s% xtra(4)                     ! extend of step function (a_ov)
    f0 = s%overshoot_f0(j)
    f2 = s%overshoot_f(j)            ! exponential decay (f_ov)

    D0 = s%overshoot_D0(j)
    Delta0 = s%overshoot_Delta0(j)

    if (f < 0.0_dp .or. f0 <= 0.0_dp .or. f2 < 0.0_dp) then
       write(*,*) 'ERROR: for extended convective penetration, must set f0 > 0, and f and f2 >= 0'
       write(*,*) 'see description of overshooting in star/defaults/control.defaults'
       ierr = -1
       return
    end if

    ! Evaluate convective boundary (_cb) parameters
    call star_eval_conv_bdy_r(s, i, r_cb, ierr)
    if (ierr /= 0) return

    call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
    if (ierr /= 0) return

    ! Evaluate overshoot boundary (_ob) parameters
    call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
    if (ierr /= 0) return

    ! Loop over cell faces, adding overshoot until D <= overshoot_D_min
    outward = s%top_conv_bdy(i)

    if (outward) then
       k_a = k_ob
       k_b = 1
       dk = -1
    else
       k_a = k_ob+1
       k_b = s%nz
       dk = 1
    endif

    if (f > 0.0_dp) then
       r_step = f*Hp_cb
    else
       r_step = 0.0_dp
    endif

    face_loop : do k = k_a, k_b, dk
       ! Evaluate the extended convective penetration factor
       r = s%r(k)
       if (outward) then
          dr = r - r_ob
       else
          dr = r_ob - r
       endif

       if (dr < r_step .AND. f > 0.0_dp) then  ! step factor
          factor = 1.0_dp
       else
          if ( f2 > 0.0_dp) then                ! exponential factor
             factor = exp(-2.0_dp*(dr-r_step)/(f2*Hp_cb))
          else
             factor = 0.0_dp
          endif
       endif

       ! Store the diffusion coefficient and velocity
       D(k) = (D0 + Delta0*D_ob)*factor
       vc(k) = (D0/D_ob + Delta0)*vc_ob*factor

       ! Check for early overshoot completion
       if (D(k) < s%overshoot_D_min) then
          k_b = k
          exit face_loop
       endif

    end do face_loop

    if (DEBUG) then
       write(*,*) 'step exponential overshoot:'
       write(*,*) '  k_a, k_b   =', k_a, k_b
       write(*,*) '  r_a, r_b   =', s%r(k_a), s%r(k_b)
       write(*,*) '  r_ob, r_cb =', r_ob, r_cb
       write(*,*) '  Hp_cb      =', Hp_cb
    end if

  end subroutine extended_convective_penetration



  subroutine dissipation_balanced_penetration(s, id)
    use eos_def
    use star_lib
    use kap_def
    type (star_info), pointer :: s
    integer, intent(in) :: id
    real(dp), parameter :: f = 0.86d0
    real(dp), parameter :: xi = 0.6d0
    integer :: k, j, ierr
    real(dp) :: Lint, V_CZ, Favg, RHS, dr, h, dLint
    real(dp) :: r_cb

    s% lxtra(1) = .true.
    V_CZ = 0d0
    Lint = 0d0

    call star_eval_conv_bdy_k(s, 1, k, ierr)
    call star_eval_conv_bdy_r(s, 1, r_cb, ierr)
    s% ixtra(2) = k
    s% xtra(5) = r_cb
    s% xtra(1) = s%m(k)
    s% xtra(6) = s%rho(k)
    h = s%scale_height(k)

    ! prescription based on Jermyn A. et al (2022)  https://arxiv.org/pdf/2203.09525.pdf
    ! Equation A1 is used here.
    ! RHS refers to right-hand side of equation A1, and Lint the integrated
    ! luminosity on either left or right side of that equation

    ! Integrate over cells that are fully in CZ
    ! r and L_conv are face values, assume they change linear within cell
    do j=s%nz,k+1,-1
       if (j .eq. s%nz) then
          dr = s%r(j)
          Lint = Lint + s%L_conv(j)*0.5 * dr
       else
          dr = s%r(j) - s%r(j+1)
          Lint = Lint + (s%L_conv(j+1) + s%L_conv(j))*0.5 * dr
       endif
    end do

    ! Take cell that is partially convective
    ! convective part of cell k
    ! L_conv goes to 0 at edge of conv zone
    dr = r_cb - s%r(k+1)
    Lint = Lint + s%L_conv(k+1)*0.5 * dr

    ! Calculate target RHS
    V_CZ = 4d0/3d0 * pi * r_cb*r_cb*r_cb
    Favg = Lint / V_CZ
    RHS = (1d0 - f) * Lint
    Lint = 0d0

    ! Integrate over RZ until we find the edge of the PZ
    ! remainder of cell k (non-convective part)
    ! Do integration explicitely, f*xi*4*pi*r^2 is moved outside of integral
    dr = s%r(k) - r_cb
    dLint = xi * f * 4d0 * pi * (pow3(s%r(k))-pow3(r_cb))/3 * Favg + (s%L(k)*0.5 * (s%grada_face(k) / s%gradr(k) - 1d0)) * dr
    Lint = dLint

    ! If remainder of cell k would already satisfy Lint > RHS
    if (Lint > RHS) then
       dr = dr*(Lint - RHS)/dLint
       s% xtra(2) =  s%rho(k) * 4d0/3d0 * pi * (pow3(r_cb+dr) - pow3(r_cb)) !s%m(k) - s%xtra(1) !only used for history output
       s% xtra(3) = dr
       s% xtra(4) = s% xtra(3) / h
       s% ixtra(1) = k
       return
    end if
    ! Else calculate dL_int for each cell, untill the total integrated L > RHS
    do j=k-1,1,-1
       dr = s%r(j) - s%r(j+1)
       dLint = xi * f * 4d0 * pi * (pow3(s%r(j))-pow3(s%r(j+1)))/3 * Favg &
            + ( (s%L(j+1)*(s%grada_face(j+1) / s%gradr(j+1) - 1d0) +(s%L(j)*(s%grada_face(j) / s%gradr(j) - 1d0)) )*0.5 * dr)

       if (Lint + dLint > RHS) then
          dr = dr*(RHS - Lint)/dLint
          s% xtra(2) = s%m(j) - s%xtra(1) !only used for history output
          s% xtra(3) = s%r(j+1)+dr - r_cb
          s% xtra(4) = s% xtra(3) / h
          s% ixtra(1) = j
          return
       end if
       Lint = Lint + dLint
    end do

  end subroutine dissipation_balanced_penetration

  subroutine mesh_delta_coeff_core_boundary(id, ierr)
    integer, intent(in) :: id
    ! real(dp), intent(in), dimension(:) :: eps_h, eps_he, eps_z
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    logical, parameter :: dbg = .false.
    integer :: k, k_max_N2comp
    real(dp) :: Hp, r_lower, r_upper_ov, r_upper_BV
    real(dp) :: xtra_dist_below, xtra_dist_above_ov, xtra_dist_above_bv, xtra_coeff_mesh
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! set defaults
    xtra_dist_below = s% x_ctrl(1)
    xtra_dist_above_ov = s% x_ctrl(2)
    xtra_dist_above_bv = s% x_ctrl(3)
    xtra_coeff_mesh = s% x_ctrl(4)

    if (xtra_coeff_mesh == 1d0 .or. xtra_coeff_mesh < 0) return
    if (s% mixing_type(s% nz) /= convective_mixing) return  ! only enable for convective cores

    r_upper_ov=-1
    r_upper_BV=-1
    k_max_N2comp = maxloc(s% brunt_n2_composition_term(:), 1)
    ! Find boundary of convective core, and go inwards by the specified distance (in Hp)
    do k = s% nz, 1, -1
       if (s% mixing_type(k) == convective_mixing) then
          continue
       else
          Hp = s% scale_height(k)  !s% P(k)/(s% rho(k)*s% grav(k))
          r_lower = max (s% r(s%nz), s% r(k) - xtra_dist_below*Hp)
          exit
       endif
    end do

    do k = s% nz, 1, -1
       ! Start increasing the mesh once closer than the given distance (in Hp) to the core boundary
       if (s%r(k) > r_lower) then
          if (xtra_coeff_mesh < s% mesh_delta_coeff_factor(k)) then
             s% mesh_delta_coeff_factor(k) = xtra_coeff_mesh
          end if
       else
          cycle
       endif

       ! Go up to the given distance past the overshoot boundary
       if (r_upper_ov<0 .and. s% mixing_type(k) /= overshoot_mixing .and. s% mixing_type(k) /= convective_mixing) then
          if (xtra_dist_above_ov > 0) then
             Hp = s% scale_height(k) !s% P(k)/(s% rho(k)*s% grav(k))
             r_upper_ov = min(s% r(1), s% r(k) + xtra_dist_above_ov*Hp)
          else
             r_upper_ov = 0
          end if
       end if

       ! Go up to the given distance past the order in magnitude decrease in BV composition term, outwards of its maximum
       if (r_upper_BV<0 .and. k < k_max_N2comp .and. (s% brunt_n2_composition_term(k)*10 < maxval(s% brunt_n2_composition_term(:))) ) then
          if (xtra_dist_above_bv > 0) then
             Hp = s% scale_height(k) !s% P(k)/(s% rho(k)*s% grav(k))
             r_upper_BV = min(s% r(1), s% r(k) + xtra_dist_above_bv*Hp)
          else
             r_upper_BV = 0
          end if
       end if

       ! Stop increasing mesh when further than the specified distance from both the overshoot boundary and BV composition peak
       if (s% r(k) > r_upper_ov .and. s% r(k) > r_upper_BV .and. r_upper_ov >=0 .and. r_upper_BV >=0) exit
    end do

  end subroutine mesh_delta_coeff_core_boundary
