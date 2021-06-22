subroutine wd_builder(id, ierr)

  use num_lib, only: look_for_brackets, safe_root

  integer, intent(in) :: id
  integer, intent(out) :: ierr
  type (star_info), pointer :: s

  integer, parameter :: wd_lipar = 1
  integer, pointer :: ipar(:)
  integer :: wd_lrpar
  real(dp), pointer :: rpar(:)

  real(dp), pointer :: xh(:,:), q(:), dq(:)

  integer :: i, j, k, nz

  real(dp) :: mstar, mstar1, rstar, T_c, rho_c, L_core

  real(dp) :: lnd, dlnd, lnd1, lnd3, y1, y3, epsx, epsy

  integer, parameter :: imax = 100

  include 'formats'

  ierr = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) return

  ! set mass variables
  mstar = s% initial_mass * Msun
  s% mstar = mstar
  s% star_mass = mstar/Msun
  s% xmstar = mstar

  ! set center variables (usual r = 0 conditions)
  s% M_center = 0
  s% L_center = 0
  s% R_center = 0
  s% v_center = 0

  ! fix central temperature
  ! needs to be such that you can have a degenerate, isothermal core
  T_c = s% x_ctrl(i_Tc)

  ! pick a luminosity for the core; will impose L(m) = Lcore * m
  ! this won't be the final L, MESA is just happier if gradT != 0
  L_core = s% x_ctrl(i_Lc) * Lsun

  ! rough guess for initial central density
  rho_c = get_rho_c(s% initial_mass)

  if (dbg) then
     write(*,*) 'wd_builder'
     write(*,*)
     write(*,1) 'T_c [K]', T_c
     write(*,1) 'L_core/Lsun', L_core/Lsun
     write(*,1) 'mstar/Msun', mstar/Msun
     write(*,*)
  end if

  wd_lrpar = 3
  allocate(rpar(wd_lrpar))
  i = 1 ! rpar(1) for mstar result
  rpar(i+1) = T_c; i = i+1
  rpar(i+1) = L_core; i = i+1

  if (i /= wd_lrpar) then
     write(*,*) 'i /= wd_lrpar', i, wd_lrpar
     write(*,*) 'wd'
     ierr = -1
     return
  end if

  allocate(ipar(wd_lipar))
  ipar(1) = id

  lnd = log(rho_c)
  dlnd = ln2

  ! given a central density, build1_wd_model creates a WD with the
  ! specified composition and temperature structure and returns the
  ! mass.  wd_f wraps this routine and returns the relative difference
  ! between the mass and the target mass, so we then do a root find of
  ! wd_f in central density of to get a WD of the specified mass.

  ! first, bracket the solution so that we can do a safe root find

  call look_for_brackets(lnd, dlnd, lnd1, lnd3, wd_f, y1, y3, &
       imax, wd_lrpar, rpar, wd_lipar, ipar, ierr)
  if (ierr /= 0) then
     if (dbg) then
        if (dbg) write(*,*) 'look_for_brackets ierr', ierr
        write(*,1) 'lnd1', lnd1
        write(*,1) 'lnd3', lnd3
        write(*,1) 'y1', y1
        write(*,1) 'y3', y3
     end if
     return
  end if

  if (dbg) then
     write(*,*)
     write(*,*) 'Bracketing'
     write(*,1) 'rho_c (initial guess)', rho_c
     write(*,1) 'm (target mass)', mstar/Msun
     write(*,1) 'rho_c (lower bracket)', exp(lnd1)
     write(*,1) 'rho_c (upper bracket)', exp(lnd3)
     write(*,*)
  end if

  ! now do the root find

  epsx = 1d-6 ! limit for variation in lnd
  epsy = 1d-6 ! limit for matching desired mass as fraction of total mass

  lnd = safe_root(wd_f, lnd1, lnd3, y1, y3, imax, epsx, epsy, &
       wd_lrpar, rpar, wd_lipar, ipar, ierr)
  if (ierr /= 0) then
     if (dbg) write(*,*) 'safe_root ierr', ierr
     return
  end if

  mstar1 = rpar(1)

  xh => s% xh
  q => s% q
  dq => s% dq
  nz = s% nz

  if (dbg) then
     write(*,*)
     write(*,*) 'Found root'
     write(*,1) 'rho_c (final value)', exp(xh(s% i_lnd,nz))
     write(*,1) 'mstar1/Msun', mstar1/Msun
     write(*,1) '(mstar-mstar1)/mstar', (mstar-mstar1)/mstar
     write(*,1) 'log10(r/Rsun)', log10(exp(xh(s% i_lnR,1))/Rsun)
     write(*,1) 'log10(L/Lsun)', log10(xh(s% i_lum,1)/Lsun)
     write(*,1) 'Tsurf (K)', exp(xh(s% i_lnT,1))
     write(*,2) 'nz', nz
     write(*,*)
  end if

  write(*,*) 'finished wd_builder model'

  ! save the data associated with final root-find model

  call star_allocate_arrays(id, ierr)
  if (ierr /= 0) then
     return
  end if

  do k=1,nz
     do j=1,s% nvar_hydro
        s% xh(j,k) = xh(j,k)
     end do
     s% q(k) = q(k)
     s% dq(k) = dq(k)
     call wd_builder_get_xa(s, s% q(k), s% xa(:, k))
  end do

  deallocate(xh, q, dq, ipar, rpar)

end subroutine wd_builder


real(dp) function wd_f(lnd, dfdx, lrpar, rpar, lipar, ipar, ierr)
  integer, intent(in) :: lrpar, lipar
  real(dp), intent(in) :: lnd
  real(dp), intent(out) :: dfdx
  integer, intent(inout), pointer :: ipar(:) ! (lipar)
  real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
  integer, intent(out) :: ierr

  type (star_info), pointer :: s
  real(dp) :: rho_c, T_c, L_core, d_log10_P
  real(dp), pointer :: xa(:)
  integer :: i, nz
  real(dp) :: mstar, mstar1

  logical, parameter :: dbg = .false.

  include 'formats'

  ierr = 0
  wd_f = 0
  if (lipar <= 0) then
     write(*,*) 'lipar', lipar
     write(*,*) 'wd f'
     ierr = -1
     return
  end if

  call star_ptr(ipar(1), s, ierr)
  if (ierr /= 0) return

  rho_c = exp(lnd)

  i = 1 ! rpar(1) for mstar result
  T_c = rpar(i+1); i = i+1
  L_core = rpar(i+1); i = i+1
  if (i > lrpar) then
     write(*,*) 'i > lrpar', i, lrpar
     write(*,*) 'wd f'
     ierr = -1
     return
  end if

  mstar = s% mstar ! desired value
  mstar1 = mstar ! to keep gfortran quiet

  if (dbg) write(*,*) 'call build1_wd_model'
  call build1_wd_model(s, T_c, rho_c, L_core, nz, mstar1, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in build1_wd_model'
     return
  end if

  s% nz = nz

  rpar(1) = mstar1 ! return the actual mass

  wd_f = (mstar - mstar1) / mstar
  dfdx = 0

  if (dbg) then
     write(*,1) 'rho_c', rho_c
     write(*,1) 'mstar1', mstar1 / msun
     write(*,1) 'wd_f', wd_f
     write(*,*)
  end if

end function wd_f


subroutine build1_wd_model(s, T_c, rho_c, L_core, nz, mstar, ierr)
  use chem_def
  use eos_def
  use kap_def
  use kap_lib
  use chem_lib
  use eos_lib, only: Radiation_Pressure
  type (star_info), pointer :: s
  real(dp), intent(in) :: T_c, rho_c, L_core
  integer, intent(out) :: nz
  real(dp), intent(out) :: mstar ! the mass of the constructed model
  integer, intent(out) :: ierr

  real(dp) :: logRho_guess
  real(dp), parameter :: LOGRHO_TOL = 1E-6_dp
  real(dp), parameter :: LOGPGAS_TOL = 1E-6_dp

  integer :: i, ii, k, i_lnd, i_lnT, i_lnR, prune
  real(dp), parameter :: &
       dlogPgas = 0.01d0, q_at_nz = 1d-5

  real(dp) :: &
       logPgas, Prad, Pgas, logPgas0, &
       P_c, logP, m, gradT, &
       r, rmid, rho, logRho, T, lnT, L, P, P0, dm, m0, L0, r0, lnT0, T0, &
       rho0, rho_mid, Pmid, chiRho0, chiRho_mid, chiT0, chiT_mid, Cp0, Cp_mid, &
       grada0, grada_mid, mmid, Tmid, Lmid

  real(dp), pointer :: xh(:,:), q(:), dq(:) ! model structure info
  real(dp), allocatable :: xa(:)

  ! composition information
  integer :: species
  real(dp) :: x, y, z, abar, zbar, z53bar, z2bar, ye, mass_correction, sumx
  real(dp) :: XC, XN, XO, XNe

  ! eos results
  real(dp), dimension(num_eos_basic_results) :: res, d_eos_dlnd, d_eos_dlnT
  real(dp), allocatable :: d_eos_dxa(:,:)
  real(dp) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
       eta, d_eta_dlnRho, d_eta_dlnT, &
       chiRho, chiT, Cp, grada

  ! kap results
  real(dp) :: opacity, kap_cond, kap_rad, dlnkap_dlnd, dlnkap_dlnT
  real(dp) :: dlnkap_dxa(s% species)
  real(dp) :: kap_fracs(num_kap_fracs)

  ! limits
  real(dp) :: logT_surf_limit, logP_surf_limit, logRho_surf_limit

  real(dp), parameter :: cgrav = standard_cgrav

  logical, parameter :: dbg = .false.

  include 'formats'

  ierr = 0

  ! set limits

  logP_surf_limit = 8d0
  if (dbg) write(*,1) 'logP_surf_limit', logP_surf_limit

  logT_surf_limit = 4d0
  if (dbg) write(*,1) 'logT_surf_limit', logT_surf_limit

  logRho_surf_limit = -12d0
  if (dbg) write(*,1) 'logRho_surf_limit', logRho_surf_limit

  i_lnd = s% i_lnd
  i_lnT = s% i_lnT
  i_lnR = s% i_lnR

  allocate(xa(s% species))
  allocate(d_eos_dxa(num_eos_d_dxa_results, s% species))


  ! set composition and evaluate EOS at center

  call wd_builder_get_xa(s, 0d0, xa)
  call set_composition_info

  call star_get_eos( &
       s% id, 0, xa, &
       rho_c, log10(rho_c), T_c, log10(T_c), &
       res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in get_eos'
     return
  end if
  call unpack_eos_results

  logPgas = res(i_lnPgas)/ln10
  Pgas = exp10(logPgas)
  P_c = Pgas + Radiation_Pressure(T_c) ! center pressure

  mstar = s% mstar ! desired total mass
  m = q_at_nz*mstar ! mass at nz
  ! pressure at innermost point using K&W 10.6
  P = P_c - 3*cgrav/(8*pi)*pow(pi4*rho_c/3,4d0/3d0)*pow(m,two_thirds)
  logP = log10(P)

  ! estimate nz from lgP
  nz = 1 + (logP - logP_surf_limit)/dlogPgas

  ! temperature at nz assuming isothermal
  lnT = log(T_c)
  T = exp(lnT)

  ! density at nz
  logRho_guess = log10(rho_c)
  call star_solve_eos_given_PgasT( &
       s% id, 0, xa, &
       lnT/ln10, log10(Pgas), logRho_guess, LOGRHO_TOL, LOGPGAS_TOL, &
       logRho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, &
       ierr)
  if (ierr /= 0) return
  rho = exp10(logRho)
  call unpack_eos_results

  r = pow(m/(pi4*rho/3),one_third) ! radius at nz

  L = L_core * (m/mstar)

  call kap_get( &
       s% kap_handle, s% species, s% chem_id, s% net_iso, xa, log10(rho), log10(T), &
       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
       eta, d_eta_dlnRho, d_eta_dlnT, &
       kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)

  if (ierr /= 0) then
     write(*,*) 'failed in kap_get'
     return
  end if

  call eval_gradT( &
       s, xa, m, &
       r, T, L, P, opacity, rho, chiRho, chiT, Cp, grada, &
       gradT, ierr)
  if (ierr /= 0) return


  if (associated(s% xh)) deallocate(s% xh)
  if (associated(s% q)) deallocate(s% q)
  if (associated(s% dq)) deallocate(s% dq)

  allocate(xh(s% nvar_hydro,nz), q(nz), dq(nz), stat=ierr)
  if (ierr /= 0) return
  s% xh => xh
  s% dq => dq
  s% q => q

  xh(i_lnd, nz) = logRho*ln10
  xh(i_lnT, nz) = lnT
  xh(i_lnR, nz) = log(r)
  xh(s% i_lum,nz) = L

  q(nz) = q_at_nz
  dq(nz) = q_at_nz

  if (dbg) write(*,*) 'nz', nz

  prune = 0
  step_loop: do k = nz-1, 1, -1

     logPgas0 = logPgas
     P0 = P
     m0 = m
     L0 = L
     r0 = r
     lnT0 = lnT
     T0 = T
     rho0 = rho
     chiRho0 = chiRho
     chiT0 = chiT
     Cp0 = Cp
     grada0 = grada
     dm = 0 ! for gfortran

     if (dbg) write(*,3) 'step', k, nz, logPgas0

     logPgas = logPgas0 - dlogPgas
     Pgas = exp10(logPgas)

     ! these inner i & ii loops follow the same approach as the pre-main-sequence builder

     do i = 1, 2

        Prad = Radiation_Pressure(T)
        P = Pgas + Prad

        rho_mid = (rho+rho0)/2

        do ii = 1, 10 ! repeat to get hydrostatic balance
           rmid = pow((r*r*r + r0*r0*r0)/2,one_third)
           mmid = (m + m0)/2
           if (ii == 10) exit
           dm = -pi4*pow4(rmid)*(P-P0)/(cgrav*mmid)
           m = m0 + dm ! mass at point k
           r = pow(r0*r0*r0 + dm/((4*pi/3)*rho_mid),one_third)
           if (dbg) write(*,2) 'r', ii, r, m, dm
        end do

        L = L0 + dm*(L_core/mstar) ! luminosity at point k
        Lmid = (L0+L)/2

        Pmid = (P+P0)/2

        chiRho_mid = (chiRho0 + chiRho)/2
        chiT_mid = (chiT0 + chiT)/2
        Cp_mid = (Cp0 + Cp)/2
        grada_mid = (grada0 + grada)/2

        call wd_builder_get_xa(s, m/mstar, xa)
        call set_composition_info

        do ii = 1, 2
           Tmid = (T+T0)/2

           call kap_get( &
                s% kap_handle, s% species, s% chem_id, s% net_iso, xa, log10(rho_mid), log10(Tmid), &
                lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                eta, d_eta_dlnRho, d_eta_dlnT, &
                kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
           if (ierr /= 0) then
              write(*,*) 'failed in kap_get'
              return
           end if

           call kap_get_elect_cond_opacity(s% kap_handle, &
                zbar, log10(rho_mid), log10(Tmid), &
                kap_cond, dlnkap_dlnd, dlnkap_dlnT, ierr)

           kap_rad = 1d0/(1d0/opacity - 1d0/kap_cond)

           if (kap_cond .lt. kap_rad) then ! suggested by Evan
              call eval_gradT( &
                   s, xa, mmid, &
                   rmid, Tmid, Lmid, Pmid, opacity, rho_mid, &
                   chiRho_mid, chiT_mid, Cp_mid, grada_mid, &
                   gradT, ierr )
              if (ierr /= 0) return
           else
              gradT = 0.235 ! suggested by Lars
           end if

           T = T0 + Tmid*gradT*(P-P0)/Pmid
           lnT = log(T)
           if (dbg) write(*,2) 'T', ii, T
        end do

        if (i == 2) exit

        logRho_guess = logRho
        call star_solve_eos_given_PgasT( &
             s% id, 0, xa, &
             lnT/ln10, logPgas, logRho_guess, LOGRHO_TOL, LOGPGAS_TOL, &
             logRho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dxa, &
             ierr)
        rho = exp10(logRho)
        if (ierr /= 0) return
        call unpack_eos_results

     end do

     ! check limits

     if (lnT/ln10 <= logT_surf_limit) then
        if (dbg) write(*,*) 'have reached logT_surf_limit', lnT/ln10, logT_surf_limit
        prune = k
        exit step_loop
     end if

     if (log10(P) <= logP_surf_limit) then
        if (dbg) write(*,1) 'have reached logP_surf limit', log10(P), logP_surf_limit
        prune = k
        exit step_loop
     end if

     if (logRho <= logRho_surf_limit) then
        if (dbg) write(*,1) 'have reached logRho_surf limit', logRho, logRho_surf_limit
        prune = k
        exit step_loop
     end if

     ! save data for cell k
     xh(i_lnd, k) = logRho*ln10
     xh(i_lnT, k) = lnT
     xh(i_lnR, k) = log(r)
     xh(s% i_lum,k) = L
     q(k) = m/mstar
     dq(k) = dm/mstar

     if (dbg) then
        write(*,2) 'xh(i_lnd, k)', k, xh(i_lnd, k)
        write(*,2) 'xh(i_lnT, k)', k, xh(i_lnT, k)
        write(*,2) 'xh(i_lnR, k)', k, xh(i_lnR, k)
        write(*,2) 'L', k, L
        write(*,2) 'q(k)', k, q(k)
        write(*,2) 'dq(k)', k, dq(k)
     end if

  end do step_loop

  ! reduce to the actual number of zones we needed
  if (prune > 0) then ! move stuff and reduce nz
     if (dbg) write(*,*) 'prune', prune
     do k=1,nz-prune
        xh(:,k) = xh(:,k+prune)
        q(k) = q(k+prune)
        dq(k) = dq(k+prune)
     end do
     m = mstar*q(1)
     nz = nz-prune
     if (dbg) write(*,*) 'final nz', nz
  end if

  mstar = m ! actual total mass

  ! clean up the qs and dqs
  call star_normalize_dqs(s% id, nz, dq, ierr)
  if (ierr /= 0) then
     if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
     return
  end if
  call star_set_qs(s% id, nz, q, dq, ierr)
  if (ierr /= 0) then
     if (s% report_ierr) write(*,*) 'set_qs failed in pre ms model'
     return
  end if

contains

  subroutine set_composition_info

    call basic_composition_info( &
         s% species, s% chem_id, xa, x, y, z, abar, zbar, z53bar, z2bar, ye, &
         mass_correction, sumx)

    xc = xa(s% net_iso(ic12))
    xn = xa(s% net_iso(in14))
    xo = xa(s% net_iso(io16))
    xne = xa(s% net_iso(ine20))

  end subroutine set_composition_info


  subroutine unpack_eos_results
    chiRho = res(i_chiRho)
    chiT = res(i_chiT)
    Cp = res(i_cp)
    grada = res(i_grad_ad)
    lnfree_e = res(i_lnfree_e)
    d_lnfree_e_dlnRho = d_eos_dlnd(i_lnfree_e)
    d_lnfree_e_dlnT = d_eos_dlnT(i_lnfree_e)
    eta = res(i_eta)
    d_eta_dlnRho = d_eos_dlnd(i_eta)
    d_eta_dlnT = d_eos_dlnT(i_eta)
  end subroutine unpack_eos_results


end subroutine build1_wd_model


! get an estimate for gradT from MLT
subroutine eval_gradT( &
     s, xa, m, &
     r, T, L, P, opacity, rho, chiRho, chiT, Cp, grada, &
     gradT, ierr)
  use chem_def, only: ih1
  use auto_diff
  use eos_lib, only: Radiation_Pressure

  type (star_info), pointer :: s
  real(dp), intent(in) :: xa(:), m, &
       r, T, L, P, opacity, rho, chiRho, chiT, Cp, grada
  real(dp), intent(out) :: gradT
  integer, intent(out) :: ierr

  ! mlt inputs
  character (len=16) :: MLT_option
  integer :: iso
  real(dp) :: X, gradr, scale_height, &
       gradL_composition_term, mixing_length_alpha

  ! mlt outputs
  integer :: mixing_type

  real(dp), parameter :: cgrav = standard_cgrav

  ierr = 0

  ! set passed values; don't care about derivatives
  gradr = P*opacity*L/(16d0*pi*clight*m*cgrav*Radiation_Pressure(T))
  scale_height = P/((cgrav*m/(r*r))*rho)

  ! set some mixing options
  ! details are irrelevant, so just pick something simple
  MLT_option = 'Henyey'
  gradL_composition_term = 0
  mixing_length_alpha = 1.8
  iso = ih1 ! doesn't matter since no thermohaline
  X = xa(ih1)

  call star_mlt_gradT(s% id, MLT_option, &
       r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
       iso, X, cgrav, m, gradL_composition_term, mixing_length_alpha, &
       gradT, mixing_type, ierr)

end subroutine eval_gradT


! estimate central density of a WD as a function of mass
! this doesn't need to be a great guess, just vaguely reasonable
real(dp) function get_rho_c(m)
  real(dp), intent(in) :: m
  real(dp), parameter, dimension(5) :: ms = (/0.15, 0.40, 0.82, 1.37, 1.43/)
  real(dp), parameter, dimension(5) :: logRhos = (/5.0, 6.0, 7.0, 8.0, 9.0/)
  integer :: k

  ! less than min ms
  if (m < ms(1)) then
     get_rho_c = exp10(logRhos(1))
     return
  end if

  ! simple log-linear interpolation
  do k = 1, 4
     if (m < ms(k+1)) then
        get_rho_c = exp10(logRhos(k) + (m - ms(k)) / (ms(k+1) - ms(k)))
        return
     end if
  end do

  ! greater than max ms
  get_rho_c = exp10(logRhos(5))

end function get_rho_c
