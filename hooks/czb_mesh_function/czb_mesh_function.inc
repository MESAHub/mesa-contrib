subroutine read_inlist_czb_mesh_function(ierr)
  integer, intent(out) :: ierr
  character (len=256) :: filename, message
  integer :: unit         
  filename = 'inlist_czb_mesh_function'

  write(*,*) 'read_inlist_czb_mesh_function'

  ! set defaults
  czb_mesh_zone_type(:) = '' !  Possible values: burn_H, burn_He, burn_Z, nonburn, any
  czb_mesh_zone_loc(:) = ''  !  Possible values: core, shell, any
  czb_mesh_bdy_loc(:) = ''   !  Possible values: bottom, top, any`

  czb_mesh_f(:) = 0d0
  czb_mesh_l(:) = 0d0

  open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
  if (ierr /= 0) then
     write(*, *) 'Failed to open control namelist file ', trim(filename)
  else
     read(unit, nml=czb_mesh_function_nml, iostat=ierr)  
     close(unit)
     if (ierr /= 0) then
        write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
        write(*, '(a)') &
             'The following runtime error message might help you find the problem'
        write(*, *) 
        open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
        read(unit, nml=czb_mesh_function_nml)
        close(unit)
     end if
  end if

end subroutine read_inlist_czb_mesh_function


subroutine czb_mesh_function (s, fcn, ierr)

  type(star_info), pointer :: s
  real(dp), intent(inout)    :: fcn(:)
  integer, intent(out)     :: ierr

  logical, parameter :: DEBUG = .TRUE.

  integer  :: i
  integer  :: j
  integer  :: k
  logical  :: is_core
  logical  :: match_zone_type
  logical  :: match_zone_loc
  logical  :: match_bdy_loc

  include 'formats'

  ! Initialize

  ierr = 0
  fcn = 1d0

  if (DEBUG) then
     write(*, 3) 'conv_bdy_mesh_function; model, n_conv_bdy=', s%model_number, s%num_conv_boundaries
  end if

  ! Loop over convective boundaries, from center to surface

  conv_bdy_loop : do i = 1, s%num_conv_boundaries

     if (DEBUG) then
        write(*,*) 'zone ', i
        write(*,*) 'burn_H', s%burn_h_conv_region(i)
        write(*,*) 'burn_He', s%burn_he_conv_region(i)
        write(*,*) 'burn_Z', s%burn_Z_conv_region(i)
        write(*,*) 'is_core', (i == 1 .AND. s%R_center == 0d0 .AND. s%top_conv_bdy(i))
        write(*,*) 'top', s%top_conv_bdy(i)
     end if

     ! Loop over criteria

     criteria_loop : do j = 1, NUM_CZB_MESH_PARAM_SETS

        ! Skip if the zone type is blank

        if (czb_mesh_zone_type(j) == '') cycle criteria_loop

        ! Check if the criteria match the current boundary

        select case (czb_mesh_zone_type(j))
        case ('burn_H')
           match_zone_type = s%burn_h_conv_region(i)
        case ('burn_He')
           match_zone_type = s%burn_he_conv_region(i)
        case ('burn_Z')
           match_zone_type = s%burn_z_conv_region(i)
        case ('nonburn')
           match_zone_type = .NOT. ( &
                s%burn_h_conv_region(i) .OR. &
                s%burn_he_conv_region(i) .OR. &
                s%burn_z_conv_region(i) )
        case ('any')
           match_zone_type = .TRUE.
        case default
           write(*,*) 'Invalid czb_mesh_zone_type: j, czb_mesh_zone_type(j)=', j, czb_mesh_zone_type(j)
           ierr = -1
           return
        end select

        is_core = (i == 1 .AND. s%R_center == 0d0 .AND. s%top_conv_bdy(i))

        select case (czb_mesh_zone_loc(j))
        case ('core')
           match_zone_loc = is_core
        case ('shell')
           match_zone_loc = .NOT. is_core
        case ('any')
           match_zone_loc = .TRUE.
        case default
           write(*,*) 'Invalid czb_mesh_zone_loc: j, czb_mesh_zone_loc(j)=', j, czb_mesh_zone_loc(j)
           ierr = -1
           return
        end select

        select case (czb_mesh_bdy_loc(j))
        case ('bottom')
           match_bdy_loc = .NOT. s%top_conv_bdy(i)
        case ('top')
           match_bdy_loc = s%top_conv_bdy(i)
        case ('any')
           match_bdy_loc = .TRUE.
        case default
           write(*,*) 'Invalid czb_mesh_bdy_loc: j, czb_mesh_bdy_loc(j)=', j, czb_mesh_bdy_loc(j)
           ierr = -1
           return
        end select

        if (.NOT. (match_zone_type .AND. match_zone_loc .AND. match_bdy_loc)) cycle criteria_loop

        if (DEBUG) then
           write(*,*) 'czb_mesh at convective boundary: i, j=', i, j
           write(*,*) '  czb_mesh_zone_type=', TRIM(czb_mesh_zone_type(j))
           write(*,*) '  czb_mesh_zone_loc=', TRIM(czb_mesh_zone_loc(j))
           write(*,*) '  czb_mesh_bdy_loc=', TRIM(czb_mesh_bdy_loc(j))
        endif

        ! Evaluate the mesh function
        call eval_czb_mesh_fcn(s, i, j, fcn, ierr)

        if (ierr /= 0) return

     end do criteria_loop

  end do conv_bdy_loop

  ! Finish
  fcn(1) = s% lnPeos(1)
  do k = 2, s% nz
     fcn(k) = fcn(k-1) + (s% lnPeos(k) - s% lnPeos(k-1)) * fcn(k)
  end do

  return

end subroutine czb_mesh_function

subroutine eval_czb_mesh_fcn (s, i, j, fcn, ierr)

  type(star_info), pointer :: s
  integer, intent(in)      :: i
  integer, intent(in)      :: j
  real(dp), intent(inout)  :: fcn(:)
  integer, intent(out)     :: ierr

  logical, parameter :: DEBUG = .FALSE.

  real(dp) :: f
  real(dp) :: l
  integer  :: k_cb
  real(dp) :: r_cb
  real(dp) :: Hp_cb
  integer  :: k
  real(dp) :: r
  real(dp) :: dr


  ierr = 0

  ! Extract parameters

  f = czb_mesh_f(j)
  l = czb_mesh_l(j)

  ! Evaluate convective boundary (_cb) parameters

  call eval_conv_bdy_k(s, i, k_cb, ierr)
  if (ierr /= 0) return

  call eval_conv_bdy_r(s, i, r_cb, ierr)
  if (ierr /= 0) return

  call eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
  if (ierr /= 0) return

  ! Loop over cells; first outward, then inward

  outward_cell_loop : do k = k_cb, 1, -1

     r = s%r(k)
     dr = r - r_cb

     fcn(k) = max(fcn(k), 1d0/f)

     if (dr > Hp_cb * l) then
        exit outward_cell_loop
     endif

  end do outward_cell_loop

  inward_cell_loop : do k = k_cb+1, s% nz, 1

     r = s%r(k)
     dr = r_cb - r

     fcn(k) = max(fcn(k), 1d0/f)

     if (dr > Hp_cb * l) then
        exit inward_cell_loop
     endif

  end do inward_cell_loop


  if (DEBUG) then
     write(*,*) 'conv_bdy_mesh:'
     write(*,*) '  f          =', f
     write(*,*) '  l          =', l
     write(*,*) '  Hp_cb      =', Hp_cb
  end if

  ! Finish

  return

end subroutine eval_czb_mesh_fcn


!****
! copy of routines from star/private/overshoot_utils.f90
!****

  subroutine eval_conv_bdy_k (s, i, k, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    integer, intent(out)     :: k
    integer, intent(out)     :: ierr

    ! Evaluate the index k of the cell containing the i'th convective
    ! boundary

    ierr = 0
    
    if (s%top_conv_bdy(i)) then
       k = s%conv_bdy_loc(i)
    else
       k = s%conv_bdy_loc(i) - 1
    endif

    if (k >= s%nz .OR. k < 1) then
       write(*,*) 'Invalid cell for convective boundary: i, k, nz=', i, k, s%nz
       ierr = -1
       return
    endif

    ! Finish

    return

  end subroutine eval_conv_bdy_k

  !****

  subroutine eval_conv_bdy_r (s, i, r, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(out)    :: r
    integer, intent(out)     :: ierr

    integer  :: k
    real(dp) :: w

    ! Evaluate the radius r at the i'th convective boundary

    ! Find the convective boundary cell

    ierr = 0

    call eval_conv_bdy_k(s, i, k, ierr)
    if (ierr /= 0) return

    ! Interpolate r based on the fact that r^3 varies linearly with q
    ! across the (constant-density) cell

    w = s%cz_bdy_dq(k)/s%dq(k)

    if (w < 0._dp .OR. w > 1._dp) then
       write(*,*) 'Invalid weight for convective boundary: i, w=', i, w
       ierr = -1
       return
    end if

    associate (k_o => k, &
               k_i => k+1)

      r = pow((      w)*s%r(k_i)*s%r(k_i)*s%r(k_i) + &
                 (1._dp-w)*s%r(k_o)*s%r(k_o)*s%r(k_o), 1._dp/3._dp)

    end associate

    ! Finish

    return

  end subroutine eval_conv_bdy_r

  !****

  subroutine eval_conv_bdy_Hp (s, i, Hp, ierr)

    use num_lib
    
    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(out)    :: Hp
    integer, intent(out)     :: ierr

    integer  :: k
    real(dp) :: r
    real(dp) :: w
    real(dp) :: x0
    real(dp) :: x1
    real(dp) :: x2
    real(dp) :: x
    real(dp) :: a0
    real(dp) :: a1
    real(dp) :: a2
    real(dp) :: P
    real(dp) :: rho
    real(dp) :: r_top
    real(dp) :: r_bot
    real(dp) :: dr

    ! Evaluate the pressure scale height Hp at the i'th convective boundary

    ! Find the convective boundary cell

    ierr = 0

    call eval_conv_bdy_k(s, i, k, ierr)
    if (ierr /= 0) return

    ! Evaluate the radius at the convective boundary

    call eval_conv_bdy_r(s, i, r, ierr)
    if (ierr /= 0) return

    ! Interpoalte the pressure and density at the boundary, using a
    ! quadratic fit across the boundary cell and its neighbors (the
    ! x's are fractional mass distances from the outer edge of cell
    ! k-1); then, evaluate the pressure scale height

    associate (k_o => k-1, &
               k_m => k, &
               k_i => k+1)

      x0 = s%dq(k_o)/2._dp
      x1 = s%dq(k_o) + s%dq(k_m)/2._dp
      x2 = s%dq(k_o) + s%dq(k_m) + s%dq(k_i)/2._dp

      x = s%dq(k_o) + s%cz_bdy_dq(k)

      call two_piece_linear_coeffs(x, x0, x1, x2, a0, a1, a2, ierr)
      if (ierr /= 0) return

      P = exp(a0*s%lnPeos(k_o) + a1*s%lnPeos(k_m) + a2*s%lnPeos(k_i))
      rho = exp(a0*s%lnd(k_o) + a1*s%lnd(k_m) + a2*s%lnd(k_i))

      ! Evaluate the pressure scale height

      Hp = P/(rho*s%cgrav(k_m)* &
           (s%M_center + s%xmstar*s%conv_bdy_q(i))/(r*r))

    end associate

    ! (Possibly) limit the scale height using the size of the
    ! convection zone

    if (s%limit_overshoot_Hp_using_size_of_convection_zone) then

       ! Determine the radial extent of the convection zone (note that
       ! r_top/r_bot don't coincide exactly with the r calculated
       ! above)

       if (s%top_conv_bdy(i)) then

          if (i == 1) then
             r_bot = s%R_center
          else
             if (s%top_conv_bdy(i-1)) then
                write(*,*) 'Double top boundary in overshoot; i=', i
                ierr = -1
                return
             end if
             r_bot = s%r(s%conv_bdy_loc(i-1))
          endif

          r_top = s%r(k)

       else

          r_bot = s%r(k+1)

          if (i == s%num_conv_boundaries) then
             r_top = s%r(1)
          else
             if (.NOT. s%top_conv_bdy(i+1)) then
                write(*,*) 'Double bottom boundary in overshoot; i=', i
                ierr = -1
                return
             endif
             r_top = s%r(s%conv_bdy_loc(i+1))
          endif

       endif
          
       dr = r_top - r_bot

       ! Apply the limit

       if (s%overshoot_alpha > 0d0) then
          if (s%overshoot_alpha*Hp > dr) Hp = dr/s%overshoot_alpha
       else
          if (s%alpha_mlt(k)*Hp > dr) Hp = dr/s%mixing_length_alpha
       end if

    end if

    ! Finish

    return

  end subroutine eval_conv_bdy_Hp
