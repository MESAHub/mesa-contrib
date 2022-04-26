subroutine detach_binary(binary_id, ierr)
   ! This routine can be called in binary_extra_finish_step to
   ! continue the evolution of the secondary as a single star, e.g.
   ! because of the explosion of the companion and the disruption of
   ! the binary.
   use binary_lib, only: binary_set_separation_eccentricity
   type (binary_info), pointer :: b
   integer, intent(in) :: binary_id
   integer, intent(out) :: ierr ! ierr = 0 means success
   call binary_ptr(binary_id, b, ierr)
   if (ierr /= 0) return

   print *, '****************************************'
   print *, '* Switching from binary to single star *'
   print *, '****************************************'

   if (b% job% evolve_both_stars .eqv. .false.) then
      ! detaching a binary with a point mass would imply
      ! evolving two detached point masses
      ! probably not what you want to do
      ierr = 1
      STOP 'Cannot detach a binary with point mass companion'
      ! you can remove this STOP statement, add a return, and check
      ! for the value of ierr where you call this routine from instead
   end if

   b% job% evolve_both_stars = .false.
   b% ignore_hard_limits_this_step = .true.
   print *, 'current mtransfer_rate', b% mtransfer_rate, 'will be set to zero'
   b% mtransfer_rate = 0d0
   b% change_factor = b% max_change_factor
   print *, 'shutting down tides, accretion of J, magnetic braking, and missing wind'
   b% do_tidal_sync = .false.
   b% do_j_accretion = .false.
   b% do_jdot_mb = .false.
   b% do_jdot_missing_wind = .false.
   print *, 'ignore RLOF from now on'
   b% ignore_rlof_flag = .true.
   ! fix Eddington mdot stuff
   b% eq_initial_bh_mass = b% s_donor% m(1)
   b% mdot_edd = 1d99
   b% mdot_edd_eta = 0d0
   ! switch donor and accretor
   b% d_i = 2
   b% a_i = 1
   ! switch pointers
   b% s_donor => b% s2
   b% s_accretor => b% s1
   ! set point mass index
   b% point_mass_i = 1
   ! update separation and eccentricity to be very large and zero respectively
   call binary_set_separation_eccentricity(binary_id, 1d99, 0d0, ierr)
   if (ierr /= 0) return
   ! for sanity checks
   ! print *, '****************************************'
   ! print *, 'new period:'
   ! print *, b% period
   ! print *, 'new separation:'
   ! print *, b% separation
   ! print *, 'new Jorb:'
   ! print *, b% angular_momentum_j
   ! print *, 'new eccentricity:'
   ! print *, b% eccentricity
   ! print *, '****************************************'
end subroutine detach_binary
