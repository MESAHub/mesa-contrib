  integer, parameter :: i_Tc = 1
  logical, parameter :: dbg = .false.

  abstract interface
     subroutine get_xa_interface(s, q, xa)
       use const_def, only: dp
       use star_def, only: star_info
       type (star_info), pointer :: s
       real(dp), intent(in) :: q
       real(dp) :: xa(:)
     end subroutine get_xa_interface
  end interface

  procedure (get_xa_interface), pointer :: wd_builder_get_xa
