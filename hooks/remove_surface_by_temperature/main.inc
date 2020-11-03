      subroutine do_remove_surface_by_temperature(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_temperature: get_star_ptr ierr', ierr
            return
         end if
         if (remove_surface_temperature <= s% T(1)) return
         do k=1,s% nz
            if (s% T(k) >= remove_surface_temperature) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_temperature
