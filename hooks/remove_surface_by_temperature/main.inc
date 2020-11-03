      subroutine remove_surface_by_temperature(id, ierr, k)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer, intent(out) :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'remove_surface_by_temperature: get_star_ptr ierr', ierr
            return
         end if
         if (remove_surface_temperature <= s% T(1)) return
         do k=1,s% nz
            if (s% T(k) >= remove_surface_temperature) then
               return
            end if
         end do
         ierr = -1
      end subroutine remove_surface_by_temperature
