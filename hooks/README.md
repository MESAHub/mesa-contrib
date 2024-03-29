# Style suggestions

## Overview

Routines and data can go anywhere in a contribution directory, including in subfolders.
You just need to be sure to point to the correct locations in your `run_star_extras.f90`.

Please include
``
This contribution is maintained by Your Name/GitHub ID
``

in your `README.md` so that your contribution has a contact and known maintainer.
You should also add your name to the top-level `CODEOWNERS` file.

## Usage

It helps to include example directories showing how your contribution is used.
Depending on which parts of MESA your contribution touches you may want to include examples for the star, rsp, and/or astero modules.

Typically an example directory contains whatever inlists are needed as well as a `src` directory containing the sample `run_star_extras.f90`.
In your `run_star_extras.f90`, around the `contains` statement you'll need

````f90
      include 'contribution_directory/any_data_to_include.inc'

      contains

      include 'contribution_directory/routines.dek'
      include 'contribution_directory/more_routines/inc'
````

If you have chosen to put your routines and data in subfolders, just amend the paths above accordingly.

Also, if your contribution makes use of any MESA hooks you'll need to point those hooks in extras_controls:

````f90
      subroutine extras_controls(id, ierr)
         ...
         call hydro_Ttau_setup(id, ierr)
         ...
         s% other_gradr_factor => hydro_Ttau_gradr_factor
         s% other_surface_PT => hydro_Ttau_surface_PT
         ...
      end subroutine extras_controls
````

Please document all such changes users need to make in order to use your contribution.

## Extra history and profile columns

Sometimes users will want to use multiple contributions in their work.
Each contribution might have routines for adding history and/or profile columns.
To make it easier for users to combine these routines, we recommend that these routines have a hook-specific name.
That avoids namespace conflicts.
Otherwise we recommend that these routines follow the call signature of the current methods for extra profile and history columns. 
If your routines follow this signature, then a user can incorporate your routines alongside other contributions as follows:

````f90
      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
         how_many_extra_history_columns = how_many_extra_history_columns + contribution1_how_many_extra_history_columns(id)
         how_many_extra_history_columns = how_many_extra_history_columns + contribution2_how_many_extra_history_columns(id)
         how_many_extra_history_columns = how_many_extra_history_columns + contribution3_how_many_extra_history_columns(id)
         ...

      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: j
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         
         j = 1
         nc = contribution1_how_many_extra_history_columns(id)
         call contribution1_data_for_extra_history_columns(id, nc, names(j:j+nc), vals(j:j+nc), ierr)

         j = j + nc
         nc = contribution2_how_many_extra_history_columns(id)
         call contribution2_data_for_extra_history_columns(id, nc, names(j:j+nc), vals(j:j+nc), ierr)

         j = j + nc
         nc = contribution3_how_many_extra_history_columns(id)
         call contribution3_data_for_extra_history_columns(id, nc, names(j:j+nc), vals(j:j+nc), ierr)
	...
	
      end subroutine data_for_extra_history_columns
````

and similarly for profile columns. Note that for profile columns you need to explicitly slice the first index of the ``vals`` array, as in ``vals(1:s%nz,j:j+nc)``.
