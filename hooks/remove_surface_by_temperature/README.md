
# ``remove_surface_by_temperature``

## Overview

These routines enable MESA to remove surface material below a fixed temperature.
This was previously implemented as part of MESA/star/private/remove_shells.f90 in the routine do_remove_surface_by_temperature.

This contribution is maintained by @adamjermyn.

## Usage

Example inlists and `run_star_extras.f90` are provided for a normal
MESA `star` run in the `star_example` subfolder.

To use these routines you need to include them in your `run_star_extras.f90`, after
the `contains` statement you'll need

````Fortran

      real(dp), parameter :: remove_surface_temperature = (desired temperature)

      contains
 
      include 'remove_surface_by_temperature/main.inc'
````

In `extras_controls` point you also need

````Fortran
      subroutine extras_controls(id, ierr)
        ...
         s% other_remove_surface => remove_surface_by_temperature
        ...
      end subroutine extras_controls
````

