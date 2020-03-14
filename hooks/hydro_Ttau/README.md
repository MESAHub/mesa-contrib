# ``hydro_Ttau``

## Overview

These routines implement atmospheric T(τ) relations and mixing-length
parameters derived from and calibrated to 3D radiation hydrodynamics
simulations of the near-surface layers of cool stars.  The T(τ)
relations are described by [Trampedach et
al. (2014a)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.442..805T/abstract)
and the calibrated MLT parameters by [Trampedach et
al. (2014b)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.445.4366T/abstract).
Their implementation in MESA is detailed by [Mosumgaard et
al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.5650M/abstract).
Two examples using the T(τ) relations were provided through the test
suite from SVN r10460 and subsequently moved into the `mesa-contrib`
repo.  If you use the `hydro_Ttau` routines, please cite these
publications.

This contribution is maintained by [@warrickball](https://github.com/warrickball).

## Usage

Example inlists and `run_star_extras.f90` are provided for both a normal
MESA `star` run and run with the `astero` submodule in the
corresponding subfolders.  These are in essence the same inlists that
were used the old test cases `hydro_Ttau_evolve` and
`hydro_Ttau_solar`, which are in turn what was used for the results in
[Mosumgaard et
al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.5650M/abstract).
Those results were computed using MESA r9575 so some controls and
parameters will have changed as MESA developed.

To use the `hydro_Ttau` relations, you need to include the definitions
in `hydro_Ttau_def.inc` and procedures in `hydro_Ttau_procs.inc` and
`624.dek` in the usual way. i.e. in your `run_star_extras.f90`, around
the `contains` statement you'll need

````Fortran
      include 'hydro_Ttau/hydro_Ttau_def.inc'

      contains

      include 'hydro_Ttau/624.dek'
      include 'hydro_Ttau/hydro_Ttau_proc.inc'
````

In `extras_controls`, call the subroutine `hydro_Ttau_setup`, which
will load the table of T(τ) relations and mixing-length parameters and
make a first call to the interpolator to put initial data in the
relevant arrays.

You'll also need to point `s% other_gradr_factor` and `s%
other_surface_PT` to the matching `hydro_Ttau` procedures, named
`hydro_Ttau_gradr_factor` and `hydro_Ttau_surface_PT`,
respectively. i.e. your `extras_controls` subroutines should look like

````Fortran
      subroutine extras_controls(id, ierr)
         ...
         call hydro_Ttau_setup(id, ierr)
         ...
         s% other_gradr_factor => hydro_Ttau_gradr_factor
         s% other_surface_PT => hydro_Ttau_surface_PT
         ...
      end subroutine extras_controls
````

Finally, you need to call `hydro_Ttau_update` to interpolate in `logg`
and `Teff` for each stellar model, which means you need to add a call like
````Fortran
         call hydro_Ttau_update(id, ierr)
````
to `extras_start_step` and `extras_check_model`.

Two quantities that `hydro_Ttau` manipulates that are worth following
in your histories are `s% mixing_length_alpha` and `s% tau_base`.  The
examples add these using the `extra_history_columns` procedures.
