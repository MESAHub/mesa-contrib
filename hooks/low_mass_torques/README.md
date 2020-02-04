# ``low_mass_torques``

## Overview

These routines implement the wind angular momentum loss schemes of [Kawaler 1988](https://ui.adsabs.harvard.edu/abs/1988ApJ...333..236K/abstract) and
[van Sanders & Pinsonneault](https://ui.adsabs.harvard.edu/abs/2013ApJ...776...67V/abstract) using the ``other_torque`` hook.
The implementation is based on the MESA 2019 summer school lab by [Marc Pinsonneault](https://zenodo.org/record/3374959#.XjmqPC2ZPUI).
These are most relevant to low-mass stars, and notably assume that the star has a surface convection zone.
If you use these routines please cite the relevant publications, as well as the summer school lab with ``doi:10.5281/zenodo.3374959``.

This contribution is maintained by @adamjermyn.

## Usage

Example inlists and `run_star_extras.f` are provided for a normal
MESA `star` run in the `star_example` subfolder.

To use these routines you need to include them in your `run_star_extras.f`, after
the `contains` statement you'll need

````Fortran

      contains

      include 'low_mass_torques/low_mass_torques.inc'
````

In `extras_controls` point you also need

````Fortran
      subroutine extras_controls(id, ierr)
        ...
        s% other_torque => low_mass_torques
        ...
      end subroutine extras_controls

      subroutine extras_startup
        ...
      	include 'low_mass_torques/low_mass_torques_extras_startup.inc'
      end subroutine extras_startup
````

Note that `star_example/run_star_extras.f` includes routines for adding extra history and profile columns which may be of interest.
