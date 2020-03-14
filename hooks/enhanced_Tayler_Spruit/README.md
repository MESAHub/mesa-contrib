# ``enhanced_Tayler_Spruit``

## Overview

These routines implement the enhanced Tayler-Spruit dynamo prescription of [Fuller et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3661F/abstract).
This prescription sets `am_nu_omega`, which is the diffusivity of angular velocity.
In principle the dynamo also results in chemical mixing, but these routines do not implement that effect.
If you use these routines please cite this publication.

This contribution is maintained by @adamjermyn.

## Usage

Example inlists and `run_star_extras.f90` are provided for a normal
MESA `star` run in the `star_example` subfolder. These are the same
as those appearing in the [Zenodo](https://zenodo.org/record/3228403#.XinEzi3MzUI) entry for [Fuller et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3661F/abstract),
but with the relevant subroutines split into `enhanced_Tayler_Spruit.inc`.

To use these routines you need to include them in your `run_star_extras.f90`, after
the `contains` statement you'll need

````Fortran

      contains

      include 'enhanced_Tayler_Spruit/enhanced_Tayler_Spruit.inc'
````

In `extras_controls` point you also need

````Fortran
      subroutine extras_controls(id, ierr)
         ...
         s% other_am_mixing => TSF
         ...
      end subroutine extras_controls
````

Note that `star_example/run_star_extras.f90` includes routines for adding extra history and profile columns which may be of interest.
