# Dissipation Balanced Convective Penetration

## Overview

The functions in `dbconvpen.inc` implement the dissipation balanced
convective penetration scheme for convective boundary mixing described
in [Johnston et al. 2024](https://ui.adsabs.harvard.edu/abs/2024ApJ...964..170J/abstract) based on the 3D simulations of [Anders et al.
2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...926..169A/abstract).

This is appropriate for convective boundary mixing due to convective
penetration in stars with convective cores with masses <40M<sub>o</sub>.

## Usage

Make sure `$MESA_CONTRIB_DIR` is set to the appropriate location and
your `run_star_extras.f90` includes this file, e.g.:

``` Fortran90
    contains
       include 'dbconvpen/dbconvpen.inc'
```
Then make sure to set the following in your inlist:

``` Fortran90
    overshoot_scheme(1) = 'other'
    overshoot_zone_type(1) = 'any'
    overshoot_zone_loc(1) = 'core'
    overshoot_bdy_loc(1) = 'top'
```
You can/should add some exponential overshooting on top of this to
include processes other than convective penetration.

N.B.: while the original implementation of [Johnston et al.
2024](https://ui.adsabs.harvard.edu/abs/2024ApJ...964..170J/abstract)
(available on [zenodo](https://doi.org/10.5281/zenodo.10286503)) used
global variables in the `run_star_extras.f90` to make some values
accessible to multiple routines, these have been replaced using
`s%xtra(1:6)`, `s%lxtra(1)`, and `s%ixtra(1:2)` in this
implementation. Therefore if using this implementation, one should not
use these extra variables elsewhere in your `run_star_extras.f90`.


## Maintainer

This contribution is maintained by Mathieu Renzo (`mathren`).

