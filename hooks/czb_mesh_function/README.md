# ``czb_mesh_function``

## Overview

This provides a mesh function that can be used to increase resolution
in the vicinity of convection zone boundaries (czb).

This contribution is maintained by Josiah Schwab (@jschwab).

## Usage

An example inlist and `run_star_extras.f90` are provided in the
`star_example` subfolder.

The czb_mesh inlist options follow the same pattern-matching approach
as the MESA overshooting controls.

```Fortran
  ! set defaults
  czb_mesh_zone_type(:) = '' !  Possible values: burn_H, burn_He, burn_Z, nonburn, any
  czb_mesh_zone_loc(:) = ''  !  Possible values: core, shell, any
  czb_mesh_bdy_loc(:) = ''   !  Possible values: bottom, top, any

  czb_mesh_f(:) = 0d0 ! target zone spacing in delta_lnP (i.e. ~= 1/f zones per scale height)
  czb_mesh_l(:) = 0d0 ! extent of enhanced resolution region (in pressure scale heights)

```

