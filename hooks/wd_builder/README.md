# ``wd_builder``

## Overview

This provides an alternate initial model builder that creates a white
dwarf of a given mass and chemical composition.  The initial thermal
structure is roughly approximated as a degenerate isothermal core and
a radiative envelope.  The initial model is not in thermal equilibrium
and so there is an initial transient phase lasting a few thermal times
that should be disregarded.

This contribution is maintained by @jschwab.

## Usage

Example inlists and `run_star_extras.f90` are provided for a normal
MESA `star` run in the `star_example` subfolder.

To use these routines you need to include them in your
`run_star_extras.f90`, after the `contains` statement you'll need

````Fortran
      include 'wd_builder/wd_builder_def.inc'

      contains

      include 'wd_builder/wd_builder.inc'
````

In `extras_controls` you also need

````Fortran
      subroutine extras_controls(id, ierr)
        ...
        s% other_build_initial_model => wd_builder
        ...
      end subroutine extras_controls

```` 

and to set `use_other_build_initial_model = .true.` in your inlist.

The user is also responsible for providing an implementation for the
subroutine `wd_builder_get_xa` that returns the chemical composition
as a function of mass coordinate.  Two example implementations are
provided: one that reads from a file and one that uses a simple
parameterization.
