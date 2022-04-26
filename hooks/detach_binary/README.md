# ``detach binary``

## Overview

This routine can be used to continue the evolution of the initially
less massive star (typically `star 2`) as an "effectively single"
star. This can be used, for example, to simulate the disruption of a
binary because of the core-collapse of the initially more massive star
([Blaauw 1961](https://ui.adsabs.harvard.edu/abs/1961BAN....15..265B/abstract),
[Renzo et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A..66R/abstract)).

When calling this routine:

 * `star 1` will be turned into a point mass, the donor pointer `b% d_i%` is
   set to `star 2`
 * the binary separation is set to something enormous (`1d99` solar
 radii)
 * tides, Roche lobe overflow, magnetic braking, and angular momentum accretion will
   be ignored in the remaining evolution

It assumes that you are running a binary with `evolve_both_stars =
.true.`, and otherwise will set `ierr=1` and hit a Fortran `STOP`. You
can comment out the `STOP` statement, exit the subroutine, and check
the returned `ierr` from where you called it

## Usage

Typically, this will be used from `extras_binary_finish_step` in
`run_binary_extras.f90`.

To use this routine you need to include it in your `run_binary_extras.f90`,
after the `contains` statement you'll need

```Fortran

      contains

        include 'detach_binary/detach_binary.inc'
```

### Example

Example inlists and `run_binary_extras.f90` (based on the
`$MESA_DIR/binary/work` template) are provided in `test_detach_binary`. 

In this example, at `model_number = 3` the binary 
is detached and the evolution of the initially 
less massive star is continued until `model_number=5`.

It has been tested with MESA version 15140 and `29bd15f`.

## Maintainer

This contribution is maintained by Mathieu Renzo (`mathren`)
