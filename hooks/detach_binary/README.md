# ``detach binary``

## Overview

This routine can be used to simulate the breaking of a binary (e.g.,
because of a SN explosion). It assumes that you are running a binary
with `evolve_both_stars = .true.` and will stop the run otherwise when
called.

## Usage

Typically this will be used from `extras_binary_finish_step` in
`run_binary_extras.f90`.

Example inlists and `run_binary_extras.f90` (based on the
`$MESA_DIR/binary/work` template) are provided in `test_detach_binary`.
In this example at `model_number = 3` the binary is detached and 
the evolution of the initially less massive star is continued until
`model_number=10`. It has been tested with MESA revision 15140.

To use this routine you need to include it in your `run_binary_extras.f90`,
after the `contains` statement you'll need

```Fortran

      contains

        include 'detach_binary/detach_binary.inc'
```

## Maintainer

This contribution is maintained by Mathieu Renzo (`mathren`)
