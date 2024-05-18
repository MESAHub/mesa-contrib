# Fuller et al. angular momentum transport

## Overview

The functions in `Fuller_AM_transport.inc` implement the modified
Tayler-Spruit dynamo by Fuller and collaborators.

Specifically `TSF_Fuller19` implements the algorithm described in
[Fuller et al.
2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3661F/abstract)
and available in the original form at
https://zenodo.org/records/3228403, and `TSF_Fuller_Lu22` implements
the algorithm described in [Fuller & Lu
2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.3951F/abstract)
and available in the original form at
https://zenodo.org/records/5778001. Here, they are only adapted for
compatibility with the current MESA version (r24.03.1). All credits go
to the original authors of the papers above, which you should cite if you
use this.

The two functions are overall similar, but the latter does not have a
minimum shear threshold to activate and has a different efficiency
value. Refer to the above papers and references therein for more
details.

## Usage

The two routines are mutually exclusive, you can only use one _or_ the
other.

Make sure `$MESA_CONTRIB_DIR` is set to the appropriate location and
your `run_star_extras.f90` includes this file, e.g.:

```Fortran

      contains

        include 'Fuller_AM/Fuller_AM_transport.inc'
```

Then, these can be used as any `other` hook in MESA, that is:

* set `use_other_am_mixing=.true.` in your inlist
* inside `run_star_extras.f90`, in the routine `extras_controls`,
  point `other_am_mixing` to one of the two functions, e.g.:
  `s% other_am_mixing => TSF_Fuller19` or `s% other_am_mixing => TSF_Fuller_lu22`

### Example

Example inlists and `run_star_extras.f90` (based on the
`$MESA_DIR/star/work` template) are provided in `Fuller_AM_example`
and by default will call the Fuller & Lu 2022 version.

It has been tested with MESA version r24.03.1.

## Maintainer

This contribution is maintained by Mathieu Renzo (`mathren`)
