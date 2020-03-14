# MESA Contributed Procedures

This repo contains contributed procedures that can be used with the
stellar evolution code [Modules for Experiments in Stellar Astrophysics
(MESA)](http://mesa.sourceforge.net/).

Although this repo is administered by the MESA team and many of the
contributions are from team members, it is *not* part of the main
trunk and the MESA team takes no responsibility for the correctness of
the procedures.  If you are having trouble with a particular
procedure, you should contact the person(s) who contributed that
procedure to see if they can help.

Many of the procedures will have been developed for a particular
scientific purpose and have an associated publication in which the
contributor used that procedure.  If so, please cite the relevant
publication(s) when using that procedure.

## Usage

The `contrib` repo can be used with the current development version of
MESA (support was added in r12477-9) by checking out this repo to some
directory (e.g. `$HOME/mesa-contrib`) and setting the environment
variable `MESA_CONTRIB_DIR` to point there. e.g.

    export MESA_CONTRIB_DIR=$HOME/mesa-contrib

Then, `$MESA_CONTRIB_DIR/hooks` will be added to the include path when
compiling a MESA `work` folder.  Use your chosen routine by adding an
`include` statements after `contains` in your `run_star_extras.f`. e.g.
to use `hooks/other_wind_grafener.f90`, you would add something like

    ...
    contains
   
    include 'other_wind_grafener.f90'       ! add this line
   
    subroutine extras_controls(id, ierr)
       integer, intent(in) :: id
    ...
