#!/usr/bin/env python3
#
# Generates the command-line argument hook for MESA, which is saved to
# `$MESA_CONTRIB_DIR/hooks/cmd_line_args.inc`.
#
# You can add or remove the parameters you'd like to control from the
# list `args`, below.
#
# To use the command line arguments in MESA, add the variable
# declarations
#
#       integer :: i_arg
#       character(len=32) :: arg
#
# to the preamble of the `extras_controls` function and
#
#       include `cmd_line_args.inc`
#
# after the call to `star_ptr` in your `run_star_extras.f90`.  In your
# work folder, run `./clean` and `./mk` and then run MESA with, e.g.
#
#     $ ./star inlist --initial-mass 2.0
#
# Added by @warrickball following the simple example by @jschwab on
# the mailing list.
#
# https://lists.mesastar.org/pipermail/mesa-users/2020-January/011068.html

import os

# Change this list to include the controls you'd like to vary from the
# command line.
args = ['initial_mass', 'initial_z', 'initial_y', 'mixing_length_alpha']

with open(os.environ['MESA_CONTRIB_DIR'] + '/hooks/cmd_line_args.inc', 'wt') as f:
    f.write('\n'.join([
        "         ! Add these two variables to extras_controls preamble",
        "         ! integer :: i_arg",
        "         ! character(len=32) :: arg",
        "         ! then include this after the star pointer is set.",
        "",
        "         i_arg = 2 ! first argument is filename of inlist",
        "         call GET_COMMAND_ARGUMENT(i_arg, arg)",
        "         do while (arg /= ' ')",
        "            select case (arg)",
        ""]))

    for arg in args:
        f.write('\n'.join([
            "            case ('--%s')" % arg.replace('_', '-'),
            "               i_arg = i_arg + 1",
            "               call GET_COMMAND_ARGUMENT(i_arg, arg)",
            "               read(arg, *) s%% %s" % arg,
            ""]))

    f.write('\n'.join([
        "            case default",
        "               stop 'invalid command-line argument: '//trim(arg)",
        "            end select",
        "",
        "            i_arg = i_arg + 1",
        "            call GET_COMMAND_ARGUMENT(i_arg, arg)",
        "         end do",
        ""]))
