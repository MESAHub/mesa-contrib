# Reporting issues

If you run into a problem with one of the hooks in ``mesa-contrib``,
you should open an issue and include any error messages and how to
reproduce the problem.  This usually includes

* your MESA version,
* your ``inlist``s,
* your ``run_star_extras.f90`` and
* anything else you've changed.

Please reduce your input files to the minimum needed to reproduce the
problem.

# Contributing

The whole point of ``mesa-contrib`` is to encourage and promote users'
implementations of useful tools and physical models, so we welcome
contributions!  A contribution might be an implementation of a new
formula for a stellar wind or a routine that writes out data in a
format that can be passed to another stellar physics program.  We
particularly encourage you to add code that you've used in
publications, so that others can easily reproduce your work and cite
it if they use or develop it.

## What to include

The code for your new hook should usually go in a folder like
``hooks/my_new_hook`` unless it's a single file, in which case it can
go straight in ``hooks`` (but must be well documented in the
file).  Your Fortran source should follow the 
[MESA style](https://docs.mesastar.org/en/latest/developing/code_style.html).

You should add a ``README.md`` file that

* describes what the code does
* gives your GitHub username and
* describes how to use the hook.

If appropriate, you should point to any publication (or other citable
object, e.g. a Zenodo entry) that users should cite.  You should also
add your name to the top-level ``CODEOWNERS`` file so that you will be
notified of changes to your contribution.

It's helpful (and users are more likely to use your hook) if you
include an example or examples.  You should design these so that they
will work with something like

```
cp -r "${MESA_DIR}"/star/work test_contrib
cd test_contrib
cp "${MESA_CONTRIB_DIR}"/my-new-hook/star_example/src/* src/
./mk
./rn
```

## How to add your code	to ``mesa-contrib``

New contributions must be proposed through GitHub's *pull request*
(PR) system.  The process is roughly:

1. fork the ``mesa-contrib`` repo (click *Fork* in the top-right of the GitHub interface),
2. clone your fork to you computer,
3. create a new branch for your additions (e.g. ``git checkout -b my-new-hook``),
4. make, commit and push your changes and
5. open a PR.

You'll usually develop your new hook with the latest release version
of MESA and should therefore make your PR against the ``release``
branch

Once the PR is opened, the repo maintainers will review your code 
and make a decision about whether to merge it into the repo, perhaps 
after some changes are made.

## Your responsibilities

After your contribution has been merged, you are committing 
to

* responding to users who come across issues and
* keeping your hook working as MESA develops.

Both might sound daunting but the maintainers and other users or
contributors will help where possible. e.g. with simple compilation
errors.  Many changes to keep up with MESA will be minor updates to
the hook interfaces.  If the ``mesa-contrib`` maintainers can simply
update your contribution to be consistent with changes in the MESA
hooks interface, in a way that makes no change to the mathematics
behind your hook, we will.

If, however, it isn't clear how your code needs to be correctly
updated (e.g. some quantity you use is removed from the star pointer),
we'll expect you to commit to fixing or removing it, though we'll try
our best to help figure out what to do.

Code that remains broken will confuse and frustrate users and
eventually be removed from ``mesa-contrib``.
