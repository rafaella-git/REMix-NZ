---
title: How to debug your REMix model
lang: en-US
---

(debugging_label)=

# How to debug your REMix model

Once you start setting up your own REMix model, you might encounter problems
writing your model to files using the
{class}`remix.framework.api.instance.Instance` class or during the solution
process with GAMS.
This part of the documentation is aiming to provide you with some hints on
common errors or infeasibilities and probable ways to fix them.
In case you do not find a solution that works for you here, you can also resort
to writing an issue or ask in the
[REMix section of the openmod forum](https://forum.openmod.org/tag/remix).

## Common model infeasibilities with REMix/GAMS

If you encounter an infeasibility during the solution process of your model, you
need to understand where it comes from and why.
This can be fairly obvious at times, at others pretty hard.
How do you know you have run into an infeasibility?
GAMS will tell you.
Even if you might oversee that, you sometimes might find that there was no new
`remix.gms` file created.
Most likely, an infeasibility is the reason for that.

In this section, we firstly introduce you to a debugging workflow
that has proven useful for us in the past. This is followed by a collection on
the most common errors and infeasibilities and how to fix them.

### Debugging your model

#### A possible workflow

When you encounter an error or infeasibility:

- 1. open the *.lst file (e.g. `run_remix.lst`) and look for details on the
error
- 2. errors are accompanied by for asterisks (`****`), so use `STRG+F` and
search for their first appearance
- 3. find the GAMS error code

If the GAMS error code is not obvious (happens...), head over to the next
sections.

#### Useful GAMS/remix.run arguments for debugging

When calling {class}`remix.framework.cli.run`, the following arguments might be
helpful for pinpointing the actual issue.

- i. `names=1` : to get understandable GAMS error/infeasibility messages
- ii. `lo=4` : writes log of REMix run in a dedicated `*.log` file (instead of
terminal)
- iii. `keep=1` : keeps temporary folder with `*.csv` files built during the
model run, might be interesting for debugging purposes; REMix default: 0
- iv. `gdx="default"` : writes out an additional `*.gdx` file containing all
symbols
- v. `iis=1` : "irreducible infeasible subsystem", writes out minimal
combination of variables in conflict

### Common issues and how to fix them aka REMix debugging FAQ

The following is a collection of more or less common issues that may occur when
creating and solving a REMix model.

#### Common mistakes using the REMix.Instance class

- years are sometimes given as integers, sometimes as strings

#### Numerical difficulties

##### issue / problem indication in `*.lst` file:

"Solution available but not proven optimal due to numerical difficulties"

##### how to get more information on the issue:

Have a look at the "Range statistics".
These should be in a range of 1e+4 and 1e-4 for numerical stability.
If you find that your LP lies outside these boundaries, this will probably be
your problem.
Turn on the GAMS equation listing by specifying `equlist=1` in
{class}`remix.framework.cli.run`.
Find the highest and/or lowest values in your model (e.g., by searching for
`E-7` in the `*.lst` file) to get an idea, where the issue may arise.

##### possible approaches

- change your output units (for large-scale problems consider using "GW(h)" as
power/energy unit, "MEUR" as cost unit and "Mt" as mass unit)
- round values before running the model (e.g. to four digits)
- have a closer look at your renewable potentials (=activity profiles); they might lead to the problem combined with small capacities
