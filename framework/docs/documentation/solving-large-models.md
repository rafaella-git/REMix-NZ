---
title: Hints on solving (very) large models
lang: en-US
---

(large_models_label)=

# Hints on solving (very) large models

There is no precise definition of what a (very) large energy system model is.
We have observed issues in solving models with a geographical resolution of 70
regions with around 100 technologies or 110 regions and 44 technologies with a
temporal resolution of hours or even less.
You might, however, encounter problems that are specific to large model
instances already earlier.
For all these cases, the following hints are thought to give some ideas on how
to possibly address the size problem before thinking of aggregation or making
some scopes smaller.

## The solver

It has shown that a mere change of solvers can affect the solvability of a model
a lot.
In one instance, e.g., for a large model (100 technologies, 70 regions) that
could not be solved even after seven days by one solver, solved with another
after six and a half hours.

## Threads

The default number of threads used in REMix is four.
Reasonable values until which performance increases seem to lie between twelve
and sixteen threads.
You can just give, e.g., `threads=16` as an argument to `remix run ...`.

## Permutation methods / "barorder"

With the solver argument "barorder", you can control which method is used to do
matrix permutations when solving the model.
REMix default is "Nested Dissection (ND)", which seems to perform best for our
models.
There are other options, e.g., for the CPLEX solver, "Approximate Minimum Degree
(AMD)" or "Approximate Minimum Fill" (AMF) or also "automatic" (cf.
[CPLEX barorder options](https://www.gams.com/latest/docs/S_CPLEX.html#CPLEXbarorder)).

```{note}
The use of the option "automatic" has lead to some performance problems for
some model sizes before, e.g. around 100 minutes for nested-dissection ordering
when using `barorder=0` ("automatic") for CPLEX compared to around 2 minutes
when using nested-dissection ordering directly with `barorder=3`.
```

## Dense columns / "barcolnz"

The solvers usually determine dense columns dynamically (e.g.
[CPLEX barcolnz](https://www.gams.com/latest/docs/S_CPLEX.html#CPLEXbarcolnz)).
This can be controlled by the user via the argument "barcolnz".
If `barcolnz=0` (default), the solver decides automatically how many entries a
column needs to have to be considered "dense".
This has shown to have influence on the solver's performance.
If you find the number of dense columns (given in the output file) suspiciously
small, you might want to try out this argument.
Reasonable values for "barcolnz" for our kinds of problems seem to lie between
60 and 80.

## Floating point operations to factor

The floating point operations needed to calculate a factor in one iteration is
given as an estimation by the solvers.
With CPLEX, e.g., this is called "Total FP ops to factor", with Gurobi it is
"Factor Ops". These can give an estimate about how long your model will need to
be solved.
Usually, the solution of a model takes 80 to 180 iterations, depending on model
size.
