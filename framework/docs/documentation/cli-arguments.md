---
title: Command-Line Interface
lang: en-US
---

(cli_label)=

# Command-Line Interface (CLI)

It is possible to run REMix models directly from the data files with a
command-line interface (CLI) without the Python API with the command `remix run`.
Most frequently required **REMix arguments** are listed in the table below,
for **GAMS core** arguments refer to the link at the end of this page.
Boolean values are to be provided with `1` (`True`) and `0` (`False`).

```{note}
In case you are using the Python API calling the
{meth}`remix.framework.api.run.run_remix` function, you can provide the same
keywords as provided in the list below without the prefixing dashes, e.g.
`run_remix(datadir="mydata")`.
```

```bash
remix run <argument>
```

## REMix-specific arguments

| argument | meaning |
|:-------- |:------- |
| `--help` | Overview on optional keywords for CLI |
| `--datadir=...` | Base data directory |
| `--scendir=...` | Path to the scenario to calculate, relative to the `--datadir` path |
| `--resultdir=...` | Directory to write the results to |
| `--resultfile=...` | Name of the result file |
| `--timestart=...` | First time step for the optimization, default 1 |
| `--timeend=...` | Last time step for the optimization, default 8760 |
| `--timeres=...` | "time resolution", e.g. timeres=24 aggregates hourly data to daily resolution, default 1 |
| `--postcalc=...` | Run post-processing of the results, default 1 |
| `--roundts=...` | Automatically round after-comma digits in large profile files where necessary to successfully read them |
| `--roundcoefs=...` | Automatically round profiles for converter activities, converter coefficients and source sink profiles to 1e-3, set all values below to 0, default 0 |
| `--equlist=...` | Write all equations with all variables in the output remix.lst file, can be helpful for model debugging, should only be used if very few time steps, otherwise GB-sized *.txt file, default 0 |
| `--method=...` | Choose method: "solve" = run the model with commercial solvers, "pareto" = solve the REMix model once and afterwards run multiple points along the pareto front, "mga" = solve the model once and afterwards maximize the distance metric of indicators, "iternodes" = run the model by iterating through all individual nodes while all connections between nodes are ignored, "pips" = build an annotated *.gdx file for PIPS-IPM++ |
| `--gdx_mipconverter=...` | Write the MIP features for converters into the results file |
| `--gdx_mipstorage=...` | Write the MIP features for storage technologies into the results file |

## Important GAMS arguments

| argument | meaning |
|:-------- |:------- |
| `gdx=...` | if argument given, e.g. gdx="all_symbols", will give out an extra *.gdx file with all existing symbols in model, not only those specified in postprocessing |
| `keep=...` | Keep the scratch directory (e.g. `/225a`) containing files used for the optimization run (usually deleted after a model run), default 0 |
| `lo=...` | 'log option' of GAMS; ensures that the output from GAMS will be visible in the terminal (`lo=3`; `lo=4` will write out a log file) |
| `profile=...` | activates GAMS profiling to show how long an assignment takes, gives out 10 most computing-expensive operations at the end of the *.lst file |

## Important solver arguments

| argument | meaning |
|:-------- |:------- |
| `--crossover=...` | Run the optimization with crossover (1) or without crossover (0, default) |
| `--datacheck=...` | datacheck=2 will make CPLEX give out warnings about disproportionate values (which lead to numerical difficulties = non-optimal solutions) at the beginning of the optimization, default 0 |
| `--iis=...` | *.gdx file with all infeasible equations (if any), default 0 |
| `--names=...` | Write out actual variable names into *.lst file in case of error, default 0 |

```{note}
Full list of available command-line arguments for

- [GAMS core parameters](https://www.gams.com/latest/docs/UG_GamsCall.html#UG_GamsCall_ListOfCommandLineParameters)
- [CPLEX solver](https://www.gams.com/latest/docs/S_CPLEX.html)
- [XPRESS solver](https://www.gams.com/latest/docs/S_XPRESS.html)
```
