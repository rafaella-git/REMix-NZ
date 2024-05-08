---
title: Introduction to REMix tutorials
lang: en-US
---

(tutorials_introduction_label)=

# REMix tutorials

This page gives a general overview of the structure of the REMix tutorials.
The REMix tutorials are Python scripts that are being translated to iPython
notebooks automatically using the `jupytext` package.

## General structure of the tutorials

The tutorials are separated into the categories "basic" and "advanced". The
numbering of the basic tutorials starts with a `1`, the advanced tutorials
with a `2`.

There is a figure for each tutorial that visualizes the energy system that
is being modeled within it.

Each REMix tutorial is split up in four different parts to simplify the learning
experience. Each of these parts is a logical unit that stands and works for itself.
The structure is as follows:

- part a: script to **build** the actual REMix data model
- part b: script to **run** the REMix data model
- part c: script with an example on how to **evaluate** REMix model results
- part d: **bonus** tasks for deeper model understanding

These parts are explained in more detail in the following paragraphs.

### Building the model (part a)

The model scope describes the fundamental dimensions (spatial and temporal) of
the model, e.g. which distinct regions and years are modeled.

It might be good to know that although REMix provides a lot of features to
build a sophisticated energy system model, it will also run if you do not
provide some features. Some features are needed for every model,
however (like setting an objective function).

In this and all other tutorials, we will use the
{class}`remix.framework.api.instance.Instance` class to store our model
parametrization in. As a result, you can find all parameters we set in that
object.

At the end of part a of the tutorials, the content of `Instance` is written out
to GAMS-readable files (either `*.csv` or `*.dat`) into the`./data` directory.

### Running the optimization (part b)

Part b describes in detail how a REMix model is run. The optimization is executed
by calling the {meth}`remix.framework.api.instance.Instance.run` method of
`Instance`. A detailed explanation is given on which (common) command line
arguments can be used to customize the optimization output.

### Evaluating model results (part c)

The output of each REMix optimization run is a `*.gdx` file. In part c, we look
at how to read that into Pandas DataFrames using the Python package `gamstransfer`
and plot some exemplary results.

### Bonus tasks (part d)

The bonus tasks try to give some extended knowledge on the model building in REMix.
In tutorial 101, for example, an introduction into error handling in REMix/GAMS is
given. Others compare the influence of implementing different technologies on the
overall optimization results.

## Tutorials

```{toctree}
:glob:
:maxdepth: 2

basic
advanced
```
