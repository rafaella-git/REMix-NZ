---
title: What can I contribute?
lang: en-US
---

(contributing_what)=

# What can I contribute?

## Feedback and bug reports

If you encounter a bug in the software, we appreciate if you open an
[issue](https://gitlab.com/dlr-ve/remix/framework/-/issues). Please
include the following information:

-  Your operating system name and version.
-  Your REMix and GAMS version.
-  A minimal example showcasing the bug or steps to reproduce the bug.

## Request new features

If you have an idea for a new feature, you can open an issue request as
well. Please

- describe what the feature aims to do and how it should work.
- try to narrow down the feature as much as possible, to make implementation
  easier.

## Actively develop REMix

You can easily participate in the active development of REMix, for
example by providing your developments to others and allow them to use
and build upon your work. For this, you need to merge the changes you
made into the framework. To make sure no accidental changes are made
this requires a so called merge request where you also get some feedback
on the code you created. There are three different entry points for your
active development:

-  documentation.
-  pre- and post-processing.
-  REMix core.

We will walk you through changing your REMix set up to make active
development possible in {ref}`this section <contributing_active_dev_label>`.
The structure of the project is as follows:

### Documentation

Improving the documentation is a very low level entry point to contribution. If
you spot any mistakes or typos in the documentation, we invite you to create a
merge request. The online documentation can be found in the `docs` subdirectory
of `framework`.

### Pre- or post-processing tools (Python)

The pre- and post-processing of any model instance is done using Python.
In the framework folder structure, the `remix/framework/` subdirectory contains
all code required for pre- and post-processing the model instance and to start
an optimization run.

### Optimization model (GAMS)

The core of REMix is written in GAMS and can be found in the same place
as the Python code in the `model` module.
