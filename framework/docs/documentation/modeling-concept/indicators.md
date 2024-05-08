---
title: Indicator Accounting
lang: en-US
---

# Indicator Accounting

## General concept

The accounting module has a system of abstraction called indicators.
Indicators are the accounting variables in REMix. They are used to associating relative
values to the decision variables in a model. One or more indicators or decision variables
can be mapped to another indicator. The user can decide the extent to which each indicator
or variable should be factored in the other indicator. The desired value of each indicator
can be specified or limited by indicator bounds.

The most fundamental application of
the accounting module is the calculation of the system cost. A model needs at
least one indicator to be optimized, namely the target variable of the
{ref}`objective function <Eq_accounting_objective>`. This has to be a quantity
to be minimized or maximized during a model optimization. Usually, we associated
this indicator to the cost of the system.

Indicators can be associated to each other in accounting networks defined in the
"perIndicator" parameter table. This allows them to influence the decision
variables. The indicators are associated to each one of the core modules by
their respective accounting parameter tables. Each core module has their own
unique form of indicator assignment. You can learn about them in their
respective sections.

## Variable indicators

One feature unique to the accounting module is the possibility of using variable
indicators. These indicators are accounted independently of the other modules
and are useful as backstop variables for other indicators. An example
application are slack variables. Sometimes models are not solvable because of
inherent conditions like not enough supply possible for a given demand. In such
cases the solver will immediately tell that the problem is infeasible and will
not try to solve it. In some cases, it is still interesting to see how such a
system would develop, so we can set a slack variable that will complement the
instability, so the model is still solvable. A practical example of this
application can be seen in the basic example {ref}`slack variables
<slack_variables_label>`.
