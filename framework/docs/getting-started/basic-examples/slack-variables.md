---
title: "Slack variables"
lang: en-US
---

(slack_variables_label)=

# Slack variables

It is common to have mathematical infeasibilities when building a model.
A way to make it solve although there are infeasibilities is the
introduction of slack variables into a model. We will show how that can be done
in this introductory tutorial.

In this example we will show how to use sources to identify when some constraint
is too tight for a system configuration.
The REMix components used in this example are one
{ref}`sink <remix_model_core_sourcesink_label>`, one
{ref}`source <remix_model_core_sourcesink_label>` and a single
{ref}`converter <remix_model_core_converter_label>`.

We set up a single-node cost-minimizing model.
Since we won't allow capacity expansion, there is initially no variables
contributing to the system cost.

```python
>>> from remix.framework.api.instance import Instance
>>> from numpy import Inf
>>> m = Instance()
>>> m.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> m.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

```

Our model consists of a converter with a fixed activity profile.
It has also a fixed 1 MW capacity and a generation coefficient of 1 MWh per
activity step.

```python
>>> m.profile.converter_activityprofile.loc[("R1_data", "2030", "WindOnshore", "upper"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
>>> m.profile.converter_activityprofile.fillna(0, inplace=True)
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "WindOnshore"), ("unitsLowerLimit", "unitsUpperLimit")] = (1, 1)
>>> m.parameter.converter_coefficient.loc[("WindOnshore", "2030", "Powergeneration", "Electricity"), ("coefficient")] = 1
>>> m.parameter.converter_techparam.loc[("WindOnshore", "2030"), "lifeTime"] = 10

```

The demand load is also fixed but with an overshot at the second time step.
Under these conditions the model is infeasible as the generation cannot satisfy
the load.
This means that the solver will not try to solve the problem.

```python
>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"],] = [-n for n in [0.1, 0.2, 0.1, 0.1, 0.1, 0.1]]
>>> m.profile.sourcesink_profile.fillna(0, inplace=True)
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Electricity"), ("usesFixedProfile")] = 1

```

Sometimes you can be interested in the results of such an infeasible model,
although it is infeasible.
For those cases we suggest including slack variables.
Using a source with infinite upper and 0 lower limits we can make the model
feasible.
In order to avoid slack from interfering with actual optimal solutions we assign
it a cost of significantly higher order of magnitude compared to the rest of the
cost variables.

```python
>>> m.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "Slack", "Electricity"), ("lower", "upper")] = (0, Inf)
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Slack", "Electricity"), ("usesUpperSum", "usesLowerSum")] = (1,1)
>>> m.parameter.accounting_sourcesinkflow.loc[("SlackPenalty", "R1_data", "2030", "Slack", "Electricity"), "perFlow"] = 1
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "SlackPenalty", "global", "2030"), "perIndicator"] = 10000

```

It should be possible to solve this model now.

```python
>>> m.set.add(["2030"], "yearssel")
>>> m.datadir = "slack_model"
>>> m.write(fileformat="dat")
>>> ();m.run(timeend=6, read_result=True);()
(...)

```

If we take a look at the electricity produced by the slack, we see that the
value is consistent with the input difference.

```python
>>> m.result["commodity_balance_annual"].loc[("global", "2030", "Slack", "Electricity", "positive"), "value"]
0.1

```