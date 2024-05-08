---
title: "Secured capacity"
lang: en-US
---

(secured_capacity_label)=

# Secured capacity

Secured capacity is defined in models to ensure that the
built capacity of a system can always satisfy the peak demand {cite}`Sasanpour2021`.

The REMix components used in this example are two
{ref}`converters <remix_model_core_converter_label>` and a single
{ref}`sink <remix_model_core_sourcesink_label>`.

For simplicity, we will only consider costs for operation and maintenance.

```python
>>> from remix.framework.api.instance import Instance
>>> from numpy import Inf
>>> n = Instance()
>>> n.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> n.parameter.accounting_perindicator.loc[("SystemCost", "OMFix", "global", "2030"), "perIndicator"] = 1
>>> n.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

```

One part of the model is a variable renewable energy generator.

```python
>>> n.profile.converter_activityprofile.loc[("R1_data", "2030", "VariableGen", "upper"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0.1, 0.2, 0.9, 0.2, 0.1, 0.2]
>>> n.profile.converter_activityprofile.fillna(0, inplace=True)
>>> n.parameter.converter_capacityparam.loc[("R1_data", "2030", "VariableGen"), ("unitsLowerLimit", "unitsUpperLimit")] = (0, 10)
>>> n.parameter.converter_coefficient.loc[("VariableGen", "2030", "Powergeneration", "Electricity"), ("coefficient")] = 1
>>> n.parameter.accounting_converterunits.loc[("OMFix", "R1_data", "VariableGen", "2030"), "perUnitBuild"] = 5
>>> n.parameter.converter_techparam.loc[("VariableGen", "2030"), "lifeTime"] = 40

```

We also include a generator that is always available, i.e. dispatchable.
Usually this kind of generator has costs associated to its production like
the ones coming from supplying fuel.
We simplify any extra costs by making the investment larger.
These are usually further constrained by carbon emissions but in this occasion
we set a hard upper limit that is less than the variable generator.

```python
>>> n.parameter.converter_capacityparam.loc[("R1_data", "2030", "Dispatchable"), ("unitsLowerLimit", "unitsUpperLimit")] = (0, 3)
>>> n.parameter.converter_coefficient.loc[("Dispatchable", "2030", "Powergeneration", "Electricity"), ("coefficient")] = 1
>>> n.parameter.accounting_converterunits.loc[("OMFix", "R1_data", "Dispatchable", "2030"), "perUnitBuild"] = 30
>>> n.parameter.converter_techparam.loc[("Dispatchable", "2030"), "lifeTime"] = 40

```

The load of our system should look similar like the one of the last example.
The secured capacity always refers to the peak load of our model, in this case
6 units (which can be MW for example).

```python
>>> n.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"],] = [-n for n in [1, 2, 6, 2, 1, 1]]
>>> n.profile.sourcesink_profile.fillna(0, inplace=True)
>>> n.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Electricity"), ("usesFixedProfile")] = 1

```

To account for the secured capacity we include a new indicator lower bound and a
mapping to the indicators associated with the limit.

```python
>>> n.parameter.accounting_indicatorbounds.loc[("R1_data", "2030", "SecuredCapacity"), ("useLower", "lowerValue")] = (1, 6)
>>> n.parameter.accounting_perindicator.loc[("SecuredCapacity", "VariableCapacity", "R1_data", "2030"), "perIndicator"] = 0.2
>>> n.parameter.accounting_perindicator.loc[("SecuredCapacity", "DispatchCapacity", "R1_data", "2030"), "perIndicator"] = 1

```

We then add the contribution of each unit to the indicators.
This is a capacity factor associated to the availability of the converter.

```python
>>> n.parameter.accounting_converterunits.loc[("DispatchCapacity", "R1_data", "Dispatchable", "2030"), "perUnitBuild"] = 1
>>> n.parameter.accounting_converterunits.loc[("VariableCapacity", "R1_data", "VariableGen", "2030"), "perUnitBuild"] = 1

```

In its current state, the model is infeasible as there is not enough
dispatchable capacity available.
To have insights of the model we will need to add a slack variable which we can
do as follows.

```python
>>> n.parameter.accounting_indicatorbounds.loc[("R1_data", "2030", "SlackCapacity"), ("isVariable", "useLower", "lowerValue")] = (1, 1, 0)
>>> n.parameter.accounting_perindicator.loc[("SecuredCapacity", "SlackCapacity", "R1_data", "2030"), "perIndicator"] = 1
>>> n.parameter.accounting_perindicator.loc[("SlackCost", "SlackCapacity", "R1_data", "2030"), "perIndicator"] = 1
>>> n.parameter.accounting_perindicator.loc[("SystemCost", "SlackCost", "global", "horizon"), "perIndicator"] = 10000

```

Now the model should be feasible.

```python
>>> n.set.add(["2030"], "yearssel")
>>> n.datadir = "slack_model_2"
>>> n.write(fileformat="dat")
>>> ();n.run(timeend=6, read_result=True);()
(...)

```

The system costs are really high because of the need for slack.

```python
>>> n.result["indicator_accounting"].loc[("global", "horizon", "SystemCost"), "value"]
10140.0

```
