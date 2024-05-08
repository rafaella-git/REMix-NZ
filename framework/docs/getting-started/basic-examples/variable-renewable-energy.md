---
title: "Variable Renewable Energy Sources"
lang: en-US
---

(vre_intro)=

# Variable Renewable Energy (VRE) Sources

## VRE background

Output from variable renewable energy (VRE) power plants is based on
time-varying input.
Examples for VRE technologies are photovoltaic power plants, onshore and
offshore wind turbines as well as run-of-river power plants.

## VRE modeling in REMix

```{figure} /img/REMix_VRE.svg
:align: center
Figure: Modeling concept for VRE.
```

In REMix, VRE technologies are usually represented by converters with
electricity as single output commodity (Option 1).
They are a special case because they can be modeled without specifying an input
commodity (as we otherwise do for converters as e.g. for conventional power
plants).
They can also be modeled as source in REMix (Option 2). This could be a valid
option for you if you e.g. have a wind park with given capacity that you do not
want to optimize.
On this page, we give an introduction to both modeling options.

We set up a small optimization model with VRE. This means creating an instance
and setting cost minimization as objective.
We are also creating an exemplary model node `"R1_data"` in this first step.

```python
>>> from remix.framework.api.instance import Instance
>>> import numpy as np
>>> import pandas as pd


>>> idx = pd.IndexSlice

>>> m1 = Instance()
>>> m1.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> m1.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

```

After this general first step, the path splits up for Option 1 and Option 2.

### Option 1: Modeling VRE as converter with one output commodity only (Electricity)

In the accounting of VRE plants, mainly investment and fixed operations and
maintenance costs are relevant, so we are setting these up as model indicators.

```python
>>> m1.parameter.accounting_perindicator.loc[("SystemCost", "Invest", "global", "2030"), "perIndicator"] = 1
>>> m1.parameter.accounting_perindicator.loc[("SystemCost", "OMFix", "global", "2030"), "perIndicator"] = 1

```

The variability of an energy source (e.g. wind) is represented by a normalized
profile in high temporal resolution (e.g. hourly).
Each value has to be in the range between 0 and 1 and represents the share of
rated capacity of all installed power plants for a specific technology in a
specific data node.
This profile represents the respective potential of each renewable energy
technology in each time step.

For the sake of conciseness, our VRE profile will only be filled with six time steps.

```python
>>> m1.profile.converter_activityprofile.loc[("R1_data", "2030", "WindOnshore", "upper"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0.1, 0.35, 0.5, 0.84, 1, 0.46]

```

Power plant capacities comprise two types of capacities:

1. existing capacities (exogenous input to optimization) representing the
available power plant park in each data node (represented by lower limits).
2. newly installable capacities on the other hand may be introduced to optimize
the capacity of each VRE technology in each node for a specific optimization
period (e.g. a year).

For this optimization step, lower and upper limits for capacities are set in
order to represent political targets or potential limitations, respectively.

```python
>>> m1.parameter.converter_capacityparam.loc[("R1_data", "2030", "WindOnshore"), ("unitsLowerLimit", "unitsUpperLimit")] = (0, np.inf)

```

Optionally, a reduction factor can be given, reducing the converter activity to
e.g. 80 % of the rated capacity.
This refers to the total installed capacity (comprising previously available and
newly installed capacities minus decommissioned units).
Decommissioning of VRE power plants is accounted for by vintage years ("2030" in
the example).

```python
>>> m1.parameter.converter_coefficient.loc[("WindOnshore", "2030", "Powergeneration", "Electricity"), ("coefficient")] = 1

```

For the cost calculation of VRE power plant expansion, additional parameters
have to be given.
These comprise lifetime used for vintaging, specific investment per installed
capacity, economic amortization time, the technology-specific interest rate and
fixed operation and maintenance costs.
We do not need to set an upper limit for the activity of a VRE technology,
because we set an activity profile already (it will be set to "1" by the
Instance class anyway, without effect in this case).

```python
>>> m1.parameter.converter_techparam.loc[("WindOnshore", "2030"), ("lifeTime", "freeDecom")] = (25, 1)
>>> m1.parameter.accounting_converterunits.loc[("Invest", "R1_data", "WindOnshore", "2030"), ("perUnitBuild", "perUnitDecom", "amorTime", "interest", "useAnnuity")] = (30000, 300, 25, 0.06, 1)
>>> m1.parameter.accounting_converterunits.loc[("OMFix", "R1_data", "WindOnshore", "2030"), ("perUnitTotal")] = 80

```

It is common in energy system modeling that amortization time and anticipated
lifetime of a technology are assumed to be equal (here: 25 years).
With `"freeDecom"`, we can specify whether the technology with the given vintage
year (=year of construction) can be decommissioned before the end of its
lifetime.
`"useAnnuity"` is used to turn depreciation in the model for the technology on
("1") or off ("0").

Because through the activity of the WindOnshore technology, the commodity
`"Electricity"` is in the system.
Therefore, a sink is needed.
Pay attention that for sinks (as for flows outside the system boundaries in
general) a negative sign (`-`) is needed.

```python
>>> m1.profile.sourcesink_profile.loc[
...    ("R1_data", "2030", "Demand", "Electricity", "fixed"),
...    ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"],
... ] = [-n * 20 for n in [0.23, 0.12, 1, 0.38, 0.69, 0.72]]
>>> m1.profile.sourcesink_profile.fillna(0, inplace=True)
>>> m1.parameter.sourcesink_config.loc[
...    ("R1_data", "2030", "Demand", "Electricity"), ("usesFixedProfile")
... ] = 1

```

We are only running the first six time steps, because our VRE time series only
has the first six time steps filled anyway.
Two more things need to be done before the model is executable: (1) defining the
years for which we want to run the optimization, and (2) setting the directory
in which to write the model data.

```python
>>> m1.set.add(["2030"], "yearssel")
>>> m1.datadir = "vre_model_option1"
>>> m1.write(fileformat="dat")
>>> ();m1.run(timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

We can check the results, for example how many `WindOnshore` units are required
or how much electricity is generated:

```python
>>> units = m1.result["converter_units"].loc[("R1_model", "2030", "WindOnshore", "2030", "build"), "value"]
>>> round(units, 1)
46.0
>>> total_production = m1.result["commodity_balance_annual"].loc[("R1_model", "2030", "WindOnshore", "Electricity", "netto"), "value"]
>>> round(total_production, 1)
62.8

```

In the post-processing, we can see the curtailment of the wind energy, i.e. the
amount of energy available in theory but not delivered by the windpower plant
due to lack of demand.
To calculate that, we can calculate the actual per-unit production from the
`commodity_balance` and compare that with the upper production limit specified
in the input data.
The time series resulting from that multiplied with the amount of converter
units will give us the total curtailed electricity.

```python
>>> actual_production = m1.result["commodity_balance"].loc[idx[:, "R1_model", "2030", "WindOnshore", "Electricity"], "value"].values
>>> per_unit_upper_limit = m1.profile.converter_activityprofile.loc[("R1_data", "2030", "WindOnshore", "upper")].dropna().values
>>> total_curtailment = (per_unit_upper_limit * units).sum() - actual_production.sum()
>>> round(total_curtailment, 1)
86.7

```

We can check whether our model was set up correctly, e.g. by comparing the total
curtailed and produced energy with the upper production limit.

```python
>>> total_curtailment + total_production == (per_unit_upper_limit * units).sum()
True

```

### Option 2: Modeling VRE as source

Besides this first and usually implemented option, it is also possible to model
renewables as source in REMix.
We need to set an upper profile to do that, which we can derive from the
`converter_activityprofile` from above.

As `sourcesink_profile` we can use the normalized `converter_activityprofile`
from above and multiply it with the wind park capacity `cap_wind`.

```python
>>> m2 = Instance()
>>> m2.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> m2.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

>>> cap_wind = units  # 46 units
>>> m2.profile.sourcesink_profile.loc[("R1_data", "2030", "WindOnshore", "Electricity", "upper"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [n * cap_wind for n in [0.1, 0.35, 0.5, 0.84, 1, 0.46]]
>>> m2.parameter.sourcesink_config.loc[("R1_data", "2030", "WindOnshore", "Electricity"), ("usesUpperProfile")] = 1

```

We do not need a `sourcesink_annualsum` in this case, because the upper limit of
electricity produced by wind turbines is already set through the
`sourcesink_profile`.
What we do need, is a sink for the model to run as intended.
We are just taking the same as for `m1`.

```python
>>> m2.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [-n * 20 for n in [0.23, 0.12, 1, 0.38, 0.69, 0.72]]
>>> m2.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Electricity"), ("usesFixedProfile")] = 1

```

Then we can run the optimization for option 2.

```python
>>> m2.set.add(["2030"], "yearssel")
>>> m2.datadir = "vre_model_option2"
>>> m2.write(fileformat="dat")
>>> ();m2.run(resultfile="vre_model_option2", timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

Again, we can have a look at the results.
Obviously, the total production should be the same.
We can also calculate the curtailment as we did in the first example, and it
should produce the same result.
Note that the upper limit is not per unit anymore but is the actual upper
production limit.
Therefore, we do not need to take the number of units into account.

```python
>>> total_production = m2.result["commodity_balance_annual"].loc[("R1_model", "2030", "WindOnshore", "Electricity", "netto"), "value"]
>>> round(total_production, 1)
62.8
>>> actual_production = m2.result["commodity_balance"].loc[idx[:, "R1_model", "2030", "WindOnshore", "Electricity"], "value"].values
>>> upper_limit = m2.profile.sourcesink_profile.loc[("R1_data", "2030", "WindOnshore", "Electricity", "upper")].dropna().values
>>> total_curtailment = upper_limit.sum() - actual_production.sum()
>>> round(total_curtailment, 1)
86.7

```

As already discussed in the introduction, a source cannot be optimized in terms
of its capacity in contrast to the first option presented.
For instance, we can change the wind capacity in a way that the demand cannot
be satisfied with the capacity available, and the model should be infeasible:

```python
>>> cap_wind /= 1.5
>>> m2.profile.sourcesink_profile.loc[("R1_data", "2030", "WindOnshore", "Electricity", "upper"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [n * cap_wind for n in [0.1, 0.35, 0.5, 0.84, 1, 0.46]]
>>> m2.write(fileformat="dat")
>>> ();m2.run(resultfile="vre_model_option2_infeasible", timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

(vre_curtailment_modeling_label)=

## Curtailment as a variable or constraint (only with Option 1)

In REMix, it is possible to include curtailment also as a variable or a
constraint inside the optimization.
This is done by including a fixed activity profile for the converter
technologies instead of an upper limit.
An additional activity is assigned, so the technology has one that produces
electricity and one that represents curtailment.
In our example we now have two activities `"Powergeneration"` and
`"Curtailment"` defined for the `"WindOnshore"` technology.
Either one of them must be used in any one time step or any linear combination
of both.

```python
>>> m1.parameter.converter_coefficient.loc[("WindOnshore", "2030", "Curtailment", "Electricity_curtailed"), ("coefficient")] = 1.0
>>> m1.profile.converter_activityprofile.loc[("R1_data", "2030", "WindOnshore", "fixed")] = m1.profile.converter_activityprofile.loc[("R1_data", "2030", "WindOnshore", "upper")]
>>> m1.profile.converter_activityprofile.drop(("R1_data", "2030", "WindOnshore", "upper"), inplace=True)

```

We also need a sink for `"Electricity_curtailed"`.

```python
>>> m1.parameter.sourcesink_annualsum.loc[
...    ("R1_data", "2030", "Curtailment", "Electricity_curtailed"), "lower"
... ] = -np.inf
>>> m1.parameter.sourcesink_config.loc[
...    ("R1_data", "2030", "Curtailment", "Electricity_curtailed"), ("usesLowerSum")
... ] = 1

```

```python
>>> m1.datadir = "vre_model_option1_curtailment"
>>> m1.write(fileformat="dat")
>>> ();m1.run(resultfile="vre_model_option1_curtailment", timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

First, we can compare the results between both variants and see that the
curtailed electricity is still the same value but can be directly read from the
results file.

```python
>>> total_production = m1.result["commodity_balance_annual"].loc[("R1_model", "2030", "WindOnshore", "Electricity", "netto"), "value"]
>>> round(total_production, 1)
62.8
>>> total_curtailment = m1.result["commodity_balance_annual"].loc[("R1_model", "2030", "WindOnshore", "Electricity_curtailed", "netto"), "value"]
>>> round(total_curtailment, 1)
86.7
>>> total_costs_without_penalty = m1.result["indicator_accounting"].loc[("R1_model", "2030", "SystemCost"), "value"]

```

The curtailment is only used for post-processing as well, so let us change that:
You might want to penalize curtailment by adding marginal cost for every unit of
electricity that is curtailed.
For this, we add an indicator our objective and connect the indicator with the
flow into the curtailment sink.
The flow of the curtailed electricity into the sink will be negative from the
model's perspective, therefore we use a negative sign for the per flow cost, for
example -10 €/MWh.

```python
>>> m1.parameter.accounting_perindicator.loc[("SystemCost", "CurtailmentCost", "global", "2030"), "perIndicator"] = 1
>>> m1.parameter.accounting_sourcesinkflow.loc[("CurtailmentCost", "global", "2030", "Curtailment", "Electricity_curtailed"), "perflow"] = -10
>>> m1.datadir = "vre_model_option1_curtailment"
>>> m1.write(fileformat="dat")
>>> ();m1.run(resultfile="vre_model_option1_curtailment", timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

After rerunning the model, we see the `CurtailmentCost` in the indicator
accounting and the change in `SystemCost`.

```python
>>> round(total_costs_without_penalty)
111633
>>> round(m1.result["indicator_accounting"].loc[("R1_model", "2030", "SystemCost"), "value"])
112500
>>> round(m1.result["indicator_accounting"].loc[("R1_model", "2030", "CurtailmentCost"), "value"], 1)
867.0
>>> round(total_curtailment * m1.parameter.accounting_sourcesinkflow.loc[("CurtailmentCost", "global", "2030", "Curtailment", "Electricity_curtailed"), "perflow"], 1)
867.0

```

## Disclaimer

There is a big heterogeneity between the VRE technologies.
For details of VRE modeling refer to
{cite}`Scholz_2012, Pfenniger_2016, Staffell_2016`.
The available renewable energy varies on seasonal and daily scales.
In order to grasp this variation, (sub-)hourly time series are used to represent
discrete historic years of this variability.
Please be aware of limitations of drawing policy-relevant conclusions if not
accounting for inter-annual variability {cite}`Collins_2018`.
