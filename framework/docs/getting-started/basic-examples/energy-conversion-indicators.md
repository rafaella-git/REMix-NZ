---
title: "Using Indicators"
lang: en-US
---

(simple_indicator_intro)=

# Using Indicators

In this example, we will have a more detailed look on indicators at the example
of tracking CO<sub>2</sub> emissions of a thermal power plant.

## Introduction

Thermal power plants, for example coal power plants, convert a (fossil) fuel
into electricity (and CO<sub>2</sub>).
The figure below shows a possible model implementation of such a plant in REMix.

```{figure} /img/REMix_CoalPP.svg
:align: center
Figure: Schema of a thermal power plant.
```

Therefore, the following REMix components are required:

- {ref}`a source and two sinks <remix_model_core_sourcesink_label>` for coal
  imports as well as electricity demand and the CO<sub>2</sub> emissions and
- {ref}`a converter <remix_model_core_converter_label>`.

In this example, the power plant will provide electricity following a demand
profile.
The table below shows an overview of the techno-economic data we want to provide
to our REMix model:

| Component   | Parameter                      |    Value | Unit    | Comment                 |
|:----------- |:------------------------------ | --------:|:------- |:----------------------- |
| Fuel        | Minimum flow value             |    0.000 | GW      | no exports              |
|             | Cost                           |    0.010 | M€/GWh  |                         |
|             | CO<sub>2</sub> emission factor |    0.360 | kt/GWh  |                         |
| Power Plant | Electric efficiency            |   40.000 | %       |                         |
|             | Capacity                       |    5.000 | GW      | relative to electricity |
|             | Availability factor            |   90.000 | %       |                         |
|             | Variable cost                  |   0.0001 | M€/GWh  | relative to electricity |
| Emissions   | Cost                           |    0.060 | M€/kt   |                         |
| Electricity | Demand profile                 |        - | -       |                         |

For the fuel we want to make sure it can only be imported.
We therefore set a lower value to the flow of 0.0 GW.
For example in case you have a networked system with links and do not want
to allow selling of fuel this is an important setting.
In this example this constraint is not strictly required, but we still introduce
it to show how it can be applied.
Importing a fuel is connected with marginal cost; marginal cost for the import
is stated as well.
On top of that, we know the specific CO<sub>2</sub> emission factor associated
with the coal.

For the power-plant model we want to set efficiency, a maximum capacity for
extension and cost-related terms, i.e. investment and operational cost.
Also, there are options to take amortization time and interest rate into account
for the economic evaluation.
The cost of emissions is provided as well.

## Tracking emissions with a commodity

### Build the REMix model

First, we have to define how we want to set up the REMix model.
Here REMix also offers a lot of flexibility, similar to what we have seen in the
{ref}`introductory example on simple energy conversion <simple_energy_conversion_intro>`.
Since we want to keep track of the CO<sub>2</sub> emissions specifically in this
example, we will first show how to model the emissions as an own commodity.

```{figure} /img/REMix_CoalPP_accounting.svg
:align: center
Figure: Schema of thermal power plant implementation in REMix with accounting.
```

As in the previous example, we only need a single region, and we will make use of
the {class}`remix.framework.api.instance.Instance` class.

```python
>>> from remix.framework.api.instance import Instance
>>> import numpy as np


>>> m = Instance()
>>> m.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""
>>> m.set.yearssel = ["2030"]

```

After creating our `Instance` object and adding the data-model lookup, we can
add the technology data from the table.
First, we add the lower profile for the fuel import, which sets a lower limit of
0.
Then we can add the cost for the coal imports.
The marginal costs are applied with the `indicators` and the respective value
per flow of the commodity.

```python
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "FuelImport", "Coal"), "usesLowerProfile"] = 1
>>> m.parameter.accounting_sourcesinkflow.loc[("FuelCost", "global", "2030", "FuelImport", "Coal"), "perFlow"] = 0.01

```

Next, we can create the thermal power plant.
We begin by specifying the number of converter units available to the model and
apply the approach with multiple units of a normed size as shown in
{ref}`this section <simple_energy_conversion_multiple_units>`, because our data
is related to the production of electricity.
For an efficiency of 40 % the plant needs to consume 2.5 units of coal for every
unit of electricity produced.
The CO<sub>2</sub> emissions are linked with the CO<sub>2</sub> emission factor:
per unit of coal, 0.36 units of CO<sub>2</sub> are emitted.
Since we have to make the definition of the CO<sub>2</sub> emissions relative to
the number of units and the activity, {math}`0.36 \cdot 2.5 = 0.9` units of
CO<sub>2</sub> are emitted per activity of the plant.

```python
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "CoalPowerPlant"), "unitsBuild"] = 5
>>> m.parameter.converter_coefficient.loc[("CoalPowerPlant", "2030", "Powergeneration", "Electricity"), "coefficient"] = 1
>>> m.parameter.converter_coefficient.loc[("CoalPowerPlant", "2030", "Powergeneration", "Coal"), "coefficient"] = -2.5
>>> m.parameter.converter_coefficient.loc[("CoalPowerPlant", "2030", "Powergeneration", "CO2"), "coefficient"] = 0.9

```

Next, we can account for the variable cost of production which is linked through
the `activity` value.
With our setup, 1 GWh of electrical energy is generated per activity per unit of
the thermal power plant.
Therefore, the value of the variable cost of operation per activity is the same
relative to the electricity production.

```python
>>> m.parameter.accounting_converteractivity.loc[("OMVar", "global", "CoalPowerPlant", "2030", "Powergeneration"), "perActivity"] = 0.001

```

Finally, we have to add a `lifeTime` value other than zero, otherwise the
converter is not available in the model.

```python
>>> m.parameter.converter_techparam.loc[("CoalPowerPlant", "2030"), "lifeTime"] = 40

```

The electricity-demand time series needs to be specified and the CO<sub>2</sub>
emission constraints have to be set.
It is a very simple ramp here, since the demand model is not of particular
interest for this introduction.

```python
>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0, -1, -2, -3, -4, -5]
>>> m.profile.sourcesink_profile.fillna(0, inplace=True)
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Electricity"), "usesFixedProfile"] = 1

```

We allow for unlimited CO<sub>2</sub> emissions but apply a variable cost to them:

```python
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Emissions", "CO2"), "usesLowerSum"] = 1
>>> m.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "Emissions", "CO2"), "lower"] = -np.inf
>>> m.parameter.accounting_sourcesinkflow.loc[("CO2Cost", "global", "2030", "Emissions", "CO2"), "perFlow"] = 0.06

```

Finally, we set up the indicator accounting to build the objective function.
We can also set the ending time step to 6 here, since we only provided demand
data for 6 time steps.

```{caution}
The factor for the CO<sub>2</sub> emissions is negative.
This is because of the setup with CO<sub>2</sub> as a separate commodity.
It leaves the system's boundaries and therefore the flow values of the emissions
to the emission sink are negative.
To not subtract the value of the CO<sub>2</sub> emission costs from the other
costs, but to add it, we therefore have to change the sign by multiplying it
with -1.
```

```python
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "OMVar", "global", "2030"), "perIndicator"] = 1
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "FuelCost", "global", "2030"), "perIndicator"] = 1
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "CO2Cost", "global", "2030"), "perIndicator"] = -1
>>> m.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1

```

### Run the model

Now we can run the model and have a look at the results.

```python
>>> m.datadir = "conversion_model"
>>> m.write()
>>> ();m.run(timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

For example, we can obtain the CO<sub>2</sub> emissions from the results' annual
commodity balance value, or the objective value for the system cost:

```python
>>> m.result["commodity_balance_annual"].loc[("global", "2030", "Emissions", "CO2", "netto"), "value"]
-13.5

>>> m.result["indicator_accounting"].loc[("global", "2030", "SystemCost"), "value"]
1.2

```

Now, let's change it up a little:
we can define negative `CO2Cost` per flow of emissions (the emissions flow will
always be negative to the system as explained in the previous section).
When doing that, we have to change the `perIndicator` factor to positive 1 and
rerun the model.

```python
>>> m.parameter.accounting_sourcesinkflow.loc[("CO2Cost", "global", "2030", "Emissions", "CO2"), "perFlow"] = -0.06
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "CO2Cost", "global", "2030"), "perIndicator"] = 1
>>> m.write()
>>> ();m.run(timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

We can see that there is no difference in any of the results:

```python
>>> m.result["commodity_balance_annual"].loc[("global", "2030", "Emissions", "CO2", "netto"), "value"]
-13.5

>>> m.result["indicator_accounting"].loc[("global", "2030", "SystemCost"), "value"]
1.2

```

```{note}
In this example, the commodity CO<sub>2</sub> is not used for purposes other
than accounting for the emission cost.
Since CO<sub>2</sub> is a model commodity here, it would however be possible to
reuse the CO<sub>2</sub> in other processes, for example, for methanation of
hydrogen or in carbon capture and storage.
The decision of modeling a component with or without commodities is always a
question of your individual modeling requirements.
```

Therefore, the next subsection presents a different approach to model the same
system without introducing CO<sub>2</sub> as a commodity but linking it to coal
imports.

## Tracking emissions with an indicator

In this section, a second approach of modeling the same system as above with an
identical outcome is presented and investigated.

### Build the REMix model

In our first model there was no real need for a CO<sub>2</sub> commodity, since
it did not have any other purpose than in the accounting.
This means we introduced more mathematical and modeling complexity than
necessary with the additional commodity balance for the CO<sub>2</sub>.
It is possible to keep track of the CO<sub>2</sub> emissions without creating a
dedicated commodity for it.
The approach for that is shown in the figure below.


```{figure} /img/REMix_CoalPP_accounting_no_CO2_commodity.svg
:align: center
Figure: Schema of the thermal power plant implementation in REMix without a
separate CO<sub>2</sub> commodity.
```

To model this, we have to modify our data accordingly.
First, we drop all data we do not need anymore.
These are the emissions caused by the power plant through the converter, the
CO<sub>2</sub> emission sink and the associated CO<sub>2</sub> cost with the
flow of CO<sub>2</sub> emissions to the sink.

```python
>>> m.parameter.converter_coefficient.drop(("CoalPowerPlant", "2030", "Powergeneration", "CO2"), inplace=True)

>>> m.parameter.sourcesink_annualsum.drop(("R1_data", "2030", "Emissions", "CO2"), inplace=True)
>>> m.parameter.sourcesink_config.drop(("R1_data", "2030", "Emissions", "CO2"), inplace=True)
>>> m.parameter.accounting_sourcesinkflow.drop(("CO2Cost", "global", "2030", "Emissions", "CO2"), inplace=True)

```

We can modify the accounting for the `FuelImport` source by adding the
associated CO<sub>2</sub> emissions.
The emissions' indicator can be weighted into the `SystemCost` with the specific
cost of the emissions.

```python
>>> m.parameter.accounting_sourcesinkflow.loc[("CO2Emissions", "global", "2030", "FuelImport", "Coal"), "perFlow"] = 0.36
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "CO2Emissions", "global", "2030"), "perIndicator"] = 0.06

```

```{note}
The sign of the cost is positive in this approach.
This is because the flow of fuel to which the emissions are linked to is always
directed into the system's boundaries in contrary to the CO<sub>2</sub>
emissions leaving it in the first example.
```

### Run the model

We can now run the model again like in the first example.

```python
>>> m.datadir = "conversion_model_without_CO2_commodity"
>>> m.write()
>>> ();m.run(timeend=6, read_result=True);()  # +doctest: ELIPSIS
(...)

```

We can verify that we get the same results.
The `CO2Emissions` are now obtained from the indicator accounting instead of the
commodity balance, since CO<sub>2</sub> is not a commodity anymore.
Note that the emission's value is positive since it is tied to the import of the
fuel with a positive value.

```python
>>> m.result["indicator_accounting"].loc[("global", "2030", "CO2Emissions"), "value"]
13.5

>>> m.result["indicator_accounting"].loc[("global", "2030", "SystemCost"), "value"]
1.2

```
