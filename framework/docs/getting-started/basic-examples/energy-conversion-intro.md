---
title: "Simple Energy Conversion"
lang: en-US
---

(simple_energy_conversion_intro)=

# Simple Energy Conversion

In this section you can learn how the converter component of REMix is working at the example of a methanation plant.
We show the flexibility the modeler has when using REMix to model energy conversion units.
There are different approaches that produce the same outcome but have implications on and result in restrictions to your model.

```{tip}
This example shows the usage of a linear model for the power plant.
REMix offers methods to create mixed-integer models as well, which can be used, for example, for discrete expansion planning or part-load modeling.
For further information on that see {ref}`MIP features <explanations_mip_label>`.
```

## Introduction

In this example a methanation unit is built, which takes hydrogen and carbon dioxide as inputs.
It produces heat and methane (and water, which is not further considered in this example) as shown in the figure below.

```{figure} /img/REMix_Methanation.svg
:align: center
Figure: Schema of a methanation plant.
```

The technical data for the methanation is listed in the table below.
Assuming a simplified and ideal process without by-products, incomplete chemical reactions or heat losses to the environment, the reaction energy balance and stoichiometry provide the energy and mass conversion data.
Per MWh of methane produced, 1.2 MWh hydrogen and 198 kg CO<sub>2</sub> are consumed.
Furthermore, about 0.2 MWh of heat are produced.
We assume a simple demand profile with a few demand values.

| Component   | Parameter                  | Value | Unit    | Comment                 |
|:------------|:---------------------------|------:|:--------|:------------------------|
| Hydrogen    | Cost                       |  50.0 | €/MWh   |                         |
| Methanation | Methane generation         | 100.0 | MWh/h   |                         |
|             | Hydrogen consumption       |   1.2 | MWh/MWh | per energy of methane   |
|             | CO<sub>2</sub> consumption | 198.0 | kg/MWh  | per energy of methane   |
|             | Heat generation            |   0.2 | MWh/MWh | per energy of methane   |
|             | Variable operation cost    |   3.0 | €/MWh   | per energy of methane   |
|             | Investment cost            | 400.0 | €/kW    | per generation capacity |
|             | Revenue for sales of heat  |  50.0 | €/MWh   | per energy of heat      |
| Methane     | Demand profile             |     - | -       |                         |

## Model a single unit

### Model setup

To model a single methanation unit with a capacity of 100 MW in REMix, we create an object of the {class}`remix.framework.api.instance.Instance` class:

```python
>>> from remix.framework.api.instance import Instance
>>> import numpy as np


>>> su = Instance()  # su: single unit

```

The next step is to create a model data/model node lookup table.
This is called `aggregatenodesmodel`. We will call our data node `R1_data` and our model node `R1_model`.

```python
>>> su.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

```

```{tip}
The `aggregatenodesmodel` is a lookup for data regions with model regions.
If you want to have a less detailed spatial representation of your mathematical model compared to the data model, the `aggregatenodesmodel` will define which model nodes are aggregations of which data nodes.
We introduce the usage of that in a later example.
```

To model the methanation unit, we can begin with the sources, i.e. hydrogen and CO<sub>2</sub>.
To create them, we can simply add the necessary information from the data table to our Instance object.
We add a configuration which limits the total import of the commodity into the system:

```python
>>> su.parameter.sourcesink_config.loc[("R1_data", "2030", "HydrogenImport", "Hydrogen"), "usesUpperSum"] = 1
>>> su.parameter.sourcesink_config.loc[("R1_data", "2030", "CO2Import", "CO2"), "usesUpperSum"] = 1

```

Then, we add a value for the maximum consumption of hydrogen and CO<sub>2</sub>, which is unlimited in the example:

```python
>>> su.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "HydrogenImport", "Hydrogen"), "upper"] = np.inf
>>> su.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "CO2Import", "CO2"), "upper"] = np.inf

```

Last, we add the cost for hydrogen like so:

```python
>>> su.parameter.accounting_sourcesinkflow.loc[("FuelCost", "global", "2030", "HydrogenImport", "Hydrogen"), "perFlow"] = 50 # €/MWh

```

The next step is to add the conversion unit.
We do that by specifying the amount of units available, which is a single unit initially:

```python
>>> su.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "unitsBuild"] = 1

```

With the unit in place, the information on consumption and generation of commodities linked to a specific activity of the unit is added.
From our data table we can see that our unit should produce 100 MWh of methane per hour as output.
As inputs 198 kg CO<sub>2</sub> per MWh of methane and 1.2 MWh of hydrogen are required.
Last, 20 MWh of heat are generated per 100 MWh of methane generation.
Both, generation and consumption data need to be specified in absolute values and are called `coefficient`:

```{note}
The sign of the `coefficient` is negative in case of consumption and positive for generation.
```

```python
>>> su.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Methane"), "coefficient"] = 100 # MW per converter unit
>>> su.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Hydrogen"), "coefficient"] = -120 # MW per 100 MWh production of methane
>>> su.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "CO2"), "coefficient"] = -19_800 # per 100 MW production of methane
>>> su.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Heat"), "coefficient"] = 20 # per 100 MW production of methane

```

```{note}
The connection between the inputs and outputs of the converter units is established by the so called `activity`.
The value of the activity is between 0 and 1, where 0 means the units are shut down and 1 means they produce at full load.
For our Methanation we therefore create an activity called `Methanegeneration`.
In our example, it produces 100 MW of the commodity `Methane` per unit available and per `activity`.
```

We can account for the variable cost of generation which is also linked through the `activity` value.
With our setup, 100 MWh of electrical energy are generated per activity per unit of the methanation per time step.
Since the variable cost given from the table is relative to the generation of a single MWh of methane, we have to transform this information to make it relative to the generation of 100 MWh by multiplying the specific cost with the total generation per activity.

```python
>>> su.parameter.accounting_converteractivity.loc[("OMVar", "global", "Methanation", "2030", "Methanegeneration"), "perActivity"] = 3 * 100 # 300 € per 100 MWh production of methane

```

In a first setup we only want to look at the operation of the methanation and do not investigate investment.
Therefore, we only need a few more settings to finalize our model.
We have to add a `lifeTime` value other than zero, otherwise the converter is not available in the model. We also have to add `noExpansion` to define that the model has not the degree of freedom to build new `Methanation` units or decomission the existing one:

```python
>>> su.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "noExpansion"] = 1 # using no expansion/decomissioning
>>> su.parameter.converter_techparam.loc[("Methanation", "2030"), "lifeTime"] = 40 # years

```

A methane-demand time series is specified, which the methanation can provide the produced methane for.
We can make the converter follow a simple ramp like so:

```python
>>> su.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Methane", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0, -20, -40, -60, -80, -100] # MWh energy of methane in every hour
>>> su.profile.sourcesink_profile.fillna(0, inplace=True)
>>> su.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Methane"), "usesFixedProfile"] = 1 # using fixed profile

```{note}

Similar to the converters, negative values indicate that a commodity is consumed. In REMix it is called a `sink` when the commodity leaves the system.
Positive values automatically indicate an import or the generation of a commodity.
```

For the heat demand we create a sink that does not follow a time series but will result in a profit per MWh of heat generated to keep it simple in this example.

```python
>>> su.parameter.sourcesink_config.loc[("R1_data", "2030", "HeatDemand", "Heat"), "usesLowerSum"] = 1 # use a sum as lower limit for the heat demand
>>> su.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "HeatDemand", "Heat"), "lower"] = -np.inf # The amount of heat which can leave the system is infinite. The default would be 0.
>>> su.parameter.accounting_sourcesinkflow.loc[("HeatRevenue", "global", "2030", "HeatDemand", "Heat"), "perFlow"] = 50 # €/MWh revenue for every flow of heat leaving the system

```

Finally, the indicator accounting defines which indicators are part of the objective function and whether the objective should be minimized or maximized.
We can also set the ending time step to 6 here, since we only provide demand data for 6 time steps.
The `indicatorbounds` for the `SystemCost` indicator is `-1`, which means that we want to run a minimization.
The `accounting_perindicator` defines which other indicators (variable costs, fuel cost, heating revenue, ...) are summed up to calculate the total system cost and what the factor is (`perIndicator`) they are multiplied with in that sum.
The `yearssel` value indicates which year should be optimized.
This value has to be specified regardless of whether there is a single year only or not.
For the variable cost, cost of fuel and the heat revenue, we have to provide the cost factor per indicator.
The heat revenue is an income, therefore its factor has to be negative.

```python
>>> su.parameter.accounting_perindicator.loc[("SystemCost", "OMVar", "global", "2030"), "perIndicator"] = 1
>>> su.parameter.accounting_perindicator.loc[("SystemCost", "FuelCost", "global", "2030"), "perIndicator"] = 1
>>> su.parameter.accounting_perindicator.loc[("SystemCost", "HeatRevenue", "global", "2030"), "perIndicator"] = -1
>>> su.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> su.set.yearssel = ["2030"]

```

Now we can run the model with first defining the directory where to write the data and take it for the model run.

```python
>>> su.datadir = "methanation-model-single-unit"
>>> su.write()
>>> ();su.run(timeend=6, read_result=True, resultdir=su.datadir);()  # +doctest: ELIPSIS
(...)

```

### Check the results

We can check in the results whether our setup yields the expected outcome.
For example, we can sum up the total methane generation and check whether the variable cost of operation of the methanation matches the specific cost multiplied by the total generation:
To do so, we extract the total generation of the methanation for the R1 region from the annual commodity balance:

```python
>>> su_methane_total = su.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "Methane", "netto"), "value"]
>>> su_methane_total
300.0

```

Then we check if the total variable cost of operation is equal to the expected value, i.e. `300 * 3 = 900`.

```python
>>> su.result["indicator_accounting"].loc[("R1_model", "2030", "OMVar"), "value"]
900.0

```

Also, we can check whether our specifications for hydrogen and CO<sub>2</sub> consumption are correct.
The annual commodity balance gives us the total consumption of these commodities, which we can then be divided by the total methane generation.
By this we find out if the specific consumption is set up correctly:

```python
>>> su_CO2_total = su.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "CO2", "brutto"), "value"]
>>> su_hydrogen_total = su.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "Hydrogen", "brutto"), "value"]
>>> su_CO2_total / su_methane_total
198.0
>>> su_hydrogen_total / su_methane_total
1.2

```

The next step is to make the same setup but use multiple units of the methanation plant with smaller capacity individually.

(simple_energy_conversion_multiple_units)=

## Multiple units

As we learned, the converter model in REMix is unit- and activity-based:
Generation and consumption of commodities are per activity of each unit.
Therefore, instead of considering a single large unit we can consider multiple units.
Having multiple units might especially make sense for large-scale and long-term energy system analysis, since it is not common to set the focus on single plants and/or locations.

### Model setup

We start with a new `Instance` object and the configuration of the sources and the demand since these are not affected by changing the number of methanation units:

```python
>>> mu = Instance()  # mu: multi-unit
>>> mu.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""

>>> mu.parameter.sourcesink_config.loc[("R1_data", "2030", "HydrogenImport", "Hydrogen"), "usesUpperSum"] = 1
>>> mu.parameter.sourcesink_config.loc[("R1_data", "2030", "CO2Import", "CO2"), "usesUpperSum"] = 1
>>> mu.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "HydrogenImport", "Hydrogen"), "upper"] = np.inf
>>> mu.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "CO2Import", "CO2"), "upper"] = np.inf
>>> mu.parameter.accounting_sourcesinkflow.loc[("FuelCost", "global", "2030", "HydrogenImport", "Hydrogen"), "perFlow"] = 50


>>> mu.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Methane", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0, -20, -40, -60, -80, -100]
>>> mu.profile.sourcesink_profile.fillna(0, inplace=True)
>>> mu.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Methane"), "usesFixedProfile"] = 1

>>> mu.parameter.sourcesink_config.loc[("R1_data", "2030", "HeatDemand", "Heat"), "usesLowerSum"] = 1
>>> mu.parameter.sourcesink_annualsum.loc[("R1_data", "2030", "HeatDemand", "Heat"), "lower"] = -np.inf
>>> mu.parameter.accounting_sourcesinkflow.loc[("HeatRevenue", "global", "2030", "HeatDemand", "Heat"), "perFlow"] = 50

```

Now, we can decide how large the individual units can get.
A clever way to do this is to assign a coefficient of `1` per activity to our main generation commodity (here: `Methane`).

```python
>>> mu.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Methane"), "coefficient"] = 1

```

```{note}
Because of this

- one unit of methanation will produce 1 MWh of methane per activity.
- one unit of methanation has a capacity of 1 MW.

Therefore

- all specifications per activity are relative to the generation of 1 MWh of methane.
- all specifications per unit are relative to 1 MW capacity of methanation.
```

That means we will need 100 units for our total capacity of 1 MW.
For the consumption of hydrogen and CO<sub>2</sub> we can simply put in our values from the overview table since these are already relative to the generation of 1 MWh of methane.
The same is true for the variable costs.

```python
>>> mu.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "unitsBuild"] = 100
>>> mu.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "noExpansion"] = 1
>>> mu.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Hydrogen"), "coefficient"] = -1.2
>>> mu.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "CO2"), "coefficient"] = -198
>>> mu.parameter.converter_coefficient.loc[("Methanation", "2030", "Methanegeneration", "Heat"), "coefficient"] = 0.2
>>> mu.parameter.accounting_converteractivity.loc[("OMVar", "global", "Methanation", "2030", "Methanegeneration"), "perActivity"] = 3
>>> mu.parameter.converter_techparam.loc[("Methanation", "2030"), "lifeTime"] = 40

```

We finish the setup and run the model:

```python
>>> mu.parameter.accounting_perindicator.loc[("SystemCost", "OMVar", "global", "2030"), "perIndicator"] = 1
>>> mu.parameter.accounting_perindicator.loc[("SystemCost", "FuelCost", "global", "2030"), "perIndicator"] = 1
>>> mu.parameter.accounting_perindicator.loc[("SystemCost", "HeatRevenue", "global", "2030"), "perIndicator"] = -1
>>> mu.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1
>>> mu.set.yearssel = ["2030"]

>>> mu.datadir = "methanation-model-multi-unit"
>>> mu.write()
>>> ();mu.run(timeend=6, read_result=True, resultdir=mu.datadir);()  # +doctest: ELIPSIS
(...)

```

### Check the results

We can compare the results of both models:

```python
>>> su_methane_total == mu.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "Methane", "netto"), "value"]
True
>>> su.result["indicator_accounting"].loc[("R1_model", "2030", "OMVar"), "value"] == mu.result["indicator_accounting"].loc[("R1_model", "2030", "OMVar"), "value"]
True
>>> su_CO2_total == mu.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "CO2", "brutto"), "value"]
True
>>> su_hydrogen_total == mu.result["commodity_balance_annual"].loc[("R1_model", "2030", "Methanation", "Hydrogen", "brutto"), "value"]
True

```

## So what is the difference?

Well, apparently there is none.
We can check what happens if we allow unit expansion.
For that, we set the number of existing units to zero and set the `noExpansion` setting to `False` in both models.
Additionally, we add the investment cost to our accounting.
Beware that the investment costs have to be relative to the amount of units, therefore the individual unit's size.
For the single unit that is 400,000 €/MW and a single unit has 100 MW capacity, therefore 400,000 * 100 €/unit.
For the multiple-units approach a single unit has 1 MW capacity, therefore 400,000 €/unit.
Then add the new `Investment` indicator to the accounting.

```python
>>> su.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "unitsBuild"] = 0
>>> su.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "noExpansion"] = False
>>> su.parameter.accounting_converterunits.loc[("Investment", "R1_data", "Methanation", "2030"), "perUnitBuild"] = 400000 * 100
>>> su.parameter.accounting_perindicator.loc[("SystemCost", "Investment", "global", "2030"), "perIndicator"] = 1

>>> su.write()
>>> ();su.run(timeend=6, read_result=True, resultdir=su.datadir);()  # +doctest: ELIPSIS
(...)

>>> mu.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "unitsBuild"] = 0
>>> mu.parameter.converter_capacityparam.loc[("R1_data", "2030", "Methanation"), "noExpansion"] = False
>>> mu.parameter.accounting_converterunits.loc[("Investment", "R1_data", "Methanation", "2030"), "perUnitBuild"] = 400000
>>> mu.parameter.accounting_perindicator.loc[("SystemCost", "Investment", "global", "2030"), "perIndicator"] = 1
>>> mu.write()
>>> ();mu.run(timeend=6, read_result=True, resultdir=mu.datadir);()  # +doctest: ELIPSIS
(...)

```

Now, we check the results for both models, specifically, the total units available in each model and the investment cost.

```python
>>> su.result["converter_units"].loc[("R1_model", "2030", "Methanation", "2030", "total"), "value"]
1.0
>>> mu.result["converter_units"].loc[("R1_model", "2030", "Methanation", "2030", "total"), "value"]
100.0

>>> su.result["indicator_accounting_detailed"].loc[("Investment", "R1_model", "2030", "Methanation"), "value"]
40000000.0
>>> mu.result["indicator_accounting_detailed"].loc[("Investment", "R1_model", "2030", "Methanation"), "value"]
40000000.0

```

It is not a surprise that we see this result.
We can modify the demand time series (effectively changing to total methanation capacity required), and see how it does affect the unit extension.

```python
>>> su.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Methane", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0, -10, -20, -30, -40, -50]
>>> su.write()
>>> ();su.run(timeend=6, read_result=True, resultdir=su.datadir);()  # +doctest: ELIPSIS
(...)

>>> mu.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Methane", "fixed"), ["t0001", "t0002", "t0003", "t0004", "t0005", "t0006"]] = [0, -10, -20, -30, -40, -50]
>>> mu.write()
>>> ();mu.run(timeend=6, read_result=True, resultdir=mu.datadir);()  # +doctest: ELIPSIS
(...)

```

```python
>>> su.result["converter_units"].loc[("R1_model", "2030", "Methanation", "2030", "total"), "value"]
0.5
>>> mu.result["converter_units"].loc[("R1_model", "2030", "Methanation", "2030", "total"), "value"]
50.0

```

We see that only half the number of units are created since our peak demand is only half as much as before.
On the model and on the result side, **there is really no difference** between the approaches shown in **this example**.
It is your choice to handle the input data in the appropriate or preferred way.

## Key Take Away

- REMix allows you to flexibly assign conversion rates and unit sizes.
- It is the modeler's task to size a unit and its conversion coefficients in the way best suited for his or her data available.
- In **linear models** the specific size of a unit can be chosen arbitrarily, if cost and coefficients are rescaled appropriately.
- In **mixed integer models** this is not the case. You can learn more about it {ref}`in this basic example <explanations_mip_label>`.
