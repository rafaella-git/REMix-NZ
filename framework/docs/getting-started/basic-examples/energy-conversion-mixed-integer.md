---
title: "Mixed-Integer Energy Conversion Units"
lang: en-US
---

(mip_example_label)=

# Mixed-Integer Energy Conversion Units

With mixed-integer programming (MIP) it is possible to allow for discrete decisions, such as on/off on the operational state of a plant, or investment into an integer number of units.
On this page some of the MIP features of REMix in terms of converters is shown using a combined cycle gas turbine plant as example as shown in the figure below.

```{figure} /img/REMix_CCGT.svg
:align: center
Figure: Combined Cycle Power Plant schema used in the examples below.
```

The power plant has a nominal production of 100 units, and - in case we model minimum part load - 30 units of electricity generation respectively.
A demand time series ranging from 0 to 100 units of electricity is employed.
The required natural gas is a result of the efficiency specified.
At nominal operation that value should be 60 % and, if modeled, the minimum efficiency in this example will be at 40 %.

```{note}
The general setup for all applications outlined on this page is identical.
Based on that, each section describes the changes necessary to modify the basic setup.
```

## Setting up the REMix model

```python
>>> import numpy as np
>>> import pandas as pd
>>> from matplotlib import pyplot as plt
>>> from remix.framework.api.instance import Instance


>>> idx = pd.IndexSlice

>>> m = Instance()
>>> m.map.aggregatenodesmodel.loc[("R1_data", "R1_model"), ""] = ""
>>> m.set.yearssel = ["2030"]

```

We are setting up our objective to minimize `FuelCost` and the `Investment`.
The investment is only a dummy value in these examples, to show the differences between mixed-integer investment and linear investment in REMix.
A demand time series is provided ranging from 0 to 100 MW with a fixed profile.
Additionally, a source for slack electricity is introduced, which will be used in case the plant cannot be operated due to minimum output power restrictions.
The slack source needs to have a lower profile (defaults to 0) to prevent exporting electricity to slack.
Finally, marginal cost for slack (very high) and import of CH<sub>4</sub> are assigned and a single unit of the combined cycle power plant is set up in order to see how the part load can be modeled.

```python
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "FuelCost", "global", "2030"), "perIndicator"] = 1
>>> m.parameter.accounting_perindicator.loc[("SystemCost", "Investment", "global", "2030"), "perIndicator"] = 1
>>> m.parameter.accounting_indicatorbounds.loc[("global", "horizon", "SystemCost"), "obj"] = -1

>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0150"] = -np.arange(1, 151)
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Demand", "Electricity"), "usesFixedProfile"] = 1

>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "FuelImport", "CH4"), "usesLowerProfile"] = 1
>>> m.parameter.sourcesink_config.loc[("R1_data", "2030", "Slack", "Electricity"), "usesLowerProfile"] = 1

>>> m.parameter.accounting_sourcesinkflow.loc[("FuelCost", "global", "2030", "FuelImport", "CH4"), "perFlow"] = 1
>>> m.parameter.accounting_sourcesinkflow.loc[("FuelCost", "global", "2030", "Slack", "Electricity"), "perFlow"] = 100000

>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "lifeTime"] = 1
>>> m.parameter.accounting_converterunits.loc[("Investment", "R1_data", "CCGT", "2030"), "perUnitBuild"] = 10

```

Furthermore, we define values for the efficiency of the power plant as well as the nominal generation and minimal generation.

```python
>>> nominal_generation = 100
>>> minimum_generation = 30
>>> maximum_efficiency = 0.6
>>> minimum_efficiency = 0.4

```

## Discrete Investment Decisions

The discrete expansion feature can be used for converter and transfer technologies.
In this example we show how to use it for energy conversion units.
The expansion of a technology can only be realized unit-wise in each region if this feature is activated.
The total installed capacity in a region therefore consists of a discrete number of units, each with a predefined capacity.

```{figure} /img/REMix_MIP_invest.svg
:align: center
Figure: Investment blocks using Mixed-Integer-Programming.
```

### Example

In the capacity expansion part, we build the linear approach first and the MIP approach second.
The solver calculates the number of units to meet the electricity demand at lowest cost in both approaches.
This results in the minimum amount of units possible being built.
The operation of the plant does not consider minimum part load, variable efficiency or other operational restrictions.

#### Linear Capacity Expansion

We define the plant's efficiency and unit size by providing the respective coefficients, add an upper limit for the number of units and run the model without further changes.

```python
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = -(nominal_generation / maximum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "coefficient"] = nominal_generation
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "CCGT"), "unitsUpperLimit"] = 3
>>> m.datadir = "linear-expansion"
>>> m.write(fileformat="dat")
>>> ();m.run(read_result=True, timeend=150, crossover=1);()  # +doctest: ELIPSIS
(...)

```

We can see that the number of units is a floating point number since we have linear expansion.

```python
>>> m.result["converter_units"].loc[("R1_model", "2030", "CCGT", "2030", "build"), "value"]
1.5

```

#### MIP Capacity Expansion

To activate the MIP expansion feature, we have to set `mipUnits` to `True`.
We run the model again and can see that since the maximum demand (150 units of electricity) cannot be met by a single unit a second full unit needs to be built.

```python
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "mipUnits"] = True
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "mipDispatch"] = True
>>> m.datadir = "mip-expansion"
>>> m.write(fileformat="dat")
>>> ();m.run(read_result=True, timeend=150, crossover=1);()  # +doctest: ELIPSIS
(...)
>>> m.result["converter_units"].loc[("R1_model", "2030", "CCGT", "2030", "build"), "value"]
2.0

```

## Modeling Part Load Operation

In this section, we have a look at the part load efficiency modeling.

With the `mipDispatch` approach it is possible to model minimum and maximum relative loads for converters.
In general, different activities can be defined for converters, each with an individual efficiency, minimum and maximum load.
The `coefficient`, `minLoad`, `maxLoad` and `constant` for each technology and activity can be set in `converter_coefficent`.
In the following example, we show three different approaches using this information to model part load operation:

1. Use `minLoad` to model part load with a constant efficiency.
2. Use `minLoad` with multiple activities to model part load with piece wise constant efficiency.
3. Use `constant` which applies an offset to generation or consumption and results in nonlinear efficiency.

Additionally, we compare the outcome of the MIP models with a full linear model.

### Examples

After the general setup, we first create the full linear model and then proceed with modeling operation considering minimum part load and

- constant efficiency,
- stepwise constant efficiency as well as
- two approaches for variable efficiency.

In the end we show a plot comparing the different approaches.
For that, the following code snippet provides functions to

- run the model and extract the results,
- sets up the plots with some annotations as well as
- a section that initializes the plot.

```python
>>> def plot(ax, col, x, y):
...     x_max = x.max()
...     x_min = x.min()
...     x_min_nonzero = x[x > 0].min()
...     y_max = y.max()
...     y_min = y.min()
...     y_min_nonzero = y[y > 0].min()
...
...     ax[0, col].plot([x_min, x_max], [y_max, y_max], "--", color="tab:red")
...     ax[0, col].plot([x_min, x_max], [y_min_nonzero, y_min_nonzero], "--", color="tab:red")
...
...     ax[0, col].plot([x_min_nonzero, x_min_nonzero], [y_min, y_max], "--", color="tab:orange")
...     ax[0, col].plot([x_max, x_max], [y_min, y_max], "--", color="tab:orange")
...
...     ax[0, col].scatter(x, y)
...
...     ax[0, col].annotate(
...         f"min consumption = {round(x_min_nonzero)}", xy=(0, 0),
...         xytext=(x_min_nonzero * 0.99, y_min_nonzero + (y_max - y_min_nonzero) / 2), textcoords='data',
...         horizontalalignment='right', verticalalignment='center', rotation=90, fontsize=8
...     )
...
...     ax[0, col].annotate(
...         f"max consumption = {round(x_max)}", xy=(0, 0),
...         xytext=(x_max * 1.01, y_min_nonzero + (y_max - y_min_nonzero) / 2), textcoords='data',
...         horizontalalignment='left', verticalalignment='center', rotation=90, fontsize=8
...     )
...
...     y = 100 * y / x
...     y[np.isnan(y)] = 0
...     y_max = y.max()
...     y_min_nonzero = y[y > 0].min()
...
...     ax[1, col].plot([x_min, x_max], [y_max, y_max], "--", color="tab:red")
...     ax[1, col].plot([x_min, x_max], [y_min_nonzero, y_min_nonzero], "--", color="tab:red")
...
...     ax[1, col].plot([x_min_nonzero, x_min_nonzero], [y_min, y_max], "--", color="tab:orange")
...     ax[1, col].plot([x_max, x_max], [y_min, y_max], "--", color="tab:orange")
...
...     ax[1, col].scatter(x, y)
...

>>> def run_model_and_get_result(m):
...     m.datadir = "mip"
...     m.write(fileformat="dat")
...     m.run(read_result=True, timeend=100, crossover=1)
...
...     x = m.result["commodity_balance"].loc[idx[:, "R1_model", "2030", "CCGT", "CH4"], "value"].abs().values
...     y = m.result["commodity_balance"].loc[idx[:, "R1_model", "2030", "CCGT", "Electricity"], "value"].values
...
...     return x, y
...

```

```python
>>> fig, ax = plt.subplots(2, 5, sharex=True, sharey="row", figsize=(15, 7.5), facecolor="#fffaeb")
>>> from matplotlib.patches import FancyBboxPatch
>>> import matplotlib

>>> ();matplotlib.artist.getp(fig.patch);()  # +doctest: ELIPSIS
(...)

>>> bb = fig.patch.get_bbox()
>>> p_bbox = FancyBboxPatch(
...     (bb.xmin, bb.ymin),
...     abs(bb.width), abs(bb.height),
...     boxstyle="round, pad=-0.0040, rounding_size=0.015",
...     ec="none", fc='#fffaeb',
...     mutation_aspect=4
... )
>>> p_bbox.set_figure(fig)

>>> _ = [(ax.grid(), ax.set_axisbelow(True), ax.set_facecolor('#fffaeb')) for ax in ax.flatten()]
>>> _ = ax[1, 0].set_ylabel("Conversion efficiency in %")
>>> _ = ax[0, 0].set_ylabel("Generation of Electricity")
>>> _ = ax[0, 0].set_title("linear model")
>>> _ = ax[0, 1].set_title("constant efficiency")
>>> _ = ax[0, 2].set_title("stepwise constant efficiency")
>>> _ = ax[0, 3].set_title("high full-load efficiency")
>>> _ = ax[0, 4].set_title("high part-load efficiency")

>>> _ = fig.text(0.5, 0.0, 'Consumption of CH4', ha='center')

```

We also prepare the model to have a fixed amount of units, i.e. a single unit by changing the parameters provided in the expansion section.

```python
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "mipUnits"] = False
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "CCGT"), "noExpansion"] = True
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "CCGT"), "unitsBuild"] = 1

```

#### Linear model

For the linear model we simply specify the coefficients for the `PowerGeneration` activity and the respective commodities.
The efficiency is constant and calculated according to the following equation:

```{math}
\eta = \frac{\text{coefficient}_\text{input}}{|\text{coefficient}_\text{output}|}
```

```python
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = -(nominal_generation / maximum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "coefficient"] = nominal_generation
>>> ();x, y = run_model_and_get_result(m);()  # +doctest: ELIPSIS
(...)
>>> _ = ax[0, 0].scatter(x, y)
>>> y = 100 * y / x
>>> _ = ax[1, 0].scatter(x, y)

```

```{note}
For the linear model, we do not add any annotations to the plot since the minimum values are 0 on both axes.
```

#### Constant efficiency with minimum load

Next, we introduce a minimum load (and constant efficiency) using the `mipDispatch` parameter and adding minimum load value, which is relative to the nominal load:

```python
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "mipDispatch"] = True
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "minLoad"] = (minimum_generation / nominal_generation)
>>> ();x, y = run_model_and_get_result(m);()  # +doctest: ELIPSIS
(...)
>>> plot(ax, 1, x, y)

```

#### Piece wise constant efficiency

Furthermore, with a MIP model it is possible to subdivide a converter into smaller units, each corresponding to a part-load range using the above approach.
For example, we distribute the loads on three different activities:

1. `PowerGeneration0` with an efficiency of 50 % ranging from 30 % to 53.3 % load (relative to electricity conversion)
2. `PowerGeneration1` with an efficiency of 58.3 % ranging from 53.3 % to 76.7 % load
3. `PowerGeneration2` with an efficiency of 66.7 % ranging from 76.7 % to 100 % load

Each of the three units has a constant efficiency as defined in the previous section.
Due to the differences in efficiency you have to either accept an overlap in the consumption curve or jumps in the generation curve.
A continuous function cannot be modeled with this concept.

```python
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = np.nan
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "coefficient"] = np.nan
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "minLoad"] = np.nan

>>> num_sections = 3
>>> load_sections = np.linspace((minimum_generation / nominal_generation), 1, num_sections + 1)
>>> efficiency_by_section = np.linspace(minimum_efficiency, maximum_efficiency, num_sections)

>>> for i in range(num_sections):
...     m.parameter.converter_coefficient.loc[("CCGT", "2030", f"PowerGeneration{i}", "CH4"), "coefficient"] = -(nominal_generation * load_sections[i + 1]) / efficiency_by_section[i]
...     m.parameter.converter_coefficient.loc[("CCGT", "2030", f"PowerGeneration{i}", "Electricity"), "coefficient"] = nominal_generation * load_sections[i + 1]
...     m.parameter.converter_coefficient.loc[("CCGT", "2030", f"PowerGeneration{i}", "Electricity"), "minLoad"] = load_sections[i] / load_sections[i + 1]

>>> ();x, y = run_model_and_get_result(m);()  # +doctest: ELIPSIS
(...)
>>> plot(ax, 2, x, y)

```

```{attention}
Remember to also set the OMVar for each activity in `accounting_converterActivity` if it should be higher than 0.
```

In the default partial-load setup the number of activities that are used by all units of a power-plant type in a region has to be equal to the number of online units.
If the detailed partial-load setup is activated (set `mipDetailedPartialLoad` to 1 in {ref}`converter_techparam <table_converter_techparam>`) the number of activities can exceed the number of online units.
This is useful for modeling processes with more than one product such as combined heat and power, since a unit can operate in several activities here.
However, beware that this might not make sense for single product processes, as the unit may operate in part load and full load at the same time.
Also, this setup increases the calculation time significantly.

#### Offset-converter with nonlinear efficiency

The last two variants use single activities again but instead the `constant` parameter is introduced to model an offset.
The two subplots to the right show the result in the third approach:
Here the minimum load point is defined with the `constant` parameter, which applies the desired efficiency at the minimum load point.
The `coefficient` parameter then determines the slope of the line in the upper plot to match the maximum load restriction.
The load-offset is the change of load relative to the constant load point.

```{math}
\eta = \frac{\text{constant}_\text{input} + \text{coefficient}_\text{input} \cdot \text{loadoffset}}
{|\text{constant}_\text{output} + \text{coefficient}_\text{output} \cdot \text{loadoffset}|}
```

First we start with a minimum efficiency of `40 %` and efficiency at nominal generation of `60 %`.
To model this, you define the minimum part load efficiency first with the constant parameter.
Then you calculate the value for the coefficient in a way that the sum of the electricity generation is equal to the desired total electricity generation at full load and the consumption of CH<sub>4</sub> for that point to match the desired full load efficiency as shown in the equation above.
We unset the specifications from the section before and add the respective data:

```python
>>> for i in range(num_sections):
...     m.parameter.converter_coefficient.drop([("CCGT", "2030", f"PowerGeneration{i}", "CH4"), ("CCGT", "2030", f"PowerGeneration{i}", "Electricity")], inplace=True)
>>> m.dataframes["sets"]["set_activities"] = []
>>> m.set.activities = []

>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = -(nominal_generation / maximum_efficiency - minimum_generation / minimum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "coefficient"] = nominal_generation - minimum_generation
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "constant"] = -(minimum_generation / minimum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "constant"] = minimum_generation

>>> ();x, y = run_model_and_get_result(m);()  # +doctest: ELIPSIS
(...)
>>> plot(ax, 3, x, y)

```

In the last approach we can reverse the efficiency curve (as observed in electrolyzers) by changing the `CH4` consumption coefficient and constant parameters.
The electricity generation will stay the same for the minimum part load and nominal generation, but the efficiencies (therefore the consumption of CH<sub>4</sub>) are changed respectively.

```python
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = -(nominal_generation / minimum_efficiency - minimum_generation / maximum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "constant"] = -(minimum_generation / maximum_efficiency)

>>> ();x, y = run_model_and_get_result(m);()  # +doctest: ELIPSIS
(...)
>>> plot(ax, 4, x, y)

>>> plt.tight_layout()
>>> fig.savefig("MIP-dispatch-efficiencies.svg")

```

### Overview

The figure below illustrates the behavior of the converter with the part-load modeling approaches applied.

```{figure} /img/REMix_MIP_partload_efficiencies.svg
:align: center
Figure: Efficiency model for part-load operation with Mixed-Integer Programming.
```

```{note}
Due to the different efficiencies given a constant operational range with respect to the generation, the consumption value range changes in the different setups.
The first three setups do have identical full load efficiency (lower image), thus their respective full load operation points are identical.
Due to the non-constant efficiency introduced from the piecewise and the offset converter the consumption at minimum load is higher in the second and third setup than in the first.
For the first and the last setup the minimum load efficiency is identical.
The full load operation in the third setup requires higher energy input due to the lower efficiency at full load and therefore higher consumption is observed at full load.
```

Of course you can combine the methods mentioned above.
Beware, however, that this can have unexpected effects, especially when modeling non-linear efficiency.

## Unit Uptime and Downtime

Additionally, it is possible to track and/or constrain the number of unit startups, the minimum uptime and minimum downtime.
By setting a minimum uptime or a minimum downtime, a unit that is switched on or off is forced to stay in this state for a minimum time.
The minimum uptime and downtime in hours can be set in `converter_techParam` as minUptime and minDowntime for each technology.

```{note}
Minimum uptime and downtime constraints only work with a minimum part load larger than zero.
```

### Examples

The following examples show how to implement the minimum uptime and minimum downtime feature for a converter.
To keep things simple, we will use the constant efficiency approach from the previous section together with 3 units.
For that, we first reset the parameters, that had been used in the latest example and adjust them in the respective way.
We also change the number of available units to three.
The minimum uptime and minimum downtime are selected as three as well.

```python
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "minLoad"] = (minimum_generation / nominal_generation)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "coefficient"] = -(nominal_generation / maximum_efficiency)
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "coefficient"] = nominal_generation
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "CH4"), "constant"] = np.nan
>>> m.parameter.converter_coefficient.loc[("CCGT", "2030", "PowerGeneration", "Electricity"), "constant"] = np.nan
>>> m.parameter.converter_capacityparam.loc[("R1_data", "2030", "CCGT"), "unitsBuild"] = 3
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "minUptime"] = 3
>>> m.parameter.converter_techparam.loc[("CCGT", "2030"), "minDowntime"] = 3

```

#### No operational violations

First, we can construct a case where no operational violations appear by modifying the demand time series:

- `t0004`: only a single unit can be operational
- `t0005`: at least two must be operational
- `t0006`: all units must be operational
- `t0007`: at least one unit must be operational
- `t0008` - `t0009`: all units need to shut down
- `t0010` - `t0012`: only a single unit can be operational

```python
>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"] = [0, 0, 0, -50, -120, -300, -100, 0, 0, -30, -30, -30, 0, 0, 0]

```

```{attention}
Beware that independently of how many time steps you are running your model on, the **time dimension is always circular!**
That means that uptime and downtime extend over the end of the time interval selected for optimization to the beginning of the time interval.
Especially when working with the minimum downtime and uptime constraints this may have unexpected effects.
```

To extract the amount of units that are operational per time step from the optimization we have to provide the keyword `gdx_mipconverter` to the solver.

```python
>>> m.datadir = "minimum-uptime"
>>> m.write(fileformat="dat")
>>> ();m.run(read_result=True, gdx_mipconverter=1, timeend=15, crossover=1);()  # +doctest: ELIPSIS
(...)

```

We can make a plot of the unit activity, where we plot the number of online units on the left y-axis and the demand on the right y-axis.

```python
>>> fig, ax = plt.subplots(3, sharex=True, sharey=True, facecolor="#fffaeb")
>>> axtwins = np.array([_.twinx() for _ in ax])
>>> _ = fig.text(0.02, 0.5, 'Online units', va='center', rotation='vertical')
>>> _ = fig.text(0.97, 0.5, 'Electricity demand', va='center', rotation='vertical')

```

We can extract the number of online units from the results and plot them:

```python
>>> online_units = m.result["converter_unitsOnline"].loc[m.result["converter_unitsOnline"].values == 1].index.get_level_values("level")
>>> demand = -m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"]

>>> l1 = ax[0].step(np.arange(1, 16), online_units, label="Online Units")
>>> l2 = axtwins[0].step(np.arange(1, 16), demand, "--", color="tab:orange", label="Demand")

>>> links = l1 + l2
>>> labels = [l.get_label() for l in links]
>>> _ = ax[0].legend(links, labels)
>>> _ = ax[0].set_ylabel("no violations")

```

#### Downtime violation

We modify the time series (so we can construct a case where a downtime violation appears) by adding additional demand on time step `t0009`:

```python
>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"] = [0, 0, 0, -50, -120, -300, -100, 0, -30, -30, -30, -30, 0, 0, 0]

```

Then we can rerun the model and add to our plot:

```python
>>> m.datadir = "minimum-uptime"
>>> m.write(fileformat="dat")
>>> ();m.run(read_result=True, gdx_mipconverter=1, timeend=15, crossover=1);()  # +doctest: ELIPSIS
(...)

>>> online_units = m.result["converter_unitsOnline"].loc[m.result["converter_unitsOnline"].values == 1].index.get_level_values("level")
>>> demand = -m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"]

>>> _ = ax[1].step(np.arange(1, 16), online_units)
>>> _ = axtwins[1].step(np.arange(1, 16), demand, "--", color="tab:orange")
>>> _ = ax[1].set_ylabel("downtime violation")

```

#### Uptime violation

Finally, we can add an uptime violation by changing the demand in a way that in time step `t0007` not all three units can be operational but only a single unit can be operational.
By changing this, it is impossible to operate two or more units in `t0005` and `t0006`, because those units can only be operated for two time steps.

```python
>>> m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"] = [0, 0, 0, -50, -120, -300, -50, 0, 0, -30, -30, -30, 0, 0, 0]

>>> m.datadir = "minimum-uptime"
>>> m.write(fileformat="dat")
>>> ();m.run(read_result=True, gdx_mipconverter=1, timeend=15, crossover=1);()  # +doctest: ELIPSIS
(...)

>>> online_units = m.result["converter_unitsOnline"].loc[m.result["converter_unitsOnline"].values == 1].index.get_level_values("level")
>>> demand = -m.profile.sourcesink_profile.loc[("R1_data", "2030", "Demand", "Electricity", "fixed"), "t0001":"t0015"]

>>> _ = ax[2].step(np.arange(1, 16), online_units)
>>> _ = axtwins[2].step(np.arange(1, 16), demand, "--", color="tab:orange")
>>> _ = ax[2].set_ylabel("uptime violation")
>>> _ = ax[2].set_xlabel("timestep")

>>> _ = [(_.set_facecolor('#fffaeb')) for _ in ax.flatten()]
>>> _ = [(_.set_facecolor('#fffaeb')) for _ in axtwins.flatten()]

>>> fig.savefig("MIP-uptime-downtime.svg")

```

```{figure} /img/REMix_MIP_uptime_downtime.svg
:align: center
Figure: Illustration of uptime and downtime violations.
```

## Key Take Away

- REMix allows modeling conversion with mixed-integer (MIP) approaches.
- You can optimize the investment into a number of units in discrete steps in MIP models, but you cannot optimize the size of each unit.
- Different approaches to part-load modeling are available for converter units.
- In case multiple units of a technology are available, the MIP approach allows determining how many integer units of a technology are operational per time step and for each unit at which part load.
