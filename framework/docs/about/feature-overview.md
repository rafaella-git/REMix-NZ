---
title: Feature overview
lang: en-US
---

(about_feature-overview)=

# Feature overview

## Usable modeling features

REMix comes with several features to address a large variety of research questions in the field of energy systems analysis.
In general, we distinguish most of these features between planning of operation (dispatch optimization) and investments into new system components (expansion planning).
Usually, the optimization is done as a minimization of costs.
However, what to optimize is rather a matter of input data you provide to REMix.

In addition, there are technical features which may help you to performantly solve your model or to scale and modify a model according to your specific needs.
In the following, we also shortly introduce those features, which are described in more detail in {ref}`Modeling concept of REMix <modeling_concept_label>`.

### Dispatch optimization

This is a basic feature of the model framework which is frequently used.
It refers to the optimal operation of a model component (e.g., a power plant).
In REMix this is basically related to the activity of the converter and storage units which can operate between defined limits (e.g. the rated power of a power plant represents an upper activity limit).
The way of how the dispatch operation is done can be modified by several additional sub-features, which are introduced in the following:

#### Storage scheduling

```{figure} /img/Overview_dispatch_storage_converter.svg
:align: center
Figure: Dispatch of a PV and gas power plant together with storage scheduling.
```

The ability to shift power generation and demand in time is another model feature that allows load-balancing.
Therefore, storage reservoirs can be used to store any type of energy carrier for a certain time period.
In context of dispatch optimization this relates to finding the optimal points in time for charging and discharging energy storage.

#### Unit commitment

```{figure} /img/REMix_MIP_mintime.svg
:align: center
Figure: Dispatch optimization with minimum up- and down-times.
```

Especially the operation of large thermal power plants is subject to technical constraints such as ramping or minimum down-times (time duration for which a power plant must be turned off).
To perform a dispatch optimization with respect to these constraints, one needs to use integer variables.
The associated unit commitment problems are a very common problem class in the field of electrical power engineering and optimal scheduling of power plant portfolios of utilities.
Solving such mixed-integer linear optimization problems is particularly useful when distinct power plant units should be explicitly modeled.
Hence, if you have a world model where Europe is a single model region, unit commitment is just overkill.

### Energy transfer

Usually, models created with REMix are multi-regional.
Thus, you define so-called nodes (or model regions) which exchange energy or other commodities with each other (load-balancing in space) over so-called transfer links.
However, defining a number of nodes (>=1) is up to your choice.

#### Optimal power flow

The optimal power flow in particular refers to power transmission problems in electricity grids (>100 kV), where the dispatch of grid-connected power plants is optimized, e.g., to minimize power production costs or losses for power transmission.
As REMix is a linear model, this feature is only implemented with linear power flow constraints, which are usually referred to as DC-power flow constraints.

### Expansion planning

By expansion planning we usually refer to an extension of the capacities of converters, transfer and storage units.
In other words, in addition to dispatch optimization you allow the increase of the upper limits that constrain the converter activity, flow and reservoir size, respectively.
This increase is usually subject to costs, i.e. investment costs for building new or extending existing energy infrastructure.
Accordingly, the result of a classical expansion planning for converters would be a least-cost power plant portfolio.

#### Integer capacity expansion

```{figure} /img/REMix_MIP_invest.svg
:align: center
Figure: In integer capacity expansion the total capacity of a power plant consists of a discrete number of units of a technology-specific size.
```

Depending on your model's scope considering investments into plants of any size might not make sense.
It does if you model large areas, where the optimized capacities represent the aggregation of many commercially available power plant units.
It also holds for decentral power plants, such as photovoltaics as long as you do not model individual buildings.
However, the higher the spatial resolution of your model, the more plausible it is to only allow investment in plants with a fixed unit size (e.g., 5 MW for wind turbines or 50 MW for open-cycle gas turbines) using integer decision variables.

#### Target systems

A very common application of expansion planning is the identification of target systems.
Therefore, strategic or political targets have to be achieved (by setting the appropriate constraints) at a certain point in the future (e.g. 100 % renewable power supply).
Thus, when optimizing target systems you usually define the targets and perform expansion planning for a single year (e.g. 2050).
However, in this simple setup, the year may only imply the techno-economic assumptions made (e.g., which decrease in costs or increase in conversion efficiencies do you expect in 2050).

#### Brown-field and green-field approach

The terms brown-field and green-field refer to the model-exogenous (prescribed) capacities of the energy system. If these capacities are 0, we consider the expansion planning to start from the green-field.
Opposed to that, if these values are >0 (e.g. based on the installed capacities of the current energy system), we consider it as brown-field expansion planning.
However, for multi-regional models, the topology of the transport networks has to be defined in any case.
In other words, with regard to networks the expansion planning is always limited to prescribed edges and edges' candidates between nodes.

```{figure} /img/Overview_brown_green_field.svg
:align: center
Figure: In a brown-field optimization preexisting capacities are available and additional capacities can be built. In a green-field optimization no preexisting capacities are available.
```

#### System transformation (myopic or path optimization)

```{figure} /img/Overview_system_transformation.svg
:align: center
Figure: A target optimization optimizes one year only. With a myopic optimization several years are calculated one after another without the information of the next years. A path optimization calculates several target years at ones with perfect foresight.
```

You may also want to know at which point in the future to deploy or decommission a particular plant if you are not only interested in the technology mix or the design of a least-cost power plant portfolio.
Furthermore, expansion planning for target systems is only plausible if the time horizon until final implementation is comparably large (e.g. greater than the lifetime of power plants commissioned today).
Thus, the path-optimization feature allows you to find optimal transformation pathways, where you usually start from a brown field and consider multiple years for making investment decisions to achieve strategic or political targets at a certain point in the future.
However, if you prefer to determine consistent system transformations without a perfect foresight in terms of investment decisions, you can also use a myopic approach, where the investments of a previous year are taken as inputs to the expansion planning of a subsequent year.

### Indicator accounting and multi-criterial assessments

REMix comes with a powerful indicator-accounting system.
This means, as long as you have the data (e.g., factors provided per installed capacity) you can directly generate indicators you are interested in or constrain your model with respect to these indicators.
For example, you may use this feature to define a minimal self-sufficiency for each region of your model.
Or you are interested in the land use caused by cost-minimal expansion planning.
Then you only provide the specific land use per technology (e.g. in m² per GW) to observe this indicator in your final results.

## Technical features

### Adjustable model scopes and resolutions

#### Space

Especially if you are setting up models that represent networks, you may be interested in aggregating the nodes of the network (e.g., to state-level or to simply reduce the model's size).
Therefore, REMix comes with two levels of the spatial dimension.
The first (nodesData) describes the input data, the second (nodesModel) the desired resolution for which the model is solved and thus the output data.
The required spatial aggregation is done within REMix based on a mapping between nodesData and nodesModel to be provided by the modeler.
In addition, one can individually enable and disable the optimization for each node.

There is also the possibility to define accounting nodes (accNodes), e.g. if you model a country like Germany with several nodes, but want to give out the overall CO<sub>2</sub> emissions of the entire country as a result.

#### Time

The temporal dimension also has two levels. The first is the operational time horizon which comprises the 8760 hours of a single year.
The second temporal level describes the investment horizon which comprises the years subject to capacity expansion (and decommissioning).
Opposed to the operational time horizon the number of years considered is up to your choice.

### Methods for solving

#### High performance computing (HPC)

REMix comes with an interface to the open-source solver PIPS-IPM++, which is suited to solve energy system optimization models on massively parallel high-performance computers.
This enables you to solve really large-scale models not solvable otherwise.

#### User-defined objectives

The most common motivation for finding "good" operational and investment decisions is related to economic objectives, i.e. the minimization of costs.
However, REMix is flexible in this regard. Accordingly, you may consider a minimization of CO<sub>2</sub> emissions to be useful to better analyze climate-neutral energy systems.
Therefore, you simply set up the indicator accounting for CO<sub>2</sub> emissions and define the total emissions as objective (instead of costs).
Moreover, if you have data on how to weight different objectives among each other, you can also perform multi-objective optimizations.

#### Modeling to generate alternatives (MGA)

In combination with an indicator that constrains the maximal deviation from a cost-minimal system, this feature allows you to perform MGA analysis.

#### Pareto fronts

Another way to solve models with multiple objectives concerns pareto fronts.
This feature helps you to identify pareto-optimal energy systems with respect to two objectives.

## Features in development

### Rolling horizon dispatch

One strong assumption of optimization models is that they implicitly have the perspective of a central planner having perfect information and thus perfect foresight in terms of events to be expected within the optimization time horizon.
For gaining insights into an ideal energy system of the future this might be okay.
However, if your model shall rather serve as a real-world model simulation, you may want to avoid decisions based on perfect foresight.
Therefore, you may use the rolling horizon dispatch which limits the (still perfect) temporal foresight horizon and optimizes several overlapping time horizons.

### Disruptive events

To plan energy systems which are robust against disruptive events you may want to consider these events also in your model.
However, due to the perfect foresight (see above), disruptions cannot be simply added by appropriately changing the input parameters of a model, especially if expansion planning is enabled. Therefore, the resilience approach can be used where simply two subsequent model runs are performed.
The first as typical expansion planning, the second as pure dispatch optimization model which is extended with an additional input, a so-called functionality time series.
This allows you to simulate a global drop of the available capacities of the energy system's components (caused by any disruptive event) at a predefined point in time and to observe the system's reaction.
