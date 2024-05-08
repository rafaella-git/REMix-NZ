---
title: Modeling Energy Systems with REMix
lang: en-US
---

(modeling_concept_label)=

# Modeling Energy Systems with REMix

The figure below illustrates the modeling concept of REMix:

```{figure} /img/JOSS_overview_remix.svg
:align: center
Figure: Illustration of the structure of a REMix model.
```

(remix_concept_label)=

## Concept

REMix is built on data nodes, typically considered as geographical regions.
Different data nodes can be interconnected by transfer links, which allow
the transport of one or multiple commodities between data nodes, for example
of gas via pipelines or electricity via AC or DC links.
The demand and supply of commodities, the conversion of one commodity to a
different commodity as well as storage are defined within each data node.
REMix provides a very flexible tool for tracking indicators over the whole
system, e.g. costs or CO<sub>2</sub> emissions, i.e. the accounting module,
which is used for the definition of the objective.

On the temporal scale REMix is built on uniform hour-long time steps for each
model year. Investment and dispatch can be optimized for multiple model years,
either with a myopic approach or using pathway optimization.

The user has access to the following building blocks to create an energy system
model:

- `sourcesink` for exogenous demand or supply of commodities,
- `converter` for the conversion and endogenous supply of commodities,
- `storage` to store commodities over a period of time,
- `transfer` to transport commodities from one region to another.

The most important ideas of these building blocks as well as the commodity and
the indicator terms are listed below.
More detailed explanations are given on the respective subpages.

### Commodities

Commodities are energy carriers or other physical flows (e.g., fuels,
electricity, heat, etc.).
They need to have dedicated sources and sinks if any exogenous demand or supply
is given.
Inside the system boundaries of a model node they can be supplied from a source,
used to meet demand, stored and converted into other commodities.
Between nodes of the system, commodities are exchanged by transfer
technologies.

### Indicators

Indicators are used for accounting purposes (e.g., costs, firm capacities,
land use, CO<sub>2</sub>, etc.).
They typically account for non-physical or non-restricted flows and are
calculated based on associated units, links, activities, etc.
Indicators can be balanced across single model nodes, globally or according to
custom regions.
Similarly, they can account for individual years or planning horizons.

### Sources and sinks

Sources and sinks describe physical commodity flows across the system boundary.
They represent, for example, fuel imports from a global market or emissions into
the environment.
These flows can be connected with monetary costs in the form of indicators.
Additionally, it is possible to limit availability.

### Converters

Converters are technologies which allow conversion processes from one commodity
to another.
Conversion activities have to be assigned to converters.
These activities describe the ratios between input commodities and output
commodities.
A converter can also have multiple different activities (e.g. generation of
power and generation of heat).

### Storage

Storage technologies can store commodities and therefore introduce a temporal
shift between supply and generation or demand and consumption.

### Transfer

Transfer links allow the transport of commodities between different model nodes:
while conversion of commodities can only happen within the same respective node
of the model (e.g. from fuel to electricity at model node A), links
can ship fuel from model node B to A, which could then be converted to
electricity.

(remix_domains)=

## Model data

In REMix the energy system model is defined by its data.
The data structure follows the modeling approach described at the beginning of
this page, i.e. the data is defined based on technological, geographical and
temporal domains.
To simplify the modeling of large-scale energy systems, the data domains to
define specific components of the system and its properties are not identical,
but change based on the information you want to give.

```{tip}
For example, to model the production of electricity from wind, the modeler
typically provides the following data:

- the existing as well as the minimum and maximum electrical capacity
- time series with a capacity factor to model the availability of wind
- cost assumptions for investment and operation of the wind turbines

In REMix this can be done using a converter.
The information on the lifetime of an investment into a wind turbine is based on

- the technology class and
- the year of investment.

The information on how much of one commodity one unit of that technology
provides (or consumes) at nominal conditions is defined by

- the technology class,
- the year of investment,
- the name of the activity and
- the name of the commodity.

The number of units that is available or its lower and upper limits for unit
expansion is defined by

- the data node,
- the modeling year and
- the technology class.

Finally, the time series for the capacity factor is defined by

- the data node,
- the modeling year,
- the technology class and
- the profile type (lower or upper value or a fixed capacity factor).
```

The same concept applies to sources, sinks, storage and transfer in similar
ways.
In order to stick to the uniform data structure, the data always has to
be provided with the respective domains, even in case there is only a single
value in a specific domain, e.g. if you model only a single year.

An overview of the specific data domains is given in the respective subsections
{ref}`model documentation <remix_model_label>`.

```{toctree}
:glob:
:maxdepth: 1
:hidden:

commodities
indicators
converter
sourcesink
storage
transfer
```
