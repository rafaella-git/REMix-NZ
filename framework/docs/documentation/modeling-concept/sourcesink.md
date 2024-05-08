---
title: Sources and Sinks
lang: en-US
---

# Sources and Sinks

## General concept

Sources and Sinks are used to model commodity flow crossing the system's boundary.
Sources can be used to provide a commodity in a specific or all considered regions.
Sinks consume commodities.
The available and consumable amount of commodities by sources and sinks can be
indicated or limited by a profile or annual sum.

For example, exogenous demand like household heating or industrial gas demand or
supply of energy like imports from non-modeled regions can be modeled
with this component.

A source or a sink is always linked to a commodity within the nodal balance as
shown in the figure below.
It cannot be directly associated with a converter, storage or transfer link.

```{figure} /img/REMix_sources-sinks.svg
:align: center
Figure: Illustration of a source and a sink with a single commodity in REMix.
```

## Distinction from converters

In contrast to converters, sources and sinks will always have one single flow
associated with a single component.
If a component acts as a source or as a sink is defined by the data and through
the commodity balance {ref}`sign definition <explanations_commodities_label>`:
a flow from a source into the modeled system has a positive value, a flow
leaving the modeled system has a negative value.

The value of the flow can have upper and lower limits, follow a predefined fixed
time series or simply have a constant value.
It is also possible to limit the total generation or consumption of the source
or the sink by constraining the annual flow.
The capacity of the component cannot be expanded as it is the case for
converters.
Therefore, the most common use case for sources and sinks in your model are
usually

- demand modeling, where the demand follows a fixed time series or has a constant
  value and
- import and export modeling, where a commodity can be imported or exported into
  the system at a limited rate.

For example, when modeling renewable electricity generation with photovoltaic
plants you do not necessarily need to use a source component.
Different approaches possible (not an exhaustive list) are shown in the figure
below:

- a) Direct electricity generation by the source, e.g. by following an upper
  generation limit profile. Since a single source is used, the generation
  capacity is fixed exogenously.
  This is therefore the simplest approach.
- b) Production of a solar radiation commodity by a source and using a converter
  to generate electricity with the solar radiation input.
  With this, the PV plant capacity can be calculated endogenously.
- c) Direct generation of electricity using a converter, that does not take any
  inputs but is restricted in its activity profile.
  This also offers the option to calculate the PV plant capacity as a result of
  the optimization but is simpler than variant b) regarding mathematical
  complexity.

That means, depending on your modeling goal, you need to choose the approach.

```{figure} /img/REMix_sources-sinks-renewable-options.svg
:align: center
Figure: Illustration of a Source and a Sink with a single commodity in REMix.
```

## General concept

The goal of an energy system model is to meet an externally specified demand for
one or more commodities based on certain criteria such as minimum cost.
Such a demand is generally referred to as a sink.
In order to meet this demand, commodities are produced, converted, stored and
imported.
The source of a commodity, in the sense of its generation, can be either inside
or outside the model scope.
If the generation is inside the system boundary, then the generation is
explicitly modeled, e.g. electricity from a wind turbine, if it is outside, then
the commodity flows are treated as imports.

```{figure} /img/SourceSink_concept.svg
:align: center
Figure: Schema of a simple ESM.
```

## Mathematical formulation

Sinks and sources can either be specified by an hourly flow or limited in their
annual sum:

```{math}
\text{SourceSink} = \sum_{timesteps}(\text{SourceSinkFlows}_\text{t}).
```
