---
title: Converters
lang: en-US
---

(explanations_converters_label)=

# Converters

In this section, you will get a general overview on the converter component of REMix. Some example applications of
converters are presented in {ref}`this section <remix_technology_examples>`. A complete overview on all model inputs
available to characterize a converter and their description can be found in the respective section of the
{ref}`technical documentation <remix_model_core_converter_label>`.

## General concept

Converters are components, that can convert commodities into different commodities within a single region of the model.
The operation of the converters is restricted by the amount of units, their size, activity limits and additional
optional features, e.g. unit commitment.
It is possible to create models with a fixed amount of converter units, and/or unit expansion can be modeled.
For example, converters are used to model power plants, photovoltaics or manufacturing processes etc.

## Converter activity

In the figure below a coal power plant is shown. The power plant consumes coal and produces electricity and emits
CO<sub>2</sub>.

```{figure} /img/REMix_CoalPP.svg
:align: center
Figure: Example converter coal power plant.
```

### Activity concept

The process of operating this power plant is defined by an `activity` in REMix.
Activities consume any amount of input commodities to generate any amount of output commodities.
For the example that means, that the electricity generation activity has to link the generation of electricity to the
emission of CO<sub>2</sub> and the consumption of coal. This is done by a respective `coefficient` {math}`c_t`
for every input commodity and every output commodity linked to this `activity`.

```{note}
It is possible to define a converter activity with no input commodity or no output commodity.
In that case, such an activity represent the behavior of source or a sink. An example for this usage can be found in
the basic examples section on {ref}`renewable energy modeling <vre_intro>`.
```

The `coefficient` specifies the amount of a commodity consumed or generated per `unit` {math}`u` of the converter.
Within a model region a converter can have multiple units, which share the same specifications.
In our example that means, that the `coefficient` has to describe the per-unit consumption of coal and the per-unit
generation of electricity and CO<sub>2</sub> at nominal load.

REMix represents the load of a converter by the `activity` value {math}`a_t`, which ranges between 0 and 1.
That means, that the total consumption or generation of a commodity {math}`g_t` can be calculated by the following
expression:

```{math}
    g_t = a_t \cdot c_t \cdot u
```

```{tip}
If a commodity is consumed or generated is decided by the coefficient's sign. A negative value indicates consumption and
a positive value indicates generation, also see {ref}`commodity <explanations_commodities_label>`.
```

### Using multiple activities

In general, converters can have multiple activities at the same time.
Each activity has the same specifications as described in the previous section.

```{tip}
For example, curtailment of renewable energy generation can be penalized with this modeling concept or you can model
combined heat and power plants with non-fixed ratio of heat and power output in a linear program.
An example implementation of curtailment penalization can be found in the respective
{ref}`VRE basic example <vre_curtailment_modeling_label>`.
```

```{math}
\sum_{i=1}^n a_t \leq 1
```

(explanations_mip_label)=

## Mixed-integer programming

Converters can also be modelled as MIP units. This is shown in a basic example:
{ref}`Mixed-Integer Energy Conversion Units <mip_example_label>`.
