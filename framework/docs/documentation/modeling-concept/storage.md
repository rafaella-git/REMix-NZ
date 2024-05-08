---
title: Storage technologies
lang: en-US
---

# Storage technologies

## General concept

In an energy system, energy demand and energy generation should always match to avoid shortages and overproduction.
In order to introduce a temporal buffer into the energy balance, energy storages are used. A storage unit in energy system models
is an entity that receives a medium, holds it for a specified duration, and then releases it subsequently. The storage
can be charged in times of abundance and discharged when needed. The available amount of storage can be pre-existing or
additionally built. The storable amount of commodity depends on the size of the storage, storage level limits and
advanced features such as degradation. Several types of losses can occur: charging and discharging as well as
self-discharging losses. They can give an indication if a storage is rather used as short-, medium- or long term
flexibility option.


## Storages in linear ESM`s

Depending on the modeling method and application, storages can be modeled in any degree of complexity. To keep model sizes within reasonable bounds, linear energy system models generally do not differ in storage characteristics for different media (e.g. electricity, heat or gas). Thus, only the energy quantities in a storage can be accounted for, but not, for example, the temperature level of heat storages. An energy storage unit can be described by the three subcomponents: the charging unit, the storage reservoir and the discharging unit with individual parameters and variables.

```{figure} /img/Storage_concept.svg
:align: center
Figure: Schema of a storage unit in ESM´s.
```

## Mathematical formulation

The storage reservoir is defined by a maximum {math}`\text{SoC}_\text{max}` and
minimum {math}`\text{SoC}_\text{min}` filling level.
The reservoir level *SOC* (state of charge) at a certain time step *t* must be
within these bounds.

```{math}
\text{SoC}_\text{min} \leq \text{SoC}_\text{t} \leq \text{SoC}_\text{max}
```

The storage filling level per time step results from the storage filling level of
the previous time step {math}`\text{SoC}_\text{t-1}` increased by the energy
stored {math}`\text{P}_\text{input}` and reduced by the energy discharged
{math}`\text{P}_\text{output}`. In addition, the storage losses
{math}`\text{Loss}_\text{t}` must be taken into account.

```{math}
\text{SoC}_\text{t} = \text{SoC}_\text{t-1}+\text{P}_\text{input}-\text{P}_\text{output} - \text{Loss}_\text{t}
```

The losses are composed of charging, discharging and self-discharging losses.
The latter can be expressed as a percentage and/or as an absolute value.
In addition, it is possible to specify the self-discharge rate discretely
{math}`\text{selfdischarge}_\text{discrete}^\text{state}` for different state of
charges or to define it linearly {math}`\text{selfdischarge}_\text{linear}` as a
function of the state of charge.

```{math}
\begin{aligned}
\text{Loss}_\text{t} &= - \text{SoC}_\text{t-1} * \text{selfdischarge}_\text{linear} \\
                       &+ selfdischargeAbs \\
                       &+ \frac{\text{chargingloss}}{(1-\text{chargingloss})}*\text{P}_\text{input}\\
                       &+ \text{dischargingloss} * \text{P}_\text{output} \\
                       &- \text{size} * \sum_{states}(\text{unitsSOC}_\text{t-1}^\text{state}*\text{SOC}^\text{state}*\text{selfdischarge}_\text{discrete}^\text{state})
\end{aligned}
```
