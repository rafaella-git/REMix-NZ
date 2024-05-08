---
title: "Electric boilers"
lang: en-US
---

# Electric boilers

## General

An electric boiler is generally used to convert electricity into usable heat (hot water or steam) for buildings,
district heating networks and/or industry by means of direct electric heating and thus with an efficiency of lower than one.

## Modeling in REMix

In REMix, electric boilers are represented by the
{ref}`converter module <remix_model_core_converter_label>` using electricity as input commodity and heat of a certain
temperature as output commodity. In doing so, electric boilers are generally characterized by their efficiency, which links
the power input to the produced heat. Additionally, in the core_accounting module, the costs of the electric boilers (OMfix, OMvar,
annuity) can be calculated based on specific investment costs as well as fixed and variable operational costs. As for
all technologies, the boiler capacities may include either or both exogenous and endogenous share. The overall capacity,
possibly reduced by a constant or time-variable availability and divided by the corresponding efficiency defines the possible
heat output in each time step.

Advanced features, such as the discrete expansion of the electric boiler size and unit commitment can be added as well.
For further information see {ref}`MIP features <explanations_mip_label>`.

In the core_accounting module, the costs of the electric boilers are
accounted for. The costs can include investment costs, fixed operational costs and variable operational costs.

The {ref}`core_balance module <remix_model_core_balance_label>` balances the commodities involved in each time step.


## Data

In this part, the data that is required for modeling electric boilers is listed and explained for the use-case of a
large electric boiler providing high-temperature (HT) heat. Further optional features can be found in
{ref}`remix model reference <remix_model_label>`.

### Source/Sink

In order for the electric boilers to run, the input and output commodities need to be available.

**set_commodities.dat**
Example: Electricity, HTHeat

### Converter

The converter represents the electric boiler. It converts electricity to high-temperature heat.

**set_converter_techs.dat**
Example: ElectricBoiler_HT_L


**set_activities.dat**
Example: HeatGen (mode in which a converter produces heat)

**converter_coefficient.dat**

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for the electricity in the “HeatGen” mode to -1.01 and the coefficient of electricity to one to indicate an efficiency of 99%. [!!! is this correct???] |

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2 to indicate that two electric boiler units are already existing. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five electric boiler units can be available in this region. This means that three units can be optimized. |

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityUpperLimit | Set to 0.98 to indicate an availability factor of 98% for the electric boiler. |
| lifeTime | Set to 20 to indicate a technical life time of 20 years. |

### Accounting

In the accounting module, the costs associated with expansion and operation of the electric boilers are calculated. For
this, the cost components have to be defined in the respective indicator set.

**set_indicators.dat**
Examples: Invest, OMFix, OMVar

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set OMVar per activity for the electric boiler to 0.001 to indicate the cost of using the boiler, i.e. 0.001M€/GWh. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 250 for “Invest” indicator to specify investment costs 250M€/GW.|
| useAnnuity | Set to one for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 20 for "Invest" indicator to use an amortization time of 20 years in the annuity calculation. This is not relevant for the technical lifetime specified in converter_tech_parameter |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 5 for "OMFix" indicator to indicate OMFix costs of 5 M€/GW rated electric power. [!!! is this correct???]|

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | You can link the different cost components, such as OMVar, and OMFix, to the total system costs by setting a value of 1 for each. |
