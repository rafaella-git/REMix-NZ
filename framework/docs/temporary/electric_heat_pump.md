---
title: "Electric heat pumps"
lang: en-US
---

# Electric heat pump technologies

## General

Like a refrigerator or an air conditioner, a heat pump is generally a machine that uses electrical energy to extract
heat from one reservoir and add heat to another. In the energy system, electric heat pumps are used to convert
electricity and ambient heat (e.g. from the air or ground) into usable heat (hot water or steam) for buildings,
district heating networks and/or industry by means of a refrigeration cycle and thus with a coefficient of performance
(COP) of usually higher than one. So far, heat pumps are used mainly in buildings and for temperatures below 100 degree
Celsius. This is increasingly being supplemented by the use of larger heat pumps, e.g. in heat networks.
High-temperature heat pumps could also be used in the future, particularly for generating process heat, but also for
storing energy.

## Modeling in REMix

In REMix, electric heat pumps are represented by the {ref}`converter module <remix_model_core_converter_label>` using
electricity as input commodity and heat of a certain temperature as output commodity. In doing so, heat pumps are
generally characterized by their COP, which links the power input to the produced heat. Additionally, in the
core_accounting module, the costs of the heat pumps (OMfix, OMvar, annuity)
can be calculated based on specific investment costs as well as fixed and variable operational costs.

As for all technologies, the heat pump capacities may include either or both exogenous and endogenous share. The
overall capacity, possibly reduced by a constant or time-variable availability and divided by the corresponding COP
defines the possible heat output in each time step. Heat pump COPs are strongly depending on the temperature difference
between heat source and heat sink. Given the considerable seasonal variations in ambient temperature, this is
particularly important in the case of air source heat pumps. For this reason, provision of regional and hourly COP
values is possible and recommended. Alternatively, a constant COP can be applied to HP technologies relying on heat
sources featuring only minor fluctuations in temperature. The
{ref}`core_balance module <remix_model_core_balance_label>` balances the commodities involved in each time step.

Advanced features, such as the discrete expansion of heat pump unit size and unit commitment can be added as well. For
further information see {ref}`MIP features <explanations_mip_label>`.

## Data

In this part, the data that is required for modeling heat pumps is listed and explained.

### Source/Sink

In order for the heat pumps to run, the input and output commodities need to be available.

**set_commodities.dat**
Example: Electricity, Hydrogen, Heat or respective aliases such as Elec and H2

### Converter

The converter represents the heat pump, e.g. an air source heat pump unit with all its ancillary equipment. It is
modeled as a black-box converting electricity to heat.

**set_converter_techs.dat**
Example: HeatPump_Air2Water


**set_activities.dat**
Example: HeatGen (mode in which a converter produces heat)

**converter_coefficient.dat**

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for electricity in the “HeatGen” mode to 1 and the coeffient of heat to 3.5 to indicate a COP of 3.5. [!!! is this correct???] |

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2.000 to indicate that two thousand heat pump units are already existing. |
| unitsUpperLimit | Set to 5.000 to indicate that a maximum of five thousand heat pump units can be available in this region. This means that three thousand units can be optimized. |

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityUpperLimit | Set to 0.98 to indicate an availability factor of 98% for the heat pump. |
| lifeTime | Set to 20 to indicate a technical life time of 20 years. |

### Accounting

In the accounting module, the costs associated with expansion and operation of the heat pumps are calculated. For this, the cost components have to be defined in the respective indicator set.

**set_indicators.dat**
Examples: Invest, OMFix, OMVar

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set variable operation and maintenance costs (OMVar) per activity for the heat pump to 0.007 to indicate the cost of its usage, i.e. 0.007 €/kWh. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 1500 for “Invest” indicator to specify investment costs 1500 €/kW rated electrical power.[!!! is this correct???]|
| useAnnuity | Set to one for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 20 for "Invest" indicator to use an amortization time of 20 years in the annuity calculation. This is not relevant for the technical lifetime specified in converter_tech_parameter |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 20 for “OMFix” indicator to indicate annual OMFix costs of 20 €/kW rated electric power. [!!! is this correct???]|

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | You can link the different cost components, such as OMVar, and OMFix, to the total system costs by setting a value of 1 for each. |
