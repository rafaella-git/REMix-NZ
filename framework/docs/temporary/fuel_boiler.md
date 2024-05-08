---
title: "Fuel boilers"
lang: en-US
---

# Fuel boilers

## General

A fuel boiler is generally used to convert a solid, liquid or gaseous fuel into usable heat (hot water or steam) for
buildings, district heating networks and/or industry. The used fuel can be of biogenic, fossil or synthetic nature.
Depending on the fuel used, CO2 emissions may occur.

## Modeling in REMix

In REMix, fuel boilers are represented by the {ref}`converter module <remix_model_core_converter_label>` using the
corresponding fuel, or also a variety of different fuels as input commodity and heat of a certain temperature as output
commodity. Depending on the fuel type considered, other output commodities such as emittants and water may be
considered. The efficiency links the input commodity to the produced heat, and may have deviating values for different
input fuels. As for all technologies, the boiler capacities may include either or both exogenous and endogenous share.
The overall capacity, possibly reduced by a constant or time-variable availability and divided by the corresponding
efficiency defines the possible heat output in each time step. Advanced features, such as the discrete capacity
expansion and unit commitment can be added to the fuel boilers as well. For further informations see
{ref}`MIP features <explanations_mip_label>`.

In the core_accounting module, the costs of the fuel boilers and the emission
of CO2 are accounted for. The costs can include investment costs, fixed operational costs and variable operational
costs.

The {ref}`core_balance module <remix_model_core_balance_label>` balances the commodities involved in each time step.

## Data

In this part, the data that is required for modeling fuel boilers is listed and explained for the use-case of a large
gas boiler providing low-temperature (LT) heat using natural gas. Further optional features can be found in
{ref}`remix model reference <remix_model_label>`.

### Source/Sink

In order to be able run the boiler, the fuel needs to be available. In our case, the fuel can be imported as natural
gas from outside the modeling scope by implementing a fuel source.

**set_commodities.dat**
Example: Electricity, LTHeat, NatGas

**set_sourcesink_techs.dat**
Example: FuelImport

**sourcesink_config.dat**

| sourcesink_config_parameter | Example |
| ------ | ------ |
| usesLowerProfile | Set to 1 to make sure that a fuel cannot be sold outside of the modeling scope, since the default value for the profile will then be zero (if no profile is further specified). |
| usesUpperSum |  Set to 1 to indicate the annual amount of fuel that is availabe. The amount can be specified in sourcesink_annualSum.dat. |

**sourcesink_annualSum.dat**

| sourcesink_annualSum_parameter | Example |
| ------ | ------ |
| upper | Set to inf to indicate that an unlimited annual amount of natural gas is available. |


### Converter

The converter represents the fuel boiler, so in our case a gas boiler. It converts the imported gas to low-temperature
heat.

**set_converter_techs.dat**
Example: NatGas_boilerLarge [!!!stimmt die Nomenklatur????]


**set_activities.dat**
Example: HeatGen (mode in which a converter produces heat)

**converter_coefficient.dat**

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for the natural gas in the “HeatGen” mode to -1.1 and the coefficient of electricity to one to indicate an efficiency of 90%. [!!! is this correct???] |

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2 to indicate that gas boiler units are already existing. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five gas boilers can be available in this region. This means that three units can be optimized. |

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityUpperLimit | Set to 0.98 to indicate an availability factor of 98% for the gas boiler. |
| lifeTime | Set to 20 to indicate a technical life time of 20 years. |

### Accounting

In the accounting module the costs of the fuel boilers are calculated. In our case, we can also account for the carbon emission as indicator here. (Another option would be to include carbon as commodity to be able to use and/or store the emitted CO2.)

**set_indicators.dat**
Examples: CO2, CarbonCost, FuelCost, Invest, OMFix, OMVar

**accounting_sourcesinkFlow.dat**

| accounting_sourcesinkFlow_parameter | Example |
| ------ | ------ |
| perFlow | Set to 0.2 for the indicators “CO2” and “FuelImport” to indicate CO2 emissions of 0.2 kt_CO2/GWh for importing gas. Set to 0.04 for the indicators “FuelCost” and “FuelImport” to indicate gas import costs of 0.04M€/GWh.. |

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set OMVar per activity for the gas boiler to 0.005 to indicate the cost of using the boiler, i.e. 0.005M€/GWh [!!!related to energy input or output???]. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 100 for "Invest" indicator to specify investment costs 100M€/GW [!!!related to energy input or output???]. |
| useAnnuity | Set to one for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 20 for "Invest" indicator to use an amortization time of 20 years. |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 5 for "OMFix" indicator to indicate OMFix costs of 5M€/GW [!!!related to energy input or output???]. |

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | Set CarbonCost for CO2 emission to 0.06 to indicate a cost of 0.06 M€/kt_CO2. Futhermore, you can link the different cost components, such as CO2 costs, fuel costs, OMVar, and OMFix, to the total system costs by setting a value of 1 for each. |
