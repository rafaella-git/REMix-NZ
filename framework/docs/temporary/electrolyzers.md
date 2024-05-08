---
title: "Electrolyzers"
lang: en-US
---

# Electrolysis technologies

## General

Electrolysis describes the process of separating water into it's chemical components hydrogen and oxygen. This can be done by different technological systems such as alkaline electrolysis (AEL) using an electrolytic fluid (usually Potassium-hydroxide - KOH), polymer exchange membrane (PEM) with solid ion transfer media or solid oxide electrolysis cells (SOEC) operating at high temperatures. This last technology can be operated exothermally with high electricity demand and endothermally with lower electricity demand. Hydrogen can then be stored in its chemical form, injected into the gas grid or further synthesized.

## Modeling in REMix

```{figure} /img/REMix_electrolysers.svg
:align: center
Figure: Schema of electrolyzer implementation in REMix.
```

Electrolyzers in REMix are implemented as converters (see
{ref}`core_converter module <remix_model_core_converter_label>`,) between the
commodities electricity and hydrogen. An efficiency has to be given as well as
techno-economic properties such as specific investments, fix and variable
operating costs.

Other than that, in the core_accounting module, the costs of the electrolyzers
(OMfix, OMvar, annuity) can be calculated.

The core_balance module balances the commodities involved in each time step.

Advanced features, such as the discrete expansion of electrolyzer unit size
and unit commitment can be added to the electrolyzer as well. For further
informations see {ref}`MIP features <explanations_mip_label>`.

## Data

In this part, the data that is required for modeling electrolyzers is listed and explained.

In order for the electrolyzer to run, the input and output commodities need to be available.

**set_commodities.dat**
Example: Electricity, Hydrogen, Heat or respective aliases such as Elec and H2

### Converter

The converter represents the electrolyzer, e.g. an AEL unit with all its ancillary equipment. It is modeled as a black-box converting electricity (and potentially heat) to hydrogen (and potentially heat).

**set_converter_techs.dat**
Example: Electrolyzer_AEL


**set_activities.dat**
Example: HydroGen (mode in which a converter produces hydrogen)

**converter_coefficient.dat**

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for electricity in the "HydroGen" mode to -2.5 and the coeffient of hydrogen to 1 to indicate an efficiency of 40%. |

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2 to indicate that two electrolyzer units are already existing. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five electrolyzer units can be available in this region. This means that three units can be optimized. |

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityUpperLimit | Set to 0.9 to indicate an availability factor of 90% for the electrolyzer. |
| lifeTime | Set to 20 to indicate a technical life time of 20 years. |

### Accounting

In the accounting module, the costs associated with expansion and operation of the electrolyzer are calculated. For this, the cost components have to be defined in the respective indicator set.

**set_indicators.dat**
Examples: Invest, OMFix, OMVar

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set variable operation and maintenance costs (OMVar) per activity for the electrolyzer plant to 0.001 to indicate the cost of using the power plant, i.e. 0.001 M€/GWh. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 1500 for "Invest" indicator to specify investment costs 1500 M€/GW rated electrical power. |
| useAnnuity | Set to one for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 40 for "Invest" indicator to use an amortization time of 40 years in the annuity calculation. This is not relevant for the technical lifetime specified in converter_tech_parameter |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 60 for "OMFix" indicator to indicate OMFix costs of 60 M€/GW rated electric power. |

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | You can link the different cost components, such as OMVar, and OMFix, to the total system costs by setting a value of 1 for each. |
