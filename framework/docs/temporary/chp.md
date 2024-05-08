---
title: "Combined Heat and Power (CHP) plants"
lang: en-US
---

# Combined Heat and Power (CHP) plants

## CHP in general

CHP power plants convert a (fossil) fuel into electricity and heat (and CO2). An example is a coal-fueled CHP power
plant. In this power plant coal is burned to generate electricity and heat. At the same time, the power plant emits CO2.

## CHP in REMix

```{figure} /img/REMix_CHP.svg
:align: center
Figure: Schema of CHP implementation in REMix.
```
The core_sourcesink module can provide the fuel needed to run the CHP power plant and the electricity and heat demand
that need to be met.

The essential part to model CHP power plants is the {ref}`core_converter module <remix_model_core_converter_label>`. In
this module commodities can be converted into other commodities. Advanced features, such as the discrete expansion of
power plants and unit commitment can be added to the CHP power plants as well. For further information see
{ref}`MIP features <explanations_mip_label>`

Other than that, in the core_accounting module, the costs of the power plants (fuel, OMfix, OMvar, annuity) and the
emission of CO2 are accounted for.

The {ref}`core balance <remix_model_core_balance_label>` module balances the commodities involved in each time step.

## Data

In this part, the data that is required for modeling CHP power plants is listed and explained for the exemplary
use-case of a coal-fueled CHP plant. Further optional features can be found in
{ref}`remix model reference <remix_model_label>`.

### Source/Sink

In order to be able run the CHP plant, a fuel needs to be available. In our case, coal can be imported from outside the
modeling scope by implementing a fuel source.

**set_commodities.dat**
Example: Electricity, Heat, Coal

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
| upper | Set to inf to indicate that an unlimited annual amount of coal is available. For other fuels such as solid biomass, an annual potential can be specified here. |

                                           lower   upper
    Albania . 2050 . FuelImport . Biomass    0.0   1.143
    Albania . 2050 . FuelImport . Coal       0.0     inf


### Converter

The converter represents the CHP power plant, so in our case a coal-fueled CHP power plant. It converts the imported coal to electricity and heat.

**set_converter_techs.dat**
Example: Coal_CHP

**set_activities.dat**
Example: Powergen (mode in which a converter produces electricity), Heatgen (mode in which converter produces heat), Combinedgen (mode in which converter produces both electricity and heat)

**converter_coefficient.dat**

Here you can indicate the efficiencies of the available activities. The
optimizer decides which activity to use. You can also use the partial load
feature to run the different units of the CHP power plant with different
activities. You can find more information on this feature
{ref}`here <explanations_mip_label>`.

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for the coal in the "Powergen" mode to -2.5 and the coefficient of electricity to one to indicate an efficiency of 40%. For the "Heatgen" mode you can set the coal coefficient to -1.25 and the coefficient of heat to one to indicate an effiency of 80%. You can set the coefficient for the coal in the "Combinedgen" mode to -3, the coefficient of electricity to one and the coefficient of heat to 1.5 to indicate an electricity effiency of 33.3% and an heat efficiency of 50%. |
                                            coefficient
    ST_Coal . 2050 . Powergen . CO2_air           0.609
    ST_Coal . 2050 . Powergen . Coal             -1.818
    ST_Coal . 2050 . Powergen . Elec              1.000


**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2 to indicate that two coal-fueled CHP power plant units are already existing. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five coal-fueled CHP power plants can be available in this region. This means that three units can be optimized. |
                                      optimize    unitsBuild    unitsLowerLimit    unitsUpperLimit
    Czech_Republic . 2050 . ST_Coal          1             0                  0                inf

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityUpperLimit | Set to 0.9 to indicate an availability factor of 90% for the CHP power plant. |
| lifeTime | Set to 40 to indicate a technical life time of 40 years. |

                    activityLowerLimit    activityUpperLimit    freeDecom    lifeTime    mipUnits
    ST_Coal . 2050                 0.0                 0.910          0.0        40.0         0.0

### Accounting

In the accounting module the costs of the CHP plant are calculated. In our case, we can also account for the carbon emission as indicator here. Another option would be to include carbon as commodity to be able to use and/or store the emitted CO2.

**set_indicators.dat**
Examples: CO2, CarbonCost, FuelCost, Invest, OMFix, OMVar

**accounting_sourcesinkFlow.dat**

| accounting_sourcesinkFlow_parameter | Example |
| ------ | ------ |
| perFlow | Set to 0.36 for the indicators "CO2" and "FuelImport" to indicate CO2 emissions of 0.36 kt_CO2/GWh for importing coal. Set to 0.01 for the indicators "FuelCost" and "FuelImport" to indicate coal import costs of 0.01 M€/GWh. |

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set OMVar per activity for the CHP power plant to 0.001 to indicate the cost of using the power plant, i.e. 0.001 M€/GWh. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 1500 for "Invest" indicator to specify investment costs 1500 M€/GW. |
| useAnnuity | Set to 1 for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 40 for "Invest" indicator to use an amortization time of 40 years. |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 60 for "OMFix" indicator to indicate OMFix costs of 60 M€/GW. |

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | Set CarbonCost for CO2 emission to 0.06 to indicate a cost of 0.06 M€/kt_CO2. Futhermore, you can link the different cost components, such as CO2 costs, fuel costs, OMVar, and OMFix, to the total system costs by setting a value of 1 for each. |
