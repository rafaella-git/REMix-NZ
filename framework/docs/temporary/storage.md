---
title: "Electric and Heat Storage"
lang: en-US
---

# Electric and Heat Storage

Energy system storage technologies are crucial sources of flexibilities in the energy system and hence essential in energy system optimization. Commonly modelled technologies are lithium-ion batteries, pumped hydro storage and hydrogen cavern or tank storage.

In the following, basic aspects and features of modelling energy storage in
REMix-MISO are explained for a linear programming problem. Additional features
like {ref}`Mixed Integer Programming (MIP) <explanations_mip_label>` that allow
advanced modelling of energy storage technologies, are not discussed here.

```{figure} /img/REMix_storage.svg
:align: center
Figure: Modeling concept for storage at the example of a Lithium Ion battery.
```

## Modelling in REMix

According to the overall structure of REMix, a storage technology can be divided into a converter part and a storage part as illustrated in the following schema:

```{figure} /img/old_files/LIB_graphic.png
:align: center
Figure: Schema of storage implementation in REMix.
```

The storage part can be understood as a reservoir storing energy. It is just the storage part that is represented by a *storage_tech* modelled in REMix. The storage technology has to be declared in the *set_storage_techs.dat* file.


To characterise storage completely, another important element to specify is the **commodity**. As commodities refer to energy-carriers, all the commodities stored and converted by the storage and converter technologies respectively must be included in the list of commodities in *set_commodities.dat*

**set_commodities.dat**
Example: Electricity, LiIon_stored

Due to these multiple requirements for characterizing storage, many electricity and thermal storage technologies can be modelled in the same way in REMix.

## Data

In this part, the data that is required for modeling storage is listed and explained for the use-case of a lithium ion battery. Further optional features can be found in the
{ref}`storage <remix_model_core_storage_label>` and the
{ref}`converter <remix_model_core_converter_label>` module documentation.

### Converter

The converter represents the charging and discharging units of the storage technology, so in our case a power electronic equipment. It converts electricity to chemical energy stored in the lithium ion cell stack.

The converter part performs the **charging** and **discharging** activities and converts energy from one form e.g. electricity to the form stored and vice versa. While modelling, it is important to specify the converter part specifically as a **converter_techs** in the *set_converter_techs.dat* even if it refers to the same technology entity in a physical sense, just as in the example of a lithium-ion battery above.

**set_converter_techs.dat**
Example: LiIonBattery

This decoupling adds a flexibility that is particularly important if multiple technologies can provide the same commodity for storage (think of a CSP plant where internal storage could be provided either by gas or by the solarthermal field) or if the charging and discharging technologies are different.
For example: In hydrogen storage the electrolyzer is associated only with charging as it converts electricity to hydrogen. The discharging is performed by the fuel cell. The activities of *charging* and *discharging* are also declared as set elements in the *set_activities.dat* file.

**set_activities.dat**
Example: Charge, Discharge

**converter_coefficient.dat**

Here you can indicate the efficiencies of activities available to the converter - charging and discharging.

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for Electricity in "charging" mode to -1 and the coefficient of LiIon_stored to 0.95 to indicate an efficiency of 95%. For the "discharging" mode you can set the Electricity coefficient to -1 and the coefficient of LiIon_stored to 0.9 to indicate a discharging effiency of 90%. |

                                                      coefficient    ...
    LiIonBattery . 2050 . Charge . Electricity              -1.00    ...
    LiIonBattery . 2050 . Charge . LiIon_stored              0.95    ...
    LiIonBattery . 2050 . Discharge . Electricity            0.90    ...
    LiIonBattery . 2050 . Discharge . LiIon_stored          -1.00    ...

The converter part of a storage technology can be limited in its activities and a lifetime can be set after which the converter is decommisioned or replaced. Those properties are assigned in *converter_techParam*.

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityLowerLimit | Set to 0 to indicate that each activity can be 0 at minimum. |
| activityUpperLimit | Set to 0.9 to indicate an upper bound of the rated charging and discharging capacity. |
| freeDecom | Allow decommissioning of the unit before the end of the technical life time. |
| lifeTime | Set to 40 to indicate a technical life time of 40 years. |

                         activityLowerLimit    activityUpperLimit    freeDecom    lifeTime
    LiIonBattery . 2050                   0                  0.98            0          25

As for power plant technologies, storage converter capacity properties are defined in *converter_capacityParam*.


For energy storage, time independent technical parameters like upper limits, decommissioning rules and lifetime can be explicitly specified in the technical parameter of each storage technology.

**storage_techParam.dat**

| storage_tech_parameter | Example |
| ------ | ------ |
| levelUpperLimit | Set to 1 to set an upper limit for storage fill level. |
| freeDecom | Allow decommissioning of the unit before the end of the technical life time. |
| lifeTime | Set to 40 to indicate a technical life time of 40 years. |

Other parameters such as rate of self discharge, and reservoir size per unit are set in *storage_sizeParam*

**storage_sizeParam.dat**

| storage_size_parameter | Example |
| ------ | ------ |
| selfdischarge | Set to 0.01 to set hourly losses of 1% of the storage level. |
| size | Storage unit size. Set to 1 for a 1 GWh / unit storage.  |

To mark the limits on charging and discharging, time series of storage level factors indicating an upper, lower or fixed bounds can be set. This is useful in problems where it is important to maintain a certain level of storage for a short or long duration.

**storage_levelprofile.dat**

| storage_levelprofile | Recommended value | Explanation |
| ------ | ------ | ------ |
| t0001 | Value between 0 and 1 | Lower or upper limit to the storage state-of-charge (SOC) in the respective time step. |
| t0002 | Value between 0 and 1 | Lower or upper limit to the storage state-of-charge (SOC) in the respective time step. |
| ... | ... | ... |

                                     t0001    t0002    ...
    DE . 2050 . LiIonBattery . upper     1        1    ...

As for all technologies, the energy capacity for a storage technology and the corresponding converter capacity for every node region modelled can be predefined in the input or determined by a lower and upper bound. Additionally, it is possible to prevent expansion of capacity completely or limit expansion upto a certain extent in a year.

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsBuild | Set to 2 to indicate that two lithium-ion battery converter units are already existing. |
| unitsLowerLimit | Set to 3 to indicate that at least three lithium-ion battery converter units have to exist. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five lithium-ion battery converter units can be available in this region. This means that five units can be optimized. |

**storage_reservoirParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 3 to indicate that at least three lithium-ion battery reservoir units have to exist. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five lithium-ion battery reservoir units can be available in this region. This means that five units can be optimized. |
| unitsBuild | Set to 2 to indicate that two lithium-ion battery reservoir units are already existing. |

```{note}
Electric Vehicles are also modelled like storage technologies. More information can be found {ref}`here <explanation_evs_label>`
```

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 2 to indicate that two coal-fueled CHP power plant units are already existing. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five coal-fueled CHP power plants can be available in this region. This means that three units can be optimized. |
                                       optimize    unitsBuild    unitsLowerLimit    unitsUpperLimit
    Czech_Republic . 2050 . ST_Coal           1             0                  0                inf

### Accounting

In the accounting module the costs of storage technologies are calculated.

**set_indicators.dat**
Examples: Invest, OMFix, OMVar

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set OMVar per activity for the lithium-ion battery unit to 0.001 to indicate the cost of using the battery unit, i.e. 0.001 M€/GWh. This has to be specified for each activity i.e. Charging and Discharging. |

**accounting_converterUnits.dat**

| accounting_converterUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 1500 for "Invest" indicator to specify investment costs 1500 M€/GW. |
| useAnnuity | Set to 1 for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 40 for "Invest" indicator to use an amortization time of 40 years. |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 60 for "OMFix" indicator to indicate OMFix costs of 60 M€/GW. |

                                            perUnitBuild  useAnnuity  amorTime  interest  perUnitTotal
    Invest . global . LiIonBattery . 2050          50000           1        25      0.07             0
    OMFix . global . LiIonBattery . 2050               0           0         0      0.07           110

**accounting_storageUnits.dat**

| accounting_storageUnits_parameter | Example |
| ------ | ------ |
| perUnitBuild | Set to 1500 for "Invest" indicator to specify investment costs 1500 M€/GWh. |
| useAnnuity | Set to 1 for "Invest" indicator to indicate that the cost should be calculated as annuities instead of absolute cost. |
| amorTime | Set to 40 for "Invest" indicator to use an amortization time of 40 years. |
| interest | Set to 0.05 for "Invest" indicator to indicate an interest rate of 5% for the calculation of the annuity factor. |
| perUnitTotal | Set to 60 for "OMFix" indicator to indicate OMFix costs of 60 M€/GWh. |
