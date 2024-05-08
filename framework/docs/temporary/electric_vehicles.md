---
title: "Electric Vehicles"
lang: en-US
---

(explanation_evs_label)=

# Electric Vehicles

## General

Electric vehicles are a major driver for future electricity demand. Their electric load can be shifted in time to some degree so it is important
to explicitly model it. Usually three modes of charging and discharging are separated: Uncontrolled charging (UC), controlled charging (CC) and
vehicle-to-grid (V2G, discharging back to the grid, thus acting as a small power plant). Types of electric vehicles that are commonly modelled in REMix are battery electric vehicles (BEVs), plug-in hybrid electric vehicles (PHEVs) and fuel cell electric vehicles (FCEVs). Additionally, different vehicle sizes have been separated, only adding to the number of vehicle type and respective profiles. For each vehicle type one set of five profiles has to be supplied.

### Electricity from the grid to battery energy
The connection of the fleet to the grid and thus the fleet's rated capacity varies across the day. This is provided in a profile of absolute fleet connection to the grid in GW defining the maximum electric flow from the grid to the vehicle fleet. It is a upper bound limit to the converter of controlled charging converting the commodity 'electricity' to the commodity 'EV energy stored'.
Part of the charged electricity cannot be shifted because parking times are shorter than required to charge the respective battery to it's full capacity. This part is provided to REMix in the uncontrolled charging profile in TWh per timestep representing a fixed electricity demand from the grid and a fixed inflow to the battery. Electric losses between grid and battery are accounted for so that the demand from the grid is higher than the electric feed-in to the battery. In REMix, uncontrolled charging is represented by a converter with a fixed activity profile.
The activity of the controlled charging converter is determined in the REMix optimization subject to the constraints described above.

### The fleet battery
The fleet battery is a model concept comprising all vehicles of a specific vehicle type. Thus, it's capacity ammounts to the the sum of all vehicle batteries of that vehicle type. Trips that follow a parking and potential plugging situation have to be possible with the battery level. This is ensured by the minimum battery level profile giving a battery level for every hour. Adversely, the battery cannot be completely filled at each time since vehicles are driving and often parking at locations without charging stations. The maximum battery level profile models the case of every vehicle of the fleet charging as soon as it connects with the full rated power until the battery is filled or the next trip is taken. These two profiles form lower and upper bound for the fleet battery level. They are given [in absolute units of TWh].
The fleet battery has a pre-defined feed-in through uncontrolled charging and a pre-defined outflow (both optimization parameters not variables) from the battery to the respective electric motors.

### Battery energy to driving demand
The converter from fleet battery to electric demand for driving translates the commodity 'EV stored energy' into 'electricity for driving'. It has a fixed time series, pre-defined by the electric drain profile for each timestep. This timeseries defines the amount of energy that is drained from the fleet battery in every timestep.


## Modeling in REMix

```{figure} /img/REMix_electric_vehicles.svg
:align: center
Figure: Schema of electric vehicle fleet implementation in REMix.
```

In REMix, electric vehicles are modeled according to the following schema using the three converters, one reservoir and one sink as well as the three commodities 'electricity', 'EV stored energy' and 'electricity for driving'. The five profiles constraining the converter activities and the storage level are shown in the schema and additionally listed in the table below.

| Profile name              | Constrainined entity          | Constraint type   | Unit  |
| ------------------------- | ----------------------------- | ----------------- | ----- |
| Plugging                  | Charging station (converter)  | Upper bound       | GW    |
| Uncontrolled charging     | Charging station (converter)  | Fixed             | TWh   |
| Maximum battery energy    | EV fleet battery (reservoir)  | Upper bound       | TWh   |
| Minimum battery energy    | EV fleet battery (reservoir)  | Lower bound       | TWh   |
| Electric drain            | Fixed (sink)                  | Fixed             | TWh   |

[Profiles constraining the charging of electric vehicles. If you want to read more on how these profiles can be calculated, visit the [documentation page of VencoPy (Vehicle Energy Consumption in Python)](https://vencopy.readthedocs.io/en/latest/index.html) which we use to calculate those profiles.]

All vehicles (in a specific vehicle class such as small BEVs) are combined to one representative vehicle fleet comprising one set of sink, storage and three converters as well as five constraining profiles.

The essential parts to model electric vehicles are the modules
{ref}`core_converter <remix_model_core_converter_label>`,
{ref}`storage <remix_model_core_storage_label>`,
{ref}`core_sourcesink  <remix_model_core_sourcesink_label>`. The converters
tranform the commodity electricity into stored energy and electricity for
driving. Two converters (uncontrolled charging and controlled charging) model
the interface between the grid and the fleet battery thus charging (controlled
and uncontrolled) and discharging (controlled, via vehicle-to-grid) the fleet
battery modeled as storage. The stored energy in turn is then drained from the
fleet battery through a fixed profile given as a sink with a temporally
resolved profile.

Other than that, in the core_accounting module, the variable and fix costs of charging and discharging (OMfix, OMvar) and optional investments in the flexibility are accounted for.

The core_balance module balances the commodities involved in each time step.

## Data

In this part, the data that is required for modeling electric vehicles is listed and explained for the use-case of a fleet of battery electric vehicles (BEV).

### Converter

The converter represents the charging and discharging units of the storage technology, so in our case the multitude of all charging points connecting the electric vehicle fleet to the grid. Thus, its rated charging power varies over time, since vehicles connect and disconnect throughout the day. The converter transforms electricity to chemical energy stored in the lithium ion cell stack. For stored energy in the fleet battery and electricity for driving, commodities have to be defined in *set_commodities*.

**set_commodities.dat**

Example: EV_stored, Driving_energy

While modelling, it is important to specify the converter part specifically as a **converter_techs** in the *set_converter_techs.dat* even if it refers to the same technology entity in a physical sense, just as in the example of a lithium-ion battery above. Additionally, activities have to be defined in *set_activities.dat* and mapped to the respective **set_converter_techs** in *converter_coefficient.dat*

**set_converter_techs.dat**

Example: EVs_CC (controlled charging), EVs_UC (uncontrolled_charging), EVs (drain from fleet battery to fleet demand)

**set_activities.dat**

Examples: Charge, Discharge, Discharge_Driving

**converter_coefficient.dat**

Here the activities available to the converter have to be defined - charging and discharging (controlled charging) and charging (uncontrolled charging).

| converter_coefficient_parameter | Example |
| ------ | ------ |
| coefficient | Set the coefficient for Electricity in "charging" mode to -1 and the coefficient of EV_CC_stored to 0.95 to indicate an efficiency of 95%. For the "discharging" mode you can set the Electricity coefficient to -1 and the coefficient of LiIon_stored to 0.9 to indicate a discharging effiency of 90%. |

Below, an exemplary file is shown with the sets converter, year, activity and commodity. Each converter is mapped in this file to a year, an activity and two commodities, one of which is negative to indicate the conversion. Adjusting the numbers incorporates efficiencies in the conversion.

                                                  coefficient
    EVs_CC . 2050 . Charge . Electricity                -1.00
    EVs_CC . 2050 . Charge . EV_Stored                   0.95
    EVs_CC . 2050 . Discharge . Electricity              0.95
    EVs_CC . 2050 . Discharge . EV_Stored               -1.00
    EVs_UC . 2050 . Charge . Electricity                -1.00
    EVs_UC . 2050 . Charge . EV_Stored                   0.95
    EVs . 2050 . Discharge_Driving . EV_Stored          -1.00
    EVs . 2050 . Discharge_Driving . Driving_Energy      1.00

The converter part of a storage technology can be limited in its activities and
a lifetime can be set after which the converter is decommisioned or replaced.
Those properties are assigned in *converter_techParam*. For more parameters see
the {ref}`core_converter reference <remix_model_core_converter_label>`.

**converter_techParam.dat**

| converter_tech_parameter | Example |
| ------ | ------ |
| activityLowerLimit | Set to 0 to indicate that each activity can be 0 at minimum. |
| activityUpperLimit | Set to 0.9 to indicate an upper bound of the rated charging and discharging capacity. |
| freeDecom | Allow decommissioning of the unit before the end of the technical life time. |
| lifeTime | Set to 40 to indicate a technical life time of 40 years. |

                    activityLowerLimit    activityUpperLimit    freeDecom    lifeTime
    EVs_CC . 2050                    0                  0.98            0          25

As for power plant technologies, electric vehicle storage capacity properties are defined in *converter_capacityParam*.

For energy storage, time independent technical parameters like upper limits, decommissioning rules and lifetime can be explicitly specified in the technical parameter of each storage technology.

**storage_techParam.dat**

| storage_tech_parameter | Example |
| ------ | ------ |
| levelUpperLimit | Set to 1 to set an upper limit for EV fleet battery fill level. |
| freeDecom | Allow decommissioning of the unit before the end of the technical life time. |
| lifeTime | Set to 40 to indicate a technical life time of 40 years. |

Other parameters such as rate of self discharge, and reservoir size per unit are set in *storage_sizeParam*

**storage_sizeParam.dat**

| storage_size_parameter | Example |
| ------ | ------ |
| selfdischarge | Set to 0.01 to set hourly losses of 1% of the storage level. |
| size | Storage unit size. Set to 1 for a 1 GWh / unit storage.  |

To mark the limits on charging and discharging, time series of storage level factors indicating an upper, lower or fixed bounds can be set. With this profiles, lower and upper boundary constraints of the fleet battery can be modeled that are caused by minimum range requirements (lower bound) and limited connectivity or charging power (upper bound). Those profiles can be calculated by electric vehicle charging models such as VencoPy.

**storage_levelprofile.dat**

| storage_levelprofile | Recommended value | Explanation |
| ------ | ------ | ------ |
| t0001 | Value between 0 and 1 | Lower or upper limit to the storage state-of-charge (SOC) in the respective time step. |
| t0002 | Value between 0 and 1 | Lower or upper limit to the storage state-of-charge (SOC) in the respective time step. |
| ... | ... | ... |

                                     t0001    t0002 ...
    DE . 2050 . EV_Battery . lower   0.213    0.215 ...
    DE . 2050 . EV_Battery . lower   0.213    0.215 ...

As for all technologies, the energy capacity for a storage technology and the corresponding converter capacity for every node region modelled can be predefined in the input or determined by a lower and upper bound. Additionally, it is possible to prevent expansion of capacity completely or limit expansion upto a certain extent in a year. For electric vehicles, optimizations of converter capacities (representing charge points and vehicles) are rarely modeled in energy system optimization since vehicle purchase decisions are driven by more than only cost factors.

**converter_capacityParam.dat**

| converter_capacity_parameter | Example |
| ------ | ------ |
| unitsBuild | Set to 2 to indicate that two lithium-ion battery converter units are already existing. |
| unitsLowerLimit | Set to 3 to indicate that at least three lithium-ion battery converter units have to exist. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five lithium-ion battery converter units can be available in this region. This means that five units can be optimized. |

**storage_reservoirParam.dat**

| storage_capacity_parameter | Example |
| ------ | ------ |
| unitsLowerLimit | Set to 3 to indicate that at least three lithium-ion battery reservoir units have to exist. |
| unitsUpperLimit | Set to 5 to indicate that a maximum of five lithium-ion battery reservoir units can be available in this region. This means that five units can be optimized. |
| unitsBuild | Set to 2 to indicate that two lithium-ion battery reservoir units are already existing. |

### Source/Sink

Lastly, driving and thus draining the electric vehicle fleet battery has to be modeled. This can be done by implementing a sink in *set_sourcesink_techs*.

**set_sourcesink_techs.dat**
Example: DemandEV

**sourcesink_config.dat**

| sourcesink_config_parameter | Example |
| ------ | ------ |
| usesLowerProfile | Set to 1 to make sure that a fuel can not be sold outside of the modeling scope, since the default value for the profile will then be zero (if no profile is further specified). |
| usesUpperSum |  Set to 1 to indicate the annual amount of fuel that is availabe. The amount can be specified in sourcesink_annualSum.dat. |

**sourcesink_annualSum.dat**

| sourcesink_annualSum_parameter | Example |
| ------ | ------ |
| usesFixedProfile | Set to 1 to indicate to use a fix value that will get translated into both a lower and upper bound in GAMS. |

                                           usesFixedProfile
    DE . 2050 . DemandEV . Driving_Energy                 1

**sourcesink_profile.dat**

| sourcesink_profile | Recommended value | Explanation |
| ------ | ------ | ------ |
| t0001 | Value between -1 and 0 | Value giving the hourly drain from the fleet battery |
| t0002 | Value between -1 and 0 | Value giving the hourly drain from the fleet battery |
| ... | ... | ... |

                                                    t0001    t0002
    DE . 2050 . DemandEV . Driving_Energy . fixed  -1.279   -0.949

### Accounting

In the accounting module the costs of the electric vehicle fleet charging and discharging are calculated.

**set_indicators.dat**
Examples: Invest, OMFix, OMVar

**accounting_converterActivity.dat**

| accounting_converterActivity_parameter | Example |
| ------ | ------ |
| perActivity | Set OMVar per activity for the converter EVs_CC and activity Charge to 5 to indicate the cost of charging, i.e. 5 M€/GWh (€/kWh). |

                                                     perActivity
    OMVar . global . EVs . 2050 . Discharge_Driving         0.00
    OMVar . global . EVs_CC . 2050 . Charge                 5.00
    OMVar . global . EVs_CC . 2050 . Discharge             23.75
    OMVar . global . EVs_UC . 2050 . Charge                 5.00

**accounting_perIndicator.dat**

| accounting_perIndicator_parameter | Example |
| ------ | ------ |
| perIndicator | You can link the different cost components, such as OMVar defined above to the total system costs by setting a value of 1 for each. |

                                        perIndicator
    SystemCost . OMVar . global . 2050             1
