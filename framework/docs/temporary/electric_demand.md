---
title: "Electricity demand"
lang: en-US
---

# Electricity demand

## General
Energy demand is a crucial component in the design of energy systems. Electricity demand volumes and characteristics drive the necessary power plant park and electricity supply.
In REMix, the major constraint to the optimization is the energy balance which enforces the supply to meet demand at every time step in every modelled region. Thus, electrictiy demand is a main *sink* in this energy balance and it has to be specified for every region and every time step.

## Modeling in REMix
Demand is modelled in REMix as a source-sink technology. This is declared in the **set_sourcesink_techs.dat** input file.
Demand is represented as a time series with **negative** values as it is the sink to the system.
*Electricity* demand is specified by using the combination of source-sink technology, *Demand* and the commodity, *Electricity*. Thus the commodity *Electricity* has to be declared in the **set_commodities.dat** file.
This representation is for flexibility and allows easy inclusions of demand for other energy carriers like heat and fuels.

## Data

With the required set declarations, REMix-specific paramters that can be used to define electricity demand are :

1. **sourcesink_config** in sourcesink_config.dat

In this file, the configuration for electricity demand is set up with the following set dependencies:

- `nodesModel`: every node with a *demand* should be specified
- `years`:  all modelled years are specified respectively
- `sourcesink_techs`: here the source-sink technology should be declared, e.g. **Demand**
- `commodity`:  specify *Electricity*

| sourcesink_config_parameter | Recommended value | Explanation |
| ------ | ------ | ------ |
| usesFixedProfile | 1 | Set to 1 to make sure that electricity demand will be a fixed value every time step and not inherited as an upper or lower limit. |
| usesLowerProfile | 0 | If set to 1, the commidity is demanded (sink) or supplied (source) above a lower limit profile in the given year for the given technology and commodity. |
| usesLowerSum | 0 | If set to 1, the commidity is demanded (sink) or supplied (source) above an annual lower limit in the given year for the given technology and commodity. |
| usesUpperProfile | 0 | If set to 1, the commidity is demanded (sink) or supplied (source) below an upper limit profile in the given year for the given technology and commodity.|
| usesUpperSum | 0 | If set to 1, the commidity is demanded (sink) or supplied (source) below an annual upper limit in the given year for the given technology and commodity.|

An example .dat file for sourcesink_config.dat can is shown below.

                                    usesFixedProfile    usesLowerProfile
    Algeria . 2050 . Demand . CH4                  1                   0
    Algeria . 2050 . Demand . Elec                 1                   0
    Algeria . 2050 . Demand . H2                   1                   0

2. **sourcesink_timeSeries** in  sourcesink_timeSeries.dat

As the name suggests, the demand for each time step is specified in this file in a tabular format with each timestep (e.g. t0001) given as a column. For each time step, the sets to be specified are:
- `nodesModel`: Every node with a *demand* should be specified
- `years`:  Every modeled year can be specified respectively
- `sourcesink_techs`: Specify **Demand**
- `commodity`:  Specify *Electricity*
- `profile_type`: Can be 'lower', 'upper' or 'fixed'
- `timeModel`: Every time step with a value

| sourcesink_timeSeries | Recommended value | Explanation |
| ------ | ------ | ------ |
| t0001 | Value between -1 and 0 | Value of the normalized profile in each time step. |
| t0002 | Value between -1 and 0 | Value of the normalized profile in each time step. |
| ... | ... | ... |

                                             t0001    t0002    ...
    Albania . 2050 . Demand . Elec . fixed  -1.955   -1.896    ...

The type of profile i.e. *fixed* as indicated above is specified here again by the sourcesink_timeSeries parameter. The given value is normalized, i.e. has to be between -1 and 1.

```{note}
The demand values should be negative to act like a sink in the energy balance!
```

For more information about the other options available for the parameters, please check out the **technical documentation** on sourcesink_techs.

## Data sources

Typical data sources for electricity demands are:

1. For profiles on different geographic and temporal resolutions:
   [Open Power System Data](https://data.open-power-system-data.org/time_series/)
2. Detailed data for Control Areas, Bidding Zones and Countries can also be retrieved from the
   [entso-e transparency platform](https://transparency.entsoe.eu/load-domain/r2/totalLoadR2/show)
3. Data from the TYNDP can be retrieved from the latest scenario report supplementaries, e.g.
   [for the TYNDP 2022](https://2022.entsos-tyndp-scenarios.eu/download/)
