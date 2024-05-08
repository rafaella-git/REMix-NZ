---
title: "Example project: FlexMex"
lang: en-US
---

(example_projects_label)=

# Example Project: FlexMex

FlexMex stands for *Model experiment for the temporally and spatially resolved investigation of future load balancing in the power system*
Within the MODEX i.e. Modelling experiments for the Energy Transition (*Modellexperimente für die Energiewende* in German), the FlexMex project invoked comparisons of nine energy system models from different research institutions.
The project was funded by the German Federal Ministry of Economic Affairs and Climate Action (BMWK), formerly Federal Ministry of Economic Affairs and Energy (BMWi). All relevant publications of data and reports are listed below.
The scope and requirements of FlexMex translate into a comprehensive model within the REMix framework. The model components (commodities, converters, etc.) and the overall structure are explained in the following.

# Model design

```{figure} /img/REMix_flexmex_example.svg
:align: center
Figure: Illustration of the FlexMex project's model design in REMix.
```

1. For **renewable energy** technologies wind and solar, one way is to consider the available resources as sourcesink technologies and converting them with the Wind turbine or solar PV converters respectively.
In such a case, the converter Solar PV for e.g. will convert the commodity solar resource to electricity. The upper limit on the electricity generation will be anyhow set extrinsically by the resource availability.
However, we simplify the mode in this case, by omitting the source-sink conversion. By setting a time-dependent upper bound to the converter's activity, we limit electricity generation by the resource availability.
Hence, the converters *Solar PV* and *Wind Turbine* have no input commodity only an outgoing commodity, electricity.

2. The **Hydro-reservoir** power plant consists of 4 major components : pump, turbine, reservoir and the run-off. The pump and turbines are connected to the grid and convert electricity to potential energy or vice-versa respectively.
The reservoir is thus a storage technology. To account for river-runoff, we have an additional converter technology that *converts* feed-in to the reservoir.
There is an extrinsically fixed feed-in and a minimum hydro-outflow that constrain the source-sink technologies. The water bypassed by the run-off and the water outlet from the turbine contribute to the overall outflow.

3. The electrolyzer fulfills **hydrogen** demand with some tank storage.

4. Gas turbines, Combined heat-power gas turbines, electric boilers, heat pumps and gas boilers are modelled as converters. The figure above shows the interplay between electricity, heat and gas(only as a source of import).

5. Li-ion battery is modelled simply with a converter and storage representation.

6. Besides the technologies illustrated above, the FlexMex project also incorporated Electric Vehicles, Demand Response and a little diverse representation of heat-based converters (differentiated on basis of capacities).
  They are not elaborated here for simplicity.

To know more about the FlexMex project, please check out the following publications: {cite}`vanOuwerkerk2022, Gils2022a, Gils2022b`.
*If you are curious about the FlexMex model in REMix, feel free to contact us!*
