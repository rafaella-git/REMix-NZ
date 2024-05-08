---
title: Power transmission
lang: en-UK
---

# Power transmission

## Data (ToDo)

In the following, we describe the setup of a model that represents the markets of the so-called CWE-countries (Central Western Europe), which are coupled to each other. Accordingly, we are dealing with a strongly spatially aggregated network for which using the simple transhipment approach is sufficient (at first glance).

### General sets
Modeling power transmission only makes sense for multi-regional models. Therefore, we need to define at least two nodesData and nodesModel

**set_nodesData.dat**
Example: Germany, France, Belgium, Netherlands, Louxembourg

With regard to nodesModel, we may want to aggregate Louxembourg to Germany, whereas the rest of the states is mapped 1:1 from nodesData to nodesModel.

**set_nodesModel.dat**
Example: GerLux, France, Belgium, Netherlands, Louxembourg

To indicate the intended aggregation the mapping has to be explicitely provided.

**map_accNodes.dat**
Example: Germany.GerLux, France.France, Belgium.Belgium, Louxemburg.GerLux, Netherlands.Netherlands

### Transfer

Now we are creating a network, which may consist of several sub-networks or sub-transmission systems. However, at least we define one grid segment

**set_gridSegments.dat**
Example: CWE

The connections for power transmission are finally provided with linksData (the aggregation according to map_accNodes.dat will be automatically performed which is why we do not need to exogenously provide linksModel).

**set_linksData.dat**
Example: GER-FR, GER-NL, GER-BE, FR-BE, FR-LU, BE-NL, BE-LU

Having both nodesData and linksData defined, we can define the network. Therefore, we make sure that each link exactly connects two nodes.

**transfer_linkStartEnd.dat**

| transfer_startEnd_parameter | Example |
| ------ | ------ |
| start | Set to 1, if the direction of a dedicated transmission line is so that it starts at the dedicated node.  |
| end | Set to 1, if the direction of a dedicated transmission line is so that it ends at the dedicated node.  |

For modeling power transmission we also define at least one technology. (Opposed to this example, in tutorial 5 the technology is defined as "HVDC".)

**set_transfer_techs.dat**
Example: cross_border_links

Each of the links may have a length provided with **transfer_lengthParam.dat**. However, for our purpose, where we do not model physical links, this is not of relevant.

Next, we define the technology's features.

**transfer_techParam.dat**

| transfer_techParam_parameter | Example |
| ------ | ------ |
| lifeTime           | Set to 30 to indicate a technical life time of the links for calculating decommissioning. |
| freeDecom          | Set to 0 to disable decommissioning of the links before the end of the technical life time. |
| mipLinks           | Set to 0 to allow any value for capacity expansion the links of this technology |
| flowUpperLimit     | Set to 0.9 to consider a security margin for the maximal power exchange |

After providing information on the technology, we specify the parameters of each individual line.

**transfer_linksParam.dat**
| transfer_linksParam | Example |
| ------ | ------ |
| linksLowerLimit | Set to the NTC published by ENTSO-E to indicate the intially possible power exchange|
| linksBuild      | Do not provide this parameter to indicate that it is equal to linksLowerLimit |
| linksUpperLimit | Set to 15 to allow a maximal capacity expansion for each line up to 15 GW |
| linksDelta      | Do not provide this parameter to indicate any delta |
| noExpansion     | Set to 1, if you like to enable capacity expansion for the dedicated line |
| circuits        | Do not provide this parameter since it is not needed for simple transport models |


## Additional power-flow contraints

For performing a LOPF with DC-power flow contraints, we have to provide the reactance of each line. In
**transfer_reactPerLength.dat** we indicate the specific reactance per length for which the typical unit of
Ohms/km.

Furthermore, it is necessary to specify **transfer_gridSegments.dat**, where we indicate for which gridSegment we want
to perform consider the LOPF by setting **useDCopf** to 1. In **transfer_linksParam.dat** the **circuits** indicate
the number of electrical circuits used for the LOPF.
