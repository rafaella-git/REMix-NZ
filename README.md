# REMix
This repository documents the modelling of New Zealand's energy system. Including the electricity, heat and transport sectors.
To set up the model, we use the file `build_instance.py` modified to fit New Zealand. In this case we are using a spatial resolution of 11 nodes. Currently, we are focused on the electricity sector and will later expand the model to include heat and transport. 

# Model of New Zealand 
The two main islands of New Zealand, the north island (NI) and south island (SI) are represented in 11 nodes corresponding to electrical subdvisions (based on study from 
University of Auckland).

| id | Node | Island |
| ------------- | ------------- |------------- |
| 1 |   NIS  | NI |
| 2 | AKL  | NI |
| 3 | WTO  |  NI |
| 4 | BOP  | NI |
| 5 | CEN  |  NI |
| 6 | HBY  | NI |
| 7 | TRN  | NI |
| 8 | WEL |  NI |
| 9 | NEL | SI |
| 10| CAN | SI |
| 11 | OTG |  SI |



### Transport network

The transport network for energy is based on the electricity transmission network. We started modelling only electricity transfers but we are working on other technologies as well. We consider the node location as the main electric substation of the region (except for Canterbury, where the node is taken as the middle point between the two most important substations, Twizel and Islington). The lenght of each conection is the eucledian distance between these points plus 20%, accounting for real lines not being straight.  

| Name of the line | Start | End |
| ------------- | ------------- | ------------- |
| NIS__AKL  | NIS  | AKL  | 
| AKL__WTO  | AKL  | WTO  |
| WTO__BOP  | WTO  | BOP  | 
| WTO__CEN  |  WTO  | CEN  |
| CEN__HBY  |  CEN  | HBY  |
| TRN__CEN  | TRN  | CEN  | 
| CEN__WEL  |  CEN | WEL  | 
| WEL__CAN  | WEL  | CAN  |
| NEL__CAN  | NEL | CAN  |
| CAN__OTG  |  CAN  |  OTG  | 
| AKL__TRN  | AKL | TRN  | 
| WTO__HBY  | WTO  |  HBY  | 

## Transporting commodities via different technologies 

The transport technology can vary (e.g. "HVDC", "Pipeline_CH4", "Pipeline_H2", "Pipeline_H2_retrofit"). Every link gets a maximum transport capacity (GW/hr/line), life time (years).

### Electricity
In electrcity, the lenght of the transmission line is directly proportional to the losses that will be experienced. 

### Hydrogen

When transporting hydrogen (H2) we start with considering shipping via trucks. Modelling this transport is simplified to treating the truck as a link with no CAPEX but high OPEX (due to high costs of fuel, paying the drivers). Regarding the hydrogen demand. We consider adding a storage technology for hydrogen 

![image](https://github.com/rafaella-git/energy-nz/assets/135769724/3eab3ebb-4d42-4593-804b-628b7811b7e2)

# Running the model
Currently the file calls these functions:  add_nodes(m), add_demand_ffe(m), add_scope(m), add_renewables(m), add_lithium_batteries(m), add_network(m), add_accounting(m).

### add_nodes:
Each node is defined in a "nodesData" column, the column "nodesModel" refers to another level of agregation for these nodes. In this case, we are interested in running the model with the 11 regios, but if we wanted to, we could aggregate them  there are 2 "nodesModel" options, north isand (NI) or south island (SI). Here we also select the years to be added to the model (currently 2020, 2025, 2030, 2035, 2040, 2045, 2050), and the years to be optimised (currently 2050).

### add_demand_ffe:
Different commodities can be defined (as shown in the table below).
Q4M (Question for Manuel): why are they named the same? why make a difference in the first place?

|Commodity | Name |
| ------------- | ------------- |
| Electricity | Elec |
|Hydrogen| H2 |
| H2-feedstock | H2 |
| Natural Gas | CH4 |
| Gas | CH4 |
| Feedstock Gas | CH4 |
| Feedstock methanol | CH3OH |
| Renewable Fuels | REfuel |

The demand profiles are drawn from a csv file and follow this structure:

| node | year | sector | carrier | t0001 | ... | t8760 |
| ---- | ---- |------- | ------- |------ | --- |------ | 
| AKL | 2020 | Wholesale|	Electricity	|0.645925 |... |0.5831185 |
| ... |  ... | ... | ... | ... | ... | ... |

The annual slack limite is set as infinite (no limit) but it has a cost of 3.000 EUR/MWh -> 3 MEUR/GWh.

### load_feedin_csv
This function loads the installable capacity for different renewable technologies (from the file "region_statistics_2012.csv" that is mapped to 2050 in the function add_renewables), as well as the generation for each timestep ("timeseries_{year}.csv"), in this case both files where supplied by Manuel and the data belongs to 2012. / _(Q4M: can/should we get more years to do a 5 year mapping from 2020 to 2050? Answer: Not necessarily now, he can provide it though)_ /

From the file region_statistics (as shown in the table below) we get for every region, the installable capacity and the maximum annual energy generation for each renewable generation technology (different types of PV and wind power).

|  region	technology	|  	installable_per_region		|  annual_energy_per_region 	| 
| ---- | ---- |------- | 
| AKL	| pv_central_fixed	|  180668.17	|  	40548476.87	|  
|  AKL |  	pv_central_track_azimuth		|  180668.17	62416678.15	|  
|  AKL		|  pv_decentral		|  2202.53		|  602355.08	|  
|  AKL	|  	wind_offshore_floating		|  34722.03	|  	16797596.68	|  
|  AKL		|  wind_offshore_foundation	|  	12564.32		|  5582488.15	|  
|  AKL	|  	wind_onshore		|  7569.17	|  	2773127.47	|  

From the file "timeseries_{year}.csv" (as shown in the table below) we get for every region and every timestep of the year (8760) the generation from each renewable technology (different types of PV and wind power). 

 | t | 	region | 	technology | 	timeseries_per_region | 
 | ---- | ---- |------- |  ---- |
 | 31/12/2011 23:30 | 	AKL	 | pv_central_fixed | 	77039.56 | 
 | 31/12/2011 23:30	 | AKL	 | pv_central_track_azimuth | 	81403.81 | 
 | 31/12/2011 23:30	 |  AKL | 	pv_decentral | 	1306.1 | 

### function add_renewables
This function adds the renewable generation technologies to the model. It reads the file region_statistics to get the installable/installed capacities (Q4M clarify). It defines as "vintage" as the years taht we can build/install new caacites. In this case we are also mapping the year 2012 (year we have data from) to 2050 (year we want to optimise). When defining technologies, we have PV and wind (csp is not included at the moment).\
Here we set the lifetime for the technologies, for instance, solar PV and wind turbines are expected to have a longer lifespan from 2030 onwards, so that value is adjusted accordingly for those years.
We also define the coefficients for each technology (Electricity is liked to power generation, and heat CSP is linked to heat generation).\
This is where the costs for the different technologies are defined. \

### function add_scope(m)
Q4M is this where we select the year to optimise? why 2045?, also I have NI and SI, but should I keep the 11 nodes?\

### function add_accounting(m)

### function add_scope(m)
#error line 232 with the model\


### Comments I am getting from the code
*"No set elements for grid segments included."
*"No set elements for storage degradation states included."
*"No set elements for SoC included." State of Charge for storage or for EVs?
*"No set elements for accounting nodes included.": what accounting for nodes?
"No set elements for nodesData to nodesAcc mapping included." wdym?
"No bounds for indicator per links included" ?
"No time-dependent coefficients for converter activities included" ?
"No electrical reactance per distance included" 
"No losses per distance for transport lines included" ?
"No grid segements for transport lines included" ?
"No grid segements for transport lines included" ? 
*** Error at line 232: Execution halted: abort$11 'Error encountered during mapping of modelNodes. See the lst file 
for a detailed description.'\

# Disclaimers
When the critical bug was discovered, the new version of REMix was pulled from the DLR's Git. Then the code did not run so some changes were made to the following files: 
* utilities.py (REMix>framework>remix>framework>tools>utilities.py) - Q4M to be verified if this is legal: copied-pasted the functions that were in the old file and that were not here anymore. \
  _Errors saying unable to import gdx function, gt and other gams related errors. Fixed by copy-pasting functions from old REMix files_
* instances.py (REMix>framework>remix>framework>api>instance.py)
  _Error saying "current_data.extend(new_data) setattr(self, f"_{label}", sorted(list(set(current_data)))) '<' not supported between instances of 'str' and 'int'" fixed with replacing lins in the rerror with "current_data.extend(new_data) current_data_strings = [str(item) for item in current_data] setattr(self, f"_{label}", sorted(set(current_data_strings)))"_\


