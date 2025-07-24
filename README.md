# REMix-NZ

REMix-NZ is a regional and hourly multi-energy system model for New Zealand. It includes detailed projections to 2050 across the power, heat, and transport sectors, supporting scenario-based analysis and planning. The model is based on the REMix framework developed by DLR (please familiarise yourself with the framework [here](https://dlr-ve.gitlab.io/esy/remix/framework/dev/about/introduction.html#about-introduction)) and it has been adapted for New Zealand. 

REMix-NZ is under continuous development. Currently it considers:  
- **Projected energy demand**: For New Zealand from 2020 to 2050
- **Sectors**: power and electrified heat and transport.  
- **Resolution**: hourly (8760 steps/year), regional (11 regions based on electricity transmission bottlenecks)
- **Energy carriers**: electricity and green hydrogen.  

## Model of New Zealand 
The model represents the two main islands of Aotearoa New Zealand, the North Island (NI) and South Island (SI), which are represented in 11 nodes corresponding to electrical subdvisions.

### Transfer network
- **Electricity**: the electriity transfer network follows **New Zealand's electricity transmission system**, with distances based on **Euclidean distance + 20%** to account for real-world deviations.

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

![image](https://github.com/rafaella-git/energy-nz/assets/135769724/3eab3ebb-4d42-4593-804b-628b7811b7e2)

### Building and Running the Model

The model setup is handled via `build_instance.py`. The key functions include:
- `add_nodes(m)`: Defines nodes, spatial resolution, and years for modeling.
- `add_demand_ffe(m)`: Loads commodity demand (electricity, hydrogen, etc.).
- `add_scope(m)`: Sets up the scenario scope (nodes, optimisation years).
- `add_renewables(m)`: Loads renewable energy potential and timeseries data.
- `add_network(m)`: Defines energy transfer infrastructure.
- `add_accounting(m)`: Implements economic and energy balance constraints.

To run the model you can use the command line or run_instance.

## Contributions
Pull requests and issues are welcome! Please use Common Commit Types.
- 
    **Common Commit Types**  
  - **feat**: New feature  
  - **fix**: Bug fix  
  - **docs**: Documentation changes  
  - **chore**: Routine tasks (e.g., updating `.gitignore`)  
  - **refactor**: Code changes without affecting behavior  
  - **ci**: Continuous integration updates  
  - **perf**: Performance improvements  
  - **test**: Adding or updating tests  
  - **style**: Code style changes (e.g., formatting, CSS updates)  
  
  **Examples**  
  ✅ **Feature:** `feat: parameterized additional runner config param settings`  
  ✅ **Fix:** `fix: Invalid value for vars parameter: vars map does not contain key …`  
  ✅ **Docs:** `docs: readme update about VPC`  
  ✅ **Chore:** `chore: update gitignore`  
  ✅ **Refactor:** `refactor: change code to parameterize properties instead of hardcoding them`  
  ✅ **CI:** `ci: added Jenkinsfile for Continuous Integration pipeline`  
  ✅ **Performance:** `perf: running multiple threads for concurrent work`  
  ✅ **Test:** `test: added test case for feature ...`  
  ✅ **Style:** `style: changing the fonts - css`  


## Funding Acknowledgement
We thank the Catalyst: Strategic Fund, administered by the Ministry of Business Innovation and Employment of New Zealand, and the German Federal Ministry of Education and Research (grant number 03SF0690) for supporting the project HINT (New Zealand-German Platform for Green Hydrogen Integration). 


## Contact
Rafaella Canessa. [Sustainable Energy Research Group (SERG)](serg.co.nz). University of Canterbury, Christchurch, Aotearoa New Zealand. rca139@uclive.ac.nz


