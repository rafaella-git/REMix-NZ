# REMix

This repository documents the under-development REMix-NZ model of New Zealand's energy system. The model is based on the REMix framework developed by DLR (please familiarise yourself with the framework [here](https://dlr-ve.gitlab.io/esy/remix/framework/dev/about/introduction.html#about-introduction)) and has been adapted for New Zealand with a **spatial resolution of 11 nodes**.

Currently, the focus is on **electricity**, with future expansion planned for heat and transport.


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



### Transfer network
- **Electricity**: The electriity transfer network follows **New Zealand's electricity transmission system**, with distances based on **Euclidean distance + 20%** to account for real-world deviations.
- **Hydrogen (H2)**: not considered (transport via **trucks** under development, modelled as a **high OPEX, low CAPEX** system).

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

## Building and Running the Model

The model setup is handled via `build_instance.py`. The key functions include:
- `add_nodes(m)`: Defines nodes, spatial resolution, and years for modeling.
- `add_demand_ffe(m)`: Loads commodity demand (electricity, hydrogen, etc.).
- `add_scope(m)`: Sets up the scenario scope (nodes, optimisation years).
- `add_renewables(m)`: Loads renewable energy potential and timeseries data.
- `add_network(m)`: Defines energy transfer infrastructure.
- `add_accounting(m)`: Implements economic and energy balance constraints.

To run the model you can use the command line or run_instance.
To evaluate the results you can use the dashboard.

## Contributions and Contact

- **Contributions:** Pull requests and issues are welcome! Please use Common Commit Types.
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

- **Contact:** Rafaella Canessa (rca139@uclive.ac.nz)


