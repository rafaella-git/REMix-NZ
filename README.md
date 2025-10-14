# REMix-NZ

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)]()
[![Framework: REMix](https://img.shields.io/badge/Framework-REMix-orange.svg)]()

**A regional, hourly multi-energy system model for Aotearoa New Zealand**

**REMix-NZ** is a regional and hourly multi-energy system model for New Zealand, built on the [REMix framework developed by DLR](https://dlr-ve.gitlab.io/esy/remix/framework/).  
It provides scenario-based projections to **2050** across the **power**, **heat**, and **transport** sectors, enabling detailed assessments of electrification and green hydrogen integration.

The model is under continuous development as part of the **HINT – NZ-German Platform for Green Hydrogen Integration** project.

---

## Current Model Scope

| Category | Description |
|-----------|--------------|
| **Projection horizon** | 2020 – 2050 |
| **Sectors** | Power, electrified heat and transport |
| **Temporal resolution** | Hourly (8760 steps per year) |
| **Spatial resolution** | 11 regions based on transmission bottlenecks |
| **Energy carriers** | Electricity and (green) hydrogen |
| **Model basis** | [REMix framework (DLR)](https://dlr-ve.gitlab.io/esy/remix/framework/) |

---

## Model of New Zealand

The country is represented as two main islands – **North Island (NI)** and **South Island (SI)** – subdivided into **11 electrical nodes** corresponding to transmission zones.

### Electricity transfer network
The inter-node network follows Transpower’s (national Grid owner and system operator) transmission system.  
Distances are approximated as **Euclidean × 1.2** to account for real-world routing.

| Line | Start | End |
|------|--------|-----|
| NIS__AKL | NIS | AKL |
| AKL__WTO | AKL | WTO |
| WTO__BOP | WTO | BOP |
| WTO__CEN | WTO | CEN |
| CEN__HBY | CEN | HBY |
| TRN__CEN | TRN | CEN |
| CEN__WEL | CEN | WEL |
| WEL__CAN | WEL | CAN |
| NEL__CAN | NEL | CAN |
| CAN__OTG | CAN | OTG |
| AKL__TRN | AKL | TRN |
| WTO__HBY | WTO | HBY |

## Overview

![image](https://github.com/rafaella-git/energy-nz/assets/135769724/3eab3ebb-4d42-4593-804b-628b7811b7e2)

### Structure layout

```
remix_nz/
├── input/                         # Raw input data organized by type
│   ├── brownfield/
│   ├── demand/                    # Demand file creates project folders
│   ├── profiles/
│   ├── shapefiles/
│   ├── technical/
│   └── xlsx/
│
├── process/                       # Core processing scripts and outputs
│   ├── build_instance.py          # Main script to build model instance
│   ├── run_scenario.py            # Scenario execution script
│   └── evaluate.py                # Evaluation script
│
├── project/                       # Project-specific configurations and inputs
│   └── project-name/case-name
│       ├── data/                  # Model data, created with build_instance.py
│       └── result/                # .gdx files, created with run_scenario.py
│
├── .gitignore                     # Git ignore rules
└── README.md                      # Project documentation
```

### Building and running the model

All model setup and execution is handled via `build_instance.py` and `run_instance.py`.

### Building the model
Key build functions:
- `add_nodes(m)` → defines nodes, spatial scope, and years  
- `add_scope(m)` → sets up scenario horizon and active optimisation years  
- `add_demand(m)` → loads electricity and hydrogen demand  
- `add_renewables(m)` → reads renewable potential and time series  
- `add_network(m)` → establishes transmission and pipeline links  
- `add_accounting(m)` → configures cost and balance equations  

### Running the Model

Run from terminal:
```bash
python src/build_instance.py
```

## Funding Acknowledgement
We thank the Catalyst: Strategic Fund, administered by the Ministry of Business Innovation and Employment of New Zealand, and the German Federal Ministry of Education and Research (grant number 03SF0690) for supporting the project HINT (New Zealand-German Platform for Green Hydrogen Integration). 
