# REMix-NZ input builder
# This script writes CSV input tables for the GAMS REMix-NZ model.
# CSVs will be exported to a data folder, REMix can be run from that folder

# REMix-NZ input builder
# This script writes CSV input tables for the GAMS REMix-NZ model.
# CSVs will be exported to a data folder, REMix can be run from that folder

import time
import numpy as np
import pandas as pd
from pathlib import Path
from remix.framework.api.instance import Instance

idx = pd.IndexSlice


# =============================================================================
# USER CHOICES
# =============================================================================

# nz_case_GP_2020-2035-2050)
group_name = "GP-NT-ELEC-BIO-H2"       # Demand folder group (output folder structure)
base_scenario = "GP"                   # Scenario name (must match Scenario column in CSV/Excel)

# Years:
# - yrs_sel: years the MODEL optimises over (subset of yrs_to_calc)
# - yrs_to_calc: years you have DATA for / want to write tables for
yrs_sel = [2020, 2050]  
# yrs_sel = [2020, 2045, 2050]                 # model run: [2020, 2025, 2030, 2035, 2040, 2045, 2050]
yrs_to_calc = list(yrs_sel)            # functions may extend this (e.g., brownfield build years)

# Case folder name (output goes to ../project/<group_name>/<case_name>/data)
years_tag = "-".join(str(y) for y in yrs_sel)
case_name = f"nz_case_{base_scenario}_{years_tag}"

# Demand inputs (GP-NT-ELEC-BIO-H2 scenarios)
demand_folder = Path("C:/Local/REMix/remix_nz/input/demand/GP-NT-ELEC-BIO-H2")
hourly_electricity_file = "hourly_electricity_GP-NT-ELEC-BIO-H2.csv"
carriers_excel_file = "data_summary.xlsx"   # must contain sheet "carriers" (used later)

# Hydro inflows (simplified hydro)
use_hydro_cascade = False  # MVP: False
simplified_hydro_inflow_file = Path("C:/Local/REMix/remix_nz/input/brownfield/hydro/inflows_remix-nz/regional-hydro-inflow.csv")

eur_per_aud = 0.62  # update as needed
# -------------------------------------------------------------------
# CSIRO fuel price medians (AUD/GJ) used for MVP
# -------------------------------------------------------------------
# Source: CSIRO GenCost fuel price spreadsheets.
# https://www.aemo.com.au/energy-systems/major-publications/integrated-system-plan-isp/2026-integrated-system-plan-isp/2025-26-inputs-assumptions-and-scenarios
# For each fuel there are multiple plant entries. We use the median value.
# Ranges (min..max) are left in comments for transparency.
csiro_fuel_aud_per_gj = {
    "Black Coal_early": 5.19,
    "Black Coal_2050": 4.09,
    "Brown Coal_early": 1.04,
    "Brown Coal_2050": 1.04,
    "Gas_early": 16.62,
    "Gas_2050": 16.65,
    "Liquid Fuel_early": 31.89,
    "Liquid Fuel_2050": 32.04,
    "Biomass_early": 0.66,
    "Biomass_2050": 0.66,
}

# Which supply blocks to include
include_renewables = True
include_geothermal = True
include_hydro_simple = True
include_thermal = True
include_gas_turbines = True
include_lithium_batteries = True
include_h2_techs = True
include_heat_transp = True


# =============================================================================
# PATHS
# =============================================================================

path_input = Path("C:/Local/REMix/remix_nz/input")
path_demand = demand_folder
path_profiles = path_input / "profiles"
path_brownfield = path_input / "brownfield"

base_dir = Path(f"../project/{group_name}/{case_name}")
data_dir = base_dir / "data"
results_dir = base_dir / "result"
data_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

if use_hydro_cascade:
    inflow_file = str(
        path_brownfield / "hydro" / "inflows_remix-nz" / "inflow_2021-to-2020_2013-to-2030_ENERGY_GW.csv"
    )
else:
    inflow_file = str(simplified_hydro_inflow_file)

print("\n--- SCENARIO SETUP COMPLETE ---")
print("case:", case_name)
print("scenario:", base_scenario)
print("years (optimise):", yrs_sel)
print("years (data):", yrs_to_calc)
print("data_dir:", data_dir.resolve())
print("demand folder:", path_demand.resolve())
print("hydro inflow:", inflow_file)

# =============================================================================
# UTILITIES
# =============================================================================

def map_region_to_remix(region_value: str) -> str:
    """
    Converts region names like  "Auckland" into REMix region codes "(AKL)".
    """
    region_to_remix = {
        "Auckland": "AKL",
        "Bay of Plenty": "BOP",
        "Canterbury": "CAN",
        "Gisborne": "HBY",
        "Hawkes Bay": "HBY",
        "Hawke's Bay": "HBY",
        "Manawatu-Whanganui": "CEN",
        "Manawatū-Whanganui": "CEN",
        "Marlborough": "NEL",
        " Marlborough": "NEL",
        "Nelson": "NEL",
        "Northland": "NIS",
        "Otago": "OTG",
        "Southland": "OTG",
        "Taranaki": "TRN",
        "Tasman": "NEL",
        "Waikato": "WTO",
        "Wellington": "WEL",
        "West Coast": "CAN",
    }
    s = str(region_value).strip()
    if "(" in s and ")" in s:
        return s.split("(")[-1].split(")")[0].strip()
    return region_to_remix.get(s, None)

def aud_per_gj_to_meur_per_gwh(aud_per_gj: float) -> float:
    """
    Convert AUD/GJ to million EUR/GWh using eur_per_aud.

      1 GJ = 0.277777... MWh  =>  1 MWh = 3.6 GJ
      AUD/MWh = AUD/GJ * 3.6
      EUR/MWh = AUD/MWh * eur_per_aud
      M€ / GWh = (EUR/MWh) * 0.001
    """
    aud_per_mwh = aud_per_gj * 3.6
    eur_per_mwh = aud_per_mwh * eur_per_aud
    meur_per_gwh = eur_per_mwh * 0.001
    return float(meur_per_gwh)

def csiro_fuel_cost_meur_per_gwh(year: int, fuel_key: str) -> float:
    """
    Return CSIRO-based fuel cost in M€ / GWh for a given year.

    Rule:
      - year < 2050  ->  use '*_early'
      - year >= 2050 ->  use '*_2050'

    fuel_key must be: "Black Coal", "Brown Coal", "Gas", "Liquid Fuel", "Biomass".
    """
    y = int(year)
    suffix = "early" if y < 2050 else "2050"
    key = f"{fuel_key}_{suffix}"
    if key not in csiro_fuel_aud_per_gj:
        raise KeyError(f"CSRIO fuel key not in table: {key}")
    aud = csiro_fuel_aud_per_gj[key]
    cost_meur_per_gwh = aud_per_gj_to_meur_per_gwh(aud)

    print(
        f"[csiro_fuel_cost] year={y}, fuel_key='{fuel_key}', suffix='{suffix}', "
        f"AUD/GJ={aud}, M€/GWh={cost_meur_per_gwh:.6f}"
    )
    return cost_meur_per_gwh

# =============================================================================
# FUNCTIONS 
# =============================================================================

def add_scope(m):
    """
    Scope: nodes mapping + set nodesdata/nodesmodel + years + yearssel.
    """
    df = pd.DataFrame(
        [
            ["AKL", "AKL", 1],
            ["BOP", "BOP", 1],
            ["NEL", "NEL", 1],
            ["NIS", "NIS", 1],
            ["OTG", "OTG", 1],
            ["TRN", "TRN", 1],
            ["WEL", "WEL", 1],
            ["WTO", "WTO", 1],
            ["CAN", "CAN", 1],
            ["CEN", "CEN", 1],
            ["HBY", "HBY", 1],
        ]
    )
    df.columns = ["nodesData", "nodesModel", "aggregate"]
    df.set_index(["nodesData", "nodesModel"], inplace=True)
    df["aggregate"] = ""
    df.columns = [""]

    m["Base"].map.add(df, "aggregatenodesmodel")

    # sets
    m["Base"].set.add(
        list(sorted(set(m["Base"].map.aggregatenodesmodel.index.get_level_values(0)))),
        "nodesdata",
    )
    m["Base"].set.add(
        list(sorted(set(m["Base"].map.aggregatenodesmodel.index.get_level_values(1)))),
        "nodesmodel",
    )

    # years
    # years = all years that data is provided for in the model
    m["Base"].set.add(sorted({int(y) for y in yrs_to_calc}), "years")
    # yearssel = years to be optimised
    m["Base"].set.add(sorted({int(y) for y in yrs_sel}), "yearssel")

def add_network(m):
    
    # First we need to set up the link connections in the data by defining the starting and ending node of each link
    link_names = [
        "NIS__AKL",
        "AKL__WTO",
        "WTO__BOP",
        "WTO__CEN",
        "CEN__HBY",
        "TRN__CEN",
        "CEN__WEL",
        "WEL__CAN",
        "NEL__CAN",
        "CAN__OTG",
        "AKL__TRN",
        "WTO__HBY",
    ]

    m["Base"].set.add(link_names, "linksdata")
    m["Base"].set.add(link_names, "linksmodel")

    link_types = ["land", "sea"]
    m["Base"].set.add(link_types, "link_types")

    # Nodes already defined by add_scope()
    nodes_data = sorted(set(m["Base"].set.nodesdata))

    # ------------------------------------------------------------------
    # transfer_linkstartend  (linksdata x nodesdata -> start/end flags)
    # ------------------------------------------------------------------
    link_connections = pd.DataFrame(
        index=pd.MultiIndex.from_product([link_names, nodes_data], names=["linksdata", "nodesdata"])
    )

    link_connections.loc[idx["NIS__AKL", "NIS"], "start"] = 1
    link_connections.loc[idx["NIS__AKL", "AKL"], "end"] = 1

    link_connections.loc[idx["AKL__WTO", "AKL"], "start"] = 1
    link_connections.loc[idx["AKL__WTO", "WTO"], "end"] = 1

    link_connections.loc[idx["WTO__BOP", "WTO"], "start"] = 1
    link_connections.loc[idx["WTO__BOP", "BOP"], "end"] = 1

    link_connections.loc[idx["WTO__CEN", "WTO"], "start"] = 1
    link_connections.loc[idx["WTO__CEN", "CEN"], "end"] = 1

    link_connections.loc[idx["CEN__HBY", "CEN"], "start"] = 1
    link_connections.loc[idx["CEN__HBY", "HBY"], "end"] = 1

    link_connections.loc[idx["TRN__CEN", "TRN"], "start"] = 1
    link_connections.loc[idx["TRN__CEN", "CEN"], "end"] = 1

    link_connections.loc[idx["CEN__WEL", "CEN"], "start"] = 1
    link_connections.loc[idx["CEN__WEL", "WEL"], "end"] = 1

    link_connections.loc[idx["WEL__CAN", "WEL"], "start"] = 1
    link_connections.loc[idx["WEL__CAN", "CAN"], "end"] = 1

    link_connections.loc[idx["NEL__CAN", "NEL"], "start"] = 1
    link_connections.loc[idx["NEL__CAN", "CAN"], "end"] = 1

    link_connections.loc[idx["CAN__OTG", "CAN"], "start"] = 1
    link_connections.loc[idx["CAN__OTG", "OTG"], "end"] = 1

    link_connections.loc[idx["AKL__TRN", "AKL"], "start"] = 1
    link_connections.loc[idx["AKL__TRN", "TRN"], "end"] = 1

    link_connections.loc[idx["WTO__HBY", "WTO"], "start"] = 1
    link_connections.loc[idx["WTO__HBY", "HBY"], "end"] = 1

    link_connections = link_connections.fillna(0)
    m["Base"].parameter.add(link_connections, "transfer_linkstartend")

    # ------------------------------------------------------------------
    # transfer_lengthparam  (linksdata x link_types -> km)
    # ------------------------------------------------------------------
    link_lengths = pd.DataFrame(
        index=pd.MultiIndex.from_product([link_names, link_types], names=["linksdata", "link_types"])
    )
    link_lengths.loc[idx["NIS__AKL", "land"], "length"] = 149.8
    link_lengths.loc[idx["AKL__WTO", "land"], "length"] = 136.0
    link_lengths.loc[idx["WTO__BOP", "land"], "length"] = 76.2
    link_lengths.loc[idx["WTO__CEN", "land"], "length"] = 154.5
    link_lengths.loc[idx["CEN__HBY", "land"], "length"] = 96.5
    link_lengths.loc[idx["TRN__CEN", "land"], "length"] = 115.1
    link_lengths.loc[idx["CEN__WEL", "land"], "length"] = 83.1
    link_lengths.loc[idx["WEL__CAN", "land"], "length"] = 316.4
    link_lengths.loc[idx["NEL__CAN", "land"], "length"] = 203.9
    link_lengths.loc[idx["CAN__OTG", "land"], "length"] = 179.5
    link_lengths.loc[idx["AKL__TRN", "land"], "length"] = 200.1
    link_lengths.loc[idx["WTO__HBY", "land"], "length"] = 114.6
    link_lengths = link_lengths.fillna(0)
    m["Base"].parameter.add(link_lengths, "transfer_lengthparam")

    # ------------------------------------------------------------------
    # Transfer tech / caps / coefficients / losses
    # ------------------------------------------------------------------
    transport_techs = ["HV"]  
    m["Base"].set.add(transport_techs, "transfer_techs")

    years_model = list(m["Base"].set.yearssel)
    commodity = ["Elec"]

    # transfer_linksparam (linksdata x years x transfer_techs)
    link_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product([link_names, years_model, transport_techs],
                                         names=["linksdata", "years", "transfer_techs"])
    )
    link_caps.loc[idx[:, :, "HV"], "linksUpperLimit"] = 100
    link_caps = link_caps.fillna(0)
    m["Base"].parameter.add(link_caps, "transfer_linksparam")

    # transfer_techparam (transfer_techs x years)
    tech_params = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, years_model], names=["transfer_techs", "years"])
    )
    tech_params["lifeTime"] = 40
    tech_params["flowUpperLimit"] = 1
    m["Base"].parameter.add(tech_params, "transfer_techparam")

    # transfer_coefficient (transfer_techs x years x commodities)
    transfer_coefficient = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, years_model, commodity],
                                         names=["transfer_techs", "years", "commodities"])
    )
    transfer_coefficient["coefficient"] = 1
    m["Base"].parameter.add(transfer_coefficient, "transfer_coefficient")

    # transfer_coefperflow
    coef_per_flow = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, years_model, commodity],
                                         names=["transfer_techs", "years", "commodities"])
    )
    coef_per_flow["coefPerFlow"] = -0.014
    m["Base"].parameter.add(coef_per_flow, "transfer_coefperflow")

    # transfer_coefperlength (transfer_techs x years x commodities x link_types)
    coef_per_dist = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, years_model, commodity, link_types],
                                         names=["transfer_techs", "years", "commodities", "link_types"])
    )
    coef_per_dist.loc[idx[:, :, :, "land"], "coefPerLength"] = -0.00004
    coef_per_dist.loc[idx[:, :, :, "sea"], "coefPerLength"] = -0.00003
    coef_per_dist = coef_per_dist.fillna(0)
    m["Base"].parameter.add(coef_per_dist, "transfer_coefperlength")

    # ------------------------------------------------------------------
    # Accounting 
    # ------------------------------------------------------------------
    cost_indicators = ["Invest", "OMFix"]
    area = ["global"]

    transfer_indicators = pd.DataFrame(
        index=pd.MultiIndex.from_product([cost_indicators, area, ["horizon"], transport_techs, years_model],
                                         names=["indicator", "regionscope", "timescope", "transfer_techs", "years"])
    )
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "perLinkBuild"] = 180
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "interest"] = 0.06
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "amorTime"] = 40
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "useAnnuity"] = 1
    transfer_indicators.loc[idx["OMFix", "global", "horizon"], "perLinkTotal"] = 1.8
    transfer_indicators = transfer_indicators.fillna(0)
    m["Base"].parameter.add(transfer_indicators, "accounting_transferlinks")

    indicators_distance = pd.DataFrame(
        index=pd.MultiIndex.from_product([cost_indicators, area, ["horizon"], transport_techs, years_model, link_types],
                                         names=["indicator", "regionscope", "timescope", "transfer_techs", "years", "link_types"])
    )

    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "perLengthBuild"] = 0.544
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "useAnnuity"] = 1
    indicators_distance.loc[idx["OMFix", "global", "horizon", :, :, "land"], "perLengthTotal"] = 0.00544

    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "perLengthBuild"] = 0.975
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "useAnnuity"] = 1
    indicators_distance.loc[idx["OMFix", "global", "horizon", :, :, "sea"], "perLengthTotal"] = 0.00975

    indicators_distance = indicators_distance.fillna(0)
    m["Base"].parameter.add(indicators_distance, "accounting_transferperlength")

    print("Network added.")

def add_accounting(m):
    # The value global uses all the regions in the system
    # whereas the value horizon takes into account all years in the set set.yearssel
    accounting_indicatorBounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], ["horizon"], ["SystemCost"]])
    )
    accounting_indicatorBounds["obj"] = -1  # minimization of system costs
    accounting_indicatorBounds["discount"] = 0.02  # social discount rate for the indicators
    m["Base"].parameter.add(accounting_indicatorBounds, "accounting_indicatorbounds")

    accounting_perIndicator = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["SystemCost"],
                ["Invest", "OMFix", "FuelCost", "FuelCost", "SlackCost"],
                ["global"],
                ["horizon"],
            ]
        )
    )
    accounting_perIndicator["perIndicator"] = 1
    m["Base"].parameter.add(accounting_perIndicator, "accounting_perindicator")

    # Price electricity slack so "SlackCost" has something to count.
    # Slack sourcesink itself is created in add_demand() (sourcesink_config + sourcesink_annualsum).
    years_model = list(m["Base"].set.yearssel)

    slack_cost = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["SlackCost"], ["global"], years_model, ["Slack"], ["Elec"]],
            names=["indicator", "regionscope", "years", "techs", "commodities"],
        )
    )
    slack_cost["perFlow"] = 1e6
    m["Base"].parameter.add(slack_cost, "accounting_sourcesinkflow")

def validate_scope(m):
    """
    Scope validation (years + nodes + a few sanity checks).

    - Confirms yearssel ⊆ years and matches yrs_sel
    - Confirms years contains all yrs_to_calc (after any function extensions)
    - Confirms nodesdata/nodesmodel exist
    """
    years = sorted(int(y) for y in m["Base"].set.years)
    ysel = sorted(int(y) for y in m["Base"].set.yearssel)

    # yearssel must be subset of years
    missing = [y for y in ysel if y not in years]
    if missing:
        raise ValueError(f"yearssel contains years not in years: {missing}")


    expected_sel = sorted(int(y) for y in yrs_sel)
    if ysel != expected_sel:
        raise ValueError(f"Model yearssel {ysel} does not match yrs_sel {expected_sel}")

    # yrs_to_calc should be contained in years set (after extensions)
    expected_years = sorted(set(int(y) for y in yrs_to_calc))
    if sorted(set(years)) != expected_years:
        # Don’t hard-fail on ordering/dup issues; fail only if missing
        missing2 = [y for y in expected_years if y not in years]
        if missing2:
            raise ValueError(f"Model years missing expected yrs_to_calc entries: {missing2}")

    nodesdata = sorted(set(m["Base"].set.nodesdata))
    nodesmodel = sorted(set(m["Base"].set.nodesmodel))
    if not nodesdata:
        raise ValueError("nodesdata is empty")
    if not nodesmodel:
        raise ValueError("nodesmodel is empty")

    print("\n--- VALIDATE SCOPE ---")
    print("nodesdata:", len(nodesdata), nodesdata)
    print("nodesmodel:", len(nodesmodel), nodesmodel)
    print("years (data):", years)
    print("yearssel (optimise):", ysel)

# conventional

# def add_thermal(m):
#     """
#     Brownfield coal / diesel / biogas plants with fuel, costs and CO2.

#     Techno-economic data & emission factors:
#       - Gulagi et al. 2025 supplementary Table S1 for CAPEX, fixed OPEX,
#         fossil fuel emission factors (coal 354.6, diesel 262.4 kgCO2/MWh). 
#       - CSIRO GenCost fuel prices (AUD/GJ) converted to M€/GWh via eur_per_aud. 
#     """
#     print("\n--- ADDING THERMAL GENERATORS ---")
#     global yrs_to_calc

#     inst_csv = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")

#     # Map primary fuel in database to model tech + model fuel commodity
#     tech_map = {
#         "Coal":   ("Thermal_Coal",   "Coal"),     # Coal steam units
#         "Diesel": ("Thermal_Diesel", "Diesel"),   # Oil / diesel engines
#         "Biogas": ("Thermal_Bio",    "Biofuel"),  # Biogas / biomass (treated as CO2‑neutral) 

#         "Biomass":  ("Thermal_Bio", "Biofuel"),
#         "Wood":     ("Thermal_Bio", "Biofuel"),
#         "Wood waste": ("Thermal_Bio", "Biofuel"),
#     }

#     techs = sorted({v[0] for v in tech_map.values()})
#     fuel_commodities = {v[0]: v[1] for v in tech_map.values()}
#     nodes = list(m["Base"].set.nodesdata)
#     activities = ["Powergen"]

#     years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
#     base_year = int(years_sel[0])

#     # ---------- FILTER BROWNFIELD UNITS ----------
#     df = inst_csv[
#         (inst_csv["Type"] == "Thermal") &
#         (inst_csv["Primary_fuel"].isin(tech_map.keys()))
#     ].copy()

#     if df.empty:
#         print("No thermal rows found in brownfield database.")
#         return

#     df["techs"] = df["Primary_fuel"].map(lambda f: tech_map[f][0])

#     df["Year_built"] = pd.to_numeric(df["Year_built"], errors="coerce")
#     df["Capacity_MW"] = pd.to_numeric(df["Capacity_MW"], errors="coerce").fillna(0.0)
#     df = df.dropna(subset=["Year_built"])
#     df["Year_built"] = df["Year_built"].astype(int)

#     # Only brownfield units up to base_year
#     df = df.loc[df["Year_built"] <= base_year].copy()

#     # Clean node names and keep only model nodes
#     df["Node"] = df["Node"].astype(str).str.strip()
#     df = df.loc[df["Node"].isin(nodes)].copy()
#     df = df.loc[df["Capacity_MW"] > 0].copy()

#     if df.empty:
#         print(f"No thermal brownfield plants with Year_built <= {base_year} in model nodes.")
#         return

#     # ---------- YEARS / CAPACITY ----------
#     years_build = sorted(df["Year_built"].unique().tolist())
#     yrs_to_calc = sorted({int(y) for y in (list(yrs_to_calc) + years_build)})
#     m["Base"].set.add(yrs_to_calc, "years")
#     years_data = list(yrs_to_calc)

#     # 1) Tech parameters: simple lifetime / activity limit 
#     techparam = pd.DataFrame(
#         index=pd.MultiIndex.from_product([techs, years_data], names=["techs", "years"])
#     )
#     techparam["lifeTime"] = 30
#     techparam["activityUpperLimit"] = 1
#     m["Base"].parameter.add(techparam, "converter_techparam")

#     # 2) Capacity: brownfield builds (GW) + generous expansion caps
#     cap = (
#         df.groupby(["Node", "Year_built", "techs"])["Capacity_MW"]
#         .sum()
#         .rename("unitsBuild")
#         .to_frame()
#         .div(1e3)  # MW → GW
#     )
#     cap.index = cap.index.set_names(["nodesdata", "years", "techs"])

#     cap_upper = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [nodes, years_data, techs],
#             names=["nodesdata", "years", "techs"],
#         )
#     )
#     cap_upper["unitsUpperLimit"] = 100.0
#     cap_upper.loc[idx[:, [base_year], :], "noExpansion"] = 1  # no expansion in base year only

#     cap_full = pd.concat([cap, cap_upper], axis=1).sort_index()
#     m["Base"].parameter.add(cap_full, "converter_capacityparam")

#     # ---------- 3) COEFFICIENTS: FUEL IN, ELECTRICITY OUT ----------
#     # Thermal efficiencies (electric output / fuel input), consistent across costs+emissions.
#     #FIXME: correct Representative values aligned with literature for steam plants. 
#     dea_eff = {
#         "Thermal_Coal":   0.35,
#         "Thermal_Diesel": 0.38,
#         "Thermal_Bio":    0.30,
#     }

#     coef = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [techs, years_data, activities,
#              ["Elec"] + sorted(set(fuel_commodities.values()))],
#             names=["techs", "years", "activities", "commodities"],
#         )
#     )
#     coef["coefficient"] = 0.0

#     # +1 GWh_ele per unit of activity
#     coef.loc[idx[:, :, "Powergen", "Elec"], "coefficient"] = 1.0

#     # Fuel input: -1/eff on fuel commodity, so 1 GWh_el requires 1/eff GWh_fuel
#     for t in techs:
#         fuel = fuel_commodities[t]
#         eff = float(dea_eff[t])
#         coef.loc[idx[t, :, "Powergen", fuel], "coefficient"] = -1.0 / eff

#     m["Base"].parameter.add(coef, "converter_coefficient")

#     # ---------- 4) ACCOUNTING: INVEST + OMFix ----------
#     # CAPEX values in M€/GW, from Gulagi et al. 2025 Table S1 (large power plants and
#     # biomass boilers; adapted from €/kW). 
#     acc_units = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [["Invest", "OMFix"], ["global"], ["horizon"], techs, years_sel],
#             names=["indicator", "regionscope", "timescope", "techs", "years"],
#         )
#     ).sort_index()

#     acc_units.loc[idx["Invest", :, :, "Thermal_Coal", :],   "perUnitBuild"] = 1600  # M€/GW
#     acc_units.loc[idx["Invest", :, :, "Thermal_Diesel", :], "perUnitBuild"] = 900   # M€/GW
#     acc_units.loc[idx["Invest", :, :, "Thermal_Bio", :],    "perUnitBuild"] = 2600  # M€/GW

#     acc_units.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
#     acc_units.loc[idx["Invest", :, :, :, :], "amorTime"]   = 30
#     acc_units.loc[idx["Invest", :, :, :, :], "interest"]   = 0.06

#     # Fixed O&M as % of CAPEX per year, approx values from Gulagi et al. 2025. 
#     invest_vals = acc_units.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
#     invest_vals.index = pd.MultiIndex.from_tuples(
#         [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
#         names=acc_units.index.names,
#     )
#     acc_units.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.033
#     m["Base"].parameter.add(acc_units, "accounting_converterunits")

#     # ---------- 5) ACCOUNTING: FuelCost + CO2_emission ----------
#     acc_act = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [["FuelCost", "CO2_emission"], ["global"], ["horizon"], techs, years_sel, activities],
#             names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
#         )
#     ).sort_index()

#     # Direct combustion emission factors from Gulagi et al. 2025 Table S1 (kgCO2/MWh_fuel). 
#     ef_kg_per_mwh_fuel = {
#         "Thermal_Coal":   354.6,   # coal
#         "Thermal_Diesel": 262.4,   # fuel oil / diesel
#         "Thermal_Bio":    0.0,     # biogenic CO2 treated as neutral 
#     }

#     # Map to CSIRO fuel cost keys (AUD/GJ) for GenCost medians.
#     tech_to_csiro_key = {
#         "Thermal_Coal":   "Black Coal",
#         "Thermal_Diesel": "Liquid Fuel",
#         "Thermal_Bio":    "Biomass",
#     }

#     for y in years_sel:
#         suffix = "early" if int(y) < 2050 else "2050"
#         for t in techs:
#             eff = float(dea_eff[t])

#             # CO2: kg/MWh_fuel → kt/GWh_el
#             ef_fuel = ef_kg_per_mwh_fuel[t]
#             kt_per_gwh_el = (ef_fuel * 1e-6) / eff   # (kg/MWh_fuel * 1e-6 kt/kg) / eff
#             acc_act.loc[idx["CO2_emission", "global", "horizon", t, y, "Powergen"],
#                         "perActivity"] = kt_per_gwh_el

#             # Fuel cost: CSIRO AUD/GJ → EUR/MWh_fuel → M€/GWh_fuel → /eff = M€/GWh_el
#             cs_key = tech_to_csiro_key[t]
#             aud_per_gj = float(csiro_fuel_aud_per_gj[f"{cs_key}_{suffix}"])  # CSIRO GenCost.
#             eur_per_mwh_fuel = (aud_per_gj * 3.6) * float(eur_per_aud)
#             meur_per_gwh_fuel = eur_per_mwh_fuel * 0.001
#             acc_act.loc[idx["FuelCost", "global", "horizon", t, y, "Powergen"],
#                         "perActivity"] = meur_per_gwh_fuel / eff

#     m["Base"].parameter.add(acc_act, "accounting_converteractivity")

#     print("yearssel:", list(m["Base"].set.yearssel))
#     print("converter_capacityparam years (unique):",
#         sorted(set(m["Base"].parameter.converter_capacityparam.index.get_level_values(1))))


#     print("Thermal added with fuel input, Gulagi 2025 costs, CSIRO fuel prices, and S1 CO2 factors.")

# def add_gas_turbines(m):
#     """
#     Add existing gas turbines + generic expansion options.

#     Techs: GT, CCGT, OCGT
#     Activity: Powergen
#     Fuel:     CH4 (Natural gas), modelled via csiro_fuel_cost_meur_per_gwh
#     Output:   Elec
#     """

#     global yrs_to_calc, path_brownfield
#     gt_inst_csv = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")

#     # -----------------------------
#     # 1) Sets / basic lists
#     # -----------------------------
#     print("\n--- ADDING GAS TURBINES ---")
#     gt_vintage    = [1950, 2030]                # parameter vintages
#     gt_techs      = ["GT", "CCGT", "OCGT"]
#     gt_nodes      = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]
#     gt_activities = ["Powergen"]

#     # -----------------------------
#     # 2) Tech parameters (per vintage)
#     # -----------------------------
#     gt_tech = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [gt_techs, gt_vintage],
#             names=["techs", "years"],
#         )
#     )
#     gt_tech["lifeTime"]          = 30
#     gt_tech["activityUpperLimit"] = 1
#     m["Base"].parameter.add(gt_tech, "converter_techparam")

#     # -----------------------------
#     # 3) Brownfield capacity + expansion caps
#     # -----------------------------
#     # Existing units from plant database (Natural gas only)
#     df = gt_inst_csv
#     filtered_df = df[df["Primary_fuel"] == "Natural gas"].copy()
#     grouped_df = (
#         filtered_df.groupby(["Node", "Year_built", "Techs"])["Capacity_MW"]
#         .sum()
#         .reset_index()
#     )
#     grouped_df["Year_built"] = grouped_df["Year_built"].astype(int)

#     gt_cap_existing = (
#         grouped_df
#         .set_index(["Node", "Year_built", "Techs"])
#         .rename(columns={"Capacity_MW": "unitsBuild"})
#         .div(1e3)  # MW → GW
#     )

#     # Upper bounds for each (node, year, tech)
#     gt_cap_upper = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [gt_nodes, yrs_to_calc, gt_techs],
#             names=["nodesdata", "years", "techs"],
#         )
#     )
#     gt_cap_upper["unitsUpperLimit"] = 100.0       # generous upper bound
#     gt_cap_upper.loc[idx[:, 2020, :], "noExpansion"] = 1  # only existing in 2020

#     gt_cap_full = pd.concat([gt_cap_existing, gt_cap_upper], axis=1)
#     m["Base"].parameter.add(gt_cap_full, "converter_capacityparam")

#     # -----------------------------
#     # 4) Conversion coefficients
#     # -----------------------------
#     gt_coef = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [gt_techs, gt_vintage, gt_activities, ["Elec", "CH4"]],
#             names=["techs", "years", "activities", "commodities"],
#         )
#     )

#     # Output: 1 GWh_e per 1 unit of activity
#     gt_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1.0

#     # Fuel input: -1 / efficiency (C. Habib)
#     # Efficiencies by vintage: [eff_1950, eff_2030]
#     eff_GT   = [0.41, 0.43]    # simple GT
#     eff_CCGT = [0.58, 0.61]    # combined cycle
#     eff_OCGT = [0.4695, 0.4695]

#     gt_coef.loc[idx["GT",   :, "Powergen", "CH4"], "coefficient"] = np.round(-1.0 / np.array(eff_GT), 3)
#     gt_coef.loc[idx["CCGT", :, "Powergen", "CH4"], "coefficient"] = np.round(-1.0 / np.array(eff_CCGT), 3)
#     gt_coef.loc[idx["OCGT", :, "Powergen", "CH4"], "coefficient"] = np.round(-1.0 / np.array(eff_OCGT), 3)

#     m["Base"].parameter.add(gt_coef, "converter_coefficient")

#     # -----------------------------
#     # 5) Accounting: converter units (CAPEX + OMFix)
#     # -----------------------------
#     gt_acc = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [["Invest", "OMFix"], ["global"], ["horizon"], gt_techs, gt_vintage],
#             names=["indicator", "regionscope", "timescope", "techs", "years"],
#         )
#     ).sort_index()

#     # CAPEX per unit [M€ / unit], per vintage (Gulagi et al. 2025)
#     gt_acc.loc[idx["Invest", "global", "horizon", "GT",   :], "perUnitBuild"] = [900, 830]
#     gt_acc.loc[idx["Invest", "global", "horizon", "CCGT", :], "perUnitBuild"] = [775, 775]
#     gt_acc.loc[idx["Invest", "global", "horizon", "OCGT", :], "perUnitBuild"] = [475, 475]

#     gt_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"] = 1
#     gt_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"]   = 30
#     gt_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"]   = 0.06

#     # OMFix = 1.93 % of CAPEX per year (LUT / Gulagi et al. 2025)
#     invest_vals = gt_acc.loc[
#         idx["Invest", "global", "horizon", :, :], "perUnitBuild"
#     ]
#     invest_vals.index = pd.MultiIndex.from_tuples(
#         [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
#         names=gt_acc.index.names,
#     )
#     gt_acc.loc[
#         idx["OMFix", "global", "horizon", :, :], "perUnitTotal"
#     ] = invest_vals * 0.0193

#     gt_acc = gt_acc.fillna(0.0)
#     m["Base"].parameter.add(gt_acc, "accounting_converterunits")

#     # -----------------------------
#     # 6) Accounting: converter activity (CO2 + fuel cost CH4)
#     # -----------------------------
#     gt_emission = pd.DataFrame(
#         index=pd.MultiIndex.from_product(
#             [["CO2_emission", "FuelCost"], ["global"], ["horizon"], gt_techs, gt_vintage, gt_activities],
#             names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
#         )
#     ).sort_index()

#     # CO2 intensity for natural gas from Gulagi et al. 2025 (kgCO2/MWh_LHV): 204.8
#     # Convert to kt CO2 / GWh_e using efficiency: (kg/MWh_fuel * 1e-6 kt/kg) / efficiency
#     kgCO2_per_MWh = 204.8

#     for tech, eff in zip(
#         ["GT", "CCGT", "OCGT"],
#         [eff_GT, eff_CCGT, eff_OCGT],
#     ):
#         eff_arr = np.array(eff)
#         co2_kt_per_GWh = (kgCO2_per_MWh * 1e-6) / eff_arr
#         gt_emission.loc[
#             idx["CO2_emission", "global", "horizon", tech, :, "Powergen"],
#             "perActivity",
#         ] = co2_kt_per_GWh

#     # FuelCost: use csiro_fuel_cost_meur_per_gwh("Gas") for each vintage year
#     # mapping vintage year to cost
#     for y in gt_vintage:
#         cost_meur_per_gwh = csiro_fuel_cost_meur_per_gwh(int(y), "Gas")
#         for tech in gt_techs:
#             gt_emission.loc[
#                 idx["FuelCost", "global", "horizon", tech, y, "Powergen"],
#                 "perActivity",
#             ] = cost_meur_per_gwh

#     gt_emission = gt_emission.fillna(0.0)
#     m["Base"].parameter.add(gt_emission, "accounting_converteractivity")


def add_thermal(m, path_brownfield, add_fuel_imports=True, debug=False):
    """
    Brownfield coal / diesel / bio steam units with DEA vintages.

    Vintages (converter_techparam years): 1950, 2020, 2030, 2050
      - 1950 uses DEA-2015 parameters (old fleet)
      - 2020 uses DEA-2020
      - 2030 uses DEA-2030
      - 2050 uses DEA-2050

    Brownfield:
      - Plants mapped to vintages by build year:
          Year_built < 2020 -> 1950
          2020 <= Year_built < 2030 -> 2020
          2030 <= Year_built < 2050 -> 2030
          Year_built >= 2050        -> 2050
      - converter_capacityparam has unitsBuild at (nodesdata, vintage_year, techs)

    Expansion:
      - 2020 is dispatch-only: noExpansion=1 (no new thermal capacity in 2020)

    Fuel:
      - Fuel commodities: Coal, Diesel, Biofuel
      - Optional FuelImport_* sourcesinks for fuel supply with FuelCost perFlow.
    """
    print("\n--- ADDING THERMAL (DEA vintages + brownfield) ---")

    idx = pd.IndexSlice
    nodes = list(m["Base"].set.nodesdata)
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)   # optimisation years (e.g. [2020, 2050])
    base_year = years_sel[0]
    activities = ["Powergen"]

    # Vintage years for converters
    vintages = [1950, 2020, 2030, 2050]
    dea_year_for_vintage_map = {1950: 2015, 2020: 2020, 2030: 2030, 2050: 2050}

    def map_build_to_vintage(y):
        if y < 2020:
            return 1950
        elif y < 2030:
            return 2020
        elif y < 2050:
            return 2030
        else:
            return 2050

    # 1) Read brownfield DB and map to techs/fuels
    inst = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")

    tech_map = {
        "Coal":       ("Thermal_Coal",   "Coal"),
        "Diesel":     ("Thermal_Diesel", "Diesel"),
        "Biogas":     ("Thermal_Bio",    "Biofuel"),
        "Biomass":    ("Thermal_Bio",    "Biofuel"),
        "Wood":       ("Thermal_Bio",    "Biofuel"),
        "Wood waste": ("Thermal_Bio",    "Biofuel"),
    }

    df = inst[
        (inst["Type"] == "Thermal")
        & (inst["Primary_fuel"].isin(tech_map.keys()))
    ].copy()
    if df.empty:
        print("Thermal: no matching rows in brownfield DB.")
        return

    df["techs"] = df["Primary_fuel"].map(lambda f: tech_map[f][0])
    df["fuel"]  = df["Primary_fuel"].map(lambda f: tech_map[f][1])

    df["Year_built"] = pd.to_numeric(df["Year_built"], errors="coerce")
    df["Capacity_MW"] = pd.to_numeric(df["Capacity_MW"], errors="coerce")
    df = df.dropna(subset=["Year_built", "Capacity_MW"]).copy()
    df["Year_built"] = df["Year_built"].astype(int)
    df["Node"] = df["Node"].astype(str).str.strip()
    df = df[df["Node"].isin(nodes)].copy()
    df = df[df["Capacity_MW"] > 0].copy()

    if df.empty:
        print("Thermal: no plants with positive capacity in model nodes.")
        return

    # Map build years to vintages
    df["vintage_year"] = df["Year_built"].apply(map_build_to_vintage)

    techs = sorted(df["techs"].unique().tolist())
    fuels = sorted(df["fuel"].unique().tolist())
    fuel_by_tech = {v[0]: v[1] for v in tech_map.values()}

    if debug:
        print("[TH] years_sel:", years_sel, "base_year:", base_year)
        print("[TH] techs:", techs, "fuels:", fuels)
        print("[TH] brownfield capacity by tech, vintage:")
        print(df.groupby(["techs", "vintage_year"])["Capacity_MW"].sum())

    # 2) DEA tables

    # efficiency (η = GWh_el / GWh_fuel)
    dea_eff_table = {
        "Thermal_Coal":   {2015: 0.46,  2020: 0.485, 2030: 0.52,  2050: 0.535},
        "Thermal_Diesel": {2015: 0.37,  2020: 0.37,  2030: 0.37,  2050: 0.37},
        "Thermal_Bio":    {2015: 0.42,  2020: 1.25,  2030: 0.45,  2050: 1.37},
    }
    # lifetime (years)
    dea_lifetime = {
        "Thermal_Coal":   {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "Thermal_Diesel": {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "Thermal_Bio":    {2015: 25, 2020: 50, 2030: 25, 2050: 50},
    }
    # CAPEX (M€/GW)
    dea_invest_meur_per_gw = {
        "Thermal_Coal":   {2015: 2182.4,      2020: 2148.4,      2030: 2103.2,      2050: 2012.8},
        "Thermal_Diesel": {2015: 372.180451,  2020: 361.546724,  2030: 361.546724,  2050: 361.546724},
        "Thermal_Bio":    {2015: 1063.4,      2020: 3136.9,      2030: 957.0,       2050: 3030.6},
    }
    # OMFix fraction
    dea_omfix_frac = {
        "Thermal_Coal":   {2015: 0.016,       2020: 0.016,       2030: 0.016,       2050: 0.016},
        "Thermal_Diesel": {2015: 0.025142857, 2020: 0.025882353, 2030: 0.024847059, 2050: 0.023811765},
        "Thermal_Bio":    {2015: 0.010,       2020: 0.012,       2030: 0.010,       2050: 0.010},
    }
    # CO2 factors (kg CO2 / MWh_fuel)
    ef_kg_per_mwh_fuel = {
        "Thermal_Coal":   354.6,
        "Thermal_Diesel": 262.4,
        "Thermal_Bio":    0.0,
    }
    # Fuel cost (M€/GWh_fuel)
    fuel_cost_meur_per_gwh_fuel = {
        "Coal":    0.03,
        "Diesel":  0.12,
        "Biofuel": 0.06,
    }

    def dea_year_for_tech_vintage(t, v):
        return dea_year_for_vintage_map[v]

    def eff_for_tech_vintage(t, v):
        dy = dea_year_for_tech_vintage(t, v)
        return float(dea_eff_table[t][dy])

    def lifetime_for_tech_vintage(t, v):
        dy = dea_year_for_tech_vintage(t, v)
        return int(dea_lifetime[t][dy])

    def capex_for_tech_year(t, y):
        candidates = [2015, 2020, 2030, 2050]
        dy = min(candidates, key=lambda d: abs(d - y))
        return float(dea_invest_meur_per_gw[t][dy])

    def omfix_for_tech_year(t, y):
        candidates = [2015, 2020, 2030, 2050]
        dy = min(candidates, key=lambda d: abs(d - y))
        capex = float(dea_invest_meur_per_gw[t][dy])
        frac = float(dea_omfix_frac[t][dy])
        return capex * frac

    # 3) converter_techparam (vintages)
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, vintages], names=["techs", "years"])
    )
    techparam["activityUpperLimit"] = 1.0
    for t in techs:
        for v in vintages:
            techparam.loc[(t, v), "lifeTime"] = lifetime_for_tech_vintage(t, v)
    m["Base"].parameter.add(techparam, "converter_techparam")

    # 4) converter_capacityparam: brownfield builds by vintage + expansion caps
    cap_build = (
        df.groupby(["Node", "vintage_year", "techs"])["Capacity_MW"]
          .sum()
          .rename("unitsBuild")
          .to_frame()
          .div(1e3)
    )
    cap_build.index = cap_build.index.set_names(["nodesdata", "years", "techs"])

    cap_bounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_sel, techs],
                                         names=["nodesdata", "years", "techs"])
    )
    cap_bounds["unitsUpperLimit"] = 100.0
    cap_bounds["noExpansion"] = 0

    # 2020 dispatch only: noExpansion=1 (no new builds)
    if base_year in years_sel:
        cap_bounds.loc[idx[:, [base_year], :], "noExpansion"] = 1

    cap_full = pd.concat([cap_build, cap_bounds], axis=1).sort_index()
    m["Base"].parameter.add(cap_full, "converter_capacityparam")

    # 5) converter_coefficient (Elec output, fuel input)
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, vintages, activities, ["Elec"] + fuels],
                                         names=["techs", "years", "activities", "commodities"])
    )
    coef["coefficient"] = 0.0
    coef.loc[idx[:, :, "Powergen", "Elec"], "coefficient"] = 1.0
    for t in techs:
        fuel = fuel_by_tech[t]
        for v in vintages:
            eff = eff_for_tech_vintage(t, v)
            coef.loc[idx[t, v, "Powergen", fuel], "coefficient"] = -1.0 / eff
    m["Base"].parameter.add(coef, "converter_coefficient")

    # 6) accounting_converterunits (Invest, OMFix)
    acc_units = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, years_sel],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()
    for t in techs:
        for y in years_sel:
            acc_units.loc[("Invest", "global", "horizon", t, y), "perUnitBuild"] = capex_for_tech_year(t, y)
            acc_units.loc[("Invest", "global", "horizon", t, y), "useAnnuity"] = 1
            acc_units.loc[("Invest", "global", "horizon", t, y), "amorTime"] = lifetime_for_tech_vintage(t, 2020)
            acc_units.loc[("Invest", "global", "horizon", t, y), "interest"] = 0.05
            acc_units.loc[("OMFix", "global", "horizon", t, y), "perUnitTotal"] = omfix_for_tech_year(t, y)
    m["Base"].parameter.add(acc_units, "accounting_converterunits")

    # 7) accounting_converteractivity (FuelCost, CO2_emission) per vintage
    acc_act = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["FuelCost", "CO2_emission"], ["global"], ["horizon"], techs, vintages, activities],
            names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
        )
    ).sort_index()
    acc_act["perActivity"] = 0.0
    for t in techs:
        fuel = fuel_by_tech[t]
        for v in vintages:
            eff = eff_for_tech_vintage(t, v)
            # Fuel cost: M€/GWh_el
            fc_fuel = float(fuel_cost_meur_per_gwh_fuel[fuel])
            acc_act.loc[("FuelCost", "global", "horizon", t, v, "Powergen"), "perActivity"] = fc_fuel / eff
            # CO2: kt/GWh_el = (kg/MWh * 1e-3)/eff
            kg_per_mwh = float(ef_kg_per_mwh_fuel[t])
            kt_per_gwh_el = (kg_per_mwh * 1e-3) / eff
            acc_act.loc[("CO2_emission", "global", "horizon", t, v, "Powergen"), "perActivity"] = kt_per_gwh_el
    m["Base"].parameter.add(acc_act, "accounting_converteractivity")

    # 8) Fuel imports: FuelImport_Coal, FuelImport_Diesel, FuelImport_Biofuel
    if add_fuel_imports:
        fuels_needed = sorted(set(fuels))  # e.g. ["Biofuel", "Coal", "Diesel"]
        # annualsum: (nodesdata, years, sourcesinks, commodities)
        rows = [(n, y, f"FuelImport_{f}", f) for n in nodes for y in years_sel for f in fuels_needed]
        ss_idx = pd.MultiIndex.from_tuples(
            rows, names=["nodesdata", "years", "sourcesinks", "commodities"]
        )
        ss = pd.DataFrame(index=ss_idx)
        ss["lower"] = 0.0
        ss["upper"] = 1e12
        m["Base"].parameter.add(ss[["lower", "upper"]], "sourcesink_annualsum")

        # config: usesLowerSum + usesUpperSum
        cfg = pd.DataFrame(index=ss_idx)
        cfg["usesLowerSum"] = 1
        cfg["usesUpperSum"] = 1
        m["Base"].parameter.add(cfg, "sourcesink_config")

        # accounting_sourcesinkflow: FuelCost perFlow for imports
        acc_rows = [
            ("FuelCost", "global", y, f"FuelImport_{f}", f)
            for y in years_sel for f in fuels_needed
        ]
        acc_idx = pd.MultiIndex.from_tuples(
            acc_rows,
            names=["indicator", "regionscope", "years", "sourcesinks", "commodities"],
        )
        acc = pd.DataFrame(index=acc_idx)
        acc["perFlow"] = [
            float(fuel_cost_meur_per_gwh_fuel[f]) for (_, _, _, _, f) in acc_rows
        ]
        m["Base"].parameter.add(acc, "accounting_sourcesinkflow")

        # ensure import techs in set
        try:
            m["Base"].set.add([f"FuelImport_{f}" for f in fuels_needed], "sourcesink_techs")
        except Exception:
            pass

    print("Thermal: added techs:", techs, "vintages:", vintages)

def add_gas_turbines(m, path_brownfield, add_fuel_imports=True, debug=False):
    """
    Gas turbines (GT, CCGT, OCGT) with DEA vintages + brownfield.

    Vintages: 1950, 2020, 2030, 2050
      - 1950 uses DEA-2015
      - 2020 uses DEA-2020
      - 2030 uses DEA-2030
      - 2050 uses DEA-2050

    Brownfield:
      - Plants mapped to vintages by build year (same mapping as thermal)
      - converter_capacityparam has unitsBuild at (nodesdata, vintage_year, techs)

    2020 is dispatch-only for GT: noExpansion=1 (no new GT builds in 2020).

    Fuel:
      - CH4 commodity
      - Optional FuelImport_CH4 sourcesinks for supply.
    """
    print("\n--- ADDING GAS TURBINES (DEA vintages + brownfield) ---")

    idx = pd.IndexSlice
    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    base_year = years_sel[0]
    activities = ["Powergen"]

    gt_techs = ["GT", "CCGT", "OCGT"]
    vintages = [1950, 2020, 2030, 2050]
    dea_year_for_vintage_map = {1950: 2015, 2020: 2020, 2030: 2030, 2050: 2050}

    def map_build_to_vintage(y):
        if y < 2020:
            return 1950
        elif y < 2030:
            return 2020
        elif y < 2050:
            return 2030
        else:
            return 2050

    inst = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")
    df = inst.copy()
    df["Primary_fuel"] = df["Primary_fuel"].astype(str).str.strip()
    df = df[df["Primary_fuel"].isin(["Natural gas", "Gas", "CH4"])].copy()
    if df.empty:
        print("Gas: no matching fuel rows.")
        return

    df["Year_built"] = pd.to_numeric(df["Year_built"], errors="coerce")
    df["Capacity_MW"] = pd.to_numeric(df["Capacity_MW"], errors="coerce").fillna(0.0)
    df = df.dropna(subset=["Year_built"]).copy()
    df["Year_built"] = df["Year_built"].astype(int)
    df = df[df["Capacity_MW"] > 0].copy()
    df["Node"] = df["Node"].astype(str).str.strip()
    df = df[df["Node"].isin(nodes)].copy()

    def map_gt_type(x):
        s = str(x).upper()
        if "CCGT" in s:
            return "CCGT"
        if "OCGT" in s:
            return "OCGT"
        if "GT" in s:
            return "GT"
        return "GT"

    if "Techs" in df.columns:
        df["techs"] = df["Techs"].apply(map_gt_type)
    else:
        df["techs"] = "GT"

    df = df[df["techs"].isin(gt_techs)].copy()
    if df.empty:
        print("Gas: no GT/CCGT/OCGT units after filters.")
        return

    # Map build years to vintages
    df["vintage_year"] = df["Year_built"].apply(map_build_to_vintage)

    if debug:
        print("[GT] years_sel:", years_sel, "base_year:", base_year)
        print("[GT] brownfield capacity by tech, vintage:")
        print(df.groupby(["techs", "vintage_year"])["Capacity_MW"].sum())

    # DEA tables
    eff_table = {
        "GT":   {2015: 0.36, 2020: 0.37, 2030: 0.39, 2050: 0.40},
        "CCGT": {2015: 0.50, 2020: 0.51, 2030: 0.53, 2050: 0.55},
        "OCGT": {2015: 0.41, 2020: 0.42, 2030: 0.43, 2050: 0.45},
    }
    lifetime_table = {
        "GT":   {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "CCGT": {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "OCGT": {2015: 25, 2020: 25, 2030: 25, 2050: 25},
    }
    invest_meur_per_gw = {
        "GT":   {2015: 797.5, 2020: 776.3, 2030: 744.4, 2050: 723.1},
        "CCGT": {2015: 1382.4, 2020: 1382.4, 2030: 1276.0, 2050: 1169.7},
        "OCGT": {2015: 499.8, 2020: 478.5, 2030: 467.9, 2050: 436.0},
    }
    omfix_frac = {
        "GT":   {2015: 0.027, 2020: 0.027, 2030: 0.027, 2050: 0.026},
        "CCGT": {2015: 0.023, 2020: 0.023, 2030: 0.023, 2050: 0.024},
        "OCGT": {2015: 0.017, 2020: 0.018, 2030: 0.018, 2050: 0.018},
    }
    kgCO2_per_MWh_fuel = 204.8
    fuel_cost_meur_per_gwh_fuel = 0.045
    fuel = "CH4"

    def eff_for_tech_vintage(t, v):
        dy = dea_year_for_vintage_map[v]
        return float(eff_table[t][dy])

    def lifetime_for_tech_vintage(t, v):
        dy = dea_year_for_vintage_map[v]
        return int(lifetime_table[t][dy])

    def capex_for_tech_year(t, y):
        candidates = [2015, 2020, 2030, 2050]
        dy = min(candidates, key=lambda d: abs(d - y))
        return float(invest_meur_per_gw[t][dy])

    def omfix_for_tech_year(t, y):
        candidates = [2015, 2020, 2030, 2050]
        dy = min(candidates, key=lambda d: abs(d - y))
        capex = float(invest_meur_per_gw[t][dy])
        frac = float(omfix_frac[t][dy])
        return capex * frac

    # converter_techparam
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product([gt_techs, vintages],
                                         names=["techs", "years"])
    )
    techparam["activityUpperLimit"] = 1.0
    for t in gt_techs:
        for v in vintages:
            techparam.loc[(t, v), "lifeTime"] = lifetime_for_tech_vintage(t, v)
    m["Base"].parameter.add(techparam, "converter_techparam")

    # converter_capacityparam
    cap_build = (
        df.groupby(["Node", "vintage_year", "techs"])["Capacity_MW"]
          .sum()
          .rename("unitsBuild")
          .to_frame()
          .div(1e3)
    )
    cap_build.index = cap_build.index.set_names(["nodesdata", "years", "techs"])

    cap_bounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_sel, gt_techs],
                                         names=["nodesdata", "years", "techs"])
    )
    cap_bounds["unitsUpperLimit"] = 100.0
    cap_bounds["noExpansion"] = 0
    if base_year in years_sel:
        cap_bounds.loc[idx[:, [base_year], :], "noExpansion"] = 1

    cap_full = pd.concat([cap_build, cap_bounds], axis=1).sort_index()
    m["Base"].parameter.add(cap_full, "converter_capacityparam")

    # converter_coefficient
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([gt_techs, vintages, activities, ["Elec", fuel]],
                                         names=["techs", "years", "activities", "commodities"])
    )
    coef["coefficient"] = 0.0
    coef.loc[idx[:, :, "Powergen", "Elec"], "coefficient"] = 1.0
    for t in gt_techs:
        for v in vintages:
            eff = eff_for_tech_vintage(t, v)
            coef.loc[idx[t, v, "Powergen", fuel], "coefficient"] = -1.0 / eff
    m["Base"].parameter.add(coef, "converter_coefficient")

    # accounting_converterunits
    acc_units = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], gt_techs, years_sel],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()
    for t in gt_techs:
        for y in years_sel:
            acc_units.loc[("Invest", "global", "horizon", t, y), "perUnitBuild"] = capex_for_tech_year(t, y)
            acc_units.loc[("Invest", "global", "horizon", t, y), "useAnnuity"] = 1
            acc_units.loc[("Invest", "global", "horizon", t, y), "amorTime"] = lifetime_for_tech_vintage(t, 2020)
            acc_units.loc[("Invest", "global", "horizon", t, y), "interest"] = 0.05
            acc_units.loc[("OMFix", "global", "horizon", t, y), "perUnitTotal"] = omfix_for_tech_year(t, y)
    m["Base"].parameter.add(acc_units, "accounting_converterunits")

    # accounting_converteractivity (FuelCost, CO2)
    acc_act = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["FuelCost", "CO2_emission"], ["global"], ["horizon"], gt_techs, vintages, activities],
            names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
        )
    ).sort_index()
    acc_act["perActivity"] = 0.0
    for t in gt_techs:
        for v in vintages:
            eff = eff_for_tech_vintage(t, v)
            # Fuel cost: M€/GWh_el
            acc_act.loc[("FuelCost", "global", "horizon", t, v, "Powergen"), "perActivity"] = fuel_cost_meur_per_gwh_fuel / eff
            # CO2: kt/GWh_el
            kt_per_gwh_el = (kgCO2_per_MWh_fuel * 1e-3) / eff
            acc_act.loc[("CO2_emission", "global", "horizon", t, v, "Powergen"), "perActivity"] = kt_per_gwh_el
    m["Base"].parameter.add(acc_act, "accounting_converteractivity")

    # Fuel imports: FuelImport_CH4
    if add_fuel_imports:
        rows = [(n, y, "FuelImport_CH4", "CH4") for n in nodes for y in years_sel]
        ss_idx = pd.MultiIndex.from_tuples(
            rows, names=["nodesdata", "years", "sourcesinks", "commodities"]
        )
        ss = pd.DataFrame(index=ss_idx)
        ss["lower"] = 0.0
        ss["upper"] = 1e12
        m["Base"].parameter.add(ss[["lower", "upper"]], "sourcesink_annualsum")

        cfg = pd.DataFrame(index=ss_idx)
        cfg["usesLowerSum"] = 1
        cfg["usesUpperSum"] = 1
        m["Base"].parameter.add(cfg, "sourcesink_config")

        acc_rows = [("FuelCost", "global", y, "FuelImport_CH4", "CH4") for y in years_sel]
        acc_idx = pd.MultiIndex.from_tuples(
            acc_rows,
            names=["indicator", "regionscope", "years", "sourcesinks", "commodities"],
        )
        acc = pd.DataFrame(index=acc_idx)
        acc["perFlow"] = [fuel_cost_meur_per_gwh_fuel for _ in acc_rows]
        m["Base"].parameter.add(acc, "accounting_sourcesinkflow")

        try:
            m["Base"].set.add(["FuelImport_CH4"], "sourcesink_techs")
        except Exception:
            pass

    print("Gas turbines: added techs:", gt_techs, "vintages:", vintages)


# renewables

def add_hydro(m):
    """
    Hydro:
      - Loads Water_in inflows (fixed profile, positive supply)
      - Adds hydro turbines + reservoirs
      - Brownfield build year = 2000
      - Capacity rows exist for all yrs_to_calc (critical for noExpansion)
      - Adds an accounting penalty for hydro spill activity
    """

    import pandas as pd
    global yrs_to_calc, yrs_sel, inflow_file, idx

    print("\n--- ADDING HYDRO INFLOWS (Water_in) ---")

    # --------------------------------------------------
    # 1) LOAD HYDRO INFLOWS  (same logic as old add_demand)
    # --------------------------------------------------
    inflow_df = pd.read_csv(inflow_file)

    # Normalise column names (permissive)
    rename = {}
    if "Region" in inflow_df.columns:
        rename["Region"] = "node"
    if "region" in inflow_df.columns:
        rename["region"] = "node"
    if "Year" in inflow_df.columns:
        rename["Year"] = "year"
    if "Sector" in inflow_df.columns:
        rename["Sector"] = "sector"
    if "Carrier" in inflow_df.columns:
        rename["Carrier"] = "commodity"
    if "carrier" in inflow_df.columns:
        rename["carrier"] = "commodity"

    if rename:
        inflow_df = inflow_df.rename(columns=rename)

    # Required base columns
    required = {"node", "year"}
    missing = required - set(inflow_df.columns)
    if missing:
        raise ValueError(
            f"Hydro inflow CSV missing required columns {missing}. "
            f"Found columns: {list(inflow_df.columns)}"
        )

    # Defaults if not present
    if "sector" not in inflow_df.columns:
        inflow_df["sector"] = "All"
    if "commodity" not in inflow_df.columns:
        inflow_df["commodity"] = "HydroInflow"

    # Map regions to REMix codes
    inflow_df["node"] = inflow_df["node"].apply(
        lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r)
    )
    inflow_df = inflow_df.dropna(subset=["node"])

    inflow_df["year"] = inflow_df["year"].astype(int)

    # Keep ONLY optimisation years (do not extend years set here)
    inflow_df = inflow_df.loc[inflow_df["year"].isin(yrs_sel)].copy()

    # Rename HydroInflow → Water_in commodity
    inflow_df["commodity"] = inflow_df["commodity"].astype(str).str.strip()
    inflow_df.loc[inflow_df["commodity"] == "HydroInflow", "commodity"] = "Water_in"

    # Hour columns t0001..t8760
    hour_cols = [c for c in inflow_df.columns if str(c).startswith("t")]
    if not hour_cols:
        raise ValueError("No hourly columns found in hydro inflow file (t0001..t8760).")

    # Build profile (positive supply)
    inflow = inflow_df.set_index(["node", "year", "sector", "commodity"])[hour_cols]

    #  inflow is negative → flip sign to make it positive supply
    inflow *= -1

    inflow["type"] = "fixed"
    inflow_fixed = inflow.set_index("type", append=True).round(3)

    m["Base"].profile.add(inflow_fixed, "sourcesink_profile")

    inflow_cfg = pd.DataFrame(index=inflow.index)
    inflow_cfg["usesFixedProfile"] = 1
    # Keep only rows that actually have nonzero inflow
    inflow_cfg = inflow_cfg.loc[inflow.select_dtypes(include="number").sum(axis=1) != 0]

    m["Base"].parameter.add(inflow_cfg, "sourcesink_config")

    # Ensure nodesdata includes all hydro inflow nodes
    m["Base"].set.add(sorted(inflow.index.get_level_values(0).unique()), "nodesdata")

    print(f"Hydro inflows loaded for years={sorted(inflow_df['year'].unique().tolist())}.")

    # --------------------------------------------------
    # 2) HYDRO PLANT (converter + reservoir)
    # --------------------------------------------------
    print("\n--- ADDING HYDRO PLANT (converter + reservoir) ---")

    hydro_vintage = [1950]
    # brownfield build year 2000 + all data years (so capacity rows exist for all)
    hydro_years = [2000] + list(yrs_to_calc)
    hydro_nodes = ["BOP", "CAN", "CEN", "HBY", "NEL", "OTG", "WTO"]
    hydro_techs = ["Hydro"]
    hydro_activities = ["Powergen", "Spill"]

    # Ocean sink for Water_out (allows water to leave the system)
    ocean_cfg = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [m["Base"].set.nodesdata, m["Base"].set.yearssel, ["Ocean"], ["Water_out"]]
        )
    )
    ocean_cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(ocean_cfg, "sourcesink_config")

    # Turbine tech parameters (lifetime, availability)
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_techs, hydro_vintage])
    )
    conv_tech["lifeTime"] = 100
    conv_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # Turbine capacity (GW) with brownfield unitsBuild in 2000
    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, hydro_techs]),
        columns=["unitsBuild", "unitsUpperLimit", "noExpansion"],
    )
    conv_cap.loc[idx["BOP", 2000, "Hydro"], "unitsBuild"] = 0.17095
    conv_cap.loc[idx["CAN", 2000, "Hydro"], "unitsBuild"] = 1.82683
    conv_cap.loc[idx["CEN", 2000, "Hydro"], "unitsBuild"] = 0.399
    conv_cap.loc[idx["HBY", 2000, "Hydro"], "unitsBuild"] = 0.1422
    conv_cap.loc[idx["NEL", 2000, "Hydro"], "unitsBuild"] = 0.0453
    conv_cap.loc[idx["OTG", 2000, "Hydro"], "unitsBuild"] = 1.664
    conv_cap.loc[idx["WTO", 2000, "Hydro"], "unitsBuild"] = 1.0873
    # No new hydro builds in any year (pure brownfield)
    conv_cap["noExpansion"] = 1
    # Large upper limit to avoid artificial cap beyond unitsBuild
    conv_cap["unitsUpperLimit"] = 100
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # Turbine coefficients (Water_in, Water_out, Elec)
    hydro_eff = 0.95
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [hydro_techs, hydro_vintage, hydro_activities, ["Water_in", "Water_out", "Elec"]]
        )
    )
    # Powergen activity: convert Water_in to Elec and Water_out
    conv_coef.loc[idx[:, :, "Powergen", "Elec"], "coefficient"] = 1.0
    conv_coef.loc[idx[:, :, "Powergen", "Water_in"], "coefficient"] = -hydro_eff
    conv_coef.loc[idx[:, :, "Powergen", "Water_out"], "coefficient"] = hydro_eff
    # Spill activity: moves water from reservoir in to reservoir out, no electricity
    conv_coef.loc[idx[:, :, "Spill", "Water_in"], "coefficient"] = -100.0
    conv_coef.loc[idx[:, :, "Spill", "Water_out"], "coefficient"] = 100.0
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # Reservoir tech: stores Water_in between periods
    stor_techs = ["Hydro_reservoir"]
    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([stor_techs, hydro_vintage])
    )
    stor_tech["lifeTime"] = 100
    stor_tech["levelUpperLimit"] = 1
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([stor_techs, hydro_vintage, ["Water_in"]])
    )
    stor_size["size"] = 1
    stor_size["selfdischarge"] = 0
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # Reservoir capacities (unitsBuild, unitsUpperLimit, noExpansion)
    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, stor_techs]),
        columns=["unitsBuild", "unitsUpperLimit", "noExpansion"],
    )
    stor_res.loc[idx["CAN", 2000, :], "unitsBuild"] = 2517.2429
    stor_res.loc[idx["HBY", 2000, :], "unitsBuild"] = 154.2635
    stor_res.loc[idx["OTG", 2000, :], "unitsBuild"] = 729.5595
    stor_res.loc[idx["WTO", 2000, :], "unitsBuild"] = 587.1371
    stor_res["noExpansion"] = 1
    stor_res["unitsUpperLimit"] = 3000
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    # --------------------------------------------------
    # 3) Spill penalty in accounting_converteractivity
    # --------------------------------------------------
    # We treat spill as an accounting indicator so the model avoids
    # unnecessary spill but can still spill when there is no better option.
    #
    # Example ranges for penalty_per_gwh (M€/GWh activity):
    #   0.001  -> very soft penalty (~1 €/MWh), mostly opportunity cost driven.
    #   0.01   -> moderate penalty (~10 €/MWh), discourages spill strongly.
    #   0.05+  -> very high penalty, spill only when absolutely unavoidable.
    penalty_per_gwh = 0.01  # adjust as needed

    years_sel_int = sorted(int(y) for y in m["Base"].set.yearssel)

    spill_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["SpillPenalty"], ["global"], ["horizon"], hydro_techs, years_sel_int, ["Spill"]],
            names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
        )
    )
    spill_acc["perActivity"] = float(penalty_per_gwh)

    m["Base"].parameter.add(spill_acc, "accounting_converteractivity")
    print(f"Hydro spill penalty added: {penalty_per_gwh} M€/GWh activity.")

    print("Hydro inflows + plant added.")

def add_geothermal(m):
    print("\n--- ADDING GEOTHERMAL ---")

    df = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))

    # Optimisation years
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    base_year = years_sel[0]

    global yrs_to_calc
    # Ensure 2000 (geothermal vintage) is part of yrs_to_calc so brownfield build is visible
    yrs_to_calc = sorted(set(yrs_to_calc) | {2000})
    m["Base"].set.add(yrs_to_calc, "years")
    years_data = list(yrs_to_calc)

    nodes = sorted(m["Base"].set.nodesdata)
    techs = ["Geothermal"]
    acts = ["Powergen"]


    # Tech params: treat geothermal as a single vintage (2000)
    vintage_year = 2000
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, [vintage_year]],
            names=["techs", "vintages"],
        )
    )
    techparam["lifeTime"] = 100
    techparam["activityUpperLimit"] = 1
    m["Base"].parameter.add(techparam, "converter_techparam")

    # Brownfield capacity: all existing geothermal aggregated and built in 2000
    df_geo = df[df["Type"] == "Geothermal"].copy()
    if df_geo.empty:
        print("No geothermal rows found in brownfield database.")
        return

    df_geo["Capacity_MW"] = pd.to_numeric(df_geo["Capacity_MW"], errors="coerce").fillna(0.0)
    df_geo = df_geo[df_geo["Capacity_MW"] > 0]

    cap_by_node = df_geo.groupby("Node")["Capacity_MW"].sum().div(1000.0)  # MW → GW

    cap_idx = pd.MultiIndex.from_product(
        [nodes, years_data, techs],
        names=["nodesdata", "years", "techs"],
    )
    cap = pd.DataFrame(index=cap_idx)
    cap["noExpansion"] = 1
    cap["unitsUpperLimit"] = 0.0

    for n in nodes:
        gw = float(cap_by_node.get(n, 0.0))
        if gw <= 0:
            continue

        # All geothermal built in vintage_year
        if vintage_year in years_data:
            cap.loc[(n, vintage_year, "Geothermal"), "unitsBuild"] = gw

        # Same usable capacity limit in every data year
        cap.loc[(n, years_data, "Geothermal"), "unitsUpperLimit"] = gw

    cap = cap.dropna(how="all")
    m["Base"].parameter.add(cap.sort_index(), "converter_capacityparam")

    # Conversion coefficients: Elec output only (no explicit fuel)
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, [vintage_year], acts, ["Elec"]],
            names=["techs", "vintages", "activities", "commodities"],
        )
    )
    coef["coefficient"] = 1.0
    m["Base"].parameter.add(coef, "converter_coefficient")

    print("Geothermal installed as brownfield in 2000 and usable in all model years.")

def load_feedin_csv(weather_year: int = 2012):
    """
    Load renewables feed-in from timeseries_2012_w_corr.csv.

    Returns
    -------
    feed : DataFrame
        Indexed by (nodesdata, techs, t_model) with one column 'value'.

    Behaviour
    ---------
    - Reads raw timeseries (t, region, technology, timeseries_per_region).
    - Parses timestamps (dayfirst, NZ style).
    - Applies +30 minutes shift (REMix convention).
    - Keeps ONLY rows whose shifted timestamp is in 'weather_year'.
    - Computes integer t_model = 1..8760.
    - For each (node, tech), reindexes to exactly 1..8760.
      * Leading missing hours (typically 1..11) are now WRAPPED from
        the last hours of the year instead of set to 0.
      * All other missing hours (internal gaps) are still set to 0.0.
    """


    ts_path = Path(path_profiles) / "timeseries_2012_w_corr.csv"

    if not ts_path.exists():
        raise FileNotFoundError(f"Timeseries file not found: {ts_path}")

    # Read raw file
    ts_raw = pd.read_csv(ts_path)

    # Check required columns
    required = {"t", "region", "technology", "timeseries_per_region"}
    missing = required - set(ts_raw.columns)
    if missing:
        raise ValueError(
            f"[load_feedin_csv] Timeseries missing columns {missing}. "
            f"Found: {list(ts_raw.columns)}"
        )

    # Parse timestamps (NZ / dayfirst)
    dt = pd.to_datetime(ts_raw["t"], dayfirst=True, errors="coerce")
    bad = int(dt.isna().sum())
    if bad:
        print("[load_feedin_csv] bad sample:",
              ts_raw.loc[dt.isna(), "t"].astype(str).head(5).tolist())
        raise ValueError("Some timestamps could not be parsed (check date format).")

    # Shift by +30 min (REMix convention)
    dt_shift = dt + pd.Timedelta(minutes=30)

    # Keep only rows with shifted year == weather_year
    mask = (dt_shift.dt.year == int(weather_year))
    kept = int(mask.sum())

    ts = ts_raw.loc[mask].copy()
    dt_shift = dt_shift.loc[mask]

    # Numeric MW values
    ts["value"] = pd.to_numeric(ts["timeseries_per_region"], errors="coerce")
    bad_val = int(ts["value"].isna().sum())
    if bad_val:
        print("[load_feedin_csv] NaN values sample:",
              ts.loc[ts["value"].isna()].head(5))
        raise ValueError("Non-numeric values in timeseries_per_region after coercion.")

    # Compute t_model = 1..8760
    day_of_year = dt_shift.dt.dayofyear  # 1..365
    hour_of_day = dt_shift.dt.hour       # 0..23
    t_model = (day_of_year - 1) * 24 + hour_of_day + 1
    ts["t_model"] = t_model.astype(int)

    # Rename to match REMix sets
    ts = ts.rename(columns={"region": "nodesdata", "technology": "techs"})
    ts = ts[["nodesdata", "techs", "t_model", "value"]]

    # Collapse duplicates
    ts = (
        ts.groupby(["nodesdata", "techs", "t_model"], as_index=False, sort=False)["value"]
        .sum()
    )

    # Enforce exactly 1..8760 hours per (node,tech)
    all_hours = pd.Index(range(1, 8761), name="t_model")
    frames = []

    # For diagnostics: track how many hours we wrap vs keep as zero
    wrap_stats = []

    for (node, tech), g in ts.groupby(["nodesdata", "techs"], sort=False):
        s = g.set_index("t_model")["value"]
        s_full = s.reindex(all_hours, fill_value=0.0)

        observed_hours = sorted(g["t_model"].unique().tolist())
        first_obs = observed_hours[0]
        last_obs = observed_hours[-1]

        # Wrap leading missing hours: copy last N hours of the year
        wrapped = s_full.copy()
        if first_obs > 1:
            n_lead = first_obs - 1
            donor = s_full.iloc[-n_lead:].values
            if len(donor) == n_lead:
                wrapped.iloc[0:n_lead] = donor
            else:
                print(f"[load_feedin_csv] WARNING: not enough donor hours "
                      f"to wrap for (node={node}, tech={tech}), leaving zeros.")
        else:
            n_lead = 0

        # Count internal zeros (inside [first_obs, last_obs]) after wrapping
        internal_mask = (wrapped.index >= first_obs) & (wrapped.index <= last_obs)
        n_internal_zeros = int((wrapped[internal_mask] == 0.0).sum())

        wrap_stats.append(
            {
                "nodesdata": node,
                "techs": tech,
                "first_obs": first_obs,
                "last_obs": last_obs,
                "n_leading_wrapped": n_lead,
                "n_internal_zeros": n_internal_zeros,
            }
        )

        df_full = wrapped.to_frame("value")
        df_full["nodesdata"] = node
        df_full["techs"] = tech

        df_full = df_full.set_index(["nodesdata", "techs"], append=True)
        df_full = df_full.reorder_levels(["nodesdata", "techs", "t_model"])
        frames.append(df_full)

    feed = pd.concat(frames).sort_index()

    # Sanity prints
    counts = feed.groupby(["nodesdata", "techs"]).size()
    # print("[load_feedin_csv] per (node,tech) hour count min/max:",
    #       counts.min(), counts.max())
    # print("[load_feedin_csv] t_model min/max:",
    #       int(feed.index.get_level_values("t_model").min()),
    #       int(feed.index.get_level_values("t_model").max()))

    # Quick summary of wrap behaviour
    wrap_df = pd.DataFrame(wrap_stats)
    # print("\n[load_feedin_csv] WRAP SUMMARY (per node,tech):")
    # print("  Leading hours wrapped (should be mostly 11):")
    # print(wrap_df["n_leading_wrapped"].describe().to_string())
    # print("  Internal zeros remaining after wrap (still potential gaps):")
    # print(wrap_df["n_internal_zeros"].describe().to_string())

    # # Optional: write detailed wrap stats for offline inspection
    # wrap_out = Path(path_profiles) / "feedin_wrap_stats.csv"
    # wrap_df.to_csv(wrap_out, index=False)
    # print(f"[load_feedin_csv] Detailed wrap stats written to: {wrap_out}")

    return feed

def add_renewables(m):
    """
    Add PV and wind renewables with:
      - Legacy tech names: pv_central_fixed, pv_decentral,
                           wind_onshore, wind_offshore_floating, wind_offshore_foundation
      - Installables from region_statistics_2012_w_corr.csv (wind buckets mapped -> legacy names)
      - Feed-in from timeseries_2012_w_corr.csv (legacy tech names already)
      - Activity profiles from weather year 2012 replicated to all yrs_to_calc
      - Investment + OMFix costs based on DEA 2022 vintage data, mapped to model years
    """

    global yrs_to_calc, yrs_sel, idx

    print("\n--- ADDING RENEWABLES (PV + WIND) ---")

    # ------------------------------------------------------------------
    # 1. Years / nodes
    # ------------------------------------------------------------------
    # yearssel: optimisation years from the instance (e.g. [2020, 2050])
    yearssel = sorted(int(y) for y in m["Base"].set.yearssel)
    # years_data: all model years we want data/costs for (global yrs_to_calc)
    years_data = sorted(int(y) for y in yrs_to_calc)
    base_year = int(yearssel[0])

    # Ensure the REMix "years" set contains all data years (numeric only)
    m["Base"].set.add(years_data, "years")

    nodes = list(m["Base"].set.nodesdata)

    # ------------------------------------------------------------------
    # 2. Tech sets 
    # ------------------------------------------------------------------
    pv_techs = ["pv_central_fixed", "pv_decentral"]
    wind_techs = ["wind_onshore", "wind_offshore_floating", "wind_offshore_foundation"]
    model_techs = pv_techs + wind_techs


    # ------------------------------------------------------------------
    # 3. Installables (region_statistics_2012_w_corr.csv)
    # ------------------------------------------------------------------
    inst_path = Path(path_profiles) / "region_statistics_2012_w_corr.csv"
    inst_raw = pd.read_csv(inst_path)

    # Required columns: region, technology, installable_per_region
    need_inst = {"region", "technology", "installable_per_region"}
    missing = need_inst - set(inst_raw.columns)
    if missing:
        raise ValueError(
            f"Installables file missing {missing}. Found {list(inst_raw.columns)}"
        )

    # Numeric installables (MW)
    inst_raw["installable_per_region"] = pd.to_numeric(
        inst_raw["installable_per_region"], errors="coerce"
    )
    if inst_raw["installable_per_region"].isna().any():
        raise ValueError("Non-numeric installable_per_region in installables file.")

    # Rename to match REMix sets: nodesdata / techs
    inst = inst_raw.rename(columns={"region": "nodesdata", "technology": "techs"})

    # Wind buckets -> legacy techs (only in statistics, not in timeseries)
    wind_bucket_map = {
        "wind_onshore_1": "wind_onshore",
        "wind_onshore_2": "wind_offshore_floating",
        "wind_onshore_4": "wind_offshore_foundation",
    }

    # PV rows: names already match model techs
    inst_pv = inst[inst["techs"].isin(pv_techs)].copy()

    # Wind rows: bucketed names mapped to legacy names
    inst_wind = inst[inst["techs"].isin(wind_bucket_map.keys())].copy()
    inst_wind["techs"] = inst_wind["techs"].map(wind_bucket_map)

    # Combine and aggregate by (nodesdata, techs)
    inst_model = pd.concat([inst_pv, inst_wind], ignore_index=True)
    inst_model = (
        inst_model.groupby(["nodesdata", "techs"], as_index=False)["installable_per_region"]
        .sum()
    )

    # Convert MW -> GW for capacity limits (unitsUpperLimit)
    inst_model["unitsUpperLimit"] = inst_model["installable_per_region"] / 1000.0
    inst_model = inst_model.set_index(["nodesdata", "techs"])[["unitsUpperLimit"]].sort_index()

    # ------------------------------------------------------------------
    # 4. converter_techparam for renewables (techs x years_data)
    # ------------------------------------------------------------------
    techparam_index = pd.MultiIndex.from_product(
        [model_techs, years_data],
        names=["techs", "years"],
    )

    techparam = pd.DataFrame(index=techparam_index)
    # Activity upper limit is 0 here; hourly availability is imposed via profile
    techparam["activityUpperLimit"] = 0.0
    techparam["lifeTime"] = np.nan

    # Lifetimes by year: PV 35/40, wind 27/30
    techparam.loc[idx[pv_techs, [y for y in years_data if y <= 2020]], "lifeTime"] = 35
    techparam.loc[idx[pv_techs, [y for y in years_data if y >  2020]], "lifeTime"] = 40

    techparam.loc[idx[wind_techs, [y for y in years_data if y <= 2020]], "lifeTime"] = 27
    techparam.loc[idx[wind_techs, [y for y in years_data if y >  2020]], "lifeTime"] = 30

    m["Base"].parameter.add(techparam, "converter_techparam")

    # ------------------------------------------------------------------
    # 5. converter_capacityparam for renewables (nodesdata, years, techs)
    # ------------------------------------------------------------------
    cap_index = pd.MultiIndex.from_product(
        [nodes, years_data, model_techs],
        names=["nodesdata", "years", "techs"],
    )
    cap_df = pd.DataFrame(index=cap_index)

    # Join installable upper limits (GW) onto capacity table
    cap_df = cap_df.join(inst_model, on=["nodesdata", "techs"], how="left")

    years_build = list(range(1989, 2021))  # brownfield build years
    re_nodes = nodes  # same node set
    re_techs = model_techs

    re_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product([re_nodes, years_build, re_techs]),
        columns=["unitsBuild"],
    )
    re_caps.index.names = ["nodesdata", "years", "techs"]

    # wind_onshore unitsBuild entries:
    re_caps.loc[idx["CAN", [2003], "wind_onshore"], "unitsBuild"] = 0.0005
    re_caps.loc[idx["CAN", [2005], "wind_onshore"], "unitsBuild"] = 0.0001

    re_caps.loc[idx["CEN", [1999], "wind_onshore"], "unitsBuild"] = 31.7  / 1000
    re_caps.loc[idx["CEN", [2004], "wind_onshore"], "unitsBuild"] = 127.05 / 1000
    re_caps.loc[idx["CEN", [2007], "wind_onshore"], "unitsBuild"] = 93    / 1000
    re_caps.loc[idx["CEN", [2011], "wind_onshore"], "unitsBuild"] = 48.5  / 1000
    re_caps.loc[idx["CEN", [2020], "wind_onshore"], "unitsBuild"] = 221.4 / 1000

    re_caps.loc[idx["OTG", [2007], "wind_onshore"], "unitsBuild"] = 58    / 1000
    re_caps.loc[idx["OTG", [2009], "wind_onshore"], "unitsBuild"] = 2.25  / 1000
    re_caps.loc[idx["OTG", [2010], "wind_onshore"], "unitsBuild"] = 0.45  / 1000
    re_caps.loc[idx["OTG", [2011], "wind_onshore"], "unitsBuild"] = 43.65 / 1000
    re_caps.loc[idx["OTG", [2015], "wind_onshore"], "unitsBuild"] = 6.8   / 1000

    re_caps.loc[idx["NEL", [2010], "wind_onshore"], "unitsBuild"] = 0.75 / 1000
    re_caps.loc[idx["NEL", [2011], "wind_onshore"], "unitsBuild"] = 1.0  / 1000
    re_caps.loc[idx["NEL", [2014], "wind_onshore"], "unitsBuild"] = 0.66 / 1000

    re_caps.loc[idx["TRN", [2020], "wind_onshore"], "unitsBuild"] = 0.1333

    re_caps.loc[idx["WEL", [1993], "wind_onshore"], "unitsBuild"] = 8.45  / 1000
    re_caps.loc[idx["WEL", [1996], "wind_onshore"], "unitsBuild"] = 8.45  / 1000
    re_caps.loc[idx["WEL", [2009], "wind_onshore"], "unitsBuild"] = 143   / 1000
    re_caps.loc[idx["WEL", [2014], "wind_onshore"], "unitsBuild"] = 71.3  / 1000
    re_caps.loc[idx["WTO", [2011], "wind_onshore"], "unitsBuild"] = 64.4  / 1000

    # align brownfield frame to cap_df index and merge in unitsBuild
    re_caps = re_caps.dropna(how="all")
    cap_df = cap_df.merge(
        re_caps,
        left_index=True,
        right_index=True,
        how="left",
    )

    # No expansion allowed in base year (e.g. 2020)
    cap_df.loc[idx[:, [base_year], :], "noExpansion"] = 1

    # Any (node,year,tech) without installable data gets unitsUpperLimit=0
    missing_ul = int(cap_df["unitsUpperLimit"].isna().sum())
    cap_df["unitsUpperLimit"] = cap_df["unitsUpperLimit"].fillna(0.0)

    # allow at least existing brownfield capacity in the base year
    existing_cap = (
        re_caps["unitsBuild"]
        .groupby([ "nodesdata", "techs" ])
        .sum()
    )
    for (n, t), cap in existing_cap.items():
        if pd.notna(cap) and cap > 0:
            if (n, base_year, t) in cap_df.index:
                if cap_df.loc[(n, base_year, t), "unitsUpperLimit"] < cap:
                    cap_df.loc[(n, base_year, t), "unitsUpperLimit"] = cap


    # Any (node,year,tech) without installable data gets unitsUpperLimit=0
    missing_ul = int(cap_df["unitsUpperLimit"].isna().sum())
    cap_df["unitsUpperLimit"] = cap_df["unitsUpperLimit"].fillna(0.0)
    m["Base"].parameter.add(cap_df, "converter_capacityparam")

    # ------------------------------------------------------------------
    # 6. converter_coefficient: Powergen -> Elec = 1.0
    # ------------------------------------------------------------------
    coef_index = pd.MultiIndex.from_product(
        [model_techs, years_data, ["Powergen"], ["Elec"]],
        names=["techs", "years", "activities", "commodities"],
    )
    coef_df = pd.DataFrame(index=coef_index)
    coef_df["coefficient"] = 1.0  # 1 unit activity -> 1 unit Elec
    m["Base"].parameter.add(coef_df, "converter_coefficient")

    # ------------------------------------------------------------------
    # 7. Activity profile: load_feedin_csv, normalise by installables,
    #    enforce 8760 hours, replicate to years_data.
    # ------------------------------------------------------------------
    feed = load_feedin_csv(weather_year=2012)  # (nodesdata, techs, t_model) -> MW

    # Keep only those techs we actually model
    available_techs = sorted(feed.index.get_level_values("techs").unique().tolist())
    keep_techs = [t for t in model_techs if t in available_techs]
    missing_ts = [t for t in model_techs if t not in available_techs]
    if missing_ts:
        print("[renewables] WARNING: techs missing in timeseries:", missing_ts)
    feed = feed.loc[idx[:, keep_techs, :], :]

    # Denominator: MW installables (convert from GW back to MW)
    denom_mw = (inst_model["unitsUpperLimit"] * 1000.0).rename("denom_mw")

    # Merge feed with denom by (nodesdata, techs)
    fr = feed.reset_index()  # columns: nodesdata, techs, t_model, value
    fr["denom_mw"] = fr.set_index(["nodesdata", "techs"]).index.map(denom_mw)

    # If any (node,tech) pair is missing installables, warn and set denom=1
    if fr["denom_mw"].isna().any():
        print("[renewables] WARNING: missing installables for some (node,tech). "
              "Filling denom_mw=1.0 there.")
        fr["denom_mw"] = fr["denom_mw"].fillna(1.0)

    # Availability = MW_profile / MW_installable, clipped to [0,1]
    fr["avail"] = (fr["value"] / fr["denom_mw"]).clip(lower=0.0)
    fr.loc[fr["avail"] > 1.0, "avail"] = 1.0
    fr["avail"] = fr["avail"].round(3)

    # Pivot to wide format: index (nodesdata, techs), columns = t_model
    avail_wide = fr.pivot_table(
        index=["nodesdata", "techs"],
        columns="t_model",
        values="avail",
        aggfunc="mean",
        fill_value=0.0,
    )

    # Ensure exactly 1..8760 hour columns
    avail_wide = avail_wide.reindex(columns=range(1, 8761), fill_value=0.0)

    # Rename columns to t0001..t8760
    avail_wide.columns = [f"t{str(t).zfill(4)}" for t in range(1, 8761)]

    # Replicate 2012 profile to each data year, add "type" = "upper"
    frames = []
    for y in years_data:
        tmp = avail_wide.copy()
        tmp["years"] = int(y)
        tmp["type"] = "upper"
        tmp = tmp.reset_index().set_index(["nodesdata", "years", "techs", "type"])
        frames.append(tmp)

    profile = pd.concat(frames).sort_index()

    m["Base"].profile.add(profile, "converter_activityprofile")

    # ------------------------------------------------------------------
    # 8. accounting_converterunits: CAPEX + OMFix (DEA-based) for PV & wind
    # ------------------------------------------------------------------
    # Index: (indicator, regionscope, timescope, techs, years)
    acc_index = pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["horizon"], model_techs, years_data],
        names=["indicator", "regionscope", "timescope", "techs", "years"],
    )
    re_acc = pd.DataFrame(index=acc_index)

    # Vintage years used in the original script to define cost trajectories
    vintages = [1950, 2030, 2040, 2050]

    # DEA 2022 investment costs per vintage (M€/GW) from old add_renewables
    pv_dec_cost = [870, 570, 460, 410]
    pv_cen_cost = [560, 380, 320, 290]
    w_on_cost   = [1330, 1040, 980, 960]
    w_off_fix   = [2120, 2287, 2168, 2130]
    w_off_float = (np.array(w_off_fix) * 1.2).tolist()  # +20% for floating

    # Helper: map each model year to the nearest cost vintage
    def nearest_vintage_idx(y):
        return min(range(len(vintages)), key=lambda i: abs(y - vintages[i]))

    # Fill Invest costs per (tech, year) using nearest vintage
    for y in years_data:
        i = nearest_vintage_idx(y)
        re_acc.loc[("Invest", "global", "horizon", "pv_decentral",            y), "perUnitBuild"] = pv_dec_cost[i]
        re_acc.loc[("Invest", "global", "horizon", "pv_central_fixed",        y), "perUnitBuild"] = pv_cen_cost[i]
        re_acc.loc[("Invest", "global", "horizon", "wind_onshore",            y), "perUnitBuild"] = w_on_cost[i]
        re_acc.loc[("Invest", "global", "horizon", "wind_offshore_foundation",y), "perUnitBuild"] = w_off_fix[i]
        re_acc.loc[("Invest", "global", "horizon", "wind_offshore_floating",  y), "perUnitBuild"] = w_off_float[i]

        # Lifetime for amortisation = same as converter_techparam logic
        if y <= 2020:
            life_pv = 35
            life_wind = 27
        else:
            life_pv = 40
            life_wind = 30

        re_acc.loc[("Invest", "global", "horizon", "pv_decentral",            y), "amorTime"] = life_pv
        re_acc.loc[("Invest", "global", "horizon", "pv_central_fixed",        y), "amorTime"] = life_pv
        re_acc.loc[("Invest", "global", "horizon", "wind_onshore",            y), "amorTime"] = life_wind
        re_acc.loc[("Invest", "global", "horizon", "wind_offshore_foundation",y), "amorTime"] = life_wind
        re_acc.loc[("Invest", "global", "horizon", "wind_offshore_floating",  y), "amorTime"] = life_wind

    # Common annuity settings for all Invest rows
    re_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
    re_acc.loc[idx["Invest", :, :, :, :], "interest"]   = 0.06

    # OMFix = 2% of investment per year (same as old script)
    invest_vals = re_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
        names=re_acc.index.names,
    )
    re_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.02

    # Drop all‑NaN rows
    re_acc = re_acc.replace(0, np.nan).dropna(how="all")
    m["Base"].parameter.add(re_acc, "accounting_converterunits")

    # ------------------------------------------------------------------
    # 9. Final check: Base.years must be numeric only
    # ------------------------------------------------------------------
    years_set = list(m["Base"].set.years)
    non_numeric = [y for y in years_set if not str(y).strip().isdigit()]
    # print("\n[renewables] Base.years sample:", years_set[:20])
    # print("[renewables] non-numeric years:", non_numeric[:20])

# elec storage

def add_lithium_batteries(m):
    """
    Add lithium-ion batteries with:
      - Tech: Battery
      - Power block (converter_* + accounting_converterunits)
      - Energy block (storage_* + accounting_storageunits)
    """

    import pandas as pd
    global idx

    print("\n--- ADDING LITHIUM BATTERIES ---")

    techs = ["Battery"]

    # Use all electricity nodes except LNG (consistent with other blocks)
    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]

    # Physics years: all model "years" in the instance (data + brownfield years)
    years = [int(y) for y in m["Base"].set.years]

    # Accounting years: optimisation years only (objective defined here)
    years_acc = [int(y) for y in m["Base"].set.yearssel]

    # Data years for cost trajectories (carry-forward stepwise)
    data_years = [2020, 2030, 2040, 2050]

    def prev_year(y):
        """
        Map any model year y to the latest data_year <= y.
        Used to carry-forward CAPEX/OMFix values.
        """
        y = int(y)
        earlier = [yy for yy in data_years if yy <= y]
        return earlier[-1] if earlier else data_years[0]

    # ------------------------------------------------------------------
    # 1. converter_techparam: Battery power block (techs x years)
    # ------------------------------------------------------------------
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, years], names=["techs", "years"])
    )
    conv_tech["lifeTime"] = 20           # inverter / power-block lifetime (years)
    conv_tech["activityUpperLimit"] = 1  # no extra availability constraint
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # ------------------------------------------------------------------
    # 2. converter_capacityparam: Battery power capacity (GW)
    # ------------------------------------------------------------------
    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [nodes, years, techs], names=["nodesData", "years", "techs"]
        )
    )
    conv_cap["unitsUpperLimit"] = 50.0   # max power capacity per node (GW_el)
    if 2020 in years:
        # Optional: disallow new power capacity in base year (pure brownfield)
        conv_cap.loc[idx[:, [2020], :], "noExpansion"] = 1
    
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # ------------------------------------------------------------------
    # 3. converter_coefficient: charge/discharge efficiencies
    # ------------------------------------------------------------------
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, years, ["Charge", "Discharge"], ["Elec", "Elec_battery"]],
            names=["techs", "years", "activities", "commodity"],
        )
    )
    # Charge: consume grid electricity, store in battery
    conv_coef.loc[idx["Battery", :, "Charge", "Elec"], "coefficient"] = -1.0
    conv_coef.loc[idx["Battery", :, "Charge", "Elec_battery"], "coefficient"] = 0.975
    # Discharge: release electricity from storage to grid
    conv_coef.loc[idx["Battery", :, "Discharge", "Elec"], "coefficient"] = 1.0
    conv_coef.loc[idx["Battery", :, "Discharge", "Elec_battery"], "coefficient"] = -1.025
    # Round-trip efficiency ≈ 0.975 / 1.025 ≈ 95%
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # ------------------------------------------------------------------
    # 4. accounting_converterunits: power part (inverter CAPEX/OMFix)
    # ------------------------------------------------------------------
    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, years_acc],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Power CAPEX (M€/GW_el) by data year, carried forward by prev_year()
    capex_power = {2020: 117, 2030: 55, 2040: 37, 2050: 30}

    for y in years_acc:
        yy = prev_year(y)
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", y], "perUnitBuild"] = capex_power[yy]
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", y], "amorTime"] = 20  # years

    conv_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "useAnnuity"] = 1
    conv_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "interest"] = 0.06

    # OMFix = 1.4% of power CAPEX per year 
    invest_vals = conv_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],
        names=conv_acc.index.names,
    )
    conv_acc.loc[idx["OMFix", "global", "horizon", "Battery", :], "perUnitTotal"] = invest_vals * 0.014
    m["Base"].parameter.add(conv_acc, "accounting_converterunits")

    # ------------------------------------------------------------------
    # 5. storage_techparam: energy block lifetime and max SOC
    # ------------------------------------------------------------------
    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, years], names=["techs", "years"])
    )
    stor_tech["lifeTime"] = 20          # energy-block lifetime (years)
    stor_tech["levelUpperLimit"] = 1.0  # SOC limited to 100% of design size
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    # ------------------------------------------------------------------
    # 6. storage_sizeparam: energy capacity per unit (4 h)
    # ------------------------------------------------------------------
    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, years, ["Elec_battery"]],
            names=["techs", "years", "commodity"],
        )
    )
    # Each Battery unit has 4 GWh of storage per GW of power (4-hour battery)
    stor_size.loc[idx["Battery", :, "Elec_battery"], "size"] = 4.0
    # No explicit self-discharge specified; default 0 in REMix
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # ------------------------------------------------------------------
    # 7. storage_reservoirparam: number of Battery units per node/year
    # ------------------------------------------------------------------
    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [nodes, years, techs],
            names=["nodesData", "years", "techs"],
        )
    )
    # Max number of Battery units per node (30 units => 30×4 GWh = 120 GWh energy)
    stor_res["unitsUpperLimit"] = 30.0
    if 2020 in years:
        # Optional: no new energy capacity in base year
        stor_res.loc[idx[:, [2020], :], "noExpansion"] = 1
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    # ------------------------------------------------------------------
    # 8. accounting_storageunits: energy part CAPEX/OMFix
    # ------------------------------------------------------------------
    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, years_acc],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Energy CAPEX (M€/GW_el of 4h storage):
    # original comment: [234, 110, 76, 61] * 4 because 4 GWh/unit
    capex_energy = {
        2020: 234 * 4,
        2030: 110 * 4,
        2040: 76  * 4,
        2050: 61  * 4,
    }

    for y in years_acc:
        yy = prev_year(y)
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", y], "perUnitBuild"] = capex_energy[yy]
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", y], "amorTime"] = 20  # years

    stor_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "useAnnuity"] = 1
    stor_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "interest"] = 0.06

    # OMFix = 1.4% of energy CAPEX per year
    invest_vals = stor_acc.loc[idx["Invest", "global", "horizon", "Battery", :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],
        names=stor_acc.index.names,
    )
    stor_acc.loc[idx["OMFix", "global", "horizon", "Battery", :], "perUnitTotal"] = invest_vals * 0.014
    m["Base"].parameter.add(stor_acc, "accounting_storageunits")

    print("Lithium batteries added (power + energy blocks).")

# hydrogen techs

def add_electrolyser(m):
    """
    Add a multi-vintage electrolyser converter:

      Tech:       Electrolyser
      Activity:   Electrolysis
      Commodities: Elec (input), H2 (output)

    Vintages: 2020, 2030, 2040, 2050
    - lifetime, efficiency and CAPEX are defined per vintage year.
    - REMix will map build years to these vintages according to its internal rules.
    """

    global yrs_to_calc

    # -----------------------------
    # 1) Define vintages and nodes
    # -----------------------------
    eltr_vintage = [2020, 2030, 2040, 2050]
    eltr_nodes = [n for n in m["Base"].set.nodesdata]

    # -----------------------------
    # 2) Tech parameters (per vintage)
    # -----------------------------
    eltr_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Electrolyser"], eltr_vintage],
            names=["techs", "years"],
        )
    )
    eltr_tech["lifeTime"] = [25, 30, 32, 35]  # years per vintage
    eltr_tech["activityUpperLimit"] = 1       # availability factor
    m["Base"].parameter.add(eltr_tech, "converter_techparam")

    # -----------------------------
    # 3) Capacity bounds (all model/data years)
    # -----------------------------
    eltr_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [eltr_nodes, yrs_to_calc, ["Electrolyser"]],
            names=["nodesdata", "years", "techs"],
        )
    )
    eltr_caps["unitsUpperLimit"] = 100.0  # GW_el (large generic upper bound)
    eltr_caps.loc[idx[:, 2020, :], "noExpansion"] = 1 
    m["Base"].parameter.add(eltr_caps, "converter_capacityparam")

    # -----------------------------
    # 4) Conversion coefficients
    # -----------------------------
    eltr_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["Electrolyser"],        # techs
                eltr_vintage,            # years (vintages)
                ["Electrolysis"],        # activities
                ["Elec", "H2"],          # commodities
            ],
            names=["techs", "years", "activities", "commodities"],
        )
    )

    # 1 GWh activity consumes 1 GWh electricity
    eltr_coef.loc[idx["Electrolyser", :, "Electrolysis", "Elec"], "coefficient"] = -1.0

    # 1 GWh activity produces H2 with vintage-specific efficiency
    eltr_coef.loc[idx["Electrolyser", :, "Electrolysis", "H2"], "coefficient"] = [
        0.665, 0.79, 0.79, 0.85
    ]
    m["Base"].parameter.add(eltr_coef, "converter_coefficient")

    # -----------------------------
    # 5) Accounting (converter units, v13 structure)
    # -----------------------------
    electrolyser_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], ["Electrolyser"], eltr_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Investment cost and annuity settings per vintage
    electrolyser_acc.loc[
        idx["Invest", "global", "horizon", "Electrolyser", :], "perUnitBuild"
    ] = [750, 570, 450, 350]  # M€ per unit
    electrolyser_acc.loc[
        idx["Invest", "global", "horizon", "Electrolyser", :], "amorTime"
    ] = [25, 30, 32, 35]      # years
    electrolyser_acc.loc[
        idx["Invest", "global", "horizon", "Electrolyser", :], "useAnnuity"
    ] = 1
    electrolyser_acc.loc[
        idx["Invest", "global", "horizon", "Electrolyser", :], "interest"
    ] = 0.06

    # Fixed O&M = 1.4 % of CAPEX per year
    invest_vals = electrolyser_acc.loc[
        idx["Invest", "global", "horizon", "Electrolyser", :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", "Electrolyser", i[-1]) for i in invest_vals.index],
        names=electrolyser_acc.index.names,
    )
    electrolyser_acc.loc[
        idx["OMFix", "global", "horizon", "Electrolyser", :], "perUnitTotal"
    ] = invest_vals * 0.014

    electrolyser_acc = electrolyser_acc.fillna(0.0)
    m["Base"].parameter.add(electrolyser_acc, "accounting_converterunits")

def add_H2_CCGT(m):
    import pandas as pd

    global yrs_to_calc

    H2_CCGT_vintage = [2030, 2035, 2040, 2045, 2050]
    H2_CCGT_nodes   = [n for n in m["Base"].set.nodesdata]

    # -----------------------------
    # 1) Tech parameters (per vintage)
    # -----------------------------
    H2_CCGT_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["H2_CCGT"], H2_CCGT_vintage],
            names=["techs", "years"],
        )
    )
    H2_CCGT_tech["lifeTime"]         = 35  # years
    H2_CCGT_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(H2_CCGT_tech, "converter_techparam")

    # -----------------------------
    # 2) Capacity bounds (all data years)
    # -----------------------------
    H2_CCGT_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [H2_CCGT_nodes, yrs_to_calc, ["H2_CCGT"]],
            names=["nodesdata", "years", "techs"],
        )
    )
    H2_CCGT_caps["unitsUpperLimit"] = 100.0
    m["Base"].parameter.add(H2_CCGT_caps, "converter_capacityparam")

    # -----------------------------
    # 3) Conversion coefficients
    # -----------------------------
    H2_CCGT_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["H2_CCGT"],           # techs
                H2_CCGT_vintage,       # years (vintages)
                ["Powergen"],          # activities
                ["Elec", "H2"],        # commodities
            ],
            names=["techs", "years", "activities", "commodities"],
        )
    )
    # 1 GWh activity consumes 1 GWh H2
    H2_CCGT_coef.loc[idx["H2_CCGT", :, "Powergen", "H2"], "coefficient"] = -1.0
    # Electrical efficiency by vintage (C. Habib)
    H2_CCGT_coef.loc[idx["H2_CCGT", :, "Powergen", "Elec"], "coefficient"] = [
        0.58, 0.59, 0.60, 0.60, 0.60
    ]
    m["Base"].parameter.add(H2_CCGT_coef, "converter_coefficient")

    # -----------------------------
    # 4) Accounting (converter units)
    # -----------------------------
    H2_CCGT_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], ["H2_CCGT"], H2_CCGT_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    H2_CCGT_acc.loc[idx["Invest", "global", "horizon", "H2_CCGT", :], "perUnitBuild"] = 100.0
    H2_CCGT_acc.loc[idx["Invest", "global", "horizon", "H2_CCGT", :], "amorTime"]     = 35
    H2_CCGT_acc.loc[idx["Invest", "global", "horizon", "H2_CCGT", :], "useAnnuity"]  = 1
    H2_CCGT_acc.loc[idx["Invest", "global", "horizon", "H2_CCGT", :], "interest"]    = 0.06

    invest_vals = H2_CCGT_acc.loc[
        idx["Invest", "global", "horizon", "H2_CCGT", :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", "H2_CCGT", i[-1]) for i in invest_vals.index],
        names=H2_CCGT_acc.index.names,
    )
    H2_CCGT_acc.loc[
        idx["OMFix", "global", "horizon", "H2_CCGT", :], "perUnitTotal"
    ] = invest_vals * 0.025  # 2.5 % fixed O&M

    H2_CCGT_acc = H2_CCGT_acc.fillna(0.0)
    m["Base"].parameter.add(H2_CCGT_acc, "accounting_converterunits")

def add_H2_FC(m):
    import pandas as pd

    global yrs_to_calc

    H2_FC_vintage = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    H2_FC_nodes   = [n for n in m["Base"].set.nodesdata]

    # -----------------------------
    # 1) Tech parameters
    # -----------------------------
    H2_FC_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["H2_FC"], H2_FC_vintage],
            names=["techs", "years"],
        )
    )
    H2_FC_tech["lifeTime"]          = 35
    H2_FC_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(H2_FC_tech, "converter_techparam")

    # -----------------------------
    # 2) Capacity bounds
    # -----------------------------
    H2_FC_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [H2_FC_nodes, yrs_to_calc, ["H2_FC"]],
            names=["nodesdata", "years", "techs"],
        )
    )
    H2_FC_caps["unitsUpperLimit"] = 100.0
    m["Base"].parameter.add(H2_FC_caps, "converter_capacityparam")

    # -----------------------------
    # 3) Conversion coefficients
    # -----------------------------
    H2_FC_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["H2_FC"],
                H2_FC_vintage,
                ["Powergen"],
                ["Elec", "H2"],
            ],
            names=["techs", "years", "activities", "commodities"],
        )
    )
    H2_FC_coef.loc[idx["H2_FC", :, "Powergen", "H2"], "coefficient"] = -1.0
    H2_FC_coef.loc[idx["H2_FC", :, "Powergen", "Elec"], "coefficient"] = [
        0.579, 0.6134, 0.6383, 0.6477, 0.6514, 0.6686, 0.6737
    ]
    m["Base"].parameter.add(H2_FC_coef, "converter_coefficient")

    # -----------------------------
    # 4) Accounting
    # -----------------------------
    H2_FC_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], ["H2_FC"], H2_FC_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    H2_FC_acc.loc[idx["Invest", "global", "horizon", "H2_FC", :], "perUnitBuild"] = 100.0
    H2_FC_acc.loc[idx["Invest", "global", "horizon", "H2_FC", :], "amorTime"]     = 35
    H2_FC_acc.loc[idx["Invest", "global", "horizon", "H2_FC", :], "useAnnuity"]  = 1
    H2_FC_acc.loc[idx["Invest", "global", "horizon", "H2_FC", :], "interest"]    = 0.06

    invest_vals = H2_FC_acc.loc[
        idx["Invest", "global", "horizon", "H2_FC", :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", "H2_FC", i[-1]) for i in invest_vals.index],
        names=H2_FC_acc.index.names,
    )
    H2_FC_acc.loc[
        idx["OMFix", "global", "horizon", "H2_FC", :], "perUnitTotal"
    ] = invest_vals * 0.05  # 5 % fixed O&M

    H2_FC_acc = H2_FC_acc.fillna(0.0)
    m["Base"].parameter.add(H2_FC_acc, "accounting_converterunits")

def add_h2_storage(m):
    import pandas as pd

    global yrs_to_calc

    # Compressor converter + storage tank
    h2_stor_techs = ["H2_storage"]
    h2_stor_nodes = [n for n in m["Base"].set.nodesdata]

    # Technology years / vintages for this block
    yrs_h2 = [2020, 2035, 2050]

    # -----------------------------
    # 1) Converter: tech parameters
    # -----------------------------
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs, yrs_h2],
            names=["techs", "years"],
        )
    )
    conv_tech["lifeTime"]          = 40
    conv_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # -----------------------------
    # 2) Converter: capacity bounds
    # -----------------------------
    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_nodes, yrs_h2, h2_stor_techs],
            names=["nodesdata", "years", "techs"],
        )
    )
    conv_cap["unitsUpperLimit"] = 50.0
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # -----------------------------
    # 3) Converter: coefficients (charge/discharge)
    # -----------------------------
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs, yrs_h2, ["Charge", "Discharge"], ["H2", "H2_stored", "Elec"]],
            names=["techs", "years", "activities", "commodities"],
        )
    )

    # Charge: consume Elec + H2, store H2_stored (incl. efficiency)
    conv_coef.loc[idx["H2_storage", :, "Charge", "Elec"], "coefficient"] = [
        -0.043346085,
        -0.035470537,
        -0.031529343,
    ]
    conv_coef.loc[idx["H2_storage", :, "Charge", "H2"], "coefficient"]        = -0.15
    conv_coef.loc[idx["H2_storage", :, "Charge", "H2_stored"], "coefficient"] = 0.15 * 0.89

    # Discharge: recover H2 from H2_stored with losses
    conv_coef.loc[idx["H2_storage", :, "Discharge", "H2"], "coefficient"]        = 0.15
    conv_coef.loc[idx["H2_storage", :, "Discharge", "H2_stored"], "coefficient"] = -0.15 * 1.11

    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # -----------------------------
    # 4) Converter accounting
    # -----------------------------
    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], h2_stor_techs, yrs_h2],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    )

    conv_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "perUnitBuild"] = 5.14
    conv_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "useAnnuity"]   = 1
    conv_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "amorTime"]     = 40
    conv_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "interest"]     = 0.06

    invest_vals = conv_acc.loc[
        idx["Invest", "global", "horizon", "H2_storage", :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", "H2_storage", i[-1]) for i in invest_vals.index],
        names=conv_acc.index.names,
    )
    conv_acc.loc[
        idx["OMFix", "global", "horizon", "H2_storage", :], "perUnitTotal"
    ] = invest_vals * 0.04  # 4 % fixed O&M

    conv_acc = conv_acc.fillna(0.0)
    m["Base"].parameter.add(conv_acc, "accounting_converterunits")

    # -----------------------------
    # 5) Storage: tech parameters
    # -----------------------------
    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs, yrs_h2],
            names=["techs", "years"],
        )
    )
    stor_tech["lifeTime"]        = 40
    stor_tech["levelUpperLimit"] = 1
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    # -----------------------------
    # 6) Storage: size per unit
    # -----------------------------
    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs, yrs_h2, ["H2_stored"]],
            names=["techs", "years", "commodities"],
        )
    )
    stor_size.loc[idx["H2_storage", :, "H2_stored"], "size"] = 6.9993  # GWh / unit
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # -----------------------------
    # 7) Storage: reservoir parameters
    # -----------------------------
    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_nodes, yrs_h2, h2_stor_techs],
            names=["nodesdata", "years", "techs"],
        )
    )
    stor_res["unitsUpperLimit"] = 50
    stor_res.loc[idx[:, 2020, :], "noExpansion"] = 1  # no expansion of existing 2020 capacity
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    # -----------------------------
    # 8) Storage accounting
    # -----------------------------
    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], h2_stor_techs, yrs_h2],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    )

    stor_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "perUnitBuild"] = 131.91
    stor_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "useAnnuity"]   = 1
    stor_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "amorTime"]     = 40
    stor_acc.loc[idx["Invest", "global", "horizon", "H2_storage", :], "interest"]     = 0.06

    invest_vals = stor_acc.loc[
        idx["Invest", "global", "horizon", "H2_storage", :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", "H2_storage", i[-1]) for i in invest_vals.index],
        names=stor_acc.index.names,
    )
    stor_acc.loc[
        idx["OMFix", "global", "horizon", "H2_storage", :], "perUnitTotal"
    ] = invest_vals * 0.03  # 3 % fixed O&M

    stor_acc = stor_acc.fillna(0.0)
    m["Base"].parameter.add(stor_acc, "accounting_storageunits")

# synfuels

def add_dac(m):
    """
    Add Direct Air Capture (DAC) technology:

    - Tech name: DAC
    - Vintages: 2020, 2030, 2040, 2050
    - Activity: "Capture"
    - Commodities:
        * Elec (input, negative coefficient)
        * CO2_feed (output, positive coefficient)
    - Accounting:
        * Investment and OMFix per installed unit (ktCO2/h or similar)
        * CO2_emission: negative (removal from atmosphere)
    """

    # ------------------------------------------------------------------
    # 1. Basic sets: techs, vintages, nodes, years
    # ------------------------------------------------------------------
    dac_techs = ["DAC"]
    dac_vintages = [2020, 2030, 2040, 2050]

    # Use all electricity nodes except LNG (same style as you showed)
    dac_nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]

    # Optimisation years (ensure int)
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)

    print("\n--- ADDING DAC ---")
    print("DAC techs:", dac_techs)
    print("DAC vintages:", dac_vintages)
    print("DAC nodes:", dac_nodes)
    print("Optimisation years:", years_sel)

    # ------------------------------------------------------------------
    # 2. converter_techparam: (techs, vintages)
    # ------------------------------------------------------------------
    dac_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [dac_techs, dac_vintages],
            names=["techs", "vintages"],
        )
    )

    # Technical lifetime (years)
    dac_tech.loc[idx[:, :], "lifeTime"] = 20

    # Availability flag (1 = can operate up to capacity; profile can further limit)
    dac_tech["activityUpperLimit"] = 1.0

    m["Base"].parameter.add(dac_tech, "converter_techparam")
    print("DAC converter_techparam rows:", len(dac_tech))

    # ------------------------------------------------------------------
    # 3. converter_capacityparam: (nodesdata, years, techs)
    # ------------------------------------------------------------------
    # Here we allow DAC to be potentially installed up to a large upper
    # limit in all optimisation years; you can refine by year later.
    cap_index = pd.MultiIndex.from_product(
        [dac_nodes, years_sel, dac_techs],
        names=["nodesdata", "years", "techs"],
    )
    dac_caps = pd.DataFrame(index=cap_index)

    # Capacity upper limit (units) – interpret as ktCO2/h capacity or similar
    dac_caps["unitsUpperLimit"] = 100.0
    dac_caps.loc[idx[:, [2020], :], "noExpansion"] = 1

    # No special brownfield builds for now; noExpansion not set.
    m["Base"].parameter.add(dac_caps, "converter_capacityparam")
    print("DAC converter_capacityparam rows:", len(dac_caps))

    # ------------------------------------------------------------------
    # 4. converter_coefficient: (techs, vintages, activities, commodities)
    # ------------------------------------------------------------------
    # Activity name "Capture":
    #   - Output: CO2_feed (coefficient > 0)
    #   - Input: Elec (coefficient < 0)
    #
    # Coefficients are in "commodity units per unit of activity".
    #FIXME: verify 1 unit activity produces 1 ktCO2_feed and consumes~1.5 GWh_el
    dac_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [dac_techs, dac_vintages, ["Capture"], ["Elec", "CO2_feed"]],
            names=["techs", "vintages", "activities", "commodities"],
        )
    )

    # CO2_feed output: 1 unit of activity -> 1 unit CO2_feed
    dac_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = 1.0

    # Electricity input: negative coefficients (kWh/GWh per unit CO2 capture)
    elec_per_kt = [-1.535, -1.458, -1.385, -1.316]  # GWh_el per ktCO2 (negative sign)
    dac_coef.loc[idx[:, 2020, :, "Elec"], "coefficient"] = elec_per_kt[0]
    dac_coef.loc[idx[:, 2030, :, "Elec"], "coefficient"] = elec_per_kt[1]
    dac_coef.loc[idx[:, 2040, :, "Elec"], "coefficient"] = elec_per_kt[2]
    dac_coef.loc[idx[:, 2050, :, "Elec"], "coefficient"] = elec_per_kt[3]

    m["Base"].parameter.add(dac_coef, "converter_coefficient")
    print("DAC converter_coefficient rows:", len(dac_coef))

    # ------------------------------------------------------------------
    # 5. accounting_converterunits: CAPEX + OMFix
    # ------------------------------------------------------------------
    # Index: (indicator, regionscope, timescope, techs, years)
    dac_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], dac_techs, dac_vintages],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Investment cost: 815 EUR/tCO2*a -> convert to M€/unit-year basis.
    # FIXME: verify
    dac_acc.loc[idx["Invest", :, :, :, :], "perUnitBuild"] = 815.0

    # Annuity parameters
    dac_acc.loc[idx["Invest", :, :, :, :], "amorTime"] = 20  # years
    dac_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
    dac_acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06

    # OMFix as fixed fraction (4%) of investment per year
    invest_vals = dac_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
        names=dac_acc.index.names,
    )
    dac_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.04

    m["Base"].parameter.add(dac_acc, "accounting_converterunits")
    print("DAC accounting_converterunits rows:", len(dac_acc))

def add_methanizer(m):
    """
    Add methanation unit:

    - Tech: Methanizer
    - Consumes: H2, CO2_feed, Elec
    - Produces: CH4
    - No direct CO2_emission accounting here; assume CO2 is internal stream.
    """

    meth_techs = ["Methanizer"]
    meth_vintages = [2020, 2030, 2040, 2050]
    meth_nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)

    print("\n--- ADDING METHANIZER ---")
    print("nodes:", meth_nodes)
    print("vintages:", meth_vintages)

    # 1. converter_techparam (techs, vintages)
    meth_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [meth_techs, meth_vintages],
            names=["techs", "vintages"],
        )
    )
    meth_tech["lifeTime"] = 30
    meth_tech["activityUpperLimit"] = 1.0
    m["Base"].parameter.add(meth_tech, "converter_techparam")

    # 2. converter_capacityparam (nodesdata, years, techs)
    cap_index = pd.MultiIndex.from_product(
        [meth_nodes, years_sel, meth_techs],
        names=["nodesdata", "years", "techs"],
    )
    meth_caps = pd.DataFrame(index=cap_index)
    meth_caps["unitsUpperLimit"] = 100.0
    meth_caps.loc[idx[:, [2020], :], "noExpansion"] = 1
    m["Base"].parameter.add(meth_caps, "converter_capacityparam")

    # 3. converter_coefficient (techs, vintages, activities, commodities)
    meth_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [meth_techs, meth_vintages, ["Synthesis"], ["Elec", "H2", "e-CH4", "CO2_feed"]],
            names=["techs", "vintages", "activities", "commodities"],
        )
    )
    # Output CH4 = 1
    meth_coef.loc[idx[:, :, :, "e-CH4"], "coefficient"] = 1.0
    # Inputs: negative coefficients (per unit CH4 output)
    meth_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.284
    meth_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.2016
    meth_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.006

    m["Base"].parameter.add(meth_coef, "converter_coefficient")

    # 4. accounting_converterunits (Invest, OMFix)
    meth_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], meth_techs, meth_vintages],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Investment per unit (you gave 4 different values)
    invest_vals_vintage = [558, 309, 251, 211]  # M€/unit
    for v, val in zip(meth_vintages, invest_vals_vintage):
        meth_acc.loc[idx["Invest", "global", "horizon", "Methanizer", v], "perUnitBuild"] = val

    meth_acc.loc[idx["Invest", :, :, :, :], "amorTime"] = [20, 25, 30, 30]
    meth_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
    meth_acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06

    # OMFix = 4% of Invest
    invest_vals = meth_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
        names=meth_acc.index.names,
    )
    meth_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.04

    m["Base"].parameter.add(meth_acc, "accounting_converterunits")

def add_methanol_syn(m):
    """
    Add methanol synthesis unit:

    - Tech: MethanolSyn
    - Consumes: H2, CO2_feed, Elec
    - Produces: CH3OH
    """


    techs = ["MethanolSyn"]
    vintages = [2020]
    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)

    print("\n--- ADDING METHANOL SYNTHESIS ---")

    # 1. converter_techparam
    tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages],
            names=["techs", "vintages"],
        )
    )
    tech["lifeTime"] = 30
    tech["activityUpperLimit"] = 1.0
    m["Base"].parameter.add(tech, "converter_techparam")

    # 2. converter_capacityparam
    cap_index = pd.MultiIndex.from_product(
        [nodes, years_sel, techs],
        names=["nodesdata", "years", "techs"],
    )
    caps = pd.DataFrame(index=cap_index)
    caps["unitsUpperLimit"] = 100.0
    caps.loc[idx[:, [2020], :], "noExpansion"] = 1
    m["Base"].parameter.add(caps, "converter_capacityparam")

    # 3. converter_coefficient
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages, ["Synthesis"], ["CH3OH", "H2", "CO2_feed", "Elec"]],
            names=["techs", "vintages", "activities", "commodities"],
        )
    )
    coef.loc[idx[:, :, :, "CH3OH"], "coefficient"] = 1.0
    coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.25
    coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.219
    coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.1

    m["Base"].parameter.add(coef, "converter_coefficient")

    # 4. accounting_converterunits
    acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, vintages],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    acc.loc[idx["Invest", :, :, :, :], "perUnitBuild"] = 835.0
    acc.loc[idx["Invest", :, :, :, :], "amorTime"] = 30
    acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
    acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06

    invest_vals = acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
        names=acc.index.names,
    )
    acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.04

    m["Base"].parameter.add(acc, "accounting_converterunits")

def add_ftropsch_syn(m):
    """
    Add Fischer-Tropsch synthesis unit:

    - Tech: FTropschSyn
    - Produces: REfuel (single synthetic liquid fuel)
    - Consumes: H2, CO2_feed, Elec
    - No direct CO2_emission accounting here (CO2 handled via CO2_feed + DAC).
    """

    techs = ["FTropschSyn"]
    vintages = [2020, 2040]  # technology vintages
    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)

    print("\n--- ADDING FTROPSCH SYNTHESIS ---")
    print("nodes:", nodes)
    print("vintages:", vintages)
    print("years_sel:", years_sel)

    # ------------------------------------------------------------------
    # 1. converter_techparam: (techs, vintages)
    # ------------------------------------------------------------------
    ft_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages],
            names=["techs", "vintages"],
        )
    )
    # Technical lifetime for each vintage (years)
    ft_tech["lifeTime"] = 30
    # Availability flag; profiles/capacities will further constrain operation
    ft_tech["activityUpperLimit"] = 1.0

    m["Base"].parameter.add(ft_tech, "converter_techparam")
    print("FTropsch converter_techparam rows:", len(ft_tech))

    # ------------------------------------------------------------------
    # 2. converter_capacityparam: (nodesdata, years, techs)
    # ------------------------------------------------------------------
    cap_index = pd.MultiIndex.from_product(
        [nodes, years_sel, techs],
        names=["nodesdata", "years", "techs"],
    )
    ft_caps = pd.DataFrame(index=cap_index)

    # Large generic upper limit (MWh/GW-equivalent units); refine later as needed
    ft_caps["unitsUpperLimit"] = 100.0
    ft_caps.loc[idx[:, [2020], :], "noExpansion"] = 1

    m["Base"].parameter.add(ft_caps, "converter_capacityparam")
    print("FTropsch converter_capacityparam rows:", len(ft_caps))

    # ------------------------------------------------------------------
    # 3. converter_coefficient: (techs, vintages, activities, commodities)
    # ------------------------------------------------------------------
    # Activity "Synthesis":
    #   Output: REfuel (coefficient > 0)
    #   Inputs: H2, CO2_feed, Elec (coefficients < 0)
    #
    # Coefficients are in "commodity units per unit of activity".
    ft_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages, ["Synthesis"], ["REfuel", "H2", "CO2_feed", "Elec"]],
            names=["techs", "vintages", "activities", "commodities"],
        )
    )

    # Output: 1 unit activity -> 1 unit REfuel
    ft_coef.loc[idx[:, :, :, "REfuel"], "coefficient"] = 1.0

    # Inputs: numbers derived from e-fuels split (65% synthesis efficiency etc.)
    # Here we aggregate to a single REfuel
    ft_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1.52    # H2 consumption per unit REfuel
    ft_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.35  # CO2_feed consumption
    ft_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -0.11      # Electricity consumption

    m["Base"].parameter.add(ft_coef, "converter_coefficient")
    print("FTropsch converter_coefficient rows:", len(ft_coef))

    # ------------------------------------------------------------------
    # 4. accounting_converterunits: CAPEX + OMFix
    # ------------------------------------------------------------------
    ft_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, vintages],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # Investment costs per vintage (M€/unit)
    invest_vals_vintage = [1017, 915]
    for v, val in zip(vintages, invest_vals_vintage):
        ft_acc.loc[idx["Invest", "global", "horizon", "FTropschSyn", v], "perUnitBuild"] = val

    # Annuity parameters (same for both vintages)
    ft_acc.loc[idx["Invest", :, :, :, :], "amorTime"] = 30
    ft_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
    ft_acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06

    # OMFix = 4% of Invest per year
    invest_vals = ft_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", "global", "horizon", i[3], i[4]) for i in invest_vals.index],
        names=ft_acc.index.names,
    )
    ft_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.04

    m["Base"].parameter.add(ft_acc, "accounting_converterunits")
    print("FTropsch accounting_converterunits rows:", len(ft_acc))


# demand

def add_demand(m):
    """
    Electricity demand + electricity slack ONLY.

    - Electricity demand: sourcesink_profile type=fixed (negative)
      sourcesink_config usesFixedProfile=1

    - Electricity slack: allow feasibility with VERY HIGH penalty
      sourcesink_config usesUpperSum=1 and usesLowerProfile=1
      sourcesink_annualsum upper = inf
      accounting_sourcesinkflow SlackCost perFlow = big

    Notes:
    - self-contained
    - robust to year types (casts to int)
    - robust to string hour columns (coerces t-cols to numeric)
    """

    print("\n--- ADDING ELECTRICITY DEMAND ---")

    # years to run optimisation (ensure int)
    years_sel = sorted([int(y) for y in list(m["Base"].set.yearssel)])
    print("yearssel (optimise):", years_sel)

    if str(group_name).strip() != "GP-NT-ELEC-BIO-H2":
        raise ValueError("This add_demand() version expects GP-NT-ELEC-BIO-H2 hourly electricity input format.")

    elec_path = Path(path_demand) / hourly_electricity_file
    if not elec_path.exists():
        raise FileNotFoundError(f"Hourly electricity file not found: {elec_path}")

    df = pd.read_csv(elec_path)

    # Normalize column names
    rename_cols = {}
    if "scenario" in df.columns and "Scenario" not in df.columns:
        rename_cols["scenario"] = "Scenario"
    if "Region" in df.columns:
        rename_cols["Region"] = "node"
    if "region" in df.columns:
        rename_cols["region"] = "node"
    if "node" not in df.columns and "Node" in df.columns:
        rename_cols["Node"] = "node"
    if "Year" in df.columns:
        rename_cols["Year"] = "year"
    if "Sector" in df.columns:
        rename_cols["Sector"] = "sector"
    if "Carrier" in df.columns:
        rename_cols["Carrier"] = "carrier"
    if rename_cols:
        df = df.rename(columns=rename_cols)

    required = {"Scenario", "node", "year", "sector", "carrier"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Hourly electricity CSV missing columns after normalisation: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    # Map regions to REMix codes
    df["node"] = df["node"].apply(lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r))
    df = df.dropna(subset=["node"])

    # Scenario filter
    df = df.loc[df["Scenario"].astype(str).str.strip() == str(base_scenario)].copy()
    df.drop(columns=["Scenario"], inplace=True)

    # Keep electricity only
    df["carrier"] = df["carrier"].astype(str).str.strip()
    df = df.loc[df["carrier"].str.lower() == "electricity"].copy()

    # Filter years (int)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"])
    df["year"] = df["year"].astype(int)
    df = df.loc[df["year"].isin(years_sel)].copy()

    # Hour columns
    hour_cols = [c for c in df.columns if str(c).startswith("t")]
    if not hour_cols:
        raise ValueError("No hourly columns found (expected t0001..t8760).")

    # Coerce hours to numeric (prevents float+str bugs)
    df[hour_cols] = df[hour_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    print("demand rows after filters:", len(df))
    print("hour_cols:", len(hour_cols), "sample:", hour_cols[:5], "...", hour_cols[-5:])

    # Aggregate to Sector=All
    df_agg = df.groupby(["node", "year", "carrier"], as_index=False)[hour_cols].sum()
    df_agg["sector"] = "All"
    df_agg["carrier"] = "Elec"

    # Build profile (negative = demand)
    ts = df_agg.set_index(["node", "year", "sector", "carrier"])[hour_cols]
    ts *= -1
    ts["type"] = "fixed"
    ts_fixed = ts.set_index("type", append=True).round(6)

    m["Base"].profile.add(ts_fixed, "sourcesink_profile")

    # Config: fixed profile
    cfg_idx = ts_fixed.index.droplevel(-1)  # drop type
    cfg = pd.DataFrame(index=cfg_idx)
    cfg["usesFixedProfile"] = 1
    m["Base"].parameter.add(cfg, "sourcesink_config")

    # Ensure nodes exist
    m["Base"].set.add(sorted(cfg.index.get_level_values(0).unique().tolist()), "nodesdata")

    print(f"Electricity demand loaded: {len(cfg)} rows (Sector='All').")

    # =========================
    # Electricity slack (strongly penalised)
    # =========================
    print("\n--- ADDING ELECTRICITY SLACK ---")

    # annualsum: allow slack production (upper=inf)
    slack_annual = cfg.copy()
    slack_annual.index = pd.MultiIndex.from_tuples(
        [(i[0], i[1], "Slack", i[3]) for i in slack_annual.index],
        names=cfg.index.names,
    )
    slack_annual = slack_annual.rename(columns={"usesFixedProfile": "upper"})
    slack_annual["upper"] = np.inf
    m["Base"].parameter.add(slack_annual[["upper"]], "sourcesink_annualsum")

    # slack config: usesUpperSum=1 + usesLowerProfile=1
    slack_cfg = pd.DataFrame(index=slack_annual.index)
    slack_cfg["usesUpperSum"] = 1
    slack_cfg["usesLowerProfile"] = 1
    m["Base"].parameter.add(slack_cfg, "sourcesink_config")

    # strong penalty
    # NOTE: keep it huge. If objective is in M€, 1e6 is “very expensive”.
    slack_cost = 1e6

    # accounting_sourcesinkflow index style: (indicator, regionscope, years, sector, commodities)
    years_data = sorted([int(y) for y in list(m["Base"].set.years)])
    slack_cost_df = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["SlackCost"], ["global"], years_data, ["Slack"], ["Elec"]],
            names=["indicator", "regionscope", "years", "sector", "commodities"],
        )
    )
    slack_cost_df["perFlow"] = float(slack_cost)
    m["Base"].parameter.add(slack_cost_df, "accounting_sourcesinkflow")

    print(f"Electricity slack added: {len(slack_cfg)} rows.")
    print("SlackCost perFlow:", slack_cost, "years_data:", len(years_data), "range:", (min(years_data), max(years_data)))

def add_fuel_demands_from_excel(m):
    """
    Fuel demand + fuel slack in ONE function, using the SAME index style as add_demand().

      - The sourcesink name lives in the 'sector' index level (NOT 'sourcesinks').
      - So we write sourcesink_annualsum and sourcesink_config with index:
            (nodesdata, years, sector, commodities)

    Demand:
      - Fixed annual sinks: lower = upper = -Demand  (cannot flip sign / cancel)

    Slack:
      - Positive-only supply slack: lower = 0, upper = +inf
      - Strong penalty SlackCost perFlow (positive)
    """

    print("\n--- ADDING FUEL DEMAND + FUEL SLACK (Heat/Transport) ---")

    global path_demand, carriers_excel_file, base_scenario, yrs_sel, map_region_to_remix

    excel_path = Path(path_demand) / carriers_excel_file
    if not excel_path.exists():
        raise FileNotFoundError(f"Carriers Excel not found: {excel_path}")

    years_sel = sorted(int(y) for y in yrs_sel)

    # -------------------------
    # STEP 1: Read + scenario filter
    # -------------------------
    df = pd.read_excel(excel_path, sheet_name="carriers")
    df.columns = [str(c).strip() for c in df.columns]

    df = df.loc[df["Scenario"].astype(str).str.strip() == str(base_scenario)].copy()
    if df.empty:
        print("[fuel demand] No rows for scenario:", base_scenario)
        return {"nodes": [], "years": years_sel, "fuels": []}

    # Map region names to REMix nodes
    df["node"] = df["Region"].apply(map_region_to_remix)
    missing = df.loc[df["node"].isna(), "Region"].unique().tolist()
    if missing:
        raise ValueError("Unmapped regions in carriers Excel: " + ", ".join(map(str, missing)))

    # -------------------------
    # STEP 2: Melt year columns
    # -------------------------
    year_cols = [c for c in df.columns if str(c).isdigit()]
    if not year_cols:
        raise ValueError("No year columns found in carriers sheet (expected 2020, 2025, ...).")

    df_long = df.melt(
        id_vars=["node", "Sector", "Carrier"],
        value_vars=year_cols,
        var_name="year",
        value_name="Demand",
    )
    df_long["year"] = pd.to_numeric(df_long["year"], errors="coerce")
    df_long = df_long.dropna(subset=["year"])
    df_long["year"] = df_long["year"].astype(int)

    df_long["Demand"] = pd.to_numeric(df_long["Demand"], errors="coerce").fillna(0.0)
    df_long = df_long.loc[df_long["year"].isin(years_sel)].copy()
    df_long = df_long.loc[df_long["Demand"] != 0.0].copy()
    if df_long.empty:
        print("[fuel demand] No non-zero demand in selected years.")
        return {"nodes": [], "years": years_sel, "fuels": []}

    # -------------------------
    # STEP 3: Carrier -> commodity mapping (robust)
    # -------------------------
    carrier_to_fuel = {
        "e-fuel (gas)": "e-CH4",
        "e-fuel (lf)": "REfuel",
        "hydrogen": "H2",

        # optional:
        "biofuel (gas)": "Bio_gas",
        "biofuel (lf)": "Bio_LF",
        "fossil (gas)": "Fossil_CH4",
        "fossil (lf)": "Fossil_LF",
    }

    df_long["carrier_key"] = df_long["Carrier"].astype(str).str.strip().str.lower()
    df_long["commodity"] = df_long["carrier_key"].map(carrier_to_fuel)

    unmapped = sorted(df_long.loc[df_long["commodity"].isna(), "Carrier"].unique().tolist())
    if unmapped:
        print("[fuel demand] WARNING unmapped carriers ignored:", unmapped)

    df_long = df_long.loc[df_long["commodity"].notna()].copy()
    if df_long.empty:
        print("[fuel demand] After mapping, no carriers remain.")
        return {"nodes": [], "years": years_sel, "fuels": []}

    # -------------------------
    # STEP 4: Heat/Transport only, group
    # -------------------------
    df_long["sector"] = df_long["Sector"].astype(str).str.strip()
    df_long = df_long.loc[df_long["sector"].isin(["Heat", "Transport"])].copy()
    if df_long.empty:
        print("[fuel demand] No Heat/Transport rows.")
        return {"nodes": [], "years": years_sel, "fuels": []}

    grp = (
        df_long.groupby(["node", "year", "sector", "commodity"], as_index=False)["Demand"]
        .sum()
    )

    nodes = sorted(grp["node"].unique().tolist())
    years = sorted(grp["year"].unique().tolist())
    fuels = sorted(grp["commodity"].unique().tolist())

    print("[fuel demand] nodes:", nodes)
    print("[fuel demand] years:", years)
    print("[fuel demand] fuels:", fuels)

    # Ensure nodes + commodities exist so they show up in balances
    m["Base"].set.add(nodes, "nodesdata")
    m["Base"].set.add(fuels, "commodities")

    # -------------------------
    # STEP 5: Demand sinks (sector names act as sourcesink names)
    # index: (node, year, sector, commodity)  <-- SAME STYLE AS add_demand()
    # -------------------------
    demand_rows = []
    for r in grp.itertuples(index=False):
        if r.sector == "Heat":
            ss_sector = f"HeatDemand_{r.commodity}"
        else:
            ss_sector = f"TranspDemand_{r.commodity}"

        fixed = -float(r.Demand)
        demand_rows.append((r.node, int(r.year), ss_sector, r.commodity, fixed, fixed))

    dem_idx = pd.MultiIndex.from_tuples(
        [(n, y, sec, c) for (n, y, sec, c, lo, up) in demand_rows],
        names=["nodesdata", "years", "sector", "commodities"],
    )
    dem = pd.DataFrame(index=dem_idx)
    dem["lower"] = [lo for (*_, lo, __) in [(n,y,sec,c,lo,up) for (n,y,sec,c,lo,up) in demand_rows]]
    dem["upper"] = [up for (*_, up) in [(n,y,sec,c,lo,up) for (n,y,sec,c,lo,up) in demand_rows]]
    m["Base"].parameter.add(dem, "sourcesink_annualsum")

    dem_cfg = pd.DataFrame(index=dem.index)
    dem_cfg["usesLowerSum"] = 1
    dem_cfg["usesUpperSum"] = 1
    m["Base"].parameter.add(dem_cfg, "sourcesink_config")

    print("[fuel demand] demand sinks written:", len(dem))

    # ---- CO2 emissions for fuel use in Heat/Transport ----
    # Units must match REMix expectations for accounting_sourcesinkflow "perFlow"
    # (commonly ktCO2 per GWh of fuel commodity flow).


    ef_kt_per_gwh = {
        # Fossil fuels (from table Gulagi et al. 2025)
        "Fossil_LF":   0.2624,   # oil / diesel / kerosene # ~262.4 kg/MWh = 0.2624 kt/GWh
        "Fossil_CH4":  0.2048,   # natural gas
        "Fossil_Coal": 0.3546,

        # E-fuels emit at combustion (count them like the equivalent fuel)
        "REfuel": 0.2624,
        "e-CH4":  0.2048,

        # Biofuels (combustion emissions not counted)
        "Bio_LF":  0.0, #0.2624,
        "Bio_gas": 0.0, #0.2048,

        # Hydrogen
        "H2": 0.0,
    }

    # accounting_sourcesinkflow index in *this* function uses:
    # (indicator, regionscope, years, sector, commodities)
    # where sector holds your sink name (HeatDemand_*, TranspDemand_*)

    co2_rows = []
    for r in grp.itertuples(index=False):
        fuel = r.commodity
        if fuel not in ef_kt_per_gwh:
            continue

        if r.sector == "Heat":
            ss_sector = f"HeatDemand_{fuel}"
        else:
            ss_sector = f"TranspDemand_{fuel}"

        co2_rows.append(("CO2_emission", "global", int(r.year), ss_sector, fuel, float(ef_kt_per_gwh[fuel])))

    cidx = pd.MultiIndex.from_tuples(
        [(ind, scope, y, sec, c) for (ind, scope, y, sec, c, v) in co2_rows],
        names=["indicator", "regionscope", "years", "sector", "commodities"],
    )
    co2 = pd.DataFrame(index=cidx)
    co2["perFlow"] = [v for (*_, v) in co2_rows]
    m["Base"].parameter.add(co2, "accounting_sourcesinkflow")

    # -------------------------
    # STEP 6: Slack sources (positive-only, penalised)
    # slack index: (node, year, sector, commodity) where sector = SlackFuel_<commodity>
    # -------------------------
    print("\n--- ADDING FUEL SLACK (positive-only, penalised) ---")

    slack_cost = 1e6

    slack_rows = []
    for n in nodes:
        for y in years:
            for f in fuels:
                slack_rows.append((n, int(y), f"SlackFuel_{f}", f, 0.0, np.inf))

    sidx = pd.MultiIndex.from_tuples(
        [(n, y, sec, c) for (n, y, sec, c, lo, up) in slack_rows],
        names=["nodesdata", "years", "sector", "commodities"],
    )
    slack = pd.DataFrame(index=sidx)
    slack["lower"] = 0.0
    slack["upper"] = np.inf
    m["Base"].parameter.add(slack, "sourcesink_annualsum")

    slack_cfg = pd.DataFrame(index=slack.index)
    slack_cfg["usesLowerSum"] = 1
    slack_cfg["usesUpperSum"] = 1
    m["Base"].parameter.add(slack_cfg, "sourcesink_config")

    cost_rows = []
    for y in years:
        for f in fuels:
            sec = f"SlackFuel_{f}"
            cost_rows.append(("SlackCost", "global", int(y), sec, f, float(slack_cost)))

    cidx = pd.MultiIndex.from_tuples(
        [(ind, scope, y, sec, c) for (ind, scope, y, sec, c, v) in cost_rows],
        names=["indicator", "regionscope", "years", "sector", "commodities"],
    )
    cost = pd.DataFrame(index=cidx)
    cost["perFlow"] = [v for (*_, v) in cost_rows]
    m["Base"].parameter.add(cost, "accounting_sourcesinkflow")

    print("[fuel demand] slack sinks written:", len(slack))
    print("[fuel demand] SlackCost perFlow:", slack_cost)

    return {"nodes": nodes, "years": years, "fuels": fuels}

def add_purchasable_fuels(m, fuels_importable=None):
    """
    Separate function: create IMPORT sourcesinks for purchasable fuels and attach CSIRO price.

    This matches the structure you described:
      - fuel demand is loaded elsewhere 
      - slack is loaded elsewhere 
      - THIS function adds supply capability via imports (purchasable fuels)

    It creates ONE import per (node,year,fuel):
      sourcesink: FuelImport_<fuel>
      commodity : <fuel>
      bounds    : 0 .. +inf (or 0..0 if not importable)
      cost      : CSIRO FuelCost (M€/GWh) via accounting_sourcesinkflow

    IMPORTANT:
      - For REfuels that are endogenous later, set them NOT importable
        (either remove from fuels_importable, or pass fuels_importable explicitly).
    """
    import numpy as np
    import pandas as pd
    global yrs_sel

    years_sel = sorted(int(y) for y in yrs_sel)
    nodes = sorted(list(m["Base"].set.nodesdata))
    all_fuels = sorted(list(m["Base"].set.commodities))

    # Default importable fuels: the four MVP fuels if present
    default_importable = {"Bio_gas", "Bio_LF", "Fossil_CH4", "Fossil_LF"}
    if fuels_importable is None:
        fuels_importable = {f for f in default_importable if f in all_fuels}
    else:
        fuels_importable = set(fuels_importable)

    # CSIRO key map 
    fuel_to_csiro = {
        "Fossil_CH4": "Gas",
        "Bio_gas": "Gas",
        "Fossil_LF": "Liquid Fuel",
        "Bio_LF": "Liquid Fuel",
    }

    fuels_for_import = [f for f in all_fuels if f in fuel_to_csiro]  # only fuels we can price
    if not fuels_for_import:
        print("\n=== add_purchasable_fuels_from_csiro ===")
        print("  No fuels with CSIRO mapping found in commodities set. Skipping imports.")
        print("=== end add_purchasable_fuels_from_csiro ===\n")
        return

    print("\n=== add_purchasable_fuels_from_csiro ===")
    print("  fuels_for_import:", fuels_for_import)
    print("  fuels_importable (allowed):", sorted(list(fuels_importable)))

    # Ensure commodities exist
    m["Base"].set.add(fuels_for_import, "commodities")

    # 1) sourcesink_annualsum for imports
    imp_rows = []
    for node in nodes:
        for year in years_sel:
            for fuel in fuels_for_import:
                ss = f"FuelImport_{fuel}"
                lower = 0.0
                upper = np.inf if fuel in fuels_importable else 0.0
                imp_rows.append((node, int(year), ss, fuel, lower, upper))

    imp_idx = pd.MultiIndex.from_tuples(
        [(n, y, ss, c) for (n, y, ss, c, _, _) in imp_rows],
        names=["nodesdata", "years", "sourcesinks", "commodities"],
    )
    imp = pd.DataFrame(index=imp_idx)
    imp["lower"] = [lo for (*_, lo, _) in imp_rows]
    imp["upper"] = [up for (*_, _, up) in imp_rows]
    m["Base"].parameter.add(imp, "sourcesink_annualsum")

    imp_cfg = pd.DataFrame(index=imp.index)
    imp_cfg["usesLowerSum"] = 1
    imp_cfg["usesUpperSum"] = 1
    m["Base"].parameter.add(imp_cfg, "sourcesink_config")

    # 2) accounting_sourcesinkflow FuelCost using CSIRO
    cost_rows = []
    for year in years_sel:
        for fuel in fuels_for_import:
            cs_key = fuel_to_csiro[fuel]
            cval = csiro_fuel_cost_meur_per_gwh(int(year), cs_key)
            ss = f"FuelImport_{fuel}"
            cost_rows.append(("FuelCost", "global", int(year), ss, fuel, float(cval)))

    cost_idx = pd.MultiIndex.from_tuples(
        [(ind, scope, y, ss, c) for (ind, scope, y, ss, c, _) in cost_rows],
        names=["indicator", "regionscope", "years", "sourcesinks", "commodities"],
    )
    cost = pd.DataFrame(index=cost_idx)
    cost["perFlow"] = [v for (*_, v) in cost_rows]
    m["Base"].parameter.add(cost, "accounting_sourcesinkflow")

    print("  wrote FuelImport_<fuel> sourcesinks on same fuel commodities")
    print("=== end add_purchasable_fuels_from_csiro ===\n")

def add_commodity_slack(m, commodities, slack_cost_meur_per_gwh=1e6):
    """
    Adds a positive slack SOURCE for each commodity:
      - sourcesink_annualsum: lower=0, upper=+inf
      - sourcesink_config: usesUpperSum=1, usesLowerSum=1
      - accounting_sourcesinkflow: SlackCost perFlow = slack_cost_meur_per_gwh

    This makes the model always feasible, but very expensive to use.
    Also forces commodities into the commodity set so they show in commodity_balance.
    """

    nodes = sorted(list(m["Base"].set.nodesdata))
    years = sorted(int(y) for y in m["Base"].set.yearssel)

    # Ensure these commodities exist in the set so they appear in results tables
    m["Base"].set.add(list(commodities), "commodities")

    rows = []
    for n in nodes:
        for y in years:
            for c in commodities:
                rows.append((n, y, f"Slack_{c}", c, 0.0, float("inf")))

    idx_ss = pd.MultiIndex.from_tuples(
        [(n, y, ss, c) for (n, y, ss, c, lo, up) in rows],
        names=["nodesdata", "years", "sourcesinks", "commodities"],
    )

    df = pd.DataFrame(index=idx_ss)
    df["lower"] = 0.0
    df["upper"] = float("inf")
    m["Base"].parameter.add(df, "sourcesink_annualsum")

    cfg = pd.DataFrame(index=idx_ss)
    cfg["usesLowerSum"] = 1
    cfg["usesUpperSum"] = 1
    m["Base"].parameter.add(cfg, "sourcesink_config")

    # SlackCost accounting (indicator, regionscope, years, sourcesinks, commodities)
    acc_rows = []
    for (n, y, ss, c, *_rest) in rows:
        acc_rows.append(("SlackCost", "global", y, ss, c, slack_cost_meur_per_gwh))

    idx_acc = pd.MultiIndex.from_tuples(
        [(ind, scope, y, ss, c) for (ind, scope, y, ss, c, v) in acc_rows],
        names=["indicator", "regionscope", "years", "sourcesinks", "commodities"],
    )
    acc = pd.DataFrame(index=idx_acc)
    acc["perFlow"] = [v for (*_, v) in acc_rows]
    m["Base"].parameter.add(acc, "accounting_sourcesinkflow")


    print("[slack] added slack for:", list(commodities), "cost:", slack_cost_meur_per_gwh)

    ss = m["Base"].parameter.sourcesink_annualsum
    cfg = m["Base"].parameter.sourcesink_config

    # show all slack rows for REfuel
    mask = ss.index.get_level_values("sourcesinks").astype(str).str.contains("Slack") & \
        (ss.index.get_level_values("commodities") == "REfuel")

    print(ss.loc[mask].head(20))
    print(cfg.loc[mask].head(20))



# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    start_time = time.time()

    print("\nBuilding REMix data for case:", case_name)
    m = {"Base": Instance(index_names=False, datadir=data_dir)}

    # build the model
    add_scope(m)
    add_demand(m)

    include_heat_transp = True
    if include_heat_transp:
        add_fuel_demands_from_excel(m)
        add_purchasable_fuels(m)
        add_methanizer(m)
        add_ftropsch_syn(m) 
        add_dac(m)
        #add_commodity_slack(m,commodities=["e-CH4", "REfuel", "H2", "CO2_feed"], slack_cost_meur_per_gwh=1e6)

    if include_h2_techs:
        add_electrolyser(m)
        add_H2_FC(m)
        add_H2_CCGT(m)
        add_h2_storage(m)

    if include_hydro_simple:
        add_hydro(m)
    
    if include_geothermal:
        add_geothermal(m)

    if include_thermal:
        add_thermal(m, path_brownfield)
        add_gas_turbines(m, path_brownfield)


    if include_lithium_batteries:
        add_lithium_batteries(m)


    if include_renewables:
        add_renewables(m)

    add_network(m)
    add_accounting(m)
    validate_scope(m)

    # write base data
    m["Base"].write(project_path=data_dir, fileformat="csv", float_format="{:.4g}".format)

    total = time.time() - start_time
    print("\nDone !!")
    print("Wrote inputs to:", data_dir)
    print(f"Total time: {total:.1f} s")