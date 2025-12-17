# %%
# REMix-NZ input builder 
#
# This script writes CSV input tables for the GAMS REMix-NZ model.

import time
import numpy as np
import pandas as pd
from pathlib import Path
from remix.framework.api.instance import Instance

idx = pd.IndexSlice

# -------------------------------------------------------------------
# USER CHOICES
# -------------------------------------------------------------------

group_name = "GP-NT-ELEC-BIO-H2" # Demand folder group (used for output folder structure)
base_scenario = "GP" # Scenario name: GP, NT, ELEC+, BIO+, H2+ (must match "Scenario" column in CSV/Excel)
yrs_to_calc = [2020, 2050] # Model years available
yrs_sel = [2050] # Model years to optimise

# Case folder name (output goes to ../project/<group_name>/<case_name>/data)
years_tag = "-".join(str(y) for y in yrs_sel)
case_name = f"nz_case_{base_scenario}_{years_tag}"


demand_folder = Path("C:/Local/REMix/remix_nz/input/demand/GP-NT-ELEC-BIO-H2") # Demand input folder (contains hourly electricity CSV and the Excel workbook)
hourly_electricity_file = "hourly_electricity_GP-NT-ELEC-BIO-H2.csv"
carriers_excel_file = "data_summary.xlsx"   # must contain sheet "carriers"

# Hydro inflows
use_hydro_cascade = False  # MVP: False
simplified_hydro_inflow_file = Path("C:/Local/REMix/remix_nz/input/brownfield/hydro/inflows_remix-nz/regional-hydro-inflow.csv")

# Currency conversion for CSIRO fuel prices (AUD/GJ) -> EUR
eur_per_aud = 0.60

# Carbon policy (MVP)
use_emission_limit_2050 = False
use_emission_slack_2050 = False
slack_cost_eur_per_tco2 = 1_000_000  # large number only used if cap is infeasible

# Which supply blocks to include
include_renewables = True
include_geothermal = True          
include_hydro_simple = True
include_lithium_batteries = True
include_thermal = True
include_gas_turbines = True
include_electrolyser = True
include_dac = True
include_ftropsch = True
include_methanizer = True
include_h2_storage = True
include_h2_ccgt = True
include_h2_fc = True

include_network = True
include_accounting = True

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------

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
    inflow_file = str(path_brownfield / "hydro" / "inflows_remix-nz" / "inflow_2021-to-2020_2013-to-2030_ENERGY_GW.csv")
else:
    inflow_file = str(simplified_hydro_inflow_file)

print("\n--- scenario setup complete ---")
print("case:", case_name)
print("scenario:", base_scenario)
print("years:", yrs_sel)
print("data_dir:", data_dir.resolve())
print("demand folder:", path_demand.resolve())
print("hydro inflow:", inflow_file)



#mapping carriers from excel 16 regions to 11 remix-nz regions
def map_region_to_remix(region_value: str) -> str:
    """
    Converts region names like "Auckland (AKL)" or "Auckland" into REMix region codes.
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


def read_carriers_excel_long(excel_path: Path, scenario_name: str) -> pd.DataFrame:
    """
    Reads the 'carriers' sheet and returns a long dataframe with columns:
    Scenario, Sector, Carrier, Region, REMixRegion, Year, Demand

    Excludes carriers that should not be treated as final energy carriers in the model:
    Electricity (hourly), Hydrogen (handled separately), Solar, Geothermal.
    """
    df = pd.read_excel(excel_path, sheet_name="carriers")
    df.columns = [str(c).strip() for c in df.columns]

    # filter scenario
    df = df.loc[df["Scenario"].astype(str).str.strip() == str(scenario_name)].copy()
    if df.empty:
        raise ValueError(f"No rows found in carriers Excel for scenario '{scenario_name}'")

    df["REMixRegion"] = df["Region"].apply(map_region_to_remix)
    missing = df.loc[df["REMixRegion"].isna(), "Region"].unique().tolist()
    if missing:
        raise ValueError("Unmapped regions in carriers Excel: " + ", ".join([str(x) for x in missing]))

    # reshape years
    year_cols = [c for c in df.columns if str(c).isdigit()]
    if not year_cols:
        raise ValueError("No year columns found in carriers sheet (expected 2020, 2025, ...).")

    df_long = df.melt(
        id_vars=["Scenario", "Sector", "Carrier", "Region", "REMixRegion"],
        value_vars=year_cols,
        var_name="Year",
        value_name="Demand",
    )
    df_long["Year"] = df_long["Year"].astype(int)
    df_long["Demand"] = pd.to_numeric(df_long["Demand"], errors="coerce").fillna(0.0)

    # keep only model years
    df_long = df_long.loc[df_long["Year"].isin(yrs_sel)].copy()

    # drop zeros
    df_long = df_long.loc[df_long["Demand"] != 0.0].copy()

    return df_long

# -------------------------------------------------------------------
# CSIRO fuel price medians (AUD/GJ) used for MVP
# -------------------------------------------------------------------
# Source: CSIRO GenCost fuel price spreadsheets. 
# https://www.aemo.com.au/energy-systems/major-publications/integrated-system-plan-isp/2026-integrated-system-plan-isp/2025-26-inputs-assumptions-and-scenarios
# For each fuel there are multiple plant entries. We use the median value.
# Ranges (min..max) are left in comments for transparency.

csiro_fuel_aud_per_gj = {
    # 2025-26 (used for "early" years like 2020 in MVP): min 2.26, median 5.19, max 7.22
    # 2050-51 (used for 2050+): min 2.26, median 4.09, max 4.38
    "Black Coal_early": 5.19,
    "Black Coal_2050": 4.09,

    # Brown coal: min 1.02, median 1.04, max 1.06 (basically constant)
    "Brown Coal_early": 1.04,
    "Brown Coal_2050": 1.04,

    # Gas 2025-26: min 13.57, median 16.62, max 18.41
    # Gas 2050-51: min 15.43, median 16.65, max 18.45
    "Gas_early": 16.62,
    "Gas_2050": 16.65,

    # Liquid fuel 2025-26: min 31.62, median 31.89, max 32.87
    # Liquid fuel 2050-51: min 31.76, median 32.04, max 33.19
    "Liquid Fuel_early": 31.89,
    "Liquid Fuel_2050": 32.04,

    # Biomass appears once in CSIRO: 0.66
    # We map this to "Wood" in the carriers list.
    "Biomass_early": 0.66,
    "Biomass_2050": 0.66,
}


def aud_per_gj_to_meur_per_gwh(aud_per_gj: float) -> float:
    """
    Converts AUD/GJ to million EUR/GWh using eur_per_aud.
    1 GJ = 0.277777... MWh, so AUD/MWh = AUD/GJ * 3.6
    EUR/MWh = AUD/MWh * eur_per_aud
    M€ / GWh = (EUR/MWh) * 0.001
    """
    aud_per_mwh = aud_per_gj * 3.6
    eur_per_mwh = aud_per_mwh * eur_per_aud
    meur_per_gwh = eur_per_mwh * 0.001
    return float(meur_per_gwh)


def csiro_fuel_cost_meur_per_gwh(year: int, fuel_key: str) -> float:
    """
    Returns a CSIRO-based fuel cost in M€ / GWh for a given year.
    MVP rule: use 2025-26 median for 2020; use 2050-51 median for 2050.
    """
    suffix = "early" if year < 2050 else "2050"
    aud = csiro_fuel_aud_per_gj[f"{fuel_key}_{suffix}"]
    return aud_per_gj_to_meur_per_gwh(aud)



# demand Elec and Water inflows

def add_scope(m):
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

    # Get the data and model regions based on the mapping
    # "set_nodesData.dat"
    m["Base"].set.add(list(sorted(set(m["Base"].map.aggregatenodesmodel.index.get_level_values(0)))), "nodesdata")
    # "set_nodesModel" & "set_nodesModelSel"
    m["Base"].set.add(list(sorted(set(m["Base"].map.aggregatenodesmodel.index.get_level_values(1)))), "nodesmodel")

    # Set the years to be considered in the model and the years to be optimized
    # "set_years"
    m["Base"].set.add(yrs_to_calc, "years")  # must include all years that data is provided for in the model
    # "set_yearsSel"
    m["Base"].set.add(yrs_sel, "yearssel")  # years to be optimised

def add_demand(m):
    """
    Demand loader.

    If group_name is "GP-NT-ELEC-BIO-H2", this uses the MVP demand folder structure:
    - reads hourly electricity CSV
    - filters to base_scenario
    - drops Scenario column
    - sums sectors into Sector="All" to keep the model small
    - does not load hourly hydrogen

    Otherwise, uses previous configuration (base_scenario.csv and base_scenario-h2.csv).

    This function accepts both naming conventions in the hourly CSV:
    - Region/Year/Sector/Carrier
    - node/year/sector/carrier  (your current file)
    """
    print("\n--- ADDING DEMAND AND INFLOWS ---")

    # ---------------------------------------------
    # Electricity demand (hourly)
    # ---------------------------------------------
    if str(group_name).strip() == "GP-NT-ELEC-BIO-H2":
        elec_path = path_demand / hourly_electricity_file
        if not elec_path.exists():
            raise FileNotFoundError(f"Hourly electricity file not found: {elec_path}")

        df = pd.read_csv(elec_path)

        # Normalise column names to a common internal schema
        # Required internal names: Scenario, node, year, sector, carrier
        rename_cols = {}
        # scenario column can be "Scenario" or "scenario"
        if "scenario" in df.columns and "Scenario" not in df.columns:
            rename_cols["scenario"] = "Scenario"
        # node can be node/Region/region
        if "Region" in df.columns:
            rename_cols["Region"] = "node"
        if "region" in df.columns:
            rename_cols["region"] = "node"
        if "node" not in df.columns and "Node" in df.columns:
            rename_cols["Node"] = "node"
        # year can be Year/year
        if "Year" in df.columns:
            rename_cols["Year"] = "year"
        # sector can be Sector/sector
        if "Sector" in df.columns:
            rename_cols["Sector"] = "sector"
        # carrier can be Carrier/carrier
        if "Carrier" in df.columns:
            rename_cols["Carrier"] = "carrier"
        if "Carrier" not in df.columns and "carrier" not in df.columns and "Carrier" in rename_cols:
            pass

        if rename_cols:
            df = df.rename(columns=rename_cols)

        required = {"Scenario", "node", "year", "sector", "carrier"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                "Hourly electricity CSV missing columns after normalisation: "
                f"{missing}. Found columns: {list(df.columns)}"
            )

        # Map node (Region) to REMix region codes if needed
        # If the value is not already an all-caps code, use map_region_to_remix
        df["node"] = df["node"].apply(
            lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r)
        )
        df = df.dropna(subset=["node"])

        # filter scenario and drop Scenario column
        df = df.loc[df["Scenario"].astype(str).str.strip() == str(base_scenario)].copy()
        df.drop(columns=["Scenario"], inplace=True)

        # keep only electricity
        df["carrier"] = df["carrier"].astype(str).str.strip()
        df = df.loc[df["carrier"].str.lower() == "electricity"].copy()

        # keep only model years
        df["year"] = df["year"].astype(int)
        df = df.loc[df["year"].isin(yrs_sel)].copy()

        hour_cols = [c for c in df.columns if str(c).startswith("t")]
        if not hour_cols:
            raise ValueError("No hourly columns found (expected t0001..t8760).")

        # sum all sectors -> All
        df_agg = df.groupby(["node", "year", "carrier"], as_index=False)[hour_cols].sum()
        df_agg["sector"] = "All"

        # map carrier to REMix commodity name
        df_agg["carrier"] = df_agg["carrier"].replace({"Electricity": "Elec", "electricity": "Elec"})

        ts = df_agg.set_index(["node", "year", "sector", "carrier"])[hour_cols]
        ts *= -1
        ts["type"] = "fixed"
        ts_fixed = ts.set_index("type", append=True).round(3)

        m["Base"].profile.add(ts_fixed, "sourcesink_profile")

        cfg = pd.DataFrame(index=ts_fixed.index.droplevel(-1))
        cfg["usesFixedProfile"] = 1
        m["Base"].parameter.add(cfg, "sourcesink_config")

        m["Base"].set.add(sorted(cfg.index.get_level_values(0).unique().tolist()), "nodesdata")
        m["Base"].set.add(sorted(cfg.index.get_level_values(1).unique().tolist()), "years")

        print(f"Electricity demand loaded: {len(cfg)} rows (Sector aggregated to 'All').")
        print("Hydrogen is not loaded hourly in MVP (annual sum handled separately).")

    else:
        # original behaviour
        file_paths = {
            "Elec": str(path_demand / f"{base_scenario}.csv"),
            "H2": str(path_demand / f"{base_scenario}-h2.csv"),
        }
        rename_commodity = {
            "Electricity": "Elec",
            "Hydrogen": "H2",
            "H2-feedstock": "H2",
            "Natural Gas": "CH4",
            "Gas": "CH4",
            "Feedstock Gas": "CH4",
            "Feedstock methanol": "CH3OH",
            "Renewable Fuels": "REfuel",
        }

        for commodity_type, file_path in file_paths.items():
            if not Path(file_path).exists():
                print(f"Skipping {commodity_type} (file not found): {file_path}")
                continue

            ts = pd.read_csv(file_path, index_col=[0, 1, 2, 3]).rename(index=rename_commodity)

            # Map long region names in index level 0 to REMix codes, if needed
            ts = ts.reset_index()
            ts.iloc[:, 0] = ts.iloc[:, 0].apply(
                lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r)
            )
            ts = ts.dropna(subset=[ts.columns[0]])
            ts = ts.set_index(list(ts.columns[:4]))

            ts *= -1
            ts["type"] = "fixed"
            ts_fixed = ts.set_index("type", append=True).round(3)
            m["Base"].profile.add(ts_fixed, "sourcesink_profile")

            ts_cfg = pd.DataFrame(index=ts.index)
            ts_cfg["usesFixedProfile"] = 1
            ts_cfg = ts_cfg.loc[ts.select_dtypes(include="number").sum(axis=1) != 0]
            m["Base"].parameter.add(ts_cfg, "sourcesink_config")

            m["Base"].set.add(list(ts.index.get_level_values(0)), "nodesdata")
            m["Base"].set.add(list(ts.index.get_level_values(1)), "years")

            print(f"{commodity_type} demand loaded: {len(ts)} rows")

    # ---------------------------------------------
    # Hydro inflows (simplified for MVP)
    # ---------------------------------------------
    print("\nProcessing hydro inflows...")
    inflow_df = pd.read_csv(inflow_file)

    inflow_rename = {}
    if "Region" in inflow_df.columns:
        inflow_rename["Region"] = "node"
    if "region" in inflow_df.columns:
        inflow_rename["region"] = "node"
    if "Year" in inflow_df.columns:
        inflow_rename["Year"] = "year"
    if "Carrier" in inflow_df.columns:
        inflow_rename["Carrier"] = "commodity"
    if "carrier" in inflow_df.columns and "commodity" not in inflow_df.columns:
        inflow_rename["carrier"] = "commodity"
    if "Sector" in inflow_df.columns:
        inflow_rename["Sector"] = "sector"

    if inflow_rename:
        inflow_df = inflow_df.rename(columns=inflow_rename)

    required_inflow = {"node", "year", "sector", "commodity"}
    missing_inflow = required_inflow - set(inflow_df.columns)
    if missing_inflow:
        raise ValueError(
            "Hydro inflow CSV missing columns after normalisation: "
            f"{missing_inflow}. Found columns: {list(inflow_df.columns)}"
        )

    # Map inflow node to REMix codes
    inflow_df["node"] = inflow_df["node"].apply(
        lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r)
    )
    inflow_df = inflow_df.dropna(subset=["node"])

    inflow = inflow_df.set_index(["node", "year", "sector", "commodity"])
    inflow = inflow.select_dtypes(include=np.number).dropna(how="all")

    # Convert HydroInflow commodity name to Waterin if needed
    inflow.index = inflow.index.set_levels(
        ["Waterin" if x == "HydroInflow" else x for x in inflow.index.levels[3]],
        level=3,
    )

    inflow *= -1
    inflow["type"] = "fixed"
    inflow_fixed = inflow.set_index("type", append=True).round(3)
    m["Base"].profile.add(inflow_fixed, "sourcesink_profile")

    inflow_cfg = pd.DataFrame(index=inflow.index)
    inflow_cfg["usesFixedProfile"] = 1
    inflow_cfg = inflow_cfg.loc[inflow.select_dtypes(include=np.number).sum(axis=1) != 0]
    m["Base"].parameter.add(inflow_cfg, "sourcesink_config")

    m["Base"].set.add(sorted(inflow.index.get_level_values(0).unique().tolist()), "nodesdata")
    m["Base"].set.add(sorted(inflow.index.get_level_values(1).unique().tolist()), "years")

    print("Hydro inflows loaded.")

# hydro
def add_hydro(m):
    # Define hydropower converters (turbines) and storage (reservoirs) for REMix NZ.
    # Flow sequence: inflow (Water_in) -> reservoir -> turbine -> electricity.
    # All quantities are in GWh-equivalent units.

    # Configuration for water balance and flow profiles
    sourcesink_cfg = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [m["Base"].set.nodesdata, m["Base"].set.yearssel, ["Ocean"], ["Water_out"]]
        )
    )
    sourcesink_cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(sourcesink_cfg, "sourcesink_config")

    hydro_vintage = [1950]
    hydro_years = [2000] + yrs_to_calc
    hydro_nodes = ["BOP", "CAN", "CEN", "HBY", "NEL", "OTG", "WTO"]
    hydro_techs = ["Hydro"]
    hydro_activities = ["Power_gen", "Spill"]

    # Turbine technology settings
    conv_tech = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_techs, hydro_vintage]))
    conv_tech["lifeTime"] = 100
    conv_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # Installed turbine capacities (GW_el)
    conv_cap = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, hydro_techs]))
    conv_cap.loc[idx[["BOP"], [2000], "Hydro"], "unitsBuild"] = 0.17095
    conv_cap.loc[idx[["CAN"], [2000], "Hydro"], "unitsBuild"] = 1.82683
    conv_cap.loc[idx[["CEN"], [2000], "Hydro"], "unitsBuild"] = 0.399
    conv_cap.loc[idx[["HBY"], [2000], "Hydro"], "unitsBuild"] = 0.1422
    conv_cap.loc[idx[["NEL"], [2000], "Hydro"], "unitsBuild"] = 0.0453
    conv_cap.loc[idx[["OTG"], [2000], "Hydro"], "unitsBuild"] = 1.664
    conv_cap.loc[idx[["WTO"], [2000], "Hydro"], "unitsBuild"] = 1.0873
    conv_cap["noExpansion"] = 1
    conv_cap["unitsUpperLimit"] = 100
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # Turbine input/output coefficients
    hydro_efficiency = 0.95
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [hydro_techs, hydro_vintage, hydro_activities, ["Water_in", "Water_out", "Elec"]]
        )
    )
    conv_coef.loc[idx[:, :, "Power_gen", "Elec"], "coefficient"] = 1.0
    conv_coef.loc[idx[:, :, "Power_gen", "Water_in"], "coefficient"] = -hydro_efficiency
    conv_coef.loc[idx[:, :, "Power_gen", "Water_out"], "coefficient"] = hydro_efficiency
    conv_coef.loc[idx[:, :, "Spill", "Water_in"], "coefficient"] = -100.0
    conv_coef.loc[idx[:, :, "Spill", "Water_out"], "coefficient"] = 100.0
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # Turbine investment and fixed O&M
    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], hydro_techs, hydro_vintage]
        )
    )
    conv_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] = 2560
    conv_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1
    conv_acc.loc[idx["Invest", :, "horizon", :, :], "amorTime"] = 20
    conv_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06

    inv = conv_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    inv.index = pd.MultiIndex.from_tuples([("OMFix", *i[1:]) for i in inv.index])
    conv_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = inv * 0.03
    m["Base"].parameter.add(conv_acc, "accounting_converterunits")

    # Reservoir technology and size
    stor_techs = ["Hydro_reservoir"]
    stor_tech = pd.DataFrame(index=pd.MultiIndex.from_product([stor_techs, hydro_vintage]))
    stor_tech["lifeTime"] = 100
    stor_tech["levelUpperLimit"] = 1
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    stor_size = pd.DataFrame(index=pd.MultiIndex.from_product([stor_techs, hydro_vintage, ["Water_in"]]))
    stor_size["size"] = 1
    stor_size["selfdischarge"] = 0
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # Reservoir installed capacities (GWh water potential)
    stor_res = pd.DataFrame(index=pd.MultiIndex.from_product([hydro_nodes, hydro_years, stor_techs]))
    stor_res.loc[idx[["CAN"], [2000], :], "unitsBuild"] = 2517.2429
    stor_res.loc[idx[["HBY"], [2000], :], "unitsBuild"] = 154.2635
    stor_res.loc[idx[["OTG"], [2000], :], "unitsBuild"] = 729.5595
    stor_res.loc[idx[["WTO"], [2000], :], "unitsBuild"] = 587.1371
    stor_res["noExpansion"] = 1
    stor_res["unitsUpperLimit"] = 3000
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    # Reservoir investment and O&M
    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], stor_techs, hydro_vintage]
        )
    )
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] = 1650
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "amorTime"] = 20
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06

    inv = stor_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    inv.index = pd.MultiIndex.from_tuples([("OMFix", *i[1:]) for i in inv.index])
    stor_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = inv * 0.03
    m["Base"].parameter.add(stor_acc, "accounting_storageunits")

def add_hydro_scheme(m):
    """
    Detailed NZ hydropower cascade using JADE data for capacities and storage,
    but keeping the original REMix cascade naming and regions.

    Flow: inflow (*_in) -> [reservoir] -> turbine (Power_gen / Spill) -> next *_in ... -> Water_out
    Units:
      - Water and storage are handled in *energy units (GWh)* via specific power (MW per cumec).
      - Mm³ from JADE files are converted to GWh with the downstream turbine's specific power.
    """
    jade_dir = Path(path_brownfield) / "hydro" / "JADE_hydro"
    f_stations = jade_dir / "hydro_stations.csv"
    f_limits   = jade_dir / "reservoir_limits.csv"
    f_resinit  = jade_dir / "reservoirs.csv"

    try:
        years = list(m["Base"].set.yearssel)
        if not years:
            years = [2000]
    except Exception:
        years = [2000]
    y0 = years[0]

    # Columns: inflow_list (comma sep), tech, outflow, scheme, region, (cap, sp filled from JADE)
    turbines = pd.DataFrame([
        # Waitaki (CAN)
        ("Tekapo_in"          , "Tekapo_A"    , "Tekapo_A_out" , "Waitaki", "CAN"),
        ("Tekapo_A_out"       , "Tekapo_B"    , "Pukaki_in"    , "Waitaki", "CAN"),
        ("Pukaki_in, Ohau_in" , "Ohau_A"      , "Ohau_A_out"   , "Waitaki", "CAN"),
        ("Ohau_A_out"         , "Ohau_B"      , "Ohau_BC_canal", "Waitaki", "CAN"),
        ("Ohau_BC_canal"      , "Ohau_C"      , "Benmore_in"   , "Waitaki", "CAN"),
        ("Benmore_in"         , "Benmore"     , "Aviemore_in"  , "Waitaki", "CAN"),
        ("Aviemore_in"        , "Aviemore"    , "Waitaki_in"   , "Waitaki", "CAN"),
        ("Waitaki_in"         , "Waitaki"     , "Water_out"    , "Waitaki", "CAN"),

        # Waikato (WTO)
        ("Aratiatia_in"       , "Aratiatia"   , "Ohakuri_in"   , "Waikato", "WTO"),
        ("Ohakuri_in"         , "Ohakuri"     , "Atiamuri_in"  , "Waikato", "WTO"),
        ("Atiamuri_in"        , "Atiamuri"    , "Whakamaru_in" , "Waikato", "WTO"),
        ("Whakamaru_in"       , "Whakamaru"   , "Maraetai_in"  , "Waikato", "WTO"),
        ("Maraetai_in"        , "Maraetai"    , "Waipapa_in"   , "Waikato", "WTO"),
        ("Waipapa_in"         , "Waipapa"     , "Arapuni_in"   , "Waikato", "WTO"),
        ("Arapuni_in"         , "Arapuni"     , "Karapiro_in"  , "Waikato", "WTO"),
        ("Karapiro_in"        , "Karapiro"    , "Water_out"    , "Waikato", "WTO"),

        # Clutha (OTG)
        ("Dunstan_in"         , "Clyde_220kV" , "Roxburgh_in"  , "Clutha" , "OTG"),
        ("Roxburgh_in"        , "Roxburgh"    , "Water_out"    , "Clutha" , "OTG"),

        # Manapouri (OTG)
        ("Manapouri_in"       , "Manapouri"   , "Water_out"    , "Manapouri","OTG"),

        # Waikaremoana (HBY)
        ("Waikaremoana_in"    , "Waikaremoana", "Water_out"    , "Waikaremoana","HBY"),

        # Singles (CEN, CAN, NEL)
        ("Rangipo_in"         , "Rangipo"     , "Water_out"    , "Rangipo",   "CEN"),
        ("Tokaanu_in"         , "Tokaanu"     , "Water_out"    , "Tokaanu",   "CEN"), #tongariro
        ("Matahina_in"        , "Matahina"    , "Water_out"    , "Matahina",  "CEN"),
        ("Mangahao_in"        , "Mangahao"    , "Water_out"    , "Mangahao",  "CEN"),
        ("Coleridge_in"       , "Coleridge"   , "Water_out"    , "Coleridge", "CAN"),
        ("Cobb_in"            , "Cobb"        , "Water_out"    , "Cobb",      "NEL"),
    ], columns=["inflow_list","tech","outflow","scheme","region"])

    # Reservoirs
    reservoirs = pd.DataFrame([
        ("Lake_Tekapo"              , "Tekapo_in"       , "CAN"),
        ("Lake_Pukaki"              , "Pukaki_in"       , "CAN"),
        ("Lake_Ohau"                , "Ohau_in"         , "CAN"),
        ("Lake_Taupo"               , "Aratiatia_in"    , "WTO"),
        ("Lake_Hawea"               , "Dunstan_in"      , "OTG"),
        ("Lakes_Manapouri_Te_Anau"  , "Manapouri_in"    , "OTG"),
        ("Lake_Waikaremoana"        , "Waikaremoana_in" , "HBY"),
    ], columns=["Reservoir","stores_commodity","Region"])

    # Ensure nodes are present
    nodes = sorted(set(turbines["region"].unique()) | set(reservoirs["Region"].unique()))
    m["Base"].set.add(nodes, "nodesdata")

    # --------- load JADE data (stations, limits, initial) ----------
    # hydro_stations: GENERATOR, CAPACITY (MW), SPECIFIC_POWER (MW per cumec)
    stations = pd.read_csv(f_stations)
    stations = stations.rename(columns=str.strip)

    # Build lookup for capacity and Specific Power by turbine name
    cap_map = stations.set_index("GENERATOR")["CAPACITY"].to_dict()
    sp_map  = stations.set_index("GENERATOR")["SPECIFIC_POWER"].to_dict()

    # Build mapping inflow commodity -> downstream turbine tech 
    inflow_to_turb = {}
    for _, r in turbines.iterrows():
        for c in [x.strip() for x in r["inflow_list"].split(",") if x.strip()]:
            inflow_to_turb.setdefault(c, r["tech"])

    # Helper: convert Mm3 to GWh given a specific power (MW per cumec)
    def mm3_to_gwh(volume_mm3, sp_mw_per_cumec):
        if sp_mw_per_cumec is None or pd.isna(sp_mw_per_cumec):
            return 0.0
        # 1 m³ = 1/3600 cumec-hour
        cumec_hours = (volume_mm3 * 1e6) / 3600.0
        return (sp_mw_per_cumec * cumec_hours) / 1000.0  # → GWh

    # reservoir_limits: take MAX of *_MAX_LEVEL columns across file
    limits = pd.read_csv(f_limits)
    limits = limits.rename(columns=str.strip)
    max_cols = [c for c in limits.columns if c.endswith("MAX_LEVEL")]
    max_by_res = {}
    for c in max_cols:
        # Column like "Lake_Tekapo MAX_LEVEL" → reservoir key "Lake_Tekapo"
        rname = c.replace(" MAX_LEVEL", "")
        max_by_res[rname] = float(limits[c].max())

    # reservoirs.csv for initial storage (INI_STATE in Mm3)
    rinit = pd.read_csv(f_resinit)
    rinit = rinit.rename(columns=str.strip)
    init_mm3_map = dict(zip(rinit["RESERVOIR"], rinit["INI_STATE"]))

    # Compute reservoir E_max_GWh and initialLevel (fraction), using downstream turbine SP
    def reservoir_sp(res_row):
        inflow_comm = res_row["stores_commodity"]  # e.g. "Tekapo_in"
        turb = inflow_to_turb.get(inflow_comm, None)  # e.g. "Tekapo_A"
        return sp_map.get(turb, np.nan)

    reservoirs["Max_Level_Mm3"]  = reservoirs["Reservoir"].map(max_by_res).fillna(0.0)
    reservoirs["Ini_Level_Mm3"]  = reservoirs["Reservoir"].map(init_mm3_map).fillna(0.0)
    reservoirs["SP_downstream"]  = reservoirs.apply(reservoir_sp, axis=1)

    reservoirs["E_max_GWh"]  = reservoirs.apply(
        lambda r: mm3_to_gwh(r["Max_Level_Mm3"], r["SP_downstream"]), axis=1
    ).fillna(0.0)
    reservoirs["E_init_GWh"] = reservoirs.apply(
        lambda r: mm3_to_gwh(r["Ini_Level_Mm3"], r["SP_downstream"]), axis=1
    ).fillna(0.0)
    reservoirs["InitFrac"]   = np.where(reservoirs["E_max_GWh"] > 0,
                                        (reservoirs["E_init_GWh"] / reservoirs["E_max_GWh"]).clip(0,1),
                                        0.0)

    # --------- converter tech params ----------
    hydro_vintage = [1950]
    hydro_eff = 1.0  # pass all water downstream
    activities = ["Power_gen", "Spill", "Power_gen_b"]  # special split for Ohau_A
    techs = turbines["tech"].tolist()

    conv_tech = pd.DataFrame(index=pd.MultiIndex.from_product([techs, hydro_vintage]))
    conv_tech["lifeTime"] = 100
    conv_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # --------- capacities (from JADE capacity, at first model year) ----------
    cap = pd.DataFrame(index=pd.MultiIndex.from_product([nodes, years, techs]))
    for _, r in turbines.iterrows():
        mw = cap_map.get(r["tech"], np.nan)
        if pd.isna(mw):
            continue
        cap.loc[(r["region"], y0, r["tech"]), "unitsBuild"] = float(mw) / 1000.0  # to GW
    cap["noExpansion"] = 1
    cap["unitsUpperLimit"] = 100
    m["Base"].parameter.add(cap.dropna(how="all"), "converter_capacityparam")

    # --------- coefficients (generation + spill) ----------
    # Collect all commodities that appear
    all_inflows = set()
    for s in turbines["inflow_list"]:
        all_inflows |= {x.strip() for x in s.split(",") if x.strip()}
    all_outflows = set(turbines["outflow"].tolist())
    commodities = sorted(all_inflows | all_outflows | {"Elec"})

    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, hydro_vintage, activities, commodities]),
        columns=["coefficient"],
        data=0.0
    )

    BIG = 1e6  # unbounded spill

    for _, r in turbines.iterrows():
        inflows = [x.strip() for x in r["inflow_list"].split(",") if x.strip()]

        if r["tech"] != "Ohau_A":
            # Power_gen: +Elec, -inflows, +outflow
            coef.loc[(r["tech"], 1950, "Power_gen", "Elec"), "coefficient"] = 1.0
            for c_in in inflows:
                coef.loc[(r["tech"], 1950, "Power_gen", c_in), "coefficient"] = -hydro_eff
            coef.loc[(r["tech"], 1950, "Power_gen", r["outflow"]), "coefficient"] = hydro_eff
        else:
            # Ohau_A special split: one activity per inflow
            for c_in in inflows:
                act = "Power_gen" if c_in == "Pukaki_in" else "Power_gen_b"
                coef.loc[(r["tech"], 1950, act, "Elec"), "coefficient"] = 1.0
                coef.loc[(r["tech"], 1950, act, c_in), "coefficient"] = -hydro_eff
                coef.loc[(r["tech"], 1950, act, r["outflow"]), "coefficient"] = hydro_eff

        # Spill: route water without energy with a very loose bound
        for c_in in inflows:
            coef.loc[(r["tech"], 1950, "Spill", c_in), "coefficient"] = -BIG
        coef.loc[(r["tech"], 1950, "Spill", r["outflow"]), "coefficient"] = BIG

    coef = coef[coef["coefficient"] != 0.0]
    m["Base"].parameter.add(coef, "converter_coefficient")

    # --------- storage: one technology per reservoir (size in GWh) ----------
    stor_names = [f"Stor_{r}" for r in reservoirs["Reservoir"]]
    stor_tech = pd.DataFrame(index=pd.MultiIndex.from_product([stor_names, hydro_vintage]))
    stor_tech["lifeTime"] = 100
    stor_tech["levelUpperLimit"] = 1.0
    stor_tech["levelLowerLimit"] = 0.0
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    # size param (tech,year,commodity) = 1, selfdischarge=0
    size_idx = [(f"Stor_{r.Reservoir}", 1950, r.stores_commodity) for _, r in reservoirs.iterrows()]
    stor_size = pd.DataFrame(index=pd.MultiIndex.from_tuples(size_idx, names=["techs","years","commodities"]))
    stor_size["size"] = 1.0
    stor_size["selfdischarge"] = 0.0
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # units (node,year,tech) carry the energy capacity in GWh, and initial level fraction
    stor_units = pd.DataFrame(index=pd.MultiIndex.from_product([nodes, years, stor_names]))
    for _, r in reservoirs.iterrows():
        tech = f"Stor_{r.Reservoir}"
        emx  = float(r.E_max_GWh)
        init = float(r.InitFrac)
        stor_units.loc[(r.Region, y0, tech), "unitsBuild"]   = emx
        stor_units.loc[(r.Region, y0, tech), "initialLevel"] = init
    stor_units["noExpansion"] = 1
    m["Base"].parameter.add(stor_units, "storage_reservoirparam")

    # --------- ocean sink  ----------
    ocean_cfg = pd.DataFrame(index=pd.MultiIndex.from_product([nodes, years, ["Ocean"], ["Water_out"]]))
    ocean_cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(ocean_cfg, "sourcesink_config")

    # --------- prints ----------
    # sanity: warn if any reservoir has zero capacity (likely SP lookup issue)
    zeros = reservoirs.loc[reservoirs["E_max_GWh"] <= 0, "Reservoir"].tolist()
    if zeros:
        print("WARNING: reservoirs with 0 GWh capacity (check SP mapping):", ", ".join(zeros))

    print("\nHydro scheme (JADE) loaded successfully.")
    print(f"  Turbines: {len(techs)} | Reservoirs: {len(reservoirs)} | Regions: {', '.join(sorted(set(nodes)))}")
    print("  Reservoir energy [GWh] (init %):")
    for _, r in reservoirs.iterrows():
        emx = float(r.E_max_GWh)
        init_pct = 100.0 * float(r.InitFrac) if emx > 0 else 0.0
        print(f"   - {r.Reservoir:<26s} {emx:6.1f} GWh  (init {init_pct:5.1f}%)")

# other renewables
    
def load_feedin_csv(year, aggregate=False, norm=True):
    inst = pd.read_csv(Path(path_profiles).joinpath("region_statistics_2012.csv"), index_col=[0, 1])
    ts = pd.read_csv(Path(path_profiles).joinpath(f"timeseries_{year}_NZT.csv"), index_col=[0, 1, 2])

    #C: we adjusted plus half an hour
    # Split year and time index
    ts["year"] = (pd.DatetimeIndex(ts.index.get_level_values(0)) + pd.Timedelta(hours=0.5)).year
    ts["t_model"] = ((pd.DatetimeIndex(ts.index.get_level_values(0)) + pd.Timedelta(hours=0.5)).dayofyear - 1) * 24 + (pd.DatetimeIndex(ts.index.get_level_values(0)) - pd.Timedelta(hours=1)).hour + 1
    ts = ts.reset_index().set_index(["region", "technology", "year", "t_model"]).drop(columns=["t"])

    if aggregate:
        inst = inst.groupby(["technology"]).sum()
        ts = ts.groupby(["technology", "year", "t_model"]).sum()

    if norm:
        ts = ts.unstack(["year", "t_model"])
        ts = ts.loc[idx[inst.index]]

        ts = pd.DataFrame(
            ts.values / pd.DataFrame(inst["installable_per_region"]).values,
            index=ts.index,
            columns=ts.columns,
        )
        ts = ts.stack(["year", "t_model"])

    ts[ts > 1] = 1
    return ts.round(3)

def add_renewables(m):
    
    re_inst_csv = pd.read_csv(Path(path_profiles).joinpath("region_statistics_2012.csv"), index_col=[0, 1])
    re_nodes = [n for n in m["Base"].set.nodesdata]

    re_vintage = [1950,2030,2040,2050]
    year_mapping = {
        2012: 2020,
        # 2011: 2025,
        # 2012: 2030,
        #2012: 2035,
        #2014: 2040,
        #2016: 2045,
        #2012: 2050 #this mapping is only used with different weather years
    }
    re_techs = list(set(re_inst_csv.index.get_level_values(1)))
    re_techs = [tech for tech in re_techs if tech != "pv_central_track_azimuth"]
    pv_techs = [i for i in re_techs if i.startswith("pv") and not i.startswith("pv_central_track_azimut") ]
    csp_techs = []#[i for i in re_techs if i.startswith("csp")]
    wind_techs = [i for i in re_techs if i.startswith("wind")]
    
    re_tech = pd.DataFrame(index=pd.MultiIndex.from_product([re_techs, re_vintage]))
    re_tech.loc[idx[:, :], "activityUpperLimit"] = 0 # so it is overwritten by the availabity from the timeseriesfiles
    re_tech.loc[idx[pv_techs, [1950]], "lifeTime"] = 35  # years
    re_tech.loc[idx[pv_techs, [2030, 2040, 2050]], "lifeTime"] = 40  # years
    re_tech.loc[idx[csp_techs + wind_techs, [1950]], "lifeTime"] = 27  # years
    re_tech.loc[idx[csp_techs + wind_techs, [2030, 2040, 2050]], "lifeTime"] = 30  # years
    m["Base"].parameter.add(re_tech, "converter_techparam")
    
    # capacities
    years_build = list(range(1989, 2021)) #brownfield

    re_caps = pd.DataFrame(index=pd.MultiIndex.from_product([re_nodes, years_build, re_techs]))
    re_caps.index.names = ["region", "years", "technology"]

    re_caps.loc[idx["CAN", [2003], "wind_onshore"], "unitsBuild"] = 0.0005 # GW_el
    re_caps.loc[idx["CAN", [2005], "wind_onshore"], "unitsBuild"] = 0.0001 # GW_el

    re_caps.loc[idx["CEN", [1999], "wind_onshore"], "unitsBuild"] = 31.7  /1000 # GW_el	
    re_caps.loc[idx["CEN", [2004], "wind_onshore"], "unitsBuild"] = 127.05 /1000 # GW_el	
    re_caps.loc[idx["CEN", [2007], "wind_onshore"], "unitsBuild"] = 93 /1000 # GW_el
    re_caps.loc[idx["CEN", [2011], "wind_onshore"], "unitsBuild"] = 48.5 /1000 # GW_el
    re_caps.loc[idx["CEN", [2020], "wind_onshore"], "unitsBuild"] = 221.4 /1000 # GW_el	

    re_caps.loc[idx["OTG", [2007], "wind_onshore"], "unitsBuild"] = 58/1000 # GW_el	
    re_caps.loc[idx["OTG", [2009], "wind_onshore"], "unitsBuild"] = 2.25/1000 # GW_el	
    re_caps.loc[idx["OTG", [2010], "wind_onshore"], "unitsBuild"] = 0.45/1000 # GW_el	
    re_caps.loc[idx["OTG", [2011], "wind_onshore"], "unitsBuild"] = 43.65/1000 # GW_el	
    re_caps.loc[idx["OTG", [2015], "wind_onshore"], "unitsBuild"] = 6.8 /1000 # GW_el
    
    re_caps.loc[idx["NEL", [2010], "wind_onshore"], "unitsBuild"] = 0.75/1000 # GW_el
    re_caps.loc[idx["NEL", [2011], "wind_onshore"], "unitsBuild"] = 1/1000 # GW_el
    re_caps.loc[idx["NEL", [2014], "wind_onshore"], "unitsBuild"] = 0.66/1000 # GW_el

    re_caps.loc[idx["TRN", [2020], "wind_onshore"], "unitsBuild"] = 0.1333 # GW_el

    re_caps.loc[idx["WEL", [1993], "wind_onshore"], "unitsBuild"] = 8.45 /1000# GW_el
    re_caps.loc[idx["WEL", [1996], "wind_onshore"], "unitsBuild"] = 8.45 /1000# GW_el
    re_caps.loc[idx["WEL", [2009], "wind_onshore"], "unitsBuild"] = 143 /1000# GW_el
    re_caps.loc[idx["WEL", [2014], "wind_onshore"], "unitsBuild"] = 71.3 /1000# GW_el	
    re_caps.loc[idx["WTO", [2011], "wind_onshore"], "unitsBuild"] = 64.4/1000 # GW_el

    re_upper = pd.DataFrame(index=pd.MultiIndex.from_product([re_nodes, yrs_to_calc, re_techs]))
    re_upper.index.names = ["region", "years", "technology"]
    re_max_inst = pd.DataFrame(re_inst_csv.div(1e3)["installable_per_region"]).rename(columns={"installable_per_region": "unitsUpperLimit"})
    re_upper = re_upper.join(re_max_inst, on=["region", "technology"], how="outer")
    if 2020 in yrs_sel:
        re_upper.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    re_caps_full = pd.concat([re_caps, re_upper], axis=1).dropna(how="all")

    m["Base"].parameter.add(re_caps_full, "converter_capacityparam")


    # activity
    # FIXME: (PART 1): it is ok while we only have 1 weather year
    re_feedin_csv = load_feedin_csv(2012).unstack("t_model").swaplevel(1, 2)  # Load data for the year 2012
    re_feedin = pd.concat([re_feedin_csv.rename(index={2012:y}) for y in yrs_to_calc])
    re_feedin["type"] = "upper"
    re_feedin = re_feedin.set_index("type", append=True)
    re_feedin = re_feedin[re_feedin >= 0.01].dropna(how="all").fillna(0)
    re_feedin = re_feedin.iloc[:, 0:8760]
    re_feedin.columns = [f"t{str(t+1).zfill(4)}" for t in range(8760)]
    re_feedin = re_feedin.sort_index(level=["region", "technology", "year"])
    m["Base"].profile.add(re_feedin, "converter_activityprofile")

    # coefficients
    re_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                re_techs,
                re_vintage,
                ["Powergen", "Heatgen"],
                ["Elec", "Heat_CSP"],
            ]
        )
    )
    re_coef.loc[idx[wind_techs + pv_techs, :, "Powergen", "Elec"], "coefficient"] = 1
    re_coef = re_coef.replace(0, np.nan).dropna(how="all")
    m["Base"].parameter.add(re_coef, "converter_coefficient")

    re_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ["horizon"], re_techs, re_vintage])
    ).sort_index()

    # TODO: Update costs for renewable technologies
    # CSP own assumptions based on: https://aip.scitation.org/doi/pdf/10.1063/5.0028883, https://elib.dlr.de/186998/1/SolarPACES_2021_Paper_Dersch_R1.pdf
    #re_acc.loc[idx["Invest", :, "horizon", "csp_parabolic_trough", :], "perUnitBuild"] = [344.5, 274.7, 230.2, 196.0]  # Child 2019 - Mio EUR per unit
    #re_acc.loc[idx["Invest", :, "horizon", "csp_solar_tower", :], "perUnitBuild"] = [482, 372, 310, 264]  # Mio EUR per unit

    re_acc.loc[idx["Invest", :, "horizon", "pv_decentral", :], "perUnitBuild"] = [870, 570, 460, 410]  # DEA2022 PV comm&indust - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "horizon", "pv_central_fixed", :], "perUnitBuild"] = [560, 380, 320, 290]  # DEA2022 utility scale - Mio EUR per unit
    #re_acc.loc[idx["Invest", :, "horizon", "pv_central_track_azimuth", :], "perUnitBuild"] = [650, 450, 380, 350]  # DEA2022 utility scale (tracking) - Mio EUR per unit

    re_acc.loc[idx["Invest", :, "horizon", "wind_onshore", :], "perUnitBuild"] = [1330, 1040, 980, 960]  # DEA2022 onshore - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "horizon", "wind_offshore_foundation", :], "perUnitBuild"] = [2120, 2287, 2168, 2130] # DEA2022 offshore - Mio EUR per unit
    re_acc.loc[idx["Invest", :, "horizon", "wind_offshore_floating", :], "perUnitBuild"] = 1.2 * np.array([2120, 2287, 2168, 2130])  # DEA2022 offshore + 20% assumption - Mio EUR per unit

    re_acc.loc[idx["Invest", :, "horizon", pv_techs, [1950]], "amorTime"] = 35  # years
    re_acc.loc[idx["Invest", :, "horizon", pv_techs, [2030, 2040, 2050]], "amorTime"] = 40  # years
    re_acc.loc[idx["Invest", :, "horizon", csp_techs + wind_techs, [1950]], "amorTime"] = 27  # years
    re_acc.loc[idx["Invest", :, "horizon", csp_techs + wind_techs, [2030, 2040, 2050]], "amorTime"] = 30  # years

    re_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1  # binary yes/no
    re_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06  # percent/100
    # re_acc.loc[idx["OMFix", :, "horizon", pv_techs + wind_techs, :], "perUnitTotal"] = (
    #     re_acc.loc[idx["Invest", :, "horizon", pv_techs + wind_techs, :], "perUnitBuild"] * 0.02
    # )  # Mio EUR per unit
    # re_acc.loc[idx["OMFix", :, "horizon", csp_techs, :], "perUnitTotal"] = (
    #     re_acc.loc[idx["Invest", :, "horizon", csp_techs, :], "perUnitBuild"] * 0.015
    # )  # Mio EUR per unit

    invest_vals = re_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    re_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.02 # 2% fixed O&M

    re_acc = re_acc.replace(0, np.nan).dropna(how="all")
    m["Base"].parameter.add(re_acc, "accounting_converterunits")  

def add_geothermal(m):
    geoth_inst_csv = pd.read_csv(f"{path_brownfield}/power-plant-nz-database.csv")
    geoth_vintage = [1950]
    geoth_techs = ["geoth"]
    #geoth_nodes = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]
    geoth_nodes = ["BOP", "NIS", "WTO"]
    geoth_activities = ["Powergen"]

    # techparam 
    geoth_tech = pd.DataFrame(index=pd.MultiIndex.from_product([geoth_techs, geoth_vintage]))
    geoth_tech.loc[idx[:, :], "activityUpperLimit"] = 0.9 # 0 for renewables 
    geoth_tech.loc[idx[:, :], "lifeTime"] = 100  # years, data from: "Financial_Technical assumptions" Ashish 2023 
    m["Base"].parameter.add(geoth_tech, "converter_techparam")
    
    df = geoth_inst_csv
    filtered_df = df[df['Type'] == 'Geothermal']
    #geoth_df = filtered_df.groupby(['Node', 'Year_built', 'Techs']).agg({'Capacity_MW': 'sum', 'Avg_Ann_Gen_GWh': 'sum'}).reset_index()            
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int)

    # capacities
    geoth_cap = (grouped_df
                 .set_index(["Node", "Year_built", "Techs"])
                 .rename(columns={"Capacity_MW": "unitsBuild"})
                 .div(1e3))

    geoth_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([geoth_nodes, yrs_to_calc, geoth_techs]))
    geoth_cap_upper.loc[idx[:, :, :], "noExpansion"] = 1  # boolean
    geoth_cap_upper.loc[idx[:, :, :], "unitsUpperLimit"] = 100 

    geoth_cap = pd.concat([geoth_cap, geoth_cap_upper]).groupby(level=[0,1,2]).sum()
    m["Base"].parameter.add(geoth_cap, "converter_capacityparam")
    
    # activity data
    geoth_inst = geoth_inst_csv[geoth_inst_csv["Type"] == "Geothermal"]

    geoth_feedin = pd.DataFrame(
        data={
            "region": geoth_inst["Node"],  # Cambiar "Node" por "region"
            "year": np.nan,
            "technology": "Geothermal",  # Agregar "technology" desde el inicio
            "type": "upper",  # Agregar "type" desde el inicio
        },
        columns=["region", "year", "technology", "type"],  # Establecer orden de columnas
    )

    geoth_inst_filtered = geoth_inst[
        (geoth_inst["Avg_Ann_Gen_GWh"] != 0) & (geoth_inst["Capacity_MW"].notnull())
    ]

    # Create the DataFrame efficiently
    data = {}
    for i in range(1, 8761):
        data["t{:04d}".format(i)] = np.nan

    # Calculate plant factors using the filtered DataFrame
    for i in range(1, 8761):
        factor = geoth_inst_filtered["Avg_Ann_Gen_GWh"] / (geoth_inst_filtered["Capacity_MW"]*8760 / 100)
        data["t{:04d}".format(i)] = factor

    geoth_feedin_temp = pd.DataFrame(data=data)
    geoth_feedin_grouped = geoth_feedin_temp.groupby(geoth_inst["Node"]).sum()
    geoth_feedin_grouped.reset_index(inplace=True)
    geoth_feedin_grouped["technology"] = "Geothermal"  
    geoth_feedin_grouped["type"] = "upper"  

    dfs = []
    for year in yrs_to_calc:
        geoth_feedin_copy = geoth_feedin_grouped.copy()
        geoth_feedin_copy["year"] = year
        dfs.append(geoth_feedin_copy)

    # Concatenate the separate dataframes into one
    geoth_feedin = pd.concat(dfs, ignore_index=True)
    geoth_feedin = geoth_feedin.rename(columns={"Node": "region"})
    geoth_feedin = geoth_feedin.set_index(["region", "year", "technology", "type"]).sort_index()  # Establecer índice ahora que las columnas existen
    m["Base"].profile.add(geoth_feedin, "converter_activityprofile")
    
    # coefficients    
    geoth_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [geoth_techs, geoth_vintage, geoth_activities, ["Elec"]]
        )
    )

    geoth_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    m["Base"].parameter.add(geoth_coef, "converter_coefficient")

    geoth_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ["horizon"],  geoth_techs, geoth_vintage])
    ).sort_index()


    geoth_acc.loc[idx["Invest",["global"], ["horizon"],  :, :], "perUnitBuild"] = 4970   #, 3610]   # data from: "Financial_Technical assumptions" Ashish 2023   - Mio EUR per unit
    geoth_acc.loc[idx["Invest",["global"], ["horizon"],  :, :], "amorTime"] = 100  # years
    geoth_acc.loc[idx["Invest",["global"], ["horizon"],  :, :], "useAnnuity"] = 1  # binary yes/no
    geoth_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "interest"] = 0.06  # percent/100
    # geoth_acc.loc[idx["OMFix", ["global"], ["horizon"],  :, :], "perUnitTotal"] = (
    #     geoth_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] * 0.016 # data from: "Financial_Technical assumptions" Ashish 2024
    # )  # Mio EUR per unit

    invest_vals = geoth_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],
        names=geoth_acc.index.names
    )
    geoth_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.016


    m["Base"].parameter.add(geoth_acc, "accounting_converterunits")

# conventional

def add_thermal(m):
    the_inst_csv = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))

    the_vintage = [1950, 2030, 2050]                   # vintages (for techparam & coefficients)
    the_techs   = ["BIO", "COAL", "DIE"]
    the_nodes   = [n for n in m["Base"].set.nodesdata]
    the_act     = ["Powergen"]
    the_years   = list(m["Base"].set.yearssel)         # use model years for accounting tables

    # converter_techparam 
    the_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([the_techs, the_vintage])
    )
    the_tech.loc[idx[:, :], "lifeTime"] = 30
    the_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m["Base"].parameter.add(the_tech, "converter_techparam")

    # converter_capacityparam
    the_cap = pd.DataFrame(index=pd.MultiIndex.from_product([the_nodes, range(1992, 2050), the_techs]))
    df = the_inst_csv
    filtered_df = df[(df['Type'] == 'Thermal') & (df['Primary_fuel'].isin(['Biogas', 'Biomass', 'Coal', 'Diesel', 'Waste heat', 'Wood', 'Wood waste']))]
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int)

    the_cap = (
        grouped_df
        .set_index(["Node", "Year_built", "Techs"])
        .rename(columns={"Capacity_MW": "unitsBuild"})
        .div(1e3)
    )

    the_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([the_nodes, yrs_to_calc, the_techs]))
    the_cap_upper.loc[idx[:, [2020], :], "noExpansion"]    = 1
    the_cap_upper.loc[idx[:, :, :],      "unitsUpperLimit"] = 100
    the_cap_full = pd.concat([the_cap, the_cap_upper], axis=1)
    m["Base"].parameter.add(the_cap_full, "converter_capacityparam")

    # converter_coefficient (OK with vintages & activities) 
    the_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([the_techs, the_vintage, the_act, ["Elec"]])
    )
    the_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1
    m["Base"].parameter.add(the_coef, "converter_coefficient")

    # accounting_converterunits (ADD timescope="horizon" and use model years)
    the_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], the_techs, the_years],
            names=["indicator", "regionscope", "timescope", "techs", "years"]
        )
    ).sort_index()

    # broadcast values to all years
    the_acc.loc[idx["Invest", "global", "horizon", "BIO", :],  "perUnitBuild"] = 2600
    the_acc.loc[idx["Invest", "global", "horizon", "COAL", :], "perUnitBuild"] = 1600
    the_acc.loc[idx["Invest", "global", "horizon", "DIE", :],  "perUnitBuild"] = 900

    the_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"] = 1
    the_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"]   = 30
    the_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"]   = 0.06

    invest_vals = the_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"

    )
    the_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.033 # 3.3% fixed O&M


    m["Base"].parameter.add(the_acc, "accounting_converterunits")

    # accounting_converteractivity (ADD timescope="horizon" and use model years)
    the_emission = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["CO2_emission", "FuelCost"], ["global"], ["horizon"], the_techs, the_years, the_act],
            names=["indicator","regionscope","timescope","techs","years","activities"]
        )
    ).sort_index()

    # Broadcast per-activity values to all years 
    # CO2 factors (kt_CO2 per unit of activity)
    the_emission.loc[idx["CO2_emission","global","horizon","BIO", :, "Powergen"], "perActivity"] = 0.0
    the_emission.loc[idx["CO2_emission","global","horizon","COAL",:, "Powergen"], "perActivity"] = 0.3406 / 0.010   # = 34.06
    the_emission.loc[idx["CO2_emission","global","horizon","DIE", :, "Powergen"], "perActivity"] = 0.2668 / 0.262   # ≈ 1.018

    # Fuel costs (m EUR per unit of activity)
    the_emission.loc[idx["FuelCost","global","horizon","BIO", :, "Powergen"], "perActivity"] = 0.03   / 0.0001  # = 300
    the_emission.loc[idx["FuelCost","global","horizon","COAL",:, "Powergen"], "perActivity"] = 0.15   / 0.010   # = 15
    the_emission.loc[idx["FuelCost","global","horizon","DIE", :, "Powergen"], "perActivity"] = 0.58   / 0.262   # ≈ 2.214

    m["Base"].parameter.add(the_emission, "accounting_converteractivity")

def add_gas_turbines(m):    
    gt_inst_csv = pd.read_csv(Path(path_brownfield).joinpath("power-plant-nz-database.csv"))

    gt_vintage = [1950,2030]

    #gt_techs = ["GT", "CCGT", "OCGT", "GT_H2", "CCGT_H2"]                     
    gt_techs = ["GT", "CCGT", "OCGT"]                       
    gt_nodes = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]
    gt_activities = ["Powergen"]

    gt_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([gt_techs, gt_vintage])
    )
    gt_tech.loc[idx[:, :], "lifeTime"] = 30
    gt_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m["Base"].parameter.add(gt_tech, "converter_techparam")

    df = gt_inst_csv
    filtered_df = df[df['Primary_fuel'] == 'Natural gas']
    grouped_df = filtered_df.groupby(['Node', 'Year_built', 'Techs'])['Capacity_MW'].sum().reset_index()
    grouped_df['Year_built'] = grouped_df['Year_built'].astype(int) 

    gt_cap = (grouped_df
              .set_index(["Node", "Year_built", "Techs"])
              .rename(columns={"Capacity_MW": "unitsBuild"})
              .div(1e3))

    gt_cap_upper = pd.DataFrame(index=pd.MultiIndex.from_product([gt_nodes, yrs_to_calc, gt_techs]))
    gt_cap_upper.loc[idx[:, :, gt_techs], "unitsUpperLimit"] = 100  # GW_el
    if 2020 in yrs_sel:
        gt_cap_upper.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean  
    gt_cap_full = pd.concat([gt_cap, gt_cap_upper], axis=1)
    m["Base"].parameter.add(gt_cap_full, "converter_capacityparam")

    # coefficients  
    gt_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [gt_techs, gt_vintage, gt_activities, ["Elec"]]
        )
    )

    gt_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    m["Base"].parameter.add(gt_coef, "converter_coefficient")
    
    
    # gt_coef = pd.DataFrame(
    #     index=pd.MultiIndex.from_product(
    #         [gt_techs, gt_vintage, gt_activities, ["CH4", "H2", "Elec"]]
    #     )
    # )

    # gt_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = 1 # GW_el
    # gt_coef.loc[idx["GT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.41, 0.43]), 3) # GW_el #-1/efficiency=coef, rounded to 3, line above output of 1 GW
    # gt_coef.loc[idx["CCGT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.58, 0.61]), 3) # GW_el
    # gt_coef.loc[idx["OCGT", :, :, "CH4"], "coefficient"] = np.round(-1 / np.array([0.4695, 0.4695]), 3) # GW_el - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    # gt_coef.loc[idx["GT_H2", :, :, "H2"], "coefficient"] = np.round(-1 / np.array([.41, 0.43]), 3) # GW_elel
    # gt_coef.loc[idx["CCGT_H2", :, :, "H2"], "coefficient"] = np.round(-1 / np.array([0.58, 0.61]), 3) # GW_el
    # m["Base"].parameter.add(gt_coef, "converter_coefficient")

    gt_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"],  gt_techs, gt_vintage]
        )
    )
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  "GT", :], "perUnitBuild"] = [900, 830]  # million EUR / unit
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  "CCGT", :], "perUnitBuild"] = [775, 775]  # million EUR / unit
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  "OCGT", :], "perUnitBuild"] = [475, 475]  # million EUR / unit  - data: LUT Breyer "Financial_Technical assumptions-newversion.docx"
    #gt_acc.loc[idx["Invest", ["global"], ["horizon"],  "GT_H2", :], "perUnitBuild"] = [900, 830]  # million EUR / unit
    #gt_acc.loc[idx["Invest", ["global"], ["horizon"],  "CCGT_H2", :], "perUnitBuild"] = [600, 560]  # million EUR / unit
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "useAnnuity"] = 1  # binary yes/no
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "amorTime"] = 30  # years
    gt_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "interest"] = 0.06  # percent/100
    # gt_acc.loc[idx["OMFix", ["global"], ["horizon"],  :, :], "perUnitTotal"] = (
    #     gt_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] * 0.0193
    # )  # Mio EUR per unit

    invest_vals = gt_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    gt_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.0193 # 
    
    
    m["Base"].parameter.add(gt_acc, "accounting_converterunits")

    # # Emit carbon from combustion
    gt_emission = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission", "FuelCost"], ["global"], ["horizon"],  ["GT", "CCGT","OCGT"], gt_vintage, gt_activities]))
    gt_emission.loc[idx["CO2_emission", ["global"], ["horizon"],  "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * 0.2016 # kt_CO2
    gt_emission.loc[idx["CO2_emission", ["global"], ["horizon"],  "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3) * 0.2016 # kt_CO2
    gt_emission.loc[idx["CO2_emission", ["global"], ["horizon"],  "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    gt_emission.loc[idx["FuelCost", ["global"], ["horizon"],  "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * np.array([0.21, 0.31]) # mEUR
    gt_emission.loc[idx["FuelCost", ["global"], ["horizon"],  "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58, 0.61]), 3) * np.array([0.21, 0.31]) # mEUR
    gt_emission.loc[idx["FuelCost", ["global"], ["horizon"],  "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.41, 0.43]), 3) * np.array([0.21, 0.31]) # mEUR
    # Emit carbon from combustion
    # gt_emission = pd.DataFrame(index=pd.MultiIndex.from_product([["CO2_emission", "FuelCost"], ["global"], ["GT", "CCGT","OCGT"], gt_vintage, gt_activities]))
    # gt_emission.loc[idx["CO2_emission", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41]), 3) * 0.2016 # kt_CO2
    # gt_emission.loc[idx["CO2_emission", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58]), 3) * 0.2016 # kt_CO2
    # gt_emission.loc[idx["CO2_emission", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695]), 3) * 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    # gt_emission.loc[idx["FuelCost", :, "GT", :, :], "perActivity"] = np.round(1 / np.array([0.41]), 3)* 0.2016 # kt_CO2
    # gt_emission.loc[idx["FuelCost", :, "CCGT", :, :], "perActivity"] = np.round(1 / np.array([0.58]), 3)* 0.2016 # kt_CO2
    # gt_emission.loc[idx["FuelCost", :, "OCGT", :, :], "perActivity"] = np.round(1 / np.array([0.4695]), 3)* 0.2016 # kt_CO2 - data: REMix Tutorial 202, number for 2020 used for 2000 and 2030
    
    m["Base"].parameter.add(gt_emission, "accounting_converteractivity")

# hydrogen

def add_electrolyser(m):
    """
    Electrolyser values from supplementary material https://doi.org/10.1016/j.energy.2025.134878 (page 9)
    capex €/kW_el: 446 (early), 291 (late)
    fixed opex: 3.5% of capex
    lifetime: 30 years
    efficiency: 70%
    """
    years = yrs_to_calc
    eltr_vintage = years
    eltr_nodes = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]
    eltr_techs = ["Electrolyser"]

    # technology
    eltr_tech = pd.DataFrame(index=pd.MultiIndex.from_product([eltr_techs, eltr_vintage]))
    eltr_tech["lifeTime"] = 30
    eltr_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(eltr_tech, "converter_techparam")

    # capacities
    eltr_caps = pd.DataFrame(index=pd.MultiIndex.from_product([eltr_nodes, years, eltr_techs]))
    eltr_caps["unitsUpperLimit"] = 100
    m["Base"].parameter.add(eltr_caps, "converter_capacityparam")

    # coefficients: 70% efficiency
    eltr_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([eltr_techs, eltr_vintage, ["Electrolysis"], ["Elec", "H2"]])
    )
    eltr_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -1.0
    eltr_coef.loc[idx[:, :, :, "H2"], "coefficient"] = 0.70
    m["Base"].parameter.add(eltr_coef, "converter_coefficient")

    # accounting
    eltr_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ["horizon"], eltr_techs, eltr_vintage])
    ).sort_index()

    # capex: map early to 2020, late to 2050 (and reuse late for any later)
    capex_by_year = {}
    for y in eltr_vintage:
        capex_by_year[y] = 446 if y < 2050 else 291

    eltr_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] = [capex_by_year[y] for y in eltr_vintage]
    eltr_acc.loc[idx["Invest", :, "horizon", :, :], "amorTime"] = 30
    eltr_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1
    eltr_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06

    invest_vals = eltr_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples([("OMFix", *i[1:]) for i in invest_vals.index])
    eltr_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.035  # 3.5% fixed O&M

    m["Base"].parameter.add(eltr_acc, "accounting_converterunits")

def add_H2_CCGT(m):
    H2_CCGT_vintage = [2030, 2035, 2040, 2045, 2050]
    H2_CCGT_nodes = [n for n in m["Base"].set.nodesdata]
    # technology
    H2_CCGT_tech = pd.DataFrame(index=pd.MultiIndex.from_product([["H2_CCGT"], H2_CCGT_vintage]))
    H2_CCGT_tech.loc[idx[:, :], ["lifeTime"]] = 35  # years
    H2_CCGT_tech["activityUpperLimit"] = 1  # availability of technology
    m["Base"].parameter.add(H2_CCGT_tech, "converter_techparam")

    # capacities
    H2_CCGT_caps = pd.DataFrame(index=pd.MultiIndex.from_product([H2_CCGT_nodes, yrs_to_calc, ["H2_CCGT"]]))
    H2_CCGT_caps["unitsUpperLimit"] = 100  # GW_el
    m["Base"].parameter.add(H2_CCGT_caps, "converter_capacityparam")

    # coefficients
    H2_CCGT_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["H2_CCGT"],
                H2_CCGT_vintage,
                ["Powergen"],
                ["Elec", "H2"],
            ]
        )
    )
    H2_CCGT_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1     
    H2_CCGT_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = [0.58, 0.59, 0.60, 0.60,  0.60] # C. Habib
    m["Base"].parameter.add(H2_CCGT_coef, "converter_coefficient")

    # accounting
    H2_CCGT_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"],["horizon"], ["H2_CCGT"], H2_CCGT_vintage])
        ).sort_index()

    H2_CCGT_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "perUnitBuild"] = 100 #853   #  
    H2_CCGT_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "amorTime"] = 35  # years
    H2_CCGT_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "useAnnuity"] = 1  # binary yes/no
    H2_CCGT_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "interest"] = 0.06  # percent/100
    # H2_CCGT_acc.loc[idx["OMFix", ["global"],["horizon"], :, :], "perUnitTotal"] = (
    #     H2_CCGT_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "perUnitBuild"] * 0.025
    #     )  # Mio EUR per unit
    invest_vals = H2_CCGT_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  
    )
    H2_CCGT_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.025 # 2.5% fixed O&M


    m["Base"].parameter.add(H2_CCGT_acc, "accounting_converterunits")

def add_H2_FC(m):
    H2_FC_vintage = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    H2_FC_nodes = [n for n in m["Base"].set.nodesdata]
    # technology
    H2_FC_tech = pd.DataFrame(index=pd.MultiIndex.from_product([["H2_FC"], H2_FC_vintage]))
    H2_FC_tech.loc[idx[:, :], ["lifeTime"]] = 35  # years
    H2_FC_tech["activityUpperLimit"] = 1  # availability of technology
    m["Base"].parameter.add(H2_FC_tech, "converter_techparam")

    # capacities
    H2_FC_caps = pd.DataFrame(index=pd.MultiIndex.from_product([H2_FC_nodes, yrs_to_calc, ["H2_FC"]]))
    H2_FC_caps["unitsUpperLimit"] = 100  # GW_el
    m["Base"].parameter.add(H2_FC_caps, "converter_capacityparam")

    # coefficients
    H2_FC_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [
                ["H2_FC"],
                H2_FC_vintage,
                ["Powergen"],
                ["Elec", "H2"],
            ]
        )
    )
    H2_FC_coef.loc[idx[:, :, :, "H2"], "coefficient"] = -1     
    H2_FC_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = [0.579, 0.6134, 0.6383, 0.6477, 0.6514, 0.6686, 0.6737] # C. Habib
    m["Base"].parameter.add(H2_FC_coef, "converter_coefficient")

    # accounting
    H2_FC_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product([["Invest", "OMFix"], ["global"], ["horizon"],  ["H2_FC"], H2_FC_vintage])
        ).sort_index()

    H2_FC_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] = [2980.992, 1468.139, 773.195, 733.604, 694.012, 654.421, 614.830]   #  
    H2_FC_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "amorTime"] = 35  # years
    H2_FC_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "useAnnuity"] = 1  # binary yes/no
    H2_FC_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "interest"] = 0.06  # percent/100
    # H2_FC_acc.loc[idx["OMFix", ["global"], ["horizon"],  :, :], "perUnitTotal"] = (
    #     H2_FC_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] * 0.05
    #     )  # Mio EUR per unit
    invest_vals = H2_FC_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    H2_FC_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.05 # 5% fixed O&M
    
    m["Base"].parameter.add(H2_FC_acc, "accounting_converterunits")

def add_dac(m):
    idx = pd.IndexSlice

    dac_vintage = [2020, 2030, 2040, 2050]
    dac_techs   = ["DAC"]
    dac_nodes   = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]

    # tech param: techs, years
    dac_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [dac_techs, dac_vintage],
            names=["techs", "years"],
        )
    )
    dac_tech["lifeTime"]          = 20      # years
    dac_tech["activityUpperLimit"] = 1      # availability
    m["Base"].parameter.add(dac_tech, "converter_techparam")

    # capacities: nodes, years, techs
    dac_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [dac_nodes, [2020, 2025, 2030, 2035, 2040, 2045, 2050], dac_techs],
            names=["nodes", "years", "techs"],
        )
    )
    dac_caps["unitsUpperLimit"] = 100  # GW_el
    m["Base"].parameter.add(dac_caps, "converter_capacityparam")

    # coefficients: techs, years, activities, commodities
    dac_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [dac_techs, dac_vintage, ["Capture"], ["Elec", "CO2_feed"]],
            names=["techs", "years", "activities", "commodities"],
        )
    )
    # Output: CO2_feed = +1 per unit activity
    dac_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = 1.0

    # Input: Elec; one value per vintage
    elec_coeff = [-1.535, -1.458, -1.385, -1.316]  # GWh_el per ktCO2
    dac_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = elec_coeff

    m["Base"].parameter.add(dac_coef, "converter_coefficient")

    # accounting (converter units): indicator, regionscope, timescope, techs, years
    dac_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], dac_techs, dac_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # CAPEX: same for all vintages, broadcast as scalar
    capex_meur_per_ktco2h = 815 * 8.76  # EUR/tCO2*a → M€/ktCO2*h
    dac_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] = capex_meur_per_ktco2h
    dac_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"]     = 20
    dac_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"]   = 1
    dac_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"]     = 0.06

    # OMFix = 4% of capex, aligned by index
    invest = dac_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    dac_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest * 0.04

    m["Base"].parameter.add(dac_acc, "accounting_converterunits")

    # remove carbon from CO2_emission indicator for DAC
    dac_activity = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["CO2_emission"], ["global"], ["horizon"], ["DAC"], dac_vintage, ["Capture"]],
            names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
        )
    )
    dac_activity["perActivity"] = -1.0
    m["Base"].parameter.add(dac_activity, "accounting_converteractivity")

    print("DAC added successfully.")


def add_methanizer(m):
    """
    Methanizer: H2 + CO2_feed -> CH4 (synthetic methane).
    """
    methanizer_vintage = [2020, 2030, 2040, 2050]
    methanizer_techs   = ["Methanizer"]
    methanizer_nodes   = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]

    #tech param
    methanizer_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanizer_techs, methanizer_vintage],
            names=["techs", "years"],
        )
    )
    methanizer_tech["lifeTime"] = [30, 30, 30, 30]
    methanizer_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(methanizer_tech, "converter_techparam")

    # caps
    methanizer_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanizer_nodes,
             [2020, 2025, 2030, 2035, 2040, 2045, 2050],
             methanizer_techs],
            names=["nodes", "years", "techs"],
        )
    )
    methanizer_caps["unitsUpperLimit"] = 100  # GW_el-equiv
    m["Base"].parameter.add(methanizer_caps, "converter_capacityparam")

    # Conversion coefficients: techs, years, activities, commodities
    methanizer_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanizer_techs, methanizer_vintage,
             methanizer_techs,  # activity name = tech name in your original
             ["Elec", "H2", "CH4", "CO2_feed"]],
            names=["techs", "years", "activities", "commodities"],
        )
    )

    methanizer_coef.loc[idx[:, :, :, "CH4"],      "coefficient"] = 1.0
    methanizer_coef.loc[idx[:, :, :, "H2"],       "coefficient"] = -1.284
    methanizer_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.2016
    methanizer_coef.loc[idx[:, :, :, "Elec"],     "coefficient"] = -0.006

    m["Base"].parameter.add(methanizer_coef, "converter_coefficient")

    # accounting (converter units): indicator, regionscope, timescope, techs, years
    methanizer_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"],
             methanizer_techs, methanizer_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # CAPEX (Mio EUR per unit) 
    methanizer_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] = [
        558, 309, 251, 211
    ]
    methanizer_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"] = [20, 25, 30, 30]
    methanizer_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"] = 1
    methanizer_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"] = 0.06

    # OMFix = 4% of capex
    invest_vals = methanizer_acc.loc[
        idx["Invest", "global", "horizon", :, :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix",) + i[1:] for i in invest_vals.index],
        names=methanizer_acc.index.names,
    )
    methanizer_acc.loc[
        idx["OMFix", "global", "horizon", :, :], "perUnitTotal"
    ] = invest_vals * 0.04

    m["Base"].parameter.add(methanizer_acc, "accounting_converterunits")

    print("Methanizer added successfully.")


def add_methanolsyn(m):
    """
    Methanol synthesis: H2 + CO2_feed -> CH3OH.
    """
    methanolsyn_vintage = [2020]
    methanolsyn_techs   = ["MethanolSyn"]
    methanol_nodes      = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]

    # tech param
    methanolsyn_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanolsyn_techs, methanolsyn_vintage],
            names=["techs", "years"],
        )
    )
    methanolsyn_tech["lifeTime"] = 30
    methanolsyn_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(methanolsyn_tech, "converter_techparam")

    # caps
    methanolsyn_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanol_nodes,
             [2020, 2025, 2030, 2035, 2040, 2045, 2050],
             methanolsyn_techs],
            names=["nodes", "years", "techs"],
        )
    )
    methanolsyn_caps["unitsUpperLimit"] = 100  # GW_el-equiv
    m["Base"].parameter.add(methanolsyn_caps, "converter_capacityparam")

    # conversion coefficients: techs, years, activities, commodities
    methanolsyn_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [methanolsyn_techs, methanolsyn_vintage,
             ["Synthesis"], ["CH3OH", "H2", "CO2_feed", "Elec"]],
            names=["techs", "years", "activities", "commodities"],
        )
    )

    methanolsyn_coef.loc[idx[:, :, :, "CH3OH"],    "coefficient"] = 1.0
    methanolsyn_coef.loc[idx[:, :, :, "H2"],       "coefficient"] = -1.25
    methanolsyn_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.219
    methanolsyn_coef.loc[idx[:, :, :, "Elec"],     "coefficient"] = -0.1

    m["Base"].parameter.add(methanolsyn_coef, "converter_coefficient")

    # accounting (converter units): indicator, regionscope, timescope, techs, years
    methanolsyn_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"],
             methanolsyn_techs, methanolsyn_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()


    # CAPEX (Mio EUR per unit) 
    methanolsyn_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] = 835
    methanolsyn_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"] = 30
    methanolsyn_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"] = 1
    methanolsyn_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"] = 0.06

    # OMFix = 4% of capex
    invest_vals = methanolsyn_acc.loc[
        idx["Invest", "global", "horizon", :, :], "perUnitBuild"
    ]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix",) + i[1:] for i in invest_vals.index],
        names=methanolsyn_acc.index.names,
    )
    methanolsyn_acc.loc[
        idx["OMFix", "global", "horizon", :, :], "perUnitTotal"
    ] = invest_vals * 0.04

    m["Base"].parameter.add(methanolsyn_acc, "accounting_converterunits")

    print("Methanol synthesis added successfully.")


def add_ftropsch_syn(m):
    """
    Fischer-Tropsch synthesis (generic liquid e-fuel 'REfuel').
    - Output: REfuel
    - Inputs: H2, CO2_feed, Elec
    """
    ftropsch_vintage = [2020, 2040, 2050]
    ftropsch_techs   = ["FTropschSyn"]
    ftropsch_nodes   = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]

    # tech param
    ftropsch_syn_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [ftropsch_techs, ftropsch_vintage],
            names=["techs", "years"],
        )
    )
    ftropsch_syn_tech["lifeTime"] = [30, 30, 30]
    ftropsch_syn_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(ftropsch_syn_tech, "converter_techparam")

    # caps
    ftropsch_syn_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [ftropsch_nodes, [2020, 2025, 2030, 2035, 2040, 2045, 2050], ftropsch_techs],
            names=["nodes", "years", "techs"],
        )
    )
    ftropsch_syn_caps["unitsUpperLimit"] = 100  # GW_el-equivalent
    m["Base"].parameter.add(ftropsch_syn_caps, "converter_capacityparam")

    # conversion coefficients: techs, years, activities, commodities
    ftropsch_syn_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [ftropsch_techs, ftropsch_vintage, ["Synthesis"], ["REfuel", "H2", "CO2_feed", "Elec"]], #, "H2O"
            names=["techs", "years", "activities", "commodities"],
        )
    )

    # Output: generic liquid e-fuel
    ftropsch_syn_coef.loc[idx[:, :, :, "REfuel"], "coefficient"] = 1.0

    # Inputs and by‑product  - data from Andi SynLink? 
    ftropsch_syn_coef.loc[idx[:, :, :, "H2"],      "coefficient"] = -1.52
    ftropsch_syn_coef.loc[idx[:, :, :, "CO2_feed"], "coefficient"] = -0.35
    ftropsch_syn_coef.loc[idx[:, :, :, "Elec"],    "coefficient"] = -0.11
    #ftropsch_syn_coef.loc[idx[:, :, :, "H2O"],     "coefficient"] = 1.0

    # Drop any pure-zero rows, if any
    ftropsch_syn_coef = ftropsch_syn_coef.replace(0, np.nan).dropna(how="all")
    m["Base"].parameter.add(ftropsch_syn_coef, "converter_coefficient")

    # accounting (converter units): indicator, regionscope, timescope, techs, years
    ftropsch_syn_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], ftropsch_techs, ftropsch_vintage],
            names=["indicator", "regionscope", "timescope", "techs", "years"],
        )
    ).sort_index()

    # CAPEX per unit (Mio EUR per unit)
    ftropsch_syn_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] = [1017, 915, 915]
    ftropsch_syn_acc.loc[idx["Invest", "global", "horizon", :, :], "amorTime"] = [30, 30, 30]
    ftropsch_syn_acc.loc[idx["Invest", "global", "horizon", :, :], "useAnnuity"] = 1
    ftropsch_syn_acc.loc[idx["Invest", "global", "horizon", :, :], "interest"] = 0.06

    # OMFix = 4% of capex
    invest_vals = ftropsch_syn_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"]
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix",) + i[1:] for i in invest_vals.index],
        names=ftropsch_syn_acc.index.names,
    )
    ftropsch_syn_acc.loc[
        idx["OMFix", "global", "horizon", :, :], "perUnitTotal"
    ] = invest_vals * 0.04

    m["Base"].parameter.add(ftropsch_syn_acc, "accounting_converterunits")
    print("Fischer-Tropsch Synthesis added successfully.")


# storage
    
def add_lithium_batteries(m):
    battery_vintage = [2020, 2030, 2040, 2050] 
    battery_techs = ["Battery"]
    battery_nodes = [n for n in m["Base"].set.nodesdata if not n.startswith("LNG")]

    # technology
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage])
    )
    conv_tech.loc[idx[:, :], "lifeTime"] = 20
    conv_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # capacities
    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [battery_nodes, yrs_sel, battery_techs]
        )
    )
    conv_cap.loc[idx[:, :, battery_techs], "unitsUpperLimit"] = 50  # GW_el
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # conversion coefficients: techs, years, activities, commodities
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [battery_techs,battery_vintage, ["Charge", "Discharge"], ["Elec", "Elec_battery"]]
        )
    )
    conv_coef.loc[idx[:, :, "Charge", "Elec"], "coefficient"] = -1 # GW_el
    conv_coef.loc[idx[:, :, "Charge", "Elec_battery"], "coefficient"] = 0.975 # GW_el
    conv_coef.loc[idx[:, :, "Discharge", "Elec"], "coefficient"] = 1 # GW_el
    conv_coef.loc[idx[:, :, "Discharge", "Elec_battery"], "coefficient"] = -1.025 # GW_el
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # accounting
    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"],["horizon"], battery_techs, battery_vintage]
        )
    ) 
    conv_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "perUnitBuild"] = [117, 55, 37, 30]  # million EUR / unit
    conv_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "useAnnuity"] = 1  # binary yes/no
    conv_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "amorTime"] = 20  # years
    conv_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "interest"] = 0.06  # percent/100
    # conv_acc.loc[idx["OMFix", ["global"],["horizon"], :, :], "perUnitTotal"] = (
    #     conv_acc.loc[idx["Invest", ["global"],["horizon"], :, :], "perUnitBuild"] * 0.014
    # )  # Mio EUR per unit
    invest_vals = conv_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    conv_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.014 # 1.4% fixed O&M

    m["Base"].parameter.add(conv_acc, "accounting_converterunits")

    # storage tech
    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage])
    )
    stor_tech.loc[idx[:, :], "lifeTime"] = 20
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1

    m["Base"].parameter.add(stor_tech, "storage_techparam")
    stor_tech


    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_techs, battery_vintage, ["Elec_battery"]])
    )
    stor_size.loc[idx["Battery", :, "Elec_battery"], "size"] = 4  # GWh_ch / unit
    m["Base"].parameter.add(stor_size, "storage_sizeparam")


    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([battery_nodes,  yrs_sel, battery_techs])
    )
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 30  # units
    if 2020 in yrs_sel:
        stor_res.loc[idx[:, ["2020"], :], "noExpansion"] = 1  # boolean
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], battery_techs, battery_vintage]
        )
    )
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] = [i * 4 for i in [234, 110, 76, 61]]  # million EUR / unit
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "amorTime"] = 20  # years
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06  # percent/100
    # stor_acc.loc[idx["OMFix", :, "horizon", :, :], "perUnitTotal"] = (
    #     stor_acc.loc[idx["Insvest", :, "horizon", :, :], "perUnitBuild"] * 0.014
    # )  # Mio EUR per unit
    invest_vals = stor_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    stor_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.014 # 1.4% fixed O&M

    m["Base"].parameter.add(stor_acc, "accounting_storageunits")

def add_h2_storage(m):

    # converter is the compressor 
    # storage is the gas tank
    h2_stor_vintage = [2020]#, 2030, 2040, 2050]
    yrs_h2=[2020, 2035, 2050]#, 2025, 2030, 2035, 2040, 2045, 2050]
    h2_stor_techs = ["H2_storage"]
    h2_stor_nodes = [n for n in m["Base"].set.nodesdata]

    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2])
    )
    conv_tech.loc[idx[:, :], "lifeTime"] = 40
    conv_tech.loc[idx[:, :], "activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_nodes, yrs_h2, h2_stor_techs]
        )
    )
    conv_cap.loc[idx[:, :, h2_stor_techs], "unitsUpperLimit"] = 50 # GW_el
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")


    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [h2_stor_techs,  yrs_h2, ["Charge", "Discharge"], ["H2", "H2_stored","Elec"]]
        )
    )
    # electrolyzer_coef.loc[idx[:, :, :, "Elec"], "coefficient"] = -1
    conv_coef.loc[idx[:, :, "Charge", "Elec"], "coefficient"] =  [ -0.043346085, -0.035470537,  -0.031529343] #[ -0.043346085]#,	-0.043346085, -0.035470537,	-0.035470537, -0.035470537, -0.031529343, -0.031529343]
    # 2020, 2025, 030, 2035, 2040, 2045, 2040...
    conv_coef.loc[idx[:, :, "Charge", "H2"], "coefficient"] = -0.15  # GW_h2
    conv_coef.loc[idx[:, :, "Charge", "H2_stored"], "coefficient"] = 0.15  * 0.89 # GW_h2
    conv_coef.loc[idx[:, :, "Discharge", "H2"], "coefficient"] = 0.15  # GW_h2
    conv_coef.loc[idx[:, :, "Discharge", "H2_stored"], "coefficient"] = -0.15 * 1.11 # GW_h2
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], h2_stor_techs, yrs_h2]
        )
    )
    #
    conv_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] = 5.14  # million EUR / unit
    conv_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "useAnnuity"] = 1  # binary yes/no
    conv_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "amorTime"] = 40  # years
    conv_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "interest"] = 0.06  # percent/100
    # conv_acc.loc[idx["OMFix", ["global"], ["horizon"],  :, :], "perUnitTotal"] = (
    #     conv_acc.loc[idx["Invest", ["global"], ["horizon"],  :, :], "perUnitBuild"] * 0.04
    # )  # Mio EUR per unit
    invest_vals = conv_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    conv_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.04 # 4% fixed O&M

    m["Base"].parameter.add(conv_acc, "accounting_converterunits")


    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2])
    )
    stor_tech.loc[idx[:, :], "lifeTime"] = 40
    stor_tech.loc[idx[:, :], "levelUpperLimit"] = 1

    m["Base"].parameter.add(stor_tech, "storage_techparam")
    stor_tech


    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_techs, yrs_h2, ["H2_stored"]])
    )
    stor_size.loc[idx["H2_storage", :, "H2_stored"], "size"] = 6.9993  # GWh / unit
    m["Base"].parameter.add(stor_size, "storage_sizeparam")


    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([h2_stor_nodes, yrs_h2, h2_stor_techs])
    )
    stor_res.loc[idx[:, :, :], "unitsUpperLimit"] = 50  # units 
    if 2020 in yrs_sel:
        stor_res.loc[idx[:, [2020], :], "noExpansion"] = 1  # boolean
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], h2_stor_techs, yrs_h2]
        )
    )
   
     #CAPEX (628.14 EUR / kgH2 * 210000 kgH2/unit ) /1M = 131.91 M EUR / unit
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] = 131.91  # million EUR / unit
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "useAnnuity"] = 1  # binary yes/no
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "amorTime"] = 40  # years 
    stor_acc.loc[idx["Invest", :, "horizon", :, :], "interest"] = 0.06  # percent/100
    # stor_acc.loc[idx["OMFix", :, "horizon", :, :], "perUnitTotal"] = (
    #     stor_acc.loc[idx["Invest", :, "horizon", :, :], "perUnitBuild"] * 0.03
    # )  # Mio EUR per unit
    invest_vals = stor_acc.loc[idx["Invest", "global", "horizon", :, :], "perUnitBuild"] 
    invest_vals.index = pd.MultiIndex.from_tuples(
        [("OMFix", *i[1:]) for i in invest_vals.index],  # replace "Invest" by "OMFix"
        #names=conv_acc.index.names
    )
    stor_acc.loc[idx["OMFix", "global", "horizon", :, :], "perUnitTotal"] = invest_vals * 0.03 # 3% fixed O&M

    m["Base"].parameter.add(stor_acc, "accounting_storageunits")

# co2 constrains

def add_emission_limit(m):
    # Add net zero restriction for 2050
    accounting_emissionLimit = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], [2050], ["CO2Emission"]])
    )
    accounting_emissionLimit["useUpper"] = 1  # minimization of system costs
    accounting_emissionLimit["upperValue"] = 0  # minimization of system costs
    m["Base"].parameter.add(accounting_emissionLimit, "accounting_indicatorbounds")

def add_emission_budget(m):
    # Add cumulativ emission budget for all years
    accounting_emissionBudget = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], ["horizon"], ["CO2Emission"]])
    )
    accounting_emissionBudget["integral"] = 1
    accounting_emissionBudget["endyear"] = 25 # length of the last year to run (2050)
    accounting_emissionBudget["useUpper"] = 1
    accounting_emissionBudget["upperValue"] = 45000  # 45 Gt CO2
    m["Base"].parameter.add(accounting_emissionBudget, "accounting_indicatorbounds")

def add_emission_slack(m, year=2050, slack_cost_eur_per_tco2=1_000_000):
    """
    Adds a high-cost slack price for CO2 emission violations.
    Follows the tutorial 202 structure: (SystemCostIndicator, indicator, regionscope, timescope).
    """
    # Convert EUR/tCO2 to M€/ktCO2
    cost_meur_per_kt = (slack_cost_eur_per_tco2 * 1000) / 1e6

    # We keep timescope = "horizon" (aggregate over all optimisation years),
    # but weight only the Slack_CO2Emission indicator.
    idx_slack = pd.MultiIndex.from_product(
        [
            ["SystemCost"],          # aggregate indicator
            ["Slack_CO2Emission"],   # component indicator
            ["global"],
            ["horizon"],
        ],
        names=["SystemCostIndicator", "indicator", "regionscope", "timescope"],
    )

    perind = pd.DataFrame(index=idx_slack)
    perind["perIndicator"] = cost_meur_per_kt

    m["Base"].parameter.add(perind, "accounting_perindicator")


# others
    
def add_network(m):
    #First we need to set up the link connections in the data by defining the starting and ending node of each link
    link_names = ["NIS__AKL",
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
                  "WTO__HBY"]
    nodes_data = set(m["Base"].set.nodesdata)
    
    link_connections = pd.DataFrame(index=pd.MultiIndex.from_product([link_names, nodes_data]))
    link_connections.loc[idx["NIS__AKL", "NIS"], ["start"]] = 1
    link_connections.loc[idx["NIS__AKL", "AKL"], ["end"]] = 1
    link_connections.loc[idx["AKL__WTO","AKL"], ["start"]] = 1
    link_connections.loc[idx["AKL__WTO", "WTO"], ["end"]] = 1
    link_connections.loc[idx["WTO__BOP", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__BOP", "BOP"], ["end"]] = 1
    link_connections.loc[idx["WTO__CEN", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__CEN", "CEN"], ["end"]] = 1
    link_connections.loc[idx["CEN__HBY", "CEN"], ["start"]] = 1
    link_connections.loc[idx["CEN__HBY", "HBY"], ["end"]] = 1
    link_connections.loc[idx["TRN__CEN", "TRN"], ["start"]] = 1
    link_connections.loc[idx["TRN__CEN", "CEN"], ["end"]] = 1
    link_connections.loc[idx["CEN__WEL", "CEN"], ["start"]] = 1
    link_connections.loc[idx["CEN__WEL", "WEL"], ["end"]] = 1
    link_connections.loc[idx["WEL__CAN", "WEL"], ["start"]] = 1
    link_connections.loc[idx["WEL__CAN", "CAN"], ["end"]] = 1
    link_connections.loc[idx["NEL__CAN", "NEL"], ["start"]] = 1
    link_connections.loc[idx["NEL__CAN", "CAN"], ["end"]] = 1
    link_connections.loc[idx["CAN__OTG", "CAN"], ["start"]] = 1
    link_connections.loc[idx["CAN__OTG", "OTG"], ["end"]] = 1
    link_connections.loc[idx["AKL__TRN", "AKL"], ["start"]] = 1
    link_connections.loc[idx["AKL__TRN", "TRN"], ["end"]] = 1
    link_connections.loc[idx["WTO__HBY", "WTO"], ["start"]] = 1
    link_connections.loc[idx["WTO__HBY", "HBY"], ["end"]] = 1
    link_connections = link_connections.dropna(how="all").fillna(0)

    m["Base"].parameter.add(link_connections, "transfer_linkstartend")
    link_connections

    #Next we define the lengths for each corridor. We can use different distance types.
    # "transport_distanceParam"
    link_types = ["land", "sea"]
 
    link_lengths = pd.DataFrame(index=pd.MultiIndex.from_product([link_names, link_types]))
    link_lengths.loc[idx["NIS__AKL", "land"], ["length"]] = 149.8
    link_lengths.loc[idx["AKL__WTO", "land"], ["length"]] = 136.0
    link_lengths.loc[idx["WTO__BOP", "land"], ["length"]] = 76.2
    link_lengths.loc[idx["WTO__CEN", "land"], ["length"]] = 154.5
    link_lengths.loc[idx["CEN__HBY", "land"], ["length"]] = 96.5
    link_lengths.loc[idx["TRN__CEN", "land"], ["length"]] = 115.1
    link_lengths.loc[idx["CEN__WEL", "land"], ["length"]] = 83.1
    link_lengths.loc[idx["WEL__CAN", "land"], ["length"]] = 316.4
    link_lengths.loc[idx["NEL__CAN", "land"], ["length"]] = 203.9
    link_lengths.loc[idx["CAN__OTG", "land"], ["length"]] = 179.5
    link_lengths.loc[idx["AKL__TRN", "land"], ["length"]] = 200.1
    link_lengths.loc[idx["WTO__HBY", "land"], ["length"]] = 114.6
    # link_lengths = link_lengths.dropna()
    m["Base"].parameter.add(link_lengths, "transfer_lengthparam")

    transfer_lengths = pd.DataFrame(
    index=pd.MultiIndex.from_product([link_names, link_types]))

    #With the line corridors now defined, we can start adding lines to be optimized to the model.
    #Note: we have to change these
    transport_techs = ["HV"] #other options: "HVDC", "Pipeline_CH4", "Pipeline_H2", "Pipeline_H2_retrofit"
    #in previous version of this function, the capacities come from a csv file
    # question: should I have a different capacity for each tech? how do I define that?
    # question: here I have transport_techs, how do I assign different transport techs to different links (to
    # where i want to have the hvDC? also should i just split the line like that and create extra virtual nodes?
    #  or how do i go about that line having only one fragment being DC)
    link_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product([link_names, m["Base"].set.yearssel, transport_techs])
    )
    link_caps.loc[idx[:,:,"HV"],["linksUpperLimit"]] = 100 # Allow to build 100 GW for all links as the upper limit
    

    m["Base"].parameter.add(link_caps, "transfer_linksparam")
    tech_params = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m["Base"].set.yearssel])
    )
    tech_params.loc[:, "lifeTime"] = 40
    tech_params.loc[:, "flowUpperLimit"] = 1

    m["Base"].parameter.add(tech_params, "transfer_techparam")

    # Define the commodity and rated capacity of the network technology
    # "transfer_coefficient"
    commodity = ["Elec"]

    transfer_coefficient = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m["Base"].set.yearssel, commodity])
    )
    transfer_coefficient["coefficient"] = 1  # GWh / h per line

    m["Base"].parameter.add(transfer_coefficient, "transfer_coefficient")
    transfer_coefficient


    
    # Define the losses for the converter stations
    # "transport_coefPerFlow"
    coef_per_flow = pd.DataFrame(
        index=pd.MultiIndex.from_product([transport_techs, m["Base"].set.yearssel, commodity])
    )
    coef_per_flow[
        "coefPerFlow"
    ] = -0.014  # electrical losses of 14 MWh/h for each flow of 1 GWh/h

    m["Base"].parameter.add(coef_per_flow, "transfer_coefperflow")

    # Define the losses for the lines per km
    # "transport_coefPerDistance"
    coef_per_dist = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, m["Base"].set.yearssel, commodity, link_types]
        )
    )
    coef_per_dist.loc[
        idx[:, :, :, "land"], idx["coefPerLength"]
    ] = (
        -0.00004
    )  # electrical losses of 40 kWh / h for each flow of 1 GWh / h and 1 km line length ~ 24 MWh / h for 600 km distance
    coef_per_dist.loc[idx[:, :, :, "sea"], idx["coefPerLength"]] = -0.00003

    m["Base"].parameter.add(coef_per_dist, "transfer_coefperlength")
    coef_per_dist

    # Define indicators for each line built (for HVDC this is an AC/DC converter station 
    # at the beginning and end of the line)
    # "accounting_transferlinks"
    cost_indicators = ["Invest", "OMFix"]
    area = ["global"]

    transfer_indicators = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, area, ["horizon"], transport_techs, m["Base"].set.yearssel]
        )
    )
    
    transfer_indicators.index.set_names(["indicator","regionscope","timescope","techs","years"], inplace=True)
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "perLinkBuild"] = 180
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "interest"] = 0.06
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "amorTime"] = 40
    transfer_indicators.loc[idx["Invest", "global", "horizon"], "useAnnuity"] = 1
    transfer_indicators.loc[idx["OMFix", "global", "horizon"], "perLinkTotal"] = 1.8
    transfer_indicators = transfer_indicators.fillna(0)

    m["Base"].parameter.add(transfer_indicators, "accounting_transferlinks")
    transfer_indicators

    # Define indicators for each line-km built 
    # (this needs the additional set for distance-type modifiers, such as land and sea)
    # ""accounting_transferperlength""
    indicators_distance = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, area, ["horizon"], transport_techs, m["Base"].set.yearssel, link_types]
        )
    )
    
    indicators_distance.index.set_names(["indicator","regionscope","timescope","techs","years","linktype"], inplace=True)
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "perLengthBuild"
    ] = 0.544
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "land"], "useAnnuity"] = 1
    indicators_distance.loc[idx["OMFix", "global", "horizon", :, :, "land"], "perLengthTotal"
    ] = 0.00544

    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "perLengthBuild"
    ] = 0.975
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "interest"] = 0.06
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "amorTime"] = 40
    indicators_distance.loc[idx["Invest", "global", "horizon", :, :, "sea"], "useAnnuity"] = 1
    indicators_distance.loc[idx["OMFix", "global", "horizon", :, :, "sea"], "perLengthTotal"
    ] = 0.00975
    indicators_distance = indicators_distance.fillna(0)

    m["Base"].parameter.add(indicators_distance, "accounting_transferperlength")
    indicators_distance

def add_accounting(m):
    # Bounds: SystemCost is the main objective over the whole horizon
    accounting_indicatorBounds = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["global"], ["horizon"], ["SystemCost"]],
            names=["regionscope", "timescope", "indicator"],
        )
    )
    accounting_indicatorBounds["obj"] = -1          # minimize SystemCost
    accounting_indicatorBounds["discount"] = 0.02   # social discount rate
    m["Base"].parameter.add(accounting_indicatorBounds, "accounting_indicatorbounds")

    # Per-indicator weights: SystemCost aggregates these sub-indicators
    # Index: (SystemCostIndicator, componentIndicator, regionscope, timescope)
    base_idx = pd.MultiIndex.from_product(
        [
            ["SystemCost"],                                  # aggregate indicator
            ["Invest", "OMFix", "FuelCost", "ImportCost", "SlackCost"],  # components
            ["global"],
            ["horizon"],                                     # all optimisation years
        ],
        names=["SystemCostIndicator", "indicator", "regionscope", "timescope"],
    )
    accounting_perIndicator = pd.DataFrame(index=base_idx)
    accounting_perIndicator["perIndicator"] = 1.0

    m["Base"].parameter.add(accounting_perIndicator, "accounting_perindicator")

def validate_scope(m):
    m["Base"].infer_set_data()
    nodes_data = set(m["Base"].set.nodesdata)
    nodes_model = set(m["Base"].set.nodesmodel)
    print(f"Not including modes nodes: {', '.join(sorted(nodes_data - nodes_model))}")

# scenario modifiers

def modify_inflow_file(m, data_dir, new_inflow_path, tag="alt-inflow", cascade_mode=True):
    """
    Replace the hydropower inflow file in an existing REMix model with an alternative one.

    This function:
      • Reads the new inflow CSV (same format as REMix input)
      • Applies sign conventions and fixed profile flags
      • Writes new sourcesink_profile and sourcesink_config files
      • Creates a new folder under data_dir/<tag> for the modified inputs
      • Prints diagnostics about the inflows

    Parameters
    ----------
    m : dict
        REMix instance dictionary (m["Base"])
    data_dir : Path
        Base data directory (e.g. Path(".../project/hadi/pypsa-cascade/data"))
    new_inflow_path : str or Path
        Path to the new inflow CSV (either _INFLOWS_cumecs.csv or _ENERGY_GW.csv)
    tag : str, optional
        Name for the output subfolder (default: 'alt-inflow')
    cascade_mode : bool, optional
        If True, keeps multiple inflow commodities (e.g. Tekapo_in, Roxburgh_in, ...).
        If False, aggregates or renames to 'Water_in'.
    """

    import pandas as pd
    import numpy as np
    from pathlib import Path

    new_inflow_path = Path(new_inflow_path)
    out_dir = Path(data_dir) / tag
    out_dir.mkdir(exist_ok=True, parents=True)

    print(f"\n--- Processing inflow modifier ({tag}) ---")
    print(f"Loading new inflow file: {new_inflow_path}")

    # --- Load inflow CSV ---
    inflow = pd.read_csv(new_inflow_path)

    # Detect and align expected column names
    rename_map = {}
    if "carrier" in inflow.columns and "commodity" not in inflow.columns:
        rename_map["carrier"] = "commodity"
    if "region" in inflow.columns and "node" not in inflow.columns:
        rename_map["region"] = "node"
    if rename_map:
        inflow.rename(columns=rename_map, inplace=True)

    # Ensure standard index format
    expected_cols = ["node", "year", "sector", "commodity"]
    if not all(c in inflow.columns for c in expected_cols):
        raise ValueError(
            f"Inflow CSV missing required columns. Expected: {expected_cols}, got: {list(inflow.columns)}"
        )

    inflow = inflow.set_index(expected_cols)
    inflow = inflow.select_dtypes(include=[np.number]).dropna(how="all")

    # Drop or rename based on mode
    if cascade_mode:
        inflow = inflow.loc[inflow.index.get_level_values(3) != "Water_in"]
    else:
        inflow.index = inflow.index.set_levels(
            ["Water_in" if x == "HydroInflow" else x for x in inflow.index.levels[3]],
            level=3,
        )

    # --- Apply REMix conventions ---
    inflow *= -1  # inflows are positive sources in REMix

    inflow["type"] = "fixed"
    inflow_fixed = inflow.set_index("type", append=True).round(6)

    # --- Write out ---
    profile_path = out_dir / "sourcesink_profile.csv"
    config_path = out_dir / "sourcesink_config.csv"

    inflow_fixed.to_csv(profile_path)

    inflow_cfg = pd.DataFrame(index=inflow.index)
    inflow_cfg["usesFixedProfile"] = 1
    inflow_cfg = inflow_cfg.loc[inflow.select_dtypes(include="number").sum(axis=1) != 0]
    inflow_cfg.to_csv(config_path)

    # --- Diagnostics ---
    inflow_names = inflow.index.get_level_values(3).unique().tolist()
    total_gwh = inflow.select_dtypes(include=[np.number]).sum().sum()

    print(f"\nInflow modification complete.")
    print(f"  Total inflow commodities: {len(inflow_names)}")
    print(f"  First inflows: {', '.join(inflow_names[:8])}{'...' if len(inflow_names) > 8 else ''}")
    print(f"  Total inflow energy (sum over all hours and nodes): {total_gwh:,.2f} GWh")
    print(f"  Files written to: {out_dir}\n")

    # --- Optionally, update model instance directly (if desired) ---
    if m is not None:
        m["Base"].profile.add(inflow_fixed, "sourcesink_profile")
        m["Base"].parameter.add(inflow_cfg, "sourcesink_config")
        m["Base"].set.add(list(inflow.index.get_level_values(0)), "nodesdata")
        m["Base"].set.add(list(inflow.index.get_level_values(1)), "years")
        print("  Updated inflow data added to current model instance (in-memory).")

    return out_dir

def modify_renewable_costs(m, data_dir, tech="wind_onshore", factor=0.8, tag="wind-cost-low"):
    """
    Modify renewable technology investment cost (e.g. wind, PV).
    Reloads renewable cost parameters from base files, applies scaling factor, and writes modified accounting files to subfolder.

    Parameters
    ----------
    m : dict
        REMix instance dictionary
    data_dir : Path
        Base data directory
    tech : str
        Technology name to modify ('wind_onshore', 'pv_central_fixed', etc.)
    factor : float
        Scaling factor (e.g., 0.8 = 20% cheaper)
    tag : str
        Folder name for scenario output
    """

    print(f"Processing modifier: renewable cost ({tech}, {factor*100:.0f}% level)...")
    from pathlib import Path

    out_dir = data_dir / tag
    out_dir.mkdir(exist_ok=True)

    # read existing accounting file
    acc_path = data_dir / "accounting_converterunits.csv"
    acc = pd.read_csv(acc_path, index_col=list(range(5)))

    # find investment rows for the given tech and apply scaling
    mask = (acc.index.get_level_values(0) == "Invest") & (acc.index.get_level_values(3) == tech)
    acc.loc[mask, "perUnitBuild"] *= factor

    acc.to_csv(out_dir / "accounting_converterunits.csv")
    print(f"  File written: {out_dir/'accounting_converterunits.csv'}")

def modify_tech_cost_by_year(m, data_dir, tech, year_costs, tag=None):
    """
    Modify investment costs for a technology for specific years.

    Parameters
    ----------
    m : dict
        REMix instance dictionary
    data_dir : Path
        Base data directory
    tech : str
        Technology to modify (e.g. 'Electrolyser', 'Battery')
    year_costs : dict[int, float]
        Mapping of year -> new perUnitBuild value (absolute, not factor)
    tag : str or None
        Folder name for scenario output (defaults to tech name)
    """

    print(f"Processing modifier: {tech} custom yearly cost changes...")

    out_dir = data_dir / (tag or f"{tech}-custom-costs")
    out_dir.mkdir(exist_ok=True)

    acc_path = data_dir / "accounting_converterunits.csv"
    acc = pd.read_csv(acc_path, index_col=list(range(5)))

    for year, new_val in year_costs.items():
        mask = (
            (acc.index.get_level_values(0) == "Invest")
            & (acc.index.get_level_values(3) == tech)
            & (acc.index.get_level_values(4).astype(int) == int(year))
        )
        acc.loc[mask, "perUnitBuild"] = new_val
        print(f"  - updated {tech} investment cost for {year} → {new_val}")

    acc.to_csv(out_dir / "accounting_converterunits.csv")
    print(f"  File written: {out_dir/'accounting_converterunits.csv'}")


#other carriers

def add_h2_annual_demand_from_excel(m, carriers_long: pd.DataFrame):
    """
    Adds annual hydrogen demand as exogenous sinks on the real H2 commodity.
    Supply must come from the model (electrolyers, etc.); no H2 imports.

    Aggregates original Excel regions (e.g. West Coast + Canterbury) into
    single REMix nodes (e.g. CAN) before writing parameters, ensuring
    one row per (node, year, sector, commodity).

    Pattern:
      - sourcesink_annualsum: upper = -Demand (GWh), lower not used
      - sourcesink_config:    usesUpperSum = 1, usesLowerSum = 0, usesUpperProfile = 1
    """

    df = carriers_long.copy()
    df["Carrier_key"] = df["Carrier"].astype(str).str.strip().str.lower()

    # Keep only hydrogen-type carriers
    h2_keys = {"hydrogen", "h2", "h2-feedstock"}
    df = df.loc[df["Carrier_key"].isin(h2_keys)].copy()
    if df.empty:
        print("No H2 carriers found in Excel; skipping H2 annual demand.")
        return

    # Map to REMix commodity H2
    df["commodity"] = "H2"

    print("\n[DEBUG] H2 demand before aggregation:")
    print(f"  Rows: {len(df)}")
    print(f"  By (REMixRegion, Year, Sector):")
    pre_agg = (
        df.groupby(["REMixRegion", "Year", "Sector"])["Demand"]
        .sum()
        .sort_index()
    )
    print(pre_agg.to_string())

    # AGGREGATE to ensure 1 row per (node, year, sector, commodity)
    # This sums original regions that map to same REMixRegion
    grp = (
        df.groupby(["REMixRegion", "Year", "Sector", "commodity"])["Demand"]
        .sum()
        .reset_index()
    )

    print(f"\n[DEBUG] H2 demand after aggregation:")
    print(f"  Aggregated rows: {len(grp)}")
    print(f"  Total demand: {grp['Demand'].sum():,.2f} GWh")
    print(grp.to_string(index=False))

    # Index: (nodesdata, years, sector, commodity)
    idx_ss = pd.MultiIndex.from_frame(
        grp[["REMixRegion", "Year", "Sector", "commodity"]]
    )
    idx_ss.names = ["nodesdata", "years", "sector", "commodity"]

    # Annual upper sums: negative for demand (GWh per year)
    ss = pd.DataFrame(index=idx_ss)
    ss["upper"] = -grp["Demand"].values
    m["Base"].parameter.add(ss, "sourcesink_annualsum")

    # Configuration: upper sum only, no lower, no profile
    cfg = pd.DataFrame(index=idx_ss)
    cfg["usesUpperSum"] = 1
    cfg["usesLowerSum"] = 0
    cfg["usesUpperProfile"] = 0
    m["Base"].parameter.add(cfg, "sourcesink_config")

    # Ensure nodes and years are in sets
    m["Base"].set.add(sorted(grp["REMixRegion"].unique().tolist()), "nodesdata")
    m["Base"].set.add(sorted(grp["Year"].unique().tolist()), "years")

    print(
        f"\n[DEBUG] Hydrogen annual demand added: {len(grp)} aggregated rows "
        f"(total {grp['Demand'].sum():,.2f} GWh)"
    )

def fuel_from_excel(m, carriers_long: pd.DataFrame):
    """
    entry point for paid fuels (fossil/bio/wood) via FuelConsumer (FuelService_* demand)
      - Optional REfuel/CH4 e-fuel demand

    Aggregates original Excel regions (e.g. West Coast + Canterbury) into
    single REMix nodes (e.g. CAN) before writing parameters, ensuring
    one row per (node, year, sector, commodity).
    """

    df = carriers_long.copy()
    df["Carrier_key"] = df["Carrier"].astype(str).str.strip().str.lower()

    # ---- 1. Classification of carriers --------------------------------
    paid_carrier_to_commodity = {
        "coal": "Coal_fossil",
        "fossil (gas)": "Gas_fossil",
        "fossil (lf)": "LF_fossil",
        "biofuel (lf)": "LF_bio",
        "biofuel (gas)": "Gas_bio",
        "wood": "Wood",
    }

    efuel_keys = {}

    # Split into paid and e-fuel groups
    df_paid = df.loc[df["Carrier_key"].isin(paid_carrier_to_commodity.keys())].copy()
    df_efuel = df.loc[df["Carrier_key"].isin(efuel_keys.keys())].copy()

    # Restrict to optimisation years
    df_paid = df_paid.loc[df_paid["Year"].isin(yrs_sel)].copy()
    df_efuel = df_efuel.loc[df_efuel["Year"].isin(yrs_sel)].copy()

    # ---- 2. Paid fuel demand via FuelConsumer -------------------------
    if df_paid.empty:
        print("No paid fuel carriers found in Excel; skipping FuelConsumer block.")
        return

    df_paid["commodity"] = df_paid["Carrier_key"].map(paid_carrier_to_commodity)
    df_paid["service_commodity"] = df_paid["commodity"].apply(lambda c: f"FuelService_{c}")

    # AGGREGATE to ensure 1 row per (node, year, sector, commodity)
    # This sums Canterbury + West Coast -> CAN totals
    grp_paid = (
        df_paid.groupby(["REMixRegion", "Year", "Sector", "service_commodity"])["Demand"]
        .sum()
        .reset_index()
    )

    if grp_paid.empty:
        print("No paid fuel demand in selected years; skipping FuelConsumer block.")
        return

    idx_ss = pd.MultiIndex.from_frame(
        grp_paid[["REMixRegion", "Year", "Sector", "service_commodity"]]
    )
    idx_ss.names = ["nodesdata", "years", "sector", "commodity"]

    ss = pd.DataFrame(index=idx_ss)
    ss["upper"] = -grp_paid["Demand"].values  # negative for demand
    m["Base"].parameter.add(ss, "sourcesink_annualsum")

    cfg = pd.DataFrame(index=idx_ss)
    cfg["usesUpperSum"] = 1
    cfg["usesLowerSum"] = 0
    cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(cfg, "sourcesink_config")

    m["Base"].set.add(sorted(grp_paid["REMixRegion"].unique().tolist()), "nodesdata")
    m["Base"].set.add(sorted(grp_paid["Year"].unique().tolist()), "years")

    # 2b. FuelConsumer converter tech (capacity non-binding, per-node where needed)
    fuel_tech = "FuelConsumer"
    vintages = yrs_to_calc

    tech = pd.DataFrame(index=pd.MultiIndex.from_product([[fuel_tech], vintages]))
    tech["lifeTime"] = 1
    tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(tech, "converter_techparam")

    nodes_with_demand = sorted(grp_paid["REMixRegion"].unique().tolist())
    cap = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes_with_demand, yrs_to_calc, [fuel_tech]])
    )
    cap["unitsBuild"] = 0.0          # start at zero
    cap["unitsUpperLimit"] = 1e6     # effectively unconstrained
    cap["noExpansion"] = 0
    m["Base"].parameter.add(cap, "converter_capacityparam")

    # 2c. Converter coefficients: one activity per fuel
    # Find all unique base commodities involved
    # We must map back from FuelService_X to X or take from original df_paid
    # simpler: just list unique commodities from the aggregated set
    # But wait, we need the BASE commodity (Wood) not service commodity (FuelService_Wood)
    # The 'service_commodity' string always starts with 'FuelService_'
    unique_services = sorted(grp_paid["service_commodity"].unique())
    unique_fuels = [s.replace("FuelService_", "") for s in unique_services]

    fuels = unique_fuels
    activities = [f"{f}_cons" for f in fuels]
    services = unique_services  # matches order of fuels
    commodities = fuels + services

    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([[fuel_tech], vintages, activities, commodities])
    )

    for f, act, serv in zip(fuels, activities, services):
        # consume 1 GWh of fuel
        coef.loc[idx[:, :, act, f], "coefficient"] = -1.0
        # produce 1 GWh of matching service
        coef.loc[idx[:, :, act, serv], "coefficient"] = 1.0

    m["Base"].parameter.add(coef, "converter_coefficient")

    # 2d. Accounting per activity: fuel cost + CO2 (unchanged)
    co2_kt_per_gwh = {
        "Coal_fossil": 0.34,
        "Gas_fossil":  0.20,
        "LF_fossil":   0.27,
        "LF_bio":      0.0,
        "Gas_bio":     0.0,
        "Wood":        0.0,
    }

    commodity_to_csiro = {
        "Coal_fossil": "Black Coal",
        "Gas_fossil":  "Gas",
        "LF_fossil":   "Liquid Fuel",
        "LF_bio":      "Liquid Fuel",
        "Gas_bio":     "Gas",
        "Wood":        "Biomass",
    }

    acc_idx = pd.MultiIndex.from_product(
        [
            ["FuelCost", "CO2_emission"],
            ["global"],
            ["horizon"],
            [fuel_tech],
            vintages,
            activities,
        ],
        names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
    )

    acc = pd.DataFrame(index=acc_idx)

    for y in vintages:
        suffix = "early" if y < 2050 else "2050"
        for f, act in zip(fuels, activities):
            fuel_key = commodity_to_csiro.get(f, None)
            if fuel_key is None:
                continue
            csiro_key = f"{fuel_key}_{suffix}"
            cost = aud_per_gj_to_meur_per_gwh(csiro_fuel_aud_per_gj[csiro_key])
            co2 = co2_kt_per_gwh.get(f, 0.0)

            acc.loc[("FuelCost", "global", "horizon", fuel_tech, y, act), "perActivity"] = cost
            acc.loc[("CO2_emission", "global", "horizon", fuel_tech, y, act), "perActivity"] = co2

    acc = acc.fillna(0.0)
    m["Base"].parameter.add(acc, "accounting_converteractivity")

    print(
        f"FuelConsumer block added for {len(fuels)} fuels, "
        f"{len(nodes_with_demand)} nodes and {len(grp_paid)} aggregated demand rows."
    )

def add_fuel_imports(m, carriers_long: pd.DataFrame):
    """
    Allow purchasable annual imports of primary fuels (Coal_fossil, Wood, Gas_fossil, Gas_bio,
    LF_fossil, LF_bio) in each node/year where there is FuelService demand.

    - sourcesink_annualsum: positive upper (supply cap) on the FUEL commodity
    - sourcesink_config:    usesUpperSum = 1, usesLowerSum = 0, usesUpperProfile = 1
    - accounting_converteractivity: FuelImportCost + CO2_emission per unit imported fuel
    """
    # same mapping used in fuel_and_efuel_from_excel
    paid_carrier_to_commodity = {
        "coal": "Coal_fossil",
        "fossil (gas)": "Gas_fossil",
        "fossil (lf)": "LF_fossil",
        "biofuel (lf)": "LF_bio",
        "biofuel (gas)": "Gas_bio",
        "wood": "Wood",
        # "e-fuel (lf)": "REfuel",
        # "e-fuel (gas)": "CH4",
    }

    df = carriers_long.copy()
    df["Carrier_key"] = df["Carrier"].astype(str).str.strip().str.lower()
    df = df.loc[df["Carrier_key"].isin(paid_carrier_to_commodity.keys())].copy()
    if df.empty:
        print("No paid fuel carriers found in Excel; skipping fuel import sinks.")
        return

    df["commodity"] = df["Carrier_key"].map(paid_carrier_to_commodity)
    df = df.loc[df["Year"].isin(yrs_sel)].copy()
    if df.empty:
        print("No paid fuel demand in selected years; skipping fuel import sinks.")
        return

    # aggregate total fuel demand per node/year/commodity to set a sensible cap
    agg = (
        df.groupby(["REMixRegion", "Year", "commodity"], as_index=False)["Demand"]
          .sum()
    )
    # import cap as factor * demand (so it is never binding, just finite)
    safety_factor = 2.0
    agg["ImportCap"] = safety_factor * agg["Demand"]

    idx_ss = pd.MultiIndex.from_arrays(
        [agg["REMixRegion"], agg["Year"], ["FuelImport"] * len(agg), agg["commodity"]],
        names=["nodesdata", "years", "sector", "commodity"],
    )

    ss = pd.DataFrame(index=idx_ss)
    # Positive upper: this is a source to the system
    ss["upper"] = agg["ImportCap"].values
    m["Base"].parameter.add(ss, "sourcesink_annualsum")

    cfg = pd.DataFrame(index=idx_ss)
    cfg["usesUpperSum"] = 1
    cfg["usesLowerSum"] = 0
    cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(cfg, "sourcesink_config")

    m["Base"].set.add(list(agg["REMixRegion"].unique()), "nodesdata")
    m["Base"].set.add(list(agg["Year"].unique()), "years")

    # Cost and CO2 per imported GWh, using CSIRO + same factors as FuelConsumer
    co2_kt_per_gwh = {
        "Coal_fossil": 0.34,
        "Gas_fossil":  0.20,
        "LF_fossil":   0.27,
        "LF_bio":      0.0,
        "Gas_bio":     0.0,
        "Wood":        0.0,
    }
    commodity_to_csiro = {
        "Coal_fossil": "Black Coal",
        "Gas_fossil":  "Gas",
        "LF_fossil":   "Liquid Fuel",
        "LF_bio":      "Liquid Fuel",
        "Gas_bio":     "Gas",
        "Wood":        "Biomass",
    }

    # Model imports as a "virtual" technology FuelImport with one activity per fuel
    fuel_import_tech = "FuelImport"
    vintages = yrs_to_calc
    activities = [f"{c}_imp" for c in sorted(co2_kt_per_gwh.keys())]

    # tech and capacity (non-binding)
    tech = pd.DataFrame(index=pd.MultiIndex.from_product([[fuel_import_tech], vintages]))
    tech["lifeTime"] = 1
    tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(tech, "converter_techparam")

    nodes = sorted(agg["REMixRegion"].unique().tolist())
    cap = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, yrs_to_calc, [fuel_import_tech]])
    )
    cap["unitsBuild"] = 0.0
    cap["unitsUpperLimit"] = 1e9
    cap["noExpansion"] = 0
    m["Base"].parameter.add(cap, "converter_capacityparam")

    # coefficients: +1 fuel commodity per unit of import activity
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [[fuel_import_tech], vintages, activities, sorted(co2_kt_per_gwh.keys())]
        )
    )
    for fuel, act in zip(sorted(co2_kt_per_gwh.keys()), activities):
        coef.loc[idx[:, :, act, fuel], "coefficient"] = 1.0
    m["Base"].parameter.add(coef, "converter_coefficient")

    # accounting per activity: FuelImportCost + CO2_emission
    acc_idx = pd.MultiIndex.from_product(
        [
            ["FuelImportCost", "CO2_emission"],
            ["global"],
            ["horizon"],
            [fuel_import_tech],
            vintages,
            activities,
        ],
        names=["indicator", "regionscope", "timescope", "techs", "years", "activities"],
    )
    acc = pd.DataFrame(index=acc_idx)

    fuels = sorted(co2_kt_per_gwh.keys())
    for y in vintages:
        suffix = "early" if y < 2050 else "2050"
        for fuel, act in zip(fuels, activities):
            fuel_key = commodity_to_csiro.get(fuel, None)
            if fuel_key is None:
                continue
            csiro_key = f"{fuel_key}_{suffix}"
            cost = aud_per_gj_to_meur_per_gwh(csiro_fuel_aud_per_gj[csiro_key])
            co2 = co2_kt_per_gwh[fuel]

            acc.loc[("FuelImportCost", "global", "horizon", fuel_import_tech, y, act), "perActivity"] = cost
            acc.loc[("CO2_emission", "global", "horizon", fuel_import_tech, y, act), "perActivity"] = co2

    acc = acc.fillna(0.0)
    m["Base"].parameter.add(acc, "accounting_converteractivity")

    fuel_list = sorted(fuels)
    node_list = sorted(nodes)
    vintage_list = sorted(vintages)

    print(
        f"Fuel import sinks added for {len(fuel_list)} fuels {fuel_list}, "
        f"{len(node_list)} nodes {node_list} and years {vintage_list}."
    )

def add_efuel_annual_demand_from_excel(m, carriers_long: pd.DataFrame):
    """
    Adds annual REfuel and CH4 demand as exogenous sinks on the real commodities.
    Supply comes from add_dac(), add_ftropsch_syn(), add_methanizer().

    Aggregates original Excel regions into REMix nodes before writing.

    Pattern:
      - sourcesink_annualsum: upper = -Demand (GWh), lower not used
      - sourcesink_config:    usesUpperSum = 1, usesLowerSum = 0, usesUpperProfile = 1
    """

    df = carriers_long.copy()
    df["Carrier_key"] = df["Carrier"].astype(str).str.strip().str.lower()

    # E-fuel carrier mapping
    efuel_carrier_to_commodity = {
        "e-fuel (lf)": "REfuel",
        "e-fuel (gas)": "CH4",
        "refuel": "REfuel",
        "ch4": "CH4",
    }

    # Filter to e-fuel carriers only
    df = df.loc[df["Carrier_key"].isin(efuel_carrier_to_commodity.keys())].copy()
    if df.empty:
        return

    # Map to commodities
    df["commodity"] = df["Carrier_key"].map(efuel_carrier_to_commodity)

    # Filter to optimisation years
    df = df.loc[df["Year"].isin(yrs_sel)].copy()
    if df.empty:
        return

    # Aggregate: ensure 1 row per (node, year, sector, commodity)
    grp = (
        df.groupby(["REMixRegion", "Year", "Sector", "commodity"])["Demand"]
        .sum()
        .reset_index()
    )

    # Index: (nodesdata, years, sector, commodity)
    idx_ss = pd.MultiIndex.from_frame(
        grp[["REMixRegion", "Year", "Sector", "commodity"]]
    )
    idx_ss.names = ["nodesdata", "years", "sector", "commodity"]

    # Annual upper sums: negative for demand
    ss = pd.DataFrame(index=idx_ss)
    ss["upper"] = -grp["Demand"].values
    m["Base"].parameter.add(ss, "sourcesink_annualsum")

    # Configuration
    cfg = pd.DataFrame(index=idx_ss)
    cfg["usesUpperSum"] = 1
    cfg["usesLowerSum"] = 0
    cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(cfg, "sourcesink_config")

    # Ensure nodes and years in sets
    m["Base"].set.add(sorted(grp["REMixRegion"].unique().tolist()), "nodesdata")
    m["Base"].set.add(sorted(grp["Year"].unique().tolist()), "years")


def add_efuel_perfect_storage(m):
    """
    Adds perfect (100% efficient, zero-cost) storage reservoirs for REfuel and CH4.
    
    Follows the H2_storage pattern but with:
      - Zero capital and operational costs
      - 100% round-trip efficiency (via converter coefficients)
      - Annual cycle (system balances over the year)
      - Effectively unlimited capacity per node
    
    Creates:
      - Converter (compressor/pump): charges/discharges storage
      - Storage reservoir: holds fuel throughout the year
      - Both with zero cost
    """
    
    efuels = ["REfuel", "CH4"]
    h2_stor_nodes = sorted(m["Base"].set.nodesdata)
    storage_suffix = "_Storage"
    converter_suffix = "_Compressor"
    
    print("\n[E-fuel Storage] Adding perfect storage for REfuel and CH4...")
    print(f"  E-fuels: {efuels}")
    print(f"  Nodes: {h2_stor_nodes}")

    # --- For each e-fuel, create storage + compressor pair ---
    for ef in efuels:
        stor_tech_name = f"{ef}{storage_suffix}"
        conv_tech_name = f"{ef}{converter_suffix}"
        stored_commodity = f"{ef}_stored"  # e.g. REfuel_stored
        
        print(f"\n  Setting up {ef}:")
        print(f"    Storage tech: {stor_tech_name}")
        print(f"    Converter tech: {conv_tech_name}")
        
        # ============================================================
        # 1. CONVERTER (compressor/pump) - charges/discharges storage
        # ============================================================
        conv_tech = pd.DataFrame(
            index=pd.MultiIndex.from_product([[conv_tech_name], yrs_to_calc])
        )
        conv_tech.index.names = ["techs", "years"]
        conv_tech.loc[:, "lifeTime"] = 40
        m["Base"].parameter.add(conv_tech, "converter_techparam")
        
        # Converter capacity (units)
        conv_cap = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [h2_stor_nodes, yrs_to_calc, [conv_tech_name]]
            )
        )
        conv_cap.index.names = ["nodesdata", "years", "techs"]
        conv_cap.loc[:, "unitsUpperLimit"] = 1000.0  # unlimited effectively
        m["Base"].parameter.add(conv_cap, "converter_capacityparam")
        
        # Converter coefficients: charge/discharge 100% efficient
        conv_coef = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [[conv_tech_name], yrs_to_calc, ["Charge", "Discharge"], 
                 [ef, stored_commodity, "dummy"]]
            )
        )
        conv_coef.index.names = ["techs", "years", "operation", "commodity"]
        
        # Charge: take 1 GWh of fuel, store as stored commodity at 100% efficiency
        conv_coef.loc[(conv_tech_name, slice(None), "Charge", ef), "coefficient"] = -1.0
        conv_coef.loc[(conv_tech_name, slice(None), "Charge", stored_commodity), "coefficient"] = 1.0
        
        # Discharge: take 1 GWh stored, release as fuel at 100% efficiency
        conv_coef.loc[(conv_tech_name, slice(None), "Discharge", stored_commodity), "coefficient"] = -1.0
        conv_coef.loc[(conv_tech_name, slice(None), "Discharge", ef), "coefficient"] = 1.0
        
        # Drop dummy rows
        conv_coef = conv_coef.dropna(subset=["coefficient"])
        m["Base"].parameter.add(conv_coef, "converter_coefficient")
        
        # Converter accounting: ZERO COST (perfect storage)
        conv_acc = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [["Invest", "OMFix"], ["global"], ["horizon"], [conv_tech_name], yrs_to_calc]
            )
        )
        conv_acc.index.names = ["indicator", "regionscope", "timescope", "techs", "years"]
        conv_acc.loc[idx["Invest", :, :, :, :], "perUnitBuild"] = 0.0
        conv_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
        conv_acc.loc[idx["Invest", :, :, :, :], "amorTime"] = 40
        conv_acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06
        conv_acc.loc[idx["OMFix", :, :, :, :], "perUnitTotal"] = 0.0
        m["Base"].parameter.add(conv_acc, "accounting_converterunits")
        
        # ============================================================
        # 2. STORAGE (reservoir)
        # ============================================================
        stor_tech = pd.DataFrame(
            index=pd.MultiIndex.from_product([[stor_tech_name], yrs_to_calc])
        )
        stor_tech.index.names = ["techs", "years"]
        stor_tech.loc[:, "lifeTime"] = 40
        stor_tech.loc[:, "levelUpperLimit"] = 1  # state of charge limit (0-1)
        m["Base"].parameter.add(stor_tech, "storage_techparam")
        
        # Storage size (capacity per unit in GWh)
        stor_size = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [[stor_tech_name], yrs_to_calc, [stored_commodity]]
            )
        )
        stor_size.index.names = ["techs", "years", "commodity"]
        stor_size.loc[:, "size"] = 1000.0  # 1000 GWh per unit (effectively unlimited)
        m["Base"].parameter.add(stor_size, "storage_sizeparam")
        
        # Storage reservoir (units per node)
        stor_res = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [h2_stor_nodes, yrs_to_calc, [stor_tech_name]]
            )
        )
        stor_res.index.names = ["nodesdata", "years", "techs"]
        stor_res.loc[:, "unitsUpperLimit"] = 1.0  # 1 unit per node (holds 1000 GWh)
        m["Base"].parameter.add(stor_res, "storage_reservoirparam")
        
        # Storage accounting: ZERO COST
        stor_acc = pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [["Invest", "OMFix"], ["global"], ["horizon"], [stor_tech_name], yrs_to_calc]
            )
        )
        stor_acc.index.names = ["indicator", "regionscope", "timescope", "techs", "years"]
        stor_acc.loc[idx["Invest", :, :, :, :], "perUnitBuild"] = 0.0
        stor_acc.loc[idx["Invest", :, :, :, :], "useAnnuity"] = 1
        stor_acc.loc[idx["Invest", :, :, :, :], "amorTime"] = 40
        stor_acc.loc[idx["Invest", :, :, :, :], "interest"] = 0.06
        stor_acc.loc[idx["OMFix", :, :, :, :], "perUnitTotal"] = 0.0
        m["Base"].parameter.add(stor_acc, "accounting_storageunits")
        
        print(f"    ✓ Storage + converter added")
    
    print(f"\n[E-fuel Storage] Perfect storage setup complete")

# #%%



# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
if __name__ == "__main__":
    start_time = time.time()

    m = {"Base": Instance(index_names=False, datadir=data_dir)}

    # scope first
    add_scope(m)

    # hourly electricity + inflows
    add_demand(m)

    # carriers from Excel
    excel_path = path_demand / carriers_excel_file
    carriers_long = read_carriers_excel_long(excel_path, base_scenario)

    # annual hydrogen demand (perfect battery)
    add_h2_annual_demand_from_excel(m, carriers_long)
    fuel_from_excel(m, carriers_long)
    add_fuel_imports(m, carriers_long)
    add_efuel_annual_demand_from_excel(m, carriers_long)
    # add_efuel_perfect_storage(m)
    # add_perfect_xstorage(m)

    # supply side blocks
    if include_renewables:
        add_renewables(m)

    if include_geothermal:
        add_geothermal(m)

    if include_hydro_simple:
        add_hydro(m)

    if include_lithium_batteries:
        add_lithium_batteries(m)

    if include_thermal:
        add_thermal(m)

    if include_gas_turbines:
        add_gas_turbines(m)

    # hydrogen + efuels chain
    if include_electrolyser:
        add_electrolyser(m)

    if include_dac:
        add_dac(m)

    if include_ftropsch:
        add_ftropsch_syn(m)

    if include_methanizer:
        add_methanizer(m)

    if include_h2_storage:
        add_h2_storage(m)

    if include_h2_ccgt:
        add_H2_CCGT(m)

    if include_h2_fc:
        add_H2_FC(m)

    if include_network:
        add_network(m)

    if include_accounting:
        add_accounting(m)

    # carbon policy
    if use_emission_limit_2050:
        add_emission_limit(m)
    if use_emission_slack_2050:
        add_emission_slack(m, year=2050, slack_cost_eur_per_tco2=slack_cost_eur_per_tco2)

    validate_scope(m)

    m["Base"].write(project_path=data_dir, fileformat="csv", float_format="{:.4g}".format)

    total = time.time() - start_time
    print("\nDone.")
    print("Wrote inputs to:", data_dir)
    print(f"Total time: {total:.1f} s")
