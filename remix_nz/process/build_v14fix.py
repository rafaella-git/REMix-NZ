from networkx import nodes
import numpy as np
import pandas as pd
import time
import shutil # warning: we use this to clean up the data folder before writing
import warnings
import remix.framework
import importlib.metadata as md
from pandas.errors import PerformanceWarning
from pathlib import Path
from remix.framework.api.instance import Instance


# ------------------------------------------------------------------------------------
# USER CHOICES
# ------------------------------------------------------------------------------------

group_name = "GP-NT-ELEC-BIO-H2"
base_scenario = "GP"
yrs_sel = [2020, 2050]

# Switch blocks on/off while debugging the full build
include_renewables = True
include_conventional_pp = True
include_lithium_batteries = True
include_h2_techs = True
include_heat_transp = True
include_emissions_constraints = True

# ------------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------------

yrs_to_calc = list(yrs_sel)
years_tag = "-".join(str(y) for y in yrs_sel)
case_name = f"nz_case_{base_scenario}_{years_tag}"

path_input = Path("C:/Local/REMix/remix_nz/input")
path_demand = path_input / "demand" / group_name
path_profiles = path_input / "profiles"
path_brownfield = path_input / "brownfield"

hourly_electricity_file = "hourly_electricity_GP-NT-ELEC-BIO-H2.csv"

base_dir = Path(f"../project/{group_name}/{case_name}")
data_dir = base_dir / "data"
results_dir = base_dir / "result"

# print("\n--- SCENARIO SETUP COMPLETE ---")
# print("  case:", case_name)
# print("  scenario:", base_scenario)
# print("  years (optimise):", yrs_sel)
# print("  years (data):", yrs_to_calc)
# print("  demand folder:", path_demand.resolve())

warnings.filterwarnings("ignore", category=PerformanceWarning) # Silence pandas performance warning from wide hourly tables
idx = pd.IndexSlice # commonly used for multiindex slicing
pd.set_option("display.float_format", lambda x: f"{x:,.4f}")
np.set_printoptions(precision=4, suppress=True)


# ------------------------------------------------------------------------------------
# UTILITIES
# ------------------------------------------------------------------------------------

def map_region_to_remix(region_value: str) -> str:
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

# Shared fuel-cost constants (used by thermal + fuel imports)
eur_per_aud = 0.62
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

# ------------------------------------------------------------------------------------
# Core functions
# ------------------------------------------------------------------------------------

def add_scope(m):
    """
    Scope: node aggregation map + core sets.
    """
    print("\n--- ADDING SCOPE ---")
    nodes = ["AKL","BOP","CAN","CEN","HBY","NEL","NIS","OTG","TRN","WEL","WTO"]

    # map_aggregatenodesmodel (schema keys: nodesData, nodesModel)
    mp = pd.DataFrame(index=pd.MultiIndex.from_product([nodes, nodes], names=["nodesData", "nodesModel"]))

    # sets_nodesdata / set_nodesmodel are derived from the map
    m["Base"].set.add(nodes, "nodesdata")
    m["Base"].set.add(nodes, "nodesmodel")

    # years + yearssel
    m["Base"].set.add([str(y) for y in sorted(set(yrs_to_calc))], "years")
    m["Base"].set.add([str(y) for y in sorted(set(yrs_sel))], "yearssel")

    print("[add_scope] nodesdata:", nodes)



def add_demand(m):
    """
    Electricity demand + slack, using the v0.14 label:
      nodesData, years, sourcesinktechs, commodity, profileTypes
    """
    print("\n--- ADDING ELECTRICITY DEMAND + SLACK ---")

    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    years_sel = [str(y) for y in years_sel]

    elec_path = Path(path_demand) / hourly_electricity_file
    if not elec_path.exists():
        raise FileNotFoundError(f"Hourly electricity file not found: {elec_path}")

    df = pd.read_csv(elec_path)

    # Normalise columns
    rename_cols = {}
    if "scenario" in df.columns and "Scenario" not in df.columns:
        rename_cols["scenario"] = "Scenario"
    if "Region" in df.columns:
        rename_cols["Region"] = "node"
    if "region" in df.columns:
        rename_cols["region"] = "node"
    if "Year" in df.columns:
        rename_cols["Year"] = "year"
    if "Carrier" in df.columns:
        rename_cols["Carrier"] = "carrier"
    if rename_cols:
        df = df.rename(columns=rename_cols)

    required = {"Scenario", "node", "year", "carrier"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found: {list(df.columns)}")

    # Map nodes
    df["node"] = df["node"].apply(
        lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r).strip()
    )
    df = df.dropna(subset=["node"]).copy()

    # Scenario filter
    df = df.loc[df["Scenario"].astype(str).str.strip() == str(base_scenario)].copy()

    # Electricity only
    df["carrier"] = df["carrier"].astype(str).str.strip().str.lower()
    df = df.loc[df["carrier"] == "electricity"].copy()

    # Years
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int).astype(str)
    df = df.loc[df["year"].isin(years_sel)].copy()

    # Hours
    hour_cols = [c for c in df.columns if str(c).startswith("t")]
    if not hour_cols:
        raise ValueError("No hourly columns found (expected t0001..t8760).")
    df[hour_cols] = df[hour_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Aggregate node-year
    df_agg = df.groupby(["node", "year"], as_index=False)[hour_cols].sum()
    print("[add_demand] rows:", len(df_agg), "unique nodes:", df_agg["node"].nunique())

    # ===== Demand profile (negative) =====
    ts = df_agg.set_index(["node", "year"])[hour_cols].copy()
    ts *= -1.0

    ts["sourcesinktechs"] = "Elec"
    ts["commodity"] = "Elec"
    ts["profileTypes"] = "fixed"

    ts = ts.set_index(["sourcesinktechs", "commodity", "profileTypes"], append=True)
    ts.index = ts.index.set_names(["nodesData", "years", "sourcesinktechs", "commodity", "profileTypes"])
    ts = ts.round(6)

    m["Base"].profile.add(ts, "sourcesink_profile")

    # Config: fixed profile
    cfg = pd.DataFrame(index=ts.index.droplevel("profileTypes").unique())
    cfg.index = cfg.index.set_names(["nodesData", "years", "sourcesinktechs", "commodity"])
    cfg["usesFixedProfile"] = 1.0
    m["Base"].parameter.add(cfg, "sourcesink_config")

    # ===== Slack (positive supply, bounded by annualsum) =====
    slack_idx = pd.MultiIndex.from_product(
        [sorted(m["Base"].set.nodesdata), years_sel, ["Slack"], ["Elec"]],
        names=["nodesData", "years", "sourcesinktechs", "commodity"],
    )

    slack_cfg = pd.DataFrame(index=slack_idx)
    slack_cfg["usesLowerSum"] = 1.0
    slack_cfg["usesUpperSum"] = 1.0
    m["Base"].parameter.add(slack_cfg, "sourcesink_config")

    slack_sum = pd.DataFrame(index=slack_idx)
    slack_sum["lower"] = 0.0
    slack_sum["upper"] = np.inf
    m["Base"].parameter.add(slack_sum, "sourcesink_annualsum")

    # Slack cost (use horizon)
    acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["SlackCost"], ["global"], ["horizon"], ["Slack"], ["Elec"]],
            names=["indicator", "accNodesData", "accYears", "sourcesinktechs", "commodity"],
        )
    )
    acc["perFlow"] = 1e6
    m["Base"].parameter.add(acc.fillna(0.0), "accounting_sourcesinkflow")

    print("[add_demand] wrote demand+slack for nodes:", len(sorted(m["Base"].set.nodesdata)), "years:", years_sel)

    # --- Sanity check: Slack must have lower=0 and usesLowerSum=1 everywhere ---
    ss = m["Base"].parameter.sourcesink_annualsum
    cfg = m["Base"].parameter.sourcesink_config

    sl = ss.loc[idx[:, years_sel, "Slack", "Elec"], :]
    if (sl["lower"] < -1e-12).any() or (sl["lower"] > 1e-12).any():
        raise ValueError("[add_demand] Slack lower bound is not exactly 0 everywhere.")

    sc = cfg.loc[idx[:, years_sel, "Slack", "Elec"], :]
    if not (sc["usesLowerSum"] == 1).all():
        raise ValueError("[add_demand] Slack usesLowerSum not 1 everywhere.")
    
    # Add hourly lower profile so Slack cannot go negative in any hour
    hour_cols = [c for c in df.columns if str(c).startswith("t")]  # use demand DF hour columns

    slack_idx = pd.MultiIndex.from_product(
        [sorted(m["Base"].set.nodesdata), [str(y) for y in m["Base"].set.yearssel], ["Slack"], ["Elec"]],
        names=["nodesData", "years", "sourcesink_techs", "commodity"],
    )

    # Turn on lower profile constraint
    slack_cfg = pd.DataFrame(index=slack_idx)
    slack_cfg["usesLowerProfile"] = 1.0
    m["Base"].parameter.add(slack_cfg, "sourcesink_config")

    # Lower profile = 0 for all hours
    slack_lower = pd.DataFrame(index=slack_idx, columns=hour_cols, data=0.0)
    slack_lower["profileTypes"] = "lower"
    slack_lower = slack_lower.set_index("profileTypes", append=True)
    slack_lower.index = slack_lower.index.set_names(
        ["nodesData", "years", "sourcesink_techs", "commodity", "profileTypes"]
    )
    m["Base"].profile.add(slack_lower, "sourcesink_profile")



def add_network(m):
    """
    Add transmission network + accounting.

    Uses v0.14 schema keys:
      * transfer_linkstartend: (linksData, nodesData) -> start/end
      * transfer_lengthparam: (linksData, link_types) -> length
      * transfer_linksparam: (linksData, years, transfer_techs) -> linksUpperLimit, etc.
      * transfer_techparam: (transfer_techs, vintage) -> lifeTime, flowUpperLimit, etc.
      * transfer_coefficient: (transfer_techs, vintage, commodity) -> coefficient
      * transfer_coefperflow: (transfer_techs, vintage, commodity) -> coefPerFlow
      * transfer_coefperlength: (transfer_techs, vintage, commodity, link_types) -> coefPerLength
      * accounting_transferlinks: (indicator, accLinksData, accYears, transfer_techs, vintage) -> perLinkBuild, ...
      * accounting_transferperlength: (indicator, accLinksData, accYears, transfer_techs, vintage, link_types) -> perLengthBuild, ...
    - Adapts automatically to any yearssel: it builds vintages from m["Base"].set.yearssel.
    """
    print("\n--- ADDING TRANSMISSION NETWORK + ACCOUNTING ---")
    # ------------------------------------------------------------------
    # Define links and types
    # ------------------------------------------------------------------
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

    # Sets (v0.14 names: linksdata/linksmodel/link_types)
    m["Base"].set.add(link_names, "linksdata")
    m["Base"].set.add(link_names, "linksmodel")

    link_types = ["land", "sea"]
    m["Base"].set.add(link_types, "link_types")

    # Use the model's nodesdata set (already added in add_scope)
    nodes_data = sorted(set(m["Base"].set.nodesdata))

    # ------------------------------------------------------------------
    # transfer_linkstartend  (linksData x nodesData -> start/end flags)
    # ------------------------------------------------------------------
    link_connections = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [link_names, nodes_data],
            names=["linksData", "nodesData"],
        )
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

    link_connections = link_connections.fillna(0.0)
    m["Base"].parameter.add(link_connections, "transfer_linkstartend")

    # ------------------------------------------------------------------
    # transfer_lengthparam  (linksData x link_types -> km)
    # ------------------------------------------------------------------
    link_lengths = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [link_names, link_types],
            names=["linksData", "link_types"],
        )
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
    link_lengths = link_lengths.fillna(0.0)
    m["Base"].parameter.add(link_lengths, "transfer_lengthparam")

    # ------------------------------------------------------------------
    # Transfer tech / caps / coefficients / losses
    # ------------------------------------------------------------------
    transport_techs = ["HV"]
    m["Base"].set.add(transport_techs, "transfer_techs")

    # Important: in v0.14 schemas, years are strings (Instance casts to str anyway).
    years_model = sorted(str(y) for y in m["Base"].set.yearssel)

    # Commodity name must match set_commodities
    commodities = ["Elec"]

    # transfer_linksparam (linksData x years x transfer_techs)
    link_caps = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [link_names, years_model, transport_techs],
            names=["linksData", "years", "transfer_techs"],
        )
    )
    link_caps.loc[idx[:, :, "HV"], "linksUpperLimit"] = 100.0
    link_caps = link_caps.fillna(0.0)
    m["Base"].parameter.add(link_caps, "transfer_linksparam")

    # transfer_techparam (transfer_techs x vintage)
    # Keep your original approach: one vintage per model year (adaptive to yearssel).
    tech_params = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, years_model],
            names=["transfer_techs", "vintage"],
        )
    )
    tech_params["lifeTime"] = 40
    tech_params["flowUpperLimit"] = 1
    m["Base"].parameter.add(tech_params, "transfer_techparam")

    # transfer_coefficient (transfer_techs x vintage x commodity)
    transfer_coefficient = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, years_model, commodities],
            names=["transfer_techs", "vintage", "commodity"],
        )
    )
    transfer_coefficient["coefficient"] = 1
    m["Base"].parameter.add(transfer_coefficient, "transfer_coefficient")

    # transfer_coefperflow (transfer_techs x vintage x commodity)
    coef_per_flow = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, years_model, commodities],
            names=["transfer_techs", "vintage", "commodity"],
        )
    )
    coef_per_flow["coefPerFlow"] = -0.014
    m["Base"].parameter.add(coef_per_flow, "transfer_coefperflow")

    # transfer_coefperlength (transfer_techs x vintage x commodity x link_types)
    coef_per_dist = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [transport_techs, years_model, commodities, link_types],
            names=["transfer_techs", "vintage", "commodity", "link_types"],
        )
    )
    coef_per_dist.loc[idx[:, :, :, "land"], "coefPerLength"] = -0.00004
    coef_per_dist.loc[idx[:, :, :, "sea"], "coefPerLength"] = -0.00003
    coef_per_dist = coef_per_dist.fillna(0.0)
    m["Base"].parameter.add(coef_per_dist, "transfer_coefperlength")

    # ------------------------------------------------------------------
    # Accounting for transfer (links and per length)
    # Use accYears = ["horizon"] for constant costs across modelling years.
    # ------------------------------------------------------------------
    cost_indicators = ["Invest", "OMFix"]
    accLinksData = ["global"]
    accYears_vals = ["horizon"]
    vintages = list(years_model)

    acc_tl = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, accLinksData, accYears_vals, transport_techs, vintages],
            names=["indicator", "accLinksData", "accYears", "transfer_techs", "vintage"],
        )
    )
    acc_tl.loc[idx["Invest", "global", :, :, :], "perLinkBuild"] = 180
    acc_tl.loc[idx["Invest", "global", :, :, :], "interest"] = 0.06
    acc_tl.loc[idx["Invest", "global", :, :, :], "amorTime"] = 40
    acc_tl.loc[idx["Invest", "global", :, :, :], "useAnnuity"] = 1
    acc_tl.loc[idx["OMFix", "global", :, :, :], "perLinkTotal"] = 1.8
    m["Base"].parameter.add(acc_tl.fillna(0.0), "accounting_transferlinks")

    acc_tpl = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [cost_indicators, accLinksData, accYears_vals, transport_techs, vintages, link_types],
            names=["indicator", "accLinksData", "accYears", "transfer_techs", "vintage", "link_types"],
        )
    )

    # land
    acc_tpl.loc[idx["Invest", "global", :, :, :, "land"], "perLengthBuild"] = 0.544
    acc_tpl.loc[idx["Invest", "global", :, :, :, "land"], "interest"] = 0.06
    acc_tpl.loc[idx["Invest", "global", :, :, :, "land"], "amorTime"] = 40
    acc_tpl.loc[idx["Invest", "global", :, :, :, "land"], "useAnnuity"] = 1
    acc_tpl.loc[idx["OMFix", "global", :, :, :, "land"], "perLengthTotal"] = 0.00544

    # sea
    acc_tpl.loc[idx["Invest", "global", :, :, :, "sea"], "perLengthBuild"] = 0.975
    acc_tpl.loc[idx["Invest", "global", :, :, :, "sea"], "interest"] = 0.06
    acc_tpl.loc[idx["Invest", "global", :, :, :, "sea"], "amorTime"] = 40
    acc_tpl.loc[idx["Invest", "global", :, :, :, "sea"], "useAnnuity"] = 1
    acc_tpl.loc[idx["OMFix", "global", :, :, :, "sea"], "perLengthTotal"] = 0.00975

    m["Base"].parameter.add(acc_tpl.fillna(0.0), "accounting_transferperlength")

    print(
        "[add_network] links:", len(link_names),
        "transfer techs:", transport_techs,
        "vintages:", vintages,
        "accYears:", accYears_vals,
    )

def add_accounting(m):
    #  The value global uses all the regions in the system
    #  Horizon takes into account all years in the set set.yearssel 
    bnd = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["global"], ["horizon"], ["SystemCost"]],
            names=["accNodesData","accYears","indicator"],
        )
    )
    bnd["obj"] = -1 # minimization of system costs
    bnd["discount"] = 0.02 # discount rate for the indicators
    m["Base"].parameter.add(bnd, "accounting_indicatorbounds")

    per = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["SystemCost"], ["Invest", "OMFix", "FuelCost", "ImportCost", "SlackCost"], ["global"], ["horizon"]],
            names=["indicator","indicator_a","accNodesData","accYears"],
        )
    )
    per["perIndicator"] = 1
    m["Base"].parameter.add(per, "accounting_perindicator")

def validate_scope(m):
    """
    Scope validation: checks that years, yearssel, nodesdata, nodesmodel match expectations.
    """
    years = sorted(int(y) for y in m["Base"].set.years)
    ysel = sorted(int(y) for y in m["Base"].set.yearssel)

    missing = [y for y in ysel if y not in years]
    if missing:
        raise ValueError(f"yearssel contains years not in years: {missing}")

    expected_sel = sorted(int(y) for y in yrs_sel)
    if ysel != expected_sel:
        raise ValueError(f"Model yearssel {ysel} does not match yrs_sel {expected_sel}")

    expected_years = sorted(set(int(y) for y in yrs_to_calc))
    if sorted(set(years)) != expected_years:
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
    print("  nodesdata:", len(nodesdata), nodesdata)
    print("  nodesmodel:", len(nodesmodel), nodesmodel)
    print("  years (data):", years)
    print("  yearssel (optimise):", ysel)

# conventional

def add_thermal(m, add_fuel_imports=True, debug=False):
    """
    Thermal brownfield:
    - unitsBuild at years = actual Year_built (extend yrs_to_calc + set years)
    - techno-economic parameters on vintage axis (1950/2020/2030/2050)
    - adds fuel import sourcesinks so plants can dispatch
    - noExpansion in base year only (2020); expansion allowed in later years (e.g. 2050)
    """
    print("\n--- ADDING THERMAL (build-year capacity + vintage parameters) ---")

    nodes = sorted(list(m["Base"].set.nodesdata))
    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    base_year = int(years_sel[0])

    # Core naming: must match your demand commodity
    elec = "Elec"
    activity = "Powergen"

    vintages = [1950, 2020, 2030, 2050]
    vintages_str = [str(v) for v in vintages]
    dea_year_for_vintage = {1950: 2015, 2020: 2020, 2030: 2030, 2050: 2050}

    def map_build_to_vintage(y):
        mode = str(globals().get("BROWNFIELDVINTAGEMODE", "A")).upper()
        if mode == "B":
            return 2020
        if y < 2020:
            return 1950
        elif y < 2030:
            return 2020
        elif y < 2050:
            return 2030
        else:
            return 2050

    inst = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")

    tech_map = {
        "Coal":       ("Thermal_Coal",   "Coal"),
        "Diesel":     ("Thermal_Diesel", "Diesel"),
        "Biogas":     ("Thermal_Bio",    "Biofuel"),
        "Biomass":    ("Thermal_Bio",    "Biofuel"),
        "Wood":       ("Thermal_Bio",    "Biofuel"),
        "Wood waste": ("Thermal_Bio",    "Biofuel"),
    }

    df = inst[(inst["Type"] == "Thermal") & (inst["Primary_fuel"].isin(tech_map.keys()))].copy()
    if df.empty:
        print("Thermal: no matching rows in brownfield DB.")
        return

    df["convertertechs"] = df["Primary_fuel"].map(lambda f: tech_map[f][0])
    df["fuel"] = df["Primary_fuel"].map(lambda f: tech_map[f][1])

    df["Year_built"] = pd.to_numeric(df["Year_built"], errors="coerce")
    df["Capacity_MW"] = pd.to_numeric(df["Capacity_MW"], errors="coerce")
    df = df.dropna(subset=["Year_built", "Capacity_MW"]).copy()
    df["Year_built"] = df["Year_built"].astype(int)

    df["nodesData"] = df["Node"].astype(str).str.strip()
    df = df[df["nodesData"].isin(nodes)].copy()
    df = df[df["Capacity_MW"] > 0].copy()
    if df.empty:
        print("Thermal: no plants with positive capacity in model nodes.")
        return

    df["vintage_bucket"] = df["Year_built"].apply(map_build_to_vintage).astype(int)

    techs = sorted(df["convertertechs"].unique().tolist())
    fuels = sorted(df["fuel"].unique().tolist())
    fuel_by_tech = {v[0]: v[1] for v in tech_map.values()}

    if debug:
        print("[add_thermal] Year_built min/max:", int(df["Year_built"].min()), int(df["Year_built"].max()))
        uniq = sorted(df["Year_built"].unique().tolist())
        print("[add_thermal] unique Year_built (first 15):", uniq[:15])
        print("[add_thermal] unique Year_built (last 15):", uniq[-15:])
        print("[add_thermal] MW mapped into vintage buckets:")
        print(df.groupby(["convertertechs", "vintage_bucket"])["Capacity_MW"].sum())

    # Extend yrs_to_calc + set years
    build_years = sorted(set(int(y) for y in df["Year_built"].unique().tolist()))
    new_years = sorted(set(int(y) for y in yrs_to_calc).union(build_years))
    if new_years != sorted(set(int(y) for y in yrs_to_calc)):
        yrs_to_calc[:] = new_years
        m["Base"].set.add([int(y) for y in new_years], "years")
        if debug:
            print("[add_thermal] extended yrs_to_calc with build years, years now:",
                  sorted(int(y) for y in m["Base"].set.years))

    # --- DEA parameter tables (as in your file) ---
    dea_eff = {
        "Thermal_Coal":   {2015: 0.46,  2020: 0.485, 2030: 0.52,  2050: 0.535},
        "Thermal_Diesel": {2015: 0.37,  2020: 0.37,  2030: 0.37,  2050: 0.37},
        "Thermal_Bio":    {2015: 0.42,  2020: 1.25,  2030: 0.45,  2050: 1.37},
    }
    dea_life = {
        "Thermal_Coal":   {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "Thermal_Diesel": {2015: 25, 2020: 25, 2030: 25, 2050: 25},
        "Thermal_Bio":    {2015: 25, 2020: 50, 2030: 25, 2050: 50},
    }
    dea_capex = {
        "Thermal_Coal":   {2015: 2182.4,      2020: 2148.4,      2030: 2103.2,      2050: 2012.8},
        "Thermal_Diesel": {2015: 372.180451,  2020: 361.546724,  2030: 361.546724,  2050: 361.546724},
        "Thermal_Bio":    {2015: 1063.4,      2020: 3136.9,      2030: 957.0,       2050: 3030.6},
    }
    dea_omfix_frac = {
        "Thermal_Coal":   {2015: 0.016,        2020: 0.016,        2030: 0.016,        2050: 0.016},
        "Thermal_Diesel": {2015: 0.025142857,  2020: 0.025882353,  2030: 0.024847059,  2050: 0.023811765},
        "Thermal_Bio":    {2015: 0.010,        2020: 0.012,        2030: 0.010,        2050: 0.010},
    }
    ef_kg_per_mwh_fuel = {"Thermal_Coal": 354.6, "Thermal_Diesel": 262.4, "Thermal_Bio": 0.0}
    fuel_cost_meur_per_gwh_fuel = {"Coal": 0.03, "Diesel": 0.12, "Biofuel": 0.06}

    def eff(t, v):
        return float(dea_eff[t][dea_year_for_vintage[int(v)]])

    def lifetime(t, v):
        return int(dea_life[t][dea_year_for_vintage[int(v)]])

    def nearest_data_year(y):
        candidates = [2015, 2020, 2030, 2050]
        return min(candidates, key=lambda d: abs(d - int(y)))

    # --- converter_techparam (vintage axis) ---
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, vintages_str], names=["convertertechs", "vintage"])
    )
    techparam["activityUpperLimit"] = 1.0
    for t in techs:
        for v in vintages:
            techparam.loc[(t, str(v)), "lifeTime"] = lifetime(t, v)
    m["Base"].parameter.add(techparam, "converter_techparam")

    # --- converter_capacityparam ---
    # Brownfield builds at actual build year
    cap_build = (
        df.groupby(["nodesData", "Year_built", "convertertechs"])["Capacity_MW"]
          .sum()
          .rename("unitsBuild")
          .to_frame()
          .div(1e3)
    )
    cap_build.index = cap_build.index.set_names(["nodesData", "years", "convertertechs"])
    cap_build.index = cap_build.index.set_levels(
        [cap_build.index.levels[0],
         cap_build.index.levels[1].astype(str),
         cap_build.index.levels[2]],
        level=[0, 1, 2]
    )

    # Bounds for ALL years in set years (not only yearssel)
    years_data_str = sorted(str(y) for y in m["Base"].set.years)
    cap_bounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_data_str, techs], names=["nodesData", "years", "convertertechs"])
    )
    cap_bounds["unitsUpperLimit"] = 100.0
    cap_bounds["noExpansion"] = 0.0
    cap_bounds.loc[idx[:, [str(base_year)], :], "noExpansion"] = 1.0  # block only in 2020

    cap_full = pd.concat([cap_build, cap_bounds], axis=1).sort_index()
    m["Base"].parameter.add(cap_full, "converter_capacityparam")

    # --- converter_coefficient (vintage axis; ensure Elec produced + fuel consumed) ---
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages_str, [activity], [elec] + fuels],
            names=["convertertechs", "vintage", "activity", "commodity"],
        )
    )
    coef["coefficient"] = 0.0
    coef.loc[idx[:, :, activity, elec], "coefficient"] = 1.0
    for t in techs:
        f = fuel_by_tech[t]
        for v in vintages:
            coef.loc[(t, str(v), activity, f), "coefficient"] = -1.0 / eff(t, v)
    m["Base"].parameter.add(coef, "converter_coefficient")

    # --- accounting tables (vintage axis, accYears="horizon") ---
    y_cost = base_year
    acc_units = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, vintages_str],
            names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage"],
        )
    ).sort_index()

    for t in techs:
        for v in vintages:
            dy = nearest_data_year(y_cost)
            capex = float(dea_capex[t][dy])
            omfix = capex * float(dea_omfix_frac[t][dy])
            acc_units.loc[("Invest", "global", "horizon", t, str(v)), "perUnitBuild"] = capex
            acc_units.loc[("Invest", "global", "horizon", t, str(v)), "useAnnuity"] = 1
            acc_units.loc[("Invest", "global", "horizon", t, str(v)), "amorTime"] = lifetime(t, 2020)
            acc_units.loc[("Invest", "global", "horizon", t, str(v)), "interest"] = 0.05
            acc_units.loc[("OMFix", "global", "horizon", t, str(v)), "perUnitTotal"] = omfix

    m["Base"].parameter.add(acc_units.fillna(0.0), "accounting_converterunits")

    acc_act = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["FuelCost", "CO2_emission"], ["global"], ["horizon"], techs, vintages_str, [activity]],
            names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage", "activity"],
        )
    ).sort_index()
    acc_act["perActivity"] = 0.0

    for t in techs:
        f = fuel_by_tech[t]
        for v in vintages:
            e = eff(t, v)
            acc_act.loc[("FuelCost", "global", "horizon", t, str(v), activity), "perActivity"] = float(fuel_cost_meur_per_gwh_fuel[f]) / e
            acc_act.loc[("CO2_emission", "global", "horizon", t, str(v), activity), "perActivity"] = (float(ef_kg_per_mwh_fuel[t]) * 1e-3) / e

    m["Base"].parameter.add(acc_act.fillna(0.0), "accounting_converteractivity")

    # --- Fuel imports sourcesinks (critical so thermal can run) ---
    if add_fuel_imports and fuels:
        # allow imports only in optimisation years (like your existing pattern)
        for f in fuels:
            ss_name = f"FuelImport_{f}"
            ss_idx = pd.MultiIndex.from_product(
                [nodes, [str(y) for y in years_sel], [ss_name], [f]],
                names=["nodesData", "years", "sourcesinks", "commodities"],
            )
            ss = pd.DataFrame(index=ss_idx)
            ss["lower"] = 0.0
            ss["upper"] = 1e12
            m["Base"].parameter.add(ss, "sourcesink_annualsum")

            cfg = pd.DataFrame(index=ss_idx)
            cfg["usesLowerSum"] = 1
            cfg["usesUpperSum"] = 1
            m["Base"].parameter.add(cfg, "sourcesink_config")

            acc_idx = pd.MultiIndex.from_product(
                [["FuelCost"], ["global"], ["horizon"], [ss_name], [f]],
                names=["indicator", "accNodesData", "accYears", "sourcesinktechs", "commodity"],
            )
            acc = pd.DataFrame(index=acc_idx)
            acc["perFlow"] = float(fuel_cost_meur_per_gwh_fuel[f])
            m["Base"].parameter.add(acc.fillna(0.0), "accounting_sourcesinkflow")

        print("[add_thermal] fuel imports added:", [f"FuelImport_{f}" for f in fuels])

    print("[add_thermal] Done.")

def add_gas_turbines(m, add_fuel_import=True, debug=False):
    """
    Gas turbines:
    - unitsBuild at years = actual Year_built
    - parameters on vintage axis (1950/2020/2030/2050)
    - adds FuelImport_CH4 sourcesink so turbines can dispatch
    - blocks new builds only in base year (2020)
    """
    print("\n--- ADDING GAS TURBINES (build-year capacity + vintage parameters) ---")

    nodes = sorted([n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")])

    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    base_year = int(years_sel[0])

    elec = "Elec"
    fuel = "CH4"
    activity = "Powergen"

    techs_allowed = ["GT", "CCGT", "OCGT"]
    vintages = [1950, 2020, 2030, 2050]
    vintages_str = [str(v) for v in vintages]
    dea_year_for_vintage = {1950: 2015, 2020: 2020, 2030: 2030, 2050: 2050}

    def map_build_to_vintage(y):
        mode = str(globals().get("BROWNFIELDVINTAGEMODE", "A")).upper()
        if mode == "B":
            return 2020
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
        print("Gas turbines: no matching fuel rows.")
        return

    df["Year_built"] = pd.to_numeric(df["Year_built"], errors="coerce")
    df["Capacity_MW"] = pd.to_numeric(df["Capacity_MW"], errors="coerce").fillna(0.0)
    df = df.dropna(subset=["Year_built"]).copy()
    df["Year_built"] = df["Year_built"].astype(int)
    df = df[df["Capacity_MW"] > 0].copy()

    df["nodesData"] = df["Node"].astype(str).str.strip()
    df = df[df["nodesData"].isin(nodes)].copy()
    if df.empty:
        print("Gas turbines: no plants with positive capacity in model nodes.")
        return

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
        df["convertertechs"] = df["Techs"].apply(map_gt_type)
    else:
        df["convertertechs"] = "GT"

    df = df[df["convertertechs"].isin(techs_allowed)].copy()
    if df.empty:
        print("Gas turbines: no GT/CCGT/OCGT units after filters.")
        return

    df["vintage_bucket"] = df["Year_built"].apply(map_build_to_vintage).astype(int)

    techs = sorted(df["convertertechs"].unique().tolist())

    if debug:
        print("[add_gas_turbines] Year_built min/max:", int(df["Year_built"].min()), int(df["Year_built"].max()))
        uniq = sorted(df["Year_built"].unique().tolist())
        print("[add_gas_turbines] unique Year_built (first 15):", uniq[:15])
        print("[add_gas_turbines] unique Year_built (last 15):", uniq[-15:])
        print("[add_gas_turbines] MW mapped into vintage buckets:")
        print(df.groupby(["convertertechs", "vintage_bucket"])["Capacity_MW"].sum())

    # Extend yrs_to_calc + set years
    build_years = sorted(set(int(y) for y in df["Year_built"].unique().tolist()))
    new_years = sorted(set(int(y) for y in yrs_to_calc).union(build_years))
    if new_years != sorted(set(int(y) for y in yrs_to_calc)):
        yrs_to_calc[:] = new_years
        m["Base"].set.add([int(y) for y in new_years], "years")
        if debug:
            print("[add_gas_turbines] extended yrs_to_calc with build years, years now:",
                  sorted(int(y) for y in m["Base"].set.years))

    # DEA tables (as in your file)
    eff_table = {
        "GT":   {2015: 0.36, 2020: 0.37, 2030: 0.39, 2050: 0.40},
        "CCGT": {2015: 0.50, 2020: 0.51, 2030: 0.53, 2050: 0.55},
        "OCGT": {2015: 0.41, 2020: 0.42, 2030: 0.43, 2050: 0.45},
    }
    lifetime_table = {t: {2015: 25, 2020: 25, 2030: 25, 2050: 25} for t in techs_allowed}
    invest_meur_per_gw = {
        "GT":   {2015: 797.5,  2020: 776.3,  2030: 744.4,  2050: 723.1},
        "CCGT": {2015: 1382.4, 2020: 1382.4, 2030: 1276.0, 2050: 1169.7},
        "OCGT": {2015: 499.8,  2020: 478.5,  2030: 467.9,  2050: 436.0},
    }
    omfix_frac = {
        "GT":   {2015: 0.027, 2020: 0.027, 2030: 0.027, 2050: 0.026},
        "CCGT": {2015: 0.023, 2020: 0.023, 2030: 0.023, 2050: 0.024},
        "OCGT": {2015: 0.017, 2020: 0.018, 2030: 0.018, 2050: 0.018},
    }

    kgCO2_per_MWh_fuel = 204.8
    fuel_cost_meur_per_gwh_fuel = 0.045

    def eff(t, v):
        return float(eff_table[t][dea_year_for_vintage[int(v)]])

    def lifetime(t, v):
        return int(lifetime_table[t][dea_year_for_vintage[int(v)]])

    def nearest_data_year(y):
        candidates = [2015, 2020, 2030, 2050]
        return min(candidates, key=lambda d: abs(d - int(y)))

    # converter_techparam (vintage)
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, vintages_str], names=["convertertechs", "vintage"])
    )
    techparam["activityUpperLimit"] = 1.0
    for t in techs:
        for v in vintages:
            techparam.loc[(t, str(v)), "lifeTime"] = lifetime(t, v)
    m["Base"].parameter.add(techparam, "converter_techparam")

    # converter_capacityparam (build-year + bounds for all years)
    cap_build = (
        df.groupby(["nodesData", "Year_built", "convertertechs"])["Capacity_MW"]
          .sum()
          .rename("unitsBuild")
          .to_frame()
          .div(1e3)
    )
    cap_build.index = cap_build.index.set_names(["nodesData", "years", "convertertechs"])
    cap_build.index = cap_build.index.set_levels(
        [cap_build.index.levels[0],
         cap_build.index.levels[1].astype(str),
         cap_build.index.levels[2]],
        level=[0, 1, 2]
    )

    years_data_str = sorted(str(y) for y in m["Base"].set.years)
    cap_bounds = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_data_str, techs], names=["nodesData", "years", "convertertechs"])
    )
    cap_bounds["unitsUpperLimit"] = 100.0
    cap_bounds["noExpansion"] = 0.0
    cap_bounds.loc[idx[:, [str(base_year)], :], "noExpansion"] = 1.0

    cap_full = pd.concat([cap_build, cap_bounds], axis=1).sort_index()
    m["Base"].parameter.add(cap_full, "converter_capacityparam")

    # converter_coefficient (vintage; produce Elec, consume CH4)
    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, vintages_str, [activity], [elec, fuel]],
            names=["convertertechs", "vintage", "activity", "commodity"],
        )
    )
    coef["coefficient"] = 0.0
    coef.loc[idx[:, :, activity, elec], "coefficient"] = 1.0
    for t in techs:
        for v in vintages:
            coef.loc[(t, str(v), activity, fuel), "coefficient"] = -1.0 / eff(t, v)
    m["Base"].parameter.add(coef, "converter_coefficient")

    # accounting (vintage, horizon)
    y_cost = base_year
    acc_units = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], ["global"], ["horizon"], techs, vintages_str],
            names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage"],
        )
    ).sort_index()

    for t in techs:
        dy = nearest_data_year(y_cost)
        for v in vintages_str:
            capex = float(invest_meur_per_gw[t][dy])
            acc_units.loc[("Invest", "global", "horizon", t, v), "perUnitBuild"] = capex
            acc_units.loc[("Invest", "global", "horizon", t, v), "useAnnuity"] = 1
            acc_units.loc[("Invest", "global", "horizon", t, v), "amorTime"] = lifetime(t, 2020)
            acc_units.loc[("Invest", "global", "horizon", t, v), "interest"] = 0.05
            acc_units.loc[("OMFix", "global", "horizon", t, v), "perUnitTotal"] = capex * float(omfix_frac[t][dy])

    m["Base"].parameter.add(acc_units.fillna(0.0), "accounting_converterunits")

    acc_act = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["FuelCost", "CO2_emission"], ["global"], ["horizon"], techs, vintages_str, [activity]],
            names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage", "activity"],
        )
    ).sort_index()
    acc_act["perActivity"] = 0.0

    for t in techs:
        for v in vintages:
            e = eff(t, v)
            acc_act.loc[("FuelCost", "global", "horizon", t, str(v), activity), "perActivity"] = float(fuel_cost_meur_per_gwh_fuel) / e
            acc_act.loc[("CO2_emission", "global", "horizon", t, str(v), activity), "perActivity"] = (float(kgCO2_per_MWh_fuel) * 1e-3) / e

    m["Base"].parameter.add(acc_act.fillna(0.0), "accounting_converteractivity")

    # Fuel import for CH4 (critical)
    if add_fuel_import:
        ss_name = "FuelImport_CH4"
        ss_idx = pd.MultiIndex.from_product(
            [nodes, [str(y) for y in years_sel], [ss_name], [fuel]],
            names=["nodesData", "years", "sourcesinks", "commodities"],
        )
        ss = pd.DataFrame(index=ss_idx)
        ss["lower"] = 0.0
        ss["upper"] = 1e12
        m["Base"].parameter.add(ss, "sourcesink_annualsum")

        cfg = pd.DataFrame(index=ss_idx)
        cfg["usesLowerSum"] = 1
        cfg["usesUpperSum"] = 1
        m["Base"].parameter.add(cfg, "sourcesink_config")

        acc_idx = pd.MultiIndex.from_product(
            [["FuelCost"], ["global"], ["horizon"], [ss_name], [fuel]],
            names=["indicator", "accNodesData", "accYears", "sourcesinktechs", "commodity"],
        )
        acc = pd.DataFrame(index=acc_idx)
        acc["perFlow"] = float(fuel_cost_meur_per_gwh_fuel)
        m["Base"].parameter.add(acc.fillna(0.0), "accounting_sourcesinkflow")

        print("[add_gas_turbines] fuel import added: FuelImport_CH4")

    print("[add_gas_turbines] Done.")

# renewables and batteries

def add_hydro(m):
    """
    Hydro:
      - Loads Water_in inflows (fixed profile, positive supply)
      - Adds hydro turbine converter + reservoir storage
      - Brownfield build year = 2000
      - No expansion
      - Spill penalty (accounting_converteractivity, horizon)
    """
    global yrs_to_calc, inflow_file, idx

    print("\n--- ADDING HYDRO INFLOWS (Water_in) ---")

    nodes = sorted(list(m["Base"].set.nodesdata))
    yearssel = sorted(int(y) for y in m["Base"].set.yearssel)

    simplified_hydro_inflow_file = Path("C:/Local/REMix/remix_nz/input/brownfield/hydro/inflows_remix-nz/regional-hydro-inflow.csv")
    inflow_df = pd.read_csv(simplified_hydro_inflow_file)

    rename = {}
    if "Region" in inflow_df.columns:  rename["Region"] = "node"
    if "region" in inflow_df.columns:  rename["region"] = "node"
    if "Year" in inflow_df.columns:    rename["Year"] = "year"
    if "Sector" in inflow_df.columns:  rename["Sector"] = "sector"
    if "Carrier" in inflow_df.columns: rename["Carrier"] = "commodity"
    if "carrier" in inflow_df.columns: rename["carrier"] = "commodity"
    if rename:
        inflow_df = inflow_df.rename(columns=rename)

    required = {"node", "year"}
    missing = required - set(inflow_df.columns)
    if missing:
        raise ValueError(f"Hydro inflow CSV missing required columns {missing}. Found {list(inflow_df.columns)}")

    if "sector" not in inflow_df.columns:
        inflow_df["sector"] = "All"
    if "commodity" not in inflow_df.columns:
        inflow_df["commodity"] = "HydroInflow"

    inflow_df["node"] = inflow_df["node"].apply(lambda r: map_region_to_remix(r) if not str(r).isupper() else str(r))
    inflow_df = inflow_df.dropna(subset=["node"])
    inflow_df["node"] = inflow_df["node"].astype(str).str.strip()
    inflow_df = inflow_df[inflow_df["node"].isin(nodes)].copy()

    inflow_df["year"] = inflow_df["year"].astype(int)
    inflow_df = inflow_df.loc[inflow_df["year"].isin(yearssel)].copy()

    inflow_df["commodity"] = inflow_df["commodity"].astype(str).str.strip()
    inflow_df.loc[inflow_df["commodity"] == "HydroInflow", "commodity"] = "Water_in"

    hour_cols = [c for c in inflow_df.columns if str(c).startswith("t")]
    if not hour_cols:
        raise ValueError("No hourly columns found in hydro inflow file (t0001..t8760).")

    inflow = inflow_df.set_index(["node", "year", "sector", "commodity"])[hour_cols]
    inflow *= -1  # CSV negative -> positive supply
    inflow["type"] = "fixed"
    inflow_fixed = inflow.set_index("type", append=True).round(3)

    m["Base"].profile.add(inflow_fixed, "sourcesink_profile")

    inflow_cfg = pd.DataFrame(index=inflow.index)
    inflow_cfg["usesFixedProfile"] = 1
    inflow_cfg = inflow_cfg.loc[inflow.select_dtypes(include="number").sum(axis=1) != 0]
    m["Base"].parameter.add(inflow_cfg, "sourcesink_config")

    print(f"Hydro inflows loaded for years={sorted(inflow_df['year'].unique().tolist())}.")

    print("\n--- ADDING HYDROPOWER PLANTS (converter + reservoir) ---")

    # Ensure 2000 exists in years (data)
    yrs_to_calc[:] = sorted(set(int(y) for y in yrs_to_calc).union({2000}))
    m["Base"].set.add([int(y) for y in yrs_to_calc], "years")
    years_data = sorted(int(y) for y in m["Base"].set.years)

    hydro_nodes = ["BOP", "CAN", "CEN", "HBY", "NEL", "OTG", "WTO"]
    hydro_nodes = [n for n in hydro_nodes if n in nodes]

    hydro_techs = ["Hydro"]
    hydro_vintage = [1950]
    hydro_years = [2000] + years_data
    hydro_acts = ["Powergen", "Spill"]

    # Ocean sink for Water_out
    ocean_idx = pd.MultiIndex.from_product(
        [nodes, [str(y) for y in yearssel], ["Ocean"], ["Water_out"]],
        names=["nodesdata", "years", "sourcesinks", "commodities"],
    )
    ocean_cfg = pd.DataFrame(index=ocean_idx)
    ocean_cfg["usesUpperProfile"] = 1
    m["Base"].parameter.add(ocean_cfg, "sourcesink_config")

    # Turbine tech parameters
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([hydro_techs, [str(y) for y in hydro_vintage]], names=["techs", "years"])
    )
    conv_tech["lifeTime"] = 100
    conv_tech["activityUpperLimit"] = 1
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # Turbine capacity
    cap_idx = pd.MultiIndex.from_product(
        [hydro_nodes, [str(y) for y in hydro_years], hydro_techs],
        names=["nodesdata", "years", "techs"],
    )
    cap = pd.DataFrame(index=cap_idx)
    cap["noExpansion"] = 1
    cap["unitsUpperLimit"] = 100.0

    cap.loc[("BOP", "2000", "Hydro"), "unitsBuild"] = 0.17095
    cap.loc[("CAN", "2000", "Hydro"), "unitsBuild"] = 1.82683
    cap.loc[("CEN", "2000", "Hydro"), "unitsBuild"] = 0.399
    cap.loc[("HBY", "2000", "Hydro"), "unitsBuild"] = 0.1422
    cap.loc[("NEL", "2000", "Hydro"), "unitsBuild"] = 0.0453
    cap.loc[("OTG", "2000", "Hydro"), "unitsBuild"] = 1.664
    cap.loc[("WTO", "2000", "Hydro"), "unitsBuild"] = 1.0873

    m["Base"].parameter.add(cap.sort_index(), "converter_capacityparam")

    # Turbine coefficients
    hydro_eff = 0.95
    coef_idx = pd.MultiIndex.from_product(
        [hydro_techs, [str(y) for y in hydro_vintage], hydro_acts, ["Water_in", "Water_out", "Elec"]],
        names=["techs", "years", "activities", "commodities"],
    )
    coef = pd.DataFrame(index=coef_idx)
    coef.loc[idx[:, :, "Powergen", "Elec"], "coefficient"] = 1.0
    coef.loc[idx[:, :, "Powergen", "Water_in"], "coefficient"] = -hydro_eff
    coef.loc[idx[:, :, "Powergen", "Water_out"], "coefficient"] = hydro_eff
    coef.loc[idx[:, :, "Spill", "Water_in"], "coefficient"] = -100.0
    coef.loc[idx[:, :, "Spill", "Water_out"], "coefficient"] = 100.0
    m["Base"].parameter.add(coef, "converter_coefficient")

    # Reservoir storage
    stortechs = ["Hydro_reservoir"]
    stortech = pd.DataFrame(
        index=pd.MultiIndex.from_product([stortechs, [str(y) for y in hydro_vintage]], names=["techs", "years"])
    )
    stortech["lifeTime"] = 100
    stortech["levelUpperLimit"] = 1
    m["Base"].parameter.add(stortech, "storage_techparam")

    storsize = pd.DataFrame(
        index=pd.MultiIndex.from_product([stortechs, [str(y) for y in hydro_vintage], ["Water_in"]],
                                         names=["techs", "years", "commodities"])
    )
    storsize["size"] = 1
    storsize["selfdischarge"] = 0
    m["Base"].parameter.add(storsize, "storage_sizeparam")

    stor_idx = pd.MultiIndex.from_product(
        [hydro_nodes, [str(y) for y in hydro_years], stortechs],
        names=["nodesdata", "years", "storage_techs"],
    )
    stor = pd.DataFrame(index=stor_idx)
    stor["noExpansion"] = 1
    stor["unitsUpperLimit"] = 3000.0
    stor.loc[("CAN", "2000", "Hydro_reservoir"), "unitsBuild"] = 2517.2429
    stor.loc[("HBY", "2000", "Hydro_reservoir"), "unitsBuild"] = 154.2635
    stor.loc[("OTG", "2000", "Hydro_reservoir"), "unitsBuild"] = 729.5595
    stor.loc[("WTO", "2000", "Hydro_reservoir"), "unitsBuild"] = 587.1371
    m["Base"].parameter.add(stor.sort_index(), "storage_reservoirparam")

    # Spill penalty (horizon)
    penalty_per_gwh = 0.01
    spill_idx = pd.MultiIndex.from_product(
        [["SpillPenalty"], ["global"], ["horizon"], hydro_techs, [str(y) for y in hydro_vintage], ["Spill"]],
        names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage", "activity"],
    )
    spill = pd.DataFrame(index=spill_idx)
    spill["perActivity"] = float(penalty_per_gwh)
    m["Base"].parameter.add(spill.fillna(0.0), "accounting_converteractivity")

    print(f"Hydro spill penalty added: {penalty_per_gwh} M€/GWh activity.")
    print("Hydro inflows + plant added.")

def add_geothermal(m):
    print("\n--- ADDING GEOTHERMAL ---")

    global yrs_to_calc
    df = pd.read_csv(Path(path_brownfield) / "power-plant-nz-database.csv")

    # Ensure 2000 exists in data years
    yrs_to_calc[:] = sorted(set(int(y) for y in yrs_to_calc).union({2000}))
    m["Base"].set.add([int(y) for y in yrs_to_calc], "years")
    years_data = sorted(int(y) for y in m["Base"].set.years)

    nodes = sorted(list(m["Base"].set.nodesdata))
    tech = "Geothermal"
    activity = "Powergen"
    vintage_year = 2000

    # Tech param
    techparam = pd.DataFrame(
        index=pd.MultiIndex.from_product([[tech], [str(vintage_year)]], names=["techs", "years"])
    )
    techparam["lifeTime"] = 100
    techparam["activityUpperLimit"] = 0.98
    m["Base"].parameter.add(techparam, "converter_techparam")

    df_geo = df[df["Type"] == "Geothermal"].copy()
    if df_geo.empty:
        print("No geothermal rows found in brownfield database.")
        return

    df_geo["Capacity_MW"] = pd.to_numeric(df_geo["Capacity_MW"], errors="coerce").fillna(0.0)
    df_geo = df_geo[df_geo["Capacity_MW"] > 0].copy()

    cap_by_node = df_geo.groupby("Node")["Capacity_MW"].sum().div(1000.0)  # MW -> GW

    cap_idx = pd.MultiIndex.from_product(
        [nodes, [str(y) for y in years_data], [tech]],
        names=["nodesdata", "years", "techs"],
    )
    cap = pd.DataFrame(index=cap_idx)
    cap["noExpansion"] = 1
    cap["unitsUpperLimit"] = 0.0

    for n in nodes:
        gw = float(cap_by_node.get(n, 0.0))
        if gw <= 0:
            continue
        cap.loc[(n, str(vintage_year), tech), "unitsBuild"] = gw
        cap.loc[idx[n, :, tech], "unitsUpperLimit"] = gw

    cap = cap.dropna(how="all")
    m["Base"].parameter.add(cap.sort_index(), "converter_capacityparam")

    coef = pd.DataFrame(
        index=pd.MultiIndex.from_product([[tech], [str(vintage_year)], [activity], ["Elec"]],
                                         names=["techs", "years", "activities", "commodities"])
    )
    coef["coefficient"] = 1.0
    m["Base"].parameter.add(coef, "converter_coefficient")

    print("Geothermal installed as brownfield in 2000 and usable in all model years.")

def load_feedin_csv(weather_year: int = 2012):
    """
    Load renewables feed-in.

    - Solar: old timeseries_2012_w_corr.csv (MW), 30-min shift, wrap, 8760h.
    - Wind: results_w_corr/timeseries_norm_2012.csv (availability 0..1), no shift,
      mapped directly to t_model=1..8760 for each (node, tech).

    Returns DataFrame indexed by (nodesdata, techs, t_model) with column 'value'.
    """

    global path_profiles

    # -----------------------------
    # Solar: keep previous behaviour
    # -----------------------------
    ts_path_solar = Path(path_profiles) / "timeseries_2012_w_corr.csv"
    if not ts_path_solar.exists():
        raise FileNotFoundError(f"[load_feedin_csv] Solar timeseries file not found: {ts_path_solar}")

    ts_raw = pd.read_csv(ts_path_solar)

    required = {"t", "region", "technology", "timeseries_per_region"}
    missing = required - set(ts_raw.columns)
    if missing:
        raise ValueError(f"[load_feedin_csv] Timeseries missing columns {missing}. Found: {list(ts_raw.columns)}")

    dt = pd.to_datetime(ts_raw["t"], dayfirst=True, errors="coerce")
    if dt.isna().any():
        bad_sample = ts_raw.loc[dt.isna(), "t"].astype(str).head(5).tolist()
        raise ValueError(f"[load_feedin_csv] Solar timestamps could not be parsed, sample: {bad_sample}")

    # keep old +30 min behaviour for solar
    dt_shift = dt + pd.Timedelta(minutes=30)
    mask = (dt_shift.dt.year == int(weather_year))
    ts = ts_raw.loc[mask].copy()
    dt_shift = dt_shift.loc[mask]

    ts["value"] = pd.to_numeric(ts["timeseries_per_region"], errors="coerce")
    if ts["value"].isna().any():
        bad = ts.loc[ts["value"].isna()].head(5)
        raise ValueError(f"[load_feedin_csv] Non-numeric solar values, sample:\n{bad}")

    day_of_year = dt_shift.dt.dayofyear
    hour_of_day = dt_shift.dt.hour
    ts["t_model"] = ((day_of_year - 1) * 24 + hour_of_day + 1).astype(int)

    ts = ts.rename(columns={"region": "nodesdata", "technology": "techs"})
    ts = ts[["nodesdata", "techs", "t_model", "value"]]

    # collapse duplicates
    ts = ts.groupby(["nodesdata", "techs", "t_model"], as_index=False, sort=False)["value"].sum()

    all_hours = pd.Index(range(1, 8761), name="t_model")
    frames = []

    # reindex + wrap leading gap from end-of-year
    for (node, tech), g in ts.groupby(["nodesdata", "techs"], sort=False):
        s = g.set_index("t_model")["value"]
        s_full = s.reindex(all_hours, fill_value=0.0)

        observed_hours = sorted(g.index.tolist())
        first_obs = observed_hours[0] if observed_hours else 1

        wrapped = s_full.copy()
        if first_obs > 1:
            n_lead = first_obs - 1
            donor = s_full.iloc[-n_lead:].values
            if len(donor) == n_lead:
                wrapped.iloc[0:n_lead] = donor

        df_full = wrapped.to_frame("value")
        df_full["nodesdata"] = node
        df_full["techs"] = tech
        df_full = df_full.set_index(["nodesdata", "techs"], append=True)
        df_full = df_full.reorder_levels(["nodesdata", "techs", "t_model"])
        frames.append(df_full)

    feed_solar = pd.concat(frames).sort_index()

    # -----------------------------
    # Wind: new file, availability
    # -----------------------------
    base = Path(path_profiles) / "results_w_corr"
    ts_file = base / "timeseries_norm_2012.csv"
    if not ts_file.exists():
        raise FileNotFoundError(f"[load_feedin_csv] Wind timeseries file not found: {ts_file}")

    df = pd.read_csv(ts_file)
    df = df.rename(columns={c: c.strip() for c in df.columns})

    need = {"technology", "t", "region", "timeseries_norm"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"[load_feedin_csv] Wind file missing {missing}. Found: {list(df.columns)}")

    df["technology"] = df["technology"].astype(str).str.strip()
    df["region"] = df["region"].astype(str).str.strip()

    # timestamps may already be ISO-like; avoid double-parsing failures
    if np.issubdtype(df["t"].dtype, np.datetime64):
        t = df["t"]
    else:
        t_raw = df["t"].astype(str).str.replace("\u00A0", " ", regex=False).str.strip()
        # first try without dayfirst (works for ISO '2012-01-13 00:30:00')
        t = pd.to_datetime(t_raw, errors="coerce")
        bad = t.isna()
        if bad.any():
            # second pass with dayfirst for any remaining weird ones
            t2 = pd.to_datetime(t_raw[bad], dayfirst=True, errors="coerce")
            t.loc[bad] = t2

    if t.isna().any():
        print("[load_feedin_csv] Wind still has NaT timestamps (showing up to 10):")
        print(df.loc[t.isna(), ["technology", "region", "t"]].head(10).to_string(index=False))
        # replace with first valid time so we don't crash during build
        first_valid = t.dropna().iloc[0]
        t = t.fillna(first_valid)

    df["t"] = t

    df["timeseries_norm"] = pd.to_numeric(df["timeseries_norm"], errors="coerce")
    nanv = df["timeseries_norm"].isna()
    if nanv.any():
        print("[load_feedin_csv] NaN timeseries_norm count:", int(nanv.sum()), "-> filling with 0.0")
        df.loc[nanv, "timeseries_norm"] = 0.0

    df.loc[df["timeseries_norm"] < 0, "timeseries_norm"] = 0.0

    df = df[df["t"].dt.year == int(weather_year)].copy()
    df["nodesdata"] = df["region"]
    df["techs"] = df["technology"]

    df = df.sort_values("t")
    day_of_year = df["t"].dt.dayofyear
    hour_of_day = df["t"].dt.hour
    df["t_model"] = ((day_of_year - 1) * 24 + hour_of_day + 1).astype(int)

    # truncate per (node,tech) to 8760 and reindex
    all_hours = pd.Index(range(1, 8761), name="t_model")
    frames_w = []
    for (node, tech), g in df.groupby(["nodesdata", "techs"], sort=False):
        g = g.sort_values("t").iloc[:8760].copy()
        s = g.set_index("t_model")["timeseries_norm"]
        s_full = s.reindex(all_hours, fill_value=0.0)
        df_full = s_full.to_frame("value")
        df_full["nodesdata"] = node
        df_full["techs"] = tech
        df_full = df_full.set_index(["nodesdata", "techs"], append=True)
        df_full = df_full.reorder_levels(["nodesdata", "techs", "t_model"])
        frames_w.append(df_full)

    feed_wind = pd.concat(frames_w).sort_index()

    feed = pd.concat([feed_solar, feed_wind]).sort_index()
    return feed

def add_renewables(m):
    """
    Renewables block (PV + Wind variants).

    - Solar: installables from old region_statistics_2012_w_corr.csv (MW -> GW),
      profiles from timeseries_2012_w_corr.csv (MW), converted to availability.
    - Wind: installables from results_w_corr/installable_per_region.csv (GW),
      annual energy from annual_energy_per_region.csv, profiles from
      results_w_corr/timeseries_norm_2012.csv (availability 0..1).
    - Brownfield wind builds are allocated into wind_onshore_1, then _2, then _3, then _4.
    """
    import numpy as np
    import pandas as pd
    from pathlib import Path

    global path_profiles, idx

    print("\n--- ADDING RENEWABLES (PV + wind variants) ---")

    years_sel = sorted(int(y) for y in m["Base"].set.yearssel)
    base_year = years_sel[0]
    years_all = sorted(int(y) for y in m["Base"].set.years)
    years_all_str = [str(y) for y in years_all]

    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]

    pv_techs = ["pv_central_fixed", "pv_decentral"]
    wind_on = [f"wind_onshore_{i}" for i in [1, 2, 3, 4]]
    wind_off = [f"wind_offshore_{i}" for i in [1, 2, 3, 4]]
    wind_techs = wind_on + wind_off
    techs = pv_techs + wind_techs

    # ----------------------------
    # 1) Installables
    # ----------------------------
    # PV + wind: keep existing PV installables from old file (MW)
    inst_path_old = Path(path_profiles) / "region_statistics_2012_w_corr.csv"
    if not inst_path_old.exists():
        raise FileNotFoundError(f"Installables file not found: {inst_path_old}")

    inst_old = pd.read_csv(inst_path_old)
    inst_old = inst_old.rename(columns={c: c.strip() for c in inst_old.columns})

    need_old = {"region", "technology", "installable_per_region"}
    missing_old = need_old - set(inst_old.columns)
    if missing_old:
        raise ValueError(f"Old installables missing {missing_old}; found {list(inst_old.columns)}")

    inst_old["installable_per_region"] = pd.to_numeric(
        inst_old["installable_per_region"], errors="coerce"
    )
    if inst_old["installable_per_region"].isna().any():
        raise ValueError("Old installables has non-numeric installable_per_region.")

    inst_old["nodesdata"] = inst_old["region"].astype(str).str.strip()
    inst_old["techs"] = inst_old["technology"].astype(str).str.strip()
    inst_old = inst_old.loc[inst_old["nodesdata"].isin(nodes)].copy()

    # wind: override with new GW-based installables
    base_new = Path(path_profiles) / "results_w_corr"
    inst_new = pd.read_csv(base_new / "installable_per_region.csv")
    inst_new = inst_new.rename(columns={c: c.strip() for c in inst_new.columns})

    need_new = {"technology", "region", "installable_per_region"}
    missing_new = need_new - set(inst_new.columns)
    if missing_new:
        raise ValueError(f"New installables missing {missing_new}; found {list(inst_new.columns)}")

    inst_new["installable_per_region"] = pd.to_numeric(
        inst_new["installable_per_region"], errors="coerce"
    ).fillna(0.0)
    inst_new["nodesdata"] = inst_new["region"].astype(str).str.strip()
    inst_new["techs"] = inst_new["technology"].astype(str).str.strip()
    inst_new = inst_new.loc[inst_new["nodesdata"].isin(nodes)].copy()
    inst_new = inst_new.loc[inst_new["techs"].isin(wind_techs)].copy()

    # old file contains MW; new wind file is GW
    inst_pv = inst_old.loc[inst_old["techs"].isin(pv_techs)].copy()
    inst_pv = (
        inst_pv.groupby(["nodesdata", "techs"], as_index=True)["installable_per_region"]
        .sum()
        .to_frame("installable_MW")
    )
    inst_pv["unitsUpperLimit"] = inst_pv["installable_MW"] / 1000.0

    inst_wind = inst_new.groupby(["nodesdata", "techs"], as_index=True)[
        "installable_per_region"
    ].sum()
    inst_wind = inst_wind.to_frame("installable_GW")
    inst_wind["unitsUpperLimit"] = inst_wind["installable_GW"]

    inst = pd.concat(
        [
            inst_pv[["unitsUpperLimit"]],
            inst_wind[["unitsUpperLimit"]],
        ],
        axis=0,
    )
    inst.index = inst.index.set_names(["nodesdata", "techs"])

    totals = inst["unitsUpperLimit"].groupby("techs").sum().sort_values(ascending=False)
    print(
        "[add_renewables] installable unitsUpperLimit (GW): "
        + ", ".join(f"{t}: {totals[t]:.1f}" for t in totals.index)
    )

    # ----------------------------
    # 2) converter_techparam
    # ----------------------------
    techparam_idx = pd.MultiIndex.from_product(
        [techs, years_all_str], names=["techs", "years"]
    )
    techparam = pd.DataFrame(index=techparam_idx)
    techparam["activityUpperLimit"] = 1.0

    for y in years_all:
        ystr = str(y)
        techparam.loc[idx[pv_techs, [ystr]], "lifeTime"] = 35 if y <= 2020 else 40
        techparam.loc[idx[wind_techs, [ystr]], "lifeTime"] = 27 if y <= 2020 else 30

    m["Base"].parameter.add(techparam, "converter_techparam")

    # ----------------------------
    # 3) converter_capacityparam
    # ----------------------------
    cap_idx = pd.MultiIndex.from_product(
        [nodes, years_all_str, techs],
        names=["nodesdata", "years", "techs"],
    )
    cap = pd.DataFrame(index=cap_idx)
    cap["unitsUpperLimit"] = 0.0
    cap["unitsBuild"] = 0.0

    for n in nodes:
        for t in techs:
            ul = float(inst["unitsUpperLimit"].get((n, t), 0.0))
            cap.loc[idx[n, :, t], "unitsUpperLimit"] = ul

    # brownfield wind builds in GW (total onshore)
    brownfield = {
        ("CAN", 2003): 0.0005,
        ("CAN", 2005): 0.0001,
        ("CEN", 1999): 31.7 / 1000,
        ("CEN", 2004): 127.05 / 1000,
        ("CEN", 2007): 93 / 1000,
        ("CEN", 2011): 48.5 / 1000,
        ("CEN", 2020): 221.4 / 1000,
        ("OTG", 2007): 58 / 1000,
        ("OTG", 2009): 2.25 / 1000,
        ("OTG", 2010): 0.45 / 1000,
        ("OTG", 2011): 43.65 / 1000,
        ("OTG", 2015): 6.8 / 1000,
        ("NEL", 2010): 0.75 / 1000,
        ("NEL", 2011): 1.0 / 1000,
        ("NEL", 2014): 0.66 / 1000,
        ("TRN", 2020): 0.1333,
        ("WEL", 1993): 8.45 / 1000,
        ("WEL", 1996): 8.45 / 1000,
        ("WEL", 2009): 143 / 1000,
        ("WEL", 2014): 71.3 / 1000,
        ("WTO", 2011): 64.4 / 1000,
    }

    total_alloc = 0.0
    for (n, y), gw in brownfield.items():
        ystr = str(int(y))
        if (n not in nodes) or (ystr not in years_all_str) or (gw <= 0):
            continue

        remaining = float(gw)

        # allocate into onshore 1..4 sequentially
        for t in wind_on:
            ul = float(cap.loc[(n, ystr, t), "unitsUpperLimit"])
            already = float(cap.loc[(n, ystr, t), "unitsBuild"])

            free = ul - already
            if free <= 0:
                continue

            take = min(remaining, free)
            if take > 0:
                cap.loc[(n, ystr, t), "unitsBuild"] = already + take
                remaining -= take
                total_alloc += take

            if remaining <= 1e-12:
                break

        if remaining > 1e-9:
            raise ValueError(
                "[add_renewables] Brownfield onshore wind cannot be allocated within installables "
                f"(node={n}, year={y}, remaining={remaining:.6f} GW)."
            )

    print(f"[add_renewables] brownfield onshore wind allocated: {total_alloc:.3f} GW")

    if str(base_year) in years_all_str:
        cap.loc[idx[:, [str(base_year)], :], "noExpansion"] = 1

    viol = cap["unitsBuild"] > (cap["unitsUpperLimit"] + 1e-9)
    if viol.any():
        raise ValueError(
            "[add_renewables] unitsBuild exceeds unitsUpperLimit:\n"
            + str(cap.loc[viol].head(30))
        )

    m["Base"].parameter.add(cap.sort_index(), "converter_capacityparam")

    # ----------------------------
    # 4) converter_coefficient
    # ----------------------------
    coef_idx = pd.MultiIndex.from_product(
        [techs, years_all_str, ["Powergen"], ["Elec"]],
        names=["techs", "years", "activities", "commodities"],
    )
    coef = pd.DataFrame(index=coef_idx)
    coef["coefficient"] = 1.0
    m["Base"].parameter.add(coef, "converter_coefficient")

    # ----------------------------
    # 5) converter_activityprofile
    # ----------------------------
    feed = load_feedin_csv(weather_year=2012)

    # solar: MW -> availability by dividing by installable MW
    inst_mw = (inst["unitsUpperLimit"] * 1000.0).rename("installable_MW")
    feed_df = feed.reset_index()
    feed_df["installable_MW"] = feed_df.set_index(["nodesdata", "techs"]).index.map(
        inst_mw
    )
    feed_df["installable_MW"] = feed_df["installable_MW"].fillna(1.0)

    # detect wind vs solar by tech name
    is_wind = feed_df["techs"].str.startswith("wind_")
    is_pv = ~is_wind

    # availability column: for PV divide MW by MW; for wind use value directly
    avail = np.zeros(len(feed_df), dtype=float)
    avail[is_pv.values] = (
        feed_df.loc[is_pv, "value"] / feed_df.loc[is_pv, "installable_MW"]
    )
    avail[is_wind.values] = feed_df.loc[is_wind, "value"]

    feed_df["avail"] = np.clip(avail, 0.0, 1.0).round(3)

    avail_wide = feed_df.pivot_table(
        index=["nodesdata", "techs"],
        columns="t_model",
        values="avail",
        aggfunc="mean",
        fill_value=0.0,
    )
    avail_wide = avail_wide.reindex(columns=range(1, 8761), fill_value=0.0)
    avail_wide.columns = [f"t{str(i).zfill(4)}" for i in range(1, 8761)]

    frames_prof = []
    for ystr in years_all_str:
        tmp = avail_wide.copy()
        tmp["years"] = ystr
        tmp["type"] = "upper"
        tmp = tmp.reset_index().set_index(["nodesdata", "years", "techs", "type"])
        frames_prof.append(tmp)

    prof = pd.concat(frames_prof).sort_index()
    m["Base"].profile.add(prof, "converter_activityprofile")

    # ----------------------------
    # 6) accounting_converterunits (unchanged)
    # ----------------------------
    acc_idx = pd.MultiIndex.from_product(
        [["Invest", "OMFix"], ["global"], ["horizon"], techs, years_all_str],
        names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage"],
    )
    acc = pd.DataFrame(index=acc_idx)

    cost_years = [1950, 2030, 2040, 2050]
    pv_dec = {1950: 870, 2030: 570, 2040: 460, 2050: 410}
    pv_cen = {1950: 560, 2030: 380, 2040: 320, 2050: 290}
    w_on_c = {1950: 1330, 2030: 1040, 2040: 980, 2050: 960}
    w_off = {1950: 2120, 2030: 2287, 2040: 2168, 2050: 2130}

    def nearest_cost(y):
        y = int(y)
        return min(cost_years, key=lambda cy: abs(cy - y))

    for v in years_all_str:
        cy = nearest_cost(v)

        acc.loc[idx["Invest", "global", "horizon", "pv_decentral", v], "perUnitBuild"] = pv_dec[cy]
        acc.loc[idx["Invest", "global", "horizon", "pv_central_fixed", v], "perUnitBuild"] = pv_cen[cy]

        for t in wind_on:
            acc.loc[idx["Invest", "global", "horizon", t, v], "perUnitBuild"] = w_on_c[cy]
        for t in wind_off:
            acc.loc[idx["Invest", "global", "horizon", t, v], "perUnitBuild"] = w_off[cy]

        acc.loc[idx["Invest", "global", "horizon", :, v], "useAnnuity"] = 1
        acc.loc[idx["Invest", "global", "horizon", :, v], "interest"] = 0.06

        life_pv = 35 if int(v) <= 2020 else 40
        life_w = 27 if int(v) <= 2020 else 30
        acc.loc[idx["Invest", "global", "horizon", "pv_decentral", v], "amorTime"] = life_pv
        acc.loc[idx["Invest", "global", "horizon", "pv_central_fixed", v], "amorTime"] = life_pv
        for t in wind_techs:
            acc.loc[idx["Invest", "global", "horizon", t, v], "amorTime"] = life_w

        for t in techs:
            capex = float(acc.loc[idx["Invest", "global", "horizon", t, v], "perUnitBuild"])
            acc.loc[idx["OMFix", "global", "horizon", t, v], "perUnitTotal"] = 0.02 * capex

    m["Base"].parameter.add(acc.fillna(0.0), "accounting_converterunits")

    print("[add_renewables] finished")

def add_lithium_batteries(m):
    """
    Add lithium-ion batteries with:
      - Tech: Battery
      - Power block: converter_* + accounting_converterunits
      - Energy block: storage_* + accounting_storageunits
    """
    import pandas as pd
    import numpy as np
    global idx

    print("\n--- ADDING LITHIUM-ION BATTERIES ---")

    techs = ["Battery"]
    nodes = [n for n in m["Base"].set.nodesdata if not str(n).startswith("LNG")]

    # "years" = data years (incl. brownfield build years, etc.)
    years = sorted(int(y) for y in m["Base"].set.years)
    years_str = [str(y) for y in years]

    # accYears handled in GAMS; use horizon keyword in accounting tables
    accNodesData = ["global"]
    accYears = ["horizon"]

    # Cost trajectory anchor years
    data_years = [2020, 2030, 2040, 2050]

    def prev_year(y):
        y = int(y)
        earlier = [yy for yy in data_years if yy <= y]
        return earlier[-1] if earlier else data_years[0]

    # ------------------------------------------------------------------
    # 1) converter_techparam (Battery power block)
    # ------------------------------------------------------------------
    conv_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, years_str], names=["techs", "years"])
    )
    conv_tech["lifeTime"] = 20
    conv_tech["activityUpperLimit"] = 1.0
    m["Base"].parameter.add(conv_tech, "converter_techparam")

    # ------------------------------------------------------------------
    # 2) converter_capacityparam (Battery power capacity, GW)
    # ------------------------------------------------------------------
    conv_cap = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_str, techs], names=["nodesdata", "years", "techs"])
    )
    conv_cap["unitsUpperLimit"] = 50.0
    if 2020 in years:
        conv_cap.loc[idx[:, [str(2020)], :], "noExpansion"] = 1
    m["Base"].parameter.add(conv_cap, "converter_capacityparam")

    # ------------------------------------------------------------------
    # 3) converter_coefficient (Charge/Discharge)
    # Use Elecbattery exactly as in paste.txt
    # ------------------------------------------------------------------
    conv_coef = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [techs, years_str, ["Charge", "Discharge"], ["Elec", "Elecbattery"]],
            names=["techs", "years", "activities", "commodities"],
        )
    )
    conv_coef.loc[idx["Battery", :, "Charge", "Elec"], "coefficient"] = -1.0
    conv_coef.loc[idx["Battery", :, "Charge", "Elecbattery"], "coefficient"] = 0.975
    conv_coef.loc[idx["Battery", :, "Discharge", "Elec"], "coefficient"] = 1.0
    conv_coef.loc[idx["Battery", :, "Discharge", "Elecbattery"], "coefficient"] = -1.025
    m["Base"].parameter.add(conv_coef, "converter_coefficient")

    # ------------------------------------------------------------------
    # 4) storage_techparam (Battery energy block)
    # ------------------------------------------------------------------
    stor_tech = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, years_str], names=["techs", "years"])
    )
    stor_tech["lifeTime"] = 20
    stor_tech["levelUpperLimit"] = 1.0
    m["Base"].parameter.add(stor_tech, "storage_techparam")

    # ------------------------------------------------------------------
    # 5) storage_sizeparam (4h battery => 4 GWh per GW power unit)
    # ------------------------------------------------------------------
    stor_size = pd.DataFrame(
        index=pd.MultiIndex.from_product([techs, years_str, ["Elecbattery"]],
                                         names=["techs", "years", "commodities"])
    )
    stor_size.loc[idx["Battery", :, "Elecbattery"], "size"] = 4.0
    stor_size.loc[idx["Battery", :, "Elecbattery"], "selfdischarge"] = 0.0
    m["Base"].parameter.add(stor_size, "storage_sizeparam")

    # ------------------------------------------------------------------
    # 6) storage_reservoirparam (unit bounds)
    # ------------------------------------------------------------------
    stor_res = pd.DataFrame(
        index=pd.MultiIndex.from_product([nodes, years_str, techs], names=["nodesdata", "years", "techs"])
    )
    stor_res["unitsUpperLimit"] = 30.0
    if 2020 in years:
        stor_res.loc[idx[:, [str(2020)], :], "noExpansion"] = 1
    m["Base"].parameter.add(stor_res, "storage_reservoirparam")

    # ------------------------------------------------------------------
    # 7) accounting_converterunits (power CAPEX/OMFix, horizon)
    # vintage axis = years_str (match your other blocks)
    # ------------------------------------------------------------------
    conv_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], accNodesData, accYears, techs, years_str],
            names=["indicator", "accNodesData", "accYears", "convertertechs", "vintage"],
        )
    ).sort_index()

    capex_power = {2020: 117, 2030: 55, 2040: 37, 2050: 30}  # M€/GW
    for v in years_str:
        yy = prev_year(v)
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "perUnitBuild"] = capex_power[yy]
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "amorTime"] = 20
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "useAnnuity"] = 1
        conv_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "interest"] = 0.06

        conv_acc.loc[idx["OMFix", "global", "horizon", "Battery", v], "perUnitTotal"] = capex_power[yy] * 0.014

    m["Base"].parameter.add(conv_acc.fillna(0.0), "accounting_converterunits")

    # ------------------------------------------------------------------
    # 8) accounting_storageunits (energy CAPEX/OMFix, horizon)
    # ------------------------------------------------------------------
    stor_acc = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["Invest", "OMFix"], accNodesData, accYears, techs, years_str],
            names=["indicator", "accNodesData", "accYears", "storagetechs", "vintage"],
        )
    ).sort_index()

    capex_energy = {2020: 234 * 4, 2030: 110 * 4, 2040: 76 * 4, 2050: 61 * 4}  # M€/GW of 4h
    for v in years_str:
        yy = prev_year(v)
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "perUnitBuild"] = capex_energy[yy]
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "amorTime"] = 20
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "useAnnuity"] = 1
        stor_acc.loc[idx["Invest", "global", "horizon", "Battery", v], "interest"] = 0.06

        stor_acc.loc[idx["OMFix", "global", "horizon", "Battery", v], "perUnitTotal"] = capex_energy[yy] * 0.014

    m["Base"].parameter.add(stor_acc.fillna(0.0), "accounting_storageunits")

    print("[add_lithium_batteries] done")

# CO2 emission limits and budgets

def add_emission_limit(m, year=2050, upper_value=0.0):
    """
    Hard cap on indicator CO2_emission in a specific model year.
    """
    df = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["global"], [str(int(year))], ["CO2_emission"]],
            names=["accNodesData", "accYears", "indicator"],
        )
    )

    df["useUpper"] = 1
    df["upperValue"] = float(upper_value)

    m["Base"].parameter.add(df, "accounting_indicatorbounds")

def add_emission_budget(m, upper_value):
    """
    Cumulative emissions budget over the optimisation horizon.

    upper_value: total budget (units must match CO2_emission accounting units)
    endyear: length of horizon in years (PREVIOUSLY we used 25 for 2020->2050)
    """
    df = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [["global"], ["horizon"], ["CO2_emission"]],
            names=["accNodesData", "accYears", "indicator"],
        )
    )
    df["integral"] = 1
    # endyear meaning is model-dependent; keep your intent (25 years) but as float/int:
    df["endyear"] = 25
    df["useUpper"] = 1
    df["upperValue"] = float(upper_value)

    m["Base"].parameter.add(df, "accounting_indicatorbounds")


if __name__ == "__main__":
    start = time.time()

    # Clean data inputs TO START fresh.
    if data_dir.exists():
        shutil.rmtree(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    m = {"Base": Instance(index_names=True, column_names=True, datadir=data_dir)}

    # fundamental builds
    add_scope(m)
    add_demand(m)
    add_network(m)
    add_accounting(m)

    # optional builds

    if include_conventional_pp:
        add_thermal(m)
        add_gas_turbines(m)

    if include_renewables:
        add_hydro(m)
        add_geothermal(m)
        add_renewables(m)
        add_lithium_batteries(m)

    # if include_heat_transp:
        

    
    if include_emissions_constraints:
        add_emission_limit(m, year=2050, upper_value=0.0)
    
    validate_scope(m)
    m["Base"].write(project_path=data_dir, fileformat="csv", float_format="{:.4g}".format)

    elapsed = time.time() - start
    print(f"\n--- [REMix {md.version('remix.framework')}] Built for {case_name} completed in {int(elapsed // 60)} m {round(elapsed % 60, 1):.0f} s ---")
