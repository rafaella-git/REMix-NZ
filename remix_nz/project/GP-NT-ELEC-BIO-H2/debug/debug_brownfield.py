import pandas as pd
import numpy as np
from pathlib import Path

idx = pd.IndexSlice

path_input = Path("C:/Local/REMix/remix_nz/input")
group_name = "GP-NT-ELEC-BIO-H2"
path_profiles = path_input / "profiles"
path_brownfield = path_input / "brownfield"

NODES = ["AKL","BOP","CAN","CEN","HBY","NEL","NIS","OTG","TRN","WEL","WTO"]
BASEYEAR = 2020

# Reference "original database" node totals in GW (brownfield criterion Year_built < 2020)
REF_GW_BY_NODE = pd.DataFrame.from_dict(
    {
        "AKL": {"Thermal_Bio": 0.0155, "CCGT": 0.0,   "Thermal_Coal": 0.112, "Thermal_Diesel": 0.0,    "GT": 0.0036, "OCGT": 0.0,   "wind": 0.0},
        "BOP": {"Thermal_Bio": 0.037,  "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.0373, "OCGT": 0.01,  "wind": 0.0},
        "CAN": {"Thermal_Bio": 0.0034, "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0047, "GT": 0.0,    "OCGT": 0.0,   "wind": 0.0006},
        "CEN": {"Thermal_Bio": 0.0,    "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.0,    "OCGT": 0.0,   "wind": 0.52165},
        "HBY": {"Thermal_Bio": 0.0128, "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.155,  "GT": 0.0,    "OCGT": 0.0,   "wind": 0.0},
        "NEL": {"Thermal_Bio": 0.0,    "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0006, "GT": 0.0,    "OCGT": 0.0,   "wind": 0.00241},
        "NIS": {"Thermal_Bio": 0.0098, "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.018,  "GT": 0.0,    "OCGT": 0.0,   "wind": 0.0},
        "OTG": {"Thermal_Bio": 0.0014, "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.0,    "OCGT": 0.0,   "wind": 0.11115},
        "TRN": {"Thermal_Bio": 0.0,    "CCGT": 0.377, "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.0786, "OCGT": 0.435, "wind": 0.1333},
        "WEL": {"Thermal_Bio": 0.0038, "CCGT": 0.0,   "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.01,   "OCGT": 0.0,   "wind": 0.22295},
        "WTO": {"Thermal_Bio": 0.0459, "CCGT": 0.385, "Thermal_Coal": 0.0,   "Thermal_Diesel": 0.0,    "GT": 0.75,   "OCGT": 0.092, "wind": 0.0644},
    },
    orient="index",
).fillna(0.0)

def _print_node_totals(df_gw_by_node, label, nodes):
    df_gw_by_node = df_gw_by_node.reindex(nodes).fillna(0.0)
    print(f"=== {label} totals by node (GW) ===")
    print(df_gw_by_node.sort_values(ascending=False).to_string())
    print("TOTAL:", float(df_gw_by_node.sum()))
    print()

def _compare_to_reference(calc_series, ref_series, label, nodes, tol=1e-6):
    calc_series = calc_series.reindex(nodes).fillna(0.0)
    ref_series = ref_series.reindex(nodes).fillna(0.0)

    diff = (calc_series - ref_series).abs()
    bad = diff > tol

    print(f"=== CHECK {label}: calc vs reference (GW) ===")
    if not bad.any():
        print("OK: all nodes match within tolerance", tol)
        print()
        return

    out = pd.DataFrame({"calc": calc_series, "ref": ref_series, "abs_diff": diff})
    out = out.loc[bad].sort_values("abs_diff", ascending=False)
    print(out.to_string())
    print()

def debug_brownfield_wind_allocation(pathprofiles, nodes, baseyear=2020):
    print("=== DEBUG brownfield wind allocation ===")
    print("pathprofiles:", Path(pathprofiles).resolve())
    print("baseyear:", baseyear)
    print("nodes:", nodes)
    print("n_nodes:", len(nodes))
    print()

    windon = [f"wind_onshore_{i}" for i in [1, 2, 3, 4]]
    windoff = [f"wind_offshore_{i}" for i in [1, 2, 3, 4]]
    windtechs = windon + windoff

    instfile = Path(pathprofiles) / "results_w_corr" / "installable_per_region.csv"
    if not instfile.exists():
        raise FileNotFoundError(f"installables file not found: {instfile}")

    instnew = pd.read_csv(instfile).rename(columns={c: c.strip() for c in pd.read_csv(instfile, nrows=0).columns})
    instnew = pd.read_csv(instfile)
    instnew = instnew.rename(columns={c: c.strip() for c in instnew.columns})

    need = {"technology", "region", "installable_per_region"}
    missing = need - set(instnew.columns)
    if missing:
        raise ValueError(f"installables missing columns {missing}, found {list(instnew.columns)}")

    instnew["technology"] = instnew["technology"].astype(str).str.strip()
    instnew["region"] = instnew["region"].astype(str).str.strip()
    instnew["installable_per_region"] = pd.to_numeric(instnew["installable_per_region"], errors="coerce").fillna(0.0)

    instnew = instnew.loc[instnew["region"].isin(nodes)].copy()
    instnew = instnew.loc[instnew["technology"].isin(windtechs)].copy()

    instwind = (
        instnew.groupby(["region", "technology"], as_index=True)["installable_per_region"]
        .sum()
        .to_frame("installableGW")
        .sort_index()
    )

    print("installables rows:", len(instwind))
    print("installables tech totals (GW):")
    print(instwind.groupby("technology")["installableGW"].sum().sort_values(ascending=False).to_string())
    print()

    yearsallstr = [str(baseyear)]
    capidx = pd.MultiIndex.from_product([nodes, yearsallstr, windtechs], names=["nodesdata", "years", "techs"])
    cap = pd.DataFrame(index=capidx, data={"unitsUpperLimit": 0.0, "unitsBuild": 0.0})

    for (n, t), row in instwind.iterrows():
        cap.loc[(n, str(baseyear), t), "unitsUpperLimit"] = float(row["installableGW"])

    brownfield_wind_gw_by_node = REF_GW_BY_NODE["wind"].to_dict()

    print("brownfield wind targets (GW) by node:")
    for n in nodes:
        v = float(brownfield_wind_gw_by_node.get(n, 0.0))
        if v != 0.0:
            print(f"  {n}: {v:.5f}")
    print("brownfield wind total (GW):", sum(float(brownfield_wind_gw_by_node.get(n, 0.0)) for n in nodes))
    print()

    for n in nodes:
        remaining = float(brownfield_wind_gw_by_node.get(n, 0.0))
        if remaining <= 0.0:
            continue

        print(f"--- node {n} target {remaining:.6f} GW ---")
        for t in windon:
            ul = float(cap.loc[(n, str(baseyear), t), "unitsUpperLimit"])
            already = float(cap.loc[(n, str(baseyear), t), "unitsBuild"])
            free = ul - already

            print(f"  tech {t}: unitsUpperLimit={ul:.6f}, unitsBuild={already:.6f}, free={free:.6f}")

            if free <= 0.0:
                continue

            take = min(remaining, free)
            if take > 0.0:
                cap.loc[(n, str(baseyear), t), "unitsBuild"] = already + take
                remaining -= take
                print(f"    take={take:.6f} -> remaining={remaining:.6f}")

            if remaining <= 1e-12:
                break

        if remaining > 1e-9:
            print(f"  ERROR: remaining {remaining:.6f} GW could not be allocated for node {n}")
            tmp = cap.loc[idx[n, str(baseyear), :], ["unitsUpperLimit", "unitsBuild"]].copy()
            tmp = tmp.reset_index().sort_values("unitsUpperLimit", ascending=False)
            print(tmp.to_string(index=False))
            raise ValueError(f"Allocation failed for node={n}, remaining={remaining:.6f} GW")

        print()

    cap2020 = cap.loc[idx[:, str(baseyear), windon], :].copy()
    by_node = cap2020.groupby("nodesdata")["unitsBuild"].sum().sort_values(ascending=False)

    _print_node_totals(by_node, "RESULT: allocated brownfield wind in 2020", nodes)
    _compare_to_reference(by_node, REF_GW_BY_NODE["wind"], "WIND", nodes)

    return cap

def _load_powerplant_db(path_brownfield):
    ppfile = Path(path_brownfield) / "power-plant-nz-database.csv"
    if not ppfile.exists():
        raise FileNotFoundError(f"power plant DB not found: {ppfile}")

    df = pd.read_csv(ppfile)
    df = df.rename(columns={c: str(c).strip() for c in df.columns})

    # Canonicalize column names from your CSV to what the debug functions expect
    colmap = {
        "Primary_fuel": "Primaryfuel",
        "Capacity_MW": "CapacityMW",
        "Year_built": "Yearbuilt",
        "Node_name": "Node",  # only if you want Node_name to override Node
    }
    for src, dst in colmap.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})

    required = {"Type", "Primaryfuel", "Yearbuilt", "CapacityMW", "Node"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"power plant DB missing columns {missing}. Found: {list(df.columns)}")

    df["Type"] = df["Type"].astype(str).str.strip()
    df["Primaryfuel"] = df["Primaryfuel"].astype(str).str.strip()
    df["Node"] = df["Node"].astype(str).str.strip()

    df["Yearbuilt"] = pd.to_numeric(df["Yearbuilt"], errors="coerce")
    df["CapacityMW"] = pd.to_numeric(df["CapacityMW"], errors="coerce")
    df = df.dropna(subset=["Yearbuilt", "CapacityMW"]).copy()

    df["Yearbuilt"] = df["Yearbuilt"].astype(int)
    df = df.loc[df["CapacityMW"] > 0].copy()

    return df

def debug_brownfield_thermal_from_ppdb(path_brownfield, nodes, baseyear=2020):
    print("=== DEBUG brownfield thermal from power-plant-nz-database.csv ===")
    print("path_brownfield:", Path(path_brownfield).resolve())
    print("brownfield criterion: Yearbuilt < baseyear")
    print("baseyear:", baseyear)
    print()

    df = _load_powerplant_db(path_brownfield)

    df = df.loc[df["Node"].isin(nodes)].copy()

    # mirror addthermalm filters/mapping [file:1]
    df = df.loc[df["Type"].eq("Thermal")].copy()

    techmap = {
        "Coal": ("ThermalCoal", "Coal"),
        "Diesel": ("ThermalDiesel", "Diesel"),
        "Biogas": ("ThermalBio", "Biofuel"),
        "Biomass": ("ThermalBio", "Biofuel"),
        "Wood": ("ThermalBio", "Biofuel"),
        "Wood waste": ("ThermalBio", "Biofuel"),
    }

    df = df.loc[df["Primaryfuel"].isin(techmap.keys())].copy()

    df["convertertechs"] = df["Primaryfuel"].map(lambda f: techmap[f][0])
    df["fuel"] = df["Primaryfuel"].map(lambda f: techmap[f][1])

    # brownfield only
    df_bf = df.loc[df["Yearbuilt"] < int(baseyear)].copy()

    # GW per node per tech
    grp = (
        df_bf.groupby(["Node", "convertertechs"], as_index=True)["CapacityMW"]
        .sum()
        .div(1000.0)
        .rename("GW")
        .sort_index()
    )

    print("rows after filters (thermal, mapped fuels, brownfield):", len(df_bf))
    print("unique Primaryfuel after mapping:", sorted(df_bf["Primaryfuel"].unique().tolist()))
    print()

    # show tech totals
    print("thermal brownfield totals by tech (GW):")
    print(grp.groupby("convertertechs").sum().sort_values(ascending=False).to_string())
    print()

    # compare vs reference (node totals)
    for ref_col, tech in [("Thermal_Bio", "ThermalBio"), ("Thermal_Coal", "ThermalCoal"), ("Thermal_Diesel", "ThermalDiesel")]:
        calc = grp.xs(tech, level="convertertechs", drop_level=False).reset_index().set_index("Node")["GW"]
        _print_node_totals(calc, f"CALC thermal {tech}", nodes)
        _compare_to_reference(calc, REF_GW_BY_NODE[ref_col], f"THERMAL {tech}", nodes)

def debug_brownfield_gas_turbines_from_ppdb(path_brownfield, nodes, baseyear=2020):
    print("=== DEBUG brownfield gas turbines from power-plant-nz-database.csv ===")
    print("path_brownfield:", Path(path_brownfield).resolve())
    print("brownfield criterion: Yearbuilt < baseyear")
    print("baseyear:", baseyear)
    print()

    df = _load_powerplant_db(path_brownfield)
    df = df.loc[df["Node"].isin(nodes)].copy()

    # mirror addgasturbinesm fuel filter [file:1]
    df = df.loc[df["Primaryfuel"].isin(["Natural gas", "Gas", "CH4"])].copy()

    # map GT type similar to addgasturbinesm [file:1]
    def mapgttype(x):
        s = str(x).upper()
        if "CCGT" in s:
            return "CCGT"
        if "OCGT" in s:
            return "OCGT"
        if "GT" in s:
            return "GT"
        return "GT"

    if "Techs" in df.columns:
        df["convertertechs"] = df["Techs"].apply(mapgttype)
    else:
        df["convertertechs"] = "GT"

    df = df.loc[df["convertertechs"].isin(["GT", "CCGT", "OCGT"])].copy()

    df_bf = df.loc[df["Yearbuilt"] < int(baseyear)].copy()

    grp = (
        df_bf.groupby(["Node", "convertertechs"], as_index=True)["CapacityMW"]
        .sum()
        .div(1000.0)
        .rename("GW")
        .sort_index()
    )

    print("rows after filters (gas fuels, mapped GT types, brownfield):", len(df_bf))
    if "Techs" in df_bf.columns:
        print("unique Techs values sample (up to 20):", df_bf["Techs"].astype(str).unique().tolist()[:20])
    print()

    print("gas turbine brownfield totals by tech (GW):")
    print(grp.groupby("convertertechs").sum().sort_values(ascending=False).to_string())
    print()

    # compare to reference
    for ref_col, tech in [("GT", "GT"), ("CCGT", "CCGT"), ("OCGT", "OCGT")]:
        calc = grp.xs(tech, level="convertertechs", drop_level=False).reset_index().set_index("Node")["GW"]
        _print_node_totals(calc, f"CALC gas {tech}", nodes)
        _compare_to_reference(calc, REF_GW_BY_NODE[ref_col], f"GAS {tech}", nodes)

if __name__ == "__main__":
    debug_brownfield_wind_allocation(pathprofiles=path_profiles, nodes=NODES, baseyear=BASEYEAR)
    debug_brownfield_thermal_from_ppdb(path_brownfield=path_brownfield, nodes=NODES, baseyear=BASEYEAR)
    debug_brownfield_gas_turbines_from_ppdb(path_brownfield=path_brownfield, nodes=NODES, baseyear=BASEYEAR)
