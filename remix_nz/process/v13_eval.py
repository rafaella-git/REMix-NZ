import numpy as np
import pandas as pd
from pathlib import Path
from remix.framework.tools.gdx import GDXEval

pd.set_option("display.width", 160)
pd.set_option("display.max_rows", 200)
pd.set_option("display.max_columns", 20)

idx = pd.IndexSlice

# ------------------------------------------------------------
# 1. Paths and GDX
# ------------------------------------------------------------
path_base = r"C:/Local/REMix"
project_name = "GP-NT-ELEC-BIO-H2"
case_name = "nz_case_GP_2020-2050"

project_dir = Path(path_base) / "remix_nz" / "project" / project_name / case_name
results_dir = project_dir / "result"
gdx_file = results_dir / f"{case_name}.gdx"

print("Project dir:", project_dir)
print("Results dir:", results_dir)
print("GDX file:", gdx_file, "exists:", gdx_file.exists())

results = GDXEval(str(gdx_file))

# ------------------------------------------------------------
# 2. Rename map
# ------------------------------------------------------------
rename_techs = {
    "H2_CCGT": "H2_CCGT",
    "pv_central_fixed": "PV",
    "pv_central_track_azimuth": "PV",
    "pv_decentral": "PV",
    "wind_onshore": "Onshore wind",
    "wind_offshore_floating": "Offshore wind",
    "wind_offshore_foundation": "Offshore wind",
    "geoth": "Geothermal",
    "COAL": "Coal",
    "DIE": "Diesel",
    "H2_storage": "H2 Storage",
    "BIO": "Biomass",
    "H2_FC": "Fuel Cell (H2)",
    "CCGT": "CCGT",
    "OCGT": "OCGT",
    "GT": "GT",
    "Thermal_Coal": "Thermal_Coal",
    "Thermal_Diesel": "Thermal_Diesel",
    "Thermal_Bio": "Thermal_Bio",
}

def safe_get(name):
    """Safely read a symbol as DataFrame with column 'value'."""
    try:
        df = results[name]
        if not isinstance(df, pd.DataFrame):
            df = df.to_frame("value")
        if "value" not in df.columns:
            df.columns = ["value"]
        print(f"\n[OK] Loaded symbol '{name}' with shape {df.shape}")
        return df
    except Exception as e:
        print(f"[WARN] Could not read symbol '{name}': {e}")
        return None

# ------------------------------------------------------------
# 3. converter_caps diagnostics
# ------------------------------------------------------------
converter_caps = safe_get("converter_caps")
if converter_caps is not None:
    converter_caps = converter_caps.rename(index=rename_techs)
    converter_caps = converter_caps[converter_caps["value"] > 0.001].dropna()

    print("\n=== converter_caps: index names ===")
    print(converter_caps.index.names)

    print("\n=== converter_caps: unique techs (after rename) ===")
    try:
        print(sorted(set(converter_caps.index.get_level_values("techs"))))
    except Exception as e:
        print("[INFO] Could not get tech level from converter_caps:", e)

    # GT / Thermal only
    try:
        tech_mask = converter_caps.index.get_level_values("techs").isin(
            ["GT", "CCGT", "OCGT", "Thermal_Coal", "Thermal_Diesel", "Thermal_Bio"]
        )
        therm_gt_caps = converter_caps[tech_mask]
        caps_elec_total = therm_gt_caps.xs(
            ("Elec", "total"), level=["commodity", "capType"], drop_level=False
        )
        print("\n=== converter_caps: GT + Thermal techs (all years, Elec/total, top 50) ===")
        print(caps_elec_total.sort_values("value", ascending=False).head(50))
    except Exception as e:
        print("[INFO] Could not slice GT/Thermal Elec/total from converter_caps:", e)

# ------------------------------------------------------------
# 4. commodity_balance_annual diagnostics
# ------------------------------------------------------------
cba = safe_get("commodity_balance_annual")
if cba is not None:
    cba = cba.rename(index=rename_techs)
    cba = cba[cba["value"].abs() > 0.01].dropna()

    print("\n=== commodity_balance_annual: index names ===")
    print(cba.index.names)

    # which Elec techs appear at all?
    try:
        elec_mask = cba.index.get_level_values("commodity") == "Elec"
        elec_techs = sorted(set(cba.index.get_level_values("techs")[elec_mask]))
        print("\n=== Elec techs present in commodity_balance_annual ===")
        print(elec_techs)
    except Exception as e:
        print("[INFO] Could not get Elec tech list:", e)

    # Global Elec net 2020
    try:
        cb_global_2020 = cba.xs(
            ("global", "2020", "Elec", "net"),
            level=["accNodesModel", "accYears", "commodity", "balanceType"],
            drop_level=False,
        )
        cb_global_2020 = (
            cb_global_2020.groupby("techs")["value"].sum().sort_values(ascending=False)
        )
        print("\n=== Global Elec net 2020 by tech (GWh) ===")
        print(cb_global_2020.round(2))
    except Exception as e:
        print("[INFO] Could not slice global/2020/Elec/net:", e)

    # Global Elec net 2050
    try:
        cb_global_2050 = cba.xs(
            ("global", "2050", "Elec", "net"),
            level=["accNodesModel", "accYears", "commodity", "balanceType"],
            drop_level=False,
        )
        cb_global_2050 = (
            cb_global_2050.groupby("techs")["value"].sum().sort_values(ascending=False)
        )
        print("\n=== Global Elec net 2050 by tech (GWh) ===")
        print(cb_global_2050.round(2))
    except Exception as e:
        print("[INFO] Could not slice global/2050/Elec/net:", e)

    # Nodal Elec net 2020 for key techs (Slack, GT, Thermal, Hydro, PV, etc.)
    try:
        cb_nodal_2020 = cba.xs(
            ("2020", "Elec", "net"),
            level=["accYears", "commodity", "balanceType"],
            drop_level=False,
        )
        focus_techs = [
            "Slack", "GT", "CCGT", "OCGT",
            "Thermal_Coal", "Thermal_Diesel", "Thermal_Bio",
            "Geothermal", "Hydro", "Battery", "PV", "Onshore wind", "Offshore wind",
        ]
        mask_focus = cb_nodal_2020.index.get_level_values("techs").isin(focus_techs)
        cb_focus = cb_nodal_2020[mask_focus]
        print("\n=== Nodal Elec net 2020 for key techs (GWh, +gen / -load) ===")
        print(
            cb_focus.groupby(["accNodesModel", "techs"])["value"]
            .sum()
            .unstack("techs")
            .round(2)
        )
    except Exception as e:
        print("[INFO] Could not build nodal 2020 Elec balance:", e)

# ------------------------------------------------------------
# 5. Hourly commodity_balance diagnostics for Elec / Slack
# ------------------------------------------------------------
cb = safe_get("commodity_balance")
if cb is not None:
    cb = cb.rename(index=rename_techs)
    print("\n=== commodity_balance: index names ===")
    print(cb.index.names)

    # Global / 2020 / Elec
    try:
        cb_elec_2020 = cb.xs(
            ("global", "2020", "Elec"),
            level=["accNodesModel", "accYears", "commodity"],
            drop_level=False,
        )
        print("\n[DEBUG] Rows for global/2020/Elec:", len(cb_elec_2020))

        gen_2020 = (
            cb_elec_2020[cb_elec_2020["value"] > 0]
            .groupby("techs")["value"]
            .sum()
            .sort_values(ascending=False)
        )
        print("\n=== Global hourly 2020 Elec: total generation by tech (GWh over year) ===")
        print(gen_2020.round(2))

        sinks_2020 = (
            cb_elec_2020[cb_elec_2020["value"] < 0]
            .groupby("techs")["value"]
            .sum()
            .sort_values()
        )
        print("\n=== Global hourly 2020 Elec: total sinks (demand+Slack) by tech (GWh over year) ===")
        print(sinks_2020.round(2))

        # Slack time series sample
        try:
            slack_2020 = cb_elec_2020.xs("Slack", level="techs", drop_level=False)
            print("\n=== First 50 rows of global/2020/Elec/Slack hourly series ===")
            print(slack_2020.head(50))
        except Exception as e2:
            print("[INFO] Could not slice Slack tech in global/2020/Elec:", e2)
    except Exception as e:
        print("[INFO] Could not slice global/2020/Elec in commodity_balance:", e)
