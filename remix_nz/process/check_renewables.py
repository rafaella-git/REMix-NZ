from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(r"C:\Local\REMix\remix_nz\input\profiles\results_w_corr")
TS_FILE = BASE / "timeseries_norm_2012.csv"
INST_FILE = BASE / "installable_per_region.csv"
ANN_FILE = BASE / "annual_energy_per_region.csv"


def read_ts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={c: c.strip() for c in df.columns})

    need = {"technology", "t", "region", "timeseries_norm"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}. Found: {list(df.columns)}")

    df["technology"] = df["technology"].astype(str).str.strip()
    df["region"] = df["region"].astype(str).str.strip()

    # If already looks like "2012-01-13 00:30:00", keep as is
    if np.issubdtype(df["t"].dtype, np.datetime64):
        t = df["t"]
    else:
        t_raw = df["t"].astype(str).str.replace("\u00A0", " ", regex=False).str.strip()
        # first try ISO-like
        t = pd.to_datetime(t_raw, errors="coerce")
        bad = t.isna()
        if bad.any():
            # second try NZ-style dayfirst
            t2 = pd.to_datetime(t_raw[bad], dayfirst=True, errors="coerce")
            t.loc[bad] = t2

    if t.isna().any():
        # log but do NOT drop whole dates like 2012-01-13
        print("\n[read_ts] still have bad timestamps (up to 20 shown):")
        print(df.loc[t.isna(), ["technology", "region", "t"]].head(20).to_string(index=False))
        # fill them with first valid timestamp just so the analyser can continue
        first_valid = t.dropna().iloc[0]
        t = t.fillna(first_valid)

    df["t"] = t

    df["timeseries_norm"] = pd.to_numeric(df["timeseries_norm"], errors="coerce")
    nanv = df["timeseries_norm"].isna()
    if nanv.any():
        print("\n[read_ts] NaN timeseries_norm count:", int(nanv.sum()), "-> filling with 0.0")
        df.loc[nanv, "timeseries_norm"] = 0.0

    df.loc[df["timeseries_norm"] < 0, "timeseries_norm"] = 0.0
    return df


def longest_run_true(mask: np.ndarray) -> int:
    if mask.size == 0:
        return 0
    x = mask.astype(np.int8)
    d = np.diff(np.r_[0, x, 0])
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    if starts.size == 0:
        return 0
    return int((ends - starts).max())



def read_installable(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={c: c.strip() for c in df.columns})

    need = {"technology", "region", "installable_per_region"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}. Found: {list(df.columns)}")

    df["technology"] = df["technology"].astype(str).str.strip()
    df["region"] = df["region"].astype(str).str.strip()
    df["installable_per_region"] = pd.to_numeric(df["installable_per_region"], errors="coerce").fillna(0.0)

    return df[["technology", "region", "installable_per_region"]]


def read_annual(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={c: c.strip() for c in df.columns})

    need = {"technology", "region", "annual_energy_per_region"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}. Found: {list(df.columns)}")

    df["technology"] = df["technology"].astype(str).str.strip()
    df["region"] = df["region"].astype(str).str.strip()
    df["annual_energy_per_region"] = pd.to_numeric(df["annual_energy_per_region"], errors="coerce").fillna(0.0)

    return df[["technology", "region", "annual_energy_per_region"]]


def truncate_to_8760(g: pd.DataFrame) -> pd.DataFrame:
    # sort then keep first 8760 points (drop extra day from leap-year)
    g = g.sort_values("t")
    if len(g) <= 8760:
        return g
    return g.iloc[:8760].copy()


def main(topn: int = 20):
    for p in [TS_FILE, INST_FILE, ANN_FILE]:
        if not p.exists():
            raise FileNotFoundError(f"Missing: {p}")

    ts = read_ts(TS_FILE)
    inst = read_installable(INST_FILE)
    ann = read_annual(ANN_FILE)

    print("\n=== TIME COVERAGE (raw) ===")
    print("min t:", ts["t"].min())
    print("max t:", ts["t"].max())
    print("unique minutes:", sorted(ts["t"].dt.minute.unique().tolist()))

    print("\n=== COUNTS ===")
    print("rows:", len(ts))
    print("regions:", ts["region"].nunique(), sorted(ts["region"].unique().tolist()))
    print("techs:", ts["technology"].nunique(), sorted(ts["technology"].unique().tolist()))

    # duplicates check
    dup = ts.duplicated(subset=["region", "technology", "t"]).sum()
    print("\nDuplicate (region,technology,t) rows:", int(dup))

    # merge installable + annual info
    meta = inst.merge(ann, on=["technology", "region"], how="outer").fillna(0.0)

    # pairs with installable=0 but annual>0 (report; you want to treat these as unusable)
    bad_meta = meta[(meta["installable_per_region"] == 0.0) & (meta["annual_energy_per_region"] > 1e-6)]
    print("\n=== installable=0 but annual_energy>0 ===")
    print("count:", len(bad_meta))
    if len(bad_meta):
        print(bad_meta.sort_values("annual_energy_per_region", ascending=False).head(50).to_string(index=False))

    # Build per-(region,tech) diagnostics after:
    # - truncating to 8760
    # - forcing availability=0 for installable=0
    meta_key = meta.set_index(["technology", "region"])

    rows = []
    missing_pairs = []

    print("\n=== PER (REGION,TECH) SERIES CHECKS (after truncate-to-8760 and installable=0->avail=0) ===")
    for (reg, tech), g in ts.groupby(["region", "technology"], sort=False):
        g = g.copy()

        # truncate to 8760
        g = truncate_to_8760(g)

        # installable lookup
        ul = float(meta_key.get("installable_per_region").get((tech, reg), np.nan)) if (tech, reg) in meta_key.index else np.nan
        ae = float(meta_key.get("annual_energy_per_region").get((tech, reg), np.nan)) if (tech, reg) in meta_key.index else np.nan
        if np.isnan(ul):
            missing_pairs.append((reg, tech))
            ul = 0.0
            ae = 0.0

        # force availability=0 if installable is 0
        if ul == 0.0:
            g["timeseries_norm"] = 0.0

        s = g.sort_values("t")["timeseries_norm"].to_numpy(dtype=float)

        # summary stats
        mn, mx, mean = float(s.min()), float(s.max()), float(s.mean())
        share_zeros = float((s == 0.0).mean())
        lz = longest_run_true(s == 0.0)

        # simple completeness check: do we have exactly 8760 points?
        n = int(len(s))
        ok_len = (n == 8760)

        # check if time step looks consistent: compute modal delta in minutes
        dt = g.sort_values("t")["t"]
        deltas = dt.diff().dropna().dt.total_seconds().div(60.0)
        modal_dt = float(deltas.value_counts().idxmax()) if len(deltas) else np.nan

        rows.append(
            (reg, tech, n, ok_len, modal_dt, ul, ae, mn, mx, mean, share_zeros, lz)
        )

    diag = pd.DataFrame(
        rows,
        columns=[
            "region", "technology",
            "n_points", "n_is_8760", "modal_dt_minutes",
            "installable_GW", "annual_energy",
            "min", "max", "mean",
            "share_zeros", "longest_zero_run_steps",
        ],
    )

    # Print key summaries
    print("\n=== LENGTH / TIME-STEP SUMMARY ===")
    print("Pairs not equal to 8760 points:", int((~diag["n_is_8760"]).sum()))
    print("Modal timestep (minutes) value counts (top 10):")
    print(diag["modal_dt_minutes"].value_counts(dropna=False).head(10).to_string())

    print("\n=== TOP by longest zero run (steps) ===")
    print(diag.sort_values(["longest_zero_run_steps", "share_zeros"], ascending=[False, False]).head(topn).to_string(index=False))

    print("\n=== TOP by share_zeros ===")
    print(diag.sort_values(["share_zeros", "longest_zero_run_steps"], ascending=[False, False]).head(topn).to_string(index=False))

    print("\n=== TOP by mean availability ===")
    print(diag.sort_values("mean", ascending=False).head(topn).to_string(index=False))

    print("\n=== Max availability > 1 (should not happen if truly normalized) ===")
    gt1 = diag[diag["max"] > 1.0]
    print("count:", len(gt1))
    if len(gt1):
        print(gt1.sort_values("max", ascending=False).head(50).to_string(index=False))

    print("\n=== Missing (region,tech) in meta (installable/annual) ===")
    print("count:", len(missing_pairs))
    if len(missing_pairs):
        print(missing_pairs[:50])

    # Cross-check: if installable>0 but mean is 0, that suggests the whole series is zeroed/missing
    weird = diag[(diag["installable_GW"] > 0) & (diag["mean"] == 0.0)]
    print("\n=== installable>0 but mean availability==0 (likely missing series) ===")
    print("count:", len(weird))
    if len(weird):
        print(weird.head(50).to_string(index=False))

    # Totals by technology (after applying your "installable=0 => avail=0" rule)
    print("\n=== BY TECHNOLOGY summary (mean of means across regions) ===")
    bytech = diag.groupby("technology").agg(
        regions=("region", "nunique"),
        avg_mean=("mean", "mean"),
        max_max=("max", "max"),
        avg_share_zeros=("share_zeros", "mean"),
        total_installable_GW=("installable_GW", "sum"),
        total_annual_energy=("annual_energy", "sum"),
    ).sort_values("total_installable_GW", ascending=False)

    print(bytech.to_string())


if __name__ == "__main__":
    main(topn=25)
