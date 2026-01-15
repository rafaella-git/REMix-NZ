# This script shows the structure of Remix schemas and how they are represented in the Instance.
# Canessa. Jan 2026

import sys
import pandas as pd
from remix.framework.api.instance import Instance
from remix.framework.schema.templates import sm

def schema_keys(label: str):
    sc = sm.get_schemas()[label]
    keys = [fk["fields"][0] for fk in sc["foreignKeys"]]
    fields = [f["name"] for f in sc["fields"]]
    value_cols = [c for c in fields if c not in keys]
    return keys, value_cols

def print_one(label: str):
    keys, value_cols = schema_keys(label)
    print(f"\n=== {label} ===")
    print("Index (foreignKeys) order:")
    for i, k in enumerate(keys):
        print(f"  {i+1}. {k}")
    print("Allowed value columns (non-index):")
    if value_cols:
        print(" ", ", ".join(value_cols))
    else:
        print("  (none)")

def main():
    print("Python:", sys.executable)
    import remix.framework
    print("remix.framework:", remix.framework.__file__)
    try:
        import importlib.metadata as md
        print("remix.framework version:", md.version("remix.framework"))
    except Exception as e:
        print("remix.framework version: (metadata not found)", e)

    schemas = sm.get_schemas()
    print("schemas loaded:", len(schemas))

    # list all schema names
    print(sorted(schemas.keys()))

    # Put any labels you care about here.
    labels = [
        # core
        "map_aggregatenodesmodel",
        "sourcesink_profile",
        "sourcesink_config",
        "sourcesink_annualsum",

        # transfer
        "transfer_linkstartend",
        "transfer_lengthparam",
        "transfer_linksparam",
        "transfer_techparam",
        "transfer_coefficient",
        "transfer_coefperflow",
        "transfer_coefperlength",

        # accounting
        "accounting_indicatorbounds",
        "accounting_perindicator",
        "accounting_sourcesinkflow",
        "accounting_transferlinks",
        "accounting_transferperlength",
        "accounting_converterunits",
        "accounting_converteractivity",
        "accounting_storageunits",
        "accounting_sourcesinkflow",
        
        # converters
        "converter_techparam",
        "converter_capacityparam",
        "converter_coefficient",
        "converter_activityprofile",# "converter_config",
        "accounting_converterunits",
        "accounting_converteractivity",
        
        # storage
        "storage_techparam",
        "storage_sizeparam",
        "storage_resevoirparam",
        # "storage_profile",
        # "storage_config",
        "accounting_storageunits",
        "accounting_storageactivity",



    ]

    for lab in labels:
        if lab not in schemas:
            print(f"\n=== {lab} ===")
            print("Not present in sm.get_schemas() on this install.")
            continue
        print_one(lab)

    # show the shapes the Instance initializes
    m = Instance(index_names=True, column_names=True)
    for lab in labels:
        short = lab.replace("map_", "").replace("set_", "")
        container = None
        if lab.startswith("map_"):
            container = m.map
        elif lab.endswith("_profile"):
            container = m.profile
        else:
            container = m.parameter

        if hasattr(container, short):
            df = getattr(container, short)
            if isinstance(df, pd.DataFrame) and isinstance(df.index, pd.MultiIndex):
                print(f"\n[instance] {short}: nlevels={df.index.nlevels}, names={df.index.names}")
            elif isinstance(df, pd.DataFrame):
                print(f"\n[instance] {short}: index type={type(df.index)}")
        else:
            # some labels may be stored with different naming; schema print above is authoritative
            pass

if __name__ == "__main__":
    main()
