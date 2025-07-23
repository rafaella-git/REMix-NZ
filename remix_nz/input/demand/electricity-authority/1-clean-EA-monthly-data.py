# This code takes the info from previously downloaded monthly data Electricity Authority - EMI (market statistics and tools) https://www.emi.ea.govt.nz/Wholesale/Reports/R_NSPL_DR?_si=v|3
# It processes the data to put them in hourly resolution (previously in trading periods, i.e. half hours)
# It verifies the original data has the correct amount of trading periods, if it doesnt

# There is 1 main outcome:
# 1. a monthly hourly demand 

import pandas as pd
from datetime import datetime
from calendar import monthrange
import os

# Path to your CSV files
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
output_path = input_path  

month_names = ["", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

for month_num in range(1, 13):
    month_label = month_names[month_num]
    input_filename = f"{input_path}Demand_trends_2019_{month_label}.csv"

    if not os.path.exists(input_filename):
        print(f"❌ File not found: {input_filename} — skipping.")
        continue

    print(f"📥 Processing: {input_filename}")

    # Load data
    df = pd.read_csv(input_filename, skiprows=11)
    df.columns = df.columns.str.strip()

    # Drop Demand ($)
    df = df.drop(columns=["Demand ($)"], errors='ignore')

    # Convert types
    df["Demand (GWh)"] = pd.to_numeric(df["Demand (GWh)"], errors="coerce")
    df["Period start"] = pd.to_datetime(df["Period start"], dayfirst=True, errors="coerce")
    df["Period end"] = pd.to_datetime(df["Period end"], dayfirst=True, errors="coerce")

    # Drop rows with missing timestamps or demand
    df = df.dropna(subset=["Period start", "Period end", "Demand (GWh)"])
    df = df.sort_values(["Region ID", "Period start"]).reset_index(drop=True)

    # Combine Half-Hours into Hourly
    df["row_number"] = df.groupby("Region ID").cumcount()
    df["hour_group"] = df["row_number"] // 2

    hourly_df = df.groupby(["Region ID", "Region", "hour_group"]).agg({
        "Period start": "first",
        "Period end": "last",
        "Demand (GWh)": "sum"
    }).reset_index()

    # Assign Hour of Year
    year = hourly_df["Period start"].dt.year.min()
    start_of_year = datetime(year, 1, 1)
    full_hours = pd.date_range(start=start_of_year, periods=8760, freq="H")
    hour_map = {ts: f"t{i+1:04d}" for i, ts in enumerate(full_hours)}
    hourly_df["Hour of Year"] = hourly_df["Period start"].map(hour_map)

    # Final formatting
    hourly_df = hourly_df.sort_values(["Region ID", "Period start"]).reset_index(drop=True)
    hourly_df = hourly_df[["Hour of Year", "Period start", "Period end", "Region ID", "Region", "Demand (GWh)"]]

    # Get summary info
    unique_hours = hourly_df[["Period start", "Period end"]].drop_duplicates().shape[0]
    region_count = hourly_df["Region ID"].nunique()
    days_in_month = monthrange(year, month_num)[1]
    expected_hours = days_in_month * 24
    missing_hours = expected_hours - unique_hours

    # Save output
    month_str = f"{month_num:02d}"
    output_filename = f"{output_path}/hourly_demand_{month_str}_{year}.csv"
    hourly_df.to_csv(output_filename, index=False)
    print(f"📁 Output saved as: {output_filename}")

    # Check for missing hours
    month_start = datetime(year, month_num, 1)
    month_end = datetime(year, month_num, days_in_month, 23, 0)
    all_month_hours = pd.date_range(start=month_start, end=month_end, freq='h')
    present_hours = hourly_df["Period start"].drop_duplicates().sort_values()
    missing_timestamps = all_month_hours.difference(present_hours)

    missing_info = []
    for ts in missing_timestamps:
        missing_regions = set(hourly_df["Region ID"].unique()) - set(
            hourly_df[hourly_df["Period start"] == ts]["Region ID"]
        )
        missing_info.append({
            "Missing Hour": ts,
            "Missing Region Count": len(missing_regions),
            "Missing Region IDs": ", ".join(sorted(missing_regions))
        })

    if missing_info:
        missing_df = pd.DataFrame(missing_info)
        report_file = f"{output_path}/missing_hours_report_{month_str}_{year}.csv"
        missing_df.to_csv(report_file, index=False)
        print(f"⚠️ Found {len(missing_timestamps)} missing hourly timestamps.")
        print(f"📄 Missing hours report saved as: {report_file}")

    completeness_msg = (
        "✅ all hours present!"
        if missing_hours == 0 else
        f"⚠️ missing {missing_hours} hours (~{round(missing_hours / 24, 2)} days)"
    )

    print(
        f"📆 {month_label.title()} summary: {unique_hours} unique hours, "
        f"{region_count} regions — {completeness_msg}"
    )
    print("-" * 50)
