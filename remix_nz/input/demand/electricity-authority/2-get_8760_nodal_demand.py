# This code takes the info from previously created hourly demand for each month and each node (originally processed from Electricity Authority - EMI (market statistics and tools) https://www.emi.ea.govt.nz/Wholesale/Reports/R_NSPL_DR?_si=v|3)
# It verifies that the original data has the correct amount of trading periods, if it doesnt, it fills them to then unify all months in a yearly file

# There is 1 main outcome:
# 1. an annual hourly demand (8760 hours for each of the nodes) for 2019
 
import pandas as pd
from datetime import datetime
from calendar import monthrange
import os

# Define your input folder
input_path = "C:/Local/REMix/remix_nz/input/demand/electricity-authority/"
year = 2019
month_list = [f"{i:02d}" for i in range(1, 13)]
monthly_data = []
missing_report = []

for month_str in month_list:
    file_path = f"{input_path}hourly_demand_{month_str}_{year}.csv"
    
    if not os.path.exists(file_path):
        print(f"❌ Missing file: {file_path}. Skipping.")
        continue

    print(f"📥 Processing file: {file_path}")
    df = pd.read_csv(file_path, parse_dates=["Period start", "Period end"])
    month_num = int(month_str)
    days = monthrange(year, month_num)[1]
    expected_hours = pd.date_range(
        start=datetime(year, month_num, 1, 0, 0),
        end=datetime(year, month_num, days, 23, 0),
        freq="H"
    )

    region_ids = df["Region ID"].unique()
    monthly_cleaned = []

    for region in region_ids:
        region_df = df[df["Region ID"] == region].copy()
        region_df = region_df.set_index("Period start").sort_index()

        full_df = pd.DataFrame(index=expected_hours)
        full_df.index.name = "Period start"
        full_df["Region ID"] = region
        full_df["Region"] = region_df["Region"].iloc[0] if not region_df.empty else "Unknown"

        full_df = full_df.merge(
            region_df[["Period end", "Demand (GWh)"]],
            left_index=True,
            right_index=True,
            how="left"
        )

        full_df["Period end"] = full_df.index + pd.Timedelta(hours=1)

        total_missing = full_df["Demand (GWh)"].isna().sum()

        ffilled = full_df["Demand (GWh)"].ffill()
        ffill_count = ffilled.notna().sum() - full_df["Demand (GWh)"].notna().sum()

        bfilled = ffilled.bfill()
        bfill_count = bfilled.notna().sum() - ffilled.notna().sum()

        zero_fill_count = total_missing - ffill_count - bfill_count

        full_df["Demand (GWh)"] = bfilled.fillna(0)

        if total_missing > 0:
            missing_report.append({
                "Region ID": region,
                "Region": full_df["Region"].iloc[0],
                "Month": month_str,
                "Missing Hours": total_missing,
                "Filled with ffill": ffill_count,
                "Filled with bfill": bfill_count,
                "Filled with 0": zero_fill_count
            })

        monthly_cleaned.append(full_df.reset_index())

    monthly_data.append(pd.concat(monthly_cleaned, ignore_index=True))

# Combine all months
full_year_df = pd.concat(monthly_data, ignore_index=True)

# Sort and assign hour of year
full_year_df = full_year_df.sort_values(["Region ID", "Period start"]).reset_index(drop=True)
full_year_df["Hour of Year"] = full_year_df.groupby("Region ID").cumcount() + 1
full_year_df["Hour of Year"] = full_year_df["Hour of Year"].apply(lambda x: f"t{x:04d}")

# Add Region Code and Region Name columns
full_year_df[["Region Code", "Region Name"]] = full_year_df["Region"].str.extract(r'^([A-Z]{3})\d+\s*-\s*(.+)$')

# Reorder columns
full_year_df = full_year_df[[
    "Hour of Year", "Period start", "Period end",
    "Region ID", "Region", "Region Code", "Region Name", "Demand (GWh)"
]]

# ✅ Validation
region_counts = full_year_df.groupby("Region ID").size()
missing_nodes = region_counts[region_counts != 8760]
if not missing_nodes.empty:
    print("⚠️ The following regions do not have 8760 hours:")
    print(missing_nodes)
else:
    print("✅ All regions have 8760 hours.")

# Missing data report
if missing_report:
    missing_report_df = pd.DataFrame(missing_report)
    print("\n📊 Missing data handling summary:")
    print(missing_report_df.sort_values(["Region ID", "Month"]).to_string(index=False))

#  Preview final dataframe
print("\n📋 Final cleaned dataframe preview:")
print(full_year_df.head(10))

# Save
full_year_df.to_csv(f"{input_path}/hourly_demand_all_{year}_cleaned.csv", index=False)
missing_report_df.to_csv(f"{input_path}/missing_fill_report_{year}.csv", index=False)
