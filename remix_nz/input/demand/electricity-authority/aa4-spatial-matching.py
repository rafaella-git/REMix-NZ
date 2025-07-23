import pandas as pd

# Define useful shortcuts
shp_attrcol = "id"
idx = pd.IndexSlice


input_path = "C:/Local/REMix/remix_nz/input"
nodal_file = f"{input_path}/demand/electricity-authority/hourly_demand_nodal_2019.csv"
sites_file = f"{input_path}/demand/electricity-authority/sites-transpower.csv"
path_geo = f"{input_path}/shapefiles"

region_type = "admin" # "remix"


if region_type == "admin":
    geojson_file = "11regionsNZ.geojson"
    shp_file = f"{path_geo}/11regionsNZ"

if region_type == "remix"
    geojson_file = "16regionsNZ.geojson"
    shp_file = f"{path_geo}/16regionsNZ"

    # Define coordinates
    nodes_lat = [-35.8758611, -36.9547333, -38.4192333, -39.3342222, -37.9867861, -39.5516472, -40.2813000, -41.1502722, -41.6735722, -43.8619833, -45.481222]
    nodes_lon = [174.4669472, 174.8625250, 175.8000111, 174.3204333, 176.8294056, 176.8208167, 175.6404750, 174.9811500, 172.8737917, 171.3427694, 169.3195194]
    nodes_coords = [(lat, lon) for lat, lon in zip(nodes_lat, nodes_lon)]

    nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG"]
    nodes_lat=[-35.8758611, -36.9547333, -38.4192333, -39.3342222, -37.9867861, -39.5516472, -40.2813000, -41.1502722, -41.6735722, -43.8619833, -45.481222]
    nodes_lon=[174.4669472, 174.8625250, 175.8000111, 174.3204333, 176.8294056, 176.8208167, 175.6404750, 174.9811500, 172.8737917, 171.3427694, 169.3195194]    
    nodes_coords=[(nodes_lat[i],nodes_lon[i]) for i in range(0,len(nodes_lon))]
    # dictionary with 'Node': [lat,long]
    coord_dict = dict.fromkeys(nodes_lst)
    for key, value in zip(coord_dict.keys(), nodes_coords):
        coord_dict[key] = value
    print(coord_dict)


# Read the demand data
df = pd.read_csv(nodal_file, parse_dates=["Period start", "Period end"])

# Extract 3-letter code (Codigo) and node name (Nombre) from 'Region'
df[['Codigo', 'Nombre']] = df['Region'].str.extract(r'^([A-Z]{3})\d+\s*-\s*(.+)$')

# === Step 2: Load site location data ===
sites_df = pd.read_csv(sites_file)

# Ensure MXLOCATION is uppercase (should already be, but just in case)
sites_df['MXLOCATION'] = sites_df['MXLOCATION'].str.upper()

# === Step 3: Merge on Codigo <-> MXLOCATION ===
merged_df = df.merge(
    sites_df[['MXLOCATION', 'X', 'Y', 'description']],
    left_on='Codigo',
    right_on='MXLOCATION',
    how='left'
)

# === Step 4: Print useful stats ===

# Total unique codes in demand data
total_codes = df['Codigo'].nunique()

# Codes successfully matched
matched_codes = merged_df['MXLOCATION'].notna().sum()

print(f"Total unique region codes in demand data: {total_codes}")
print(f"Rows matched with site coordinates: {matched_codes} / {len(df)}")
print("\nSample of matched rows:")
print(merged_df[['Region', 'Codigo', 'Nombre', 'description', 'X', 'Y']].dropna().head())

# Show unmatched codes
unmatched = merged_df[merged_df['MXLOCATION'].isna()]['Codigo'].unique()
print("\nUnmatched codes:")
print(unmatched)

# Optional: check if coordinates are in New Zealand bounds (NZTM approx)
valid_coords = merged_df['X'].between(1e6, 2.5e6) & merged_df['Y'].between(5e6, 6.5e6)
print(f"\nRows with valid NZTM coordinates: {valid_coords.sum()} / {len(merged_df)}")

