# This code takes the info from transpower sites NSP https://www.emi.ea.govt.nz/Wholesale/Datasets/MappingsAndGeospatial/NetworkSupplyPointsTable   
# And shapefiles (in this case remix and pypsa regions)
# So we can asign these transpower location to the regions we use

# There are 3 main outcomes:
# 1. Plotting the region boundaries and where the sites are located
# 2. Sorts the sites and assigns them a region (or provides a list of all the sites in each region, depending on what you want)
# 3. Rename the shape file from pypsa regions into words and export them

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  
import pandas as pd

# --- 1. File paths ---
shape_path = "C:/Local/REMix/remix_nz/input/shapefiles"
sites_fp = f"{shape_path}/Sites.geojson"
regions_fp_11 = f"{shape_path}/11regionsNZ.geojson"
regions_fp_19 = f"{shape_path}/19districtsnz.geojson"

# --- 2. Load data ---
sites = gpd.read_file(sites_fp)
regions_11 = gpd.read_file(regions_fp_11)
regions_19 = gpd.read_file(regions_fp_19)

# --- 3. Ensure CRS match ---
regions_11 = regions_11.to_crs(sites.crs)
regions_19 = regions_19.to_crs(sites.crs)

# --- 4. Spatial join (within) ---
joined_11 = gpd.sjoin(sites, regions_11[['geometry', 'id']], how="left", predicate="within")
joined_11 = joined_11.drop(columns=["index_right"])

joined_both = gpd.sjoin(joined_11, regions_19[['geometry', 'GADM_ID']], how="left", predicate="within")

# --- 5. Assign sites not within to closest region (by boundary distance) ---
def assign_nearest_region(sites_df, region_df, region_col):
    unmatched = sites_df[sites_df[region_col].isna()].copy()
    matched = sites_df[sites_df[region_col].notna()].copy()

    if unmatched.empty:
        return sites_df

    region_df = region_df[[region_col, 'geometry']].copy()

    # Assign based on minimum distance to region boundaries
    def find_closest_region(site_geom):
        distances = region_df.geometry.distance(site_geom)
        nearest_idx = distances.idxmin()
        return region_df.loc[nearest_idx, region_col]

    unmatched[region_col] = unmatched.geometry.apply(find_closest_region)
    return pd.concat([matched, unmatched], ignore_index=True)

# Apply to both region types
joined_both = assign_nearest_region(joined_both, regions_11, 'id')
joined_both = assign_nearest_region(joined_both, regions_19, 'GADM_ID')

# --- 6. Group into dictionaries and counts ---
def group_sites_by_region(df, region_col):
    df_valid = df.dropna(subset=[region_col, 'MXLOCATION'])
    grouped = df_valid.groupby(region_col)['MXLOCATION'].apply(list).to_dict()
    counts = df_valid.groupby(region_col)['MXLOCATION'].count().to_dict()
    return grouped, counts

region_11_dict, region_11_counts = group_sites_by_region(joined_both, 'id')
region_19_dict, region_19_counts = group_sites_by_region(joined_both, 'GADM_ID')

# --- 7. Export to CSV ---
site_region_df = joined_both[['OBJECTID', 'MXLOCATION', 'id', 'GADM_ID']].drop_duplicates()
site_region_df.to_csv("site_region_assignment.csv", index=False)
print(f"✅ Exported {len(site_region_df)} site-region assignments to 'site_region_assignment.csv'")

new_site_region_df = site_region_df.rename(columns={"GADM_ID": "Administrative Region"})

# Define mapping from GADM_ID to region names ---
region_name_map = {
    "NZ.1_1": "Auckland",
    "NZ.2_1": "Bay of Plenty",
    "NZ.3_1": "Canterbury",
    "NZ.5_1": "Gisborne",
    "NZ.6_1": "Hawkes Bay",
    "NZ.7_1": "Manawatu-Whanganui",
    "NZ.8_1": "Marlborough",
    "NZ.9_1": "Nelson",
    "NZ.11_1": "Northland",
    "NZ.12_1": "Otago",
    "NZ.14_1": "Southland",
    "NZ.15_1": "Taranaki",
    "NZ.16_1": "Tasman",
    "NZ.17_1": "Waikato",
    "NZ.18_1": "Wellington",
    "NZ.19_1": "West Coast"
}

#  Replace region codes with readable names ---
new_site_region_df["Administrative Region"] = new_site_region_df["Administrative Region"].map(region_name_map).fillna(new_site_region_df["Administrative Region"])

#  Save to new CSV ---
new_site_region_df.to_csv("site_region_assignment_named.csv", index=False)



# --- 8. Plotting ---
def get_site_bounds_with_padding(pad=0.06):
    bounds = sites.total_bounds  # [minx, miny, maxx, maxy]
    dx = (bounds[2] - bounds[0]) * pad
    dy = (bounds[3] - bounds[1]) * pad
    return (bounds[0] - dx, bounds[2] + dx), (bounds[1] - dy, bounds[3] + dy)


def save_figure(fig, name_base):
    fig.savefig(f"{name_base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{name_base}.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def plot_side_by_side():
    xlim, ylim = get_site_bounds_with_padding(pad=0.12)
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # --- 11 Regions Plot ---
    regions_11.plot(ax=axes[0], edgecolor="blue", facecolor="#add8e6", alpha=0.5)
    sites.plot(ax=axes[0], color="black", markersize=10)
    axes[0].set_xlim(xlim)
    axes[0].set_ylim(ylim)
    axes[0].axis("off")

    # --- 19 Districts Plot ---
    regions_19.plot(ax=axes[1], edgecolor="red", facecolor="#f08080", alpha=0.5)
    sites.plot(ax=axes[1], color="black", markersize=10)
    axes[1].set_xlim(xlim)
    axes[1].set_ylim(ylim)
    axes[1].axis("off")

    # Custom legend
    legend_elements = [
        Line2D([0], [0], color="blue", lw=2, label="11 Regions"),
        Line2D([0], [0], color="red", lw=2, label="19 Districts"),
        Line2D([0], [0], marker="o", color="black", markersize=6, linestyle="None", label="Sites")
    ]
    axes[1].legend(handles=legend_elements, loc="lower right", fontsize="medium")

    plt.tight_layout()
    save_figure(fig, "site_maps_side_by_side")


def plot_overlay_map():
    xlim, ylim = get_site_bounds_with_padding(pad=0.09)
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot both region types
    regions_11.plot(ax=ax, edgecolor="blue", facecolor="white", alpha=0.5)
    regions_19.plot(ax=ax, edgecolor="red", facecolor="white", alpha=0.5)
    sites.plot(ax=ax, color="black", markersize=10)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.axis("off")

    # Custom legend
    legend_elements = [
        Line2D([0], [0], color="blue", lw=2, label="11 Regions"),
        Line2D([0], [0], color="red", lw=2, label="19 Districts"),
        Line2D([0], [0], marker="o", color="black", markersize=6, linestyle="None", label="Sites")
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize="medium")

    plt.tight_layout()
    save_figure(fig, "site_overlay_map")


# --- 9. Run ---
if __name__ == "__main__":
    plot_side_by_side()
    plot_overlay_map()

    print("\n📘 FULL 11-region dictionary (id ➝ MXLOCATION list):")
    for region, site_list in region_11_dict.items():
        print(f"{region} ({region_11_counts[region]} sites): {site_list}")

    print("\n📕 FULL 19-district dictionary (GADM_ID ➝ MXLOCATION list):")
    for region, site_list in region_19_dict.items():
        print(f"{region} ({region_19_counts[region]} sites): {site_list}")
