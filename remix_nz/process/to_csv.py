# %% [markdown]
# import dependencies
from pathlib import Path
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import geopandas as gpd
import warnings

# ## Define cases
# Global path variables

will_elec = ["01-battery-distributed", "02-battery-overnight", "03-battery-recharging", "04-battery-solar"]
will_h2 = ["01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-h2-solar"]
sdewes_ap = ["base", "high", "low", "med", "ev-med"]
mbie=["base","h2"]

scenario_dict = {       
    "will": [will_elec, [2020, 2035, 2050]],
    "sdewes-ap": [sdewes_ap, [2020, 2030, 2040, 2050]],
    "mbie": [mbie, [2020, 2030, 2040, 2050]]
}
group_name="mbie"
indx=0

files_lst = scenario_dict[group_name][0]
yrs_sel = scenario_dict[group_name][1] 
yrs_str='-'.join([str(item) for item in yrs_sel])
yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
# Define paths as global variables
path_base = "C:/Local/REMix"
path_input = f"{path_base}/remix_nz/input"
path_demand = f"{path_input}/demand/{group_name}"
path_profiles = f"{path_input}/profiles"      # renewables
path_brownfield = f"{path_input}/brownfield"  # info hydro and existing power plants database
path_geo = f"{path_input}/shapefiles"         # geojson
geofile="11regionsNZ.geojson"
demand_file=files_lst[indx] 
case_name=f"{demand_file}_{yrs_str}"
path_output = f"{path_base}/remix_nz/output/{group_name}"
data_dir = Path(f"{path_output}/{case_name}/data")
data_dir.mkdir(parents=True, exist_ok=True)
results_dir = Path(f"{path_output}/{case_name}/result")
results_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# Export results csv
def results_to_csv(demand_file):
    case_name=f"{demand_file}_{yrs_str}"
    output_dir = Path(f"{path_output}/{case_name}")
    data_dir = Path(f"{output_dir}/data")
    result_dir = Path(f"{output_dir}/result")
    
    # read in the output `*.gdx` file from the optimization in GAMS
    results = GDXEval(f"{result_dir}/{case_name}.gdx")
    
    # Converter capacities
    converter_caps = results["converter_caps"]   # convert converter capacities to a Pandas DataFrame
    converter_caps = converter_caps[converter_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    converter_caps.to_csv(f"{result_dir}/{case_name}_converter_caps.csv")
    # Example of slicing:
    # accNodesModel=all, accYears=2050 only, techs=all, commodit=only Elec, capType=total
    # caps.loc[idx[:, "2050", :, "Elec", "total"], :].round(2) 

    # Storage capacities
    storage_caps = results["storage_caps"]   # convert storage capacities to a Pandas DataFrame
    storage_caps = storage_caps[storage_caps > 0.01].dropna()  # remove all capacities with less than 10 MW
    storage_caps.to_csv(f"{result_dir}/{case_name}_storage_caps.csv")

    # Commodity balance annual
    commodity_balance_annual = results["commodity_balance_annual"]   # convert commodity_balance to a Pandas DataFrame
    commodity_balance_annual.to_csv(f"{result_dir}/{case_name}_commodity_balance_annual.csv")

    # Marginals
    marginals_sourcesink_profile = results["marginals_sourcesink_profile"]   # convert marginals to a Pandas DataFrame
    marginals_sourcesink_profile.to_csv(f"{result_dir}/{case_name}_marginals_sourcesink_profile.csv")

    # Accounting indicators
    indicator_accounting = results["indicator_accounting"]   # convert  to a Pandas DataFrame
    indicator_accounting.to_csv(f"{result_dir}/{case_name}_indicator_accounting.csv")


results_to_csv(demand_file)