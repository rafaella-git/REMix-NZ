# %% [markdown]
# import dependencies
#from will_build_instance import *
from build_instance_2020 import *
from remix.framework.tools.gdx import GDXEval
from remix.framework.tools.plots import plot_network, plot_choropleth
import geopandas as gpd
import warnings


# ## Define cases
# Global path variables
path_base = "C:/Local/REMix"
path_input, path_demand, path_profiles, path_geo,  path_brownfield, path_output, yrs_sel, yrs_to_calc, group_name, files_lst, data_dir, results_dir, geofile = [None] * 13
def set_up(name):
    global group_name, files_lst, yrs_sel, yrs_to_calc, path_input, path_demand, path_profiles, path_geo, path_brownfield, path_output


    will_elec = ["01-battery-distributed", "02-battery-overnight", "03-battery-recharging", "04-battery-solar"]
    will_h2 = ["01-h2-distributed", "02-h2-overnight", "03-h2-recharging", "04-h2-solar"]
    sdewes_ap = ["base", "high", "low", "med", "med_ev"]
        
    scenario_dict = {       
        "will": [will_h2, [2020, 2035, 2050]],
        "sdewes-ap": [sdewes_ap, [2020, 2030, 2040, 2050]]
    }
    group_name="sdewes-ap"

    # Check if the name exists as a key in the scenario_dict
    if name not in scenario_dict:
        print("Error: You need to create the case dictionary")
        return None

    group_name = name
    files_lst = scenario_dict[name][0]
    yrs_sel = scenario_dict[group_name][1] 
    yrs_to_calc = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
    # Define paths as global variables
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

    return group_name, files_lst, yrs_to_calc , yrs_sel
#%%
# Customise: set up the scenarios to run
set_up("sdewes-ap")

# %%
# Build data
def build_data(demand_file):
    nodes_lst=["NIS","AKL","WTO","TRN","BOP","HBY","CEN","WEL","NEL","CAN","OTG" ]
    yrs_str='-'.join([str(item) for item in yrs_sel])
    case_name=f"{demand_file}_{yrs_str}"
    geofile="11regionsNZ.geojson"

    # output
    data_dir = Path(f"{path_output}/{case_name}/data")
    data_dir.mkdir(parents=True, exist_ok=True)
    results_dir = Path(f"{path_output}/{case_name}/result")
    results_dir.mkdir(parents=True, exist_ok=True)
    print(f"----------Creating data for {case_name}")

    if __name__ == "__main__":
        # Create instance
        s1 = time.perf_counter()
        m = Instance(datadir=data_dir)

        add_scope(m)
        add_demand(m)

        # renewables
        add_renewables(m)
        add_geothermal(m)
        add_hydro(m)

        # batteries
        add_lithium_batteries(m)

        # conventional
        add_thermal(m)
        add_gas_turbines(m)

        # hydrogen
        if "h2" in demand_file:
            add_electrolyser(m)
            add_h2_storage(m)

        #add_methanizer(m)
        #add_methanol_syn(m)
        #add_ftropsch_syn(m)

        #carbon capture
        #add_dac(m)

        #others
        add_network(m)
        add_accounting(m)
        validate_scope(m)

        e1 = time.perf_counter()
        d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
        print(f"Dataframe management took {d1}.")

        # Create data
        s2 = int(time.perf_counter())
        m.write(output_path=f"{path_output}/{case_name}/data", fileformat="csv")
        e2 = time.perf_counter()
        d2=time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
        print(f"Writing dataset took {d2}.")

s0 = int(time.perf_counter())
for file in files_lst:
    s1 = int(time.perf_counter())
    build_data(file)
    e1 = time.perf_counter()
    d1 = time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
    print(f"Driting dataset took {d1}.")

# %% [markdown]
# Run the model
def case_run(demand_file):
    yrs_str='-'.join([str(item) for item in yrs_sel])
    case_name=f"{demand_file}_{yrs_str}"
    output_dir = Path(f"{path_output}/{case_name}")
    data_dir = Path(f"{output_dir}/data")
    result_dir = Path(f"{output_dir}/result")
    if not data_dir.exists():
        raise IOError("You need to build the data!")

    s1 = time.perf_counter()
    # running GAMS from Python script
    m = Instance(datadir=data_dir)
    m.run(
        resultdir=results_dir,
        resultfile=case_name,
        threads=12,
        lo=4,
        timeres=1,
        names=1,
        roundts=1,
        iis=1,
        profile=1,
        gdx="default",
        pathopt="myopic"
    )
    print(os. getcwd())

    e1 = time.perf_counter()
    d1=time.strftime("%Hh %Mm %Ss", time.gmtime(e1-s1))
    print(f"------------- Running {demand_file} took {d1}.")
for file in files_lst:
    s2 = int(time.perf_counter())
    case_run(file)
    e2 = time.perf_counter()
    d2 = time.strftime("%Hh %Mm %Ss", time.gmtime(e2-s2))
    print(f"Running the model took {d2}.")


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
for file in files_lst:
    s3 = int(time.perf_counter())
    results_to_csv(file)
    e3 = time.perf_counter()
    d3 = time.strftime("%Hh %Mm %Ss", time.gmtime(e3-s3))
    print(f"Writing dataset took {d1}.")
    print(f"Running the model took {d2}.")
    print(f"Exporting to csv files took {d3}.")
e4 = time.perf_counter()
d4 = time.strftime("%Hh %Mm %Ss", time.gmtime(e4-s0))
print(f"----------- Total time for {group_name} runs: {d4}.")
# %%
